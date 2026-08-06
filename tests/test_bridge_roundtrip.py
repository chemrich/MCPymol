"""End-to-end tests over a real TCP socket.

Every other test in the suite mocks either the socket or ``send_request``,
which means the wire format itself — half-close framing, chunked reads, the
incremental JSON parse on both sides — is only ever checked against a mock
that was written to match the implementation.  These tests run the real
bridge against the real plugin listener on an ephemeral port, so a framing
regression fails here even when the mocks still agree with each other.
"""

import json
import socket
import threading
import time
from unittest.mock import patch

import pytest

import mcpymol.bridge as bridge_module
from mcpymol.bridge import send_request
from mcpymol.plugin import PyMOLSocketServer


@pytest.fixture
def live_bridge():
    """Start a real listener on an ephemeral port and point send_request at it.

    Yields the running PyMOLSocketServer.
    """
    server = PyMOLSocketServer(port=0)
    server.start()

    # _listen_loop binds on its own thread; wait for a real port to appear.
    deadline = time.monotonic() + 5.0
    port = 0
    while time.monotonic() < deadline:
        sock = server.server_socket
        if sock is not None:
            try:
                port = sock.getsockname()[1]
            except OSError:
                port = 0
            if port:
                break
        time.sleep(0.01)
    assert port, "listener never bound"

    # send_request reads PORT out of its own module globals, so that is the
    # binding that has to move — not the re-export on mcpymol.server.
    with patch.object(bridge_module, "PORT", port):
        yield server

    server.stop()


def test_round_trip_success(live_bridge):
    """A normal command crosses the wire and comes back decoded."""
    with patch.object(live_bridge, "handle_request") as handler:
        handler.return_value = {"status": "success", "result": ["A", "B"]}
        res = send_request("get_chains", args=["all"])

    assert res == {"status": "success", "result": ["A", "B"]}
    payload = handler.call_args[0][0]
    assert payload["action"] == "get_chains"
    assert payload["args"] == ["all"]
    assert payload["kwargs"] == {}


def test_round_trip_large_response(live_bridge):
    """Regression for the 8 KB truncation that corrupted get_fastastr.

    A response several chunks long must arrive intact, not clipped at the
    first recv().
    """
    big = "M" * 500_000
    with patch.object(live_bridge, "handle_request") as handler:
        handler.return_value = {"status": "success", "result": big}
        res = send_request("get_fastastr", args=["1abc and chain A"])

    assert res["status"] == "success"
    assert res["result"] == big


def test_round_trip_large_request(live_bridge):
    """conservation_view sends a batched script that can be tens of KB."""
    script = "\n".join(f"alter resi {i}, b={i}" for i in range(20_000))
    with patch.object(live_bridge, "handle_request") as handler:
        handler.return_value = {"status": "success", "result": "OK"}
        send_request("do", args=[script])

    assert handler.call_args[0][0]["args"] == [script]


def test_round_trip_error_response(live_bridge):
    with patch.object(live_bridge, "handle_request") as handler:
        handler.return_value = {"status": "error", "error": "no such object"}
        res = send_request("zoom", args=["nope"])

    assert res == {"status": "error", "error": "no such object"}


def test_round_trip_unicode(live_bridge):
    """Multi-byte characters must survive being split across recv() chunks."""
    text = "α-helix → β-sheet ± 2.5 Å " * 5000
    with patch.object(live_bridge, "handle_request") as handler:
        handler.return_value = {"status": "success", "result": text}
        res = send_request("help", args=["cartoon"])

    assert res["result"] == text


def test_oversized_request_is_rejected_not_absorbed(live_bridge):
    """The size cap must produce a real error response, not a dead connection."""
    with (
        patch("mcpymol.plugin.MAX_REQUEST_BYTES", 4096),
        patch.object(live_bridge, "handle_request") as handler,
    ):
        res = send_request("do", args=["x" * 100_000])

    assert res["status"] == "error"
    assert "Request rejected" in res["error"]
    handler.assert_not_called()


def test_stalled_client_does_not_wedge_the_listener(live_bridge):
    """The core regression: a peer that connects and never finishes its request
    must not block the next command for the rest of the session."""
    # Must be patched before the stalled peer is accepted — serve_connection
    # reads RECV_TIMEOUT once, at the moment it takes the connection.
    with patch("mcpymol.plugin.RECV_TIMEOUT", 0.5):
        stalled = socket.create_connection(("127.0.0.1", bridge_module.PORT))
        try:
            # A deliberately incomplete JSON document, and never half-closed.
            stalled.sendall(b'{"action": "refresh"')

            # Give the listener a moment to pick the connection up and block on it.
            time.sleep(0.2)

            with patch.object(live_bridge, "handle_request") as handler:
                handler.return_value = {"status": "success", "result": "OK"}
                res = send_request("refresh", timeout=15.0)

            assert res == {"status": "success", "result": "OK"}
        finally:
            stalled.close()


def test_concurrent_clients_are_serialized_not_dropped(live_bridge):
    """pymol.cmd is not thread-safe, so requests are handled one at a time —
    but every one of them must still get its own correct answer."""
    seen: list[str] = []
    lock = threading.Lock()

    def handler(payload):
        with lock:
            seen.append(payload["args"][0])
        return {"status": "success", "result": payload["args"][0]}

    results: dict[int, dict] = {}

    def worker(i):
        results[i] = send_request("color", args=[f"sel{i}"], timeout=15.0)

    with patch.object(live_bridge, "handle_request", side_effect=handler):
        threads = [threading.Thread(target=worker, args=(i,)) for i in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=20)

    assert len(results) == 8
    for i, res in results.items():
        assert res["result"] == f"sel{i}", f"client {i} got another client's answer"
    assert sorted(seen) == sorted(f"sel{i}" for i in range(8))


def test_non_json_garbage_gets_an_error_response(live_bridge):
    """A client speaking the wrong protocol gets told so, and the listener
    stays up for the next real request."""
    with socket.create_connection(("127.0.0.1", bridge_module.PORT)) as raw:
        raw.sendall(b"GET / HTTP/1.1\r\nHost: localhost\r\n\r\n")
        raw.shutdown(socket.SHUT_WR)
        reply = b""
        while True:
            chunk = raw.recv(4096)
            if not chunk:
                break
            reply += chunk

    decoded = json.loads(reply)
    assert decoded["status"] == "error"
    assert "Invalid JSON payload" in decoded["error"]

    # Listener survived.
    with patch.object(live_bridge, "handle_request") as handler:
        handler.return_value = {"status": "success", "result": "OK"}
        assert send_request("refresh")["status"] == "success"
