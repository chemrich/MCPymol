import json
import socket

# We need to mock pymol.cmd before we can even import plugin.py
import sys
from unittest.mock import MagicMock

import pytest

mock_pymol = MagicMock()
sys.modules["pymol"] = mock_pymol
sys.modules["pymol.cmd"] = mock_pymol.cmd

# Now we can import the plugin
import mcpymol.plugin as plugin_module
from mcpymol.plugin import (
    MAX_REQUEST_BYTES,
    PyMOLSocketServer,
    RequestTooLarge,
    _recv_all,
)


@pytest.fixture
def plugin_server():
    """Initializes a PyMOLSocketServer instance ready for simulated payloads."""
    server = PyMOLSocketServer()
    # Reset the mock cmd calls before every test so they are isolated
    mock_pymol.cmd.reset_mock()
    return server


def test_handle_load_structure(plugin_server):
    """Test standard file loading mapped correctly to cmd.load."""
    payload = json.dumps({"action": "load", "args": ["1ubq.pdb"]})

    res = plugin_server.handle_request(payload)

    assert res["status"] == "success"
    assert "Loaded 1ubq.pdb" in res["result"]
    mock_pymol.cmd.load.assert_called_once_with("1ubq.pdb")


def test_handle_get_chains(plugin_server):
    """Test custom get_chains implementation."""
    payload = json.dumps({"action": "get_chains", "args": ["all"]})
    mock_pymol.cmd.get_chains.return_value = ["A", "B", "X"]

    res = plugin_server.handle_request(payload)

    assert res["status"] == "success"
    assert res["result"] == ["A", "B", "X"]
    mock_pymol.cmd.get_chains.assert_called_once_with("all")


def test_dynamic_method_resolution(plugin_server):
    """Test that arbitrary PyMOL cmd functions (e.g., cmd.hide) are dynamically invoked."""

    # We add a fake callable attribute to the mock for testing dynamic resolution
    def fake_hide(*args, **kwargs):
        pass

    mock_pymol.cmd.hide = fake_hide

    payload = json.dumps({"action": "hide", "args": ["everything", "all"]})
    res = plugin_server.handle_request(payload)

    assert res["status"] == "success"
    assert "Executed 'hide' successfully" in res["result"]


def test_unknown_action(plugin_server):
    """Test behavior when an invalid action is requested."""
    payload = json.dumps({"action": "nonexistent_pymol_command"})

    # We implicitly mock missing methods by explicitly deleting them on the MagicMock
    # to force `hasattr(cmd, action)` to return False
    del mock_pymol.cmd.nonexistent_pymol_command

    res = plugin_server.handle_request(payload)

    assert res["status"] == "error"
    assert "Unknown action or method not found on cmd" in res["error"]


def test_malformed_json(plugin_server):
    """Test handling of invalid JSON payloads from a buggy client."""
    payload = '{"action": "load", "args": ["1ubq.pdb"]'  # Missing closing brace

    res = plugin_server.handle_request(payload)

    assert res["status"] == "error"
    assert "Invalid JSON" in res["error"]


def test_pymol_internal_exception(plugin_server):
    """Test handling of PyMOL execution errors."""

    mock_pymol.cmd.color.side_effect = Exception("Invalid color name 'bleargh'")

    payload = json.dumps({"action": "color", "args": ["bleargh", "all"]})
    res = plugin_server.handle_request(payload)

    assert res["status"] == "error"
    assert "Invalid color name 'bleargh'" in res["error"]


@pytest.mark.parametrize("payload", ["[1, 2, 3]", '"just a string"', "42"])
def test_non_object_payload_rejected(plugin_server, payload):
    """A JSON document that isn't an object must not blow up on .get()."""
    res = plugin_server.handle_request(payload)

    assert res["status"] == "error"
    assert "must be a JSON object" in res["error"]


# ── Import must not bind a port ──────────────────────────────────────────────


def test_import_does_not_autostart():
    """Importing the module as a library must not open a listening socket.

    `run plugin.py` inside PyMOL still auto-starts (different __name__), but a
    plain import — the test suite, tooling — would otherwise race a real PyMOL
    session already holding the port.
    """
    assert plugin_module.__name__ == "mcpymol.plugin"
    assert plugin_module.mcp_bridge_plugin is None


# ── _recv_all framing ────────────────────────────────────────────────────────


def _fake_sock(*chunks):
    """A socket whose recv() yields the given chunks, then EOF forever."""
    sock = MagicMock()
    sock.recv.side_effect = list(chunks) + [b""] * 10
    return sock


def test_recv_all_stops_at_complete_json_without_half_close():
    """A well-formed request is served even if the peer never half-closes."""
    body = json.dumps({"action": "refresh"}).encode()
    sock = _fake_sock(body)

    data = _recv_all(sock)

    assert json.loads(data) == {"action": "refresh"}
    # Stopped on the complete document rather than reading on to EOF.
    assert sock.recv.call_count == 1


def test_recv_all_reassembles_split_json():
    body = json.dumps({"action": "color", "args": ["red", "all"]}).encode()
    sock = _fake_sock(body[:10], body[10:])

    assert json.loads(_recv_all(sock)) == {"action": "color", "args": ["red", "all"]}


def test_recv_all_returns_empty_on_immediate_eof():
    assert _recv_all(_fake_sock()) == b""


def test_recv_all_rejects_oversized_request():
    """A peer streaming junk must be cut off, not allowed to exhaust memory."""
    # Never parses as JSON, so only the size cap can stop it.
    sock = MagicMock()
    sock.recv.side_effect = [b"x" * 4096] * 10_000

    with pytest.raises(RequestTooLarge):
        _recv_all(sock, max_bytes=64 * 1024)


def test_recv_all_default_cap_is_the_module_constant():
    sock = MagicMock()
    sock.recv.side_effect = [b"x" * (MAX_REQUEST_BYTES + 1)]

    with pytest.raises(RequestTooLarge):
        _recv_all(sock)


def test_recv_all_propagates_timeout():
    sock = MagicMock()
    sock.recv.side_effect = TimeoutError("timed out")

    with pytest.raises(TimeoutError):
        _recv_all(sock)


# ── serve_connection: every failure path answers ─────────────────────────────


def _served(plugin_server, sock):
    """Run serve_connection and return the decoded response (or None)."""
    plugin_server.serve_connection(sock)
    if not sock.sendall.call_args_list:
        return None
    return json.loads(sock.sendall.call_args[0][0])


def test_serve_connection_happy_path(plugin_server):
    mock_pymol.cmd.get_chains.return_value = ["A"]
    sock = _fake_sock(json.dumps({"action": "get_chains", "args": ["all"]}).encode())

    res = _served(plugin_server, sock)

    assert res == {"status": "success", "result": ["A"]}
    sock.settimeout.assert_called_once_with(plugin_module.RECV_TIMEOUT)


def test_serve_connection_answers_on_timeout(plugin_server):
    """A stalled peer gets an error back instead of wedging the accept loop."""
    sock = MagicMock()
    sock.recv.side_effect = TimeoutError("timed out")

    res = _served(plugin_server, sock)

    assert res["status"] == "error"
    assert "Timed out" in res["error"]


def test_serve_connection_answers_on_oversized_request(plugin_server):
    sock = MagicMock()
    sock.recv.side_effect = [b"x" * 4096] * 10_000

    res = _served(plugin_server, sock)

    assert res["status"] == "error"
    assert "Request rejected" in res["error"]


def test_serve_connection_answers_on_bad_json(plugin_server):
    sock = _fake_sock(b'{"action": "load"')  # truncated, then EOF

    res = _served(plugin_server, sock)

    assert res["status"] == "error"
    assert "Invalid JSON payload" in res["error"]


def test_serve_connection_silent_on_empty_request(plugin_server):
    """A peer that connects and closes without sending gets no reply."""
    assert _served(plugin_server, _fake_sock()) is None


def test_serve_connection_survives_dead_peer(plugin_server):
    """sendall failing on a vanished peer must not raise into the accept loop."""
    mock_pymol.cmd.refresh.return_value = None
    sock = _fake_sock(json.dumps({"action": "refresh"}).encode())
    sock.sendall.side_effect = BrokenPipeError("peer gone")

    plugin_server.serve_connection(sock)  # must not raise


def test_serve_connection_survives_unsettable_timeout(plugin_server):
    """settimeout can fail on an already-closed socket; keep going."""
    mock_pymol.cmd.refresh.return_value = None
    sock = _fake_sock(json.dumps({"action": "refresh"}).encode())
    sock.settimeout.side_effect = OSError("bad fd")

    res = _served(plugin_server, sock)

    assert res["status"] == "success"


def test_serve_connection_read_error_is_swallowed(plugin_server):
    """A hard read failure leaves nobody to answer; don't raise."""
    sock = MagicMock()
    sock.recv.side_effect = ConnectionResetError("peer reset")

    assert _served(plugin_server, sock) is None


def test_serve_connection_handles_unserializable_result(plugin_server):
    """A cmd function returning an opaque object must still produce a reply."""

    class Opaque:
        def __repr__(self):
            return "<chempy model>"

    mock_pymol.cmd.get_model.return_value = Opaque()
    sock = _fake_sock(json.dumps({"action": "get_model", "args": ["all"]}).encode())

    res = _served(plugin_server, sock)

    assert res["status"] == "success"
    assert res["result"] == "<chempy model>"


def test_timeout_error_is_an_oserror_subclass():
    """Guards the except-ordering in serve_connection: TimeoutError must be
    caught before the generic OSError arm or stalled peers get no reply."""
    assert issubclass(TimeoutError, OSError)
    assert socket.timeout is TimeoutError
