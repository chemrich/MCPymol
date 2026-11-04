import json
import os
import queue
import re
import socket
import sys
import threading
import time
from types import SimpleNamespace
from unittest.mock import MagicMock, call, patch

import pytest

# conftest installs the pymol mock before any test module is imported; reuse
# that exact object — a fresh MagicMock here would not be the one the plugin
# bound at import time, so every assertion on it would silently pass on an
# object nothing calls.
mock_pymol = sys.modules["pymol"]

import mcpymol.bridge as bridge_module
import mcpymol.plugin as plugin_module
from mcpymol.bridge import send_request
from mcpymol.plugin import (
    MAX_REQUEST_BYTES,
    MainThreadTimeout,
    PyMOLSocketServer,
    RequestTooLarge,
    _AppNapGuard,
    _do_prefix,
    _OutputCapture,
    _problem_lines,
    _recv_all,
    run_do,
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


# ── iterate_to_list ─────────────────────────────────────────────────────────


def test_iterate_to_list_reads_per_atom_properties(plugin_server):
    """The only route to atom-level properties over the wire.

    cmd.get_model returns a chempy object that does not survive JSON, so
    occupancy and altloc are otherwise unreachable by a client.
    """

    def fake_iterate(selection, expression, space=None):
        # Mimic PyMOL: evaluate the expression once per atom, in `space`.
        for chain, resi, q in (("A", "1", 1.0), ("A", "2", 0.6)):
            eval(expression, {"chain": chain, "resi": resi, "q": q, **(space or {})})

    mock_pymol.cmd.iterate = fake_iterate

    payload = json.dumps({"action": "iterate_to_list", "args": ["polymer", "chain, resi, q"]})
    res = plugin_server.handle_request(payload)

    assert res["status"] == "success"
    assert res["result"] == [("A", "1", 1.0), ("A", "2", 0.6)]


def test_iterate_to_list_single_field_still_yields_tuples(plugin_server):
    """Arity follows the expression, so callers never special-case one field."""

    def fake_iterate(selection, expression, space=None):
        for q in (1.0, 0.5):
            eval(expression, {"q": q, **(space or {})})

    mock_pymol.cmd.iterate = fake_iterate

    res = plugin_server.handle_request(
        json.dumps({"action": "iterate_to_list", "args": ["all", "q"]})
    )

    assert res["result"] == [(1.0,), (0.5,)]


def test_iterate_to_list_accepts_kwargs(plugin_server):
    def fake_iterate(selection, expression, space=None):
        assert selection == "chain B"
        eval(expression, {"resi": "7", **(space or {})})

    mock_pymol.cmd.iterate = fake_iterate

    res = plugin_server.handle_request(
        json.dumps(
            {
                "action": "iterate_to_list",
                "kwargs": {"selection": "chain B", "expression": "resi"},
            }
        )
    )

    assert res["status"] == "success"
    assert res["result"] == [("7",)]


def test_iterate_to_list_without_an_expression_is_an_error(plugin_server):
    res = plugin_server.handle_request(json.dumps({"action": "iterate_to_list"}))

    assert res["status"] == "error"
    assert "needs an expression" in res["error"]


def test_iterate_to_list_empty_selection_returns_an_empty_list(plugin_server):
    """No atoms is a legitimate answer, not an error — the caller decides."""

    def fake_iterate(selection, expression, space=None):
        return None  # PyMOL calls the expression zero times

    mock_pymol.cmd.iterate = fake_iterate

    res = plugin_server.handle_request(
        json.dumps({"action": "iterate_to_list", "args": ["none", "q"]})
    )

    assert res["status"] == "success"
    assert res["result"] == []


def test_iterate_to_list_propagates_pymol_errors(plugin_server):
    def boom(*a, **k):
        raise RuntimeError("Invalid selection")

    mock_pymol.cmd.iterate = boom

    res = plugin_server.handle_request(
        json.dumps({"action": "iterate_to_list", "args": ["bogus(", "q"]})
    )

    assert res["status"] == "error"
    assert "Invalid selection" in res["error"]


# ── numpy coercion in the response path ─────────────────────────────────────


class _FakeArray:
    """Stands in for a numpy array: not JSON-able, but has .tolist()."""

    def __init__(self, data):
        self._data = data

    def tolist(self):
        return self._data

    def __repr__(self):
        return "array(...)"


def test_numpy_result_is_converted_not_reprd(plugin_server):
    """cmd.get_coords returns a numpy array.

    Before coercion this hit the repr fallback and reached the client as
    ``status: success`` carrying the *string* "array(...)" — a reply that
    looks like data and is not. Regression guard for that.
    """
    mock_pymol.cmd.get_coords.return_value = _FakeArray([[1.0, 2.0, 3.0]])
    sock = _fake_sock(json.dumps({"action": "get_coords", "args": ["all"]}).encode())

    res = _served(plugin_server, sock)

    assert res["status"] == "success"
    assert res["result"] == [[1.0, 2.0, 3.0]]
    assert res["result"] != "array(...)"


def test_nested_numpy_is_converted(plugin_server):
    mock_pymol.cmd.get_extent.return_value = {"box": _FakeArray([1.0, 2.0])}
    sock = _fake_sock(json.dumps({"action": "get_extent", "args": ["all"]}).encode())

    res = _served(plugin_server, sock)

    assert res["result"] == {"box": [1.0, 2.0]}


def test_set_result_is_converted_to_a_list(plugin_server):
    mock_pymol.cmd.get_names.return_value = {"a", "b"}
    sock = _fake_sock(json.dumps({"action": "get_names", "args": []}).encode())

    res = _served(plugin_server, sock)

    assert sorted(res["result"]) == ["a", "b"]


def test_opaque_object_still_falls_back_to_repr(plugin_server):
    """Coercion must not remove the existing safety net."""

    class Opaque:
        def __repr__(self):
            return "<chempy model>"

    mock_pymol.cmd.get_model.return_value = Opaque()
    sock = _fake_sock(json.dumps({"action": "get_model", "args": ["all"]}).encode())

    res = _served(plugin_server, sock)

    assert res["status"] == "success"
    assert res["result"] == "<chempy model>"


def test_tolist_that_raises_falls_back_to_repr(plugin_server):
    """A broken .tolist() must not take down the reply."""

    class Hostile:
        def tolist(self):
            raise ValueError("nope")

        def __repr__(self):
            return "<hostile>"

    mock_pymol.cmd.get_coords.return_value = Hostile()
    sock = _fake_sock(json.dumps({"action": "get_coords", "args": ["all"]}).encode())

    res = _served(plugin_server, sock)

    assert res["status"] == "success"
    assert res["result"] == "<hostile>"


# ── `do`: PyMOL's verdict on every line, after the lines have run ───────────
#
# The old handler queued the text with cmd.do and replied "success" at once:
# from the plugin's thread cmd.do only queues a line for PyMOL's main thread,
# so the reply came before the command ran and could never carry its error.


class FakeParser:
    """Stands in for pymol.parser.Parser.

    PyMOL's parser prints its errors rather than raising, and returns 0 for a
    line that failed and 1 for one that ran.  Each line here gets a scripted
    verdict and prints what PyMOL would have.  ``states`` sets the parser's top
    layer after a line: an open python block (``embed_sentinel``) or a
    trailing-backslash continuation (``cont``).
    """

    def __init__(self, verdicts, prints=None, raises=None, states=None):
        self.verdicts = list(verdicts)
        self.prints = dict(prints or {})
        self.raises = dict(raises or {})
        self.states = dict(states or {})
        self.layer = {1: SimpleNamespace(embed_sentinel=None, cont="")}
        self.lines: list[str] = []

    def parse(self, line, secure=0):
        self.lines.append(line)
        n = len(self.lines)
        if n in self.states:
            self.layer[1].embed_sentinel, self.layer[1].cont = self.states[n]
        if n in self.prints:
            print(self.prints[n])
        if n in self.raises:
            raise self.raises[n]
        return self.verdicts[n - 1]


@pytest.fixture
def fake_parser(monkeypatch):
    """Make run_do parse with a FakeParser.  Tests run on the main thread, so
    run_do runs the lines inline, as it does when PyMOL's main thread is the
    caller."""
    mock_pymol.cmd.reset_mock()

    def install(*args, **kwargs):
        parser = FakeParser(*args, **kwargs)
        monkeypatch.setattr(plugin_module, "_new_parser", lambda: parser)
        return parser

    return install


def test_do_reports_the_rejected_line_and_runs_the_rest(fake_parser):
    """The bug: a bad line in a `do` read as "Executed command" because the
    reply never looked at what PyMOL made of it."""
    parser = fake_parser([1, 0, 1], prints={2: ' Error: Invalid selection name "x".\nx<--'})

    res = run_do("set ray_shadow, 0\ncolor red, x\nshow sticks")

    assert res["status"] == "error"
    assert "PyMOL rejected line 2 of 3: color red, x" in res["error"]
    assert 'Invalid selection name "x"' in res["error"]
    assert "The other 2 line(s) ran." in res["error"]
    assert res["failures"][0]["line"] == 2
    assert 'Invalid selection name "x"' in res["output"]
    assert "result" not in res
    # Like cmd.do, the lines after a failure still run.
    assert parser.lines == ["set ray_shadow, 0", "color red, x", "show sticks"]


def test_do_handler_replies_with_pymol_verdict(plugin_server, fake_parser):
    """The handler itself: it called cmd.do and returned success, always."""
    fake_parser([0], prints={1: ' Error: Invalid selection name "nosuch".\nnosuch<--'})

    res = plugin_server.handle_request({"action": "do", "args": ["color red, nosuch"]})

    assert res["status"] == "error"
    assert res["error"].startswith("PyMOL rejected the command: color red, nosuch\n")
    assert 'Invalid selection name "nosuch"' in res["error"]
    assert mock_pymol.cmd.do.call_count == 0  # parsed here, not queued blind
    assert not plugin_server._in_do


def test_do_single_line_rejection_names_the_command(fake_parser):
    fake_parser([0], prints={1: " Error: unknown Setting: 'cartoon_tube_radiusx'"})

    res = run_do("set cartoon_tube_radiusx, 0.5")

    assert res["error"].startswith("PyMOL rejected the command: set cartoon_tube_radiusx, 0.5\n")
    assert "unknown Setting" in res["error"]
    assert "line(s) ran" not in res["error"]


def test_do_success_carries_output_and_warnings(fake_parser):
    """What PyMOL printed comes back, and warning lines from lines that still
    ran are listed: a level-mismatched set "succeeds" with only a warning."""
    fake_parser(
        [1, 1],
        prints={
            1: ' Setting: cartoon_tube_radius set for 43 atoms in object "ub".\n'
            " Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting",
            2: "### cut below here and paste into script ###\nset_view (\\",
        },
    )

    res = run_do("set cartoon_tube_radius, 0.5, ub and resi 1-5\nget_view")

    assert res["status"] == "success"
    assert (
        res["result"] == "Executed command: set cartoon_tube_radius, 0.5, ub and resi 1-5\nget_view"
    )
    assert res["warnings"] == [
        "Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting"
    ]
    assert "set_view" in res["output"]


def test_do_reply_always_has_output(fake_parser):
    """The bridge reads a missing `output` as a plugin too old to report
    errors at all, so a quiet success must still carry an empty one."""
    fake_parser([1])

    res = run_do("hide everything")

    assert res == {"status": "success", "result": "Executed command: hide everything", "output": ""}


def test_empty_do_does_not_wait_on_the_main_thread(monkeypatch):
    def must_not_parse():
        raise AssertionError("an empty request went to the parser")

    monkeypatch.setattr(plugin_module, "_new_parser", must_not_parse)
    monkeypatch.setattr(plugin_module, "_on_main_thread", must_not_parse)

    assert run_do("") == {"status": "success", "result": "Executed command: ", "output": ""}


def test_parser_crash_reports_the_error_it_was_printing(fake_parser, monkeypatch):
    """`alter all, name='CA'; rebuild`: alter takes the rest of the line
    literally, so 'rebuild' is evaluated per atom and raises NameError.  Under
    Python 3.14 PyMOL's error printer then fails with a TypeError of its own,
    which escapes parse(); the NameError is at the root of its context."""
    monkeypatch.setattr(mock_pymol, "parsing", SimpleNamespace(LITERAL=20), raising=False)
    monkeypatch.setattr(mock_pymol.cmd, "keyword", {"alter": [None, 0, 0, "", 21]}, raising=False)
    printer_failure = TypeError("'wrapper.Wrapper' object is not iterable")
    printer_failure.__context__ = NameError("name 'rebuild' is not defined")
    fake_parser([0], raises={1: printer_failure})

    res = run_do("alter all, name='CA'; rebuild")

    assert res["status"] == "error"
    assert "NameError: name 'rebuild' is not defined" in res["error"]
    assert "alter takes the rest of the line literally" in res["error"]
    assert "Wrapper" not in res["error"]


@pytest.mark.parametrize(
    "raw, expected",
    [
        # "_ x" runs x without an echo, but still logs it.
        ("_ set ray_shadow, 1", ("set ray_shadow, 1", None, True)),
        # A pasted prompt is stripped before parsing; the raw line is echoed.
        ("PyMOL> zoom", ("zoom", "PyMOL> zoom", True)),
        ("PyMOL>zoom", ("zoom", "PyMOL>zoom", True)),
        # Internal call-backs: neither echoed nor logged.
        ("_cmd.refresh()", ("_cmd.refresh()", None, False)),
        ("cmd._mcpymol_run_job('x')", ("cmd._mcpymol_run_job('x')", None, False)),
        ("color red, all", ("color red, all", "color red, all", True)),
    ],
)
def test_do_prefix_mirrors_c_cmddo(raw, expected):
    assert _do_prefix(raw) == expected


def test_do_echoes_and_logs_each_line_as_cmd_do_did(fake_parser, capsys):
    """The console and a .pml log must read as they did when cmd.do ran the
    text: echo per line, `_ ` lines logged but not echoed, internals neither."""
    parser = fake_parser([1, 1, 1])

    run_do("_ set ray_shadow, 1\nPyMOL>set ray_shadow, 0\ncmd._mcpymol_noop()")

    out = capsys.readouterr().out
    assert "PyMOL>PyMOL>set ray_shadow, 0" in out
    assert "PyMOL>_ set" not in out and "PyMOL>set ray_shadow, 1" not in out
    assert "cmd._mcpymol_noop" not in out
    assert mock_pymol.cmd.log.call_args_list == [
        call("set ray_shadow, 1\n"),
        call("set ray_shadow, 0\n"),
    ]
    assert parser.lines == ["set ray_shadow, 1", "set ray_shadow, 0", "cmd._mcpymol_noop()"]


def test_multi_line_do_defers_updates_and_restores_them(fake_parser, monkeypatch):
    fake_parser([1, 1])
    monkeypatch.setattr(mock_pymol.cmd, "get_setting_int", MagicMock(return_value=0))

    run_do("show sticks\nhide lines")

    assert mock_pymol.cmd.set.call_args_list == [
        call("defer_updates", 1),
        call("defer_updates", 0),
    ]


def test_unclosed_python_block_is_a_failure_not_a_trap(fake_parser):
    """Each request gets its own parser, so an unclosed block cannot swallow
    the user's next command; it is reported instead of silently doing nothing."""
    fake_parser([1, 1], states={1: ("python end", "")})

    res = run_do("python\nprint('never runs')")

    assert res["status"] == "error"
    assert "PyMOL rejected line 1 of 2: python" in res["error"]
    assert "never closed" in res["error"] and "'python end'" in res["error"]
    assert "No line ran." in res["error"]


def test_trailing_continuation_is_a_failure(fake_parser):
    fake_parser([1, 1], states={2: (None, "color red, ")})

    res = run_do("set ray_shadow, 0\ncolor red, \\")

    assert "PyMOL rejected line 2 of 2" in res["error"]
    assert "continuation" in res["error"]
    assert "The other 1 line(s) ran." in res["error"]


def test_do_output_keeps_head_and_tail(fake_parser, monkeypatch):
    monkeypatch.setattr(plugin_module, "MAX_OUTPUT_CHARS", 400)
    fake_parser([1], prints={1: "\n".join(f"atom {i}" for i in range(5000))})

    res = run_do("iterate all, print(name)")

    assert len(res["output"]) < 500
    assert res["output"].startswith("atom 0\n")
    assert res["output"].rstrip().endswith("atom 4999")
    assert "lines omitted" in res["output"]


def test_runaway_warnings_and_rejections_are_capped(fake_parser):
    """A command can warn once per atom, and a long script can fail on every
    line; the reply lists twenty of each and counts the rest."""
    warns = "\n".join(f' Selector-Error: Invalid selection name "s{i}".' for i in range(100))
    fake_parser([1], prints={1: warns + "\n" + warns})  # repeats are listed once
    res = run_do("zoom many")
    assert len(res["warnings"]) == 21
    assert res["warnings"][-1] == "... and 80 more"

    fake_parser([0] * 30)
    res = run_do("\n".join(f"color red, bad{i}" for i in range(30)))
    assert len(res["failures"]) == 20
    assert res["error"].count("color red, bad") == 21  # the first, then 20 more
    assert "; and 9 more lines" in res["error"]


def test_a_huge_rejected_line_is_quoted_short(fake_parser):
    fake_parser([0], prints={1: "SyntaxError: invalid syntax"})

    res = run_do("x" * 5000)

    assert len(res["error"]) < 500
    assert "[4700 more characters]" in res["error"]


def test_problem_lines_picks_pymol_warnings_and_errors():
    text = (
        " Setting: ray_shadow set to off.\n"
        " Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting\n"
        ' Selector-Error: Invalid selection name "nosuch".\n'
        "nosuch<--\n"
        "Setting-Error: unknown color 'nosuchcolor'\n"
        ' ScenePNG-Error: error writing "/nonexistent_dir/x.png"! Please check directory...\n'
        " Executive: Colored 660 atoms and 1 object.\n"
    )

    assert _problem_lines(text) == [
        "Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting",
        'Selector-Error: Invalid selection name "nosuch". nosuch<--',
        "Setting-Error: unknown color 'nosuchcolor'",
        'ScenePNG-Error: error writing "/nonexistent_dir/x.png"! Please check directory...',
    ]


# ── _OutputCapture: C feedback on fd 1, Python text on sys.stdout ───────────


def test_capture_sees_c_and_python_output_and_replays_both(tmp_path, capsys):
    """C-level feedback (a level warning, a selector error) is written straight
    to fd 1 and never passes through Python; parser errors are print()ed.  The
    capture has to see both, and the terminal must still get both afterwards.
    Three captures in a row exercise the reused spool file."""
    sink = tmp_path / "fd1.txt"
    fd = os.open(sink, os.O_RDWR | os.O_CREAT | os.O_TRUNC)
    saved = os.dup(1)
    os.dup2(fd, 1)
    seen = []
    try:
        for i in range(3):
            with _OutputCapture() as captured:
                os.write(1, f" Setting-Warning: c{i}\n".encode())
                print(f"py{i}")
            seen.append((captured.c_text, captured.python_text))
    finally:
        os.dup2(saved, 1)
        os.close(saved)
        os.close(fd)

    assert seen == [(f" Setting-Warning: c{i}\n", f"py{i}\n") for i in range(3)]
    assert sink.read_text() == "".join(f" Setting-Warning: c{i}\n" for i in range(3))
    assert capsys.readouterr().out == "py0\npy1\npy2\n"


def test_capture_with_fd1_closed_leaves_it_closed():
    """Under pythonw fd 1 does not exist.  dup2 onto it would *create* fd 1,
    pointing at the spool, and leave it there; the capture must fall back to
    the Python-level text alone."""
    saved = os.dup(1)
    os.close(1)
    try:
        with _OutputCapture() as captured:
            print("still seen")
        with pytest.raises(OSError):
            os.fstat(1)
    finally:
        os.dup2(saved, 1)
        os.close(saved)

    assert captured.python_text == "still seen\n"
    assert captured.c_text == ""


def test_capture_open_across_a_plugin_rerun_releases_its_own_lock(monkeypatch):
    """A `do` of "run plugin.py" re-executes this file while a capture made
    from the old class is open around that very command, rebinding the global
    name to a new class with a new lock.  The capture used to release the new
    class's lock — RuntimeError: release unlocked lock — so the request failed,
    and the old lock stayed held for good."""
    old_class = plugin_module._OutputCapture
    namespace = {k: v for k, v in vars(old_class).items() if k not in ("__dict__", "__weakref__")}
    namespace.update(_spool=None, _spool_lock=threading.Lock())
    rerun_class = type("_OutputCapture", (), namespace)

    with old_class():
        monkeypatch.setattr(plugin_module, "_OutputCapture", rerun_class)
        print("run plugin.py")

    assert old_class._spool_lock.acquire(blocking=False), "the old capture kept its lock"
    old_class._spool_lock.release()
    assert not rerun_class._spool_lock.locked()


# ── generic actions: what PyMOL printed reaches the reply ───────────────────


class CmdException(Exception):
    """Named like pymol.CmdException, which _error_text keys on."""


def test_generic_action_success_carries_pymol_warnings(plugin_server, monkeypatch):
    """cmd.set returns None and raises nothing for a value stored at a level
    PyMOL never reads; the C-level warning on fd 1 is the only trace."""

    def fake_set(*args, **kwargs):
        os.write(1, b' Setting: cartoon_tube_radius set for 43 atoms in object "ub".\n')
        os.write(1, b" Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting\n")

    monkeypatch.setattr(mock_pymol.cmd, "set", fake_set)

    res = plugin_server.handle_request(
        {"action": "set", "args": ["cartoon_tube_radius", "0.5", "ub and resi 1-5"]}
    )

    assert res["status"] == "success"
    assert res["warnings"] == [
        "Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting"
    ]


def test_generic_action_error_gets_the_reason_pymol_printed(plugin_server, monkeypatch):
    """super/align on a bad selection raise CmdException(' Error: ') and print
    the reason only as C feedback: the reply used to be 'PyMOL execution
    error:  Error: ' and nothing else."""

    def fake_super(*args, **kwargs):
        os.write(1, b' Selector-Error: Invalid selection name "nosuch1".\nnosuch1<--\n')
        raise CmdException(" Error: ")

    monkeypatch.setattr(mock_pymol.cmd, "super", fake_super)

    res = plugin_server.handle_request({"action": "super", "args": ["nosuch1", "ub"]})

    assert res["status"] == "error"
    assert 'Selector-Error: Invalid selection name "nosuch1". nosuch1<--' in res["error"]


def test_generic_action_quiet_success_has_no_warnings_key(plugin_server, monkeypatch):
    monkeypatch.setattr(mock_pymol.cmd, "zoom", lambda *a, **k: None)

    res = plugin_server.handle_request({"action": "zoom", "args": ["ub"]})

    assert res == {"status": "success", "result": "Executed 'zoom' successfully."}


@pytest.mark.parametrize("action", ["alter", "iterate_to_list"])
def test_failing_traceback_printer_still_answers(plugin_server, monkeypatch, action):
    """alter('all', 'b=undefined_name') and iterate_to_list with a NameError
    expression returned 'Empty response from PyMOL plugin.': under Python 3.14
    traceback.print_exc() itself raised TypeError ('wrapper.Wrapper' object
    is not iterable) inside the error handler, and nothing was sent back."""

    def name_error(*args, **kwargs):
        raise NameError("name 'undefined_name' is not defined")

    def broken_print_exc(*args, **kwargs):
        raise TypeError("'wrapper.Wrapper' object is not iterable")

    monkeypatch.setattr(mock_pymol.cmd, "alter", name_error)
    monkeypatch.setattr(mock_pymol.cmd, "iterate", name_error)
    monkeypatch.setattr(plugin_module.traceback, "print_exc", broken_print_exc)

    res = plugin_server.handle_request({"action": action, "args": ["all", "b=undefined_name"]})

    assert res["status"] == "error"
    assert "NameError: name 'undefined_name' is not defined" in res["error"]


# ── running a `do` on PyMOL's main thread, and waiting for it ───────────────


@pytest.fixture
def pymol_queue(pymol_api):
    """PyMOL's command queue, drained by this test's thread standing in for
    PyMOL's main thread: from any other thread a line is only queued."""
    api = pymol_api()

    def run_next(timeout=5.0):
        line = api.lines.get(timeout=timeout)
        token = re.fullmatch(r"cmd\._mcpymol_run_job\('(\w+)'\)", line).group(1)
        mock_pymol.cmd._mcpymol_run_job(token)  # what PyMOL's parser does with it
        return token

    api.run_next = run_next
    return api


class FakeCmdApi:
    """Stands in for pymol._cmd, the C layer under cmd.do: do() queues a line
    for the main thread, and refuses it during a modal draw, as CmdDo does."""

    def __init__(self, modal_for=0.0, fail=None):
        self.lines: queue.Queue = queue.Queue()
        self.modal_until = time.monotonic() + modal_for
        self.fail = fail
        self.do_calls = 0
        self.run_next = None

    def _modal(self):
        return time.monotonic() < self.modal_until

    def do(self, cob, line, log, echo):
        self.do_calls += 1
        if self.fail is not None:
            raise self.fail
        if self._modal():
            raise CmdException(" Error: APIEnterNotModal(G)")
        self.lines.put(line)

    def get_modal_draw(self, cob):
        return int(self._modal())


@pytest.fixture
def pymol_api(monkeypatch):
    """Install a FakeCmdApi as cmd._cmd, with a real lock as PyMOL's API lock."""

    def install(**kwargs):
        api = FakeCmdApi(**kwargs)
        monkeypatch.setattr(mock_pymol.cmd, "_cmd", api, raising=False)
        monkeypatch.setattr(mock_pymol.cmd, "lock_api", threading.RLock(), raising=False)
        return api

    return install


def _in_thread(func, *args, **kwargs):
    """Call ``func`` on a thread of its own, as the plugin's socket thread."""
    box: dict = {}

    def target():
        try:
            box["value"] = func(*args, **kwargs)
        except BaseException as e:
            box["error"] = e

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    return box, thread


def test_do_work_runs_on_pymol_main_thread_and_is_waited_for(pymol_queue):
    """Lines run where cmd.do always ran them — PyMOL's main thread; in the
    GUI, 'window size' run anywhere else killed PyMOL — and the reply waits."""
    seen = {}

    def work():
        seen["thread"] = threading.current_thread()
        return 42

    box, worker = _in_thread(plugin_module._on_main_thread, work)
    pymol_queue.run_next()
    worker.join(5)

    assert box == {"value": 42}
    assert seen["thread"] is threading.main_thread()


def test_job_the_main_thread_never_starts_is_cancelled_and_never_runs(pymol_queue):
    """A main thread that is stuck (a modal dialog, a hung command) must not
    hold the request forever.  Once the reply has said nothing was run,
    nothing may run when PyMOL finally reaches the queued line."""
    ran = []

    box, worker = _in_thread(
        plugin_module._on_main_thread, lambda: ran.append(1), start_timeout=0.2
    )
    worker.join(5)

    assert isinstance(box.get("error"), MainThreadTimeout)
    assert "did not pick up the command within 0.2 s" in str(box["error"])
    assert "nothing was run" in str(box["error"])

    token = pymol_queue.run_next()  # PyMOL gets there in the end
    assert ran == []
    assert token not in plugin_module._job_registry()


def test_started_job_is_waited_for_past_the_start_timeout(pymol_queue):
    """The timeout bounds only the wait to *start*: a long command that began
    in time is waited on to the end, or its reply would precede its effect."""

    def slow():
        time.sleep(0.6)
        return "done"

    box, worker = _in_thread(plugin_module._on_main_thread, slow, start_timeout=0.2)
    pymol_queue.run_next()
    worker.join(5)

    assert box == {"value": "done"}


@pytest.mark.parametrize(
    "raised, expected", [(ValueError("bad"), ValueError), (SystemExit(3), RuntimeError)]
)
def test_job_errors_reach_the_waiting_thread_not_pymol(pymol_queue, raised, expected):
    def boom():
        raise raised

    box, worker = _in_thread(plugin_module._on_main_thread, boom)
    pymol_queue.run_next()  # must not raise on PyMOL's main thread
    worker.join(5)

    assert isinstance(box["error"], expected)


def test_on_the_main_thread_already_runs_inline(pymol_api):
    api = pymol_api(fail=AssertionError("queued"))

    assert plugin_module._on_main_thread(lambda: "inline") == "inline"
    assert api.do_calls == 0


def test_cancelled_job_does_nothing_when_it_starts_late():
    """The main thread can pop a job off the registry just as the plugin gives
    up on it; the cancelled flag, checked under the job's lock, closes that gap."""
    ran = []
    job = plugin_module._MainThreadJob(lambda: ran.append(1))

    assert job.cancel_unless_started()
    job.run()

    assert ran == [] and not job.started.is_set()


def test_started_job_cannot_be_cancelled():
    job = plugin_module._MainThreadJob(lambda: "ran")
    job.run()

    assert not job.cancel_unless_started()
    assert job.value == "ran" and job.done.is_set()


def test_start_timeout_becomes_an_error_reply(monkeypatch):
    def stuck(func, start_timeout=None):
        raise MainThreadTimeout("PyMOL's main thread did not pick up the command within 120 s")

    monkeypatch.setattr(plugin_module, "_on_main_thread", stuck)

    res = run_do("color red, all")

    assert res == {
        "status": "error",
        "error": "PyMOL's main thread did not pick up the command within 120 s",
        "output": "",
    }


def test_job_waits_out_a_modal_draw(pymol_api):
    """C CmdDo refuses every line while PyMOL renders a movie with a modal
    draw loop (mpng modal=1): cmd.do raised 'APIEnterNotModal(G)' at once."""
    api = pymol_api(modal_for=0.3)
    started = time.monotonic()

    plugin_module._queue_job("tok", time.monotonic() + 5, 5)

    assert time.monotonic() - started >= 0.3
    assert api.lines.get_nowait() == "cmd._mcpymol_run_job('tok')"


def test_modal_draw_that_outlasts_the_budget_is_reported(pymol_api):
    api = pymol_api(modal_for=60)

    with pytest.raises(MainThreadTimeout, match="rendering a movie"):
        plugin_module._queue_job("tok", time.monotonic() + 0.3, 0.3)
    assert api.lines.empty()


def test_queueing_failure_without_a_modal_draw_is_raised(pymol_api):
    pymol_api(fail=RuntimeError("boom"))

    with pytest.raises(RuntimeError, match="boom"):
        plugin_module._queue_job("tok", time.monotonic() + 5, 5)


def _hold_api_lock(seconds):
    """Hold PyMOL's API lock on another thread, as PyMOL's main thread does
    while it runs a queued command."""
    holding = threading.Event()

    def hold():
        with mock_pymol.cmd.lock_api:
            holding.set()
            time.sleep(seconds)

    threading.Thread(target=hold, daemon=True).start()
    assert holding.wait(5)


def test_a_main_thread_holding_the_api_lock_cannot_hold_the_request(pymol_api):
    """cmd.do takes PyMOL's API lock with no timeout, and the main thread held
    it for the whole of a 4 s command (measured): the start budget elapsed
    unseen inside that wait, and a hung main thread would have held the
    request, and the plugin, forever."""
    api = pymol_api()
    _hold_api_lock(3.0)
    started = time.monotonic()

    with pytest.raises(MainThreadTimeout, match=r"PyMOL was busy for 0\.3 s"):
        plugin_module._queue_job("tok", time.monotonic() + 0.3, 0.3)

    assert time.monotonic() - started < 1.5
    assert api.do_calls == 0


def test_a_lock_freed_within_the_budget_is_taken(pymol_api):
    api = pymol_api()
    _hold_api_lock(0.3)

    plugin_module._queue_job("tok", time.monotonic() + 5, 5)

    assert api.lines.get_nowait() == "cmd._mcpymol_run_job('tok')"
    # From another thread: the lock is re-entrant, so this one would always get it.
    free, _ = _in_thread(lambda: mock_pymol.cmd.lock_api.acquire(timeout=1.0))
    _.join(5)
    assert free == {"value": True}, "the API lock was left held"


# ── replacing the listener while it waits on a `do` ─────────────────────────


def _wait_for(predicate, timeout=5.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return True
        time.sleep(0.02)
    return False


def _listening(port):
    try:
        with socket.create_connection(("127.0.0.1", port), timeout=0.2):
            return True
    except OSError:
        return False


@pytest.fixture
def held_do(monkeypatch):
    """Make every `do` wait until the test releases it — a request the
    listener is still answering while PyMOL's main thread does something
    else, such as re-running plugin.py."""
    entered, release = threading.Event(), threading.Event()

    def run_do_until_released(commands):
        entered.set()
        release.wait(10)
        return {"status": "success", "result": "first listener", "output": ""}

    monkeypatch.setattr(plugin_module, "run_do", run_do_until_released)
    yield SimpleNamespace(entered=entered, release=release)
    release.set()


def test_rerun_through_do_hands_the_port_to_the_new_listener(held_do, monkeypatch):
    """A `do` of "run plugin.py" replaces the server while the old listener
    waits for that very `do`.  The new listener bound at once, failed with
    'address in use', and once the old one had replied nothing listened at
    all.  And stop() spent 2 s joining a listener that could not finish
    before it returned."""
    old = PyMOLSocketServer(port=0)
    old.start()
    assert _wait_for(lambda: old.port != 0)
    port = old.port
    new = None
    try:
        with patch.object(bridge_module, "PORT", port):
            box, client = _in_thread(send_request, "do", args=["run plugin.py"], timeout=10)
            assert held_do.entered.wait(5)

            # What start_plugin() does, on PyMOL's main thread (this one).
            started = time.monotonic()
            old.stop()
            assert time.monotonic() - started < 1.0
            new = PyMOLSocketServer(port=port, replaces=old.thread)
            new.start()

            held_do.release.set()  # the `do` finishes; the old listener replies
            client.join(5)
            assert box["value"]["result"] == "first listener"

            monkeypatch.setattr(
                plugin_module,
                "run_do",
                lambda c: {"status": "success", "result": "second listener", "output": ""},
            )
            assert _wait_for(lambda: not old.thread.is_alive())
            assert _wait_for(lambda: _listening(port)), "no listener after the handover"
            assert send_request("do", args=["x"])["result"] == "second listener"
    finally:
        held_do.release.set()
        old.stop()
        if new is not None:
            new.stop()


def test_stop_then_start_during_a_do_leaves_one_listener(held_do):
    """stop_mcp then start_mcp in one `do`: the old listener finishes its
    request and exits instead of serving on beside its successor."""
    server = PyMOLSocketServer(port=0)
    server.start()
    assert _wait_for(lambda: server.port != 0)
    port = server.port
    try:
        with patch.object(bridge_module, "PORT", port):
            box, client = _in_thread(send_request, "do", args=["stop_mcp\nstart_mcp"], timeout=10)
            assert held_do.entered.wait(5)
            first = server.thread
            server.stop()
            server.start()
            held_do.release.set()
            client.join(5)

            assert box["value"]["status"] == "success"
            assert _wait_for(lambda: not first.is_alive())
            assert server.thread is not first and server.thread.is_alive()
            assert _wait_for(lambda: _listening(port))
    finally:
        held_do.release.set()
        server.stop()


def test_failed_bind_leaves_the_server_restartable():
    """A server whose bind failed stayed `running`, so start_mcp was a no-op
    and only re-running the plugin could bring it back."""
    blocker = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    blocker.bind(("127.0.0.1", 0))
    blocker.listen(1)
    port = blocker.getsockname()[1]
    server = PyMOLSocketServer(port=port)
    try:
        server.start()
        server.thread.join(5)
        assert not server.running

        blocker.close()
        server.start()
        assert _wait_for(lambda: _listening(port))
    finally:
        blocker.close()
        server.stop()


def test_stop_from_its_own_listener_does_not_join_itself():
    """The generic action 'stop_mcp' runs on the listener thread itself."""
    server = PyMOLSocketServer(port=0)
    server.running = True
    server.thread = threading.current_thread()

    server.stop()  # used to raise RuntimeError: cannot join current thread

    assert not server.running


# ── the App Nap guard ────────────────────────────────────────────────────────


class FakeObjC:
    """Stands in for _ObjCRuntime: records activities begun and ended."""

    def __init__(self, fail=False):
        self.begun: list[int] = []
        self.ended: list[int] = []
        self.fail = fail

    def begin_activity(self, options, reason):
        if self.fail:
            raise OSError("no Objective-C runtime here")
        self.begun.append(options)
        return 1000 + len(self.begun)

    def end_activity(self, token):
        self.ended.append(token)


class FakeClock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


@pytest.fixture
def macos_gui(monkeypatch):
    """A Qt GUI PyMOL on macOS, as far as the guard can tell."""
    monkeypatch.setattr(plugin_module.sys, "platform", "darwin")
    monkeypatch.setattr(plugin_module, "_headless", lambda: False)
    monkeypatch.setattr(plugin_module, "_qt_window_exists", lambda: True)


def test_app_nap_guard_holds_one_activity_and_lets_go_when_idle(macos_gui):
    """Taken on the first request and held across requests — per request was
    too late: a napped thread took 2-76 s to reach the call that un-naps it —
    then released after idle_s without one, so an idle PyMOL naps as usual."""
    objc, clock = FakeObjC(), FakeClock()
    guard = _AppNapGuard(idle_s=900, runtime=objc, clock=clock)

    for _ in range(50):
        guard.touch()
        clock.now += 1
    assert objc.begun == [0x00EFFFFF]

    clock.now += 800  # 801 s since the last request
    guard.tick()
    assert objc.ended == [] and guard.held

    clock.now += 100
    guard.tick()
    assert objc.ended == [1001] and not guard.held
    guard.tick()
    assert objc.ended == [1001]

    guard.touch()  # the next request takes it again
    assert objc.begun == [0x00EFFFFF, 0x00EFFFFF]


def test_app_nap_options_still_let_the_mac_sleep():
    """NSActivityUserInitiated (0x00FFFFFF) includes IdleSystemSleepDisabled:
    `pmset -g assertions` showed PreventUserIdleSystemSleep for PyMOL with it,
    which would keep the user's Mac awake whenever PyMOL is open."""
    assert plugin_module._NS_ACTIVITY_NO_APP_NAP == 0x00EFFFFF
    assert not plugin_module._NS_ACTIVITY_NO_APP_NAP & (1 << 20)


@pytest.mark.parametrize(
    "platform, headless, idle_s",
    [("linux", False, 900.0), ("darwin", True, 900.0), ("darwin", False, 0.0)],
    ids=["not-macos", "headless", "disabled-by-zero"],
)
def test_app_nap_guard_stays_off(monkeypatch, platform, headless, idle_s):
    monkeypatch.setattr(plugin_module.sys, "platform", platform)
    monkeypatch.setattr(plugin_module, "_headless", lambda: headless)
    monkeypatch.setattr(plugin_module, "_qt_window_exists", lambda: True)
    objc = FakeObjC()
    guard = _AppNapGuard(idle_s=idle_s, runtime=objc, clock=FakeClock())

    guard.touch()
    guard.touch()

    assert objc.begun == []
    assert guard.enabled is False


def test_app_nap_guard_waits_for_the_qt_window(macos_gui, monkeypatch):
    """The first request can arrive before PyMOL's window exists; that is not
    a reason to decide the session is headless for good."""
    window = {"up": False}
    monkeypatch.setattr(plugin_module, "_qt_window_exists", lambda: window["up"])
    objc = FakeObjC()
    guard = _AppNapGuard(idle_s=900, runtime=objc, clock=FakeClock())

    guard.touch()
    assert objc.begun == [] and guard.enabled is None

    window["up"] = True
    guard.touch()
    assert objc.begun == [0x00EFFFFF]


def test_app_nap_guard_failure_turns_it_off_with_one_line(macos_gui, capsys):
    guard = _AppNapGuard(idle_s=900, runtime=FakeObjC(fail=True), clock=FakeClock())

    guard.touch()
    guard.touch()
    guard.tick()
    guard.release()

    assert capsys.readouterr().out.count("App Nap guard off") == 1
    assert guard.enabled is False and not guard.held


def test_requests_are_served_when_the_app_nap_guard_fails(macos_gui):
    server = PyMOLSocketServer(port=0)
    server.app_nap = _AppNapGuard(idle_s=900, runtime=FakeObjC(fail=True), clock=FakeClock())
    mock_pymol.cmd.get_chains.return_value = ["A"]
    server.start()
    try:
        assert _wait_for(lambda: server.port != 0)
        with patch.object(bridge_module, "PORT", server.port):
            assert send_request("get_chains", args=["all"]) == {
                "status": "success",
                "result": ["A"],
            }
    finally:
        server.stop()


def test_app_nap_guard_is_released_on_stop(macos_gui):
    objc = FakeObjC()
    server = PyMOLSocketServer(port=0)
    server.app_nap = _AppNapGuard(idle_s=900, runtime=objc, clock=FakeClock())
    server.app_nap.touch()
    server.running = True

    server.stop()

    assert objc.ended == [1001]
    assert not server.app_nap.held


@pytest.mark.skipif(sys.platform != "darwin", reason="NSProcessInfo is macOS-only")
def test_objc_runtime_begins_and_ends_a_real_activity():
    """The ctypes prototypes against the real Objective-C runtime.  A wrong
    objc_msgSend signature is a crash no try/except can catch — in the user's
    PyMOL — so it is exercised here, in a process that can afford one."""
    runtime = plugin_module._ObjCRuntime()

    token = runtime.begin_activity(plugin_module._NS_ACTIVITY_NO_APP_NAP, "mcpymol unit test")
    try:
        assert token
    finally:
        runtime.end_activity(token)
