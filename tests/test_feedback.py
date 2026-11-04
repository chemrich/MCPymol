"""What PyMOL says reaching the caller, and honest words when it says nothing.

Covers the tool side of the feedback fixes (the plugin side is in
test_plugin.py): the bridge's errors for a refused connection versus a reply
that never came, sleep measured across a request, PyMOL's warnings appended to
tool results, execute_pymol_command and set_setting reporting what PyMOL made
of their command, and render saying when the computer slept through it.

tests/test_live_feedback.py checks the same contracts against a real PyMOL.
"""

import asyncio
import json
import queue
import re
import socket
import sys
import threading
import time
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

mock_pymol = sys.modules["pymol"]  # conftest's mock: the one plugin.py bound

from mcp.server.fastmcp import Image

import mcpymol.bridge as bridge_module
import mcpymol.plugin as plugin_module
import mcpymol.server  # noqa: F401  (registers every tool on the app)
from mcpymol.app import mcp
from mcpymol.bridge import _SLOW_OP_TIMEOUT, _call, _clock_with_sleep, send_request
from mcpymol.cli import plugin_path
from mcpymol.plugin import PyMOLSocketServer
from mcpymol.primitives import execute_pymol_command, intra_rms, rms_cur, set_setting
from mcpymol.rendering import render, turntable
from mcpymol.wiggles.port import PortError, SendRequestPort

# ── the bridge: what went wrong, in words that point the right way ──────────


@pytest.fixture
def fake_socket():
    with patch("socket.socket") as cls:
        inst = MagicMock()
        cls.return_value.__enter__.return_value = inst
        yield inst


def _replies(*documents):
    """recv() side effect: each JSON document, then EOF."""
    return [json.dumps(d).encode() for d in documents] + [b""]


def test_timeout_after_connecting_does_not_blame_the_plugin(fake_socket):
    """Every timeout read 'Socket connection failed: timed out. Is the PyMOL
    plugin running?'.  But a render that outlives its budget was accepted and
    may still be running — re-sending it only queues the work twice."""
    fake_socket.recv.side_effect = TimeoutError("timed out")

    res = send_request("png", args=["x.png"], timeout=600)

    assert res["status"] == "error"
    assert "PyMOL accepted 'png' but sent no reply within 600 s" in res["error"]
    assert "may still be running it" in res["error"]
    assert "was not cancelled" in res["error"]
    assert "Is the PyMOL plugin running?" not in res["error"]


def test_refused_connection_says_nothing_is_listening(fake_socket):
    fake_socket.connect.side_effect = ConnectionRefusedError("[Errno 61] Connection refused")

    res = send_request("count_atoms")

    assert res["error"].startswith(
        f"Socket connection failed: nothing is listening on {bridge_module.HOST}:"
        f"{bridge_module.PORT}"
    )
    assert "Is the PyMOL plugin running?" in res["error"]


def test_timeout_while_connecting_is_still_a_connection_failure(fake_socket):
    fake_socket.connect.side_effect = TimeoutError("timed out")

    res = send_request("count_atoms")

    assert res["error"].startswith("Socket connection failed: timed out connecting to")


def _clocks(monkeypatch, awake, total):
    """Make the bridge's two clocks read the given values in turn: ``awake``
    is time.monotonic() (stops while the machine sleeps), ``total`` the clock
    that keeps counting."""
    awake_reads, total_reads = iter(awake), iter(total)
    monkeypatch.setattr(bridge_module, "time", SimpleNamespace(monotonic=lambda: next(awake_reads)))
    monkeypatch.setattr(bridge_module, "_clock_with_sleep", lambda: next(total_reads))


def test_a_sleep_during_a_request_is_measured(fake_socket, monkeypatch):
    """The two renders that 'stalled for ten minutes' had the Mac asleep for
    579 s and 716 s of their wait; nothing said so."""
    _clocks(monkeypatch, awake=[100.0, 110.0], total=[1000.0, 1600.0])
    fake_socket.recv.side_effect = _replies({"status": "success", "result": None})

    res = send_request("png", args=["x.png"])

    assert res["asleep_s"] == 590.0


def test_a_timeout_that_spanned_a_sleep_says_so(fake_socket, monkeypatch):
    _clocks(monkeypatch, awake=[100.0, 110.0], total=[1000.0, 1600.0])
    fake_socket.recv.side_effect = TimeoutError("timed out")

    res = send_request("png", args=["x.png"], timeout=10)

    assert res["error"].endswith(" The computer was asleep for 590 s of this wait.")


def test_a_short_sleep_is_not_worth_mentioning(fake_socket, monkeypatch):
    _clocks(monkeypatch, awake=[100.0, 110.0], total=[1000.0, 1013.0])
    fake_socket.recv.side_effect = _replies({"status": "success", "result": 1})

    assert send_request("count_atoms") == {"status": "success", "result": 1}


def test_no_sleep_clock_means_no_sleep_reported(fake_socket, monkeypatch):
    """Windows has no clock that counts sleep; nothing is measured there."""
    monkeypatch.setattr(bridge_module, "_clock_with_sleep", lambda: None)
    fake_socket.recv.side_effect = _replies({"status": "success", "result": 1})

    assert send_request("count_atoms") == {"status": "success", "result": 1}


@pytest.mark.parametrize(
    "platform, clocks, expected",
    [
        ("darwin", {"CLOCK_MONOTONIC": 1, "CLOCK_BOOTTIME": 7}, 1),
        ("linux", {"CLOCK_MONOTONIC": 1, "CLOCK_BOOTTIME": 7}, 7),
        ("win32", {}, None),
    ],
)
def test_clock_with_sleep_is_the_one_that_counts_sleep(monkeypatch, platform, clocks, expected):
    """macOS time.monotonic() is mach_absolute_time, which stops while the
    machine sleeps; CLOCK_MONOTONIC there keeps counting (measured: 218 h
    against 671 h since boot).  On Linux the counting one is CLOCK_BOOTTIME."""
    read = []
    fake_time = SimpleNamespace(clock_gettime=lambda clk: read.append(clk) or 5.0, **clocks)
    monkeypatch.setattr(bridge_module, "time", fake_time)
    monkeypatch.setattr(bridge_module.sys, "platform", platform)

    value = _clock_with_sleep()

    assert read == ([] if expected is None else [expected])
    assert value == (None if expected is None else 5.0)


def test_a_reply_that_is_not_an_object_is_an_error(fake_socket):
    fake_socket.recv.side_effect = [b"null", b""]

    res = send_request("count_atoms")

    assert res["status"] == "error"
    assert "expected a JSON object, got NoneType" in res["error"]


@patch("mcpymol.bridge.send_request")
def test_call_appends_what_pymol_warned(mock_sr):
    """zoom on a selection that matches nothing 'succeeds'; the selector error
    PyMOL printed was the only sign, and it reached only PyMOL's console."""
    mock_sr.return_value = {
        "status": "success",
        "result": None,
        "warnings": ['Selector-Error: Invalid selection name "nosuch". nosuch<--'],
    }

    assert _call("zoom", selection="nosuch") == (
        'Executed zoom successfully.\nPyMOL warned: Selector-Error: Invalid selection name "nosuch". nosuch<--'
    )


@patch("mcpymol.bridge.send_request")
def test_call_lists_several_warnings_one_per_line(mock_sr):
    mock_sr.return_value = {
        "status": "success",
        "result": 2.5,
        "warnings": ["A-Warning: a", "B-Error: b"],
    }

    out = _call("get_distance", atom1="a", atom2="b", _measures=("Distance", "A"))

    assert out == "Distance: 2.50 A.\nPyMOL warned:\n  A-Warning: a\n  B-Error: b"


@patch("mcpymol.bridge.send_request")
def test_call_without_warnings_is_unchanged(mock_sr):
    mock_sr.return_value = {"status": "success", "result": None}

    assert _call("zoom", selection="all") == "Executed zoom successfully."


def test_unit_tests_cannot_reach_a_real_pymol():
    """tests/conftest.py refuses connections to PyMOL's port outside live
    tests. Without the guard this raises TypeError instead (the C method
    rejects a non-socket), so the test cannot pass by accident."""
    with pytest.raises(ConnectionRefusedError, match="real PyMOL"):
        socket.socket.connect(object(), ("127.0.0.1", 9876))


# ── selections PyMOL crashes on are vetted first ─────────────────────────────


def _pymol_with(counts):
    """send_request stand-in: count_atoms answers from ``counts`` (an error
    for a name it lacks, as PyMOL does for an invalid selection)."""

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "count_atoms":
            if args[0] not in counts:
                return {
                    "status": "error",
                    "error": f'PyMOL execution error:  Error: Invalid selection name "{args[0]}".',
                }
            return {"status": "success", "result": counts[args[0]]}
        return {"status": "success", "result": 0.25}

    return fake


@patch("mcpymol.bridge.send_request")
def test_rms_cur_never_sends_a_selection_pymol_would_crash_on(mock_sr):
    """`rms_cur ub, nosuch` aborted PyMOL itself — the user's whole session —
    with a double free in ExecutiveRMS; the tool answered 'Empty response'."""
    mock_sr.side_effect = _pymol_with({"ub": 660})

    out = rms_cur(mobile="ub", target="nosuch")

    assert out.startswith("Error: rms_cur was not run — the target selection 'nosuch' is not valid")
    assert 'Invalid selection name "nosuch"' in out
    assert "rms_cur" not in [c.args[0] for c in mock_sr.call_args_list]


@patch("mcpymol.bridge.send_request")
def test_rms_cur_refuses_an_empty_selection(mock_sr):
    mock_sr.side_effect = _pymol_with({"ub": 660, "ub and resi 999": 0})

    out = rms_cur(mobile="ub and resi 999", target="ub")

    assert out == (
        "Error: rms_cur was not run — the mobile selection 'ub and resi 999' matches no atoms."
    )


@patch("mcpymol.bridge.send_request")
def test_rms_cur_runs_once_both_selections_resolve(mock_sr):
    mock_sr.side_effect = _pymol_with({"a": 76, "b": 76})

    assert rms_cur(mobile="a", target="b") == "RMSD (as positioned): 0.25 A."
    assert [c.args[0] for c in mock_sr.call_args_list] == ["count_atoms", "count_atoms", "rms_cur"]


@pytest.mark.parametrize(
    "error",
    [
        "PyMOL accepted 'count_atoms' but sent no reply within 10 s. It may still be running it.",
        "Socket connection failed: nothing is listening on 127.0.0.1:9876 (refused). "
        "Is the PyMOL plugin running?",
    ],
)
@patch("mcpymol.bridge.send_request")
def test_rms_cur_guard_does_not_call_a_busy_pymol_an_invalid_selection(mock_sr, error):
    """With PyMOL busy or down, the guard's own check fails; calling the
    selection invalid would send a model off to "fix" a correct one."""
    mock_sr.return_value = {"status": "error", "error": error}

    out = rms_cur(mobile="1ubq and name CA", target="1ubq_b and name CA")

    assert "is not valid" not in out
    assert out == (
        "Error: rms_cur was not run — the mobile selection '1ubq and name CA' could not be "
        f"checked first (PyMOL crashes on an invalid selection here): {error}"
    )
    assert [c.args[0] for c in mock_sr.call_args_list] == ["count_atoms"]


@patch("mcpymol.bridge.send_request")
def test_intra_rms_never_sends_a_selection_pymol_would_crash_on(mock_sr):
    """`intra_rms nosuch` segfaulted PyMOL in CmdIntraFit."""
    mock_sr.side_effect = _pymol_with({})

    out = intra_rms(selection="nosuch")

    assert out.startswith("Error: intra_rms was not run — the selection 'nosuch' is not valid")
    assert [c.args[0] for c in mock_sr.call_args_list] == ["count_atoms"]


# ── execute_pymol_command: PyMOL's words, after the command has run ─────────


@patch("mcpymol.primitives.send_request")
def test_execute_returns_what_pymol_printed_and_warned(mock_sr):
    mock_sr.return_value = {
        "status": "success",
        "result": "Executed command: get_view\nzoom nosuch",
        "output": "set_view (\\\n  1.0, 0.0 )\n",
        "warnings": ['Selector-Error: Invalid selection name "nosuch". nosuch<--'],
    }

    out = execute_pymol_command("get_view\nzoom nosuch")

    assert out == (
        "Executed command: get_view\nzoom nosuch\n"
        "PyMOL output:\nset_view (\\\n  1.0, 0.0 )\n"
        'PyMOL warned:\n  Selector-Error: Invalid selection name "nosuch". nosuch<--'
    )


@patch("mcpymol.primitives.send_request")
def test_execute_waits_as_long_as_a_slow_command_needs(mock_sr):
    """`do` now replies once the command has run, so a ray-traced png through
    here needs the slow budget, not the 10 s default."""
    mock_sr.return_value = {"status": "success", "result": "x", "output": ""}

    execute_pymol_command("png x.png, ray=1")

    assert mock_sr.call_args.kwargs["timeout"] == _SLOW_OP_TIMEOUT


@patch("mcpymol.primitives.send_request")
def test_execute_passes_pymol_rejection_through(mock_sr):
    error = (
        "PyMOL rejected line 2 of 3: color red, nosuch\n"
        ' Error: Invalid selection name "nosuch".\nnosuch<--\nThe other 2 line(s) ran.'
    )
    mock_sr.return_value = {"status": "error", "error": error, "output": ""}

    assert execute_pymol_command("a\ncolor red, nosuch\nb") == error


@patch("mcpymol.primitives.send_request")
def test_execute_keeps_output_and_warnings_when_a_line_fails(mock_sr):
    """The lines that ran may have printed what was asked for, or warned that
    a value they set is never used; a failure elsewhere used to drop both, and
    'The other 2 line(s) ran' then read as a clean bill for them."""
    error = (
        "PyMOL rejected line 3 of 3: color red, typo_sel\n"
        ' Error: Invalid selection name "typo_sel".\ntypo_sel<--\nThe other 2 line(s) ran.'
    )
    warning = "Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting"
    mock_sr.return_value = {
        "status": "error",
        "error": error,
        "output": "PyMOL>iterate p1 and name CA, print(resn)\nMET\nGLN\n",
        "warnings": [warning],
    }

    out = execute_pymol_command("iterate ...\nset ...\ncolor red, typo_sel")

    assert out.startswith(error)
    assert "PyMOL output:\nPyMOL>iterate p1 and name CA, print(resn)\nMET\nGLN" in out
    assert out.endswith(f"PyMOL warned:\n  {warning}")
    assert "predates error reporting" not in out


@patch("mcpymol.primitives.send_request")
def test_execute_on_an_old_plugin_says_errors_went_unseen(mock_sr):
    """A plugin from before error reporting replies success to everything
    without running it first; its replies have no `output`."""
    mock_sr.return_value = {"status": "success", "result": "Executed command: color red, x"}

    out = execute_pymol_command("color red, x")

    assert out.startswith("Executed command: color red, x\n")
    assert "predates error reporting" in out
    # The path of this install's plugin, not a checkout-relative one.
    assert f"run {plugin_path()}" in out


@patch("mcpymol.primitives.send_request")
def test_execute_quiet_success_is_just_the_command(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "x", "output": ""}

    assert execute_pymol_command("hide everything") == "Executed command: hide everything"


# ── set_setting: a value PyMOL stores but never uses is not a success ───────


@pytest.mark.parametrize(
    "warning, advice",
    [
        (
            "Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting",
            "Pass an object name to set it for a whole object, or omit the selection",
        ),
        ("Setting-Warning: 'ray_shadow' is a global-level setting", "Omit the selection"),
    ],
)
@patch("mcpymol.primitives.send_request")
def test_set_setting_reports_a_value_pymol_will_never_use(mock_sr, warning, advice):
    """`set cartoon_tube_radius, 0.5, ub and resi 1-5` is stored per atom and
    never drawn; PyMOL only warns, and the tool said 'Executed set
    successfully.'"""
    mock_sr.return_value = {"status": "success", "result": None, "warnings": [warning]}
    name = re.search(r"'(\w+)'", warning).group(1)

    out = set_setting(name, "0.5", "ub and resi 1-5")

    assert out.startswith(f"Error: PyMOL did not apply {name} = 0.5 to 'ub and resi 1-5': ")
    assert "-level setting" in out
    assert advice in out


@patch("mcpymol.primitives.send_request")
def test_set_setting_success_and_other_warnings(mock_sr):
    mock_sr.return_value = {"status": "success", "result": None}
    assert set_setting("cartoon_tube_radius", "0.5", "ub") == "Executed set successfully."
    assert mock_sr.call_args.kwargs["args"] == ["cartoon_tube_radius", "0.5", "ub"]

    set_setting("ray_shadow", "0")
    assert mock_sr.call_args.kwargs["args"] == ["ray_shadow", "0"]

    mock_sr.return_value = {"status": "success", "result": None, "warnings": ["X-Warning: y"]}
    assert (
        set_setting("ray_shadow", "0") == "Executed set successfully.\nPyMOL warned: X-Warning: y"
    )


@patch("mcpymol.primitives.send_request")
def test_set_setting_error_passes_through(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "PyMOL execution error: unknown setting"}

    assert set_setting("nosuch", "1") == "PyMOL execution error: unknown setting"


# ── render: when the computer slept through it, say so ──────────────────────


def _png():
    import struct
    import zlib

    def chunk(tag, data):
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", zlib.crc32(tag + data))

    raw = b"".join(
        b"\x00" + b"".join(bytes(((x * 7) % 256, (y * 11) % 256, 90)) for x in range(8))
        for y in range(8)
    )
    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", 8, 8, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw))
        + chunk(b"IEND", b"")
    )


def _png_writer(asleep_s=None, status="success"):
    """send_request stand-in: writes a PNG where asked, and reports a sleep."""

    def fake(action, args=None, kwargs=None, **_ignored):
        res: dict = {"status": status, "result": None}
        if status == "error":
            res["error"] = "ray tracing failed"
        elif action == "png":
            with open(args[0], "wb") as fh:
                fh.write(_png())
        if asleep_s is not None:
            res["asleep_s"] = asleep_s
        return res

    return fake


def test_render_says_when_the_computer_slept_through_it():
    """'The render stalled until the next command' was the Mac asleep: awake
    28 of 607 s.  The image comes back with a note saying so."""
    with patch("mcpymol.rendering.send_request", side_effect=_png_writer(asleep_s=700.0)):
        result = render(width=8, height=8)

    assert isinstance(result, list) and len(result) == 2
    image, note = result
    assert isinstance(image, Image)
    assert "the computer was asleep for 700 s of it" in note
    assert "PyMOL only works while the machine is awake" in note


def test_render_without_a_sleep_is_just_the_image():
    with patch("mcpymol.rendering.send_request", side_effect=_png_writer()):
        assert isinstance(render(width=8, height=8), Image)


def test_render_error_after_a_sleep_says_so():
    with patch("mcpymol.rendering.send_request", side_effect=_png_writer(700.0, "error")):
        result = render(width=8, height=8)

    assert result.startswith("Error rendering: ray tracing failed")
    assert "asleep for 700 s" in result


def test_render_note_reaches_the_client_as_text_after_the_image():
    """FastMCP must turn [Image, note] into an image block and a text block;
    a list it could not convert would reach the client as a JSON dump."""
    with patch("mcpymol.rendering.send_request", side_effect=_png_writer(asleep_s=700.0)):
        content = asyncio.run(mcp.call_tool("render", {"width": 8, "height": 8}))

    assert [block.type for block in content] == ["image", "text"]
    assert "asleep for 700 s" in content[1].text


def test_turntable_says_when_the_computer_slept(tmp_path):
    with patch("mcpymol.rendering.send_request", side_effect=_png_writer(asleep_s=300.0)):
        out = turntable(frames=2, out_dir=str(tmp_path))

    assert out.startswith("Wrote 2 ray-traced frames")
    assert "asleep for 600 s" in out


# ── wiggles: a rejected line is an exception, as its views expect ───────────


def test_wiggles_do_raises_on_a_line_pymol_rejected():
    port = SendRequestPort(
        lambda action, args=None, kwargs=None, **kw: {
            "status": "error",
            "error": 'PyMOL rejected the command: color red, nosuch_wiggles\n Error: Invalid selection name "nosuch_wiggles".',
            "output": "",
        }
    )

    with pytest.raises(PortError, match="do failed: PyMOL rejected the command"):
        port.do("color red, nosuch_wiggles")


def test_wiggles_do_returns_the_result_of_a_new_style_success():
    port = SendRequestPort(
        lambda action, args=None, kwargs=None, **kw: {
            "status": "success",
            "result": "Executed command: stored.wiggles_b = []",
            "output": "",
        }
    )

    assert port.do("stored.wiggles_b = []") == "Executed command: stored.wiggles_b = []"


# ── end to end: plugin, wire, bridge and tool ────────────────────────────────


class _OneLineFails:
    """A parser that rejects line 2 the way PyMOL's does: print, return 0."""

    def __init__(self):
        self.layer = {1: SimpleNamespace(embed_sentinel=None, cont="")}

    def parse(self, line, secure=0):
        if "nosuch" in line:
            print(' Error: Invalid selection name "nosuch".\nnosuch<--')
            return 0
        return 1


def test_execute_reports_a_rejected_line_over_the_real_wire(monkeypatch):
    """The whole chain: the plugin's listener queues the lines for PyMOL's
    main thread (here a thread draining a queue, as PyMOL's does), waits for
    them to run, and the tool returns PyMOL's verdict on each."""
    monkeypatch.setattr(plugin_module, "_new_parser", _OneLineFails)
    queued: queue.Queue = queue.Queue()
    monkeypatch.setattr(
        mock_pymol.cmd,
        "_cmd",
        SimpleNamespace(
            do=lambda cob, line, log, echo: queued.put(line), get_modal_draw=lambda cob: 0
        ),
        raising=False,
    )
    monkeypatch.setattr(mock_pymol.cmd, "lock_api", threading.RLock(), raising=False)

    def pymol_main_loop():
        while (line := queued.get()) is not None:
            token = re.fullmatch(r"cmd\._mcpymol_run_job\('(\w+)'\)", line).group(1)
            mock_pymol.cmd._mcpymol_run_job(token)

    main_loop = threading.Thread(target=pymol_main_loop, daemon=True)
    main_loop.start()
    server = PyMOLSocketServer(port=0)
    server.start()
    try:
        deadline = time.monotonic() + 5
        while server.port == 0 and time.monotonic() < deadline:
            time.sleep(0.01)
        with patch.object(bridge_module, "PORT", server.port):
            out = execute_pymol_command("set ray_shadow, 0\ncolor red, nosuch\nshow sticks")
    finally:
        server.stop()
        queued.put(None)
        main_loop.join(5)

    assert out == (
        "PyMOL rejected line 2 of 3: color red, nosuch\n"
        ' Error: Invalid selection name "nosuch".\nnosuch<--\n'
        "The other 2 line(s) ran."
    )
