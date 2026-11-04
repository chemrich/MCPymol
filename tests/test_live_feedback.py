"""PyMOL's feedback reaching the caller, against a real PyMOL.

    MCPYMOL_PORT=<port> pytest -m live tests/test_live_feedback.py

WARNING: like tests/test_live.py, this clears the session of the PyMOL it
connects to, and one test re-runs this checkout's plugin inside it.

The mocked suite (test_plugin.py, test_feedback.py) pins the logic; this pins
the facts it rests on, which only PyMOL can confirm: that `do` replies after
the lines have run, that PyMOL's printed errors and C-level warnings are
caught, that a stuck main thread is reported rather than waited on forever, and
that re-running the plugin through `do` leaves one working listener.
"""

import os
import re
import socket
import time

import pytest

import mcpymol.bridge as bridge
from mcpymol.bridge import send_request

pytestmark = pytest.mark.live

OBJ = "_fb_probe"
SEQUENCE = "ACDEFGHIKLMNPQRSTVWY"


def _reachable() -> bool:
    try:
        with socket.create_connection((bridge.HOST, bridge.PORT), timeout=2):
            return True
    except OSError:
        return False


def _do(text: str, timeout: float = 60.0) -> dict:
    return send_request("do", args=[text], timeout=timeout)


def _rebuild_probe() -> None:
    send_request("delete", args=[OBJ])
    send_request("fab", args=[SEQUENCE, OBJ])
    send_request("show", args=["cartoon", OBJ])


@pytest.fixture(scope="module", autouse=True)
def probe():
    if not _reachable():
        pytest.skip(
            f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT} — start PyMOL and "
            f"`run src/mcpymol/plugin.py`",
            allow_module_level=True,
        )
    if "output" not in _do("print('mcpymol feedback probe')", timeout=30):
        pytest.skip(
            "the running PyMOL plugin predates error reporting; load this checkout's "
            "src/mcpymol/plugin.py into it",
            allow_module_level=True,
        )
    _do("reinitialize")
    _rebuild_probe()
    yield
    _do("reinitialize")


@pytest.fixture
def fresh_probe():
    yield
    _rebuild_probe()


def _plugin_value(name: str):
    """A module global of the running plugin, read through `do`.

    PyMOL runs plugin.py in the namespace `do`'s Python lines execute in, so
    its globals are reachable there."""
    res = _do(f"/print(repr(globals().get({name!r}, '<absent>')))")
    return res.get("output", "").strip()


# ── what PyMOL made of a `do` ────────────────────────────────────────────────


def test_a_semicolon_after_alter_is_explained(fresh_probe):
    """alter takes the rest of the line literally, so `; rebuild` became part
    of its expression; the tool used to answer 'Executed command'."""
    from mcpymol.primitives import execute_pymol_command

    out = execute_pymol_command(f"alter {OBJ} and resi 1 and name N, name='N'; rebuild")

    assert "NameError" in out and "rebuild" in out
    assert "takes the rest of the line literally" in out


def test_a_bad_line_is_named_and_the_rest_still_run(fresh_probe):
    from mcpymol.primitives import execute_pymol_command

    out = execute_pymol_command(
        f"set ray_shadow, 0\ncolor red, nosuch_fb\nshow sticks, {OBJ} and resi 1-3"
    )

    assert "PyMOL rejected line 2 of 3: color red, nosuch_fb" in out
    assert 'Invalid selection name "nosuch_fb"' in out
    assert "The other 2 line(s) ran." in out
    shown = send_request("count_atoms", args=[f"{OBJ} and rep sticks"])["result"]
    assert shown > 0, "the line after the rejected one did not run"


def test_what_pymol_printed_comes_back():
    from mcpymol.primitives import execute_pymol_command

    out = execute_pymol_command("get_view")

    assert "PyMOL output:" in out and "set_view" in out


def test_do_replies_after_the_command_has_run(tmp_path):
    """From the plugin's thread cmd.do only queued the line and returned after
    ~0.13 s, finished or not: the reply could precede the work it reported."""
    from mcpymol.primitives import execute_pymol_command

    started = time.monotonic()
    execute_pymol_command("/import time; time.sleep(0.5)")
    assert time.monotonic() - started >= 0.5

    png = tmp_path / "fb.png"
    execute_pymol_command(f"png {png}, width=100, height=80, ray=1")
    assert png.exists() and png.stat().st_size > 0, "the reply came before the file"


def test_pml_log_reads_as_cmd_do_wrote_it(tmp_path):
    """Echo and logging moved from C CmdDo into the plugin; a .pml log must
    still record what ran, with the `_ ` and `PyMOL>` prefixes removed."""
    log = tmp_path / "fb.pml"
    _do(f"log_open {log}")
    _do("_ set ray_shadow, 1\nPyMOL>set ray_shadow, 0")
    _do("log_close")

    assert log.read_text().startswith("set ray_shadow, 1\nset ray_shadow, 0\n")


# ── what PyMOL printed while a tool's command "succeeded" ───────────────────


def test_set_setting_refuses_a_value_pymol_never_uses(fresh_probe):
    from mcpymol.primitives import set_setting

    on_atoms = set_setting("cartoon_tube_radius", "0.5", f"{OBJ} and resi 1-5")
    assert on_atoms.startswith("Error: PyMOL did not apply cartoon_tube_radius = 0.5")
    assert "object-state-level" in on_atoms

    on_object = set_setting("ray_shadow", "0", OBJ)
    assert "global-level" in on_object and "Omit the selection" in on_object

    assert set_setting("cartoon_tube_radius", "0.5", OBJ) == "Executed set successfully."


def test_tools_carry_pymol_warnings_and_reasons(fresh_probe):
    from mcpymol.primitives import align, alter, zoom

    assert "Selector-Error" in zoom(selection="nosuch_fb")
    assert 'Invalid selection name "nosuch_fb"' in align(mobile="nosuch_fb", target=OBJ)

    # A NameError in an alter expression came back as 'Empty response from
    # PyMOL plugin.': under Python 3.14 the traceback printer itself raised.
    out = alter(selection=OBJ, expression="b=undefined_name")
    assert "NameError" in out
    assert "Empty response" not in out


def test_internal_do_callers_are_all_accepted(monkeypatch):
    """Every view still sends lines with `do`; now that PyMOL's verdict comes
    back, none of them may be one PyMOL rejects."""
    import sys

    from mcpymol.structures import fetch_structure

    _do("reinitialize")
    if "Successfully fetched" not in fetch_structure(pdb_code="1hsg"):
        pytest.skip("could not fetch 1hsg")

    real = bridge.send_request
    sent: list[tuple[str, dict]] = []

    def recording(action, args=None, kwargs=None, timeout=bridge._DEFAULT_TIMEOUT):
        res = real(action, args=args, kwargs=kwargs, timeout=timeout)
        if action == "do":
            sent.append(((args or [""])[0], res))
        return res

    for name, module in list(sys.modules.items()):
        if name.startswith("mcpymol") and getattr(module, "send_request", None) is real:
            monkeypatch.setattr(module, "send_request", recording)

    from mcpymol.comparison import superposition_view
    from mcpymol.views import (
        bfactor_view,
        crosslink_view,
        interface_view,
        ligand_view,
        mutation_view,
        pocket_view,
        putty_view,
        textbook_view,
    )

    try:
        ligand_view(obj_name="1hsg", ligand_resn="MK1")
        pocket_view(obj_name="1hsg", resn="MK1")
        bfactor_view(obj_name="1hsg")
        putty_view(obj_name="1hsg")
        crosslink_view(obj_name="1hsg")
        mutation_view(obj_name="1hsg", mutations="A25N")
        textbook_view(obj_name="1hsg")
        interface_view(obj_name="1hsg", chain_a="A", chain_b="B")
        real("do", args=["create _fb_a, 1hsg and chain A and polymer"])
        real("do", args=["create _fb_b, 1hsg and chain B and polymer"])
        superposition_view(mobile="_fb_a", target="_fb_b")
    finally:
        _do("reinitialize")
        _rebuild_probe()

    assert len(sent) > 20, f"the views sent only {len(sent)} `do` requests"
    rejected = [(line, res.get("error")) for line, res in sent if res.get("status") != "success"]
    assert not rejected, rejected


# ── the plugin itself: never wedged ──────────────────────────────────────────


def test_a_main_thread_that_never_gets_there_is_reported_not_waited_on():
    """While PyMOL's main thread is busy, a `do` waits only
    MCPYMOL_DO_START_TIMEOUT for it; then it says nothing was run, nothing
    does run later, and the plugin keeps answering.

    Headless, the main thread held PyMOL's API lock for the whole busy line,
    so the refusal comes from the bounded wait for that lock ('PyMOL was busy
    ...'); where it does not hold it, from the wait for the queued line to
    start ('did not pick up the command ...').  Either way nothing ran."""
    before = _plugin_value("DO_START_TIMEOUT")
    if before in ("", "'<absent>'"):
        pytest.skip("the plugin's globals are not reachable through `do`")
    assert _do("set ray_shadow, 0")["status"] == "success"
    assert _do("/DO_START_TIMEOUT = 1.0")["status"] == "success"
    busy_until = None
    try:
        # A helper thread queues a 4 s line, which PyMOL's main thread then runs.
        _do(
            "/import threading, time; threading.Thread(target=lambda: (time.sleep(0.2), "
            "cmd.do('/import time; time.sleep(4)')), daemon=True).start()"
        )
        busy_until = time.monotonic() + 4.2
        time.sleep(0.8)

        started = time.monotonic()
        stuck = _do("set ray_shadow, 1", timeout=30)
        elapsed = time.monotonic() - started

        assert stuck["status"] == "error", stuck
        assert "within 1 s" in stuck["error"] or "busy for 1 s" in stuck["error"]
        assert "nothing was run" in stuck["error"]
        assert elapsed < 3.0, f"waited {elapsed:.1f} s on a 1 s budget"

        time.sleep(max(0.0, busy_until - time.monotonic()) + 0.5)
        assert send_request("get_setting_int", args=["ray_shadow"])["result"] == 0
        assert _do("print('still answering')")["output"] == "still answering\n"
    finally:
        # Until the busy line ends every `do`, this restore included, is refused.
        if busy_until is not None:
            time.sleep(max(0.0, busy_until - time.monotonic()) + 0.5)
        for _ in range(10):
            if _do(f"/DO_START_TIMEOUT = {before}")["status"] == "success":
                break
            time.sleep(1.0)


def test_a_movie_rendering_in_a_modal_draw_is_waited_out(tmp_path):
    """PyMOL refuses every queued line during an mpng with modal=1; `do`
    failed at once with 'APIEnterNotModal(G)'.  It now waits for the movie."""
    setup = _do(f"mset 1x20\nshow cartoon, {OBJ}")
    assert setup["status"] == "success", setup
    try:
        res = send_request(
            "mpng",
            args=[str(tmp_path / "frame")],
            kwargs={"modal": 1, "width": 120, "height": 90},
            timeout=120,
        )
        assert res["status"] == "success", res
        after = _do("print('after the movie')", timeout=120)

        assert after["status"] == "success", after
        assert after["output"] == "after the movie\n"
        assert len(list(tmp_path.glob("frame*.png"))) == 20
    finally:
        _do("mclear\nmset")


def test_rerunning_the_plugin_through_do_keeps_a_listener():
    """`execute_pymol_command("run .../plugin.py")` re-runs the plugin while
    its listener waits for that very `do`: the new listener could not bind, and
    once the old one replied nothing listened at all."""
    from mcpymol.cli import plugin_path
    from mcpymol.primitives import execute_pymol_command

    out = execute_pymol_command(f"run {plugin_path()}")

    assert out.startswith(f"Executed command: run {plugin_path()}"), out
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        res = send_request("count_atoms", args=[OBJ], timeout=5)
        if res.get("status") == "success":
            break
        time.sleep(0.2)
    assert res.get("status") == "success", res
    threads = _do(
        "/import threading; print(sum(t.name.endswith('(_listen_loop)') for t in threading.enumerate()))"
    )
    assert threads["output"].strip() == "1"


# ── GUI PyMOL only ───────────────────────────────────────────────────────────


@pytest.fixture
def gui():
    if "False" not in _do("/print(bool(cmd._pymol.invocation.options.no_gui))")["output"]:
        pytest.skip("needs a GUI PyMOL (pymol -kq plugin.py)")


def test_gui_window_lines_run_on_the_gui_thread(gui):
    """'window size' run off the GUI thread killed PyMOL outright."""
    from mcpymol.primitives import execute_pymol_command

    execute_pymol_command("window size, 640, 480")

    assert send_request("count_atoms", args=[OBJ])["status"] == "success"


def test_gui_rapid_two_line_requests_never_hang(gui):
    """PyMOL's own MainThreadCaller lost a wake-up and hung by the 11th
    request; the command queue must not."""
    for i in range(300):
        res = _do(f"set ray_shadow, {i % 2}\ncount_atoms {OBJ}", timeout=30)
        assert res["status"] == "success", (i, res)


@pytest.mark.skipif(os.uname().sysname != "Darwin", reason="App Nap is macOS-only")
def test_gui_holds_off_app_nap_while_driven(gui):
    held = _do("/print(mcp_bridge_plugin.app_nap.held)")
    if "True" not in held["output"] and "False" not in held["output"]:
        pytest.skip("the plugin's globals are not reachable through `do`")
    if _plugin_value("APP_NAP_IDLE_S") in ("0.0", "0"):
        pytest.skip("MCPYMOL_APP_NAP_IDLE_S=0 turned the guard off")

    assert held["output"].strip() == "True"
    assert re.fullmatch(r"\d+", _do("/print(mcp_bridge_plugin.app_nap._token)")["output"].strip())
