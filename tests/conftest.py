"""Shared test setup.

``mcpymol.plugin`` does ``from pymol import cmd`` at import time and binds the
result as a module global, so the mock has to be installed before *any* test
module imports it — and it has to be the *same* mock for every test file, or a
test asserting on its own mock silently watches an object the plugin never
uses.  conftest is imported before the test modules, which makes this the one
place that can guarantee both.
"""

import os
import socket
import sys
from unittest.mock import MagicMock

import pytest

mock_pymol = MagicMock()
sys.modules.setdefault("pymol", mock_pymol)
sys.modules.setdefault("pymol.cmd", mock_pymol.cmd)


# ── Live runs: only against a port named for testing, never a vacuous green ──
#
# The live suites clear the session they reach. "Not 9876" is not a safe rule:
# the README tells people with a port clash to move their everyday PyMOL with
# MCPYMOL_PORT, so the port in a shell can be a real session too. Live tests
# therefore run only against MCPYMOL_TEST_PORT, a port someone named for
# testing. It is copied into MCPYMOL_PORT here, before any test module imports
# the bridge (which reads the port once, at import).
TEST_PORT = "MCPYMOL_TEST_PORT"
MIN_PASSED = "MCPYMOL_LIVE_MIN_PASSED"
DEFAULT_PORT = 9876
_live_run = {"selected": 0, "passed": 0, "required": 1}


def _port_number(value: str | None) -> int | None:
    try:
        return int(value) if value else None
    except ValueError:
        return None


_SHELL_PORT = os.environ.get("MCPYMOL_PORT")
if _port_number(os.environ.get(TEST_PORT)) is not None:
    # A malformed value is left for the refusal below to name, rather than
    # crashing the bridge's import with a bare ValueError.
    os.environ["MCPYMOL_PORT"] = os.environ[TEST_PORT]


def _pymol_ports() -> set[int]:
    """The ports a real PyMOL plugin may be listening on, on this machine:
    the default, the shell's MCPYMOL_PORT, and the test port. Unit tests may
    reach none of them."""
    ports = {DEFAULT_PORT}
    for value in (_SHELL_PORT, os.environ.get(TEST_PORT)):
        port = _port_number(value)
        if port is not None:
            ports.add(port)
    return ports


def live_run_refusal(selected: int, test_port: str | None, min_passed: str | None) -> str | None:
    """Why a run selecting ``selected`` live tests must not start, or None."""
    if not selected:
        return None
    if not test_port:
        return (
            f"{selected} live tests would clear the session of whatever PyMOL they reach, "
            f"so they run only against a port named for testing. Start one with "
            f"`MCPYMOL_PORT=9898 pymol -ckqK src/mcpymol/plugin.py` and run "
            f"`{TEST_PORT}=9898 pytest -m live`."
        )
    if _port_number(test_port) is None:
        return f"{TEST_PORT}={test_port!r} is not a port number."
    if min_passed and _port_number(min_passed) is None:
        return f"{MIN_PASSED}={min_passed!r} is not a number."
    return None


def live_shortfall(selected: int, passed: int, required: int, port: str | None) -> str | None:
    """Why a run that selected live tests failed for want of passes, or None."""
    if not selected or passed >= required:
        return None
    return (
        f"only {passed} of {selected} live tests passed; at least {required} must "
        f"(is a PyMOL with the plugin listening on port {port}?)"
    )


@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config, items):
    """Refuse a live run that has no test port, before anything runs.

    trylast, so it sees what ``-m`` left selected: a default run deselects
    every live test and is never affected. It fails the run rather than
    skipping, because a skipped live suite reads as a green one, and it fires
    before any fixture, including the one that sends ``cd`` to the live PyMOL.
    ``--collect-only`` runs nothing, so it is let through.
    """
    if config.option.collectonly:
        return
    selected = sum(1 for item in items if item.get_closest_marker("live"))
    refusal = live_run_refusal(selected, os.environ.get(TEST_PORT), os.environ.get(MIN_PASSED))
    if refusal:
        raise pytest.UsageError(refusal)
    required = _port_number(os.environ.get(MIN_PASSED)) or 1
    _live_run.update(selected=selected, passed=0, required=required)


def pytest_runtest_logreport(report):
    if report.when == "call" and report.passed and "live" in report.keywords:
        _live_run["passed"] += 1


def pytest_sessionfinish(session, exitstatus):
    """Fail a run that asked for live tests and got too few passes.

    The live suites skip themselves when nothing answers on the port, so a run
    aimed at the wrong port, or at no PyMOL at all, used to end green with
    every live test skipped. MCPYMOL_LIVE_MIN_PASSED raises the bar (a release
    gate sets it near the full count); the default is one. Only a run that
    would otherwise pass is changed: an interrupted or failing run keeps its
    own exit code and its own explanation.
    """
    if exitstatus != pytest.ExitCode.OK:
        return
    shortfall = live_shortfall(
        _live_run["selected"],
        _live_run["passed"],
        _live_run["required"],
        os.environ.get(TEST_PORT),
    )
    if shortfall is None:
        return
    session.exitstatus = pytest.ExitCode.TESTS_FAILED
    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    if reporter is not None:
        reporter.write_line(f"ERROR: {shortfall}", red=True, bold=True)


@pytest.fixture(autouse=True)
def _keep_unit_tests_off_a_real_pymol(request, monkeypatch):
    """Refuse every connection a unit test makes to a real PyMOL's port.

    Unit tests talk to fakes, but nothing stopped one that slipped past its
    mock from reaching whatever PyMOL was listening, and the default port is
    the session the developer has open. That happened while this guard was
    being written: a mutation-testing harness ran old function bodies outside
    the module the tests patch, and their requests ran in a developer's
    session — ``hide everything``, ``show sticks``, a 22 s ray trace written
    into the checkout. Here the request fails as if no PyMOL were running;
    the test itself fails only if the code checks the result. (Some view tests
    leave mcpymol.style.set_background unpatched, and it ignores its result:
    this guard is what keeps those off a real session.) Tests marked ``live``
    are exempt: reaching PyMOL is their whole point.
    """
    if request.node.get_closest_marker("live"):
        return
    blocked = _pymol_ports()
    real_connect = socket.socket.connect

    def connect(sock, address):
        if isinstance(address, tuple) and len(address) >= 2 and address[1] in blocked:
            raise ConnectionRefusedError(
                f"unit test tried to reach a real PyMOL at {address[0]}:{address[1]}"
            )
        return real_connect(sock, address)

    monkeypatch.setattr(socket.socket, "connect", connect)


@pytest.fixture(scope="session")
def _pymol_working_in_a_temp_dir(tmp_path_factory):
    """Point the live PyMOL's working directory at a temp dir for the run.

    PyMOL writes fetched .cif files and mpng frames into its own working
    directory: wherever it was launched, which is usually the checkout. Live
    runs left eleven .cif files and _live_sym*.png frames there. A setting
    such as fetch_path would not hold, because the suites reinitialize all
    the time; the process working directory survives that. Put back after.
    """
    from mcpymol.bridge import send_request

    res = send_request("do", args=["pwd"], timeout=10)
    if res.get("status") != "success":
        yield None  # no PyMOL listening; the live tests skip themselves
        return
    lines = str(res.get("output") or "").strip().splitlines()
    original = lines[-1].strip() if lines else None
    work = tmp_path_factory.mktemp("pymol_cwd")
    send_request("do", args=[f"cd {work}"], timeout=10)
    yield work
    if original:
        send_request("do", args=[f"cd {original}"], timeout=10)


@pytest.fixture(autouse=True)
def _live_tests_leave_the_checkout_alone(request):
    if request.node.get_closest_marker("live"):
        request.getfixturevalue("_pymol_working_in_a_temp_dir")
