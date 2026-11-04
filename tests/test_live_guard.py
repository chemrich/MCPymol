"""The live-run guard in conftest: no default port by accident, no vacuous green.

Each case runs a real pytest in a subprocess against a copy of this suite's
conftest, so what is tested is the hooks themselves, not helpers they call.
Nothing here reaches a PyMOL: the stand-in live tests don't use the bridge,
and the port they are pointed at has nothing listening.
"""

import socket
from pathlib import Path

import pytest

pytest_plugins = ["pytester"]

CONFTEST = (Path(__file__).parent / "conftest.py").read_text()

LIVE_TESTS = """
import pytest

# The conftest's autouse fixture sends `cd` to the live PyMOL. These stand-ins
# never need it: override it, so nothing here ever opens a connection.
@pytest.fixture(autouse=True)
def _live_tests_leave_the_checkout_alone():
    yield

@pytest.mark.live
def test_live_passes():
    pass

@pytest.mark.live
def test_live_skips():
    pytest.skip("no PyMOL listening")

def test_unit():
    pass
"""


def _free_port() -> int:
    """A port with nothing listening on it, for the cases that must not reach 9876."""
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


@pytest.fixture
def suite(pytester, monkeypatch):
    for var in ("MCPYMOL_PORT", "MCPYMOL_TEST_PORT", "MCPYMOL_LIVE_MIN_PASSED"):
        monkeypatch.delenv(var, raising=False)
    pytester.makeconftest(CONFTEST)
    pytester.makeini(
        "[pytest]\nmarkers =\n    live: drives a running PyMOL\n    network: hits the net\n"
    )
    pytester.makepyfile(test_stand_in=LIVE_TESTS)
    return pytester


def test_a_live_run_with_no_test_port_is_refused_before_anything_runs(suite):
    result = suite.runpytest_subprocess("-m", "live")
    assert result.ret == pytest.ExitCode.USAGE_ERROR
    result.stderr.fnmatch_lines(["*would clear the session*MCPYMOL_TEST_PORT=9898 pytest -m live*"])
    result.stdout.no_fnmatch_line("*passed*")


def test_the_shells_mcpymol_port_is_not_a_test_port(suite, monkeypatch):
    # The README tells people with a port clash to move their everyday PyMOL
    # with MCPYMOL_PORT, so a port in the shell can be a real session.
    monkeypatch.setenv("MCPYMOL_PORT", "9867")
    assert suite.runpytest_subprocess("-m", "live").ret == pytest.ExitCode.USAGE_ERROR


def test_a_run_that_deselects_live_tests_is_untouched(suite):
    # The project's default addopts deselect live tests; nothing should change.
    result = suite.runpytest_subprocess("-m", "not live")
    assert result.ret == pytest.ExitCode.OK
    result.assert_outcomes(passed=1, deselected=2)


def test_collect_only_is_let_through(suite):
    result = suite.runpytest_subprocess("-m", "live", "--collect-only", "-q")
    assert result.ret == pytest.ExitCode.OK
    result.stdout.no_fnmatch_line("*only 0 of*")


def test_a_test_port_runs_normally(suite, monkeypatch):
    monkeypatch.setenv("MCPYMOL_TEST_PORT", str(_free_port()))
    result = suite.runpytest_subprocess("-m", "live")
    assert result.ret == pytest.ExitCode.OK
    result.assert_outcomes(passed=1, skipped=1, deselected=1)


def test_the_test_port_is_what_the_bridge_reaches(suite, monkeypatch):
    port = _free_port()
    monkeypatch.setenv("MCPYMOL_PORT", "9867")
    monkeypatch.setenv("MCPYMOL_TEST_PORT", str(port))
    suite.makepyfile(
        test_bridge_port=f"""
import pytest

@pytest.fixture(autouse=True)
def _live_tests_leave_the_checkout_alone():
    yield

@pytest.mark.live
def test_port():
    from mcpymol import bridge
    assert bridge.PORT == {port}
"""
    )
    result = suite.runpytest_subprocess("-m", "live", "test_bridge_port.py")
    assert result.ret == pytest.ExitCode.OK, result.stdout.str()


def test_a_live_run_where_every_live_test_skipped_fails(suite, monkeypatch):
    port = _free_port()
    monkeypatch.setenv("MCPYMOL_TEST_PORT", str(port))
    result = suite.runpytest_subprocess("-m", "live", "-k", "skips")
    assert result.ret == pytest.ExitCode.TESTS_FAILED
    result.stdout.fnmatch_lines([f"*only 0 of 1 live tests passed*port {port}*"])


def test_the_minimum_pass_count_is_enforced(suite, monkeypatch):
    monkeypatch.setenv("MCPYMOL_TEST_PORT", str(_free_port()))
    monkeypatch.setenv("MCPYMOL_LIVE_MIN_PASSED", "2")
    result = suite.runpytest_subprocess("-m", "live")
    assert result.ret == pytest.ExitCode.TESTS_FAILED
    result.stdout.fnmatch_lines(["*only 1 of 2 live tests passed; at least 2 must*"])


def test_a_failing_run_keeps_its_own_exit_code_and_story(suite, monkeypatch):
    monkeypatch.setenv("MCPYMOL_TEST_PORT", str(_free_port()))
    suite.makepyfile(test_broken="import pytest\n\nraise RuntimeError('broken module')\n")
    result = suite.runpytest_subprocess("-m", "live")
    assert result.ret == pytest.ExitCode.INTERRUPTED
    result.stdout.no_fnmatch_line("*live tests passed*")


@pytest.mark.parametrize(
    ("var", "value"), [("MCPYMOL_TEST_PORT", "abc"), ("MCPYMOL_LIVE_MIN_PASSED", "all")]
)
def test_a_malformed_setting_is_a_usage_error_not_a_traceback(suite, monkeypatch, var, value):
    monkeypatch.setenv("MCPYMOL_TEST_PORT", str(_free_port()))
    monkeypatch.setenv(var, value)
    result = suite.runpytest_subprocess("-m", "live")
    assert result.ret == pytest.ExitCode.USAGE_ERROR
    result.stderr.fnmatch_lines([f"*{var}=*{value}*"])


def test_a_malformed_minimum_does_not_break_a_unit_run(suite, monkeypatch):
    monkeypatch.setenv("MCPYMOL_LIVE_MIN_PASSED", "all")
    assert suite.runpytest_subprocess("-m", "not live").ret == pytest.ExitCode.OK
