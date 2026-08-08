"""Tests for the command-line entry points.

The plugin half of MCPymol runs inside PyMOL, whose Python is a separate
interpreter that cannot import this package — so it has to be loaded from a
file path, and that path lives somewhere unhelpful inside a pip or uvx
install. These commands exist to find and register it.
"""

import pathlib

import pytest

from mcpymol.cli import (
    BLOCK_END,
    BLOCK_START,
    build_parser,
    install_plugin,
    plugin_path,
    run,
    uninstall_plugin,
)


def test_plugin_path_points_at_a_real_file():
    """It is handed to PyMOL as `run <path>`, so it has to exist."""
    path = plugin_path()

    assert path.is_file()
    assert path.name == "plugin.py"
    assert path.is_absolute()
    assert "PyMOLSocketServer" in path.read_text()


# ── install ──────────────────────────────────────────────────────────────────


def test_install_writes_a_runnable_line(tmp_path):
    rc = tmp_path / ".pymolrc.py"

    message = install_plugin(rc)

    body = rc.read_text()
    assert f'cmd.do("run {plugin_path()}")' in body
    assert "from pymol import cmd" in body
    assert "Added" in message


def test_install_creates_the_file_when_absent(tmp_path):
    rc = tmp_path / "nested" / ".pymolrc.py"

    install_plugin(rc)

    assert rc.is_file()


def test_install_preserves_existing_settings(tmp_path):
    """A pymolrc is somebody's accumulated configuration. Appending to it must
    not disturb what is already there."""
    rc = tmp_path / ".pymolrc.py"
    rc.write_text("set ray_trace_mode, 1\nset ambient, 0.3\n")

    install_plugin(rc)

    body = rc.read_text()
    assert "set ray_trace_mode, 1" in body
    assert "set ambient, 0.3" in body
    assert BLOCK_START in body


def test_install_twice_updates_rather_than_duplicates(tmp_path):
    """The embedded path points into one installation and goes stale on
    upgrade, so re-running is the fix — and must not stack up blocks."""
    rc = tmp_path / ".pymolrc.py"

    install_plugin(rc)
    message = install_plugin(rc)

    body = rc.read_text()
    assert body.count(BLOCK_START) == 1
    assert body.count(BLOCK_END) == 1
    assert "Updated" in message


def test_install_updates_a_stale_path(tmp_path):
    """Simulates an upgrade: the old block names a path that no longer exists."""
    rc = tmp_path / ".pymolrc.py"
    rc.write_text(
        f"set ambient, 0.3\n{BLOCK_START}\n"
        f'cmd.do("run /old/version/mcpymol/plugin.py")\n{BLOCK_END}\n'
    )

    install_plugin(rc)

    body = rc.read_text()
    assert "/old/version/mcpymol/plugin.py" not in body
    assert str(plugin_path()) in body
    assert "set ambient, 0.3" in body


def test_install_says_how_to_confirm_it_worked(tmp_path):
    message = install_plugin(tmp_path / ".pymolrc.py")

    assert "listening on 127.0.0.1:9876" in message
    assert "upgrading" in message


# ── uninstall ────────────────────────────────────────────────────────────────


def test_uninstall_removes_only_the_managed_block(tmp_path):
    rc = tmp_path / ".pymolrc.py"
    rc.write_text("set ambient, 0.3\n")
    install_plugin(rc)

    message = uninstall_plugin(rc)

    body = rc.read_text()
    assert BLOCK_START not in body
    assert "set ambient, 0.3" in body
    assert "untouched" in message


def test_uninstall_is_harmless_when_nothing_is_installed(tmp_path):
    rc = tmp_path / ".pymolrc.py"
    rc.write_text("set ambient, 0.3\n")

    assert "Nothing to remove" in uninstall_plugin(rc)
    assert rc.read_text() == "set ambient, 0.3\n"


def test_uninstall_is_harmless_when_the_file_is_absent(tmp_path):
    assert "does not exist" in uninstall_plugin(tmp_path / "nope.py")


def test_install_and_uninstall_round_trip(tmp_path):
    """The file should come back to what it was."""
    rc = tmp_path / ".pymolrc.py"
    original = "set ray_trace_mode, 1\nset ambient, 0.3\n"
    rc.write_text(original)

    install_plugin(rc)
    uninstall_plugin(rc)

    assert rc.read_text() == original


# ── dispatch ─────────────────────────────────────────────────────────────────


def test_plugin_path_flag_prints_the_path(tmp_path, capsys):
    assert run(["--plugin-path"]) == 0
    assert capsys.readouterr().out.strip() == str(plugin_path())


def test_install_flag_writes_the_file(tmp_path, capsys):
    rc = tmp_path / ".pymolrc.py"

    assert run(["--install-plugin", "--pymolrc", str(rc)]) == 0

    assert BLOCK_START in rc.read_text()
    assert "Added" in capsys.readouterr().out


def test_uninstall_flag_removes_it(tmp_path, capsys):
    rc = tmp_path / ".pymolrc.py"
    run(["--install-plugin", "--pymolrc", str(rc)])

    assert run(["--uninstall-plugin", "--pymolrc", str(rc)]) == 0
    assert BLOCK_START not in rc.read_text()


def test_the_management_flags_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["--install-plugin", "--uninstall-plugin"])


def test_bare_invocation_serves_rather_than_managing(monkeypatch):
    """An MCP client launches `mcpymol` with no arguments; that must still
    start the server and not be swallowed by argument parsing."""
    served = []
    import mcpymol.app

    monkeypatch.setattr(mcpymol.app.mcp, "run", lambda *a, **k: served.append(True))

    assert run([]) == 0
    assert served == [True]


def test_install_reports_a_missing_plugin(tmp_path, monkeypatch):
    monkeypatch.setattr("mcpymol.cli.plugin_path", lambda: pathlib.Path("/nowhere/plugin.py"))

    message = install_plugin(tmp_path / ".pymolrc.py")

    assert message.startswith("Error")
    assert run(["--install-plugin", "--pymolrc", str(tmp_path / "x.py")]) == 1
