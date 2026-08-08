"""Command-line entry points.

Bare ``mcpymol`` runs the MCP server over stdio, which is what an MCP client
launches. The flags exist for the other half of the install: the plugin runs
*inside* PyMOL, whose Python is a separate interpreter that cannot import this
package, so it has to be loaded from a file path. Finding that path inside a
pip or uvx installation is unpleasant, and these do it for you.
"""

import argparse
import pathlib
import sys

# The block is delimited so it can be rewritten in place. The path embedded in
# it points into a specific installation and goes stale when the package is
# upgraded, so re-running the command has to update it rather than append a
# second copy.
BLOCK_START = "# >>> mcpymol >>>"
BLOCK_END = "# <<< mcpymol <<<"

DEFAULT_PYMOLRC = pathlib.Path.home() / ".pymolrc.py"


def plugin_path() -> pathlib.Path:
    """Absolute path to the plugin that runs inside PyMOL."""
    return (pathlib.Path(__file__).parent / "plugin.py").resolve()


def _block(path: pathlib.Path) -> str:
    return (
        f"{BLOCK_START}\n"
        f"# Added by `mcpymol --install-plugin`.\n"
        f"# The path points into one installation of mcpymol, so re-run that\n"
        f"# command after upgrading — this block is rewritten, not duplicated.\n"
        f"from pymol import cmd\n"
        f'cmd.do("run {path}")\n'
        f"{BLOCK_END}\n"
    )


def _without_block(text: str) -> str:
    """Strip a previously managed block, leaving everything else untouched."""
    if BLOCK_START not in text:
        return text
    before, _, rest = text.partition(BLOCK_START)
    _, _, after = rest.partition(BLOCK_END)
    return (before.rstrip("\n") + "\n" + after.lstrip("\n")).strip("\n")


def install_plugin(rc_path: pathlib.Path | None = None) -> str:
    """Point PyMOL's startup file at this installation's plugin."""
    rc = pathlib.Path(rc_path) if rc_path else DEFAULT_PYMOLRC
    plugin = plugin_path()
    if not plugin.is_file():
        return f"Error: no plugin found at {plugin}. Is the installation complete?"

    existing = rc.read_text() if rc.is_file() else ""
    had_block = BLOCK_START in existing
    kept = _without_block(existing)

    rc.parent.mkdir(parents=True, exist_ok=True)
    body = (kept + "\n\n" if kept else "") + _block(plugin)
    rc.write_text(body)

    verb = "Updated" if had_block else "Added"
    return (
        f"{verb} the MCPymol block in {rc}, pointing at:\n"
        f"  {plugin}\n\n"
        f"Restart PyMOL (or run that line now) and you should see:\n"
        f"  MCPymol Native Plugin listening on 127.0.0.1:9876\n\n"
        f"Re-run `mcpymol --install-plugin` after upgrading mcpymol; the path "
        f"changes with the version."
    )


def uninstall_plugin(rc_path: pathlib.Path | None = None) -> str:
    """Remove the managed block, leaving the rest of the file alone."""
    rc = pathlib.Path(rc_path) if rc_path else DEFAULT_PYMOLRC
    if not rc.is_file():
        return f"Nothing to remove: {rc} does not exist."
    existing = rc.read_text()
    if BLOCK_START not in existing:
        return f"Nothing to remove: no MCPymol block in {rc}."

    remaining = _without_block(existing)
    rc.write_text(remaining + "\n" if remaining else "")
    return f"Removed the MCPymol block from {rc}. Your other settings are untouched."


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mcpymol",
        description=(
            "MCP server for PyMOL. With no arguments it serves over stdio, which "
            "is how an MCP client launches it."
        ),
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--plugin-path",
        action="store_true",
        help="print the path of the plugin to load inside PyMOL, and exit",
    )
    group.add_argument(
        "--install-plugin",
        action="store_true",
        help="add the plugin to ~/.pymolrc.py so PyMOL loads it on startup",
    )
    group.add_argument(
        "--uninstall-plugin",
        action="store_true",
        help="remove the block that --install-plugin added",
    )
    parser.add_argument(
        "--pymolrc",
        metavar="PATH",
        help="operate on this file instead of ~/.pymolrc.py",
    )
    return parser


def run(argv: list[str] | None = None) -> int:
    """Dispatch the CLI. Returns a process exit code."""
    args = build_parser().parse_args(argv)
    rc = pathlib.Path(args.pymolrc) if args.pymolrc else None

    if args.plugin_path:
        print(plugin_path())
        return 0
    if args.install_plugin:
        message = install_plugin(rc)
        print(message)
        return 1 if message.startswith("Error") else 0
    if args.uninstall_plugin:
        print(uninstall_plugin(rc))
        return 0

    # No flags: serve. Imported here so the flags stay fast and do not need a
    # working MCP stack just to print a path.
    import mcpymol.server  # noqa: F401  (registers every tool)
    from mcpymol.app import mcp

    mcp.run()
    return 0


def main() -> None:
    sys.exit(run())
