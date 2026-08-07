"""The package split must not change what the outside world sees.

mcpymol.server is the documented entry point and the historical import
location for every tool, so it has to keep re-exporting the whole surface even
though the code now lives in focused modules.
"""

import importlib
import pkgutil

import pytest

import mcpymol
import mcpymol.server as server

TOOL_MODULES = [
    "mcpymol.analysis",
    "mcpymol.bridge",
    "mcpymol.comparison",
    "mcpymol.pdbtext",
    "mcpymol.conservation",
    "mcpymol.primitives",
    "mcpymol.printing",
    "mcpymol.rendering",
    "mcpymol.structures",
    "mcpymol.style",
    "mcpymol.views",
]


def _public_names(mod):
    """Top-level functions and constants a module owns (not its imports)."""
    return {
        name
        for name, obj in vars(mod).items()
        if not name.startswith("__")
        and getattr(obj, "__module__", mod.__name__) == mod.__name__
        and not isinstance(obj, type(importlib))
    }


@pytest.mark.parametrize("mod_name", TOOL_MODULES)
def test_server_reexports_every_module_name(mod_name):
    """A new tool added to a module must also appear on the facade."""
    mod = importlib.import_module(mod_name)
    missing = sorted(n for n in _public_names(mod) if not hasattr(server, n))
    assert not missing, f"{mod_name} defines {missing}, not re-exported by mcpymol.server"


@pytest.mark.parametrize("mod_name", TOOL_MODULES)
def test_reexports_are_the_same_objects(mod_name):
    """Re-export, not re-definition — patching identity has to hold."""
    mod = importlib.import_module(mod_name)
    for name in _public_names(mod):
        if hasattr(server, name):
            assert getattr(server, name) is getattr(mod, name), f"{name} diverged"


def test_all_matches_what_is_actually_exported():
    for name in server.__all__:
        assert hasattr(server, name), f"__all__ lists {name}, which is not present"


def test_every_module_registers_against_one_app():
    """Tools must land on a single FastMCP instance or half of them vanish."""
    from mcpymol.app import mcp

    assert server.mcp is mcp
    for mod_name in TOOL_MODULES:
        mod = importlib.import_module(mod_name)
        if hasattr(mod, "mcp"):
            assert mod.mcp is mcp, f"{mod_name} registers against a different app"


def test_importing_server_pulls_in_every_tool_module():
    """Registration happens on import, so a module nobody imports is a tool
    silently missing from the MCP client's list."""
    for mod_name in TOOL_MODULES:
        assert mod_name in __import__("sys").modules, f"{mod_name} not imported by mcpymol.server"


def test_no_module_was_left_behind():
    """Every module in the package is either imported by the facade or a
    deliberate standalone (the in-PyMOL plugin, the __main__ shim)."""
    standalone = {"mcpymol.plugin", "mcpymol.__main__", "mcpymol.app", "mcpymol.server"}
    found = {f"mcpymol.{m.name}" for m in pkgutil.iter_modules(mcpymol.__path__)}
    unreferenced = found - set(TOOL_MODULES) - standalone
    assert not unreferenced, f"module(s) not wired into the server: {sorted(unreferenced)}"


# ── tool schema quality ──────────────────────────────────────────────────────


def _all_tools():
    import asyncio

    from mcpymol.app import mcp

    return asyncio.run(mcp.list_tools())


def test_every_parameter_has_a_schema_description():
    """Parameter docs must live on the parameter, not only in the docstring.

    FastMCP builds the JSON schema from the signature, so an Args: block never
    reaches inputSchema.properties[].description — the model would be left
    inferring arguments from names. Descriptions come from
    Annotated[..., Field(description=...)].
    """
    undocumented = [
        f"{tool.name}.{param}"
        for tool in _all_tools()
        for param, spec in (tool.inputSchema.get("properties") or {}).items()
        if not spec.get("description")
    ]

    assert not undocumented, f"parameters missing a schema description: {undocumented}"


def test_every_tool_has_a_description():
    missing = [t.name for t in _all_tools() if not (t.description or "").strip()]

    assert not missing, f"tools with no description: {missing}"


def test_no_tool_still_carries_an_args_block():
    """Args: in a docstring is now duplicated state — the signature owns
    parameter docs, so a leftover block will drift out of sync."""
    stale = [t.name for t in _all_tools() if "Args:" in (t.description or "")]

    assert not stale, f"tools with a leftover Args: block: {stale}"


def test_descriptions_are_substantive():
    """A one-word description is not documentation."""
    thin = [
        f"{tool.name}.{param}"
        for tool in _all_tools()
        for param, spec in (tool.inputSchema.get("properties") or {}).items()
        if len(spec.get("description", "")) < 15
    ]

    assert not thin, f"parameters with a too-short description: {thin}"
