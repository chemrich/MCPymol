"""What mcpymol.app adds on top of mcp's MCPServer.

mcp 2.x changed two things this server depends on, so both are pinned here
through the same path a client takes (``mcp.call_tool``), not by calling the
tool functions directly.
"""

import asyncio
import sys
import threading
import time
from importlib.metadata import version

import pytest
from mcp.server.mcpserver.exceptions import ToolError, UnexpectedToolError

import mcpymol.server  # noqa: F401  (registers every tool)
from mcpymol.app import mcp


def _patch_send_request(monkeypatch, fake):
    # Each tool module imported send_request by name, so patch every copy.
    for name, module in list(sys.modules.items()):
        if name.startswith("mcpymol") and hasattr(module, "send_request"):
            monkeypatch.setattr(module, "send_request", fake)


def test_parallel_tool_calls_run_one_at_a_time(monkeypatch):
    """mcp 2.x runs synchronous tools on worker threads, so two calls a client
    issues together would interleave inside the one PyMOL session. Without
    SESSION_LOCK this records two calls inside PyMOL at once."""
    inside, peak, guard = 0, 0, threading.Lock()

    def slow_send_request(*args, **kwargs):
        nonlocal inside, peak
        with guard:
            inside += 1
            peak = max(peak, inside)
        time.sleep(0.05)
        with guard:
            inside -= 1
        return {"status": "success", "result": "OK"}

    _patch_send_request(monkeypatch, slow_send_request)

    async def two_at_once():
        await asyncio.gather(mcp.call_tool("zoom", {}), mcp.call_tool("orient", {}))

    asyncio.run(two_at_once())
    assert peak == 1, f"{peak} tool calls were inside PyMOL at the same time"


def test_an_unexpected_exception_reaches_the_client_with_its_text(monkeypatch):
    """mcp 2.x reduces an unexpected exception to "Error executing tool X". The
    model needs the cause, as it had under 1.x, so the wrapper passes it on."""

    def explode(*_a, **_kw):
        raise RuntimeError("the socket layer exploded")

    _patch_send_request(monkeypatch, explode)

    with pytest.raises(ToolError) as caught:
        asyncio.run(mcp.call_tool("list_objects", {}))

    assert type(caught.value) is ToolError, "reduced to UnexpectedToolError: the text was dropped"
    assert not isinstance(caught.value, UnexpectedToolError)
    assert str(caught.value) == "Error executing tool list_objects: the socket layer exploded"


def test_the_server_reports_mcpymols_version():
    """mcp 2.x reports an empty version unless the server passes its own."""
    assert mcp.version == version("mcpymol")


def test_modules_keep_the_plain_function():
    """Only calls arriving through MCP take the lock. A tool that uses another
    tool's function directly must not wait on a lock its caller holds."""
    from mcpymol.structures import list_objects

    assert not hasattr(list_objects, "__wrapped__")


def test_every_tool_has_explicit_annotations():
    """Every registered tool carries complete, conservative MCP ToolAnnotations."""
    from mcpymol.app import (
        DESTRUCTIVE_TOOLS,
        NON_IDEMPOTENT_TOOLS,
        OPEN_WORLD_TOOLS,
        READ_ONLY_TOOLS,
    )

    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}

    # Every set in app.py must refer only to real registered tools.
    assert names >= READ_ONLY_TOOLS
    assert names >= DESTRUCTIVE_TOOLS
    assert names >= OPEN_WORLD_TOOLS
    assert names >= NON_IDEMPOTENT_TOOLS
    assert READ_ONLY_TOOLS.isdisjoint(DESTRUCTIVE_TOOLS)
    assert READ_ONLY_TOOLS.isdisjoint(NON_IDEMPOTENT_TOOLS)

    # Tools that evaluate Python or temporarily mutate session state must never
    # be marked read-only.
    for non_read_only in (
        "atom_properties",
        "interface_report",
        "execute_pymol_command",
        "distance",
        "angle",
        "dihedral",
    ):
        assert non_read_only not in READ_ONLY_TOOLS

    for tool in tools:
        ann = tool.annotations
        assert ann is not None, f"{tool.name} has no ToolAnnotations"
        assert isinstance(ann.read_only_hint, bool)
        assert isinstance(ann.destructive_hint, bool)
        assert isinstance(ann.idempotent_hint, bool)
        assert isinstance(ann.open_world_hint, bool)
        assert ann.read_only_hint is (tool.name in READ_ONLY_TOOLS)
        if ann.read_only_hint:
            assert ann.destructive_hint is False
            assert ann.idempotent_hint is True
