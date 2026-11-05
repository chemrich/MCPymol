"""The shared MCP application object.

Kept in its own module so every tool module can register against the same
instance without importing each other — the tool modules import ``mcp`` from
here, and ``mcpymol.server`` imports the tool modules.
"""

import functools
import threading
from collections.abc import Callable
from importlib.metadata import PackageNotFoundError, version
from typing import Any, TypeVar

from mcp import MCPError
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

_F = TypeVar("_F", bound=Callable[..., Any])

try:
    _VERSION = version("mcpymol")
except PackageNotFoundError:  # a source tree that was never installed
    _VERSION = ""

# Synchronizes tool execution across mcp 2.x worker threads so multi-step
# operations on the single shared PyMOL session cannot interleave.
SESSION_LOCK = threading.RLock()

# Provably read-only tools: these only query PyMOL or read local headers and
# never mutate objects, representations, settings, coordinates, or files.
# Note: `atom_properties` is intentionally excluded because `cmd.iterate`
# evaluates its expression as Python per atom; `interface_report` is excluded
# because it creates/deletes temporary objects during SASA calculation.
READ_ONLY_TOOLS: frozenset[str] = frozenset(
    {
        "contact_report",
        "count_atoms",
        "get_sequence",
        "help",
        "intra_rms",
        "list_chains",
        "list_ligands",
        "list_objects",
        "ls",
        "map_info",
        "pwd",
        "rms_cur",
        "sasa",
        "structure_info",
    }
)

# Tools that delete/replace session objects, mutate coordinates or atom state,
# write files to disk, or execute arbitrary commands.
DESTRUCTIVE_TOOLS: frozenset[str] = frozenset(
    {
        "align",
        "alter",
        "alter_state",
        "atom_properties",
        "bond",
        "delete",
        "execute_pymol_command",
        "fetch_alphafold",
        "fetch_structure",
        "h_add",
        "h_fill",
        "intra_fit",
        "load_ensemble",
        "load_map",
        "load_session",
        "load_structure",
        "morph_states",
        "move",
        "mpng",
        "png",
        "print_export",
        "remove",
        "render",
        "reset",
        "restore_bfactors",
        "save",
        "save_session",
        "sculpt_iterate",
        "set_symmetry",
        "super",
        "system",
        "turntable",
        "unbond",
    }
)

# Tools that reach external network APIs or arbitrary OS/Python execution.
OPEN_WORLD_TOOLS: frozenset[str] = frozenset(
    {
        "conservation_view",
        "execute_pymol_command",
        "fetch_alphafold",
        "fetch_structure",
        "structure_info",
        "system",
    }
)

# Tools whose effect is relative/incremental rather than idempotent.
NON_IDEMPOTENT_TOOLS: frozenset[str] = frozenset(
    {
        "alter",
        "alter_state",
        "atom_properties",
        "backward",
        "execute_pymol_command",
        "fab",
        "forward",
        "fragment",
        "move",
        "sculpt_iterate",
        "symexp",
        "system",
        "turn",
    }
)


def _tool_annotations(tool_name: str) -> ToolAnnotations:
    """Return conservative MCP ToolAnnotations for ``tool_name``."""
    read_only = tool_name in READ_ONLY_TOOLS
    return ToolAnnotations(
        read_only_hint=read_only,
        destructive_hint=False if read_only else (tool_name in DESTRUCTIVE_TOOLS),
        idempotent_hint=True if read_only else (tool_name not in NON_IDEMPOTENT_TOOLS),
        open_world_hint=tool_name in OPEN_WORLD_TOOLS,
    )


class _MCPymolServer(MCPServer):
    def tool(self, *args: Any, **kwargs: Any) -> Callable[[_F], _F]:
        def decorator(fn: _F) -> _F:
            tool_name = kwargs.get("name") or (
                args[0] if args and isinstance(args[0], str) else fn.__name__
            )
            call_kwargs = dict(kwargs)
            if "annotations" not in call_kwargs:
                call_kwargs["annotations"] = _tool_annotations(tool_name)
            register = super(_MCPymolServer, self).tool(*args, **call_kwargs)

            @functools.wraps(fn)
            def serialized(*a: Any, **kw: Any) -> Any:
                with SESSION_LOCK:
                    try:
                        return fn(*a, **kw)
                    except (ToolError, MCPError):
                        raise
                    except Exception as exc:
                        raise ToolError(str(exc)) from exc

            register(serialized)
            return fn

        return decorator


mcp = _MCPymolServer("MCPymol", version=_VERSION)
