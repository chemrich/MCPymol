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

_F = TypeVar("_F", bound=Callable[..., Any])

try:
    _VERSION = version("mcpymol")
except PackageNotFoundError:  # a source tree that was never installed
    _VERSION = ""

# One PyMOL session sits behind every tool, and many tools are several
# requests long: print_export hides, saves and restores; the views style, then
# zoom, then render. mcp 1.x ran synchronous tools inline on the event loop,
# which serialised them by accident. mcp 2.x runs them on worker threads, so
# tool calls a client issues in parallel would interleave inside the one
# session. This lock puts the serialisation back on purpose: one tool call at
# a time, start to finish. Re-entrant, so a tool that calls another tool's
# function through MCP machinery cannot deadlock itself.
SESSION_LOCK = threading.RLock()


class _MCPymolServer(MCPServer):
    def tool(self, *args: Any, **kwargs: Any) -> Callable[[_F], _F]:
        register = super().tool(*args, **kwargs)

        def decorator(fn: _F) -> _F:
            @functools.wraps(fn)  # keeps the name, docstring and signature the schema is built from
            def serialized(*a: Any, **kw: Any) -> Any:
                with SESSION_LOCK:
                    try:
                        return fn(*a, **kw)
                    except (ToolError, MCPError):
                        raise
                    except Exception as exc:
                        # mcp 2.x keeps an unexpected exception's text on the
                        # server and tells the client only "Error executing tool
                        # X". Right for a shared HTTP server; wrong for this one,
                        # a local stdio tool whose user is its developer, where
                        # that text is how the model sees what went wrong. A
                        # ToolError is passed through, so this restores 1.x's
                        # message exactly.
                        raise ToolError(str(exc)) from exc

            register(serialized)
            # The module keeps the plain function. Only calls arriving through
            # MCP take the lock; direct calls (tests, one tool using another's
            # function) do not, so they can never wait on themselves.
            return fn

        return decorator


mcp = _MCPymolServer("MCPymol", version=_VERSION)
