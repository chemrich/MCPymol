"""
MCPymol's native PyMOL plugin.

This script runs *inside* PyMOL.  It opens a tiny localhost TCP server that
accepts JSON-encoded command requests from the external MCPymol bridge
(``mcpymol.server``), invokes the corresponding ``pymol.cmd`` (or
``pymol.util``) function, and returns the result.

Usage inside PyMOL::

    run /path/to/MCPymol/src/mcpymol/plugin.py

The listening port defaults to 9876 and can be overridden with the
``MCPYMOL_PORT`` environment variable *before* starting PyMOL.  The plugin
exposes ``start_mcp`` / ``stop_mcp`` PyMOL commands so you can toggle the
bridge at runtime.
"""
from __future__ import annotations

import json
import os
import socket
import threading
import traceback

# Read once at import so the whole session is consistent.
MCP_PORT = int(os.environ.get("MCPYMOL_PORT", 9876))

# Cap on a single recv()'d response chunk.  We loop until the peer half-closes,
# so this is just a memory upper bound on a *single* call, not a total cap.
RECV_CHUNK = 65536

try:
    from pymol import cmd
    try:
        # Optional: gives us access to ``util.cbc``, ``util.chainbow`` etc.
        from pymol import util as pymol_util
    except ImportError:  # pragma: no cover — only inside PyMOL
        pymol_util = None
except ImportError:  # pragma: no cover — only when imported outside PyMOL
    print("Warning: pymol.cmd not found. This script must be run inside PyMOL.")
    cmd = None
    pymol_util = None


def _recv_all(sock: socket.socket) -> bytes:
    """Read from ``sock`` until the peer half-closes the connection.

    The bridge sends one request per connection and then half-closes after
    sendall(), so this loop terminates promptly.
    """
    chunks: list[bytes] = []
    while True:
        chunk = sock.recv(RECV_CHUNK)
        if not chunk:
            break
        chunks.append(chunk)
    return b"".join(chunks)


def _sendall(sock: socket.socket, payload: bytes) -> None:
    """Send the whole response, then half-close the write side so the client
    knows the message is complete (mirrors the framing the bridge uses)."""
    sock.sendall(payload)
    try:
        sock.shutdown(socket.SHUT_WR)
    except OSError:
        pass


def _resolve_dotted(name: str):
    """Resolve a dotted action name like ``util.cbc`` to a callable.

    Falls back through ``pymol.util`` for ``util.*`` and through ``cmd`` for
    everything else.  Returns ``None`` if the symbol can't be found.
    """
    if cmd is None:
        return None
    if "." not in name:
        attr = getattr(cmd, name, None)
        return attr if callable(attr) else None

    head, _, tail = name.partition(".")
    if head == "util" and pymol_util is not None:
        attr = getattr(pymol_util, tail, None)
        return attr if callable(attr) else None

    # Last-ditch: walk attributes (e.g. cmd.foo.bar)
    obj = cmd
    for part in name.split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return None
    return obj if callable(obj) else None


class PyMOLSocketServer:
    """Tiny TCP/JSON server that dispatches requests to ``pymol.cmd``."""

    def __init__(self, host: str = "127.0.0.1", port: int | None = None) -> None:
        self.host = host
        self.port = port if port is not None else MCP_PORT
        self.running = False
        self.thread: threading.Thread | None = None
        self.server_socket: socket.socket | None = None

    # ── Request dispatch ────────────────────────────────────────────────
    def handle_request(self, payload):
        """Execute the requested PyMOL command and return a JSON-able dict."""
        if cmd is None:
            return {"status": "error", "error": "PyMOL cmd module not available."}

        if isinstance(payload, (str, bytes)):
            try:
                payload = json.loads(payload)
            except json.JSONDecodeError as e:
                return {"status": "error", "error": f"Invalid JSON payload: {e}"}

        action = payload.get("action")
        args = payload.get("args", []) or []
        kwargs = payload.get("kwargs", {}) or {}

        try:
            # ── Special, custom handlers ─────────────────────────────────
            if action == "do":
                command_str = args[0] if args else kwargs.get("command", "")
                cmd.do(command_str)
                return {"status": "success", "result": f"Executed command: {command_str}"}

            if action == "fetch":
                cmd.fetch(*args, **kwargs)
                return {"status": "success", "result": f"Fetched {args[0] if args else 'structure'}"}

            if action == "load":
                cmd.load(*args, **kwargs)
                return {"status": "success", "result": f"Loaded {args[0] if args else 'structure'}"}

            if action == "get_chains":
                selection = args[0] if args else kwargs.get("selection", "all")
                return {"status": "success", "result": cmd.get_chains(selection)}

            # ── Dynamic resolution (incl. dotted names like util.cbc) ───
            method = _resolve_dotted(action) if action else None
            if method is not None:
                rv = method(*args, **kwargs)
                return {
                    "status": "success",
                    "result": rv if rv is not None else f"Executed '{action}' successfully.",
                }

            return {"status": "error", "error": f"Unknown action or method not found on cmd: {action}"}

        except Exception as e:
            err_msg = str(e)
            print(f"MCPymol Plugin Error executing {action}: {err_msg}")
            traceback.print_exc()
            return {"status": "error", "error": f"PyMOL execution error: {err_msg}"}

    # ── Socket lifecycle ────────────────────────────────────────────────
    def _listen_loop(self) -> None:
        self.server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

        try:
            self.server_socket.bind((self.host, self.port))
            self.server_socket.listen(5)
            self.server_socket.settimeout(1.0)
            print(f"MCPymol Native Plugin listening on {self.host}:{self.port}")

            while self.running:
                try:
                    client, _addr = self.server_socket.accept()
                except socket.timeout:
                    continue
                except OSError:  # socket closed during shutdown
                    break

                try:
                    data = _recv_all(client)
                    if not data:
                        continue
                    try:
                        payload = json.loads(data.decode("utf-8"))
                        response = self.handle_request(payload)
                    except json.JSONDecodeError:
                        response = {"status": "error", "error": "Invalid JSON payload"}
                    _sendall(client, json.dumps(response).encode("utf-8"))
                except Exception as e:
                    print(f"MCPymol Socket Server Error: {e}")
                finally:
                    try:
                        client.close()
                    except OSError:
                        pass

        except Exception as e:
            print(f"MCPymol Failed to start socket server: {e}")
        finally:
            if self.server_socket is not None:
                try:
                    self.server_socket.close()
                except OSError:
                    pass

    def start(self) -> None:
        if self.running:
            return
        self.running = True
        self.thread = threading.Thread(target=self._listen_loop, daemon=True)
        self.thread.start()
        print("MCPymol Bridge started.")

    def stop(self) -> None:
        if not self.running:
            return
        self.running = False
        if self.thread is not None:
            self.thread.join(timeout=2.0)
        print("MCPymol Bridge stopped.")


# ── Auto-start singleton (idempotent across `run plugin.py` re-runs) ────
try:
    if "mcp_bridge_plugin" in globals():
        mcp_bridge_plugin.stop()  # type: ignore[name-defined]
    mcp_bridge_plugin = PyMOLSocketServer()
    mcp_bridge_plugin.start()
    if cmd is not None:
        cmd.extend("stop_mcp", mcp_bridge_plugin.stop)
        cmd.extend("start_mcp", mcp_bridge_plugin.start)
except Exception as e:
    print(f"Failed to initialize MCPymol plugin: {e}")
