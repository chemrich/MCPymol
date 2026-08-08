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

# How long a single connection may take to deliver its request.  The accept
# loop is deliberately single-threaded (``pymol.cmd`` is not thread-safe, so
# serializing requests is correct), which means a peer that connects and then
# neither completes a JSON document nor half-closes would otherwise wedge the
# bridge for the rest of the PyMOL session.  Requests come from localhost and
# are a few KB at most, so this is generous.
RECV_TIMEOUT = float(os.environ.get("MCPYMOL_RECV_TIMEOUT", 30.0))

# Upper bound on one request.  Real requests are small; the largest in normal
# use is conservation_view's batched `do` script (tens of KB for a long chain).
MAX_REQUEST_BYTES = int(os.environ.get("MCPYMOL_MAX_REQUEST_BYTES", 8 * 1024 * 1024))


# How often the accept loop checks whether it has been asked to stop.
_ACCEPT_POLL_INTERVAL = 0.25


class RequestTooLarge(Exception):
    """Raised when a peer sends more than ``MAX_REQUEST_BYTES`` in one request."""


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


def _recv_all(sock: socket.socket, max_bytes: int | None = None) -> bytes:
    """Read one request from ``sock``.

    Returns as soon as the accumulated bytes parse as a JSON document, or when
    the peer half-closes.  The bridge sends one request per connection and then
    half-closes after sendall(), so either exit fires promptly — but a peer that
    does neither must not be able to block the accept loop forever, so the
    caller sets a socket timeout and this raises ``TimeoutError`` in that case.

    Args:
        sock: Connected socket to read one request from.
        max_bytes: Size cap; defaults to the current ``MAX_REQUEST_BYTES``.
            Resolved at call time rather than bound as a default argument, so
            the module constant stays genuinely overridable.

    Raises:
        RequestTooLarge: the peer sent more than ``max_bytes``.
        TimeoutError: the socket timeout elapsed mid-request.
    """
    limit = MAX_REQUEST_BYTES if max_bytes is None else max_bytes
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = sock.recv(RECV_CHUNK)
        if not chunk:
            break
        chunks.append(chunk)
        total += len(chunk)
        if total > limit:
            raise RequestTooLarge(f"request exceeded {limit} bytes")
        # JSON is self-delimiting once we have the whole document.  Checking
        # after every chunk means a well-formed request is served immediately
        # even if the peer never half-closes (and keeps mock sockets simple).
        try:
            json.loads(b"".join(chunks).decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        break
    return b"".join(chunks)


def _sendall(sock: socket.socket, payload: bytes) -> None:
    """Send the whole response, then half-close the write side so the client
    knows the message is complete (mirrors the framing the bridge uses).

    Best-effort: a peer that has already gone away is not worth raising over,
    since there is nobody left to receive the error.
    """
    try:
        sock.sendall(payload)
    except OSError:
        return
    try:
        sock.shutdown(socket.SHUT_WR)
    except OSError:
        pass


def _error_bytes(message: str) -> bytes:
    """Encode an error response in the bridge's wire format."""
    return json.dumps({"status": "error", "error": message}).encode("utf-8")


def _jsonable(value):
    """Best-effort conversion of a PyMOL return value into JSON-able data.

    Several ``pymol.cmd`` functions return numpy arrays — ``get_coords`` and
    ``get_atom_coords`` among them.  Without this they fall through to the
    repr fallback below and reach the client as a *string* that merely looks
    like data: a ``success`` response carrying something no caller can use.
    Converting here keeps the repr fallback for genuinely opaque objects
    (chempy models, handles) where a repr really is the best available.
    """
    if value is None or isinstance(value, (str, bytes, bool, int, float)):
        return value
    tolist = getattr(value, "tolist", None)  # numpy arrays and scalars
    if callable(tolist):
        try:
            return tolist()
        except Exception:
            return value
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(v) for v in value]
    return value


def _dumps_response(response: dict) -> bytes:
    """Serialize a response, degrading gracefully rather than failing to reply.

    ``pymol.cmd`` functions are free to return objects json knows nothing about
    (chempy models, numpy arrays, opaque handles).  Letting ``json.dumps``
    raise would send the client nothing at all, which it can only report as an
    empty response — so coerce what is coercible, then fall back to a repr of
    the result, and only then to an explicit error.
    """
    try:
        return json.dumps(response).encode("utf-8")
    except (TypeError, ValueError):
        pass

    coerced = dict(response)
    try:
        coerced["result"] = _jsonable(coerced.get("result"))
        return json.dumps(coerced).encode("utf-8")
    except (TypeError, ValueError):
        pass

    coerced["result"] = repr(response.get("result"))
    try:
        return json.dumps(coerced).encode("utf-8")
    except (TypeError, ValueError):
        return _error_bytes("PyMOL returned a value that could not be serialized to JSON.")


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

        if not isinstance(payload, dict):
            return {
                "status": "error",
                "error": f"Request must be a JSON object, got {type(payload).__name__}",
            }

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
                return {
                    "status": "success",
                    "result": f"Fetched {args[0] if args else 'structure'}",
                }

            if action == "load":
                cmd.load(*args, **kwargs)
                return {"status": "success", "result": f"Loaded {args[0] if args else 'structure'}"}

            if action == "get_chains":
                selection = args[0] if args else kwargs.get("selection", "all")
                return {"status": "success", "result": cmd.get_chains(selection)}

            if action == "iterate_to_list":
                # Read per-atom properties back to the caller.  Nothing else in
                # the protocol can: the properties that live on atoms rather
                # than on objects — occupancy, altloc, per-atom b — are only
                # reachable through cmd.iterate or cmd.get_model, and get_model
                # returns a chempy object that does not survive the JSON wire.
                #
                # Always returns a list of tuples, one per atom, with one
                # element per comma-separated field in the expression, so a
                # single-field query is a list of 1-tuples rather than a list
                # of scalars.  Callers can then unpack without special-casing.
                selection = args[0] if args else kwargs.get("selection", "all")
                expression = args[1] if len(args) > 1 else kwargs.get("expression")
                if not expression:
                    return {
                        "status": "error",
                        "error": (
                            "iterate_to_list needs an expression, "
                            "e.g. iterate_to_list('polymer', 'chain, resi, q')"
                        ),
                    }
                space: dict = {"_rows": []}
                cmd.iterate(selection, f"_rows.append(({expression},))", space=space)
                return {"status": "success", "result": space["_rows"]}

            # ── Dynamic resolution (incl. dotted names like util.cbc) ───
            method = _resolve_dotted(action) if action else None
            if method is not None:
                rv = method(*args, **kwargs)
                return {
                    "status": "success",
                    "result": rv if rv is not None else f"Executed '{action}' successfully.",
                }

            return {
                "status": "error",
                "error": f"Unknown action or method not found on cmd: {action}",
            }

        except Exception as e:
            err_msg = str(e)
            print(f"MCPymol Plugin Error executing {action}: {err_msg}")
            traceback.print_exc()
            return {"status": "error", "error": f"PyMOL execution error: {err_msg}"}

    # ── Per-connection handling ─────────────────────────────────────────
    def serve_connection(self, client: socket.socket) -> None:
        """Read one request off ``client``, run it, and write the response.

        Every failure path still writes something back, so the bridge sees a
        real error instead of an empty response it has to guess about.
        """
        try:
            client.settimeout(RECV_TIMEOUT)
        except OSError:
            pass

        try:
            data = _recv_all(client)
        except RequestTooLarge as e:
            _sendall(client, _error_bytes(f"Request rejected: {e}."))
            return
        except TimeoutError:
            # Must precede OSError: TimeoutError is a subclass of it.
            _sendall(
                client,
                _error_bytes(
                    f"Timed out after {RECV_TIMEOUT}s waiting for a complete request. "
                    "The client must send one JSON object and then half-close."
                ),
            )
            return
        except OSError as e:
            print(f"MCPymol Socket Server: read failed: {e}")
            return

        if not data:
            return

        try:
            payload = json.loads(data.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            _sendall(client, _error_bytes(f"Invalid JSON payload: {e}"))
            return

        response = self.handle_request(payload)
        _sendall(client, _dumps_response(response))

    # ── Socket lifecycle ────────────────────────────────────────────────
    def _listen_loop(self) -> None:
        self.server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

        try:
            self.server_socket.bind((self.host, self.port))
            self.server_socket.listen(5)
            # Bounds how long stop() waits: the loop only notices `running` went
            # false when accept() returns. Closing the socket from the stopping
            # thread would be faster but does not reliably wake accept() on
            # Linux, and risks an fd-reuse race — a short poll is simpler and
            # costs a few no-op wakeups per second on an idle session.
            self.server_socket.settimeout(_ACCEPT_POLL_INTERVAL)
            # Report the port actually bound, which differs from self.port when
            # 0 was requested (ephemeral).
            self.port = self.server_socket.getsockname()[1]
            print(f"MCPymol Native Plugin listening on {self.host}:{self.port}")

            while self.running:
                try:
                    client, _addr = self.server_socket.accept()
                except TimeoutError:
                    continue
                except OSError:  # socket closed during shutdown
                    break

                try:
                    self.serve_connection(client)
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


def start_plugin() -> PyMOLSocketServer | None:
    """Start the bridge singleton and register the PyMOL commands.

    Idempotent across ``run plugin.py`` re-runs: a server left over from a
    previous run is stopped before the new one binds.
    """
    global mcp_bridge_plugin
    try:
        existing = globals().get("mcp_bridge_plugin")
        if existing is not None:
            existing.stop()
        mcp_bridge_plugin = PyMOLSocketServer()
        mcp_bridge_plugin.start()
        if cmd is not None:
            cmd.extend("stop_mcp", mcp_bridge_plugin.stop)
            cmd.extend("start_mcp", mcp_bridge_plugin.start)
        return mcp_bridge_plugin
    except Exception as e:
        print(f"Failed to initialize MCPymol plugin: {e}")
        return None


# ── Auto-start ──────────────────────────────────────────────────────────
#
# PyMOL's `run plugin.py` execs this file as a script, so __name__ is whatever
# namespace PyMOL ran it in — never "mcpymol.plugin".  Importing it as a
# library module (the test suite, tooling, anything doing
# `import mcpymol.plugin`) must NOT bind a port: doing so races a real PyMOL
# session that already owns 9876.  That distinction is exactly __name__, so
# gate on it rather than on an opt-in flag nobody would remember to set.
# NB: only seed the global when it is absent.  PyMOL's `run` re-executes this
# file in a namespace that persists across runs, and the previous run's server
# must survive long enough for start_plugin() to stop it and free the port.
if "mcp_bridge_plugin" not in globals():
    mcp_bridge_plugin = None

if __name__ != "mcpymol.plugin":
    start_plugin()
