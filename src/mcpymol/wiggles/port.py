"""The entire coupling surface between Wiggles and a running PyMOL.

SPEC.md §0: Wiggles talks to PyMOL through one small protocol and nothing else.
Everything above this module is pure and testable with :class:`FakePort` and no
PyMOL installed.

Two operations:

``do(command)``
    Run a PyMOL command line. Fire-and-forget; returns a status string.

``query(action, *args, **kwargs)``
    Call something on the PyMOL side and get a JSON-able value back. The shape
    is deliberately **structured rather than a string to be parsed** — it maps
    one-to-one onto the wire protocol's ``{action, args, kwargs}``, so no
    adapter has to reverse-engineer an expression.

Two adapters implement that protocol, and inside MCPymol only one of them is
used. :class:`SendRequestPort` wraps :func:`mcpymol.bridge.send_request`, which
already speaks the wire — reimplementing it here would be two sources of truth
for one socket.

:class:`BridgePort` is the other, and it stays. It speaks the protocol directly
over the standard library without importing ``mcpymol``, which is what makes
this package runnable on its own, and it is what the Wiggles live-fire harness
drives. :func:`interpret` is shared by both, so a response translates
identically whichever one is in play — asserted by a test, because if they
diverge the behaviour depends on how the code happens to be installed.

Per-atom queries need the ``iterate_to_list`` plugin action. It is on
``mcpymol.plugin`` as of PR #41 and is load-bearing for MCPymol's own
``atom_properties`` tool as well as for everything here — see
:data:`ITERATE_TO_LIST`.
"""

from __future__ import annotations

import json
import os
import socket
from typing import Any, Protocol, runtime_checkable

# Loopback, with no environment override — deliberately matching
# ``mcpymol.bridge``. Standalone, this adapter read a host-override variable
# that MCPymol itself does not honour, which would have been a new public knob
# arriving by accident. The constructor still takes ``host``, which is what the
# live-fire harness passes.
HOST = "127.0.0.1"
PORT = int(os.environ.get("MCPYMOL_PORT", 9876))

DEFAULT_TIMEOUT = float(os.environ.get("MCPYMOL_WIGGLES_TIMEOUT", 10.0))
_RECV_CHUNK = 65536

#: Reading per-atom properties (occupancy, altloc, per-atom b) has no other
#: route: they live on atoms, and ``cmd.get_model`` returns a chempy object
#: that does not survive JSON. Added for these tools in PR #41 and now used by
#: MCPymol's own ``atom_properties`` as well.
ITERATE_TO_LIST = "iterate_to_list"

#: Actions whose *return value* Wiggles actually consumes. A test that forgets
#: to stub one of these has a hole in it, so FakePort raises. Everything else
#: is a command whose result is ignored, and stubbing each one would be noise
#: that hides the queries that matter.
DATA_ACTIONS = frozenset({ITERATE_TO_LIST, "count_states", "get_coords", "get_names"})


class PortError(RuntimeError):
    """A command or query failed on the PyMOL side, or never got there."""


def call(port: PymolPort, action: str, *args: Any, **kwargs: Any) -> Any:
    """Invoke a ``pymol.cmd`` function, raising if PyMOL rejects it.

    **Prefer this over** :meth:`PymolPort.do`. ``cmd.do`` runs a command *line*
    and reports success for any string at all — verified against a live
    session, it returns ``"Executed command: ..."`` for
    ``this_is_not_a_pymol_command foo, bar``. Errors go to PyMOL's log, not the
    return value, so a command line sent that way is unverifiable: a live test
    would pass whether the string was right or wrong.

    A structured action makes the plugin call the function directly, so PyMOL's
    own exception comes back:

        PortError: isomesh failed: PyMOL execution error: Error: Map object
        "no_such_map" not found.

    ``do`` remains correct for the few things that are not a ``cmd`` call —
    assigning into PyMOL's ``stored`` namespace, for instance.
    """
    return port.query(action, *args, **kwargs)


@runtime_checkable
class PymolPort(Protocol):
    """What Wiggles needs from a PyMOL session. Nothing more."""

    def do(self, command: str) -> str:
        """Run a PyMOL command line. Returns a status string."""
        ...

    def query(self, action: str, *args: Any, **kwargs: Any) -> Any:
        """Call ``action`` on the PyMOL side, returning a JSON-able value."""
        ...


class FakePort:
    """Recording port for tests: no PyMOL, no sockets, no I/O.

    Commands land in :attr:`commands` in order. Queries are answered from
    ``responses``, keyed by **action name**. A response may be a plain value or
    a callable, which is invoked with the query's arguments — that is how a
    test varies an answer by state without threading state through the fixture.

    An unanticipated query raises :class:`KeyError` rather than returning
    ``None``, so a view asking for something the test did not set up fails
    loudly instead of silently rendering nothing.
    """

    def __init__(self, responses: dict[str, Any] | None = None) -> None:
        self.commands: list[str] = []
        self.queries: list[tuple[str, tuple, dict]] = []
        self.responses: dict[str, Any] = dict(responses or {})

    def do(self, command: str) -> str:
        self.commands.append(command)
        return "OK"

    def query(self, action: str, *args: Any, **kwargs: Any) -> Any:
        self.queries.append((action, args, kwargs))
        if action in self.responses:
            response = self.responses[action]
            return response(*args, **kwargs) if callable(response) else response
        if action in DATA_ACTIONS:
            raise KeyError(
                f"FakePort has no response for data action {action!r}, whose "
                f"return value the caller uses. Known: {sorted(self.responses)}"
            )
        return "OK"  # a command; its result is not consumed

    # -- test conveniences -------------------------------------------------

    def commands_matching(self, needle: str) -> list[str]:
        """Every recorded command containing ``needle``."""
        return [c for c in self.commands if needle in c]

    def ran(self, needle: str) -> bool:
        """Did any recorded command contain ``needle``?"""
        return any(needle in c for c in self.commands)

    def queried(self, action: str) -> bool:
        """Was ``action`` asked for?"""
        return any(a == action for a, _, _ in self.queries)

    def calls(self, action: str) -> list[tuple[tuple, dict]]:
        """Every (args, kwargs) this port saw for ``action``."""
        return [(a, k) for name, a, k in self.queries if name == action]

    def called(self, action: str, *args: Any, **kwargs: Any) -> bool:
        """Was ``action`` invoked with exactly these arguments?"""
        return (tuple(args), dict(kwargs)) in self.calls(action)

    @property
    def call_log(self) -> str:
        """Every call, one per line — for assertion messages."""
        return "\n".join(
            f"{name}({', '.join([*map(repr, a), *[f'{k}={v!r}' for k, v in kw.items()]])})"
            for name, a, kw in self.queries
        )

    @property
    def transcript(self) -> str:
        """All commands, newline-joined — useful in assertion messages."""
        return "\n".join(self.commands)


def interpret(action: str, response: Any) -> Any:
    """Turn one plugin response into a value, or raise.

    Shared by every adapter. MCPymol's ``send_request`` returns an error
    *dict*; Wiggles' views are written against exceptions, because an error
    dict that looks like a result is the failure mode this project keeps
    finding. This is the one place that translation happens.
    """
    if not isinstance(response, dict):
        raise PortError(f"expected a JSON object for {action!r}, got {response!r}")

    if response.get("status") != "success":
        error = response.get("error", "no error message given")
        if action == ITERATE_TO_LIST and "unknown action" in str(error).lower():
            raise PortError(
                f"the plugin does not support {ITERATE_TO_LIST!r}. It reads "
                f"per-atom properties (occupancy, altloc, per-atom b) and "
                f"nothing else can. The running plugin predates PR #41 — "
                f"reload it from this install."
            )
        raise PortError(f"{action} failed: {error}")

    return response.get("result")


class SendRequestPort:
    """Adapter over a ``send_request``-shaped callable.

    This is what Wiggles uses when it lives *inside* MCPymol as
    ``mcpymol.wiggles``: the socket protocol is already implemented there, and
    reimplementing it in-process would be two sources of truth for one wire.

        from mcpymol.bridge import send_request
        port = SendRequestPort(send_request)

    The callable must accept ``(action, args=..., kwargs=...)`` and return the
    plugin's response dict. Anything matching that shape works, which is also
    what makes this trivial to fake.
    """

    def __init__(self, send_request: Any, timeout: float | None = None) -> None:
        self._send_request = send_request
        self._timeout = timeout

    def _call(self, action: str, args: list, kwargs: dict) -> Any:
        extra = {} if self._timeout is None else {"timeout": self._timeout}
        try:
            response = self._send_request(action, args=args, kwargs=kwargs, **extra)
        except Exception as exc:
            raise PortError(f"{action} could not be sent: {type(exc).__name__}: {exc}") from exc
        return interpret(action, response)

    def do(self, command: str) -> str:
        result = self._call("do", [command], {})
        return str(result) if result is not None else "OK"

    def query(self, action: str, *args: Any, **kwargs: Any) -> Any:
        return self._call(action, list(args), dict(kwargs))


class BridgePort:
    """Talks to the in-PyMOL plugin over MCPymol's socket protocol.

    One request per connection: write a JSON object, half-close the write side
    so the plugin sees EOF, then read until the plugin half-closes. Framing is
    MCPymol's; this is a reimplementation over the standard library rather than
    an import, so Wiggles stays installable on its own.

    Unlike MCPymol's ``send_request``, which returns an error *dict*, every
    failure here raises :class:`PortError`. Wiggles' views are written against
    exceptions, and an error dict that looks like a result is the failure mode
    this project keeps finding in other people's code.
    """

    def __init__(
        self,
        host: str = HOST,
        port: int = PORT,
        timeout: float = DEFAULT_TIMEOUT,
    ) -> None:
        self.host = host
        self.port = port
        self.timeout = timeout

    def _send(
        self,
        action: str,
        args: list | None = None,
        kwargs: dict | None = None,
        timeout: float | None = None,
    ) -> Any:
        payload = json.dumps({"action": action, "args": args or [], "kwargs": kwargs or {}}).encode(
            "utf-8"
        )

        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                sock.settimeout(self.timeout if timeout is None else timeout)
                sock.connect((self.host, self.port))
                sock.sendall(payload)
                try:
                    sock.shutdown(socket.SHUT_WR)
                except OSError:
                    pass
                raw = self._drain(sock)
        except (TimeoutError, OSError) as exc:
            raise PortError(
                f"could not reach the PyMOL plugin at {self.host}:{self.port} "
                f"({type(exc).__name__}: {exc}). Is PyMOL running with the "
                f"MCPymol plugin loaded?"
            ) from exc

        if not raw:
            raise PortError(f"empty response from the PyMOL plugin for {action!r}")

        try:
            response = json.loads(raw.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise PortError(f"malformed response for {action!r}: {exc}") from exc

        return interpret(action, response)

    @staticmethod
    def _drain(sock: socket.socket) -> bytes:
        """Read until the response parses as JSON, or the peer half-closes.

        Parsing after each chunk means a well-formed reply is taken as soon as
        it is complete, without depending on the peer to close first.
        """
        chunks: list[bytes] = []
        while True:
            chunk = sock.recv(_RECV_CHUNK)
            if not chunk:
                break
            chunks.append(chunk)
            try:
                json.loads(b"".join(chunks).decode("utf-8"))
            except (json.JSONDecodeError, UnicodeDecodeError):
                continue
            break
        return b"".join(chunks)

    # -- PymolPort ---------------------------------------------------------

    def do(self, command: str) -> str:
        result = self._send("do", args=[command])
        return str(result) if result is not None else "OK"

    def query(self, action: str, *args: Any, **kwargs: Any) -> Any:
        return self._send(action, args=list(args), kwargs=dict(kwargs))
