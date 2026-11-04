"""The wire to PyMOL.

Everything here is about getting a request to the in-PyMOL plugin and a
response back; nothing here knows what a protein is.  See ``mcpymol.plugin``
for the other half of the protocol.
"""

import json
import os
import socket
import sys
import time
from typing import Any

HOST = "127.0.0.1"

# Port can be overridden via environment variable, e.g.: MCPYMOL_PORT=9867 uv run mcpymol
PORT = int(os.environ.get("MCPYMOL_PORT", 9876))

# Operations that legitimately run for minutes: ray-tracing a large scene,
# writing a movie frame series, saving a big assembly.  The 10 s default is
# right for the interactive commands but guarantees a spurious timeout here.
_SLOW_OP_TIMEOUT = float(os.environ.get("MCPYMOL_SLOW_OP_TIMEOUT", 600.0))

# A single recv() chunk size.  We loop until the plugin half-closes, so this
# is just a memory hint, not a cap on the total response size.
_RECV_CHUNK = 65536

# Default socket budget for an interactive command.  Slow operations override
# it with _SLOW_OP_TIMEOUT.
_DEFAULT_TIMEOUT = 10.0

# Sleep shorter than this is not worth mentioning: clock jitter, a lid opened
# and closed.  The two renders that "stalled" for 10 minutes had the Mac
# asleep for 579 s and 716 s of their wait.
_ASLEEP_REPORT_S = 5.0


def _clock_with_sleep() -> float | None:
    """A clock that keeps counting while the computer sleeps, or None.

    ``time.monotonic()`` stops during sleep — on macOS it is
    mach_absolute_time — so a request that spanned a sleep looks exactly as
    long as the time the machine was awake, and socket timeouts are counted
    the same way.  Measured on this Mac: since boot, monotonic() read 218 h
    and CLOCK_MONOTONIC 671 h.  The difference between the two over a request
    is how long the machine slept during it.  Linux spells the sleep-counting
    clock CLOCK_BOOTTIME; Windows has none, so nothing is measured there.
    """
    if sys.platform == "darwin":
        return time.clock_gettime(time.CLOCK_MONOTONIC)
    boottime = getattr(time, "CLOCK_BOOTTIME", None)
    if boottime is not None:
        return time.clock_gettime(boottime)
    return None


def _asleep_since(awake_start: float, total_start: float | None) -> float:
    """Seconds the computer slept since the two clock readings were taken."""
    if total_start is None:
        return 0.0
    total_now = _clock_with_sleep()
    if total_now is None:
        return 0.0
    return max(0.0, (total_now - total_start) - (time.monotonic() - awake_start))


def send_request(
    action: str,
    args: list | None = None,
    kwargs: dict | None = None,
    timeout: float | None = None,
) -> dict:
    """Send a JSON request to the PyMOL plugin socket server.

    The framing is one request per TCP connection: we write the JSON payload,
    half-close our write side so the plugin sees EOF, then drain the response
    until the plugin half-closes its write side.  This avoids the 8 KB
    response truncation that bit ``get_fastastr`` on long chains.

    Args:
        action: The PyMOL command or custom action to execute.
        args: Positional arguments for the action.
        kwargs: Keyword arguments for the action.
        timeout: Socket timeout in seconds. Defaults to 10.0, except for ``do``,
            which defaults to the slow-operation timeout: the plugin answers a
            ``do`` only once the command has finished, so loading a large map
            or a long batched ``alter`` can legitimately take minutes. With
            10 s, the caller saw a timeout, moved on, and its next requests
            queued behind the unfinished command.

    Returns:
        A dict with at least a ``status`` key (``success`` or ``error``).
        ``error`` responses always include a human-readable ``error`` field.
        When the computer slept for more than a few seconds of the wait, a
        response from PyMOL also carries ``asleep_s``, and an error the bridge
        raised itself says so: PyMOL only works while the machine is awake.
    """
    if timeout is None:
        timeout = _SLOW_OP_TIMEOUT if action == "do" else _DEFAULT_TIMEOUT
    payload = {
        "action": action,
        "args": args or [],
        "kwargs": kwargs or {},
    }
    awake_start, total_start = time.monotonic(), _clock_with_sleep()
    connected = False
    chunks: list[bytes] = []
    res = None
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(timeout)
            s.connect((HOST, PORT))
            connected = True
            s.sendall(json.dumps(payload).encode("utf-8"))
            try:
                s.shutdown(socket.SHUT_WR)
            except OSError:
                pass
            while True:
                chunk = s.recv(_RECV_CHUNK)
                if not chunk:
                    break
                chunks.append(chunk)
                # JSON is self-delimiting once the whole document has arrived,
                # so a parse attempt doubles as an end-of-message test — but
                # only attempt it on a short read.
                #
                # A full-size chunk means the socket had at least that much
                # waiting, so the message is very unlikely to end there, and
                # parsing anyway costs a scan of everything received so far.
                # Doing that per chunk is quadratic: an 8 MB response spent
                # ~476 ms on parses guaranteed to fail. The plugin half-closes
                # when it is done, which is the normal exit; this check only
                # exists so a peer that does not (a test mock) still works.
                if len(chunk) < _RECV_CHUNK:
                    try:
                        res = json.loads(b"".join(chunks).decode("utf-8"))
                        break
                    except (json.JSONDecodeError, UnicodeDecodeError):
                        continue
    except ConnectionRefusedError as e:
        return {
            "status": "error",
            "error": f"Socket connection failed: nothing is listening on {HOST}:{PORT} ({e}). "
            f"Is the PyMOL plugin running?",
        }
    except TimeoutError as e:
        if connected:
            # The kernel completes the connection before the plugin accepts
            # it, so "connected" does not mean PyMOL has started on this
            # request: it may still be busy with an earlier one.
            error = (
                f"PyMOL accepted '{action}' but sent no reply within {timeout:g} s. It may "
                f"still be running it, or still be busy with an earlier request; it was not "
                f"cancelled, so check what PyMOL has done before re-sending."
            )
        else:
            error = (
                f"Socket connection failed: {e} connecting to {HOST}:{PORT}. The plugin may be "
                f"busy with a long request (it serves one at a time) or PyMOL may be hung."
            )
        return {
            "status": "error",
            "error": error + _asleep_sentence(_asleep_since(awake_start, total_start)),
        }
    except Exception as e:
        return {
            "status": "error",
            "error": f"Socket connection failed: {e}. Is the PyMOL plugin running?",
        }

    asleep = _asleep_since(awake_start, total_start)
    if not chunks:
        return {
            "status": "error",
            "error": "Empty response from PyMOL plugin." + _asleep_sentence(asleep),
        }
    if res is None:
        try:
            res = json.loads(b"".join(chunks).decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            return {
                "status": "error",
                "error": f"Malformed response from PyMOL plugin: {e}" + _asleep_sentence(asleep),
            }
    if not isinstance(res, dict):
        return {
            "status": "error",
            "error": f"Malformed response from PyMOL plugin: expected a JSON object, got "
            f"{type(res).__name__}." + _asleep_sentence(asleep),
        }
    if asleep > _ASLEEP_REPORT_S:
        res["asleep_s"] = round(asleep, 1)
    return res


def _asleep_sentence(asleep: float) -> str:
    """' The computer was asleep for N s of this wait.' — or nothing, for a
    sleep too short to explain anything."""
    if asleep <= _ASLEEP_REPORT_S:
        return ""
    return f" The computer was asleep for {asleep:.0f} s of this wait."


def _warnings_note(res: dict) -> str:
    """'\\nPyMOL warned: ...' for a response that carries PyMOL's warnings.

    The plugin collects the warning and error lines PyMOL printed while a
    command it still counted as a success ran — a setting stored at a level it
    is never read from, a selection that matched nothing.  Without this they
    reached only PyMOL's console, and the tool reported plain success.
    """
    warnings = [str(w) for w in (res.get("warnings") or [])]
    if not warnings:
        return ""
    if len(warnings) == 1:
        return f"\nPyMOL warned: {warnings[0]}"
    return "\nPyMOL warned:\n" + "\n".join(f"  {w}" for w in warnings)


def format_measurement(
    label: str, value: Any, unit: str, name: str | None = None, signed: bool = False
) -> str:
    """Render a numeric PyMOL result as an answer rather than an acknowledgement.

    ``cmd.distance`` / ``angle`` / ``dihedral`` / ``get_area`` all *return* the
    quantity they measure, and the wrappers used to discard it — so asking
    MCPymol to measure something got you a drawing and no number.

    PyMOL signals "nothing matched" by returning -1, which is worth saying
    plainly rather than reporting as a measurement of -1.  That inference is
    only valid for quantities that cannot be negative, though: set ``signed``
    for dihedrals, where -57.8 deg is an ordinary alpha-helical phi and not a
    failure at all.
    """
    stored = f" (stored as '{name}')" if name else ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return f"{label} was performed{stored}, but PyMOL returned no numeric value ({value!r})."

    if not signed and number < 0:
        return (
            f"{label} failed{stored}: PyMOL returned {number:g}, which means the "
            f"selections matched no atoms, or no pair fell within the cutoff."
        )
    return f"{label}: {number:.2f} {unit}{stored}."


def _call(
    _action: str,
    /,
    *,
    _timeout: float = _DEFAULT_TIMEOUT,
    _measures: tuple[str, str] | tuple[str, str, bool] | None = None,
    **values: str | None,
) -> str:
    """Forward a primitive PyMOL command and report the outcome.

    This is the shared body of every thin ``pymol.cmd`` wrapper below.  Each
    wrapper is still a real ``def`` — the signature and docstring are what the
    MCP client sees, and keeping them written out means mypy, ruff and grep all
    still work — but the mechanical part (drop unset arguments, forward,
    translate the response) lives here so a fix lands once instead of 77 times.

    ``values`` is an ordered name → value mapping matching the wrapper's
    parameters.  Unset (``None``) values are dropped from the end, because
    PyMOL takes these positionally: ``ray(height="600")`` with no width cannot
    send just the height without PyMOL reading it *as* the width.  The wrappers
    used to do exactly that, silently.  A gap now returns an error naming the
    parameters instead of quietly producing a wrong render.

    ``_measures`` is ``(label, unit)`` — or ``(label, unit, signed)`` — for
    commands that return a quantity; with it, the result is reported as a
    number instead of "Executed ... successfully".

    Warnings PyMOL printed while the command ran are appended either way.
    """
    names = list(values)
    args = [values[n] for n in names]

    while args and args[-1] is None:
        args.pop()
        names.pop()

    gaps = [n for n, v in zip(names, args, strict=True) if v is None]
    if gaps:
        given = names[-1]
        return (
            f"Error: {_action} was given '{given}' without "
            f"{' and '.join(repr(g) for g in gaps)}. PyMOL takes these arguments "
            f"positionally, so earlier ones cannot be skipped — pass them all, "
            f"or drop the later one."
        )

    res = send_request(_action, args=args, timeout=_timeout)
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    if _measures is not None:
        label, unit = _measures[0], _measures[1]
        signed = bool(_measures[2]) if len(_measures) > 2 else False
        return format_measurement(
            label, res.get("result"), unit, values.get("name"), signed
        ) + _warnings_note(res)
    return f"Executed {_action} successfully." + _warnings_note(res)
