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
from collections.abc import Callable
from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    _BRIDGE_VERSION = version("mcpymol")
except PackageNotFoundError:
    _BRIDGE_VERSION = ""

_WARNED_VERSIONS: set[str] = set()

HOST = "127.0.0.1"

# Port can be overridden via environment variable, e.g.: MCPYMOL_PORT=9867 uv run mcpymol
PORT = int(os.environ.get("MCPYMOL_PORT", 9876))

# Timeout for long-running operations (ray-tracing, movie frames, large saves).
_SLOW_OP_TIMEOUT = float(os.environ.get("MCPYMOL_SLOW_OP_TIMEOUT", 600.0))

# Single recv() chunk size; send_request loops until the plugin half-closes.
_RECV_CHUNK = 65536

# Default socket budget for interactive commands.
_DEFAULT_TIMEOUT = 10.0

# Minimum sleep duration (seconds) worth reporting on a slow or timed-out call.
_ASLEEP_REPORT_S = 5.0


def _clock_with_sleep() -> float | None:
    """Return a monotonic clock reading that continues counting during system sleep, or None."""
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
                # Only attempt an incremental JSON parse on a short read so
                # multi-chunk responses aren't re-parsed on every chunk.
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
    _check_plugin_version(res)
    return res


def _check_plugin_version(res: dict) -> None:
    """Warn once per session on stderr if the PyMOL plugin version differs from the bridge."""
    if _WARNED_VERSIONS or not _BRIDGE_VERSION:
        return
    plugin_ver = res.get("plugin_version")
    if plugin_ver == _BRIDGE_VERSION:
        return
    reported = str(plugin_ver) if plugin_ver else "<=1.6.0 (no plugin_version)"
    _WARNED_VERSIONS.add(reported)
    print(
        f"MCPymol warning: bridge is v{_BRIDGE_VERSION} but the PyMOL plugin is {reported}. "
        f"Run `mcpymol --install-plugin` and restart PyMOL so both match.",
        file=sys.stderr,
    )


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


def format_superposition(result: Any) -> str | None:
    """align's and super's 7-number answer, in words; None if it is not one.

    PyMOL returns (RMSD after refinement, atoms in it, refinement cycles, RMSD
    before refinement, atoms before, raw score, residues aligned).
    """
    if not isinstance(result, (list, tuple)) or len(result) < 7:
        return None
    try:
        rms, atoms, cycles, rms0, atoms0 = (
            float(result[0]),
            int(result[1]),
            int(result[2]),
            float(result[3]),
            int(result[4]),
        )
        residues = int(result[6])
    except (TypeError, ValueError):
        return None
    return (
        f"RMSD {rms:.2f} A over {atoms} atoms after {cycles} refinement "
        f"cycle{'s' if cycles != 1 else ''} (before refinement: {rms0:.2f} A over "
        f"{atoms0} atoms); {residues} residues aligned."
    )


def format_per_state_rmsd(result: Any, moved: bool = False) -> str | None:
    """intra_rms / intra_fit's per-state RMSDs, in words; None if not that shape.

    PyMOL answers one value per state, with -1.0 at the reference state.
    """
    if not isinstance(result, (list, tuple)) or not result:
        return None
    try:
        values = [float(v) for v in result]
    except (TypeError, ValueError):
        return None
    reference = next((i for i, v in enumerate(values) if v == -1.0), None)
    others = [(i + 1, v) for i, v in enumerate(values) if i != reference]
    if reference is None or not others:
        return None
    listed = ", ".join(f"state {n}: {v:.2f} A" for n, v in others)
    text = f"RMSD of each state against state {reference + 1}: {listed}."
    if moved:
        text += (
            f" States were fitted onto state {reference + 1}, so their coordinates in "
            f"the session have changed."
        )
    return text


def _call(
    _action: str,
    /,
    *,
    _timeout: float = _DEFAULT_TIMEOUT,
    _measures: tuple[str, str] | tuple[str, str, bool] | None = None,
    _format: Callable[[Any], str | None] | None = None,
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
    number instead of "Executed ... successfully". ``_format`` does the same
    for a structured answer (align's seven numbers, intra_rms's list); it
    returns None for anything it does not recognise, which falls back to the
    plain wording rather than misreporting.

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
    if _format is not None:
        text = _format(res.get("result"))
        if text is not None:
            return text + _warnings_note(res)
    if _measures is not None:
        label, unit = _measures[0], _measures[1]
        signed = bool(_measures[2]) if len(_measures) > 2 else False
        return format_measurement(
            label, res.get("result"), unit, values.get("name"), signed
        ) + _warnings_note(res)
    return f"Executed {_action} successfully." + _warnings_note(res)
