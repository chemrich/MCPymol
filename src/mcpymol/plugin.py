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

import io
import json
import os
import re
import socket
import sys
import tempfile
import threading
import time
import traceback
import uuid
from typing import IO, TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable

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

# Cap on the PyMOL output one `do` sends back.  "iterate all, print(name)"
# prints a line per atom; the head and the tail are kept, the middle counted.
MAX_OUTPUT_CHARS = int(os.environ.get("MCPYMOL_MAX_OUTPUT_CHARS", 8000))

# How long a `do` waits for PyMOL's main thread to *start* it.  Once started it
# is waited on to the end, however long that takes; this only bounds how long a
# main thread that is stuck (a modal dialog, a hung command) can hold a request
# that has not run at all.
DO_START_TIMEOUT = float(os.environ.get("MCPYMOL_DO_START_TIMEOUT", 120.0))

# macOS only: seconds without a request after which PyMOL may App Nap again.
# 0 turns the guard off.  See _AppNapGuard.
APP_NAP_IDLE_S = float(os.environ.get("MCPYMOL_APP_NAP_IDLE_S", 900.0))


# How often the accept loop checks whether it has been asked to stop.
_ACCEPT_POLL_INTERVAL = 0.25


class RequestTooLarge(Exception):
    """Raised when a peer sends more than ``MAX_REQUEST_BYTES`` in one request."""


class MainThreadTimeout(Exception):
    """PyMOL's main thread did not start a queued job in time; nothing was run."""


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


def _jsonable(value: object) -> object:
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


def _resolve_dotted(name: str) -> Callable[..., Any] | None:
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


def _error_text(exc: BaseException) -> str:
    """How an exception reads in an error response.

    PyMOL's own exceptions carry a message written for the user, so they stand
    alone.  Anything else is a Python error, and its type is half of what it
    says: "name 'x' is not defined" is a NameError, "'x'" alone is a KeyError.
    """
    if type(exc).__name__ in ("CmdException", "QuietException") or type(exc) is Exception:
        return str(exc)
    return f"{type(exc).__name__}: {exc}"


def _print_traceback(exc: BaseException) -> None:
    """Print the traceback of the exception being handled, whatever happens.

    Under Python 3.14, formatting a NameError raised inside an alter or
    iterate expression raises TypeError ('wrapper.Wrapper' object is not
    iterable): the traceback module looks through the frame's locals for a
    name to suggest, and PyMOL evaluates those expressions with a Wrapper as
    locals.  That TypeError used to escape the request handler, so the client
    got an empty response instead of the NameError.
    """
    try:
        traceback.print_exc()
    except Exception:
        print(f"{type(exc).__name__}: {exc}")


# ── Seeing what PyMOL said ──────────────────────────────────────────────────
#
# PyMOL reports most problems by printing them, not by raising:
#
# * Parser-level failures (a bad selection, an unknown setting name, a Python
#   error in an alter or iterate expression) are printed from Python, and the
#   parser returns 0 instead of 1.
# * Some problems exist only as C-level feedback, written straight to file
#   descriptor 1 (and, in the GUI, to the console's queue).  The command still
#   "succeeds", so this text is their only trace:
#     " Setting-Warning: 'cartoon_tube_radius' is a object-state-level setting"
#     " Selector-Error: Invalid selection name "x"."   (zoom, orient)
#     "Setting-Error: unknown color 'x'"                 (bg_color)
#     " ScenePNG-Error: error writing "/no/dir/x.png"!"  (png, ray=0)

# "<Module>-Error:", "<Module>-Warning:", " Error:", " Warning:".  Not Python's
# "NameError:" — a failed line is reported through the parser's verdict.
_PROBLEM_LINE = re.compile(r"^\s*(?:[A-Za-z][\w.]*-)?(?:Error|Warning)s?\b")

# A selector error is followed by the offending expression and a caret.
_CARET_LINE = re.compile(r"<--\s*$")

# How many warnings, or rejected lines, one reply spells out.  A command can
# warn once per atom, and a 20,000-line script can fail on every line; the
# rest are counted, and their text is still in `output` (within its own cap).
_MAX_LISTED = 20

# How much of a rejected line an error message quotes.
_MAX_QUOTED_LINE = 300


def _problem_lines(text: str) -> list[str]:
    """The warning and error lines in ``text``, each with its caret line folded in."""
    found = []
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if _PROBLEM_LINE.match(line):
            message = line.strip()
            if i + 1 < len(lines) and _CARET_LINE.search(lines[i + 1]):
                message += " " + lines[i + 1].strip()
            found.append(message)
    return found


def _listed(problems: list[str]) -> list[str]:
    """``problems`` without repeats, and at most ``_MAX_LISTED`` of them plus a
    count of the rest."""
    unique = list(dict.fromkeys(problems))
    if len(unique) <= _MAX_LISTED:
        return unique
    return [*unique[:_MAX_LISTED], f"... and {len(unique) - _MAX_LISTED} more"]


def _quoted(line: str) -> str:
    """A command line as an error message quotes it: whole unless it is huge."""
    if len(line) <= _MAX_QUOTED_LINE:
        return line
    return f"{line[:_MAX_QUOTED_LINE]}... [{len(line) - _MAX_QUOTED_LINE} more characters]"


def _truncate(text: str, limit: int | None = None) -> str:
    """``text`` cut to about ``limit`` characters, keeping its head and tail."""
    limit = MAX_OUTPUT_CHARS if limit is None else limit
    if len(text) <= limit:
        return text
    keep = max(limit, 0) // 2
    middle = text[keep : len(text) - keep]
    lines = middle.count("\n")
    omitted = f"{lines} lines" if lines else f"{len(middle)} characters"
    return f"{text[:keep]}\n... [{omitted} omitted] ...\n{text[len(text) - keep :]}"


def _write_fd(fd: int, data: bytes) -> None:
    """os.write all of ``data``; a short write is legal on a pipe."""
    view = memoryview(data)
    while view:
        try:
            written = os.write(fd, view)
        except OSError:
            return
        view = view[written:]


class _Tap(io.StringIO):
    """Stands in for sys.stdout or sys.stderr while output is captured.

    Records each write tagged with the stream it was meant for, so the capture
    reads as one text in the order it was printed and can still be replayed to
    the stream each part belongs to.
    """

    def __init__(self, record: list[tuple[str, str]], stream_name: str, original: IO[str]) -> None:
        super().__init__()
        self._record = record
        self._stream_name = stream_name
        self._original = original

    def write(self, s: str) -> int:
        self._record.append((self._stream_name, s))
        return len(s)

    # A StringIO has no file descriptor and no encoding, so code that asks the
    # stream for them (subprocess.run(stdout=sys.stdout), faulthandler, some C
    # extensions) failed inside a capture. Answer for the real stream: fd 1 is
    # pointed at the spool meanwhile, so what is written there is still caught.
    def fileno(self) -> int:
        return self._original.fileno()

    @property
    def encoding(self) -> str:  # type: ignore[override]
        return getattr(self._original, "encoding", None) or "utf-8"

    def isatty(self) -> bool:
        return False


class _OutputCapture:
    """Capture what PyMOL prints while a command runs, without hiding it.

    Python-level text (parser errors, tracebacks, print) is caught by swapping
    ``sys.stdout``/``sys.stderr``.  C-level feedback never passes through
    Python, so file descriptor 1 is pointed at a spool file for the duration;
    PyMOL flushes after every feedback line (captures were byte-identical with
    and without an explicit fflush), so nothing is left in a C buffer.

    Both are replayed on exit, so the terminal still shows everything.  The GUI
    console reads PyMOL's own feedback queue, which this never touches; draining
    that queue instead would steal the lines from the console.

    One spool file is reused (truncated per capture): creating a temporary file
    per capture cost ~1.5 ms a request against ~20 us.  A capture nested inside
    another gets its own.  fd capture is best-effort — where fd 1 cannot be
    duplicated (closed, as under pythonw), only the Python-level text is caught,
    and fd 1 is left exactly as it was.

    The spool and its lock are reached through ``type(self)``, never the global
    name: a `do` of "run plugin.py" re-executes this file while a capture made
    from the old class is open around that very command, and the global then
    names a new class whose lock this capture never took.
    """

    _spool: IO[bytes] | None = None
    _spool_lock = threading.Lock()

    def __enter__(self) -> _OutputCapture:
        self.c_text = ""
        self.python_text = ""
        self._record: list[tuple[str, str]] = []
        self._streams = {"stdout": sys.stdout, "stderr": sys.stderr}
        # Anything already buffered (an echoed command line) belongs before
        # the output it introduces, not after the replay.
        for stream in self._streams.values():
            try:
                stream.flush()
            except Exception:
                pass
        self._saved_fd: int | None = None
        self._spool_file: IO[bytes] | None = None
        self._owns_spool = False
        self._start_fd_capture()
        sys.stdout = _Tap(self._record, "stdout", self._streams["stdout"])
        sys.stderr = _Tap(self._record, "stderr", self._streams["stderr"])
        return self

    def _start_fd_capture(self) -> None:
        cls = type(self)
        try:
            # Duplicate first: if fd 1 is closed, dup2 below would *create* it
            # pointing at the spool, and leave it there.
            saved = os.dup(1)
        except OSError:
            return
        spool = None
        try:
            # Both outlive any with-block on purpose: the shared spool lives as
            # long as PyMOL, a private one until __exit__ closes it.
            if cls._spool_lock.acquire(blocking=False):
                self._owns_spool = True
                if cls._spool is None:
                    cls._spool = tempfile.TemporaryFile()  # noqa: SIM115
                spool = cls._spool
            else:
                spool = tempfile.TemporaryFile()  # noqa: SIM115
            fd = spool.fileno()
            os.ftruncate(fd, 0)
            os.lseek(fd, 0, os.SEEK_SET)
            os.dup2(fd, 1)
        except OSError:
            os.close(saved)
            self._spool_file = spool
            self._release_spool()
            return
        self._saved_fd = saved
        self._spool_file = spool

    def _release_spool(self) -> None:
        if self._owns_spool:
            self._owns_spool = False
            type(self)._spool_lock.release()
        elif self._spool_file is not None:
            try:
                self._spool_file.close()
            except OSError:
                pass
        self._spool_file = None

    def _read_spool(self) -> bytes:
        if self._spool_file is None:
            return b""
        fd = self._spool_file.fileno()
        chunks = []
        try:
            os.lseek(fd, 0, os.SEEK_SET)
            while True:
                chunk = os.read(fd, 65536)
                if not chunk:
                    break
                chunks.append(chunk)
        except OSError:
            pass
        return b"".join(chunks)

    def __exit__(self, *exc_info: object) -> None:
        sys.stdout, sys.stderr = self._streams["stdout"], self._streams["stderr"]
        raw = b""
        if self._saved_fd is not None:
            try:
                os.dup2(self._saved_fd, 1)
            finally:
                os.close(self._saved_fd)
                self._saved_fd = None
            raw = self._read_spool()
            self._release_spool()
            if raw:
                _write_fd(1, raw)
        self.c_text = raw.decode("utf-8", "replace")
        self.python_text = "".join(text for _, text in self._record)
        self._replay()

    def _replay(self) -> None:
        """Write the captured Python text to the streams it was meant for."""
        merged: list[list[str]] = []
        for name, text in self._record:
            if merged and merged[-1][0] == name:
                merged[-1][1] += text
            else:
                merged.append([name, text])
        for name, text in merged:
            stream = self._streams[name]
            try:
                stream.write(text)
                # Per chunk: a block-buffered stdout would otherwise land after
                # the stderr text that followed it.
                stream.flush()
            except Exception:
                pass

    @property
    def text(self) -> str:
        return self.c_text + self.python_text


# ── Reading a `do` line the way PyMOL's C CmdDo does ────────────────────────


def _do_prefix(raw: str) -> tuple[str, str | None, bool]:
    """Mirror C ``CmdDo``: the line to parse, the text to echo, and whether to log it.

    ``_ cmd`` runs without an echo but is still logged; ``_cmd.x`` and
    ``cmd._x`` are internal and neither echoed nor logged; a pasted ``PyMOL>``
    prompt is stripped before parsing (the raw line is what gets echoed).
    """
    if raw.startswith("_"):
        if raw[1:2] == " ":
            return raw[2:], None, True
        return raw, None, False
    if raw.startswith("cmd._"):
        return raw, None, False
    line = raw
    if line.startswith("PyMOL>"):
        line = line[7:] if line[6:7] == " " else line[6:]
    return line, raw, True


def _literal_command(line: str) -> str | None:
    """The command name if ``line`` starts with one that takes the rest of the
    line literally (alter, iterate, label, ...), so a ``;`` does not split it."""
    word = line.strip().split(" ", 1)[0].split(",", 1)[0]
    try:
        from pymol import parsing

        name = word
        if name not in cmd.keyword:
            expanded = cmd.kwhash.interpret(name)
            if isinstance(expanded, str):
                name = expanded
        entry = cmd.keyword.get(name)
        if entry and len(entry) > 4 and entry[4] >= parsing.LITERAL:
            return name
    except Exception:
        pass
    return None


def _original_error(exc: BaseException) -> str:
    """The error worth reporting when the parser itself raised.

    parse() catches everything a command raises, so what escapes it is its own
    error printer failing.  Under Python 3.14, PyMOL 3.1's printer raises
    TypeError ('wrapper.Wrapper' object is not iterable) while formatting a
    NameError from an alter, iterate or label expression; the error it was
    printing is at the root of the context chain.
    """
    seen = set()
    while exc.__context__ is not None and id(exc) not in seen:
        seen.add(id(exc))
        exc = exc.__context__
    return f"{type(exc).__name__}: {exc}"


def _failure_message(text: str) -> str:
    """PyMOL's account of a failed line: what it printed, minus echoes, and for
    a traceback only the last frame and the exception."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("PyMOL>")]
    if any(ln.startswith("Traceback (most recent call last)") for ln in lines):
        frames = [i for i, ln in enumerate(lines) if ln.startswith("  File ")]
        if frames:
            lines = lines[frames[-1] :]
            lines[0] = lines[0].strip()
    return _truncate("\n".join(lines))


# ── Running `do` on PyMOL's main thread, and waiting for it ─────────────────
#
# cmd.do from this (socket) thread does not run anything: it queues the line
# for PyMOL's main thread and returns after at most ~0.13 s, finished or not.
# So the old reply could arrive before the command ran, and could never carry
# its error.  Parsing on this thread instead is not an option either: in the
# Qt GUI, lines that touch widgets ("window size, 640, 480") killed PyMOL when
# run off the main thread, and holding the API lock while waiting for the main
# thread deadlocks png(ray=0), whose GUI-thread half needs that lock.
#
# So a request still reaches the main thread through PyMOL's own command
# queue, the path cmd.do always used — but as one Python call this thread can
# wait on: "cmd._mcpymol_run_job('<token>')".  PyMOL neither echoes nor logs a
# "cmd._" line, and runs it on the main thread under the API lock, where the
# job parses every line of the request itself, keeping the parser's verdict and
# what each line printed, then sets an Event.  PyMOL's cmd._call_in_gui_thread
# was the obvious alternative and was measured hanging after ~11 requests: its
# QWaitCondition loses the wake-up when the call finishes first.


class _MainThreadJob:
    """One piece of work for PyMOL's main thread, as the plugin thread sees it."""

    def __init__(self, func: Callable[[], Any]) -> None:
        self._func = func
        self._lock = threading.Lock()
        self.started = threading.Event()
        self.done = threading.Event()
        self.cancelled = False
        self.value = None
        self.error: BaseException | None = None

    def run(self) -> None:
        """Called on the main thread.  A job the plugin gave up on does nothing:
        its caller has already said that nothing was run."""
        with self._lock:
            if self.cancelled:
                return
            self.started.set()
        try:
            self.value = self._func()
        except BaseException as e:
            self.error = e
        finally:
            self.done.set()

    def cancel_unless_started(self) -> bool:
        """Cancel the job if the main thread has not picked it up; True if cancelled."""
        with self._lock:
            if self.started.is_set():
                return False
            self.cancelled = True
            return True


def _job_registry() -> dict:
    """Queued jobs by token.  Kept on the cmd module rather than here, so a
    re-run of plugin.py (a fresh module namespace) still finds them."""
    jobs = getattr(cmd, "_mcpymol_jobs", None)
    if not isinstance(jobs, dict):
        jobs = {}
        cmd._mcpymol_jobs = jobs
    return jobs


def _run_queued_job(token: str) -> None:
    """The main-thread end of the line _on_main_thread queues."""
    job = _job_registry().pop(token, None)
    if job is not None:
        job.run()


def _on_pymol_main_thread() -> bool:
    """True on the thread that runs PyMOL's command queue."""
    if threading.current_thread() is threading.main_thread():
        return True
    gui_thread = getattr(getattr(cmd, "_pymol", None), "glutThread", None)
    return isinstance(gui_thread, int) and gui_thread == threading.get_ident()


def _modal_draw(api: Any) -> bool:
    """True while PyMOL renders frames in a modal draw loop (``mpng`` with
    ``modal=1``), during which it refuses every queued command.  Called with
    the API lock held."""
    try:
        value = api.get_modal_draw(cmd._COb)
    except Exception:
        return False
    return isinstance(value, int) and bool(value)


def _queue_job(token: str, deadline: float, timeout: float) -> None:
    """Put the job's line on PyMOL's command queue, within the start budget.

    This is what cmd.do does for one line, with its two open-ended waits
    bounded by the time left for the job to start:

    * PyMOL's API lock.  cmd.do takes it with no timeout, and PyMOL's main
      thread can hold it for as long as a command runs — measured headless: a
      4 s line queued by another thread held it the whole 4 s.  A main thread
      stuck in a command would have held this request, and so the plugin,
      forever, and a deadline counted around that wait expired unseen.
    * A modal draw (mpng with modal=1).  C CmdDo refuses every line until it
      ends — cmd.do raised 'APIEnterNotModal(G)' at once (measured) — and the
      render ends by itself, so wait for it.
    """
    api = cmd._cmd
    lock = cmd.lock_api
    line = f"cmd._mcpymol_run_job('{token}')"
    while True:
        if not lock.acquire(timeout=max(0.0, deadline - time.monotonic())):
            raise MainThreadTimeout(
                f"PyMOL was busy for {timeout:g} s — its main thread held PyMOL's API lock "
                f"the whole time: a long command, an open dialog, or a hang — so nothing was "
                f"run. Check the PyMOL window before re-sending. (MCPYMOL_DO_START_TIMEOUT "
                f"sets this wait.)"
            )
        modal = False
        try:
            api.do(cmd._COb, line, 0, 0)
        except Exception:
            modal = _modal_draw(api)
            if not modal:
                raise
        finally:
            lock.release()
        if not modal:
            return
        if time.monotonic() >= deadline:
            # Seen in the GUI: with PyMOL's window minimized, a modal mpng drew
            # no frames at all, and stayed modal until the window came back.
            raise MainThreadTimeout(
                f"PyMOL was still rendering a movie (a modal draw loop, as mpng with modal=1 "
                f"uses) after {timeout:g} s, and takes no commands until it finishes, so "
                f"nothing was run. Re-send once the movie is done; in the GUI it renders only "
                f"while the PyMOL window is visible, not minimized. (MCPYMOL_DO_START_TIMEOUT "
                f"sets this wait.)"
            )
        time.sleep(0.1)


def _on_main_thread(func: Callable[[], Any], start_timeout: float | None = None) -> Any:
    """Run ``func`` on PyMOL's main thread and return what it returns.

    Waits up to ``start_timeout`` (default ``DO_START_TIMEOUT``) for the main
    thread to start it, then however long it takes to finish.  Never holds
    the API lock while waiting.  Already on the main thread, it just calls it.

    Raises:
        MainThreadTimeout: the main thread never started it; it never will.
    """
    if _on_pymol_main_thread():
        return func()
    timeout = DO_START_TIMEOUT if start_timeout is None else start_timeout
    deadline = time.monotonic() + timeout
    job = _MainThreadJob(func)
    token = uuid.uuid4().hex
    jobs = _job_registry()
    jobs[token] = job
    # Set on every call: a re-run of plugin.py must route jobs to its own code.
    cmd._mcpymol_run_job = _run_queued_job
    try:
        _queue_job(token, deadline, timeout)
    except BaseException:
        jobs.pop(token, None)
        raise
    remaining = max(0.0, deadline - time.monotonic())
    if not job.started.wait(remaining) and job.cancel_unless_started():
        jobs.pop(token, None)
        raise MainThreadTimeout(
            f"PyMOL's main thread did not pick up the command within {timeout:g} s, so "
            f"nothing was run. It may be busy with a long command, held by an open "
            f"dialog, or hung; check the PyMOL window before re-sending. "
            f"(MCPYMOL_DO_START_TIMEOUT sets this wait.)"
        )
    job.done.wait()
    if job.error is not None:
        if isinstance(job.error, Exception):
            raise job.error
        # SystemExit and friends belong to the main thread, not to this one.
        raise RuntimeError(_error_text(job.error)) from job.error
    return job.value


def _new_parser() -> Any:
    """A private command parser for one request.

    Never the console's cmd._parser: its state (an open python block, a
    trailing-backslash continuation) is shared with whatever the user types, so
    a request that left one open would swallow the user's next command.
    """
    from pymol.parser import Parser

    pymol_module = cmd._pymol
    script = getattr(pymol_module, "__script__", None)
    parser = Parser(cmd)
    if script is not None:
        # Parser() resets this to top level; a script the main thread is in the
        # middle of (this job can run nested inside one) still owns it.
        pymol_module.__script__ = script
    return parser


def _run_lines(lines: list[str]) -> dict:
    """Parse ``lines`` one at a time and record PyMOL's verdict on each.

    Runs on PyMOL's main thread.  Does with each line what C CmdDo does —
    echo, log, parse — except that it keeps the parser's return value (1 ran,
    0 failed) and the text the line printed, where CmdDo throws both away.
    Lines keep running after one fails, as they do with cmd.do.
    """
    parser = _new_parser()
    failures: list[dict] = []
    warnings: list[str] = []
    outputs: list[str] = []
    not_run: set[int] = set()
    block_start = None  # the line that opened a python/skip block still open
    continued_from = None  # the line a trailing-backslash continuation began on
    defer = None
    if len(lines) > 1:
        # As cmd.do does: no redraws between the lines of one request.
        defer = cmd.get_setting_int("defer_updates")
        cmd.set("defer_updates", 1)
    try:
        for number, raw in enumerate(lines, 1):
            line, echo, logged = _do_prefix(raw)
            if echo is not None:
                print("PyMOL>" + echo)
            if logged and line.strip().split(" ", 1)[0].lower() != "quit":
                try:
                    cmd.log(line + "\n")
                except Exception:
                    pass
            raised = None
            with _OutputCapture() as captured:
                try:
                    status = parser.parse(line, 0)
                except Exception as e:
                    status, raised = 0, e
                    print("Error: " + _original_error(e))
            text = captured.text
            if text:
                outputs.append(text)

            top = getattr(parser, "layer", {}).get(1)
            if getattr(top, "embed_sentinel", None) is None:
                block_start = None
            elif block_start is None:
                block_start = number
            if not getattr(top, "cont", ""):
                continued_from = None
            elif continued_from is None:
                continued_from = number

            if status == 0:
                message = _original_error(raised) if raised is not None else _failure_message(text)
                literal = _literal_command(line)
                if literal and ";" in line:
                    message += (
                        f"\nNote: {literal} takes the rest of the line literally, so ';' does "
                        f"not start a new command there. Put each command on its own line."
                    )
                failures.append(
                    {"line": number, "command": raw, "message": message or "PyMOL gave no reason."}
                )
                not_run.add(number)
            else:
                warnings.extend(_problem_lines(text))
    finally:
        if defer is not None:
            cmd.set("defer_updates", defer)

    top = getattr(parser, "layer", {}).get(1)
    if block_start is not None:
        sentinel = getattr(top, "embed_sentinel", None)
        failures.append(
            {
                "line": block_start,
                "command": lines[block_start - 1],
                "message": f"This opens a block that is never closed (PyMOL is still waiting "
                f"for '{sentinel}'), so nothing after it was run.",
            }
        )
        not_run.update(range(block_start, len(lines) + 1))
    elif continued_from is not None:
        failures.append(
            {
                "line": continued_from,
                "command": lines[continued_from - 1],
                "message": "The text ends in a '\\' continuation with nothing after it, so "
                "this command was never run.",
            }
        )
        not_run.update(range(continued_from, len(lines) + 1))
    failures.sort(key=lambda f: f["line"])
    return {
        "failures": failures,
        "warnings": warnings,
        "output": "".join(outputs),
        "ran": len(lines) - len(not_run),
    }


def _do_response(commands: object, lines: list[str], report: dict) -> dict:
    """The reply to a `do`.  ``output`` is always present: the bridge reads its
    absence as a plugin too old to report errors at all."""
    response: dict = {
        "status": "success",
        "result": f"Executed command: {commands}",
        "output": _truncate(report["output"]),
    }
    if report["warnings"]:
        response["warnings"] = _listed(report["warnings"])
    failures = report["failures"]
    if not failures:
        return response
    first = failures[0]
    where = f"line {first['line']} of {len(lines)}" if len(lines) > 1 else "the command"
    error = f"PyMOL rejected {where}: {_quoted(first['command'])}\n{first['message']}"
    if len(failures) > 1:
        others = failures[1 : 1 + _MAX_LISTED]
        error += "\nAlso rejected: " + "; ".join(
            f"line {f['line']} ({_quoted(f['command'])})" for f in others
        )
        if len(failures) > 1 + len(others):
            error += f"; and {len(failures) - 1 - len(others)} more lines"
    if len(lines) > 1:
        error += f"\nThe other {report['ran']} line(s) ran." if report["ran"] else "\nNo line ran."
    response.update(status="error", error=error, failures=failures[:_MAX_LISTED])
    del response["result"]
    return response


def run_do(commands: object) -> dict:
    """Run PyMOL command-language text and report what PyMOL made of it.

    Splits ``commands`` into lines the way cmd.do does, runs them on PyMOL's
    main thread through a private parser, and replies only once they have run:

    * ``status`` is ``error`` when PyMOL rejected any line, and ``error`` says
      which, in PyMOL's own words; ``failures`` lists each one;
    * ``output`` is what the lines printed (capped at ``MAX_OUTPUT_CHARS``);
    * ``warnings`` lists PyMOL's warning and error lines from lines that ran
      anyway — a level-mismatched ``set``, ``zoom`` on a bad selection.
    """
    command_list = commands if isinstance(commands, (list, tuple)) else [commands]
    lines = [line for c in command_list for line in str(c).splitlines() if line]
    if not lines:
        # Nothing for the main thread to do; cmd.do("") is a no-op too.
        return _do_response(commands, lines, {"failures": [], "warnings": [], "output": ""})
    try:
        report = _on_main_thread(lambda: _run_lines(lines))
    except MainThreadTimeout as e:
        return {"status": "error", "error": str(e), "output": ""}
    return _do_response(commands, lines, report)


# ── Keeping macOS App Nap off a PyMOL that is being driven ──────────────────

# NSActivityUserInitiatedAllowingIdleSystemSleep.  NSActivityUserInitiated is
# 0x00FFFFFF: one more bit, NSActivityIdleSystemSleepDisabled (1 << 20), which
# would keep the Mac awake as long as PyMOL is open.
_NS_ACTIVITY_NO_APP_NAP = 0x00EFFFFF


class _ObjCRuntime:
    """Just enough of the Objective-C runtime, through ctypes, to begin and end
    an NSProcessInfo activity."""

    def __init__(self) -> None:
        import ctypes
        import ctypes.util

        path = ctypes.util.find_library("objc")
        if not path:
            raise OSError("libobjc not found")
        objc = ctypes.CDLL(path)
        # NSProcessInfo lives in Foundation.  A Qt GUI has it loaded already;
        # loading it here makes that an observation rather than an assumption.
        ctypes.CDLL("/System/Library/Frameworks/Foundation.framework/Foundation")
        objc.objc_getClass.restype = ctypes.c_void_p
        objc.objc_getClass.argtypes = [ctypes.c_char_p]
        objc.sel_registerName.restype = ctypes.c_void_p
        objc.sel_registerName.argtypes = [ctypes.c_char_p]
        objc.objc_autoreleasePoolPush.restype = ctypes.c_void_p
        objc.objc_autoreleasePoolPush.argtypes = []
        objc.objc_autoreleasePoolPop.restype = None
        objc.objc_autoreleasePoolPop.argtypes = [ctypes.c_void_p]
        send = ctypes.cast(objc.objc_msgSend, ctypes.c_void_p).value
        if not send:
            raise OSError("objc_msgSend not found")

        def method(restype: Any, *argtypes: Any) -> Any:
            # objc_msgSend has to be called through a prototype that matches
            # the method exactly.  On arm64, variadic and fixed arguments are
            # passed differently, so the generic signature scrambles them — a
            # crash try/except cannot catch.
            return ctypes.CFUNCTYPE(restype, ctypes.c_void_p, ctypes.c_void_p, *argtypes)(send)

        self._objc = objc
        self._id = method(ctypes.c_void_p)
        self._id_cstr = method(ctypes.c_void_p, ctypes.c_char_p)
        self._id_u64_id = method(ctypes.c_void_p, ctypes.c_uint64, ctypes.c_void_p)
        self._void_id = method(None, ctypes.c_void_p)
        self._void = method(None)

    def _process_info(self) -> int | None:
        cls = self._objc.objc_getClass(b"NSProcessInfo")
        return self._id(cls, self._objc.sel_registerName(b"processInfo"))

    def begin_activity(self, options: int, reason: str) -> int:
        pool = self._objc.objc_autoreleasePoolPush()
        try:
            sel = self._objc.sel_registerName
            text = self._id_cstr(
                self._objc.objc_getClass(b"NSString"),
                sel(b"stringWithUTF8String:"),
                reason.encode("utf-8"),
            )
            token = self._id_u64_id(
                self._process_info(), sel(b"beginActivityWithOptions:reason:"), options, text
            )
            if not token:
                raise OSError("beginActivityWithOptions:reason: returned nil")
            self._id(token, sel(b"retain"))  # it is autoreleased; outlive the pool
            return token
        finally:
            self._objc.objc_autoreleasePoolPop(pool)

    def end_activity(self, token: int) -> None:
        pool = self._objc.objc_autoreleasePoolPush()
        try:
            sel = self._objc.sel_registerName
            self._void_id(self._process_info(), sel(b"endActivity:"), token)
            self._void(token, sel(b"release"))
        finally:
            self._objc.objc_autoreleasePoolPop(pool)


def _headless() -> bool:
    """True for ``pymol -c``: no window, so nothing for App Nap to throttle."""
    try:
        return bool(cmd._pymol.invocation.options.no_gui)
    except Exception:
        return True


def _qt_window_exists() -> bool:
    """True once PyMOL's Qt main window exists.  Looked up, never imported:
    importing pmg_qt from a PyMOL without it would load Qt for nothing."""
    module = sys.modules.get("pmg_qt.pymol_qt_gui")
    return getattr(module, "window", None) is not None


class _AppNapGuard:
    """Hold off macOS App Nap while an MCP client is using this PyMOL.

    A GUI PyMOL whose window is in the background gets App-Napped: every
    thread, the plugin's socket thread and the ray tracer included, drops to
    Mach priority 4, and the Qt idle timer that runs queued commands is
    throttled.  Measured under load: a 400x300 ray trace took 74-349 s napped
    against 0.85-1.5 s with an NSProcessInfo activity held.

    The activity is taken on the first request and released once no request
    has arrived for ``idle_s`` seconds, so a PyMOL nobody is driving naps as
    usual.  Taking it per request is not enough: a napped thread took 2-76 s
    just to reach the call that would have un-napped it.

    macOS Qt GUI sessions only.  Anything going wrong turns the guard off with
    one printed line; it never affects request handling.
    """

    def __init__(
        self,
        idle_s: float | None = None,
        runtime: Any = None,
        clock: Callable[[], float] | None = None,
    ) -> None:
        self.idle_s = APP_NAP_IDLE_S if idle_s is None else idle_s
        self._runtime = runtime
        self._clock = clock or time.monotonic
        self._lock = threading.Lock()
        self._token = None
        self._last = 0.0
        self.enabled: bool | None = None  # decided on the first request

    @property
    def held(self) -> bool:
        return self._token is not None

    def _disable(self, what: str, exc: Exception) -> None:
        self.enabled = False
        self._token = None
        print(f"MCPymol: App Nap guard off — {what} failed ({type(exc).__name__}: {exc}).")

    def touch(self) -> None:
        """A request arrived: hold the activity, and restart the idle clock."""
        try:
            with self._lock:
                self._last = self._clock()
                if self._token is not None or self.enabled is False:
                    return
                if self.enabled is None:
                    if sys.platform != "darwin" or self.idle_s <= 0 or _headless():
                        self.enabled = False
                        return
                    if not _qt_window_exists():
                        return  # not up yet, or not the Qt GUI: look again next time
                    self.enabled = True
                try:
                    if self._runtime is None:
                        self._runtime = _ObjCRuntime()
                    self._token = self._runtime.begin_activity(
                        _NS_ACTIVITY_NO_APP_NAP, "MCPymol is serving an MCP client"
                    )
                except Exception as e:
                    self._disable("beginning the activity", e)
        except Exception as e:  # never let the guard break a request
            self.enabled = False
            print(f"MCPymol: App Nap guard off ({type(e).__name__}: {e}).")

    def tick(self) -> None:
        """Called on the accept loop's idle poll: release once idle long enough."""
        if self._token is not None and self._clock() - self._last > self.idle_s:
            self.release()

    def release(self) -> None:
        """End the activity, if one is held.  PyMOL may nap again."""
        try:
            with self._lock:
                token, self._token = self._token, None
                if token is None:
                    return
                try:
                    self._runtime.end_activity(token)
                except Exception as e:
                    self._disable("ending the activity", e)
        except Exception as e:
            self.enabled = False
            print(f"MCPymol: App Nap guard off ({type(e).__name__}: {e}).")


# PyMOL 3.1 takes its whole process down, and the user's session with it, when
# these get a selection it cannot resolve: SIGABRT in ExecutiveRMS for rms,
# rms_cur and fit; SIGSEGV in CmdIntraFit for intra_rms and intra_rms_cur.
# Valid selections are fine (an empty one returns garbage, refused below too),
# and intra_fit and pair_fit report the same bad selection as a plain error (all checked against PyMOL 3.1.0, one fresh
# process per call). PyMOL's bug, but a typo must not be able to set it off.
# Checked here, in the one place every action passes through, so the tools,
# wiggles' port and any future caller are all covered.
_CRASHES_ON_UNRESOLVED_SELECTION: dict[str, tuple[str, ...]] = {
    "rms": ("mobile", "target"),
    "rms_cur": ("mobile", "target"),
    "fit": ("mobile", "target"),
    "intra_rms": ("selection",),
    "intra_rms_cur": ("selection",),
}


def _refuse_crashing_selection(action: Any, args: list[Any], kwargs: dict[str, Any]) -> dict | None:
    """The error to answer instead of running ``action``, or None to run it."""
    params = _CRASHES_ON_UNRESOLVED_SELECTION.get(action)
    if params is None:
        return None
    for position, name in enumerate(params):
        selection = args[position] if position < len(args) else kwargs.get(name)
        if selection is None:
            continue  # PyMOL's own default, which resolves
        which = "selection" if name == "selection" else f"{name} selection"
        try:
            atoms = cmd.count_atoms(str(selection))
        except Exception as e:
            return {
                "status": "error",
                "error": (
                    f"Error: {action} was not run — the {which} '{selection}' is not valid: "
                    f"{_error_text(e)} (Checked first, because PyMOL crashes on an invalid "
                    "selection here.)"
                ),
            }
        if atoms == 0:
            # Does not crash, but answers with a meaningless number: an "RMSD"
            # of 1.4e-45 between two empty selections.
            return {
                "status": "error",
                "error": f"Error: {action} was not run — the {which} '{selection}' matches no atoms.",
            }
    return None


class PyMOLSocketServer:
    """Tiny TCP/JSON server that dispatches requests to ``pymol.cmd``."""

    def __init__(
        self,
        host: str = "127.0.0.1",
        port: int | None = None,
        replaces: threading.Thread | None = None,
    ) -> None:
        self.host = host
        self.port = port if port is not None else MCP_PORT
        self.running = False
        self.thread: threading.Thread | None = None
        self.server_socket: socket.socket | None = None
        self.app_nap = _AppNapGuard()
        # Which start() a listener belongs to.  One left over from before a
        # stop/start pair finishes the request it is on and exits, rather than
        # serving on beside its successor.
        self._generation = 0
        # True while the listener is waiting for PyMOL's main thread to run a
        # `do` — which is when stop() must not wait for it (see stop).
        self._in_do = False
        # The listener of the server this one replaces (see start_plugin): it
        # may still hold the port, answering its last request.
        self._replaces = replaces

    # ── Request dispatch ────────────────────────────────────────────────
    def handle_request(self, payload: object) -> dict:
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

        if action == "do":
            command_str = args[0] if args else kwargs.get("command", "")
            self._in_do = True
            try:
                return run_do(command_str)
            except Exception as e:
                print(f"MCPymol Plugin Error executing do: {e}")
                _print_traceback(e)
                return {
                    "status": "error",
                    "error": f"PyMOL execution error: {_error_text(e)}",
                    "output": "",
                }
            finally:
                self._in_do = False

        with _OutputCapture() as captured:
            response = self._dispatch(action, args, kwargs)
        problems = _listed(_problem_lines(captured.text))
        if problems:
            if response.get("status") == "success":
                response["warnings"] = problems
            else:
                # align/super on a bad selection raise CmdException(" Error: ")
                # and print the reason only as C-level feedback.
                error = str(response.get("error", ""))
                extra = [p for p in problems if p not in error]
                if extra:
                    response["error"] = error.rstrip() + "\n" + "\n".join(extra)
        return response

    def _dispatch(self, action: Any, args: list[Any], kwargs: dict[str, Any]) -> dict:
        """Run one action other than `do`.  Exceptions become error responses."""
        refused = _refuse_crashing_selection(action, args, kwargs)
        if refused is not None:
            return refused
        try:
            # ── Special, custom handlers ─────────────────────────────────
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
            print(f"MCPymol Plugin Error executing {action}: {e}")
            _print_traceback(e)
            return {"status": "error", "error": f"PyMOL execution error: {_error_text(e)}"}

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
    def _serving(self, generation: int | None) -> bool:
        return self.running and (generation is None or generation == self._generation)

    def _listen_loop(
        self, generation: int | None = None, previous: threading.Thread | None = None
    ) -> None:
        if previous is not None and previous.is_alive():
            # The listener this one takes over from is still answering a
            # request, and holds the port until it has replied.  Often that
            # request is the very `do` of "run plugin.py" (or "stop_mcp" then
            # "start_mcp") that is replacing it: it runs on PyMOL's main thread
            # while that listener waits for it, so the port cannot be free
            # before the replacement starts.  Binding now failed, and left no
            # listener at all once the old one exited.
            print(
                f"MCPymol: waiting for the previous listener to finish its request "
                f"before listening on {self.host}:{self.port}."
            )
            previous.join()
        if not self._serving(generation):
            return

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.server_socket = sock
        try:
            try:
                sock.bind((self.host, self.port))
                sock.listen(5)
            except OSError as e:
                print(f"MCPymol Failed to start socket server: {e}")
                # Not left looking "running": start_mcp must be able to retry.
                if self._serving(generation):
                    self.running = False
                return
            # Bounds how long stop() waits: the loop only notices `running` went
            # false when accept() returns. Closing the socket from the stopping
            # thread would be faster but does not reliably wake accept() on
            # Linux, and risks an fd-reuse race — a short poll is simpler and
            # costs a few no-op wakeups per second on an idle session.
            sock.settimeout(_ACCEPT_POLL_INTERVAL)
            # Report the port actually bound, which differs from self.port when
            # 0 was requested (ephemeral).
            self.port = sock.getsockname()[1]
            print(f"MCPymol Native Plugin listening on {self.host}:{self.port}")

            while self._serving(generation):
                try:
                    client, _addr = sock.accept()
                except TimeoutError:
                    self.app_nap.tick()
                    continue
                except OSError:  # socket closed during shutdown
                    break

                # Before serving, so the request itself runs un-napped.
                self.app_nap.touch()
                try:
                    self.serve_connection(client)
                except Exception as e:
                    print(f"MCPymol Socket Server Error: {e}")
                finally:
                    try:
                        client.close()
                    except OSError:
                        pass
                    # Idle time counts from the end of the last request, not
                    # its start: a 20-minute render is not 20 minutes of idle.
                    self.app_nap.touch()

        except Exception as e:
            print(f"MCPymol Socket Server Error: {e}")
        finally:
            self.app_nap.release()
            try:
                sock.close()
            except OSError:
                pass

    def start(self) -> None:
        if self.running:
            return
        self.running = True
        self._generation += 1
        previous, self._replaces = self._replaces, None
        if self.thread is not None and self.thread.is_alive():
            # stop() then start() while the old listener is still finishing a
            # request: it exits once it has replied, and hands over the port.
            previous = self.thread
        self.thread = threading.Thread(
            target=self._listen_loop, args=(self._generation, previous), daemon=True
        )
        self.thread.start()
        print("MCPymol Bridge started.")

    def stop(self) -> None:
        if not self.running:
            return
        self.running = False
        thread = self.thread
        # Waiting would be pointless when this is PyMOL's main thread running a
        # `do` the listener is itself waiting on — a "stop_mcp", or a re-run of
        # plugin.py, sent through the bridge: the listener can only finish once
        # this returns, so a join just stalls PyMOL for its whole timeout.  It
        # exits by itself after replying.  (And a listener cannot join itself.)
        waiting_on_us = self._in_do and _on_pymol_main_thread()
        if thread is not None and thread is not threading.current_thread() and not waiting_on_us:
            thread.join(timeout=2.0)
        self.app_nap.release()
        print("MCPymol Bridge stopped.")


def start_plugin() -> PyMOLSocketServer | None:
    """Start the bridge singleton and register the PyMOL commands.

    Idempotent across ``run plugin.py`` re-runs: a server left over from a
    previous run is stopped, and the new one takes the port over once the old
    listener has finished whatever request it is on.
    """
    global mcp_bridge_plugin
    try:
        existing = globals().get("mcp_bridge_plugin")
        replaces = None
        if existing is not None:
            existing.stop()
            replaces = getattr(existing, "thread", None)
        mcp_bridge_plugin = PyMOLSocketServer(replaces=replaces)
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
