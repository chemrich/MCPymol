"""Rendering the scene into something the model can actually look at.

The primitives in :mod:`mcpymol.primitives` expose PyMOL's ``ray`` and ``png``
separately, which leaves the caller holding a filename it cannot see.  These
tools close that loop: they render and hand the image straight back as MCP
image content.
"""

import os
import tempfile
import time
from typing import Annotated

from mcp.server.fastmcp import Image
from pydantic import Field

from mcpymol.app import mcp
from mcpymol.bridge import _SLOW_OP_TIMEOUT, send_request

# Ray-tracing is not cheap and the result is inlined into the conversation, so
# the default is a size that reads well without dominating the context.  A
# 1000x750 trace of a typical structure lands around 300-700 KB.
DEFAULT_RENDER_WIDTH = 1000
DEFAULT_RENDER_HEIGHT = 750

# Images are base64-encoded into the response, inflating them by a third.
# Beyond this, hand back the path instead of the pixels.
MAX_INLINE_IMAGE_BYTES = int(os.environ.get("MCPYMOL_MAX_IMAGE_BYTES", 5_000_000))

# How long to wait for the PNG to appear on disk after PyMOL reports success.
_PNG_APPEAR_TIMEOUT = 60.0

# Every PNG ends with a zero-length IEND chunk; its presence means the file is
# complete rather than still being written.
_PNG_EOF = b"IEND\xae\x42\x60\x82"

# PyMOL's OpenGL frame grab (cmd.png with ray=0) needs the GUI thread, and the
# plugin dispatches on a socket worker thread — so over this bridge it writes a
# blank image rather than failing. The unshaded fast path is therefore not
# reachable here at all, whatever the caller asks for.
_OPENGL_UNAVAILABLE_NOTE = (
    "Note: the fast unshaded path (ray_trace=False) writes a blank image over "
    "the MCPymol bridge — PyMOL's OpenGL capture needs its GUI thread, and the "
    "plugin runs off it. Ray-traced instead."
)

# Below this many distinct bytes in the decompressed pixel stream, an image is
# effectively one flat colour. A heuristic, but a decisive one: real renders
# are nowhere near it.
_BLANK_IMAGE_MAX_DISTINCT_BYTES = 8


def _looks_blank(png: bytes) -> bool:
    """True if a PNG is essentially a single flat colour.

    Catches the two ways a render comes back empty: the OpenGL path writing
    nothing, and a scene where everything is hidden or the selection matched
    no atoms. Counts distinct bytes in the decompressed pixel stream — for a
    uniform image the filtered rows are uniform too, so the count stays tiny.
    """
    import struct
    import zlib

    idat, offset = b"", 8
    try:
        while offset + 8 <= len(png):
            length = struct.unpack(">I", png[offset : offset + 4])[0]
            if png[offset + 4 : offset + 8] == b"IDAT":
                idat += png[offset + 8 : offset + 8 + length]
            offset += 12 + length
        if not idat:
            return False
        raw = zlib.decompress(idat)
    except (struct.error, zlib.error):
        return False
    return len(set(raw)) <= _BLANK_IMAGE_MAX_DISTINCT_BYTES


def _read_complete_png(path: str, timeout: float | None = None) -> bytes | None:
    """Return the PNG's bytes once it is fully written, or None on timeout.

    PyMOL writes the file from its own thread, so it can lag the socket
    response that reports success.  Waiting for the IEND terminator rather than
    mere existence avoids handing back a truncated image.

    ``timeout`` defaults to the current ``_PNG_APPEAR_TIMEOUT``, resolved at
    call time so the constant stays overridable.
    """
    deadline = time.monotonic() + (_PNG_APPEAR_TIMEOUT if timeout is None else timeout)
    while time.monotonic() < deadline:
        try:
            with open(path, "rb") as fh:
                data = fh.read()
        except OSError:
            data = b""
        if data.endswith(_PNG_EOF):
            return data
        time.sleep(0.05)
    return None


# structured_output=False: the return is either image content or an error
# string, which is not a shape pydantic can build an output model for.
@mcp.tool(structured_output=False)
def render(
    width: Annotated[
        int, Field(description="Image width in pixels. Larger costs render time and context.")
    ] = DEFAULT_RENDER_WIDTH,
    height: Annotated[int, Field(description="Image height in pixels.")] = DEFAULT_RENDER_HEIGHT,
    ray_trace: Annotated[
        bool,
        Field(
            description="Ray-trace for publication quality (default), or take a fast unshaded viewport grab."
        ),
    ] = True,
    filename: Annotated[
        str | None,
        Field(
            description="Optional path to also keep the PNG at. Without it the render goes to a temporary file that is cleaned up afterwards."
        ),
    ] = None,
) -> Image | str:
    """
    Renders the current scene and returns the image, so you can see it.

    Use this instead of calling ``ray`` and ``png`` separately — those leave
    you holding a filename you cannot look at. Call this after setting up a
    view to check what it actually looks like, and iterate.

    Ray-tracing gives shadows and smooth surfaces but takes seconds to minutes
    on a large assembly; set ``ray_trace=False`` for a fast, flat OpenGL
    snapshot when you only need to confirm a selection or orientation.
    """
    if width < 1 or height < 1:
        return f"Error: width and height must be positive, got {width}x{height}."

    note = "" if ray_trace else f" {_OPENGL_UNAVAILABLE_NOTE}"
    ray_trace = True  # see _OPENGL_UNAVAILABLE_NOTE

    if filename is not None:
        path = os.path.abspath(os.path.expanduser(filename))
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    else:
        fd, path = tempfile.mkstemp(suffix=".png", prefix="mcpymol_render_")
        os.close(fd)
        # PyMOL refuses to overwrite in some builds; give it a clean slate.
        os.unlink(path)
    keep = filename is not None

    try:
        res = send_request(
            "png",
            args=[path],
            kwargs={"width": width, "height": height, "ray": 1 if ray_trace else 0},
            timeout=_SLOW_OP_TIMEOUT,
        )
        if res.get("status") == "error":
            return f"Error rendering: {res.get('error')}"

        data = _read_complete_png(path)
        if data is None:
            return (
                f"PyMOL reported success but no complete PNG appeared at {path} "
                f"within {_PNG_APPEAR_TIMEOUT:.0f}s. If PyMOL is running on a "
                f"different machine than this bridge, they cannot exchange files."
            )

        if _looks_blank(data):
            return (
                f"Rendered {width}x{height}, but the image is a single flat colour — "
                f"nothing was visible. Check that an object is loaded and enabled "
                f"(list_objects), and that your selection matches atoms "
                f"(count_atoms).{note}"
            )

        if len(data) > MAX_INLINE_IMAGE_BYTES:
            saved = path if keep else "a temporary file"
            return (
                f"Rendered {width}x{height} ({len(data) / 1e6:.1f} MB) — too large to "
                f"return inline (limit {MAX_INLINE_IMAGE_BYTES / 1e6:.1f} MB). It is at "
                f"{saved}. Re-run with a smaller width/height to see it here."
            )

        return Image(data=data, format="png")
    finally:
        if not keep:
            try:
                os.unlink(path)
            except OSError:
                pass


@mcp.tool()
def turntable(
    obj_name: Annotated[
        str, Field(description="Object to centre the rotation on (default: the whole scene).")
    ] = "all",
    frames: Annotated[
        int, Field(description="Number of frames over the full 360°. 36 gives 10° steps.")
    ] = 36,
    out_dir: Annotated[str, Field(description="Directory to write the PNG sequence into.")] = ".",
    prefix: Annotated[
        str, Field(description="Filename prefix; frames are ``<prefix>_0000.png`` and up.")
    ] = "turntable",
    width: Annotated[int, Field(description="Frame width in pixels.")] = 800,
    height: Annotated[int, Field(description="Frame height in pixels.")] = 600,
    ray_trace: Annotated[
        bool,
        Field(
            description="Kept for compatibility and ignored — every frame is "
            "ray-traced, because PyMOL's unshaded OpenGL capture does not work "
            "over this bridge."
        ),
    ] = True,
) -> str:
    """
    Renders a full 360° rotation as a numbered PNG sequence.

    Spins the camera around the vertical axis in equal steps and writes one
    frame per step, ready to assemble into a GIF or MP4 (e.g.
    ``ffmpeg -i turntable_0000.png out.mp4``).

    Every frame is ray-traced, which is slow — 36 frames of a large assembly
    can take an hour, so start with few frames and a small width/height. The
    unshaded OpenGL renderer would be far faster but does not work over this
    bridge (it needs PyMOL's GUI thread), and silently produced blank frames
    until this was found.
    """
    if frames < 2:
        return f"Error: frames must be at least 2, got {frames}."

    note = "" if ray_trace else f" {_OPENGL_UNAVAILABLE_NOTE}"
    ray_trace = True  # see _OPENGL_UNAVAILABLE_NOTE

    out_dir = os.path.abspath(os.path.expanduser(out_dir))
    os.makedirs(out_dir, exist_ok=True)

    # Rotate about the object's own centre, or the turntable wobbles.
    send_request("do", args=[f"origin {obj_name}"])
    send_request("center", args=[obj_name])

    step = 360.0 / frames
    written = []
    for i in range(frames):
        if i:
            turn = send_request("turn", args=["y", str(step)])
            if turn.get("status") == "error":
                return f"Error rotating at frame {i}: {turn.get('error')}"

        path = os.path.join(out_dir, f"{prefix}_{i:04d}.png")
        res = send_request(
            "png",
            args=[path],
            kwargs={"width": width, "height": height, "ray": 1 if ray_trace else 0},
            timeout=_SLOW_OP_TIMEOUT,
        )
        if res.get("status") == "error":
            return f"Error writing frame {i}: {res.get('error')}"
        written.append(path)

    return (
        f"Wrote {len(written)} ray-traced frames ({step:.1f}° apart) to {out_dir} as "
        f"{prefix}_0000.png … {prefix}_{frames - 1:04d}.png. Assemble with:\n"
        f"  ffmpeg -framerate 24 -i {os.path.join(out_dir, prefix)}_%04d.png "
        f"-pix_fmt yuv420p {prefix}.mp4{note}"
    )
