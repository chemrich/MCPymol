"""Tests for render() and turntable()."""

import os
import struct
import zlib
from unittest.mock import patch

import pytest
from mcp.server.fastmcp import Image

from mcpymol.rendering import _PNG_EOF, _read_complete_png, render, turntable


def _png_bytes(width=2, height=2):
    """A real, minimal, valid PNG — not a stub, so the IEND check is exercised."""

    def chunk(tag, data):
        return (
            struct.pack(">I", len(data))
            + tag
            + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
        )

    raw = b"".join(b"\x00" + b"\xff\x00\x00" * width for _ in range(height))
    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw))
        + chunk(b"IEND", b"")
    )


def _writes_png(path_holder, data=None):
    """send_request side effect that writes a PNG where PyMOL was asked to."""

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "png" and args:
            path_holder.append(args[0])
            with open(args[0], "wb") as fh:
                fh.write(data if data is not None else _png_bytes())
        return {"status": "success", "result": "OK"}

    return fake


# ── _read_complete_png ───────────────────────────────────────────────────────


def test_read_complete_png_returns_bytes(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(_png_bytes())

    assert _read_complete_png(str(p)) == _png_bytes()


def test_read_complete_png_rejects_a_truncated_file(tmp_path):
    """A half-written PNG must not be handed back as an image."""
    p = tmp_path / "x.png"
    p.write_bytes(_png_bytes()[:-20])  # no IEND yet

    assert _read_complete_png(str(p), timeout=0.2) is None


def test_read_complete_png_times_out_when_absent(tmp_path):
    assert _read_complete_png(str(tmp_path / "never.png"), timeout=0.2) is None


def test_png_eof_constant_matches_a_real_png():
    assert _png_bytes().endswith(_PNG_EOF)


# ── render ───────────────────────────────────────────────────────────────────


@patch("mcpymol.rendering.send_request")
def test_render_returns_image_content(mock_sr):
    """The whole point: the model gets pixels back, not a filename."""
    mock_sr.side_effect = _writes_png([])

    result = render()

    assert isinstance(result, Image)
    assert result.data == _png_bytes()
    assert result._mime_type == "image/png"


@patch("mcpymol.rendering.send_request")
def test_render_asks_pymol_to_ray_trace_at_the_requested_size(mock_sr):
    mock_sr.side_effect = _writes_png([])

    render(width=640, height=480)

    kwargs = mock_sr.call_args.kwargs["kwargs"]
    assert kwargs == {"width": 640, "height": 480, "ray": 1}
    assert mock_sr.call_args.args[0] == "png"


@patch("mcpymol.rendering.send_request")
def test_render_can_skip_ray_tracing(mock_sr):
    mock_sr.side_effect = _writes_png([])

    render(ray_trace=False)

    assert mock_sr.call_args.kwargs["kwargs"]["ray"] == 0


@patch("mcpymol.rendering.send_request")
def test_render_uses_the_slow_timeout(mock_sr):
    """A 1920x1080 trace does not finish inside the 10 s interactive budget."""
    from mcpymol.bridge import _SLOW_OP_TIMEOUT

    mock_sr.side_effect = _writes_png([])

    render()

    assert mock_sr.call_args.kwargs["timeout"] == _SLOW_OP_TIMEOUT


@patch("mcpymol.rendering.send_request")
def test_render_cleans_up_its_temp_file(mock_sr):
    seen: list[str] = []
    mock_sr.side_effect = _writes_png(seen)

    render()

    assert seen, "png was never requested"
    assert not os.path.exists(seen[0]), "temporary render was left behind"


@patch("mcpymol.rendering.send_request")
def test_render_keeps_a_named_file(mock_sr, tmp_path):
    target = tmp_path / "nested" / "shot.png"
    mock_sr.side_effect = _writes_png([])

    result = render(filename=str(target))

    assert isinstance(result, Image)
    assert target.exists(), "named output should be kept"


@patch("mcpymol.rendering.send_request")
def test_render_reports_pymol_errors(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no objects to render"}

    result = render()

    assert isinstance(result, str)
    assert "no objects to render" in result


@patch("mcpymol.rendering.send_request")
def test_render_reports_a_missing_file(mock_sr):
    """PyMOL claiming success without writing means the halves are split
    across machines — say so rather than returning a broken image."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    with patch("mcpymol.rendering._PNG_APPEAR_TIMEOUT", 0.2):
        result = render()

    assert isinstance(result, str)
    assert "no complete PNG appeared" in result


@patch("mcpymol.rendering.send_request")
def test_render_refuses_to_inline_a_huge_image(mock_sr, tmp_path):
    """Base64 inflates by a third; a 40 MB render would swamp the context."""
    big = _png_bytes(width=2, height=2)
    mock_sr.side_effect = _writes_png([], data=big)

    with patch("mcpymol.rendering.MAX_INLINE_IMAGE_BYTES", 10):
        result = render(filename=str(tmp_path / "big.png"))

    assert isinstance(result, str)
    assert "too large to return inline" in result
    assert "smaller width/height" in result


@pytest.mark.parametrize("w,h", [(0, 100), (100, 0), (-5, 100)])
@patch("mcpymol.rendering.send_request")
def test_render_rejects_nonsense_dimensions(mock_sr, w, h):
    result = render(width=w, height=h)

    assert isinstance(result, str)
    assert "must be positive" in result
    mock_sr.assert_not_called()


# ── turntable ────────────────────────────────────────────────────────────────


@patch("mcpymol.rendering.send_request")
def test_turntable_writes_one_frame_per_step(mock_sr, tmp_path):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = turntable(obj_name="1abc", frames=4, out_dir=str(tmp_path))

    png_calls = [c for c in mock_sr.call_args_list if c.args[0] == "png"]
    assert len(png_calls) == 4
    assert [os.path.basename(c.kwargs["args"][0]) for c in png_calls] == [
        "turntable_0000.png",
        "turntable_0001.png",
        "turntable_0002.png",
        "turntable_0003.png",
    ]
    assert "Wrote 4 OpenGL frames" in result


@patch("mcpymol.rendering.send_request")
def test_turntable_rotates_a_full_circle(mock_sr, tmp_path):
    """N frames must span exactly 360°, with no rotation before the first."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    turntable(frames=8, out_dir=str(tmp_path))

    turns = [c for c in mock_sr.call_args_list if c.args[0] == "turn"]
    assert len(turns) == 7, "the first frame is the starting orientation"
    steps = [float(c.kwargs["args"][1]) for c in turns]
    assert all(s == 45.0 for s in steps)
    assert sum(steps) == 315.0  # the 8th step would close the loop


@patch("mcpymol.rendering.send_request")
def test_turntable_sets_the_rotation_origin(mock_sr, tmp_path):
    """Without this the model orbits the scene centre and appears to wobble."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    turntable(obj_name="1abc", frames=2, out_dir=str(tmp_path))

    dos = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert "origin 1abc" in dos


@patch("mcpymol.rendering.send_request")
def test_turntable_defaults_to_the_fast_renderer(mock_sr, tmp_path):
    """36 ray-traced frames of a big assembly can take an hour."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    turntable(frames=3, out_dir=str(tmp_path))

    for c in mock_sr.call_args_list:
        if c.args[0] == "png":
            assert c.kwargs["kwargs"]["ray"] == 0


@patch("mcpymol.rendering.send_request")
def test_turntable_can_ray_trace(mock_sr, tmp_path):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = turntable(frames=3, out_dir=str(tmp_path), ray_trace=True)

    assert all(c.kwargs["kwargs"]["ray"] == 1 for c in mock_sr.call_args_list if c.args[0] == "png")
    assert "ray-traced" in result


@patch("mcpymol.rendering.send_request")
def test_turntable_stops_on_a_write_error(mock_sr, tmp_path):
    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "png":
            return {"status": "error", "error": "disk full"}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = turntable(frames=5, out_dir=str(tmp_path))

    assert "Error writing frame 0" in result
    assert "disk full" in result


@patch("mcpymol.rendering.send_request")
def test_turntable_rejects_too_few_frames(mock_sr):
    assert "at least 2" in turntable(frames=1)
    mock_sr.assert_not_called()


@patch("mcpymol.rendering.send_request")
def test_turntable_creates_the_output_directory(mock_sr, tmp_path):
    mock_sr.return_value = {"status": "success", "result": "OK"}
    target = tmp_path / "frames" / "run1"

    turntable(frames=2, out_dir=str(target))

    assert target.is_dir()
