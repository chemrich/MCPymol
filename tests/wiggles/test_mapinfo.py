"""Tests for wiggles.mapinfo.

Fixtures are synthesised rather than checked in as binary blobs: a 1024-byte
header is small enough to build honestly, and building it makes the field
offsets explicit, which is the part most likely to be wrong.
"""

from __future__ import annotations

import gzip
import struct

import pytest

from mcpymol.wiggles.mapinfo import ISOTROPY_RTOL, map_info, read_map_header

HEADER_SIZE = 1024


def make_header(
    *,
    nx: int = 256,
    ny: int = 256,
    nz: int = 256,
    mode: int = 2,
    mx: int | None = None,
    my: int | None = None,
    mz: int | None = None,
    cella: tuple[float, float, float] | None = None,
    mapc: int = 1,
    mapr: int = 2,
    maps: int = 3,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    magic: bytes = b"MAP ",
    machst: bytes = b"\x44\x44\x00\x00",
    endian: str = "<",
    voxel: float = 1.0,
    nxstart: int = 0,
    nystart: int = 0,
    nzstart: int = 0,
    dmin: float = -1.0,
    dmax: float = 1.0,
    dmean: float = 0.0,
    rms: float = 0.5,
) -> bytes:
    """Build a 1024-byte MRC2014 header.

    By default the grid sampling matches the extent and the cell is sized so
    that ``cella/m == voxel`` exactly.

    The density statistics are parameterised because a *local-resolution* map
    stores Ångström rather than density, so its dmin/dmax/dmean/rms are the
    input to the ramp conversion rather than incidental — see
    ``test_localres.py``.
    """
    mx = nx if mx is None else mx
    my = ny if my is None else my
    mz = nz if mz is None else mz
    if cella is None:
        cella = (mx * voxel, my * voxel, mz * voxel)

    buf = bytearray(HEADER_SIZE)

    def put(fmt: str, off: int, *vals: object) -> None:
        struct.pack_into(endian + fmt, buf, off, *vals)

    put("3i", 0, nx, ny, nz)
    put("i", 12, mode)
    put("3i", 16, nxstart, nystart, nzstart)
    put("3i", 28, mx, my, mz)
    put("3f", 40, *cella)
    put("3f", 52, 90.0, 90.0, 90.0)  # cellb
    put("3i", 64, mapc, mapr, maps)
    put("3f", 76, dmin, dmax, dmean)
    put("i", 88, 1)  # ispg
    put("i", 92, 0)  # nsymbt
    put("i", 108, 20140)  # nversion
    put("3f", 196, *origin)
    buf[208:212] = magic
    buf[212:216] = machst
    put("f", 216, rms)
    put("i", 220, 0)  # nlabl
    return bytes(buf)


def write_map(tmp_path, name="test.mrc", *, gz: bool = False, **kw):
    data = make_header(**kw) + b"\x00" * 64  # header + a little dummy payload
    path = tmp_path / (name + (".gz" if gz else ""))
    if gz:
        path.write_bytes(gzip.compress(data))
    else:
        path.write_bytes(data)
    return path


# -- the central correctness claim ----------------------------------------


def test_voxel_size_uses_grid_sampling_not_extent(tmp_path):
    """cella/m, not cella/n. This is the whole point of the module.

    A cropped map: 100 columns stored, but the grid it was sampled on is 200,
    and the cell is 200 Å. The true spacing is 200/200 = 1.0 Å. Dividing by
    the extent would give 200/100 = 2.0 Å — plausible, and wrong by 2x.
    """
    p = write_map(
        tmp_path, nx=100, ny=100, nz=100, mx=200, my=200, mz=200, cella=(200.0, 200.0, 200.0)
    )
    h = read_map_header(p)
    assert h.voxel_size == pytest.approx((1.0, 1.0, 1.0))
    # the wrong answer, stated explicitly so the test documents the trap
    assert h.cella[0] / h.nx == pytest.approx(2.0)


def test_cropped_map_is_flagged(tmp_path):
    p = write_map(
        tmp_path, nx=100, ny=100, nz=100, mx=200, my=200, mz=200, cella=(200.0, 200.0, 200.0)
    )
    warnings = read_map_header(p).warnings()
    assert any("boxed or cropped" in w for w in warnings)


# -- isotropy, with the EMD-30913 regression -------------------------------


def test_isotropic_map_has_no_anisotropy_warning(tmp_path):
    p = write_map(tmp_path, voxel=1.0)
    h = read_map_header(p)
    assert h.is_isotropic is True
    assert not any("ANISOTROPIC" in w for w in h.warnings())


def test_float_noise_spacing_is_still_isotropic(tmp_path):
    """EMD-30913 reports 0.7999967, not 0.8. Exact equality would call this
    anisotropic; it is not."""
    spacing = 0.7999967
    n = 110
    p = write_map(tmp_path, nx=n, ny=n, nz=n, cella=(n * spacing, n * spacing, n * 0.8))
    h = read_map_header(p)
    assert h.is_isotropic is True, f"spacings {h.voxel_size} wrongly flagged"
    assert not any("ANISOTROPIC" in w for w in h.warnings())


def test_genuinely_anisotropic_is_flagged_loudly(tmp_path):
    p = write_map(tmp_path, nx=100, ny=100, nz=100, cella=(100.0, 100.0, 150.0))  # Z spacing 1.5x
    h = read_map_header(p)
    assert h.is_isotropic is False
    assert any("ANISOTROPIC" in w for w in h.warnings())
    assert "ANISOTROPIC" in map_info(p)


def test_isotropy_boundary(tmp_path):
    """Just inside tolerance passes; comfortably outside fails."""
    n = 100
    inside = 1.0 + ISOTROPY_RTOL / 2
    p = write_map(tmp_path, nx=n, ny=n, nz=n, cella=(100.0, 100.0, n * inside))
    assert read_map_header(p).is_isotropic is True

    outside = 1.0 + ISOTROPY_RTOL * 100
    q = write_map(tmp_path, "b.mrc", nx=n, ny=n, nz=n, cella=(100.0, 100.0, n * outside))
    assert read_map_header(q).is_isotropic is False


# -- malformed headers should degrade, not explode -------------------------


def test_zero_cell_does_not_divide_by_zero(tmp_path):
    p = write_map(tmp_path, cella=(0.0, 0.0, 0.0))
    h = read_map_header(p)
    assert h.voxel_size == (None, None, None)
    assert h.is_isotropic is None
    assert any("undefined" in w for w in h.warnings())
    assert "unknown" in map_info(p)  # renders without raising


def test_zero_sampling_does_not_divide_by_zero(tmp_path):
    p = write_map(tmp_path, mx=0, my=0, mz=0, cella=(100.0, 100.0, 100.0))
    assert read_map_header(p).voxel_size == (None, None, None)


def test_bad_magic_is_warned_not_fatal(tmp_path):
    p = write_map(tmp_path, magic=b"JUNK")
    h = read_map_header(p)
    assert any("magic" in w for w in h.warnings())
    assert h.voxel_size == pytest.approx((1.0, 1.0, 1.0))  # still readable


def test_truncated_file_raises(tmp_path):
    p = tmp_path / "short.mrc"
    p.write_bytes(b"\x00" * 100)
    with pytest.raises(ValueError, match="too short"):
        read_map_header(p)


def test_nonsense_header_raises(tmp_path):
    p = tmp_path / "junk.mrc"
    p.write_bytes(b"\xff" * HEADER_SIZE)
    with pytest.raises(ValueError, match=r"nonsense|too short"):
        read_map_header(p)


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_map_header(tmp_path / "nope.mrc")


# -- axis order ------------------------------------------------------------


def test_default_axis_order_is_not_warned(tmp_path):
    p = write_map(tmp_path, mapc=1, mapr=2, maps=3)
    h = read_map_header(p)
    assert h.axis_mapping == ("X", "Y", "Z")
    assert not any("axis order" in w for w in h.warnings())


def test_permuted_axis_order_is_warned(tmp_path):
    p = write_map(tmp_path, mapc=3, mapr=1, maps=2)
    h = read_map_header(p)
    assert h.axis_mapping == ("Z", "X", "Y")
    assert any("non-default axis order" in w for w in h.warnings())


def test_malformed_axis_order_is_warned(tmp_path):
    p = write_map(tmp_path, mapc=1, mapr=1, maps=1)
    assert any("malformed" in w for w in read_map_header(p).warnings())


# -- byte order and compression -------------------------------------------


def test_big_endian_header(tmp_path):
    p = write_map(tmp_path, endian=">", machst=b"\x11\x11\x00\x00", nx=64, ny=64, nz=64, voxel=2.0)
    h = read_map_header(p)
    assert h.byte_order == "big"
    assert (h.nx, h.ny, h.nz) == (64, 64, 64)
    assert h.voxel_size == pytest.approx((2.0, 2.0, 2.0))


def test_junk_machine_stamp_falls_back_to_inference(tmp_path):
    """Real files carry junk stamps; dimensions should still parse."""
    p = write_map(tmp_path, machst=b"\x00\x00\x00\x00", nx=64, ny=64, nz=64)
    h = read_map_header(p)
    assert (h.nx, h.ny, h.nz) == (64, 64, 64)


def test_gzipped_map_is_read_transparently(tmp_path):
    """EMDB ships .map.gz."""
    p = write_map(tmp_path, "emd.map", gz=True, nx=64, ny=64, nz=64, voxel=1.5)
    h = read_map_header(p)
    assert h.voxel_size == pytest.approx((1.5, 1.5, 1.5))


# -- the report ------------------------------------------------------------


def test_report_leads_with_voxel_size(tmp_path):
    p = write_map(tmp_path, voxel=0.5332)
    out = map_info(p)
    assert out.index("Voxel size") < out.index("Grid")
    assert "0.5332" in out


def test_report_states_no_warnings_when_clean(tmp_path):
    assert "Warnings: none" in map_info(write_map(tmp_path))


def test_emd_11638_like_apoferritin(tmp_path):
    """The 1.22 Å apoferritin fixture from the compendium: 0.5332 Å, 256³."""
    p = write_map(tmp_path, nx=256, ny=256, nz=256, voxel=0.5332)
    h = read_map_header(p)
    assert h.voxel_size == pytest.approx((0.5332, 0.5332, 0.5332), rel=1e-6)
    assert h.is_isotropic is True
    assert h.warnings() == []
