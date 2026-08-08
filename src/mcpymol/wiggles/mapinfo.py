"""Read an MRC/CCP4 map header and report what it actually says.

The point of this module is one number that nothing displays: **voxel size**.
It is not stored in the header — it is derived — and the derivation has a trap
that yields a plausible-looking wrong answer on any cropped map. Meanwhile the
nominal value it encodes is only expected to be accurate to ±5–15% in the first
place (Wu, Lander & Herzik, *J Struct Biol X* 4:100020, 2020), which at 1.2 Å
is a systematic stretch of every distance in the model that looks like a
slightly strained structure rather than an error.

So: parse the header, compute the spacing correctly, and say it out loud.

No dependencies, no network, no PyMOL. Gzipped maps are handled because that is
how EMDB ships them.

Header layout is MRC2014 (Cheng et al., *J Struct Biol* 192(2):146–150, 2015).
"""

from __future__ import annotations

import gzip
import struct
from dataclasses import dataclass
from pathlib import Path

HEADER_SIZE = 1024

# Relative tolerance for calling two axis spacings equal. EMD-30913 reports
# 0.7999967 rather than 0.8; exact equality would flag a perfectly isotropic
# map as anisotropic.
ISOTROPY_RTOL = 1e-4

# MRC mode -> human label. Unknown modes are reported rather than rejected;
# a mode we do not recognise does not stop the geometry from being readable.
_MODES = {
    0: "int8",
    1: "int16",
    2: "float32",
    3: "complex int16",
    4: "complex float32",
    6: "uint16",
    12: "float16",
    101: "4-bit",
}

_AXIS_NAMES = {1: "X", 2: "Y", 3: "Z"}


@dataclass(frozen=True)
class MapHeader:
    """Everything from an MRC/CCP4 header that bears on geometry.

    Field names follow the format specification rather than being prettified,
    so that comparing this against the spec — or against another parser — is
    a direct read.
    """

    path: str
    byte_order: str  # "little" | "big"
    # Words 1-3: extent of the stored array, in columns/rows/sections.
    nx: int
    ny: int
    nz: int
    mode: int
    # Words 5-7: location of the first column/row/section.
    nxstart: int
    nystart: int
    nzstart: int
    # Words 8-10: grid sampling along the crystallographic X/Y/Z axes.
    mx: int
    my: int
    mz: int
    # Words 11-13: cell dimensions in Angstrom, along X/Y/Z.
    cella: tuple[float, float, float]
    cellb: tuple[float, float, float]
    # Words 17-19: which crystallographic axis each of column/row/section is.
    mapc: int
    mapr: int
    maps: int
    dmin: float
    dmax: float
    dmean: float
    ispg: int
    nsymbt: int
    exttyp: str
    nversion: int
    origin: tuple[float, float, float]
    magic: str
    rms: float
    nlabl: int
    labels: tuple[str, ...]

    # -- derived -----------------------------------------------------------

    @property
    def voxel_size(self) -> tuple[float | None, float | None, float | None]:
        """Voxel size in Angstrom along crystallographic X, Y, Z.

        ``cella / m``, **not** ``cella / n``. ``mx/my/mz`` are the grid
        sampling; ``nx/ny/nz`` are the extent of the stored array. On a map
        that has been boxed or cropped these differ, and dividing by ``n``
        silently returns a wrong spacing that is the right order of magnitude
        — the worst kind of wrong.

        ``None`` for an axis whose sampling or cell dimension is zero, rather
        than raising: a header can be malformed in that way and the rest of it
        is still worth reporting.
        """
        return tuple(  # type: ignore[return-value]
            (cell / m) if (m > 0 and cell > 0) else None
            for cell, m in zip(self.cella, (self.mx, self.my, self.mz), strict=True)
        )

    @property
    def axis_mapping(self) -> tuple[str, str, str]:
        """Which crystallographic axis each of column/row/section corresponds to."""
        return tuple(  # type: ignore[return-value]
            _AXIS_NAMES.get(a, f"?({a})") for a in (self.mapc, self.mapr, self.maps)
        )

    @property
    def is_isotropic(self) -> bool | None:
        """True/False, or None if any axis spacing is unknown."""
        known = [s for s in self.voxel_size if s is not None]
        if len(known) != 3:
            return None
        lo, hi = min(known), max(known)
        return (hi - lo) <= ISOTROPY_RTOL * max(abs(hi), abs(lo))

    def warnings(self) -> list[str]:
        """Things a user should be told, loudly. Empty list means nothing odd."""
        out: list[str] = []

        if self.magic != "MAP ":
            out.append(
                f"magic string is {self.magic!r}, expected 'MAP ' — file may not be "
                f"a conformant MRC2014 map"
            )

        spacing = self.voxel_size
        unknown = [ax for ax, s in zip("XYZ", spacing, strict=True) if s is None]
        if unknown:
            out.append(
                f"voxel size undefined along {', '.join(unknown)} "
                f"(cell dimension or grid sampling is zero) — geometry is unusable"
            )
        elif self.is_isotropic is False:
            vals = ", ".join(f"{ax}={s:.6g}" for ax, s in zip("XYZ", spacing, strict=True))
            out.append(
                f"ANISOTROPIC VOXELS: {vals} Å. This is usually a header or "
                f"processing error. Every distance measured in this map is "
                f"scaled differently along each axis."
            )

        if (self.mx, self.my, self.mz) != (self.nx, self.ny, self.nz):
            out.append(
                f"grid sampling (mx,my,mz)=({self.mx},{self.my},{self.mz}) differs "
                f"from extent (nx,ny,nz)=({self.nx},{self.ny},{self.nz}) — map is "
                f"boxed or cropped. Voxel size uses m, which is correct; dividing "
                f"by n here would give a different and wrong answer"
            )

        if sorted((self.mapc, self.mapr, self.maps)) != [1, 2, 3]:
            out.append(
                f"axis mapping (mapc,mapr,maps)=({self.mapc},{self.mapr},{self.maps}) "
                f"is not a permutation of 1,2,3 — header is malformed"
            )
        elif (self.mapc, self.mapr, self.maps) != (1, 2, 3):
            out.append(
                f"non-default axis order: column/row/section map to "
                f"{'/'.join(self.axis_mapping)}. Anything assuming C-order XYZ will "
                f"transpose this map"
            )

        if self.mode not in _MODES:
            out.append(f"unrecognised data mode {self.mode}")

        return out


def _read_header_bytes(path: Path) -> bytes:
    """Read the first 1024 bytes, transparently handling gzip."""
    with open(path, "rb") as fh:
        head = fh.read(2)
    opener = gzip.open if head == b"\x1f\x8b" else open
    with opener(path, "rb") as fh:  # type: ignore[operator]
        raw = fh.read(HEADER_SIZE)
    if len(raw) < HEADER_SIZE:
        raise ValueError(
            f"{path}: file is {len(raw)} bytes, too short for a {HEADER_SIZE}-byte MRC header"
        )
    return raw


def _plausible(raw: bytes, endian: str) -> bool:
    """Do the dimensions parse as something sane under this byte order?"""
    nx, ny, nz, mode = struct.unpack_from(f"{endian}4i", raw, 0)
    return all(0 < n < 100_000 for n in (nx, ny, nz)) and -1 < mode < 200


def read_map_header(path: str | Path) -> MapHeader:
    """Parse the header of an MRC/CCP4 map.

    Byte order is taken from the machine-stamp when it is one of the documented
    values, and otherwise inferred by checking which interpretation yields sane
    dimensions — real files in the wild do carry junk stamps.

    Raises:
        FileNotFoundError: no such file.
        ValueError: too short, or the header parses as nonsense under both
            byte orders.
    """
    path = Path(path)
    raw = _read_header_bytes(path)

    machst = raw[212:216]
    if machst[:2] == b"\x44\x44" or machst[:2] == b"\x44\x41":
        endian, order = "<", "little"
    elif machst[:2] == b"\x11\x11":
        endian, order = ">", "big"
    elif _plausible(raw, "<"):
        endian, order = "<", "little"
    elif _plausible(raw, ">"):
        endian, order = ">", "big"
    else:
        raise ValueError(
            f"{path}: header parses as nonsense under both byte orders — not an MRC/CCP4 map?"
        )

    # Trust the stamp only if it also produces sane numbers.
    if not _plausible(raw, endian):
        other = ">" if endian == "<" else "<"
        if _plausible(raw, other):
            endian, order = other, ("big" if other == ">" else "little")

    u = lambda fmt, off: struct.unpack_from(endian + fmt, raw, off)  # noqa: E731

    nx, ny, nz, mode = u("4i", 0)
    nxstart, nystart, nzstart = u("3i", 16)
    mx, my, mz = u("3i", 28)
    cella = u("3f", 40)
    cellb = u("3f", 52)
    mapc, mapr, maps = u("3i", 64)
    dmin, dmax, dmean = u("3f", 76)
    (ispg,) = u("i", 88)
    (nsymbt,) = u("i", 92)
    exttyp = raw[104:108].decode("ascii", "replace").strip()
    (nversion,) = u("i", 108)
    origin = u("3f", 196)
    magic = raw[208:212].decode("ascii", "replace")
    (rms,) = u("f", 216)
    (nlabl,) = u("i", 220)
    # Words 57-256: ten 80-character labels. Programs stamp their name
    # here and EMDB stamps the accession, so this is the only provenance
    # evidence carried inside the file itself.
    labels = tuple(
        raw[224 + 80 * i : 224 + 80 * (i + 1)].decode("ascii", "replace").rstrip("\x00 ")
        for i in range(max(0, min(nlabl, 10)))
    )

    return MapHeader(
        path=str(path),
        byte_order=order,
        nx=nx,
        ny=ny,
        nz=nz,
        mode=mode,
        nxstart=nxstart,
        nystart=nystart,
        nzstart=nzstart,
        mx=mx,
        my=my,
        mz=mz,
        cella=cella,
        cellb=cellb,
        mapc=mapc,
        mapr=mapr,
        maps=maps,
        dmin=dmin,
        dmax=dmax,
        dmean=dmean,
        ispg=ispg,
        nsymbt=nsymbt,
        exttyp=exttyp,
        nversion=nversion,
        origin=origin,
        magic=magic,
        rms=rms,
        nlabl=nlabl,
        labels=labels,
    )


def _fmt(value: float | None, spec: str = ".6g") -> str:
    return "unknown" if value is None else format(value, spec)


def map_info(path: str | Path) -> str:
    """Human-readable geometry report for an MRC/CCP4 map.

    Leads with voxel size because that is the number this tool exists to
    surface, and ends with warnings because those are what a user needs to act
    on. Does not load the map data and does not touch the network.

    Args:
        path: Path to a ``.mrc``/``.map``/``.ccp4`` file, optionally gzipped.

    Returns:
        A multi-line report.
    """
    h = read_map_header(path)
    vx, vy, vz = h.voxel_size

    iso = h.is_isotropic
    iso_note = {
        True: "isotropic",
        False: "ANISOTROPIC — see warnings",
        None: "indeterminate",
    }[iso]

    lines = [
        f"{Path(h.path).name}",
        "",
        "Voxel size (Å)",
        f"  X {_fmt(vx)}   Y {_fmt(vy)}   Z {_fmt(vz)}      [{iso_note}]",
        "  computed as cella/m — NOT cella/n; see warnings if these differ",
        "",
        "Grid",
        f"  extent   (nx,ny,nz) = ({h.nx}, {h.ny}, {h.nz})   columns/rows/sections",
        f"  sampling (mx,my,mz) = ({h.mx}, {h.my}, {h.mz})   along X/Y/Z",
        f"  start    = ({h.nxstart}, {h.nystart}, {h.nzstart})",
        "",
        "Cell",
        f"  dimensions (Å) = ({h.cella[0]:.6g}, {h.cella[1]:.6g}, {h.cella[2]:.6g})",
        f"  angles     (°) = ({h.cellb[0]:.6g}, {h.cellb[1]:.6g}, {h.cellb[2]:.6g})",
        f"  origin         = ({h.origin[0]:.6g}, {h.origin[1]:.6g}, {h.origin[2]:.6g})",
        "",
        "Layout",
        f"  axis order  = column→{h.axis_mapping[0]}, row→{h.axis_mapping[1]}, "
        f"section→{h.axis_mapping[2]}  (mapc,mapr,maps = {h.mapc},{h.mapr},{h.maps})",
        f"  data mode   = {h.mode} ({_MODES.get(h.mode, 'unrecognised')})",
        f"  byte order  = {h.byte_order}-endian",
        f"  magic       = {h.magic!r}",
        f"  labels      = {'; '.join(h.labels) if h.labels else '(none)'}",
        f"  ext. header = {h.nsymbt} bytes" + (f" ({h.exttyp})" if h.exttyp else ""),
        "",
        "Density statistics (from header, not recomputed)",
        f"  min {h.dmin:.6g}   max {h.dmax:.6g}   mean {h.dmean:.6g}   rms {h.rms:.6g}",
    ]

    warnings = h.warnings()
    lines.append("")
    if warnings:
        lines.append(f"Warnings ({len(warnings)})")
        for w in warnings:
            lines.append(f"  ! {w}")
    else:
        lines.append("Warnings: none")

    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """CLI: ``python -m mcpymol.wiggles.mapinfo MAP [MAP ...]``."""
    import sys

    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print("usage: python -m mcpymol.wiggles.mapinfo MAP [MAP ...]")
        print("\nReport MRC/CCP4 geometry, especially voxel size. Gzip is handled.")
        print("Reads only the 1024-byte header — the map data is never loaded.")
        return 0 if args else 2

    status = 0
    for i, arg in enumerate(args):
        if i:
            print("\n" + "─" * 60 + "\n")
        try:
            print(map_info(arg))
        except (OSError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    raise SystemExit(main())
