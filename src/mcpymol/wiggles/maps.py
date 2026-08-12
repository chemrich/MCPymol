"""Loading volumes into PyMOL, with the geometry and the provenance stated.

``map_info`` in :mod:`mcpymol.wiggles.mapinfo` reads a header and touches nothing.
This module is the PyMOL-facing half: it loads the volume, records where it
came from (SPEC invariant **I1**), and reports the voxel size at load time
rather than leaving it to be discovered later — see
the Wiggles compendium entry `voxel-size` for why a silently wrong voxel size is
a systematic stretch of every distance in the model.

Loader discipline, from MCPymol issue #15: verify the object actually arrived
before doing anything else, scope any cleanup to the object created here, and
return an error on an empty load rather than a success message.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mcpymol.wiggles.mapinfo import MapHeader, read_map_header
from mcpymol.wiggles.port import PortError, PymolPort, call
from mcpymol.wiggles.provenance import (
    Evidence,
    Provenance,
    declare,
    gather_evidence,
    provenance_banner,
)


@dataclass(frozen=True)
class LoadedMap:
    """What we know about a volume that is in the session.

    Kept because a contour level cannot be reported honestly without it: PyMOL
    contours in sigma, EMDB publishes absolute levels, and converting between
    them needs the header's ``dmean`` and ``rms``.
    """

    obj: str
    path: Path
    header: MapHeader
    evidence: Evidence


# Process-global. Correct for a single-session MCP server, which is what
# MCPymol is today; wrong the moment one process serves two sessions. Keyed
# by object name, so two sessions with the same object name would collide.
# Flagged rather than solved — see MOVING.md.
_LOADED: dict[str, LoadedMap] = {}


def loaded_map(obj: str) -> LoadedMap | None:
    """The record for ``obj``, or None if it was not loaded through load_map."""
    return _LOADED.get(obj)


def forget_map(obj: str | None = None) -> None:
    """Drop the record for ``obj``, or all of them."""
    if obj is None:
        _LOADED.clear()
    else:
        _LOADED.pop(obj, None)


def _default_name(path: Path) -> str:
    """A PyMOL-safe object name from a filename."""
    stem = path.name
    for suffix in (".gz", ".mrc", ".map", ".ccp4"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
    return "".join(c if (c.isalnum() or c == "_") else "_" for c in stem) or "map"


def load_map(
    port: PymolPort,
    path: str | Path,
    name: str | None = None,
    *,
    provenance: Provenance | str = Provenance.UNKNOWN,
) -> str:
    """Load an MRC/CCP4 volume, recording its provenance and reporting geometry.

    The header is parsed **before** PyMOL is touched, so a malformed file fails
    without leaving a half-loaded object in the session.

    Args:
        port: A live or fake PyMOL port.
        path: Path to a ``.mrc`` / ``.map`` / ``.ccp4``, optionally gzipped.
        name: PyMOL object name. Defaults to a sanitised filename stem.
        provenance: How this volume came to exist. **Defaults to UNKNOWN and is
            never inferred** — adopting a guess would risk asserting that a
            generated volume was measured. The report shows what the file says
            about itself so you can declare it.

    Returns:
        A report: voxel size, geometry warnings, and the provenance banner.

    Raises:
        PortError: the object did not arrive in the session.
        ValueError / OSError: the file is not a readable MRC/CCP4 map.
    """
    path = Path(path)
    header = read_map_header(path)  # fails here, before PyMOL sees anything
    obj = name or _default_name(path)

    if isinstance(provenance, str):
        provenance = Provenance(provenance)

    call(port, "load", str(path), obj)

    # Issue #15 discipline: confirm it arrived rather than assuming.
    names = port.query("get_names", "objects")
    if not isinstance(names, (list, tuple)) or obj not in names:
        # Scope any cleanup to the object we tried to create — never the session.
        call(port, "delete", obj)
        raise PortError(
            f"loading {path} produced no object named {obj!r} "
            f"(session has: {sorted(names) if isinstance(names, (list, tuple)) else names!r}). "
            f"Nothing else in the session was touched."
        )

    declare(obj, provenance)
    evidence = gather_evidence(header, path)
    _LOADED[obj] = LoadedMap(obj=obj, path=path, header=header, evidence=evidence)

    vx, vy, vz = header.voxel_size
    iso = header.is_isotropic
    iso_note = {True: "isotropic", False: "ANISOTROPIC", None: "indeterminate"}[iso]

    def fmt(v: float | None) -> str:
        return "unknown" if v is None else f"{v:.6g}"

    lines = [
        f"load_map({path.name} -> {obj})",
        "",
        f"  Voxel size (Å): X {fmt(vx)}  Y {fmt(vy)}  Z {fmt(vz)}   [{iso_note}]",
        f"  Grid: {header.nx}x{header.ny}x{header.nz}   axis order {'/'.join(header.axis_mapping)}",
        "",
        provenance_banner(obj),
    ]

    if provenance is Provenance.UNKNOWN:
        lines += ["", "  What the file says about itself:"]
        lines += [f"    - {reason}" for reason in evidence.reasons]
        if evidence.suggested is not Provenance.UNKNOWN:
            lines.append(
                f"    Suggested: {evidence.suggested.value}. Not adopted — pass "
                f"provenance={evidence.suggested.value!r} to declare it."
            )

    warnings = header.warnings()
    if warnings:
        lines += ["", f"  Geometry warnings ({len(warnings)})"]
        lines += [f"    ! {w}" for w in warnings]

    return "\n".join(lines)


def normalisation_state(port: PymolPort) -> bool | None:
    """Is ``normalize_ccp4_maps`` on? None when PyMOL will not say.

    Read now, not at load time, which is the honest limitation: a session that
    had it off when the map was loaded and on now would report the wrong thing.

    **What a real PyMOL returns, checked rather than assumed:** ``'1'`` — the
    *string* ``'1'``, not ``'on'`` and not the integer ``1`` (open-source PyMOL,
    2026-08-09, via ``tools/livefire.py``). The parse stays tolerant of the
    other spellings because that answer is one build's, and an older plugin may
    not expose ``get`` at all, which is unknown rather than an error.

    Public because ``tools/livefire.py`` prints the raw answer next to this
    verdict: what a real PyMOL returns here is the one thing in this module
    that ``FakePort`` cannot establish.
    """
    try:
        raw = port.query("get", "normalize_ccp4_maps")
    except PortError:
        return None
    text = str(raw).strip().lower()
    if text in ("on", "1", "1.0", "true", "yes"):
        return True
    if text in ("off", "0", "0.0", "false", "no"):
        return False
    return None
