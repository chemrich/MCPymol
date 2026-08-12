"""Colouring density by a local-resolution volume, with the correspondence checked.

A single global FSC number misrepresents almost every map: a rigid core at
2.5 Å and a flexible periphery at 5 Å live in the same volume, and one
isosurface at one contour level invites the user to believe all of it is
equally real. Local resolution estimation produces a per-voxel resolution
field, and nearly every viewer throws it away. See the Wiggles compendium entry
`local-resolution`.

Putting that field on the surface is the same move MCPymol already makes with
``conservation_view`` and ``bfactor_view``: take a per-position scalar nobody
looks at and paint it on. Two things make it harder than it sounds, and this
module exists for both.

**1. The two maps must share a grid.** Colouring surface A by volume B samples
B at A's coordinates. If the grids differ — different extent, different voxel
size, different origin, a permuted axis order — the sample comes from the wrong
place, and the result is not a broken picture but a *plausible* one: smoothly
varying colour that means nothing. This is the same failure as morphing an
ensemble whose states differ in atom count, so it gets the same treatment. The
correspondence is checked, and where it does not hold the tool refuses and says
what differs.

**2. The ramp breakpoints are in sigma, not Ångström.** PyMOL normalises
MRC/CCP4 maps on load (``normalize_ccp4_maps``, on by default), so a
resolution map's stored values are no longer resolutions — 3.2 Å has become
some number of standard deviations about that map's own mean. ``ramp_new``
reads the stored values, so breakpoints given in Å contour the wrong thing
entirely. The conversion is :func:`mcpymol.wiggles.density.to_sigma` against **the
resolution map's** header, not the main map's:

======================  ==========================  ==========================
Breakpoint              Against the main map        Against the resolution map
======================  ==========================  ==========================
3.2 Å                   meaningless — wrong rms     the value ramp_new needs
======================  ==========================  ==========================

Two maps, two sets of header statistics, two sigma scales. Mixing them is the
trap here, exactly as absolute-versus-sigma is the trap in
:mod:`mcpymol.wiggles.density`, and the report states every breakpoint in both units so
the mistake has nowhere to hide.
"""

from __future__ import annotations

from itertools import pairwise

from mcpymol.wiggles.density import (
    DEFAULT_CARVE,
    DEFAULT_SIGMA,
    to_absolute,
    to_sigma,
    usable_rms,
)
from mcpymol.wiggles.mapinfo import MapHeader
from mcpymol.wiggles.maps import LoadedMap, loaded_map, normalisation_state
from mcpymol.wiggles.port import PortError, PymolPort, call
from mcpymol.wiggles.provenance import provenance_banner

#: Relative tolerance for calling two voxel spacings the same, for the same
#: reason :data:`mcpymol.wiggles.mapinfo.ISOTROPY_RTOL` exists: EMD-30913 reports
#: 0.7999967 rather than 0.8, and exact equality would refuse a matched pair.
GRID_RTOL = 1e-4

#: Absolute tolerance (Å) for comparing map origins. Origins are written as
#: float32, so a round trip through a processing package perturbs the last bits.
ORIGIN_ATOL = 1e-3

#: Blue is the *lowest* number of Ångström — the best-resolved density. The
#: scale runs opposite to the usual "red means more" reading, which is why the
#: report says so in words rather than trusting the colours to explain
#: themselves.
DEFAULT_PALETTE = ("blue", "cyan", "green", "yellow", "red")

LOCAL_RESOLUTION_LEGEND = (
    "Local resolution is an ESTIMATE, not a measurement. Estimators disagree with "
    "each other on the same map, and the value at a voxel depends on the window "
    "size and the mask as much as on the data. Read this as a relative picture of "
    "where the reconstruction is better and worse resolved — not as a calibrated "
    "per-voxel resolution, and not as a reason to trust or distrust any individual "
    "atom."
)


def grid_differences(a: MapHeader, b: MapHeader) -> list[str]:
    """Every way two maps fail to share a voxel grid. Empty means they match.

    Each entry is a sentence naming the field and both values, because "the
    grids differ" is not actionable and "the origins differ by 12 Å along Z"
    is.
    """
    out: list[str] = []

    if (a.nx, a.ny, a.nz) != (b.nx, b.ny, b.nz):
        out.append(
            f"extent differs: (nx,ny,nz) = ({a.nx},{a.ny},{a.nz}) vs "
            f"({b.nx},{b.ny},{b.nz}) — the two volumes do not cover the same "
            f"number of voxels"
        )

    if (a.mx, a.my, a.mz) != (b.mx, b.my, b.mz):
        out.append(
            f"grid sampling differs: (mx,my,mz) = ({a.mx},{a.my},{a.mz}) vs "
            f"({b.mx},{b.my},{b.mz}) — one map is boxed or cropped relative to "
            f"the other"
        )

    for axis, va, vb in zip("XYZ", a.voxel_size, b.voxel_size, strict=True):
        if va is None or vb is None:
            out.append(
                f"voxel size along {axis} is undefined in at least one map "
                f"(cell dimension or grid sampling is zero), so the grids "
                f"cannot be compared at all"
            )
        elif abs(va - vb) > GRID_RTOL * max(abs(va), abs(vb)):
            out.append(
                f"voxel size along {axis} differs: {va:.6g} vs {vb:.6g} Å — "
                f"every sampled position drifts, and the drift grows with "
                f"distance from the origin"
            )

    if (a.nxstart, a.nystart, a.nzstart) != (b.nxstart, b.nystart, b.nzstart):
        out.append(
            f"start position differs: ({a.nxstart},{a.nystart},{a.nzstart}) vs "
            f"({b.nxstart},{b.nystart},{b.nzstart}) — the volumes are offset by "
            f"whole voxels"
        )

    for axis, oa, ob in zip("XYZ", a.origin, b.origin, strict=True):
        if abs(oa - ob) > ORIGIN_ATOL:
            out.append(
                f"origin along {axis} differs: {oa:.6g} vs {ob:.6g} Å — the "
                f"volumes are placed differently in space, so the colour at a "
                f"point comes from somewhere else in the resolution map"
            )

    if (a.mapc, a.mapr, a.maps) != (b.mapc, b.mapr, b.maps):
        out.append(
            f"axis order differs: column/row/section map to "
            f"{'/'.join(a.axis_mapping)} vs {'/'.join(b.axis_mapping)} — one "
            f"volume is a transpose of the other"
        )

    return out


def _breakpoints(header: MapHeader, breaks: list[float] | None) -> tuple[list[float], str]:
    """The ramp breakpoints in Å, and a sentence saying where they came from.

    Defaults span the resolution map's own observed range, taken from the
    header's ``dmin``/``dmax`` rather than recomputed — the whole map does not
    need reading to know what it contains.

    Raises:
        ValueError: the breakpoints are unusable, or the header carries no
            positive values for a default to be derived from.
    """
    if breaks is not None:
        if len(breaks) < 2:
            raise ValueError(f"need at least two breakpoints, got {len(breaks)}")
        if any(b <= 0 for b in breaks):
            raise ValueError(
                f"breakpoints are resolutions in Ångström and must be positive, got {breaks}"
            )
        if any(x >= y for x, y in pairwise(breaks)):
            raise ValueError(
                f"breakpoints must ascend — ramp_new requires it, and a descending "
                f"list silently inverts the colour scale. Got {breaks}"
            )
        return list(breaks), "breakpoints as given"

    best, worst = header.dmin, header.dmax
    if worst <= 0:
        raise ValueError(
            f"the resolution map's header reports dmax={worst:.6g}, so it holds no "
            f"positive values. A local-resolution volume stores resolutions in "
            f"Ångström; this one does not look like a resolution field at all."
        )
    if best <= 0:
        # Estimators write 0 outside the mask. Zero is not a resolution, and
        # taking it as the best-resolved end would compress the whole ramp.
        best = worst / 10 if worst > 0 else worst
        note = (
            f"breakpoints span {best:.3g}–{worst:.3g} Å; the header's dmin is "
            f"{header.dmin:.3g}, which is outside-the-mask padding rather than a "
            f"resolution, so the low end is estimated instead"
        )
    else:
        note = f"breakpoints span the header's own range, {best:.3g}–{worst:.3g} Å"

    n = len(DEFAULT_PALETTE)
    step = (worst - best) / (n - 1)
    return [best + step * i for i in range(n)], note


def _require_loaded(port: PymolPort, obj: str, role: str) -> LoadedMap:
    record = loaded_map(obj, port)
    if record is None:
        raise PortError(
            f"{obj!r} (the {role}) was not loaded through load_map, so its header "
            f"is unavailable. Without both headers this tool cannot check that "
            f"the two volumes share a grid, and colouring one map by another that "
            f"does not match it produces a picture that looks fine and means "
            f"nothing. Load it with load_map first."
        )
    return record


def local_resolution_view(
    port: PymolPort,
    map_obj: str,
    res_obj: str,
    *,
    level: float | None = None,
    units: str = "sigma",
    breaks: list[float] | None = None,
    palette: list[str] | tuple[str, ...] | None = None,
    selection: str | None = None,
    carve: float = DEFAULT_CARVE,
    name: str | None = None,
    ramp_name: str | None = None,
    validate_only: bool = False,
) -> str:
    """Draw ``map_obj``'s isosurface coloured by the resolution field in ``res_obj``.

    Both volumes must have been loaded through :func:`mcpymol.wiggles.maps.load_map`:
    the headers are what make the grid check and the Å-to-sigma conversion
    possible, and without them this tool would be guessing at both.

    Args:
        port: A live or fake PyMOL port.
        map_obj: The volume to draw. Its contour level is in *its* sigma.
        res_obj: A local-resolution volume, values in Ångström.
        level: Contour level for the isosurface. ``None`` uses
            :data:`mcpymol.wiggles.density.DEFAULT_SIGMA` and says so.
        units: ``"sigma"`` or ``"absolute"``, for ``level`` only. Breakpoints
            are always in Ångström.
        breaks: Ascending resolution breakpoints in Å. Defaults to the
            resolution map's own observed range.
        palette: One colour per breakpoint. Defaults to
            :data:`DEFAULT_PALETTE`, blue (best) through red (worst).
        selection: Restrict the surface to a carve around this selection.
        carve: Carve radius in Å. Ignored without ``selection``.
        name: Surface object name. Defaults to ``<map_obj>_localres``.
        ramp_name: Ramp object name. Defaults to ``<res_obj>_ramp``. The ramp
            must outlive the call — deleting it un-colours the surface.
        validate_only: Run the grid check and report, creating nothing.

    Returns:
        A report: the grid check, every breakpoint in both Å and sigma, both
        provenance banners, and the legend. Refusals name what differs.

    Raises:
        PortError: either object was not loaded through ``load_map``.
        ValueError: bad ``units``, bad breakpoints, or a palette that does not
            match the breakpoints.
    """
    if units not in ("sigma", "absolute"):
        raise ValueError(f"units must be 'sigma' or 'absolute', got {units!r}")

    if map_obj == res_obj:
        return "\n".join(
            [
                f"local_resolution_view({map_obj} by {res_obj})",
                "",
                "  REFUSED: the same object was given as both the density and the",
                "  resolution field. Colouring a map by itself paints contour level,",
                "  not resolution, and the picture is indistinguishable from a real",
                "  one. Load the local-resolution volume as a separate object.",
            ]
        )

    main = _require_loaded(port, map_obj, "density map")
    res = _require_loaded(port, res_obj, "resolution map")

    differences = grid_differences(main.header, res.header)
    if differences:
        return "\n".join(
            [
                f"local_resolution_view({map_obj} by {res_obj})",
                "",
                "  REFUSED: the two volumes do not share a voxel grid. A point on",
                f"  {map_obj}'s surface does not correspond to the same point in {res_obj},",
                "  so every colour would be sampled from the wrong place.",
                "",
                *(f"    - {d}" for d in differences),
                "",
                "  This is not a rendering that comes out visibly broken. It comes",
                "  out smooth, plausible, and wrong.",
                "",
                "  A local-resolution map is normally written on the grid of the map",
                "  it was computed from. If these came from different reconstructions",
                "  or one has been resampled, regenerate the resolution volume on this",
                "  map's grid. `map_info` on both files shows the full geometry.",
            ]
        )

    try:
        points, breaks_note = _breakpoints(res.header, breaks)
    except ValueError as exc:
        if breaks is not None:
            raise
        return "\n".join(
            [
                f"local_resolution_view({map_obj} by {res_obj})",
                "",
                f"  REFUSED: {exc}",
                "",
                "  Pass breaks=[...] explicitly if the header statistics are stale but",
                "  the data are good.",
            ]
        )

    colours = list(palette) if palette is not None else list(DEFAULT_PALETTE)
    if len(colours) != len(points):
        raise ValueError(
            f"{len(points)} breakpoints but {len(colours)} colours — ramp_new pairs "
            f"them one to one, so they must be the same length"
        )

    normalised = normalisation_state(port)
    if normalised is False:
        # Values are stored as written, so the breakpoints are already right.
        sigmas = list(points)
    else:
        try:
            sigmas = [to_sigma(res.header, p) for p in points]
        except ValueError as exc:
            return "\n".join(
                [
                    f"local_resolution_view({map_obj} by {res_obj})",
                    "",
                    f"  REFUSED: {exc}",
                    "",
                    "  PyMOL normalised this volume on load, so its stored values are in",
                    "  sigma and the breakpoints have to be converted — which needs a",
                    "  usable rms. Reload with `set normalize_ccp4_maps, off` to keep the",
                    "  Ångström values, or repair the header statistics.",
                ]
            )

    # The main map's contour level, in the main map's sigma — a different scale
    # from the breakpoints above, which are in the resolution map's sigma.
    used_default = level is None
    if used_default:
        contour = DEFAULT_SIGMA
    elif units == "absolute":
        contour = to_sigma(main.header, float(level))  # type: ignore[arg-type]
    else:
        contour = float(level)  # type: ignore[arg-type]
    contour_absolute = to_absolute(main.header, contour) if usable_rms(main.header) else None

    surface = name or f"{map_obj}_localres"
    ramp = ramp_name or f"{res_obj}_ramp"

    lines = [
        f"local_resolution_view({map_obj} by {res_obj})",
        "",
        f"  Grid check passed: {main.header.nx}x{main.header.ny}x{main.header.nz} voxels, "
        f"{_voxel_text(main.header)}, same origin and axis order.",
        f"  {map_obj} and {res_obj} sample the same points in space.",
        "",
    ]

    if validate_only:
        lines += ["  validate_only=True — nothing created."]
        lines += ["", _ramp_table(points, sigmas, colours, normalised)]
        lines += ["", provenance_banner(map_obj), "", provenance_banner(res_obj)]
        lines += ["", LOCAL_RESOLUTION_LEGEND]
        return "\n".join(lines)

    if selection:
        call(port, "isosurface", surface, map_obj, contour, selection, carve=carve)
    else:
        call(port, "isosurface", surface, map_obj, contour)
    call(port, "ramp_new", ramp, res_obj, sigmas, colours)
    call(port, "set", "surface_color", ramp, surface)

    absolute_text = f"{contour_absolute:.6g}" if contour_absolute is not None else "unknown (rms=0)"
    lines += [
        f"  Surface `{surface}` at {contour:.3g} sigma  =  {absolute_text} absolute"
        + (f", carved {carve:g} Å around {selection}." if selection else "."),
        f"  Coloured by ramp `{ramp}` over {res_obj}. Deleting the ramp un-colours",
        "  the surface — it is not baked in.",
        "",
        _ramp_table(points, sigmas, colours, normalised),
        "",
        f"  {breaks_note}.",
    ]

    if used_default:
        lines += [
            f"  No contour level given, so {DEFAULT_SIGMA} sigma was used — a generic",
            "  starting point, not a recommendation for this map. density_view reports",
            "  the author-recommended contour for EMDB depositions.",
        ]

    lines += [
        "",
        "  Two sigma scales are in play and they are not interchangeable: the",
        f"  contour above is in {map_obj}'s sigma, the breakpoints in {res_obj}'s.",
        "",
        provenance_banner(map_obj),
        "",
        provenance_banner(res_obj),
    ]

    for label, header in ((map_obj, main.header), (res_obj, res.header)):
        warnings = header.warnings()
        if warnings:
            lines += ["", f"  Geometry warnings, {label} ({len(warnings)})"]
            lines += [f"    ! {w}" for w in warnings]

    lines += ["", LOCAL_RESOLUTION_LEGEND]
    return "\n".join(lines)


def _voxel_text(header: MapHeader) -> str:
    vx, vy, vz = header.voxel_size
    parts = ["unknown" if v is None else f"{v:.4g}" for v in (vx, vy, vz)]
    return f"voxel {'/'.join(parts)} Å"


def _ramp_table(
    points: list[float],
    sigmas: list[float],
    colours: list[str],
    normalised: bool | None,
) -> str:
    """Every breakpoint in both units, plus what the colours mean."""
    rows = [
        f"    {p:>6.2f} Å   ->  {s:>8.3g}   {c}"
        for p, s, c in zip(points, sigmas, colours, strict=True)
    ]
    unit_header = "as sent" if normalised is False else "sigma"
    head = [
        "  Resolution ramp",
        f"    {'Å':>6}         {unit_header:>8}   colour",
        *rows,
        "",
        f"    {colours[0]} is {points[0]:.3g} Å — the BEST-resolved density. "
        f"{colours[-1]} is {points[-1]:.3g} Å,",
        "    the worst. Low numbers are good, so the scale runs opposite to the",
        "    usual reading of a colour ramp.",
    ]

    if normalised is False:
        head += [
            "",
            "    normalize_ccp4_maps is OFF, so the volume keeps its Ångström values",
            "    and the breakpoints were sent unconverted.",
        ]
    elif normalised is True:
        head += [
            "",
            "    normalize_ccp4_maps is ON, so PyMOL rescaled the volume on load and",
            "    the Ångström breakpoints were converted to sigma against the",
            "    resolution map's own mean and rms. Sending Å directly would ramp",
            "    over a range the data never reaches.",
        ]
    else:
        head += [
            "",
            "    PyMOL would not report normalize_ccp4_maps, so its default (ON) was",
            "    assumed and the breakpoints were converted to sigma. If the map was",
            "    loaded with normalisation off, the ramp is wrong — check the setting.",
        ]

    head += [
        "",
        "    The setting is read now, not as it was at load time. A session that",
        "    changed it in between would be reported wrongly.",
    ]
    return "\n".join(head)
