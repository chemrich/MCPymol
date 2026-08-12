"""Showing density around a selection, with the contour level stated honestly.

Every viewer makes users hunt for a contour threshold. EMDB has been publishing
the depositor's own recommended level all along — ``map.contour_list[].level``
with ``source: AUTHOR`` — and nothing surfaces it. See
``REPOSITORIES.md``.

**The unit trap this module exists to close.** PyMOL normalises MRC/CCP4 maps
on load (``normalize_ccp4_maps``, on by default), so an ``isomesh`` level is in
**sigma**. EMDB's author contour is an **absolute** map value. They are not
interchangeable, and confusing them is not a subtle error:

======================  =================  ==================================
Entry                   Author contour     As a PyMOL level if used directly
======================  =================  ==================================
EMD-30913               0.05 absolute      0.05 sigma — contours the noise
EMD-11638               0.116 absolute     0.116 sigma — contours the noise
======================  =================  ==================================

Converted properly those are **3.16 sigma** and **4.42 sigma**. The conversion
needs the map's mean and RMS, which the MRC header carries, which is why
:func:`density_view` requires the map to have been loaded through
:func:`mcpymol.wiggles.maps.load_map`: without the header it cannot state a level in
both units, and stating one without saying which is how the mistake happens.
"""

from __future__ import annotations

from mcpymol.wiggles.mapinfo import MapHeader
from mcpymol.wiggles.maps import loaded_map, normalisation_state
from mcpymol.wiggles.port import PortError, PymolPort, call
from mcpymol.wiggles.provenance import provenance_banner

#: Used when no level is given and no author contour is available. A generic
#: starting point, not a recommendation for any particular map — the report
#: says so rather than letting it read as considered.
DEFAULT_SIGMA = 1.5

#: Default carve radius (Å) around the selection.
DEFAULT_CARVE = 2.0


def to_sigma(header: MapHeader, absolute: float) -> float:
    """Convert an absolute map value to sigma units.

    Raises:
        ValueError: the header reports zero RMS, so sigma is undefined and no
            conversion is possible.
    """
    if not header.rms:
        raise ValueError(
            "header reports rms=0, so sigma is undefined and an absolute level "
            "cannot be converted. The map may be unnormalised or its header "
            "statistics stale."
        )
    return (absolute - header.dmean) / header.rms


def to_absolute(header: MapHeader, sigma: float) -> float:
    """Convert a sigma level back to an absolute map value."""
    return header.dmean + sigma * header.rms


def _normalisation_note(normalised: bool | None, level_sent: float) -> str:
    """One line saying which number reached ``isomesh``, and why.

    The sigma/absolute pair above it is a statement about the map; this is a
    statement about the command that was actually run. They differ exactly
    when PyMOL is not normalising, which is the case the user most needs
    told — silently swapping units is how a correct-looking report ends up
    describing a mesh that was drawn somewhere else.
    """
    if normalised is False:
        return (
            f"  normalize_ccp4_maps is OFF, so isomesh was given {level_sent:.6g} "
            f"in the map's own units."
        )
    if normalised is None:
        return (
            f"  PyMOL would not report normalize_ccp4_maps; assumed on, so isomesh "
            f"was given {level_sent:.3g} sigma. If it is off, this contour is wrong."
        )
    return f"  normalize_ccp4_maps is on, so isomesh was given {level_sent:.3g} sigma."


def density_view(
    port: PymolPort,
    map_obj: str,
    selection: str,
    *,
    level: float | None = None,
    units: str = "sigma",
    carve: float = DEFAULT_CARVE,
    name: str | None = None,
) -> str:
    """Draw an isomesh around ``selection``, reporting the level in both units.

    Args:
        port: A live or fake PyMOL port.
        map_obj: A volume loaded through :func:`mcpymol.wiggles.maps.load_map`.
        selection: What to carve the mesh around.
        level: Contour level. ``None`` uses :data:`DEFAULT_SIGMA` and says so.
        units: ``"sigma"`` (PyMOL's units) or ``"absolute"`` (EMDB's). An
            absolute level is converted, and both values are reported.
        carve: Carve radius in Å around the selection.
        name: Mesh object name. Defaults to ``<map_obj>_mesh``.

    Returns:
        A report: the level in sigma *and* absolute units, the sigma value
        itself, the provenance banner, and — where the map is an EMDB
        deposition and no level was given — a pointer to the author contour.

    Raises:
        PortError: ``map_obj`` was not loaded through ``load_map``, so its
            header is unavailable and the level cannot be stated in both units.
        ValueError: ``units`` is not recognised, or conversion is impossible.
    """
    if units not in ("sigma", "absolute"):
        raise ValueError(f"units must be 'sigma' or 'absolute', got {units!r}")

    record = loaded_map(map_obj)
    if record is None:
        raise PortError(
            f"{map_obj!r} was not loaded through load_map, so its header is "
            f"unavailable. Without it a contour level cannot be stated in both "
            f"sigma and absolute units, and stating one without saying which is "
            f"exactly how an EMDB author contour gets used as a sigma level. "
            f"Load it with load_map first."
        )

    header = record.header
    used_default = level is None

    if used_default:
        sigma = DEFAULT_SIGMA
        absolute: float | None = None
    elif units == "absolute":
        sigma = to_sigma(header, float(level))  # type: ignore[arg-type]
        absolute = float(level)  # type: ignore[arg-type]
    else:
        sigma = float(level)  # type: ignore[arg-type]
        absolute = None

    if absolute is None and header.rms:
        absolute = to_absolute(header, sigma)

    # Which number isomesh wants depends on whether PyMOL normalised the volume
    # on load. With normalize_ccp4_maps on (the default) a level is in sigma;
    # with it off the map keeps its stored values and the level is read as a
    # raw density. Sending sigma to an un-normalised map is the failure this
    # module exists to prevent: EMD-30913's published 0.05 is 3.16 sigma, and
    # 3.16 as a raw value contours nothing at all — an empty mesh under a
    # report claiming the depositor's own level was applied.
    normalised = normalisation_state(port)
    if normalised is False:
        if absolute is None:
            raise PortError(
                f"{map_obj!r} was loaded with normalize_ccp4_maps off, so a contour "
                f"level is read in the map's own units — but its header reports "
                f"rms={header.rms:g}, so the sigma level given cannot be converted "
                f"to one. Give the level in absolute units, or repair the header "
                f"statistics. Contouring on the unconverted number would draw a "
                f"surface at the wrong density and say nothing was wrong."
            )
        level_sent = absolute
    else:
        level_sent = sigma

    mesh = name or f"{map_obj}_mesh"
    call(port, "isomesh", mesh, map_obj, level_sent, selection, carve=carve)

    absolute_text = f"{absolute:.6g}" if absolute is not None else "unknown (rms=0)"
    lines = [
        f"density_view({map_obj} around {selection})",
        "",
        f"  Contour: {sigma:.3g} sigma  =  {absolute_text} absolute",
        f"  Map sigma (header rms): {header.rms:.6g}   mean: {header.dmean:.6g}",
        f"  Mesh `{mesh}`, carved {carve:g} Å around the selection.",
        _normalisation_note(normalised, level_sent),
        "",
    ]

    if used_default:
        lines += [
            f"  No level given, so {DEFAULT_SIGMA} sigma was used. That is a generic",
            "  starting point, not a recommendation for this map.",
        ]
        accession = record.evidence.emdb_accession
        if accession:
            lines += [
                "",
                f"  {accession} has an author-recommended contour published by the",
                "  depositor. It is an ABSOLUTE value, so pass it as:",
                f"      density_view(port, {map_obj!r}, {selection!r}, "
                f"level=<value>, units='absolute')",
                "  Retrieve it from:",
                f"      https://www.ebi.ac.uk/emdb/api/entry/{accession}"
                "  ->  map.contour_list.contour[].level",
            ]
        lines.append("")
    elif units == "absolute":
        lines += [
            "  Level was given in absolute units and converted to sigma, which is",
            "  what PyMOL contours in. Passing an absolute EMDB contour straight",
            "  to PyMOL would contour near zero and show mostly noise.",
            "",
        ]

    lines.append(provenance_banner(map_obj))

    warnings = header.warnings()
    if warnings:
        lines += ["", f"  Geometry warnings ({len(warnings)})"]
        lines += [f"    ! {w}" for w in warnings]

    return "\n".join(lines)
