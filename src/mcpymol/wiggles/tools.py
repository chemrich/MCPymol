"""The MCP boundary: ten tools, wrapping the port-taking functions.

Every view in this package takes a :class:`~mcpymol.wiggles.port.PymolPort` as
its first argument. That is right for testing — the whole suite runs against
``FakePort`` with no PyMOL installed — and wrong for an MCP tool, because the
model cannot supply a port. So each tool is a thin wrapper: it builds the port,
calls the function, and turns the exceptions into the error strings the rest of
MCPymol returns.

The port-taking function stays the testable core and is where the behaviour
lives. Nothing here should grow logic; if a wrapper starts making decisions,
that decision belongs in the module it wraps, where a test can reach it without
a socket.

**Errors come back as text, not exceptions.** The views raise — ``PortError``
when PyMOL rejects something, ``ValueError`` when an argument is unusable —
because that is the right shape for a library. An MCP client wants a string it
can show, and MCPymol's other tools already answer that way.
"""

from collections.abc import Callable
from typing import Annotated, Any

from pydantic import Field

from mcpymol.app import mcp
from mcpymol.bridge import _SLOW_OP_TIMEOUT, send_request
from mcpymol.wiggles import bfactors as _bfactors
from mcpymol.wiggles import density as _density
from mcpymol.wiggles import ensembles as _ensembles
from mcpymol.wiggles import localres as _localres
from mcpymol.wiggles import mapinfo as _mapinfo
from mcpymol.wiggles import maps as _maps
from mcpymol.wiggles import occupancy as _occupancy
from mcpymol.wiggles import qscore as _qscore
from mcpymol.wiggles.port import PortError, SendRequestPort


def _port() -> SendRequestPort:
    """A port over MCPymol's own bridge.

    ``_SLOW_OP_TIMEOUT`` rather than the default ten seconds: these tools read
    every atom of a structure or push a whole volume across the socket, and a
    map load is measured in tens of megabytes.
    """
    return SendRequestPort(send_request, timeout=_SLOW_OP_TIMEOUT)


def _report(fn: Callable[..., str], *args: Any, **kwargs: Any) -> str:
    """Run a view, returning its report or a readable error.

    ``OSError`` is caught alongside the two the views raise deliberately
    because the map tools open files the model named, and a missing path is an
    ordinary answer here rather than a crash.
    """
    try:
        return fn(*args, **kwargs)
    except (PortError, ValueError, OSError) as exc:
        return f"Error: {exc}"


# ── tier 1: model only ───────────────────────────────────────────────────────


@mcp.tool()
def occupancy_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1ejg")')],
    preserve_bfactors: Annotated[
        bool,
        Field(
            description="Stash the original B-factors before overwriting them so restore_bfactors can put them back (default True)"
        ),
    ] = True,
) -> str:
    """Colour and scale a structure by per-atom crystallographic occupancy (q).

    Alternates at q<1 are visually de-emphasised in proportion, so a
    half-occupied conformer looks half-there instead of looking like a
    confident model.

    This is occupancy in the crystallographic sense — the fraction of unit
    cells in which an atom sits at this position. It is NOT the fraction of
    imaged particles containing a subunit; a model can be q=1.0 everywhere
    while the subunit is present in half the particles. Both are true and they
    are different questions, so this tool never reports the other sense.
    """
    return _report(
        _occupancy.occupancy_view, _port(), obj_name, preserve_bfactors=preserve_bfactors
    )


@mcp.tool()
def altloc_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1ejg")')],
    label: Annotated[
        bool,
        Field(description="Label each alternate location with its occupancy value (default True)"),
    ] = True,
) -> str:
    """Show every altloc group at once, one colour per group, occupancies labelled.

    A multiconformer model says "at this site, these discrete alternatives, in
    these proportions". PyMOL shows one of them by default, which quietly turns
    a statement about heterogeneity into a single structure.
    """
    return _report(_occupancy.altloc_view, _port(), obj_name, label=label)


@mcp.tool()
def ensemble_spread_view(
    obj_name: Annotated[
        str, Field(description="Multi-state PyMOL object name (e.g. an NMR or EM ensemble)")
    ],
    as_putty: Annotated[
        bool,
        Field(
            description="Render as a putty cartoon so spread reads as tube width as well as colour (default True)"
        ),
    ] = True,
) -> str:
    """Colour and thicken a multi-state object by how much its states disagree.

    Per-residue RMS deviation across states, pushed into the B-factor column
    and spectrum-coloured blue (rigid) to red (variable).

    Spread is a description of how much the deposited members differ. It is NOT
    a calibrated uncertainty and NOT an error bar: the number of members and
    the refinement protocol both shape it, so a wider tube means these models
    disagree more, not that the true position is less well determined.
    """
    return _report(_ensembles.ensemble_spread_view, _port(), obj_name, as_putty=as_putty)


@mcp.tool()
def morph_states(
    obj_name: Annotated[str, Field(description="Multi-state PyMOL object name to interpolate")],
    name: Annotated[
        str | None,
        Field(description='Name for the morph object (defaults to "<obj_name>_morph")'),
    ] = None,
    steps: Annotated[
        int, Field(description="Number of interpolated frames to generate (default 30)")
    ] = 30,
    validate_only: Annotated[
        bool,
        Field(description="Run the topology check and report, without creating a morph"),
    ] = False,
) -> str:
    """Interpolate between states, refusing when interpolation is ill-posed.

    The value here is the refusal, not the interpolation — PyMOL morphs
    natively. A morph is only meaningful when every state shares a topology, so
    that atom i in state 1 is the same atom as atom i in state 2. That holds
    for a deformation-model ensemble (3DFlex, DynaMight). It does not hold for
    volumes reconstructed and modelled independently, and a morph across those
    animates a correspondence nobody established.

    Note that cmd.morph is Incentive-only. On open-source PyMOL the topology
    check still runs and reports; only the interpolation is unavailable.
    """
    return _report(
        _ensembles.morph_states,
        _port(),
        obj_name,
        name=name,
        steps=steps,
        validate_only=validate_only,
    )


@mcp.tool()
def qscore_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "9c0k")')],
    validation_path: Annotated[
        str,
        Field(
            description="Path to a wwPDB validation report (mmCIF or XML, optionally gzipped) for this entry"
        ),
    ],
    preserve_bfactors: Annotated[
        bool,
        Field(description="Stash the original B-factors before overwriting them (default True)"),
    ] = True,
) -> str:
    """Colour a model by per-residue Q-score parsed from its wwPDB validation report.

    No network, no map, and no computation — Q-scores are already published,
    and this reads the numbers rather than re-deriving them.

    Two things to expect from real data: entries deposited before the
    September 2023 validation rollout carry no Q-scores at all, which is the
    common case for older EM structures; and real Q-scores go negative, so the
    published 0–1 framing is the intended range and not the observed one.
    """
    return _report(
        _qscore.qscore_view,
        _port(),
        obj_name,
        validation_path,
        preserve_bfactors=preserve_bfactors,
    )


@mcp.tool()
def restore_bfactors(
    obj_name: Annotated[
        str, Field(description="PyMOL object whose original B-factors should be put back")
    ],
) -> str:
    """Restore the B-factors a Wiggles view overwrote.

    occupancy_view, ensemble_spread_view and qscore_view all push their values
    into the B-factor column so PyMOL can spectrum-colour by them. The originals
    are stashed first; this puts them back in one pass.
    """
    return _report(_bfactors.restore_bfactors, _port(), obj_name)


# ── tier 2: maps ─────────────────────────────────────────────────────────────


@mcp.tool()
def map_info(
    path: Annotated[
        str,
        Field(description="Path to an MRC/CCP4 map file (.mrc/.map/.ccp4, optionally gzipped)"),
    ],
) -> str:
    """Report an MRC/CCP4 map's geometry, above all its voxel size.

    Reads only the 1024-byte header — the map data is never loaded, PyMOL is
    never touched, and nothing goes over the network.

    Voxel size is not stored in the header, it is derived, and the derivation
    has a trap: it is cella/m (the grid sampling), not cella/n (the stored
    extent). Those differ on any boxed or cropped map, and dividing by n gives
    a wrong answer of the right order of magnitude. The nominal value is only
    expected to be accurate to ±5–15% in the first place, which at 1.2 Å is a
    systematic stretch of every distance in the model. Anisotropy, cropping and
    axis permutation are all flagged loudly.
    """
    return _report(_mapinfo.map_info, path)


@mcp.tool()
def load_map(
    path: Annotated[
        str,
        Field(description="Path to an MRC/CCP4 volume (.mrc/.map/.ccp4, optionally gzipped)"),
    ],
    name: Annotated[
        str | None,
        Field(description="PyMOL object name (defaults to a sanitised filename stem)"),
    ] = None,
    provenance: Annotated[
        str,
        Field(
            description='How this volume came to exist: "measured", "sharpened", "nn_enhanced", "generated", or "unknown". Defaults to "unknown" and is NEVER inferred.'
        ),
    ] = "unknown",
) -> str:
    """Load a volume into PyMOL, recording where it came from and reporting its geometry.

    The header is parsed before PyMOL is touched, so a malformed file fails
    without leaving a half-loaded object in the session, and the load is
    confirmed rather than assumed.

    Provenance defaults to unknown and is never guessed. A measured
    reconstruction, a sharpened map, a network-enhanced map and a decoder
    output all render identically once they are an isosurface, so defaulting to
    "measured" would quietly assert that somebody observed a generated volume.
    The report shows what the file says about itself — MRC labels, filename
    tokens — so the caller can declare it.
    """
    return _report(_maps.load_map, _port(), path, name, provenance=provenance)


@mcp.tool()
def density_view(
    map_obj: Annotated[str, Field(description="A volume already loaded through the load_map tool")],
    selection: Annotated[
        str, Field(description='PyMOL selection to carve the mesh around (e.g. "chain A")')
    ],
    level: Annotated[
        float | None,
        Field(description="Contour level. If omitted, 1.5 sigma is used and labelled as generic."),
    ] = None,
    units: Annotated[
        str,
        Field(
            description='Units of `level`: "sigma" (what PyMOL contours in) or "absolute" (what EMDB publishes). An absolute level is converted.'
        ),
    ] = "sigma",
    carve: Annotated[
        float, Field(description="Carve radius in Ångström around the selection (default 2.0)")
    ] = 2.0,
    name: Annotated[
        str | None, Field(description='Mesh object name (defaults to "<map_obj>_mesh")')
    ] = None,
) -> str:
    """Draw an isomesh around a selection, stating the contour level in both units.

    PyMOL normalises MRC/CCP4 maps on load, so an isomesh level is in sigma.
    EMDB's author-recommended contour is an absolute map value. They are not
    interchangeable and the difference is not subtle: EMD-30913 publishes 0.05,
    which is 3.16 sigma — used directly as a sigma level it contours noise.

    So the level is reported in sigma AND absolute units every time, an
    absolute level is converted before it reaches PyMOL, and EMDB depositions
    get a pointer to their published contour with the units named. Requires the
    map to have been loaded through load_map, because the conversion needs the
    header's mean and RMS.
    """
    return _report(
        _density.density_view,
        _port(),
        map_obj,
        selection,
        level=level,
        units=units,
        carve=carve,
        name=name,
    )


@mcp.tool()
def local_resolution_view(
    map_obj: Annotated[
        str, Field(description="The density volume to draw, loaded through the load_map tool")
    ],
    res_obj: Annotated[
        str,
        Field(
            description="A local-resolution volume (values in Ångström), loaded through the load_map tool"
        ),
    ],
    level: Annotated[
        float | None,
        Field(description="Contour level for the isosurface. If omitted, 1.5 sigma is used."),
    ] = None,
    units: Annotated[
        str,
        Field(
            description='Units of `level` only: "sigma" or "absolute". Ramp breakpoints are always in Ångström.'
        ),
    ] = "sigma",
    breaks: Annotated[
        list[float] | None,
        Field(
            description="Ascending resolution breakpoints in Ångström. Defaults to the resolution map's own observed range."
        ),
    ] = None,
    palette: Annotated[
        list[str] | None,
        Field(
            description="One PyMOL colour per breakpoint. Defaults to blue (best) through red (worst)."
        ),
    ] = None,
    selection: Annotated[
        str | None,
        Field(description="Restrict the surface to a carve around this PyMOL selection"),
    ] = None,
    carve: Annotated[
        float, Field(description="Carve radius in Ångström; ignored without a selection")
    ] = 2.0,
    name: Annotated[
        str | None, Field(description='Surface object name (defaults to "<map_obj>_localres")')
    ] = None,
    ramp_name: Annotated[
        str | None, Field(description='Ramp object name (defaults to "<res_obj>_ramp")')
    ] = None,
    validate_only: Annotated[
        bool, Field(description="Run the grid check and report, creating nothing")
    ] = False,
) -> str:
    """Colour a map's isosurface by a local-resolution volume instead of by chain.

    A single global resolution number misrepresents almost every map: a rigid
    core at 2.5 Å and a flexible periphery at 5 Å live in the same volume, and
    one isosurface invites the reader to believe all of it is equally real.

    Two things make this harder than it looks, and both are handled. The
    volumes must share a voxel grid — colouring one map by another samples it
    at the first map's coordinates, so a mismatch in extent, spacing, origin or
    axis order draws colour from the wrong place and renders smooth, plausible
    and wrong. This tool refuses and names what differs. And PyMOL normalises
    maps on load, so ramp breakpoints given in Ångström must be converted to
    sigma against the RESOLUTION map's own header — a different sigma scale
    from the contour level, which comes from the density map's. Every
    breakpoint is reported in both units.

    Local resolution is an estimate. Estimators disagree with each other on the
    same map, and the value at a voxel depends on the window and the mask as
    much as on the data.
    """
    return _report(
        _localres.local_resolution_view,
        _port(),
        map_obj,
        res_obj,
        level=level,
        units=units,
        breaks=breaks,
        palette=palette,
        selection=selection,
        carve=carve,
        name=name,
        ramp_name=ramp_name,
        validate_only=validate_only,
    )
