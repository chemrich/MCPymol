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
from mcpymol.wiggles import composition as _composition
from mcpymol.wiggles import deformation as _deformation
from mcpymol.wiggles import density as _density
from mcpymol.wiggles import ensembles as _ensembles
from mcpymol.wiggles import heterogeneity as _heterogeneity
from mcpymol.wiggles import latent as _latent
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


# ── tier 3: ensembles ────────────────────────────────────────────────────────


@mcp.tool()
def load_ensemble(
    directory: Annotated[
        str,
        Field(
            description="Path to a heterogeneity job directory (cryoDRGN, 3DVA, RECOVAR, 3DFlex, DynaMight)"
        ),
    ],
    name: Annotated[
        str | None,
        Field(description="Prefix for the loaded frame objects (defaults to the directory name)"),
    ] = None,
    method: Annotated[
        str | None,
        Field(
            description='Declare the generating method instead of detecting it: "cryodrgn", "3dva", "recovar", "3dflex", "dynamight", "cryospire"'
        ),
    ] = None,
    provenance: Annotated[
        str,
        Field(
            description='How the volumes came to exist. Defaults to "generated" because heterogeneity output is decoder output or subspace interpolation.'
        ),
    ] = "generated",
    max_volumes: Annotated[
        int, Field(description="Ceiling on frames loaded; any excess is reported, not hidden")
    ] = 50,
    trust_pickle: Annotated[
        bool,
        Field(
            description="Permit reading a .pkl latent table. Off by default because unpickling runs arbitrary code from the file."
        ),
    ] = False,
) -> str:
    """Load a cryo-EM heterogeneity job's volumes in trajectory order.

    Every method in this space writes the same thing: a directory of maps plus a
    table of latent coordinates. This reads them in trajectory order — naturally,
    so frame 10 does not sort between 1 and 2 and silently reorder the motion.

    The generating method is identified from documented on-disk markers and is
    never guessed. A directory that matches nothing loads its volumes fine and
    stays unidentified, which makes latent views refuse to render: an unlabelled
    latent plot is what SPEC invariant I2 forbids, because each method's caveat
    is different and they are not interchangeable. Pass method= if you know.

    Provenance defaults to "generated" rather than "unknown" — nothing observed a
    decoder's output — and the report says it declared that and on what evidence.
    """
    return _report(
        _heterogeneity.load_ensemble,
        _port(),
        directory,
        name,
        method=method,
        provenance=provenance,
        max_volumes=max_volumes,
        trust_pickle=trust_pickle,
    )


@mcp.tool()
def latent_traverse_view(
    ensemble_name: Annotated[
        str, Field(description="An ensemble previously loaded through the load_ensemble tool")
    ],
    level: Annotated[
        float | None,
        Field(description="Contour level. If omitted, 1.5 sigma against the first frame is used."),
    ] = None,
    units: Annotated[
        str,
        Field(description='Units of `level`: "sigma" (read against frame 1) or "absolute"'),
    ] = "sigma",
    name: Annotated[
        str | None, Field(description='Prefix for the isosurfaces (defaults to "<ensemble>_surf")')
    ] = None,
    color: Annotated[
        str,
        Field(
            description="One colour for every frame; a per-frame spectrum would encode frame index as if measured"
        ),
    ] = "skyblue",
    build_movie: Annotated[
        bool,
        Field(
            description="Wire the frames to PyMOL's movie timeline so they can be stepped and played"
        ),
    ] = True,
) -> str:
    """Step through an ensemble's conformations as a contoured trajectory.

    Refuses when the generating method was never identified. That is SPEC
    invariant I2 working, not a failure: cryoDRGN's latent density can bear no
    relation to the truth while 3DVA's frames are linear interpolations no
    particle occupied, and rendering both the same way asserts something for one
    that holds only for the other.

    The contour is held at a constant ABSOLUTE value and converted to each
    frame's own sigma. PyMOL normalises every map independently on load, so a
    fixed sigma level would contour each frame at its own scale and flatten away
    the density change the traversal exists to show.

    No latent scatter and no density estimate are drawn, deliberately. In a
    41-submission blind challenge most methods missed a genuinely present middle
    state, so an empty region of latent space is a region of unknown occupancy —
    never evidence that a conformation does not occur.
    """
    return _report(
        _latent.latent_traverse_view,
        _port(),
        ensemble_name,
        level=level,
        units=units,
        name=name,
        color=color,
        build_movie=build_movie,
    )


@mcp.tool()
def deformation_view(
    obj_name: Annotated[
        str, Field(description="A multi-state object whose states are conformations")
    ],
    start_state: Annotated[int, Field(description="The state to measure displacement from")] = 1,
    end_state: Annotated[
        int | None, Field(description="The state to measure displacement to (defaults to the last)")
    ] = None,
    arrows: Annotated[
        bool, Field(description="Draw CGO arrows from each residue's start to its end position")
    ] = True,
    arrow_scale: Annotated[
        float,
        Field(
            description="Multiply arrow length; above 1.0 this is an exaggeration and is labelled as one"
        ),
    ] = 1.0,
    max_arrows: Annotated[
        int, Field(description="Ceiling on arrows drawn; anything dropped is reported")
    ] = 60,
    as_putty: Annotated[
        bool, Field(description="Also scale cartoon tube width by displacement")
    ] = False,
    uncertainty_path: Annotated[
        str | None,
        Field(
            description="Path to a 'chain resi value' table of half-set uncertainty; colours the arrows when supplied"
        ),
    ] = None,
    preserve_bfactors: Annotated[
        bool,
        Field(description="Stash the original B-factors before overwriting them (default True)"),
    ] = True,
) -> str:
    """Colour and arrow a model by how far each residue moved between two states.

    This is the well-supported tool in the tier. A 41-submission blind challenge
    found the molecular motions these methods recover resemble both each other
    and the ground truth; it was relative populations that fell apart. So the
    arrows mean: density moved this way. They do not mean: this fraction of
    particles moved this way, and no population can be read out of them.

    Refuses when the two states differ in atom count — atoms cannot be paired
    across them, so a displacement field would be meaningless rather than
    approximate.

    Where a method estimates deformation on independent half-sets, pass that
    table and the arrows are coloured by the disagreement. Without it every arrow
    is drawn with identical confidence and they do not have identical confidence,
    which the report says rather than leaving the picture to imply otherwise.
    """
    return _report(
        _deformation.deformation_view,
        _port(),
        obj_name,
        start_state=start_state,
        end_state=end_state,
        arrows=arrows,
        arrow_scale=arrow_scale,
        max_arrows=max_arrows,
        as_putty=as_putty,
        uncertainty=uncertainty_path,
        preserve_bfactors=preserve_bfactors,
    )


@mcp.tool()
def composition_view(
    obj_name: Annotated[str, Field(description="The object the table's selections apply to")],
    table: Annotated[
        str,
        Field(
            description='Presence fractions, inline as "chain A=0.4, chain B=1.0" or a path to a file of selection<TAB>fraction lines'
        ),
    ],
    transparency: Annotated[
        bool,
        Field(
            description="Also make rarely-present parts transparent in proportion, so half-present reads as half there"
        ),
    ] = True,
    label: Annotated[
        bool, Field(description="Write each part's presence fraction next to it")
    ] = True,
) -> str:
    """Colour parts of a structure by the fraction of particles containing them.

    This is occupancy in SENSE 2 — compositional. A ligand at 40% here means 40%
    of imaged complexes had it bound. It is not the per-atom crystallographic
    occupancy q, and it is never derived from q: a model can be q=1.0 at every
    atom while the subunit it belongs to is present in half the particles. Both
    statements are true and they answer different questions, so the table must be
    supplied and this tool never reads an atom property at all. For sense 1, use
    occupancy_view.

    A selection matching no atoms is refused rather than skipped — PyMOL would
    accept it silently, colour nothing, and leave a render that looks like a
    fully-present structure.
    """
    return _report(
        _composition.composition_view,
        _port(),
        obj_name,
        table,
        transparency=transparency,
        label=label,
    )
