"""Stepping through a latent trajectory, without claiming what it cannot.

This is the tool the compendium is most careful about, and the care is not
decorative. The Flatiron Institute blind challenge put 41 submissions on the
same thyroglobulin data and found the field recovers **motion** well and
**populations** badly: the simulated ground truth had three modes, most
submissions recovered two, and only three of the 41 found all three. The mode
they missed was the middle one. `limits` gives the information-theoretic reason
— measurement noise sets a floor on how finely conformation space can be
resolved, and a state sandwiched between two others is what that floor eats
first.

Two SPEC invariants fall out of that, and both are enforced here.

**I2 — no unlabelled latent plot.** Every rendering names its generating method
and carries that method's caveat. The caveats are not interchangeable: cryoDRGN
latent density can bear no relation to the truth, while 3DVA frames are linear
interpolations no particle ever occupied. When
:func:`~mcpymol.wiggles.heterogeneity.load_ensemble` could not identify the
method, this **refuses to render** rather than draw an unlabelled trajectory.

**I3 — a gap is not an absence.** Nothing here may state or imply that an
unpopulated region of latent space is unvisited. :data:`ABSENCE_CLAIMS` lists
the phrasings that would breach that, and a test asserts no output contains
one. There is deliberately **no scatter plot and no density estimate** in this
tool: the single most requested latent rendering is the one the evidence says
is most often wrong, and the honest version of it is not a better-drawn plot,
it is not drawing it.

**The trap this tool exists to close.** PyMOL normalises each map independently
on load, so contouring every frame at "1.5 sigma" contours each frame at *its
own* sigma. A traversal is precisely the case where the density genuinely
changes between frames — that is the signal — and per-frame normalisation
rescales it away. Frames drawn that way look reassuringly similar no matter
what the ensemble does. So the level is fixed in **absolute** terms and
converted to each frame's own sigma, and the spread of those sigmas is reported,
because that spread *is* the density change.
"""

from __future__ import annotations

from mcpymol.wiggles.density import DEFAULT_SIGMA, to_absolute, to_sigma, usable_rms
from mcpymol.wiggles.heterogeneity import Ensemble, Method, loaded_ensemble
from mcpymol.wiggles.port import PortError, PymolPort, call
from mcpymol.wiggles.provenance import provenance_banner

#: Phrasings that would turn a gap into an absence claim. I3 forbids all of
#: them, in any output this package produces about latent space. Checked by
#: :func:`contains_absence_claim` and asserted over every tier-3 report.
ABSENCE_CLAIMS: tuple[str, ...] = (
    "unvisited",
    "never visited",
    "unpopulated",
    "not populated",
    "no molecules",
    "no particles here",
    "forbidden state",
    "forbidden region",
    "does not occur",
    "never occupied",
    "empty region means",
    "absent conformation",
)

GAP_LEGEND = (
    "A GAP IS NOT AN ABSENCE. An empty region between frames is a region of "
    "UNKNOWN OCCUPANCY, not a state the molecule avoids. In a blind challenge on "
    "41 submissions the commonest failure was missing an intermediate state that "
    "was genuinely present — usually the middle one of three — and there is an "
    "information-theoretic floor underneath that result, not just method error. "
    "Read this trajectory as: these conformations are supported. Read nothing "
    "into what lies between them."
)

POPULATION_LEGEND = (
    "NO POPULATION IS SHOWN, deliberately. Cryo-EM heterogeneity methods recover "
    "motion reliably and relative populations unreliably; those are different "
    "claims and a viewer conflates them constantly. This tool renders the motion. "
    "It draws no latent scatter and estimates no density, because the rendering "
    "users most want from latent space is the one the evidence says is most often "
    "wrong."
)


def contains_absence_claim(text: str) -> str | None:
    """The first I3-breaching phrase in ``text``, or None.

    Public so the invariant can be asserted over any report this package grows
    later, not only the ones that exist today.
    """
    lowered = text.lower()
    return next((phrase for phrase in ABSENCE_CLAIMS if phrase in lowered), None)


def frame_levels(ensemble: Ensemble, absolute: float) -> list[float | None]:
    """One sigma level per frame for a shared absolute contour.

    ``None`` for a frame whose header reports zero RMS, where sigma is
    undefined — reported rather than silently substituted, because a
    substituted level would draw a frame at a contour nobody chose.
    """
    levels: list[float | None] = []
    for header in ensemble.headers:
        try:
            levels.append(to_sigma(header, absolute))
        except ValueError:
            levels.append(None)
    return levels


def latent_traverse_view(
    port: PymolPort,
    ensemble_name: str,
    *,
    level: float | None = None,
    units: str = "sigma",
    name: str | None = None,
    color: str = "skyblue",
    build_movie: bool = True,
) -> str:
    """Build a steppable isosurface trajectory through an ensemble's frames.

    Args:
        port: A live or fake PyMOL port.
        ensemble_name: An ensemble loaded through ``load_ensemble``.
        level: Contour level. ``None`` uses :data:`DEFAULT_SIGMA` against the
            first frame whose header can define a sigma scale, then holds that
            **absolute** value across all frames.
        units: ``"sigma"`` (interpreted against the first frame with a usable
            RMS) or ``"absolute"`` (no conversion needed, so no anchor).
        name: Prefix for the isosurfaces. Defaults to ``<ensemble>_surf``.
        color: Colour for every frame. One colour on purpose — a per-frame
            spectrum would encode frame index as if it were a measured quantity.
        build_movie: Wire the frames to PyMOL's movie timeline so they can be
            stepped and played.

    Returns:
        A report: the method and its caveat, the shared absolute level with the
        per-frame sigmas it corresponds to, and the I2/I3 legends.

    Raises:
        PortError: no such ensemble, or its method was never identified (I2).
        ValueError: ``units`` is not recognised, or the level cannot be
            converted.
    """
    if units not in ("sigma", "absolute"):
        raise ValueError(f"units must be 'sigma' or 'absolute', got {units!r}")

    ensemble = loaded_ensemble(ensemble_name)
    if ensemble is None:
        raise PortError(
            f"{ensemble_name!r} was not loaded through load_ensemble, so its method, "
            f"frame order and headers are unavailable. Every one of those is needed "
            f"here: the method to satisfy invariant I2, the order to make the "
            f"trajectory mean anything, and the headers to hold one contour level "
            f"across frames PyMOL normalised separately."
        )

    # I2, the refusal. An unlabelled latent rendering is the thing the
    # invariant exists to prevent, so this is the tool working, not failing.
    if ensemble.method is Method.UNKNOWN:
        return "\n".join(
            [
                f"latent_traverse_view({ensemble_name})",
                "",
                "  REFUSED: this ensemble's generating method was never identified,",
                "  so no interpretive caveat can be attached to what would be drawn.",
                "",
                "  SPEC invariant I2 — no unlabelled latent plot. The caveats are not",
                "  interchangeable and they are not decoration: cryoDRGN's latent",
                "  density can bear no relation to the truth, while 3DVA's frames are",
                "  linear interpolations no particle ever occupied. Rendering these",
                "  identically would assert something for one that holds only for the",
                "  other, and the picture looks equally convincing either way.",
                "",
                "  " + "\n  ".join(ensemble.evidence),
                "",
                "  If you know what produced this directory, say so and reload:",
                f"      load_ensemble(path, {ensemble_name!r}, method='cryodrgn')",
                "  Accepted: " + ", ".join(m.value for m in Method if m is not Method.UNKNOWN),
            ]
        )

    no_usable_rms = PortError(
        f"no frame of {ensemble_name!r} has a usable RMS in its header, so an "
        f"absolute contour level cannot be converted to the sigma PyMOL "
        f"contours in. The headers may be stale or the maps unnormalised."
    )

    # Anchor the absolute level on the first header that can define a sigma
    # scale, not simply on frame 0. A header carrying rms=0, or MRC's rms=-1
    # "statistics not computed" marker, converts nothing — and taking it
    # anyway used to yield an anchor of dmean, a number unrelated to the level
    # asked for and then applied to every other frame. A level already given in
    # absolute units needs no anchor at all.
    anchor = next((h for h in ensemble.headers if usable_rms(h)), None)

    if level is None:
        if anchor is None:
            raise no_usable_rms
        absolute = to_absolute(anchor, DEFAULT_SIGMA)
        used_default = True
    elif units == "absolute":
        absolute = float(level)
        used_default = False
    else:
        if anchor is None:
            raise no_usable_rms
        absolute = to_absolute(anchor, float(level))
        used_default = False

    levels = frame_levels(ensemble, absolute)
    usable = [
        (obj, sigma)
        for obj, sigma in zip(ensemble.objects, levels, strict=True)
        if sigma is not None
    ]
    if not usable:
        raise no_usable_rms

    prefix = name or f"{ensemble_name}_surf"
    width = max(2, len(str(ensemble.n_frames)))
    surfaces: list[str] = []
    for index, (obj, sigma) in enumerate(usable, start=1):
        surface = f"{prefix}_{index:0{width}d}"
        call(port, "isosurface", surface, obj, sigma)
        call(port, "color", color, surface)
        surfaces.append(surface)

    if build_movie and len(surfaces) > 1:
        call(port, "mset", f"1 x{len(surfaces)}")
        for index, surface in enumerate(surfaces, start=1):
            # mdo attaches a command to a frame. It is a cmd function, so this
            # is still a structured call — the do() prohibition is about
            # cmd.do, which reports success for any string at all.
            call(port, "mdo", index, f"disable {prefix}_*; enable {surface}")
        call(port, "frame", 1)

    return _report(ensemble, prefix, surfaces, absolute, levels, used_default, build_movie)


def _report(
    ensemble: Ensemble,
    prefix: str,
    surfaces: list[str],
    absolute: float,
    levels: list[float | None],
    used_default: bool,
    build_movie: bool,
) -> str:
    known = [s for s in levels if s is not None]
    lo, hi = min(known), max(known)
    skipped = len(levels) - len(known)

    lines = [
        f"latent_traverse_view({ensemble.name})",
        "",
        f"  Method: {ensemble.method.label}",
        f"  {ensemble.method.caveat}",
        "",
        f"  {len(surfaces)} frames as `{surfaces[0]}` … `{surfaces[-1]}`, one colour.",
    ]
    if build_movie and len(surfaces) > 1:
        lines += [
            "  Wired to the movie timeline: `mplay` to run it, `frame N` to step.",
        ]

    lines += [
        "",
        f"  Contour: {absolute:.6g} absolute, held constant across every frame.",
        f"  That is {lo:.3g}–{hi:.3g} sigma depending on the frame.",
    ]

    spread = (hi - lo) / max(abs(hi), abs(lo), 1e-12)
    if spread > 0.05:
        lines += [
            "",
            f"  Those sigmas differ by {spread * 100:.0f}% across the trajectory, and",
            "  that difference IS the density change this traversal shows. PyMOL",
            "  normalises every map independently on load, so contouring each frame",
            "  at a fixed sigma would have contoured each at its own scale and",
            "  flattened exactly this signal — the frames would look reassuringly",
            "  alike whatever the ensemble does.",
        ]
    else:
        lines += [
            "",
            "  Those sigmas barely differ, so the frames carry near-identical",
            "  density statistics. The level is still held in absolute terms, which",
            "  is what makes that a finding rather than an artefact of normalising",
            "  each frame separately.",
        ]

    if used_default:
        lines += [
            "",
            f"  No level given, so {DEFAULT_SIGMA} sigma against the FIRST frame was",
            "  used and converted to an absolute value. A generic starting point,",
            "  not a recommendation for this ensemble.",
        ]

    if skipped:
        lines += [
            "",
            f"  {skipped} frame(s) were skipped: their headers report rms=0, so no",
            "  sigma level exists to draw them at. They are loaded but not contoured.",
        ]

    lines += [
        "",
        provenance_banner(ensemble.objects[0]),
        "",
        "  " + POPULATION_LEGEND,
        "",
        "  " + GAP_LEGEND,
    ]
    return "\n".join(lines)
