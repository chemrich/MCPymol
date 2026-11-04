"""3D-print export: watertight, multi-colour STL out of a PyMOL scene.

PyMOL's open-source build cannot write STL. Its OBJ exporter dumps the whole
visible scene as non-manifold triangle soup, in camera coordinates, and leaves
out everything PyMOL draws as spheres and cylinders: sticks, spheres, lines,
distance dashes. This module isolates each colour group, exports it in PyMOL's
model frame, takes the sticks and spheres from a POV-Ray export when asked to,
and rebuilds one watertight solid per group, so every group lands in the same
frame whatever the camera was doing.
"""

import functools
import math
import os
import re
import shutil
import tempfile
import typing
from collections import Counter
from collections.abc import Iterable
from typing import Annotated, Any, NamedTuple, cast

from pydantic import Field

from mcpymol.app import mcp
from mcpymol.bridge import _SLOW_OP_TIMEOUT, send_request

if typing.TYPE_CHECKING:
    import numpy as np
    import trimesh

# ── 3D printing export ───────────────────────────────────────────────────────
#
# PyMOL's open-source build cannot write STL, and its OBJ exporter writes the
# whole *visible* scene (it ignores the selection argument). The mesh it does
# produce is a non-manifold surface soup — unprintable as-is. This tool works
# around all three: it isolates each colour group, exports it, then rebuilds a
# watertight solid per group in PyMOL's model coordinates, so the groups
# assemble correctly as multi-material parts in a slicer.

_PRINT_DEPS_HINT = (
    "3D-print export needs extra libraries. Install them with:\n"
    "  uv sync --extra print     (from a MCPymol checkout)\n"
    "  uv pip install 'mcpymol[print]'\n"
    "Provides: trimesh, pymeshlab, scipy, scikit-image, manifold3d."
)

_REPRESENTATIONS = ("surface", "cartoon", "displayed")
_METHODS = ("auto", "poisson", "voxel")
_BODY_POLICIES = ("all", "largest")

# ── the export frame ──
#
# PyMOL's OBJ and POV-Ray writers are ray-tracer back ends and write camera
# space, with R = reshape(view[:9], (3, 3)) and x the model coordinates that
# get_coords reports:
#   OBJ  v = (x - view[12:15]) @ R + (view[9], view[10], 0)
#   POV  p = (x - view[12:15]) @ R + view[9:12]
# Turning the camera 90 degrees between two exports swapped x and z in the
# split-GFP build. An identity rotation with the camera and the origin of
# rotation at 0 makes the OBJ map the identity, so the file comes out in model
# coordinates: measured exact to 5e-6 Å, under geometry_export_mode 0 and 1
# alike (a camera left at z=-150 was 150 Å off under mode 1). Clipping,
# orthoscopic, field of view and the viewport do not enter, so the user's
# view[15:18] are kept; view[17] *sets* orthoscopic and field_of_view.
#
# The POV-Ray file keeps the camera at the user's distance instead. Dashes,
# lines, ribbon, mesh and dots are sized in pixels (their *_radius settings
# default to 0), and PyMOL turns pixels into Å from the camera: with it at
# the origin every such rod came out with radius NaN (perspective) or 0
# (orthoscopic). At the user's distance they come out exactly as the user's
# own ray trace sizes them (0.0975 Å for a 2.5-pixel dash, either way), and
# the file is model coordinates moved by that distance along z, which
# _parse_pov takes back out. That holds under geometry_export_mode 0 only
# (mode 1 left them 67 Å off), so the POV is written under mode 0.
_FRAME_ATTEMPTS = 3

# ── bodies ──
#
# PyMOL geometry is legitimately several bodies: chains that do not touch
# (4hhb's two dimers in print_ribbon_view), separate ligands (its four hemes),
# and the inward-wound cavity shells a default surface encloses. Keeping only
# the body with the most faces dropped one dimer, or three of four hemes, and
# still said "OK"; which dimer survived flipped with the voxel pitch.
_MIN_BODY_FACES = 8  # PyMOL debris was 1-2 faces; marching-cubes slivers 2 faces, no volume
_DEBRIS_VOXELS = 8  # a voxel body this small cannot stand for anything real
_LARGEST_DROP_WARNING = 0.01  # bodies="largest" warns when it drops more than this share

# ── voxel growth ──
#
# trimesh's voxeliser marks every cell a surface touches, and marching cubes
# then puts the surface half-way to the next empty cell, so every surface came
# out 0.48-0.71 voxel (mean 0.62) outside PyMOL's: a 0.33 Å rod printed 0.38 Å
# at 0.1 Å pitch, 1ubq's cartoon came out at 4x its volume at 0.7, and each
# re-voxelisation grew it again. Moving every vertex back in along its normal
# by the mean undoes it. Where a feature is thinner than that the inset is
# halved, up to 8 times, on the vertices of any face it would turn inside out;
# that left 0 flipped faces on 1ubq, 4hhb and an 843k-face ribbon.
_VOXEL_GROWTH = 0.62
_INSET_ROUNDS = 8
# A flat face moved along its own normal never turns over, so that check
# missed thin parts whose two sides passed through each other: a 0.3 Å slab
# came out inside out, and at 0.7 Å the inset took 3% of 1ubq's cartoon
# vertices past half their part's thickness. So no vertex moves more than half
# the part's thickness along its normal (read from the filled voxel grid),
# less this much, which leaves at least 0.1 voxel.
_INSET_MARGIN = 0.05  # voxels
# Gaps narrower than about this many voxels fuse while voxelising: the inset
# cannot reopen them (a 0.53 Å gap fused at 0.2 Å pitch, 0.78 Å at 0.35).
_GAP_FUSE_VOXELS = 2.5

# ── sticks and spheres ──
#
# Each sphere and cylinder becomes a polyhedron that encloses the true shape
# (circumscribed), so a stick never prints thinner than PyMOL drew it. An
# "open" cylinder is a half-bond that PyMOL caps with a sphere, or meets the
# other half face to face at the bond midpoint. Those halves are joined into
# one rod (_merge_continuing_rods); every open end is also pushed this much
# further, so that an end meeting other geometry face to face overlaps it
# rather than only touching it (before the halves were joined, touching ends
# split MK1's sticks into 6 pieces; with them joined, no case tried needed it).
_PRIMITIVE_SEGMENTS = 32
_OPEN_END_OVERLAP = 0.02  # Å, at each open end
# The exact union is fast for ligands and side chains (MK1: 0.8 s; GFP's
# cartoon plus 3,425 side-chain primitives: 20 s, 1.9 GB) and grows faster
# than linearly: 8,000 spheres took 10 s and 2.9 GB, and 1ubq's 38,567 dots
# ten minutes and a 1 GB STL. Dots and mesh are refused outright in displayed
# mode; this caps everything else (thousands of side chains as sticks).
_MAX_PRIMITIVES = 10_000

# ── thin features ──
#
# Under 3 voxels across, a feature prints lumpy (its radius varies by about
# +-0.35 voxel); under 1 voxel it comes out about one voxel thick whatever its
# real size. 1ubq's default loops are 0.4 Å across: 0.6 voxel at 0.7 Å.
_THIN_VOXELS = 3.0
# Rods and spheres thinner than two widths of a 0.4 mm nozzle snap or vanish.
_THIN_ROD_MM = 0.8

# Slicers guess an STL's units from the volume it encloses, one file at a time,
# and scale a file below these volumes (mm^3) up as if it were in inches (25.4x)
# or metres (1000x). PrusaSlicer 2.9 uses 9.0 and 0.001, Bambu Studio and
# OrcaSlicer 8.0 and 0.008; the larger of each pair, so one warning covers all
# three. PrusaSlicer and OrcaSlicer ask first; Bambu Studio 2.8.4 does not.
_SLICER_INCHES_BELOW_MM3 = 9.0
_SLICER_METRES_BELOW_MM3 = 0.008

# Representations PyMOL draws with ray-tracer spheres and cylinders. The OBJ
# writer drops the cylinders and writes each sphere as one degenerate triangle
# at its centre, so only representation="displayed", which reads them from the
# POV-Ray export, carries them.
_PRIMITIVE_REPS = (
    "sticks",
    "spheres",
    "nb_spheres",
    "lines",
    "nonbonded",
    "ribbon",
    "dots",
    "mesh",
)
_CENSUS_REPS = ("cartoon", "surface", *_PRIMITIVE_REPS, "ellipsoids", "labels")

# Thickness PyMOL draws for each cartoon type, as (setting, factor, feature),
# measured across exported OBJs: a strand is sqrt(2) x cartoon_rect_width
# thick, a helix ribbon 2 x cartoon_oval_width, a loop 2 x cartoon_loop_radius.
_CARTOON_THICKNESS = {
    1: (("cartoon_loop_radius", 2.0, "loop tubes"),),
    2: (("cartoon_rect_width", math.sqrt(2), "strands"),),
    3: (("cartoon_oval_width", 2.0, "helix ribbons"),),
    4: (("cartoon_tube_radius", 2.0, "tubes"),),
    5: (("cartoon_rect_width", math.sqrt(2), "strand arrows"),),
    6: (
        ("cartoon_dumbbell_width", math.sqrt(2), "dumbbell helices"),
        ("cartoon_dumbbell_radius", 2.0, "dumbbell edges"),
    ),
    7: (("cartoon_putty_radius", 2.0, "putty tubes"),),
}

# PyMOL's POV-Ray writer: one primitive per line, coordinates first.
_POV_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_POV_VEC = rf"<\s*({_POV_NUM})\s*,\s*({_POV_NUM})\s*,\s*({_POV_NUM})\s*>"
_POV_SPHERE_RE = re.compile(rf"^sphere\{{{_POV_VEC},\s*({_POV_NUM})", re.M)
_POV_CYLINDER_RE = re.compile(
    rf"^cylinder\{{{_POV_VEC},\s*{_POV_VEC},\s*({_POV_NUM})\s*(open)?", re.M
)
_POV_KIND_RE = re.compile(r"^([a-z_][a-z_0-9]*)\s*\{", re.M)
# Everything else PyMOL's writer emits at the start of a line. Triangles come
# from the OBJ export, which carries the same ones, so they are only named
# here; a kind outside this set is geometry print_export cannot read.
_POV_KNOWN_KINDS = frozenset(
    {
        "camera",
        "light_source",
        "plane",
        "sphere",
        "cylinder",
        "mesh2",
        "smooth_color_triangle",
        "smooth_triangle",
        "triangle",
        "texture_list",
        "face_indices",
        "pigment",
        "texture",
    }
)
_OBJ_VERTEX_RE = re.compile(r"^v\s+(\S+)\s+(\S+)\s+(\S+)", re.M)
_OBJ_FACE_RE = re.compile(r"^f\s+(\d+)\S*\s+(\d+)\S*\s+(\d+)", re.M)


# A label names the group's STL file, so it is kept to characters every file
# system takes. Labels are compared ignoring case: 'A' and 'a' are one file on
# macOS's default file system, so the second group silently overwrote the
# first (and they were one PyMOL object name too, while temp objects were
# named after labels).
_LABEL = re.compile(r"[A-Za-z0-9_.+-]+")


def _parse_groups(groups: str) -> "list[tuple[str, str]]":
    """Parse ``"label=selection; label2=selection2"`` into ordered pairs."""
    pairs = []
    seen: dict[str, str] = {}
    for chunk in groups.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"bad group spec {chunk!r}; expected 'label=selection'")
        label, sel = chunk.split("=", 1)
        label, sel = label.strip(), sel.strip()
        if not label or not sel:
            raise ValueError(f"bad group spec {chunk!r}")
        if not _LABEL.fullmatch(label):
            raise ValueError(
                f"group label {label!r} may use only letters, digits and _ . + - "
                f"(it becomes part of the STL file name), e.g. {re.sub(r'[^A-Za-z0-9_.+-]+', '_', label)!r}"
            )
        if label.casefold() in seen:
            raise ValueError(
                f"group labels {seen[label.casefold()]!r} and {label!r} would write the same "
                f"STL file (file names ignore case on macOS); give each group its own label"
            )
        seen[label.casefold()] = label
        pairs.append((label, sel))
    if not pairs:
        raise ValueError("no groups given")
    return pairs


# ── result text ──────────────────────────────────────────────────────────────


def _format_dims(lo: Iterable[float], hi: Iterable[float]) -> str:
    """``"25.1 x 14.2 x 18.0"``: the extent along x, y and z."""
    return " x ".join(f"{top - bottom:.1f}" for bottom, top in zip(lo, hi, strict=True))


def _format_volume(volume: float) -> str:
    """mm³ below 1000, cm³ from there up: ``"240 mm³"``, ``"86.29 cm³"``."""
    return f"{volume / 1000:,.2f} cm³" if volume >= 1000 else f"{volume:.3g} mm³"


def _describe_size(bounds: Iterable[Iterable[float]], volume: float) -> str:
    """``"25.1 x 14.2 x 18.0 mm, 1.23 cm³"`` for a group's line in the result."""
    return f"{_format_dims(*bounds)} mm, {_format_volume(volume)}"


def _unit_guess_warning(label: str, volume: float, mm_per_angstrom: float) -> "str | None":
    """Warn when a slicer would take this part's STL for inches or metres.

    Slicers judge each imported file on its own, so a small part loaded with a
    large one is the case that goes wrong: only the small one is enlarged.
    """
    if volume >= _SLICER_INCHES_BELOW_MM3:
        return None
    if volume < _SLICER_METRES_BELOW_MM3:
        guess = (
            "as metres and scale it 1000x (Bambu Studio and OrcaSlicer below "
            "0.008 mm³, PrusaSlicer below 0.001 mm³)"
        )
    else:
        guess = (
            "as inches and scale it 25.4x (PrusaSlicer below 9 mm³, Bambu Studio "
            "and OrcaSlicer below 8 mm³)"
        )
    text = (
        f"WARNING: '{label}' encloses only {volume:.3g} mm³. Slicers read an STL "
        f"this small {guess}, judging each file on its own, so it would load out "
        f"of scale with the other parts. PrusaSlicer and OrcaSlicer ask first; "
        f"Bambu Studio's 2.8.4 beta does it without asking."
    )
    if volume > 0:
        # Volume grows with the cube of the scale.
        needed = mm_per_angstrom * (_SLICER_INCHES_BELOW_MM3 / volume) ** (1 / 3)
        text += (
            f" Re-export with mm_per_angstrom >= {math.ceil(needed * 100) / 100:g}, "
            f"or check the part's size after import and undo the conversion."
        )
    return text


def _describe_bodies(info: dict) -> "list[str]":
    """The pieces a group prints as, and what the body policy filled or removed."""
    lines = []
    kept = info["bodies"]
    if len(kept) > 1:
        shown = ", ".join(_format_volume(v) for v in kept[:6])
        more = f" and {len(kept) - 6} more" if len(kept) > 6 else ""
        lines.append(
            f"prints as {len(kept)} separate pieces: {shown}{more} "
            f"(bodies='largest' keeps only the first)"
        )
    if info["dropped"]:
        dropped = info["dropped"]
        share = sum(dropped) / (sum(dropped) + info["volume"])
        text = (
            f"kept the largest piece and dropped {len(dropped)} more "
            f"({_format_volume(sum(dropped))}, {share:.1%} of the group)"
        )
        if share > _LARGEST_DROP_WARNING:
            text = (
                f"WARNING: {text}. They are missing from the STL; "
                f"bodies='all' keeps them as separate pieces."
            )
        lines.append(text)
    if info["cavities"]:
        n = len(info["cavities"])
        lines.append(
            f"filled {n} enclosed {'cavity' if n == 1 else 'cavities'} "
            f"({_format_volume(sum(info['cavities']))})"
        )
    if info["debris"]:
        # Mostly zero-volume slivers (a boolean union's, marching cubes'), so
        # say how little they held.
        n = len(info["debris"])
        lines.append(
            f"removed {n} stray {'fragment' if n == 1 else 'fragments'} of debris "
            f"({_format_volume(sum(info['debris']))} in all)"
        )
    return lines


# ── talking to PyMOL ─────────────────────────────────────────────────────────


def _result_list(res: dict) -> list:
    """A list result, or ``[]``.

    The plugin turns a ``None`` return into the *string* "Executed
    'get_object_list' successfully.", and iterating that enabled its letters
    one by one: PyMOL matches unique name prefixes, so ``enable e`` enabled an
    unrelated object, which was then exported and reported OK.
    """
    value = res.get("result") if res.get("status") == "success" else None
    return value if isinstance(value, list) else []


def _numeric_result(res: dict, default: float = 0.0) -> float:
    """A numeric result, or ``default`` for an error or a non-number."""
    value = res.get("result") if res.get("status") == "success" else None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return default
    return float(value)


def _setting_float(name: str, obj: str) -> "float | None":
    res = send_request("get_setting_float", args=[name, obj])
    value = _numeric_result(res, math.nan)
    return None if math.isnan(value) else value


def _setting_on(name: str, obj: str) -> bool:
    # cmd.get answers "on"/"off" for a boolean; get_setting_boolean answers a bool.
    res = send_request("get_setting_boolean", args=[name, obj])
    return res.get("status") == "success" and res.get("result") is True


def _enabled_names() -> "list[str]":
    """Every enabled object, group and selection, to put back exactly."""
    res = send_request("get_names", args=["all", 1])
    if res.get("status") == "error":
        raise RuntimeError(f"could not read which objects are enabled: {res.get('error')}")
    return [str(name) for name in _result_list(res)]


def _restore_enabled(names: "list[str]") -> None:
    # Exactly the set that was enabled: the old `enable all` also re-enabled
    # every object the user had switched off. Groups do not cascade either
    # way, so enabling each name restores the set.
    send_request("disable", args=["all"])
    for name in names:
        send_request("enable", args=[name])


def _show_only(names: "list[str]") -> None:
    """Make ``names`` the only enabled objects: PyMOL exports the whole visible scene."""
    send_request("disable", args=["all"])
    for name in names:
        # parents=1: an object inside a PyMOL group draws, and exports, nothing
        # while its group is disabled.
        send_request("enable", args=[name], kwargs={"parents": 1})
    # A grouped object's visibility updates lazily: without this, the export
    # straight after a toggle came out stale in both directions.
    send_request("refresh")


def _pause_motion() -> "dict[str, bool]":
    """Switch stereo and rock off for the export; returns what to switch back on.

    Stereo exports one eye's view (1.2 degrees and up to 0.9 Å off), and with
    the camera at the origin it writes NaN. Rock turns the camera every tick.
    """
    paused = {}
    for setting, off in (("stereo", "off"), ("rock", 0)):
        res = send_request("get", args=[setting])
        on = res.get("status") == "success" and str(res.get("result")).strip().lower() == "on"
        if on:
            send_request(setting, args=[off])
        paused[setting] = on
    return paused


def _resume_motion(paused: "dict[str, bool]") -> None:
    if paused.get("rock"):
        send_request("rock", args=[1])
    if paused.get("stereo"):
        send_request("stereo", args=["on"])


def _model_frame_view(view: "list[float]", camera_z: float = 0.0) -> "list[float]":
    """Identity rotation, origin of rotation at 0, the camera ``camera_z`` along z.

    The user's clipping and lens (view[15:18]) are kept.
    """
    return [
        *(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        *(0.0, 0.0, float(camera_z)),
        *(0.0, 0.0, 0.0),
        *view[15:18],
    ]


def _read_view() -> "list[float]":
    res = send_request("get_view")
    view = res.get("result") if res.get("status") == "success" else None
    if not isinstance(view, list) or len(view) != 18:
        raise RuntimeError(f"could not read the camera view: {res.get('error', view)!r}")
    return [float(v) for v in view]


def _save_in_model_frame(paths: "list[str]") -> float:
    """Write the visible scene to each path (.obj, .pov) in model coordinates (Å).

    Returns how far along z the POV-Ray file is moved (the user's camera
    distance, see the frame notes above), for _parse_pov to take back out.
    The user's view is put back straight after the saves, even when one
    fails, so the camera is only disturbed while PyMOL writes and never during
    the slow repair. If the camera moved meanwhile (a movie, a drag in the
    GUI), the frame cannot be trusted: try again, then give up with an error.
    """
    mode = 0.0
    if any(path.endswith(".pov") for path in paths):
        mode = _numeric_result(send_request("get_setting_float", args=["geometry_export_mode"]))
        if mode:
            send_request("set", args=["geometry_export_mode", 0])
    try:
        for _ in range(_FRAME_ATTEMPTS):
            view = _read_view()
            moved = False
            try:
                for path in paths:
                    pov = path.endswith(".pov")
                    target = _model_frame_view(view, view[11] if pov else 0.0)
                    res = send_request("set_view", args=[target])
                    if res.get("status") == "error":
                        raise RuntimeError(f"could not set the export view: {res.get('error')}")
                    # No selection: PyMOL's OBJ and POV writers ignore it and
                    # write the whole visible scene.
                    res = send_request("save", args=[path], timeout=_SLOW_OP_TIMEOUT)
                    if res.get("status") == "error":
                        raise RuntimeError(
                            f"PyMOL could not write {os.path.basename(path)}: {res.get('error')}"
                        )
                    after = _read_view()
                    if any(abs(a - b) > 1e-6 for a, b in zip(after[:15], target[:15], strict=True)):
                        moved = True
                        break
            finally:
                send_request("set_view", args=[view])
            if not moved:
                return view[11]
    finally:
        if mode:
            send_request("set", args=["geometry_export_mode", int(mode)])
    raise RuntimeError(
        "the camera moved while PyMOL was writing the geometry (rock, a movie "
        "playing, or the view being dragged?), so its frame can't be trusted. "
        "Stop it and export again."
    )


# A name in a selection expression: anything between spaces, parentheses and
# the operators ! & |.
_SELECTION_TOKEN = r"[^\s()!&|]+"


def _atoms_only(sel: str, names: Iterable[str]) -> str:
    """``sel`` with each of ``names`` (objects without atoms) swapped for "none".

    Named in a selection, a distance or CGO object makes PyMOL's selector
    reject the whole expression ("Invalid selection name").
    """
    for name in names:
        sel = re.sub(rf"(?<![^\s()!&|]){re.escape(name)}(?![^\s()!&|])", "none", sel)
    return sel


def _resolve_group(sel: str) -> "tuple[list[str], list[tuple[str, str]]]":
    """The molecular objects a group covers, and the other objects it names.

    get_object_list only reports objects with atoms in the selection, so a
    distance or CGO object has to be named in the group to be found. Left in
    the selection, it made PyMOL reject the whole expression, and "(lig or
    d1)" exported the distance without the ligand.
    """
    names = {str(name) for name in _result_list(send_request("get_names", args=["objects"]))}
    others = []
    for token in dict.fromkeys(re.findall(_SELECTION_TOKEN, sel)):
        if token not in names:
            continue
        kind = str(send_request("get_type", args=[token]).get("result", ""))
        if kind.startswith("object:") and kind not in ("object:molecule", "object:group"):
            others.append((token, kind.removeprefix("object:")))
    atoms = _atoms_only(sel, [name for name, _ in others])
    res = send_request("get_object_list", args=[f"({atoms})"])
    if res.get("status") == "error":
        raise RuntimeError(str(res.get("error")))
    molecules = [str(obj) for obj in _result_list(res)]
    # A typo in the selection: PyMOL prints a Selector-Error and the newer
    # plugin passes it on. Only that counts. Any other warning here was
    # printed by a redraw and landed in this request (KNOWN-ISSUES §5): taken
    # as a rejection, an Extrude-Warning would have failed a valid export.
    rejected = [str(w) for w in res.get("warnings") or [] if "Selector-Error" in str(w)]
    if not molecules and rejected:
        raise RuntimeError(f"PyMOL rejected the selection: {'; '.join(rejected)}")
    return molecules, others


# Bit of each census representation in an atom's `reps` mask: PyMOL's own
# cmd.repmasks (as of 3.1.0). The bridge cannot import pymol, so
# the table is written out; tests/test_live_printing.py checks it against a
# real PyMOL.
_REP_BITS = {
    "sticks": 1,
    "spheres": 2,
    "surface": 4,
    "labels": 8,
    "nb_spheres": 16,
    "cartoon": 32,
    "ribbon": 64,
    "lines": 128,
    "mesh": 256,
    "dots": 512,
    "nonbonded": 2048,
    "ellipsoids": 524288,
}


def _census(objs: "list[str]") -> "dict[str, dict[str, int]]":
    """How many atoms of each object show each representation.

    One request per object, reading each atom's ``reps`` mask and decoding it
    here. It was one count_atoms per object per representation: 13 sequential
    round trips for every object, before the export began.
    """
    census = {}
    for obj in objs:
        res = send_request("iterate_to_list", args=[f"({obj})", "reps"])
        rows = res.get("result") if res.get("status") == "success" else None
        masks = [int(row[0]) for row in rows] if isinstance(rows, list) else []
        census[obj] = {rep: sum(1 for m in masks if m & _REP_BITS[rep]) for rep in _CENSUS_REPS}
    return census


def _shown_reps(counts: "dict[str, int]", reps: Iterable[str]) -> str:
    return ", ".join(f"{rep} ({counts[rep]:,} atoms)" for rep in reps if counts.get(rep))


def _solid_dashes(obj: str) -> "tuple[str, float, bool]":
    """dash_gap 0 on a distance object for the export: one solid rod, not dashes.

    Returns what _restore_dashes needs to put it back.
    """
    old = _setting_float("dash_gap", obj)
    default = _numeric_result(send_request("get_setting_float", args=["dash_gap"]), math.nan)
    send_request("set", args=["dash_gap", 0, obj])
    value = old if old is not None else default
    return obj, value, old is not None and old != default


def _restore_dashes(saved: "tuple[str, float, bool]") -> None:
    obj, value, own = saved
    if own:
        send_request("set", args=["dash_gap", value, obj])
    else:
        # It followed the global value; an object-level copy would stop it doing so.
        send_request("unset", args=["dash_gap", obj])


def _thinnest_cartoon(objs: "list[str]") -> "tuple[float, str, float, str, str] | None":
    """The thinnest cartoon feature PyMOL draws for these objects.

    Returns (thickness Å, setting, its value, object, feature). Read from the
    cartoon type and secondary structure of each shown guide atom, so it costs
    a few round trips, not an export.
    """
    best = None
    for obj in objs:
        res = send_request(
            "iterate_to_list",
            args=[f"({obj}) and guide and rep cartoon and present", "cartoon, ss"],
        )
        kinds = {(int(row[0]), str(row[1])) for row in _result_list(res) if len(row) == 2}
        if not kinds:
            continue
        fancy = _setting_on("cartoon_fancy_helices", obj)
        cylinders = _setting_on("cartoon_cylindrical_helices", obj)
        types = set()
        for ctype, ss in kinds:
            if ctype == 0:  # automatic: by secondary structure
                if ss == "H":
                    ctype = 9 if cylinders else 6 if fancy else 3
                else:
                    ctype = 2 if ss == "S" else 1
            types.add(ctype)
        for ctype in sorted(types):
            for setting, factor, feature in _CARTOON_THICKNESS.get(ctype, ()):
                value = _setting_float(setting, obj)
                if value is None:
                    continue
                thickness = factor * value
                if best is None or thickness < best[0]:
                    best = (thickness, setting, value, obj, feature)
    return best


def _thin_cartoon_warning(
    thin: "tuple[float, str, float, str, str] | None", voxel_pitch: float
) -> "str | None":
    if thin is None:
        return None
    thickness, setting, value, obj, feature = thin
    if thickness >= _THIN_VOXELS * voxel_pitch:
        return None
    return (
        f"WARNING: the thinnest feature PyMOL draws here is {thickness:.2f} Å across "
        f"({feature} of '{obj}', {setting} {value:g}): {thickness / voxel_pitch:.1f} "
        f"voxels at voxel_pitch {voxel_pitch:g} Å. Under {_THIN_VOXELS:g} voxels it "
        f"prints lumpy, and under 1 it comes out about a voxel thick whatever its "
        f"real size. Use voxel_pitch <= {math.floor(thickness / _THIN_VOXELS * 100) / 100:g} "
        f"(slower), or thicken it with `set {setting}, <value>, {obj}`."
    )


def _thin_rod_warning(
    primitives: dict, census: "dict[str, dict[str, int]]", valence: bool, dashes: bool, s: float
) -> "str | None":
    """Name rods and spheres too thin to print at this scale, and what draws them."""
    import numpy as np

    radii = np.concatenate([primitives["cylinders"][:, 6], primitives["spheres"][:, 3]])
    radii = radii[radii > 0]
    thin = radii[2 * radii * s < _THIN_ROD_MM]
    if not len(thin):
        return None
    totals: Counter[str] = Counter()
    for counts in census.values():
        totals.update(counts)
    causes = []
    if totals["sticks"]:
        causes.append(
            "sticks (raise stick_radius"
            + (
                "; double bonds are drawn as thinner valence rods: `set valence, 0`)"
                if valence
                else ")"
            )
        )
    # Lines, ribbon, dots, mesh and dashes are sized in pixels unless their
    # *_radius setting is given in Å, so their thickness follows the zoom.
    if totals["lines"] or totals["nonbonded"]:
        causes.append(
            "lines and nonbonded crosses, sized in pixels (show them as sticks, or set line_radius)"
        )
    if totals["ribbon"]:
        causes.append("ribbon, sized in pixels (set ribbon_radius, or show cartoon)")
    if totals["dots"] or totals["mesh"]:
        causes.append("dots and mesh (hide them, or set dot_radius / mesh_radius)")
    if totals["spheres"] or totals["nb_spheres"]:
        causes.append("spheres (raise sphere_scale)")
    if dashes:
        causes.append("distance rods (set dash_radius; 0 sizes them in pixels)")
    smallest = float(thin.min())
    needed = math.ceil(_THIN_ROD_MM / (2 * smallest) * 100) / 100
    return (
        f"WARNING: {len(thin):,} of {len(radii):,} rods and spheres are under "
        f"{_THIN_ROD_MM:g} mm across at {s:g} mm per Å (the thinnest is "
        f"{2 * smallest * s:.2g} mm, radius {smallest:.2g} Å), too thin to print "
        f"reliably. They come from: {'; '.join(causes) or 'the displayed representations'}. "
        f"Or print larger: mm_per_angstrom >= {needed:g}."
    )


def _describe_empty_group(
    label: str,
    sel: str,
    objs: "list[str]",
    others: "list[tuple[str, str]]",
    census: "dict[str, dict[str, int]]",
    representation: str,
) -> str:
    """Why a group exported nothing, from what its objects show."""
    lines = [f"Error: group '{label}' ({sel}) exported no printable geometry."]
    for obj in objs:
        shown = _shown_reps(census[obj], _CENSUS_REPS) or "nothing (every representation is hidden)"
        lines.append(f"  {obj}: molecule; shows {shown}")
    for name, kind in others:
        lines.append(f"  {name}: {kind} object")
    totals: Counter[str] = Counter()
    for counts in census.values():
        totals.update(counts)
    hints = []
    primitive = [rep for rep in _PRIMITIVE_REPS if totals[rep]]
    primitive += [f"the {kind} object {name}" for name, kind in others]
    if representation == "cartoon" and primitive:
        hints.append(
            f"cartoon mode exports only triangles (cartoon, surface), and PyMOL "
            f"draws {', '.join(primitive)} as spheres and rods, which "
            f"representation='displayed' includes"
        )
    if totals["labels"]:
        hints.append("labels have no 3D geometry in any PyMOL export")
    if any(kind == "cgo" for _, kind in others):
        hints.append("PyMOL's exporters drop CGO cones and text")
    if not any(totals.values()) and not others:
        hints.append("show cartoon, sticks or surface for it first")
    if hints:
        lines.append("  Note: " + "; ".join(hints) + ".")
    return "\n".join(lines)


# ── reading PyMOL's exports ──────────────────────────────────────────────────


def _read_obj_triangles(path: str) -> "tuple[Any, int]":
    """The triangles of a PyMOL OBJ, (T, 3, 3) in Å, and how many sphere markers it held.

    PyMOL writes each sphere as one degenerate triangle at its centre. Those
    are dropped: voxelised they made a spheres-only group a single 0.2 Å³
    voxel reported "OK", where the real spheres enclosed 495 Å³.
    """
    import numpy as np

    with open(path) as fh:
        text = fh.read()
    vertices = np.array(_OBJ_VERTEX_RE.findall(text), dtype=float).reshape(-1, 3)
    if np.isnan(vertices).any():
        raise RuntimeError("PyMOL wrote NaN coordinates (is stereo on?)")
    faces = np.array(_OBJ_FACE_RE.findall(text), dtype=np.int64).reshape(-1, 3) - 1
    triangles = vertices[faces]
    marker = (triangles[:, 0] == triangles[:, 1]).all(axis=1) & (
        triangles[:, 1] == triangles[:, 2]
    ).all(axis=1)
    return triangles[~marker], int(marker.sum())


def _write_obj(path: str, triangles: "np.ndarray") -> None:
    import numpy as np

    vertices = triangles.reshape(-1, 3)
    with open(path, "w") as fh:
        np.savetxt(fh, vertices, fmt="v %.6f %.6f %.6f")
        np.savetxt(fh, np.arange(1, len(vertices) + 1).reshape(-1, 3), fmt="f %d %d %d")


def _parse_pov(path: str, z_shift: float = 0.0) -> dict:
    """The spheres and cylinders of a PyMOL POV-Ray export, in Å.

    ``z_shift`` is what _save_in_model_frame returned: the file's z is that
    far off model coordinates. Returns ``spheres`` (S, 4: x y z r),
    ``cylinders`` (K, 7: two ends and r), ``cyl_open`` (K: no flat caps, i.e.
    a half-bond PyMOL caps with spheres), and ``unknown``: counts of any
    primitive kind print_export cannot read.
    """
    import numpy as np

    with open(path) as fh:
        text = fh.read()
    spheres = np.array(_POV_SPHERE_RE.findall(text), dtype=float).reshape(-1, 4)
    found = _POV_CYLINDER_RE.findall(text)
    cylinders = np.array([c[:7] for c in found], dtype=float).reshape(-1, 7)
    kinds = Counter(_POV_KIND_RE.findall(text))
    # A NaN coordinate or radius does not match the number pattern at all.
    unusable = kinds["sphere"] - len(spheres) + kinds["cylinder"] - len(cylinders)
    unusable += int((spheres[:, 3] <= 0).sum() + (cylinders[:, 6] <= 0).sum())
    if unusable:
        raise RuntimeError(
            f"PyMOL wrote {unusable} spheres or rods with no usable position or size "
            f"(NaN or 0). Lines, dashes, ribbon, mesh and dots are sized in pixels "
            f"unless line_radius, dash_radius, ribbon_radius, mesh_radius or "
            f"dot_radius is set, and with the camera at the origin there is no pixel "
            f"size; set those radii in Å, or zoom out, and export again"
        )
    spheres[:, 2] -= z_shift
    cylinders[:, [2, 5]] -= z_shift
    return {
        "spheres": spheres,
        "cylinders": cylinders,
        "cyl_open": np.array([c[7] == "open" for c in found], dtype=bool),
        "unknown": {kind: n for kind, n in kinds.items() if kind not in _POV_KNOWN_KINDS},
    }


# ── mesh repair ──────────────────────────────────────────────────────────────


def _select_bodies(mesh: Any, pitch: "float | None", keep: str) -> "tuple[list, dict]":
    """Split a mesh into bodies and apply the body policy.

    Call before fix_normals: PyMOL and Poisson wind enclosed cavity shells
    inward, and that sign is what tells a cavity from a body. Debris (under 8
    faces, or under 8 voxels when a pitch applies) is dropped; an inward shell
    inside a solid's bounds is a cavity, dropped so it prints filled; an
    inward shell outside every solid is a solid wound the wrong way, so it is
    turned out and kept. Solids are ranked by volume, not face count.

    Returns (kept bodies, stats), stats holding the volumes (Å³) of the
    ``dropped`` solids, filled ``cavities`` and removed ``debris``.
    """
    floor = _DEBRIS_VOXELS * pitch**3 if pitch else 0.0
    debris, cavities, solids = [], [], []
    for body in mesh.split(only_watertight=False):
        volume = float(body.volume)
        if len(body.faces) < _MIN_BODY_FACES or abs(volume) <= floor:
            debris.append(abs(volume))
        elif volume < 0 and body.is_watertight:
            cavities.append(body)
        else:
            solids.append(body)
    filled = []
    for cavity in cavities:
        lo, hi = cavity.bounds
        if any((lo >= s.bounds[0]).all() and (hi <= s.bounds[1]).all() for s in solids):
            filled.append(-float(cavity.volume))
        else:
            cavity.invert()
            solids.append(cavity)
    solids.sort(key=lambda b: abs(float(b.volume)), reverse=True)
    kept = solids[:1] if keep == "largest" else solids
    dropped = [abs(float(b.volume)) for b in solids[len(kept) :]]
    return kept, {"dropped": dropped, "cavities": filled, "debris": debris}


def _closed_apart_from_debris(mesh: "trimesh.Trimesh") -> bool:
    """Watertight once stray fragments are set aside (1hsg's surface went to
    Poisson only because of three loose faces)."""
    if mesh.is_watertight:
        return True
    bodies = [b for b in mesh.split(only_watertight=False) if len(b.faces) >= _MIN_BODY_FACES]
    return bool(bodies) and all(b.is_watertight for b in bodies)


def _room_to_inset(
    vox: "trimesh.voxel.VoxelGrid", body: "trimesh.Trimesh", pitch: float
) -> "np.ndarray":
    """How far each vertex of a voxel body may move in: half the part's thickness.

    The thickness along the vertex's inward normal is where a walk from the
    vertex, in 0.05-voxel steps, first reaches an empty cell of the filled
    grid. Beyond twice the inset (plus a little) it cannot bind, so the walk
    stops there.
    """
    import numpy as np

    normals = np.asarray(body.vertex_normals, dtype=float)
    vertices = np.asarray(body.vertices, dtype=float)
    to_index = np.linalg.inv(np.asarray(vox.transform, dtype=float))
    filled = np.asarray(vox.matrix, dtype=bool)
    shape = np.array(filled.shape)
    thickness = np.full(len(vertices), np.inf)
    walking = np.ones(len(vertices), dtype=bool)
    for step in np.arange(1, int(2 * _VOXEL_GROWTH / 0.05) + 4) * 0.05:
        here = vertices - normals * (step * pitch)
        ijk = np.rint(here @ to_index[:3, :3].T + to_index[:3, 3]).astype(np.int64)
        inside = ((ijk >= 0) & (ijk < shape)).all(axis=1)
        inside[inside] = filled[ijk[inside, 0], ijk[inside, 1], ijk[inside, 2]]
        thickness[walking & ~inside] = step * pitch
        walking &= inside
    return np.maximum(thickness / 2 - _INSET_MARGIN * pitch, 0.0)


def _inset(body: "trimesh.Trimesh", distance: float, room: "np.ndarray | None" = None) -> int:
    """Move every vertex ``distance`` in along its normal, without folding.

    ``room`` caps each vertex's move (see _room_to_inset). On the vertices of
    any face the move would turn inside out, the inset is halved, up to
    _INSET_ROUNDS times, then dropped. Returns how many vertices moved less
    than ``distance``.
    """
    import numpy as np

    normals = np.array(body.vertex_normals, dtype=float)
    start = np.array(body.vertices, dtype=float)
    before = np.array(body.face_normals, dtype=float)
    real = np.linalg.norm(before, axis=1) > 0
    amount = np.full(len(start), float(distance))
    if room is not None:
        amount = np.minimum(amount, room)

    def folded() -> "np.ndarray":
        body.vertices = start - amount[:, None] * normals
        return real & (np.einsum("ij,ij->i", before, body.face_normals) <= 0)

    bad = folded()
    for _ in range(_INSET_ROUNDS):
        if not bad.any():
            break
        amount[np.unique(body.faces[bad])] *= 0.5
        bad = folded()
    while bad.any():
        # Ends: a face whose vertices have all stopped is back where it was.
        amount[np.unique(body.faces[bad])] = 0.0
        bad = folded()
    return int((amount < distance).sum())


@functools.lru_cache(maxsize=4)
def _sphere_circumscribe(segments: int) -> float:
    """The scale that makes manifold3d's geodesic sphere enclose the true one.

    Its vertices lie on the sphere and its faces cut inside it, by 1.3% of
    the radius at 32 segments.
    """
    import manifold3d as mf
    import numpy as np

    mesh = mf.Manifold.sphere(1.0, segments).to_mesh()
    vertices = np.asarray(mesh.vert_properties, dtype=float)[:, :3]
    tri = vertices[np.asarray(mesh.tri_verts, dtype=np.int64)]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    return 1.0 / float(np.abs(np.einsum("ij,ij->i", normals, tri[:, 0])).min())


def _merge_continuing_rods(
    cylinders: "np.ndarray", cyl_open: "np.ndarray", tol: float = 1e-4
) -> "tuple[np.ndarray, np.ndarray]":
    """Join open cylinders that continue each other into one, and drop duplicates.

    PyMOL draws a bond as two half-cylinders, one per atom colour, that meet
    end to end. Extended to overlap, the halves share their side walls, and a
    boolean union of coincident surfaces leaves zero-area slivers whose
    vertices coincide: in memory the mesh was closed, but the STL, re-merged
    by position as a slicer reads it, had 627 non-manifold edges and fell into
    573 pieces for MK1. One cylinder per bond has no coincident walls. A joint
    is merged only where exactly two ends meet, in a straight line, at the
    same radius.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    cyl = np.array(cylinders, dtype=float).reshape(-1, 7)
    opened = np.array(cyl_open, dtype=bool).reshape(-1)
    if not len(cyl):
        return cyl, opened
    # The same rod twice (two objects drawing one bond) coincides everywhere.
    # Compare them end-order independent: lower end first, lexicographically.
    swap = np.array([tuple(row[3:6]) < tuple(row[0:3]) for row in cyl], dtype=bool)
    canon = np.round(np.where(swap[:, None], cyl[:, [3, 4, 5, 0, 1, 2, 6]], cyl) / tol)
    _, first = np.unique(np.column_stack([canon, opened]), axis=0, return_index=True)
    keep = np.sort(first)
    cyl, opened = cyl[keep], opened[keep]
    while len(cyl) > 1:
        n = len(cyl)
        ends = np.concatenate([cyl[:, 0:3], cyl[:, 3:6]])
        tree = cKDTree(ends)
        degree = np.array([len(hits) for hits in tree.query_ball_point(ends, tol)])
        gone = np.zeros(n, dtype=bool)
        # One merge per rod per pass: a merged rod is re-oriented, so its old
        # end indices no longer say which end is which. Chains of more than
        # two take more passes.
        touched = np.zeros(n, dtype=bool)
        for p, q in tree.query_pairs(tol):
            i, j = p % n, q % n
            if i == j or touched[i] or touched[j] or degree[p] != 2 or degree[q] != 2:
                continue
            if not (opened[i] and opened[j]) or abs(cyl[i, 6] - cyl[j, 6]) > 1e-6:
                continue
            joint = ends[p]
            far_i = cyl[i, 3:6] if p < n else cyl[i, 0:3]
            far_j = cyl[j, 3:6] if q < n else cyl[j, 0:3]
            into, out = joint - far_i, far_j - joint
            if into @ out < (1 - 1e-6) * np.linalg.norm(into) * np.linalg.norm(out):
                continue
            cyl[i, 0:3], cyl[i, 3:6] = far_i.copy(), far_j.copy()
            gone[j] = touched[i] = touched[j] = True
        if not gone.any():
            break
        cyl, opened = cyl[~gone], opened[~gone]
    return cyl, opened


def _union_primitives(primitives: dict, solid: Any = None) -> "trimesh.Trimesh | None":
    """The exact boolean union of spheres, cylinders and an optional solid.

    Returns a trimesh in Å with its vertices merged by position, as an STL is
    read back, or None when ``solid`` is not a closed manifold manifold3d
    accepts. Voxelising primitives instead grew MK1's sticks by 62% at 0.1 Å
    (30% at 0.05, in 27 s); the union is within 0.4% in 0.8 s.
    """
    import manifold3d as mf
    import numpy as np
    import trimesh

    segments = _PRIMITIVE_SEGMENTS
    k_sphere = _sphere_circumscribe(segments)
    k_cylinder = 1.0 / math.cos(math.pi / segments)
    parts = []
    balls = np.array(primitives["spheres"], dtype=float).reshape(-1, 4)
    if len(balls):
        # Identical spheres (a joint drawn by two objects) coincide everywhere.
        balls = balls[np.sort(np.unique(np.round(balls / 1e-4), axis=0, return_index=True)[1])]
    for x, y, z, r in balls:
        if r > 0:
            ball = mf.Manifold.sphere(float(r) * k_sphere, segments)
            parts.append(ball.translate((float(x), float(y), float(z))))
    rods, rod_open = _merge_continuing_rods(primitives["cylinders"], primitives["cyl_open"])
    for row, is_open in zip(rods, rod_open, strict=True):
        start = np.array(row[0:3], dtype=float)
        axis = np.array(row[3:6], dtype=float) - start
        length = float(np.linalg.norm(axis))
        radius = float(row[6])
        if length < 1e-9 or radius <= 0:
            continue
        z_axis = axis / length
        if is_open:
            start = start - z_axis * _OPEN_END_OVERLAP
            length += 2 * _OPEN_END_OVERLAP
        x_axis = np.cross(z_axis, [1.0, 0.0, 0.0] if abs(z_axis[0]) < 0.9 else [0.0, 1.0, 0.0])
        x_axis /= np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)  # right-handed: a mirror would turn it inside out
        rod = mf.Manifold.cylinder(length, radius * k_cylinder, radius * k_cylinder, segments)
        parts.append(rod.transform(np.column_stack([x_axis, y_axis, z_axis, start])))
    if solid is not None and len(solid.faces):
        mesh = mf.Mesh(
            vert_properties=np.ascontiguousarray(solid.vertices, dtype=np.float32),
            tri_verts=np.ascontiguousarray(solid.faces, dtype=np.uint32),
        )
        mesh.merge()
        body = mf.Manifold(mesh)
        if body.status() != mf.Error.NoError:
            return None
        parts.append(body)
    union = mf.Manifold.batch_boolean(parts, mf.OpType.Add)
    out = union.to_mesh()
    # process=False: trimesh's processing deletes the slivers a boolean leaves
    # and opens the mesh. Vertices are still merged by position, which is all
    # an STL keeps, so watertightness is judged on what the file will hold.
    mesh = trimesh.Trimesh(
        np.asarray(out.vert_properties, dtype=float)[:, :3],
        np.asarray(out.tri_verts, dtype=np.int64),
        process=False,
    )
    mesh.merge_vertices()
    return mesh


# PoissonRecon, as pymeshlab runs it, ends the whole process with exit(0) when
# it gives up ("Failed to close loop"): 4 of 14 runs on one print-ribbon export
# did, at random, with no exception to catch, so in the server it took the MCP
# connection down mid-call. The code before this did the same. Run in a child
# process, a failure is only a missing file, which auto turns into voxel.
_POISSON_CHILD = (
    "import sys\n"
    # `python -c` puts the working directory first on sys.path, so a project
    # file named like a stdlib module (calendar.py, secrets.py) would be
    # imported, and run, in place of it: Poisson then failed every time. The
    # child also runs in its own scratch directory; this drops the entry anyway.
    "sys.path[:] = [p for p in sys.path if p]\n"
    "from mcpymol.printing import _poisson_reconstruct\n"
    "_poisson_reconstruct(sys.argv[1], sys.argv[2], int(sys.argv[3]))\n"
)


def _poisson_reconstruct(src: str, dst: str, depth: int) -> None:
    """Screened-Poisson rebuild of the mesh ``src`` into the STL ``dst``, in this process."""
    import pymeshlab

    def apply_first(ms: Any, names: list[str], **kw: Any) -> None:
        last = None
        for n in names:
            try:
                ms.apply_filter(n, **kw)
                return
            except Exception as e:
                last = e
        raise RuntimeError(f"none of {names} worked: {last}")

    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(src)
    apply_first(ms, ["meshing_remove_duplicate_vertices", "remove_duplicate_vertices"])
    apply_first(ms, ["meshing_remove_null_faces", "remove_zero_area_faces"])
    apply_first(ms, ["meshing_remove_unreferenced_vertices", "remove_unreferenced_vertices"])
    apply_first(ms, ["compute_normal_per_vertex", "re_orient_all_faces_coherently"])
    apply_first(
        ms,
        [
            "generate_surface_reconstruction_screened_poisson",
            "surface_reconstruction_screened_poisson",
        ],
        depth=depth,
        preclean=True,
    )
    ms.save_current_mesh(dst, binary=True)


def _run_poisson(src: str, dst: str, depth: int) -> None:
    """_poisson_reconstruct in a child process; RuntimeError when it wrote nothing."""
    import importlib.util
    import subprocess
    import sys

    if importlib.util.find_spec("pymeshlab") is None:
        raise ImportError("pymeshlab is not installed")
    if os.path.exists(dst):
        os.remove(dst)
    # The child imports mcpymol the way this process did.
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(p for p in sys.path if p)}
    try:
        proc = subprocess.run(
            [sys.executable, "-c", _POISSON_CHILD, src, dst, str(depth)],
            capture_output=True,
            text=True,
            env=env,
            timeout=_SLOW_OP_TIMEOUT,
            # The private directory the export writes its intermediates to,
            # never the server's own working directory (see _POISSON_CHILD).
            cwd=os.path.dirname(os.path.abspath(dst)),
            stdin=subprocess.DEVNULL,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"Poisson reconstruction took longer than {_SLOW_OP_TIMEOUT:g} s"
        ) from None
    if proc.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 0:
        return
    said = [line.strip() for line in (proc.stderr + proc.stdout).splitlines() if line.strip()]
    errors = [line for line in said if "error" in line.lower()]
    reason = (errors or said or [f"exit status {proc.returncode}"])[-1]
    raise RuntimeError(f"Poisson reconstruction failed and wrote nothing ({reason})")


def _repair_to_stl(
    src_obj: "str | None",
    dst_stl: str,
    method: str,
    voxel_pitch: float,
    poisson_depth: int,
    *,
    bodies: str = "all",
    undo_voxel_growth: bool = True,
    mm_per_angstrom: float = 1.0,
    primitives: "dict | None" = None,
) -> dict:
    """Rebuild PyMOL geometry into a watertight STL; returns what it did.

    ``src_obj`` is an OBJ of triangles (None if there are none) and
    ``primitives`` the spheres and cylinders of a POV-Ray export (see
    _parse_pov). ``method``: ``poisson`` keeps surface detail (good for bulky
    chains); ``voxel`` is robust for thin tubular geometry; ``auto`` uses the
    cheapest path that works: compact structures (e.g. a GFP barrel) often
    export already-watertight, so a light cleanup beats Poisson, which can
    degrade an already-closed surface. Primitives are added by an exact
    boolean union, never through the light path (overlapping closed shells
    pass as watertight and it kept one sphere of 45) or Poisson.

    Repair runs in Ångström, so ``voxel_pitch`` means the same thing at every
    print size. ``mm_per_angstrom`` is applied afterwards, about the origin
    rather than the part's own centre, so every group keeps its place.

    The info dict: ``method`` (with the voxel pitch and inset), ``faces``,
    ``watertight``, ``bounds`` and ``volume`` in mm, ``bodies`` (each kept
    piece, mm³, largest first), ``dropped``, ``cavities`` and ``debris``
    (volumes in mm³), ``inset`` (Å) and ``primitives`` (rods, spheres).
    """
    import trimesh

    if method not in _METHODS:
        raise ValueError(f"unknown method {method!r}")
    n_rods = len(primitives["cylinders"]) if primitives else 0
    n_balls = len(primitives["spheres"]) if primitives else 0
    if (n_rods or n_balls) and method == "poisson":
        raise ValueError(
            "Poisson reconstruction cannot merge sticks and spheres (it lost 78% of "
            "MK1's sticks); use method='auto' or 'voxel'"
        )
    if src_obj is None and not (n_rods or n_balls):
        raise ValueError("nothing to rebuild: no triangles and no spheres or rods")
    # With primitives, the policy waits for the union: a stick can join two
    # bodies that the triangles alone leave apart.
    policy = "all" if (n_rods or n_balls) else bodies

    def _finish(body: "trimesh.Trimesh") -> "trimesh.Trimesh":
        body.merge_vertices()
        body.update_faces(body.unique_faces())
        body.update_faces(body.nondegenerate_faces())
        body.remove_unreferenced_vertices()
        trimesh.repair.fix_normals(body)
        trimesh.repair.fill_holes(body)
        return body

    def _light(raw: "trimesh.Trimesh") -> "tuple[list[trimesh.Trimesh], dict]":
        kept, stats = _select_bodies(raw, None, policy)
        return [_finish(b) for b in kept], stats

    def _voxel(
        src: Any, pitch: float
    ) -> "tuple[list[trimesh.Trimesh], dict, trimesh.voxel.VoxelGrid]":
        # force="mesh" collapses a scene into a single Trimesh, but the
        # declared return type is the Geometry base class, so narrow it.
        m = cast("trimesh.Trimesh", trimesh.load(src, force="mesh"))
        vox = m.voxelized(pitch=pitch).fill()
        out = vox.marching_cubes
        # marching_cubes is in voxel-index space; map back to world coords.
        out.apply_transform(vox.transform)
        kept, stats = _select_bodies(out, pitch, policy)
        # The grid too: it says how thick each part is, for the inset.
        return [_finish(b) for b in kept], stats, vox

    def _poisson(src: Any, depth: int) -> "tuple[list[trimesh.Trimesh], dict]":
        # Written beside the source, not as the STL: the body policy still has
        # to fill the cavity voids Poisson keeps.
        rebuilt = f"{src}.poisson.stl"
        try:
            _run_poisson(src, rebuilt, depth)
            mesh = trimesh.load(rebuilt, force="mesh")
        finally:
            if os.path.exists(rebuilt):
                os.remove(rebuilt)
        kept, stats = _select_bodies(mesh, None, policy)
        return [_finish(b) for b in kept], stats

    def _triangles(how: str) -> "tuple[list[trimesh.Trimesh], dict, str, float]":
        """Repair the OBJ: (pieces, stats, method text, inset)."""
        if how == "auto":
            raw = cast("trimesh.Trimesh", trimesh.load(src_obj, force="mesh", process=True))
            if _closed_apart_from_debris(raw):
                return (*_light(raw), "light (already watertight)", 0.0)
            try:
                return (*_poisson(src_obj, poisson_depth), "poisson", 0.0)
            except Exception as e:
                # Say why: a voxel solid differs enough from a Poisson one
                # that a silent swap is worth explaining.
                pieces, stats, vox = _voxel(src_obj, voxel_pitch)
                why = str(e) if len(str(e)) <= 160 else f"{str(e)[:157]}..."
                used = f"voxel (poisson fallback: {why})"
        elif how == "poisson":
            return (*_poisson(src_obj, poisson_depth), "poisson", 0.0)
        else:
            pieces, stats, vox = _voxel(src_obj, voxel_pitch)
            used = "voxel"
        used += f", {voxel_pitch:g} Å pitch"
        growth = _VOXEL_GROWTH * voxel_pitch
        if not undo_voxel_growth:
            return pieces, stats, f"{used}, growth left in (about {growth:.2f} Å per surface)", 0.0
        held = sum(_inset(body, growth, _room_to_inset(vox, body, voxel_pitch)) for body in pieces)
        total = sum(len(body.vertices) for body in pieces)
        used += f", surfaces moved in {growth:.2f} Å to undo voxel growth"
        if held and total:
            used += f" (less on {held / total:.1%} of vertices, where a part is too thin for it)"
        return pieces, stats, used, growth

    pieces: list = []
    stats: dict = {"dropped": [], "cavities": [], "debris": []}
    used, inset = "", 0.0
    if src_obj is not None:
        pieces, stats, used, inset = _triangles(method)
    if n_rods or n_balls:
        solid = None
        if pieces:
            solid = pieces[0] if len(pieces) == 1 else trimesh.util.concatenate(pieces)
        union = _union_primitives(primitives or {}, solid)
        if union is None and src_obj is not None and not used.startswith("voxel"):
            # Light and Poisson output is not always a closed manifold; the
            # voxel path always is.
            pieces, stats, used, inset = _triangles("voxel")
            solid = pieces[0] if len(pieces) == 1 else trimesh.util.concatenate(pieces)
            union = _union_primitives(primitives or {}, solid)
        if union is None:
            raise RuntimeError("the triangle part is not a closed manifold, so it cannot be merged")
        pieces, union_stats = _select_bodies(union, None, bodies)
        for key, values in union_stats.items():
            stats[key] = stats[key] + values
        joined = f"exact union of {n_rods:,} rods and {n_balls:,} spheres"
        used = f"{used} + {joined}" if used else joined
    if not pieces:
        raise RuntimeError("nothing printable was left once stray fragments were removed")

    volumes = sorted((abs(float(b.volume)) for b in pieces), reverse=True)
    mesh = cast(
        "trimesh.Trimesh", pieces[0] if len(pieces) == 1 else trimesh.util.concatenate(pieces)
    )
    if mm_per_angstrom != 1.0:
        # A scale matrix with no translation: about (0, 0, 0), never the
        # part's own centre, or two groups would drift apart as they grew.
        mesh.apply_scale(mm_per_angstrom)
    mesh.export(dst_stl)
    lo, hi = (tuple(float(c) for c in corner) for corner in mesh.bounds)
    cube = mm_per_angstrom**3
    return {
        "method": used,
        "faces": len(mesh.faces),
        "watertight": bool(mesh.is_watertight),
        # Millimetres once scaled: what a slicer measures. A slicer flips an
        # inside-out mesh before it reads the volume, so report the magnitude.
        "bounds": (lo, hi),
        "volume": sum(volumes) * cube,
        "bodies": [v * cube for v in volumes],
        "dropped": [v * cube for v in stats["dropped"]],
        "cavities": [v * cube for v in stats["cavities"]],
        "debris": [v * cube for v in stats["debris"]],
        "inset": inset,
        "primitives": (n_rods, n_balls),
    }


# ── print_ribbon_view ────────────────────────────────────────────────────────
#
# The look is set on the object and its spine only. It used to be set with a
# bare `set`, which changes the global value, so every other object in the
# session, and every one loaded afterwards, took the chunky strands too, and
# two objects could not be styled differently. PyMOL reads each of these per
# object, and set per object they export byte-identical OBJs to the old
# global style (1ubq, 4hhb and 2b3p, main object and spine alike).
_RIBBON_STYLE = {
    "cartoon_loop_radius": "1.3",
    "cartoon_oval_width": "1.3",
    "cartoon_oval_length": "2.0",
    # Inert unless cartoon_cylindrical_helices is on; kept so a session with
    # it on still draws what it did.
    "cartoon_helix_radius": "1.3",
    "cartoon_fancy_sheets": "1",  # keep 3D arrowheads (strand direction)
    "cartoon_flat_sheets": "1",  # well-defined arrow plane
    "cartoon_smooth_loops": "1",
}
# PyMOL puts a strand's corners at +-cos 45° x cartoon_rect_width and
# +-sin 45° x cartoon_rect_length, so it draws strands sqrt(2) times both
# settings: the faces of the old "2.8 Å wide, 2.0 Å thick" strands measured
# 4.0 and 2.8 Å in the export. GFP's neighbouring strands have their Cα about
# 4.4 Å apart, so the barrel wall closed up: 8% of it was see-through.
_STRAND_DRAWN = math.sqrt(2)
# Spine residues kept at each end of a strand when spine_in_strands=False,
# where the loops fuse to the strands. With 1, 10 of GFP's 22 strand ends
# lost their junction and the barrel came apart into 6 pieces.
_MIN_STRAND_OVERLAP = 2
_RIBBON_VOXEL_PITCH = 0.2  # Å, for the export the result suggests
# The Cα-like atoms of a residue with no CA: GFP's chromophore (CRO, made
# from Thr65, Tyr66 and Gly67) has CA1, CA2 and CA3.
_CA_LIKE = re.compile(r"CA\d+")


def _ribbon_call(action: str, *args: object, notes: "list[str]") -> Any:
    """Send one PyMOL call and return its result; RuntimeError if PyMOL refused it.

    Every call is direct, and none goes through `do`: older plugins queue a
    `do` line onto PyMOL's main thread and answer before it has run, so its
    errors were lost (print_ribbon_view on a missing object reported
    success) and it could run after the calls sent behind it. Warnings
    PyMOL printed for a call it still ran are added to ``notes``.
    """
    res = send_request(action, args=list(args))
    warnings = [str(w) for w in res.get("warnings") or []]
    # A bad selection in iterate or orient is only a printed Selector-Error.
    # Other printed errors are not taken as this call's failure: PyMOL prints
    # some while it redraws (ObjectSurfaceUpdate-Error, Extrude-Warning), and
    # those can land in whichever request is running at the time.
    selector = [w for w in warnings if "Selector-Error" in w]
    if res.get("status") != "success" or selector:
        error = res.get("error") or "; ".join(selector)
        shown = ", ".join(str(a) for a in args)
        raise RuntimeError(f"`{action} {shown}` failed: {error}")
    notes.extend(warnings)
    return res.get("result")


def _by_index(obj: str, indices: Iterable[int]) -> str:
    """``(obj) and index 1+2+3``: atoms picked by index, so no chain or segment
    identifier ever has to be quoted into a selection."""
    return f"({obj}) and index {'+'.join(str(i) for i in indices)}"


class _CalessLink(NamedTuple):
    """A residue with no CA atom, peptide-linked on both sides, and how to bridge it.

    ``through`` is the atoms to trace across it (its CA-like atoms, else the
    atoms carrying its peptide links), ``before`` and ``after`` the Cα of the
    residues either side, and ``why_not`` why it cannot be bridged, or "".
    All are atom indices in the spine.
    """

    label: str
    before: "list[int]"
    through: "list[int]"
    after: "list[int]"
    why_not: str


def _caless_links(spine: str, notes: "list[str]") -> "list[_CalessLink]":
    """Residues with no CA atom that are peptide-linked on both sides.

    PyMOL joins two Cα guide atoms only when they are exactly three bonds
    apart (CA-C-N-CA), so a residue with no atom named CA breaks the cartoon:
    in GFP, CA64 and CA68 sit nine bonds apart through the chromophore, and
    the spine broke there into two pieces and a few loose dashes. Everything
    is found by atom index, so no chain or segment identifier ever has to be
    quoted into a selection.
    """
    no_ca = f"({spine}) and not (byres (({spine}) and name CA))"
    fields = "index, chain, resn, resi, segi"
    entries = _ribbon_call(
        "iterate_to_list",
        f"{no_ca} and (neighbor (({spine}) and polymer and name C))",
        fields,
        notes=notes,
    )
    exits = {
        tuple(row[1:]): int(row[0])
        for row in _ribbon_call(
            "iterate_to_list",
            f"{no_ca} and (neighbor (({spine}) and polymer and name N))",
            fields,
            notes=notes,
        )
    }
    found: dict[tuple, _CalessLink] = {}
    for row in entries:
        key = tuple(row[1:])
        if key in found or key not in exits:
            continue  # linked on one side only: a cap such as ACE or NH2
        chain, resn, resi, _segi = key
        residue = f"byres ({_by_index(spine, [row[0]])})"
        atoms = {
            int(index): str(name)
            for index, name in _ribbon_call(
                "iterate_to_list", f"({residue}) and present", "index, name", notes=notes
            )
        }
        ca_like = sorted(index for index, name in atoms.items() if _CA_LIKE.fullmatch(name))
        through = ca_like or sorted({int(row[0]), exits[key]} & set(atoms))
        before, after = (
            sorted(
                int(ca[0])
                for ca in _ribbon_call(
                    "iterate_to_list",
                    f"({spine}) and name CA and present and byres "
                    f"((({spine}) and polymer and name {link}) and neighbor ({residue}))",
                    "index",
                    notes=notes,
                )
            )
            for link in ("C", "N")
        )
        why_not = ""
        if not through:
            why_not = "it has no coordinates"
        elif not before or not after:
            why_not = "a neighbouring residue has no Cα coordinates"
        elif not before[-1] < through[0] <= through[-1] < after[0]:
            # PyMOL scans guide atoms in atom order, so a trace through atoms
            # stored elsewhere would join the wrong residues.
            why_not = "its atoms are not stored between its neighbours'"
        found[key] = _CalessLink(f"{resn} {chain}{resi}", before, through, after, why_not)
    return list(found.values())


def _traced_already(spine: str, notes: "list[str]") -> int:
    """How many of the spine's atoms PyMOL already draws as a trace.

    That is any atom with cartoon_trace_atoms on (create copies a per-atom
    setting from the source object) and any carbon Cα with no bonds, which
    PyMOL traces on its own: 1grl has 413 of those. A bridge needs
    trace_atoms_mode 0 on the spine, which joins every traced atom to the
    shown Cα before it, across chain breaks too.

    Only carbon counts: PyMOL takes an atom named CA as a guide atom only if
    it is a carbon, and a calcium ion is also named CA. Counting the four
    unbonded Ca2+ of a calcium sensor (GCaMP2, 3EVR) refused its chromophore
    bridge, though the ions are never traced.
    """
    rows = _ribbon_call(
        "iterate_to_list", f"({spine}) and present", "s.cartoon_trace_atoms", notes=notes
    )
    unbonded = _ribbon_call(
        "count_atoms",
        f"(({spine}) and name CA and elem C and present) and not (bound_to ({spine}))",
        notes=notes,
    )
    return sum(1 for (on,) in rows if on) + int(unbonded or 0)


def _strand_interiors(
    spine: str, overlap: int, breaks: "set[int]", notes: "list[str]"
) -> "list[list[int]]":
    """Cα indices of each β-strand's interior: its residues less ``overlap`` per end.

    A strand is a run of residues with ss S, each following the last (resv
    + 1, or a new insertion code) in one chain and segment. A run is cut
    before any Cα in ``breaks``, the residue after one with no CA: a bridge
    must stay shown on both sides, since trace_atoms_mode 0 joins it to the
    last *shown* Cα, wherever that is.
    """
    rows = _ribbon_call(
        "iterate_to_list",
        f"({spine}) and polymer and name CA and present",
        "chain, segi, resv, resi, ss, index",
        notes=notes,
    )
    runs: list[list[int]] = []
    run: list[int] = []
    last = None
    for chain, segi, resv, resi, ss, index in rows:
        here = (chain, segi, str(resi))
        if last is not None and here == last[0]:
            continue  # an alternate location of the same residue
        follows = (
            last is not None
            and here[:2] == last[0][:2]
            and (resv == last[1] + 1 or (resv == last[1] and here[2] != last[0][2]))
            and int(index) not in breaks
        )
        strand = str(ss).upper() == "S"
        if strand and follows and run:
            run.append(int(index))
        else:
            if run:
                runs.append(run)
            run = [int(index)] if strand else []
        last = (here, resv)
    if run:
        runs.append(run)
    return [run[overlap : len(run) - overlap] for run in runs if len(run) > 2 * overlap]


@mcp.tool()
def print_ribbon_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc").')],
    spine_radius: Annotated[
        float,
        Field(
            description="Radius (Å) of the continuous backbone tube. Larger values give more internal reinforcement. Default 0.9."
        ),
    ] = 0.9,
    strand_width: Annotated[
        float,
        Field(
            description="β-strand width setting (PyMOL's cartoon_rect_length). PyMOL draws the strand √2 times this wide: the default 2.8 is about 4.0 Å. 2.0 (2.8 Å) opens the gaps in a β-barrel."
        ),
    ] = 2.8,
    strand_thickness: Annotated[
        float,
        Field(
            description="β-strand thickness setting (PyMOL's cartoon_rect_width). PyMOL draws the strand √2 times this thick: the default 2.0 is about 2.8 Å, 1.0 about 1.4 Å."
        ),
    ] = 2.0,
    spine_in_strands: Annotated[
        bool,
        Field(
            description="True (default) runs the spine through the strands as well, where it shows as bumps along thin ones. False takes it out of each strand's interior, keeping strand_overlap residues at each end so the loops still fuse to the strands."
        ),
    ] = True,
    strand_overlap: Annotated[
        int,
        Field(
            description="With spine_in_strands=False, the residues of spine kept at each end of a strand, where it fuses with the loop. At least 2: with 1, 10 of GFP's 22 strand ends lost their junction."
        ),
    ] = 2,
    bridge_modified: Annotated[
        bool,
        Field(
            description="Carry the spine through residues that are peptide-linked on both sides but have no CA atom, such as GFP's chromophore (CRO), which otherwise break it. Done on the spine only, with per-atom cartoon_trace_atoms and PyMOL's no_smooth flag on the Cα either side; no atom is renamed or renumbered. Default true."
        ),
    ] = True,
) -> str:
    """
    Chunky β-arrow ribbons plus a continuous backbone "spine", tuned for
    rigid, gap-free 3D printing.

    Configures the look developed for FDM printing: thick β-strand arrows
    and a fat helix on the main object with loop cartoon hidden, plus a
    separate ``<obj>_spine`` object showing PyMOL's ``cartoon tube`` (which
    ignores secondary structure) running unbroken through the whole
    backbone. The spine threads through the strand bodies, so when the two
    objects are exported together the voxel step fuses them into one
    watertight solid with no strand→loop discontinuity. The spine also acts
    as internal rebar, reinforcing the thin junctions for print rigidity.

    The style is set on the object and its spine only; other objects and
    PyMOL's global settings are left alone. The spine is carried through
    residues with no CA atom, such as GFP's chromophore, without renaming
    anything.

    PyMOL draws a strand √2 times the strand_width and strand_thickness
    settings, so the defaults (2.8, 2.0) give strands about 4.0 x 2.8 Å, and
    a β-barrel's wall closes up. To see into a barrel, use
    ``strand_width=2.0, strand_thickness=1.0, spine_in_strands=False``.
    Measured on GFP (2b3p): the see-through share of the barrel wall goes
    from 8% to 29%; thinner strands alone leave 399 Å³ of spine bulging out
    of them, which spine_in_strands=False takes to 1 Å³; and it still
    exports as one solid (unbridged at the chromophore, it split 69:31).

    After calling this, export the fused solid with::

        print_export(obj_name="<obj>",
                     groups="<obj>=(<obj> or <obj>_spine)",
                     representation="cartoon",
                     method="voxel", voxel_pitch=0.2)
    """
    for name, size in (
        ("spine_radius", spine_radius),
        ("strand_width", strand_width),
        ("strand_thickness", strand_thickness),
    ):
        if not (math.isfinite(size) and size > 0):
            return f"Error: {name} must be a positive number of Å, got {size}."
    if strand_overlap < _MIN_STRAND_OVERLAP:
        return (
            f"Error: strand_overlap must be at least {_MIN_STRAND_OVERLAP} (with 1, 10 of "
            f"GFP's 22 strand ends lost their junction), got {strand_overlap}."
        )

    res = send_request("get_type", args=[obj_name])
    if res.get("status") != "success" and "not found" not in str(res.get("error", "")):
        # The plugin is not listening, or is still busy with an earlier
        # request: say that, rather than that the object does not exist, which
        # invites re-loading a structure that is already there.
        return f"Error: could not check '{obj_name}' in PyMOL: {res.get('error')}"
    kind = str(res.get("result")) if res.get("status") == "success" else ""
    if kind != "object:molecule":
        what = f"a {kind.removeprefix('object:')}" if kind else "not an object"
        return (
            f"Error: '{obj_name}' is {what}. print_ribbon_view styles a whole molecular "
            f"object; to style part of one, copy it out first (create <name>, <selection>)."
        )

    spine = f"{obj_name}_spine"
    notes: list[str] = []
    style = {
        **_RIBBON_STYLE,
        "cartoon_rect_length": str(strand_width),
        "cartoon_rect_width": str(strand_thickness),
    }

    def call(action: str, *args: object) -> Any:
        return _ribbon_call(action, *args, notes=notes)

    try:
        # Main object: chunky β-arrows + fat helix.
        call("hide", "everything", obj_name)
        call("show", "cartoon", obj_name)
        call("cartoon", "automatic", obj_name)
        # Loops are represented by the spine only - drop their fat cartoon.
        call("hide", "cartoon", f"({obj_name}) and not (ss H or ss S)")

        # Spine: one continuous tube through the entire backbone (no SS gaps).
        call("delete", spine)
        call("create", spine, obj_name)
        call("hide", "everything", spine)
        call("show", "cartoon", spine)
        call("cartoon", "tube", spine)
        # The spine takes the style too (smooth loops, flat sheets): create
        # copies per-atom settings but not the object's own.
        for target in (obj_name, spine):
            for setting, value in style.items():
                call("set", setting, value, target)
        call("set", "cartoon_tube_radius", str(spine_radius), spine)

        links = _caless_links(spine, notes)
        bridgeable = [link for link in links if not link.why_not]
        bridged: list[str] = []
        # Checked whether or not a bridge was asked for, so the result never
        # promises one that cannot be made.
        traced = _traced_already(spine, notes) if bridgeable else 0
        if bridge_modified and bridgeable and not traced:
            # Trace mode joins the CA-like atoms in atom order, whatever
            # their names, numbering or HETATM flag; on the spine only.
            trace = [i for link in bridgeable for i in link.through + link.after]
            call("set", "cartoon_trace_atoms", "1", _by_index(spine, trace))
            call("set", "trace_atoms_mode", "0", spine)
            # Loop smoothing averages each Cα with its neighbours, and the
            # chromophore's CA-like atoms sit close together, so it pulled
            # the tube 1.6 Å off Leu64's Cα (2.2 Å on a 0.3 Å spine) and
            # print_export flagged a gap there. no_smooth holds the Cα on
            # either side in place; the rest of the spine is smoothed as
            # before.
            pins = [i for link in bridgeable for i in link.before + link.after]
            call("flag", "no_smooth", _by_index(spine, pins), "set")
            # A flag does not invalidate a cartoon PyMOL has already drawn.
            call("rebuild", spine)
            bridged = [link.label for link in bridgeable]

        interiors: list[list[int]] = []
        if not spine_in_strands:
            breaks = {i for link in links for i in link.after}
            interiors = _strand_interiors(spine, strand_overlap, breaks, notes)
            hidden = [i for run in interiors for i in run]
            if hidden:
                call("hide", "cartoon", f"byres ({_by_index(spine, hidden)})")

        call("util.chainbow", obj_name)
        call("util.chainbow", spine)
        call("orient", obj_name)
    except RuntimeError as e:
        # Not rolled back, but every run starts from scratch (hide everything,
        # recreate the spine, set every setting), so a re-run repairs it.
        return (
            f"Error: {e}\n{obj_name} may be left partly restyled. print_ribbon_view "
            f"starts from scratch each time, so re-running it once the cause is fixed "
            f"puts it right."
        )

    wide, thick = _STRAND_DRAWN * strand_width, _STRAND_DRAWN * strand_thickness
    across = 2 * spine_radius
    lines = [
        f"Print-ribbon view applied to {obj_name} (spine object '{spine}' created). Style "
        f"set on {obj_name} and {spine} only; other objects and the global settings are "
        f"unchanged.",
        f"β-strands: cartoon_rect_length {style['cartoon_rect_length']} x cartoon_rect_width "
        f"{style['cartoon_rect_width']}, which PyMOL draws √2 larger: about {wide:.1f} Å "
        f"wide x {thick:.1f} Å thick. Loops are the spine only.",
    ]
    spine_line = f"Spine: a {spine_radius:g} Å-radius tube ({across:.1f} Å across)"
    if spine_in_strands:
        spine_line += " through the whole backbone, strands included."
        if across > thick:
            spine_line += (
                f" It is thicker than the {thick:.1f} Å strands and will show as bumps "
                f"along them; spine_in_strands=False takes it out of strand interiors."
            )
    elif interiors:
        strands = "1 strand" if len(interiors) == 1 else f"{len(interiors)} strands"
        spine_line += (
            f", taken out of the interior of {strands} ({strand_overlap} residues kept at "
            f"each end so the loops still fuse to the strands)."
        )
    else:
        spine_line += (
            f" through the whole backbone: no strand is longer than {2 * strand_overlap} "
            f"residues, so spine_in_strands=False had nothing to take it out of."
        )
    lines.append(spine_line)
    labels = ", ".join(link.label for link in bridgeable)
    if bridged:
        residues = "1 residue" if len(bridged) == 1 else f"{len(bridged)} residues"
        lines.append(
            f"Spine bridged through {residues} with no CA atom: {labels} (traced through "
            f"the Cα-like atoms with trace_atoms_mode 0, the Cα either side held against "
            f"loop smoothing, on {spine} only; nothing renamed)."
        )
    elif bridgeable and not bridge_modified and not traced:
        lines.append(f"Spine breaks at {labels} (no CA atom); bridge_modified=True joins it.")
    elif bridgeable:
        lines.append(
            f"Spine breaks at {labels} (no CA atom) and was not bridged: {spine} already "
            f"traces {traced} atom(s) (cartoon_trace_atoms, or Cα with no bonds, which PyMOL "
            f"traces itself), and the trace_atoms_mode 0 a bridge needs would join those "
            f"across chain breaks."
        )
    for link in links:
        if link.why_not:
            lines.append(f"Spine breaks at {link.label} (no CA atom): not bridged, {link.why_not}.")
    if notes:
        lines.append("PyMOL warned: " + "; ".join(dict.fromkeys(notes)))

    # Fine enough for print_export's own thin-feature check, 3 voxels across
    # the thinnest thing drawn, compared the way it compares them: in binary,
    # a 0.6 Å spine is under 3 x 0.2.
    thinnest = min(thick, across, 2 * float(_RIBBON_STYLE["cartoon_oval_width"]))
    pitch = min(_RIBBON_VOXEL_PITCH, math.floor(thinnest / _THIN_VOXELS * 100) / 100)
    if thinnest < _THIN_VOXELS * pitch:
        pitch = round(pitch - 0.01, 2)
    lines.append(
        f"Export the two fused into one solid (chains that do not touch come out as "
        f"separate pieces, which print_export reports) with:\n"
        f'  print_export(obj_name="{obj_name}", groups="{obj_name}=({obj_name} or {spine})", '
        f'representation="cartoon", method="voxel", voxel_pitch={pitch:g})'
    )
    return "\n".join(lines)


# ── residue coverage ─────────────────────────────────────────────────────────
#
# A shown residue whose Cα ends up outside the printed mesh is a gap nothing
# reported: in the split-GFP build, Leu64's Cα sat 1.84 Å outside the white
# tube because loop smoothing bends the tube off the Cα next to a hidden
# residue, so the part meant to attach there would not have touched. Normal
# smoothed residues sit up to 1.9 Å off too, so a residue whose sequence
# neighbour is hidden (a junction) gets a tighter tolerance.
_JUNCTION_TOLERANCE = 1.0  # Å
_SMOOTHED_TOLERANCE = 2.0  # Å, elsewhere, while cartoon_smooth_loops or flat_sheets is on
_COVERAGE_LISTED = 10
# Indexing costs ~2 s per million faces and a 20-million-face mesh took
# minutes; past this the check is skipped, and the result says so.
_COVERAGE_MAX_FACES = 4_000_000


class _SurfaceIndex:
    """Signed distance from points to a closed mesh (positive outside).

    trimesh's proximity queries and ``contains`` need rtree, which the print
    extra does not install, so this uses scipy: a cKDTree on the centroids of
    the triangles, subdivided until no edge is longer than ``max_edge``, finds
    candidates; trimesh.triangles.closest_point picks the nearest; the sign
    comes from the vertex normals interpolated at that point.
    """

    def __init__(self, mesh: "trimesh.Trimesh", max_edge: float) -> None:
        import numpy as np
        import trimesh
        from scipy.spatial import cKDTree

        tri = np.asarray(mesh.vertices, dtype=float)[np.asarray(mesh.faces)]
        corner_normals = np.asarray(mesh.vertex_normals, dtype=float)[np.asarray(mesh.faces)]
        area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
        tri, corner_normals = tri[area > 1e-12], corner_normals[area > 1e-12]
        flat = tri.reshape(-1, 3)
        vertices, faces, parent = trimesh.remesh.subdivide_to_size(
            flat,
            np.arange(len(flat)).reshape(-1, 3),
            max_edge=max_edge,
            max_iter=12,
            return_index=True,
        )
        self.triangles = vertices[faces]
        bary = trimesh.triangles.points_to_barycentric(
            np.repeat(tri[parent], 3, axis=0), self.triangles.reshape(-1, 3)
        )
        self.normals = np.einsum(
            "ij,ijk->ik", bary, np.repeat(corner_normals[parent], 3, axis=0)
        ).reshape(-1, 3, 3)
        centres = self.triangles.mean(axis=1)
        self.reach = float(np.linalg.norm(self.triangles - centres[:, None, :], axis=2).max())
        self.tree = cKDTree(centres)

    def signed_distance(self, points: "np.ndarray", k: int = 16) -> "np.ndarray":
        import numpy as np
        import trimesh

        points = np.asarray(points, dtype=float).reshape(-1, 3)
        k = min(k, len(self.triangles))
        near, idx = self.tree.query(points, k=k)
        near, idx = near.reshape(len(points), -1), idx.reshape(len(points), -1)
        closest = trimesh.triangles.closest_point(
            self.triangles[idx.ravel()], np.repeat(points, idx.shape[1], axis=0)
        ).reshape(len(points), -1, 3)
        dist = np.linalg.norm(closest - points[:, None, :], axis=2)
        pick = dist.argmin(axis=1)
        rows = np.arange(len(points))
        best, face, where = dist[rows, pick], idx[rows, pick], closest[rows, pick]
        # A triangle whose centroid is further than the k-th one can still be
        # nearer, by up to its reach: check those points exhaustively.
        for i in np.nonzero(near[:, -1] - self.reach < best)[0]:
            cand = np.asarray(self.tree.query_ball_point(points[i], best[i] + self.reach))
            if not len(cand):
                continue
            pts = trimesh.triangles.closest_point(
                self.triangles[cand], np.repeat(points[i][None], len(cand), axis=0)
            )
            d = np.linalg.norm(pts - points[i], axis=1)
            j = int(d.argmin())
            if d[j] < best[i]:
                best[i], face[i], where[i] = d[j], cand[j], pts[j]
        bary = trimesh.triangles.points_to_barycentric(self.triangles[face], where)
        normal = np.einsum("ij,ijk->ik", bary, self.normals[face])
        sign = np.sign(np.einsum("ij,ij->i", points - where, normal))
        sign[sign == 0] = 1.0
        return sign * best


def _junction_selection(obj: str) -> str:
    """Shown guide atoms whose backbone is bonded to a residue drawn without cartoon.

    That covers a hidden neighbour and one with no guide atom at all (GFP's
    chromophore has CA1-3 but no CA).
    """
    shown = f"(byres (({obj}) and guide and rep cartoon))"
    return (
        f"(({obj}) and guide and rep cartoon and present) and byres "
        f"((({obj}) and name C+N) and neighbor (({obj}) and not {shown}))"
    )


def _coverage_warnings(objs: "list[str]", stl_path: str, mm_per_angstrom: float) -> "list[str]":
    """WARNING lines for shown residues whose Cα is not inside the final STL."""
    rows: list[tuple[str, str, str, str]] = []
    coords: list = []
    junction: list[bool] = []
    smooth: list[bool] = []
    for obj in objs:
        # "present": PyMOL 3 loads unobserved residues as atoms without
        # coordinates, which iterate lists and get_coords skips; without it
        # 2b3p gave 243 rows and 228 coordinates, and every label slipped.
        sel = f"({obj}) and guide and rep cartoon and present"
        found = _result_list(send_request("iterate_to_list", args=[sel, "chain, resi, resn"]))
        if not found:
            continue
        xyz = _result_list(send_request("get_coords", args=[sel]))
        if len(xyz) != len(found):
            continue
        joins = {
            (str(chain), str(resi))
            for chain, resi in _result_list(
                send_request("iterate_to_list", args=[_junction_selection(obj), "chain, resi"])
            )
        }
        smoothed = _setting_on("cartoon_smooth_loops", obj) or _setting_on(
            "cartoon_flat_sheets", obj
        )
        for (chain, resi, resn), point in zip(found, xyz, strict=True):
            rows.append((obj, str(chain), str(resi), str(resn)))
            coords.append(point)
            junction.append((str(chain), str(resi)) in joins)
            smooth.append(smoothed)
    if not rows:
        return []

    import numpy as np
    import trimesh

    mesh = cast("trimesh.Trimesh", trimesh.load(stl_path, force="mesh"))
    index = _SurfaceIndex(mesh, max_edge=0.8 * mm_per_angstrom)
    distance = index.signed_distance(np.asarray(coords, dtype=float) * mm_per_angstrom)
    distance = distance / mm_per_angstrom
    tolerance = np.where(
        junction, _JUNCTION_TOLERANCE, np.where(smooth, _SMOOTHED_TOLERANCE, _JUNCTION_TOLERANCE)
    )
    out = distance > tolerance
    warnings = []

    def _name(i: int) -> str:
        obj, chain, resi, resn = rows[i]
        return f"{obj} {chain}/{resn}{resi}"

    at_junction = [int(i) for i in np.nonzero(out)[0] if junction[i]]
    if at_junction:
        listed = ", ".join(
            f"{_name(i)} ({distance[i]:.2f} Å out)" for i in at_junction[:_COVERAGE_LISTED]
        )
        more = (
            f" and {len(at_junction) - _COVERAGE_LISTED} more"
            if len(at_junction) > _COVERAGE_LISTED
            else ""
        )
        objs_hit = sorted({rows[i][0] for i in at_junction})
        warnings.append(
            f"WARNING: {len(at_junction)} shown residue(s) next to a hidden one have their "
            f"Cα outside the mesh: {listed}{more}. A part meant to attach there will not "
            f"touch. Loop smoothing bends the tube off the Cα where a neighbour is hidden: "
            + "; ".join(f"`set cartoon_smooth_loops, 0, {obj}`" for obj in objs_hit)
            + ", or show the neighbour."
        )
    elsewhere = [int(i) for i in np.nonzero(out)[0] if not junction[i]]
    if elsewhere:
        chains = Counter(f"{rows[i][0]} chain {rows[i][1]}" for i in elsewhere)
        listed = ", ".join(_name(i) for i in elsewhere[:_COVERAGE_LISTED])
        more = (
            f" and {len(elsewhere) - _COVERAGE_LISTED} more"
            if len(elsewhere) > _COVERAGE_LISTED
            else ""
        )
        warnings.append(
            f"WARNING: {len(elsewhere)} residue(s) PyMOL shows as cartoon are missing from "
            f"the mesh, up to {max(distance[i] for i in elsewhere):.1f} Å away "
            f"({', '.join(f'{c}: {n}' for c, n in chains.most_common())}): {listed}{more}. "
            f"A dropped piece (bodies='largest'), or a feature too thin for the voxel pitch?"
        )
    return warnings


@mcp.tool()
def print_export(
    obj_name: Annotated[
        str,
        Field(
            description='PyMOL object the groups are taken from (e.g. "1abc"); the STL files are named after it.'
        ),
    ],
    groups: Annotated[
        str,
        Field(
            description='Semicolon-separated ``label=selection`` pairs, one per colour. Example: "protein=polymer.protein; nucleic=polymer.nucleic". In cartoon and displayed modes each group names whole objects, e.g. "barrel=(barrel or barrel_spine)".'
        ),
    ],
    out_dir: Annotated[
        str, Field(description="Directory for the STL files (default: current directory).")
    ] = ".",
    method: Annotated[
        str,
        Field(
            description='"auto" (light cleanup if the export is already watertight apart from stray fragments, else poisson with voxel fallback), "poisson" (keeps detail, best for bulky surfaces; cannot merge sticks and spheres), or "voxel" (robust for thin cartoons, tubes and nucleic acids).'
        ),
    ] = "auto",
    voxel_pitch: Annotated[
        float,
        Field(
            description="Voxel size in Ångström for the voxel method; it stays in Å whatever mm_per_angstrom is. Use about one third of the thinnest feature you need to keep: 0.7 keeps a ~10 Å helix, print_ribbon_view ribbons want 0.2, thin rods 0.1. The result warns when a displayed cartoon feature is under 3 voxels across. Halving it costs about 5x the time."
        ),
    ] = 0.7,
    poisson_depth: Annotated[
        int, Field(description="Screened-Poisson octree depth; higher = more detail.")
    ] = 10,
    representation: Annotated[
        str,
        Field(
            description='What geometry to export. "surface" (default) isolates each group as its own temporary object and exports its molecular surface. "cartoon" exports the cartoon and surface triangles the group\'s objects *display*, keeping per-residue rep flags and per-object cartoon type (e.g. a `cartoon tube` spine from print_ribbon_view); it leaves out sticks, spheres, lines and dashes and names what it left out. "displayed" exports everything PyMOL draws for the group\'s objects: those triangles plus sticks, spheres, lines, ribbon and distance dashes (as solid rods), merged exactly. In cartoon and displayed modes groups are isolated by object visibility, so each group names whole object(s), one colour per object (e.g. "1abc or 1abc_spine"); name measurement and CGO objects in the group too.'
        ),
    ] = "surface",
    bodies: Annotated[
        str,
        Field(
            description='"all" (default) keeps every separate piece of a group, e.g. both dimers of 4hhb\'s ribbon or all four hemes, and the result says it prints as N pieces. "largest" keeps only the biggest piece and reports what it dropped, with a WARNING when that is over 1% of the volume. Either way enclosed cavities are filled and stray fragments under 8 faces removed.'
        ),
    ] = "all",
    undo_voxel_growth: Annotated[
        bool,
        Field(
            description="Voxel method only: move every surface back in by 0.62 x voxel_pitch, the amount voxelising grows it (a 0.33 Å rod came out 0.38 Å at 0.1 Å pitch). Default true. Gaps narrower than about 2.5 x voxel_pitch fuse while voxelising and cannot be reopened, so where a clearance matters use voxel_pitch <= gap / 2.5."
        ),
    ] = True,
    mm_per_angstrom: Annotated[
        float,
        Field(
            description="Print scale in millimetres per Ångström. 1.0 (default) prints 1 Å as 1 mm; 3.0 makes a 40 Å protein 120 mm across. Applied after repair to every group about the same origin, so the parts stay aligned; voxel_pitch stays in Å. Raise it when the result warns that a part is small enough for a slicer to mistake it for inches."
        ),
    ] = 1.0,
) -> str:
    """
    Exports a structure as watertight STL files ready for multi-colour 3D printing.

    Each colour group becomes one STL file. PyMOL's exporters write the whole
    visible scene, so each group is isolated on its own before export, then
    rebuilt into a watertight solid. Every STL is written in PyMOL model
    coordinates (Å, the numbers get_coords reports) whatever the camera view,
    so groups from one call, or from separate calls, line up; the camera is
    left as it was. A slicer can load them as aligned multi-material parts
    (e.g. add the second STL as a "part" of the first in Bambu Studio).

    representation="displayed" prints what PyMOL shows, sticks and spheres
    included; "cartoon" keeps to cartoon and surface triangles and names what
    it left out. A group keeps every separate piece (two chains that do not
    touch print as two pieces) unless bodies="largest". The voxel method's
    growth of about 0.62 voxel per surface is undone unless
    undo_voxel_growth=False, but gaps narrower than about 2.5 voxels fuse
    while voxelising and cannot be reopened: a clearance needs
    voxel_pitch <= gap / 2.5.

    STL has no units and slicers read it as millimetres, so ``mm_per_angstrom``
    sets the print size. The result gives each part's size in mm and warns about
    any part small enough for a slicer to take for inches and enlarge 25.4x.

    Requires the optional ``print`` extra (trimesh, pymeshlab, manifold3d); see
    the install hint returned if the libraries are missing.
    """
    try:
        import trimesh  # noqa: F401
    except ImportError:
        return _PRINT_DEPS_HINT

    try:
        pairs = _parse_groups(groups)
    except ValueError as e:
        return f"Error: {e}"

    representation = representation.lower()
    if representation not in _REPRESENTATIONS:
        return (
            f"Error: representation must be 'surface', 'cartoon' or 'displayed', "
            f"got {representation!r}"
        )
    method = method.lower()
    if method not in _METHODS:
        return f"Error: method must be 'auto', 'poisson' or 'voxel', got {method!r}"
    bodies = bodies.lower()
    if bodies not in _BODY_POLICIES:
        return f"Error: bodies must be 'all' or 'largest', got {bodies!r}"
    if not (math.isfinite(voxel_pitch) and voxel_pitch > 0):
        return f"Error: voxel_pitch must be a positive number of Å, got {voxel_pitch}."
    # Zero collapses the model and a negative scale mirrors it, which swaps
    # every stereocentre: an L-protein would print as its D enantiomer.
    if not (math.isfinite(mm_per_angstrom) and mm_per_angstrom > 0):
        return f"Error: mm_per_angstrom must be a positive number, got {mm_per_angstrom}."
    displayed = representation == "displayed"
    if displayed:
        try:
            import manifold3d  # noqa: F401
        except ImportError:
            return _PRINT_DEPS_HINT

    try:
        os.makedirs(out_dir, exist_ok=True)
    except OSError as e:
        return f"Error: could not create the output directory {out_dir!r}: {e}"
    # Intermediates go to a directory of their own, removed afterwards: they
    # used to land beside the STLs, in another frame and unscaled, so a folder
    # dragged into a slicer brought them along.
    workdir = tempfile.mkdtemp(prefix="mcpymol-print-")
    results: list[tuple[str, str, str, dict, list[str]]] = []

    def _stopped(message: str) -> str:
        """``message``, plus the groups already written when a later one fails."""
        if not results:
            return message
        written = ", ".join(os.path.basename(r[2]) for r in results)
        return (
            f"{message}\nThe export stopped there. Already written to {out_dir} before it "
            f"did: {written}. Each is complete for its own group, but the set is not; "
            f"re-run once the cause is fixed rather than printing these alone."
        )

    try:
        for i, (label, sel) in enumerate(pairs):
            stl_path = os.path.join(out_dir, f"{obj_name}_{label}.stl")
            obj_path = os.path.join(workdir, f"group{i}.obj")
            pov_path = os.path.join(workdir, f"group{i}.pov")
            objs: list[str] = []
            others: list[tuple[str, str]] = []
            census: dict[str, dict[str, int]] = {}
            primitives = None
            thin = None
            extra: list[str] = []
            # The session is disturbed only while this group is written out:
            # visibility, rock, stereo and the temporary object are put back
            # before the repair, which can take minutes, so a server that dies
            # mid-repair leaves the session as the user had it, and a change
            # the user makes meanwhile is not undone afterwards.
            try:
                enabled = _enabled_names()
            except RuntimeError as e:
                return _stopped(f"Error: {e}")
            paused: dict[str, bool] = {}
            tmp = None
            try:
                paused = _pause_motion()
                if representation == "surface":
                    # Named by position, never by label: a label with a space
                    # made `delete` read a list of names and delete the user's
                    # objects ('chain A' deleted an object called A).
                    tmp = f"_mcpymol_print_{i}"
                    send_request("delete", args=[tmp])
                    # zoom=0: by default create re-zooms and re-centres the
                    # user's camera.
                    res = send_request(
                        "create", args=[tmp, f"({obj_name}) and ({sel})"], kwargs={"zoom": 0}
                    )
                    if res.get("status") == "error":
                        return _stopped(f"Error isolating group '{label}': {res.get('error')}")
                    if not _numeric_result(send_request("count_atoms", args=[tmp])):
                        return _stopped(
                            f"Error: group '{label}' selection '{sel}' matched no atoms "
                            f"in '{obj_name}'."
                        )
                    # Isolate: only this temp object visible, only as a surface.
                    send_request("hide", args=["everything", tmp])
                    send_request("show", args=["surface", tmp])
                    _show_only([tmp])
                    _save_in_model_frame([obj_path])
                else:
                    # Export what the real objects display, so per-residue rep
                    # flags (hidden loops) and per-object cartoon type (a
                    # `cartoon tube` spine) are kept exactly. PyMOL's exporters
                    # write the whole visible scene, so isolate a group by
                    # enabling only its objects. One colour per object.
                    objs, others = _resolve_group(sel)
                    if not objs and not others:
                        return _stopped(
                            f"Error: {representation}-mode group '{label}' selection "
                            f"'{sel}' matched no objects. Each group must name whole "
                            f"object(s), e.g. '1abc or 1abc_spine'."
                        )
                    census = _census(objs)
                    if displayed:
                        pellets = [
                            f"{obj} ({census[obj]['dots'] + census[obj]['mesh']:,} atoms)"
                            for obj in objs
                            if census[obj].get("dots") or census[obj].get("mesh")
                        ]
                        if pellets:
                            # Refused before anything is written: 1ubq with dots
                            # shown was 38,567 spheres, ten minutes and a 1 GB
                            # STL of loose pellets.
                            return _stopped(
                                f"Error: group '{label}' shows dots or mesh on "
                                f"{', '.join(pellets)}. PyMOL draws those as a tiny sphere "
                                f"or rod per point, which would print as a cloud of loose "
                                f"pellets. Hide them (hide dots; hide mesh), or use "
                                f"representation='cartoon' or 'surface'."
                            )
                    _show_only(objs + [name for name, _ in others])
                    thin = _thinnest_cartoon(objs)
                    measurements = [name for name, kind in others if kind == "measurement"]
                    saved = []
                    try:
                        if displayed and measurements:
                            # One solid rod per distance: dashes would print as
                            # loose pellets.
                            saved = [_solid_dashes(name) for name in measurements]
                            send_request("refresh")
                        pov_shift = _save_in_model_frame(
                            [obj_path, pov_path] if displayed else [obj_path]
                        )
                    finally:
                        for item in saved:
                            _restore_dashes(item)
                    if saved:
                        extra.append(
                            f"NOTE: the distance dashes of {', '.join(measurements)} were "
                            f"exported as solid rods (dash_gap 0 during the export only)."
                        )
                    if displayed:
                        primitives = _parse_pov(pov_path, pov_shift)
                triangles, markers = _read_obj_triangles(obj_path)
            except (RuntimeError, OSError) as e:
                return _stopped(f"Error exporting group '{label}': {e}")
            finally:
                if tmp is not None:
                    send_request("delete", args=[tmp])
                _restore_enabled(enabled)
                _resume_motion(paused)

            n_prims = len(primitives["cylinders"]) + len(primitives["spheres"]) if primitives else 0
            if not len(triangles) and not n_prims:
                if representation == "surface":
                    return _stopped(
                        f"Error: group '{label}' ({sel}) exported no surface triangles."
                    )
                return _stopped(
                    _describe_empty_group(label, sel, objs, others, census, representation)
                )
            if n_prims > _MAX_PRIMITIVES:
                return _stopped(
                    f"Error: group '{label}' draws {n_prims:,} spheres and rods, more than "
                    f"print_export merges ({_MAX_PRIMITIVES:,}); dots and the mesh "
                    f"representation write one per point. Hide them, or use "
                    f"representation='cartoon' or 'surface'."
                )
            if n_prims and method == "poisson":
                return _stopped(
                    f"Error: group '{label}' shows sticks or spheres, and Poisson "
                    f"reconstruction cannot merge them (it lost 78% of a ligand's "
                    f"sticks). Use method='auto' or 'voxel'."
                )
            src: str | None = obj_path
            if not len(triangles):
                src = None
            elif markers:
                # The sphere markers would voxelise into specks; the spheres
                # themselves come from the POV-Ray export.
                src = os.path.join(workdir, f"group{i}_triangles.obj")
                _write_obj(src, triangles)

            if representation == "cartoon":
                left_out = [
                    f"{obj}: {_shown_reps(census[obj], _PRIMITIVE_REPS)}"
                    for obj in objs
                    if _shown_reps(census[obj], _PRIMITIVE_REPS)
                ]
                left_out += [f"{name} ({kind} object)" for name, kind in others]
                if left_out:
                    extra.append(
                        f"NOTE: left out {'; '.join(left_out)}. PyMOL draws these as "
                        f"spheres and rods, which cartoon mode does not export; "
                        f"representation='displayed' includes them."
                    )
            if objs:
                atoms = _atoms_only(sel, [name for name, _ in others])
                covered = _numeric_result(send_request("count_atoms", args=[f"({atoms})"]))
                total = _numeric_result(send_request("count_atoms", args=[" or ".join(objs)]))
                if 0 < covered < total:
                    extra.append(
                        f"NOTE: {representation} mode exports whole objects. '{sel}' "
                        f"covers {covered:,.0f} of the {total:,.0f} atoms of "
                        f"{', '.join(objs)}, but everything they show went into this "
                        f"STL. To print part of an object, copy it out first "
                        f"(create <name>, <selection>)."
                    )
            if primitives and primitives["unknown"]:
                kinds = ", ".join(f"{n} '{k}'" for k, n in primitives["unknown"].items())
                extra.append(
                    f"WARNING: PyMOL's POV-Ray export held {kinds} primitive(s) that "
                    f"print_export cannot read; they are missing from the STL."
                )

            try:
                info = _repair_to_stl(
                    src,
                    stl_path,
                    method,
                    voxel_pitch,
                    poisson_depth,
                    bodies=bodies,
                    undo_voxel_growth=undo_voxel_growth,
                    mm_per_angstrom=mm_per_angstrom,
                    primitives=primitives if n_prims else None,
                )
            except ImportError:
                return _stopped(_PRINT_DEPS_HINT)
            except Exception as e:
                return _stopped(f"Error rebuilding group '{label}': {type(e).__name__}: {e}")

            if "voxel" in info["method"]:
                warning = _thin_cartoon_warning(thin, voxel_pitch)
                if warning:
                    extra.append(warning)
            if primitives is not None and n_prims:
                valence = any(
                    census[obj].get("sticks") and _setting_on("valence", obj) for obj in objs
                )
                warning = _thin_rod_warning(
                    primitives,
                    census,
                    valence,
                    any(kind == "measurement" for _, kind in others),
                    mm_per_angstrom,
                )
                if warning:
                    extra.append(warning)
            if objs and info["faces"] > _COVERAGE_MAX_FACES:
                extra.append(
                    f"NOTE: the residue coverage check was skipped: the mesh has "
                    f"{info['faces']:,} faces, over {_COVERAGE_MAX_FACES:,}."
                )
            elif objs:
                try:
                    extra.extend(_coverage_warnings(objs, stl_path, mm_per_angstrom))
                except Exception as e:
                    extra.append(
                        f"NOTE: the residue coverage check could not run "
                        f"({type(e).__name__}: {e}); the STL itself is written."
                    )
            warning = _unit_guess_warning(label, info["volume"], mm_per_angstrom)
            if warning:
                extra.append(warning)
            results.append((label, sel, stl_path, info, extra))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    s = mm_per_angstrom
    lines = [f"Exported {len(results)} group(s) for 3D printing to {out_dir}, at {s:g} mm per Å:"]
    for label, sel, stl_path, info, extra in results:
        flag = "watertight" if info["watertight"] else "NOT watertight - check the mesh"
        lines.append(
            f"  {label} ({sel}) -> {os.path.basename(stl_path)}  "
            f"[{info['method']}, {info['faces']:,} faces, {flag}]"
        )
        size = _describe_size(info["bounds"], info["volume"])
        lines.append(f"      {size}, one piece" if len(info["bodies"]) == 1 else f"      {size}")
        lines.extend(f"      {line}" for line in _describe_bodies(info) + extra)
    if len(results) > 1:
        lo = [min(r[3]["bounds"][0][k] for r in results) for k in range(3)]
        hi = [max(r[3]["bounds"][1][k] for r in results) for k in range(3)]
        lines.append(f"Assembled: {_format_dims(lo, hi)} mm.")
    lines.append(
        f"Coordinates: PyMOL model coordinates (what get_coords reports) whatever the "
        f"camera view, scaled by {s:g} about the origin (0, 0, 0), so STLs from this "
        f"and any other print_export call at the same scale line up. Your view was "
        f"left as it was."
    )
    if any("voxel" in r[3]["method"] for r in results):
        gap = _GAP_FUSE_VOXELS * voxel_pitch
        lines.append(
            f"Voxel pitch {voxel_pitch:g} Å is {voxel_pitch * s:.2g} mm in the print. Gaps "
            f"narrower than about {gap:.2g} Å ({gap * s:.2g} mm) fused while voxelising "
            f"and cannot be reopened: for a clearance, use voxel_pitch <= gap / "
            f"{_GAP_FUSE_VOXELS:g}."
        )
    lines.append(
        "In your slicer, load the first STL, then add the others as parts (don't "
        "re-centre) and assign a filament per part."
    )
    return "\n".join(lines)
