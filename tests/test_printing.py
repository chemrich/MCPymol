"""3D-print tools, against a mocked PyMOL.

These pin what the tools send. Whether PyMOL then draws the right geometry is
for tests/test_live_printing.py, against a real one. The print_export tests
also run the mesh code itself on small synthetic meshes where trimesh and
manifold3d are installed; CI's test job has neither, and those tests skip.
"""

import inspect
import math
import os
import re
import sys
from unittest.mock import MagicMock, patch

import pytest

from mcpymol.printing import _REP_BITS
from mcpymol.server import (
    _atoms_only,
    _coverage_warnings,
    _describe_bodies,
    _inset,
    _merge_continuing_rods,
    _parse_pov,
    _read_obj_triangles,
    _repair_to_stl,
    _run_poisson,
    _select_bodies,
    _strand_interiors,
    _SurfaceIndex,
    _thin_cartoon_warning,
    _thin_rod_warning,
    _thinnest_cartoon,
    _union_primitives,
    _unit_guess_warning,
    print_export,
    print_ribbon_view,
)

# ── print_export ─────────────────────────────────────────────────────────────
#
# A fake PyMOL that keeps the state print_export reads back: a camera that
# remembers set_view, the set of enabled objects, and the settings it changes.
# Mesh repair is patched out where a test is about what reaches PyMOL; the
# tests that need real meshes import trimesh/manifold3d or skip (CI's test job
# runs without the print extra).

# A rotated, shifted user view; print_export must not care what it is.
_USER_VIEW = [
    *(0.5, -0.5, 0.7071, 0.7071, 0.7071, 0.0, -0.5, 0.5, 0.7071),
    *(5.0, -3.0, -110.0),
    *(30.0, 28.0, 15.0),
    *(84.0, 135.0, -20.0),
]
_MODEL_FRAME = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


class _PyMOL:
    """A send_request stand-in that answers the way PyMOL does.

    ``objects`` maps object name to type ("molecule", "measurement", ...).
    ``counts`` maps a selection substring to count_atoms' answer (first
    match wins). ``drift`` turns the camera a little on every save, as rock
    or a playing movie would.
    """

    def __init__(
        self,
        objects=None,
        enabled=None,
        counts=None,
        rows=None,
        settings=None,
        state=None,
        save_error=None,
        drift=0.0,
        object_list=None,
        warnings=None,
    ):
        self.objects = dict(objects or {"ub": "molecule"})
        self.enabled = list(self.objects if enabled is None else enabled)
        self.counts = dict(counts or {})
        self.rows = dict(rows or {})
        self.settings = {"dash_gap": 0.45, "geometry_export_mode": 0.0, **(settings or {})}
        self.state = {"stereo": "off", "rock": "off", **(state or {})}
        self.view = list(_USER_VIEW)
        self.save_error = save_error
        self.drift = drift
        self.object_list = object_list
        self.warnings = warnings
        self.calls = []
        self.saved = []

    def __call__(self, action, args=None, kwargs=None, **_):
        args, kwargs = list(args or []), dict(kwargs or {})
        self.calls.append((action, args, kwargs))
        handler = getattr(self, f"_{action}", None)
        result = handler(args, kwargs) if handler else "OK"
        if isinstance(result, dict):
            return result
        if result is None:
            # What the plugin sends for a None return.
            result = f"Executed '{action}' successfully."
        return {"status": "success", "result": result}

    def sent(self, action):
        return [args for name, args, _ in self.calls if name == action]

    def _get_names(self, args, kwargs):
        if args[:2] == ["all", 1]:
            return list(self.enabled)
        return list(self.objects)

    def _get_type(self, args, kwargs):
        kind = self.objects.get(args[0])
        if kind is None:
            return {"status": "error", "error": f"object {args[0]} not found"}
        return f"object:{kind}"

    def _get_object_list(self, args, kwargs):
        if self.object_list is not None:
            return self.object_list
        response = {"status": "success", "result": None}
        named = [n for n, kind in self.objects.items() if kind == "molecule" and n in args[0]]
        response["result"] = named or "Executed 'get_object_list' successfully."
        if self.warnings:
            response["warnings"] = list(self.warnings)
        return response

    def _get_view(self, args, kwargs):
        return list(self.view)

    def _set_view(self, args, kwargs):
        self.view = [float(v) for v in args[0]]

    def _save(self, args, kwargs):
        if self.save_error:
            return {"status": "error", "error": self.save_error}
        self.saved.append((args, list(self.view), self.settings["geometry_export_mode"]))
        self.view[0] += self.drift

    def _get(self, args, kwargs):
        return self.state.get(args[0], "off")

    def _stereo(self, args, kwargs):
        self.state["stereo"] = args[0]

    def _rock(self, args, kwargs):
        self.state["rock"] = "on" if args[0] else "off"

    def _disable(self, args, kwargs):
        self.enabled = [] if args[0] == "all" else [n for n in self.enabled if n != args[0]]

    def _enable(self, args, kwargs):
        if args[0] not in self.enabled:
            self.enabled.append(args[0])

    def _count_atoms(self, args, kwargs):
        for key, value in self.counts.items():
            if key in args[0]:
                return value
        return 0 if " rep " in args[0] else 100

    def _iterate_to_list(self, args, kwargs):
        if args[1:2] == ["reps"]:
            # The census reads each atom's reps mask. Answer with one atom per
            # count given as "rep <name>", so `counts` drives both forms.
            masks = []
            for key, n in self.counts.items():
                m = re.search(r"\brep (\w+)", key)
                if m and m.group(1) in _REP_BITS:
                    masks += [(_REP_BITS[m.group(1)],)] * n
            return masks
        for key, value in self.rows.items():
            if key in args[0]:
                return value
        return []

    def _get_setting_float(self, args, kwargs):
        return self.settings.get(args[0], 0.0)

    def _get_setting_boolean(self, args, kwargs):
        return bool(self.settings.get(args[0], False))

    def _set(self, args, kwargs):
        if args[0] == "geometry_export_mode":
            self.settings["geometry_export_mode"] = float(args[1])


def _repaired(volume=48000.0, lo=(0.0, 0.0, 0.0), hi=(40.0, 30.0, 40.0), method="voxel"):
    """What _repair_to_stl reports for one group, sizes in mm."""
    return {
        "method": method,
        "faces": 100,
        "watertight": True,
        "bounds": (lo, hi),
        "volume": volume,
        "bodies": [volume],
        "dropped": [],
        "cavities": [],
        "debris": [],
        "inset": 0.0,
        "primitives": (0, 0),
    }


def _export(pymol, tmp_path, triangles=12, repaired=None, **kwargs):
    """print_export against ``pymol``, with mesh reading and repair patched out."""
    kwargs = {"obj_name": "ub", "groups": "ub=ub", "out_dir": str(tmp_path), **kwargs}
    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock(), "manifold3d": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * triangles, 0)),
        patch("mcpymol.printing._repair_to_stl", return_value=repaired or _repaired()) as repair,
    ):
        result = print_export(**kwargs)
    return result, repair


def _index(calls, action, start=0):
    return next(i for i in range(start, len(calls)) if calls[i][0] == action)


# ── the export frame ──


def test_print_export_saves_in_the_model_frame_and_puts_the_view_back(tmp_path):
    """The OBJ used to come out in camera space: turning the view 90 degrees
    between two exports swapped x and z in the split-GFP build."""
    pymol = _PyMOL()

    result, _ = _export(pymol, tmp_path, representation="cartoon", method="voxel")

    assert "Error" not in result, result
    (save,) = [i for i, call in enumerate(pymol.calls) if call[0] == "save"]
    # No selection: PyMOL ignores it and writes the whole visible scene.
    assert len(pymol.calls[save][1]) == 1
    view_before = [c for c in pymol.calls[:save] if c[0] == "set_view"][-1][1][0]
    assert view_before == [*_MODEL_FRAME, 84.0, 135.0, -20.0]
    view_after = next(c for c in pymol.calls[save:] if c[0] == "set_view")[1][0]
    assert view_after == _USER_VIEW
    assert pymol.view == _USER_VIEW
    assert "PyMOL model coordinates (what get_coords reports) whatever the camera view" in result
    assert "Your view was left as it was." in result


def test_print_export_writes_the_pov_with_the_camera_back_and_shifts_it_back(tmp_path):
    """Dashes, lines and ribbon are sized in pixels; with the camera at the
    origin PyMOL wrote their radius as NaN, so the POV-Ray file is written
    with the camera at the user's distance and moved back along z. Only
    geometry_export_mode 0 writes it that way, so the mode is held at 0."""
    pymol = _PyMOL(settings={"geometry_export_mode": 1.0})
    primitives = {"spheres": [], "cylinders": [[0, 0, 0, 1, 0, 0, 0.25]], "cyl_open": [True]}
    with (
        patch("mcpymol.printing._parse_pov", return_value={**primitives, "unknown": {}}) as parse,
        patch("mcpymol.printing._thin_rod_warning", return_value=None),
    ):
        result, _ = _export(pymol, tmp_path, representation="displayed")

    assert "Error" not in result, result
    (obj_save, obj_view, obj_mode), (pov_save, pov_view, pov_mode) = pymol.saved
    assert obj_save[0].endswith(".obj") and pov_save[0].endswith(".pov")
    assert obj_view[:15] == _MODEL_FRAME
    assert pov_view[:15] == [*_MODEL_FRAME[:11], -110.0, 0.0, 0.0, 0.0]
    assert (obj_mode, pov_mode) == (0.0, 0.0)
    parse.assert_called_once_with(pov_save[0], -110.0)
    # Both the user's view and their export mode come back.
    assert pymol.view == _USER_VIEW
    assert pymol.settings["geometry_export_mode"] == 1.0


def test_print_export_puts_the_view_back_when_the_save_fails(tmp_path):
    pymol = _PyMOL(save_error="disk full")

    result, repair = _export(pymol, tmp_path, representation="cartoon")

    assert result.startswith("Error exporting group 'ub'")
    assert "disk full" in result
    assert pymol.view == _USER_VIEW
    repair.assert_not_called()


def test_print_export_puts_the_view_back_when_the_repair_raises(tmp_path):
    pymol = _PyMOL(state={"rock": "on", "stereo": "on"})
    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * 12, 0)),
        patch("mcpymol.printing._repair_to_stl", side_effect=RuntimeError("mesh exploded")),
    ):
        result = print_export(obj_name="ub", groups="ub=ub", out_dir=str(tmp_path))

    assert result.startswith("Error rebuilding group 'ub': RuntimeError: mesh exploded")
    assert pymol.view == _USER_VIEW
    assert pymol.state == {"stereo": "on", "rock": "on"}


def test_print_export_refuses_a_camera_that_keeps_moving(tmp_path):
    """Rock or a playing movie turns the camera between set_view and save:
    the frame cannot be trusted, so three tries and an error, never an STL
    in some other frame."""
    pymol = _PyMOL(drift=0.01)

    result, repair = _export(pymol, tmp_path, representation="cartoon")

    assert "camera moved" in result, result
    assert len(pymol.saved) == 3
    repair.assert_not_called()
    calls = pymol.calls
    for i, call in enumerate(calls):
        if call[0] == "save":
            restore = calls[_index(calls, "set_view", i)]
            assert restore[1][0][1:] == _USER_VIEW[1:]
    assert pymol.view == _USER_VIEW


def test_print_export_pauses_stereo_and_rock_while_it_saves(tmp_path):
    """Stereo exports one eye's view (1.2 degrees off) and NaN with the camera
    at the origin; rock turns the camera every tick."""
    pymol = _PyMOL(state={"rock": "on", "stereo": "on"})

    result, _ = _export(pymol, tmp_path, groups="a=ub; b=ub", representation="cartoon")

    assert "Error" not in result, result
    acts = [(name, args) for name, args, _ in pymol.calls]
    first_save = acts.index(("save", pymol.saved[0][0]))
    last_save = len(acts) - 1 - acts[::-1].index(("save", pymol.saved[-1][0]))
    assert ("stereo", ["off"]) in acts[:first_save]
    assert ("rock", [0]) in acts[:first_save]
    assert ("rock", [1]) in acts[last_save:]
    assert ("stereo", ["on"]) in acts[last_save:]
    assert pymol.state == {"stereo": "on", "rock": "on"}


def test_surface_mode_builds_its_objects_without_moving_the_camera(tmp_path):
    """create re-zooms and re-centres the camera unless told zoom=0."""
    pymol = _PyMOL()

    result, _ = _export(pymol, tmp_path, groups="protein=polymer.protein")

    assert "Error" not in result, result
    creates = [(args, kwargs) for name, args, kwargs in pymol.calls if name == "create"]
    assert creates == [(["_mcpymol_print_0", "(ub) and (polymer.protein)"], {"zoom": 0})]
    assert pymol.view == _USER_VIEW


def test_print_export_gives_the_session_back_before_each_repair(tmp_path):
    """Visibility, rock and stereo come back as soon as a group is written out,
    not after every group's repair. A repair can take minutes, and a server
    killed meanwhile left the user with their objects hidden, rock and stereo
    off, and temporary objects in the list."""
    pymol = _PyMOL(
        objects={"ub": "molecule", "lig": "molecule", "bystander": "molecule"},
        state={"rock": "on", "stereo": "on"},
    )
    seen = []

    def repair(*args, **kwargs):
        seen.append((sorted(pymol.enabled), dict(pymol.state)))
        return _repaired()

    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock(), "manifold3d": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * 12, 0)),
        patch("mcpymol.printing._repair_to_stl", side_effect=repair),
    ):
        result = print_export(
            obj_name="ub", groups="a=ub; b=lig", out_dir=str(tmp_path), representation="cartoon"
        )

    assert "Error" not in result, result
    assert seen == [(["bystander", "lig", "ub"], {"stereo": "on", "rock": "on"})] * 2


def test_surface_mode_temp_objects_are_named_by_position_not_label(tmp_path):
    """A temp object named after its label let a label with a space make
    `delete` read a list and delete the user's objects: 'chain A' deleted an
    object called A. Every temp object is deleted again, whatever happens."""
    pymol = _PyMOL(objects={"ub": "molecule", "A": "molecule"})

    result, _ = _export(pymol, tmp_path, groups="chain_A=chain A; rest=not chain A")

    assert "Error" not in result, result
    created = [args[0] for args in pymol.sent("create")]
    assert created == ["_mcpymol_print_0", "_mcpymol_print_1"]
    deleted = [args[0] for args in pymol.sent("delete")]
    assert deleted.count("_mcpymol_print_0") == 2  # cleared first, removed after
    assert deleted.count("_mcpymol_print_1") == 2
    assert "A" not in deleted


@pytest.mark.parametrize(
    ("groups", "reason"),
    [
        ("chain A=all", "may use only letters"),
        ("a/b=all", "may use only letters"),
        ("a=chain A; a=chain B", "would write the same STL file"),
        ("A=chain A; a=chain B", "would write the same STL file"),
    ],
)
def test_print_export_refuses_a_label_before_touching_pymol(tmp_path, groups, reason):
    """Labels name the STL files: a space or slash leaked temp objects (and
    deleted the user's), and labels equal ignoring case wrote one file holding
    the last group while the result listed both."""
    pymol = _PyMOL()

    result, repair = _export(pymol, tmp_path, groups=groups)

    assert result.startswith("Error: "), result
    assert reason in result
    assert pymol.calls == []
    repair.assert_not_called()


def test_displayed_mode_refuses_dots_before_writing_anything(tmp_path):
    """1ubq with dots shown was 38,567 spheres: ten minutes and a 1 GB STL of
    loose pellets, under the old 50,000-primitive cap."""
    pymol = _PyMOL(counts={"rep dots": 38567})

    result, repair = _export(pymol, tmp_path, representation="displayed")

    assert result.startswith("Error: group 'ub' shows dots or mesh on ub (38,567 atoms)")
    assert "hide dots" in result
    assert pymol.saved == []
    repair.assert_not_called()


# ── isolating a group ──


def test_print_export_reads_a_none_object_list_as_no_objects(tmp_path):
    """A None result reaches us as the string "Executed 'get_object_list'
    successfully.", and iterating it enabled its letters one by one: PyMOL
    matches unique name prefixes, so an unrelated object was exported "OK"."""
    pymol = _PyMOL(object_list="Executed 'get_object_list' successfully.")

    result, repair = _export(pymol, tmp_path, groups="hb=d1", representation="cartoon")

    assert result.startswith("Error: cartoon-mode group 'hb' selection 'd1' matched no objects")
    assert not [args for args in pymol.sent("enable") if len(args[0]) == 1]
    repair.assert_not_called()


def test_print_export_puts_back_exactly_what_was_enabled(tmp_path):
    """`enable all` afterwards re-enabled every object the user had switched off."""
    pymol = _PyMOL(
        objects={"ub": "molecule", "off": "molecule", "grp": "group"}, enabled=["ub", "grp"]
    )

    result, _ = _export(pymol, tmp_path, representation="cartoon")

    assert "Error" not in result, result
    assert sorted(pymol.enabled) == ["grp", "ub"]
    assert "do" not in [name for name, _, _ in pymol.calls]
    last_disable = max(i for i, c in enumerate(pymol.calls) if c[:2] == ("disable", ["all"]))
    assert ["off"] not in [args for name, args, _ in pymol.calls[last_disable:] if name == "enable"]


def test_print_export_enables_a_grouped_object_with_its_group(tmp_path):
    """An object inside a disabled PyMOL group draws, and exports, nothing;
    its visibility also updates lazily, so a refresh follows."""
    pymol = _PyMOL()

    _export(pymol, tmp_path, representation="cartoon")

    names = [name for name, _, _ in pymol.calls]
    enable = names.index("enable")
    assert pymol.calls[enable][1:] == (["ub"], {"parents": 1})
    assert "refresh" in names[enable : names.index("save")]


def test_print_export_keeps_its_intermediates_out_of_the_output_folder(tmp_path):
    """The OBJ used to be left beside the STLs, in another frame and unscaled,
    so a folder dragged into a slicer brought it along."""
    pymol = _PyMOL()

    _export(pymol, tmp_path / "out", representation="cartoon")

    (path,) = [args[0] for args, _, _ in pymol.saved]
    assert not path.startswith(str(tmp_path))
    assert not os.path.exists(os.path.dirname(path))
    assert os.listdir(tmp_path / "out") == []


def test_print_export_takes_measurement_objects_out_of_the_selection(tmp_path):
    """Named in a selection, a distance object makes PyMOL reject the whole
    expression, and "(lig or d1)" exported the distance without the ligand."""
    pymol = _PyMOL(objects={"lig": "molecule", "d1": "measurement"})

    result, _ = _export(pymol, tmp_path, groups="lig=(lig or d1)", representation="cartoon")

    assert pymol.sent("get_object_list") == [["((lig or none))"]]
    assert not [args for args in pymol.sent("count_atoms") if "d1" in args[0]]
    assert {args[0] for args in pymol.sent("enable")} >= {"lig", "d1"}
    assert "left out d1 (measurement object)" in result, result


def test_atoms_only_swaps_whole_names():
    assert _atoms_only("(d10 or d1) and not d1x", ["d1"]) == "(d10 or none) and not d1x"
    assert _atoms_only("d1|lig&!d1", ["d1"]) == "none|lig&!none"


def test_print_export_passes_on_why_pymol_rejected_a_selection(tmp_path):
    pymol = _PyMOL(warnings=['Selector-Error: Invalid selection name "lgi".'])

    result, _ = _export(pymol, tmp_path, groups="x=lgi", representation="cartoon")

    assert result == (
        "Error exporting group 'x': PyMOL rejected the selection: "
        'Selector-Error: Invalid selection name "lgi".'
    )


def test_a_redraw_warning_is_not_a_rejected_selection(tmp_path):
    """PyMOL prints some warnings while it redraws, and they land in whichever
    request is running (KNOWN-ISSUES §5). Only a Selector-Error means the
    selection was refused; anything else here once failed a valid export."""
    pymol = _PyMOL(warnings=["Extrude-Warning: no cartoon data"])

    result, _ = _export(pymol, tmp_path, groups="x=lig", representation="cartoon")

    assert "rejected the selection" not in result, result


def test_print_export_cannot_create_its_output_folder(tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("")
    pymol = _PyMOL()

    result, _ = _export(pymol, blocker / "out")

    assert result.startswith("Error: could not create the output directory")
    assert pymol.calls == []


# ── what a group shows ──


def test_print_export_explains_an_empty_group(tmp_path):
    """A group with nothing printable used to raise 'zero-size array to
    reduction operation minimum which has no identity'."""
    pymol = _PyMOL(objects={"lig3": "molecule"}, counts={"rep labels": 45})

    result, repair = _export(
        pymol, tmp_path, triangles=0, groups="lig=lig3", representation="cartoon"
    )

    assert result.splitlines() == [
        "Error: group 'lig' (lig3) exported no printable geometry.",
        "  lig3: molecule; shows labels (45 atoms)",
        "  Note: labels have no 3D geometry in any PyMOL export.",
    ]
    repair.assert_not_called()


def test_print_export_says_what_cartoon_mode_left_out(tmp_path):
    """Cartoon mode keeps to triangles, so sticks and spheres vanished
    without a word; now it names them and the mode that keeps them."""
    pymol = _PyMOL(counts={"rep sticks": 20, "rep cartoon": 76})

    result, _ = _export(pymol, tmp_path, representation="cartoon")

    assert (
        "NOTE: left out ub: sticks (20 atoms). PyMOL draws these as spheres and rods, "
        "which cartoon mode does not export; representation='displayed' includes them."
    ) in result


def test_print_export_notes_a_selection_that_covers_part_of_an_object(tmp_path):
    pymol = _PyMOL(object_list=["ub"], counts={" rep ": 0, "(chain A)": 40, "ub": 100})

    result, _ = _export(pymol, tmp_path, groups="a=chain A", representation="cartoon")

    assert "NOTE: cartoon mode exports whole objects. 'chain A' covers 40 of the 100" in result


def test_print_export_makes_distance_dashes_solid_and_puts_them_back(tmp_path):
    """Dashes would print as loose pellets: dash_gap 0 makes one rod, for the
    export only, and the object goes back to following the global value."""
    pymol = _PyMOL(objects={"lig": "molecule", "d1": "measurement"})
    rod = {"spheres": [], "cylinders": [[0, 0, 0, 3, 0, 0, 0.1]], "cyl_open": [True]}
    with (
        patch("mcpymol.printing._parse_pov", return_value={**rod, "unknown": {}}),
        patch("mcpymol.printing._thin_rod_warning", return_value=None),
    ):
        result, _ = _export(pymol, tmp_path, groups="lig=(lig or d1)", representation="displayed")

    names = [(name, args) for name, args, _ in pymol.calls]
    solid = names.index(("set", ["dash_gap", 0, "d1"]))
    save = [i for i, (name, _) in enumerate(names) if name == "save"]
    assert solid < save[0]
    assert names.index(("unset", ["dash_gap", "d1"])) > save[-1]
    assert "the distance dashes of d1 were exported as solid rods" in result


def test_print_export_refuses_poisson_for_sticks(tmp_path):
    """Poisson lost 78% of a ligand's sticks; there is no point trying."""
    pymol = _PyMOL()
    rod = {"spheres": [], "cylinders": [[0, 0, 0, 3, 0, 0, 0.25]], "cyl_open": [True]}
    with patch("mcpymol.printing._parse_pov", return_value={**rod, "unknown": {}}):
        result, repair = _export(pymol, tmp_path, representation="displayed", method="poisson")

    assert result.startswith("Error: group 'ub' shows sticks or spheres, and Poisson")
    repair.assert_not_called()


def test_print_export_needs_manifold3d_for_displayed(tmp_path):
    with (
        patch("mcpymol.printing.send_request") as mock_sr,
        patch.dict(sys.modules, {"trimesh": MagicMock(), "manifold3d": None}),
    ):
        result = print_export(
            obj_name="ub", groups="ub=ub", out_dir=str(tmp_path), representation="displayed"
        )

    assert "uv sync --extra print" in result and "manifold3d" in result
    mock_sr.assert_not_called()


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"voxel_pitch": 0.0}, "voxel_pitch must be a positive number of Å"),
        ({"voxel_pitch": float("nan")}, "voxel_pitch must be a positive number of Å"),
        ({"bodies": "biggest"}, "bodies must be 'all' or 'largest'"),
        ({"method": "sculpt"}, "method must be 'auto', 'poisson' or 'voxel'"),
        ({"mm_per_angstrom": 0.0}, "mm_per_angstrom must be a positive number"),
        ({"mm_per_angstrom": -1.0}, "mm_per_angstrom must be a positive number"),
        ({"mm_per_angstrom": float("nan")}, "mm_per_angstrom must be a positive number"),
        ({"mm_per_angstrom": float("inf")}, "mm_per_angstrom must be a positive number"),
    ],
)
def test_print_export_checks_its_arguments_before_touching_pymol(tmp_path, kwargs, message):
    """A negative scale would mirror the model: an L-protein printed as its D enantiomer."""
    with (
        patch("mcpymol.printing.send_request") as mock_sr,
        patch.dict(sys.modules, {"trimesh": MagicMock()}),
    ):
        result = print_export(obj_name="ub", groups="ub=ub", out_dir=str(tmp_path), **kwargs)

    assert result.startswith("Error:") and message in result, result
    mock_sr.assert_not_called()


# ── scale and the report ──


def test_print_export_scale_defaults_to_one_mm_per_angstrom():
    """The default keeps every earlier export the same size."""
    parameters = inspect.signature(print_export).parameters
    assert parameters["mm_per_angstrom"].default == 1.0
    assert list(parameters)[-1] == "mm_per_angstrom"


def test_print_export_gives_every_group_the_same_scale(tmp_path):
    pymol = _PyMOL()

    result, repair = _export(pymol, tmp_path, groups="a=ub; b=ub", mm_per_angstrom=2.5)

    assert [c.kwargs["mm_per_angstrom"] for c in repair.call_args_list] == [2.5, 2.5]
    assert "at 2.5 mm per Å:" in result.splitlines()[0]
    assert "scaled by 2.5 about the origin (0, 0, 0)" in result


def test_print_export_reports_each_part_in_millimetres(tmp_path):
    pymol = _PyMOL()
    sizes = [
        _repaired(86290.0, (0.0, 0.0, 0.0), (103.2, 94.9, 143.3)),
        _repaired(240.0, (10.0, 20.0, 30.0), (21.7, 29.6, 56.1)),
    ]
    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * 12, 0)),
        patch("mcpymol.printing._repair_to_stl", side_effect=sizes),
    ):
        result = print_export(obj_name="x", groups="a=ub; b=ub", out_dir=str(tmp_path))

    assert "      103.2 x 94.9 x 143.3 mm, 86.29 cm³, one piece" in result
    assert "      11.7 x 9.6 x 26.1 mm, 240 mm³, one piece" in result
    assert "Assembled: 103.2 x 94.9 x 143.3 mm." in result
    assert "WARNING" not in result


def test_print_export_warns_about_a_part_a_slicer_would_read_as_inches(tmp_path):
    """Bambu Studio scales any part under 8 mm³ by 25.4, judging each file on
    its own: the split-GFP kit's Arg96 side chain was 8.9 mm³."""
    pymol = _PyMOL()
    sizes = [_repaired(3195.8), _repaired(8.9)]
    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * 12, 0)),
        patch("mcpymol.printing._repair_to_stl", side_effect=sizes),
    ):
        result = print_export(obj_name="x", groups="barrel=ub; r96=ub", out_dir=str(tmp_path))

    warnings = [line for line in result.splitlines() if "WARNING" in line]
    assert len(warnings) == 1
    assert "'r96' encloses only 8.9 mm³" in warnings[0]
    assert "25.4x" in warnings[0] and "mm_per_angstrom >= 1.01" in warnings[0]


def test_unit_guess_warning_boundaries():
    """Slicers test volume < threshold: 9 mm³ is safe, just under is not."""
    assert _unit_guess_warning("a", 9.0, 1.0) is None
    assert "as inches and scale it 25.4x" in _unit_guess_warning("a", 8.99, 1.0)
    assert "as metres and scale it 1000x" in _unit_guess_warning("a", 0.0079, 1.0)
    # Volume grows with the cube of the scale: 0.5 * (9 / 7) ** (1/3) = 0.544.
    assert "mm_per_angstrom >= 0.55" in _unit_guess_warning("a", 7.0, 0.5)
    empty = _unit_guess_warning("a", 0.0, 1.0)
    assert "as metres" in empty and "mm_per_angstrom >=" not in empty


def test_describe_bodies_reports_pieces_and_what_was_dropped():
    info = {
        "bodies": [7624.5, 7570.1],
        "dropped": [],
        "cavities": [3.0, 1.5],
        "debris": [0.0, 0.0],
        "volume": 15194.6,
    }
    assert _describe_bodies(info) == [
        "prints as 2 separate pieces: 7.62 cm³, 7.57 cm³ (bodies='largest' keeps only the first)",
        "filled 2 enclosed cavities (4.5 mm³)",
        "removed 2 stray fragments of debris (0 mm³ in all)",
    ]
    largest = {
        "bodies": [7624.5],
        "dropped": [7570.1],
        "cavities": [],
        "debris": [],
        "volume": 7624.5,
    }
    assert _describe_bodies(largest) == [
        "WARNING: kept the largest piece and dropped 1 more (7.57 cm³, 49.8% of the group). "
        "They are missing from the STL; bodies='all' keeps them as separate pieces."
    ]
    sliver = {**largest, "dropped": [10.0], "volume": 7624.5}
    assert _describe_bodies(sliver) == [
        "kept the largest piece and dropped 1 more (10 mm³, 0.1% of the group)"
    ]


# ── thin features ──


def test_thin_check_names_the_setting_to_change():
    """A 0.66 Å rung was one voxel at the default 0.7 Å pitch, and nothing said so."""
    rows = {"guide and rep cartoon": [[0, "H"], [0, "L"], [4, ""]]}
    settings = {"cartoon_oval_width": 0.25, "cartoon_loop_radius": 0.2, "cartoon_tube_radius": 0.5}
    pymol = _PyMOL(rows=rows, settings=settings)
    with patch("mcpymol.printing.send_request", side_effect=pymol):
        thin = _thinnest_cartoon(["ub"])

    assert thin == (pytest.approx(0.4), "cartoon_loop_radius", 0.2, "ub", "loop tubes")
    warning = _thin_cartoon_warning(thin, 0.7)
    assert warning.startswith(
        "WARNING: the thinnest feature PyMOL draws here is 0.40 Å across "
        "(loop tubes of 'ub', cartoon_loop_radius 0.2): 0.6 voxels at voxel_pitch 0.7 Å."
    )
    assert "Use voxel_pitch <= 0.13" in warning
    assert "`set cartoon_loop_radius, <value>, ub`" in warning
    assert _thin_cartoon_warning(thin, 0.13) is None


def test_thin_check_reads_fancy_helices_as_dumbbells():
    rows = {"guide and rep cartoon": [[0, "H"]]}
    settings = {
        "cartoon_fancy_helices": True,
        "cartoon_dumbbell_width": 0.17,
        "cartoon_dumbbell_radius": 0.16,
    }
    pymol = _PyMOL(rows=rows, settings=settings)
    with patch("mcpymol.printing.send_request", side_effect=pymol):
        thin = _thinnest_cartoon(["ub"])

    assert thin[1:] == ("cartoon_dumbbell_width", 0.17, "ub", "dumbbell helices")


def test_thin_rod_warning_names_valence_rods():
    np = pytest.importorskip("numpy")
    primitives = {
        "cylinders": np.array([[0, 0, 0, 1, 0, 0, 0.25], [0, 0, 0, 1, 0, 0, 0.1]]),
        "spheres": np.array([[0, 0, 0, 0.25]]),
    }
    census = {"lig": {"sticks": 20}}

    warning = _thin_rod_warning(primitives, census, True, False, 1.0)

    assert warning.startswith("WARNING: 3 of 3 rods and spheres are under 0.8 mm across")
    assert "`set valence, 0`" in warning
    assert "mm_per_angstrom >= 4" in warning
    assert _thin_rod_warning(primitives, census, True, False, 4.0) is None


# ── reading PyMOL's files ──

_POV = """\
camera {direction<0.0,0.0,  -2.835>
 location <0.0 , 0.0 , 0.0>
 right   -1.33*x up y
 }
light_source{<4000.0,4000.0,-10000.0> rgb<1.0,1.0,1.0>}
cylinder{<1.0,2.0,-97.0>,
<1.5,2.0,-97.0>,
 0.25 open
pigment{color rgb<0.2,0.2,0.9>}}
cylinder{<0.0,0.0,-100.0>,
<0.0,3.0,-100.0>,
 0.1
pigment{color rgb<0.2,0.2,0.9>}}
sphere{<1.0,2.0,-97.0>, 0.25
pigment{color rgb<0.2,0.2,0.9>}}
mesh2 { vertex_vectors { 3, <0,0,-100>,
<1,0,-100>,
<0,1,-100>}
 face_indices { 1, <0,1,2>}}
smooth_color_triangle{<0,0,-100>,<0,0,1>,<1,1,1>,<1,0,-100>,<0,0,1>,<1,1,1>,<0,1,-100>,<0,0,1>,<1,1,1>}
disc{<0,0,-100>,<0,0,1>,1.0}
"""


def test_parse_pov_reads_spheres_and_rods_and_names_the_rest(tmp_path):
    np = pytest.importorskip("numpy")
    path = tmp_path / "scene.pov"
    path.write_text(_POV)

    geometry = _parse_pov(str(path), -100.0)

    np.testing.assert_allclose(geometry["spheres"], [[1.0, 2.0, 3.0, 0.25]])
    np.testing.assert_allclose(
        geometry["cylinders"], [[1.0, 2.0, 3.0, 1.5, 2.0, 3.0, 0.25], [0, 0, 0, 0, 3, 0, 0.1]]
    )
    assert geometry["cyl_open"].tolist() == [True, False]
    assert geometry["unknown"] == {"disc": 1}


@pytest.mark.parametrize("radius", ["nan", "0.0"])
def test_parse_pov_refuses_rods_without_a_size(tmp_path, radius):
    """With the camera at the origin PyMOL wrote pixel-sized dashes with
    radius NaN (perspective) or 0 (orthoscopic); they must not vanish quietly."""
    pytest.importorskip("numpy")
    path = tmp_path / "scene.pov"
    path.write_text(f"cylinder{{<0,0,0>,\n<1,0,0>,\n {radius} open\npigment{{}}}}\n")

    with pytest.raises(RuntimeError, match="no usable position or size"):
        _parse_pov(str(path))


def test_read_obj_triangles_drops_sphere_markers(tmp_path):
    """PyMOL writes each sphere as one degenerate triangle; voxelised, a
    spheres-only group came out a 0.2 Å³ speck reported OK."""
    pytest.importorskip("numpy")
    path = tmp_path / "scene.obj"
    path.write_text(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 5 5 5\nv 5 5 5\nv 5 5 5\n"
        "vn 0 0 1\nf 1//1 2//1 3//1\nf 4//1 5//1 6//1\n"
    )

    triangles, markers = _read_obj_triangles(str(path))

    assert triangles.shape == (1, 3, 3) and markers == 1
    path.write_text("v nan 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n")
    with pytest.raises(RuntimeError, match="NaN"):
        _read_obj_triangles(str(path))


# ── sticks and spheres ──


def test_merge_continuing_rods_joins_half_bonds():
    """PyMOL draws a bond as two half-cylinders, one per atom colour."""
    np = pytest.importorskip("numpy")
    pytest.importorskip("scipy")
    halves = [[0, 0, 0, 1, 0, 0, 0.25], [1, 0, 0, 2, 0, 0, 0.25], [2, 0, 0, 3, 0, 0, 0.25]]

    rods, is_open = _merge_continuing_rods(halves, [True, True, True])

    np.testing.assert_allclose(rods, [[0, 0, 0, 3, 0, 0, 0.25]])
    assert is_open.tolist() == [True]
    # Different radii, a bend, a three-way joint, a capped end: all left alone.
    for rods_in, opened in (
        ([[0, 0, 0, 1, 0, 0, 0.25], [1, 0, 0, 2, 0, 0, 0.1]], [True, True]),
        ([[0, 0, 0, 1, 0, 0, 0.25], [1, 0, 0, 1, 1, 0, 0.25]], [True, True]),
        (
            [[0, 0, 0, 1, 0, 0, 0.25], [1, 0, 0, 2, 0, 0, 0.25], [1, 0, 0, 1, 1, 0, 0.25]],
            [True, True, True],
        ),
        ([[0, 0, 0, 1, 0, 0, 0.25], [1, 0, 0, 2, 0, 0, 0.25]], [True, False]),
    ):
        assert len(_merge_continuing_rods(rods_in, opened)[0]) == len(rods_in)
    # The same rod drawn twice (two objects) is one rod.
    twice = [[0, 0, 0, 1, 0, 0, 0.25], [1, 0, 0, 0, 0, 0, 0.25]]
    assert len(_merge_continuing_rods(twice, [True, True])[0]) == 1


def _stl_bodies(mesh, tmp_path):
    """Write and read back an STL, as a slicer reads it: vertices merged by position."""
    import trimesh

    path = tmp_path / "round_trip.stl"
    mesh.export(path)
    back = trimesh.load(path, force="mesh")
    return back, back.split(only_watertight=False)


def _half_bond_chain(n=12):
    """Sticks along a helix of n atoms away from the origin, as PyMOL writes
    them: two open half-bonds per bond, float32 coordinates (6 decimals)."""
    import numpy as np

    t = np.arange(n) * 1.9
    atoms = np.column_stack([1.2 * np.cos(t) + 13.37, 1.2 * np.sin(t) - 7.21, 0.55 * t + 21.9])
    rods = []
    for i in range(n - 1):
        a, b = atoms[i], atoms[i + 1]
        middle = (a + b) / 2
        rods += [[*a, *middle, 0.25], [*middle, *b, 0.25]]
    return {
        "spheres": np.round(np.column_stack([atoms, np.full(n, 0.25)]), 6),
        "cylinders": np.round(np.array(rods), 6),
        "cyl_open": np.ones(len(rods), dtype=bool),
    }


def test_union_of_half_bonds_survives_an_stl_round_trip(tmp_path):
    """Half-bonds extended to overlap share their side walls, and the union
    left coincident slivers: closed in memory, but MK1's STL fell into 573
    pieces once a slicer merged its vertices by position."""
    np = pytest.importorskip("numpy")
    pytest.importorskip("trimesh")
    pytest.importorskip("manifold3d")
    chain = _half_bond_chain()

    back, bodies = _stl_bodies(_union_primitives(chain), tmp_path)
    assert back.is_watertight and len(bodies) == 1

    # The test can tell: union the halves as they come and the file falls apart.
    def as_they_come(cylinders, cyl_open, tol=1e-4):
        return np.asarray(cylinders, dtype=float), np.asarray(cyl_open, dtype=bool)

    with patch("mcpymol.printing._merge_continuing_rods", side_effect=as_they_come):
        back, bodies = _stl_bodies(_union_primitives(chain), tmp_path)
    assert not back.is_watertight or len(bodies) > 1


def test_union_encloses_the_true_spheres():
    """Circumscribed: a stick never prints thinner than PyMOL drew it."""
    np = pytest.importorskip("numpy")
    pytest.importorskip("trimesh")
    pytest.importorskip("manifold3d")
    # Two unit spheres one radius apart: 2 x 4/3 pi - pi (4 + 1)(2 - 1)^2 / 12.
    exact = 8 * math.pi / 3 - 5 * math.pi / 12
    spheres = {
        "spheres": np.array([[0, 0, 0, 1.0], [1, 0, 0, 1.0]]),
        "cylinders": np.zeros((0, 7)),
        "cyl_open": np.zeros(0, dtype=bool),
    }

    union = _union_primitives(spheres)

    assert union.is_watertight and len(union.split(only_watertight=False)) == 1
    assert exact <= union.volume <= exact * 1.04


def test_repair_unions_primitives_and_never_sends_them_through_poisson(tmp_path):
    pytest.importorskip("trimesh")
    pytest.importorskip("manifold3d")

    chain = _half_bond_chain()

    info = _repair_to_stl(None, str(tmp_path / "chain.stl"), "auto", 0.7, 10, primitives=chain)

    assert info["method"] == "exact union of 22 rods and 12 spheres"
    assert info["watertight"] and len(info["bodies"]) == 1
    with pytest.raises(ValueError, match="Poisson reconstruction cannot merge sticks"):
        _repair_to_stl(None, str(tmp_path / "x.stl"), "poisson", 0.7, 10, primitives=chain)


# ── bodies ──


def _obj_file(path, *meshes):
    import trimesh

    trimesh.util.concatenate(list(meshes)).export(path)
    return str(path)


def test_voxel_keeps_every_separate_body(tmp_path):
    """Keeping the body with the most faces dropped one of 4hhb's two dimers,
    or three of four hemes, and still said OK."""
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    spheres = [
        trimesh.creation.icosphere(subdivisions=3, radius=3.0),
        trimesh.creation.icosphere(subdivisions=3, radius=3.0).apply_translation([12, 0, 0]),
        trimesh.creation.icosphere(subdivisions=3, radius=2.0).apply_translation([0, 12, 0]),
    ]
    src = _obj_file(tmp_path / "three.obj", *spheres)

    info = _repair_to_stl(src, str(tmp_path / "three.stl"), "voxel", 0.5, 10)

    assert len(info["bodies"]) == 3
    assert len(trimesh.load(tmp_path / "three.stl").split(only_watertight=False)) == 3


def test_light_fills_a_cavity_and_keeps_separate_pieces(tmp_path):
    """PyMOL winds an enclosed cavity's shell inward; it prints filled, while
    a separate ligand stays as its own piece."""
    trimesh = pytest.importorskip("trimesh")
    outer = trimesh.creation.icosphere(subdivisions=3, radius=10.0)
    cavity = trimesh.creation.icosphere(subdivisions=3, radius=3.0)
    cavity.invert()
    ligand = trimesh.creation.icosphere(subdivisions=3, radius=2.0).apply_translation([30, 0, 0])
    src = _obj_file(tmp_path / "cavity.obj", outer, cavity, ligand)

    info = _repair_to_stl(src, str(tmp_path / "cavity.stl"), "auto", 0.7, 10)

    assert info["method"] == "light (already watertight)"
    assert len(info["bodies"]) == 2
    assert info["cavities"] == [pytest.approx(-cavity.volume)]
    stl = trimesh.load(tmp_path / "cavity.stl")
    assert stl.is_watertight and len(stl.split(only_watertight=False)) == 2
    assert stl.volume == pytest.approx(outer.volume + ligand.volume)


def test_debris_is_dropped_and_counted(tmp_path):
    trimesh = pytest.importorskip("trimesh")
    sphere = trimesh.creation.icosphere(subdivisions=2, radius=3.0)
    stray = trimesh.Trimesh([[20, 0, 0], [21, 0, 0], [20, 1, 0]], [[0, 1, 2]])
    src = _obj_file(tmp_path / "debris.obj", sphere, stray)

    info = _repair_to_stl(src, str(tmp_path / "debris.stl"), "auto", 0.7, 10)

    assert info["method"] == "light (already watertight)"
    assert len(info["debris"]) == 1 and len(info["bodies"]) == 1
    assert trimesh.load(tmp_path / "debris.stl").is_watertight


def test_bodies_largest_keeps_one_and_says_what_it_dropped(tmp_path):
    trimesh = pytest.importorskip("trimesh")
    a = trimesh.creation.icosphere(subdivisions=2, radius=3.0)
    b = trimesh.creation.icosphere(subdivisions=2, radius=2.9).apply_translation([10, 0, 0])
    src = _obj_file(tmp_path / "two.obj", a, b)

    info = _repair_to_stl(src, str(tmp_path / "two.stl"), "auto", 0.7, 10, bodies="largest")

    assert info["bodies"] == [pytest.approx(a.volume)]
    assert info["dropped"] == [pytest.approx(b.volume)]
    assert _describe_bodies(info)[0].startswith(
        "WARNING: kept the largest piece and dropped 1 more"
    )
    assert len(trimesh.load(tmp_path / "two.stl").split(only_watertight=False)) == 1


def test_select_bodies_turns_an_unenclosed_inward_shell_out():
    """An inward shell inside no solid is a solid wound the wrong way, not a cavity."""
    trimesh = pytest.importorskip("trimesh")
    lone = trimesh.creation.icosphere(subdivisions=2, radius=2.0)
    lone.invert()

    kept, stats = _select_bodies(lone, None, "all")

    assert len(kept) == 1 and kept[0].volume > 0
    assert stats == {"dropped": [], "cavities": [], "debris": []}


def _mean_radius(mesh):
    import numpy as np

    points, _ = __import__("trimesh").sample.sample_surface(mesh, 20000, seed=0)
    return float(np.linalg.norm(points - mesh.centroid, axis=1).mean())


def test_voxel_growth_is_undone(tmp_path):
    """Voxelising grew every surface by about 0.62 voxel: a 0.33 Å rod printed
    0.38 Å at 0.1 Å pitch, and every re-voxelisation grew it again."""
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    sphere = trimesh.creation.icosphere(subdivisions=4, radius=2.0)
    src = _obj_file(tmp_path / "sphere.obj", sphere)
    pitch = 0.2

    _repair_to_stl(src, str(tmp_path / "in.stl"), "voxel", pitch, 10)
    _repair_to_stl(src, str(tmp_path / "grown.stl"), "voxel", pitch, 10, undo_voxel_growth=False)

    corrected = trimesh.load(tmp_path / "in.stl")
    grown = trimesh.load(tmp_path / "grown.stl")
    assert abs(_mean_radius(corrected) - 2.0) < 0.1 * pitch
    assert _mean_radius(grown) - 2.0 > 0.4 * pitch
    assert corrected.is_watertight and corrected.volume > 0

    # Re-voxelising the corrected part leaves it where it is.
    again = _obj_file(tmp_path / "again.obj", corrected)
    _repair_to_stl(again, str(tmp_path / "again.stl"), "voxel", pitch, 10)
    assert abs(_mean_radius(trimesh.load(tmp_path / "again.stl")) - _mean_radius(corrected)) < 0.02


def test_voxel_inset_never_turns_a_thin_part_inside_out(tmp_path):
    """A flat face moved along its own normal never turns over, so the fold
    check missed a thin part whose two sides passed through each other: a
    0.2 Å slab, one voxel thick at 0.3 Å, came out inside out."""
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    slab = trimesh.creation.box(extents=[6.0, 6.0, 0.2]).apply_translation([0.13, -0.27, 0.31])
    src = _obj_file(tmp_path / "slab.obj", slab)

    _repair_to_stl(src, str(tmp_path / "slab.stl"), "voxel", 0.3, 10)

    out = trimesh.load(tmp_path / "slab.stl")
    assert out.is_watertight and out.volume > 0
    # The test can tell: without the thickness cap the slab comes out inside out.
    with patch("mcpymol.printing._room_to_inset", return_value=None):
        _repair_to_stl(src, str(tmp_path / "flipped.stl"), "voxel", 0.3, 10)
    assert trimesh.load(tmp_path / "flipped.stl").volume < 0


def test_inset_backs_off_where_a_face_would_fold():
    """Moved in further than its radius, a thin rod's faces turn over."""
    np = pytest.importorskip("numpy")
    trimesh = pytest.importorskip("trimesh")
    rod = trimesh.creation.cylinder(radius=0.2, height=4.0, sections=16)
    before = rod.face_normals.copy()

    held = _inset(rod, 0.3)

    assert held > 0
    assert (np.einsum("ij,ij->i", before, rod.face_normals) > 0).all()
    assert rod.volume > 0


# ── residue coverage ──


def test_surface_index_needs_no_rtree():
    """trimesh's own proximity queries need rtree, which the print extra lacks."""
    np = pytest.importorskip("numpy")
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    sphere = trimesh.creation.icosphere(subdivisions=4, radius=5.0)
    points = np.array([[0, 0, 0], [0, 0, 4.0], [0, 6.0, 0], [8.0, 0, 0]])

    with patch.dict(sys.modules, {"rtree": None}):
        distance = _SurfaceIndex(sphere, max_edge=0.8).signed_distance(points)

    np.testing.assert_allclose(distance, [-5.0, -1.0, 1.0, 3.0], atol=0.05)


def _coverage_pymol(tmp_path):
    """A 20 Å rod of radius 1 along x, and four Cα: on the axis, 1.5 Å out
    where loops are smoothed, 1.5 Å out next to a hidden residue, 5 Å out."""
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    rod = trimesh.creation.cylinder(radius=1.0, height=20.0, sections=48)
    rod.apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [0, 1, 0]))
    stl = tmp_path / "rod.stl"
    rod.export(stl)
    rows = {
        "byres": [["A", "3"]],
        "guide and rep cartoon and present": [
            ["A", "1", "GLY"],
            ["A", "2", "ALA"],
            ["A", "3", "LEU"],
            ["A", "4", "SER"],
        ],
    }
    pymol = _PyMOL(rows=rows, settings={"cartoon_smooth_loops": True})
    pymol._get_coords = lambda args, kwargs: [[-6, 0, 0], [-2, 0, 2.5], [2, 0, 2.5], [6, 0, 6]]
    return pymol, str(stl)


def test_coverage_flags_residues_outside_the_mesh(tmp_path):
    """Leu64's Cα sat 1.84 Å outside the split-GFP tube, beside a hidden
    residue, and nothing said so: the part meant to attach there floated."""
    pymol, stl = _coverage_pymol(tmp_path)
    with patch("mcpymol.printing.send_request", side_effect=pymol):
        warnings = _coverage_warnings(["ub"], stl, 1.0)

    assert len(warnings) == 2
    assert warnings[0].startswith(
        "WARNING: 1 shown residue(s) next to a hidden one have their Cα outside the mesh: "
        "ub A/LEU3 (1.50 Å out)"
    )
    assert "`set cartoon_smooth_loops, 0, ub`" in warnings[0]
    assert warnings[1].startswith(
        "WARNING: 1 residue(s) PyMOL shows as cartoon are missing from the mesh, up to 5.0 Å away"
    )
    assert "ub A/SER4" in warnings[1]


def test_coverage_selections_skip_atoms_without_coordinates(tmp_path):
    """PyMOL 3 loads unobserved residues as atoms without coordinates: 2b3p
    listed 243 residues and 228 coordinates, and every label slipped."""
    pymol, stl = _coverage_pymol(tmp_path)
    with patch("mcpymol.printing.send_request", side_effect=pymol):
        _coverage_warnings(["ub"], stl, 1.0)

    selections = [
        args[0] for name, args, _ in pymol.calls if name in ("iterate_to_list", "get_coords")
    ]
    assert selections and all("present" in sel for sel in selections)


# ── Poisson in a child process ──


def test_a_poisson_crash_does_not_take_the_server_down(tmp_path):
    """PoissonRecon ends its process with exit(0) when it gives up: 4 of 14
    runs on one print-ribbon export did, and in the server that dropped the
    MCP connection mid-call. In a child process it is an ordinary error."""
    with (
        patch("mcpymol.printing._POISSON_CHILD", "import os; os._exit(0)"),
        patch("importlib.util.find_spec", return_value=object()),
        pytest.raises(RuntimeError, match="Poisson reconstruction failed and wrote nothing"),
    ):
        _run_poisson(str(tmp_path / "in.obj"), str(tmp_path / "out.stl"), 10)


def test_auto_falls_back_to_voxel_when_poisson_dies(tmp_path):
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    # An open surface, so auto tries Poisson first.
    sphere = trimesh.creation.icosphere(subdivisions=2, radius=3.0)
    sphere.update_faces(range(len(sphere.faces) - 4))
    src = _obj_file(tmp_path / "open.obj", sphere)

    with (
        patch("mcpymol.printing._POISSON_CHILD", "import os; os._exit(0)"),
        patch("importlib.util.find_spec", return_value=object()),
    ):
        info = _repair_to_stl(src, str(tmp_path / "open.stl"), "auto", 0.5, 10)

    # The reason rides along, and the pitch after it.
    assert info["method"].startswith("voxel (poisson fallback: ")
    assert "Poisson reconstruction failed" in info["method"]
    assert "), 0.5 Å pitch" in info["method"]


def test_poisson_child_never_imports_from_the_servers_working_directory(tmp_path, monkeypatch):
    """`python -c` puts the working directory first on sys.path, so a project
    file named like a stdlib module was imported, and run, in the child:
    Poisson failed every time, and auto fell back to voxel without saying why."""
    trimesh = pytest.importorskip("trimesh")
    pytest.importorskip("pymeshlab")
    project, work = tmp_path / "project", tmp_path / "work"
    project.mkdir()
    work.mkdir()
    ran = tmp_path / "the project calendar.py ran"
    (project / "calendar.py").write_text(
        f"open({str(ran)!r}, 'w').close()\nraise ImportError('the project calendar.py')\n"
    )
    monkeypatch.chdir(project)
    src = _obj_file(work / "ball.obj", trimesh.creation.icosphere(subdivisions=2, radius=3.0))

    _run_poisson(src, str(work / "ball.stl"), 6)

    assert (work / "ball.stl").stat().st_size > 0
    assert not ran.exists()


def test_scale_is_applied_about_the_origin(tmp_path):
    """About the origin, never a part's own centre, so every group keeps its place."""
    np = pytest.importorskip("numpy")
    trimesh = pytest.importorskip("trimesh")
    cube = trimesh.creation.box(extents=[2.0, 2.0, 2.0]).apply_translation([10, -5, 3])
    src = _obj_file(tmp_path / "cube.obj", cube)

    info = _repair_to_stl(src, str(tmp_path / "cube.stl"), "auto", 0.7, 10, mm_per_angstrom=2.0)

    np.testing.assert_allclose(
        trimesh.load(tmp_path / "cube.stl").bounds, [[18, -12, 4], [22, -8, 8]]
    )
    np.testing.assert_allclose(info["bounds"], [[18, -12, 4], [22, -8, 8]])
    assert info["volume"] == pytest.approx(64.0)


# ── print_ribbon_view ────────────────────────────────────────────────────────

# What the default style sets, on the object and on its spine alike.
_DEFAULT_STYLE = {
    "cartoon_loop_radius": "1.3",
    "cartoon_rect_width": "2.0",
    "cartoon_rect_length": "2.8",
    "cartoon_oval_width": "1.3",
    "cartoon_oval_length": "2.0",
    "cartoon_helix_radius": "1.3",
    "cartoon_fancy_sheets": "1",
    "cartoon_flat_sheets": "1",
    "cartoon_smooth_loops": "1",
}

# GFP's chromophore as 2b3p loads into PyMOL: CRO A66 has CA1, CA2 and CA3
# and no CA; its N1 bonds Leu64's C and its C3 bonds Val68's N.
_CRO = {
    "entry": [483, "A", "CRO", "66", "A"],
    "exit": [487, "A", "CRO", "66", "A"],
    "atoms": [[468, "CA1"], [469, "CA2"], [470, "CA3"], [483, "N1"], [487, "C3"]],
    "before": 461,
    "after": 491,
}


def _ca_rows(strands, last=20, first_index=101, step=1):
    """iterate_to_list rows (chain, segi, resv, resi, ss, index) for residues 1..last."""
    return [
        [
            "A",
            "",
            resv,
            str(resv),
            "S" if any(a <= resv <= b for a, b in strands) else "L",
            first_index + step * (resv - 1),
        ]
        for resv in range(1, last + 1)
    ]


def _ribbon_pymol(ca_rows=(), link=None, traced=0, unbonded=0, fail=None, warn=None, redraw=None):
    """A send_request stand-in that answers print_ribbon_view's queries.

    ``link`` is a residue with no CA (see _CRO); ``traced`` atoms already have
    cartoon_trace_atoms on; ``unbonded`` Cα have no bonds; ``fail`` is an
    action that errors, ``warn`` one that succeeds with a Selector-Error, and
    ``redraw`` one that succeeds carrying an error PyMOL printed while redrawing.
    """

    def ok(result):
        return {"status": "success", "result": result}

    def fake(action, args=None, kwargs=None, **_):
        args = args or []
        if action == fail:
            return {"status": "error", "error": f"{action} went wrong"}
        if action == warn:
            return {
                "status": "success",
                "result": [],
                "warnings": ['Selector-Error: Invalid selection name "x".'],
            }
        if action == redraw:
            return {
                "status": "success",
                "result": "OK",
                "warnings": ["ObjectSurfaceUpdate-Error: surface of 'other' is out of date"],
            }
        if action == "get_type":
            return ok("object:molecule")
        if action == "count_atoms":
            return ok(unbonded)
        if action != "iterate_to_list":
            return ok("OK")
        selection, expression = args
        if expression == "chain, segi, resv, resi, ss, index":
            return ok([list(row) for row in ca_rows])
        if expression == "s.cartoon_trace_atoms":
            return ok([[1]] * traced + [[0]] * 40)
        if link is None:
            return ok([])
        if expression == "index, chain, resn, resi, segi":
            return ok([link["entry"] if "name C))" in selection else link["exit"]])
        if expression == "index, name":
            return ok(link["atoms"])
        if expression == "index":
            return ok([[link["before"] if "name C)" in selection else link["after"]]])
        return ok([])

    return fake


def _sent(mock_sr, action):
    return [c.kwargs["args"] for c in mock_sr.call_args_list if c.args[0] == action]


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_sets_style_per_object_not_globally(mock_sr):
    """The style used to go out as `set name, value` with no object: the
    global value, so every other object in the session, and every one loaded
    later, took the chunky strands too."""
    mock_sr.side_effect = _ribbon_pymol()

    result = print_ribbon_view(obj_name="1UBQ")

    sets = _sent(mock_sr, "set")
    assert sets
    assert all(len(args) == 3 for args in sets), sets
    assert {args[2] for args in sets} == {"1UBQ", "1UBQ_spine"}
    assert "other objects and the global settings are unchanged" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_defaults_match_todays_style(mock_sr):
    """The defaults keep the look, value for value: set per object, they
    export byte-identical OBJs to the old global style (1ubq, 4hhb, 2b3p)."""
    mock_sr.side_effect = _ribbon_pymol()

    print_ribbon_view(obj_name="1UBQ")

    by_target: dict[str, dict[str, str]] = {}
    for setting, value, target in _sent(mock_sr, "set"):
        by_target.setdefault(target, {})[setting] = value
    assert by_target["1UBQ"] == _DEFAULT_STYLE
    assert by_target["1UBQ_spine"] == {**_DEFAULT_STYLE, "cartoon_tube_radius": "0.9"}
    # Nothing to bridge and the spine left in the strands: no trace, no pins,
    # nothing of the spine hidden.
    acts = [c.args[0] for c in mock_sr.call_args_list]
    assert "flag" not in acts and "rebuild" not in acts
    assert not [a for a in _sent(mock_sr, "hide") if a[0] == "cartoon" and "1UBQ_spine" in a[1]]


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_strand_dimensions(mock_sr):
    """Strands 2.8 Å wide and 2.0 Å thick in name were drawn 4.0 x 2.8 Å and
    closed GFP's barrel; there was no way to make them thinner."""
    mock_sr.side_effect = _ribbon_pymol()

    result = print_ribbon_view(obj_name="1UBQ", strand_width=2.0, strand_thickness=1.0)

    for target in ("1UBQ", "1UBQ_spine"):
        assert ["cartoon_rect_length", "2.0", target] in _sent(mock_sr, "set")
        assert ["cartoon_rect_width", "1.0", target] in _sent(mock_sr, "set")
    # The size PyMOL actually draws, not the setting.
    assert "about 2.8 Å wide x 1.4 Å thick" in result
    # A 1.8 Å spine through 1.4 Å strands bulges out of them.
    assert "spine_in_strands=False" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_hides_spine_on_strand_interiors(mock_sr):
    """spine_in_strands=False takes the spine out of each strand but keeps two
    residues at each end, where the loops fuse to the strand."""
    # Residues 1-20 are Cα indices 101-120; 4-12 is a strand, 15-18 a short one.
    mock_sr.side_effect = _ribbon_pymol(ca_rows=_ca_rows([(4, 12), (15, 18)]))

    result = print_ribbon_view(obj_name="1UBQ", spine_in_strands=False)

    spine_hides = [args for args in _sent(mock_sr, "hide") if "1UBQ_spine" in args[1]]
    spine_hides = [args for args in spine_hides if args[0] == "cartoon"]
    # Residues 6-10 of 4-12; the four-residue strand keeps all of its spine.
    assert spine_hides == [["cartoon", "byres ((1UBQ_spine) and index 106+107+108+109+110)"]]
    assert "interior of 1 strand (2 residues kept at each end" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_rejects_an_overlap_below_two(mock_sr):
    """With one residue kept per end, 10 of GFP's 22 strand ends lost their
    junction and the barrel printed as 6 pieces."""
    result = print_ribbon_view(obj_name="1UBQ", spine_in_strands=False, strand_overlap=1)

    assert result.startswith("Error: strand_overlap must be at least 2")
    mock_sr.assert_not_called()


@pytest.mark.parametrize(
    "kwargs", [{"strand_thickness": 0.0}, {"strand_width": -1.0}, {"spine_radius": float("nan")}]
)
@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_rejects_sizes_that_are_not_positive(mock_sr, kwargs):
    result = print_ribbon_view(obj_name="1UBQ", **kwargs)

    assert result.startswith("Error:") and "positive" in result
    mock_sr.assert_not_called()


def test_strand_interiors_reads_alt_locs_and_insertion_codes_as_one_strand():
    """A residue with two alternate locations is one residue, and 52, 52A, 53
    follow one another; either used to cut the run in two."""
    rows = [
        ["A", "", 50, "50", "S", 1],
        ["A", "", 51, "51", "S", 2],
        ["A", "", 52, "52", "S", 3],
        ["A", "", 52, "52", "S", 4],  # alternate location of 52
        ["A", "", 52, "52A", "S", 5],
        ["A", "", 53, "53", "S", 6],
        ["A", "", 54, "54", "S", 7],
        ["A", "", 55, "55", "S", 8],
        ["B", "", 56, "56", "S", 9],  # another chain
    ]
    with patch("mcpymol.printing.send_request", return_value={"status": "success", "result": rows}):
        assert _strand_interiors("x", 2, set(), []) == [[3, 5, 6]]


def test_strand_interiors_cuts_a_strand_at_a_bridge():
    """A bridge must never sit next to a hidden residue: trace_atoms_mode 0
    joins it to the last *shown* Cα, wherever that is. Insertion codes can
    number the residues either side of one consecutively."""
    rows = _ca_rows([(1, 20)])
    with patch("mcpymol.printing.send_request", return_value={"status": "success", "result": rows}):
        # 111 (residue 11) follows a residue with no CA.
        runs = _strand_interiors("x", 2, {111}, [])
    assert runs == [list(range(103, 109)), list(range(113, 119))]


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_bridges_a_residue_with_no_ca(mock_sr):
    """GFP's chromophore has CA1-3 but no CA, so the spine broke at 64/68 and
    thin prints lost 30% of the barrel. The bridge traces it, on the spine
    only, without renaming anything."""
    mock_sr.side_effect = _ribbon_pymol(link=_CRO)

    result = print_ribbon_view(obj_name="1UBQ")

    sets = _sent(mock_sr, "set")
    assert ["cartoon_trace_atoms", "1", "(1UBQ_spine) and index 468+469+470+491"] in sets
    assert ["trace_atoms_mode", "0", "1UBQ_spine"] in sets
    assert not [
        s
        for s in sets
        if s[0] in ("cartoon_trace_atoms", "trace_atoms_mode") and "1UBQ_spine" not in s[2]
    ]
    # Loop smoothing pulled the tube 1.6 Å off Leu64's Cα; the Cα either side is held.
    assert _sent(mock_sr, "flag") == [["no_smooth", "(1UBQ_spine) and index 461+491", "set"]]
    assert _sent(mock_sr, "rebuild") == [["1UBQ_spine"]]
    assert not _sent(mock_sr, "alter")
    assert "Spine bridged through 1 residue with no CA atom: CRO A66" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_reports_a_break_it_was_told_not_to_bridge(mock_sr):
    mock_sr.side_effect = _ribbon_pymol(link=_CRO)

    result = print_ribbon_view(obj_name="1UBQ", bridge_modified=False)

    assert not [
        s for s in _sent(mock_sr, "set") if s[0] in ("cartoon_trace_atoms", "trace_atoms_mode")
    ]
    assert not _sent(mock_sr, "flag")
    assert "Spine breaks at CRO A66 (no CA atom); bridge_modified=True joins it." in result


@pytest.mark.parametrize("already", [{"traced": 2}, {"unbonded": 413}])
@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_will_not_bridge_a_spine_that_already_traces(mock_sr, already):
    """trace_atoms_mode 0 joins every traced atom to the shown Cα before it,
    across chain breaks too, so a spine with traced atoms of its own (set by
    the user, or CA-only residues like 1grl's 413) is left unbridged."""
    mock_sr.side_effect = _ribbon_pymol(link=_CRO, **already)

    result = print_ribbon_view(obj_name="1UBQ")

    assert not [
        s for s in _sent(mock_sr, "set") if s[0] in ("cartoon_trace_atoms", "trace_atoms_mode")
    ]
    assert "Spine breaks at CRO A66 (no CA atom) and was not bridged" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_leaves_a_cap_alone(mock_sr):
    """A residue linked on one side only, such as an N-terminal ACE or a
    C-terminal NH2 cap, ends its chain: there is nothing to bridge."""
    caps = {**_CRO, "entry": [700, "A", "NH2", "77", "A"], "exit": [10, "A", "ACE", "0", "A"]}
    mock_sr.side_effect = _ribbon_pymol(link=caps)

    result = print_ribbon_view(obj_name="1UBQ")

    assert not [s for s in _sent(mock_sr, "set") if s[0] == "cartoon_trace_atoms"]
    assert "NH2" not in result and "ACE" not in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_will_not_trace_atoms_stored_out_of_order(mock_sr):
    """PyMOL scans guide atoms in atom order: a residue whose atoms sit after
    its neighbour's would be joined to the wrong residues."""
    late = {**_CRO, "atoms": [[900, "CA1"], [901, "CA2"], [902, "CA3"]]}
    mock_sr.side_effect = _ribbon_pymol(link=late)

    result = print_ribbon_view(obj_name="1UBQ")

    assert not [s for s in _sent(mock_sr, "set") if s[0] == "cartoon_trace_atoms"]
    assert "Spine breaks at CRO A66 (no CA atom): not bridged" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_keeps_the_spine_shown_either_side_of_a_bridge(mock_sr):
    """With strands on both sides of the bridge, the residues next to it are
    strand ends and keep their spine."""
    # Residue n's Cα is atom 10n. The chromophore's atoms (112-114) sit
    # between residue 11 and residue 12, numbered as if nothing were between.
    rows = _ca_rows([(1, 30)], last=30, first_index=10, step=10)
    link = {
        **_CRO,
        "atoms": [[112, "CA1"], [113, "CA2"], [114, "CA3"]],
        "before": 110,
        "after": 120,
    }
    mock_sr.side_effect = _ribbon_pymol(ca_rows=rows, link=link)

    result = print_ribbon_view(obj_name="1UBQ", spine_in_strands=False)

    assert "Spine bridged through 1 residue" in result
    hidden = [args[1] for args in _sent(mock_sr, "hide") if "1UBQ_spine) and index" in args[1]]
    assert len(hidden) == 1
    indices = {int(i) for i in hidden[0].split("index ")[1].rstrip(")").split("+")}
    assert not indices & {100, 110, 120, 130}
    assert indices == set(range(30, 100, 10)) | set(range(140, 290, 10))


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_sends_no_do(mock_sr):
    """`do` is queued onto PyMOL's main thread by older plugins, which answer
    before it runs: its errors were lost, and it could run out of order."""
    mock_sr.side_effect = _ribbon_pymol(ca_rows=_ca_rows([(4, 12)]), link=_CRO)

    print_ribbon_view(obj_name="1UBQ", spine_in_strands=False)

    assert "do" not in [c.args[0] for c in mock_sr.call_args_list]
    assert ["util.chainbow", ["1UBQ"]] in [
        [c.args[0], c.kwargs["args"]] for c in mock_sr.call_args_list
    ]


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_fails_on_a_missing_object(mock_sr):
    """It used to return its full success text for an object that does not exist."""
    mock_sr.return_value = {"status": "error", "error": "object not found"}

    result = print_ribbon_view(obj_name="no_such_object")

    assert result.startswith("Error: 'no_such_object' is not an object")
    assert [c.args[0] for c in mock_sr.call_args_list] == ["get_type"]


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_refuses_a_selection(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "selection"}

    result = print_ribbon_view(obj_name="sele")

    assert result.startswith("Error: 'sele' is a selection")


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_stops_at_the_first_failed_call(mock_sr):
    mock_sr.side_effect = _ribbon_pymol(fail="create")

    result = print_ribbon_view(obj_name="1UBQ")

    assert result.startswith("Error: `create 1UBQ_spine, 1UBQ` failed: create went wrong")
    # Not rolled back, so it says so, and that a re-run puts it right.
    assert result.endswith(
        "1UBQ may be left partly restyled. print_ribbon_view starts from scratch each "
        "time, so re-running it once the cause is fixed puts it right."
    )
    acts = [c.args[0] for c in mock_sr.call_args_list]
    assert acts[-1] == "create"
    assert "util.chainbow" not in acts


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_treats_a_selector_error_as_a_failure(mock_sr):
    """iterate and orient report a bad selection only as a printed
    Selector-Error, with an empty, successful result."""
    mock_sr.side_effect = _ribbon_pymol(warn="iterate_to_list")

    result = print_ribbon_view(obj_name="1UBQ")

    assert result.startswith("Error:") and "Selector-Error" in result


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_counts_only_carbon_as_an_untraceable_ca(mock_sr):
    """A calcium ion is also named CA and PyMOL leaves it unbonded, but PyMOL
    only takes a carbon as a CA guide atom: counting the four Ca2+ of a
    calcium sensor (GCaMP2, 3EVR) refused its chromophore bridge."""
    mock_sr.side_effect = _ribbon_pymol(link=_CRO)

    print_ribbon_view(obj_name="1UBQ")

    counted = [args[0] for args in _sent(mock_sr, "count_atoms") if "bound_to" in args[0]]
    assert counted and all("name CA and elem C" in sel for sel in counted)


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_promises_no_bridge_it_could_not_make(mock_sr):
    """With bridge_modified=False on a spine that already traces atoms, the
    old hint 'bridge_modified=True joins it' was false."""
    mock_sr.side_effect = _ribbon_pymol(link=_CRO, traced=2)

    result = print_ribbon_view(obj_name="1UBQ", bridge_modified=False)

    assert "bridge_modified=True joins it" not in result
    assert "was not bridged" in result


@pytest.mark.parametrize(
    "error",
    [
        "Socket connection failed: nothing is listening on 127.0.0.1:9876 (refused). "
        "Is the PyMOL plugin running?",
        "PyMOL accepted 'get_type' but sent no reply within 10 s. It may still be running it.",
    ],
)
@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_passes_on_why_it_could_not_check(mock_sr, error):
    """PyMOL down or busy used to read as "'1UBQ' is not an object", inviting an
    agent to re-load a structure that is already there."""
    mock_sr.return_value = {"status": "error", "error": error}

    result = print_ribbon_view(obj_name="1UBQ")

    assert result == f"Error: could not check '1UBQ' in PyMOL: {error}"


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_reports_a_redraw_error_without_failing(mock_sr):
    """PyMOL prints some errors while it redraws, and they land in whatever
    request is running; only a Selector-Error is the call's own failure."""
    mock_sr.side_effect = _ribbon_pymol(redraw="orient")

    result = print_ribbon_view(obj_name="1UBQ")

    assert result.startswith("Print-ribbon view applied to 1UBQ")
    assert "PyMOL warned: ObjectSurfaceUpdate-Error" in result


@pytest.mark.parametrize(("radius", "pitch"), [(0.9, "0.2"), (0.3, "0.19"), (0.25, "0.16")])
@patch("mcpymol.printing.send_request")
def test_print_ribbon_view_suggests_a_pitch_print_export_accepts(mock_sr, radius, pitch):
    """The suggested export must not draw print_export's own thin-feature
    warning: 3 voxels across the thinnest thing drawn."""
    mock_sr.side_effect = _ribbon_pymol()

    result = print_ribbon_view(obj_name="1UBQ", spine_radius=radius)

    assert f"voxel_pitch={pitch})" in result
    assert 'groups="1UBQ=(1UBQ or 1UBQ_spine)", representation="cartoon", method="voxel"' in result


def test_print_export_names_the_files_it_wrote_before_a_later_group_failed(tmp_path):
    """It returned only the later group's error, while the earlier STLs sat in
    out_dir unreported, ready to be dragged into a slicer as a complete set."""
    pymol = _PyMOL(objects={"pa": "molecule", "pb": "molecule"})
    repairs = iter([_repaired(), RuntimeError("mesh exploded")])

    def repair(*args, **kwargs):
        outcome = next(repairs)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock(), "manifold3d": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * 12, 0)),
        patch("mcpymol.printing._repair_to_stl", side_effect=repair),
    ):
        result = print_export(
            obj_name="ub",
            groups="a=pa; b=pb",
            out_dir=str(tmp_path),
            representation="cartoon",
            method="voxel",
        )

    first, *rest = result.splitlines()
    assert first == "Error rebuilding group 'b': RuntimeError: mesh exploded"
    assert rest == [
        f"The export stopped there. Already written to {tmp_path} before it did: ub_a.stl. "
        "Each is complete for its own group, but the set is not; re-run once the cause "
        "is fixed rather than printing these alone."
    ]


def test_print_export_says_nothing_extra_when_the_first_group_fails(tmp_path):
    pymol = _PyMOL(objects={"pa": "molecule", "pb": "molecule"})
    with (
        patch("mcpymol.printing.send_request", side_effect=pymol),
        patch.dict(sys.modules, {"trimesh": MagicMock(), "manifold3d": MagicMock()}),
        patch("mcpymol.printing._read_obj_triangles", return_value=([0] * 12, 0)),
        patch("mcpymol.printing._repair_to_stl", side_effect=RuntimeError("mesh exploded")),
    ):
        result = print_export(
            obj_name="ub",
            groups="a=pa",
            out_dir=str(tmp_path),
            representation="cartoon",
            method="voxel",
        )

    assert result == "Error rebuilding group 'a': RuntimeError: mesh exploded"
