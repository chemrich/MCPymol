"""3D-print tools against a real PyMOL: the geometry, not the payload.

    pytest -m live tests/test_live_printing.py

WARNING: this clears the PyMOL session it connects to. Save anything you care
about first.

The mocked tests in test_printing.py pin what the tools send. Only PyMOL can
say whether the result is right: whether a setting really stays on its object,
whether a cartoon really runs through a residue with no CA, whether an export
really comes out as one piece. The GFP tests fetch 2b3p and skip when offline.
"""

import socket

import pytest

import mcpymol.bridge as bridge
from mcpymol.bridge import send_request
from mcpymol.printing import _CENSUS_REPS, _census, print_export, print_ribbon_view

pytestmark = pytest.mark.live


def _bridge_reachable() -> bool:
    try:
        with socket.create_connection((bridge.HOST, bridge.PORT), timeout=2):
            return True
    except OSError:
        return False


@pytest.fixture(autouse=True)
def _clean_session():
    """A fresh session per test; skip the module unless a PyMOL plugin is listening."""
    if not _bridge_reachable():
        pytest.skip(
            f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT} — start PyMOL and "
            f"`run src/mcpymol/plugin.py`"
        )
    send_request("reinitialize")
    yield
    send_request("reinitialize")


def _get(setting: str, obj: str = "") -> str:
    res = send_request("get", args=[setting, obj] if obj else [setting])
    assert res.get("status") == "success", res
    return str(res["result"])


def _rows(selection: str, expression: str) -> list:
    res = send_request("iterate_to_list", args=[selection, expression])
    assert res.get("status") == "success", res
    return res["result"]


def _count(selection: str) -> int:
    res = send_request("count_atoms", args=[selection])
    assert res.get("status") == "success", res
    return int(res["result"])


@pytest.fixture
def gfp():
    """GFP (2b3p): one chain whose chromophore, CRO 66, has CA1-3 and no CA."""
    from mcpymol.structures import fetch_structure

    result = fetch_structure(pdb_code="2b3p")
    if "Successfully fetched" not in result:
        pytest.skip(f"could not fetch 2b3p: {result}")
    return "2b3p"


# ── print_export ─────────────────────────────────────────────────────────────
#
# Built with fab (or a few CA atoms written to a file), so these need no
# network; pitches are coarse, so the whole section runs in seconds.


def _pe_call(action: str, *args, **kwargs):
    res = send_request(action, args=list(args), kwargs=kwargs)
    assert res.get("status") == "success", (action, args, res)
    return res.get("result")


def _pe_peptide(name: str = "pep", sequence: str = "ACDEFGHIKLMNPQRSTVWY", **kwargs):
    """A fab peptide showing only its cartoon."""
    _pe_call("fab", sequence, name, **kwargs)
    _pe_call("hide", "everything", name)
    _pe_call("show", "cartoon", name)
    return name


def _pe_xyz(selection: str):
    import numpy as np

    return np.array(_pe_call("get_coords", selection), dtype=float)


def _pe_mesh(path):
    import trimesh

    return trimesh.load(path, force="mesh")


@pytest.fixture
def _pe_mesh_libs():
    """The print extra, which print_export needs to build an STL at all."""
    pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    pytest.importorskip("manifold3d")


def test_census_reads_reps_masks_the_way_pymol_counts_them():
    """_census decodes each atom's `reps` mask with a bit table written out in
    printing.py. Each figure must match PyMOL's own `rep <name>` selection, or
    the table has drifted from cmd.repmasks."""
    for name in ("pep", "lig"):
        assert send_request("fab", args=["ACDEFGHIK", name]).get("status") == "success"
    shows = [
        ("cartoon", "pep and resi 1-5"),
        ("sticks", "pep and resi 4-9"),
        ("spheres", "pep and name CA"),
        ("lines", "lig"),
        ("dots", "lig and resi 2"),
        ("mesh", "pep and resi 9"),
        ("labels", "lig and name CA"),
        ("nb_spheres", "lig"),
        ("ribbon", "lig and resi 1-3"),
    ]
    send_request("hide", args=["everything"])
    for rep, sel in shows:
        assert send_request("show", args=[rep, sel]).get("status") == "success"

    census = _census(["pep", "lig"])

    for obj in ("pep", "lig"):
        for rep in _CENSUS_REPS:
            assert census[obj][rep] == _count(f"({obj}) and rep {rep}"), (obj, rep)
    assert census["pep"]["sticks"] > 0 and census["lig"]["lines"] > 0  # not vacuous


def test_journey_export_frame_ignores_the_camera(_pe_mesh_libs, tmp_path):
    """The STL came out in camera space: turning the view 90 degrees between
    two exports swapped x and z in the split-GFP build. Now both views give
    the same file, in the coordinates get_coords reports, and neither the
    camera nor what is enabled is touched."""
    import numpy as np

    from mcpymol.printing import _SurfaceIndex

    pep = _pe_peptide()
    _pe_call("fab", "GG", "bystander")
    _pe_call("disable", "bystander")
    _pe_call("orient", pep)
    enabled = _pe_call("get_names", "all", 1)
    stls = []
    for turn in (0, 90):
        _pe_call("turn", "y", turn)
        _pe_call("move", "x", 7)
        view = _pe_call("get_view")
        out = tmp_path / f"turn{turn}"
        result = print_export(
            obj_name=pep,
            groups=f"chain={pep}",
            out_dir=str(out),
            representation="cartoon",
            method="voxel",
        )
        assert "Error" not in result, result
        assert _pe_call("get_view") == view
        assert _pe_call("get_names", "all", 1) == enabled
        stls.append(_pe_mesh(out / f"{pep}_chain.stl"))
    a, b = stls
    np.testing.assert_allclose(a.bounds, b.bounds, atol=1e-6)
    np.testing.assert_allclose(np.sort(a.vertices, axis=0), np.sort(b.vertices, axis=0), atol=1e-6)
    # In model coordinates: every Cα is on the tube (camera space put them
    # tens of Å away).
    distance = _SurfaceIndex(a, max_edge=0.8).signed_distance(_pe_xyz(f"{pep} and name CA"))
    assert (abs(distance) < 1.5).all(), distance


def test_journey_surface_export_holds_every_atom(_pe_mesh_libs, tmp_path):
    """Under a rotated, shifted view the surface used to come out rotated
    with it; every atom must lie inside the exported surface."""
    from mcpymol.printing import _SurfaceIndex

    pep = _pe_peptide(sequence="ACDEFGHIK")
    _pe_call("turn", "x", 40)
    _pe_call("move", "y", 5)
    view = _pe_call("get_view")

    result = print_export(obj_name=pep, groups="all=all", out_dir=str(tmp_path))

    assert "Error" not in result, result
    assert _pe_call("get_view") == view
    mesh = _pe_mesh(tmp_path / f"{pep}_all.stl")
    inside = _SurfaceIndex(mesh, max_edge=0.8).signed_distance(_pe_xyz(pep))
    assert (inside < 0).all(), inside.max()


def test_journey_sticks_print_as_one_watertight_piece(_pe_mesh_libs, tmp_path):
    """A sticks-only group exported no triangles and crashed with 'zero-size
    array to reduction operation minimum which has no identity'."""
    import numpy as np

    lig = "lig"
    _pe_call("fab", "WYF", lig)
    # Hydrogens are drawn at stick_h_scale (0.4 x the radius): leave them out
    # so every atom has the same stick radius.
    _pe_call("remove", f"{lig} and hydro")
    _pe_call("hide", "everything", lig)
    _pe_call("show", "sticks", lig)
    _pe_call("turn", "z", 30)

    result = print_export(
        obj_name=lig, groups=f"sticks={lig}", out_dir=str(tmp_path), representation="displayed"
    )

    assert "exact union of" in result and "one piece" in result, result
    mesh = _pe_mesh(tmp_path / f"{lig}_sticks.stl")
    assert mesh.is_watertight
    assert len(mesh.split(only_watertight=False)) == 1
    # Within a stick radius of the atoms (circumscribed: never thinner).
    radius = float(_pe_call("get", "stick_radius"))
    xyz = _pe_xyz(lig)
    np.testing.assert_allclose(
        mesh.bounds, [xyz.min(axis=0) - radius, xyz.max(axis=0) + radius], atol=0.05
    )


def test_journey_separate_chains_print_as_separate_pieces(_pe_mesh_libs, tmp_path):
    """Keeping the body with the most faces dropped one of 4hhb's two dimers
    and still said OK. Both pieces stay, or bodies='largest' says what it
    dropped, and the residue check notices the chain that is gone."""
    # fab replaces an object it is given, so build the chains apart and join them.
    _pe_call("fab", "ACDEFGHIKLM", "_pa", chain="A")
    _pe_call("fab", "ACDEFGHIKLM", "_pb", chain="B")
    _pe_call("translate", [0.0, 0.0, 20.0], "_pb", camera=0)
    _pe_call("create", "pair", "_pa or _pb")
    _pe_call("delete", "_pa or _pb")
    _pe_call("hide", "everything", "pair")
    _pe_call("show", "cartoon", "pair")

    both = print_export(
        obj_name="pair",
        groups="both=pair",
        out_dir=str(tmp_path / "all"),
        representation="cartoon",
        method="voxel",
    )
    one = print_export(
        obj_name="pair",
        groups="both=pair",
        out_dir=str(tmp_path / "largest"),
        representation="cartoon",
        method="voxel",
        bodies="largest",
    )

    assert "prints as 2 separate pieces" in both, both
    assert len(_pe_mesh(tmp_path / "all" / "pair_both.stl").split(only_watertight=False)) == 2
    assert "WARNING: kept the largest piece and dropped 1 more" in one, one
    assert "residue(s) PyMOL shows as cartoon are missing from the mesh" in one, one
    assert len(_pe_mesh(tmp_path / "largest" / "pair_both.stl").split(only_watertight=False)) == 1


def test_journey_scale_multiplies_about_the_origin(_pe_mesh_libs, tmp_path):
    """mm_per_angstrom scales every group about the origin, so parts from
    separate calls at one scale still line up."""
    import numpy as np

    pep = _pe_peptide()
    groups = f"head={pep} and resi 1-10; tail={pep} and resi 11-20"
    sizes = {}
    for scale in (1.0, 2.0):
        result = print_export(
            obj_name=pep,
            groups=groups,
            out_dir=str(tmp_path / str(scale)),
            mm_per_angstrom=scale,
        )
        assert f"at {scale:g} mm per Å" in result, result
        sizes[scale] = [
            _pe_mesh(tmp_path / str(scale) / f"{pep}_{g}.stl") for g in ("head", "tail")
        ]
    for small, large in zip(sizes[1.0], sizes[2.0], strict=True):
        np.testing.assert_allclose(large.bounds, 2 * small.bounds, atol=1e-3)
        assert large.volume == pytest.approx(8 * small.volume, rel=1e-4)


def test_journey_distance_prints_as_a_solid_rod(_pe_mesh_libs, tmp_path):
    """Distance dashes printed as loose pellets, or not at all; named in the
    group, a distance object made PyMOL reject the whole selection."""
    import numpy as np

    from mcpymol.printing import _SurfaceIndex

    _pe_call("fab", "WWW", "lig")
    _pe_call("hide", "everything", "lig")
    _pe_call("show", "sticks", "lig")
    _pe_call("distance", "d1", "lig and resi 1 and name CA", "lig and resi 3 and name CA")
    _pe_call("set", "dash_radius", 0.2, "d1")

    result = print_export(
        obj_name="lig",
        groups="lig=(lig or d1)",
        out_dir=str(tmp_path),
        representation="displayed",
    )

    assert "the distance dashes of d1 were exported as solid rods" in result, result
    assert "one piece" in result, result
    assert float(_pe_call("get", "dash_gap", "d1")) > 0
    a, b = _pe_xyz("lig and resi 1 and name CA")[0], _pe_xyz("lig and resi 3 and name CA")[0]
    along = a + np.outer(np.linspace(0.2, 0.8, 7), b - a)
    mesh = _pe_mesh(tmp_path / "lig_lig.stl")
    assert (_SurfaceIndex(mesh, max_edge=0.8).signed_distance(along) < 0).all()


def test_journey_pixel_sized_lines_export_at_their_drawn_size(_pe_mesh_libs, tmp_path):
    """With the camera at the origin PyMOL wrote every line and dash it sizes
    in pixels with radius NaN."""
    pep = _pe_peptide(sequence="ACDEF")
    _pe_call("hide", "everything", pep)
    _pe_call("show", "lines", pep)

    result = print_export(
        obj_name=pep, groups=f"lines={pep}", out_dir=str(tmp_path), representation="displayed"
    )

    assert "exact union of" in result, result
    assert "lines and nonbonded crosses, sized in pixels" in result


def test_journey_an_empty_group_is_explained(_pe_mesh_libs, tmp_path):
    _pe_call("fab", "WYF", "lig")
    _pe_call("hide", "everything", "lig")
    _pe_call("label", "lig and name CA", "resn")

    result = print_export(
        obj_name="lig", groups="lig=lig", out_dir=str(tmp_path), representation="displayed"
    )

    assert result.startswith("Error: group 'lig' (lig) exported no printable geometry."), result
    assert "labels (3 atoms)" in result and "labels have no 3D geometry" in result


def test_journey_cartoon_mode_names_what_it_left_out(_pe_mesh_libs, tmp_path):
    pep = _pe_peptide()
    _pe_call("show", "sticks", f"{pep} and resi 5-7 and not name N+C+O")

    result = print_export(
        obj_name=pep,
        groups=f"pep={pep}",
        out_dir=str(tmp_path),
        representation="cartoon",
        method="voxel",
    )

    assert f"NOTE: left out {pep}: sticks" in result, result
    assert "representation='displayed' includes them" in result


def test_journey_rock_and_groups_survive_an_export(_pe_mesh_libs, tmp_path):
    """Rock turned the camera mid-export; an object inside a disabled PyMOL
    group exported nothing; and the group came back enabled afterwards."""
    pep = _pe_peptide(sequence="ACDEFGHIK")
    _pe_call("group", "grp", pep)
    _pe_call("disable", "grp")
    _pe_call("rock", 1)

    result = print_export(
        obj_name=pep,
        groups=f"pep={pep}",
        out_dir=str(tmp_path),
        representation="cartoon",
        method="voxel",
    )

    assert "Error" not in result, result
    assert _pe_call("get", "rock") == "on"
    assert "grp" not in _pe_call("get_names", "all", 1)
    _pe_call("rock", 0)


def test_journey_voxel_growth_is_undone_on_a_rod(_pe_mesh_libs, tmp_path):
    """A tube PyMOL draws 0.48 Å across on average came out 0.61 Å at 0.2 Å
    pitch: voxelising grows every surface by about 0.62 voxel."""
    import numpy as np

    # Seven Cα along an oblique line off the voxel grid, traced as a tube.
    start = np.array([3.13, -2.71, 5.37])
    axis = np.array([1.0, 0.7, 0.3]) / np.linalg.norm([1.0, 0.7, 0.3])
    lines = [
        f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        for i, (x, y, z) in enumerate(start + np.outer(np.arange(7) * 3.8, axis))
    ]
    pdb = tmp_path / "rod.pdb"
    pdb.write_text("\n".join(lines) + "\nEND\n")
    _pe_call("load", str(pdb), "rod")
    _pe_call("hide", "everything", "rod")
    _pe_call("set", "cartoon_trace_atoms", 1)
    _pe_call("cartoon", "tube", "rod")
    _pe_call("set", "cartoon_tube_radius", 0.5)
    _pe_call("show", "cartoon", "rod")

    radii = {}
    for undo in (True, False):
        out = tmp_path / str(undo)
        result = print_export(
            obj_name="rod",
            groups="rod=rod",
            out_dir=str(out),
            representation="cartoon",
            method="voxel",
            voxel_pitch=0.2,
            undo_voxel_growth=undo,
        )
        assert "Error" not in result, result
        points = _pe_mesh(out / "rod_rod.stl").sample(20000, seed=0) - start
        along = points @ axis
        middle = points[(along > 7) & (along < 16)]
        across = middle - np.outer(middle @ axis, axis)
        radii[undo] = float(np.linalg.norm(across, axis=1).mean())
    # PyMOL's own polygonal tube averages 0.48 Å from its axis.
    assert radii[True] == pytest.approx(0.48, abs=0.03), radii
    assert radii[False] > radii[True] + 0.08, radii


def test_journey_thin_cartoon_is_warned_about(_pe_mesh_libs, tmp_path):
    """A 0.66 Å rung was a single voxel at the default 0.7 Å pitch, and nothing said so."""
    pep = _pe_peptide(ss=1)

    coarse = print_export(
        obj_name=pep,
        groups=f"pep={pep}",
        out_dir=str(tmp_path),
        representation="cartoon",
        method="voxel",
    )

    assert "WARNING: the thinnest feature PyMOL draws here is" in coarse, coarse
    assert f"`set cartoon_oval_width, <value>, {pep}`" in coarse or "cartoon_loop_radius" in coarse


# ── print_ribbon_view ────────────────────────────────────────────────────────


def test_journey_ribbon_style_stays_on_its_objects():
    """The style used to be set globally: every other object in the session,
    and every one loaded later, took the chunky strands too."""
    send_request("fab", args=["ACDEFGHIKLMNPQRSTVWY", "_ribbon_a"])
    send_request("create", args=["_ribbon_b", "_ribbon_a"])
    default = _get("cartoon_rect_width")

    result = print_ribbon_view(obj_name="_ribbon_a")

    assert not result.startswith("Error"), result
    assert _get("cartoon_rect_width") == default
    assert _get("cartoon_rect_width", "_ribbon_a") == "2.00000"
    assert _get("cartoon_rect_width", "_ribbon_a_spine") == "2.00000"
    assert _get("cartoon_rect_width", "_ribbon_b") == default


def test_journey_ribbon_reports_a_missing_object():
    """It used to return its full success text, and create nothing."""
    result = print_ribbon_view(obj_name="_no_such_object")

    assert result.startswith("Error:"), result
    names = send_request("get_names", args=["objects"]).get("result") or []
    assert "_no_such_object_spine" not in names


def test_journey_ribbon_bridges_a_chromophore_beside_a_calcium_ion(gfp):
    """A calcium ion is named CA too, and PyMOL leaves it unbonded: counting
    it as a traced Cα refused the bridge for every calcium sensor (GCaMP2,
    3EVR). PyMOL only takes a carbon as a CA guide atom."""
    ion = f"first ({gfp} and resn CD)"
    assert _count(ion) == 1, "2b3p's cadmium ions are the stand-in for calcium here"
    _pe_call("alter", ion, "name='CA'; resn='CA'; elem='CA'")

    result = print_ribbon_view(obj_name=gfp)

    assert "Spine bridged through 1 residue with no CA atom: CRO A66" in result, result


def test_journey_ribbon_spine_bridges_the_gfp_chromophore(gfp, tmp_path):
    """CRO has CA1-3 but no CA, so the spine broke at 64/68: three pieces at
    a thin radius. The bridge traces it on the spine only and renames nothing."""
    result = print_ribbon_view(obj_name=gfp, spine_radius=0.3)

    assert "Spine bridged through 1 residue with no CA atom: CRO A66" in result, result
    spine = f"{gfp}_spine"
    traced = [row[:3] for row in _rows(spine, "resn, resi, name, s.cartoon_trace_atoms") if row[3]]
    assert traced == [
        ["CRO", "66", "CA1"],
        ["CRO", "66", "CA2"],
        ["CRO", "66", "CA3"],
        ["VAL", "68", "CA"],
    ]
    assert _get("trace_atoms_mode", spine) == "0"
    assert _get("trace_atoms_mode", gfp) == "5"
    # The Cα either side is held against loop smoothing, on the spine only.
    assert _rows(f"{spine} and flag 26", "resn, resi, name") == [
        ["LEU", "64", "CA"],
        ["VAL", "68", "CA"],
    ]
    assert _count(f"{gfp} and flag 26") == 0
    # Nothing renamed or renumbered, in either object.
    assert _count(f"({gfp} or {spine}) and resn CRO and name CA") == 0
    assert _count(f"{gfp} and resn CRO and resi 66") == _count(f"{spine} and resn CRO and resi 66")

    pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    bridged = print_export(
        obj_name=gfp,
        groups=f"spine={spine}",
        out_dir=str(tmp_path / "bridged"),
        representation="cartoon",
        method="voxel",
        voxel_pitch=0.15,
    )
    assert "one piece" in bridged, bridged
    # Without the bridge the same spine comes apart: the test can tell.
    print_ribbon_view(obj_name=gfp, spine_radius=0.3, bridge_modified=False)
    broken = print_export(
        obj_name=gfp,
        groups=f"spine={spine}",
        out_dir=str(tmp_path / "broken"),
        representation="cartoon",
        method="voxel",
        voxel_pitch=0.15,
    )
    assert "separate pieces" in broken, broken


def test_journey_ribbon_spine_leaves_strand_interiors(gfp, tmp_path):
    """The open-barrel recipe: thinner strands, no spine inside them, and the
    barrel still prints as one piece with the spine through Leu64's Cα."""
    result = print_ribbon_view(
        obj_name=gfp, strand_width=2.0, strand_thickness=1.0, spine_in_strands=False
    )

    assert "interior of 11 strands" in result, result
    spine = f"{gfp}_spine"
    # Strand 25-36 keeps two residues of spine at each end.
    assert _count(f"{spine} and name CA and rep cartoon and resi 27-34") == 0
    assert _count(f"{spine} and name CA and rep cartoon and resi 25+26+35+36") == 4
    # A bridge never sits next to a hidden residue.
    assert _count(f"{spine} and name CA and rep cartoon and resi 62-64+68-70") == 6
    assert _get("cartoon_rect_width", gfp) == "1.00000"

    pytest.importorskip("trimesh")
    pytest.importorskip("skimage")
    exported = print_export(
        obj_name=gfp,
        groups=f"{gfp}=({gfp} or {spine})",
        out_dir=str(tmp_path),
        representation="cartoon",
        method="voxel",
        voxel_pitch=0.25,
    )
    assert "one piece" in exported, exported
    # Loop smoothing used to pull the tube 1.5 Å off Leu64's Cα next to the
    # chromophore, which print_export reports as a gap.
    assert "LEU64" not in exported, exported
