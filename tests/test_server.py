import inspect
import json
import subprocess
import sys
from unittest.mock import MagicMock, patch

import pytest

from mcpymol.server import (
    _DEFAULT_TIMEOUT,
    _SLOW_OP_TIMEOUT,
    DEFAULT_MULTIMER_CUTOFF,
    _apply_multimer_heuristic,
    _call,
    _parse_groups,
    _repair_to_stl,
    cinematic_view,
    color,
    crosslink_view,
    draw,
    electrostatic_view,
    fetch_structure,
    hydrophobic_surface_view,
    interface_view,
    ligand_view,
    list_chains,
    list_ligands,
    list_objects,
    load_structure,
    mpng,
    mutation_view,
    pharmacophore_view,
    png,
    pocket_view,
    pointillist_view,
    poisson_boltzmann_view,
    print_export,
    print_ribbon_view,
    putty_view,
    ray,
    save,
    select,
    send_request,
    show,
    textbook_view,
)

# ── Helpers ───────────────────────────────────────────────────────────────────


def _sr_mock(**action_results):
    """Return a side_effect for patching send_request.

    Calls whose action key is in action_results get that result value; all
    other calls return a generic success response.  If a value is callable
    it is invoked with (action, args, kwargs) to produce the result.
    """

    def fake(action, args=None, kwargs=None, **_ignored):
        if action in action_results:
            val = action_results[action]
            result = val(action, args, kwargs) if callable(val) else val
            return {"status": "success", "result": result}
        # Loaders now verify that atoms actually arrived before doing anything
        # destructive, so the default has to look like a structure that loaded.
        if action == "count_atoms":
            return {"status": "success", "result": 1200}
        if action == "get_object_list":
            return {"status": "success", "result": []}
        return {"status": "success", "result": "OK"}

    return fake


def _actions(mock_sr):
    """Return the ordered list of action names sent via a patched send_request."""
    return [c.args[0] for c in mock_sr.call_args_list]


# ── socket.socket fixture ─────────────────────────────────────────────────────


@pytest.fixture
def mock_socket():
    """Mock the socket context manager; default recv returns a success payload."""
    with patch("socket.socket") as mock_cls:
        inst = MagicMock()
        mock_cls.return_value.__enter__.return_value = inst
        inst.recv.return_value = json.dumps(
            {"status": "success", "result": "Mocked execution"}
        ).encode()
        yield inst


# ── send_request unit tests ───────────────────────────────────────────────────


def test_send_request_success(mock_socket):
    """Low-level wire format and return value."""
    res = send_request("test_method", args=["arg1"], kwargs={"kw1": "val1"})
    assert res == {"status": "success", "result": "Mocked execution"}

    mock_socket.sendall.assert_called_once()
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "test_method"
    assert payload["args"] == ["arg1"]
    assert payload["kwargs"] == {"kw1": "val1"}


def test_send_request_connection_refused():
    """ConnectionRefusedError returns a structured error dict."""
    with patch("socket.socket") as mock_cls:
        mock_cls.return_value.__enter__.return_value.connect.side_effect = ConnectionRefusedError(
            "Connection refused"
        )
        res = send_request("test")
    assert res["status"] == "error"
    assert "Socket connection failed" in res["error"]
    assert "Connection refused" in res["error"]


def test_send_request_timeout():
    """Socket timeout returns a structured error dict."""
    with patch("socket.socket") as mock_cls:
        mock_cls.return_value.__enter__.return_value.connect.side_effect = TimeoutError("Timed out")
        res = send_request("test")
    assert res["status"] == "error"
    assert "Timed out" in res["error"]


# ── Primitive tool wrappers ───────────────────────────────────────────────────


def test_tool_show(mock_socket):
    result = show(representation="cartoon", selection="chain A")
    assert "Showing cartoon" in result
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "show"
    assert payload["args"] == ["cartoon", "chain A"]


def test_tool_color(mock_socket):
    result = color(color_name="red", selection="all")
    assert "Colored selection 'all' with red" in result
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "color"
    assert payload["args"] == ["red", "all"]


def test_tool_select(mock_socket):
    result = select(name="my_selection", selection="resi 1-10")
    assert "Created named selection" in result
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "select"
    assert payload["args"] == ["my_selection", "resi 1-10"]


def test_tool_error_propagation(mock_socket):
    """Plugin errors surface back to the MCP caller."""
    mock_socket.recv.return_value = json.dumps(
        {
            "status": "error",
            "error": "PyMOL encountered a problem: invalid selection",
        }
    ).encode()
    result = show(representation="spheres")
    assert "invalid selection" in result
    assert result.startswith("PyMOL encountered")


# ── fetch_structure ───────────────────────────────────────────────────────────


@patch("mcpymol.structures.send_request")
def test_fetch_structure_multimer(mock_sr):
    """Chains A/B/C → remove non-proximal chains and apply ghost heart."""
    mock_sr.side_effect = _sr_mock(get_chains=["A", "B", "C"])

    result = fetch_structure(pdb_code="1ubq")

    assert "Successfully fetched 1ubq" in result
    acts = _actions(mock_sr)
    # A scoped delete of the object being replaced — nothing session-wide runs
    # before the fetch is known to have produced atoms.
    assert acts[0] == "delete"
    assert "fetch" in acts
    assert "get_chains" in acts
    assert "remove" in acts  # multimer cleanup
    assert "hide" in acts  # solvent hidden


@patch("mcpymol.structures.send_request")
def test_fetch_structure_single_chain(mock_sr):
    """Single chain still applies the keep-selection + ghost heart."""
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    result = fetch_structure(pdb_code="1ubq")

    assert "Successfully fetched 1ubq" in result
    acts = _actions(mock_sr)
    assert "remove" in acts
    assert "hide" in acts


@patch("mcpymol.structures.send_request")
def test_fetch_structure_no_chains(mock_sr):
    """Empty chain list falls through to the fallback return message."""
    mock_sr.side_effect = _sr_mock(get_chains=[])

    result = fetch_structure(pdb_code="1abc")

    assert "Successfully fetched 1abc" in result
    acts = _actions(mock_sr)
    assert "fetch" in acts
    assert "remove" not in acts  # no cleanup attempted


@patch("mcpymol.structures.send_request")
def test_fetch_structure_error(mock_sr):
    """Fetch failure is propagated as a readable error string."""

    def fake(action, args=None, kwargs=None):
        if action == "fetch":
            return {"status": "error", "error": "PDB ID not found"}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = fetch_structure(pdb_code="XXXX")

    assert "Error fetching XXXX" in result
    assert "PDB ID not found" in result


@patch("mcpymol.structures.send_request")
def test_fetch_structure_custom_obj_name(mock_sr):
    """Custom obj_name is threaded through all downstream send_request calls."""
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    fetch_structure(pdb_code="1ubq", obj_name="my_protein")

    for c in mock_sr.call_args_list:
        if c.args[0] == "fetch":
            assert c.kwargs["args"] == ["1ubq", "my_protein"]
            break
    else:
        pytest.fail("fetch action was never called")


# ── load_structure ────────────────────────────────────────────────────────────


@patch("mcpymol.structures.send_request")
def test_load_structure_success(mock_sr):
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    result = load_structure(file_path="/tmp/test.pdb", obj_name="test")

    assert "Successfully loaded" in result
    acts = _actions(mock_sr)
    assert "load" in acts
    assert "remove" in acts
    assert "hide" in acts


@patch("mcpymol.structures.send_request")
def test_load_structure_error(mock_sr):
    def fake(action, args=None, kwargs=None):
        if action == "load":
            return {"status": "error", "error": "File not found"}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = load_structure(file_path="/bad/path.pdb", obj_name="test")

    assert "Error loading" in result
    assert "File not found" in result


# ── View functions ────────────────────────────────────────────────────────────


@patch("mcpymol.views.send_request")
def test_ligand_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = ligand_view(obj_name="1ATP", ligand_resn="ATP")

    assert "ATP" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    assert "zoom" in acts
    # H-bonds drawn via a 'do distance ...' call
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("distance" in a for a in do_args)


@patch("mcpymol.views.send_request")
def test_interface_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = interface_view(obj_name="1BRS", chain_a="A", chain_b="D")

    assert "chain A" in result and "chain D" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("distance" in a for a in do_args)


@patch("mcpymol.views.send_request")
def test_putty_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = putty_view(obj_name="1UBQ")

    assert "Putty view" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts
    assert "set" in acts


@patch("mcpymol.printing.send_request")
def test_print_ribbon_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = print_ribbon_view(obj_name="1UBQ", spine_radius=1.1)

    assert "Print-ribbon view" in result
    # Tells the user how to export the fused solid, in cartoon mode.
    assert "1UBQ_spine" in result
    assert "print_export" in result
    assert 'representation="cartoon"' in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "create" in acts
    assert "set" in acts
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("cartoon tube, 1UBQ_spine" in a for a in do_args)


@patch("mcpymol.views.send_request")
def test_hydrophobic_surface_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = hydrophobic_surface_view(obj_name="1TCA")

    assert "Hydrophobic" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts


@patch("mcpymol.views.send_request")
def test_electrostatic_view_atomic(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = electrostatic_view(obj_name="1LYZ")

    assert "Electrostatic" in result
    assert "atomic" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts


@patch("mcpymol.views.send_request")
def test_electrostatic_view_residue(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = electrostatic_view(obj_name="1LYZ", mode="residue")

    assert "residue" in result
    assert "Electrostatic" in result


@patch("mcpymol.views.send_request")
def test_crosslink_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = crosslink_view(obj_name="1CEL")

    assert "Crosslink" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("distance" in a for a in do_args)


@patch("mcpymol.views.send_request")
def test_pocket_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = pocket_view(obj_name="1HSG", resn="MK1")

    assert "MK1" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    assert "zoom" in acts


@patch("mcpymol.views.send_request")
def test_pharmacophore_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = pharmacophore_view(obj_name="1HSG", resn="MK1")

    assert "Pharmacophore" in result
    assert "MK1" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "zoom" in acts


@patch("mcpymol.views.send_request")
def test_mutation_view_valid(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = mutation_view(obj_name="4HHB", mutations="E6V,K16E")

    assert "Mutation view" in result
    assert "E6V" in result or "K16E" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    assert "zoom" in acts


@patch("mcpymol.views.send_request")
def test_textbook_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = textbook_view(obj_name="1ABC")

    assert "Textbook Illustration" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "set" in acts


@patch("mcpymol.views.send_request")
def test_cinematic_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = cinematic_view(obj_name="1ABC")

    assert "Cinematic view" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "set" in acts
    # The background moved into style.set_background, so it is no longer
    # visible on this module's send_request; test_rendering covers it.


@patch("mcpymol.views.send_request")
def test_pointillist_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = pointillist_view(obj_name="1ABC")

    assert "Pointillist/Starfield" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts
    assert "set" in acts
    assert "color" in acts


@patch("mcpymol.views.send_request")
def test_mutation_view_invalid_input(mock_sr):
    """Non-parseable mutation strings return an error before any PyMOL calls."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = mutation_view(obj_name="4HHB", mutations="bad_input_no_digits")

    assert "No valid mutations" in result
    assert mock_sr.call_count == 0


# ── Scene introspection ───────────────────────────────────────────────────────


@patch("mcpymol.structures.send_request")
def test_list_objects_with_objects(mock_sr):
    mock_sr.return_value = {"status": "success", "result": ["1ubq", "lig"]}
    result = list_objects()
    assert "1ubq" in result and "lig" in result
    assert mock_sr.call_args.args[0] == "get_object_list"


@patch("mcpymol.structures.send_request")
def test_list_objects_empty(mock_sr):
    mock_sr.return_value = {"status": "success", "result": []}
    assert list_objects() == "No objects are loaded."


@patch("mcpymol.structures.send_request")
def test_list_chains(mock_sr):
    mock_sr.return_value = {"status": "success", "result": ["A", "B", "C"]}
    result = list_chains("1abc")
    assert "A" in result and "B" in result and "C" in result
    assert mock_sr.call_args.args[0] == "get_chains"
    assert mock_sr.call_args.kwargs["args"] == ["1abc"]


@patch("mcpymol.structures.send_request")
def test_list_chains_empty(mock_sr):
    mock_sr.return_value = {"status": "success", "result": []}
    assert "No chains" in list_chains("empty")


@patch("mcpymol.structures.send_request")
def test_list_ligands_parses_pdb(mock_sr):
    """list_ligands parses HETATM resn columns out of the PDB dump."""
    pdb = (
        "HETATM    1  C1  ATP A 200      11.111  22.222  33.333  1.00 20.00           C\n"
        "HETATM    2  N1  ATP A 200      11.111  22.222  33.333  1.00 20.00           N\n"
        "HETATM    3  MG  MG  A 201      14.000  25.000  36.000  1.00 20.00          MG\n"
        "END\n"
    )
    mock_sr.return_value = {"status": "success", "result": pdb}
    result = list_ligands("1atp")
    assert "ATP" in result and "MG" in result
    assert "1atp" in result


@patch("mcpymol.structures.send_request")
def test_list_ligands_no_organic(mock_sr):
    mock_sr.return_value = {"status": "success", "result": ""}
    assert "No organic ligands" in list_ligands("1ubq")


# ── print_export ──────────────────────────────────────────────────────────────


def test_parse_groups_valid():
    pairs = _parse_groups("protein=polymer.protein; nucleic=chain N+R+T")
    assert pairs == [("protein", "polymer.protein"), ("nucleic", "chain N+R+T")]


def test_parse_groups_trailing_semicolon_and_spaces():
    assert _parse_groups("  a = chain A ; ") == [("a", "chain A")]


@pytest.mark.parametrize("bad", ["", "noequals", "label=", "=selection"])
def test_parse_groups_invalid(bad):
    with pytest.raises(ValueError):
        _parse_groups(bad)


@patch("mcpymol.printing.send_request")
def test_print_export_missing_deps(mock_sr):
    """Without the optional 'print' extra, returns an install hint, no PyMOL calls."""
    with patch.dict(sys.modules, {"trimesh": None}):
        result = print_export(obj_name="1MSW", groups="protein=polymer.protein")
    assert "uv sync --extra print" in result
    assert mock_sr.call_count == 0


@patch("mcpymol.printing._repair_to_stl")
@patch("mcpymol.printing.send_request")
def test_print_export_happy_path(mock_sr, mock_repair, tmp_path):
    mock_sr.return_value = {"status": "success", "result": "OK"}
    mock_repair.return_value = {"method": "poisson", "faces": 144290, "watertight": True}

    with patch.dict(sys.modules, {"trimesh": MagicMock()}):
        result = print_export(
            obj_name="1MSW",
            groups="protein=polymer.protein; nucleic=polymer.nucleic",
            out_dir=str(tmp_path),
        )

    assert "1MSW_protein.stl" in result
    assert "1MSW_nucleic.stl" in result
    assert "OK" in result  # watertight flag
    acts = _actions(mock_sr)
    # Each group is isolated (create) and exported (save), then cleaned up.
    assert acts.count("create") == 2
    assert acts.count("save") == 2
    assert acts.count("delete") >= 2
    assert mock_repair.call_count == 2


@patch("mcpymol.printing._repair_to_stl")
@patch("mcpymol.printing.send_request")
def test_print_export_surface_mode_shows_surface(mock_sr, mock_repair, tmp_path):
    """Regression: default surface mode must actually show a surface."""
    mock_sr.return_value = {"status": "success", "result": "OK"}
    mock_repair.return_value = {"method": "voxel", "faces": 1000, "watertight": True}

    with patch.dict(sys.modules, {"trimesh": MagicMock()}):
        print_export(obj_name="1MSW", groups="protein=polymer.protein", out_dir=str(tmp_path))

    show_args = [c.kwargs["args"] for c in mock_sr.call_args_list if c.args[0] == "show"]
    assert any(a[0] == "surface" for a in show_args)


@patch("mcpymol.printing._repair_to_stl")
@patch("mcpymol.printing.send_request")
def test_print_export_cartoon_mode(mock_sr, mock_repair, tmp_path):
    """Cartoon mode exports the displayed cartoon of the real objects:

    no temp object (create), no forced surface; isolate via get_object_list
    + enable, then save.
    """

    def _sr(action, args=None, **kw):
        if action == "get_object_list":
            return {"status": "success", "result": ["1ema", "1ema_spine"]}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = _sr
    mock_repair.return_value = {"method": "voxel", "faces": 634188, "watertight": True}

    with patch.dict(sys.modules, {"trimesh": MagicMock()}):
        result = print_export(
            obj_name="1ema",
            groups="1ema=(1ema or 1ema_spine)",
            out_dir=str(tmp_path),
            representation="cartoon",
            method="voxel",
            voxel_pitch=0.2,
        )

    assert "1ema_1ema.stl" in result
    assert "OK" in result
    acts = _actions(mock_sr)
    assert "get_object_list" in acts
    assert "save" in acts
    assert "enable" in acts
    # The bug being fixed: cartoon mode must NOT force a surface or build
    # a throwaway temp object.
    assert "create" not in acts
    show_args = [c.kwargs["args"] for c in mock_sr.call_args_list if c.args[0] == "show"]
    assert not any(a and a[0] == "surface" for a in show_args)
    mock_repair.assert_called_once()


@patch("mcpymol.printing.send_request")
def test_print_export_bad_representation(mock_sr):
    with patch.dict(sys.modules, {"trimesh": MagicMock()}):
        result = print_export(obj_name="1ema", groups="x=1ema", representation="ribbon")
    assert result.startswith("Error:")
    assert "surface" in result and "cartoon" in result


@patch("mcpymol.printing.send_request")
def test_print_export_bad_group(mock_sr):
    with patch.dict(sys.modules, {"trimesh": MagicMock()}):
        result = print_export(obj_name="1MSW", groups="garbage")
    assert result.startswith("Error:")
    assert mock_sr.call_count == 0


@patch("mcpymol.printing._repair_to_stl")
@patch("mcpymol.printing.send_request")
def test_print_export_save_error_reported(mock_sr, mock_repair, tmp_path):
    mock_sr.return_value = {"status": "error", "error": "disk full"}
    with patch.dict(sys.modules, {"trimesh": MagicMock()}):
        result = print_export(
            obj_name="1MSW", groups="protein=polymer.protein", out_dir=str(tmp_path)
        )
    assert "disk full" in result
    mock_repair.assert_not_called()


def test_repair_auto_light_path_when_already_watertight():
    """auto: an already-watertight export uses light cleanup, not Poisson.

    Regression for the sfGFP finding — Poisson degraded an already-closed
    compact barrel surface, so auto must skip it when not needed.
    """
    fake_tm = MagicMock()
    shell = MagicMock()
    shell.faces = list(range(2000))
    shell.is_watertight = True
    raw = MagicMock()
    raw.is_watertight = True
    raw.split.return_value = [shell]
    fake_tm.load.return_value = raw

    # pymeshlab set to None so any Poisson attempt would raise ImportError.
    with patch.dict(sys.modules, {"trimesh": fake_tm, "pymeshlab": None}):
        info = _repair_to_stl("in.obj", "out.stl", "auto", voxel_pitch=0.7, poisson_depth=10)

    assert info["method"] == "light (already watertight)"
    assert info["watertight"] is True
    assert info["faces"] == 2000
    shell.export.assert_called_once_with("out.stl")
    fake_tm.repair.fix_normals.assert_called_once()


def test_repair_voxel_consolidates_to_single_body():
    """voxel: fragmented marching-cubes output is reduced to the largest
    watertight body and hole-filled (regression for the cartoon export
    coming out as 19 loose, non-watertight shells)."""
    fake_tm = MagicMock()
    m = MagicMock()
    fake_tm.load.return_value = m
    vox = MagicMock()
    m.voxelized.return_value.fill.return_value = vox

    out = MagicMock()
    vox.marching_cubes = out
    small = MagicMock()
    small.faces = list(range(10))
    big = MagicMock()
    big.faces = list(range(900))
    big.is_watertight = True
    out.split.return_value = [small, big]

    with patch.dict(sys.modules, {"trimesh": fake_tm}):
        info = _repair_to_stl("in.obj", "out.stl", "voxel", voxel_pitch=0.2, poisson_depth=10)

    # Largest body kept, closed, and exported — not the raw soup.
    assert info["method"] == "voxel"
    assert info["watertight"] is True
    assert info["faces"] == 900
    big.export.assert_called_once_with("out.stl")
    fake_tm.repair.fill_holes.assert_called_once_with(big)


# ── poisson_boltzmann_view: external-tool failure handling ────────────────────


def _pb_send_request(save_status="success"):
    """Fake send_request that actually writes the PDB the solver will read."""

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "do" and args and args[0].startswith("save "):
            if save_status == "error":
                return {"status": "error", "error": "session is busy"}
            path = args[0].split("save ", 1)[1].split(",")[0].strip()
            with open(path, "w") as fh:
                fh.write("ATOM      1  N   MET A   1       0.0   0.0   0.0\n")
        return {"status": "success", "result": "OK"}

    return fake


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_reports_save_failure(mock_sr, mock_run):
    """A failed save must be reported, not left to surface as a PDB2PQR error
    about a file that was never written."""
    mock_sr.side_effect = _pb_send_request(save_status="error")

    result = poisson_boltzmann_view(obj_name="1LYZ")

    assert "Error saving 1LYZ" in result
    assert "session is busy" in result
    mock_run.assert_not_called()


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_reports_missing_pdb(mock_sr, mock_run):
    """PyMOL claiming success without producing the file means the two halves
    aren't sharing a filesystem — say so instead of running the solver."""
    mock_sr.side_effect = _sr_mock()  # never writes anything

    result = poisson_boltzmann_view(obj_name="1LYZ")

    assert "did not write" in result
    assert "same machine" in result
    mock_run.assert_not_called()


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_reports_missing_pdb2pqr(mock_sr, mock_run):
    mock_sr.side_effect = _pb_send_request()
    mock_run.side_effect = FileNotFoundError("pdb2pqr")

    result = poisson_boltzmann_view(obj_name="1LYZ")

    assert "PDB2PQR not found on PATH" in result
    assert "pip install pdb2pqr" in result


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_reports_pdb2pqr_timeout(mock_sr, mock_run):
    """A wedged solver must not hang the MCP server forever."""
    mock_sr.side_effect = _pb_send_request()
    mock_run.side_effect = subprocess.TimeoutExpired(cmd="pdb2pqr", timeout=600.0)

    result = poisson_boltzmann_view(obj_name="1LYZ")

    assert "PDB2PQR timed out" in result


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_reports_missing_apbs(mock_sr, mock_run):
    mock_sr.side_effect = _pb_send_request()
    mock_run.side_effect = [MagicMock(returncode=0, stderr=""), FileNotFoundError("apbs")]

    result = poisson_boltzmann_view(obj_name="1LYZ")

    assert "APBS not found on PATH" in result
    assert "brewsci/bio/apbs" in result


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_reports_apbs_timeout(mock_sr, mock_run):
    mock_sr.side_effect = _pb_send_request()
    mock_run.side_effect = [
        MagicMock(returncode=0, stderr=""),
        subprocess.TimeoutExpired(cmd="apbs", timeout=600.0),
    ]

    result = poisson_boltzmann_view(obj_name="1LYZ")

    assert "APBS timed out" in result
    assert "MCPYMOL_PB_TIMEOUT" in result


@patch("subprocess.run")
@patch("mcpymol.views.send_request")
def test_pb_view_passes_a_timeout_to_both_solvers(mock_sr, mock_run):
    """Regression: neither subprocess may run unbounded."""
    mock_sr.side_effect = _pb_send_request()
    mock_run.return_value = MagicMock(returncode=1, stderr="boom")

    poisson_boltzmann_view(obj_name="1LYZ")

    assert mock_run.call_args_list, "solver was never invoked"
    for call in mock_run.call_args_list:
        assert call.kwargs.get("timeout"), f"no timeout on {call.args[0]!r}"


# ── slow operations get a longer socket budget ───────────────────────────────


@pytest.mark.parametrize(
    "func,args,action",
    [
        (ray, ("1920", "1080"), "ray"),
        (draw, ("1920", "1080"), "draw"),
        (mpng, ("frame_",), "mpng"),
        (png, ("out.png",), "png"),
        (save, ("out.pdb",), "save"),
    ],
)
@patch("mcpymol.bridge.send_request")
def test_slow_ops_use_the_long_timeout(mock_sr, func, args, action):
    """Ray-tracing a large scene takes minutes; the 10 s default guaranteed a
    spurious 'Socket connection failed' on exactly the calls worth waiting for."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    func(*args)

    assert mock_sr.call_args.kwargs["timeout"] == _SLOW_OP_TIMEOUT
    assert mock_sr.call_args.args[0] == action


# ── _repair_to_stl: Poisson and the auto fallback chain ──────────────────────


def _poisson_stack(reconstructed_faces=5000, watertight=True):
    """Return (fake_trimesh, fake_pymeshlab, reconstructed_mesh)."""
    fake_tm = MagicMock()
    mesh = MagicMock()
    mesh.faces = list(range(reconstructed_faces))
    mesh.is_watertight = watertight
    fake_tm.load.return_value = mesh

    fake_ml = MagicMock()
    return fake_tm, fake_ml, mesh


def test_repair_poisson_runs_the_reconstruction_pipeline():
    fake_tm, fake_ml, _mesh = _poisson_stack()

    with patch.dict(sys.modules, {"trimesh": fake_tm, "pymeshlab": fake_ml}):
        info = _repair_to_stl("in.obj", "out.stl", "poisson", voxel_pitch=0.7, poisson_depth=9)

    assert info["method"] == "poisson"
    assert info["faces"] == 5000
    assert info["watertight"] is True

    ms = fake_ml.MeshSet.return_value
    ms.load_new_mesh.assert_called_once_with("in.obj")
    ms.save_current_mesh.assert_called_once_with("out.stl", binary=True)

    # The screened-Poisson filter must receive the requested octree depth.
    depths = [c.kwargs.get("depth") for c in ms.apply_filter.call_args_list]
    assert 9 in depths


def test_repair_poisson_tries_legacy_filter_names():
    """pymeshlab renamed its filters between releases; the shim must fall
    through to the older name instead of failing the whole export."""
    fake_tm, fake_ml, _ = _poisson_stack()
    ms = fake_ml.MeshSet.return_value

    # Every modern name raises; only the legacy alias succeeds.
    modern = {
        "meshing_remove_duplicate_vertices",
        "meshing_remove_null_faces",
        "meshing_remove_unreferenced_vertices",
        "compute_normal_per_vertex",
        "generate_surface_reconstruction_screened_poisson",
    }

    def apply_filter(name, **kw):
        if name in modern:
            raise RuntimeError(f"unknown filter {name}")

    ms.apply_filter.side_effect = apply_filter

    with patch.dict(sys.modules, {"trimesh": fake_tm, "pymeshlab": fake_ml}):
        info = _repair_to_stl("in.obj", "out.stl", "poisson", voxel_pitch=0.7, poisson_depth=10)

    assert info["method"] == "poisson"
    tried = [c.args[0] for c in ms.apply_filter.call_args_list]
    assert "surface_reconstruction_screened_poisson" in tried


def test_repair_poisson_raises_when_no_filter_name_works():
    fake_tm, fake_ml, _ = _poisson_stack()
    fake_ml.MeshSet.return_value.apply_filter.side_effect = RuntimeError("nope")

    with (
        patch.dict(sys.modules, {"trimesh": fake_tm, "pymeshlab": fake_ml}),
        pytest.raises(RuntimeError, match="none of"),
    ):
        _repair_to_stl("in.obj", "out.stl", "poisson", voxel_pitch=0.7, poisson_depth=10)


def test_repair_auto_uses_poisson_when_not_watertight():
    """auto on an open surface: Poisson, not the light cleanup."""
    fake_tm, fake_ml, mesh = _poisson_stack()
    raw = MagicMock()
    raw.is_watertight = False
    # First load() is the watertight probe; later loads return the rebuilt mesh.
    fake_tm.load.side_effect = [raw] + [mesh] * 5

    with patch.dict(sys.modules, {"trimesh": fake_tm, "pymeshlab": fake_ml}):
        info = _repair_to_stl("in.obj", "out.stl", "auto", voxel_pitch=0.7, poisson_depth=10)

    assert info["method"] == "poisson"


def test_repair_auto_falls_back_to_voxel_when_poisson_fails():
    """Regression: a Poisson blow-up must degrade to voxel, not kill the export."""
    fake_tm = MagicMock()
    raw = MagicMock()
    raw.is_watertight = False

    voxel_body = MagicMock()
    voxel_body.faces = list(range(777))
    voxel_body.is_watertight = True
    voxel_src = MagicMock()
    voxel_src.voxelized.return_value.fill.return_value.marching_cubes.split.return_value = [
        voxel_body
    ]
    fake_tm.load.side_effect = [raw, voxel_src]

    fake_ml = MagicMock()
    fake_ml.MeshSet.side_effect = RuntimeError("poisson exploded")

    with patch.dict(sys.modules, {"trimesh": fake_tm, "pymeshlab": fake_ml}):
        info = _repair_to_stl("in.obj", "out.stl", "auto", voxel_pitch=0.7, poisson_depth=10)

    assert info["method"] == "voxel (poisson fallback)"
    assert info["faces"] == 777
    voxel_body.export.assert_called_once_with("out.stl")


def test_repair_rejects_unknown_method():
    with (
        patch.dict(sys.modules, {"trimesh": MagicMock()}),
        pytest.raises(ValueError, match="unknown method"),
    ):
        _repair_to_stl("in.obj", "out.stl", "sculpt", voxel_pitch=0.7, poisson_depth=10)


# ── multimer cutoff: one default, documented ─────────────────────────────────


@pytest.mark.parametrize("func", [fetch_structure, load_structure])
def test_multimer_cutoff_default_is_the_documented_constant(func):
    """The helper's own default used to be 5.0 while both callers passed 8.0,
    so the signature contradicted the README."""
    default = inspect.signature(func).parameters["multimer_cutoff"].default
    assert default == DEFAULT_MULTIMER_CUTOFF
    assert inspect.signature(_apply_multimer_heuristic).parameters["cutoff"].default == default


@patch("mcpymol.structures.send_request")
def test_fetch_structure_passes_cutoff_through_to_the_heuristic(mock_sr):
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    fetch_structure(pdb_code="1abc", multimer_cutoff=4.5)

    around = [
        c.kwargs["args"][0]
        for c in mock_sr.call_args_list
        if c.args[0] == "get_chains" and "around" in str(c.kwargs.get("args"))
    ]
    assert around and "around 4.5" in around[0]


# ── _call: the shared body of every primitive wrapper ────────────────────────


@patch("mcpymol.bridge.send_request")
def test_call_forwards_values_in_declaration_order(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    assert _call("zoom", selection="1abc", buffer="5") == "Executed zoom successfully."
    assert mock_sr.call_args.args[0] == "zoom"
    assert mock_sr.call_args.kwargs["args"] == ["1abc", "5"]


@patch("mcpymol.bridge.send_request")
def test_call_drops_unset_trailing_arguments(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    _call("ray", width="1920", height=None)

    assert mock_sr.call_args.kwargs["args"] == ["1920"]


@patch("mcpymol.bridge.send_request")
def test_call_with_no_arguments_sends_none(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    _call("deselect")

    assert mock_sr.call_args.kwargs["args"] == []


@patch("mcpymol.bridge.send_request")
def test_call_rejects_a_gap_instead_of_mislabelling_it(mock_sr):
    """The old wrappers filtered out None anywhere in the list, so a height
    with no width was sent as the *width* — a silently wrong render."""
    result = _call("ray", width=None, height="1080")

    assert result.startswith("Error: ray was given 'height' without 'width'")
    mock_sr.assert_not_called()


@patch("mcpymol.bridge.send_request")
def test_call_names_every_missing_argument(mock_sr):
    result = _call("symexp", prefix=None, selection=None, cutoff="20")

    assert "'prefix'" in result and "'selection'" in result
    mock_sr.assert_not_called()


@patch("mcpymol.bridge.send_request")
def test_call_propagates_plugin_errors(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no such object: nope"}

    assert _call("zoom", selection="nope") == "no such object: nope"


@patch("mcpymol.bridge.send_request")
def test_call_falls_back_when_error_field_is_missing(mock_sr):
    mock_sr.return_value = {"status": "error"}

    assert _call("zoom", selection="x") == "Unknown error"


@patch("mcpymol.bridge.send_request")
def test_call_uses_the_default_timeout(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    _call("zoom", selection="all")

    assert mock_sr.call_args.kwargs["timeout"] == _DEFAULT_TIMEOUT


@patch("mcpymol.primitives.send_request")
def test_wrapper_surfaces_the_gap_error(mock_sr):
    """End-to-end through a real tool, not just the helper."""
    assert ray(height="1080").startswith("Error: ray was given 'height' without 'width'")
    mock_sr.assert_not_called()


# ── tool registry ────────────────────────────────────────────────────────────


def _registered_tools():
    import asyncio

    from mcpymol.server import mcp

    return {t.name: t for t in asyncio.run(mcp.list_tools())}


def test_every_tool_has_a_description_and_schema():
    """Guards the wrapper refactor: the bodies are shared, but each tool must
    still present its own real signature and docstring to the MCP client."""
    tools = _registered_tools()

    assert len(tools) > 100
    for name, tool in tools.items():
        assert (tool.description or "").strip(), f"{name} has no description"
        assert tool.inputSchema.get("type") == "object", f"{name} has no object schema"


@pytest.mark.parametrize(
    "name,required,optional",
    [
        ("zoom", [], ["selection", "buffer"]),
        ("util_cbc", [], ["selection"]),
        ("deselect", [], []),
        ("set", ["setting", "value"], ["selection"]),  # tool-name override
        ("as", ["representation"], ["selection"]),
        ("super", ["mobile"], ["target", "options"]),
        ("symexp", ["prefix", "obj_name", "selection"], ["cutoff", "segi"]),
    ],
)
def test_wrapper_schemas_match_their_signatures(name, required, optional):
    """Including the three tools whose MCP name differs from the Python name."""
    schema = _registered_tools()[name].inputSchema

    assert sorted(schema.get("required", [])) == sorted(required)
    assert sorted(schema["properties"]) == sorted(required + optional)
