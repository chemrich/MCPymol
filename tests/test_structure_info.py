"""Tests for the shared PDB parser, structure_info and get_sequence."""

import urllib.error
from unittest.mock import patch

import pytest

from mcpymol.pdbtext import Atom, parse_atoms, residue_order
from mcpymol.structures import (
    _rcsb_metadata,
    atom_properties,
    fetch_structure,
    get_sequence,
    load_structure,
    structure_info,
)


def _atom_line(resi=1, chain="A", name="CA", resn="ALA", x=1.0, y=2.0, z=3.0, b=50.0, het=False):
    record = "HETATM" if het else "ATOM  "
    element = name.lstrip("0123456789")[:1]
    return (
        f"{record}{resi:>5}  {name:<3}{resn:>4} {chain}{resi:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00{b:>6.2f}          {element:>2}"
    )


# ── parse_atoms ──────────────────────────────────────────────────────────────


def test_parse_atoms_reads_every_field():
    (atom,) = parse_atoms(_atom_line(resi=27, chain="B", name="NZ", resn="LYS", b=31.5))

    assert atom == Atom(
        name="NZ",
        resn="LYS",
        chain="B",
        resi="27",
        x=1.0,
        y=2.0,
        z=3.0,
        occupancy=1.0,
        bfactor=31.5,
        element="N",
        hetatm=False,
    )


def test_parse_atoms_marks_hetatm():
    (atom,) = parse_atoms(_atom_line(resn="HOH", name="O", het=True))
    assert atom.hetatm is True


def test_parse_atoms_ca_only_filter():
    text = "\n".join([_atom_line(name="CA"), _atom_line(name="CB"), _atom_line(name="N")])

    assert [a.name for a in parse_atoms(text, ca_only=True)] == ["CA"]


def test_parse_atoms_skips_headers_and_terminators():
    text = "HEADER    HYDROLASE\n" + _atom_line() + "\nTER\nEND\n"

    assert len(parse_atoms(text)) == 1


def test_parse_atoms_skips_truncated_lines():
    """PyMOL dumps are not always well-formed; one bad line must not lose the
    other few thousand good ones."""
    text = "ATOM      1  CA  ALA A   1\n" + _atom_line(resi=2)

    assert [a.resi for a in parse_atoms(text)] == ["2"]


def test_parse_atoms_skips_unparseable_coordinates():
    text = "ATOM      1  CA  ALA A   1       xx.x    2.000   3.000  1.00 50.00           C"

    assert parse_atoms(text) == []


def test_parse_atoms_infers_element_when_the_column_is_missing():
    """Older files omit columns 77-78 entirely."""
    line = _atom_line(name="CB")[:76]

    (atom,) = parse_atoms(line)
    assert atom.element == "C"


def test_parse_atoms_infers_element_past_a_leading_digit():
    line = _atom_line(name="1HB")[:76]

    (atom,) = parse_atoms(line)
    assert atom.element == "H"


def test_atom_labels_read_naturally():
    (atom,) = parse_atoms(_atom_line(resi=27, chain="B", resn="LYS"))

    assert atom.label == "B/LYS27"
    assert atom.residue_key == ("B", "27")


@pytest.mark.parametrize(
    "a,b", [("9", "10"), ("10", "10A"), ("10A", "11"), ("-5", "1"), ("99", "100")]
)
def test_residue_order_sorts_numerically_with_insertion_codes(a, b):
    assert residue_order(a) < residue_order(b)


def test_residue_order_survives_nonsense():
    assert residue_order("") == (0, "")
    assert residue_order("XYZ") == (0, "XYZ")


# ── RCSB metadata ────────────────────────────────────────────────────────────


ENTRY_JSON = {
    "struct": {"title": "HIV-1 PROTEASE WITH INHIBITOR"},
    "exptl": [{"method": "X-RAY DIFFRACTION"}],
    "rcsb_entry_info": {"resolution_combined": [2.0]},
    "rcsb_accession_info": {"initial_release_date": "1995-01-26T00:00:00Z"},
    "rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1"]},
}
ENTITY_JSON = {
    "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Human immunodeficiency virus 1"}]
}


@patch("mcpymol.structures._rcsb_get")
def test_rcsb_metadata_extracts_the_useful_fields(mock_get):
    mock_get.side_effect = [ENTRY_JSON, ENTITY_JSON]

    meta = _rcsb_metadata("1hsg")

    assert meta["title"] == "HIV-1 PROTEASE WITH INHIBITOR"
    assert meta["method"] == "X-RAY DIFFRACTION"
    assert meta["resolution"] == 2.0
    assert meta["released"] == "1995-01-26"
    assert meta["organisms"] == ["Human immunodeficiency virus 1"]


@patch("mcpymol.structures._rcsb_get")
def test_rcsb_metadata_is_empty_when_the_api_fails(mock_get):
    mock_get.return_value = None

    assert _rcsb_metadata("1hsg") == {}


@patch("mcpymol.structures.urllib.request.urlopen")
def test_rcsb_get_swallows_network_failures(mock_open):
    """Metadata is a bonus on top of what PyMOL knows; an unreachable API must
    degrade to 'no metadata', never fail the tool."""
    from mcpymol.structures import _rcsb_get

    mock_open.side_effect = urllib.error.URLError("offline")

    assert _rcsb_get("entry/1hsg") is None


# ── structure_info ───────────────────────────────────────────────────────────


def _sr_structure(chains=("A", "B"), counts=None, pdb="", states=1, symmetry=None, ligands=""):
    counts = counts or {}

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_chains":
            return {"status": "success", "result": list(chains)}
        if action == "count_atoms":
            sel = args[0]
            for key, value in counts.items():
                if key in sel:
                    return {"status": "success", "result": value}
            # Default to "something loaded" so the loaders' emptiness guard
            # does not fire in tests that are about something else.
            return {"status": "success", "result": 1200}
        if action == "get_object_list":
            return {"status": "success", "result": []}
        if action == "count_states":
            return {"status": "success", "result": states}
        if action == "get_symmetry":
            return {"status": "success", "result": symmetry}
        if action == "get_pdbstr":
            return {"status": "success", "result": ligands if "organic" in args[0] else pdb}
        return {"status": "success", "result": "OK"}

    return fake


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_combines_pymol_and_rcsb(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_structure(
        chains=("A", "B"),
        counts={"and solvent": 120, "name CA": 198, "(1hsg)": 1846},
        ligands=_atom_line(resn="MK1", het=True),
    )
    mock_meta.return_value = {
        "title": "HIV-1 PROTEASE",
        "method": "X-RAY DIFFRACTION",
        "resolution": 2.0,
        "released": "1995-01-26",
        "organisms": ["Human immunodeficiency virus 1"],
    }

    result = structure_info(obj_name="1hsg")

    assert "HIV-1 PROTEASE" in result
    assert "X-RAY DIFFRACTION" in result
    assert "2.00 A resolution" in result
    assert "Human immunodeficiency virus 1" in result
    assert "Chains (2): A, B" in result
    assert "MK1" in result


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_works_without_metadata(mock_sr, mock_meta):
    """A locally loaded or renamed structure still gets the PyMOL-side facts."""
    mock_sr.side_effect = _sr_structure(chains=("A",))
    mock_meta.return_value = {}

    result = structure_info(obj_name="my_model")

    assert "my_model" in result
    assert "Chains (1): A" in result


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_skips_the_api_for_non_pdb_names(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_structure(chains=("A",))

    structure_info(obj_name="my_local_model")

    mock_meta.assert_not_called()


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_accepts_an_explicit_pdb_id(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_structure(chains=("A",))
    mock_meta.return_value = {"title": "T"}

    structure_info(obj_name="renamed_thing", pdb_id="1hsg")

    mock_meta.assert_called_once_with("1hsg")


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_flags_a_probable_prediction(mock_sr, mock_meta):
    """pLDDT in the B-factor column changes which tools are meaningful."""
    mock_sr.side_effect = _sr_structure(
        chains=("A",),
        pdb="\n".join(_atom_line(resi=i, b=b) for i, b in enumerate([95.0, 88.0, 42.0], start=1)),
    )
    mock_meta.return_value = {}

    result = structure_info(obj_name="AF_P69905")

    assert "pLDDT" in result
    assert "plddt_view" in result


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_reports_nmr_ensembles(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_structure(chains=("A",), states=20)
    mock_meta.return_value = {}

    assert "20 states" in structure_info(obj_name="2kkk")


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_structure_info_reports_the_space_group(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_structure(
        chains=("A",), symmetry=[50.8, 40.7, 60.1, 90.0, 90.0, 90.0, "P 21 21 21"]
    )
    mock_meta.return_value = {}

    assert "P 21 21 21" in structure_info(obj_name="1ubq")


@patch("mcpymol.structures.send_request")
def test_structure_info_reports_a_missing_object(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "Invalid selection name"}

    assert "Error inspecting 'nope'" in structure_info(obj_name="nope")


# ── get_sequence ─────────────────────────────────────────────────────────────


def _sr_sequence(fasta, ca_lines):
    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_fastastr":
            return {"status": "success", "result": fasta}
        if action == "get_pdbstr":
            return {"status": "success", "result": "\n".join(ca_lines)}
        return {"status": "success", "result": "OK"}

    return fake


@patch("mcpymol.structures.send_request")
def test_get_sequence_returns_fasta(mock_sr):
    mock_sr.side_effect = _sr_sequence(
        ">1hsg_A\nPQITLW\n", [_atom_line(resi=i) for i in range(1, 7)]
    )

    result = get_sequence(obj_name="1hsg", chain="A")

    assert ">1hsg_A" in result
    assert "PQITLW" in result
    assert "6 modelled residues" in result


@patch("mcpymol.structures.send_request")
def test_get_sequence_reports_the_numbering_offset(mock_sr):
    """PDB numbering rarely starts at 1, so 'residue 50' in a paper and
    position 50 in the sequence are usually different residues."""
    mock_sr.side_effect = _sr_sequence(">x\nAAAA\n", [_atom_line(resi=i) for i in range(21, 25)])

    result = get_sequence(obj_name="x")

    assert "numbered 21-24" in result
    assert "offset by 20" in result


@patch("mcpymol.structures.send_request")
def test_get_sequence_reports_chain_breaks(mock_sr):
    """Unmodelled loops are invisible in the sequence but very much present
    in the structure."""
    mock_sr.side_effect = _sr_sequence(
        ">x\nAAAA\n",
        [_atom_line(resi=i) for i in (1, 2, 9, 10)],
    )

    result = get_sequence(obj_name="x")

    assert "Chain breaks" in result
    assert "2->9 (6 missing)" in result


@patch("mcpymol.structures.send_request")
def test_get_sequence_handles_multiple_chains(mock_sr):
    mock_sr.side_effect = _sr_sequence(
        ">x_A\nAA\n>x_B\nAA\n",
        [
            _atom_line(resi=1, chain="A"),
            _atom_line(resi=2, chain="A"),
            _atom_line(resi=1, chain="B"),
            _atom_line(resi=2, chain="B"),
        ],
    )

    result = get_sequence(obj_name="x")

    assert "Chain A:" in result
    assert "Chain B:" in result


@patch("mcpymol.structures.send_request")
def test_get_sequence_reports_an_empty_object(mock_sr):
    mock_sr.side_effect = _sr_sequence("", [])

    assert "No protein sequence found" in get_sequence(obj_name="x")


@patch("mcpymol.structures.send_request")
def test_get_sequence_propagates_errors(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no such object"}

    assert "Error reading sequence" in get_sequence(obj_name="nope")


@pytest.mark.parametrize("name", ["structure_info", "get_sequence"])
def test_new_tools_are_registered(name):
    import asyncio

    from mcpymol.server import mcp

    assert name in {t.name for t in asyncio.run(mcp.list_tools())}


# ── loading more than one structure ──────────────────────────────────────────


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_fetch_structure_clears_the_session_by_default(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_structure(chains=("A",))
    mock_meta.return_value = {}

    fetch_structure(pdb_code="1ubq")

    dos = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    # Settings are reset; objects are cleared by name, after the fetch is known
    # to have worked. A blanket "reinitialize" would run before it.
    assert "reinitialize settings" in dos


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_fetch_structure_can_add_to_the_session(mock_sr, mock_meta):
    """Without this, comparing two structures is impossible: the second fetch
    wipes the first, so superposition_view has nothing to compare."""
    mock_sr.side_effect = _sr_structure(chains=("A",))
    mock_meta.return_value = {}

    fetch_structure(pdb_code="4ake", replace=False)

    dos = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert "reinitialize" not in dos


@patch("mcpymol.structures.send_request")
def test_load_structure_can_add_to_the_session(mock_sr):
    mock_sr.side_effect = _sr_structure(chains=("A",))

    load_structure(file_path="/tmp/model.pdb", obj_name="model", replace=False)

    dos = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert "reinitialize" not in dos


# ── a failed load must not cost you the session ──────────────────────────────
#
# Reported in #15: behind a VPN blocking the RCSB, a fetch returned a success
# message naming the object and the session was gone. cmd.fetch does not raise
# when the download fails, and the plugin reports success regardless, so the
# only defence is checking that atoms actually arrived.


def _sr_empty_fetch(existing=("1ubq", "4hhb")):
    """A fetch that 'succeeds' while producing nothing, as a blocked download
    does."""
    deleted = []

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "count_atoms":
            return {"status": "success", "result": 0}
        if action == "get_object_list":
            return {"status": "success", "result": list(existing)}
        if action == "delete":
            deleted.append(args[0])
        return {"status": "success", "result": "OK"}

    fake.deleted = deleted
    return fake


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_empty_fetch_reports_an_error_not_success(mock_sr, mock_meta):
    mock_sr.side_effect = _sr_empty_fetch()
    mock_meta.return_value = {}

    result = fetch_structure(pdb_code="1abc")

    assert "Loaded nothing for '1abc'" in result
    assert "Successfully fetched" not in result
    assert "Nothing else in the session was touched" in result


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_empty_fetch_leaves_other_objects_alone(mock_sr, mock_meta):
    """The destructive part of #15: the session was cleared before the fetch,
    so a failed download took unrelated structures with it."""
    fake = _sr_empty_fetch(existing=("1ubq", "4hhb"))
    mock_sr.side_effect = fake
    mock_meta.return_value = {}

    fetch_structure(pdb_code="1abc")

    assert "1ubq" not in fake.deleted
    assert "4hhb" not in fake.deleted
    dos = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert not any("reinitialize" in d for d in dos), "session was cleared on a failed fetch"


@patch("mcpymol.structures.send_request")
def test_empty_load_reports_an_error(mock_sr):
    mock_sr.side_effect = _sr_empty_fetch()

    result = load_structure(file_path="/tmp/empty.pdb", obj_name="model")

    assert "Loaded nothing for '/tmp/empty.pdb'" in result


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_successful_fetch_clears_only_other_objects(mock_sr, mock_meta):
    """replace=True still gives a clean session — but by name, after the fetch
    is known to have worked, rather than by wiping everything beforehand."""

    deleted = []

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "count_atoms":
            return {"status": "success", "result": 1200}
        if action == "get_object_list":
            return {"status": "success", "result": ["1ubq", "4hhb", "1abc"]}
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        if action == "delete":
            deleted.append(args[0])
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake
    mock_meta.return_value = {}

    result = fetch_structure(pdb_code="1abc")

    assert "Successfully fetched 1abc" in result
    assert "1ubq" in deleted and "4hhb" in deleted
    assert deleted.count("1abc") == 1, "the fetched object must survive the cleanup"


@patch("mcpymol.structures._rcsb_metadata")
@patch("mcpymol.structures.send_request")
def test_settings_are_reset_without_deleting_objects(mock_sr, mock_meta):
    """A previous preset's fog or ray_trace_mode would otherwise leak into the
    new scene, which is what the blanket reinitialize was buying."""

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "count_atoms":
            return {"status": "success", "result": 1200}
        if action == "get_object_list":
            return {"status": "success", "result": ["1abc"]}
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake
    mock_meta.return_value = {}

    fetch_structure(pdb_code="1abc")

    dos = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert "reinitialize settings" in dos


# ── atom_properties ──────────────────────────────────────────────────────────
#
# The only route to properties that live on atoms rather than residues or
# objects: occupancy, altloc, per-atom B-factor. cmd.get_model returns a chempy
# object that does not survive JSON, so the plugin's iterate_to_list action is
# what makes them reachable at all.


def _sr_atoms(rows, status="success"):
    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "iterate_to_list":
            if status == "error":
                return {"status": "error", "error": "invalid selection"}
            return {"status": "success", "result": rows}
        return {"status": "success", "result": "OK"}

    return fake


@patch("mcpymol.structures.send_request")
def test_atom_properties_renders_a_table(mock_sr):
    mock_sr.side_effect = _sr_atoms(
        [["A", "25", "ASP", "OD1", 18.42, 1.0], ["A", "25", "ASP", "OD2", 21.07, 0.5]]
    )

    result = atom_properties(selection="1hsg and resi 25")

    assert "2 atoms" in result
    assert "OD1" in result and "OD2" in result
    assert "18.42" in result  # floats rendered to 2 dp
    assert "0.50" in result  # partial occupancy is the point of the tool


@patch("mcpymol.structures.send_request")
def test_atom_properties_passes_the_requested_fields(mock_sr):
    mock_sr.side_effect = _sr_atoms([["A", 0.5]])

    atom_properties(selection="x", properties="chain, q")

    call = next(c for c in mock_sr.call_args_list if c.args[0] == "iterate_to_list")
    assert call.kwargs["args"][1] == "chain, q"


@patch("mcpymol.structures.send_request")
def test_atom_properties_tolerates_untidy_field_lists(mock_sr):
    mock_sr.side_effect = _sr_atoms([["A", 0.5]])

    atom_properties(selection="x", properties=" chain , q ,")

    call = next(c for c in mock_sr.call_args_list if c.args[0] == "iterate_to_list")
    assert call.kwargs["args"][1] == "chain, q"


@patch("mcpymol.structures.send_request")
def test_atom_properties_shows_blanks_rather_than_swallowing_them(mock_sr):
    """An empty chain or altloc is information, not absence."""
    mock_sr.side_effect = _sr_atoms([["", "1", "ALA", "CA", 0.0, 1.0]])

    result = atom_properties(selection="x")

    assert "-" in result


@patch("mcpymol.structures.send_request")
def test_atom_properties_caps_output_but_says_so(mock_sr):
    mock_sr.side_effect = _sr_atoms([["A", str(i)] for i in range(100)])

    result = atom_properties(selection="x", properties="chain, resi", max_atoms=5)

    assert "100 atoms" in result
    assert "and 95 more atoms" in result


@patch("mcpymol.structures.send_request")
def test_atom_properties_reports_an_empty_selection(mock_sr):
    mock_sr.side_effect = _sr_atoms([])

    result = atom_properties(selection="chain Z")

    assert "No atoms matched" in result
    assert "count_atoms" in result


@patch("mcpymol.structures.send_request")
def test_atom_properties_propagates_errors(mock_sr):
    mock_sr.side_effect = _sr_atoms([], status="error")

    assert "Error reading" in atom_properties(selection="nope")


@pytest.mark.parametrize("bad", ["", "  ", " , , "])
@patch("mcpymol.structures.send_request")
def test_atom_properties_rejects_an_empty_property_list(mock_sr, bad):
    assert "no properties requested" in atom_properties(selection="x", properties=bad)
    mock_sr.assert_not_called()


@patch("mcpymol.structures.send_request")
def test_atom_properties_rejects_a_nonsense_cap(mock_sr):
    assert "at least 1" in atom_properties(selection="x", max_atoms=0)
    mock_sr.assert_not_called()


def test_atom_properties_is_registered():
    import asyncio

    from mcpymol.server import mcp

    assert "atom_properties" in {t.name for t in asyncio.run(mcp.list_tools())}
