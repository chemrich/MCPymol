"""Tests for superposition_view and its residue-pairing helpers."""

import json
from unittest.mock import patch

import pytest

from mcpymol.comparison import (
    _distance,
    _match_residues,
    _parse_ca_coords,
    _resi_sort_key,
    superposition_view,
)


def _atom(resi, chain, x, y, z, name="CA"):
    return (
        f"ATOM  {resi:>5}  {name:<3} ALA {chain}{resi:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C"
    )


def _pdb(atoms):
    return "\n".join(atoms) + "\nEND\n"


def _sr(mobile_pdb, target_pdb, fit_result=(1.25, 300, 5, 2.0, 320, 0.0, 150)):
    """send_request stand-in serving two structures and a superposition result."""
    calls = {"n": 0}

    def fake(action, args=None, kwargs=None, **_ignored):
        if action in ("super", "align"):
            return {"status": "success", "result": list(fit_result)}
        if action == "get_pdbstr":
            calls["n"] += 1
            return {
                "status": "success",
                "result": mobile_pdb if args[0] == "_sup_mobile" else target_pdb,
            }
        return {"status": "success", "result": "OK"}

    return fake


# ── coordinate parsing ───────────────────────────────────────────────────────


def test_parse_ca_coords_reads_fixed_width_columns():
    text = _pdb([_atom(10, "A", 1.0, 2.0, 3.0), _atom(11, "A", 4.5, 5.5, 6.5)])

    assert _parse_ca_coords(text) == {
        ("A", "10"): (1.0, 2.0, 3.0),
        ("A", "11"): (4.5, 5.5, 6.5),
    }


def test_parse_ca_coords_ignores_non_ca_atoms():
    """Only alpha carbons — otherwise a residue contributes many points."""
    text = _pdb([_atom(1, "A", 0, 0, 0, name="CB"), _atom(1, "A", 9, 9, 9, name="CA")])

    assert _parse_ca_coords(text) == {("A", "1"): (9.0, 9.0, 9.0)}


def test_parse_ca_coords_skips_malformed_lines():
    text = "ATOM      1  CA  ALA A   1      bad  coords  here\n" + _atom(2, "A", 1, 1, 1)

    assert list(_parse_ca_coords(text)) == [("A", "2")]


# ── pairing ──────────────────────────────────────────────────────────────────


def test_match_residues_pairs_on_chain_and_number():
    mobile = {("A", "1"): (0.0, 0.0, 0.0), ("A", "2"): (0.0, 0.0, 0.0)}
    target = {("A", "1"): (3.0, 4.0, 0.0), ("A", "2"): (0.0, 0.0, 0.0)}

    assert _match_residues(mobile, target) == [("A", "1", 5.0), ("A", "2", 0.0)]


def test_match_residues_ignores_residues_present_in_only_one():
    mobile = {("A", "1"): (0.0, 0.0, 0.0), ("A", "99"): (0.0, 0.0, 0.0)}
    target = {("A", "1"): (0.0, 0.0, 1.0)}

    assert _match_residues(mobile, target) == [("A", "1", 1.0)]


def test_match_residues_falls_back_to_residue_number():
    """Same protein, different chain IDs — strict matching finds nothing, so
    the fallback keeps the comparison useful instead of returning empty."""
    mobile = {("A", "5"): (0.0, 0.0, 0.0)}
    target = {("B", "5"): (0.0, 0.0, 2.0)}

    assert _match_residues(mobile, target) == [("", "5", 2.0)]


def test_match_residues_returns_empty_when_nothing_lines_up():
    assert _match_residues({("A", "1"): (0.0, 0.0, 0.0)}, {("B", "900"): (0.0, 0.0, 0.0)}) == []


def test_match_residues_is_sorted_by_position():
    mobile = {("A", str(i)): (0.0, 0.0, 0.0) for i in (2, 10, 1)}
    target = dict(mobile)

    assert [r[1] for r in _match_residues(mobile, target)] == ["1", "2", "10"]


@pytest.mark.parametrize("a,b,expected", [((0, 0, 0), (0, 0, 0), 0.0), ((0, 0, 0), (3, 4, 0), 5.0)])
def test_distance(a, b, expected):
    assert _distance(a, b) == expected


def test_resi_sort_key_handles_insertion_codes():
    assert _resi_sort_key("10") < _resi_sort_key("10A") < _resi_sort_key("11")


def test_resi_sort_key_survives_garbage():
    assert _resi_sort_key("") == (0, "")


# ── superposition_view ───────────────────────────────────────────────────────


@patch("mcpymol.comparison.send_request")
def test_superposition_view_reports_rmsd_and_shifts(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0), _atom(2, "A", 0, 0, 0), _atom(3, "A", 0, 0, 0)])
    target = _pdb([_atom(1, "A", 0, 0, 0), _atom(2, "A", 0, 0, 6), _atom(3, "A", 0, 0, 3)])
    mock_sr.side_effect = _sr(mobile, target)

    result = superposition_view(mobile="1ake", target="4ake")

    assert "RMSD 1.25 A over 300 atoms" in result
    assert "Compared 3 residues" in result
    assert "mean shift 3.00 A" in result
    assert "max 6.00 A" in result
    # Worst first, so the interesting residue is named up front.
    assert result.index("2 (6.0 A)") < result.index("3 (3.0 A)")


@patch("mcpymol.comparison.send_request")
def test_superposition_view_writes_deviations_into_b_factors(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0), _atom(2, "A", 0, 0, 0)])
    target = _pdb([_atom(1, "A", 0, 0, 0), _atom(2, "A", 0, 0, 4)])
    mock_sr.side_effect = _sr(mobile, target)

    superposition_view(mobile="m", target="t")

    scripts = [
        c.kwargs["args"][0]
        for c in mock_sr.call_args_list
        if c.args[0] == "do" and "stored.dev" in str(c.kwargs.get("args"))
    ]
    assert len(scripts) == 1, "deviations should go over in one batched script"
    assert json.dumps({"1": 0.0, "2": 4.0}) in scripts[0]
    assert "alter (m), b=0" in scripts[0]


@patch("mcpymol.comparison.send_request")
def test_superposition_view_scales_colour_to_the_largest_shift(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    target = _pdb([_atom(1, "A", 0, 0, 7)])
    mock_sr.side_effect = _sr(mobile, target)

    superposition_view(mobile="m", target="t")

    spectra = [
        c.kwargs["args"][0]
        for c in mock_sr.call_args_list
        if c.args[0] == "do" and "spectrum" in str(c.kwargs.get("args"))
    ]
    assert "maximum=7.0" in spectra[0]


@patch("mcpymol.comparison.send_request")
def test_superposition_view_honours_an_explicit_ceiling(mock_sr):
    """Needed to put two different comparisons on one colour scale."""
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    target = _pdb([_atom(1, "A", 0, 0, 7)])
    mock_sr.side_effect = _sr(mobile, target)

    result = superposition_view(mobile="m", target="t", max_deviation=3.0)

    spectra = [
        c.kwargs["args"][0]
        for c in mock_sr.call_args_list
        if c.args[0] == "do" and "spectrum" in str(c.kwargs.get("args"))
    ]
    assert "maximum=3.0" in spectra[0]
    assert "shifted >= 3.0 A" in result


@patch("mcpymol.comparison.send_request")
def test_superposition_view_copies_objects_before_reading_coordinates(mock_sr):
    """super/align leave the fit in the object matrix, so raw coordinates can
    be stale — copying bakes the transform in. Regression guard."""
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    mock_sr.side_effect = _sr(mobile, mobile)

    superposition_view(mobile="m", target="t")

    order = [c.args[0] for c in mock_sr.call_args_list]
    assert order.index("create") < order.index("get_pdbstr")
    created = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "create"]
    assert created == ["_sup_mobile", "_sup_target"]


@patch("mcpymol.comparison.send_request")
def test_superposition_view_cleans_up_its_temp_objects(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    mock_sr.side_effect = _sr(mobile, mobile)

    superposition_view(mobile="m", target="t")

    deleted = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "delete"]
    assert deleted.count("_sup_mobile") >= 1
    assert deleted.count("_sup_target") >= 1


@patch("mcpymol.comparison.send_request")
def test_superposition_view_cleans_up_after_a_read_failure(mock_sr):
    def fake(action, args=None, kwargs=None, **_ignored):
        if action in ("super", "align"):
            return {"status": "success", "result": [1.0, 10]}
        if action == "get_pdbstr":
            return {"status": "error", "error": "object not found"}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = superposition_view(mobile="m", target="t")

    assert "Error reading coordinates" in result
    deleted = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "delete"]
    assert "_sup_mobile" in deleted and "_sup_target" in deleted


@patch("mcpymol.comparison.send_request")
def test_superposition_view_reports_a_failed_fit(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no atoms in selection"}

    result = superposition_view(mobile="m", target="t")

    assert "Error superposing m onto t" in result
    assert "no atoms in selection" in result


@patch("mcpymol.comparison.send_request")
def test_superposition_view_explains_when_nothing_pairs(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    target = _pdb([_atom(900, "B", 0, 0, 0)])
    mock_sr.side_effect = _sr(mobile, target)

    result = superposition_view(mobile="m", target="t")

    assert "no residues could be paired" in result
    assert "RMSD 1.25 A" in result  # the fit itself still succeeded


@patch("mcpymol.comparison.send_request")
def test_superposition_view_survives_an_unparseable_fit_result(mock_sr):
    """Some builds return something other than the 7-tuple; the per-residue
    analysis is still the valuable part."""
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    mock_sr.side_effect = _sr(mobile, mobile, fit_result=())

    result = superposition_view(mobile="m", target="t")

    assert "RMSD" not in result
    assert "Compared 1 residues" in result


@patch("mcpymol.comparison.send_request")
def test_superposition_view_can_use_align(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    mock_sr.side_effect = _sr(mobile, mobile)

    superposition_view(mobile="m", target="t", method="align")

    assert "align" in [c.args[0] for c in mock_sr.call_args_list]


@patch("mcpymol.comparison.send_request")
def test_superposition_view_rejects_a_bad_method(mock_sr):
    assert "must be 'super' or 'align'" in superposition_view(mobile="m", target="t", method="fit")
    mock_sr.assert_not_called()


@patch("mcpymol.comparison.send_request")
def test_superposition_view_rejects_comparing_an_object_to_itself(mock_sr):
    assert "same object" in superposition_view(mobile="m", target="m")
    mock_sr.assert_not_called()


@patch("mcpymol.comparison.send_request")
def test_superposition_view_leaves_the_target_as_a_grey_reference(mock_sr):
    mobile = _pdb([_atom(1, "A", 0, 0, 0)])
    mock_sr.side_effect = _sr(mobile, mobile)

    superposition_view(mobile="m", target="t")

    colors = [c.kwargs["args"] for c in mock_sr.call_args_list if c.args[0] == "color"]
    assert ["grey70", "t"] in colors


def test_superposition_view_is_registered_as_a_tool():
    import asyncio

    from mcpymol.server import mcp

    assert "superposition_view" in {t.name for t in asyncio.run(mcp.list_tools())}
