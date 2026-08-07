"""Tests for interface_report and buried-surface-area bookkeeping.

The arithmetic is the whole tool: free SASA minus bound SASA, halved to get
the per-side interface area people quote. These tests drive it with known
areas so the numbers are checked rather than mocked.
"""

from unittest.mock import patch

import pytest

from mcpymol.analysis import (
    _CRYSTAL_CONTACT_MAX_AREA,
    _LARGE_INTERFACE_MIN_AREA,
    _residue_class,
    _residue_sasa,
    interface_report,
)


def _atom_line(resi, chain, resn="ALA", name="CA", sasa=0.0):
    """A CA record carrying per-atom SASA in the B-factor column, which is
    where get_area(load_b=1) puts it."""
    return (
        f"ATOM  {resi:>5}  {name:<3}{resn:>4} {chain}{resi:>4}    "
        f"   1.000   2.000   3.000  1.00{sasa:>6.2f}           C"
    )


def _sr_interface(free, bound, fail_area=False, fail_create=False):
    """Serve per-chain SASA dumps.

    ``free`` / ``bound`` map chain -> {resi: (resn, area)}.
    """

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "create" and fail_create:
            return {"status": "error", "error": "no such chain"}
        if action == "get_area":
            if fail_area:
                return {"status": "error", "error": "surface calculation failed"}
            return {"status": "success", "result": 100.0}
        if action == "get_pdbstr":
            sel = args[0]
            source = bound if sel.startswith("_iface_bound") else free
            chain = sel[-1] if "chain " in sel else ("A" if sel.endswith("free_a") else "B")
            rows = source.get(chain, {})
            return {
                "status": "success",
                "result": "\n".join(
                    _atom_line(resi, chain, resn=resn, sasa=area)
                    for resi, (resn, area) in rows.items()
                ),
            }
        return {"status": "success", "result": "OK"}

    return fake


# ── helpers ──────────────────────────────────────────────────────────────────


@patch("mcpymol.analysis.send_request")
def test_residue_sasa_sums_atoms_per_residue(mock_sr):
    """get_area writes per-atom SASA; a residue's area is the sum over its
    atoms."""

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_pdbstr":
            return {
                "status": "success",
                "result": "\n".join(
                    [
                        _atom_line(10, "A", name="CA", sasa=12.0),
                        _atom_line(10, "A", name="CB", sasa=8.0),
                        _atom_line(11, "A", name="CA", sasa=5.0),
                    ]
                ),
            }
        return {"status": "success", "result": 1.0}

    mock_sr.side_effect = fake

    result = _residue_sasa("x")

    assert result == {("A", "10"): ("ALA", 20.0), ("A", "11"): ("ALA", 5.0)}


@patch("mcpymol.analysis.send_request")
def test_residue_sasa_requests_load_b(mock_sr):
    """Without load_b the B-factor column keeps temperature factors and the
    numbers would be silently meaningless."""
    mock_sr.return_value = {"status": "success", "result": ""}

    _residue_sasa("x")

    area_call = next(c for c in mock_sr.call_args_list if c.args[0] == "get_area")
    assert area_call.kwargs["kwargs"]["load_b"] == 1


@patch("mcpymol.analysis.send_request")
def test_residue_sasa_returns_none_on_error(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "nope"}

    assert _residue_sasa("x") is None


@pytest.mark.parametrize(
    "resn,expected",
    [("ASP", "charged"), ("ARG", "charged"), ("LEU", "hydrophobic"), ("SER", "polar")],
)
def test_residue_class(resn, expected):
    assert _residue_class(resn) == expected


# ── interface_report ─────────────────────────────────────────────────────────


@patch("mcpymol.analysis.send_request")
def test_interface_report_computes_buried_area(mock_sr):
    """Chain A buries 300 A^2 and chain B 300 A^2, so the total is 600 and the
    per-side figure people quote is 300."""
    free = {
        "A": {1: ("LEU", 200.0), 2: ("ASP", 200.0)},
        "B": {1: ("ARG", 200.0), 2: ("SER", 200.0)},
    }
    bound = {
        "A": {1: ("LEU", 50.0), 2: ("ASP", 50.0)},
        "B": {1: ("ARG", 50.0), 2: ("SER", 50.0)},
    }
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="1brs", chain_a="A", chain_b="B")

    assert "300 A^2 per side" in result
    assert "600 A^2 total" in result
    assert "2 in chain A, 2 in chain B" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_ranks_hot_spots_by_burial(mock_sr):
    free = {"A": {1: ("LEU", 200.0), 2: ("ASP", 100.0)}, "B": {1: ("ARG", 100.0)}}
    bound = {"A": {1: ("LEU", 20.0), 2: ("ASP", 90.0)}, "B": {1: ("ARG", 50.0)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    # LEU1 buries 180, ASP2 only 10 — LEU must come first.
    assert result.index("LEU1") < result.index("ASP2")
    assert "LEU1 (180 A^2)" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_calls_a_small_interface_a_crystal_contact(mock_sr):
    free = {"A": {1: ("LEU", 100.0)}, "B": {1: ("LEU", 100.0)}}
    bound = {"A": {1: ("LEU", 90.0)}, "B": {1: ("LEU", 90.0)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    assert "crystal packing contact" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_calls_a_large_interface_specific(mock_sr):
    big = _LARGE_INTERFACE_MIN_AREA + 500
    free = {"A": {1: ("LEU", big * 2)}, "B": {1: ("LEU", big * 2)}}
    bound = {"A": {1: ("LEU", big)}, "B": {1: ("LEU", big)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    assert "large interface" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_uses_the_middle_verdict_between_thresholds(mock_sr):
    mid = (_CRYSTAL_CONTACT_MAX_AREA + _LARGE_INTERFACE_MIN_AREA) / 2
    free = {"A": {1: ("LEU", mid * 2)}, "B": {1: ("LEU", mid * 2)}}
    bound = {"A": {1: ("LEU", mid)}, "B": {1: ("LEU", mid)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    assert "transient" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_breaks_down_composition(mock_sr):
    free = {"A": {1: ("LEU", 200.0)}, "B": {1: ("ASP", 200.0)}}
    bound = {"A": {1: ("LEU", 100.0)}, "B": {1: ("ASP", 100.0)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    assert "50% hydrophobic" in result
    assert "50% charged" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_detects_no_contact(mock_sr):
    """Two chains that bury nothing are simply not in contact."""
    free = {"A": {1: ("LEU", 100.0)}, "B": {1: ("LEU", 100.0)}}
    mock_sr.side_effect = _sr_interface(free, free)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    assert "not in contact" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_sets_the_solvent_dot_mode(mock_sr):
    """dot_solvent=0 measures the molecular surface, not the accessible one —
    the numbers would be wrong but plausible-looking."""
    free = {"A": {1: ("LEU", 100.0)}, "B": {1: ("LEU", 100.0)}}
    bound = {"A": {1: ("LEU", 10.0)}, "B": {1: ("LEU", 10.0)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    interface_report(obj_name="x", chain_a="A", chain_b="B")

    settings = {tuple(c.kwargs["args"]) for c in mock_sr.call_args_list if c.args[0] == "set"}
    assert ("dot_solvent", "1") in settings


@patch("mcpymol.analysis.send_request")
def test_interface_report_cleans_up_temp_objects(mock_sr):
    free = {"A": {1: ("LEU", 100.0)}, "B": {1: ("LEU", 100.0)}}
    bound = {"A": {1: ("LEU", 10.0)}, "B": {1: ("LEU", 10.0)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    interface_report(obj_name="x", chain_a="A", chain_b="B")

    deleted = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "delete"]
    for tmp in ("_iface_free_a", "_iface_free_b", "_iface_bound"):
        assert deleted.count(tmp) >= 2  # cleared before use and after


@patch("mcpymol.analysis.send_request")
def test_interface_report_cleans_up_after_a_failure(mock_sr):
    mock_sr.side_effect = _sr_interface({}, {}, fail_create=True)

    result = interface_report(obj_name="x", chain_a="A", chain_b="Z")

    assert "Error isolating" in result
    deleted = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "delete"]
    assert "_iface_bound" in deleted


@patch("mcpymol.analysis.send_request")
def test_interface_report_reports_a_failed_area_calculation(mock_sr):
    mock_sr.side_effect = _sr_interface({}, {}, fail_area=True)

    assert "Error computing surface area" in interface_report(
        obj_name="x", chain_a="A", chain_b="B"
    )


@patch("mcpymol.analysis.send_request")
def test_interface_report_flags_an_empty_chain(mock_sr):
    mock_sr.side_effect = _sr_interface({"A": {}, "B": {}}, {"A": {}, "B": {}})

    result = interface_report(obj_name="x", chain_a="A", chain_b="B")

    assert "no polymer atoms" in result
    assert "list_chains" in result


@patch("mcpymol.analysis.send_request")
def test_interface_report_rejects_identical_chains(mock_sr):
    assert "pick two different chains" in interface_report(obj_name="x", chain_a="A", chain_b="A")
    mock_sr.assert_not_called()


@patch("mcpymol.analysis.send_request")
def test_interface_report_caps_the_hot_spot_list_but_says_so(mock_sr):
    free = {"A": {i: ("LEU", 200.0) for i in range(1, 21)}, "B": {1: ("LEU", 200.0)}}
    bound = {"A": {i: ("LEU", 100.0) for i in range(1, 21)}, "B": {1: ("LEU", 100.0)}}
    mock_sr.side_effect = _sr_interface(free, bound)

    result = interface_report(obj_name="x", chain_a="A", chain_b="B", max_residues=5)

    assert "and 15 more" in result


def test_interface_report_is_registered():
    import asyncio

    from mcpymol.server import mcp

    assert "interface_report" in {t.name for t in asyncio.run(mcp.list_tools())}
