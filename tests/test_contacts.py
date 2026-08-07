"""Tests for contact_report and its interaction geometry.

The classification criteria are the point of this tool, so these tests build
atoms at known separations and assert the call, rather than mocking the answer.
"""

import math
from unittest.mock import patch

import pytest

from mcpymol.analysis import (
    HBOND_MAX,
    HYDROPHOBIC_MAX,
    SALT_BRIDGE_MAX,
    _classify,
    _neighbour_pairs,
    _plane_normal,
    _ring_systems,
    _stacking_interactions,
    contact_report,
)
from mcpymol.pdbtext import Atom


def _atom(name="CA", resn="ALA", chain="A", resi="1", x=0.0, y=0.0, z=0.0, element=None, het=False):
    return Atom(
        name=name,
        resn=resn,
        chain=chain,
        resi=resi,
        x=x,
        y=y,
        z=z,
        occupancy=1.0,
        bfactor=0.0,
        element=(element or name.lstrip("0123456789")[:1]).upper(),
        hetatm=het,
    )


def _pdb_line(atom):
    record = "HETATM" if atom.hetatm else "ATOM  "
    return (
        f"{record}    1  {atom.name:<3}{atom.resn:>4} {atom.chain}{atom.resi:>4}    "
        f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}  1.00  0.00          "
        f"{atom.element:>2}"
    )


def _sr(atoms1, atoms2):
    """Serve two selections as PDB text, distinguished by call order."""
    calls = []

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_pdbstr":
            calls.append(args[0])
            group = atoms1 if len(calls) == 1 else atoms2
            return {"status": "success", "result": "\n".join(_pdb_line(a) for a in group)}
        return {"status": "success", "result": "OK"}

    return fake


# ── neighbour search ─────────────────────────────────────────────────────────


def test_neighbour_pairs_finds_close_atoms():
    a = [_atom(x=0.0)]
    b = [_atom(x=2.0, resi="2"), _atom(x=50.0, resi="3")]

    pairs = _neighbour_pairs(a, b, cutoff=4.0)

    assert len(pairs) == 1
    assert pairs[0][2] == pytest.approx(2.0)


def test_neighbour_pairs_respects_the_cutoff_exactly():
    a = [_atom(x=0.0)]
    b = [_atom(x=4.0, resi="2")]

    assert len(_neighbour_pairs(a, b, cutoff=4.0)) == 1
    assert _neighbour_pairs(a, b, cutoff=3.99) == []


def test_neighbour_pairs_spans_grid_cell_boundaries():
    """A pair straddling a cell boundary must still be found — the classic
    bug in grid-based neighbour searches."""
    a = [_atom(x=3.99)]
    b = [_atom(x=4.01, resi="2")]

    assert len(_neighbour_pairs(a, b, cutoff=4.0)) == 1


def test_neighbour_pairs_handles_negative_coordinates():
    """Floor division on negatives is where naive grid indexing breaks."""
    a = [_atom(x=-0.1, y=-0.1, z=-0.1)]
    b = [_atom(x=0.1, y=0.1, z=0.1, resi="2")]

    assert len(_neighbour_pairs(a, b, cutoff=4.0)) == 1


def test_neighbour_pairs_is_empty_for_empty_input():
    assert _neighbour_pairs([], [_atom()], 4.0) == []
    assert _neighbour_pairs([_atom()], [], 4.0) == []


def test_neighbour_pairs_finds_everything_a_brute_force_scan_would():
    """The grid is an optimisation; it must not change the answer."""
    import random

    rng = random.Random(11)
    a = [
        _atom(resi=str(i), x=rng.uniform(-20, 20), y=rng.uniform(-20, 20), z=rng.uniform(-20, 20))
        for i in range(120)
    ]
    b = [
        _atom(resi=str(i), x=rng.uniform(-20, 20), y=rng.uniform(-20, 20), z=rng.uniform(-20, 20))
        for i in range(120)
    ]

    grid = {(p[0].resi, p[1].resi) for p in _neighbour_pairs(a, b, 4.0)}
    brute = {
        (x.resi, y.resi) for x in a for y in b if math.dist((x.x, x.y, x.z), (y.x, y.y, y.z)) <= 4.0
    }

    assert grid == brute


# ── classification ───────────────────────────────────────────────────────────


def test_classify_salt_bridge():
    asp = _atom(name="OD2", resn="ASP")
    lys = _atom(name="NZ", resn="LYS", resi="2")

    assert _classify(asp, lys, 3.2) == "salt bridge"
    assert _classify(lys, asp, 3.2) == "salt bridge"


def test_classify_salt_bridge_respects_its_cutoff():
    asp = _atom(name="OD2", resn="ASP")
    lys = _atom(name="NZ", resn="LYS", resi="2")

    # Beyond 4.0 the charge pair is no longer a salt bridge, but the two are
    # still N and O, so it degrades to a polar contact rather than vanishing.
    assert _classify(asp, lys, SALT_BRIDGE_MAX + 0.1) == "polar contact"


def test_classify_hydrogen_bond():
    donor = _atom(name="ND1", resn="HIS")
    acceptor = _atom(name="O", resn="ALA", resi="2")

    assert _classify(donor, acceptor, HBOND_MAX - 0.1) == "hydrogen bond"


def test_classify_polar_beyond_hbond_range():
    n = _atom(name="N", resn="ALA")
    o = _atom(name="O", resn="GLY", resi="2")

    assert _classify(n, o, HBOND_MAX + 0.2) == "polar contact"


def test_classify_hydrophobic():
    a = _atom(name="CB", resn="LEU")
    b = _atom(name="CG", resn="VAL", resi="2")

    assert _classify(a, b, HYDROPHOBIC_MAX - 0.5) == "hydrophobic"


def test_classify_sulfur_counts_as_hydrophobic():
    a = _atom(name="SD", resn="MET", element="S")
    b = _atom(name="CB", resn="LEU", resi="2")

    assert _classify(a, b, 4.0) == "hydrophobic"


def test_classify_falls_back_to_plain_contact():
    a = _atom(name="CB", resn="LEU")
    b = _atom(name="CG", resn="VAL", resi="2")

    assert _classify(a, b, HYDROPHOBIC_MAX + 1.0) == "contact"


def test_histidine_counts_as_cationic():
    """HIS is protonated often enough at physiological pH that its salt
    bridges are conventionally reported."""
    his = _atom(name="NE2", resn="HIS")
    glu = _atom(name="OE1", resn="GLU", resi="2")

    assert _classify(his, glu, 3.5) == "salt bridge"


# ── aromatic rings ───────────────────────────────────────────────────────────


def _phe_ring(chain="A", resi="1", centre=(0.0, 0.0, 0.0), plane="xy"):
    """A hexagonal PHE ring of radius 1.4 A in the requested plane."""
    names = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
    atoms = []
    for i, name in enumerate(names):
        angle = 2 * math.pi * i / 6
        dx, dy = 1.4 * math.cos(angle), 1.4 * math.sin(angle)
        # "xz" is the same ring rotated 90 degrees, for T-shaped geometry.
        offset = (dx, dy, 0.0) if plane == "xy" else (dx, 0.0, dy)
        atoms.append(
            _atom(
                name=name,
                resn="PHE",
                chain=chain,
                resi=resi,
                x=centre[0] + offset[0],
                y=centre[1] + offset[1],
                z=centre[2] + offset[2],
            )
        )
    return atoms


def test_plane_normal_is_a_unit_vector():
    n = _plane_normal([(0, 0, 0), (1, 0, 0), (0, 1, 0)])

    assert n is not None
    assert math.isclose(math.sqrt(sum(c * c for c in n)), 1.0)
    assert abs(n[2]) == pytest.approx(1.0)  # xy-plane normal points along z


def test_plane_normal_rejects_collinear_points():
    assert _plane_normal([(0, 0, 0), (1, 0, 0), (2, 0, 0)]) is None


def test_ring_systems_finds_the_phenylalanine_ring():
    rings = _ring_systems(_phe_ring())

    assert len(rings) == 1
    _key, label, centroid, _normal = rings[0]
    assert label == "PHE1"
    assert centroid == pytest.approx((0.0, 0.0, 0.0), abs=1e-6)


def test_ring_systems_skips_incomplete_sidechains():
    """Disordered aromatics are common; a partial ring must not produce a
    bogus centroid."""
    partial = _phe_ring()[:2]

    assert _ring_systems(partial) == []


def test_ring_systems_finds_both_tryptophan_rings():
    """The indole is two fused rings, and both stack — so both are reported."""
    # A planar (xy) approximation of indole: a 6-ring and a fused 5-ring.
    six = ["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"]
    five = ["CG", "CD1", "NE1", "CE2", "CD2"]
    coords = {}
    for i, name in enumerate(six):
        angle = 2 * math.pi * i / 6
        coords[name] = (1.4 * math.cos(angle), 1.4 * math.sin(angle))
    for i, name in enumerate(five):
        if name in coords:
            continue
        angle = 2 * math.pi * i / 5
        coords[name] = (-2.2 + 1.2 * math.cos(angle), 1.2 * math.sin(angle))

    atoms = [_atom(name=n, resn="TRP", x=xy[0], y=xy[1]) for n, xy in coords.items()]

    assert len(_ring_systems(atoms)) == 2


def test_ring_systems_ignores_non_aromatic_residues():
    assert _ring_systems([_atom(name="CB", resn="LEU")]) == []


def test_stacking_detects_parallel_rings():
    a = _phe_ring(chain="A", resi="1", centre=(0, 0, 0), plane="xy")
    b = _phe_ring(chain="B", resi="2", centre=(0, 0, 4.0), plane="xy")

    found = _stacking_interactions(a, b)

    assert len(found) == 1
    kind, gap = next(iter(found.values()))
    assert kind == "pi-stacking (parallel)"
    assert gap == pytest.approx(4.0)


def test_stacking_detects_t_shaped_rings():
    a = _phe_ring(chain="A", resi="1", centre=(0, 0, 0), plane="xy")
    b = _phe_ring(chain="B", resi="2", centre=(0, 0, 5.0), plane="xz")

    kind, _gap = next(iter(_stacking_interactions(a, b).values()))

    assert kind == "pi-stacking (T-shaped)"


def test_stacking_ignores_distant_rings():
    a = _phe_ring(chain="A", resi="1", centre=(0, 0, 0))
    b = _phe_ring(chain="B", resi="2", centre=(0, 0, 20.0))

    assert _stacking_interactions(a, b) == {}


def test_stacking_ignores_a_ring_against_itself():
    ring = _phe_ring(chain="A", resi="1")

    assert _stacking_interactions(ring, ring) == {}


# ── contact_report ───────────────────────────────────────────────────────────


@patch("mcpymol.analysis.send_request")
def test_contact_report_lists_pairs_closest_first(mock_sr):
    ligand = [_atom(name="O1", resn="LIG", chain="L", resi="1", x=0.0, het=True)]
    protein = [
        _atom(name="NZ", resn="LYS", chain="A", resi="27", x=2.6),
        _atom(name="CB", resn="LEU", chain="A", resi="30", x=3.9),
    ]
    mock_sr.side_effect = _sr(ligand, protein)

    result = contact_report(selection1="lig", selection2="prot")

    assert "2 residue pairs in contact" in result
    assert result.index("LYS27") < result.index("LEU30")
    assert "2.60 A" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_names_the_interaction_types(mock_sr):
    asp = [_atom(name="OD2", resn="ASP", chain="A", resi="30", x=0.0)]
    lys = [_atom(name="NZ", resn="LYS", chain="B", resi="27", x=3.0)]
    mock_sr.side_effect = _sr(asp, lys)

    result = contact_report(selection1="a", selection2="b")

    assert "salt bridge" in result
    assert "Interaction types across all pairs" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_aggregates_atoms_into_residue_pairs(mock_sr):
    """A residue pair touching at four atoms is one interaction, not four."""
    a = [_atom(name=f"C{i}", resn="LIG", chain="L", resi="1", x=0.0, y=i * 0.3) for i in range(4)]
    b = [_atom(name="CB", resn="LEU", chain="A", resi="30", x=3.0)]
    mock_sr.side_effect = _sr(a, b)

    result = contact_report(selection1="a", selection2="b")

    assert "1 residue pairs" in result
    assert "4 atom contacts" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_reports_stacking(mock_sr):
    a = _phe_ring(chain="A", resi="10", centre=(0, 0, 0))
    b = _phe_ring(chain="B", resi="20", centre=(0, 0, 4.0))
    mock_sr.side_effect = _sr(a, b)

    result = contact_report(selection1="a", selection2="b", cutoff=5.0)

    assert "pi-stacking (parallel)" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_caps_output_but_says_so(mock_sr):
    """Silent truncation would read as 'these are all the contacts'."""
    a = [_atom(name="CA", chain="L", resi=str(i), x=i * 0.1) for i in range(30)]
    b = [_atom(name="CB", chain="A", resi=str(100 + i), x=i * 0.1 + 1.0) for i in range(30)]
    mock_sr.side_effect = _sr(a, b)

    result = contact_report(selection1="a", selection2="b", max_pairs=5)

    assert "more pairs (raise max_pairs to see)" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_reports_no_contacts_helpfully(mock_sr):
    mock_sr.side_effect = _sr([_atom(x=0.0)], [_atom(x=100.0, resi="2")])

    result = contact_report(selection1="a", selection2="b")

    assert "No contacts within" in result
    assert "count_atoms" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_flags_an_empty_selection(mock_sr):
    mock_sr.side_effect = _sr([], [_atom()])

    assert "matched no atoms" in contact_report(selection1="nothing", selection2="b")


@patch("mcpymol.analysis.send_request")
def test_contact_report_excludes_water_by_default(mock_sr):
    mock_sr.side_effect = _sr([_atom()], [_atom(resi="2")])

    contact_report(selection1="a", selection2="b")

    queries = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "get_pdbstr"]
    assert all("not solvent" in q for q in queries)


@patch("mcpymol.analysis.send_request")
def test_contact_report_can_include_water(mock_sr):
    mock_sr.side_effect = _sr([_atom()], [_atom(resi="2")])

    contact_report(selection1="a", selection2="b", include_water=True)

    queries = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "get_pdbstr"]
    assert all("not solvent" not in q for q in queries)


@patch("mcpymol.analysis.send_request")
def test_contact_report_ignores_hydrogens(mock_sr):
    """Heavy-atom criteria would be double-counted by explicit hydrogens."""
    a = [_atom(name="H1", element="H", x=0.0), _atom(name="CA", resi="1", x=0.0)]
    b = [_atom(name="CB", chain="B", resi="9", x=3.0)]
    mock_sr.side_effect = _sr(a, b)

    result = contact_report(selection1="a", selection2="b")

    assert "1 residue pairs" in result


@patch("mcpymol.analysis.send_request")
def test_contact_report_skips_a_residue_contacting_itself(mock_sr):
    """Overlapping selections would otherwise report every residue against
    itself at distance zero."""
    same = [_atom(name="CA", chain="A", resi="5", x=0.0)]
    mock_sr.side_effect = _sr(same, same)

    assert "No contacts" in contact_report(selection1="x", selection2="x")


@patch("mcpymol.analysis.send_request")
def test_contact_report_rejects_a_bad_cutoff(mock_sr):
    assert "cutoff must be positive" in contact_report(selection1="a", selection2="b", cutoff=0)
    mock_sr.assert_not_called()


@patch("mcpymol.analysis.send_request")
def test_contact_report_propagates_read_errors(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no such object"}

    assert "Error reading" in contact_report(selection1="nope", selection2="b")


def test_contact_report_is_registered():
    import asyncio

    from mcpymol.server import mcp

    assert "contact_report" in {t.name for t in asyncio.run(mcp.list_tools())}
