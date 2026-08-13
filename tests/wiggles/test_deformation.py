"""Tests for deformation_view — the well-supported tool, and its one boundary.

Motion survives blind testing; populations do not. So the assertions here are
about getting the motion right, plus one that matters more than it looks: the
report must never let an arrow be read as a population. An arrow is a statement
about where density moved, not about how many particles moved that way, and
nothing in the picture distinguishes those.
"""

from __future__ import annotations

import pytest

from mcpymol.wiggles.deformation import (
    deformation_view,
    per_atom_displacement,
    per_residue_displacement,
    read_uncertainty_table,
)
from mcpymol.wiggles.port import FakePort, PortError

# chain, resi, resn, name, alt, q, b
ROWS = [
    ("A", "1", "ALA", "N", "", 1.0, 20.0),
    ("A", "1", "ALA", "CA", "", 1.0, 20.0),
    ("A", "2", "GLY", "CA", "", 1.0, 30.0),
    ("A", "3", "SER", "CA", "", 1.0, 40.0),
]


def _coords(shifts):
    """State-1 coords at the origin, state-2 shifted along x by ``shifts``."""
    start = [(0.0, float(i), 0.0) for i in range(len(ROWS))]
    end = [(dx, float(i), 0.0) for i, dx in enumerate(shifts)]
    return start, end


def _port(shifts, *, n_states=2):
    start, end = _coords(shifts)
    return FakePort(
        {
            "iterate_to_list": list(ROWS),
            "count_states": n_states,
            "get_coords": lambda obj, state=1: start if state == 1 else end,
        }
    )


# ── the arithmetic ───────────────────────────────────────────────────────────


def test_displacement_is_the_distance_moved():
    got = per_atom_displacement([(0.0, 0.0, 0.0)], [(3.0, 4.0, 0.0)])
    assert got == pytest.approx([5.0])


def test_mismatched_atom_counts_are_refused_not_approximated():
    """The same refusal morph_states makes, for the same reason: unpaired
    atoms make the field meaningless rather than merely imprecise."""
    with pytest.raises(PortError, match="cannot be paired"):
        per_atom_displacement([(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])


def test_the_refusal_points_at_the_tool_that_needs_no_pairing():
    with pytest.raises(PortError, match="ensemble_spread_view"):
        per_atom_displacement([(0.0, 0.0, 0.0)], [])


def test_per_residue_is_the_mean_over_its_atoms():
    atoms = _port([0.0] * 4).query("iterate_to_list", "x", "y")
    from mcpymol.wiggles.atoms import fetch_atoms

    got = per_residue_displacement(fetch_atoms(_port([0.0] * 4), "obj"), [2.0, 4.0, 1.0, 9.0])
    assert got[("A", "1")] == pytest.approx(3.0)
    assert got[("A", "2")] == pytest.approx(1.0)
    assert atoms  # the fixture really did hand back rows


def test_diverged_lengths_are_refused():
    from mcpymol.wiggles.atoms import fetch_atoms

    atoms = fetch_atoms(_port([0.0] * 4), "obj")
    with pytest.raises(PortError, match="orders have diverged"):
        per_residue_displacement(atoms, [1.0])


# ── refusals ─────────────────────────────────────────────────────────────────


def test_a_single_state_object_is_refused():
    out = deformation_view(_port([0.0] * 4, n_states=1), "obj")
    assert "REFUSED" in out
    assert "nothing here to difference" in out


def test_an_out_of_range_state_is_refused():
    with pytest.raises(PortError, match=r"outside the object's 1\.\.2 states"):
        deformation_view(_port([1.0] * 4), "obj", end_state=7)


def test_a_state_against_itself_is_refused():
    with pytest.raises(PortError, match="zero everywhere"):
        deformation_view(_port([1.0] * 4), "obj", start_state=1, end_state=1)


# ── the view ─────────────────────────────────────────────────────────────────


def test_displacement_reaches_the_bfactor_column():
    port = _port([0.0, 0.0, 3.0, 6.0])
    deformation_view(port, "obj", arrows=False)

    altered = {args[0]: args[1] for args, _ in port.calls("alter")}
    assert "b=3.0000" in altered['(obj) and chain "A" and resi "2"']
    assert "b=6.0000" in altered['(obj) and chain "A" and resi "3"']


def test_the_report_names_the_biggest_movers():
    out = deformation_view(_port([0.0, 0.0, 3.0, 6.0]), "obj", arrows=False)
    assert "A/3 6.00Å" in out
    assert "Moved most" in out


def test_bfactors_are_stashed_before_being_overwritten():
    out = deformation_view(_port([1.0] * 4), "obj", arrows=False)
    assert "restore_bfactors" in out


def test_arrows_are_drawn_as_cgo():
    port = _port([0.0, 0.0, 3.0, 6.0])
    deformation_view(port, "obj")

    calls = port.calls("load_cgo")
    assert len(calls) == 1, port.call_log
    cgo, name = calls[0][0]
    assert name == "obj_arrows"
    assert cgo[0] == 9.0, "a cylinder opcode starts the shaft"
    assert 27.0 in cgo, "a cone opcode makes the head"


def test_only_ca_atoms_get_arrows():
    """One arrow per residue, at its CA — an arrow per atom is a thicket."""
    port = _port([5.0, 5.0, 5.0, 5.0])
    out = deformation_view(port, "obj")
    assert "3 drawn" in out, out


def test_motions_below_the_noise_floor_get_no_arrow():
    port = _port([0.0, 0.01, 0.01, 5.0])
    out = deformation_view(port, "obj")

    assert "1 drawn" in out
    assert "moved less than" in out
    assert "assert a motion nobody measured" in out


def test_the_arrow_ceiling_is_reported_never_silent():
    port = _port([0.0, 5.0, 6.0, 7.0])
    out = deformation_view(port, "obj", max_arrows=1)

    assert "2 further residue(s)" in out
    assert "truncation of the largest movers" in out


def test_exaggeration_is_labelled():
    out = deformation_view(_port([0.0, 0.0, 3.0, 6.0]), "obj", arrow_scale=3.0)
    assert "EXAGGERATED 3x" in out
    assert "no longer" in out


def test_unexaggerated_output_says_nothing_about_exaggeration():
    out = deformation_view(_port([0.0, 0.0, 3.0, 6.0]), "obj")
    assert "EXAGGERATED" not in out


# ── half-set uncertainty ─────────────────────────────────────────────────────


def test_missing_uncertainty_is_stated_loudly():
    """A uniformly orange field is not evidence of uniform reliability."""
    out = deformation_view(_port([0.0, 0.0, 3.0, 6.0]), "obj")

    assert "NO UNCERTAINTY SUPPLIED" in out
    assert "do not have identical confidence" in out


def test_an_uncertainty_table_is_read_and_reported(tmp_path):
    path = tmp_path / "unc.txt"
    path.write_text("# chain resi value\nA 2 0.1\nA 3 0.9\n")

    out = deformation_view(_port([0.0, 0.0, 3.0, 6.0]), "obj", uncertainty=str(path))
    assert "coloured by half-set uncertainty" in out
    assert "2 residues" in out
    assert "NO UNCERTAINTY SUPPLIED" not in out


def test_uncertainty_changes_the_arrow_colours(tmp_path):
    """Not just "a table was accepted" — the geometry must be identical and the
    colours must not be, or the table is being read and thrown away."""
    path = tmp_path / "unc.txt"
    path.write_text("A 2 0.0\nA 3 1.0\n")

    plain = _port([0.0, 0.0, 3.0, 6.0])
    deformation_view(plain, "obj")
    coloured = _port([0.0, 0.0, 3.0, 6.0])
    deformation_view(coloured, "obj", uncertainty=str(path))

    a = plain.calls("load_cgo")[0][0][0]
    b = coloured.calls("load_cgo")[0][0][0]

    assert len(a) == len(b), "same arrows, so the same CGO length"
    assert a != b, "the uncertainty table should have changed the colours"
    # The shaft endpoints are geometry and must be untouched; CYLINDER lays out
    # opcode + 6 coordinates before the radius and colours begin.
    assert a[1:7] == b[1:7], "geometry must not depend on the uncertainty table"


def test_separators_are_flexible(tmp_path):
    path = tmp_path / "unc.csv"
    path.write_text("A,2,0.4\nA\t3\t0.6\n")
    assert read_uncertainty_table(path) == {("A", "2"): 0.4, ("A", "3"): 0.6}


def test_an_unparseable_table_says_what_it_wanted(tmp_path):
    path = tmp_path / "unc.txt"
    path.write_text("nonsense\n")
    with pytest.raises(ValueError, match="chain resi value"):
        read_uncertainty_table(path)


# ── the claim boundary ───────────────────────────────────────────────────────


def test_the_report_forbids_reading_a_population_off_an_arrow():
    out = deformation_view(_port([0.0, 0.0, 3.0, 6.0]), "obj")

    assert "MOTION IS THE WELL-SUPPORTED CLAIM" in out
    assert "Do NOT read them as: this fraction of particles moved" in out


def test_blank_chain_and_negative_resi_stay_scoped():
    """The same selection defect qscore and ensembles had, at this call site.

    A blank chain leaves `chain  and resi 2`, where PyMOL takes `and` as the
    chain name and the selection stops being scoped to the object — verified
    against PyMOL 3.1.0, where it matched more atoms than the named object
    held. A negative residue number is read as a range. Both are routine in
    real depositions, and this module writes B-factors, so an unscoped alter
    destroys every other loaded structure's values.
    """
    rows = [
        ("", "-3", "ALA", "CA", "", 1.0, 20.0),
        ("", "-2", "GLY", "CA", "", 1.0, 30.0),
    ]
    start = [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    end = [(1.0, 0.0, 0.0), (2.0, 1.0, 0.0)]
    port = FakePort(
        {
            "iterate_to_list": rows,
            "count_states": 2,
            "get_coords": lambda obj, state=1: start if state == 1 else end,
        }
    )

    deformation_view(port, "obj")

    selections = [args[0] for args, _ in port.calls("alter") if args]
    assert selections, "no alter issued"
    for sel in selections:
        assert "chain  and" not in sel, f"unscoped selection: {sel}"
    assert any('chain "" and resi "\\-3"' in s for s in selections), selections
    for sel in selections:
        assert 'resi "-3"' not in sel, f"quoted but not escaped, still a range: {sel}"
