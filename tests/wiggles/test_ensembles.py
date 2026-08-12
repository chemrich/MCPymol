"""Tests for multi-state ensemble views."""

from __future__ import annotations

import pytest

from mcpymol.wiggles.atoms import Atom
from mcpymol.wiggles.bfactors import has_stash
from mcpymol.wiggles.ensembles import (
    SPREAD_LEGEND,
    ensemble_spread_view,
    morph_states,
    per_atom_spread,
    per_residue_spread,
)
from mcpymol.wiggles.port import FakePort, PortError


def make_port(n_states, coords_by_state, atom_rows):
    """FakePort answering count_states, iterate_to_list, and get_coords.

    get_coords is a callable response so each state returns its own coords.
    """
    return FakePort(
        {
            "count_states": n_states,
            "iterate_to_list": atom_rows,
            # callable response: varies the answer by the state argument
            "get_coords": lambda _sel, state: coords_by_state[state - 1],
        }
    )


ROWS_2 = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "ALA", "CA", "", 1.0, 20.0),
]


# -- the maths -------------------------------------------------------------


def test_no_motion_gives_zero_spread():
    spread = per_atom_spread([[(0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0)]])
    assert spread == [0.0]


def test_spread_is_rms_about_the_centroid():
    """Two states at ±1 along x: centroid 0, RMS deviation 1."""
    spread = per_atom_spread([[(-1.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]])
    assert spread[0] == pytest.approx(1.0)


def test_spread_needs_two_states():
    with pytest.raises(PortError, match="at least two states"):
        per_atom_spread([[(0.0, 0.0, 0.0)]])


def test_differing_atom_counts_make_spread_undefined():
    with pytest.raises(PortError, match="differing atom counts"):
        per_atom_spread([[(0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]])


def test_residue_spread_averages_its_atoms():
    atoms = [
        Atom("A", "1", "MET", "N", "", 1.0, 0.0),
        Atom("A", "1", "MET", "CA", "", 1.0, 0.0),
        Atom("A", "2", "ALA", "CA", "", 1.0, 0.0),
    ]
    out = per_residue_spread(atoms, [1.0, 3.0, 5.0])
    assert out[("A", "1")] == pytest.approx(2.0)
    assert out[("A", "2")] == pytest.approx(5.0)


def test_misaligned_atoms_and_spread_is_an_error():
    """iterate order and coordinate order diverging is silent corruption."""
    with pytest.raises(PortError, match="orders have diverged"):
        per_residue_spread([Atom("A", "1", "MET", "CA", "", 1.0, 0.0)], [1.0, 2.0])


# -- ensemble_spread_view --------------------------------------------------


def test_single_state_object_declines_and_redirects():
    port = FakePort({"count_states": 1})
    out = ensemble_spread_view(port, "obj")
    assert "needs at least 2" in out
    assert "altloc_view" in out  # points at the right tool for the other case
    assert not port.queried("spectrum")


def test_spread_view_colours_and_putties():
    port = make_port(
        2,
        [[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]],
        ROWS_2,
    )
    out = ensemble_spread_view(port, "obj")
    assert "2 states, 2 atoms, 2 residues" in out
    assert port.queried("spectrum")
    assert port.called("cartoon", "putty", "obj"), port.call_log
    assert SPREAD_LEGEND in out


def test_spread_legend_denies_being_an_error_bar():
    port = make_port(
        2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], [("A", "1", "MET", "CA", "", 1.0, 20.0)]
    )
    out = ensemble_spread_view(port, "obj")
    assert "NOT a calibrated uncertainty" in out
    assert "not mean the true position is less well determined" in out


def test_spread_view_preserves_bfactors():
    port = make_port(
        2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], [("A", "1", "MET", "CA", "", 1.0, 20.0)]
    )
    out = ensemble_spread_view(port, "obj")
    assert has_stash("obj")
    assert "restore_bfactors" in out


def test_putty_can_be_disabled():
    port = make_port(
        2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], [("A", "1", "MET", "CA", "", 1.0, 20.0)]
    )
    ensemble_spread_view(port, "obj", as_putty=False)
    assert not port.queried("cartoon")


# -- morph_states — the refusal is the point -------------------------------


def test_morph_refuses_differing_topology_and_explains():
    port = make_port(2, [[(0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]], ROWS_2)
    out = morph_states(port, "obj")
    assert "REFUSED" in out
    assert "state 1: 1" in out and "state 2: 2" in out
    assert "ensemble_spread_view" in out  # offers the tool that does work
    assert not port.queried("morph"), port.call_log


def test_morph_refuses_single_state():
    port = FakePort({"count_states": 1})
    out = morph_states(port, "obj")
    assert "REFUSED" in out
    assert not port.queried("morph")


def test_morph_proceeds_on_shared_topology():
    port = make_port(2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ROWS_2)
    out = morph_states(port, "obj", steps=10)
    assert "Topology check passed" in out
    assert port.called("morph", "obj_morph", "obj", refinement=0, steps=10), port.call_log


def test_morph_states_says_interpolated_frames_are_invented():
    """SPEC invariant I1 in spirit: generated content must be labelled."""
    port = make_port(2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ROWS_2)
    out = morph_states(port, "obj")
    assert "generated, not observed" in out


def test_validate_only_creates_nothing():
    port = make_port(2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ROWS_2)
    out = morph_states(port, "obj", validate_only=True)
    assert "no morph created" in out
    assert not port.queried("morph")


def test_morph_honours_custom_name():
    port = make_port(2, [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)]], ROWS_2)
    morph_states(port, "obj", name="swing")
    assert port.calls("morph")[0][0][:2] == ("swing", "obj"), port.call_log


def test_morph_degrades_when_pymol_is_open_source():
    """cmd.morph is Incentive-only. Most users are on open-source PyMOL, so
    for most users the interpolation half of this tool does not exist. The
    validation half does, and it is the part carrying a judgement."""

    class RefusingPort(FakePort):
        def query(self, action, *args, **kwargs):
            if action == "morph":
                raise PortError(
                    "morph failed: PyMOL execution error: Incentive-Only-Error: "
                    '"morph" is not available in Open-Source PyMOL'
                )
            return super().query(action, *args, **kwargs)

    port = RefusingPort(
        {
            "count_states": 2,
            "iterate_to_list": ROWS_2,
            "get_coords": lambda _s, state: [(0.0, 0.0, 0.0)],
        }
    )
    out = morph_states(port, "obj")

    assert "Topology check passed" in out
    assert "Incentive-only" in out
    assert "set all_states, on" in out  # a way to see the motion regardless
    assert "ensemble_spread_view" in out  # and the tool that needs no licence


def test_other_morph_errors_still_raise():
    """Only the licence case degrades — a real failure must not be swallowed."""

    class BrokenPort(FakePort):
        def query(self, action, *args, **kwargs):
            if action == "morph":
                raise PortError("morph failed: Error: Invalid selection name")
            return super().query(action, *args, **kwargs)

    port = BrokenPort(
        {
            "count_states": 2,
            "iterate_to_list": ROWS_2,
            "get_coords": lambda _s, state: [(0.0, 0.0, 0.0)],
        }
    )
    with pytest.raises(PortError, match="Invalid selection"):
        morph_states(port, "obj")


def test_blank_chain_and_negative_resi_stay_scoped():
    """The same selection defect qscore had, at the ensembles call site.

    The reported fix named qscore.py only; ensembles.py builds the identical
    string and would have kept writing one residue's spread across a range —
    or, with a blank chain, across every atom in the session.
    """
    rows = [
        ("", "-3", "MET", "CA", "", 1.0, 20.0),
        ("", "-2", "ALA", "CA", "", 1.0, 20.0),
    ]
    port = make_port(
        n_states=2,
        atom_rows=rows,
        coords_by_state=[[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0)]],
    )

    ensemble_spread_view(port, "obj")

    selections = [args[0] for args, _ in port.calls("alter") if args]
    assert selections, "no alter issued"
    for sel in selections:
        assert "chain  and" not in sel, f"unscoped selection: {sel}"
    assert any('chain "" and resi "-3"' in s for s in selections), selections
