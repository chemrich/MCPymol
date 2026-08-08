"""Tests for the atom-reading layer shared by every tier-1 view."""

from __future__ import annotations

import pytest

from mcpymol.wiggles.atoms import (
    Atom,
    altloc_groups,
    count_states,
    fetch_atoms,
    fetch_state_coords,
    group_by_residue,
)
from mcpymol.wiggles.port import FakePort, PortError


def test_fetch_atoms_parses_rows():
    port = FakePort({"iterate_to_list": [("A", "1", "MET", "CA", "", 1.0, 20.0)]})
    (atom,) = fetch_atoms(port, "obj")
    assert atom == Atom("A", "1", "MET", "CA", "", 1.0, 20.0)
    assert atom.residue == ("A", "1")


def test_fetch_atoms_coerces_numeric_strings():
    """PyMOL round-trips through JSON; numbers may arrive as strings."""
    port = FakePort({"iterate_to_list": [("A", 1, "MET", "CA", "", "0.5", "20")]})
    (atom,) = fetch_atoms(port, "obj")
    assert atom.resi == "1"
    assert atom.q == pytest.approx(0.5)
    assert atom.b == pytest.approx(20.0)


@pytest.mark.parametrize(
    "response,match",
    [
        (None, "returned nothing"),
        ("not a list", "expected a list"),
        ([], "matched no atoms"),
    ],
)
def test_fetch_atoms_failure_modes(response, match):
    with pytest.raises(PortError, match=match):
        fetch_atoms(FakePort({"iterate_to_list": response}), "obj")


def test_wrong_field_count_names_the_expected_layout():
    with pytest.raises(PortError, match="expected 7 fields"):
        fetch_atoms(FakePort({"iterate_to_list": [("A", "1", "MET")]}), "obj")


def test_unconvertible_occupancy_is_an_error():
    port = FakePort({"iterate_to_list": [("A", "1", "MET", "CA", "", "not-a-number", 20.0)]})
    with pytest.raises(PortError, match="malformed atom row"):
        fetch_atoms(port, "obj")


# -- states ----------------------------------------------------------------


def test_count_states():
    assert count_states(FakePort({"count_states": 7}), "obj") == 7


def test_count_states_accepts_a_numeric_string():
    assert count_states(FakePort({"count_states": "3"}), "obj") == 3


def test_count_states_rejects_nonsense():
    with pytest.raises(PortError, match="count_states"):
        count_states(FakePort({"count_states": "many"}), "obj")


def test_fetch_state_coords():
    port = FakePort({"get_coords": [(1.0, 2.0, 3.0)]})
    assert fetch_state_coords(port, "obj", 1) == [(1.0, 2.0, 3.0)]


def test_fetch_state_coords_rejects_empty():
    with pytest.raises(PortError, match="no coordinates"):
        fetch_state_coords(FakePort({"get_coords": []}), "obj", 1)


def test_fetch_state_coords_rejects_malformed_point():
    with pytest.raises(PortError, match="malformed coordinate"):
        fetch_state_coords(FakePort({"get_coords": [(1.0, 2.0)]}), "obj", 1)


# -- grouping --------------------------------------------------------------


def test_group_by_residue_preserves_first_seen_order():
    atoms = [
        Atom("A", "2", "SER", "CA", "", 1.0, 0.0),
        Atom("A", "1", "MET", "CA", "", 1.0, 0.0),
        Atom("A", "2", "SER", "CB", "", 1.0, 0.0),
    ]
    grouped = group_by_residue(atoms)
    assert list(grouped) == [("A", "2"), ("A", "1")]
    assert len(grouped[("A", "2")]) == 2


def test_altloc_groups_excludes_blank_and_whitespace():
    atoms = [
        Atom("A", "1", "MET", "CA", "", 1.0, 0.0),
        Atom("A", "2", "SER", "CA", " ", 1.0, 0.0),
        Atom("A", "3", "SER", "CA", "B", 1.0, 0.0),
        Atom("A", "4", "SER", "CA", "A", 1.0, 0.0),
    ]
    assert altloc_groups(atoms) == ["A", "B"]
