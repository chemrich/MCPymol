"""Tests for the sense-1 occupancy views.

Every one of these runs against FakePort. No PyMOL, no sockets, no network —
which is the point of SPEC.md §0.
"""

from __future__ import annotations

import json

import pytest

from mcpymol.wiggles.bfactors import has_stash, restore_bfactors
from mcpymol.wiggles.occupancy import SENSE_1_LEGEND, altloc_view, occupancy_view
from mcpymol.wiggles.port import FakePort, PortError


def atoms(*rows):
    """Build an iterate response. Row: (chain, resi, resn, name, alt, q, b)."""
    return {"iterate_to_list": list(rows)}


FULL = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "ALA", "CA", "", 1.0, 22.0),
]

PARTIAL = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "SER", "CA", "A", 0.6, 25.0),
    ("A", "2", "SER", "CA", "B", 0.4, 31.0),
]


# -- the invariant that matters most --------------------------------------


@pytest.mark.parametrize(
    "view,args",
    [(occupancy_view, ()), (altloc_view, ())],
    ids=["occupancy_view", "altloc_view"],
)
@pytest.mark.parametrize("rows", [FULL, PARTIAL], ids=["full", "partial"])
def test_every_report_names_the_sense(view, args, rows):
    """SPEC: occupancy_view and composition_view stay separate forever, and a
    reader must never have to guess which sense a render shows."""
    port = FakePort(atoms(*rows))
    out = view(port, "obj", *args)
    assert SENSE_1_LEGEND in out
    assert "SENSE 1" in out


@pytest.mark.parametrize("rows", [FULL, PARTIAL], ids=["full", "partial"])
def test_never_claims_compositional_occupancy(rows):
    """No sense-2 language may leak into a sense-1 view."""
    port = FakePort(atoms(*rows))
    out = occupancy_view(port, "obj").lower()
    # The legend explicitly disclaims these; check they only appear as denial.
    assert "fraction of imaged particles" in out  # the disclaimer
    assert "particles containing" not in out.replace("particles containing this subunit", "")


# -- occupancy_view --------------------------------------------------------


def test_fully_occupied_model_is_not_coloured():
    """A spectrum across a constant is a rainbow that means nothing."""
    port = FakePort(atoms(*FULL))
    out = occupancy_view(port, "obj")
    assert "All 2 atoms are fully occupied" in out
    assert not port.queried("spectrum"), port.call_log
    assert not port.queried("alter"), port.call_log


def test_partial_occupancy_is_coloured_on_a_fixed_scale():
    port = FakePort(atoms(*PARTIAL))
    out = occupancy_view(port, "obj")
    assert "2 of 3 atoms partially occupied" in out
    assert port.called("alter", "obj", "b=q"), port.call_log
    # Fixed 0-1 scale: occupancy is absolute, so autoscaling would make a
    # 0.9-occupied model look as partial as a 0.1-occupied one.
    assert port.calls("spectrum")[0][1] == {"minimum": 0, "maximum": 1}, port.call_log


def test_bfactors_are_preserved_by_default():
    """Assert the values actually come back, not merely that a command ran."""
    port = FakePort(atoms(*PARTIAL))
    out = occupancy_view(port, "obj")

    assert has_stash("obj")
    assert "restore_bfactors" in out

    restore = FakePort()
    restore_bfactors(restore, "obj")
    payload = next(c for c in restore.commands if "stored.wiggles_b" in c)
    saved = json.loads(payload.split("= ", 1)[1])
    # PARTIAL's original B-factors, not the occupancies written over them.
    assert sorted(saved.values()) == [20.0, 25.0, 31.0]


def test_bfactor_destruction_is_warned_when_opted_into():
    port = FakePort(atoms(*PARTIAL))
    out = occupancy_view(port, "obj", preserve_bfactors=False)
    assert not has_stash("obj")
    assert "WARNING" in out and "not preserved" in out


def test_qfit_floor_is_flagged():
    rows = [("A", "1", "SER", "CA", "A", 0.05, 20.0), ("A", "1", "SER", "CA", "B", 0.95, 20.0)]
    out = occupancy_view(FakePort(atoms(*rows)), "obj")
    assert "below q=0.1" in out
    assert "qFit removes" in out


def test_float_noise_counts_as_fully_occupied():
    """Refined occupancies are float32; exact 1.0 comparison mislabels them."""
    rows = [("A", "1", "MET", "CA", "", 0.9999, 20.0)]
    out = occupancy_view(FakePort(atoms(*rows)), "obj")
    assert "fully occupied" in out


# -- altloc_view -----------------------------------------------------------


def test_no_altlocs_says_so_without_overclaiming():
    port = FakePort(atoms(*FULL))
    out = altloc_view(port, "obj")
    assert "No alternate conformations" in out
    # Must not imply the molecule is rigid — only that this file doesn't model it.
    assert "says nothing about" in out
    assert not port.queried("color"), port.call_log


def test_altloc_groups_get_distinct_colours():
    port = FakePort(atoms(*PARTIAL))
    out = altloc_view(port, "obj")
    assert "2 altloc group(s)" in out
    shown = [a[1] for a, _ in port.calls("show")]
    assert any("alt A" in s for s in shown) and any("alt B" in s for s in shown), port.call_log
    used = {a[0] for a, _ in port.calls("color") if a[0] != "grey70"}
    assert len(used) == 2, f"expected 2 distinct colours, got {used}"


def test_altloc_view_reports_occupancies_not_just_colours():
    """Drawing every alternate at equal weight discards the measurement."""
    out = altloc_view(FakePort(atoms(*PARTIAL)), "obj")
    assert "0.6" in out and "0.4" in out
    assert "would discard the measurement" in out


def test_labels_can_be_suppressed():
    port = FakePort(atoms(*PARTIAL))
    altloc_view(port, "obj", label=False)
    assert not port.queried("label")


def test_blank_altloc_is_not_a_group():
    """A blank altloc is the absence of grouping, not a group named ''."""
    rows = [
        ("A", "1", "MET", "CA", "", 1.0, 20.0),
        ("A", "2", "SER", "CA", "A", 0.5, 20.0),
        ("A", "2", "SER", "CA", "B", 0.5, 20.0),
    ]
    out = altloc_view(FakePort(atoms(*rows)), "obj")
    assert "2 altloc group(s)" in out


# -- failure behaviour -----------------------------------------------------


def test_empty_selection_is_an_error_not_an_empty_view():
    """MCPymol issue #15: an empty result reported as success is the bug."""
    with pytest.raises(PortError, match="matched no atoms"):
        occupancy_view(FakePort({"iterate_to_list": []}), "obj")


def test_malformed_atom_row_is_an_error():
    with pytest.raises(PortError, match="malformed atom row"):
        occupancy_view(FakePort({"iterate_to_list": [("A", "1")]}), "obj")


def test_unanticipated_query_raises_rather_than_returning_none():
    with pytest.raises(KeyError):
        occupancy_view(FakePort({}), "obj")
