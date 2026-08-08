"""Tests for B-factor save/restore around destructive views."""

from __future__ import annotations

import json

import pytest

from mcpymol.wiggles.atoms import Atom
from mcpymol.wiggles.bfactors import (
    clear_stash,
    has_stash,
    preservation_note,
    restore_bfactors,
    stash_bfactors,
)
from mcpymol.wiggles.port import FakePort, PortError

ATOMS = [
    Atom("A", "1", "MET", "CA", "", 1.0, 20.5),
    Atom("A", "1", "MET", "CB", "", 1.0, 22.0),
    Atom("A", "2", "SER", "CA", "A", 0.6, 31.25),
]


def test_stash_records_every_atom():
    assert stash_bfactors("obj", ATOMS) == 3
    assert has_stash("obj")


def test_stash_is_per_object():
    stash_bfactors("a", ATOMS)
    assert has_stash("a")
    assert not has_stash("b")


def test_clear_one_and_clear_all():
    stash_bfactors("a", ATOMS)
    stash_bfactors("b", ATOMS)
    clear_stash("a")
    assert not has_stash("a") and has_stash("b")
    clear_stash()
    assert not has_stash("b")


def test_restore_pushes_values_back_in_one_batch():
    """One alter, not one per atom — a 10k-atom structure would otherwise
    mean 10k round trips."""
    stash_bfactors("obj", ATOMS)
    port = FakePort()

    out = restore_bfactors(port, "obj")

    assert "Restored 3 B-factors" in out
    assert len(port.calls("alter")) == 1, port.call_log

    payload = next(c for c in port.commands if "stored.wiggles_b" in c)
    values = json.loads(payload.split("= ", 1)[1])
    # Key is chain|resi|name|alt — residue name is deliberately not part of it.
    assert values["A|1|CA|"] == 20.5
    assert values["A|1|CB|"] == 22.0
    assert values["A|2|CA|A"] == 31.25


def test_restore_key_includes_altloc():
    """Two atoms differing only by altloc must not collide."""
    atoms = [
        Atom("A", "2", "SER", "CA", "A", 0.6, 10.0),
        Atom("A", "2", "SER", "CA", "B", 0.4, 90.0),
    ]
    stash_bfactors("obj", atoms)
    port = FakePort()
    restore_bfactors(port, "obj")

    payload = next(c for c in port.commands if "stored.wiggles_b" in c)
    values = json.loads(payload.split("= ", 1)[1])
    assert len(values) == 2
    assert values["A|2|CA|A"] == 10.0
    assert values["A|2|CA|B"] == 90.0


def test_restore_falls_back_to_current_value_for_unknown_atoms():
    """The alter expression must default to `b`, not 0 — an atom added since
    the stash was taken should keep its value, not be zeroed."""
    stash_bfactors("obj", ATOMS)
    port = FakePort()
    restore_bfactors(port, "obj")

    (_, expression), _ = port.calls("alter")[0]
    assert expression.rstrip().endswith(", b)"), expression


def test_restore_without_a_stash_is_an_error():
    """Restoring from an empty stash would silently zero the column."""
    with pytest.raises(PortError, match="no saved B-factors"):
        restore_bfactors(FakePort(), "obj")


def test_restore_after_clear_is_an_error():
    stash_bfactors("obj", ATOMS)
    clear_stash("obj")
    with pytest.raises(PortError, match="no saved B-factors"):
        restore_bfactors(FakePort(), "obj")


def test_preservation_note_is_actionable():
    """A claim the user cannot act on is not worth making."""
    note = preservation_note("myobj", 42)
    assert "42" in note
    assert "restore_bfactors" in note
    assert "myobj" in note
