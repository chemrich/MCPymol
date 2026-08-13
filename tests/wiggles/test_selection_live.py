"""What a residue selection actually matches, against a running PyMOL.

    pytest -m live tests/wiggles/test_selection_live.py

WARNING: like the rest of the live suite, this drives the PyMOL session you
have open. It only adds and deletes objects under the ``_wsel_`` prefix and
touches nothing else, but ``pytest -m live`` as a whole clears the session.

Why this file exists
--------------------
``test_atoms.py`` asserted that ``residue_selection("obj", "A", "-3")``
contained ``resi "-3"`` — that the string had been quoted. It said nothing
about what PyMOL does with that string, and the answer turned out to be "reads
it as the range 1-3 anyway". The bug the quoting was added to fix was still
there, under a passing test, because the test and the bug shared a mental
model.

So this asserts **atom counts from a live session**: build an object whose
residue numbering makes the range reading and the literal reading return
different answers, then select through the real code path and count.

The structure has chain A residues -3, 1 and 2, two atoms each. Selecting
residue -3 must match 2. Under the range reading it matches all 6.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from mcpymol.wiggles.atoms import residue_clause, residue_selection
from mcpymol.wiggles.port import BridgePort, call

pytestmark = pytest.mark.live

OBJ = "_wsel_probe"

# (chain, resseq, atom name) -- two atoms per residue, and a blank-chain
# residue so the blank-chain half of quote() is covered by the same file.
ROWS = [
    ("A", -3, "N"),
    ("A", -3, "CA"),
    ("A", 1, "N"),
    ("A", 1, "CA"),
    ("A", 2, "N"),
    ("A", 2, "CA"),
    ("", 5, "N"),
    ("", 5, "CA"),
]


def _pdb() -> str:
    lines = []
    for serial, (chain, resseq, name) in enumerate(ROWS, start=1):
        name_field = f" {name:<3}"
        lines.append(
            f"{'ATOM':<6}{serial:>5} {name_field:<4}{' ':1}{'ALA':>3} "
            f"{chain:1}{resseq:>4}{' ':1}   "
            f"{float(serial):>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}          {'C':>2}"
        )
    return "\n".join(lines) + "\nEND\n"


@pytest.fixture
def probe():
    """The object, loaded into the live session and removed afterwards."""
    port = BridgePort()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "wsel.pdb"
        path.write_text(_pdb())
        call(port, "delete", OBJ)
        call(port, "load", str(path), OBJ)
        try:
            yield port
        finally:
            call(port, "delete", OBJ)


def test_the_probe_is_numbered_the_way_the_assertions_assume(probe):
    """If PyMOL ever stops reading these residue numbers as written, every
    other assertion here would pass or fail for the wrong reason."""
    assert call(probe, "count_atoms", f"({OBJ})") == 8
    assert call(probe, "count_atoms", f'({OBJ}) and chain "A"') == 6


def test_a_negative_residue_matches_only_itself(probe):
    """The whole point. Two atoms, not the six that residues -3, 1 and 2 hold
    between them."""
    sel = residue_selection(OBJ, "A", "-3")

    assert call(probe, "count_atoms", sel) == 2


def test_an_unescaped_negative_residue_really_does_over_match(probe):
    """The bug this guards against, asserted directly — so the test above
    cannot quietly start passing for some unrelated reason.

    This is what the previous implementation emitted."""
    assert call(probe, "count_atoms", f'({OBJ}) and chain "A" and resi "-3"') == 6
    assert call(probe, "count_atoms", f'({OBJ}) and chain "A" and resi -3') == 6


def test_ordinary_residues_are_unaffected(probe):
    for resi in ("1", "2"):
        assert call(probe, "count_atoms", residue_selection(OBJ, "A", resi)) == 2


def test_a_blank_chain_stays_scoped_to_the_object(probe):
    """The other half of quote(): `chain  and resi 5` would take `and` as the
    chain name and stop being scoped."""
    sel = residue_selection(OBJ, "", "5")

    assert call(probe, "count_atoms", sel) == 2


def test_the_clause_form_matches_the_same_atoms(probe):
    """residue_clause is the unscoped twin; a divergence between them is how a
    fix lands at one call site and not the other."""
    scoped = f"({OBJ}) and {residue_clause('A', '-3')}"

    assert call(probe, "count_atoms", scoped) == 2
