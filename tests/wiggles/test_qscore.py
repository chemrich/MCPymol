"""Tests for Q-score parsing and the view built on it."""

from __future__ import annotations

import gzip

import pytest

from mcpymol.wiggles.bfactors import has_stash
from mcpymol.wiggles.port import ITERATE_TO_LIST, FakePort, PortError
from mcpymol.wiggles.qscore import NO_DATA_COLOUR, QSCORE_LEGEND, parse_validation_xml, qscore_view

VALIDATION = """<?xml version="1.0" encoding="UTF-8"?>
<wwPDB-validation-information>
  <Entry pdbid="7abc"/>
  <ModelledSubgroup chain="A" resnum="1" resname="MET" Q_score="0.92"/>
  <ModelledSubgroup chain="A" resnum="2" resname="SER" Q_score="0.41"/>
  <ModelledSubgroup chain="B" resnum="1" resname="GLY" Q_score="0.77"/>
  <ModelledSubgroup chain="A" resnum="3" resname="HOH"/>
</wwPDB-validation-information>
"""

ROWS = [
    ("A", "1", "MET", "CA", "", 1.0, 20.0),
    ("A", "2", "SER", "CA", "", 1.0, 30.0),
    ("B", "1", "GLY", "CA", "", 1.0, 25.0),
]


@pytest.fixture
def validation(tmp_path):
    p = tmp_path / "7abc_validation.xml"
    p.write_text(VALIDATION)
    return p


# -- parsing ---------------------------------------------------------------


def test_parses_per_residue_scores(validation):
    scores = parse_validation_xml(validation)
    assert scores[("A", "1")] == pytest.approx(0.92)
    assert scores[("A", "2")] == pytest.approx(0.41)
    assert scores[("B", "1")] == pytest.approx(0.77)


def test_residue_without_a_score_is_omitted_not_defaulted(validation):
    """Absence must stay distinguishable from a low score."""
    scores = parse_validation_xml(validation)
    assert ("A", "3") not in scores
    assert len(scores) == 3


def test_gzipped_validation_file(tmp_path):
    p = tmp_path / "v.xml.gz"
    p.write_bytes(gzip.compress(VALIDATION.encode()))
    assert len(parse_validation_xml(p)) == 3


def test_alternate_attribute_spellings(tmp_path):
    p = tmp_path / "v.xml"
    p.write_text('<r><ModelledSubgroup chain="A" resnum="1" q_score="0.5"/></r>')
    assert parse_validation_xml(p)[("A", "1")] == pytest.approx(0.5)


def test_xray_validation_file_gives_a_useful_error(tmp_path):
    """No Q-score is the normal case for non-3DEM entries; say why."""
    p = tmp_path / "v.xml"
    p.write_text('<r><ModelledSubgroup chain="A" resnum="1" resname="MET"/></r>')
    with pytest.raises(ValueError, match="3DEM"):
        parse_validation_xml(p)


def test_malformed_xml_raises(tmp_path):
    p = tmp_path / "v.xml"
    p.write_text("<not closed")
    with pytest.raises(ValueError, match="not valid XML"):
        parse_validation_xml(p)


# -- the view --------------------------------------------------------------


def test_colours_scored_residues_and_reports_stats(validation):
    port = FakePort({"iterate_to_list": ROWS})
    out = qscore_view(port, "obj", validation)
    assert "3 of 3 residues scored" in out
    assert "0.41" in out and "0.92" in out
    assert port.queried("spectrum")
    assert port.calls("spectrum")[0][1] == {"minimum": 0, "maximum": 1}
    assert QSCORE_LEGEND in out


def test_unscored_residues_are_greyed_not_zeroed(validation):
    """The correctness point: absent data is not a bad score."""
    rows = [*ROWS, ("A", "99", "HOH", "O", "", 1.0, 40.0)]
    port = FakePort({"iterate_to_list": rows})
    out = qscore_view(port, "obj", validation)

    assert "1 unscored" in out
    assert f"{NO_DATA_COLOUR}" in out
    assert "deliberately off the scale rather than scored 0" in out
    # the unscored residue must be coloured, not left on the spectrum
    assert any(a[0] == NO_DATA_COLOUR for a, _ in port.calls("color")), port.call_log
    # and it must not appear in the spectrum selection
    # spectrum(expression, palette, selection) — selection is the third arg
    spectrum_sel = port.calls("spectrum")[0][0][2]
    assert "resi 99" not in spectrum_sel


def test_no_matching_residues_is_a_clear_error(validation):
    port = FakePort({"iterate_to_list": [("Z", "500", "MET", "CA", "", 1.0, 20.0)]})
    with pytest.raises(PortError, match="author numbering"):
        qscore_view(port, "obj", validation)


def test_bfactors_preserved_by_default(validation):
    port = FakePort({"iterate_to_list": ROWS})
    out = qscore_view(port, "obj", validation)
    assert has_stash("obj")
    assert "restore_bfactors" in out


def test_legend_states_it_is_a_map_model_metric(validation):
    out = qscore_view(FakePort({"iterate_to_list": ROWS}), "obj", validation)
    assert "map–MODEL agreement metric" in out
    assert "voxel size is wrong" in out  # links back to map_info


def test_worst_residues_are_surfaced(validation):
    out = qscore_view(FakePort({"iterate_to_list": ROWS}), "obj", validation)
    assert "Least resolvable" in out
    assert "A/2 0.41" in out  # the lowest scorer leads


BLANK_CHAIN_VALIDATION = """<?xml version="1.0" encoding="UTF-8"?>
<wwPDB-validation-information>
  <Entry pdbid="7abc"/>
  <ModelledSubgroup resnum="-3" resname="MET" Q_score="0.92"/>
</wwPDB-validation-information>
"""


def test_blank_chain_does_not_escape_the_object(tmp_path):
    """A validation entry with no chain attribute made every alter unscoped.

    qscore.py defaults chain to "" when neither `chain` nor `said` is present,
    and the raw interpolation `chain {chain} and resi {resi}` then left
    `chain  and resi -3`, where PyMOL takes `and` as the chain name. Verified
    on PyMOL 3.1.0: the selection matched more atoms than the named object
    held. Since restore_bfactors only ever restores the object it was given,
    every other structure in the session lost its B-factors with no way back.
    """
    path = tmp_path / "blank.xml"
    path.write_text(BLANK_CHAIN_VALIDATION)
    port = FakePort({ITERATE_TO_LIST: [("", "-3", "MET", "CA", "", 1.0, 20.0)]})

    qscore_view(port, "obj", str(path))

    selections = [args[0] for args, _ in port.calls("alter") if args]
    assert selections, "no alter issued"
    for sel in selections:
        assert "chain  and" not in sel, f"unscoped selection: {sel}"
        assert "resi -3" not in sel + " ", f"negative resi read as a range: {sel}"
    assert any('chain "" and resi "-3"' in s for s in selections), selections
