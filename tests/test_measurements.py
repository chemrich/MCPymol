"""Measurement tools must report the number, not just acknowledge the request.

cmd.distance / angle / dihedral / get_area all return the quantity they
measure. The wrappers used to discard it, so asking MCPymol to measure
something produced a drawing and no answer.
"""

from unittest.mock import patch

import pytest

from mcpymol.bridge import format_measurement
from mcpymol.primitives import angle, count_atoms, dihedral, distance, rms_cur, sasa


def _returns(value):
    def fake(action, args=None, kwargs=None, **_ignored):
        return {"status": "success", "result": value}

    return fake


# ── format_measurement ───────────────────────────────────────────────────────


def test_format_measurement_reports_the_number():
    assert format_measurement("Distance", 3.14159, "A") == "Distance: 3.14 A."


def test_format_measurement_names_the_object():
    out = format_measurement("Distance", 2.0, "A", "d1")
    assert "stored as 'd1'" in out
    assert "2.00 A" in out


def test_format_measurement_accepts_a_numeric_string():
    assert "1.50 A" in format_measurement("Distance", "1.5", "A")


def test_format_measurement_explains_pymols_negative_sentinel():
    """PyMOL returns -1 when nothing matched; reporting that as a measurement
    of -1 A would be actively misleading."""
    out = format_measurement("Distance", -1.0, "A")

    assert "failed" in out
    assert "matched no atoms" in out
    assert "-1.00 A" not in out


def test_format_measurement_handles_a_non_numeric_result():
    out = format_measurement("Angle", "Executed angle successfully.", "degrees")

    assert "no numeric value" in out


def test_format_measurement_handles_none():
    assert "no numeric value" in format_measurement("Angle", None, "degrees")


# ── the tools ────────────────────────────────────────────────────────────────


@patch("mcpymol.primitives.send_request")
def test_distance_returns_the_distance(mock_sr):
    mock_sr.side_effect = _returns(3.42)

    result = distance(name="d1", selection1="A/10/CA", selection2="A/20/CA")

    assert "3.42 A" in result
    assert "stored as 'd1'" in result


@patch("mcpymol.primitives.send_request")
def test_distance_reports_a_failed_measurement(mock_sr):
    mock_sr.side_effect = _returns(-1.0)

    assert "failed" in distance(name="d1", selection1="nope", selection2="also_nope")


@patch("mcpymol.primitives.send_request")
def test_distance_propagates_errors(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "invalid selection"}

    assert distance(name="d1", selection1="x", selection2="y") == "invalid selection"


@patch("mcpymol.bridge.send_request")
def test_angle_returns_degrees(mock_sr):
    mock_sr.side_effect = _returns(109.47)

    result = angle(name="a1", selection1="s1", selection2="s2", selection3="s3")

    assert "109.47 degrees" in result
    assert "stored as 'a1'" in result


@patch("mcpymol.bridge.send_request")
def test_dihedral_returns_a_negative_torsion_as_a_measurement(mock_sr):
    """-57.8 deg is an ordinary alpha-helical phi, not PyMOL's -1 failure
    sentinel. Signed quantities must not be swallowed by the failure check."""
    mock_sr.side_effect = _returns(-57.8)

    result = dihedral(name="phi", selection1="a", selection2="b", selection3="c", selection4="d")

    assert "-57.80 degrees" in result
    assert "failed" not in result


def test_format_measurement_signed_allows_negatives():
    assert "-57.80 degrees" in format_measurement("Dihedral", -57.8, "degrees", signed=True)


def test_format_measurement_unsigned_still_catches_the_sentinel():
    """Distances and areas cannot be negative, so -1 is unambiguous there."""
    assert "failed" in format_measurement("Distance", -1.0, "A", signed=False)


@patch("mcpymol.bridge.send_request")
def test_sasa_returns_area(mock_sr):
    mock_sr.side_effect = _returns(14523.75)

    result = sasa(selection="1brs and chain A")

    assert "14523.75 A^2" in result
    assert mock_sr.call_args.args[0] == "get_area"


@patch("mcpymol.bridge.send_request")
def test_rms_cur_returns_rmsd_without_moving_anything(mock_sr):
    mock_sr.side_effect = _returns(1.85)

    result = rms_cur(mobile="a", target="b")

    assert "1.85 A" in result
    assert mock_sr.call_args.args[0] == "rms_cur"


@patch("mcpymol.primitives.send_request")
def test_count_atoms_reports_the_count(mock_sr):
    mock_sr.side_effect = _returns(1231)

    assert "1,231 atoms" in count_atoms(selection="polymer")


@patch("mcpymol.primitives.send_request")
def test_count_atoms_flags_an_empty_selection(mock_sr):
    """An empty selection is otherwise invisible until the picture comes out
    blank, which is the single most common way a scene silently goes wrong."""
    mock_sr.side_effect = _returns(0)

    result = count_atoms(selection="chain Z")

    assert "matches no atoms" in result
    assert "Check the object name" in result


@patch("mcpymol.primitives.send_request")
def test_count_atoms_handles_a_non_numeric_result(mock_sr):
    mock_sr.side_effect = _returns("banana")

    assert "unexpected value" in count_atoms()


@pytest.mark.parametrize("name", ["sasa", "rms_cur", "count_atoms"])
def test_new_measurement_tools_are_registered(name):
    import asyncio

    from mcpymol.server import mcp

    assert name in {t.name for t in asyncio.run(mcp.list_tools())}
