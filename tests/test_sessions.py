"""Tests for save_session / load_session."""

from unittest.mock import patch

import pytest

from mcpymol.structures import load_session, save_session


def _writes_session(size=1024):
    """send_request stand-in that actually creates the .pse PyMOL was asked for."""

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "save" and args:
            with open(args[0], "wb") as fh:
                fh.write(b"\x00" * size)
        return {"status": "success", "result": "OK"}

    return fake


# ── save_session ─────────────────────────────────────────────────────────────


@patch("mcpymol.structures.send_request")
def test_save_session_writes_and_reports_size(mock_sr, tmp_path):
    mock_sr.side_effect = _writes_session(size=2_500_000)
    target = tmp_path / "figure.pse"

    result = save_session(filename=str(target))

    assert str(target) in result
    assert "2.5 MB" in result
    assert mock_sr.call_args.args[0] == "save"


@patch("mcpymol.structures.send_request")
def test_save_session_adds_the_extension(mock_sr, tmp_path):
    """Without .pse PyMOL writes bare coordinates, silently losing the scene."""
    mock_sr.side_effect = _writes_session()

    result = save_session(filename=str(tmp_path / "figure"))

    assert mock_sr.call_args.kwargs["args"][0].endswith(".pse")
    assert result.startswith("Saved session")


@pytest.mark.parametrize("name", ["scene.pse", "scene.PSE", "scene.pse.gz"])
@patch("mcpymol.structures.send_request")
def test_save_session_keeps_a_valid_extension(mock_sr, tmp_path, name):
    mock_sr.side_effect = _writes_session()

    save_session(filename=str(tmp_path / name))

    assert mock_sr.call_args.kwargs["args"][0].endswith(name)


@patch("mcpymol.structures.send_request")
def test_save_session_creates_missing_directories(mock_sr, tmp_path):
    mock_sr.side_effect = _writes_session()
    target = tmp_path / "figures" / "final" / "x.pse"

    save_session(filename=str(target))

    assert target.parent.is_dir()


@patch("mcpymol.structures.send_request")
def test_save_session_reports_pymol_errors(mock_sr, tmp_path):
    mock_sr.return_value = {"status": "error", "error": "permission denied"}

    result = save_session(filename=str(tmp_path / "x.pse"))

    assert "Error saving session" in result
    assert "permission denied" in result


@patch("mcpymol.structures.send_request")
def test_save_session_notices_a_file_that_never_appeared(mock_sr, tmp_path):
    """Success with no file means the two halves are on different machines."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = save_session(filename=str(tmp_path / "x.pse"))

    assert "no file appeared" in result
    assert "different machine" in result


# ── load_session ─────────────────────────────────────────────────────────────


@patch("mcpymol.structures.send_request")
def test_load_session_restores_and_lists_objects(mock_sr, tmp_path):
    session = tmp_path / "figure.pse"
    session.write_bytes(b"\x00")

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_object_list":
            return {"status": "success", "result": ["1ubq", "1ubq_spine"]}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = load_session(filename=str(session))

    assert "Loaded session" in result
    assert "1ubq, 1ubq_spine" in result
    load_calls = [c for c in mock_sr.call_args_list if c.args[0] == "load"]
    assert load_calls[0].kwargs["kwargs"] == {}


@patch("mcpymol.structures.send_request")
def test_load_session_can_merge(mock_sr, tmp_path):
    session = tmp_path / "figure.pse"
    session.write_bytes(b"\x00")
    mock_sr.return_value = {"status": "success", "result": []}

    result = load_session(filename=str(session), merge=True)

    load_calls = [c for c in mock_sr.call_args_list if c.args[0] == "load"]
    assert load_calls[0].kwargs["kwargs"] == {"partial": 1}
    assert "Merged session" in result


@patch("mcpymol.structures.send_request")
def test_load_session_rejects_a_missing_file(mock_sr, tmp_path):
    result = load_session(filename=str(tmp_path / "nope.pse"))

    assert "no such file" in result
    mock_sr.assert_not_called()


@patch("mcpymol.structures.send_request")
def test_load_session_rejects_a_coordinate_file(mock_sr, tmp_path):
    """Loading a PDB here would work but silently do the wrong thing —
    load_structure is what applies the prep."""
    pdb = tmp_path / "1ubq.pdb"
    pdb.write_text("ATOM\n")

    result = load_session(filename=str(pdb))

    assert "not a PyMOL session" in result
    assert "load_structure" in result
    mock_sr.assert_not_called()


@patch("mcpymol.structures.send_request")
def test_load_session_reports_pymol_errors(mock_sr, tmp_path):
    session = tmp_path / "broken.pse"
    session.write_bytes(b"\x00")
    mock_sr.return_value = {"status": "error", "error": "corrupt session file"}

    result = load_session(filename=str(session))

    assert "Error loading session" in result
    assert "corrupt session file" in result


@patch("mcpymol.structures.send_request")
def test_load_session_survives_an_object_list_failure(mock_sr, tmp_path):
    session = tmp_path / "figure.pse"
    session.write_bytes(b"\x00")

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_object_list":
            return {"status": "error", "error": "busy"}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = load_session(filename=str(session))

    assert "Loaded session" in result
    assert "no objects" in result


def test_session_tools_are_registered():
    import asyncio

    from mcpymol.server import mcp

    assert {"save_session", "load_session"} <= {t.name for t in asyncio.run(mcp.list_tools())}
