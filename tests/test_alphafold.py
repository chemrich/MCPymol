"""Tests for AlphaFold DB fetching and pLDDT confidence colouring."""

import urllib.error
from unittest.mock import patch

import pytest

from mcpymol.structures import (
    _alphafold_accession,
    _download_alphafold,
    fetch_alphafold,
    fetch_structure,
)
from mcpymol.views import _PLDDT_BANDS, _read_ca_bfactors, plddt_view


def _sr_ok(action, args=None, kwargs=None, **_ignored):
    return {"status": "success", "result": "OK"}


def _pdb_line(resi, bfactor):
    """A CA ATOM record with the B-factor in columns 61-66."""
    return (
        f"ATOM  {resi:>5}  CA  ALA A{resi:>4}      "
        f"11.000  12.000  13.000  1.00{bfactor:>6.2f}           C"
    )


def _sr_with_bfactors(values):
    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "get_pdbstr":
            body = "\n".join(_pdb_line(i + 1, v) for i, v in enumerate(values))
            return {"status": "success", "result": body}
        return {"status": "success", "result": "OK"}

    return fake


# ── identifier routing ───────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "ident,expected",
    [
        ("P69905", "P69905"),
        ("p69905", "P69905"),
        ("af-P69905", "P69905"),
        ("AF-P69905", "P69905"),
        ("AF-P69905-F1", "P69905"),
        ("AF-P69905-F1-model_v4", "P69905"),
        ("Q9Y6K9", "Q9Y6K9"),
        ("A0A022YWF9", "A0A022YWF9"),  # 10-character accession
    ],
)
def test_recognises_alphafold_identifiers(ident, expected):
    assert _alphafold_accession(ident) == expected


@pytest.mark.parametrize("ident", ["1abc", "1UBQ", "4hhb", "", "   ", "not-an-id", "AF-nonsense"])
def test_leaves_pdb_codes_and_junk_alone(ident):
    """A 4-character PDB code must never be mistaken for a UniProt accession."""
    assert _alphafold_accession(ident) is None


@patch("mcpymol.structures.fetch_alphafold")
@patch("mcpymol.structures.send_request")
def test_fetch_structure_routes_uniprot_to_alphafold(mock_sr, mock_af):
    mock_af.return_value = "fetched"

    result = fetch_structure(pdb_code="P69905")

    assert result == "fetched"
    mock_af.assert_called_once_with(uniprot_id="P69905", obj_name=None)
    mock_sr.assert_not_called()  # never touched the RCSB


@patch("mcpymol.structures.fetch_alphafold")
@patch("mcpymol.structures.send_request")
def test_fetch_structure_still_uses_the_pdb_for_pdb_codes(mock_sr, mock_af):
    mock_sr.side_effect = _sr_ok

    fetch_structure(pdb_code="1ubq")

    mock_af.assert_not_called()
    assert "fetch" in [c.args[0] for c in mock_sr.call_args_list]


# ── download ─────────────────────────────────────────────────────────────────


@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_writes_a_temp_file(mock_open):
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"

    path, err = _download_alphafold("P69905", 4)

    assert err == ""
    assert path is not None
    with open(path, "rb") as fh:
        assert fh.read() == b"data_AF\n"
    import os

    os.unlink(path)
    assert "P69905" in mock_open.call_args[0][0]
    assert "model_v4" in mock_open.call_args[0][0]


@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_explains_a_404(mock_open):
    """The overwhelmingly common failure: a valid-looking accession AlphaFold
    has no model for."""
    mock_open.side_effect = urllib.error.HTTPError("u", 404, "Not Found", {}, None)

    path, err = _download_alphafold("P00000", 4)

    assert path is None
    assert "No AlphaFold model for 'P00000'" in err
    assert "model_version" in err


@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_reports_other_http_errors(mock_open):
    mock_open.side_effect = urllib.error.HTTPError("u", 503, "Service Unavailable", {}, None)

    path, err = _download_alphafold("P69905", 4)

    assert path is None
    assert "HTTP 503" in err


@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_reports_unreachable_db(mock_open):
    mock_open.side_effect = urllib.error.URLError("no route to host")

    path, err = _download_alphafold("P69905", 4)

    assert path is None
    assert "Could not reach AlphaFold DB" in err


@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_rejects_an_empty_body(mock_open):
    mock_open.return_value.__enter__.return_value.read.return_value = b""

    path, err = _download_alphafold("P69905", 4)

    assert path is None
    assert "empty file" in err


# ── fetch_alphafold ──────────────────────────────────────────────────────────


@patch("mcpymol.views.send_request")
@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_loads_and_colours_by_plddt(mock_open, mock_sr, mock_views_sr):
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"
    mock_sr.side_effect = _sr_ok
    mock_views_sr.side_effect = _sr_with_bfactors([95.0, 92.0, 60.0])

    result = fetch_alphafold(uniprot_id="P69905")

    assert "Fetched AlphaFold model AF-P69905-F1 (v4) as 'AF_P69905'" in result
    assert "pLDDT confidence" in result  # plddt_view ran
    assert "load" in [c.args[0] for c in mock_sr.call_args_list]


@patch("mcpymol.views.send_request")
@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_cleans_up_the_download(mock_open, mock_sr, mock_views_sr):
    import os

    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"
    mock_sr.side_effect = _sr_ok
    mock_views_sr.side_effect = _sr_with_bfactors([90.0])

    fetch_alphafold(uniprot_id="P69905")

    loaded = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "load"]
    assert loaded and not os.path.exists(loaded[0]), "temp CIF left on disk"


@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_propagates_download_failure(mock_open, mock_sr):
    mock_open.side_effect = urllib.error.HTTPError("u", 404, "Not Found", {}, None)

    result = fetch_alphafold(uniprot_id="P00000")

    assert "No AlphaFold model" in result
    mock_sr.assert_not_called()


@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_reports_a_load_failure(mock_open, mock_sr):
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"

    def fake(action, args=None, kwargs=None, **_ignored):
        if action == "load":
            return {"status": "error", "error": "unrecognized file format"}
        return {"status": "success", "result": "OK"}

    mock_sr.side_effect = fake

    result = fetch_alphafold(uniprot_id="P69905")

    assert "Error loading AlphaFold model for P69905" in result
    assert "unrecognized file format" in result


@patch("mcpymol.views.send_request")
@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_honours_version_and_fragment(mock_open, mock_sr, mock_views_sr):
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"
    mock_sr.side_effect = _sr_ok
    mock_views_sr.side_effect = _sr_with_bfactors([90.0])

    result = fetch_alphafold(uniprot_id="P69905", model_version=3, fragment=2)

    url = mock_open.call_args[0][0]
    assert "F2" in url and "model_v3" in url
    assert "AF-P69905-F2 (v3)" in result


# ── _read_ca_bfactors ────────────────────────────────────────────────────────


@patch("mcpymol.views.send_request")
def test_read_ca_bfactors_parses_the_column(mock_sr):
    mock_sr.side_effect = _sr_with_bfactors([95.5, 42.25, 8.0])

    assert _read_ca_bfactors("AF_P69905") == [95.5, 42.25, 8.0]


@patch("mcpymol.views.send_request")
def test_read_ca_bfactors_returns_none_on_error(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no such object"}

    assert _read_ca_bfactors("nope") is None


@patch("mcpymol.views.send_request")
def test_read_ca_bfactors_ignores_non_atom_lines(mock_sr):
    mock_sr.return_value = {
        "status": "success",
        "result": "HEADER something\n" + _pdb_line(1, 77.0) + "\nEND\n",
    }

    assert _read_ca_bfactors("x") == [77.0]


# ── plddt_view ───────────────────────────────────────────────────────────────


def test_plddt_bands_use_the_official_alphafold_thresholds():
    assert [b[1] for b in _PLDDT_BANDS] == [0.0, 50.0, 70.0, 90.0]
    # Very high is AlphaFold's dark blue, very low its orange.
    assert _PLDDT_BANDS[-1][2] == (0.051, 0.341, 0.827)
    assert _PLDDT_BANDS[0][2] == (1.000, 0.490, 0.271)


@patch("mcpymol.views.send_request")
def test_plddt_view_colours_every_band(mock_sr):
    mock_sr.side_effect = _sr_with_bfactors([95.0, 80.0, 60.0, 30.0])

    plddt_view(obj_name="AF_P69905")

    colored = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "color"]
    assert colored == [b[0] for b in _PLDDT_BANDS]


@patch("mcpymol.views.send_request")
def test_plddt_view_bands_do_not_overlap(mock_sr):
    """Each selection must be bounded above by the next band's floor, or the
    later colours would overpaint the earlier ones."""
    mock_sr.side_effect = _sr_with_bfactors([95.0])

    plddt_view(obj_name="AF_X")

    sels = [c.kwargs["args"][1] for c in mock_sr.call_args_list if c.args[0] == "color"]
    assert "b > 0.0 and b < 50.0" in sels[0]
    assert "b > 50.0 and b < 70.0" in sels[1]
    assert "b > 70.0 and b < 90.0" in sels[2]
    assert sels[3].endswith("b > 90.0")  # open-ended top band


@patch("mcpymol.views.send_request")
def test_plddt_view_reports_the_confidence_breakdown(mock_sr):
    # 2 very high, 1 confident, 1 low  → 50 / 25 / 25 / 0
    mock_sr.side_effect = _sr_with_bfactors([95.0, 92.0, 75.0, 55.0])

    result = plddt_view(obj_name="AF_X")

    assert "4 residues" in result
    assert "very high (>90): 50%" in result
    assert "confident (70-90): 25%" in result
    assert "low (50-70): 25%" in result
    assert "very low (<50): 0%" in result
    assert "mean 79.2" in result


@patch("mcpymol.views.send_request")
def test_plddt_view_warns_on_a_non_alphafold_structure(mock_sr):
    """Crystallographic B-factors routinely exceed 100 — colouring them as
    confidence would be nonsense, so say so."""
    mock_sr.side_effect = _sr_with_bfactors([12.0, 340.0])

    result = plddt_view(obj_name="1ubq")

    assert "do not look like pLDDT" in result
    assert "bfactor_view" in result


@patch("mcpymol.views.send_request")
def test_plddt_view_handles_an_empty_object(mock_sr):
    mock_sr.side_effect = _sr_with_bfactors([])

    assert "no residues found" in plddt_view(obj_name="empty")


@patch("mcpymol.views.send_request")
def test_plddt_view_reports_a_read_failure(mock_sr):
    mock_sr.return_value = {"status": "error", "error": "no such object"}

    assert "could not read B-factors" in plddt_view(obj_name="nope")


def test_plddt_view_is_registered_as_a_tool():
    import asyncio

    from mcpymol.server import mcp

    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {"plddt_view", "fetch_alphafold"} <= names
