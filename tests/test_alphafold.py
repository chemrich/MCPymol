"""Tests for AlphaFold DB fetching and pLDDT confidence colouring."""

import json
import os
import urllib.error
import urllib.request
from unittest.mock import MagicMock, patch

import pytest

from mcpymol.structures import (
    _alphafold_accession,
    _download_alphafold,
    _resolve_alphafold_url,
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
    mock_af.assert_called_once_with(uniprot_id="P69905", obj_name=None, replace=True)
    mock_sr.assert_not_called()  # never touched the RCSB


@patch("mcpymol.structures.fetch_alphafold")
@patch("mcpymol.structures.send_request")
def test_fetch_structure_still_uses_the_pdb_for_pdb_codes(mock_sr, mock_af):
    mock_sr.side_effect = _sr_ok

    fetch_structure(pdb_code="1ubq")

    mock_af.assert_not_called()
    assert "fetch" in [c.args[0] for c in mock_sr.call_args_list]


# ── URL resolution ───────────────────────────────────────────────────────────
#
# The filename is asked for, not constructed. AlphaFold DB renumbers its model
# version periodically and *removes* the old files (v4 was current in 2024, v6
# by 2026), and some entries — the SARS-CoV-2 proteome among them — are keyed
# by an internal numeric ID rather than the accession. A hardcoded
# AF-{acc}-F1-model_v4 shipped in v1.3.0 and resolved for nothing at all.


def _api_response(*urls):
    """urlopen stand-in returning an AlphaFold prediction API payload."""
    payload = json.dumps([{"cifUrl": u} for u in urls]).encode()
    resp = MagicMock()
    resp.read.return_value = payload
    ctx = MagicMock()
    ctx.__enter__.return_value = resp
    return ctx


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_asks_the_api_for_the_url(mock_open):
    mock_open.return_value = _api_response(
        "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif"
    )

    url, err = _resolve_alphafold_url("P69905")

    assert err == ""
    assert url.endswith("AF-P69905-F1-model_v6.cif")
    assert "api/prediction/P69905" in mock_open.call_args[0][0]


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_handles_a_non_accession_entry_id(mock_open):
    """SARS-CoV-2 spike (P0DTC2) is served under an internal numeric ID, so no
    accession-based filename exists for it in any version."""
    mock_open.return_value = _api_response(
        "https://alphafold.ebi.ac.uk/files/AF-0000000365840314-model_v1.cif"
    )

    url, err = _resolve_alphafold_url("P0DTC2")

    assert err == ""
    assert "0000000365840314" in url


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_picks_the_requested_fragment(mock_open):
    mock_open.return_value = _api_response("first.cif", "second.cif", "third.cif")

    url, _err = _resolve_alphafold_url("P12345", fragment=2)

    assert url == "second.cif"


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_reports_a_missing_fragment(mock_open):
    mock_open.return_value = _api_response("only.cif")

    url, err = _resolve_alphafold_url("P12345", fragment=3)

    assert url is None
    assert "1 fragment(s)" in err


@pytest.mark.parametrize("code", [400, 404, 422])
@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_explains_an_unknown_accession(mock_open, code):
    """400 is a malformed accession, 404/422 a valid one with no model; neither
    is worth surfacing as a raw HTTP code."""
    mock_open.side_effect = urllib.error.HTTPError("u", code, "nope", {}, None)

    url, err = _resolve_alphafold_url("NOPE")

    assert url is None
    assert "no prediction for 'NOPE'" in err


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_reports_other_http_errors(mock_open):
    mock_open.side_effect = urllib.error.HTTPError("u", 503, "Service Unavailable", {}, None)

    url, err = _resolve_alphafold_url("P69905")

    assert url is None
    assert "HTTP 503" in err


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_reports_an_unreachable_api(mock_open):
    mock_open.side_effect = urllib.error.URLError("offline")

    url, err = _resolve_alphafold_url("P69905")

    assert url is None
    assert "Could not reach AlphaFold DB" in err


@patch("mcpymol.structures.urllib.request.urlopen")
def test_resolve_handles_an_empty_prediction_list(mock_open):
    mock_open.return_value = _api_response()

    url, err = _resolve_alphafold_url("P69905")

    assert url is None
    assert "no prediction" in err


# ── download ─────────────────────────────────────────────────────────────────


@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_writes_a_temp_file(mock_open, mock_resolve):
    mock_resolve.return_value = ("https://example/AF-P69905-F1-model_v6.cif", "")
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"

    path, err = _download_alphafold("P69905")

    assert err == ""
    with open(path, "rb") as fh:
        assert fh.read() == b"data_AF\n"
    os.unlink(path)


@patch("mcpymol.structures._resolve_alphafold_url")
def test_download_propagates_a_resolution_failure(mock_resolve):
    mock_resolve.return_value = (None, "AlphaFold DB has no prediction for 'X'.")

    path, err = _download_alphafold("X")

    assert path is None
    assert "no prediction" in err


@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_skips_the_api_when_a_version_is_pinned(mock_open, mock_resolve):
    """An explicit model_version is an override, so it must not be silently
    replaced by whatever the database currently serves."""
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"

    path, err = _download_alphafold("P69905", version=4)

    assert err == ""
    mock_resolve.assert_not_called()
    assert "model_v4" in mock_open.call_args[0][0]
    os.unlink(path)


@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_explains_a_retired_version(mock_open, mock_resolve):
    mock_resolve.return_value = ("https://example/model.cif", "")
    mock_open.side_effect = urllib.error.HTTPError("u", 404, "Not Found", {}, None)

    path, err = _download_alphafold("P69905")

    assert path is None
    assert "retired" in err


@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_download_rejects_an_empty_body(mock_open, mock_resolve):
    mock_resolve.return_value = ("https://example/model.cif", "")
    mock_open.return_value.__enter__.return_value.read.return_value = b""

    path, err = _download_alphafold("P69905")

    assert path is None
    assert "empty file" in err


# ── fetch_alphafold ──────────────────────────────────────────────────────────


@patch("mcpymol.views.send_request")
@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_loads_and_colours_by_plddt(
    mock_open, mock_resolve, mock_sr, mock_views_sr
):
    mock_resolve.return_value = ("https://example/model.cif", "")
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"
    mock_sr.side_effect = _sr_ok
    mock_views_sr.side_effect = _sr_with_bfactors([95.0, 92.0, 60.0])

    result = fetch_alphafold(uniprot_id="P69905")

    assert "Fetched AlphaFold model for P69905, fragment 1, as 'AF_P69905'" in result
    assert "pLDDT confidence" in result  # plddt_view ran
    assert "load" in [c.args[0] for c in mock_sr.call_args_list]


@patch("mcpymol.views.send_request")
@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_cleans_up_the_download(mock_open, mock_resolve, mock_sr, mock_views_sr):
    mock_resolve.return_value = ("https://example/model.cif", "")
    import os

    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"
    mock_sr.side_effect = _sr_ok
    mock_views_sr.side_effect = _sr_with_bfactors([90.0])

    fetch_alphafold(uniprot_id="P69905")

    loaded = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "load"]
    assert loaded and not os.path.exists(loaded[0]), "temp CIF left on disk"


@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_propagates_download_failure(mock_open, mock_resolve, mock_sr):
    mock_resolve.return_value = (None, "AlphaFold DB has no prediction for 'P00000'.")

    result = fetch_alphafold(uniprot_id="P00000")

    assert "no prediction for 'P00000'" in result
    mock_sr.assert_not_called()


@patch("mcpymol.structures.send_request")
@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_reports_a_load_failure(mock_open, mock_resolve, mock_sr):
    mock_resolve.return_value = ("https://example/model.cif", "")
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
@patch("mcpymol.structures._resolve_alphafold_url")
@patch("mcpymol.structures.urllib.request.urlopen")
def test_fetch_alphafold_honours_version_and_fragment(
    mock_open, mock_resolve, mock_sr, mock_views_sr
):
    mock_open.return_value.__enter__.return_value.read.return_value = b"data_AF\n"
    mock_sr.side_effect = _sr_ok
    mock_views_sr.side_effect = _sr_with_bfactors([90.0])

    result = fetch_alphafold(uniprot_id="P69905", model_version=3, fragment=2)

    url = mock_open.call_args[0][0]
    assert "F2" in url and "model_v3" in url
    assert "fragment 2 (pinned v3)" in result
    mock_resolve.assert_not_called()  # an explicit pin bypasses the API


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


# ── live contract ────────────────────────────────────────────────────────────


@pytest.mark.network
def test_alphafold_api_contract_is_still_what_we_assume():
    """Hits AlphaFold DB for real.

    Every other test here mocks urlopen, which is why v1.3.0 shipped a
    hardcoded model_v4 URL that resolved for nothing: the mocks proved the URL
    was *built* correctly, never that it *existed*. Opt-in via
    `pytest -m network` so CI stays offline-safe, but run it before a release.
    """
    url, err = _resolve_alphafold_url("P69905")

    assert err == "", err
    assert url and url.startswith("https://alphafold.ebi.ac.uk/files/")
    assert url.endswith((".cif", ".pdb"))

    with urllib.request.urlopen(url, timeout=60) as resp:
        head = resp.read(512)
    assert head, "model file was empty"
