"""Tests for evolutionary conservation analysis (conservation_view + helpers)."""

import json
import math
import pytest
from unittest.mock import patch, MagicMock

import mcpymol.server as server_module
from mcpymol.server import (
    _parse_a3m,
    _compute_shannon_entropy,
    _run_mmseqs2,
    conservation_view,
    _AA_ALPHABET,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _sr_mock(**action_results):
    """Return a side_effect for patching send_request."""
    def fake(action, args=None, kwargs=None):
        if action in action_results:
            val = action_results[action]
            result = val(action, args, kwargs) if callable(val) else val
            return {"status": "success", "result": result}
        return {"status": "success", "result": "OK"}
    return fake


def _actions(mock_sr):
    """Return the ordered list of action names sent via a patched send_request."""
    return [c.args[0] for c in mock_sr.call_args_list]


# ── Sample A3M data ─────────────────────────────────────────────────────────

SAMPLE_A3M = """\
>query
MKFLILLFNILCR
>hit1
MKFLILLFNILCR
>hit2
MKFLILLFNILCR
>hit3
MKFLILLFNILCR
>hit4
MKYLILLFNILCR
>hit5
MKFLVLLFNILCR
"""

# Same query but with A3M lowercase insertions in hits
SAMPLE_A3M_WITH_INSERTIONS = """\
>query
ACDEF
>hit1
ACDaEF
>hit2
AcCDEF
>hit3
ACDEY
"""


# ── _parse_a3m tests ────────────────────────────────────────────────────────

def test_parse_a3m_basic():
    """Parses sequences and strips lowercase insertions."""
    msa = _parse_a3m(SAMPLE_A3M)
    assert len(msa) == 6  # query + 5 hits
    # Query sequence preserved exactly
    assert "".join(msa[0]) == "MKFLILLFNILCR"
    # All sequences should be the same length as query (insertions stripped)
    query_len = len(msa[0])
    for seq in msa:
        assert len(seq) == query_len


def test_parse_a3m_strips_insertions():
    """Lowercase characters (insertions) are removed."""
    msa = _parse_a3m(SAMPLE_A3M_WITH_INSERTIONS)
    assert len(msa) == 4
    # Query is ACDEF (5 residues)
    assert "".join(msa[0]) == "ACDEF"
    # hit1 "ACDaEF" → strip 'a' → "ACDEF" (5 residues)
    assert "".join(msa[1]) == "ACDEF"
    # hit2 "AcCDEF" → strip 'c' → "ACDEF" (5 residues)
    assert "".join(msa[2]) == "ACDEF"
    # hit3 "ACDEY" → no insertions → "ACDEY"
    assert "".join(msa[3]) == "ACDEY"


def test_parse_a3m_empty():
    """Empty input returns empty list."""
    assert _parse_a3m("") == []


def test_parse_a3m_single_sequence():
    """Single sequence (query only) is handled."""
    msa = _parse_a3m(">query\nACDEF\n")
    assert len(msa) == 1
    assert "".join(msa[0]) == "ACDEF"


# ── _compute_shannon_entropy tests ──────────────────────────────────────────

def test_entropy_perfectly_conserved():
    """A column where every sequence has the same residue → entropy 0."""
    # 10 sequences, all with 'A' at position 0
    msa = [["A"] for _ in range(10)]
    entropies = _compute_shannon_entropy(msa)
    assert len(entropies) == 1
    assert entropies[0] == pytest.approx(0.0, abs=1e-6)


def test_entropy_maximally_variable():
    """A column with uniform distribution over all 20 AAs → entropy ~1.0."""
    aas = sorted(list(_AA_ALPHABET))
    # Each AA appears exactly once
    msa = [[aa] for aa in aas]
    entropies = _compute_shannon_entropy(msa)
    assert len(entropies) == 1
    assert entropies[0] == pytest.approx(1.0, abs=0.01)


def test_entropy_two_residues_equal():
    """50/50 split between two AAs should give known entropy."""
    msa = [["A"], ["A"], ["G"], ["G"]]
    entropies = _compute_shannon_entropy(msa)
    expected = 1.0 / math.log2(20)  # log2(2) / log2(20)
    assert entropies[0] == pytest.approx(expected, abs=1e-4)


def test_entropy_multiple_positions():
    """Multiple columns produce one entropy value each."""
    msa = [
        ["A", "A"],
        ["A", "G"],
        ["A", "C"],
        ["A", "D"],
    ]
    entropies = _compute_shannon_entropy(msa)
    assert len(entropies) == 2
    # Position 0: all A → entropy 0
    assert entropies[0] == pytest.approx(0.0, abs=1e-6)
    # Position 1: 4 different AAs → some entropy
    assert entropies[1] > 0


def test_entropy_gaps_ignored():
    """Gap characters ('-') are not counted as amino acids."""
    msa = [
        ["A"],
        ["-"],
        ["-"],
        ["A"],
    ]
    entropies = _compute_shannon_entropy(msa)
    # Only A is counted → perfectly conserved
    assert entropies[0] == pytest.approx(0.0, abs=1e-6)


def test_entropy_empty_msa():
    """Empty MSA returns empty list."""
    assert _compute_shannon_entropy([]) == []


# ── _run_mmseqs2 tests ──────────────────────────────────────────────────────

@patch("mcpymol.server.urllib.request.urlopen")
def test_run_mmseqs2_success(mock_urlopen):
    """Successful submit → poll → download cycle."""
    import io
    import tarfile

    # Create a mock tar.gz with an .a3m file
    a3m_content = ">query\nACDEF\n>hit1\nACDEY\n"
    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode="w:gz") as tar:
        info = tarfile.TarInfo(name="result.a3m")
        data = a3m_content.encode()
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    tar_bytes = tar_buf.getvalue()

    # Mock three calls: submit, poll, download
    submit_resp = MagicMock()
    submit_resp.read.return_value = json.dumps({"id": "test-ticket-123"}).encode()
    submit_resp.__enter__ = MagicMock(return_value=submit_resp)
    submit_resp.__exit__ = MagicMock(return_value=False)

    poll_resp = MagicMock()
    poll_resp.read.return_value = json.dumps({"status": "COMPLETE"}).encode()
    poll_resp.__enter__ = MagicMock(return_value=poll_resp)
    poll_resp.__exit__ = MagicMock(return_value=False)

    dl_resp = MagicMock()
    dl_resp.read.return_value = tar_bytes
    dl_resp.__enter__ = MagicMock(return_value=dl_resp)
    dl_resp.__exit__ = MagicMock(return_value=False)

    mock_urlopen.side_effect = [submit_resp, poll_resp, dl_resp]

    result = _run_mmseqs2("ACDEF", server_url="https://test.example.com")
    assert ">query" in result
    assert "ACDEF" in result


@patch("mcpymol.server.urllib.request.urlopen")
def test_run_mmseqs2_error_status(mock_urlopen):
    """MMseqs2 returning ERROR status raises RuntimeError."""
    submit_resp = MagicMock()
    submit_resp.read.return_value = json.dumps({"id": "ticket-456"}).encode()
    submit_resp.__enter__ = MagicMock(return_value=submit_resp)
    submit_resp.__exit__ = MagicMock(return_value=False)

    error_resp = MagicMock()
    error_resp.read.return_value = json.dumps({"status": "ERROR", "message": "DB not found"}).encode()
    error_resp.__enter__ = MagicMock(return_value=error_resp)
    error_resp.__exit__ = MagicMock(return_value=False)

    mock_urlopen.side_effect = [submit_resp, error_resp]

    with pytest.raises(RuntimeError, match="MMseqs2 search failed"):
        _run_mmseqs2("ACDEF")


# ── conservation_view integration tests ─────────────────────────────────────

@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_success(mock_sr, mock_mmseqs):
    """Full pipeline: get chain → get FASTA → run mmseqs2 → color."""
    # Mock send_request: return chain A and a FASTA sequence
    def fake_sr(action, args=None, kwargs=None):
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        if action == "get_fastastr":
            return {"status": "success", "result": ">chain_A\nMKFLILLFNILCRGSG\n"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake_sr

    # Mock mmseqs2 to return a small MSA
    mock_mmseqs.return_value = (
        ">query\nMKFLILLFNILCRGSG\n"
        ">hit1\nMKFLILLFNILCRGSG\n"
        ">hit2\nMKYLILLFNILCRGSG\n"
        ">hit3\nMKFLVLLFNILCRGSG\n"
    )

    result = conservation_view(obj_name="1ubq")

    assert "conservation view" in result.lower()
    assert "chain A" in result
    assert "4 sequences" in result
    assert "relative scale" in result  # default is relative
    assert "Entropy range" in result

    # Verify key actions were performed
    acts = _actions(mock_sr)
    assert "get_chains" in acts
    assert "get_fastastr" in acts
    assert "do" in acts  # B-factor alteration, spectrum, etc.


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_specific_chain(mock_sr, mock_mmseqs):
    """When chain is specified, skip get_chains and use that chain directly."""
    def fake_sr(action, args=None, kwargs=None):
        if action == "get_fastastr":
            return {"status": "success", "result": ">chain_B\nACDEFGHIKL\n"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake_sr

    mock_mmseqs.return_value = ">query\nACDEFGHIKL\n>hit1\nACDEFGHIKL\n>hit2\nACYEFGHIKL\n"

    result = conservation_view(obj_name="1ubq", chain="B")

    assert "chain B" in result
    # get_chains should NOT have been called
    acts = _actions(mock_sr)
    assert acts[0] != "get_chains"


@patch("mcpymol.server.send_request")
def test_conservation_view_no_chains(mock_sr):
    """Error when no protein chains found."""
    mock_sr.side_effect = _sr_mock(get_chains=[])

    result = conservation_view(obj_name="empty")

    assert "Error" in result
    assert "chains" in result.lower()


@patch("mcpymol.server.send_request")
def test_conservation_view_short_sequence(mock_sr):
    """Sequence shorter than 10 residues is rejected."""
    def fake_sr(action, args=None, kwargs=None):
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        if action == "get_fastastr":
            return {"status": "success", "result": ">chain_A\nACDE\n"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake_sr

    result = conservation_view(obj_name="tiny")

    assert "too short" in result.lower()


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_mmseqs_error(mock_sr, mock_mmseqs):
    """MMseqs2 failure is propagated as a user-readable error."""
    def fake_sr(action, args=None, kwargs=None):
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        if action == "get_fastastr":
            return {"status": "success", "result": ">chain_A\nMKFLILLFNILCRGSG\n"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake_sr

    mock_mmseqs.side_effect = RuntimeError("Server unreachable")

    result = conservation_view(obj_name="1ubq")

    assert "Error running MMseqs2" in result
    assert "Server unreachable" in result


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_insufficient_msa(mock_sr, mock_mmseqs):
    """MSA with only the query sequence produces a warning."""
    def fake_sr(action, args=None, kwargs=None):
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        if action == "get_fastastr":
            return {"status": "success", "result": ">chain_A\nMKFLILLFNILCRGSG\n"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake_sr

    mock_mmseqs.return_value = ">query\nMKFLILLFNILCRGSG\n"

    result = conservation_view(obj_name="1ubq")

    assert "Warning" in result or "1 sequence" in result


# ── Caching tests ───────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def clear_conservation_cache():
    """Clear the in-memory cache before each test to ensure isolation."""
    server_module._conservation_cache.clear()
    yield
    server_module._conservation_cache.clear()


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_cache_hit_skips_api(mock_sr, mock_mmseqs):
    """Second call for the same sequence does not call MMseqs2 again."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    # First call — should hit the API
    conservation_view(obj_name="1ubq")
    assert mock_mmseqs.call_count == 1

    # Reset send_request mock for second call
    mock_sr.side_effect = _conservation_sr_mock()

    # Second call with same protein — should use cache
    conservation_view(obj_name="1ubq")
    assert mock_mmseqs.call_count == 1  # still 1, not called again


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_cache_hit_reported_in_message(mock_sr, mock_mmseqs):
    """Return message says 'cached' on a cache hit."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    conservation_view(obj_name="1ubq")
    mock_sr.side_effect = _conservation_sr_mock()
    result = conservation_view(obj_name="1ubq")

    assert "cached" in result


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_force_refresh_bypasses_cache(mock_sr, mock_mmseqs):
    """force_refresh=True re-fetches even when the cache has an entry."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    conservation_view(obj_name="1ubq")
    assert mock_mmseqs.call_count == 1

    mock_sr.side_effect = _conservation_sr_mock()
    conservation_view(obj_name="1ubq", force_refresh=True)
    assert mock_mmseqs.call_count == 2  # called again


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_cache_keyed_by_sequence_not_name(mock_sr, mock_mmseqs):
    """Same sequence loaded under a different object name hits the cache."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    conservation_view(obj_name="1ubq")
    assert mock_mmseqs.call_count == 1

    mock_sr.side_effect = _conservation_sr_mock()
    conservation_view(obj_name="ubiquitin")  # different name, same sequence
    assert mock_mmseqs.call_count == 1  # cache hit


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_scale_change_uses_cache(mock_sr, mock_mmseqs):
    """Changing scale from relative to absolute reuses cached entropies."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    conservation_view(obj_name="1ubq", scale="relative")
    assert mock_mmseqs.call_count == 1

    mock_sr.side_effect = _conservation_sr_mock()
    result = conservation_view(obj_name="1ubq", scale="absolute")
    assert mock_mmseqs.call_count == 1  # no second API call
    assert "absolute scale" in result
    assert "cached" in result


# ── Scaling mode tests ──────────────────────────────────────────────────────

def _conservation_sr_mock():
    """Shared send_request side_effect for scaling tests."""
    def fake_sr(action, args=None, kwargs=None):
        if action == "get_chains":
            return {"status": "success", "result": ["A"]}
        if action == "get_fastastr":
            return {"status": "success", "result": ">chain_A\nMKFLILLFNILCRGSG\n"}
        return {"status": "success", "result": "OK"}
    return fake_sr


# MSA where positions 1-2 are perfectly conserved, position 3 varies
_SCALING_MSA = (
    ">query\nMKFLILLFNILCRGSG\n"
    ">hit1\nMKFLILLFNILCRGSG\n"
    ">hit2\nMKYLILLFNILCRGSG\n"
    ">hit3\nMKFLVLLFNILCRGSG\n"
)


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_relative_scale_b_factors(mock_sr, mock_mmseqs):
    """Relative scaling: most conserved position → 100, most variable → 0."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    result = conservation_view(obj_name="1ubq", scale="relative")

    assert "relative scale" in result

    # Collect all B-factor alter calls to inspect the scores
    b_scores = {}
    for call in mock_sr.call_args_list:
        if call.args[0] == "do":
            cmd_str = call.kwargs.get("args", call.args[1] if len(call.args) > 1 else [""])[0]
            if "name CA, b=" in cmd_str:
                # Extract resi and b-value
                import re
                m = re.search(r"resi (\d+) and name CA, b=([\d.]+)", cmd_str)
                if m:
                    b_scores[int(m.group(1))] = float(m.group(2))

    if b_scores:
        scores = list(b_scores.values())
        # In relative mode, the most conserved should be 100.0 and most variable 0.0
        assert max(scores) == pytest.approx(100.0, abs=0.1)
        assert min(scores) == pytest.approx(0.0, abs=0.1)


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_absolute_scale_b_factors(mock_sr, mock_mmseqs):
    """Absolute scaling: scores are (1 - entropy) * 100 without rescaling."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    result = conservation_view(obj_name="1ubq", scale="absolute")

    assert "absolute scale" in result

    # Collect B-factor alter calls
    b_scores = {}
    for call in mock_sr.call_args_list:
        if call.args[0] == "do":
            cmd_str = call.kwargs.get("args", call.args[1] if len(call.args) > 1 else [""])[0]
            if "name CA, b=" in cmd_str:
                import re
                m = re.search(r"resi (\d+) and name CA, b=([\d.]+)", cmd_str)
                if m:
                    b_scores[int(m.group(1))] = float(m.group(2))

    if b_scores:
        scores = list(b_scores.values())
        # In absolute mode, perfectly conserved positions should be 100.0
        # but variable positions should NOT be 0.0 (they'd only be 0 at max theoretical entropy)
        assert max(scores) == pytest.approx(100.0, abs=0.1)
        # The most variable position has some entropy but nowhere near log2(20),
        # so its absolute score should be well above 0
        assert min(scores) > 10.0


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_default_is_relative(mock_sr, mock_mmseqs):
    """Default scale parameter is relative."""
    mock_sr.side_effect = _conservation_sr_mock()
    mock_mmseqs.return_value = _SCALING_MSA

    # Call without specifying scale
    result = conservation_view(obj_name="1ubq")

    assert "relative scale" in result


@patch("mcpymol.server._run_mmseqs2")
@patch("mcpymol.server.send_request")
def test_conservation_view_uniform_entropy_falls_back(mock_sr, mock_mmseqs):
    """When all positions have identical entropy, relative mode degrades gracefully."""
    mock_sr.side_effect = _conservation_sr_mock()
    # All sequences identical → all entropies are 0 → entropy_range is 0
    mock_mmseqs.return_value = (
        ">query\nMKFLILLFNILCRGSG\n"
        ">hit1\nMKFLILLFNILCRGSG\n"
        ">hit2\nMKFLILLFNILCRGSG\n"
    )

    result = conservation_view(obj_name="1ubq", scale="relative")

    # Should not crash; falls back to absolute when range is 0
    assert "conservation view" in result.lower()
    assert "relative scale" in result
