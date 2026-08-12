"""Tests for SPEC invariant I1 — provenance is never dropped, never guessed."""

from __future__ import annotations

import pytest
from test_mapinfo import write_map

from mcpymol.wiggles.mapinfo import read_map_header
from mcpymol.wiggles.maps import load_map
from mcpymol.wiggles.port import FakePort, PortError
from mcpymol.wiggles.provenance import (
    Provenance,
    declare,
    forget,
    gather_evidence,
    provenance_banner,
    provenance_of,
)


def loaded_port(obj="test"):
    """A port whose session already contains `obj` after a load."""
    return FakePort({"get_names": [obj]})


# ── the invariant ───────────────────────────────────────────────────────────


def test_unknown_is_the_default():
    """Defaulting to 'measured' would assert a generated volume was observed."""
    assert provenance_of("never-declared") is Provenance.UNKNOWN


def test_unknown_never_raises():
    """A volume nobody declared is unknown — a real answer, not an error."""
    assert provenance_of("anything") is Provenance.UNKNOWN


@pytest.mark.parametrize("p", list(Provenance))
def test_every_provenance_has_a_caveat(p):
    assert p.caveat.strip()


def test_banner_states_the_category_and_its_caveat():
    declare("obj", Provenance.GENERATED)
    banner = provenance_banner("obj")
    assert "GENERATED" in banner
    assert "No particle was reconstructed" in banner


def test_unknown_banner_says_do_not_assume_measured():
    banner = provenance_banner("undeclared")
    assert "UNKNOWN" in banner
    assert "Do not" in banner and "assume it was measured" in banner


def test_only_measured_and_sharpened_count_as_observed():
    assert Provenance.MEASURED.is_observed
    assert Provenance.SHARPENED.is_observed
    assert not Provenance.NN_ENHANCED.is_observed
    assert not Provenance.GENERATED.is_observed
    assert not Provenance.UNKNOWN.is_observed


def test_forget_one_and_forget_all():
    declare("a", Provenance.MEASURED)
    declare("b", Provenance.GENERATED)
    forget("a")
    assert provenance_of("a") is Provenance.UNKNOWN
    assert provenance_of("b") is Provenance.GENERATED
    forget()
    assert provenance_of("b") is Provenance.UNKNOWN


# ── evidence is a suggestion, not a verdict ─────────────────────────────────


def test_emdb_label_is_recognised(tmp_path):
    p = write_map(tmp_path)
    header = read_map_header(p)
    header = type(header)(**{**header.__dict__, "labels": ("::::EMDATABANK.org::::EMD-30913::::",)})
    evidence = gather_evidence(header, p)
    assert evidence.emdb_accession == "EMD-30913"
    assert any("EMD-30913" in r for r in evidence.reasons)


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("run_deepemhancer.mrc", Provenance.NN_ENHANCED),
        ("cryodrgn_vol_000.mrc", Provenance.GENERATED),
        ("job_postprocess.mrc", Provenance.SHARPENED),
        ("run_half1_unfil.mrc", Provenance.MEASURED),
    ],
)
def test_filename_tokens_suggest_a_category(tmp_path, filename, expected):
    p = write_map(tmp_path, filename)
    evidence = gather_evidence(read_map_header(p), p)
    assert evidence.suggested is expected


def test_no_signal_suggests_nothing(tmp_path):
    p = write_map(tmp_path, "aaa.mrc")
    evidence = gather_evidence(read_map_header(p), p)
    assert evidence.suggested is Provenance.UNKNOWN
    assert any("nothing in the filename" in r for r in evidence.reasons)


# ── load_map ────────────────────────────────────────────────────────────────


def test_load_map_reports_voxel_size_and_provenance(tmp_path):
    p = write_map(tmp_path, "m.mrc", voxel=0.5332)
    port = loaded_port("m")
    out = load_map(port, p, "m", provenance=Provenance.MEASURED)

    assert "0.5332" in out
    assert "MEASURED" in out
    assert port.called("load", str(p), "m"), port.call_log


def test_load_map_defaults_to_unknown_and_shows_the_evidence(tmp_path):
    """A suggestion must never be silently adopted."""
    p = write_map(tmp_path, "cryodrgn_vol_007.mrc")
    port = loaded_port("cryodrgn_vol_007")
    out = load_map(port, p)

    assert provenance_of("cryodrgn_vol_007") is Provenance.UNKNOWN
    assert "UNKNOWN" in out
    assert "Suggested: generated" in out
    assert "Not adopted" in out


def test_declared_provenance_suppresses_the_evidence_block(tmp_path):
    p = write_map(tmp_path, "cryodrgn_vol_007.mrc")
    port = loaded_port("cryodrgn_vol_007")
    out = load_map(port, p, provenance=Provenance.GENERATED)

    assert "GENERATED" in out
    assert "Suggested:" not in out


def test_provenance_accepts_a_plain_string(tmp_path):
    p = write_map(tmp_path, "m.mrc")
    out = load_map(loaded_port("m"), p, "m", provenance="nn_enhanced")
    assert "NN_ENHANCED" in out
    assert provenance_of("m") is Provenance.NN_ENHANCED


def test_geometry_warnings_survive_into_the_report(tmp_path):
    p = write_map(tmp_path, "m.mrc", nx=100, ny=100, nz=100, cella=(100.0, 100.0, 150.0))
    out = load_map(loaded_port("m"), p, "m")
    assert "ANISOTROPIC" in out


def test_bad_file_fails_before_pymol_is_touched(tmp_path):
    """A malformed map must not leave a half-loaded object behind."""
    p = tmp_path / "junk.mrc"
    p.write_bytes(b"\x00" * 100)
    port = FakePort()
    with pytest.raises(ValueError):
        load_map(port, p)
    assert port.queries == [], port.call_log


def test_object_that_does_not_arrive_is_an_error(tmp_path):
    """Issue #15: an empty load reported as success is the bug."""
    p = write_map(tmp_path, "m.mrc")
    port = FakePort({"get_names": ["something_else"]})

    with pytest.raises(PortError, match="produced no object"):
        load_map(port, p, "m")


def test_failed_load_scopes_its_cleanup_to_its_own_object(tmp_path):
    """Issue #15: cleanup must never touch the rest of the session."""
    p = write_map(tmp_path, "m.mrc")
    port = FakePort({"get_names": ["unrelated_structure"]})

    with pytest.raises(PortError, match="Nothing else in the session"):
        load_map(port, p, "m")

    assert port.calls("delete") == [(("m",), {})], port.call_log
    assert not any(a[0] == "unrelated_structure" for a, _ in port.calls("delete"))
    assert not port.queried("reinitialize")


def test_default_object_name_strips_compound_suffixes(tmp_path):
    p = write_map(tmp_path, "emd_1234.map", gz=True)
    port = loaded_port("emd_1234")
    out = load_map(port, p)
    assert "-> emd_1234" in out


def test_unsharpened_is_not_read_as_sharpened(tmp_path):
    """'sharp' is a substring of 'unsharp', and SHARPENED is declared first,
    so the first-match scan inverted the one thing the depositor had troubled
    to record in the filename."""
    from test_mapinfo import write_map

    path = write_map(tmp_path, "emd_1234_unsharpened.mrc")
    evidence = gather_evidence(read_map_header(path), path)

    assert evidence.suggested is Provenance.MEASURED, evidence.reasons
    assert any("unsharp" in r for r in evidence.reasons), evidence.reasons


def test_a_genuinely_sharpened_name_still_reads_as_sharpened(tmp_path):
    """The longest-token rule must not break the case it was guarding."""
    from test_mapinfo import write_map

    path = write_map(tmp_path, "emd_1234_sharpened.mrc")
    evidence = gather_evidence(read_map_header(path), path)

    assert evidence.suggested is Provenance.SHARPENED, evidence.reasons
