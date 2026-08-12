"""Tests for density_view — mostly about the sigma/absolute unit trap."""

from __future__ import annotations

import pytest
from test_mapinfo import write_map

from mcpymol.wiggles.density import DEFAULT_SIGMA, density_view, to_absolute, to_sigma
from mcpymol.wiggles.mapinfo import read_map_header
from mcpymol.wiggles.maps import forget_map, load_map
from mcpymol.wiggles.port import FakePort, PortError
from mcpymol.wiggles.provenance import Provenance


@pytest.fixture(autouse=True)
def _clean_maps():
    forget_map()
    yield
    forget_map()


@pytest.fixture
def loaded(tmp_path):
    """A map in the session, with known header statistics.

    write_map stamps dmin/dmax/dmean = -1/1/0 and rms = 0.5, so sigma
    conversions here are exact and easy to reason about.
    """

    def _load(filename="m.mrc", **kw):
        path = write_map(tmp_path, filename, **kw)
        obj = filename.split(".")[0]
        port = FakePort({"get_names": [obj]})
        load_map(port, path, obj, provenance=Provenance.MEASURED)
        return FakePort({"get_names": [obj]}), obj, read_map_header(path)

    return _load


# ── the conversion ──────────────────────────────────────────────────────────


def test_sigma_and_absolute_are_inverses(loaded):
    _, _, header = loaded()
    for sigma in (0.0, 1.5, 3.16, -2.0):
        assert to_sigma(header, to_absolute(header, sigma)) == pytest.approx(sigma)


def test_absolute_to_sigma_uses_mean_and_rms(loaded):
    """dmean=0, rms=0.5 in the fixture, so 1.0 absolute is 2 sigma."""
    _, _, header = loaded()
    assert to_sigma(header, 1.0) == pytest.approx(2.0)


@pytest.mark.parametrize("rms", [0.0, -1.0])
def test_unusable_rms_cannot_be_converted(tmp_path, rms):
    """Neither zero nor negative RMS can define a sigma scale.

    -1 is MRC2014's marker for "statistics not computed" — mrcfile writes it,
    with dmean=-2, whenever a map is saved without update_header_stats(), and
    plenty of processing pipelines never rewrite them. It is the more
    dangerous of the two because it divides cleanly: `if not header.rms`
    passed it straight through, silently inverting the sign of every
    conversion.
    """
    path = write_map(tmp_path, "z.mrc")
    header = read_map_header(path)
    header = type(header)(**{**header.__dict__, "rms": rms})
    with pytest.raises(ValueError, match="cannot define a sigma scale"):
        to_sigma(header, 1.0)


def test_real_emdb_contour_converts_to_a_sensible_sigma(loaded):
    """Regression for the trap: EMD-30913 publishes 0.05 absolute, which is
    3.16 sigma against its real header — not 0.05 sigma."""
    _, _, header = loaded()
    header = type(header)(**{**header.__dict__, "dmean": 0.000921512, "rms": 0.0155224})
    assert to_sigma(header, 0.05) == pytest.approx(3.16, abs=0.01)


# ── the view ────────────────────────────────────────────────────────────────


def test_reports_the_level_in_both_units(loaded):
    port, obj, _ = loaded()
    out = density_view(port, obj, "chain A", level=2.0)

    assert "2 sigma" in out
    assert "absolute" in out
    assert port.calls("isomesh")[0][0][:3] == (f"{obj}_mesh", obj, 2.0), port.call_log


def test_absolute_level_is_converted_before_reaching_pymol(loaded):
    """The whole point: PyMOL contours in sigma."""
    port, obj, _ = loaded()
    out = density_view(port, obj, "chain A", level=1.0, units="absolute")

    assert port.calls("isomesh")[0][0][2] == pytest.approx(2.0), port.call_log
    assert "would contour near zero and show mostly noise" in out


def test_default_level_is_labelled_as_generic(loaded):
    port, obj, _ = loaded()
    out = density_view(port, obj, "chain A")

    assert f"{DEFAULT_SIGMA} sigma was used" in out
    assert "not a recommendation for this map" in out


def test_emdb_map_points_at_the_author_contour(tmp_path):
    """The finding this tool exists for — and it must say the level is
    absolute, or the advice would cause the very bug it warns about."""
    path = write_map(tmp_path, "emd_30913.mrc")
    port = FakePort({"get_names": ["emd_30913"]})
    load_map(port, path, "emd_30913")

    port = FakePort({"get_names": ["emd_30913"]})
    out = density_view(port, "emd_30913", "chain A")

    assert "EMD-30913 has an author-recommended contour" in out
    assert "ABSOLUTE value" in out
    assert "units='absolute'" in out
    assert "ebi.ac.uk/emdb/api/entry/EMD-30913" in out


def test_non_emdb_map_does_not_invent_an_author_contour(loaded):
    port, obj, _ = loaded("plain.mrc")
    out = density_view(port, obj, "chain A")
    assert "author-recommended" not in out


def test_carries_the_provenance_banner(loaded):
    """I1: this renders a volume, so the readout must say where it came from."""
    port, obj, _ = loaded()
    out = density_view(port, obj, "chain A")
    assert "Provenance: MEASURED" in out


def test_unloaded_map_is_refused_with_the_reason(loaded):
    port = FakePort()
    with pytest.raises(PortError, match="not loaded through load_map"):
        density_view(port, "never_loaded", "chain A")


def test_refusal_explains_the_unit_hazard(loaded):
    """The error should teach the trap, not just decline."""
    with pytest.raises(PortError, match="author contour gets used as a sigma level"):
        density_view(FakePort(), "never_loaded", "chain A")


def test_bad_units_are_rejected(loaded):
    port, obj, _ = loaded()
    with pytest.raises(ValueError, match="must be 'sigma' or 'absolute'"):
        density_view(port, obj, "chain A", level=1.0, units="angstrom")


def test_carve_radius_is_passed_through(loaded):
    port, obj, _ = loaded()
    density_view(port, obj, "chain A", level=1.0, carve=3.5)
    assert port.calls("isomesh")[0][1]["carve"] == 3.5


def test_mesh_name_can_be_overridden(loaded):
    port, obj, _ = loaded()
    density_view(port, obj, "chain A", level=1.0, name="pocket")
    assert port.calls("isomesh")[0][0][0] == "pocket"


def test_geometry_warnings_reach_the_report(loaded):
    port, obj, _ = loaded("aniso.mrc", nx=100, ny=100, nz=100, cella=(100.0, 100.0, 150.0))
    out = density_view(port, obj, "chain A", level=1.0)
    assert "ANISOTROPIC" in out


# ── which number actually reaches isomesh ───────────────────────────────────


def _isomesh_level(port):
    args, _ = port.calls("isomesh")[0]
    return args[2]


def test_normalised_map_is_contoured_in_sigma(loaded):
    port, obj, _header = loaded()
    port.responses["get"] = "1"  # normalize_ccp4_maps on

    report = density_view(port, obj, "chain A", level=1.0, units="absolute")

    # dmean=0, rms=0.5, so 1.0 absolute is 2 sigma.
    assert _isomesh_level(port) == pytest.approx(2.0)
    assert "is on" in report


def test_unnormalised_map_is_contoured_in_absolute_units(loaded):
    """With normalize_ccp4_maps off, PyMOL reads the level as a raw map value.

    Sending sigma there is the failure this module exists to prevent:
    EMD-30913's published 0.05 is 3.16 sigma, and 3.16 as a raw density
    contours nothing at all — an empty mesh under a report claiming the
    depositor's own level was applied.
    """
    port, obj, _header = loaded()
    port.responses["get"] = "0"  # normalize_ccp4_maps off

    report = density_view(port, obj, "chain A", level=1.0, units="absolute")

    assert _isomesh_level(port) == pytest.approx(1.0), (
        "sent sigma to an un-normalised map, so the mesh is at the wrong density"
    )
    assert "OFF" in report
    # The sigma/absolute pair still describes the map honestly.
    assert "2 sigma" in report


def test_the_report_says_when_normalisation_could_not_be_determined(loaded):
    """An older plugin may not expose `get`. Assuming on is the right default,
    but the report has to say the contour rests on that assumption."""
    port, obj, _header = loaded()
    port.responses["get"] = "something unexpected"

    report = density_view(port, obj, "chain A", level=1.0, units="absolute")

    assert _isomesh_level(port) == pytest.approx(2.0)
    assert "would not report" in report


def test_a_deleted_map_is_not_contoured_from_its_stale_header(loaded):
    """The registry is keyed by object name and nothing evicts from it.

    Without the existence check, deleting a map leaves its header behind, and
    the next density_view converts levels with the statistics of a volume that
    is no longer loaded — plus a provenance banner asserting a measurement for
    whatever now holds the name.
    """
    port, obj, _header = loaded()
    port.responses["get_names"] = []  # the object is gone from the session

    with pytest.raises(PortError, match="was not loaded through load_map"):
        density_view(port, obj, "chain A", level=1.0, units="absolute")


def test_the_report_names_the_file_the_header_came_from(loaded):
    """A map deleted and replaced under the same name still passes the
    existence check, because PyMOL does not expose an object's source file.
    Printing the path is what makes that substitution visible to a reader."""
    port, obj, _header = loaded()
    port.responses["get"] = "1"

    report = density_view(port, obj, "chain A", level=1.0, units="absolute")

    assert "Header read from:" in report
    assert "m.mrc" in report
