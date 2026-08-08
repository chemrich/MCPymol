"""Tests for local_resolution_view.

Two things carry the judgement in this tool and both are asserted here: the
refusal when the volumes do not share a grid, and the Ångström-to-sigma
conversion of the ramp breakpoints against the *resolution* map's header.

The resolution-map fixture stamps dmin/dmax/dmean/rms = 2/6/4/1, so a
breakpoint in Å converts to sigma by subtracting 4 — 2 Å is -2 sigma, 6 Å is
+2 sigma. The default five-colour palette then lands on exactly
[-2, -1, 0, 1, 2], which makes a wrong conversion impossible to miss.
"""

from __future__ import annotations

import pytest
from test_mapinfo import write_map

from mcpymol.wiggles.density import DEFAULT_SIGMA
from mcpymol.wiggles.localres import (
    DEFAULT_PALETTE,
    grid_differences,
    local_resolution_view,
)
from mcpymol.wiggles.mapinfo import read_map_header
from mcpymol.wiggles.maps import forget_map, load_map
from mcpymol.wiggles.port import FakePort, PortError
from mcpymol.wiggles.provenance import Provenance

#: A resolution field runs 2–6 Å. These are resolutions, not densities, which
#: is the whole reason the conversion needs this map's own statistics.
RES_STATS = {"dmin": 2.0, "dmax": 6.0, "dmean": 4.0, "rms": 1.0}


@pytest.fixture(autouse=True)
def _clean_maps():
    forget_map()
    yield
    forget_map()


@pytest.fixture
def session(tmp_path):
    """Two matched volumes in the session: a density map and a resolution map.

    ``get`` is what ``normalize_ccp4_maps`` comes back as. The default "OK" is
    what FakePort returns for an unstubbed command, i.e. PyMOL declining to
    say — which is a path the tool has to handle, so it is the default here.
    """

    def _make(*, main_kw=None, res_kw=None, get="OK"):
        main_path = write_map(tmp_path, "main.mrc", **(main_kw or {}))
        res_path = write_map(tmp_path, "locres.mrc", **{**RES_STATS, **(res_kw or {})})

        loader = FakePort({"get_names": ["main", "locres"]})
        load_map(loader, main_path, "main", provenance=Provenance.MEASURED)
        load_map(loader, res_path, "locres")  # UNKNOWN, the realistic default

        port = FakePort({"get_names": ["main", "locres"], "get": get})
        return port, read_map_header(main_path), read_map_header(res_path)

    return _make


# ── the grid check ──────────────────────────────────────────────────────────


def test_identical_grids_have_no_differences(tmp_path):
    a = read_map_header(write_map(tmp_path, "a.mrc"))
    b = read_map_header(write_map(tmp_path, "b.mrc", **RES_STATS))
    assert grid_differences(a, b) == []


def test_voxel_size_is_compared_with_a_tolerance(tmp_path):
    """EMD-30913 reports 0.7999967 where another map reports 0.8. Exact
    equality would refuse a pair that matches."""
    a = read_map_header(write_map(tmp_path, "a.mrc", voxel=0.8))
    b = read_map_header(write_map(tmp_path, "b.mrc", voxel=0.7999967))
    assert grid_differences(a, b) == []


@pytest.mark.parametrize(
    ("kw", "expected"),
    [
        ({"nx": 128, "mx": 256, "cella": (256.0, 256.0, 256.0)}, "extent differs"),
        ({"mx": 512, "cella": (512.0, 256.0, 256.0)}, "grid sampling differs"),
        ({"voxel": 1.2}, "voxel size along X differs"),
        ({"origin": (0.0, 0.0, 12.0)}, "origin along Z differs"),
        ({"nzstart": 5}, "start position differs"),
        ({"mapc": 3, "mapr": 1, "maps": 2}, "axis order differs"),
    ],
)
def test_each_way_the_grids_can_differ_is_named(tmp_path, kw, expected):
    a = read_map_header(write_map(tmp_path, "a.mrc"))
    b = read_map_header(write_map(tmp_path, "b.mrc", **kw))
    differences = grid_differences(a, b)
    assert any(expected in d for d in differences), differences


def test_undefined_voxel_size_is_a_difference_not_a_crash(tmp_path):
    a = read_map_header(write_map(tmp_path, "a.mrc"))
    b = read_map_header(write_map(tmp_path, "b.mrc", cella=(0.0, 256.0, 256.0)))
    assert any("undefined" in d for d in grid_differences(a, b))


def test_mismatched_grids_are_refused_and_nothing_is_drawn(session):
    port, _, _ = session(res_kw={"origin": (0.0, 0.0, 12.0)})
    out = local_resolution_view(port, "main", "locres")

    assert "REFUSED" in out
    assert "origin along Z differs" in out
    assert not port.calls("isosurface"), port.call_log
    assert not port.calls("ramp_new"), port.call_log


def test_the_refusal_explains_why_a_wrong_grid_is_dangerous(session):
    """A mismatch does not render visibly broken — that is the point."""
    port, _, _ = session(res_kw={"voxel": 1.2})
    out = local_resolution_view(port, "main", "locres")
    assert "smooth, plausible, and wrong" in out


def test_the_same_object_twice_is_refused(session):
    port, _, _ = session()
    out = local_resolution_view(port, "main", "main")
    assert "REFUSED" in out
    assert "Colouring a map by itself" in out
    assert not port.calls("isosurface"), port.call_log


# ── the conversion, which is the other half ─────────────────────────────────


def test_default_breakpoints_span_the_headers_range(session):
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")
    assert "2–6 Å" in out or "2-6 Å" in out


def test_breakpoints_reach_pymol_in_the_resolution_maps_sigma(session):
    """The central claim: PyMOL normalised the volume, so ramp_new needs sigma
    against *this* map's mean and rms, not Ångström and not the main map's."""
    port, _, _ = session()
    local_resolution_view(port, "main", "locres")

    args, _ = port.calls("ramp_new")[0]
    assert args[0] == "locres_ramp"
    assert args[1] == "locres"
    assert args[2] == pytest.approx([-2.0, -1.0, 0.0, 1.0, 2.0]), port.call_log
    assert args[3] == list(DEFAULT_PALETTE)


def test_both_units_appear_in_the_report(session):
    """Stating one unit without the other is how the two get confused."""
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")

    assert "Resolution ramp" in out
    rows = [line for line in out.splitlines() if "->" in line and "Å" in line]
    assert len(rows) == len(DEFAULT_PALETTE), out
    assert "2.00 Å" in rows[0] and "-2" in rows[0], rows[0]
    assert "6.00 Å" in rows[-1] and "2" in rows[-1], rows[-1]


def test_a_pymol_that_will_not_answer_the_setting_query_is_handled(session):
    """Older plugins may not expose `get` at all. An unanswerable question is
    not a failure — it is the unknown case, which is already handled."""

    def _refuse(*args, **kwargs):
        raise PortError("unknown action 'get'")

    port, _, _ = session(get=_refuse)
    out = local_resolution_view(port, "main", "locres")

    assert "would not report normalize_ccp4_maps" in out
    assert port.calls("ramp_new")[0][0][2] == pytest.approx([-2.0, -1.0, 0.0, 1.0, 2.0])


def test_normalisation_off_sends_angstrom_unconverted(session):
    """With normalize_ccp4_maps off the stored values are still resolutions,
    so converting would be the error."""
    port, _, _ = session(get="off")
    local_resolution_view(port, "main", "locres")

    args, _ = port.calls("ramp_new")[0]
    assert args[2] == pytest.approx([2.0, 3.0, 4.0, 5.0, 6.0]), port.call_log


def test_unknown_normalisation_assumes_the_pymol_default_and_says_so(session):
    port, _, _ = session()  # FakePort answers "OK", which is not a setting
    out = local_resolution_view(port, "main", "locres")

    assert "would not report normalize_ccp4_maps" in out
    assert port.calls("ramp_new")[0][0][2] == pytest.approx([-2.0, -1.0, 0.0, 1.0, 2.0])


def test_the_setting_is_read_now_not_at_load_time(session):
    """A limitation worth stating rather than papering over."""
    port, _, _ = session(get="on")
    out = local_resolution_view(port, "main", "locres")
    assert "read now, not as it was at load time" in out


def test_zero_rms_refuses_rather_than_dividing(session):
    port, _, _ = session(res_kw={"rms": 0.0})
    out = local_resolution_view(port, "main", "locres")

    assert "REFUSED" in out
    assert "sigma is undefined" in out
    assert "normalize_ccp4_maps, off" in out
    assert not port.calls("ramp_new"), port.call_log


def test_a_map_with_no_positive_values_is_not_a_resolution_field(session):
    port, _, _ = session(res_kw={"dmin": -1.0, "dmax": 0.0, "dmean": -0.5})
    out = local_resolution_view(port, "main", "locres")

    assert "REFUSED" in out
    assert "does not look like a resolution field" in out


def test_zero_padding_outside_the_mask_does_not_anchor_the_ramp(session):
    """Estimators write 0 outside the mask. Taking that as the best-resolved
    end would compress every real value into the top of the scale."""
    port, _, _ = session(res_kw={"dmin": 0.0})
    out = local_resolution_view(port, "main", "locres")

    assert "outside-the-mask padding" in out
    assert port.calls("ramp_new")[0][0][2][0] > -4.0, port.call_log


# ── breakpoints and palette given explicitly ────────────────────────────────


def test_explicit_breakpoints_are_used(session):
    port, _, _ = session()
    local_resolution_view(port, "main", "locres", breaks=[3.0, 5.0], palette=["blue", "red"])
    assert port.calls("ramp_new")[0][0][2] == pytest.approx([-1.0, 1.0])


def test_descending_breakpoints_are_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="must ascend"):
        local_resolution_view(port, "main", "locres", breaks=[5.0, 3.0], palette=["blue", "red"])


def test_non_positive_breakpoints_are_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="must be positive"):
        local_resolution_view(port, "main", "locres", breaks=[0.0, 3.0], palette=["blue", "red"])


def test_a_single_breakpoint_is_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="at least two breakpoints"):
        local_resolution_view(port, "main", "locres", breaks=[3.0], palette=["blue"])


def test_palette_must_match_the_breakpoints(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="one to one"):
        local_resolution_view(port, "main", "locres", breaks=[2.0, 4.0, 6.0])


def test_explicit_breakpoints_raise_rather_than_refusing(session):
    """A bad default is the map's fault and gets a report; a bad argument is
    the caller's and gets an exception."""
    port, _, _ = session(res_kw={"dmax": 0.0, "dmin": -1.0, "dmean": -0.5})
    with pytest.raises(ValueError, match="must be positive"):
        local_resolution_view(port, "main", "locres", breaks=[-1.0, 2.0], palette=["blue", "red"])


# ── the surface itself ──────────────────────────────────────────────────────


def test_the_surface_is_contoured_in_the_main_maps_sigma(session):
    port, _, _ = session()
    local_resolution_view(port, "main", "locres", level=2.0)

    args, _ = port.calls("isosurface")[0]
    assert args == ("main_localres", "main", 2.0), port.call_log


def test_an_absolute_contour_is_converted_against_the_main_map(session):
    """main.mrc has dmean=0, rms=0.5, so 1.0 absolute is 2 sigma — and the
    resolution map's statistics must not be the ones used."""
    port, _, _ = session()
    local_resolution_view(port, "main", "locres", level=1.0, units="absolute")
    assert port.calls("isosurface")[0][0][2] == pytest.approx(2.0), port.call_log


def test_default_contour_is_labelled_generic(session):
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")

    assert port.calls("isosurface")[0][0][2] == DEFAULT_SIGMA
    assert "not a recommendation for this map" in out


def test_bad_units_are_rejected(session):
    port, _, _ = session()
    with pytest.raises(ValueError, match="must be 'sigma' or 'absolute'"):
        local_resolution_view(port, "main", "locres", level=1.0, units="angstrom")


def test_the_ramp_is_attached_to_the_surface(session):
    port, _, _ = session()
    local_resolution_view(port, "main", "locres")
    assert port.called("set", "surface_color", "locres_ramp", "main_localres"), port.call_log


def test_a_selection_carves_the_surface(session):
    port, _, _ = session()
    local_resolution_view(port, "main", "locres", selection="chain A", carve=3.5)

    args, kwargs = port.calls("isosurface")[0]
    assert args[3] == "chain A"
    assert kwargs["carve"] == 3.5


def test_names_can_be_overridden(session):
    port, _, _ = session()
    local_resolution_view(port, "main", "locres", name="surf", ramp_name="rr")

    assert port.calls("isosurface")[0][0][0] == "surf"
    assert port.calls("ramp_new")[0][0][0] == "rr"
    assert port.called("set", "surface_color", "rr", "surf")


def test_validate_only_creates_nothing(session):
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres", validate_only=True)

    assert "Grid check passed" in out
    assert "nothing created" in out
    assert not port.calls("isosurface"), port.call_log
    assert not port.calls("ramp_new"), port.call_log


# ── what the report has to say ──────────────────────────────────────────────


def test_both_maps_get_a_provenance_banner(session):
    """I1: this renders a volume and colours it by a second one, so both
    origins have to be stated."""
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")

    assert "Provenance: MEASURED" in out
    assert "Provenance: UNKNOWN" in out


def test_the_colour_direction_is_stated_in_words(session):
    """Low Å is good, so the ramp runs opposite to the usual reading."""
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")

    assert "BEST-resolved" in out
    assert "Low numbers are good" in out


def test_the_two_sigma_scales_are_distinguished(session):
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")
    assert "not interchangeable" in out


def test_the_legend_says_the_field_is_an_estimate(session):
    port, _, _ = session()
    out = local_resolution_view(port, "main", "locres")
    assert "ESTIMATE, not a measurement" in out


def test_geometry_warnings_from_both_maps_are_labelled(session):
    """Both are anisotropic in the same way, so the grids still match and the
    warnings are what is left to report."""
    port, _, _ = session(
        main_kw={"cella": (256.0, 256.0, 384.0)},
        res_kw={"cella": (256.0, 256.0, 384.0)},
    )
    out = local_resolution_view(port, "main", "locres")

    assert "Geometry warnings, main" in out
    assert "Geometry warnings, locres" in out
    assert "ANISOTROPIC" in out


# ── refusing to guess ───────────────────────────────────────────────────────


def test_an_unloaded_density_map_is_refused(session):
    port, _, _ = session()
    with pytest.raises(PortError, match="the density map"):
        local_resolution_view(port, "never_loaded", "locres")


def test_an_unloaded_resolution_map_is_refused(session):
    port, _, _ = session()
    with pytest.raises(PortError, match="the resolution map"):
        local_resolution_view(port, "main", "never_loaded")


def test_the_refusal_explains_what_the_header_is_needed_for(session):
    port, _, _ = session()
    with pytest.raises(PortError, match="share a grid"):
        local_resolution_view(port, "main", "never_loaded")
