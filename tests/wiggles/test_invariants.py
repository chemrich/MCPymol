"""The SPEC invariants, asserted across the whole tier rather than per tool.

Each tier-3 module tests its own behaviour. This file tests the three claims
that have to hold no matter which tool produced the text, because an invariant
enforced in one function and forgotten in the next is not an invariant.

**I1** — every volume readout states provenance.
**I2** — no latent rendering goes out without its method and that method's caveat.
**I3** — no output turns a gap into an absence.

The parametrisation is deliberately over *reports*, not over functions: what
these invariants constrain is what the package says to a user.
"""

from __future__ import annotations

import pytest
from test_mapinfo import write_map

from mcpymol.wiggles.composition import composition_view
from mcpymol.wiggles.deformation import deformation_view
from mcpymol.wiggles.heterogeneity import forget_ensemble, load_ensemble
from mcpymol.wiggles.latent import contains_absence_claim, latent_traverse_view
from mcpymol.wiggles.port import FakePort

ROWS = [
    # One partial occupancy so occupancy_view takes its real path rather than
    # the "everything is fully occupied, nothing to show" shortcut.
    ("A", "1", "ALA", "CA", "", 0.5, 20.0),
    ("A", "2", "GLY", "CA", "", 1.0, 30.0),
    ("A", "3", "SER", "CA", "", 1.0, 40.0),
]


@pytest.fixture(autouse=True)
def _clean():
    forget_ensemble()
    yield
    forget_ensemble()


def _ensemble(tmp_path, name, *, marker="z.pkl", rms=(0.5, 0.8, 1.2)):
    root = tmp_path / name
    root.mkdir()
    if marker:
        (root / marker).write_text("x")
    for i, r in enumerate(rms, start=1):
        write_map(root, f"vol_{i}.mrc", rms=r, dmean=0.0)
    names = [f"{name}_f{i:02d}" for i in range(1, len(rms) + 1)]
    load_ensemble(FakePort({"get_names": names}), root, name)
    return FakePort({"get_names": names}), name


def _model_port():
    start = [(0.0, float(i), 0.0) for i in range(len(ROWS))]
    end = [(float(i) * 2.0, float(i), 0.0) for i in range(len(ROWS))]
    return FakePort(
        {
            "iterate_to_list": list(ROWS),
            "count_states": 2,
            "get_coords": lambda obj, state=1: start if state == 1 else end,
        }
    )


def _all_tier3_reports(tmp_path) -> dict[str, str]:
    """One report from every tier-3 tool, including its refusal paths."""
    port, name = _ensemble(tmp_path, "ens")
    unknown_port, unknown_name = _ensemble(tmp_path, "bare", marker=None)
    comp_port = FakePort({"count_atoms": lambda sel: 10})

    return {
        "load_ensemble": load_ensemble(
            FakePort({"get_names": [f"ens2_f{i:02d}" for i in (1, 2, 3)]}),
            (tmp_path / "ens"),
            "ens2",
        ),
        "latent_traverse_view": latent_traverse_view(port, name),
        "latent_traverse_view/refused": latent_traverse_view(unknown_port, unknown_name),
        "deformation_view": deformation_view(_model_port(), "obj"),
        "deformation_view/no_arrows": deformation_view(_model_port(), "obj", arrows=False),
        "composition_view": composition_view(comp_port, "obj", {"chain A": 0.2, "chain B": 1.0}),
    }


# ── I3: a gap is not an absence ──────────────────────────────────────────────


def test_no_tier3_report_contains_an_absence_claim(tmp_path):
    """The invariant is about what the package says, so it is asserted over
    every report the tier can produce — refusals included, since a refusal is
    still text a user reads."""
    offenders = {
        label: contains_absence_claim(report)
        for label, report in _all_tier3_reports(tmp_path).items()
        if contains_absence_claim(report)
    }
    assert not offenders, f"I3 breached: {offenders}"


def test_the_detector_would_catch_a_breach_if_one_were_introduced(tmp_path):
    """Guards the guard. A detector that never fires proves nothing about the
    reports it is pointed at."""
    reports = _all_tier3_reports(tmp_path)
    doctored = reports["latent_traverse_view"] + "\n  This region is unvisited."
    assert contains_absence_claim(doctored) == "unvisited"


# ── I2: no unlabelled latent rendering ───────────────────────────────────────


def test_every_rendered_latent_view_names_its_method(tmp_path):
    port, name = _ensemble(tmp_path, "ens")
    out = latent_traverse_view(port, name)
    assert "cryoDRGN" in out


def test_a_latent_view_without_a_method_renders_nothing(tmp_path):
    port, name = _ensemble(tmp_path, "bare", marker=None)
    out = latent_traverse_view(port, name)

    assert "REFUSED" in out
    assert not port.calls("isosurface"), port.call_log
    assert not port.calls("mset"), port.call_log


# ── I1: provenance is never dropped ──────────────────────────────────────────


def test_every_volume_report_states_provenance(tmp_path):
    volume_reports = {
        label: report
        for label, report in _all_tier3_reports(tmp_path).items()
        if label.startswith(("load_ensemble", "latent_traverse_view")) and "REFUSED" not in report
    }
    assert volume_reports, "the fixture should produce at least one volume report"

    for label, report in volume_reports.items():
        assert "Provenance:" in report, f"{label} dropped the provenance banner"


# ── the sense split, across both tools that touch occupancy ──────────────────


def test_both_occupancy_tools_name_their_sense(tmp_path):
    """The compendium's single most important design decision: a legend that
    quietly omits which sense it means makes the two tools interchangeable,
    and they are not."""
    from mcpymol.wiggles.occupancy import occupancy_view

    sense1 = occupancy_view(FakePort({"iterate_to_list": list(ROWS)}), "obj")
    sense2 = composition_view(FakePort({"count_atoms": lambda sel: 10}), "obj", {"chain A": 0.4})

    # Each report must *declare* its own sense, not merely mention the word.
    assert "Occupancy shown is SENSE 1" in sense1
    assert "THIS IS OCCUPANCY IN SENSE 2" in sense2

    # And each must disclaim the other's quantity by name, which is the part
    # that makes the declaration useful — "sense 1" means nothing to a reader
    # who does not already know there are two.
    assert "NOT compositional occupancy" in sense1
    assert "NOT the per-atom crystallographic occupancy q" in sense2
