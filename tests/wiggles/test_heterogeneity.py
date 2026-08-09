"""Tests for load_ensemble — ordering, method detection, and what is not guessed.

Two things here are load-bearing rather than convenient. Frame order *is* the
trajectory, so natural sorting is asserted directly. And method detection is
what invariant I2 hangs off: a directory that matches no documented marker must
come back UNKNOWN rather than plausible, because a wrong method label attaches
the wrong caveat to a rendering and the rendering looks the same either way.
"""

from __future__ import annotations

import pytest
from test_mapinfo import write_map

from mcpymol.wiggles.heterogeneity import (
    Method,
    detect_method,
    find_volumes,
    forget_ensemble,
    load_ensemble,
    loaded_ensemble,
    read_latent_table,
)
from mcpymol.wiggles.port import FakePort, PortError
from mcpymol.wiggles.provenance import Provenance, provenance_of


@pytest.fixture(autouse=True)
def _clean():
    forget_ensemble()
    yield
    forget_ensemble()


@pytest.fixture
def job(tmp_path):
    """A heterogeneity job directory with ``n`` frames and optional markers."""

    def _make(n=4, *, marker=None, subdir=None, rms=0.5, extra=()):
        root = tmp_path / "job042"
        target = root / subdir if subdir else root
        target.mkdir(parents=True, exist_ok=True)
        for i in range(1, n + 1):
            write_map(target, f"vol_{i}.mrc", rms=rms, dmean=0.0)
        if marker:
            (root / marker).write_text("x")
        for name in extra:
            write_map(root, name)
        return root

    return _make


def _port(root, n, subdir=None):
    names = [f"job042_f{i:02d}" for i in range(1, n + 1)]
    return FakePort({"get_names": names})


# ── ordering: the trajectory is the order ────────────────────────────────────


def test_frames_sort_naturally_not_lexicographically(tmp_path):
    """vol_10 after vol_9. Lexicographic order would put it between 1 and 2 and
    silently reorder the motion — the animation still plays, and it is wrong."""
    d = tmp_path / "j"
    d.mkdir()
    for i in (1, 2, 9, 10, 11):
        write_map(d, f"vol_{i}.mrc")

    order = [p.name for p in find_volumes(d)]
    assert order == ["vol_1.mrc", "vol_2.mrc", "vol_9.mrc", "vol_10.mrc", "vol_11.mrc"], order


def test_half_maps_and_masks_are_not_frames(tmp_path):
    d = tmp_path / "j"
    d.mkdir()
    write_map(d, "vol_1.mrc")
    write_map(d, "run_half_map_1.mrc")
    write_map(d, "emd_30913_msk_1.mrc")

    assert [p.name for p in find_volumes(d)] == ["vol_1.mrc"]


def test_volumes_are_found_in_subdirectories(job):
    """cryoDRGN puts sampled volumes under analyze.N/kmeans20/, not at the top."""
    root = job(3, subdir="analyze.9/kmeans20")
    assert len(find_volumes(root)) == 3


# ── method detection: evidence, never a guess ────────────────────────────────


@pytest.mark.parametrize(
    ("marker", "expected"),
    [
        ("z.pkl", Method.CRYODRGN),
        ("config.yaml", Method.CRYODRGN),
        ("deformations.star", Method.DYNAMIGHT),
        ("eigenvolumes", Method.RECOVAR),
        ("part_segmentation", Method.CRYOSPIRE),
        ("particles.cs", Method.THREE_DVA),
    ],
)
def test_documented_markers_identify_their_method(job, marker, expected):
    root = job(2, marker=marker)
    method, evidence = detect_method(root)

    assert method is expected
    assert evidence and marker.split(".")[0] in " ".join(evidence).lower()


def test_an_unrecognised_directory_stays_unknown(job):
    """The honest answer, and the one that makes I2 bite."""
    root = job(2)
    method, evidence = detect_method(root)

    assert method is Method.UNKNOWN
    assert "no documented marker found" in evidence[0]
    assert "z.pkl" in evidence[0], "the report should say what it looked for"


def test_the_unknown_evidence_points_at_the_way_out(job):
    _, evidence = detect_method(job(2))
    assert any("method=" in line for line in evidence)


# ── loading ──────────────────────────────────────────────────────────────────


def test_frames_load_in_order_and_register(job):
    root = job(3, marker="z.pkl")
    port = _port(root, 3)
    out = load_ensemble(port, root, "job042")

    ensemble = loaded_ensemble("job042")
    assert ensemble is not None
    assert ensemble.n_frames == 3
    assert ensemble.method is Method.CRYODRGN
    assert ensemble.objects == ("job042_f01", "job042_f02", "job042_f03")
    assert "cryoDRGN" in out
    loaded = [args[1] for args, _ in port.calls("load")]
    assert loaded == list(ensemble.objects), port.call_log


def test_provenance_is_declared_generated_and_says_so(job):
    """Not the inference I1 forbids — I1 exists to stop a generated volume being
    called measured, and this asserts the conservative direction."""
    root = job(2, marker="z.pkl")
    out = load_ensemble(_port(root, 2), root, "job042")

    assert provenance_of("job042_f01") is Provenance.GENERATED
    assert "GENERATED" in out
    assert "Declared rather than inferred" in out


def test_provenance_can_be_overridden(job):
    root = job(2, marker="z.pkl")
    load_ensemble(_port(root, 2), root, "job042", provenance="measured")
    assert provenance_of("job042_f01") is Provenance.MEASURED


def test_a_declared_method_is_marked_as_declared(job):
    root = job(2)
    out = load_ensemble(_port(root, 2), root, "job042", method="3dva")

    assert loaded_ensemble("job042").method is Method.THREE_DVA
    assert "[declared]" in out
    assert "not detected" in out


def test_an_unknown_method_warns_that_latent_views_will_refuse(job):
    root = job(2)
    out = load_ensemble(_port(root, 2), root, "job042")

    assert "could not be identified" in out
    assert "REFUSE" in out
    assert "I2" in out


def test_truncation_is_reported_never_silent(job):
    root = job(6, marker="z.pkl")
    out = load_ensemble(_port(root, 6), root, "job042", max_volumes=3)

    assert loaded_ensemble("job042").n_frames == 3
    assert "3 further volume(s) were NOT loaded" in out
    assert "not the whole job" in out


def test_a_directory_with_no_volumes_is_refused(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(PortError, match="no MRC/CCP4 volumes"):
        load_ensemble(FakePort(), empty, "x")


def test_a_missing_directory_is_refused(tmp_path):
    with pytest.raises(PortError, match="not a directory"):
        load_ensemble(FakePort(), tmp_path / "nope", "x")


def test_a_frame_that_does_not_arrive_cleans_up_after_itself(job):
    """Issue #15 discipline: scope cleanup to what this call created."""
    root = job(2, marker="z.pkl")
    port = FakePort({"get_names": ["job042_f01"]})  # second never arrives

    with pytest.raises(PortError, match="did not arrive"):
        load_ensemble(port, root, "job042")

    deleted = [args[0] for args, _ in port.calls("delete")]
    assert deleted == ["job042_f01", "job042_f02"], port.call_log
    assert loaded_ensemble("job042") is None


# ── latent tables ────────────────────────────────────────────────────────────


def test_a_text_latent_table_is_read(job):
    root = job(2, marker="z.pkl")
    (root / "z_values.txt").write_text("0.1 0.2\n0.3 0.4\n0.5 0.6\n")

    table = read_latent_table(root)
    assert table.available
    assert table.dimensions == 2
    assert len(table.rows) == 3


def test_star_headers_are_skipped(job):
    root = job(2)
    (root / "latent.star").write_text(
        "data_particles\nloop_\n_rlnZ1 #1\n_rlnZ2 #2\n1.0 2.0\n3.0 4.0\n"
    )

    table = read_latent_table(root)
    assert table.rows == ((1.0, 2.0), (3.0, 4.0))


def test_a_pickle_is_found_but_not_read_by_default(job):
    """Unpickling runs arbitrary code. Opt-in, and the reason is stated."""
    root = job(2, marker="z.pkl")
    table = read_latent_table(root)

    assert not table.available
    assert "not read" in table.unread_reason
    assert "arbitrary code" in table.unread_reason
    assert "trust_pickle=True" in table.unread_reason


def test_a_pickle_is_read_when_trusted(job):
    import pickle

    root = job(2)
    (root / "z.pkl").write_bytes(pickle.dumps([[1.0, 2.0], [3.0, 4.0]]))

    table = read_latent_table(root, trust_pickle=True)
    assert table.rows == ((1.0, 2.0), (3.0, 4.0))


def test_a_plain_npy_is_read_without_numpy_and_without_unpickling(job):
    import struct

    root = job(2)
    rows = [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]
    header = b"{'descr': '<f8', 'fortran_order': False, 'shape': (3, 2), }"
    header += b" " * ((64 - (10 + len(header)) % 64) % 64) + b"\n"
    body = b"".join(struct.pack("<2d", *r) for r in rows)
    (root / "z_latent.npy").write_bytes(
        b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header + body
    )

    table = read_latent_table(root)
    assert table.rows == tuple(rows)


def test_an_object_array_is_refused_rather_than_unpickled(job):
    import struct

    root = job(2)
    header = b"{'descr': '|O', 'fortran_order': False, 'shape': (2, 2), }\n"
    (root / "z_latent.npy").write_bytes(
        b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header
    )

    table = read_latent_table(root)
    assert not table.available
    assert "unpickling" in table.unread_reason


def test_a_missing_latent_table_is_reported_not_fatal(job):
    root = job(3, marker="config.yaml")
    out = load_ensemble(_port(root, 3), root, "job042")

    assert "Latent table: not read" in out
    assert "Nothing here depends on it" in out
    assert loaded_ensemble("job042").n_frames == 3
