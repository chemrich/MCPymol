"""Live regression against real EMDB maps.

Skipped unless ``WIGGLES_LIVE=1``, because it downloads and needs network. The
synthetic tests in ``test_mapinfo.py`` prove the parser is self-consistent;
only this file proves it agrees with reality.

    WIGGLES_LIVE=1 pytest tests/test_mapinfo_live.py -v

Each case asserts against the value **EMDB's own API reports**, fetched in the
same run, so the test cannot drift from the archive.
"""

from __future__ import annotations

import gzip
import json
import os
import ssl
import urllib.error
import urllib.request
import zlib

import pytest

from mcpymol.wiggles.mapinfo import read_map_header

pytestmark = pytest.mark.skipif(
    os.environ.get("WIGGLES_LIVE") != "1",
    reason="set WIGGLES_LIVE=1 to run live EMDB downloads",
)

API = "https://www.ebi.ac.uk/emdb/api/entry/{acc}"
FTP = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{acc}/map/{lower}.map.gz"

# (accession, why this one)
CASES = [
    ("EMD-30913", "reports 0.7999967 — the float-noise isotropy case"),
    ("EMD-11638", "1.22 A apoferritin — the compendium's worked example"),
]


def _fetch(url: str, timeout: int, limit: int | None = None) -> bytes:
    """GET, converting transport failures into skips.

    A TLS or DNS failure is an environment problem, not a parser bug, and this
    file exists to test the parser. Skipping keeps the distinction visible:
    a real disagreement with EMDB fails, an unreachable EMDB skips.

    Behind a TLS-inspecting proxy, ``uv run`` is the answer: uv-managed Python
    resolves its CA store to the macOS system bundle
    (``/private/etc/ssl/cert.pem``) and trusts the proxy root. A python.org
    interpreter ships its own bundle, reports ``cafile = None``, and skips
    every test here with CERTIFICATE_VERIFY_FAILED. If you are stuck on one,
    export the keychain and set ``SSL_CERT_FILE`` to it.
    """
    try:
        with urllib.request.urlopen(url, timeout=timeout) as fh:
            return fh.read(limit) if limit else fh.read()
    except (urllib.error.URLError, ssl.SSLError, TimeoutError, OSError) as exc:
        pytest.skip(f"cannot reach {url}: {type(exc).__name__}: {exc}")


def _emdb_metadata(acc: str) -> dict:
    return json.loads(_fetch(API.format(acc=acc), timeout=60))


def _download_header(acc: str, tmp_path):
    """Fetch just enough of the gzipped map to recover a 1024-byte header.

    Downloads a bounded prefix rather than the whole map — EMD-11638 is 62 MB
    and we need the first kilobyte. ``zlib.decompressobj`` is used directly
    because ``gzip.decompress`` rejects a truncated stream, which is exactly
    what a bounded prefix is.
    """
    url = FTP.format(acc=acc, lower=acc.replace("-", "_").lower())
    compressed = _fetch(url, timeout=300, limit=1024 * 1024)  # 1 MiB covers a 1 KiB header

    obj = zlib.decompressobj(16 + zlib.MAX_WBITS)  # 16 => gzip wrapper
    raw = obj.decompress(compressed, 8192)
    assert len(raw) >= 1024, f"{acc}: recovered only {len(raw)} bytes, need 1024"

    dest = tmp_path / f"{acc}.map.gz"
    dest.write_bytes(gzip.compress(raw[:4096]))
    return dest


@pytest.mark.parametrize("acc,why", CASES, ids=[c[0] for c in CASES])
def test_parser_agrees_with_emdb(acc, why, tmp_path):
    meta = _emdb_metadata(acc)["map"]
    expected = tuple(float(meta["pixel_spacing"][a]["valueOf_"]) for a in "xyz")
    expected_dims = (
        meta["dimensions"]["col"],
        meta["dimensions"]["row"],
        meta["dimensions"]["sec"],
    )

    h = read_map_header(_download_header(acc, tmp_path))

    assert (h.nx, h.ny, h.nz) == expected_dims, f"{acc}: dimensions disagree"
    for axis, got, want in zip("XYZ", h.voxel_size, expected, strict=True):
        assert got == pytest.approx(want, rel=1e-6), (
            f"{acc}: {axis} spacing {got} != EMDB's {want} ({why})"
        )


def test_emd_30913_float_noise_is_isotropic(tmp_path):
    """The real file, not a synthetic one: 0.7999967 must not read anisotropic."""
    h = read_map_header(_download_header("EMD-30913", tmp_path))
    assert h.is_isotropic is True
    assert not any("ANISOTROPIC" in w for w in h.warnings())
