"""Where a volume came from — SPEC.md invariant **I1**.

A latent-traversal frame, a neural-decoder volume, a U-Net-enhanced map and a
measured reconstruction all render identically once they are an isosurface.
Provenance is the one thing an isosurface cannot show, so Wiggles carries it
alongside the volume and puts it in every readout.

**The default is UNKNOWN, and that is the whole design.** Defaulting to
"measured" would be the dangerous direction to be wrong in: it would take a
generated volume and quietly assert it was observed. So nothing here classifies
a map automatically. :func:`gather_evidence` reports what the file says about
itself and may *suggest* a category; only the caller can declare one.

See the Wiggles compendium entry `benchmarks` and
the Wiggles compendium entry `local-resolution` for why this matters.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from mcpymol.wiggles.mapinfo import MapHeader


class Provenance(Enum):
    """How a volume came to exist. Ordered from observed to invented."""

    MEASURED = "measured"
    """Reconstructed from images. Sharpening not applied, or not known to be."""

    SHARPENED = "sharpened"
    """B-factor sharpened or locally scaled. Still derived from the data."""

    NN_ENHANCED = "nn_enhanced"
    """Passed through a trained network — DeepEMhancer, EMReady, LocScale-FEM.
    These produce plausible density where the data are weak, by design."""

    GENERATED = "generated"
    """Produced by a model rather than reconstructed: a latent traversal frame,
    a decoder output, a 3DVA component volume. Nothing observed this."""

    UNKNOWN = "unknown"
    """Not declared and not inferable. The honest default."""

    @property
    def is_observed(self) -> bool:
        """Did measurement, rather than a model, put the density there?"""
        return self in (Provenance.MEASURED, Provenance.SHARPENED)

    @property
    def caveat(self) -> str:
        """The one-line warning that belongs in any readout of this volume."""
        return {
            Provenance.MEASURED: "Reconstructed from images.",
            Provenance.SHARPENED: (
                "Sharpened or locally scaled — contrast has been altered, and "
                "features are easier to over-read than in the unsharpened map."
            ),
            Provenance.NN_ENHANCED: (
                "NETWORK-ENHANCED. A model trained on good density produces good-"
                "looking density, including where the data do not support it. "
                "Treat this map as a hypothesis, not a measurement."
            ),
            Provenance.GENERATED: (
                "GENERATED, not observed. This volume was produced by a model. "
                "No particle was reconstructed into exactly this density."
            ),
            Provenance.UNKNOWN: (
                "PROVENANCE UNKNOWN. Nothing here establishes whether this map "
                "was measured, sharpened, network-enhanced or generated. Do not "
                "assume it was measured."
            ),
        }[self]


# Filename and label tokens that are evidence of a category. Deliberately
# suggestive only — see gather_evidence.
_TOKENS: tuple[tuple[Provenance, tuple[str, ...]], ...] = (
    (Provenance.NN_ENHANCED, ("deepemhancer", "emready", "locscale", "cryolvm")),
    (
        Provenance.GENERATED,
        (
            "cryodrgn",
            "3dva",
            "3dflex",
            "dynamight",
            "recovar",
            "cryospire",
            "kmeans",
            "latent",
            "traverse",
        ),
    ),
    (Provenance.SHARPENED, ("sharp", "postprocess", "_pp", "locres", "bfac")),
    (Provenance.MEASURED, ("unfil", "unsharp", "half1", "half2", "raw")),
)

# EMDB writes "EMD-30913" into the MRC label but names the downloaded file
# "emd_30913.map", so both separators have to be accepted.
_EMDB = re.compile(r"EMD[-_](\d+)", re.I)


@dataclass
class Evidence:
    """What a file says about its own origin. Not a verdict."""

    emdb_accession: str | None = None
    labels: tuple[str, ...] = ()
    suggested: Provenance = Provenance.UNKNOWN
    reasons: list[str] = field(default_factory=list)

    @property
    def has_any(self) -> bool:
        return bool(self.reasons)


def gather_evidence(header: MapHeader, path: str | Path | None = None) -> Evidence:
    """Read what the file says about where it came from.

    Looks at the MRC label records — the only provenance information carried
    *inside* a map — and at the filename. Returns a suggestion, never a
    decision: `load_map` will not adopt it without the caller saying so.
    """
    name = Path(path or header.path).name
    evidence = Evidence(labels=header.labels)

    joined = " ".join(header.labels)
    match = _EMDB.search(joined) or _EMDB.search(name)
    if match:
        evidence.emdb_accession = f"EMD-{match.group(1)}"
        evidence.reasons.append(
            f"MRC label identifies this as EMDB deposition {evidence.emdb_accession}. "
            f"EMDB primary maps are usually post-processed, but the deposition "
            f"itself does not say which — check the entry."
        )

    haystack = f"{name} {joined}".lower()

    # Longest matching token wins, not the first in declaration order. These
    # tokens are substrings of one another: 'sharp' occurs inside 'unsharp',
    # and SHARPENED is declared first, so a file named emd_1234_unsharpened.mrc
    # was suggested as *sharpened* — inverting the one thing its depositor had
    # taken the trouble to record. Sorting by length makes the more specific
    # token win whatever order the table is written in; ties keep declaration
    # order, since sorted() is stable.
    candidates = [
        (provenance, token)
        for provenance, tokens in _TOKENS
        for token in tokens
        if token in haystack
    ]
    if candidates:
        provenance, hit = sorted(candidates, key=lambda pair: -len(pair[1]))[0]
        evidence.suggested = provenance
        evidence.reasons.append(
            f"the name or labels contain {hit!r}, which suggests {provenance.value}"
        )

    if not evidence.reasons:
        evidence.reasons.append("nothing in the filename or MRC labels indicates an origin")
    return evidence


# object name -> declared provenance
# Process-global. Correct for a single-session MCP server, which is what
# MCPymol is today; wrong the moment one process serves two sessions. Keyed
# by object name, so two sessions with the same object name would collide.
# Flagged rather than solved — see MOVING.md.
_REGISTRY: dict[str, Provenance] = {}


def declare(obj: str, provenance: Provenance) -> None:
    """Record how ``obj`` came to exist."""
    _REGISTRY[obj] = provenance


def provenance_of(obj: str) -> Provenance:
    """What was declared for ``obj``. UNKNOWN if nothing was.

    Never raises: a volume nobody declared is *unknown*, which is a real answer
    and the one I1 requires be shown rather than glossed over.
    """
    return _REGISTRY.get(obj, Provenance.UNKNOWN)


def forget(obj: str | None = None) -> None:
    """Drop the record for ``obj``, or all of them."""
    if obj is None:
        _REGISTRY.clear()
    else:
        _REGISTRY.pop(obj, None)


def provenance_banner(obj: str) -> str:
    """The block every volume readout carries. Invariant I1 in one function."""
    provenance = provenance_of(obj)
    return f"  Provenance: {provenance.value.upper()}\n  {provenance.caveat}"
