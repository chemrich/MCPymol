"""Q-score — per-atom resolvability, read from a wwPDB validation file.

Q-score compares the map values around each atom against a Gaussian-like
reference for a well-resolved atom: 1 is a perfect match, values near 0 mean
the local density does not look like a resolved atom at all (Pintilie *et al.*,
*Nat Methods* 17:328–334, 2020).

The reason this is the cheapest useful tool in tier 1 is that **the values are
already published**. wwPDB validation of 3DEM entries with both a model and a
map has included Q-score since September 2023, with per-residue and per-chain
averages written into the validation mmCIF and XML. So this needs no map, no
network and no computation — it parses a file the archive already distributes
and puts it on the structure.

A global resolution number is one claim about a whole map, and readers apply it
uniformly to every atom they see. Q-score breaks that apart: it shows which
parts of the model the density actually supports.

See the Wiggles compendium entry `qscore`.
"""

from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from pathlib import Path

from mcpymol.wiggles.atoms import (
    fetch_atoms,
    group_by_residue,
    residue_clause,
    residue_selection,
)
from mcpymol.wiggles.bfactors import preservation_note, stash_bfactors
from mcpymol.wiggles.port import PortError, PymolPort, call

# wwPDB has spelled this differently across releases; accept what we have seen
# rather than failing on a file that plainly contains the data.
_Q_ATTRS = ("Q_score", "Q-score", "qscore", "q_score")

# Colour for residues the validation file says nothing about. Deliberately NOT
# on the Q-score scale: absent data is not a low score, and colouring it red
# would assert something the file does not say.
NO_DATA_COLOUR = "grey50"

QSCORE_LEGEND = (
    "Q-score is a map–MODEL agreement metric, not a property of the density "
    "alone. A wrongly-placed atom sitting in someone else's density can score "
    "well, and a correctly-placed atom in a poorly-sharpened region can score "
    "badly. It also assumes the map is correctly sampled — if the voxel size is "
    "wrong, so is this (see map_info)."
)


def _open_maybe_gzip(path: Path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")


def parse_validation_xml(path: str | Path) -> dict[tuple[str, str], float]:
    """Read per-residue Q-scores from a wwPDB validation XML file.

    Handles ``*_validation.xml`` and ``*.xml.gz``. Residues carrying no
    Q-score attribute are omitted from the result rather than defaulted —
    absence is information, and a caller must be able to tell it from a low
    score.

    Returns:
        ``{(chain, resnum): q_score}``.

    Raises:
        ValueError: not parseable as XML, or contains no Q-score at all.
    """
    path = Path(path)
    try:
        with _open_maybe_gzip(path) as fh:
            tree = ET.parse(fh)
    except ET.ParseError as exc:
        raise ValueError(f"{path}: not valid XML: {exc}") from exc

    scores: dict[tuple[str, str], float] = {}
    subgroups = 0
    for element in tree.getroot().iter("ModelledSubgroup"):
        subgroups += 1
        raw = next((element.get(a) for a in _Q_ATTRS if element.get(a) is not None), None)
        if raw is None:
            continue
        chain = element.get("chain") or element.get("said") or ""
        resnum = element.get("resnum")
        if resnum is None:
            continue
        try:
            scores[(str(chain), str(resnum))] = float(raw)
        except ValueError:
            continue

    if not scores:
        raise ValueError(
            f"{path}: no Q-score found in {subgroups} residue record(s).\n"
            f"Two reasons this happens, both normal:\n"
            f"  - the entry is not 3DEM. X-ray and NMR reports never carry it.\n"
            f"  - the entry is 3DEM but was validated before the September 2023\n"
            f"    rollout, and older reports were not regenerated. Confirmed:\n"
            f"    7A4M (EM, 2020) has none; 9C0K (EM, 2024) has 1123.\n"
            f"There is no way to compute Q-score from this file — it needs the\n"
            f"map and the model together (MapQ, or Phenix)."
        )
    return scores


def qscore_view(
    port: PymolPort,
    obj: str,
    validation_path: str | Path,
    *,
    preserve_bfactors: bool = True,
) -> str:
    """Colour ``obj`` by per-residue Q-score from a wwPDB validation file.

    Residues with no Q-score in the file are coloured ``grey50`` and counted
    separately. They are *not* given a score of zero — absent data is not a bad
    score, and putting it on the same colour scale would assert something the
    archive never said.

    Args:
        port: A live or fake PyMOL port.
        obj: Object or selection name.
        validation_path: A wwPDB ``*_validation.xml`` (optionally gzipped).
        preserve_bfactors: Save B-factors before overwriting, so
            ``restore_bfactors`` can put them back.

    Returns:
        A report, always ending with the Q-score legend.
    """
    scores = parse_validation_xml(validation_path)
    atoms = fetch_atoms(port, obj)
    residues = group_by_residue(atoms)

    matched = {key: scores[key] for key in residues if key in scores}
    missing = [key for key in residues if key not in scores]

    if not matched:
        raise PortError(
            f"no residue in {obj!r} matched any of the {len(scores)} scored "
            f"residues in {Path(validation_path).name}. Chain or numbering "
            f"probably differs — validation files use author numbering."
        )

    stashed = stash_bfactors(obj, atoms) if preserve_bfactors else 0
    call(port, "alter", obj, "b=0")
    for (chain, resi), value in matched.items():
        call(port, "alter", residue_selection(obj, chain, resi), f"b={value:.4f}")

    scored_sel = " or ".join(residue_clause(c, r) for c, r in matched)
    call(
        port,
        "spectrum",
        "b",
        "red_yellow_green",
        f"({obj}) and ({scored_sel})",
        minimum=0,
        maximum=1,
    )
    if missing:
        missing_sel = " or ".join(residue_clause(c, r) for c, r in missing)
        call(port, "color", NO_DATA_COLOUR, f"({obj}) and ({missing_sel})")

    values = sorted(matched.values())
    mean = sum(values) / len(values)
    median = values[len(values) // 2]
    worst = sorted(matched.items(), key=lambda kv: kv[1])[:5]
    worst_text = ", ".join(f"{c}/{r} {v:.2f}" for (c, r), v in worst)

    lines = [
        f"qscore_view({obj})",
        "",
        f"  {len(matched)} of {len(residues)} residues scored"
        + (f"; {len(missing)} unscored" if missing else ""),
        f"  Q-score: {values[0]:.2f} – {values[-1]:.2f}  (mean {mean:.2f}, median {median:.2f})",
        f"  Least resolvable: {worst_text}",
        "",
        "  Coloured red (Q=0, unresolvable) → yellow → green (Q=1), fixed scale.",
    ]
    if missing:
        lines += [
            f"  {len(missing)} residues have no Q-score and are {NO_DATA_COLOUR} —",
            "  absent data, deliberately off the scale rather than scored 0.",
        ]

    negative = [key for key, value in matched.items() if value < 0]
    if negative:
        lines += [
            f"  {len(negative)} residues score BELOW ZERO. Q-score is a similarity to a",
            "  Gaussian, so a negative value means the density around those atoms is",
            "  anti-correlated with what a resolved atom looks like — worse than",
            "  absent. They clamp to the bottom of the colour scale; check whether",
            "  the model is placed correctly there.",
        ]
    if preserve_bfactors:
        lines.append(preservation_note(obj, stashed))
    lines += ["", QSCORE_LEGEND]
    return "\n".join(lines)
