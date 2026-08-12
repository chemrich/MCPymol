"""Compositional occupancy — sense 2 — from an explicit table, never inferred.

"Occupancy" names two incompatible quantities, and this module exists because
conflating them produces a viewer that lies while looking fine.

**Sense 1** is the per-atom crystallographic ``q`` in the coordinate file: what
fraction of copies have that atom at that position. That is
:mod:`mcpymol.wiggles.occupancy`, and it is already addressable in PyMOL.

**Sense 2**, here, is the fraction of *imaged particles* in which a subunit,
domain or ligand is present at all. A ligand at 40% in this sense means 40% of
complexes had it bound. It lives in method output — class populations, part
presence probabilities, a normalised local map filter — and **not in the
coordinate file.** Nothing in PyMOL represents it today.

A model can be ``q = 1.0`` at every atom, a perfectly ordinary single-conformer
deposition, while the subunit it belongs to is present in half the particles.
Sense 1 says "fully occupied". Sense 2 says "half there". Both are correct, and
they answer different questions.

**So this module never reads ``q``, and that is not an oversight to be tidied up
later.** Deriving sense 2 from sense 1 would take a number that means one thing
and present it as another, and the render would look identical either way. The
table is supplied by the caller or there is no view. A test asserts no atom
property is ever fetched here.

The two tools stay separately named forever. Every legend either one draws
states which sense it is showing — see the compendium entry
`occupancy-two-senses`.
"""

from __future__ import annotations

import re
from pathlib import Path

from mcpymol.wiggles.port import PortError, PymolPort, call

SENSE_2_LEGEND = (
    "THIS IS OCCUPANCY IN SENSE 2 — COMPOSITIONAL. Each value is the fraction of "
    "imaged particles in which that part is present at all. It is NOT the "
    "per-atom crystallographic occupancy q, it was NOT derived from q, and the "
    "two are not interchangeable: a model can be q=1.0 everywhere while a subunit "
    "is present in half the particles. Both statements are true and they answer "
    "different questions. For sense 1, use occupancy_view."
)

PROVENANCE_NOTE = (
    "These numbers came from the table you supplied, and this tool has no way to "
    "check them against the structure. Compositional occupancy is estimated — by "
    "class population, part-presence probability, or a normalised local map "
    "filter — and the estimate carries the uncertainty of whichever method made "
    "it. Nothing here narrows that."
)


def parse_composition_table(
    source: dict[str, float] | str | Path,
) -> dict[str, float]:
    """Read a ``selection -> presence fraction`` table.

    Accepts a dict, an inline ``"chain A=0.4, chain B=1.0"`` string, or a path
    to a file of ``selection<TAB>fraction`` lines. Selections are PyMOL
    selection expressions, so they may contain spaces — which is why the
    delimiter is ``=`` inline and a tab in a file, rather than whitespace.

    Raises:
        ValueError: nothing parseable, or a fraction outside [0, 1].
    """
    if isinstance(source, dict):
        table = {str(k): float(v) for k, v in source.items()}
    else:
        text = Path(source).read_text(errors="replace") if _looks_like_path(source) else str(source)
        table = _parse_text(text)

    if not table:
        raise ValueError(
            "no entries found. Give a dict, an inline 'chain A=0.4, chain B=1.0' "
            "string, or a file of 'selection<TAB>fraction' lines."
        )

    bad = {sel: v for sel, v in table.items() if not 0.0 <= v <= 1.0}
    if bad:
        raise ValueError(
            f"compositional occupancy is a fraction of particles and must lie in "
            f"[0, 1]; got {bad}. A value above 1 usually means a raw particle count "
            f"or a percentage — divide by the total, or by 100."
        )
    return table


def _looks_like_path(source: str | Path) -> bool:
    if isinstance(source, Path):
        return True
    return "=" not in source and "\n" not in source and Path(source).exists()


def _parse_text(text: str) -> dict[str, float]:
    table: dict[str, float] = {}
    for chunk in re.split(r"[\n,]", text):
        stripped = chunk.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "=" in stripped:
            selection, _, raw = stripped.partition("=")
        elif "\t" in stripped:
            selection, _, raw = stripped.partition("\t")
        else:
            selection, _, raw = stripped.rpartition(" ")
        selection, raw = selection.strip(), raw.strip()
        if not selection or not raw:
            continue
        try:
            table[selection] = float(raw.rstrip("%")) / (100.0 if raw.endswith("%") else 1.0)
        except ValueError:
            continue
    return table


def composition_view(
    port: PymolPort,
    obj: str,
    table: dict[str, float] | str | Path,
    *,
    transparency: bool = True,
    label: bool = True,
    palette: tuple[str, str] = ("red", "skyblue"),
) -> str:
    """Colour parts of ``obj`` by how often they are present across particles.

    Args:
        port: A live or fake PyMOL port.
        obj: The object the table's selections apply to.
        table: ``selection -> fraction``, as a dict, an inline string, or a path.
        transparency: Also make rarely-present parts transparent in proportion,
            so a half-present subunit reads as half there.
        label: Write the fraction next to each part.
        palette: Colours for (rare, always present). The default runs red for
            rarely-present to blue for always-present, which is the opposite of
            the usual "red = more" reading, so the legend states the direction.

    Returns:
        A report naming sense 2 explicitly, per-part fractions, and any
        selection that matched no atoms.

    Raises:
        PortError: a selection matched nothing, so its colour never landed.
        ValueError: the table is unparseable or holds a value outside [0, 1].
    """
    parsed = parse_composition_table(table)

    # An empty selection is the dangerous case: PyMOL accepts it, nothing is
    # coloured, and the render looks like a structure that is simply all
    # present. Check before drawing rather than discovering it visually.
    empty: list[str] = []
    counts: dict[str, int] = {}
    for selection in parsed:
        scoped = f"({obj}) and ({selection})"
        result = port.query("count_atoms", scoped)
        count = result if isinstance(result, int) else 0
        counts[selection] = count
        if count == 0:
            empty.append(selection)

    if empty:
        raise PortError(
            f"{len(empty)} selection(s) in the table matched no atoms in {obj!r}: "
            f"{empty}. PyMOL would have accepted these silently, coloured nothing, "
            f"and left a render that looks like a fully-present structure. Check the "
            f"chain identifiers or the object name."
        )

    rare, present = palette
    for selection, fraction in parsed.items():
        scoped = f"({obj}) and ({selection})"
        call(port, "set_color", f"_wg_comp_{_slug(selection)}", _blend(rare, present, fraction))
        call(port, "color", f"_wg_comp_{_slug(selection)}", scoped)
        if transparency:
            # A part present in 40% of particles is drawn 60% transparent.
            call(port, "set", "cartoon_transparency", round(1.0 - fraction, 3), scoped)
            call(port, "set", "transparency", round(1.0 - fraction, 3), scoped)
        if label:
            call(port, "label", f"({scoped}) and name CA and rank 0", f'"{fraction:.0%}"')

    return _report(obj, parsed, counts, transparency)


def _slug(selection: str) -> str:
    return re.sub(r"\W+", "_", selection).strip("_").lower() or "sel"


def _blend(rare: str, present: str, fraction: float) -> list[float]:
    """Interpolate between two named colours. Named colours are resolved by
    PyMOL, so this uses their standard RGB rather than asking — the palette is
    a default, and a caller wanting exact control passes their own colours."""
    known = {
        "red": (1.0, 0.0, 0.0),
        "skyblue": (0.34, 0.63, 0.83),
        "blue": (0.0, 0.0, 1.0),
        "white": (1.0, 1.0, 1.0),
        "grey": (0.5, 0.5, 0.5),
        "gray": (0.5, 0.5, 0.5),
        "yellow": (1.0, 1.0, 0.0),
        "orange": (1.0, 0.5, 0.0),
        "green": (0.0, 1.0, 0.0),
    }
    a = known.get(rare, (1.0, 0.0, 0.0))
    b = known.get(present, (0.34, 0.63, 0.83))
    return [round(x + (y - x) * fraction, 4) for x, y in zip(a, b, strict=True)]


def _report(
    obj: str,
    table: dict[str, float],
    counts: dict[str, int],
    transparency: bool,
) -> str:
    ordered = sorted(table.items(), key=lambda kv: kv[1])
    rows = [
        f"    {fraction:>6.0%}  {counts[selection]:>6} atoms  {selection}"
        for selection, fraction in ordered
    ]
    lowest, highest = ordered[0], ordered[-1]

    lines = [
        f"composition_view({obj})",
        "",
        f"  {len(table)} part(s), by fraction of particles containing them:",
        *rows,
        "",
        f"  Least present: {lowest[0]} at {lowest[1]:.0%}.  "
        f"Most: {highest[0]} at {highest[1]:.0%}.",
        "",
        "  Coloured red (rarely present) → blue (always present). Low values are",
        "  the interesting ones here, so the scale runs opposite to the usual",
        "  reading of red as 'more'.",
    ]
    if transparency:
        lines.append("  Transparency tracks absence too: a part in 40% of particles is drawn")
        lines.append("  60% transparent, so it reads as partly there rather than solid.")

    lines += ["", "  " + PROVENANCE_NOTE, "", SENSE_2_LEGEND]
    return "\n".join(lines)
