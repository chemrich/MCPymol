"""Per-atom displacement between conformations: the well-supported claim.

Of everything cryo-EM heterogeneity methods produce, **motion is the part that
survives blind testing.** The Flatiron challenge put 41 submissions on the same
data and reported that the molecular motions they identified resembled both each
other and the ground truth; it was the *population distributions* that fell
apart. So this tool is built with confidence, and its legend says what it is
showing rather than hedging it into uselessness.

That confidence is bounded in one specific way. A displacement field is a
statement about **where density moved**, not about how many particles moved that
way. Nothing here reports a population, and the legend says so — the honesty
problem in this corner of the field is not that motion is overclaimed, it is
that motion and population get drawn with the same confidence.

**Half-set uncertainty.** DynaMight estimates deformations twice, on independent
half-sets, and the disagreement between them is the only per-position signal of
how much any individual arrow can be trusted. When that table is supplied the
arrows are coloured by it. When it is not, the report says so loudly rather than
drawing uniformly confident arrows over a field that has none — an arrow looks
exactly as convincing either way, which is the whole problem.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

from mcpymol.wiggles.atoms import (
    Atom,
    count_states,
    fetch_atoms,
    fetch_state_coords,
    residue_selection,
)
from mcpymol.wiggles.bfactors import preservation_note, stash_bfactors
from mcpymol.wiggles.port import PortError, PymolPort, call

#: PyMOL CGO opcodes. Numbers rather than an import because this package does
#: not import from PyMOL — it talks to it over a socket, and the plugin at the
#: other end is what holds a `pymol` module.
_CGO_CYLINDER = 9.0
_CGO_CONE = 27.0

#: Ceiling on arrows drawn. One per residue over a large assembly is an
#: unreadable thicket, and drawing every one would also make the biggest
#: displacements impossible to pick out. Whatever is dropped is stated.
DEFAULT_MAX_ARROWS = 60

#: Below this, a displacement is within the noise of most refinements and an
#: arrow would assert a motion nobody measured.
MIN_ARROW_LENGTH = 0.25

MOTION_LEGEND = (
    "MOTION IS THE WELL-SUPPORTED CLAIM HERE. In a 41-submission blind challenge "
    "the recovered molecular motions resembled both each other and the ground "
    "truth; it was relative populations that most submissions got wrong. So read "
    "these arrows as: density moved this way. Do NOT read them as: this fraction "
    "of particles moved this way. No population is shown, and none can be "
    "inferred from an arrow's length or colour."
)


def _displacement(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))


def per_atom_displacement(
    start: list[tuple[float, float, float]],
    end: list[tuple[float, float, float]],
) -> list[float]:
    """Distance each atom moved between two states.

    Raises:
        PortError: the two states disagree on atom count, which makes pairing
            undefined — the same refusal ``morph_states`` makes, for the same
            reason. An unpaired displacement field is not approximate, it is
            meaningless.
    """
    if len(start) != len(end):
        raise PortError(
            f"the two states hold {len(start)} and {len(end)} atoms, so atoms cannot "
            f"be paired across them and displacement has no defined meaning. If these "
            f"came from independently reconstructed volumes, that is expected — "
            f"compare them with ensemble_spread_view, which needs no pairing."
        )
    return [_displacement(a, b) for a, b in zip(start, end, strict=True)]


def per_residue_displacement(
    atoms: list[Atom], displacement: list[float]
) -> dict[tuple[str, str], float]:
    """Mean displacement within each residue, keyed by (chain, resi)."""
    if len(atoms) != len(displacement):
        raise PortError(
            f"{len(atoms)} atoms but {len(displacement)} displacements — the iterate "
            f"and coordinate orders have diverged"
        )
    acc: dict[tuple[str, str], list[float]] = {}
    for atom, value in zip(atoms, displacement, strict=True):
        acc.setdefault(atom.residue, []).append(value)
    return {key: sum(v) / len(v) for key, v in acc.items()}


def read_uncertainty_table(path: str | Path) -> dict[tuple[str, str], float]:
    """Read a per-residue uncertainty table: ``chain resi value`` per line.

    Deliberately a plain text format rather than a reader for any one method's
    output. DynaMight's half-set estimates are the motivating case, but its
    on-disk conventions are RELION's and change; a three-column table is
    something a user can produce from any method that estimates uncertainty at
    all, and this package does not have to track six formats to accept it.

    Raises:
        ValueError: nothing parseable was found.
    """
    path = Path(path)
    table: dict[tuple[str, str], float] = {}
    for line in path.read_text(errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = re.split(r"[,\t ]+", stripped)
        if len(parts) < 3:
            continue
        try:
            table[(parts[0], parts[1])] = float(parts[2])
        except ValueError:
            continue  # a header row
    if not table:
        raise ValueError(
            f"{path}: no rows of the form 'chain resi value' found. Comments start "
            f"with #; separators may be spaces, tabs or commas."
        )
    return table


def _arrow_cgo(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    rgb: tuple[float, float, float],
    radius: float,
) -> list[float]:
    """One arrow: a cylinder for the shaft, a cone for the head."""
    shaft_end = tuple(s + (e - s) * 0.75 for s, e in zip(start, end, strict=True))
    r, g, b = rgb
    return [
        _CGO_CYLINDER,
        *start,
        *shaft_end,
        radius,
        r,
        g,
        b,
        r,
        g,
        b,
        _CGO_CONE,
        *shaft_end,
        *end,
        radius * 2.2,
        0.0,
        r,
        g,
        b,
        r,
        g,
        b,
        1.0,
        1.0,
    ]


def _uncertainty_colour(value: float, worst: float) -> tuple[float, float, float]:
    """Blue where the half-sets agree, red where they do not."""
    t = 0.0 if worst <= 0 else max(0.0, min(1.0, value / worst))
    return (t, 0.35 * (1 - t), 1.0 - t)


def deformation_view(
    port: PymolPort,
    obj: str,
    *,
    start_state: int = 1,
    end_state: int | None = None,
    arrows: bool = True,
    arrow_scale: float = 1.0,
    max_arrows: int = DEFAULT_MAX_ARROWS,
    as_putty: bool = False,
    uncertainty: dict[tuple[str, str], float] | str | Path | None = None,
    preserve_bfactors: bool = True,
) -> str:
    """Colour a model by how far each residue moved, and draw the motion.

    Args:
        port: A live or fake PyMOL port.
        obj: A multi-state object — an ensemble, a morph, or any model whose
            states are conformations.
        start_state: The state to measure from.
        end_state: The state to measure to. Defaults to the last.
        arrows: Draw CGO arrows from start to end position.
        arrow_scale: Multiply arrow length. Above 1.0 this is an exaggeration
            and the report says so, because an exaggerated arrow is a claim
            about magnitude that the data did not make.
        max_arrows: Ceiling on arrows drawn; the rest are reported, not hidden.
        as_putty: Also scale cartoon width by displacement.
        uncertainty: Per-residue half-set disagreement, as a table or a path to
            one. Colours the arrows when present.
        preserve_bfactors: Stash the original B-factors before overwriting.

    Returns:
        A report: the displacement range, the residues that moved most, what
        the uncertainty is or is not, and the motion legend.

    Raises:
        PortError: fewer than two states, an out-of-range state, or states that
            cannot be paired.
        ValueError: the uncertainty table is unreadable.
    """
    n_states = count_states(port, obj)
    if n_states < 2:
        return (
            f"deformation_view({obj})\n\n"
            f"  REFUSED: object has {n_states} state; a displacement needs two.\n"
            f"  A deformation is a difference — there is nothing here to difference."
        )

    end_state = n_states if end_state is None else end_state
    for label, state in (("start_state", start_state), ("end_state", end_state)):
        if not 1 <= state <= n_states:
            raise PortError(f"{label}={state} is outside the object's 1..{n_states} states")
    if start_state == end_state:
        raise PortError(
            f"start_state and end_state are both {start_state}; a displacement "
            f"between a state and itself is zero everywhere."
        )

    if isinstance(uncertainty, (str, Path)):
        uncertainty = read_uncertainty_table(uncertainty)

    atoms = fetch_atoms(port, obj)
    start_coords = fetch_state_coords(port, obj, start_state)
    end_coords = fetch_state_coords(port, obj, end_state)

    atom_shift = per_atom_displacement(start_coords, end_coords)
    residue_shift = per_residue_displacement(atoms, atom_shift)

    values = sorted(residue_shift.values())
    lo, hi, median = values[0], values[-1], values[len(values) // 2]

    stashed = stash_bfactors(obj, atoms) if preserve_bfactors else 0
    for (chain, resi), value in residue_shift.items():
        call(port, "alter", residue_selection(obj, chain, resi), f"b={value:.4f}")
    call(port, "spectrum", "b", "blue_white_red", obj, minimum=0, maximum=round(hi, 4))
    if as_putty:
        call(port, "show", "cartoon", obj)
        call(port, "set", "cartoon_putty_scale_min", 0.4)
        call(port, "set", "cartoon_putty_scale_max", 3.0)
        call(port, "cartoon", "putty", obj)

    drawn, skipped_short, dropped = 0, 0, 0
    arrow_name = f"{obj}_arrows"
    if arrows:
        drawn, skipped_short, dropped = _draw_arrows(
            port,
            obj,
            arrow_name,
            atoms,
            start_coords,
            end_coords,
            arrow_scale,
            max_arrows,
            uncertainty,
        )

    return _report(
        obj=obj,
        start_state=start_state,
        end_state=end_state,
        n_residues=len(residue_shift),
        n_atoms=len(atoms),
        lo=lo,
        hi=hi,
        median=median,
        residue_shift=residue_shift,
        stashed=stashed,
        as_putty=as_putty,
        arrows=arrows,
        arrow_name=arrow_name,
        drawn=drawn,
        skipped_short=skipped_short,
        dropped=dropped,
        arrow_scale=arrow_scale,
        uncertainty=uncertainty,
    )


def _draw_arrows(
    port: PymolPort,
    obj: str,
    arrow_name: str,
    atoms: list[Atom],
    start_coords: list[tuple[float, float, float]],
    end_coords: list[tuple[float, float, float]],
    arrow_scale: float,
    max_arrows: int,
    uncertainty: dict[tuple[str, str], float] | None,
) -> tuple[int, int, int]:
    """One arrow per residue, at its CA. Returns (drawn, too_short, dropped)."""
    candidates: list[tuple[float, tuple[str, str], tuple, tuple]] = []
    for atom, start, end in zip(atoms, start_coords, end_coords, strict=True):
        if atom.name != "CA":
            continue
        length = _displacement(start, end)
        candidates.append((length, atom.residue, start, end))

    if not candidates:  # not a protein, or no CAs — fall back to every atom
        for atom, start, end in zip(atoms, start_coords, end_coords, strict=True):
            candidates.append((_displacement(start, end), atom.residue, start, end))

    long_enough = [c for c in candidates if c[0] >= MIN_ARROW_LENGTH]
    skipped_short = len(candidates) - len(long_enough)

    long_enough.sort(key=lambda c: -c[0])
    dropped = max(0, len(long_enough) - max_arrows)
    chosen = long_enough[:max_arrows]

    worst = max((uncertainty or {}).values(), default=0.0)
    cgo: list[float] = []
    for length, residue, start, end in chosen:
        tip = (
            start[0] + (end[0] - start[0]) * arrow_scale,
            start[1] + (end[1] - start[1]) * arrow_scale,
            start[2] + (end[2] - start[2]) * arrow_scale,
        )
        if uncertainty is not None:
            rgb = _uncertainty_colour(uncertainty.get(residue, worst), worst)
        else:
            rgb = (1.0, 0.55, 0.0)
        cgo.extend(_arrow_cgo(start, tip, rgb, radius=0.12 + 0.04 * min(length, 5.0)))

    if cgo:
        call(port, "delete", arrow_name)
        call(port, "load_cgo", cgo, arrow_name)
    return len(chosen), skipped_short, dropped


def _report(**kw) -> str:
    obj, hi = kw["obj"], kw["hi"]
    biggest = sorted(kw["residue_shift"].items(), key=lambda kv: -kv[1])[:5]
    biggest_text = ", ".join(f"{c}/{r} {v:.2f}Å" for (c, r), v in biggest)

    lines = [
        f"deformation_view({obj}, state {kw['start_state']} -> {kw['end_state']})",
        "",
        f"  {kw['n_atoms']} atoms, {kw['n_residues']} residues.",
        f"  Per-residue displacement: {kw['lo']:.2f} – {hi:.2f} Å (median {kw['median']:.2f} Å)",
        f"  Moved most: {biggest_text}",
        "",
        f"  Coloured blue (still) → red ({hi:.2f} Å), scaled to this object.",
        preservation_note(obj, kw["stashed"]),
    ]
    if kw["as_putty"]:
        lines.append("  Tube width also tracks displacement (putty).")

    if kw["arrows"]:
        lines += ["", f"  Arrows: {kw['drawn']} drawn as `{kw['arrow_name']}`."]
        if kw["skipped_short"]:
            lines.append(
                f"    {kw['skipped_short']} residue(s) moved less than "
                f"{MIN_ARROW_LENGTH} Å and got no arrow — that is within the noise of "
                f"most refinements, and an arrow there would assert a motion nobody "
                f"measured."
            )
        if kw["dropped"]:
            lines.append(
                f"    {kw['dropped']} further residue(s) moved enough to warrant an "
                f"arrow but were NOT drawn — the ceiling is {kw['drawn']}. This is a "
                f"truncation of the largest movers, not the whole field; raise "
                f"max_arrows to see the rest."
            )
        if kw["arrow_scale"] != 1.0:
            lines.append(
                f"    EXAGGERATED {kw['arrow_scale']:g}x. Arrow lengths are no longer "
                f"the measured displacement; the colours and the numbers above still "
                f"are."
            )

    lines.append("")
    if kw["uncertainty"] is not None:
        lines += [
            f"  Arrows coloured by half-set uncertainty over {len(kw['uncertainty'])} residues:",
            "  blue where the two half-set estimates agree, red where they disagree.",
            "  A red arrow is a motion one half of the data does not confirm.",
        ]
    else:
        lines += [
            "  NO UNCERTAINTY SUPPLIED, so every arrow is drawn with identical",
            "  confidence — and they do not have identical confidence. Methods that",
            "  estimate deformation on independent half-sets (DynaMight does) can",
            "  tell you which arrows the data support; without that table this view",
            "  cannot, and a uniformly orange field is not evidence that the field is",
            "  uniformly reliable. Pass uncertainty= with a 'chain resi value' table.",
        ]

    lines += ["", MOTION_LEGEND]
    return "\n".join(lines)
