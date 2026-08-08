"""Occupancy views — **sense 1 only**, per-atom crystallographic ``q``.

These tools read the ``occupancy`` column of a coordinate file and nothing
else. They never infer, estimate or display compositional occupancy — what
fraction of imaged particles contain a subunit. That is sense 2, it lives in
method output rather than the model, and it belongs to ``composition_view`` in
tier 3.

The distinction is the single most important design decision in SPEC.md, and
the reason is that the failure is invisible: a model can be ``q = 1.0``
everywhere while its subunit is present in half the particles. Both statements
are true, they answer different questions, and a render that conflates them
looks entirely normal. So every legend these functions emit names the sense.

See the Wiggles compendium entry `occupancy-two-senses`.
"""

from __future__ import annotations

from mcpymol.wiggles.atoms import Atom, altloc_groups, fetch_atoms, group_by_residue
from mcpymol.wiggles.bfactors import preservation_note, stash_bfactors
from mcpymol.wiggles.port import PymolPort, call

# Occupancy at or above this is "full" for reporting purposes. Occupancies are
# stored as float32 and refined, so exact 1.0 comparison mislabels ordinary
# fully-occupied atoms.
FULL_OCCUPANCY = 0.999

# qFit removes conformers below 10% occupancy during refinement, so a
# multiconformer model's conformer count is partly a threshold artifact. Worth
# saying when low-occupancy atoms are present.
QFIT_FLOOR = 0.10

# Distinct colours for altloc groups, in assignment order.
_ALTLOC_COLOURS = (
    "skyblue",
    "salmon",
    "palegreen",
    "wheat",
    "lightpink",
    "paleyellow",
    "lightblue",
    "lightorange",
)

SENSE_1_LEGEND = (
    "Occupancy shown is SENSE 1: per-atom crystallographic occupancy (q) read "
    "from the coordinate file. It is NOT compositional occupancy — the fraction "
    "of imaged particles containing this subunit. A model can be fully occupied "
    "here and half-present in that sense, simultaneously and correctly."
)


def _stats(atoms: list[Atom]) -> dict[str, float | int]:
    qs = [a.q for a in atoms]
    partial = [a for a in atoms if a.q < FULL_OCCUPANCY]
    return {
        "n_atoms": len(atoms),
        "n_partial": len(partial),
        "n_below_qfit_floor": sum(1 for a in atoms if 0.0 < a.q < QFIT_FLOOR),
        "q_min": min(qs),
        "q_max": max(qs),
    }


def occupancy_view(
    port: PymolPort,
    obj: str,
    *,
    preserve_bfactors: bool = True,
) -> str:
    """Colour ``obj`` by per-atom occupancy, de-emphasising partial atoms.

    Occupancy is pushed into the B-factor column so PyMOL's ``spectrum`` can
    colour by it — the same mechanism MCPymol's ``conservation_view`` uses.
    That is destructive to B-factors, so by default the originals are saved
    (see :mod:`mcpymol.wiggles.bfactors`) and the report says how to restore them.

    If nothing in the object is partially occupied, this says so plainly and
    does not colour. A spectrum across a constant value is a rainbow that means
    nothing, and a viewer that renders one is inviting a false reading.

    Args:
        port: A live or fake PyMOL port.
        obj: Object or selection name.
        preserve_bfactors: Save B-factors before overwriting, so
            ``restore_bfactors`` can put them back.

    Returns:
        A report, always ending with the sense-1 legend.
    """
    atoms = fetch_atoms(port, obj)
    st = _stats(atoms)
    lines = [f"occupancy_view({obj})", ""]

    if st["n_partial"] == 0:
        lines += [
            f"  All {st['n_atoms']} atoms are fully occupied (q >= {FULL_OCCUPANCY}).",
            "  Nothing to show: no colouring applied.",
            "",
            "  A spectrum over a constant value would be a rainbow that means",
            "  nothing. If you expected partial occupancy here, check that the",
            "  model actually carries alternate conformations — see altloc_view.",
            "",
            SENSE_1_LEGEND,
        ]
        return "\n".join(lines)

    stashed = stash_bfactors(obj, atoms) if preserve_bfactors else 0
    call(port, "alter", obj, "b=q")
    call(port, "spectrum", "b", "red_white_blue", obj, minimum=0, maximum=1)
    # Partial atoms shown as sticks so they read as "modelled alternative"
    # rather than as the single truth. Note this must NOT exclude protein:
    # side-chain alternates are the main thing this view exists to show.
    call(port, "show", "sticks", f"({obj}) and q<{FULL_OCCUPANCY}")

    lines += [
        f"  {st['n_partial']} of {st['n_atoms']} atoms partially occupied (q < {FULL_OCCUPANCY}).",
        f"  q range: {st['q_min']:.3f} – {st['q_max']:.3f}",
        "  Coloured red (q=0) → white (q=0.5) → blue (q=1), fixed scale.",
        "",
    ]
    if st["n_below_qfit_floor"]:
        lines += [
            f"  {st['n_below_qfit_floor']} atoms sit below q={QFIT_FLOOR:g}. qFit removes",
            "  conformers under that during refinement, so their presence here",
            "  suggests a different refinement protocol — worth checking.",
            "",
        ]
    if preserve_bfactors:
        lines += [preservation_note(obj, stashed), ""]
    else:
        lines += ["  WARNING: B-factors overwritten and not preserved.", ""]

    lines.append(SENSE_1_LEGEND)
    return "\n".join(lines)


def altloc_view(port: PymolPort, obj: str, *, label: bool = True) -> str:
    """Show every alternate conformation at once, one colour per altloc group.

    Occupancies go into the labels rather than being implied by the colouring,
    because the occupancy *is* the measurement — drawing all altlocs in equal
    weight discards it, which is the commonest way this view is got wrong.

    Args:
        port: A live or fake PyMOL port.
        obj: Object or selection name.
        label: Label one atom per alternate residue with its altloc and occupancy.

    Returns:
        A report, always ending with the sense-1 legend.
    """
    atoms = fetch_atoms(port, obj)
    groups = altloc_groups(atoms)
    lines = [f"altloc_view({obj})", ""]

    if not groups:
        lines += [
            f"  No alternate conformations: all {len(atoms)} atoms have a blank altloc.",
            "  Nothing to show.",
            "",
            "  This is an ordinary single-conformer model. It says nothing about",
            "  whether the molecule is heterogeneous — only that this file does",
            "  not model it.",
            "",
            SENSE_1_LEGEND,
        ]
        return "\n".join(lines)

    call(port, "hide", "everything", obj)
    call(port, "show", "cartoon", obj)
    call(port, "color", "grey70", obj)

    per_group: list[str] = []
    for i, alt in enumerate(groups):
        colour = _ALTLOC_COLOURS[i % len(_ALTLOC_COLOURS)]
        sel = f"({obj}) and alt {alt}"
        call(port, "show", "sticks", sel)
        call(port, "color", colour, sel)
        members = [a for a in atoms if a.alt == alt]
        occ = sorted({round(a.q, 3) for a in members})
        occ_text = ", ".join(f"{q:g}" for q in occ[:4]) + (" …" if len(occ) > 4 else "")
        per_group.append(
            f"  alt {alt}: {colour:<12} {len(members):>5} atoms   occupancy {occ_text}"
        )

    if label:
        # One label per alternate residue — labelling every atom is unreadable.
        for alt in groups:
            call(
                port,
                "label",
                f"({obj}) and alt {alt} and name CA",
                f'"{alt} %.2f" % q',
            )

    residues = group_by_residue([a for a in atoms if a.alt.strip()])
    lines += [
        f"  {len(groups)} altloc group(s) across {len(residues)} residue(s).",
        "",
        *per_group,
        "",
        "  Backbone in grey; alternates as coloured sticks.",
    ]
    if label:
        lines.append("  CA atoms labelled with altloc and occupancy.")
    lines += [
        "",
        "  Occupancies are shown rather than encoded in the colouring: drawing",
        "  every alternate at equal weight would discard the measurement.",
        "",
        SENSE_1_LEGEND,
    ]
    return "\n".join(lines)
