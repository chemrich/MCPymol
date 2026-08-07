"""Comparing two structures: superposition and per-residue deviation.

``align``/``super`` in :mod:`mcpymol.primitives` give you a single RMSD number,
which says a complex moved but not *where*.  ``superposition_view`` superposes
and then colours the mobile structure by how far each residue actually shifted,
so a hinge motion or a loop rearrangement is visible at a glance.
"""

import json
from typing import Annotated

from pydantic import Field

from mcpymol.app import mcp
from mcpymol.bridge import send_request
from mcpymol.pdbtext import distance3d, parse_atoms, residue_order

# Deviations above this are almost always a domain-scale motion rather than
# local jitter; used only to pick a sensible default colour ceiling.
_DEFAULT_DEVIATION_CEILING = 5.0

# How many of the worst-shifted residues to name in the summary.
_TOP_N_SHIFTED = 8


def _parse_ca_coords(pdb_text: str) -> dict[tuple[str, str], tuple[float, float, float]]:
    """Map (chain, resi) -> (x, y, z) for every alpha carbon."""
    return {a.residue_key: (a.x, a.y, a.z) for a in parse_atoms(pdb_text, ca_only=True)}


def _match_residues(
    mobile: dict[tuple[str, str], tuple[float, float, float]],
    target: dict[tuple[str, str], tuple[float, float, float]],
) -> list[tuple[str, str, float]]:
    """Pair residues between two coordinate maps and measure the gap.

    Matches on (chain, resi) first.  Two structures of the same protein
    deposited with different chain naming — very common, e.g. an apo monomer
    labelled A against a complex where it is chain B — would otherwise match
    nothing, so fall back to residue number alone when the strict pass is
    empty and the chains are unambiguous.
    """
    keys = mobile.keys() & target.keys()
    if not keys:
        by_resi_m: dict[str, tuple[float, float, float]] = {}
        by_resi_t: dict[str, tuple[float, float, float]] = {}
        for (_chain, resi), xyz in mobile.items():
            by_resi_m.setdefault(resi, xyz)
        for (_chain, resi), xyz in target.items():
            by_resi_t.setdefault(resi, xyz)
        shared = by_resi_m.keys() & by_resi_t.keys()
        return sorted(
            (("", resi, _distance(by_resi_m[resi], by_resi_t[resi])) for resi in shared),
            key=lambda r: _resi_sort_key(r[1]),
        )

    return sorted(
        (
            (chain, resi, _distance(mobile[(chain, resi)], target[(chain, resi)]))
            for chain, resi in keys
        ),
        key=lambda r: (r[0], _resi_sort_key(r[1])),
    )


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    """Kept as a named helper because the pairing code reads better with it."""
    return distance3d(a, b)


def _resi_sort_key(resi: str) -> tuple[int, str]:
    """Sort 10A after 10 and before 11, without choking on insertion codes."""
    return residue_order(resi)


@mcp.tool()
def superposition_view(
    mobile: Annotated[str, Field(description='Object to move and color (e.g. "1ake").')],
    target: Annotated[
        str, Field(description='Object to superpose onto and leave in place (e.g. "4ake").')
    ],
    method: Annotated[
        str,
        Field(
            description='"super" (default) is sequence-independent and handles low identity or different folds; "align" uses a sequence alignment first and is better for near-identical sequences.'
        ),
    ] = "super",
    max_deviation: Annotated[
        float | None,
        Field(
            description="Angstrom value mapped to full red. Defaults to the largest observed shift, which maximises contrast; set it explicitly to compare two different pairs on one scale."
        ),
    ] = None,
) -> str:
    """
    Superposes two structures and colors the mobile one by per-residue shift.

    An RMSD alone tells you that something moved, not where. This superposes
    ``mobile`` onto ``target``, measures how far each residue's CA ended up
    from its counterpart, and colors the mobile structure blue (unchanged)
    through white to red (most shifted). The target is left as a grey
    reference cartoon.

    Best on two states of the same protein — apo vs holo, open vs closed, a
    mutant against wild type. Residues are paired by chain and residue number,
    falling back to residue number alone when the two use different chain IDs.
    """
    if method not in ("super", "align"):
        return f"Error: method must be 'super' or 'align', got {method!r}."
    if mobile == target:
        return "Error: mobile and target are the same object; nothing to compare."

    fit = send_request(method, args=[mobile, target], timeout=120.0)
    if fit.get("status") == "error":
        return f"Error superposing {mobile} onto {target}: {fit.get('error')}"

    result = fit.get("result")
    rmsd = None
    n_atoms = None
    if isinstance(result, (list, tuple)) and result:
        try:
            rmsd = float(result[0])
            n_atoms = int(result[1]) if len(result) > 1 else None
        except (TypeError, ValueError):
            rmsd = None

    # `super`/`align` put the fit in the object's transformation matrix, so the
    # stored coordinates can still be pre-superposition. Copying each object
    # bakes the current transform into fresh coordinates, which puts both in
    # one frame and makes the numbers below trustworthy.
    tmp_mobile, tmp_target = "_sup_mobile", "_sup_target"
    try:
        for tmp, src in ((tmp_mobile, mobile), (tmp_target, target)):
            send_request("delete", args=[tmp])
            send_request("create", args=[tmp, f"({src}) and polymer and name CA"])

        pdbs = {}
        for tmp in (tmp_mobile, tmp_target):
            res = send_request("get_pdbstr", args=[tmp], timeout=120.0)
            if res.get("status") == "error":
                return f"Error reading coordinates from {tmp}: {res.get('error')}"
            pdbs[tmp] = res.get("result") or ""

        pairs = _match_residues(
            _parse_ca_coords(pdbs[tmp_mobile]), _parse_ca_coords(pdbs[tmp_target])
        )
    finally:
        send_request("delete", args=[tmp_mobile])
        send_request("delete", args=[tmp_target])

    if not pairs:
        return (
            f"Superposed {mobile} onto {target}"
            + (f" (RMSD {rmsd:.2f} A)" if rmsd is not None else "")
            + ", but no residues could be paired for a per-residue comparison. "
            "The two structures share no residue numbering — check they are the "
            "same protein, or compare them with align/super directly."
        )

    deviations = [d for _chain, _resi, d in pairs]
    observed_max = max(deviations)
    ceiling = max_deviation if max_deviation is not None else observed_max
    if ceiling <= 0:
        ceiling = _DEFAULT_DEVIATION_CEILING

    # Push the deviations into the B-factor column in one batched script, the
    # same round-trip-collapsing trick conservation_view uses.
    per_chain: dict[str, dict[str, float]] = {}
    for chain, resi, dev in pairs:
        per_chain.setdefault(chain, {})[resi] = round(dev, 3)

    script = [f"alter ({mobile}), b=0"]
    for chain, mapping in per_chain.items():
        sel = f"({mobile}) and chain {chain}" if chain else f"({mobile})"
        script.append(f"stored.dev = {json.dumps(mapping)}")
        script.append(f"alter {sel}, b=stored.dev.get(resi, 0.0)")
    script.append(f"rebuild {mobile}")
    send_request("do", args=["\n".join(script)], timeout=120.0)

    send_request("hide", args=["everything", mobile])
    send_request("hide", args=["everything", target])
    send_request("show", args=["cartoon", mobile])
    send_request("show", args=["cartoon", target])
    send_request("color", args=["grey70", target])
    send_request("set", args=["cartoon_transparency", "0.6", target])
    send_request(
        "do",
        args=[f"spectrum b, blue_white_red, ({mobile}) and polymer, minimum=0, maximum={ceiling}"],
    )
    send_request("do", args=["bg_color black"])
    send_request("zoom", args=[f"({mobile}) or ({target})"])

    worst = sorted(pairs, key=lambda r: r[2], reverse=True)[:_TOP_N_SHIFTED]
    worst_txt = ", ".join(
        f"{chain + '/' if chain else ''}{resi} ({dev:.1f} A)" for chain, resi, dev in worst
    )
    mean_dev = sum(deviations) / len(deviations)

    header = f"Superposed {mobile} onto {target} with {method}"
    if rmsd is not None:
        header += f": RMSD {rmsd:.2f} A"
        if n_atoms:
            header += f" over {n_atoms} atoms"

    return (
        f"{header}. Compared {len(pairs)} residues — mean shift {mean_dev:.2f} A, "
        f"max {observed_max:.2f} A. {mobile} is coloured blue (unchanged) to red "
        f"(shifted >= {ceiling:.1f} A) with {target} as a grey reference.\n"
        f"Most-shifted residues: {worst_txt}."
    )
