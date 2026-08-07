"""Quantitative analysis: what touches what, and how tightly.

The view presets in :mod:`mcpymol.views` *draw* interactions; these tools
*report* them. A picture of a binding site tells you there are contacts; this
tells you which residues, at what distance, and of what kind — the numbers
that go in a figure legend or a methods section.

All geometry here uses heavy-atom distance criteria, because crystal
structures usually have no hydrogens. That makes the classifications
conservative and reproducible, but it also means an "H-bond" here is a
donor-acceptor pair with plausible geometry, not one verified against an
actual hydrogen position.
"""

import math
from typing import Annotated

from pydantic import Field

from mcpymol.app import mcp
from mcpymol.bridge import send_request
from mcpymol.pdbtext import Atom, distance3d, parse_atoms, residue_order

# ── Interaction criteria ─────────────────────────────────────────────────────
#
# Heavy-atom cutoffs in Angstrom, chosen to match what the structural biology
# literature conventionally reports:
#   hydrogen bond   N/O/S donor-acceptor pair    <= 3.5   (Baker & Hubbard)
#   salt bridge     charged N ... charged O      <= 4.0   (Barlow & Thornton)
#   hydrophobic     C/S ... C/S                  <= 4.5
#   pi-stacking     aromatic ring centroids      <= 5.5
HBOND_MAX = 3.5
SALT_BRIDGE_MAX = 4.0
HYDROPHOBIC_MAX = 4.5
AROMATIC_CENTROID_MAX = 5.5

# Interplanar angle below this reads as parallel (sandwich) stacking, above
# its complement as edge-to-face (T-shaped).
STACKING_PARALLEL_MAX_ANGLE = 30.0
STACKING_TSHAPED_MIN_ANGLE = 60.0

_POLAR_ELEMENTS = {"N", "O"}
_HYDROPHOBIC_ELEMENTS = {"C", "S"}

# Formally charged sidechain tips.
_ANIONIC_ATOMS = {
    ("ASP", "OD1"),
    ("ASP", "OD2"),
    ("GLU", "OE1"),
    ("GLU", "OE2"),
}
_CATIONIC_ATOMS = {
    ("LYS", "NZ"),
    ("ARG", "NE"),
    ("ARG", "NH1"),
    ("ARG", "NH2"),
    ("HIS", "ND1"),
    ("HIS", "NE2"),
}

# Aromatic ring atom sets, per residue. TRP has two fused rings.
_AROMATIC_RINGS: dict[str, list[tuple[str, ...]]] = {
    "PHE": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "TYR": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "HIS": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "TRP": [
        ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
        ("CG", "CD1", "NE1", "CE2", "CD2"),
    ],
}

# Ranked most to least specific: a contact is reported as the strongest
# interaction its geometry supports.
_INTERACTION_ORDER = [
    "salt bridge",
    "hydrogen bond",
    "hydrophobic",
    "polar contact",
    "contact",
]


def _atom_distance(a: Atom, b: Atom) -> float:
    return distance3d((a.x, a.y, a.z), (b.x, b.y, b.z))


def _neighbour_pairs(
    atoms_a: list[Atom], atoms_b: list[Atom], cutoff: float
) -> list[tuple[Atom, Atom, float]]:
    """Every cross-set atom pair within ``cutoff``, via a uniform grid.

    Bucketing by ``cutoff``-sized cells means each atom only tests the 27
    surrounding cells rather than the whole opposite set, which keeps this
    linear in atom count instead of quadratic — the difference between a
    ligand pocket and a whole ribosome interface being feasible.
    """
    if not atoms_a or not atoms_b:
        return []

    cell = max(cutoff, 1e-6)
    grid: dict[tuple[int, int, int], list[Atom]] = {}
    for atom in atoms_b:
        key = (int(atom.x // cell), int(atom.y // cell), int(atom.z // cell))
        grid.setdefault(key, []).append(atom)

    pairs: list[tuple[Atom, Atom, float]] = []
    offsets = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for a in atoms_a:
        base = (int(a.x // cell), int(a.y // cell), int(a.z // cell))
        for di, dj, dk in offsets:
            for b in grid.get((base[0] + di, base[1] + dj, base[2] + dk), ()):
                d = _atom_distance(a, b)
                if d <= cutoff:
                    pairs.append((a, b, d))
    return pairs


def _is_anionic(atom: Atom) -> bool:
    return (atom.resn, atom.name) in _ANIONIC_ATOMS


def _is_cationic(atom: Atom) -> bool:
    return (atom.resn, atom.name) in _CATIONIC_ATOMS


def _classify(a: Atom, b: Atom, distance: float) -> str:
    """Name the strongest interaction this atom pair's geometry supports."""
    charged = (_is_anionic(a) and _is_cationic(b)) or (_is_cationic(a) and _is_anionic(b))
    if charged and distance <= SALT_BRIDGE_MAX:
        return "salt bridge"

    both_polar = a.element in _POLAR_ELEMENTS and b.element in _POLAR_ELEMENTS
    if both_polar and distance <= HBOND_MAX:
        return "hydrogen bond"

    if (
        a.element in _HYDROPHOBIC_ELEMENTS
        and b.element in _HYDROPHOBIC_ELEMENTS
        and distance <= HYDROPHOBIC_MAX
    ):
        return "hydrophobic"

    if both_polar:
        return "polar contact"
    return "contact"


def _ring_systems(atoms: list[Atom]) -> list[tuple[tuple[str, str], str, tuple, tuple]]:
    """Aromatic rings as (residue_key, label, centroid, unit normal).

    Only the standard aromatic amino acids are recognised. Ring perception for
    an arbitrary ligand needs bond orders, which a PDB dump does not carry —
    so ligand aromatics are not detected here and their contacts are reported
    as hydrophobic instead.
    """
    by_residue: dict[tuple[str, str], list[Atom]] = {}
    for atom in atoms:
        if atom.resn in _AROMATIC_RINGS:
            by_residue.setdefault(atom.residue_key, []).append(atom)

    rings = []
    for key, residue_atoms in by_residue.items():
        positions = {a.name: (a.x, a.y, a.z) for a in residue_atoms}
        resn = residue_atoms[0].resn
        for ring_names in _AROMATIC_RINGS[resn]:
            points = [positions[n] for n in ring_names if n in positions]
            if len(points) < 3:
                continue  # incomplete sidechain
            centroid = tuple(sum(p[i] for p in points) / len(points) for i in range(3))
            normal = _plane_normal(points)
            if normal is None:
                continue
            rings.append((key, f"{resn}{key[1]}", centroid, normal))
    return rings


def _plane_normal(points: list[tuple[float, float, float]]) -> tuple[float, float, float] | None:
    """Unit normal of the plane through the first three points."""
    (ax, ay, az), (bx, by, bz), (cx, cy, cz) = points[:3]
    u = (bx - ax, by - ay, bz - az)
    v = (cx - ax, cy - ay, cz - az)
    n = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    length = math.sqrt(sum(c * c for c in n))
    if length < 1e-9:  # collinear
        return None
    return (n[0] / length, n[1] / length, n[2] / length)


def _stacking_interactions(
    atoms_a: list[Atom], atoms_b: list[Atom]
) -> dict[tuple[tuple[str, str], tuple[str, str]], tuple[str, float]]:
    """Pi-stacking between aromatic residues, keyed by residue pair."""
    found: dict[tuple[tuple[str, str], tuple[str, str]], tuple[str, float]] = {}
    rings_b = _ring_systems(atoms_b)
    if not rings_b:
        return found

    for key_a, _label_a, centroid_a, normal_a in _ring_systems(atoms_a):
        for key_b, _label_b, centroid_b, normal_b in rings_b:
            if key_a == key_b:
                continue
            gap = math.sqrt(sum((centroid_a[i] - centroid_b[i]) ** 2 for i in range(3)))
            if gap > AROMATIC_CENTROID_MAX:
                continue

            dot = abs(sum(normal_a[i] * normal_b[i] for i in range(3)))
            angle = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
            if angle <= STACKING_PARALLEL_MAX_ANGLE:
                kind = "pi-stacking (parallel)"
            elif angle >= STACKING_TSHAPED_MIN_ANGLE:
                kind = "pi-stacking (T-shaped)"
            else:
                continue  # geometry too oblique to call

            pair = (key_a, key_b)
            if pair not in found or gap < found[pair][1]:
                found[pair] = (kind, gap)
    return found


def _fetch_atoms(selection: str, include_water: bool, include_hydrogen: bool) -> list[Atom] | None:
    """Pull a selection out of PyMOL as parsed atom records."""
    query = f"({selection})"
    if not include_water:
        query += " and not solvent"
    res = send_request("get_pdbstr", args=[query], timeout=120.0)
    if res.get("status") == "error":
        return None
    atoms = parse_atoms(res.get("result") or "")
    if not include_hydrogen:
        atoms = [a for a in atoms if a.element != "H"]
    return atoms


@mcp.tool()
def contact_report(
    selection1: Annotated[str, Field(description='One side, e.g. "1hsg and resn MK1" (a ligand).')],
    selection2: Annotated[str, Field(description='The other side, e.g. "1hsg and polymer".')],
    cutoff: Annotated[
        float,
        Field(
            description="Maximum heavy-atom separation to count as a contact, in Angstrom. 4.0 captures the interactions above; raise toward 5.0 for a looser survey."
        ),
    ] = 4.0,
    max_pairs: Annotated[
        int,
        Field(
            description="How many residue pairs to list, closest first. The count of any omitted pairs is always reported."
        ),
    ] = 40,
    include_water: Annotated[
        bool,
        Field(
            description="Include waters, for water-mediated contacts. Off by default, as ordered waters otherwise dominate the list."
        ),
    ] = False,
) -> str:
    """
    Lists the residues in contact across two selections, with distances and types.

    This is the numeric counterpart to ligand_view and interface_view: instead
    of drawing the interactions, it reports them — which residue pairs touch,
    how close they get, and whether the contact is a salt bridge, hydrogen
    bond, hydrophobic packing or pi-stacking. Use it to answer "what holds
    this ligand in the pocket" or "which residues form this interface".

    Classification uses heavy-atom distance criteria, since crystal structures
    usually have no hydrogens: salt bridge <= 4.0 A between formally charged
    sidechain tips, hydrogen bond <= 3.5 A between N/O pairs, hydrophobic
    <= 4.5 A between C/S pairs, pi-stacking <= 5.5 A between aromatic ring
    centroids (classified parallel or T-shaped by interplanar angle). A
    reported hydrogen bond is therefore a donor-acceptor pair with plausible
    geometry, not one verified against a hydrogen position.
    """
    if cutoff <= 0:
        return f"Error: cutoff must be positive, got {cutoff}."

    atoms1 = _fetch_atoms(selection1, include_water, include_hydrogen=False)
    if atoms1 is None:
        return f"Error reading '{selection1}' from PyMOL."
    atoms2 = _fetch_atoms(selection2, include_water, include_hydrogen=False)
    if atoms2 is None:
        return f"Error reading '{selection2}' from PyMOL."

    if not atoms1:
        return f"Selection '{selection1}' matched no atoms."
    if not atoms2:
        return f"Selection '{selection2}' matched no atoms."

    pairs = _neighbour_pairs(atoms1, atoms2, cutoff)

    # Aggregate atom-level contacts up to residue pairs, which is the level
    # people actually reason and write about.
    residues: dict[tuple[tuple[str, str], tuple[str, str]], dict] = {}
    for a, b, d in pairs:
        if a.residue_key == b.residue_key and a.chain == b.chain:
            continue  # a residue contacting itself, when the selections overlap
        key = (a.residue_key, b.residue_key)
        entry = residues.setdefault(
            key,
            {"labels": (a.label, b.label), "min": d, "count": 0, "kinds": set()},
        )
        entry["count"] += 1
        entry["min"] = min(entry["min"], d)
        entry["kinds"].add(_classify(a, b, d))

    for pair, (kind, gap) in _stacking_interactions(atoms1, atoms2).items():
        if pair in residues:
            residues[pair]["kinds"].add(kind)
            residues[pair]["min"] = min(residues[pair]["min"], gap)

    if not residues:
        return (
            f"No contacts within {cutoff:.1f} A between '{selection1}' and "
            f"'{selection2}'. They may be far apart, or the selections may "
            f"overlap — check with count_atoms."
        )

    ordered = sorted(
        residues.items(),
        key=lambda item: (item[1]["min"], residue_order(item[0][0][1])),
    )

    lines = [
        f"{len(ordered)} residue pairs in contact within {cutoff:.1f} A "
        f"between '{selection1}' and '{selection2}':",
        "",
    ]
    for _key, entry in ordered[:max_pairs]:
        kinds = sorted(entry["kinds"], key=_kind_rank)
        left, right = entry["labels"]
        lines.append(
            f"  {left:<12} -- {right:<12} {entry['min']:5.2f} A  "
            f"{', '.join(kinds)} ({entry['count']} atom contacts)"
        )

    if len(ordered) > max_pairs:
        lines.append(f"  ... and {len(ordered) - max_pairs} more pairs (raise max_pairs to see).")

    tally: dict[str, int] = {}
    for _key, entry in ordered:
        for kind in entry["kinds"]:
            tally[kind] = tally.get(kind, 0) + 1
    summary = ", ".join(f"{n} {kind}" for kind, n in sorted(tally.items(), key=lambda kv: -kv[1]))
    lines += ["", f"Interaction types across all pairs: {summary}."]
    return "\n".join(lines)


def _kind_rank(kind: str) -> int:
    if kind.startswith("pi-stacking"):
        return _INTERACTION_ORDER.index("hydrogen bond")  # rank alongside specific interactions
    try:
        return _INTERACTION_ORDER.index(kind)
    except ValueError:
        return len(_INTERACTION_ORDER)


# ── Interface burial ─────────────────────────────────────────────────────────
#
# Buried surface area is the standard measure of how much of a complex is
# actually complex. Convention: total buried = SASA(A alone) + SASA(B alone) -
# SASA(AB), and the "interface area" quoted in papers is half of that, i.e.
# the area contributed by one side.

# Interpretation thresholds for a *per-side* interface area, from surveys of
# the PDB (Janin et al.; Krissinel & Henrick). These are guidance, not a
# verdict — small biological interfaces and large crystal contacts both exist.
_CRYSTAL_CONTACT_MAX_AREA = 400.0
_LARGE_INTERFACE_MIN_AREA = 1000.0

# A residue is counted as part of the interface once it buries this much.
_INTERFACE_RESIDUE_MIN_DELTA = 1.0

_HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}
_CHARGED_RESIDUES = {"ASP", "GLU", "LYS", "ARG"}


def _residue_sasa(selection: str) -> dict[tuple[str, str], tuple[str, float]] | None:
    """Per-residue SASA for a selection, as (chain, resi) -> (resn, area).

    ``get_area(load_b=1)`` writes each atom's SASA into the B-factor column,
    so one area calculation plus one PDB dump yields the whole per-residue
    breakdown — rather than one round trip per residue.
    """
    area = send_request(
        "get_area", args=[selection], kwargs={"state": 1, "load_b": 1}, timeout=180.0
    )
    if area.get("status") == "error":
        return None

    dump = send_request("get_pdbstr", args=[selection], timeout=180.0)
    if dump.get("status") == "error":
        return None

    per_residue: dict[tuple[str, str], tuple[str, float]] = {}
    for atom in parse_atoms(dump.get("result") or ""):
        resn, total = per_residue.get(atom.residue_key, (atom.resn, 0.0))
        per_residue[atom.residue_key] = (resn, total + atom.bfactor)
    return per_residue


def _residue_class(resn: str) -> str:
    if resn in _CHARGED_RESIDUES:
        return "charged"
    if resn in _HYDROPHOBIC_RESIDUES:
        return "hydrophobic"
    return "polar"


@mcp.tool()
def interface_report(
    obj_name: Annotated[str, Field(description='PyMOL object holding the complex (e.g. "1brs").')],
    chain_a: Annotated[str, Field(description='First chain ID (e.g. "A").')],
    chain_b: Annotated[str, Field(description='Second chain ID (e.g. "D").')],
    max_residues: Annotated[
        int, Field(description="How many of the most-buried residues to list per chain.")
    ] = 15,
) -> str:
    """
    Measures how large a protein-protein interface is, and which residues form it.

    Reports buried surface area — the standard measure of how much of a
    complex is actually complex — by comparing each chain's solvent-accessible
    area free and bound. Also ranks the residues by how much surface each
    buries, and breaks the interface down by residue chemistry.

    Interpretation: a per-side area under ~400 A^2 is usually a crystal packing
    contact rather than a biological interface, while over ~1000 A^2 indicates
    a substantial, likely specific association. These are guides from PDB-wide
    surveys, not a verdict — small biological interfaces exist.

    For the interactions themselves — which pairs hydrogen bond, which form
    salt bridges — use ``contact_report`` on the same two chains.
    """
    if chain_a == chain_b:
        return f"Error: chain_a and chain_b are both '{chain_a}'; pick two different chains."

    # SASA is meaningless without the solvent-accessible dot mode, and the
    # default density is too coarse for per-residue numbers.
    send_request("set", args=["dot_solvent", "1"])
    send_request("set", args=["dot_density", "3"])

    free_a, free_b, bound = "_iface_free_a", "_iface_free_b", "_iface_bound"
    try:
        for tmp, sel in (
            (free_a, f"({obj_name}) and chain {chain_a} and polymer"),
            (free_b, f"({obj_name}) and chain {chain_b} and polymer"),
            (bound, f"({obj_name}) and chain {chain_a}+{chain_b} and polymer"),
        ):
            send_request("delete", args=[tmp])
            res = send_request("create", args=[tmp, sel])
            if res.get("status") == "error":
                return f"Error isolating '{sel}': {res.get('error')}"

        unbound = {
            chain_a: _residue_sasa(free_a),
            chain_b: _residue_sasa(free_b),
        }
        complexed = {
            chain_a: _residue_sasa(f"{bound} and chain {chain_a}"),
            chain_b: _residue_sasa(f"{bound} and chain {chain_b}"),
        }
    finally:
        for tmp in (free_a, free_b, bound):
            send_request("delete", args=[tmp])

    for chain, values in list(unbound.items()) + list(complexed.items()):
        if values is None:
            return f"Error computing surface area for chain {chain} of '{obj_name}'."
        if not values:
            return (
                f"Chain {chain} of '{obj_name}' has no polymer atoms. "
                f"Check the chain IDs with list_chains."
            )

    buried: dict[str, list[tuple[str, str, float]]] = {}
    totals: dict[str, float] = {}
    for chain in (chain_a, chain_b):
        rows = []
        for key, (resn, free_area) in unbound[chain].items():  # type: ignore[union-attr]
            bound_area = (complexed[chain] or {}).get(key, (resn, free_area))[1]
            delta = free_area - bound_area
            if delta >= _INTERFACE_RESIDUE_MIN_DELTA:
                rows.append((resn, key[1], delta))
        rows.sort(key=lambda r: r[2], reverse=True)
        buried[chain] = rows
        totals[chain] = sum(r[2] for r in rows)

    total_buried = totals[chain_a] + totals[chain_b]
    per_side = total_buried / 2.0

    if per_side < 1.0:
        return (
            f"Chains {chain_a} and {chain_b} of '{obj_name}' bury no measurable "
            f"surface — they are not in contact."
        )

    if per_side < _CRYSTAL_CONTACT_MAX_AREA:
        verdict = "small enough to be a crystal packing contact rather than a biological interface"
    elif per_side >= _LARGE_INTERFACE_MIN_AREA:
        verdict = "a large interface, typical of a specific and possibly obligate complex"
    else:
        verdict = "a typical size for a specific but transient protein-protein interface"

    lines = [
        f"Interface between chains {chain_a} and {chain_b} of {obj_name}:",
        f"  Buried surface area: {per_side:,.0f} A^2 per side ({total_buried:,.0f} A^2 total).",
        f"  That is {verdict}.",
        f"  Interface residues: {len(buried[chain_a])} in chain {chain_a}, "
        f"{len(buried[chain_b])} in chain {chain_b}.",
    ]

    composition: dict[str, float] = {}
    for chain in (chain_a, chain_b):
        for resn, _resi, delta in buried[chain]:
            kind = _residue_class(resn)
            composition[kind] = composition.get(kind, 0.0) + delta
    if composition and total_buried:
        parts = [
            f"{100.0 * v / total_buried:.0f}% {k}"
            for k, v in sorted(composition.items(), key=lambda kv: -kv[1])
        ]
        lines.append(f"  Composition by buried area: {', '.join(parts)}.")

    for chain in (chain_a, chain_b):
        rows = buried[chain][:max_residues]
        if not rows:
            continue
        listed = ", ".join(f"{resn}{resi} ({delta:.0f} A^2)" for resn, resi, delta in rows)
        omitted = len(buried[chain]) - len(rows)
        suffix = f", and {omitted} more" if omitted > 0 else ""
        lines.append(f"  Chain {chain} hot spots: {listed}{suffix}.")

    return "\n".join(lines)
