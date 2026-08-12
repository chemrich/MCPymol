"""Reading per-atom properties out of PyMOL, and the residue key they group by.

Shared by every tier-1 view. Kept separate from the views so that the one
place which knows how an atom is fetched is not also the place that decides
how it is coloured.
"""

from __future__ import annotations

from dataclasses import dataclass

from mcpymol.wiggles.port import ITERATE_TO_LIST, PortError, PymolPort

# The iterate expression tier-1 views read. Order matters: it defines the
# tuple layout that _to_atom unpacks.
ATOM_EXPR = "chain, resi, resn, name, alt, q, b"


@dataclass(frozen=True)
class Atom:
    """One atom's identity and the two scalars tier 1 cares about."""

    chain: str
    resi: str
    resn: str
    name: str
    alt: str
    q: float  # occupancy — SENSE 1, per-atom. Never compositional.
    b: float

    @property
    def residue(self) -> tuple[str, str]:
        """Chain + residue number. The grouping key for per-residue views."""
        return (self.chain, self.resi)


def quote(value: str) -> str:
    """Quote an identifier read off an atom for use in a selection.

    PyMOL's selection language parses bare identifiers, so interpolating a raw
    ``chain`` or ``resi`` is not safe for the two values that occur constantly
    in real depositions:

    * A **blank chain** — every atom in a file with no chain ID — leaves
      ``chain  and resi 2``, where ``and`` is consumed as the chain name. The
      selection stops being scoped to the object and matches the whole session.
    * A **negative residue number**, the remnant of an expression tag in most
      NMR and EM entries, is read as a *range*: ``resi -3`` means 1-3, so one
      residue's value gets written across three.

    Both were verified against PyMOL 3.1.0, in both directions. Quoting fixes
    both and is a no-op for ordinary values, so every call site quotes rather
    than testing whether it needs to.

    Raises:
        ValueError: ``value`` contains a double quote, which would close the
            quoting and hand the rest of the string to the parser. Nothing in
            the PDB or mmCIF grammar puts one in a chain or residue
            identifier, so this is a corrupt file rather than a case to
            support — and guessing at it is how the blank-chain bug behaved.
    """
    if '"' in value:
        raise ValueError(
            f"identifier {value!r} contains a double quote, which cannot appear "
            f"in a chain or residue identifier — the file is probably corrupt"
        )
    return f'"{value}"'


def residue_selection(obj: str, chain: str, resi: str) -> str:
    """Selection matching exactly one residue of ``obj``.

    The single place that knows how to name a residue safely; see :func:`quote`
    for what goes wrong without it.
    """
    return f"({obj}) and chain {quote(chain)} and resi {quote(resi)}"


def residue_clause(chain: str, resi: str) -> str:
    """Like :func:`residue_selection` but unscoped, for OR-ing into a group."""
    return f"(chain {quote(chain)} and resi {quote(resi)})"


def _to_atom(row: object) -> Atom:
    if not isinstance(row, (list, tuple)) or len(row) != 7:
        raise PortError(
            f"malformed atom row {row!r}: expected 7 fields ({ATOM_EXPR}), "
            f"got {len(row) if isinstance(row, (list, tuple)) else type(row).__name__}"
        )
    try:
        return Atom(
            chain=str(row[0]),
            resi=str(row[1]),
            resn=str(row[2]),
            name=str(row[3]),
            alt=str(row[4]),
            q=float(row[5]),
            b=float(row[6]),
        )
    except (TypeError, ValueError) as exc:
        raise PortError(f"malformed atom row {row!r}: {exc}") from exc


def fetch_atoms(port: PymolPort, selection: str) -> list[Atom]:
    """Read every atom in ``selection``.

    Raises:
        PortError: the selection is empty, or the port returned something
            that is not a list of atom rows.

    An empty selection is an error rather than an empty list on purpose. Every
    caller here is about to render a view, and a view of nothing is worse than
    a message saying the selection matched nothing — see MCPymol issue #15,
    where an empty result was reported as success.
    """
    raw = port.query(ITERATE_TO_LIST, selection, ATOM_EXPR)
    if raw is None:
        raise PortError(f"selection {selection!r} returned nothing")
    if not isinstance(raw, (list, tuple)):
        raise PortError(f"expected a list of atom rows, got {type(raw).__name__}")
    if not raw:
        raise PortError(f"selection {selection!r} matched no atoms")
    return [_to_atom(row) for row in raw]


def count_states(port: PymolPort, obj: str) -> int:
    """How many states (models) the object has."""
    raw = port.query("count_states", obj)
    try:
        return int(raw)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise PortError(f"count_states({obj!r}) returned {raw!r}") from exc


def fetch_state_coords(
    port: PymolPort, selection: str, state: int
) -> list[tuple[float, float, float]]:
    """Coordinates of ``selection`` in one state, in atom order."""
    raw = port.query("get_coords", selection, state)
    if not isinstance(raw, (list, tuple)) or not raw:
        raise PortError(f"no coordinates for {selection!r} in state {state}")
    out: list[tuple[float, float, float]] = []
    for point in raw:
        try:
            x, y, z = point  # type: ignore[misc]
        except (TypeError, ValueError) as exc:
            raise PortError(f"malformed coordinate {point!r}") from exc
        out.append((float(x), float(y), float(z)))
    return out


def group_by_residue(atoms: list[Atom]) -> dict[tuple[str, str], list[Atom]]:
    """Group atoms by (chain, resi), preserving first-seen order."""
    out: dict[tuple[str, str], list[Atom]] = {}
    for atom in atoms:
        out.setdefault(atom.residue, []).append(atom)
    return out


def altloc_groups(atoms: list[Atom]) -> list[str]:
    """Distinct non-blank altloc identifiers, sorted.

    A blank altloc means "no alternate here" and is deliberately excluded —
    it is not a group, it is the absence of grouping.
    """
    return sorted({a.alt for a in atoms if a.alt.strip()})
