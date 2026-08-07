"""Parsing PDB records that PyMOL hands back as text.

PyMOL has no "send me a variable" command, so the way to get data *out* of a
session is to ask for a PDB dump and read the fixed-width columns.  Three
callers now do this — B-factors for pLDDT, CA coordinates for superposition,
full atom records for contact analysis — so the column arithmetic lives here
once.

Column positions follow the PDB format specification (1-indexed in the spec,
0-indexed here):

===========  =======  ==========================================
Spec cols    Slice    Field
===========  =======  ==========================================
1-6          0:6      record name (``ATOM``/``HETATM``)
13-16        12:16    atom name
18-20        17:20    residue name
22           21       chain identifier
23-27        22:27    residue sequence number + insertion code
31-38        30:38    x
39-46        38:46    y
47-54        46:54    z
55-60        54:60    occupancy
61-66        60:66    temperature factor (B-factor / pLDDT)
77-78        76:78    element symbol
===========  =======  ==========================================
"""

from typing import NamedTuple


class Atom(NamedTuple):
    """One ATOM/HETATM record."""

    name: str
    resn: str
    chain: str
    resi: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    hetatm: bool

    @property
    def residue_key(self) -> tuple[str, str]:
        """(chain, resi) — the identity a residue is addressed by."""
        return (self.chain, self.resi)

    @property
    def label(self) -> str:
        """Human-readable residue label, e.g. ``A/LYS27``."""
        chain = f"{self.chain}/" if self.chain else ""
        return f"{chain}{self.resn}{self.resi}"


def parse_atoms(pdb_text: str, ca_only: bool = False) -> list[Atom]:
    """Parse ATOM/HETATM records, skipping anything malformed.

    A truncated or non-coordinate line is skipped rather than raised on:
    PyMOL emits headers, CONECT and END records in the same dump, and one bad
    line should not lose the other few thousand good ones.

    Args:
        pdb_text: Raw text from ``get_pdbstr``.
        ca_only: Keep only alpha carbons — one point per residue.
    """
    atoms: list[Atom] = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:  # too short to hold coordinates
            continue

        name = line[12:16].strip()
        if ca_only and name != "CA":
            continue

        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue

        atoms.append(
            Atom(
                name=name,
                resn=line[17:20].strip(),
                chain=line[21].strip(),
                resi=line[22:27].strip(),
                x=x,
                y=y,
                z=z,
                occupancy=_float_or(line[54:60], 1.0),
                bfactor=_float_or(line[60:66], 0.0),
                # The element column is optional in older files; fall back to
                # the leading alphabetic characters of the atom name.
                element=(line[76:78].strip() or _element_from_name(name)).upper(),
                hetatm=line.startswith("HETATM"),
            )
        )
    return atoms


def _float_or(text: str, default: float) -> float:
    try:
        return float(text)
    except ValueError:
        return default


def _element_from_name(name: str) -> str:
    """Guess an element from a PDB atom name.

    Atom names are right-justified in a way that usually puts the element
    first (``CA``, ``NZ``, ``OD1``), but hydrogens may be prefixed with a
    digit (``1HB``), so leading digits are stripped.
    """
    stripped = name.lstrip("0123456789")
    return stripped[:1] if stripped else ""


def residue_order(resi: str) -> tuple[int, str]:
    """Sort key placing ``10A`` after ``10`` and before ``11``.

    Insertion codes are real in antibody and protease numbering, and plain
    ``int(resi)`` raises on them.
    """
    digits = ""
    for ch in resi:
        if ch.isdigit() or (ch == "-" and not digits):
            digits += ch
        else:
            break
    try:
        return (int(digits), resi)
    except ValueError:
        return (0, resi)
