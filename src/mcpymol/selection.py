"""Safe pieces of PyMOL selection syntax, for values read off atoms or given by a user.

PyMOL parses selections as a language, so a chain or residue identifier pasted
in raw is syntax, not data. The two values that break it occur constantly in
real files: a blank chain, and a negative residue number. ``quote`` handles
both; ``chain_clause`` builds the ``chain ...`` part every module needs.

Born in ``mcpymol.wiggles.atoms``, where the blank-chain case made three views
overwrite every atom in the session; moved here when the same raw
``chain {chain}`` turned up at a dozen sites outside wiggles.
"""

from __future__ import annotations

from collections.abc import Iterable


def quote(value: str) -> str:
    """Quote an identifier read off an atom for use in a selection.

    PyMOL's selection language parses bare identifiers, so interpolating a raw
    ``chain`` or ``resi`` is not safe for the two values that occur constantly
    in real depositions:

    * A **blank chain** — every atom in a file with no chain ID — leaves
      ``chain  and resi 2``, where ``and`` is consumed as the chain name. The
      selection stops being scoped to the object and matches the whole session.
      Quoting fixes this: ``chain ""`` stays scoped.
    * A **negative residue number**, the remnant of an expression tag in most
      NMR and EM entries, is read as a *range*: ``resi -3`` means 1-3, so one
      residue's value gets written across three.

    **Quoting does not fix the second one.** ``resi "-3"`` is still read as the
    range — checked against PyMOL 3.1.0 on an object holding residues -3, 1 and
    2, where it matched all six atoms instead of residue -3's two. The escape
    the grammar actually honours is a **backslash**, and it composes with the
    quoting (``resi "\\-3"`` matches exactly residue -3) and with ``+`` lists
    (``resi \\-3+1`` matches two residues, not four).

    So the value is quoted *and* a leading minus is escaped. The escape is a
    no-op for every other shape — plain numbers, zero, insertion codes like
    ``52A``, blank and lowercase chains were all checked, and a chain literally
    named ``-`` matches under both forms.

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
    # A leading '-' is the range operator to PyMOL's parser even inside quotes.
    escaped = f"\\{value}" if value.startswith("-") else value
    return f'"{escaped}"'


def chain_clause(chains: str | Iterable[str]) -> str:
    """``chain "A"+"B"`` for one chain ID or several, each one quoted.

    A string is split on ``+``, so a caller's ``"A+B"`` still means two
    chains. A blank chain becomes ``chain ""``, which stays scoped to its
    object; unquoted, ``(m) and chain  and polymer.protein`` matched 82 atoms
    where the object held 40, the rest from another object (PyMOL 3.1.0).
    Quoted members compose with ``+`` and ``not``: ``chain ""+"A"`` matched
    exactly those two chains.
    """
    ids = chains.split("+") if isinstance(chains, str) else list(chains)
    return "chain " + "+".join(quote(c.strip()) for c in ids)
