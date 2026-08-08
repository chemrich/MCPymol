"""Saving and restoring B-factors around views that overwrite them.

Several views colour by pushing a scalar into the B-factor column, because
that is what PyMOL's ``spectrum`` reads. It is destructive, so anything doing
it owes the user a way back.

**The values are held here, in Python, not in PyMOL.** They come free: every
view already calls :func:`mcpymol.wiggles.atoms.fetch_atoms`, and that reads ``b``.
Keeping them on this side means a restore does not depend on some PyMOL field
surviving whatever the user does to the session in between.

Restoring pushes them back in **one** command via PyMOL's ``stored`` namespace
and a single ``alter``, rather than one ``alter`` per atom — a 10,000-atom
structure would otherwise mean 10,000 round trips. This mirrors what MCPymol's
``conservation_view`` already does in production, which is the reason to
prefer it over the ``custom`` field: ``custom`` is a *string* property, so a
float pushed into it would stringify, and a view claiming "B-factors preserved"
on that basis would be asserting something it had not established.
"""

from __future__ import annotations

import json

from mcpymol.wiggles.atoms import Atom
from mcpymol.wiggles.port import PortError, PymolPort, call

# obj name -> {(chain, resi, name, alt): b}
# Process-global. Correct for a single-session MCP server, which is what
# MCPymol is today; wrong the moment one process serves two sessions. Keyed
# by object name, so two sessions with the same object name would collide.
# Flagged rather than solved — see MOVING.md.
_STASH: dict[str, dict[tuple[str, str, str, str], float]] = {}


def _key(atom: Atom) -> tuple[str, str, str, str]:
    """Atom identity for restore purposes.

    Chain + residue number + atom name + altloc. This does not include the
    PDB insertion code, which PyMOL exposes separately — two atoms differing
    only by insertion code would collide. Rare, and the failure is a restored
    B-factor being wrong rather than lost, but it is a real limit.
    """
    return (atom.chain, atom.resi, atom.name, atom.alt)


def stash_bfactors(obj: str, atoms: list[Atom]) -> int:
    """Record the B-factors of ``atoms`` under ``obj``. Returns how many.

    Takes atoms already fetched rather than querying again — the caller has
    them, and a second read would be a second round trip for data we hold.
    """
    _STASH[obj] = {_key(a): a.b for a in atoms}
    return len(_STASH[obj])


def has_stash(obj: str) -> bool:
    """Is there anything saved for ``obj``?"""
    return obj in _STASH


def clear_stash(obj: str | None = None) -> None:
    """Forget the stash for ``obj``, or all of them."""
    if obj is None:
        _STASH.clear()
    else:
        _STASH.pop(obj, None)


def restore_bfactors(port: PymolPort, obj: str) -> str:
    """Put ``obj``'s original B-factors back.

    Raises:
        PortError: nothing was stashed for ``obj`` — restoring from an empty
            stash would silently zero the column, which is worse than failing.
    """
    saved = _STASH.get(obj)
    if not saved:
        raise PortError(
            f"no saved B-factors for {obj!r}. Either no view has overwritten "
            f"them in this session, or the stash was cleared."
        )

    # JSON round-trip so the dict survives as a PyMOL command-line literal.
    # Tuple keys are not JSON-able, so key on a joined string and rebuild the
    # same string in the alter expression.
    flat = {"|".join(k): v for k, v in saved.items()}
    port.do(f"stored.wiggles_b = {json.dumps(flat)}")
    call(port, "alter", obj, "b=stored.wiggles_b.get('|'.join((chain, resi, name, alt)), b)")
    call(port, "rebuild", obj)
    return (
        f"Restored {len(saved)} B-factors on {obj}. Atoms not present when the "
        f"stash was taken keep their current value."
    )


def preservation_note(obj: str, stashed: int) -> str:
    """The line a view prints about what it did to the B-factor column.

    Deliberately states where the values are and how to get them back, rather
    than the bare word "preserved" — a claim the user cannot act on is not
    worth making.
    """
    return (
        f"  B-factor column overwritten. The original {stashed} values are held "
        f"in this session;\n  call restore_bfactors(port, '{obj}') to put them "
        f"back."
    )
