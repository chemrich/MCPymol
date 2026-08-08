"""Shared test fixtures."""

from __future__ import annotations

import pytest

from mcpymol.wiggles.bfactors import clear_stash
from mcpymol.wiggles.provenance import forget


@pytest.fixture(autouse=True)
def _isolate_module_state():
    """Both the B-factor stash and the provenance registry are module-level.

    Without this, a view in one test leaves state that changes what another
    test sees — the sort of coupling that makes a suite pass in one order and
    fail in another. Provenance especially: a leaked declaration would make
    "UNKNOWN by default" pass for the wrong reason.
    """
    clear_stash()
    forget()
    yield
    clear_stash()
    forget()
