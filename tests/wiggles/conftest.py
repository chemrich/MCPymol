"""Shared test fixtures."""

from __future__ import annotations

import socket

import pytest

import mcpymol.bridge as bridge
from mcpymol.wiggles.bfactors import clear_stash
from mcpymol.wiggles.provenance import forget


def _bridge_reachable() -> bool:
    """Is a PyMOL plugin listening? Mirrors ``tests/test_live.py``."""
    try:
        with socket.create_connection((bridge.HOST, bridge.PORT), timeout=2):
            return True
    except OSError:
        return False


@pytest.fixture(autouse=True)
def _require_live_bridge(request):
    """Skip ``live`` tests when no PyMOL is listening, rather than erroring.

    ``tests/test_live.py`` has had this guard all along, but as a module-level
    fixture it does not reach ``tests/wiggles/`` — so the first live test added
    here turned `pytest -m live` from "171 skipped" into six errors for anyone
    without a viewer running. CI is unaffected because it deselects ``live``,
    which is exactly why nobody would notice: it breaks only the workflow
    CONTRIBUTING documents.

    Lives in conftest rather than in the one file that needed it, because the
    defect was that the guard was module-local. Fixing only the file would
    recreate the gap the next time a live test is added here.
    """
    if request.node.get_closest_marker("live") and not _bridge_reachable():
        pytest.skip(
            f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT} — start PyMOL and "
            f"`run src/mcpymol/plugin.py`"
        )


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
