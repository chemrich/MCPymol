"""Shared test setup.

``mcpymol.plugin`` does ``from pymol import cmd`` at import time and binds the
result as a module global, so the mock has to be installed before *any* test
module imports it — and it has to be the *same* mock for every test file, or a
test asserting on its own mock silently watches an object the plugin never
uses.  conftest is imported before the test modules, which makes this the one
place that can guarantee both.
"""

import os
import socket
import sys
from unittest.mock import MagicMock

import pytest

mock_pymol = MagicMock()
sys.modules.setdefault("pymol", mock_pymol)
sys.modules.setdefault("pymol.cmd", mock_pymol.cmd)


def _pymol_ports() -> set[int]:
    """The ports a real PyMOL plugin may be listening on, on this machine."""
    ports = {9876}
    try:
        ports.add(int(os.environ.get("MCPYMOL_PORT", 9876)))
    except ValueError:
        pass
    return ports


@pytest.fixture(autouse=True)
def _keep_unit_tests_off_a_real_pymol(request, monkeypatch):
    """Refuse every connection a unit test makes to a real PyMOL's port.

    Unit tests talk to fakes, but nothing stopped one that slipped past its
    mock from reaching whatever PyMOL was listening, and the default port is
    the session the developer has open. That happened while this guard was
    being written: a mutation-testing harness ran old function bodies outside
    the module the tests patch, and their requests ran in a developer's
    session — ``hide everything``, ``show sticks``, a 22 s ray trace written
    into the checkout. Here the request fails as if no PyMOL were running;
    the test itself fails only if the code checks the result. (Some view tests
    leave mcpymol.style.set_background unpatched, and it ignores its result:
    this guard is what keeps those off a real session.) Tests marked ``live``
    are exempt: reaching PyMOL is their whole point.
    """
    if request.node.get_closest_marker("live"):
        return
    blocked = _pymol_ports()
    real_connect = socket.socket.connect

    def connect(sock, address):
        if isinstance(address, tuple) and len(address) >= 2 and address[1] in blocked:
            raise ConnectionRefusedError(
                f"unit test tried to reach a real PyMOL at {address[0]}:{address[1]}"
            )
        return real_connect(sock, address)

    monkeypatch.setattr(socket.socket, "connect", connect)
