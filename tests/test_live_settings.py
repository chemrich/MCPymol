"""Tools must leave the user's settings and working directory as they found them.

    MCPYMOL_TEST_PORT=9898 pytest -m live tests/test_live_settings.py

WARNING: this clears the PyMOL session it connects to.
"""

import os
import socket

import pytest

import mcpymol.bridge as bridge
from mcpymol.bridge import send_request

pytestmark = pytest.mark.live


@pytest.fixture(autouse=True)
def _clean_session():
    try:
        socket.create_connection((bridge.HOST, bridge.PORT), timeout=2).close()
    except OSError:
        pytest.skip(f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT}")
    send_request("reinitialize")
    yield
    send_request("reinitialize")


def _do(line: str) -> None:
    res = send_request("do", args=[line])
    assert res.get("status") == "success", (line, res)


def _setting(name: str) -> str:
    res = send_request("get", args=[name])
    assert res.get("status") == "success", res
    return str(res["result"])


def test_interface_report_puts_the_dot_settings_back():
    """It set dot_solvent and dot_density globally and left them, which also
    changes how every dot representation in the session is drawn."""
    from mcpymol.analysis import interface_report

    _do("fab AAAAAAAA, two, chain=A")
    _do("fab GGGGGGGG, b_side, chain=B")
    _do("alter b_side, chain='B'")
    _do("translate [0, 4, 0], b_side")
    _do("create two, two or b_side")
    _do("set dot_density, 1")
    before = (_setting("dot_solvent"), _setting("dot_density"))

    interface_report(obj_name="two", chain_a="A", chain_b="B")

    assert (_setting("dot_solvent"), _setting("dot_density")) == before


def test_a_long_ligand_name_is_listed_whole():
    """Read from PDB text columns 18-20, a ligand named ABCDE was listed as ABC."""
    from mcpymol.structures import list_ligands

    _do("pseudoatom lig1, resn=ABCDE, hetatm=1, elem=C")
    _do("pseudoatom lig2, resn=ATP, hetatm=1, elem=C")
    _do("create both, lig1 or lig2")
    assert list_ligands("both") == "Ligands in 'both': ABCDE, ATP"


@pytest.mark.network
def test_fetch_leaves_pymols_working_directory_alone():
    """fetch wrote 1ubq.cif into PyMOL's working directory: wherever it was
    launched, which is how eleven .cif files built up in this repo's root."""
    from mcpymol.structures import FETCH_CACHE, fetch_structure

    cwd = str(send_request("do", args=["pwd"]).get("output", "")).strip().splitlines()[-1]
    before = set(os.listdir(cwd))

    out = fetch_structure(pdb_code="1ubq")

    assert out.startswith("Successfully fetched"), out
    assert set(os.listdir(cwd)) == before, sorted(set(os.listdir(cwd)) - before)
    assert any(name.lower().startswith("1ubq") for name in os.listdir(FETCH_CACHE))
