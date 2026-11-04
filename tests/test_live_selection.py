"""A blank chain, against a real PyMOL: selections stay inside their object.

    pytest -m live tests/test_live_selection.py

WARNING: this clears the PyMOL session it connects to.

get_chains returns [''] for a file with no chain IDs, and a raw
``chain {chain}`` then reads ``chain  and polymer.protein``: PyMOL takes
``and`` as the chain name and the selection escapes the object (82 atoms
matched where the object held 40, PyMOL 3.1.0). Loading such a file colours
each chain, so with replace=False it recoloured every other object's protein.
"""

import socket

import pytest

import mcpymol.bridge as bridge
from mcpymol.bridge import send_request
from mcpymol.structures import load_structure

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


def _colours(selection: str) -> set:
    res = send_request("iterate_to_list", args=[selection, "color"])
    assert res.get("status") == "success", res
    return {row[0] for row in res["result"]}


def test_loading_a_blank_chain_file_leaves_other_objects_alone(tmp_path):
    send_request("fab", args=["AAAA", "blank_src"])
    send_request("alter", args=["blank_src", "chain=''"])
    path = tmp_path / "blank.pdb"
    assert send_request("save", args=[str(path), "blank_src"]).get("status") == "success"
    send_request("delete", args=["blank_src"])
    assert send_request("get_chains", args=["all"]).get("result") in ([], None)

    send_request("fab", args=["GGGGGG", "bystander"])
    send_request("alter", args=["bystander", "chain='B'"])
    send_request("color", args=["red", "bystander"])
    red = _colours("bystander")

    result = load_structure(file_path=str(path), obj_name="newcomer", replace=False)

    assert "Error" not in result, result
    assert send_request("get_chains", args=["newcomer"])["result"] == [""]
    assert _colours("bystander") == red, "loading the blank-chain file recoloured another object"
    assert send_request("count_atoms", args=["bystander"])["result"] > 0
