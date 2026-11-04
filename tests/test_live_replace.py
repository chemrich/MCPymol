"""replace=True against a real PyMOL: "the only structure loaded" means every object.

    pytest -m live tests/test_live_replace.py

WARNING: this clears the PyMOL session it connects to.

The loaders clear the rest of the session by name once the new structure is
known to have loaded. They listed names with get_object_list, which is
selector-backed and returns molecules only, so a previous structure's map,
ramp, measurement and group survived a "replace". The mocked tests can only
say which call is sent; this checks what PyMOL actually lists afterwards.
"""

import socket

import pytest

import mcpymol.bridge as bridge
from mcpymol.bridge import send_request
from mcpymol.structures import load_structure

pytestmark = pytest.mark.live


def _bridge_reachable() -> bool:
    try:
        with socket.create_connection((bridge.HOST, bridge.PORT), timeout=2):
            return True
    except OSError:
        return False


@pytest.fixture(autouse=True)
def _clean_session():
    if not _bridge_reachable():
        pytest.skip(f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT}")
    send_request("reinitialize")
    yield
    send_request("reinitialize")


def _objects() -> list[str]:
    res = send_request("get_names", args=["objects"])
    assert res.get("status") == "success", res
    return list(res["result"])


def _do(line: str) -> None:
    res = send_request("do", args=[line])
    assert res.get("status") == "success", (line, res)


def test_replace_leaves_only_the_new_structure_whatever_was_loaded(tmp_path):
    new_file = tmp_path / "newcomer.pdb"
    _do("fab GGGG, newcomer_src")
    assert send_request("save", args=[str(new_file), "newcomer_src"]).get("status") == "success"
    _do("delete newcomer_src")

    # The previous structure, and everything a view might have hung off it.
    _do("fab AAAA, oldmol")
    _do("map_new oldmap, gaussian, 1.0, oldmol")
    _do("ramp_new oldramp, oldmap, [0, 1]")
    _do("distance old_d, oldmol and resi 1 and name CA, oldmol and resi 3 and name CA")
    _do("pseudoatom oldps")
    _do("group oldgroup, oldps")
    before = set(_objects())
    assert {"oldmol", "oldmap", "oldramp", "old_d", "oldgroup"} <= before, before

    result = load_structure(file_path=str(new_file), obj_name="newcomer", replace=True)

    assert "Error" not in result, result
    assert _objects() == ["newcomer"], _objects()


def test_add_mode_keeps_everything(tmp_path):
    new_file = tmp_path / "newcomer.pdb"
    _do("fab GGGG, newcomer_src")
    assert send_request("save", args=[str(new_file), "newcomer_src"]).get("status") == "success"
    _do("delete newcomer_src")
    _do("fab AAAA, oldmol")
    _do("map_new oldmap, gaussian, 1.0, oldmol")

    load_structure(file_path=str(new_file), obj_name="newcomer", replace=False)

    assert {"oldmol", "oldmap", "newcomer"} <= set(_objects())


def test_list_objects_names_every_kind_of_object():
    """get_object_list lists molecules only: a session holding just a map
    read "No objects are loaded", and a structure's map, measurement and
    group went unlisted (PyMOL 3.1.0)."""
    from mcpymol.structures import list_objects

    _do("fab AAAA, mol1")
    _do("map_new themap, gaussian, 1.0, mol1")
    _do("distance d1, mol1 and name CA and resi 1, mol1 and name CA and resi 2")
    _do("group grp, mol1")
    out = list_objects()
    assert out.startswith("Loaded objects: "), out
    listed = out.removeprefix("Loaded objects: ").split(", ")
    assert sorted(listed) == sorted(["mol1", "themap (map)", "d1 (measurement)", "grp (group)"]), (
        out
    )

    _do("delete mol1 or d1 or grp")
    assert list_objects() == "Loaded objects: themap (map)"


# ── A name must never delete what it does not name ─────────────────────────────


def _bystanders() -> tuple[list[str], int]:
    """Two objects the call under test must leave alone, one of them recoloured."""
    _do("fab AAAAA, keepme1")
    _do("fab GGGGG, keepme2")
    _do("color red, keepme2")
    return _objects(), _red_atoms()


def _red_atoms() -> int:
    res = send_request("count_atoms", args=["keepme2 and color red"])
    return int(res.get("result") or 0)


@pytest.mark.parametrize("name", ["*", "all", "keep*", "pep two", "a;b"])
def test_a_name_that_is_a_pattern_is_refused_and_deletes_nothing(tmp_path, name):
    """`delete` takes patterns and lists. Given `*`, `all` or `keep*`,
    load_structure deleted every other object and reported success, even
    with replace=False; `pep two` loaded as `pep_two` and then said nothing
    had loaded. Seen on PyMOL 3.1.0 before the fix."""
    source = tmp_path / "pep.pdb"
    _do("fab CCC, src")
    assert send_request("save", args=[str(source), "src"]).get("status") == "success"
    _do("delete src")
    before = _bystanders()

    out = load_structure(file_path=str(source), obj_name=name, replace=False)

    assert out.startswith("Error:") and "can't be used as an object name" in out, out
    assert (_objects(), _red_atoms()) == before


def test_a_preset_given_a_selection_does_not_delete_its_object():
    """crosslink_view names its distance objects after obj_name and deletes the
    old ones first. Given a selection, `delete keep1 and chain A_disulfides`
    deleted the object keep1 itself (PyMOL 3.1.0)."""
    from mcpymol.views import crosslink_view

    _do("fab CACAC, keep1, chain=A")
    _do("fab GGG, keep2")

    out = crosslink_view(obj_name="keep1 and chain A")

    assert not out.startswith("Error"), out
    names = _objects()
    assert {"keep1", "keep2"} <= set(names), names
    assert "keep1_and_chain_A_disulfides" in names, names
