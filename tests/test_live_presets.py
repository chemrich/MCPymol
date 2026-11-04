"""Every view preset against a real PyMOL: nothing it sends is rejected.

    pytest -m live tests/test_live_presets.py

WARNING: this clears the PyMOL session it connects to. It fetches 1HSG (HIV
protease with its inhibitor MK1) and skips when offline.

Presets send dozens of commands and discard the answers, so a call PyMOL
rejects — an unknown setting, a selection keyword that does not exist — never
reaches the tool's result: the preset reports success over a step that did
nothing. Four such calls sat in the presets until a request log showed them
(KNOWN-ISSUES §4). This records every response while each preset runs.
"""

import socket
import sys

import pytest

import mcpymol.bridge as bridge
import mcpymol.server  # noqa: F401  (registers every tool)
from mcpymol.bridge import send_request
from mcpymol.structures import load_structure
from mcpymol.views import (
    bfactor_view,
    cinematic_view,
    crosslink_view,
    hydrophobic_surface_view,
    interface_view,
    ligand_view,
    mutation_view,
    pharmacophore_view,
    pocket_view,
    pointillist_view,
    putty_view,
    textbook_view,
)

pytestmark = pytest.mark.live

OBJ = "1hsg"

# (preset, the arguments it needs beyond obj_name). Left out: the presets that
# shell out to APBS or call a remote MSA service, and plddt_view, which wants
# an AlphaFold model.
PRESETS = [
    (ligand_view, {"ligand_resn": "MK1"}),
    (bfactor_view, {}),
    (interface_view, {"chain_a": "A", "chain_b": "B"}),
    (putty_view, {}),
    (hydrophobic_surface_view, {}),
    (crosslink_view, {}),
    (pocket_view, {"resn": "MK1"}),
    (pharmacophore_view, {"resn": "MK1"}),
    (mutation_view, {"mutations": "L10G,I50V"}),
    (textbook_view, {}),
    (cinematic_view, {}),
    (pointillist_view, {}),
]


def _bridge_reachable() -> bool:
    try:
        with socket.create_connection((bridge.HOST, bridge.PORT), timeout=2):
            return True
    except OSError:
        return False


@pytest.fixture
def structure():
    if not _bridge_reachable():
        pytest.skip(f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT}")
    send_request("reinitialize")
    if send_request("fetch", args=[OBJ, OBJ], timeout=60).get("status") != "success":
        pytest.skip("could not fetch 1HSG (offline?)")
    yield OBJ
    send_request("reinitialize")


@pytest.fixture
def rejected(monkeypatch):
    """Every response, during the test, that PyMOL rejected or warned about."""
    real = bridge.send_request
    found: list[str] = []

    def recording(action, args=None, kwargs=None, **kw):
        res = real(action, args=args, kwargs=kwargs, **kw)
        warned = [
            w
            for w in res.get("warnings") or []
            if any(k in w for k in ("Error", "Setting-Warning"))
        ]
        if res.get("status") == "error" or warned:
            found.append(f"{action} {args}: {res.get('error') or '; '.join(warned)}")
        return res

    for name, module in list(sys.modules.items()):
        if name.startswith("mcpymol") and getattr(module, "send_request", None) is real:
            monkeypatch.setattr(module, "send_request", recording)
    return found


@pytest.mark.parametrize("preset, extra", PRESETS, ids=[p.__name__ for p, _ in PRESETS])
def test_preset_sends_nothing_pymol_rejects(structure, rejected, preset, extra):
    result = preset(obj_name=structure, **extra)

    assert not result.startswith("Error"), result
    assert rejected == [], f"{preset.__name__} sent calls PyMOL rejected:\n" + "\n".join(rejected)


@pytest.mark.parametrize("preset, extra", PRESETS, ids=[p.__name__ for p, _ in PRESETS])
def test_preset_refuses_an_object_that_is_not_there(structure, preset, extra):
    """A typo used to style an empty selection and report success."""
    result = preset(obj_name="nosuch_obj", **extra)

    # PyMOL answers an unknown name as an invalid selection, an empty one as
    # zero atoms; either way the preset must refuse before styling anything.
    assert result.startswith("Error: object 'nosuch_obj' "), result
    assert "is not a valid selection" in result or "matches no atoms" in result, result


def test_ligand_view_refuses_a_ligand_that_is_not_there(structure):
    assert ligand_view(obj_name=structure, ligand_resn="ZZZ") == (
        "Error: ligand 'ZZZ' in 1hsg matches no atoms."
    )


def test_loading_nucleic_acid_sets_its_colours_where_pymol_reads_them(rejected, tmp_path):
    """cartoon_nucleic_acid_color given an atom selection was stored per atom
    and never drawn; PyMOL only warned. It is set on the object now."""
    if not _bridge_reachable():
        pytest.skip(f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT}")
    send_request("reinitialize")
    send_request("do", args=["fnab ACGTACGT, name=dna_src, mode=DNA, form=B, dbl_helix=1"])
    path = tmp_path / "dna.pdb"
    assert send_request("save", args=[str(path), "dna_src"]).get("status") == "success"
    send_request("delete", args=["dna_src"])

    result = load_structure(file_path=str(path), obj_name="dna")

    assert "Error" not in result, result
    assert rejected == [], "\n".join(rejected)
    assert (
        send_request("get", args=["cartoon_nucleic_acid_color", "dna"])["result"] == "brightorange"
    )
    send_request("reinitialize")
