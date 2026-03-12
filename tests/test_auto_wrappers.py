import pytest
import json
import socket
from unittest.mock import patch, MagicMock

# Import all generic wrapper functions
from mcpymol.server import (
    hide, remove, distance, execute_pymol_command,
    as_tool, set_setting, cartoon, spectrum, label,
    angle, dihedral, center, orient, zoom, reset, turn, move, clip,
    save, png, deselect, create, extract, delete, align, super_tool,
    intra_fit, intra_rms, alter, alter_state, h_add, h_fill, bond, unbond,
    rebuild, refresh, util_cbc, util_cbaw, util_cbag, util_cbac, util_cbam,
    util_cbay, util_cbas, util_cbab, util_cbao, util_cbap, util_cbak,
    util_chainbow, util_rainbow, util_ss, util_color_by_element, util_color_secondary,
    spheroid, isomesh, isosurface, sculpt_activate, sculpt_deactivate, sculpt_iterate,
    scene, scene_order, mset, mplay, mstop, frame, forward, backward, rock,
    ray, draw, mpng, symexp, set_symmetry, fab, fragment, full_screen, viewport,
    cd, pwd, ls, system, help
)

@pytest.fixture
def mock_socket():
    with patch("socket.socket") as mock_cls:
        inst = MagicMock()
        mock_cls.return_value.__enter__.return_value = inst
        inst.recv.return_value = json.dumps(
            {"status": "success", "result": "Mocked execution"}
        ).encode()
        yield inst


def _check_send_request(mock_socket, expected_action, expected_args):
    """Helper to verify socket payload."""
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == expected_action
    
    # FastMCP tools often drop optional None arguments or pass them as str. 
    # Check that whatever args were passed match expectation.
    assert payload["args"] == expected_args

# Define (function, tuple of args, expected_action_name, expected_args_list)
WRAPPER_TESTS = [
    (hide, ("sticks", "1abc"), "hide", ["sticks", "1abc"]),
    (remove, ("solvent",), "remove", ["solvent"]),
    (distance, ("dist1", "A", "B"), "distance", ["dist1", "A", "B"]),
    (execute_pymol_command, ("bg_color white",), "do", ["bg_color white"]),
    (as_tool, ("cartoon", "all"), "as", ["cartoon", "all"]),
    (set_setting, ("transparency", "0.5", "1abc"), "set", ["transparency", "0.5", "1abc"]),
    (cartoon, ("putty", "all"), "cartoon", ["putty", "all"]),
    (spectrum, ("b", "rainbow", "all"), "spectrum", ["b", "rainbow", "all"]),
    (label, ("CA", "name"), "label", ["CA", "name"]),
    (angle, ("ang1", "A", "B", "C"), "angle", ["ang1", "A", "B", "C"]),
    (dihedral, ("dih1", "A", "B", "C", "D"), "dihedral", ["dih1", "A", "B", "C", "D"]),
    (center, ("1abc",), "center", ["1abc"]),
    (orient, ("1abc",), "orient", ["1abc"]),
    (zoom, ("1abc", "5"), "zoom", ["1abc", "5"]),
    (reset, ("1abc",), "reset", ["1abc"]),
    (turn, ("x", "90"), "turn", ["x", "90"]),
    (move, ("x", "10"), "move", ["x", "10"]),
    (clip, ("near", "10"), "clip", ["near", "10"]),
    (save, ("test.pdb", "all", "-1"), "save", ["test.pdb", "all", "-1"]),
    (png, ("test.png", "100"), "png", ["test.png", "100"]),
    (deselect, (), "deselect", []),
    (create, ("new_obj", "old_obj", "1"), "create", ["new_obj", "old_obj", "1"]),
    (extract, ("new_obj", "old_obj"), "extract", ["new_obj", "old_obj"]),
    (delete, ("old_obj",), "delete", ["old_obj"]),
    (align, ("mobile", "target", "10"), "align", ["mobile", "target", "10"]),
    (super_tool, ("mobile", "target", "10"), "super", ["mobile", "target", "10"]),
    (intra_fit, ("1abc",), "intra_fit", ["1abc"]),
    (intra_rms, ("1abc",), "intra_rms", ["1abc"]),
    (alter, ("1abc", "b=10"), "alter", ["1abc", "b=10"]),
    (alter_state, ("1", "1abc", "b=10"), "alter_state", ["1", "1abc", "b=10"]),
    (h_add, ("all",), "h_add", ["all"]),
    (h_fill, ("all",), "h_fill", ["all"]),
    (bond, ("A", "B", "2"), "bond", ["A", "B", "2"]),
    (unbond, ("A", "B"), "unbond", ["A", "B"]),
    (rebuild, ("all",), "rebuild", ["all"]),
    (refresh, (), "refresh", []),
    (util_cbc, ("all",), "util.cbc", ["all"]),
    (util_cbaw, ("all",), "util.cbaw", ["all"]),
    (util_chainbow, ("all",), "util.chainbow", ["all"]),
    (util_rainbow, ("all",), "util.rainbow", ["all"]),
    (util_ss, ("all",), "util.ss", ["all"]),
    (util_color_by_element, ("all",), "util.color_by_element", ["all"]),
    (util_color_secondary, ("all",), "util.color_secondary", ["all"]),
    (spheroid, ("all",), "spheroid", ["all"]),
    (isomesh, ("mesh1", "map1", "1.0", "all"), "isomesh", ["mesh1", "map1", "1.0", "all"]),
    (isosurface, ("surf1", "map1", "1.0", "all"), "isosurface", ["surf1", "map1", "1.0", "all"]),
    (sculpt_activate, ("1abc",), "sculpt_activate", ["1abc"]),
    (sculpt_deactivate, ("1abc",), "sculpt_deactivate", ["1abc"]),
    (sculpt_iterate, ("10", "1abc"), "sculpt_iterate", ["10", "1abc"]),
    (scene, ("F1", "store"), "scene", ["F1", "store"]),
    (scene_order, ("F1 F2 F3",), "scene_order", ["F1 F2 F3"]),
    (mset, ("1 x100",), "mset", ["1 x100"]),
    (mplay, (), "mplay", []),
    (mstop, (), "mstop", []),
    (frame, ("10",), "frame", ["10"]),
    (forward, (), "forward", []),
    (backward, (), "backward", []),
    (rock, (), "rock", []),
    (ray, ("800", "600"), "ray", ["800", "600"]),
    (draw, ("800", "600"), "draw", ["800", "600"]),
    (mpng, ("frame_",), "mpng", ["frame_"]),
    (symexp, ("sym_", "1abc", "20", "A"), "symexp", ["sym_", "1abc", "20", "A"]),
    (set_symmetry, ("1abc", "10", "10", "10", "90", "90", "90"), "set_symmetry", ["1abc", "10", "10", "10", "90", "90", "90"]),
    (fab, ("ACDEFGH", "1"), "fab", ["ACDEFGH", "1"]),
    (fragment, ("benzene",), "fragment", ["benzene"]),
    (full_screen, (), "full_screen", []),
    (viewport, ("800", "600"), "viewport", ["800", "600"]),
    (cd, ("/tmp",), "cd", ["/tmp"]),
    (pwd, (), "pwd", []),
    (ls, ("/tmp",), "ls", ["/tmp"]),
    (system, ("echo hello",), "system", ["echo hello"]),
    (help, ("show",), "help", ["show"])
]

@pytest.mark.parametrize("func,args,expected_action,expected_args", WRAPPER_TESTS)
def test_auto_wrappers(mock_socket, func, args, expected_action, expected_args):
    """Test that all boilerplate generic wrappers correctly pass arguments to send_request."""
    result = func(*args)
    assert "Executed" in result or "successfully" in result or result.startswith(expected_action) or "Mocked execution" in result or "Executed command" in result or "Hiding" in result or "Removed" in result or "Measured" in result
    _check_send_request(mock_socket, expected_action, expected_args)

