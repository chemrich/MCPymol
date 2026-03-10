import pytest
import json
import socket
from unittest.mock import patch, MagicMock
from mcpymol.server import (
    send_request, fetch_structure, load_structure,
    show, color, select, distance,
    ligand_view, interface_view, putty_view,
    hydrophobic_surface_view, electrostatic_view,
    crosslink_view, pocket_view, pharmacophore_view,
    mutation_view, textbook_view, cinematic_view, pointillist_view
)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _sr_mock(**action_results):
    """Return a side_effect for patching send_request.

    Calls whose action key is in action_results get that result value; all
    other calls return a generic success response.  If a value is callable
    it is invoked with (action, args, kwargs) to produce the result.
    """
    def fake(action, args=None, kwargs=None):
        if action in action_results:
            val = action_results[action]
            result = val(action, args, kwargs) if callable(val) else val
            return {"status": "success", "result": result}
        return {"status": "success", "result": "OK"}
    return fake


def _actions(mock_sr):
    """Return the ordered list of action names sent via a patched send_request."""
    return [c.args[0] for c in mock_sr.call_args_list]


# ── socket.socket fixture ─────────────────────────────────────────────────────

@pytest.fixture
def mock_socket():
    """Mock the socket context manager; default recv returns a success payload."""
    with patch("socket.socket") as mock_cls:
        inst = MagicMock()
        mock_cls.return_value.__enter__.return_value = inst
        inst.recv.return_value = json.dumps(
            {"status": "success", "result": "Mocked execution"}
        ).encode()
        yield inst


# ── send_request unit tests ───────────────────────────────────────────────────

def test_send_request_success(mock_socket):
    """Low-level wire format and return value."""
    res = send_request("test_method", args=["arg1"], kwargs={"kw1": "val1"})
    assert res == {"status": "success", "result": "Mocked execution"}

    mock_socket.sendall.assert_called_once()
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "test_method"
    assert payload["args"] == ["arg1"]
    assert payload["kwargs"] == {"kw1": "val1"}


def test_send_request_connection_refused():
    """ConnectionRefusedError returns a structured error dict."""
    with patch("socket.socket") as mock_cls:
        mock_cls.return_value.__enter__.return_value.connect.side_effect = (
            ConnectionRefusedError("Connection refused")
        )
        res = send_request("test")
    assert res["status"] == "error"
    assert "Socket connection failed" in res["error"]
    assert "Connection refused" in res["error"]


def test_send_request_timeout():
    """Socket timeout returns a structured error dict."""
    with patch("socket.socket") as mock_cls:
        mock_cls.return_value.__enter__.return_value.connect.side_effect = (
            socket.timeout("Timed out")
        )
        res = send_request("test")
    assert res["status"] == "error"
    assert "Timed out" in res["error"]


# ── Primitive tool wrappers ───────────────────────────────────────────────────

def test_tool_show(mock_socket):
    result = show(representation="cartoon", selection="chain A")
    assert "Showing cartoon" in result
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "show"
    assert payload["args"] == ["cartoon", "chain A"]


def test_tool_color(mock_socket):
    result = color(color_name="red", selection="all")
    assert "Colored selection 'all' with red" in result
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "color"
    assert payload["args"] == ["red", "all"]


def test_tool_select(mock_socket):
    result = select(name="my_selection", selection="resi 1-10")
    assert "Created named selection" in result
    assert mock_socket.sendall.call_count == 1
    payload = json.loads(mock_socket.sendall.call_args[0][0])
    assert payload["action"] == "select"
    assert payload["args"] == ["my_selection", "resi 1-10"]


def test_tool_error_propagation(mock_socket):
    """Plugin errors surface back to the MCP caller."""
    mock_socket.recv.return_value = json.dumps({
        "status": "error",
        "error": "PyMOL encountered a problem: invalid selection",
    }).encode()
    result = show(representation="spheres")
    assert "invalid selection" in result
    assert result.startswith("PyMOL encountered")


# ── fetch_structure ───────────────────────────────────────────────────────────

@patch("mcpymol.server.send_request")
def test_fetch_structure_multimer(mock_sr):
    """Chains A/B/C → remove non-proximal chains and apply ghost heart."""
    mock_sr.side_effect = _sr_mock(get_chains=["A", "B", "C"])

    result = fetch_structure(pdb_code="1ubq")

    assert "Successfully fetched 1ubq" in result
    acts = _actions(mock_sr)
    assert acts[0] == "do"          # reinitialize
    assert "fetch" in acts
    assert "get_chains" in acts
    assert "remove" in acts         # multimer cleanup
    assert "hide" in acts           # solvent hidden


@patch("mcpymol.server.send_request")
def test_fetch_structure_single_chain(mock_sr):
    """Single chain still applies the keep-selection + ghost heart."""
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    result = fetch_structure(pdb_code="1ubq")

    assert "Successfully fetched 1ubq" in result
    acts = _actions(mock_sr)
    assert "remove" in acts
    assert "hide" in acts


@patch("mcpymol.server.send_request")
def test_fetch_structure_no_chains(mock_sr):
    """Empty chain list falls through to the fallback return message."""
    mock_sr.side_effect = _sr_mock(get_chains=[])

    result = fetch_structure(pdb_code="1abc")

    assert "no chains were found" in result
    acts = _actions(mock_sr)
    assert "fetch" in acts
    assert "remove" not in acts     # no cleanup attempted


@patch("mcpymol.server.send_request")
def test_fetch_structure_error(mock_sr):
    """Fetch failure is propagated as a readable error string."""
    def fake(action, args=None, kwargs=None):
        if action == "fetch":
            return {"status": "error", "error": "PDB ID not found"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake

    result = fetch_structure(pdb_code="XXXX")

    assert "Error fetching XXXX" in result
    assert "PDB ID not found" in result


@patch("mcpymol.server.send_request")
def test_fetch_structure_custom_obj_name(mock_sr):
    """Custom obj_name is threaded through all downstream send_request calls."""
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    fetch_structure(pdb_code="1ubq", obj_name="my_protein")

    for c in mock_sr.call_args_list:
        if c.args[0] == "fetch":
            assert c.kwargs["args"] == ["1ubq", "my_protein"]
            break
    else:
        pytest.fail("fetch action was never called")


# ── load_structure ────────────────────────────────────────────────────────────

@patch("mcpymol.server.send_request")
def test_load_structure_success(mock_sr):
    mock_sr.side_effect = _sr_mock(get_chains=["A"])

    result = load_structure(file_path="/tmp/test.pdb", obj_name="test")

    assert "Successfully loaded" in result
    acts = _actions(mock_sr)
    assert "load" in acts
    assert "remove" in acts
    assert "hide" in acts


@patch("mcpymol.server.send_request")
def test_load_structure_error(mock_sr):
    def fake(action, args=None, kwargs=None):
        if action == "load":
            return {"status": "error", "error": "File not found"}
        return {"status": "success", "result": "OK"}
    mock_sr.side_effect = fake

    result = load_structure(file_path="/bad/path.pdb", obj_name="test")

    assert "Error loading" in result
    assert "File not found" in result


# ── View functions ────────────────────────────────────────────────────────────

@patch("mcpymol.server.send_request")
def test_ligand_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = ligand_view(obj_name="1ATP", ligand_resn="ATP")

    assert "ATP" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    assert "zoom" in acts
    # H-bonds drawn via a 'do distance ...' call
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("distance" in a for a in do_args)


@patch("mcpymol.server.send_request")
def test_interface_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = interface_view(obj_name="1BRS", chain_a="A", chain_b="D")

    assert "chain A" in result and "chain D" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("distance" in a for a in do_args)


@patch("mcpymol.server.send_request")
def test_putty_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = putty_view(obj_name="1UBQ")

    assert "Putty view" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts
    assert "set" in acts


@patch("mcpymol.server.send_request")
def test_hydrophobic_surface_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = hydrophobic_surface_view(obj_name="1TCA")

    assert "Hydrophobic" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts


@patch("mcpymol.server.send_request")
def test_electrostatic_view_atomic(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = electrostatic_view(obj_name="1LYZ")

    assert "Electrostatic" in result
    assert "atomic" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts


@patch("mcpymol.server.send_request")
def test_electrostatic_view_residue(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = electrostatic_view(obj_name="1LYZ", mode="residue")

    assert "residue" in result
    assert "Electrostatic" in result


@patch("mcpymol.server.send_request")
def test_crosslink_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = crosslink_view(obj_name="1CEL")

    assert "Crosslink" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    do_args = [c.kwargs["args"][0] for c in mock_sr.call_args_list if c.args[0] == "do"]
    assert any("distance" in a for a in do_args)


@patch("mcpymol.server.send_request")
def test_pocket_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = pocket_view(obj_name="1HSG", resn="MK1")

    assert "MK1" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    assert "zoom" in acts


@patch("mcpymol.server.send_request")
def test_pharmacophore_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = pharmacophore_view(obj_name="1HSG", resn="MK1")

    assert "Pharmacophore" in result
    assert "MK1" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "zoom" in acts


@patch("mcpymol.server.send_request")
def test_mutation_view_valid(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = mutation_view(obj_name="4HHB", mutations="E6V,K16E")

    assert "Mutation view" in result
    assert "E6V" in result or "K16E" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "label" in acts
    assert "zoom" in acts


@patch("mcpymol.server.send_request")
def test_textbook_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = textbook_view(obj_name="1ABC")

    assert "Textbook Illustration" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "color" in acts
    assert "set" in acts


@patch("mcpymol.server.send_request")
def test_cinematic_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = cinematic_view(obj_name="1ABC")

    assert "Cinematic view" in result
    acts = _actions(mock_sr)
    assert "show" in acts
    assert "do" in acts
    assert "set" in acts


@patch("mcpymol.server.send_request")
def test_pointillist_view(mock_sr):
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = pointillist_view(obj_name="1ABC")

    assert "Pointillist/Starfield" in result
    acts = _actions(mock_sr)
    assert "hide" in acts
    assert "show" in acts
    assert "set" in acts
    assert "color" in acts


@patch("mcpymol.server.send_request")
def test_mutation_view_invalid_input(mock_sr):
    """Non-parseable mutation strings return an error before any PyMOL calls."""
    mock_sr.return_value = {"status": "success", "result": "OK"}

    result = mutation_view(obj_name="4HHB", mutations="bad_input_no_digits")

    assert "No valid mutations" in result
    assert mock_sr.call_count == 0
