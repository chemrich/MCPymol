import pytest
import json
import socket
from unittest.mock import patch, MagicMock
from mcpymol.server import send_request, fetch_structure, show, color, select, distance

@pytest.fixture
def mock_socket():
    """Fixture to mock the socket connection and provide a predictable response."""
    with patch("socket.socket") as mock_socket_class:
        mock_sock_instance = MagicMock()
        
        # When entering the context manager `with socket.socket() as s:`
        mock_socket_class.return_value.__enter__.return_value = mock_sock_instance
        
        # Default mock response: just a basic success JSON
        mock_sock_instance.recv.return_value = json.dumps({"status": "success", "result": "Mocked execution"}).encode('utf-8')
        yield mock_sock_instance

def test_send_request_success(mock_socket):
    """Test the low-level send_request wrapper function."""
    res = send_request("test_method", args=["arg1"], kwargs={"kw1": "val1"})
    
    assert res == {"status": "success", "result": "Mocked execution"}
    
    # Verify what was sent over the wire
    mock_socket.sendall.assert_called_once()
    sent_data = mock_socket.sendall.call_args[0][0]
    payload = json.loads(sent_data.decode('utf-8'))
    
    assert payload["action"] == "test_method"
    assert payload["args"] == ["arg1"]
    assert payload["kwargs"] == {"kw1": "val1"}

def test_send_request_connection_refused():
    """Test that send_request gracefully handles a downed PyMOL plugin socket."""
    with patch("socket.socket") as mock_socket_class:
        mock_socket_class.return_value.__enter__.return_value.connect.side_effect = ConnectionRefusedError("Connection refused")
        res = send_request("test")
        
        assert res["status"] == "error"
        assert "Socket connection failed" in res["error"]
        assert "Connection refused" in res["error"]

def test_send_request_timeout():
    """Test that socket timeouts return clear error messages."""
    with patch("socket.socket") as mock_socket_class:
        mock_socket_class.return_value.__enter__.return_value.connect.side_effect = socket.timeout("Timed out")
        res = send_request("test")
        
        assert res["status"] == "error"
        assert "Timed out" in res["error"]


# -- Tool Integration Tests --

def test_tool_show(mock_socket):
    """Test the 'show' MCP tool wrapper."""
    result = show(representation="cartoon", selection="chain A")
    
    assert "Showing cartoon" in result
    
    # Analyze the JSON payload sent by the tool
    assert mock_socket.sendall.call_count == 1
    sent_data = mock_socket.sendall.call_args[0][0]
    payload = json.loads(sent_data.decode('utf-8'))
    
    assert payload["action"] == "show"
    assert payload["args"] == ["cartoon", "chain A"]

def test_tool_color(mock_socket):
    """Test the 'color' MCP tool wrapper."""
    result = color(color_name="red", selection="all")
    assert "Colored selection 'all' with red" in result
    
    assert mock_socket.sendall.call_count == 1
    sent_data = mock_socket.sendall.call_args[0][0]
    payload = json.loads(sent_data.decode('utf-8'))
    assert payload["action"] == "color"
    assert payload["args"] == ["red", "all"]

def test_tool_select(mock_socket):
    """Test the 'select' MCP tool wrapper."""
    result = select(name="my_selection", selection="resi 1-10")
    assert "Created named selection" in result
    
    assert mock_socket.sendall.call_count == 1
    sent_data = mock_socket.sendall.call_args[0][0]
    payload = json.loads(sent_data.decode('utf-8'))
    assert payload["action"] == "select"
    assert payload["args"] == ["my_selection", "resi 1-10"]

def test_tool_fetch_structure_multimer(mock_socket):
    """Test the complex fetch_structure logic which loops internally."""
    
    # We need to simulate multiple responses for fetch_structure
    # 1. First fetch command succeeds
    # 2. get_chains returns multiple chains
    # 3. distance command simulates heuristics
    # 4. remove commands simulate cleanup
    
    responses = [
        json.dumps({"status": "success", "result": "Fetched"}).encode('utf-8'),
        json.dumps({"status": "success", "result": ["A", "B", "C"]}).encode('utf-8'),
        json.dumps({"status": "success", "result": "Removed"}).encode('utf-8'), # Remove B
        json.dumps({"status": "success", "result": "Removed"}).encode('utf-8'), # Remove C
        json.dumps({"status": "success", "result": "Hidden"}).encode('utf-8'),  # Hide Solvents
    ]
    mock_socket.recv.side_effect = responses

    result = fetch_structure(pdb_code="1ubq")
    
    assert "Successfully fetched 1ubq" in result
    
    # We should have sent 4 distinct socket payloads (fetch, get_chains, remove, hide_solvents)
    assert mock_socket.sendall.call_count == 4
    
    calls = mock_socket.sendall.call_args_list
    # The get_chains call
    chains_payload = json.loads(calls[1][0][0].decode('utf-8'))
    assert chains_payload["action"] == "get_chains"
    
def test_tool_error_propagation(mock_socket):
    """Test that errors from the PyMOL plugin are surfaced to the MCP client."""
    
    mock_socket.recv.return_value = json.dumps({
        "status": "error", 
        "error": "PyMOL encountered a problem: invalid selection"
    }).encode('utf-8')
    
    result = show(representation="spheres")
    
    assert "invalid selection" in result
    assert result.startswith("PyMOL encountered")
