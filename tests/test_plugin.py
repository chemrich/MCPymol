import pytest
import json
from unittest.mock import patch, MagicMock

# We need to mock pymol.cmd before we can even import plugin.py
import sys
mock_pymol = MagicMock()
sys.modules['pymol'] = mock_pymol
sys.modules['pymol.cmd'] = mock_pymol.cmd

# Now we can import the plugin
from mcpymol.plugin import PyMOLSocketServer

@pytest.fixture
def plugin_server():
    """Initializes a PyMOLSocketServer instance ready for simulated payloads."""
    server = PyMOLSocketServer()
    # Reset the mock cmd calls before every test so they are isolated
    mock_pymol.cmd.reset_mock()
    return server

def test_handle_load_structure(plugin_server):
    """Test standard file loading mapped correctly to cmd.load."""
    payload = json.dumps({"action": "load", "args": ["1ubq.pdb"]})
    
    res = plugin_server.handle_request(payload)
    
    assert res["status"] == "success"
    assert "Loaded 1ubq.pdb" in res["result"]
    mock_pymol.cmd.load.assert_called_once_with("1ubq.pdb")

def test_handle_get_chains(plugin_server):
    """Test custom get_chains implementation."""
    payload = json.dumps({"action": "get_chains", "args": ["all"]})
    mock_pymol.cmd.get_chains.return_value = ["A", "B", "X"]
    
    res = plugin_server.handle_request(payload)
    
    assert res["status"] == "success"
    assert res["result"] == ["A", "B", "X"]
    mock_pymol.cmd.get_chains.assert_called_once_with("all")

def test_dynamic_method_resolution(plugin_server):
    """Test that arbitrary PyMOL cmd functions (e.g., cmd.hide) are dynamically invoked."""
    
    # We add a fake callable attribute to the mock for testing dynamic resolution
    def fake_hide(*args, **kwargs):
        pass
    mock_pymol.cmd.hide = fake_hide
    
    payload = json.dumps({"action": "hide", "args": ["everything", "all"]})
    res = plugin_server.handle_request(payload)
    
    assert res["status"] == "success"
    assert "Executed 'hide' successfully" in res["result"]

def test_unknown_action(plugin_server):
    """Test behavior when an invalid action is requested."""
    payload = json.dumps({"action": "nonexistent_pymol_command"})
    
    # We implicitly mock missing methods by explicitly deleting them on the MagicMock
    # to force `hasattr(cmd, action)` to return False
    del mock_pymol.cmd.nonexistent_pymol_command
    
    res = plugin_server.handle_request(payload)
    
    assert res["status"] == "error"
    assert "Unknown action or method not found on cmd" in res["error"]

def test_malformed_json(plugin_server):
    """Test handling of invalid JSON payloads from a buggy client."""
    payload = '{"action": "load", "args": ["1ubq.pdb"]' # Missing closing brace
    
    res = plugin_server.handle_request(payload)
    
    assert res["status"] == "error"
    assert "Invalid JSON" in res["error"]

def test_pymol_internal_exception(plugin_server):
    """Test handling of PyMOL execution errors."""
    
    mock_pymol.cmd.color.side_effect = Exception("Invalid color name 'bleargh'")
    
    payload = json.dumps({"action": "color", "args": ["bleargh", "all"]})
    res = plugin_server.handle_request(payload)
    
    assert res["status"] == "error"
    assert "Invalid color name 'bleargh'" in res["error"]
    assert "Invalid color name 'bleargh'" in res["error"]
