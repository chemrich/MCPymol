import json
import socket
from typing import Optional
from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("MCPymol")

HOST = '127.0.0.1'
PORT = 9876

def send_request(action: str, args: list = None, kwargs: dict = None) -> dict:
    """Send a JSON request to the PyMOL plugin socket server."""
    payload = {
        "action": action,
        "args": args or [],
        "kwargs": kwargs or {}
    }
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(2.0)
            s.connect((HOST, PORT))
            s.sendall(json.dumps(payload).encode('utf-8'))
            data = s.recv(8192).decode('utf-8')
            return json.loads(data)
    except Exception as e:
        return {"status": "error", "error": f"Socket connection failed: {e}. Is the PyMOL plugin running?"}

@mcp.tool()
def fetch_structure(pdb_code: str, obj_name: Optional[str] = None) -> str:
    """
    Fetches a protein structure from the PDB by its 4-letter code, and applies the multimer heuristic.
    The heuristic automatically removes all chains except the first chain and any other chains
    that are within 5 angstroms of the first chain (preserving associated multimers).
    """
    name = obj_name if obj_name else pdb_code
    
    res = send_request("fetch", args=[pdb_code, name])
    if res.get("status") == "error":
        return f"Error fetching {pdb_code}: {res.get('error')}"
        
    # Apply heuristic
    chains_res = send_request("get_chains", args=[f"({name})"])
    if chains_res.get("status") == "success":
        chains = chains_res.get("result", [])
        if isinstance(chains, list) and chains:
            first = chains[0]
            keep_sel = f"({name} and chain {first}) or bychain ({name} and chain {first} around 5)"
            send_request("remove", args=[f"({name}) and not ({keep_sel})"])
            return f"Successfully fetched {pdb_code} as '{name}' and applied multimer chain heuristic (kept chain {first} and nearby chains)."
    
    return f"Successfully fetched {pdb_code} as '{name}', but no chains were found to apply heuristic."

@mcp.tool()
def load_structure(file_path: str, obj_name: str) -> str:
    """
    Loads a structure from a local file path and applies the multimer heuristic.
    """
    res = send_request("load", args=[file_path, obj_name])
    if res.get("status") == "error":
        return f"Error loading {file_path}: {res.get('error')}"
        
    # Apply heuristic
    chains_res = send_request("get_chains", args=[f"({obj_name})"])
    if chains_res.get("status") == "success":
        chains = chains_res.get("result", [])
        if isinstance(chains, list) and chains:
            first = chains[0]
            keep_sel = f"({obj_name} and chain {first}) or bychain ({obj_name} and chain {first} around 5)"
            send_request("remove", args=[f"({obj_name}) and not ({keep_sel})"])
            return f"Successfully loaded {file_path} as '{obj_name}' and applied multimer chain heuristic."
    
    return f"Loaded {file_path} successfully."

@mcp.tool()
def show(representation: str, selection: str = "all") -> str:
    """Shows a representation for a given selection."""
    res = send_request("show", args=[representation, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Showing {representation} for selection '{selection}'"

@mcp.tool()
def hide(representation: str, selection: str = "all") -> str:
    """Hides a representation for a given selection."""
    res = send_request("hide", args=[representation, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Hiding {representation} for selection '{selection}'"

@mcp.tool()
def color(color_name: str, selection: str = "all") -> str:
    """Sets the color for the specified selection."""
    res = send_request("color", args=[color_name, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Colored selection '{selection}' with {color_name}"
        
@mcp.tool()
def select(name: str, selection: str) -> str:
    """Creates a named selection for subsequent use."""
    res = send_request("select", args=[name, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Created named selection '{name}' for '{selection}'"

@mcp.tool()
def remove(selection: str) -> str:
    """Removes atoms or structures matching the selection."""
    res = send_request("remove", args=[selection])
    if res.get("status") == "error": return res.get("error")
    return f"Removed selection '{selection}'"

@mcp.tool()
def distance(name: str, selection1: str, selection2: str) -> str:
    """Measures and displays the distance between two selections."""
    res = send_request("distance", args=[name, selection1, selection2])
    if res.get("status") == "error": return res.get("error")
    return f"Measured distance between '{selection1}' and '{selection2}' as '{name}'"

@mcp.tool()
def execute_pymol_command(command: str) -> str:
    """Executes a raw PyMOL command string. Use this only for commands not covered by primitive tools."""
    res = send_request("do", args=[command])
    if res.get("status") == "error": return res.get("error")
    return f"Executed command: {command}"

def main():
    mcp.run()

if __name__ == "__main__":
    main()
