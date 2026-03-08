"""
PyMOL Native Socket Plugin Bridge for MCPymol
Run this script inside PyMOL:
    run /path/to/MCPymol/src/mcpymol/plugin.py
"""
import socket
import json
import threading
import traceback

try:
    from pymol import cmd
except ImportError:
    print("Warning: pymol.cmd not found. This script must be run inside PyMOL.")
    cmd = None

class PyMOLSocketServer:
    def __init__(self, host='127.0.0.1', port=9876):
        self.host = host
        self.port = port
        self.running = False
        self.thread = None
        self.server_socket = None

    def handle_request(self, payload):
        """Execute the requested PyMOL command and return the result."""
        if cmd is None:
            return {"status": "error", "error": "PyMOL cmd module not available."}
            
        action = payload.get("action")
        args = payload.get("args", [])
        kwargs = payload.get("kwargs", {})

        try:
            if action == "do":
                # Execute arbitrary pymol command string
                command_str = args[0] if args else kwargs.get("command", "")
                cmd.do(command_str)
                return {"status": "success", "result": f"Executed command: {command_str}"}
                
            elif action == "fetch":
                # cmd.fetch(code, name, ...)
                cmd.fetch(*args, **kwargs)
                return {"status": "success", "result": f"Fetched {args[0] if args else 'structure'}"}
                
            elif action == "load":
                # cmd.load(filename, object, ...)
                cmd.load(*args, **kwargs)
                return {"status": "success", "result": f"Loaded {args[0] if args else 'structure'}"}
                
            elif action == "get_chains":
                # Custom handler to return a list of chains for a selection
                selection = args[0] if args else kwargs.get("selection", "all")
                chains = cmd.get_chains(selection)
                return {"status": "success", "result": chains}
                
            elif action == "remove":
                cmd.remove(*args, **kwargs)
                return {"status": "success", "result": "Removed selection"}
                
            elif action == "show":
                cmd.show(*args, **kwargs)
                return {"status": "success", "result": "Show executed"}
                
            elif action == "hide":
                cmd.hide(*args, **kwargs)
                return {"status": "success", "result": "Hide executed"}
                
            elif action == "color":
                cmd.color(*args, **kwargs)
                return {"status": "success", "result": "Color executed"}
                
            elif action == "select":
                cmd.select(*args, **kwargs)
                return {"status": "success", "result": "Selection created"}
                
            elif action == "distance":
                cmd.distance(*args, **kwargs)
                return {"status": "success", "result": "Distance measured"}
                
            else:
                return {"status": "error", "error": f"Unknown action: {action}"}
                
        except Exception as e:
            err_msg = str(e)
            print(f"MCPymol Plugin Error executing {action}: {err_msg}")
            traceback.print_exc()
            return {"status": "error", "error": err_msg}

    def _listen_loop(self):
        self.server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        
        try:
            self.server_socket.bind((self.host, self.port))
            self.server_socket.listen(5)
            self.server_socket.settimeout(1.0)
            print(f"MCPymol Native Plugin listening on {self.host}:{self.port}")
            
            while self.running:
                try:
                    client, addr = self.server_socket.accept()
                    data = client.recv(8192).decode('utf-8')
                    if data:
                        try:
                            payload = json.loads(data)
                            response = self.handle_request(payload)
                        except json.JSONDecodeError:
                            response = {"status": "error", "error": "Invalid JSON payload"}
                        
                        client.sendall(json.dumps(response).encode('utf-8'))
                    client.close()
                except socket.timeout:
                    continue
                except Exception as e:
                    print(f"MCPymol Socket Server Error: {e}")
                    
        except Exception as e:
            print(f"MCPymol Failed to start socket server: {e}")
        finally:
            if self.server_socket:
                self.server_socket.close()

    def start(self):
        if not self.running:
            self.running = True
            self.thread = threading.Thread(target=self._listen_loop, daemon=True)
            self.thread.start()
            print("MCPymol Bridge started.")

    def stop(self):
        if self.running:
            self.running = False
            if self.thread:
                self.thread.join(timeout=2.0)
            print("MCPymol Bridge stopped.")

# Auto-start singleton
try:
    if 'mcp_bridge_plugin' in globals():
        mcp_bridge_plugin.stop()
    mcp_bridge_plugin = PyMOLSocketServer()
    mcp_bridge_plugin.start()
    
    # Expose stop function to PyMOL CLI
    cmd.extend("stop_mcp", mcp_bridge_plugin.stop)
    cmd.extend("start_mcp", mcp_bridge_plugin.start)
except Exception as e:
    print(f"Failed to initialize MCPymol plugin: {e}")
