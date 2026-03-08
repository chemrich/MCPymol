# MCPymol: An MCP Server for PyMOL

**MCPymol** is a Model Context Protocol (MCP) server that provides a conversational interface for viewing and analyzing protein structures using PyMOL. It exposes PyMOL's powerful molecular visualization capabilities to AI assistants (like Claude), allowing you to seamlessly load structures, manipulate views, and explore proteins using natural language.

## 🧬 What it Does

MCPymol acts as a bridge between an AI assistant and a running PyMOL desktop instance. It provides:
- **50+ Auto-Generated PyMOL Commands**: Claude has direct access to PyMOL primitives like `show`, `hide`, `color`, `distance`, `get_chains`, `select`, and many more.
- **Smart Multimer & Solvent Heuristics**: When fetching or loading structures, MCPymol automatically applies heuristics to isolate the primary chain (and any chains within 5Å of it) while hiding waters, non-standard crystallization additives, and solvents. This ensures a clean, relevant view of the protein structure immediately upon loading.
- **Dual-Process Architecture**: To circumvent PyMOL's internal Python dependency limitations, MCPymol runs via a two-part bridge. A native PyMOL script runs a lightweight background socket listener inside your PyMOL App, while a standalone FastMCP server handles communication with the AI assistant.

## 🛠️ Intended Use

MCPymol is designed for structural biologists, bioinformaticians, and anyone interested in protein structures who wants a natural language conversational partner for PyMOL. You can ask Claude to:
- "Fetch ubiquitin (1ubq) and show it as a cartoon."
- "Color the alpha helices red and the beta sheets blue."
- "Measure the distance between residue 10 and residue 20."
- "Highlight the active site."

## 💾 Installation

This project relies on `uv` for fast dependency management. 

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/MCPymol.git
   cd MCPymol
   ```

2. **Install dependencies:**
   Ensure you have `uv` installed, then run:
   ```bash
   uv sync
   ```

## 🚀 Usage

Because PyMOL requires its own isolated Python environment to run its GUI and rendering loops, using MCPymol requires starting both the PyMOL plugin and the external MCP server.

### 1. Start the Native PyMOL Plugin
1. Open your standard PyMOL desktop application.
2. In the PyMOL command line, manually initialize the background listener script. Adjust the path to where you cloned the repository:
   ```pymol
   run /path/to/MCPymol/src/mcpymol/plugin.py
   ```
   *You should see a message in the PyMOL console indicating the plugin is listening on `127.0.0.1:9876`.*

**💡 Pro Tip: Auto-Start the Plugin**
The standard way to run initialization scripts in PyMOL is through its `pymolrc` resource file. To automatically run this plugin every time you launch PyMOL, add the following to `~/.pymolrc.py`:
```python
import os, pymol
pymol.cmd.do("run /absolute/path/to/MCPymol/src/mcpymol/plugin.py")
```

**🔌 Changing the Port**
The default port is `9876`. If you need to run multiple MCP servers simultaneously, override it with the `MCPYMOL_PORT` environment variable before launching PyMOL **and** the bridge:
```bash
# PyMOL (macOS example)
MCPYMOL_PORT=9867 open -a PyMOL

# MCP bridge
MCPYMOL_PORT=9867 uv run mcpymol
```

### 2. Start the MCP Server Bridge
From your terminal, in the root of the cloned `MCPymol` directory, start the server:
```bash
uv run mcpymol
```

### 3. Configure Claude Desktop
To interact with the server via Claude Desktop, add the following to your `claude_desktop_config.json` file:

```json
{
  "mcpServers": {
    "mcpymol": {
      "command": "uv",
      "args": [
        "--directory",
        "/absolute/path/to/MCPymol",
        "run",
        "mcpymol"
      ]
    }
  }
}
```
Restart Claude Desktop, ensure the PyMOL plugin is running, and you're ready to start asking conversational questions about protein structures!

## 🧪 Running Tests
The repository includes a rigorous, "Google Engineer" grade `pytest` suite testing both the socket payload generation and simulated PyMOL API execution boundaries.
To run the automated tests:
```bash
PYTHONPATH=src uv run pytest tests/
```
