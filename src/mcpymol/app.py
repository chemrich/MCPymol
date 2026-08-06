"""The shared FastMCP application object.

Kept in its own module so every tool module can register against the same
instance without importing each other — the tool modules import ``mcp`` from
here, and ``mcpymol.server`` imports the tool modules.
"""

from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("MCPymol")
