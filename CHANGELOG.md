# Changelog

All notable changes to MCPymol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-03-30

### Added
- 50+ auto-generated PyMOL commands exposed as MCP tools (`show`, `hide`, `color`, `distance`, `get_chains`, `select`, and more)
- Biological assembly fetching with automatic BFS-based multimer heuristic — isolates the functional multimer while discarding crystallographic copies
- Automatic solvent/water hiding for clean, relevant views
- Dual-process socket bridge architecture to work around PyMOL's internal Python environment constraints
- `MCPYMOL_PORT` environment variable for running multiple instances simultaneously
- Auto-start support via `~/.pymolrc.py`
- Tested and supported with Claude Code and Gemini CLI
