# Changelog

All notable changes to MCPymol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2026-05-15

### Fixed
- **Bridge framing.** Both the in-PyMOL plugin and the external bridge now drain TCP responses to EOF (with incremental JSON parsing as a fallback for mock-style peers), instead of truncating at the first 8 KB chunk. Long PyMOL responses — `get_fastastr` on multi-hundred-residue chains, `get_chains` on large assemblies, error tracebacks — no longer corrupt the JSON.
- **`util.*` tool dispatch.** The plugin previously checked `hasattr(cmd, action)` for every action, so dotted names like `util.cbc` / `util.cbaw` / `util.chainbow` silently failed even though the tools were registered. Plugin now resolves dotted names through `pymol.util` (and falls back to a general attribute walk).
- **`conservation_view` residue mapping.** Previous version assumed `resi == i + 1` along the FASTA, which silently misaligned scores in structures with non-contiguous residue numbering (gaps, modified termini). Now walks the actual CA `resi` values from PyMOL and maps via a stored dict.

### Performance
- **`conservation_view`** alter loop collapsed from O(2N) socket round-trips to a single batched `cmd.do` script. For a 300-residue chain this drops from ~10 s of TCP overhead alone to one round-trip.

### Added
- `list_objects`, `list_chains(obj_name)`, `list_ligands(obj_name)` introspection tools so models can ground themselves in actual session state instead of guessing object names, chain IDs, or 3-letter ligand codes.
- `python -m mcpymol` entry point via `__main__.py`.
- GitHub Actions CI running pytest on Python 3.10–3.13.
- `print_export` tool — exports a structure as watertight, manifold STL files for multi-colour 3D printing. Per-colour-group isolation works around PyMOL's whole-scene OBJ export. Adaptive mesh repair: `auto` does a light cleanup when the export is already watertight (compact barrels like GFP — keeps the largest body, drops internal cavity shells), otherwise screened-Poisson reconstruction with a voxel-remesh fallback (robust for thin nucleic acids); all groups stay in one coordinate frame for slicer assembly. Optional `print` extra (trimesh, pymeshlab, scipy, scikit-image, networkx); degrades gracefully with an install hint when the libraries are absent.

### Changed
- Tool descriptions for `show`, `hide`, `color`, `select`, `remove`, `distance`, `execute_pymol_command` now enumerate valid argument vocabularies (representation names, color names) and include a brief PyMOL selection-syntax primer. Stronger guardrail on `execute_pymol_command` so models reach for it less.
- `pyproject.toml` enriched with authors, urls, classifiers, keywords, and a `[tool.pytest.ini_options]` block so `pytest` Just Works from the repo root.
- `.gitignore` extended to cover `venv/`, `refresh/`, `.vscode/`, `build/`, `dist/`.
- README rewritten: fixed duplicate "Option C", fixed `yourusername/MCPymol` placeholder, added missing views (`bfactor_view`, `textbook_view`, `cinematic_view`, `pointillist_view`, `conservation_view`), added a how-it-talks architecture diagram, a troubleshooting table, and a "Try it" prompt list.
- New `CONTRIBUTING.md`.

## [1.1.0] - 2026-03-31

### Added
- `conservation_view` tool — evolutionary conservation visualization using Shannon entropy
- MMseqs2 integration via ColabFold public API (configurable for local servers via `MCPYMOL_MMSEQS_URL` env var)
- Full pipeline: sequence extraction → MSA generation → entropy scoring → B-factor mapping → spectrum coloring
- A3M parser with insertion-stripping for clean MSA alignment
- Per-residue Shannon entropy calculation normalized to [0, 1]
- 18 new tests covering A3M parsing, entropy math, API mocking, and end-to-end conservation_view

## [1.0.0] - 2026-03-30

### Added
- 50+ auto-generated PyMOL commands exposed as MCP tools (`show`, `hide`, `color`, `distance`, `get_chains`, `select`, and more)
- Biological assembly fetching with automatic BFS-based multimer heuristic — isolates the functional multimer while discarding crystallographic copies
- Automatic solvent/water hiding for clean, relevant views
- Dual-process socket bridge architecture to work around PyMOL's internal Python environment constraints
- `MCPYMOL_PORT` environment variable for running multiple instances simultaneously
- Auto-start support via `~/.pymolrc.py`
- Tested and supported with Claude Code and Gemini CLI
