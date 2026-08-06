# Changelog

All notable changes to MCPymol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.0] - 2026-08-07

### Added
- **`render`** — ray-traces the current scene and returns the image as MCP image content, so the model can actually see what it built and iterate. Previously it had to call `ray`, then `png`, and was left holding a filename it could not inspect. One round trip (`cmd.png(..., ray=1)`), waits for the PNG's IEND terminator rather than mere existence (PyMOL writes from its own thread and can lag the response that reported success), and declines to inline anything over `MCPYMOL_MAX_IMAGE_BYTES` (5 MB) since base64 inflates by a third.
- **`turntable`** — a 360° PNG sequence plus the `ffmpeg` line to assemble it. Sets the rotation origin to the object first, or the camera orbits the scene centre and the model appears to wobble. Defaults to the OpenGL renderer.
- **`fetch_alphafold`** and AlphaFold routing in `fetch_structure` — a UniProt accession or `AF-` prefixed identifier (`P69905`, `af-P69905`, `AF-P69905-F1-model_v4`) fetches from AlphaFold DB. PyMOL's `cmd.fetch` only knows the RCSB, so the model is downloaded over HTTP and loaded from a temp file. The anchored UniProt regex cannot match a 4-character PDB code, so the namespaces stay separate.
- **`plddt_view`** — colours a predicted model by pLDDT in AlphaFold's official palette (dark blue >90, light blue 70–90, yellow 50–70, orange <50) and reports the confidence breakdown. This matters because pLDDT rides in the B-factor column, so `bfactor_view` and `putty_view` read these models exactly backwards — they treat low values as *rigid*, whereas low pLDDT means the prediction is unreliable. Warns when the B-factors do not look like pLDDT.
- **`superposition_view`** — superposes two structures and colours the mobile one by per-residue CA shift (blue unchanged → red most shifted), naming the worst-shifted residues. An RMSD says a structure moved; this says where. Copies both objects before reading coordinates, because `super`/`align` leave the fit in the object's transformation matrix and the stored coordinates can still be pre-superposition.
- **`save_session` / `load_session`** — `.pse` round-trips the whole scene (objects, selections, representations, colours, scenes, camera), unlike `save`, which writes bare coordinates. `load_session(merge=True)` adds to the current session instead of replacing it.

### Fixed
- **A stalled client could wedge PyMOL for the rest of the session.** The plugin's accept loop is single-threaded on purpose (`pymol.cmd` is not thread-safe), but nothing bounded one connection: `_recv_all` read until the peer half-closed, with no timeout on the *accepted* socket — only the listening one. A client that connected and stalled blocked every subsequent command. Accepted sockets now get `MCPYMOL_RECV_TIMEOUT` (30 s), requests are capped at `MCPYMOL_MAX_REQUEST_BYTES` (8 MB), and every failure path answers instead of leaving the bridge to guess at an empty response.
- **A non-JSON-serializable `cmd` return killed the response.** `json.dumps` raised and the client got nothing back; it now degrades to a repr.
- **`ray`, `draw`, `mpng`, `png` and `save` used the 10 s interactive socket timeout**, which a 1920×1080 ray trace cannot possibly meet — they now use `MCPYMOL_SLOW_OP_TIMEOUT` (600 s). These were exactly the calls worth waiting for.
- **`pdb2pqr`/`apbs` ran without a timeout**, so a wedged solver hung the whole MCP server. Bounded by `MCPYMOL_PB_TIMEOUT` (600 s); a missing binary now reports its install command.
- **`poisson_boltzmann_view` ignored the result of the save it depends on**, so a failure surfaced as a baffling PDB2PQR error about a file that was never written.
- **Optional arguments could be silently mislabelled.** Every thin wrapper dropped *all* `None` values, collapsing gaps in PyMOL's positional argument list — `ray(height="1080")` with no width sent `["1080"]`, and PyMOL read the height as the width. Unset *trailing* arguments are still dropped; a gap now returns an error naming the missing parameter.
- **Importing `mcpymol.plugin` bound TCP port 9876**, racing any real PyMOL session. `run plugin.py` inside PyMOL still auto-starts; a library import (the test suite, tooling) no longer does.
- The multimer cutoff disagreed with itself: the helper defaulted to 5.0 Å while both callers passed 8.0 Å, and the README documented 5 Å. Unified on `DEFAULT_MULTIMER_CUTOFF = 8.0`.

### Changed
- **`server.py` split into a package.** It had reached 3,207 lines mixing the socket protocol, structure loading, 77 command wrappers, 14 scene presets, an MMseqs2 pipeline and mesh repair. Now `app`, `bridge`, `structures`, `primitives`, `views`, `rendering`, `comparison`, `conservation` and `printing`, with `server.py` as the entry point and re-export facade — `from mcpymol.server import ligand_view` still works. Verified by diffing `mcp.list_tools()` across the change: every tool kept a byte-identical name, description and input schema.
- **77 duplicated wrapper bodies collapsed into one shared helper.** Each wrapper is still a real `def` — the signature and docstring are the tool's public interface — but the mechanical forwarding lives in `_call`, so a fix lands once instead of 77 times.
- The `typecheck` CI job now installs the `print` extra. Without it trimesh resolved to `Any` and mypy silently skipped the mesh-repair code entirely; enabling it surfaced two real errors.
- Accept-loop poll interval cut from 1.0 s to 0.25 s, so `stop_mcp` returns promptly.

### Tests
- 155 → 346. New coverage for the Poisson/voxel repair paths, the plugin's framing and every `serve_connection` failure path, the solver error handling, and all six new tools.
- New `tests/test_bridge_roundtrip.py` runs the real bridge against the real plugin listener over TCP on an ephemeral port. Everything else mocks one side or the other, so the wire format was only ever checked against a mock written to match the implementation. It immediately found two bugs (an unoverridable size cap, and the listener reporting the requested port rather than the bound one).
- New `tests/test_package_layout.py` enforces that a tool added to a module is re-exported, that nothing registers against a second FastMCP instance, and that no module is left unimported (which would silently drop its tools).

## [1.2.1] - 2026-05-18

### Added
- `.claude/skills/mcpymol-guide/SKILL.md` — a committed Claude Code skill documenting how to drive the mcpymol MCP server: structure prep and the multimer heuristic, the `*_view` preset catalogue, PyMOL selection syntax, rendering, and the 3D-print STL workflow (including the `representation="cartoon"` vs surface trap, verifying the STL with trimesh rather than the viewport, and the `/mcp` reconnect requirement after `server.py` edits). `.gitignore` now tracks `.claude/skills/` while the rest of `.claude/` stays local.

### Fixed
- `uv.lock` was not regenerated when the version was bumped to 1.2.0, so CI's `uv sync --locked` failed on every branch. The lockfile is back in sync with `pyproject.toml`.

## [1.2.0] - 2026-05-18

### Added
- `print_ribbon_view` tool — a 3D-print preset that pairs chunky β-strand arrows (and a fat helix) with a continuous backbone "spine" object (`<obj>_spine`, PyMOL `cartoon tube`). Because the tube ignores secondary structure, it runs unbroken through every strand→loop junction; exported together with the chunky cartoon via `print_export(representation="cartoon", method="voxel", voxel_pitch=0.2)` the voxel step fuses them into one watertight solid with no junction discontinuity, and the spine doubles as internal rebar for print rigidity. Configurable `spine_radius`.
- `print_export` `representation` parameter. `"surface"` (default) is the existing behaviour, unchanged. `"cartoon"` exports the *currently displayed* cartoon geometry of the real objects — preserving per-residue rep flags (hidden loops) and per-object cartoon type (the `cartoon tube` spine) — instead of recreating a temp object and forcing a molecular surface. Groups are isolated by toggling object visibility (one colour per object).

### Fixed
- `print_export` always exported the molecular **surface**, even when the scene was set up as a cartoon (e.g. by `print_ribbon_view`): it hid all reps on a throwaway temp object and forced `show surface`. The new `representation="cartoon"` path exports the actual displayed ribbon/tube geometry, so `print_ribbon_view` now produces a ribbon STL rather than a surface blob.
- `_repair_to_stl` voxel method produced fragmented, **non-watertight** output (e.g. a cartoon export came out as 19 loose shells) because it never consolidated the marching-cubes result. It now keeps the largest body and fills holes — the same consolidation the `light` path already does — yielding one watertight, printable solid. Mesh volume is unchanged (no inflation of fine features).

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
