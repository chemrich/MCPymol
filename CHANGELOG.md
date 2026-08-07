# Changelog

All notable changes to MCPymol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **`CONTRIBUTING.md` recommended a test pattern that does not work.** It told contributors to mock `mcpymol.server.send_request`; after the v1.3.0 package split that rebinds only the facade's name, leaving every module still calling the real socket — the exact trap that broke 74 tests during the split. It now shows the correct target and explains why.
- Docs still pointed at `server.py` for the bridge protocol and for "after editing", both of which moved when the package was split.
- README claimed 340+ tests (actually 469) and ~60 primitives (actually ~85, of 120 tools).
- `load_structure` was documented nowhere, so opening a local PDB/mmCIF — a built model, a docking pose, an MD frame — was undiscoverable. `fetch_alphafold` and `count_atoms` were likewise unnamed.

### Added
- README figures for `conservation_view` (lysozyme against 5,507 homologues), `print_ribbon_view` (GFP, spine visible through the strands), `textbook_view` (haemoglobin), `cinematic_view` (GroEL down its seven-fold axis), `plddt_view` (AlphaFold spike) and `superposition_view` (adenylate kinase), plus an animated `turntable` GIF. Every documented view now carries an image except `bfactor_view` and `pointillist_view`, which are deliberate: `putty_view` already illustrates the B-factor palette, and the pointillist style is decorative.
- Real output for `contact_report`, `interface_report` and `structure_info` in the README. These were the headline features of the release and it showed none of what they return. The contact and interface examples use the same systems the `ligand_view` and `interface_view` figures show, so the picture and the numbers sit side by side.
- All figures were produced through the tools themselves, which is how three of the bugs in this release were found.
- **Configuration section in the README.** Six of nine `MCPYMOL_*` environment variables were undocumented, including every timeout introduced in v1.3.0. All nine now have a table entry saying what they do and what the default is.
- Troubleshooting rows for the failure modes introduced since v1.2.1: files not appearing when PyMOL and the bridge are on different machines, oversized renders, AlphaFold entries with no model at a given version, blank views from empty selections, and slow-operation timeouts.
- `tests/test_docs.py` keeps these true. An environment variable added in code but not the README, a headline tool missing from the README or the skill, a doc that sends readers to `server.py` for implementation, or a `CONTRIBUTING` that stops warning about facade patching — each now fails CI. Written after this audit found six undocumented variables and three unmentioned tools by hand; the point is not to do that by hand again.

## [Unreleased]

### Fixed
- **`render(ray_trace=False)` and `turntable`'s default produced blank images.** PyMOL's unshaded OpenGL capture needs its GUI thread, and the plugin dispatches on a socket worker thread, so `cmd.png(..., ray=0)` writes a flat empty frame rather than failing. `turntable` defaulted to exactly that path, so every frame of every animation came out blank. Both tools now always ray-trace and say why if asked not to. `render` additionally detects a flat single-colour result and points at `list_objects`/`count_atoms`, which catches the other way a render comes back empty: a scene where nothing is visible.
- **Every view preset rendered with a transparent background.** They set `bg_color` but never `opaque_background`, so ray-traced output kept an alpha channel: correct in the viewport, transparent in the PNG, which appears white wherever the image is used. 16 of 17 presets were affected, `textbook_view` included (its background is white rather than black, and had the identical bug). Backgrounds now go through `style.set_background()`, and a test fails if any module reaches for `bg_color` directly.
- **`superposition_view` was unreachable through the documented workflow.** `fetch_structure` calls `reinitialize`, so fetching a second structure wiped the first — leaving nothing to compare, for a tool whose only purpose is comparing two structures. All three loaders (`fetch_structure`, `load_structure`, `fetch_alphafold`) take `replace=False` to add to the session instead of clearing it.
- **`fetch_alphafold` did not work for any accession.** It built the model URL from a hardcoded `AF-{accession}-F1-model_v4.cif` template. AlphaFold DB has since moved to v6 and *removes* retired versions, so every one of those URLs now 404s — the feature shipped in v1.3.0 and resolved for nothing. It also assumed the filename is keyed by accession, which is untrue for some entries: SARS-CoV-2 spike (P0DTC2) is served as `AF-0000000365840314-model_v1.cif`, so no accession-based filename exists for it in any version.

  The file URL is now asked for rather than constructed, via AlphaFold DB's prediction API, which handles version drift, non-accession entry IDs and multi-fragment proteins. `model_version` becomes an optional pin that bypasses the lookup, and defaults to unset.

  Every test for this code mocked `urlopen`, so they proved the URL was *built* correctly and never that it *existed*. A `network`-marked test now resolves and fetches a real model; it is deselected by default so CI stays offline-safe (`pytest -m network` to run it, worth doing before a release).
- A malformed accession returned a raw `HTTP 400 Bad Request`; it now says which accession failed and what one looks like.

## [1.4.0] - 2026-08-07

MCPymol could make pictures but could not answer questions with numbers. This
release closes that gap.

### Added
- **`contact_report`** — lists contacting residue pairs across two selections, closest first, with the minimum heavy-atom distance, the number of atom contacts, and a classification: salt bridge, hydrogen bond, hydrophobic, polar contact, or π-stacking (parallel vs T-shaped by interplanar angle). Criteria are heavy-atom distances, documented and pinned by tests, since crystal structures usually have no hydrogens — so a reported hydrogen bond is a donor–acceptor pair with plausible geometry, not one verified against a hydrogen position. Neighbour search uses a uniform grid, keeping cost linear in atom count; a property test asserts it finds exactly what a brute-force scan does. Ring perception needs bond orders that a PDB dump does not carry, so aromatics cover the standard aromatic amino acids only and ligand rings are reported as hydrophobic rather than mis-called as stacking.
- **`interface_report`** — buried surface area from ΔSASA (free minus bound, halved for the per-side figure papers quote), residues ranked by how much surface each buries, and a breakdown by residue chemistry. Interprets the number against PDB-wide survey thresholds: under ~400 Å² per side is usually crystal packing, over ~1000 Å² a substantial and likely specific association. Uses `get_area(load_b=1)`, which writes per-atom SASA into the B-factor column, so the whole per-residue breakdown costs one area calculation and one dump rather than hundreds of round trips — and sets `dot_solvent=1` first, without which the numbers are wrong but plausible-looking.
- **`structure_info`** — what a structure actually is, in one call: chains, atom/residue/water counts, ligands, states and space group from PyMOL, plus title, method, resolution, release date and source organism from the RCSB data API. Best-effort metadata: skipped for objects not named after a PDB entry, degrading to "no metadata" when the network is unavailable. Flags a probable predicted model when the B-factor column looks like pLDDT.
- **`get_sequence`** — there was previously no way to get a sequence out of MCPymol at all. Returns FASTA plus the two things the sequence hides and that routinely cause silent mistakes: the numbering offset (PDB numbering rarely starts at 1, so "residue 50" in a paper and position 50 in the sequence are usually different residues) and chain breaks where loops went unmodelled.
- **`sasa`**, **`rms_cur`**, **`count_atoms`** — value-returning primitives. `count_atoms` in particular catches an empty selection, which is otherwise invisible until the render comes out blank.
- New `mcpymol.pdbtext` module for parsing the PDB dumps that are the only way to read data back out of PyMOL, and `mcpymol.analysis` for the reporting tools.

### Fixed
- **Measurements discarded their results.** `cmd.distance`, `angle` and `dihedral` all return the quantity they measure; the wrappers threw it away and returned "Measured distance between 'X' and 'Y' as 'd1'" or "Executed angle successfully." So you could ask MCPymol to measure something and it would draw the measurement without telling you the number. All three now report the value. PyMOL's -1 "nothing matched" sentinel is reported as a failure rather than as a measurement of -1 — but only for quantities that cannot be negative, since -57.8° is an ordinary α-helical phi.

### Changed
- **Every tool parameter is now documented in the JSON schema.** All 222 parameters across all 120 tools carry `Annotated[..., Field(description=...)]`, so `inputSchema.properties[].description` is populated and the model no longer infers arguments from parameter names — it can see what `spectrum`'s `expression` accepts, what `clip`'s `mode` means, and what syntax `mset` takes. The parameter docs previously lived only in `Args:` docstring blocks, which FastMCP does not read into the schema; the signature is now the single source of truth and the redundant blocks are gone. Verified by diffing every schema across the change: descriptions added, types/defaults/`required` untouched. Four new tests keep it that way — a tool shipping an undocumented, empty or one-word parameter description now fails CI.
- README reorganised around the question being asked rather than the tool being called, with a "What can I ask?" index and a new Analysis section. The mcpymol-guide skill leads with report-first-then-draw.
- `views._read_ca_bfactors` and `comparison._parse_ca_coords` now delegate to `pdbtext.parse_atoms`, which additionally handles missing element columns, insertion codes and truncated lines. The two `_distance` implementations are unified on `pdbtext.distance3d`.

### Note on dependencies
The optional `analysis` extra was considered and deliberately not added. For this scope nothing earned it: per-residue SASA comes from PyMOL in two round trips via `get_area(load_b=1)`, and contact search is a spatial-grid problem where numpy would add a second code path for no measurable gain. The core stays dependency-free.

### Developer experience
- **Pre-commit hooks that actually cover the gates.** The config had a single hook (the `uv.lock` guard), `pre-commit` was not a declared dependency so nothing installed the binary `CONTRIBUTING.md` told you to run, and the hook was consequently never installed in practice — every check happened only in CI, after a push. Commits now run ruff, ruff-format, mypy, the lockfile guard, and file hygiene (trailing whitespace, final newline, line endings, YAML/TOML validity, merge markers, stray `breakpoint()`, oversized files); pushes run the full test suite, which is too slow at ~16 s to pay per commit. `pre-commit` is in the dev group, and `default_install_hook_types` wires both stages from one `pre-commit install`.
- The ruff and mypy hooks run via `uv run` rather than pre-commit's own tool repos, so they use exactly the versions in `uv.lock` — the versions CI installs — instead of introducing a second place for tool versions to drift.

### Tests
- 346 → 469.

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
