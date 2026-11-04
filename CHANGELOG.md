# Changelog

All notable changes to MCPymol will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- `print_export` no longer fails a valid group because a warning PyMOL printed
  while redrawing landed in its selection check; only a `Selector-Error`
  counts as a rejected selection.
- **A PyMOL whose Python is too old now says so.** The plugin runs in PyMOL's
  own Python and needs 3.10 or newer. On 3.9 (older Incentive builds) it
  loaded, then its listener ended at the first idle poll, so PyMOL looked
  ready and never answered. It now refuses at startup and says why, and the
  `.pymolrc.py` block `--install-plugin` writes checks first, for Pythons
  where the plugin file wouldn't even parse. Re-run `mcpymol --install-plugin`
  to get the new block. The README says older Incentive 2.x builds aren't
  supported.
- **`--install-plugin` writes a block that survives any path.** The plugin's
  path was pasted into a Python string as is, so a Windows path
  (`C:\Users\...`) or one containing a quote was a SyntaxError that stopped
  PyMOL reading the rest of `.pymolrc.py`, the user's own settings included.
- **Tools leave your settings and folders as they found them.**
  `interface_report` set `dot_solvent` and `dot_density` globally and left
  them, which also changed how every dot representation in the session was
  drawn; it now puts them back. `fetch_structure` wrote each downloaded file
  into PyMOL's working directory (wherever PyMOL was launched); it now uses a
  private per-user cache directory, unless you have set PyMOL's
  `fetch_path` yourself.
- **`list_ligands` lists five-character ligand names whole.** It read residue
  names from PDB text columns, which hold three characters, so `A1AAA` was
  listed as `A1A`; the CCD has issued five-character codes since the
  three-character ones ran out.
- **A name can no longer delete what it doesn't name.** The loaders delete the
  object they replace by name, and PyMOL's `delete` takes patterns and lists:
  `load_structure`, `fetch_structure` or `fetch_alphafold` given `*`, `all` or
  `keep*` as the object name deleted every other object and reported success,
  even with `replace=False`, and `pep two` loaded as `pep_two` while the reply
  said nothing had loaded. A name PyMOL would change (its own
  `get_legal_name`) is now refused before anything is deleted. The same class
  was in the helper objects presets and wiggles name after the caller's object:
  `crosslink_view(obj_name="1ubq and chain A")` ran `delete 1ubq and chain
  A_disulfides`, which deleted `1ubq` itself. `crosslink_view`, `pocket_view`,
  `pharmacophore_view`, `poisson_boltzmann_view`, `deformation_view`,
  `morph_states`, `load_ensemble` and `latent_traverse_view` now always name
  helpers as one plain name (`1ubq_and_chain_A_disulfides`); ordinary names
  are unchanged. `load_map` given `all` loaded it as `all_` and then, finding
  no `all`, cleaned up with `delete all`; it now refuses such a name too.
- `fetch_alphafold` refuses an identifier that isn't an accession; one with a
  `/` in it raised a bare `FileNotFoundError`.
- **`mcpymol --install-plugin` and the other plugin flags no longer need a
  working MCP stack.** The console script pointed at `mcpymol.server`, which
  imports every tool and mcp itself, so when an mcp release broke that import
  (as with 1.5.0), finding, installing or removing the plugin broke too. The
  script and `python -m mcpymol` now start in the CLI, which imports the
  server only to serve.
- **`conservation_view` stops waiting when the answer won't change.** Its
  poll of the MMseqs2 server skipped every failed status check and kept going
  on any state but COMPLETE or ERROR, for up to ten minutes while every other
  tool queued behind it. A job the server has forgotten (it answers UNKNOWN
  once a job expires) now ends the wait after 30 s, a refused status check
  or an unexpected state at once, and a rate-limited submission says so
  instead of a bare "HTTP Error 429"; network blips are tolerated six times
  in a row. The User-Agent carries the real version, and
  the README says plainly that the chain's sequence goes to the server.
- **`list_objects` lists every object, not only molecules.** It used
  `get_object_list`, which PyMOL backs with the selector, so a session holding
  just a map read "No objects are loaded", and a structure's maps,
  measurements and groups went unlisted. It now names them all, a
  non-molecule with its type (`emd_1234 (map)`); molecules read as before.
  `load_session`'s summary of what the session now holds had the same gap.
- **`pytest -m live` no longer drives your own PyMOL.** The live suites clear
  the session they reach, and with nothing set that was port 9876, a
  developer's everyday PyMOL. They now run only against `MCPYMOL_TEST_PORT`, a
  port named for testing, and refuse to start without one. A run that selects
  live tests and passes none of them now fails instead of reporting a row of
  skips, and `MCPYMOL_LIVE_MIN_PASSED` sets the bar for a release gate.

## [1.6.0] - 2026-11-04

The release that drains a long Unreleased section: the cryo-EM tools in
`mcpymol.wiggles`, a rebuilt 3D-print export, and a sweep of defects found by
checking the tools against a real PyMOL rather than mocks.

**Two breaking changes.**
- **Python 3.11 or newer.** Python 3.10 reached end of life in October 2026.
  On 3.10, pip and uv keep resolving to 1.5.2, which works.
- **mcp 2.x.** The tools a client sees are unchanged, byte for byte: the
  `tools/list` JSON for all 135 tools was diffed before and after.

**Upgrading.** `uv tool upgrade mcpymol`, then re-run `mcpymol --install-plugin`
so PyMOL loads the matching plugin; the bridge and the plugin should be the
same version. On Linux, the `print` extra now needs the system OpenGL library
(`sudo apt-get install libopengl0`).

### Added

- **`qscore_view` reads mmCIF validation reports**, as its schema already claimed. `*_validation.cif(.gz)` files keep Q-scores in `pdbx_vrpt_model_instance_map_fitting`, keyed by an instance id that `pdbx_vrpt_model_instance` resolves to a residue. The reader joins the two and keys the result by author chain and number, exactly as the XML reader does: on 10EC, 10LJ and 10ND all 1262 residues read identically in both formats. The format is decided by content, not extension, because downloaded reports get renamed. The mmCIF reader uses only the standard library, with a small tokenizer for the STAR grammar these files use (quoted values, `;` text fields, comments).
- **`MCPYMOL_AUTOSTART=0` loads the plugin without opening its port** (issue #15). Type `start_mcp` in PyMOL to start listening, and `stop_mcp` to stop. The default stays on: every install's `~/.pymolrc.py` runs the plugin and expects PyMOL to listen, and quietly not doing so after an upgrade would look like a broken plugin.
- **A new `Fresh install` workflow installs the built wheel without the lockfile and imports it.** It runs on every PR and every Monday, because an upstream release can break a fresh install with no commit here. It covers the ends of the supported range, 3.11 and 3.14.
- **`print_export` prints sticks, spheres and distances** (`representation="displayed"`). PyMOL's OBJ writer emits triangles only, so everything drawn as a cylinder or sphere was dropped: a sticks-only group crashed with *"zero-size array to reduction operation minimum which has no identity"*, and a spheres-only group exported 0.2 Å³ of degenerate triangles and reported OK. Displayed mode reads the cylinders and spheres from PyMOL's POV-Ray export and merges them with the triangles by an exact boolean union (`manifold3d`, new in the `print` extra). Half-bonds are fused into one rod, since abutting rods only touch face to face and split the STL into hundreds of open pieces, and distance dashes print as solid rods. 1hsg's ligand as sticks exports as one watertight piece within 0.003 Å of atoms ± stick radius. Cartoon mode stays triangle-only and now names what it left out; a group with nothing printable says what it does show instead of crashing. Dots and mesh are refused in displayed mode: PyMOL draws them as a sphere or rod per point, and 1ubq's dots came to 38,567 pellets, a 1 GB STL after ten minutes.
- **`print_export(mm_per_angstrom=…)`**, the print scale. It is applied to every group about the origin after repair, so parts stay registered while `voxel_pitch` stays in Å. The result now gives each part's size in mm, its volume and the assembled size, and warns when a part is small enough for a slicer to take it for inches and enlarge it 25.4× on its own (under 9 mm³ in PrusaSlicer 2.9.6, 8 in Bambu Studio 2.8 and OrcaSlicer 2.4; Bambu's 2.8.4 beta does it without asking), with the scale that avoids it.
- **`print_export(bodies=…)`** and **`print_export(undo_voxel_growth=…)`** — see Changed.
- **`print_export` warns about what will not print**: a cartoon feature under 3 voxels across (naming the setting and object), rods and spheres under 0.8 mm at the chosen scale (naming what draws them, e.g. valence rods), and shown residues whose Cα lies outside the printed mesh — the gap at Leu64 beside GFP's hidden chromophore that went unreported while a lesson model was being built.
- **`print_ribbon_view(strand_width, strand_thickness, spine_in_strands, strand_overlap, bridge_modified)`**. PyMOL draws a strand √2 times its `cartoon_rect_*` settings, so the old fixed strands were 4.0 × 2.8 Å and a β-barrel's wall closed up; `strand_width=2.0, strand_thickness=1.0, spine_in_strands=False` takes GFP's barrel from 8% to 29% see-through and removes the spine bulging out of thin strands (399 Å³ to 1 Å³), still as one piece.
- **`execute_pymol_command` returns what PyMOL printed**, and when PyMOL rejects a line, which line and PyMOL's own message — with a note when `alter`, `iterate` or `label` swallowed a `;` meant to start the next command. It waits for the command to finish (`MCPYMOL_SLOW_OP_TIMEOUT`). Plugin knobs `MCPYMOL_MAX_OUTPUT_CHARS` and `MCPYMOL_DO_START_TIMEOUT`.
- **macOS: the plugin holds off App Nap while it is serving requests** (`MCPYMOL_APP_NAP_IDLE_S`, default 900 s after the last request; 0 turns it off). A backgrounded PyMOL drops every thread to the lowest priority, and under load a 1 s ray trace was measured at 74–349 s. It uses the activity that still allows idle sleep, so it never keeps the Mac awake.
- **`mcpymol.wiggles` tier 3: ensembles.** `load_ensemble`, `latent_traverse_view`, `deformation_view`, `composition_view`. One blind-challenge result governs all four: put on the same thyroglobulin data, 41 submissions recovered molecular motions that resembled both each other and the ground truth, but the simulated truth had three population modes, most found two, and only three of 41 found all three — the missed one being the middle. So these tools draw motion and decline to draw populations. `deformation_view` is built with confidence; `latent_traverse_view` renders no latent scatter and estimates no density at all.
- **Two invariants, enforced by tests rather than by convention.** **I2**: no latent rendering goes out unlabelled, so `latent_traverse_view` *refuses* when `load_ensemble` could not identify the generating method — cryoDRGN's latent density can bear no relation to the truth while 3DVA's frames are linear interpolations no particle occupied, and drawing them identically asserts something for one that holds only for the other. **I3**: no output may turn a gap into an absence; a phrase list is checked against every tier-3 report, refusals included.
- **A per-frame normalisation trap, closed.** PyMOL normalises each map independently on load, so contouring a trajectory at a fixed sigma contours every frame at *its own* sigma — and a traversal is exactly the case where density genuinely changes between frames. That is the signal, and per-frame normalisation rescales it away, leaving frames that look reassuringly alike whatever the ensemble does. The level is now held in absolute terms and converted per frame, and the spread of the resulting sigmas is reported, because that spread is the change.
- **Method detection that refuses to guess.** Heterogeneity jobs are identified from markers their own documentation describes. A directory matching nothing loads its volumes and stays `UNKNOWN`, which is what makes I2 bite. Latent tables are read from text, STAR and plain `.npy`; a `.pkl` is found but **not** read unless `trust_pickle=True`, because unpickling runs arbitrary code from a file that may have come from a collaborator or a shared scratch volume.

### Changed

- **MCPymol runs on mcp 2.x** (`mcp>=2.3.0,<3`; 1.x pinned it below 2). What a client sees is unchanged, and checked: the `tools/list` wire JSON for all 135 tools is identical before and after, and so are call results apart from two fields the 2.x protocol adds. Three 2.x behaviours needed handling. **Synchronous tools now run on worker threads**, so tool calls a client issues in parallel would have interleaved inside the one PyMOL session; 1.x serialised them by accident, and `app.SESSION_LOCK` now does it on purpose, one tool call at a time, start to finish. **An unexpected exception's text no longer reached the client**, only "Error executing tool X"; the server passes it on as 1.x did, because for a local tool that text is how the model learns what went wrong. **An unversioned server reports an empty version**, so the server now reports mcpymol's own.
- **Python 3.11 is now the minimum, and 3.14 is supported.** Python 3.10 reached end of life in October 2026, and its dependencies can no longer be kept clean: pillow's current security fixes are in pillow 13, which requires 3.11. On 3.10, pip and uv keep resolving to 1.5.2, the last release that supports it, which works. This is the bridge's Python only. `plugin.py` runs inside PyMOL's own interpreter, which MCPymol never imports into, and its syntax stays 3.10-compatible (pinned with ruff's `per-file-target-version`).
- **STLs are written in PyMOL model coordinates whatever the camera is doing.** PyMOL's OBJ writer emits camera space, so turning the view between two exports swapped x and z, and "all groups share one frame" held only within one call. The export now saves under an identity view and restores yours exactly; separate calls at the same scale line up.
- **`print_export` keeps every separate piece by default** (`bodies="all"`) and says *"prints as N separate pieces"*; `bodies="largest"` keeps one and warns about what it dropped. Enclosed cavity shells are filled and stray fragments removed in every method, and counted.
- **Voxel exports are no longer fat** (`undo_voxel_growth=True`). Voxelising grew every surface by about 0.62 voxel — 0.43 Å at the default pitch, +68% volume on a 4hhb print ribbon — and each re-voxelisation grew it again. Surfaces are moved back in, never past half a thin part's thickness, so voxel parts come out about that much thinner than before. Gaps under about 2.5 voxels still fuse while voxelising; the result says so.
- **Poisson reconstruction runs in a child process**, and intermediate OBJ files no longer land in `out_dir` beside the STLs. When `auto` falls back from Poisson to voxel, the result says why.
- **`print_export` gives the session back group by group.** Each group's objects are shown alone, and rock and stereo paused, only while that group is written out; all three are back before its repair, which can take minutes. A server killed mid-repair used to leave the user's objects hidden and temporary objects in the list.
- **`print_ribbon_view` styles only its object and spine.** It used to change the global cartoon settings, so every other object in the session, and every one loaded later, took the chunky strands. Output at the defaults is byte-identical.
- **`execute_pymol_command` and every internal `do` are synchronous.** `cmd.do` from the plugin's socket thread only queues the line for PyMOL's main thread and returns within about 0.13 s, so the old reply could arrive before the command ran and could never carry its error. Lines still run on the main thread (some, like `window size`, kill the GUI PyMOL anywhere else), through a private parser whose verdict is kept.
- **Tools carry PyMOL's warnings.** Problems PyMOL only prints — `Selector-Error` from `zoom`, `Setting-Error: unknown color` — used to reach the PyMOL console alone while the tool said it succeeded; they now end the tool's result as *"PyMOL warned: …"*.
- **`print_export` counts what each object shows in one request**, reading each atom's `reps` mask, instead of thirteen `count_atoms` round trips per object before every export.
- **`execute_pymol_command` formats PyMOL's warnings like every other tool**: a single warning on the same line, several as a list.
- **CI tests what it used to assume.** A `Minimum versions (py3.11)` job installs every dependency at the floor `pyproject.toml` declares and runs the unit suite. Its first run also exposed a latent break: pytest 9.2 undoes `monkeypatch.setattr` on a `MagicMock` child by deleting it, which left the shared `pymol` mock missing `get_setting_int` for every later test, so the next Dependabot bump would have turned 13 tests red. Those patches now use `mocker.patch.object`. The print extra's real-mesh tests run in CI for the first time (they skipped wherever trimesh was absent, which was everywhere). The dependency audit covers the extras: pillow alone carried 52 open advisories that it never saw. Every function is now annotated, mypy runs with `disallow_untyped_defs` and `warn_unused_ignores`, and three stale `# type: ignore`s are gone. Live test runs no longer leave `.cif` files and `_live_sym*.png` frames in the checkout: PyMOL's working directory points at a temp dir for the run. `.coverage` is no longer tracked.
- **The live suite now catches two classes it was blind to.** `WIRING_FAILURES` listed only Python binding errors, so "PyMOL accepted the signature and rejected what the argument meant" — `Object not found`, `Selector-Error`, `Invalid selection` — passed silently; those signatures are now included. Separately, the sweep supplies most arguments from its own table and so never calls a tool the way its schema says it can be called, which is why `spheroid` passed every run despite an unusable default. `test_defaults_are_usable` calls every fully-defaulted tool with no arguments; it found the `frame` defect above on its first run.
- **`tests/test_docs.py` walks the package recursively.** The glob stopped at the top level, so tools defined in a subpackage would have escaped the checks that every tool is named in the README and in the agent skill — the exact rot this file exists to prevent, and it would have let ten of them through silently.
- **`tests/test_package_layout.py` covers subpackages.** `mcpymol.wiggles` and `mcpymol.wiggles.tools` are both listed: the first is what the module walk sees, the second is where the tools are defined and so the one the re-export checks have to look at.

### Fixed

- **A structure with no chain IDs made selections escape their object, outside wiggles too.** `get_chains` returns `['']` for such a file, and a raw `chain {chain}` then reads `chain  and polymer.protein`, where PyMOL takes `and` as the chain name. Checked on PyMOL 3.1.0: it matched 82 atoms where the object held 40. Loading colours each chain, so loading such a file with `replace=False` recoloured every other object's protein (reproduced live). The multimer heuristic built its neighbour query and its `remove` the same way. About a dozen sites in `structures`, `analysis`, `comparison`, `conservation` and `views` now go through `mcpymol.selection.chain_clause`, which quotes each chain. `quote` moved there from `mcpymol.wiggles.atoms`, where this bug was first found and fixed. A user's `A+B` still means two chains. In `superposition_view`, a blank chain no longer widens to the whole object, which wrote that chain's deviations onto every chain with the same residue numbers.
- **Four presets sent a command PyMOL rejected, and reported success anyway.** `pharmacophore_view` never coloured ring carbons: `inring` is not a selection keyword, and the ring operator is `byring`. `cinematic_view` set a `fog_color` that does not exist; fog takes the background colour. `pointillist_view`'s grey ribbon was a malformed selection, and is now the `ribbon_color` setting, because colouring the atoms would recolour the dots. Loading set the nucleic-acid cartoon colours on an atom selection, where PyMOL stores them and never draws them; they are now set on the object.
- **A preset given a name that matches nothing reports an error** instead of styling an empty selection and claiming success, as "Showing ligand view ... H-bonds stored" did over a ligand that was not there. Every view preset checks its object, and the ligand, pocket, pharmacophore and interface presets also check their ligand or chains. Only a definite answer refuses (zero atoms, or PyMOL rejecting the selection); if the count can't be read, the preset runs as before.
- **`tests/test_live_presets.py` runs every preset against a real PyMOL** and records each response, failing on anything PyMOL rejects or warns about.
- **`local_resolution_view` let two maps with different cell angles through its grid check.** It compared extent, sampling, voxel size, start, origin and axis order, but never `cellb`, so a map pair whose grid axes are not parallel passed the check the module exists to enforce. Cell angles are now compared too.
- **`map_info` called a transposed map cropped.** `nx/ny/nz` count columns, rows and sections in file order, while `mx/my/mz` run along X, Y and Z. They were compared as written, so any map with a non-default axis order whose extents are not all equal looked boxed or cropped. The extents are now mapped through `mapc/mapr/maps` first.
- **An unstubbed `get` in a wiggles test silently chose an answer.** `FakePort` answers a data action it has no response for by raising, but `get` was missing from that list, so it answered "OK" and `density_view` took its "normalisation on" branch without the test saying so. `get` is now a data action, and the density tests state PyMOL's default explicitly.
- **Every wiggles tool is checked to forward every argument it accepts.** `tests/wiggles/test_tools.py` gives each parameter a distinctive value and asserts it reaches the view. A parameter the schema advertises but the wrapper drops is otherwise invisible: the view runs on its default and reports success.
- **`align`, `super`, `intra_rms` and `intra_fit` report their numbers** instead of "Executed X successfully." (issue #15). `align` and `super` give the RMSD after refinement, the atoms in it, the cycles, the RMSD before refinement, and the residues aligned. `intra_rms` and `intra_fit` give one RMSD per state against the reference state, and `intra_fit` adds that it moved coordinates in the session. An answer of any other shape falls back to the old wording rather than being misreported.
- **Replacing a structure left the old one's maps, meshes, ramps, measurements and groups behind.** `fetch_structure`, `load_structure` and `fetch_alphafold` with `replace=True` promise "the only structure loaded", but cleared the rest by listing `get_object_list`, which returns molecules only. They now list every object. Checked live: with the old call, a map, ramp, distance and group all survived.
- **The plugin answered success for a fetch that downloaded nothing.** `cmd.fetch` returns -1 rather than raising when a download fails. The loaders' own atom count caught it, but any other caller saw "Fetched". It is now an error.
- **An internal `do` could time out and run behind the next request.** `do` now replies once the command has finished rather than when it was queued, and about twenty internal callers still had the 10 s default, so a view loading a large APBS map or applying a long batched `alter` timed out, moved on, and built on data not there yet. `do` now defaults to the slow-operation timeout.
- **`frame()` with no argument never said which frame.** The query answered "Executed get_frame successfully." It now returns the frame number.
- **`print_export` dropped what it had already written when a later group failed.** It returned only that group's error, while the earlier STLs sat in `out_dir` unreported. It now names them, and says the set is incomplete.
- **`print_ribbon_view` did not say a failed run had left the object partly restyled.** It is not rolled back, but every run starts from scratch, so the error now says a re-run puts it right.
- **Code run inside a captured request could not reach the real stream.** The plugin's stand-in for `sys.stdout` was a bare `StringIO`, so `fileno()` raised and `encoding` was `None`, which broke `subprocess.run(stdout=sys.stdout)`, `faulthandler` and some C extensions inside a `do`.
- **`set` reported success for a value PyMOL never uses.** An object-level setting such as `cartoon_tube_radius` given a residue selection is stored per atom and ignored, and PyMOL only warns; `set` now returns an error saying what to do instead.
- **`alter` or `iterate` with a Python error in its expression returned *"Empty response from PyMOL plugin."*** Under Python 3.14, PyMOL's traceback printer raises while formatting the NameError, and that escaped the plugin's request handler.
- **A mistyped selection could crash PyMOL** — the whole process and the user's session with it. On PyMOL 3.1.0, `rms`, `rms_cur` and `fit` abort (SIGABRT in ExecutiveRMS) and `intra_rms` and `intra_rms_cur` segfault (SIGSEGV in CmdIntraFit) on a selection PyMOL cannot resolve; `intra_fit` and `pair_fit` report the same selection as a plain error. Each case was checked in a fresh PyMOL process. The plugin now vets these five in its dispatcher, the one place every request passes through, so the tools, wiggles' port and any future caller are covered. An empty selection is refused too: it does not crash, but `rms_cur` over nothing answered 1.4e-45. A raw command string through `execute_pymol_command` is not vetted (Known issues §7).
- **A failed Poisson reconstruction ended the MCP server mid-call.** PoissonRecon, run in-process through pymeshlab, calls `exit()` when it gives up — 4 of 14 runs on one 1ubq ribbon — and `method="auto"` uses Poisson for any open cartoon. It now runs in a child process, and `auto` falls back to voxel.
- **A `print_export` group label with a space could delete your objects.** Surface mode named its temporary objects after the labels, and `delete` reads a name with a space as a list: `chain A=all` deleted an object called `A`, and a label ending in the exported object's own name deleted that object. Temporary objects are now named by position, and labels, which also name the STL files, are limited to letters, digits and `_ . + -`. Two labels that differ only in case are refused too: macOS file names ignore case, so the second group's STL silently replaced the first while the result listed both.
- **`print_export` re-enabled objects you had disabled** (`enable all` at the end), exported nothing for an object inside a disabled PyMOL group, re-zoomed your camera in surface mode, and — when PyMOL's object list came back empty — enabled unrelated objects by name prefix and exported one of them as OK. A distance or CGO object named in a group made PyMOL reject the whole selection.
- **`print_ribbon_view` reported success for an object that does not exist**, and its spine broke at residues with no CA atom: GFP's chromophore leaves Cα64 and Cα68 nine bonds apart, so a thin spine came out in three pieces and the old keep-largest step dropped about 30% of a thin-strand barrel. The spine is now traced through such residues on the spine copy only; nothing in your object is renamed. A calcium ion, also named CA, does not count as a traced Cα, so calcium sensors (GCaMP) are bridged too.
- **A timeout blamed the plugin.** A request PyMOL accepted and has not answered now says it may still be running and was not cancelled, rather than *"Is the PyMOL plugin running?"*; when the computer slept during a call that timed out, the error says for how long, and so does a render that finished. The two renders reported as "stalling until the next command" had the Mac asleep for 579 s and 716 s of their wait.
- **Unit tests can no longer reach a real PyMOL.** Nothing stopped a test that slipped past its mock from driving whatever PyMOL was listening on the default port — the session a developer has open. `tests/conftest.py` now refuses those connections outside tests marked `live`.
- **The live suite was wiping its own map fixtures before any test ran.** `_require_a_healthy_viewport()` calls `reinitialize`, and the volumes were loaded *before* it. Because map headers live in a process-global registry, a deleted object left the registry intact: `density_view` and `local_resolution_view` proceeded to contour something that no longer existed, failed with an honest domain error, passed the wiring check, and proved nothing. Both now run against volumes that are really there, and the fixture asserts its own setup rather than trusting wrappers that turn failures into strings. `map_info` and `load_map` were unaffected — they take a path and load fresh.
- **`count_atoms` joins `DATA_ACTIONS`.** `composition_view` consumes its return value, so a test that forgot to stub it would have seen `"OK"`, read zero atoms, and exercised the empty-selection refusal instead of the tool.

- **`mcpymol.wiggles`: ten cryo-EM tools for occupancy, ensembles, maps and local resolution.** `occupancy_view`, `altloc_view`, `ensemble_spread_view`, `morph_states`, `qscore_view`, `restore_bfactors`, `map_info`, `load_map`, `density_view`, `local_resolution_view`. Each exists because a render that looked fine was hiding something, and they share one rule: say what is being shown, and refuse when the picture would be meaningful-looking and wrong. `morph_states` declines to interpolate states that do not share a topology, because a morph across independently reconstructed volumes animates a correspondence nobody established. `local_resolution_view` declines to colour one map by another that does not share its voxel grid, because sampling the resolution field at the wrong coordinates renders smooth, plausible and wrong rather than visibly broken. `load_map` records provenance and never infers it — a measured reconstruction and a network-enhanced one are the same isosurface once drawn, so defaulting to "measured" would assert that somebody observed a generated volume.
- **Two unit traps closed, both of the kind that yields a wrong answer of the right order of magnitude.** `map_info` computes voxel size as `cella/m`, not `cella/n` — the two differ on any boxed or cropped map — and the nominal figure it reports is itself only accurate to ±5–15%, which at 1.2 Å is a systematic stretch of every distance in the model. `density_view` states its contour level in both σ and absolute units, because PyMOL contours in σ while EMDB publishes absolute levels: EMD-30913's published `0.05` is 3.16 σ, and used directly it contours noise. `local_resolution_view` hits the same trap on a second axis — its ramp breakpoints are converted to σ against the *resolution* map's header, a different scale from the contour level's.
- **The live suite now builds its own cryo-EM fixtures.** A 16 KB synthesised MRC volume and a minimal wwPDB validation report, both written by the session fixture, so all ten tools are swept against a real PyMOL without a network round trip or a checked-in binary. Nothing was added to `EXCLUDED`.
- **A residue with no chain ID made three views overwrite every atom in the session.** `qscore_view` defaults a missing chain to `""`, and the selection was built by interpolation, so `chain {chain} and resi {resi}` became `chain  and resi 2` — where PyMOL takes `and` as the chain name and the selection stops being scoped to the object. Checked on PyMOL 3.1.0: in a session holding a 7-atom object and a 10-atom one, that selection matched 10. `qscore_view` on a validation file whose entries carry no chain attribute therefore ran `alter <every atom loaded>, b=<one residue's Q-score>`, and since `restore_bfactors` only ever restores the object it is given, every other structure's B-factors were gone with no way back. `ensemble_spread_view` built the same string.
- **A negative residue number was read as a range.** `resi -3` means 1-3 to PyMOL, not residue minus three, so one residue's value was written across three. Negative auth numbering is routine in NMR and EM entries, where it marks expression-tag remnants. Reported against `qscore.py`; `ensembles.py` had it identically and is fixed with it. **Quoting did not fix this** — see below; the escape that does is a backslash.
- **The negative-residue fix above was incomplete, and the correction is in this same unreleased section.** `resi "-3"` is still read as the range: quoting an identifier does not escape PyMOL's range operator. Checked against PyMOL 3.1.0 on an object holding residues -3, 1 and 2, where the quoted form matched all six atoms instead of residue -3's two. The escape the grammar honours is a **backslash**, and it composes with the quoting (`resi "\-3"`) and with `+` lists. `quote` now emits both, and `tests/wiggles/test_selection_live.py` asserts atom counts from a live session rather than the spelling of the string — which is how the original fix passed a test while the bug was still there.
- Neither shape can be reintroduced at a call site: `residue_selection` and `residue_clause` in `mcpymol.wiggles.atoms` are the only places that name a residue, and both go through `quote`, which quotes (for the blank chain) and escapes a leading minus (for the range). A quote *inside* an identifier is refused rather than guessed at, since nothing in the PDB or mmCIF grammar puts one there.
- **Re-stashing B-factors destroyed the originals it existed to protect.** `stash_bfactors` replaced its entry unconditionally, but every caller reads `b` *after* an earlier view may already have overwritten it — so occupancy_view followed by qscore_view saved the occupancies as though they were the user's data, and `restore_bfactors` then wrote those back and reported success. The first stash now wins, and a restore clears it so the next view can take a fresh baseline. `has_stash` existed for exactly this and had no caller.
- **`density_view` contoured in the wrong units when PyMOL was not normalising.** It converted levels to σ on the assumption that `normalize_ccp4_maps` is on, without ever checking — while `local_resolution_view`, which defines the function that checks, consulted it for its colour ramp but not for the contour it computes twenty lines later. With normalisation off a level is read as a raw density, so EMD-30913's published `0.05` went in as `3.16` and drew an empty mesh under a report stating the depositor's own contour had been applied. The level sent now follows the setting, the report says which units reached `isomesh`, and a map whose header cannot supply the conversion is refused rather than contoured on the unconverted number.

- **A negative RMS inverted every conversion it touched.** The guard was `if not header.rms`, which rejects only exactly zero — but MRC2014 writes `rms = -1` (with `dmean = -2`) to mean *statistics not computed*, which is what `mrcfile` leaves behind whenever a map is saved without `update_header_stats()`, and plenty of pipelines never rewrite them. A negative RMS divides cleanly, so it sign-flipped every sigma conversion: `to_sigma(0.05)` returned −2.05, reported as a contour beside a "map sigma: −1", and in `local_resolution_view` it turned ascending Ångström breakpoints into a descending ramp — binding blue to the *worst*-resolved density under a legend stating blue was the best, the exact inversion the module refuses when a user supplies it. `usable_rms` now gates both conversions.
- **`ensemble_spread_view` measured rigid-body drift as flexibility.** Spread was the RMS of each atom about its own position across states, with no superposition step anywhere in the module — and PyMOL loads a multi-model PDB's models as states without aligning them. An MD trajectory that drifted across its box, or two independently refined ensemble members 2 Å apart, reported that offset at *every* residue and painted an internally rigid molecule red end to end. States are now fitted onto the first before any coordinate is read; `superpose=False` opts out, and the report says which happened, because `intra_fit` moves coordinates in the session.
- **The map and provenance registries were never invalidated.** Both are keyed by PyMOL object name and nothing evicts from them, so deleting a map left its header behind: the next contour was converted with the statistics of a volume no longer loaded, under a provenance banner asserting a measurement for whatever now held the name. `loaded_map` now takes the port and drops records whose object has gone. A map deleted and *replaced* under the same name still passes — PyMOL does not expose an object's source file — so `density_view` prints the path its header came from, which is the one thing that makes a substitution visible.
- **`unsharpened` was read as sharpened.** Provenance tokens were scanned in declaration order and `sharp` occurs inside `unsharp`, so a file whose name recorded that it had *not* been sharpened was suggested as sharpened — inverting the one thing its depositor had troubled to write down. A token contained in another match is now discarded, so the longer reading wins whatever order the table is written in.
- **Longest-token-wins alone lost the category priority, and that is now restored.** The declaration-order scan it replaced was carrying one: `Provenance` is ordered from observed to invented, and the warnings are not interchangeable. `postprocess` (11 characters) is longer than `emready` (7), so `postprocess_emready.mrc` — the ordinary name for running EMReady on a RELION postprocess map — was suggested as SHARPENED, offering "features are easier to over-read" in place of "treat this map as a hypothesis, not a measurement". The same happened in reverse wherever the network token was longer: `kmeans20_emready.mrc` and `3dflex_deepemhancer.mrc` lost GENERATED, which is the stronger claim still — no particle was reconstructed into that density at all. Six filenames change classification; every case the old rule got right is unchanged. Priority is read from the `Provenance` enum rather than from the token table, which lists NN_ENHANCED before GENERATED and would invert that pair.
- **`spheroid`'s default was the value its own docstring says fails.** `cmd.spheroid` resolves its argument as an object name and answers "Object not found" for `"all"`, while its own default is the empty string meaning *every object*. Calling the tool the way its schema says you may — with no arguments — therefore always failed.
- **`frame()` could not do the thing its schema advertised.** The parameter description says to omit it to query the current frame, but `cmd.frame` takes a required argument, so the no-argument call came back "missing 1 required positional argument". The query half now goes to `cmd.get_frame`. Found by the new defaults test below, not by review.
- **`atom_properties` split commas inside expressions.** `properties` is an expression list evaluated by `iterate`, so a comma inside brackets or quotes belongs to the expression; splitting on every comma tore `resi in ('1','2')` into fragments handed to PyMOL as field names.

## [1.5.2] - 2026-11-03

A dependency-only hotfix: 1.5.0 and 1.5.1 have failed to start on a fresh
install since 2026-08-31. No tool changed. If yours already runs, this changes
nothing about how it behaves.

### Fixed

- **A fresh install of 1.5.0 or 1.5.1 no longer starts.** `mcp` had a floor and no ceiling. mcp 2.0 removed `mcp.server.fastmcp`, which `app.py` imports. What kept mcp 2 out was an accident: `fastmcp`, a dependency nothing in `src/` imports, pinned `mcp<2` throughout its 3.x line. fastmcp 4.0 (2026-08-31) lifted that pin, and since then `pip install mcpymol` has resolved to fastmcp 4 and mcp 2.3, where `import mcpymol.server` raises `ModuleNotFoundError`. Existing installs kept working until something upgraded them. CI never saw it because every job installs from `uv.lock`.
- **The dependencies are now what the code imports.** `fastmcp` is gone, along with the 39 packages that came only with it (authlib, keyring, websockets, opentelemetry and others). `mcp` is capped at `<2`, and its floor rises to 1.30.0, which carries every security fix released for 1.x. `pydantic` was imported by nine modules without being declared, and now is.
- **The `print` extra's bounds now match what works.** These ranges were checked by running the real mesh repair, not just the mocked suite. `pymeshlab>=2023.12` allowed a release with no macOS arm64 wheel; the floor is now 2023.12.post2. `trimesh>=5.0.0` was only Dependabot raising the floor, and trimesh 4.0 repairs the same meshes, so it is back to `>=4.0`. `networkx` was redundant, since scikit-image already requires it.
- **The lockfile picks up the security fixes now within reach:** mcp 1.30.0, pyjwt 2.15.1, python-multipart, starlette, cryptography, pillow and others. The advisories were all in HTTP, SSE and auth code that MCPymol's stdio transport never reaches. Pillow on Python 3.10 cannot be fixed: its fixes are in pillow 13, which requires 3.11. 1.5.x is the last line to support Python 3.10.

## [1.5.1] - 2026-08-09

Everything v1.5.0 got wrong about installing and upgrading itself, plus two
rendering bugs found in the same review. Most of this release is documentation,
which is the part of a newly-published package people actually run.

Note that the PyPI project page for 1.5.0 renders its figures broken and always
will — a published version is immutable. This release is what fixes that page.

### Added

- **README: a Requirements section.** The README never said that PyMOL itself is a prerequisite — the one thing MCPymol cannot install for you — nor which Python version the bridge needs, nor that the two halves must share a machine because several tools hand files between them through the filesystem. Optional per-feature dependencies (APBS, the `print` extra, network access) are now listed in one table instead of being mentioned wherever the feature happens to be described.
- **README: an Upgrading section.** A failure mode that only exists now that the two halves install separately: both the plugin and the bridge keep running the code they loaded at startup, so upgrading the package on disk changes neither until PyMOL is restarted and the MCP server reconnected. It presents as a tool that is broken or a fix that did not take, rather than as a stale process. Troubleshooting gains rows for that and its neighbours.
- **Release: the tagged commit is now tested before it is published.** `ci.yml` runs on push-to-main and pull_request, neither of which a tag push matches, so nothing between `git tag` and an immutable PyPI upload had ever run the test suite against that exact ref — `twine check` validates metadata, not whether the code works. `release.yml` now runs the full matrix, lint and mypy against the tagged commit, and the build job depends on them.
- **README: the two opt-in test suites.** `pytest -m live` and `pytest -m network` were undocumented despite being the safety net for everything the mocked suite structurally cannot see, including the warning that `-m live` clears the PyMOL session it connects to.

### Fixed

- **The README documented two flags that do the opposite of what it said.** `render` was described as taking `ray_trace=False` for "a fast unshaded snapshot", and `turntable` as defaulting to "the fast OpenGL renderer". Both paths were disabled in v1.4.0 — PyMOL's OpenGL frame grab needs its GUI thread, which the plugin does not run on, and it silently wrote blank images — and `rendering.py` forces ray-tracing back on in both tools.
- **`render`'s tool schema still advertised the fast path v1.4.0 removed.** That release updated `turntable`'s parameter description to say the flag is ignored but left `render`'s promising "a fast unshaded viewport grab", and its docstring recommending the flag for quick checks. The schema is the only one of these an agent ever reads, so the correction had been made everywhere except where it mattered: a model trying to keep a mid-iteration check cheap would ask for the fast path, get a full ray-trace, and on a large assembly hit the socket timeout instead of a quick preview.
- **`assets/pharmacophore_view.png` was referenced but never committed**, so it rendered as a broken image on GitHub — the only one of twenty image references that did not resolve. Generated through the tools, as the other figures were.
- **Stale counts.** 87 primitives and 121 tools, not "~85" and 120; 555 tests by default plus 122 opt-in, not 469. Re-derived from the registered tool list and a test collection rather than counted by hand. The module layout table was missing `analysis.py`, `pdbtext.py`, `style.py` and `cli.py`, which have existed since the v1.3.0 split and v1.4.0.
- **An oversized render deleted the file it told you about.** Above `MCPYMOL_MAX_IMAGE_BYTES`, `render` returns a message instead of the image — but when no `filename` was given it unlinked the temporary PNG on the way out while reporting that the render "is at a temporary file". A 4000×3000 publication render therefore cost minutes and left nothing behind. The file is now kept whenever that branch is taken, and the message names its real path.
- **The README's images are absolute URLs**, because the README is also the PyPI long description and PyPI does not rewrite relative links — all 17 figures, including the hero image, rendered as broken-image placeholders on the project page. They point at `raw.githubusercontent.com` now. 1.5.0's page keeps its broken images permanently, since a published version is immutable — this release is the fix.
- **Install instructions assumed a PATH that the recommended installers do not guarantee.** `uv tool install` and `pipx` put `mcpymol` in `~/.local/bin`, which is not universally on PATH — and an MCP client launched by the OS inherits no PATH at all — so registering a bare `mcpymol` failed at client launch with `spawn mcpymol ENOENT` rather than at install time. Every registration command now uses the full path, with a `which mcpymol` check up front.
- **The `uvx` path documented a bridge with no plugin.** It dropped the clone that used to supply `plugin.py` while still pointing at step 2's `mcpymol --install-plugin`, a command a reader on that path does not have. It now gives `uvx mcpymol --install-plugin` and states the cache-path trade-off plainly.
- **The `print` extra pointed at the wrong environment.** A bare `pip install 'mcpymol[print]'` resolves to whatever `pip` is first on PATH, which on a Homebrew or Debian Python fails outright and otherwise installs a second copy of MCPymol somewhere the tool venv never sees. Now given per install method.
- **Upgrade advice was `uv`-only**, including in the row a reader reaches from the "no `uv`" install path, where the single prescribed remedy is the one tool their environment forbids.
- **The pinned `actions/checkout` SHA was labelled `# v4` but is v7.0.1** — in `ci.yml` since it was written, and copied into `release.yml`. The comment is what Dependabot reads to decide what "within v4" means, so the label being three majors behind the pin invites a silent downgrade.
- The test command no longer sets `PYTHONPATH=src`, which `pyproject.toml`'s `pythonpath` setting has made redundant.

## [1.5.0] - 2026-08-08

The first release published to PyPI. Every install path used to begin with
`git clone`, which asked people to take on a checkout in order to run a tool
they only wanted to use — and left them tracking a branch rather than a
version.

### Added

- **MCPymol is on PyPI**, so the bridge installs with `uvx mcpymol` and needs
  no checkout. This is what makes the plugin the interesting half of the
  install: the bridge is now a one-line registration, while the plugin still
  has to be found on disk, which is what `--install-plugin` below is for.
  Publishing uses [Trusted Publishing](https://docs.pypi.org/trusted-publishers/)
  over OIDC rather than a stored API token, so there is no long-lived
  credential in the repository to leak or rotate. `.github/workflows/release.yml`
  fires on a `v*` tag and refuses to upload if the tag disagrees with the
  version it just built — a published version is immutable, and a wrong one
  can only be superseded, never replaced.
- **`mcpymol --install-plugin`**, which adds the plugin to `~/.pymolrc.py` so PyMOL loads it on startup, plus `--plugin-path` and `--uninstall-plugin`. The plugin runs inside PyMOL, whose Python is a separate interpreter that cannot import this package, so it has to be loaded from a file path — and finding that path inside a pip or `uvx` installation is unpleasant. The block is delimited and rewritten in place rather than appended, because the path it embeds points into one installation and goes stale on upgrade; re-running the command is the fix. Existing `.pymolrc.py` content is preserved, and uninstall removes only the managed block.

## [1.4.0] - 2026-08-08

MCPymol could make pictures but could not answer questions with numbers, and
a good deal of it did not work at all. This release closes both gaps.

### Added

- README figures for `conservation_view` (lysozyme against 5,507 homologues), `print_ribbon_view` (GFP, spine visible through the strands), `textbook_view` (haemoglobin), `cinematic_view` (GroEL down its seven-fold axis), `plddt_view` (AlphaFold spike) and `superposition_view` (adenylate kinase), plus an animated `turntable` GIF. Every documented view now carries an image except `bfactor_view` and `pointillist_view`, which are deliberate: `putty_view` already illustrates the B-factor palette, and the pointillist style is decorative.
- Real output for `contact_report`, `interface_report` and `structure_info` in the README. These were the headline features of the release and it showed none of what they return. The contact and interface examples use the same systems the `ligand_view` and `interface_view` figures show, so the picture and the numbers sit side by side.
- All figures were produced through the tools themselves, which is how three of the bugs in this release were found.
- **Configuration section in the README.** Six of nine `MCPYMOL_*` environment variables were undocumented, including every timeout introduced in v1.3.0. All nine now have a table entry saying what they do and what the default is.
- Troubleshooting rows for the failure modes introduced since v1.2.1: files not appearing when PyMOL and the bridge are on different machines, oversized renders, AlphaFold entries with no model at a given version, blank views from empty selections, and slow-operation timeouts.
- `tests/test_docs.py` keeps these true. An environment variable added in code but not the README, a headline tool missing from the README or the skill, a doc that sends readers to `server.py` for implementation, or a `CONTRIBUTING` that stops warning about facade patching — each now fails CI. Written after this audit found six undocumented variables and three unmentioned tools by hand; the point is not to do that by hand again.
- **A live test suite: `pytest -m live`.** Calls every registered tool once against a running PyMOL. It exists because the mocked suite asserts the payload we *send*, which cannot tell whether PyMOL will accept it — and every wiring bug shipped so far lived in that gap. It asserts only that a tool is *wired* (the action resolves, the arguments bind), not that it succeeds, since many tools error honestly without the right context. Arguments are generated from a table keyed on parameter name, so a new tool is swept automatically, and anything unsuitable must be excluded with a stated reason. It skips cleanly when no PyMOL is listening, and refuses to run against a PyMOL whose camera has stopped responding rather than reporting a hundred confusing failures.
- **`contact_report`** — lists contacting residue pairs across two selections, closest first, with the minimum heavy-atom distance, the number of atom contacts, and a classification: salt bridge, hydrogen bond, hydrophobic, polar contact, or π-stacking (parallel vs T-shaped by interplanar angle). Criteria are heavy-atom distances, documented and pinned by tests, since crystal structures usually have no hydrogens — so a reported hydrogen bond is a donor–acceptor pair with plausible geometry, not one verified against a hydrogen position. Neighbour search uses a uniform grid, keeping cost linear in atom count; a property test asserts it finds exactly what a brute-force scan does. Ring perception needs bond orders that a PDB dump does not carry, so aromatics cover the standard aromatic amino acids only and ligand rings are reported as hydrophobic rather than mis-called as stacking.
- **`interface_report`** — buried surface area from ΔSASA (free minus bound, halved for the per-side figure papers quote), residues ranked by how much surface each buries, and a breakdown by residue chemistry. Interprets the number against PDB-wide survey thresholds: under ~400 Å² per side is usually crystal packing, over ~1000 Å² a substantial and likely specific association. Uses `get_area(load_b=1)`, which writes per-atom SASA into the B-factor column, so the whole per-residue breakdown costs one area calculation and one dump rather than hundreds of round trips — and sets `dot_solvent=1` first, without which the numbers are wrong but plausible-looking.
- **`structure_info`** — what a structure actually is, in one call: chains, atom/residue/water counts, ligands, states and space group from PyMOL, plus title, method, resolution, release date and source organism from the RCSB data API. Best-effort metadata: skipped for objects not named after a PDB entry, degrading to "no metadata" when the network is unavailable. Flags a probable predicted model when the B-factor column looks like pLDDT.
- **`get_sequence`** — there was previously no way to get a sequence out of MCPymol at all. Returns FASTA plus the two things the sequence hides and that routinely cause silent mistakes: the numbering offset (PDB numbering rarely starts at 1, so "residue 50" in a paper and position 50 in the sequence are usually different residues) and chain breaks where loops went unmodelled.
- **`sasa`**, **`rms_cur`**, **`count_atoms`** — value-returning primitives. `count_atoms` in particular catches an empty selection, which is otherwise invisible until the render comes out blank.
- New `mcpymol.pdbtext` module for parsing the PDB dumps that are the only way to read data back out of PyMOL, and `mcpymol.analysis` for the reporting tools.

### Fixed

- **`CONTRIBUTING.md` recommended a test pattern that does not work.** It told contributors to mock `mcpymol.server.send_request`; after the v1.3.0 package split that rebinds only the facade's name, leaving every module still calling the real socket — the exact trap that broke 74 tests during the split. It now shows the correct target and explains why.
- Docs still pointed at `server.py` for the bridge protocol and for "after editing", both of which moved when the package was split.
- README claimed 340+ tests (actually 469) and ~60 primitives (actually ~85, of 120 tools).
- `load_structure` was documented nowhere, so opening a local PDB/mmCIF — a built model, a docking pose, an MD frame — was undiscoverable. `fetch_alphafold` and `count_atoms` were likewise unnamed.
- **A failed fetch could destroy the session** ([#15](https://github.com/chemrich/MCPymol/issues/15)). `fetch_structure` ran `reinitialize` *before* the fetch, and `cmd.fetch` does not raise when a download fails — it simply produces nothing, which the plugin reports as success. Behind a proxy or VPN blocking the RCSB, fetching one PDB code therefore destroyed every unrelated structure already loaded and returned "Successfully fetched". All three loaders now delete only the object being replaced, verify atoms actually arrived, and clear the rest of the session by name afterwards; settings still get a clean slate via `reinitialize settings`, which does not touch objects.
- **Seven tools failed on every call** ([#15](https://github.com/chemrich/MCPymol/issues/15)), all from wrappers pointing at functions PyMOL does not have or binding arguments in the wrong order: `as` (a Python keyword, so `getattr(cmd, "as")` can never resolve — the real name is `cmd.show_as`), `util_color_by_element` (→ `util.cnc`), `util_color_secondary` (→ `util.cbss`), `spheroid` (needs an object name, not a selection), `h_fill` (has no selection parameter; the selection was landing in `quiet` and being `int()`d), `sculpt_iterate` (arguments reversed) and `symexp` (needs an object *and* a selection). Names and argument orders were verified by introspecting `cmd` and `pymol.util` in a live PyMOL 3.1, and are now pinned by tests — the existing suite could not catch any of them, because it asserts the payload we send and we were faithfully sending the wrong thing.
- **`render(ray_trace=False)` and `turntable`'s default produced blank images.** PyMOL's unshaded OpenGL capture needs its GUI thread, and the plugin dispatches on a socket worker thread, so `cmd.png(..., ray=0)` writes a flat empty frame rather than failing. `turntable` defaulted to exactly that path, so every frame of every animation came out blank. Both tools now always ray-trace and say why if asked not to. `render` additionally detects a flat single-colour result and points at `list_objects`/`count_atoms`, which catches the other way a render comes back empty: a scene where nothing is visible.
- **Every view preset rendered with a transparent background.** They set `bg_color` but never `opaque_background`, so ray-traced output kept an alpha channel: correct in the viewport, transparent in the PNG, which appears white wherever the image is used. 16 of 17 presets were affected, `textbook_view` included (its background is white rather than black, and had the identical bug). Backgrounds now go through `style.set_background()`, and a test fails if any module reaches for `bg_color` directly.
- **`superposition_view` was unreachable through the documented workflow.** `fetch_structure` calls `reinitialize`, so fetching a second structure wiped the first — leaving nothing to compare, for a tool whose only purpose is comparing two structures. All three loaders (`fetch_structure`, `load_structure`, `fetch_alphafold`) take `replace=False` to add to the session instead of clearing it.
- **`fetch_alphafold` did not work for any accession.** It built the model URL from a hardcoded `AF-{accession}-F1-model_v4.cif` template. AlphaFold DB has since moved to v6 and *removes* retired versions, so every one of those URLs now 404s — the feature shipped in v1.3.0 and resolved for nothing. It also assumed the filename is keyed by accession, which is untrue for some entries: SARS-CoV-2 spike (P0DTC2) is served as `AF-0000000365840314-model_v1.cif`, so no accession-based filename exists for it in any version.
  The file URL is now asked for rather than constructed, via AlphaFold DB's prediction API, which handles version drift, non-accession entry IDs and multi-fragment proteins. `model_version` becomes an optional pin that bypasses the lookup, and defaults to unset.
  Every test for this code mocked `urlopen`, so they proved the URL was *built* correctly and never that it *existed*. A `network`-marked test now resolves and fetches a real model; it is deselected by default so CI stays offline-safe (`pytest -m network` to run it, worth doing before a release).
- A malformed accession returned a raw `HTTP 400 Bad Request`; it now says which accession failed and what one looks like.
- **Measurements discarded their results.** `cmd.distance`, `angle` and `dihedral` all return the quantity they measure; the wrappers threw it away and returned "Measured distance between 'X' and 'Y' as 'd1'" or "Executed angle successfully." So you could ask MCPymol to measure something and it would draw the measurement without telling you the number. All three now report the value. PyMOL's -1 "nothing matched" sentinel is reported as a failure rather than as a measurement of -1 — but only for quantities that cannot be negative, since -57.8° is an ordinary α-helical phi.

### Changed

- **Every tool parameter is now documented in the JSON schema.** All 222 parameters across all 120 tools carry `Annotated[..., Field(description=...)]`, so `inputSchema.properties[].description` is populated and the model no longer infers arguments from parameter names — it can see what `spectrum`'s `expression` accepts, what `clip`'s `mode` means, and what syntax `mset` takes. The parameter docs previously lived only in `Args:` docstring blocks, which FastMCP does not read into the schema; the signature is now the single source of truth and the redundant blocks are gone. Verified by diffing every schema across the change: descriptions added, types/defaults/`required` untouched. Four new tests keep it that way — a tool shipping an undocumented, empty or one-word parameter description now fails CI.
- README reorganised around the question being asked rather than the tool being called, with a "What can I ask?" index and a new Analysis section. The mcpymol-guide skill leads with report-first-then-draw.
- `views._read_ca_bfactors` and `comparison._parse_ca_coords` now delegate to `pdbtext.parse_atoms`, which additionally handles missing element columns, insertion codes and truncated lines. The two `_distance` implementations are unified on `pdbtext.distance3d`.

### Developer experience

- **Pre-commit hooks that actually cover the gates.** The config had a single hook (the `uv.lock` guard), `pre-commit` was not a declared dependency so nothing installed the binary `CONTRIBUTING.md` told you to run, and the hook was consequently never installed in practice — every check happened only in CI, after a push. Commits now run ruff, ruff-format, mypy, the lockfile guard, and file hygiene (trailing whitespace, final newline, line endings, YAML/TOML validity, merge markers, stray `breakpoint()`, oversized files); pushes run the full test suite, which is too slow at ~16 s to pay per commit. `pre-commit` is in the dev group, and `default_install_hook_types` wires both stages from one `pre-commit install`.
- The ruff and mypy hooks run via `uv run` rather than pre-commit's own tool repos, so they use exactly the versions in `uv.lock` — the versions CI installs — instead of introducing a second place for tool versions to drift.

### Note on dependencies

The optional `analysis` extra was considered and deliberately not added. For this scope nothing earned it: per-residue SASA comes from PyMOL in two round trips via `get_area(load_b=1)`, and contact search is a spatial-grid problem where numpy would add a second code path for no measurable gain. The core stays dependency-free.

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
