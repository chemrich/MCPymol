---
name: mcpymol-guide
description: >-
  How to drive the mcpymol MCP server well — structure prep, the view-preset
  catalogue, PyMOL selection syntax, rendering, and 3D-print STL export.
  TRIGGER when the user asks to visualize/fetch/render a protein structure
  through the mcpymol tools, apply a *_view preset, or produce a 3D-printable
  model/STL. SKIP for unrelated PyMOL-less work.
---

# Driving the mcpymol MCP server

The `mcp__mcpymol__*` tools talk to a running PyMOL over a socket. Prefer the
high-level tools (`fetch_structure`, the `*_view` presets, `show`/`color`/
`select`) over `execute_pymol_command` — the presets do coloring,
transparency, H-bonds and labels in one call. Reach for
`execute_pymol_command` only when no dedicated tool covers the need.

## Report first, then draw

The presets make pictures; the report tools answer questions with numbers.
When the user asks *why* or *how much* — what holds this ligand, how big is
this interface, where do these differ — run the report and quote the numbers,
then apply the matching view if a picture helps. Answering a quantitative
question with only a rendering leaves the user to measure it themselves.

| The question | Reach for |
| --- | --- |
| What is this structure? | `structure_info` |
| What is the sequence / numbering / where are the gaps? | `get_sequence` |
| What contacts what, and how tightly? | `contact_report` -> `ligand_view` |
| How big is this interface, which residues matter? | `interface_report` -> `interface_view` |
| Where do two structures differ? | `superposition_view` |
| How far / what angle / how much area? | `distance`, `angle`, `dihedral`, `sasa`, `rms_cur` |
| Does my selection match anything? | `count_atoms` |
| Occupancy / altloc / per-atom B-factor? | `atom_properties` |

`contact_report` classification is heavy-atom geometry (no hydrogens assumed):
salt bridge <= 4.0 A, H-bond <= 3.5 A between N/O, hydrophobic <= 4.5 A,
ring centroids <= 5.5 A. Report it as such — a "hydrogen bond" here is a
donor-acceptor pair with plausible geometry, not one verified against a
hydrogen. Aromatic detection covers the standard aromatic amino acids only;
ligand rings appear as hydrophobic contacts.

`interface_report` gives buried surface area per side. Under ~400 A^2 usually
means crystal packing rather than a real interface; over ~1000 A^2 means a
substantial, likely specific association. Say which regime the number is in.

## Core workflow

0. **Predicted structures.** `fetch_structure` also accepts a UniProt
   accession (`P69905`, `af-P69905`) and routes it to AlphaFold DB, colouring
   by pLDDT confidence; `fetch_alphafold` is the same thing called directly,
   and takes `model_version` / `fragment` for older or split entries. For predicted models use `plddt_view`, never
   `bfactor_view`/`putty_view` — pLDDT lives in the B-factor column and those
   presets read it backwards (they assume low = rigid; low pLDDT = unreliable).
0b. **Local files.** `load_structure(file_path, obj_name)` takes a PDB/mmCIF
   off disk — a built model, a docking pose, an MD frame — and applies the
   same multimer heuristic and styling as a fetch. Use it whenever the user
   has a file rather than an accession.
1. **`fetch_structure(pdb_code)`** — grabs the biological assembly when one
   exists, runs a BFS chain-contact heuristic (`multimer_cutoff`, default
   8.0 Å) so functional multimers stay whole while crystallographic copies
   are dropped, hides waters/additives, and applies the default style. The
   object is named after the PDB code unless `obj_name` is given.
2. **Ground yourself before composing selections.** Use `list_objects`,
   `list_chains(obj_name)`, `list_ligands(obj_name)` instead of guessing
   object names, chain IDs, or 3-letter ligand codes.
3. **Apply a view preset** (below) or build the scene manually.
4. **Check a selection with `count_atoms`** before building a scene on it —
   an empty selection is invisible until the render comes out blank.
5. **Render with `render()`** — it ray-traces and returns the image directly,
   so you can see the result and iterate. Prefer it over `ray` + `png`, which
   only leave a file on disk. Pass `ray_trace=False` for a fast unshaded
   check of a selection or orientation, and keep the width/height modest —
   the image is inlined into the conversation.

## View-preset catalogue

All take `obj_name`. Each is a one-call scene.

- **Binding / functional sites:** `ligand_view`, `pocket_view`,
  `pharmacophore_view`, `interface_view` (chain–chain), `crosslink_view`
  (disulfides & metal coordination), `mutation_view`.
- **Per-residue scalar coloring:** `bfactor_view`, `putty_view` (B-factor as
  tube width), `conservation_view` (evolutionary; first call hits the
  ColabFold MMseqs2 API and is slow, then cached).
- **Surface chemistry / electrostatics:** `hydrophobic_surface_view`,
  `electrostatic_view` (approximate), `poisson_boltzmann_view` (true PB;
  needs `apbs` + `pdb2pqr`).
- **Illustrative:** `textbook_view` (cel-shaded; needs `ray`),
  `cinematic_view` (fog/shadows; needs `ray`), `pointillist_view`.
- **Predicted-model confidence:** `plddt_view` (AlphaFold pLDDT in the
  official palette, plus a confidence breakdown).
- **Comparing two structures:** `superposition_view(mobile, target)` —
  superposes and colours the mobile structure by per-residue shift, so you
  see *where* it moved rather than just an RMSD. Names the worst-shifted
  residues.
- **3D printing:** `print_ribbon_view` — see below.

## Cryo-EM: occupancy, ensembles and maps

Ten tools for heterogeneity, model quality and volumes. They are opinionated in
one specific way: **they refuse rather than draw something meaningful-looking
and wrong.** A refusal here is the tool working, not an error to route around —
report it and say what differed.

- **Is the atom really there:** `occupancy_view` (per-atom crystallographic `q`;
  alternates de-emphasised in proportion), `altloc_view` (every alternate at
  once, occupancies labelled). Note that "occupancy" names two incompatible
  quantities — `q`, and the fraction of *particles* containing a subunit. These
  tools report the first only, and never infer the second from it.
- **Ensembles:** `ensemble_spread_view` (per-residue RMS deviation across
  states, as colour and putty width — a description of disagreement, not an
  error bar), `morph_states` (interpolates, but **refuses** states that do not
  share a topology; `cmd.morph` is Incentive-only, so on open-source PyMOL the
  check reports and the morph does not happen).
- **Undo the B-factor overwrite:** `restore_bfactors`. The three views above
  push their values into the B-factor column so PyMOL can spectrum-colour by
  them; the originals are stashed first.
- **Model quality:** `qscore_view(obj_name, validation_path)` — reads
  per-residue Q-scores out of a wwPDB validation report. No network, no map, no
  compute. Expect entries before Sept 2023 to have none, and expect some real
  scores to be negative.
- **Map geometry:** `map_info(path)` — header only, no PyMOL, no network. Leads
  with voxel size, which is derived rather than stored and is a classic source
  of a silently wrong answer.
- **Volumes:** `load_map(path, name, provenance)` then `density_view(map_obj,
  selection, level, units)`. **Always pass `units`.** PyMOL contours in σ and
  EMDB publishes absolute contours; EMD-30913's `0.05` is 3.16 σ, and passing it
  as a σ level contours noise. Set `provenance` when you know it — it defaults
  to `unknown` and is deliberately never guessed, because a generated volume and
  a measured one are the same isosurface once drawn.
- **Where the map is good:** `local_resolution_view(map_obj, res_obj)` — colours
  the isosurface by a local-resolution volume. Refuses when the two do not share
  a voxel grid. Breakpoints are given in Å and reported in both Å and σ.

## Rendering, movies and sessions

- `render(width, height, ray_trace, filename)` — the one to reach for. Returns
  the image so you can look at it.
- `turntable(obj_name, frames, out_dir)` — a 360° PNG sequence plus the
  ffmpeg line to assemble it. Defaults to the fast renderer; `ray_trace=True`
  is for finals only (36 traced frames of a big assembly can take an hour).
- `save_session(filename)` / `load_session(filename)` — `.pse` round-trips
  the whole scene (objects, colours, camera, scenes). Save before
  experimenting with a scene that took effort to build.

## Selection syntax (PyMOL)

Selections are an algebra, not just residue lists:
`chain A`, `resi 10-50`, `resn ATP`, `polymer.protein`, `organic`,
`ss H` / `ss S` (helix/strand), `byres (... around 5)`, combined with
`and`/`or`/`not` and parentheses. Cross-object selections are fine:
`(1abc or 1abc_spine)`.

## 3D-print STL export

Two-step: configure the scene, then export.

```
print_ribbon_view(obj_name="1ema")
print_export(obj_name="1ema", groups="1ema=(1ema or 1ema_spine)",
             representation="cartoon", method="voxel", voxel_pitch=0.2)
```

`print_ribbon_view` builds chunky β-arrows + a fat helix and a separate
`<obj>_spine` (`cartoon tube`, ignores secondary structure) that threads
unbroken through every strand→loop junction — fused on export into one
watertight solid, with the spine acting as internal rebar for rigidity.

`print_export` key arguments:

- **`representation`** — `"surface"` (default) exports the molecular
  surface; **`"cartoon"`** exports the *currently displayed* cartoon/tube
  geometry (what `print_ribbon_view` sets up). Using the default with a
  ribbon scene silently gives you a surface blob.
- **`method`** — `"voxel"` (robust for thin tubular ribbon geometry;
  `voxel_pitch` ~0.2 keeps ribbon detail), `"poisson"` (bulky chains),
  `"auto"`.
- **`groups`** — `label=selection` pairs, one STL per colour. In cartoon
  mode each group should name whole object(s).

### Hard-won verification rules

- **Verify the STL itself, not the viewport.** `ray`/`png` show the live
  PyMOL scene, which can look like ribbons while the exported mesh is a
  surface. PyMOL cannot load STL back, so inspect the file with the
  project's own `trimesh` (already a `print`-extra dependency):
  check `is_watertight`, `len(mesh.split())` (want **1** component), and
  to tell ribbon from surface compare bbox-fill fraction and
  area/volume — a ribbon is thin (low fill, high area/volume), a surface
  is a fat shell.
- **A printable result is watertight AND a single component.** Fragmented,
  non-watertight output won't slice.
- The STL is in the same coordinate frame across groups — in a slicer,
  load the first then add others as *parts* (don't re-centre).

### After editing anything under `src/mcpymol/`

The MCP server imports the package at startup. Code changes are **not live**
until the user reconnects it (`/mcp` → reconnect mcpymol, or restart Claude
Code) — you cannot restart it yourself. After a reconnect, confirm the new
behavior from tool output, don't assume it loaded (a stale server returns
the old result).

## Troubleshooting

- *"Socket connection failed"* on every tool → PyMOL plugin not running /
  server disconnected. Ask the user to reconnect mcpymol via `/mcp`.
- Long `get_fastastr`/`get_chains` JSON parse errors → pre-2026-05 build.
- `poisson_boltzmann_view` fails → missing `apbs`/`pdb2pqr`.
- 3D-print tools say the `print` extra is missing → `uv sync --extra print`.
- A tool reports a file "did not appear" → PyMOL and the bridge are not on the
  same machine; they exchange files through the filesystem.
- `render` returns a size complaint instead of an image → ask for a smaller
  width/height, or raise `MCPYMOL_MAX_IMAGE_BYTES`.
- A scene comes out blank → check the selection with `count_atoms` before
  assuming the view preset failed.
