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

## Core workflow

1. **`fetch_structure(pdb_code)`** — grabs the biological assembly when one
   exists, runs a BFS chain-contact heuristic (`multimer_cutoff`, default
   8.0 Å) so functional multimers stay whole while crystallographic copies
   are dropped, hides waters/additives, and applies the default style. The
   object is named after the PDB code unless `obj_name` is given.
2. **Ground yourself before composing selections.** Use `list_objects`,
   `list_chains(obj_name)`, `list_ligands(obj_name)` instead of guessing
   object names, chain IDs, or 3-letter ligand codes.
3. **Apply a view preset** (below) or build the scene manually.
4. **Render**: `ray(width, height)` then `png(filename)`, then `Read` the PNG
   to actually see it. A bare `png` without `ray` is fast but unshaded.

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
- **3D printing:** `print_ribbon_view` — see below.

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

### After editing `src/mcpymol/server.py`

The MCP server imports `server.py` at startup. Code changes are **not live**
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
