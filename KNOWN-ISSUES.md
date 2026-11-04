# Known issues

Open defects with enough diagnosis attached to act on. Fix one, delete its
section. Anything here that turns out to be wrong should be corrected in place
rather than left to mislead.

---

## 1. `fetch_structure` can wipe the session — **fixed**

**Was [#15](https://github.com/chemrich/MCPymol/issues/15).**

The mechanism was not the post-fetch cleanup. `_apply_multimer_heuristic`
already returns early on an object with no chains, so it was innocent. The real
chain was:

1. `reinitialize` ran **before** the fetch, unconditionally — that is what
   cleared the session.
2. `cmd.fetch` does not raise when a download fails; it just produces nothing.
3. The plugin's fetch handler returns `{"status": "success"}` without checking
   that an object was created.

So the destruction happened before the failure, and the success report hid it.

All three of the proposed guards were right and are now in place, in
`fetch_structure`, `load_structure` and `fetch_alphafold` alike:

- Only the object being replaced is deleted up front.
- `_verify_loaded` counts atoms after loading and returns an error naming the
  likely cause (a proxy or VPN in front of the RCSB) when there are none.
- Clearing the rest of the session happens **after** that check, by object
  name, via `_clear_other_objects`. Settings still get a clean slate through
  `reinitialize settings`, which resets a previous preset's fog or
  `ray_trace_mode` without touching objects.

Still open, deliberately: the plugin's `fetch` handler reports success without
verifying anything arrived. The bridge-side check makes that harmless, but the
handler is worth hardening — left alone here only to avoid colliding with the
`_jsonable` work on `feat/wiggles-query`.

---

## 2. Seven tools error on every input — **fixed**

**Was also [#15](https://github.com/chemrich/MCPymol/issues/15).** Every
diagnosis in the original report was correct, verified by introspecting
`cmd` and `pymol.util` in a live PyMOL 3.1.

| Tool | Was | Now |
|---|---|---|
| `as` | `cmd.as` — a Python keyword, unresolvable | `cmd.show_as` |
| `util_color_by_element` | `util.color_by_element` — does not exist | `util.cnc` |
| `util_color_secondary` | `util.color_secondary` — does not exist | `util.cbss` |
| `spheroid` | selection in an object-name slot | takes `obj` |
| `h_fill` | selection landed in `quiet` and was `int()`d | takes no arguments |
| `sculpt_iterate` | `(iterations, obj)` | `(obj, state, cycles)` |
| `symexp` | one selection for two parameters | `(prefix, obj_name, selection, cutoff, segi)` |

One correction to the original report: **`spheroid` did not error on every
input.** It works when given an object name and fails only on a selection
expression, so it was a parameter-contract mismatch rather than a dead tool.

`h_fill` is now argument-free and honest about being unusable from here:
`cmd.h_fill` acts on PyMOL's editor pick, which cannot be set over the bridge.
It returns `Editor not active` unless an atom is picked in the GUI. `h_add` is
the selection-based alternative.

Why the suite missed all seven: `tests/test_auto_wrappers.py` asserts the
payload we *send*, and we were faithfully sending the wrong thing. The fix adds
tests that pin each action name and argument order against the real PyMOL API.

---

## 3. numpy results reach clients as a repr string — **fixed**

`_dumps_response` fell back to `repr()` whenever `json.dumps` raised, so the
numpy arrays returned by `get_coords` and `get_atom_coords` arrived as strings
that looked like data and could not be used as data.

A `_jsonable()` pass now runs before the repr fallback, converting anything
exposing `.tolist()` and recursing through dicts and sequences. The repr
fallback is preserved for genuinely opaque objects. Verified against a live
PyMOL: `get_coords` returns `[-0.0009, 0.0637, -0.4905]` rather than
`"array([[...]])"`.

The same branch added the `iterate_to_list` action — the only route to
properties that live on atoms rather than residues or objects — and it is now
exposed to clients as the `atom_properties` tool, so occupancy, altloc and
per-atom B-factor are reachable end to end.

---

## 4. Four presets send commands PyMOL rejects — open

Found while making `execute_pymol_command` report PyMOL's errors: a request log
run through the real tools showed four calls PyMOL refuses or ignores. Each sits
in a direct `send_request` whose result the preset discards, so the preset still
reports success. Verified on PyMOL 3.1.0.

| Where | Call | What PyMOL says |
|---|---|---|
| `views.py` `pharmacophore_view` | `color violet, <lig> and elem C and inring` (and `not inring`) | `Invalid selection name "inring"` — ring carbons are never coloured |
| `views.py` `cinematic_view` | `set fog_color, black` | `unknown Setting: 'fog_color'` |
| `views.py` `pointillist_view` | `color grey30, (<obj>) and polymer.protein and ribbon` | `Malformed selection` |
| `structures.py` loading | `set cartoon_nucleic_acid_color, brightorange, (<obj>) and polymer.nucleic` | `Setting-Warning: … is a object-state-level setting` — stored per atom, never drawn |

The plugin now returns these as errors or warnings in the response; the fix is
at each call site (and a test that reads the response).

---

## 5. GUI-only behaviour the headless live suite cannot see — open

The live suites run against `pymol -ckqK`. These were seen only in the Qt GUI:

- **A modal movie render stalls while the window is minimized.** `mpng` with
  `modal=1` draws no frames until the window is visible again, and PyMOL takes
  no commands meanwhile; `execute_pymol_command` reports that after
  `MCPYMOL_DO_START_TIMEOUT`.
- **Generic actions still run on the plugin's socket thread.** A raw
  `send_request("window", …)` kills a GUI PyMOL (no tool exposes it).
  `execute_pymol_command("window size, 640, 480")` is safe: `do` lines run on
  the main thread.
- **Console output can lag a long request.** The plugin captures what PyMOL
  prints while a request runs (that is how it sees warnings), so Python-level
  console output from the same moment appears when the request ends.
- **`print_export` shows the export view for the length of each save**, then
  restores yours exactly. For each group it also shows only that group's objects
  and pauses rock and stereo, and puts both back as soon as the group is written
  out, before its repair.
- **A warning can be reported one call late.** The plugin captures what PyMOL
  prints while a request runs, but PyMOL prints some messages while it redraws
  after a command has returned (`Extrude-Warning` from `cartoon putty` on
  uniform B-factors, surface update errors). Those land in the next request's
  warnings, under a tool that did not cause them (reproduced headless, 5 of 5).
  Running every action on the main thread would fix it, at the cost of freezing
  the window during ray traces.
- **Re-running the plugin drops a request sent at that moment.** `run
  plugin.py` (through `execute_pymol_command`, or `stop_mcp`/`start_mcp`)
  closes the old listening socket before the new one binds, so a request
  already waiting in its backlog is reset and never run.
- **Plugin and bridge should be the same version.** A released bridge (1.5.x)
  driving this plugin sends `do` with a 10 s budget, but `do` now replies only
  after the command has run, so a long command through `execute_pymol_command`
  times out there although it still runs. `uv tool install` and
  `mcpymol --install-plugin` keep the two together.

---

## 6. Print export limits — open

- **Poisson on open cartoon soup** leaves small fragments (1ubq print ribbon:
  one 857 mm³ piece plus nine of 0.004–2.6 mm³). `bodies="all"` now lists them
  as separate pieces where the old keep-largest step hid them; use
  `method="voxel"` for cartoons.
- **The light path** (`auto` on an already-closed export) counts two closed
  shells that overlap as separate pieces, though they print as one.
- **Grid-aligned geometry is slightly over-corrected** by the voxel-growth
  inset: an axis-aligned rod through the origin grew 0.41 voxel rather than
  0.62, so it comes back 0.05 Å thin at 0.2 Å pitch. Off-grid rods come back
  within 0.012 Å.
- **`print_ribbon_view`'s bridge** is verified on CRO in mmCIF (2b3p, 1ema)
  only. Other chromophores (GYS, CR2, NRQ) and PDB files with LINK records are
  untested, and two CA-less residues in a row are neither bridged nor reported.
- **The thin-feature warning** has no entry for cylindrical helices or the
  nucleic-acid cartoon modes, so those never warn.

---

## Not yet triaged

Nothing currently. New reports should land here with a symptom and a
reproduction before being promoted to a numbered section.
