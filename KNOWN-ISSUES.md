# Known issues

Open defects and limitations with diagnosis attached. Fixed items are recorded in `CHANGELOG.md` and removed from here.

---

## 1. GUI-only behaviour the headless live suite cannot see — open

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
  warnings, under a tool that did not cause them.
- **Re-running the plugin drops a request sent at that moment.** `run
  plugin.py` (through `execute_pymol_command`, or `stop_mcp`/`start_mcp`)
  closes the old listening socket before the new one binds, so a request
  already waiting in its backlog is reset and never run.
- **Plugin and bridge should be the same version.** A pre-1.6.0 bridge driving
  this plugin sends `do` with a 10 s budget, whereas `do` now replies only after
  the command has run, so a long command through `execute_pymol_command` times
  out on the old bridge while continuing to run. `uv tool install` and
  `mcpymol --install-plugin` keep the two together, and the bridge now warns on
  a `plugin_version` mismatch.

---

## 2. Print export limits — open

- **Poisson on open cartoon soup** leaves small fragments (1ubq print ribbon:
  one 857 mm³ piece plus nine of 0.004–2.6 mm³). `bodies="all"` lists them as
  separate pieces; use `method="voxel"` for cartoons.
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

## 3. A raw command can still crash PyMOL — open

PyMOL 3.1.0 aborts or segfaults when `rms`, `rms_cur`, `fit`, `intra_rms` or
`intra_rms_cur` get a selection it cannot resolve. The plugin vets those five
when they arrive as actions, which covers every tool and wiggles' port. A
command string sent through `execute_pymol_command`, such as `rms_cur nosuch,
1ubq`, goes to PyMOL's parser as typed and is not vetted, so it can still take
PyMOL down. Vetting it requires guarding the plugin's private parser.

---

## Not yet triaged

Nothing currently. New reports should land here with a symptom and a
reproduction before being promoted to a numbered section.
