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

## 3. numpy results reach clients as a repr string

`_dumps_response` in `plugin.py` falls back to `repr()` when `json.dumps`
raises, so the client always gets *a* reply. But several `pymol.cmd` functions
return **numpy arrays** — `get_coords` and `get_atom_coords` among them — and
those are not JSON-serialisable. Today they arrive as:

```json
{"status": "success", "result": "array([[1., 2., 3.]])"}
```

A success response carrying a *string* that looks like data and cannot be used
as data. Same failure shape as issue 1: reporting success while returning
nothing usable.

**Fix written, on branch `feat/wiggles-query`.** A `_jsonable()` pass before the
repr fallback converts anything exposing `.tolist()`, recursing through dicts
and sequences. The repr fallback is preserved for genuinely opaque objects
(chempy models, handles), and a `.tolist()` that raises falls through to it
rather than taking down the reply.

That branch also adds an `iterate_to_list` action, which is the only way to
read **per-atom** properties over the wire — occupancy, altloc and per-atom
b-factor live on atoms, reachable through `cmd.iterate` or `cmd.get_model`, and
`get_model` returns a chempy object that does not survive JSON.

Review note: `iterate_to_list` evaluates a Python expression inside PyMOL. That
is real code execution, but it is no wider than the existing `do` action or
`_resolve_dotted`, and the listener is bound to localhost.

---

## Not yet triaged

Nothing currently. New reports should land here with a symptom and a
reproduction before being promoted to a numbered section.
