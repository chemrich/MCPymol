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

## Not yet triaged

Nothing currently. New reports should land here with a symptom and a
reproduction before being promoted to a numbered section.
