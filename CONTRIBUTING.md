# Contributing to MCPymol

Thanks for the interest! This is a small, opinionated tool; PRs are welcome.

## Dev setup

```bash
git clone https://github.com/chemrich/MCPymol.git
cd MCPymol
uv sync --all-groups
uv run pre-commit install   # installs both the commit and push hooks
```

`pre-commit` is in the dev dependency group, so `uv sync --all-groups` provides
it — there is nothing to install globally.

## The hooks

They mirror CI, so a red pipeline is something you find in a second rather than
after a push. `pre-commit install` wires up both stages in one go.

**On commit** (fast — under a second on a small change):

| Hook | What it does |
| --- | --- |
| `ruff check --fix` | lint, auto-fixing what it can |
| `ruff format` | formatting |
| `mypy` | type check over `src/` |
| `uv lock --locked` | blocks a lockfile that has drifted from `pyproject.toml` |
| hygiene | trailing whitespace, final newline, line endings, YAML/TOML validity, merge markers, stray `breakpoint()`, files over 512 KB |

**On push:** the full `pytest` suite. It takes about 16 seconds — too slow to
pay on every commit, unremarkable on a push.

Two things worth knowing:

- The `ruff` and `mypy` hooks run through `uv run`, so they use exactly the
  versions in `uv.lock` — the same ones CI installs. Pinning tool versions
  separately in `.pre-commit-config.yaml` would give the project two sources of
  truth that drift until a hook passes locally and CI fails.
- The auto-fixing hooks **fail the commit when they change a file**. That is
  intended: inspect what changed, `git add` it, and commit again.

CI enforces every one of these independently, so Dependabot PRs and anyone who
skipped `pre-commit install` are still covered. Use `git commit --no-verify` to
bypass the hooks in a pinch — CI will still catch it.

## Run the tests

```bash
uv run pytest tests/
```

The suite mocks PyMOL and the socket layer, so you don't need PyMOL installed to run it. If you're touching the bridge protocol (`send_request` framing in `bridge.py`, the `PyMOLSocketServer` loop in `plugin.py`), please also test against a real PyMOL — and note `tests/test_bridge_roundtrip.py`, which runs the real bridge against the real listener over a socket rather than a mock.

## When you add a new MCP tool

- Give it a docstring an LLM can use. State expected argument vocabularies (valid representation names, palette names, etc.) inline. Mention compatible/related tools when relevant.
- If the tool exists to drive a *view*, add an entry to the view table in `README.md` and include an example prompt.
- Put it in the module that owns its concern: `views.py` for a `*_view` preset, `analysis.py` for something that reports numbers, `primitives.py` for a thin `pymol.cmd` wrapper, `structures.py` for loading and introspection. `server.py` is a re-export facade — add the new name to its import block and `__all__`, which `tests/test_package_layout.py` enforces.
- Document every parameter with `Annotated[str, Field(description="...")]`. FastMCP builds the JSON schema from the signature, so a description in the docstring never reaches the client; tests fail if a parameter has none.
- Add a test in the file matching the module (`tests/test_contacts.py`, `tests/test_rendering.py`, and so on).

**Patch `send_request` where the code looks it up, not on the facade.** Each
module does `from mcpymol.bridge import send_request`, which binds its own
name, so:

```python
@patch("mcpymol.views.send_request")        # correct — the module under test
@patch("mcpymol.server.send_request")       # WRONG — rebinds only the facade
```

Patching `mcpymol.server` leaves the real function in place for every module
that actually calls it, so the test hits a live socket instead of the mock.

## Style

- We prefer prose over bullet lists in docstrings.
- We prefer the highest-level tool that does the job — composing five primitives at the LLM layer to do what one well-written view tool already handles is friction we'd like to avoid.

## Protocol notes

The bridge speaks JSON over plain TCP, one request per connection. After sending the request payload the bridge half-closes its write side; the plugin drains, dispatches, writes back, and half-closes. Either side parses JSON incrementally so neither needs to send a length prefix. Keep new fields backward-compatible with `{"action": str, "args": list, "kwargs": dict}`.
