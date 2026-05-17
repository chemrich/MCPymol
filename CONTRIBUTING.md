# Contributing to MCPymol

Thanks for the interest! This is a small, opinionated tool; PRs are welcome.

## Dev setup

```bash
git clone https://github.com/chemrich/MCPymol.git
cd MCPymol
uv sync --all-groups
```

## Run the tests

```bash
uv run pytest tests/
```

The suite mocks PyMOL and the socket layer, so you don't need PyMOL installed to run it. If you're touching the bridge protocol (`send_request` framing in `server.py`, the `PyMOLSocketServer` loop in `plugin.py`), please also test against a real PyMOL.

## When you add a new MCP tool

- Give it a docstring an LLM can use. State expected argument vocabularies (valid representation names, palette names, etc.) inline. Mention compatible/related tools when relevant.
- If the tool exists to drive a *view*, add an entry to the view table in `README.md` and include an example prompt.
- Add a test in `tests/test_server.py` (or `test_auto_wrappers.py` for one-shot primitives). Mocking `mcpymol.server.send_request` is the typical pattern.

## Style

- We prefer prose over bullet lists in docstrings.
- We prefer the highest-level tool that does the job — composing five primitives at the LLM layer to do what one well-written view tool already handles is friction we'd like to avoid.

## Protocol notes

The bridge speaks JSON over plain TCP, one request per connection. After sending the request payload the bridge half-closes its write side; the plugin drains, dispatches, writes back, and half-closes. Either side parses JSON incrementally so neither needs to send a length prefix. Keep new fields backward-compatible with `{"action": str, "args": list, "kwargs": dict}`.
