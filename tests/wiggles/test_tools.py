"""The MCP boundary of mcpymol.wiggles: every tool forwards what it was given.

tools.py is thin by design: each tool hands its arguments to a view through
``_report``. Thin code still breaks in one quiet way — a parameter the schema
advertises, the model fills in, and the wrapper never passes on. The view then
runs on its default and reports success. Until now only the live sweep called
these bodies, and it supplies its own arguments, so it could not tell.
"""

from __future__ import annotations

import inspect
import typing
from unittest.mock import MagicMock

import pytest

from mcpymol.wiggles import tools

TOOLS = sorted(
    name
    for name, fn in vars(tools).items()
    if inspect.isfunction(fn) and fn.__module__ == tools.__name__ and not name.startswith("_")
)


def _distinct_value(name: str, annotation: object, index: int) -> object:
    """A value for ``name`` that cannot be mistaken for any other argument or default."""
    base = typing.get_origin(annotation) or annotation
    args = typing.get_args(annotation)
    if base is typing.Union or str(base) in ("types.UnionType", "<class 'types.UnionType'>"):
        non_none = [a for a in args if a is not type(None)]
        annotation = non_none[0] if non_none else str
        base = typing.get_origin(annotation) or annotation
    if base is typing.Annotated:
        return _distinct_value(name, args[0], index)
    if base is bool:
        return "flip"  # resolved against the default below
    if base is int:
        return 7000 + index
    if base is float:
        return 7000.5 + index
    if base is list:
        return [7000.25 + index, 7001.25 + index]
    return f"<{name}>"


def _arguments(fn: typing.Callable[..., str]) -> dict[str, object]:
    hints = typing.get_type_hints(fn, include_extras=False)
    out: dict[str, object] = {}
    for index, (name, param) in enumerate(inspect.signature(fn).parameters.items()):
        value = _distinct_value(name, hints.get(name, str), index)
        if value == "flip":
            value = not param.default if param.default is not inspect.Parameter.empty else True
        out[name] = value
    return out


def test_every_wiggles_tool_is_covered():
    assert len(TOOLS) >= 14, TOOLS


@pytest.mark.parametrize("tool_name", TOOLS)
def test_tool_forwards_every_argument_it_accepts(tool_name, monkeypatch):
    port = object()
    report = MagicMock(return_value="report")
    monkeypatch.setattr(tools, "_port", lambda: port)
    monkeypatch.setattr(tools, "_report", report)
    fn = getattr(tools, tool_name)
    given = _arguments(fn)

    assert fn(**given) == "report"

    report.assert_called_once()
    view, *positional = report.call_args.args
    forwarded = [*positional, *report.call_args.kwargs.values()]
    if "port" in inspect.signature(view).parameters:  # map_info reads a file, no session
        assert port in forwarded, f"{tool_name} did not pass the port"
    dropped = {
        name: value
        for name, value in given.items()
        if not any(value is f or value == f for f in forwarded)
    }
    assert not dropped, f"{tool_name} accepted but never forwarded: {dropped}"
    assert callable(view)
