"""Static code and protocol invariants checked on every test run.

Enforces five structural properties across ``src/mcpymol/``:
1. Imported packages are stdlib, ``mcpymol``, declared in ``pyproject.toml``,
   or ``pymol`` (inside ``plugin.py``, which itself imports only stdlib + ``pymol``).
2. ``get_object_list`` is called only inside ``structures._describe_objects``.
3. No f-string interpolates a bare ``chain {x}`` without ``chain_clause`` / ``quote``.
4. Discarded ``send_request(...)`` calls never exceed the recorded baseline.
5. Every bridge action in ``src/mcpymol/`` and every ``primitives.py`` wrapper
   matches the captured PyMOL 3.1.0 ``cmd`` / ``util`` signatures, and
   ``READ_ONLY_TOOLS`` only invoke read-only bridge actions.
"""

from __future__ import annotations

import ast
import collections
import json
import pathlib
import re
import sys
import tomllib

REPO = pathlib.Path(__file__).resolve().parent.parent
SRC = REPO / "src" / "mcpymol"
FIXTURES = REPO / "tests" / "fixtures"

# PyPI distribution name -> top-level import name where they differ.
_DIST_TO_IMPORT = {
    "pillow": "PIL",
    "scikit-image": "skimage",
}

_CUSTOM_PLUGIN_ACTIONS = frozenset({"do", "iterate_to_list"})

_READ_ONLY_BRIDGE_ACTIONS = frozenset(
    {
        "count_atoms",
        "get_area",
        "get_chains",
        "get_fastastr",
        "get_frame",
        "get_model",
        "get_names",
        "get_object_list",
        "get_pdbstr",
        "get_setting_float",
        "get_setting_int",
        "get_setting_text",
        "get_symmetry",
        "get_type",
        "get_version",
        "get_view",
        "help",
        "intra_rms",
        "iterate_to_list",
        "ls",
        "pwd",
        "rms_cur",
    }
)


def _declared_imports() -> tuple[set[str], set[str]]:
    """Return (core_imports, optional_imports) declared in pyproject.toml."""
    data = tomllib.loads((REPO / "pyproject.toml").read_text())
    project = data["project"]

    def _parse_req(req: str) -> str:
        dist = re.split(r"[\s\[>=<!~;]", req, maxsplit=1)[0].strip().lower()
        return _DIST_TO_IMPORT.get(dist, dist.replace("-", "_"))

    core = {_parse_req(r) for r in project.get("dependencies", [])}
    optional = {
        _parse_req(r) for group in project.get("optional-dependencies", {}).values() for r in group
    }
    return core, optional


def test_every_import_in_src_is_declared_or_stdlib():
    """No undeclared third-party imports in src/mcpymol/, and plugin.py stays
    self-contained (stdlib + pymol only) so it runs inside PyMOL's interpreter."""
    core, optional = _declared_imports()
    allowed = set(sys.stdlib_module_names) | {"mcpymol"} | core | optional
    plugin_allowed = set(sys.stdlib_module_names) | {"pymol"}

    undeclared: list[str] = []
    plugin_violations: list[str] = []

    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(REPO).as_posix()
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            tops: list[str] = []
            if isinstance(node, ast.Import):
                tops = [alias.name.split(".", 1)[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                tops = [node.module.split(".", 1)[0]]

            for top in tops:
                if path.name == "plugin.py":
                    if top not in plugin_allowed:
                        plugin_violations.append(f"{rel}:{node.lineno} imports {top!r}")
                elif top not in allowed:
                    undeclared.append(f"{rel}:{node.lineno} imports {top!r}")

    assert not undeclared, f"Undeclared imports in src/mcpymol: {undeclared}"
    assert not plugin_violations, f"plugin.py must only import stdlib or pymol: {plugin_violations}"


# Allowed callers of get_object_list:
# - structures._describe_objects: fallback when get_names("public_objects") is unavailable
# - printing._resolve_group: queries which molecule objects contain atoms matching a selection
_ALLOWED_GET_OBJECT_LIST_CALLERS = frozenset(
    {
        ("structures.py", "_describe_objects"),
        ("printing.py", "_resolve_group"),
    }
)


def test_get_object_list_is_confined_to_allowed_callers():
    """Session-wide object enumeration must go through get_names / _describe_objects,
    not raw get_object_list calls."""
    violations: list[str] = []
    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(REPO).as_posix()
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for sub in ast.walk(node):
                    if (
                        isinstance(sub, ast.Constant)
                        and sub.value == "get_object_list"
                        and (path.name, node.name) not in _ALLOWED_GET_OBJECT_LIST_CALLERS
                    ):
                        violations.append(f"{rel}:{sub.lineno} in {node.name}")

    assert not violations, f"get_object_list called outside allowed helpers: {violations}"


def test_no_unquoted_chain_interpolation_in_fstrings():
    """Selection f-strings must not interpolate a raw variable after ``and chain ``
    or ``(chain ``; use ``chain_clause(chain)``, ``residue_selection(...)``, or
    ``quote(chain)``."""
    selection_chain = re.compile(r"(?:\band|\()\s*chain\s+$")
    violations: list[str] = []

    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(REPO).as_posix()
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.JoinedStr) or not node.values:
                continue
            first = node.values[0]
            if (
                isinstance(first, ast.Constant)
                and isinstance(first.value, str)
                and first.value.lstrip()[:1].isupper()
            ):
                # English prose message ("Showing ...", "Error ...", "Chain ...")
                continue
            for idx, part in enumerate(node.values[:-1]):
                if not (
                    isinstance(part, ast.Constant)
                    and isinstance(part.value, str)
                    and selection_chain.search(part.value)
                ):
                    continue
                next_part = node.values[idx + 1]
                if not isinstance(next_part, ast.FormattedValue):
                    continue
                expr = next_part.value
                is_quoted = (
                    isinstance(expr, ast.Call)
                    and isinstance(expr.func, ast.Name)
                    and expr.func.id in ("quote", "chain_clause", "residue_selection")
                )
                if not is_quoted:
                    violations.append(f"{rel}:{node.lineno}")

    assert not violations, (
        f"Unquoted 'chain {{...}}' selection f-string interpolation (use chain_clause/quote): {violations}"
    )


def _collect_discarded_send_requests() -> dict[str, int]:
    counts: collections.Counter[str] = collections.Counter()

    def visit_stmts(stmts: list[ast.stmt], rel: str, fn_name: str) -> None:
        for stmt in stmts:
            if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                nested = f"{fn_name}.{stmt.name}" if fn_name else stmt.name
                visit_stmts(stmt.body, rel, nested)
                continue
            if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call):
                call = stmt.value
                fn = call.func
                name = (
                    fn.id
                    if isinstance(fn, ast.Name)
                    else (fn.attr if isinstance(fn, ast.Attribute) else None)
                )
                if name == "send_request" and call.args and isinstance(call.args[0], ast.Constant):
                    action = str(call.args[0].value)
                    scope = fn_name or "<module>"
                    counts[f"{rel}::{scope}::{action}"] += 1
            for field in ("body", "orelse", "finalbody", "handlers", "cases"):
                children = getattr(stmt, field, None)
                if isinstance(children, list):
                    for child in children:
                        if isinstance(child, (ast.ExceptHandler, ast.match_case)):
                            visit_stmts(child.body, rel, fn_name)
                        elif isinstance(child, ast.stmt):
                            visit_stmts([child], rel, fn_name)

    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(REPO).as_posix()
        tree = ast.parse(path.read_text())
        visit_stmts(tree.body, rel, "")

    return dict(sorted(counts.items()))


def test_discarded_send_request_calls_do_not_exceed_baseline():
    """New code must inspect send_request() responses rather than discarding them."""
    baseline_path = FIXTURES / "discarded_send_request_baseline.json"
    baseline: dict[str, int] = json.loads(baseline_path.read_text())
    current = _collect_discarded_send_requests()

    new_sites = sorted(set(current) - set(baseline))
    increased = {
        k: (baseline[k], current[k])
        for k in set(current) & set(baseline)
        if current[k] > baseline[k]
    }

    assert not new_sites, f"New discarded send_request sites added: {new_sites}"
    assert not increased, f"Discarded send_request count increased (baseline, current): {increased}"


def test_all_bridge_actions_and_primitives_match_pymol_signatures():
    """Every action sent over the bridge exists in PyMOL 3.1.0, and every
    wrapper in primitives.py passes a valid positional argument count."""
    sigs: dict[str, list[dict[str, object]]] = json.loads(
        (FIXTURES / "pymol_3_1_signatures.json").read_text()
    )

    unknown_actions: list[str] = []
    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(REPO).as_posix()
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                fn = node.func
                name = (
                    fn.id
                    if isinstance(fn, ast.Name)
                    else (fn.attr if isinstance(fn, ast.Attribute) else None)
                )
                if (
                    name in ("send_request", "_call", "request")
                    and node.args
                    and isinstance(node.args[0], ast.Constant)
                    and isinstance(node.args[0].value, str)
                ):
                    act = node.args[0].value
                    if act not in _CUSTOM_PLUGIN_ACTIONS and act not in sigs:
                        unknown_actions.append(f"{rel}:{node.lineno} -> {act!r}")

    assert not unknown_actions, (
        f"Bridge actions missing from pymol_3_1_signatures.json: {unknown_actions}"
    )

    prim_tree = ast.parse((SRC / "primitives.py").read_text())
    arity_errors: list[str] = []
    for node in prim_tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        params_ast = node.args.args
        defaults_ast = node.args.defaults
        offset = len(params_ast) - len(defaults_ast)
        param_defaults: dict[str, object] = {}
        for i, p in enumerate(params_ast):
            if i >= offset:
                d = defaults_ast[i - offset]
                param_defaults[p.arg] = d.value if isinstance(d, ast.Constant) else "NON_NONE"
            else:
                param_defaults[p.arg] = "REQUIRED"

        for sub in ast.walk(node):
            if not isinstance(sub, ast.Call):
                continue
            fn = sub.func
            name = fn.id if isinstance(fn, ast.Name) else None
            if name == "_call" and sub.args and isinstance(sub.args[0], ast.Constant):
                act = str(sub.args[0].value)
                kw_names = [kw.arg for kw in sub.keywords if kw.arg and not kw.arg.startswith("_")]
                sent = list(kw_names)
                while sent and param_defaults.get(sent[-1]) is None and node.name != "frame":
                    sent.pop()

                pymol_params = sigs[act]
                pos_params = [
                    p
                    for p in pymol_params
                    if p["kind"] in ("POSITIONAL_ONLY", "POSITIONAL_OR_KEYWORD")
                ]
                has_varargs = any(p["kind"] == "VAR_POSITIONAL" for p in pymol_params)
                req_count = sum(1 for p in pos_params if p["required"])

                if not has_varargs and len(kw_names) > len(pos_params):
                    arity_errors.append(
                        f"{node.name} ({act}): passes up to {len(kw_names)} args, "
                        f"PyMOL accepts {len(pos_params)}"
                    )
                if len(sent) < req_count:
                    arity_errors.append(
                        f"{node.name} ({act}): sends as few as {len(sent)} args, "
                        f"PyMOL requires {req_count}"
                    )

    assert not arity_errors, f"Primitive signature mismatches: {arity_errors}"


def test_read_only_tools_only_invoke_read_only_bridge_actions():
    """Tools marked readOnlyHint=True in app.py must not call mutating actions."""
    from mcpymol.app import READ_ONLY_TOOLS

    tool_actions: dict[str, set[str]] = collections.defaultdict(set)
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text())
        for node in tree.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            tool_name = node.name
            for dec in node.decorator_list:
                if isinstance(dec, ast.Call) and getattr(dec.func, "attr", None) == "tool":
                    for kw in dec.keywords:
                        if kw.arg == "name" and isinstance(kw.value, ast.Constant):
                            tool_name = str(kw.value.value)
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call):
                    fn = sub.func
                    name = (
                        fn.id
                        if isinstance(fn, ast.Name)
                        else (fn.attr if isinstance(fn, ast.Attribute) else None)
                    )
                    if (
                        name in ("send_request", "_call", "request")
                        and sub.args
                        and isinstance(sub.args[0], ast.Constant)
                        and isinstance(sub.args[0].value, str)
                    ):
                        tool_actions[tool_name].add(sub.args[0].value)

    violations: dict[str, list[str]] = {}
    for tool in sorted(READ_ONLY_TOOLS):
        non_ro = sorted(tool_actions.get(tool, set()) - _READ_ONLY_BRIDGE_ACTIONS)
        if non_ro:
            violations[tool] = non_ro

    assert not violations, (
        f"READ_ONLY_TOOLS invoke non-read-only bridge actions directly: {violations}"
    )
