"""Documentation claims that must stay true.

These are the facts that rot silently: an environment variable added in code
and never documented, or a headline tool that ships without a mention in the
README. Both happened between v1.2.1 and v1.4.0 — six of nine env vars were
undocumented and `fetch_alphafold` was never named — so they are checked here
rather than re-audited by hand each release.

Deliberately *not* checked: the test count and other prose numbers. Requiring
a README edit for every added test would be friction with no payoff.
"""

import ast
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent
README = REPO / "README.md"
SKILL = REPO / ".claude" / "skills" / "mcpymol-guide" / "SKILL.md"
CONTRIBUTING = REPO / "CONTRIBUTING.md"
SRC = REPO / "src" / "mcpymol"

# Tools in primitives.py are thin one-per-pymol.cmd wrappers, covered by the
# README's summary line rather than individually — documenting 87 of them
# would bury the tools that need explaining.
SUMMARISED_MODULE = "primitives.py"


def _env_vars_in_code() -> set[str]:
    return {
        match
        for path in SRC.glob("*.py")
        for match in re.findall(r"MCPYMOL_[A-Z_]+", path.read_text())
    }


def _headline_tools() -> dict[str, str]:
    """Tool name -> module, for every module except the summarised primitives."""
    tools = {}
    for path in sorted(SRC.glob("*.py")):
        if path.name == SUMMARISED_MODULE:
            continue
        for node in ast.parse(path.read_text()).body:
            if not isinstance(node, ast.FunctionDef):
                continue
            if any(
                (isinstance(d, ast.Call) and getattr(d.func, "attr", None) == "tool")
                or (isinstance(d, ast.Attribute) and d.attr == "tool")
                for d in node.decorator_list
            ):
                tools[node.name] = path.name
    return tools


def test_every_env_var_is_documented_in_the_readme():
    """A knob nobody can discover is a knob that does not exist."""
    undocumented = sorted(
        _env_vars_in_code() - set(re.findall(r"MCPYMOL_[A-Z_]+", README.read_text()))
    )

    assert not undocumented, f"env vars in code but not in README: {undocumented}"


def test_the_readme_documents_no_phantom_env_vars():
    """The reverse: a documented variable the code no longer reads."""
    in_code = _env_vars_in_code()
    phantom = sorted(set(re.findall(r"MCPYMOL_[A-Z_]+", README.read_text())) - in_code)

    assert not phantom, f"README documents env vars the code does not read: {phantom}"


def test_every_headline_tool_is_named_in_the_readme():
    """View presets, analysis and rendering tools each need a mention. The
    thin pymol.cmd primitives are covered by the README's summary line."""
    readme = README.read_text()
    missing = sorted(name for name in _headline_tools() if name not in readme)

    assert not missing, f"tools with no mention in README.md: {missing}"


def test_every_headline_tool_is_named_in_the_skill():
    """The skill is what an agent reads before driving the server; a tool it
    does not list is a tool that will not get used."""
    skill = SKILL.read_text()
    missing = sorted(name for name in _headline_tools() if name not in skill)

    assert not missing, f"tools with no mention in the mcpymol-guide skill: {missing}"


@pytest.mark.parametrize("doc", [README, SKILL, CONTRIBUTING])
def test_docs_do_not_describe_server_py_as_holding_the_implementation(doc):
    """server.py became a re-export facade in v1.3.0. Docs that still send a
    reader there to find the bridge or a view send them to the wrong file."""
    text = doc.read_text()
    stale = [
        line.strip()
        for line in text.splitlines()
        if "server.py" in line
        and any(word in line.lower() for word in ("framing", "bridge protocol", "editing "))
    ]

    assert not stale, f"{doc.name} still points at server.py for implementation: {stale}"


def test_contributing_warns_against_patching_the_facade():
    """Patching mcpymol.server.send_request rebinds only the facade's name, so
    the module under test keeps calling the real socket. CONTRIBUTING used to
    recommend exactly that."""
    text = CONTRIBUTING.read_text()

    assert "mcpymol.views.send_request" in text or "module under test" in text, (
        "CONTRIBUTING must explain that send_request is patched on the owning "
        "module, not on mcpymol.server"
    )
