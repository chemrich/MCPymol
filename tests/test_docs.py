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
        for path in SRC.rglob("*.py")
        for match in re.findall(r"MCPYMOL_[A-Z_]+", path.read_text())
    }


def _headline_tools() -> dict[str, str]:
    """Tool name -> module, for every module except the summarised primitives.

    Recursive because tools live in subpackages now (``mcpymol.wiggles``), and
    a glob that stopped at the top level would let a whole family of them ship
    undocumented — the exact failure this file exists to catch.
    """
    tools = {}
    for path in sorted(SRC.rglob("*.py")):
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


def test_readme_tool_counts_match_registered_tools():
    """The total tool count and primitives summary count in README.md must match
    the tools actually registered on the MCP server."""
    import asyncio

    import mcpymol.server  # noqa: F401
    from mcpymol.app import mcp

    total = len(asyncio.run(mcp.list_tools()))
    primitives = total - len(_headline_tools())
    readme = README.read_text()

    assert f"{total} tools in total" in readme, (
        f"README.md does not state the current total tool count '{total} tools in total'"
    )
    assert f"{primitives} PyMOL primitives" in readme, (
        f"README.md does not state the current primitive count '{primitives} PyMOL primitives'"
    )


def test_readme_env_var_defaults_match_code():
    """Numeric defaults in the README Configuration table must match the code."""
    import mcpymol.bridge as bridge
    import mcpymol.plugin as plugin
    import mcpymol.rendering as rendering
    import mcpymol.views as views
    import mcpymol.wiggles.port as wiggles_port

    expected = {
        "MCPYMOL_PORT": "9876",
        "MCPYMOL_AUTOSTART": "1",
        "MCPYMOL_RECV_TIMEOUT": str(int(plugin.RECV_TIMEOUT)),
        "MCPYMOL_MAX_REQUEST_BYTES": str(int(plugin.MAX_REQUEST_BYTES)),
        "MCPYMOL_MAX_OUTPUT_CHARS": str(int(plugin.MAX_OUTPUT_CHARS)),
        "MCPYMOL_DO_START_TIMEOUT": str(int(plugin.DO_START_TIMEOUT)),
        "MCPYMOL_APP_NAP_IDLE_S": str(int(plugin.APP_NAP_IDLE_S)),
        "MCPYMOL_SLOW_OP_TIMEOUT": str(int(bridge._SLOW_OP_TIMEOUT)),
        "MCPYMOL_MAX_IMAGE_BYTES": str(int(rendering.MAX_INLINE_IMAGE_BYTES)),
        "MCPYMOL_PB_TIMEOUT": str(int(views._PB_SUBPROCESS_TIMEOUT)),
        "MCPYMOL_WIGGLES_TIMEOUT": str(int(wiggles_port.DEFAULT_TIMEOUT)),
    }
    readme = README.read_text()
    for var, default in expected.items():
        pattern = rf"\|\s*`{var}`\s*\|\s*`{re.escape(default)}`"
        assert re.search(pattern, readme), (
            f"README Configuration table does not list default `{default}` for `{var}`"
        )


def test_readme_contains_no_relative_image_urls():
    """Relative image links render on GitHub but break on PyPI."""
    readme = README.read_text()
    md_images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", readme)
    html_images = re.findall(r'<img\b[^>]*\bsrc=["\']([^"\']+)["\']', readme, flags=re.I)
    relative = [
        url for url in [*md_images, *html_images] if not url.startswith(("http://", "https://"))
    ]
    assert not relative, f"Relative image URLs in README.md (break on PyPI): {relative}"
