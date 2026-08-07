"""Every tool, called once against a real PyMOL.

    pytest -m live

WARNING: this clears the PyMOL session it connects to. Save anything you care
about first.

Why this exists
---------------
The mocked suite asserts the payload we *send*. That is worth having, but it
cannot see whether PyMOL will accept it — and every wiring bug this project has
shipped lived in exactly that blind spot:

* seven tools pointing at ``cmd`` functions that do not exist, or binding
  arguments in the wrong order (#15);
* ``render(ray_trace=False)`` writing a blank image;
* ``turntable`` defaulting to that same blank path;
* ``fetch_alphafold`` building a URL that stopped resolving in 2024.

Each one passed the whole mocked suite and failed on first contact with PyMOL.

What it asserts
---------------
Deliberately **not** "no tool errors". Plenty of tools error honestly without
the right context: ``isomesh`` needs a map loaded, ``symexp`` needs crystal
symmetry, ``h_fill`` needs an atom picked in the GUI. Failing on those would
make the suite noise, and noise gets switched off.

Instead it fails only on *wiring* errors — the response signatures that can
only mean the tool is misconfigured rather than misapplied:

* ``Unknown action or method not found`` — the name does not resolve
* ``invalid literal for int()`` — an argument landed in a numeric slot
* ``missing N required positional argument`` — arguments did not reach it
* ``unexpected keyword argument`` / ``takes N positional arguments``

A domain error passes. A misconfigured wrapper cannot.
"""

import socket

import pytest

import mcpymol.bridge as bridge

# Importing mcpymol.server is what registers the tools: mcpymol.app only holds
# the empty FastMCP instance, so importing it alone yields a sweep of nothing.
import mcpymol.server  # noqa: F401
from mcpymol.app import mcp
from mcpymol.bridge import send_request

pytestmark = pytest.mark.live

# One object every smoke check operates on, built rather than fetched so the
# sweep needs no network and is byte-identical every run.
LIVE_OBJECT = "_live_probe"
LIVE_OBJECT_2 = "_live_probe_2"
LIVE_SEQUENCE = "ACDEFGHIKLMNPQRSTVWY"


# ── the signatures that can only mean a misconfigured wrapper ────────────────

WIRING_FAILURES = (
    "unknown action or method not found",
    "invalid literal for int()",
    "required positional argument",
    "unexpected keyword argument",
    "positional arguments but",
    "takes no arguments",
)


def assert_wired(tool: str, result: str) -> None:
    lowered = (result or "").lower()
    for signature in WIRING_FAILURES:
        if signature in lowered:
            pytest.fail(
                f"{tool} is misconfigured, not misapplied — PyMOL said:\n  {result}\n"
                f"This is the class of bug the mocked suite cannot see: the payload "
                f"is sent faithfully, PyMOL just cannot accept it."
            )


# ── connection ───────────────────────────────────────────────────────────────


def _bridge_reachable() -> bool:
    try:
        with socket.create_connection((bridge.HOST, bridge.PORT), timeout=2):
            return True
    except OSError:
        return False


@pytest.fixture(scope="session", autouse=True)
def live_session():
    """Skip the whole module unless a PyMOL plugin is listening."""
    if not _bridge_reachable():
        pytest.skip(
            f"no PyMOL plugin on {bridge.HOST}:{bridge.PORT} — start PyMOL and "
            f"`run src/mcpymol/plugin.py`",
            allow_module_level=True,
        )
    send_request("do", args=["reinitialize"])
    # fab builds a real peptide with chains, residues and CA atoms, without a
    # network round trip or a PDB mirror to depend on.
    send_request("fab", args=[LIVE_SEQUENCE, LIVE_OBJECT])
    send_request("create", args=[LIVE_OBJECT_2, LIVE_OBJECT])
    send_request("do", args=["set_symmetry " + LIVE_OBJECT + ", 50, 50, 50, 90, 90, 90, P 1"])
    _require_a_healthy_viewport()
    yield
    send_request("do", args=["reinitialize"])


def _require_a_healthy_viewport() -> None:
    """Refuse to run against a PyMOL that cannot frame an object.

    A session that has been driven hard can end up with a camera that no zoom
    recovers, and every render then comes back as one dark pixel. Failing forty
    tests on that is noise; saying so once, and pointing at the fix, is not.
    """
    from mcpymol.rendering import _looks_blank, render

    send_request("do", args=["reinitialize"])
    send_request("fab", args=[LIVE_SEQUENCE, LIVE_OBJECT])
    send_request("show", args=["cartoon", LIVE_OBJECT])
    send_request("zoom", args=[LIVE_OBJECT])

    probe = render(width=200, height=150)
    # A framed 200x150 cartoon runs to several KB. A flat frame, or a structure
    # rendered as one distant pixel, comes back around 1 KB — _looks_blank
    # catches the first, this catches the second.
    unusable = isinstance(probe, str) or _looks_blank(probe.data) or len(probe.data) < 4000
    if unusable:
        pytest.skip(
            "PyMOL is loaded but renders an empty frame — its camera is not "
            "responding to zoom. Restart PyMOL, re-run the plugin, and try "
            "again. (Seen after a session has been driven very hard.)",
            allow_module_level=True,
        )
    send_request("create", args=[LIVE_OBJECT_2, LIVE_OBJECT])


@pytest.fixture(autouse=True)
def _keep_the_probe_alive():
    """Some tools delete or mangle the probe; rebuild it if one did."""
    yield
    res = send_request("count_atoms", args=[f"({LIVE_OBJECT})"])
    try:
        alive = int(res.get("result"))
    except (TypeError, ValueError):
        alive = 0
    if not alive:
        send_request("delete", args=[LIVE_OBJECT])
        send_request("fab", args=[LIVE_SEQUENCE, LIVE_OBJECT])


# ── argument defaults, by parameter name ─────────────────────────────────────
#
# Keyed on the parameter name so a newly added tool is swept automatically —
# the completeness test below fails if one appears that this cannot supply.

DEFAULTS: dict[str, str] = {
    "selection": LIVE_OBJECT,
    "selection1": f"{LIVE_OBJECT} and resi 1",
    "selection2": f"{LIVE_OBJECT} and resi 2",
    "selection3": f"{LIVE_OBJECT} and resi 3",
    "selection4": f"{LIVE_OBJECT} and resi 4",
    "obj_name": LIVE_OBJECT,
    "obj": LIVE_OBJECT,
    "mobile": LIVE_OBJECT,
    "target": LIVE_OBJECT_2,
    "name": "_live_made",
    "representation": "cartoon",
    "color_name": "red",
    "setting": "sphere_scale",
    "value": "0.5",
    "item_type": "automatic",
    "expression": "b=1.0",
    "palette": "rainbow",
    "axis": "y",
    "angle": "5",
    "distance": "1",
    "mode": "near",
    "state": "1",
    "source_state": "1",
    "buffer": "5",
    "width": "120",
    "height": "90",
    "level": "1.0",
    "cutoff": "5",
    "order": "1",
    "iterations": "1",
    "key": "F1",
    "action": "store",
    "scene_list": "F1",
    "specification": "1 x2",
    "frame_number": "1",
    "sequence": "AG",
    "command": "echo mcpymol-live-probe",
    "atom1": f"{LIVE_OBJECT} and resi 1 and name CA",
    "atom2": f"{LIVE_OBJECT} and resi 2 and name CA",
    "a": "50",
    "b": "50",
    "c": "50",
    "alpha": "90",
    "beta": "90",
    "gamma": "90",
    "chain_a": "A",
    "chain_b": "A",
    "pdb_code": "1ubq",
    "uniprot_id": "P69905",
    "file_path": "/nonexistent/live-probe.pdb",
    "map_object": "_live_missing_map",
    "prefix": "_live_sym",
    "segi": None,
    "options": None,
    "path": None,
    "filename": None,  # filled per-test with a tmp path
    "out_dir": None,
    "merge": False,
    "replace": False,
    "include_water": False,
    "ray_trace": True,
    "force_refresh": False,
    "use_env": False,
    "scale": "relative",
    "method": "super",
    "max_deviation": None,
    "max_pairs": 3,
    "max_residues": 3,
    "max_residues_listed": 3,
    "frames": 2,
    "spine_radius": 0.9,
    "multimer_cutoff": 8.0,
    "model_version": None,
    "fragment": 1,
    "chain": None,
    "server_url": None,
    "pdb_id": None,
    "groups": f"probe={LIVE_OBJECT}",
    "voxel_pitch": 0.7,
    "poisson_depth": 6,
    "mutations": "A1G",
    "ligand_resn": "ALA",
    "resn": "ALA",
}

# Tools left out of the sweep, each for a stated reason. The completeness test
# treats this as the only permitted way to not be covered.
EXCLUDED: dict[str, str] = {
    "conservation_view": "submits an MSA job to ColabFold; minutes, and hammers a public API",
    "poisson_boltzmann_view": "shells out to apbs/pdb2pqr, which may not be installed",
    "print_export": "mesh reconstruction takes minutes and needs the print extra",
    "fetch_alphafold": "network; covered by the network-marked contract test",
    "fetch_structure": "network; exercised by the journey tests below",
    "full_screen": "toggles the PyMOL window into fullscreen, which is not recoverable headlessly",
    "cd": "changes PyMOL's working directory, which would misdirect later file writes",
    "turntable": "ray-traces a frame sequence; exercised by the journey tests",
}


def _tool_functions() -> dict[str, object]:
    """Registered tool name -> the callable behind it."""
    import asyncio

    return {t.name: mcp._tool_manager.get_tool(t.name).fn for t in asyncio.run(mcp.list_tools())}


def _base_type(annotation):
    """The concrete type behind `int`, `str | None`, `Annotated[float, ...]`."""
    import typing

    for candidate in typing.get_args(annotation) or (annotation,):
        if candidate is not type(None):
            return candidate
    return str


def _coerce(value, annotation):
    """Match the declared type. The tools are typed now — handing "120" to an
    int parameter fails inside the tool, which would look like a wiring bug and
    is really a harness bug."""
    target = _base_type(annotation)
    if target in (int, float, bool) and isinstance(value, str):
        return target(value)
    if target is str and not isinstance(value, str):
        return str(value)
    return value


def _arguments_for(fn, tmp_path) -> dict:
    """Build a plausible argument set from the parameter names and types."""
    import inspect
    import typing

    hints = typing.get_type_hints(fn, include_extras=False)
    args = {}
    for name, param in inspect.signature(fn).parameters.items():
        if name in ("filename",):
            args[name] = str(tmp_path / "live_probe.png")
        elif name in ("out_dir", "path"):
            args[name] = str(tmp_path)
        elif name in DEFAULTS:
            value = DEFAULTS[name]
            if value is not None:
                args[name] = _coerce(value, hints.get(name, str))
        elif param.default is not inspect.Parameter.empty:
            continue  # optional and unknown: let the default stand
        else:
            pytest.fail(
                f"no live-test argument for required parameter '{name}'. Add one to "
                f"DEFAULTS in this file, or exclude the tool with a reason."
            )
    return args


# ── completeness ─────────────────────────────────────────────────────────────


def test_every_tool_is_swept_or_excluded():
    """A tool added without either is a tool nobody ever calls for real."""
    uncovered = sorted(set(_tool_functions()) - set(EXCLUDED))
    # Everything not excluded is parametrized below; this asserts the split is
    # exhaustive rather than that a particular tool is present.
    assert uncovered, "the sweep is empty — tool discovery is broken"
    assert not (set(EXCLUDED) - set(_tool_functions())), (
        f"EXCLUDED names tools that no longer exist: "
        f"{sorted(set(EXCLUDED) - set(_tool_functions()))}"
    )


# ── the sweep ────────────────────────────────────────────────────────────────


def _sweepable():
    import asyncio

    names = sorted(t.name for t in asyncio.run(mcp.list_tools()))
    return [n for n in names if n not in EXCLUDED]


@pytest.mark.parametrize("tool_name", _sweepable())
def test_tool_is_wired_correctly(tool_name, tmp_path):
    """Call the tool for real and reject only misconfiguration."""
    fn = _tool_functions()[tool_name]
    result = fn(**_arguments_for(fn, tmp_path))
    assert_wired(tool_name, str(result))


# ── journeys: real workflows, real numbers ───────────────────────────────────
#
# The sweep proves every tool is reachable. These prove a few of them are
# right, by checking results a wrong implementation could not produce.


@pytest.fixture
def hiv_protease():
    """A clean session with 1HSG framed.

    The sweep leaves the camera wherever its last `turn`/`clip`/`zoom` put it,
    so the journeys reset rather than inherit it.
    """
    from mcpymol.structures import fetch_structure

    send_request("do", args=["reinitialize"])
    result = fetch_structure(pdb_code="1hsg")
    if "Successfully fetched" not in result:
        pytest.skip(f"could not fetch 1hsg: {result}")
    return "1hsg"


def test_journey_render_returns_a_real_image(hiv_protease):
    """The whole point of render: pixels come back, and they are not blank."""
    from mcp.server.fastmcp import Image

    from mcpymol.rendering import _looks_blank, render

    send_request("orient", args=[hiv_protease])
    send_request("zoom", args=[hiv_protease])
    result = render(width=400, height=300)

    assert isinstance(result, Image), result
    assert not _looks_blank(result.data)


def test_journey_contacts_find_the_catalytic_aspartates(hiv_protease):
    """Both copies of Asp25 hydrogen bond the inhibitor. A contact finder that
    cannot see that is not working, whatever it returns."""
    from mcpymol.analysis import contact_report

    result = contact_report(
        selection1="1hsg and resn MK1", selection2="1hsg and polymer", max_pairs=40
    )

    assert "ASP25" in result
    assert "hydrogen bond" in result


def test_journey_structure_info_reads_the_entry(hiv_protease):
    from mcpymol.structures import structure_info

    result = structure_info(obj_name="1hsg")

    assert "MK1" in result
    assert "Chains (2)" in result


def test_journey_sequence_reports_numbering(hiv_protease):
    from mcpymol.structures import get_sequence

    result = get_sequence(obj_name="1hsg", chain="A")

    assert "modelled residues" in result
    assert "PQITLW" in result.upper().replace("\n", "")


def test_journey_interface_area_is_physically_plausible():
    """Barnase-barstar buries ~780 A^2 per side. Anything wildly off means the
    free/bound bookkeeping is wrong, which no mocked test can tell."""
    from mcpymol.analysis import interface_report
    from mcpymol.structures import fetch_structure

    if "Successfully fetched" not in fetch_structure(pdb_code="1brs"):
        pytest.skip("could not fetch 1brs")

    result = interface_report(obj_name="1brs", chain_a="A", chain_b="D")

    area = float(result.split("Buried surface area: ")[1].split(" A^2")[0].replace(",", ""))
    assert 500 < area < 1100, f"implausible interface area: {area} A^2\n{result}"


def test_journey_superposition_separates_core_from_lid():
    """Adenylate kinase open vs closed: the core fits tightly while the LID
    swings ~20 A. One number cannot show that; this view exists to."""
    from mcpymol.comparison import superposition_view
    from mcpymol.structures import fetch_structure

    if "Successfully fetched" not in fetch_structure(pdb_code="4ake"):
        pytest.skip("could not fetch 4ake")
    if "Successfully fetched" not in fetch_structure(pdb_code="1ake", replace=False):
        pytest.skip("could not fetch 1ake")

    send_request("do", args=["create open_a, 4ake and chain A and polymer"])
    send_request("do", args=["create closed_a, 1ake and chain A and polymer"])

    result = superposition_view(mobile="open_a", target="closed_a")

    rmsd = float(result.split("RMSD ")[1].split(" A")[0])
    largest = float(result.split("max ")[1].split(" A")[0])
    assert rmsd < 4.0, f"core should fit tightly, got {rmsd} A\n{result}"
    assert largest > 12.0, f"the LID should move a long way, got {largest} A\n{result}"


def test_journey_a_failed_fetch_keeps_the_session(hiv_protease):
    """#15: a fetch that returns nothing used to take the session with it."""
    from mcpymol.structures import fetch_structure, list_objects

    before = list_objects()
    result = fetch_structure(pdb_code="zzzz")

    assert "Successfully fetched" not in result
    assert "1hsg" in list_objects(), f"the session was cleared by a failed fetch\n{before}"
