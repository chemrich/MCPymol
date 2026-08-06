"""Loading structures and grounding the model in what is loaded.

``fetch_structure`` / ``load_structure`` do the standard prep — biological
assembly, the BFS multimer heuristic, the default style — and the ``list_*``
tools let the model check the session state instead of guessing object names,
chain IDs or ligand codes.
"""

from mcpymol.app import mcp
from mcpymol.bridge import send_request

# Chain–chain contact radius for the BFS multimer heuristic.  8 Å keeps
# sprawling functional assemblies (CRP pentamer, ferritin cage) whole while
# still dropping crystallographic neighbours.  Single source of truth: the
# helper's own default used to be 5.0 while every caller passed 8.0, so the
# signature disagreed with the documented behaviour.
DEFAULT_MULTIMER_CUTOFF = 8.0

# Ordered green shades used by the ghost heart style
_GHOST_HEART_GREENS = ["forest", "limegreen", "chartreuse", "palegreen", "lime", "tv_green"]


def _apply_ghost_heart(name: str):
    """Applies the ghost heart visualization style to an object.

    Ghost heart = cartoon + semi-transparent surface, chains colored in
    shades of green, black background.
    """
    send_request("show", args=["cartoon", name])
    send_request("show", args=["surface", name])
    chains_res = send_request("get_chains", args=[name])
    if chains_res.get("status") == "success":
        chains = chains_res.get("result", [])
        for i, chain in enumerate(chains):
            green = _GHOST_HEART_GREENS[i % len(_GHOST_HEART_GREENS)]
            send_request("color", args=[green, f"{name} and chain {chain} and polymer.protein"])
    send_request("set", args=["transparency", "0.6", name])
    send_request("do", args=["bg_color black"])
    send_request("set", args=["opaque_background", "1"])

    # Organic cofactors/ligands: sticks, colored by atom with lightblue carbons
    send_request("show", args=["sticks", f"({name}) and organic"])
    send_request("do", args=[f"util.cbaw ({name}) and organic"])
    send_request("color", args=["lightblue", f"({name}) and organic and elem C"])

    # Inorganic ions and metals: spheres, colored by chemical element
    send_request("show", args=["spheres", f"({name}) and inorganic"])
    send_request("color", args=["atomic", f"({name}) and inorganic"])
    send_request("set", args=["sphere_scale", "0.3", f"({name}) and inorganic"])

    # Nucleic acids (DNA/RNA): brightorange backbone, deepteal ladders
    na_sel = f"({name}) and polymer.nucleic"
    send_request("set", args=["cartoon_nucleic_acid_color", "brightorange", na_sel])
    send_request("set", args=["cartoon_ladder_color", "deepteal", na_sel])

    # Center view and set rotation pivot to structure center
    send_request("center", args=[name])
    send_request("do", args=[f"origin {name}"])


def _apply_multimer_heuristic(name: str, cutoff: float = DEFAULT_MULTIMER_CUTOFF):
    """BFS expansion to find all connected chains in a multimer."""
    # 1. Get initial chains
    res = send_request("get_chains", args=[name])
    if res.get("status") != "success" or not res.get("result"):
        return

    all_chains = res.get("result", [])  # guard above guarantees this is populated
    kept_chains = {all_chains[0]}

    # 2. Expand until stable
    while True:
        chain_sel = "+".join(list(kept_chains))
        # Find chains nearby the current set
        nearby_res = send_request(
            "get_chains",
            args=[
                f"({name} and not chain {chain_sel}) and bychain (({name} and chain {chain_sel}) around {cutoff})"
            ],
        )

        if nearby_res.get("status") == "success":
            new_chains = [c for c in nearby_res.get("result", []) if c in all_chains]
            if new_chains and not set(new_chains).issubset(kept_chains):
                kept_chains.update(new_chains)
                continue
        break

    # 3. Apply the removal
    final_sel = "+".join(list(kept_chains))
    send_request("remove", args=[f"({name}) and not chain {final_sel}"])
    send_request("hide", args=["everything", f"({name}) and solvent"])


@mcp.tool()
def fetch_structure(
    pdb_code: str, obj_name: str | None = None, multimer_cutoff: float = DEFAULT_MULTIMER_CUTOFF
) -> str:
    """
    Fetches a protein structure from the PDB.
    By default, it attempts to fetch the first biological assembly (multimer),
    and removes any unrelated chains/states that are not part of the primary multimer.

    Args:
        pdb_code: 4-letter PDB code (e.g. "1abc")
        obj_name: Optional custom name for the object in PyMOL
        multimer_cutoff: Distance (A) between chains to keep them in the same multimer.
                         Default 8.0A is suitable for most functional assemblies.
    """
    name = obj_name if obj_name else pdb_code

    send_request("do", args=["reinitialize"])
    send_request("set", args=["mouse_wheel_scale", "0.1"])
    send_request("delete", args=[name])

    # Use standard fetch
    res = send_request("fetch", args=[pdb_code, name])
    if res.get("status") == "error":
        return f"Error fetching {pdb_code}: {res.get('error')}"

    _apply_multimer_heuristic(name, multimer_cutoff)
    _apply_ghost_heart(name)
    send_request("zoom", args=[name])
    return f"Successfully fetched {pdb_code} as '{name}' with ghost heart style and BFS multimer heuristic (cutoff={multimer_cutoff}A)."


@mcp.tool()
def load_structure(
    file_path: str, obj_name: str, multimer_cutoff: float = DEFAULT_MULTIMER_CUTOFF
) -> str:
    """
    Loads a structure from a local file path and applies the BFS multimer heuristic.

    Args:
        file_path: Path to the structure file (PDB, MMCIF, etc.)
        obj_name: Name for the object in PyMOL
        multimer_cutoff: Distance (A) between chains to keep them in the same multimer.
                         Default 8.0A is suitable for most functional assemblies.
    """
    send_request("do", args=["reinitialize"])
    send_request("set", args=["mouse_wheel_scale", "0.1"])
    send_request("delete", args=[obj_name])
    res = send_request("load", args=[file_path, obj_name])
    if res.get("status") == "error":
        return f"Error loading {file_path}: {res.get('error')}"

    _apply_multimer_heuristic(obj_name, multimer_cutoff)
    _apply_ghost_heart(obj_name)
    send_request("zoom", args=[obj_name])
    return f"Successfully loaded {file_path} as '{obj_name}' with ghost heart style and BFS multimer heuristic (cutoff={multimer_cutoff}A)."


# ── Scene introspection ──────────────────────────────────────────────────────


@mcp.tool()
def list_objects() -> str:
    """Lists all loaded PyMOL objects.

    Call this when you don't know what's currently in the session — for
    example before composing a selection or running a view tool that needs
    an ``obj_name``.
    """
    res = send_request("get_object_list", args=["all"])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    objs = res.get("result") or []
    if not objs:
        return "No objects are loaded."
    return "Loaded objects: " + ", ".join(objs)


@mcp.tool()
def list_chains(obj_name: str = "all") -> str:
    """Lists the chain IDs present in an object (or in all objects).

    Useful before calling :func:`interface_view`, :func:`conservation_view`
    or any tool that needs a specific chain ID.
    """
    res = send_request("get_chains", args=[obj_name])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    chains = res.get("result") or []
    if not chains:
        return f"No chains found in '{obj_name}'."
    return f"Chains in '{obj_name}': " + ", ".join(chains)


@mcp.tool()
def list_ligands(obj_name: str) -> str:
    """Lists the small-molecule (organic) ligand residue names in an object.

    Call this before :func:`ligand_view`, :func:`pocket_view`, or
    :func:`pharmacophore_view` when you don't already know the ligand's
    3-letter residue name.
    """
    # We could use `iterate (...) and organic, stored.ligs.add(resn)` and then
    # read `stored.ligs` back, but PyMOL doesn't have a built-in "send me a
    # variable" command, so we'd need a side channel. Cheaper: ask PyMOL to
    # dump the organic atoms as PDB text and parse the resn column.
    fetch = send_request("get_pdbstr", args=[f"({obj_name}) and organic"])
    if fetch.get("status") == "error":
        return fetch.get("error", "Unknown error")
    pdb = fetch.get("result") or ""
    ligs: set[str] = set()
    for line in pdb.splitlines():
        if line.startswith(("HETATM", "ATOM  ")):
            resn = line[17:20].strip()
            if resn:
                ligs.add(resn)
    if not ligs:
        return f"No organic ligands found in '{obj_name}'."
    return f"Ligands in '{obj_name}': " + ", ".join(sorted(ligs))
