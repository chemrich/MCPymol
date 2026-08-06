"""Loading structures and grounding the model in what is loaded.

``fetch_structure`` / ``load_structure`` do the standard prep — biological
assembly, the BFS multimer heuristic, the default style — and the ``list_*``
tools let the model check the session state instead of guessing object names,
chain IDs or ligand codes.
"""

import os
import re
import tempfile
import urllib.error
import urllib.request

from mcpymol.app import mcp
from mcpymol.bridge import send_request

# AlphaFold DB serves predicted models by UniProt accession.  PyMOL's own
# cmd.fetch only knows the RCSB, so these are downloaded here and loaded from
# disk.
ALPHAFOLD_URL = os.environ.get(
    "MCPYMOL_ALPHAFOLD_URL", "https://alphafold.ebi.ac.uk/files/AF-{acc}-F{frag}-model_v{ver}.cif"
)
DEFAULT_ALPHAFOLD_VERSION = 4

# UniProt accession format (the official regex from uniprot.org). Deliberately
# anchored: a 4-character PDB code can never match, so the two namespaces stay
# distinguishable when routing a bare identifier.
_UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$", re.IGNORECASE
)

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


def _alphafold_accession(identifier: str) -> str | None:
    """Return the UniProt accession if ``identifier`` names an AlphaFold model.

    Accepts what people actually type: ``af-P12345``, ``AF-P12345-F1``, the
    full ``AF-P12345-F1-model_v4`` filename, or a bare accession.  Returns None
    for anything else (notably 4-character PDB codes), so the caller can fall
    through to the RCSB.
    """
    ident = identifier.strip()
    if not ident:
        return None

    upper = ident.upper()
    if upper.startswith("AF-") or upper.startswith("AF_"):
        # AF-P12345-F1-model_v4 → P12345
        acc = re.split(r"[-_]", ident[3:])[0]
        return acc.upper() if _UNIPROT_RE.match(acc) else None

    # A bare accession. PDB codes are exactly four characters and can never
    # match the UniProt pattern, so this does not shadow them.
    return upper if _UNIPROT_RE.match(ident) else None


def _download_alphafold(accession: str, version: int, fragment: int = 1) -> tuple[str | None, str]:
    """Download an AlphaFold model to a temp file. Returns (path, message)."""
    url = ALPHAFOLD_URL.format(acc=accession, frag=fragment, ver=version)
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            data = resp.read()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None, (
                f"No AlphaFold model for '{accession}' (model_v{version}). Check the "
                f"UniProt accession, or try another model_version — the DB is "
                f"versioned and older entries are not always rebuilt."
            )
        return None, f"AlphaFold DB returned HTTP {e.code} for {accession}: {e.reason}"
    except (urllib.error.URLError, TimeoutError) as e:
        return None, f"Could not reach AlphaFold DB ({url}): {e}"

    if not data:
        return None, f"AlphaFold DB returned an empty file for {accession}."

    fd, path = tempfile.mkstemp(suffix=".cif", prefix=f"AF-{accession}-")
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
    return path, ""


@mcp.tool()
def fetch_structure(
    pdb_code: str, obj_name: str | None = None, multimer_cutoff: float = DEFAULT_MULTIMER_CUTOFF
) -> str:
    """
    Fetches a protein structure from the PDB, or a predicted model from AlphaFold DB.

    By default, it attempts to fetch the first biological assembly (multimer),
    and removes any unrelated chains/states that are not part of the primary multimer.

    A UniProt accession or an ``AF-`` prefixed identifier routes to AlphaFold DB
    instead and is coloured by pLDDT confidence — see :func:`fetch_alphafold`.

    Args:
        pdb_code: 4-letter PDB code (e.g. "1abc"), or an AlphaFold identifier
            (e.g. "P69905", "af-P69905").
        obj_name: Optional custom name for the object in PyMOL
        multimer_cutoff: Distance (A) between chains to keep them in the same multimer.
                         Default 8.0A is suitable for most functional assemblies.
    """
    accession = _alphafold_accession(pdb_code)
    if accession is not None:
        return fetch_alphafold(uniprot_id=accession, obj_name=obj_name)

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


@mcp.tool()
def fetch_alphafold(
    uniprot_id: str,
    obj_name: str | None = None,
    model_version: int = DEFAULT_ALPHAFOLD_VERSION,
    fragment: int = 1,
) -> str:
    """
    Fetches a predicted structure from AlphaFold DB by UniProt accession.

    These are *predictions*, not experimental structures, so the model is
    coloured by pLDDT confidence rather than the usual style — dark blue is
    reliable, orange is essentially unmodelled. Read the orange and yellow
    regions as "probably disordered or wrong", not as flexible loops.

    Note that pLDDT rides in the B-factor column, so ``bfactor_view`` and
    ``putty_view`` will mis-colour these models (they assume low = rigid,
    which is backwards for confidence). Use ``plddt_view`` instead.

    Args:
        uniprot_id: UniProt accession (e.g. "P69905" for human haemoglobin
            alpha). An "AF-" prefix is accepted and stripped.
        obj_name: Optional custom name for the object in PyMOL.
        model_version: AlphaFold DB model version. 4 is current.
        fragment: Fragment number for long proteins split across models
            (F1, F2, …). Most entries only have F1.
    """
    accession = _alphafold_accession(uniprot_id) or uniprot_id.strip().upper()
    name = obj_name if obj_name else f"AF_{accession}"

    path, error = _download_alphafold(accession, model_version, fragment)
    if path is None:
        return error

    try:
        send_request("do", args=["reinitialize"])
        send_request("set", args=["mouse_wheel_scale", "0.1"])
        send_request("delete", args=[name])

        res = send_request("load", args=[path, name], timeout=120.0)
        if res.get("status") == "error":
            return f"Error loading AlphaFold model for {accession}: {res.get('error')}"
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass

    # A predicted monomer — the multimer heuristic has nothing to do here.
    from mcpymol.views import plddt_view

    summary = plddt_view(obj_name=name)
    return f"Fetched AlphaFold model AF-{accession}-F{fragment} (v{model_version}) as '{name}'.\n{summary}"


# ── Sessions ─────────────────────────────────────────────────────────────────

# PyMOL decides it is writing a session from the extension, so these are the
# only two that round-trip a whole scene rather than bare coordinates.
_SESSION_SUFFIXES = (".pse", ".pse.gz")


@mcp.tool()
def save_session(filename: str) -> str:
    """
    Saves the entire PyMOL session to a .pse file.

    A session captures everything — every object, selection, representation,
    colour, scene and the camera — so the work can be reopened exactly as it
    was. Use this before experimenting with a scene you would not want to
    rebuild, and to hand a finished figure to a colleague.

    Unlike ``save``, which writes bare coordinates, this preserves the whole
    visual state.

    Args:
        filename: Path to write. A ``.pse`` extension is added if missing.
    """
    path = os.path.abspath(os.path.expanduser(filename.strip()))
    if not path.lower().endswith(_SESSION_SUFFIXES):
        path += ".pse"

    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    res = send_request("save", args=[path], timeout=300.0)
    if res.get("status") == "error":
        return f"Error saving session: {res.get('error')}"

    if not os.path.exists(path):
        return (
            f"PyMOL reported success but no file appeared at {path}. If PyMOL is "
            f"running on a different machine than this bridge, they cannot "
            f"exchange files."
        )
    return f"Saved session to {path} ({os.path.getsize(path) / 1e6:.1f} MB)."


@mcp.tool()
def load_session(filename: str, merge: bool = False) -> str:
    """
    Restores a PyMOL session from a .pse file.

    By default this replaces everything currently loaded, exactly as opening
    the file in PyMOL would. Set ``merge=True`` to add its objects to the
    current session instead, which is how you get two saved scenes side by
    side — though note that objects with the same name will collide.

    Args:
        filename: Path to the .pse file to open.
        merge: Add to the current session rather than replacing it.
    """
    path = os.path.abspath(os.path.expanduser(filename.strip()))
    if not os.path.exists(path):
        return f"Error: no such file: {path}"
    if not path.lower().endswith(_SESSION_SUFFIXES):
        return (
            f"Error: {path} is not a PyMOL session (.pse). Use load_structure "
            f"for coordinate files such as PDB or mmCIF."
        )

    kwargs = {"partial": 1} if merge else {}
    res = send_request("load", args=[path], kwargs=kwargs, timeout=300.0)
    if res.get("status") == "error":
        return f"Error loading session: {res.get('error')}"

    listed = send_request("get_object_list", args=["all"])
    objects = listed.get("result") or [] if listed.get("status") == "success" else []
    how = "Merged" if merge else "Loaded"
    summary = ", ".join(objects) if objects else "no objects"
    return f"{how} session {path}. Objects now in the session: {summary}."


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
