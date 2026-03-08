import json
import os
import socket
from typing import Optional
from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("MCPymol")

HOST = '127.0.0.1'
# Port can be overridden via environment variable, e.g.: MCPYMOL_PORT=9867 uv run mcpymol
PORT = int(os.environ.get("MCPYMOL_PORT", 9876))

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

    # Organic cofactors/ligands: sticks, colored by atom with lightblue carbons
    send_request("show", args=["sticks", f"({name}) and organic"])
    send_request("do", args=[f"util.cbaw ({name}) and organic"])
    send_request("color", args=["lightblue", f"({name}) and organic and elem C"])

    # Inorganic ions and metals: spheres, colored by chemical element
    send_request("show", args=["spheres", f"({name}) and inorganic"])
    send_request("color", args=["atomic", f"({name}) and inorganic"])
    send_request("set", args=["sphere_scale", "0.4", f"({name}) and inorganic"])

    # DNA: brightorange backbone, deepteal ladders
    dna_sel = f"({name}) and resn DA+DC+DG+DT"
    send_request("set", args=["cartoon_nucleic_acid_color", "brightorange", dna_sel])
    send_request("set", args=["cartoon_ladder_color", "deepteal", dna_sel])

    # Center view and set rotation pivot to structure center
    send_request("center", args=[name])
    send_request("do", args=[f"origin {name}"])

def send_request(action: str, args: list = None, kwargs: dict = None) -> dict:
    """Send a JSON request to the PyMOL plugin socket server."""
    payload = {
        "action": action,
        "args": args or [],
        "kwargs": kwargs or {}
    }
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(2.0)
            s.connect((HOST, PORT))
            s.sendall(json.dumps(payload).encode('utf-8'))
            data = s.recv(8192).decode('utf-8')
            return json.loads(data)
    except Exception as e:
        return {"status": "error", "error": f"Socket connection failed: {e}. Is the PyMOL plugin running?"}

@mcp.tool()
def fetch_structure(pdb_code: str, obj_name: Optional[str] = None) -> str:
    """
    Fetches a protein structure from the PDB by its 4-letter code, and applies the multimer heuristic.
    The heuristic automatically removes all chains except the first chain and any other chains
    that are within 5 angstroms of the first chain (preserving associated multimers).
    """
    name = obj_name if obj_name else pdb_code

    send_request("do", args=["reinitialize"])
    send_request("set", args=["mouse_wheel_scale", "0.2"])
    send_request("delete", args=[name])
    res = send_request("fetch", args=[pdb_code, name])
    if res.get("status") == "error":
        return f"Error fetching {pdb_code}: {res.get('error')}"
        
    # Apply heuristic
    chains_res = send_request("get_chains", args=[f"({name})"])
    if chains_res.get("status") == "success":
        chains = chains_res.get("result", [])
        if isinstance(chains, list) and chains:
            first = chains[0]
            keep_sel = f"({name} and chain {first}) or bychain ({name} and chain {first} around 5)"
            send_request("remove", args=[f"({name}) and not ({keep_sel})"])
            send_request("hide", args=["everything", f"({name}) and solvent"])
            _apply_ghost_heart(name)
            return f"Successfully fetched {pdb_code} as '{name}' and applied multimer, solvent, and ghost heart heuristics."

    _apply_ghost_heart(name)
    return f"Successfully fetched {pdb_code} as '{name}' with ghost heart style, but no chains were found to apply heuristic."

@mcp.tool()
def load_structure(file_path: str, obj_name: str) -> str:
    """
    Loads a structure from a local file path and applies the multimer heuristic.
    """
    send_request("do", args=["reinitialize"])
    send_request("set", args=["mouse_wheel_scale", "0.2"])
    send_request("delete", args=[obj_name])
    res = send_request("load", args=[file_path, obj_name])
    if res.get("status") == "error":
        return f"Error loading {file_path}: {res.get('error')}"
        
    # Apply heuristic
    chains_res = send_request("get_chains", args=[f"({obj_name})"])
    if chains_res.get("status") == "success":
        chains = chains_res.get("result", [])
        if isinstance(chains, list) and chains:
            first = chains[0]
            keep_sel = f"({obj_name} and chain {first}) or bychain ({obj_name} and chain {first} around 5)"
            send_request("remove", args=[f"({obj_name}) and not ({keep_sel})"])
            send_request("hide", args=["everything", f"({obj_name}) and solvent"])
            _apply_ghost_heart(obj_name)
            return f"Successfully loaded {file_path} as '{obj_name}' and applied multimer, solvent, and ghost heart heuristics."

    _apply_ghost_heart(obj_name)
    return f"Loaded {file_path} as '{obj_name}' with ghost heart style."

@mcp.tool()
def show(representation: str, selection: str = "all") -> str:
    """Shows a representation for a given selection."""
    res = send_request("show", args=[representation, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Showing {representation} for selection '{selection}'"

@mcp.tool()
def hide(representation: str, selection: str = "all") -> str:
    """Hides a representation for a given selection."""
    res = send_request("hide", args=[representation, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Hiding {representation} for selection '{selection}'"

@mcp.tool()
def color(color_name: str, selection: str = "all") -> str:
    """Sets the color for the specified selection."""
    res = send_request("color", args=[color_name, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Colored selection '{selection}' with {color_name}"
        
@mcp.tool()
def select(name: str, selection: str) -> str:
    """Creates a named selection for subsequent use."""
    res = send_request("select", args=[name, selection])
    if res.get("status") == "error": return res.get("error")
    return f"Created named selection '{name}' for '{selection}'"

@mcp.tool()
def remove(selection: str) -> str:
    """Removes atoms or structures matching the selection."""
    res = send_request("remove", args=[selection])
    if res.get("status") == "error": return res.get("error")
    return f"Removed selection '{selection}'"

@mcp.tool()
def distance(name: str, selection1: str, selection2: str) -> str:
    """Measures and displays the distance between two selections."""
    res = send_request("distance", args=[name, selection1, selection2])
    if res.get("status") == "error": return res.get("error")
    return f"Measured distance between '{selection1}' and '{selection2}' as '{name}'"

@mcp.tool()
def execute_pymol_command(command: str) -> str:
    """Executes a raw PyMOL command string. Use this only for commands not covered by primitive tools."""
    res = send_request("do", args=[command])
    if res.get("status") == "error": return res.get("error")
    return f"Executed command: {command}"


@mcp.tool()
def ligand_view(obj_name: str, ligand_resn: str) -> str:
    """
    Shows a binding-site view focused on a ligand.

    Protein rendered as a semi-transparent cartoon. Pocket residues (within 5Å
    of the ligand) shown as sticks with element coloring and lightblue carbons.
    Ligand shown as thick sticks with yellow carbons. H-bonds drawn as yellow
    dashes. Pocket residues labeled. View zooms to the ligand.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
        ligand_resn: 3-letter residue name of the ligand (e.g. "ATP", "HEM", "LIG")
    """
    lig_sel = f"({obj_name}) and resn {ligand_resn}"
    pocket_sel = f"byres (({obj_name}) and polymer.protein and ({lig_sel} around 5))"

    send_request("hide", args=["everything", obj_name])
    send_request("do", args=["delete hbonds"])

    # Protein as semi-transparent cartoon
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("color", args=["lightblue", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["cartoon_transparency", "0.5", f"({obj_name}) and polymer.protein"])

    # Pocket residues as sticks, element-colored with lightblue carbons
    send_request("show", args=["sticks", pocket_sel])
    send_request("do", args=[f"util.cbaw {pocket_sel}"])
    send_request("color", args=["lightblue", f"({pocket_sel}) and elem C"])

    # Ligand as thick sticks with yellow carbons
    send_request("show", args=["sticks", lig_sel])
    send_request("set", args=["stick_radius", "0.25", lig_sel])
    send_request("do", args=[f"util.cbaw {lig_sel}"])
    send_request("color", args=["yellow", f"({lig_sel}) and elem C"])

    # H-bonds between ligand and pocket (mode=2: polar contacts by geometry)
    send_request("do", args=[f"distance hbonds, ({lig_sel}), ({pocket_sel}), 3.5, 2"])
    send_request("color", args=["yellow", "hbonds"])
    send_request("hide", args=["labels", "hbonds"])
    send_request("set", args=["dash_gap", "0.3", "hbonds"])
    send_request("set", args=["dash_width", "3", "hbonds"])

    # Label pocket residues at CA (one label per residue)
    send_request("label", args=[f"({pocket_sel}) and name CA", '"%s%s" % (resn, resi)'])
    send_request("set", args=["label_color", "white"])
    send_request("set", args=["label_size", "14"])

    send_request("do", args=["bg_color black"])
    send_request("zoom", args=[lig_sel, "8"])
    send_request("do", args=[f"origin {lig_sel}"])

    return f"Showing ligand view for {ligand_resn} in {obj_name}. H-bonds stored as 'hbonds' object."


@mcp.tool()
def bfactor_view(obj_name: str) -> str:
    """
    Colors the structure by crystallographic B-factor (temperature factor).

    Blue = rigid/ordered (low B), white = intermediate, red = flexible/disordered
    (high B). Useful for identifying dynamic loops, disordered termini, and
    rigid structural cores. Shown as cartoon on black background.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", obj_name])
    send_request("do", args=[f"spectrum b, blue_white_red, {obj_name}"])
    send_request("do", args=["bg_color black"])
    send_request("center", args=[obj_name])
    send_request("do", args=[f"origin {obj_name}"])

    return f"Showing B-factor view for {obj_name}: blue=rigid, red=flexible."


@mcp.tool()
def interface_view(obj_name: str, chain_a: str, chain_b: str) -> str:
    """
    Highlights the protein-protein binding interface between two chains.

    Chain A shown in marine blue, chain B in salmon. Interface residues (within
    4Å of the partner chain) shown as a solid surface patch with sticks.
    H-bonds across the interface drawn as yellow dashes.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
        chain_a: First chain ID (e.g. "A")
        chain_b: Second chain ID (e.g. "B")
    """
    sel_a = f"({obj_name}) and chain {chain_a} and polymer.protein"
    sel_b = f"({obj_name}) and chain {chain_b} and polymer.protein"
    iface_a = f"byres ({sel_a} and ({sel_b} around 4))"
    iface_b = f"byres ({sel_b} and ({sel_a} around 4))"

    send_request("hide", args=["everything", obj_name])

    # Both chains as semi-transparent cartoon
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("color", args=["marine", sel_a])
    send_request("color", args=["salmon", sel_b])
    send_request("set", args=["cartoon_transparency", "0.3", f"({obj_name}) and polymer.protein"])

    # Interface surface patches
    send_request("show", args=["surface", iface_a])
    send_request("show", args=["surface", iface_b])
    send_request("color", args=["tv_blue", iface_a])
    send_request("color", args=["tv_red", iface_b])
    send_request("set", args=["transparency", "0.1", iface_a])
    send_request("set", args=["transparency", "0.1", iface_b])

    # Interface residues as sticks — sidechain + CA only (no backbone N, C, O)
    iface_sticks = f"({iface_a} or {iface_b}) and not name N+C+O"
    send_request("show", args=["sticks", iface_sticks])
    send_request("do", args=[f"util.cbaw {iface_a}"])
    send_request("do", args=[f"util.cbaw {iface_b}"])
    send_request("color", args=["tv_blue", f"({iface_a}) and elem C"])
    send_request("color", args=["tv_red", f"({iface_b}) and elem C"])

    # Labels at CA — one per interface residue
    iface_ca = f"({iface_a} or {iface_b}) and name CA"
    send_request("label", args=[iface_ca, '"%s%s" % (resn, resi)'])
    send_request("set", args=["label_color", "white"])
    send_request("set", args=["label_size", "14"])

    # H-bonds across the interface
    send_request("do", args=["delete iface_hbonds"])
    send_request("do", args=[f"distance iface_hbonds, ({sel_a}), ({sel_b}), 3.5, 2"])
    send_request("color", args=["yellow", "iface_hbonds"])
    send_request("hide", args=["labels", "iface_hbonds"])
    send_request("set", args=["dash_gap", "0.3", "iface_hbonds"])
    send_request("set", args=["dash_width", "3", "iface_hbonds"])

    send_request("do", args=["bg_color black"])
    send_request("center", args=[obj_name])
    send_request("do", args=[f"origin {obj_name}"])

    return (
        f"Showing interface between chain {chain_a} (marine/blue) and chain {chain_b} "
        f"(salmon/red) in {obj_name}. Cross-chain H-bonds stored as 'iface_hbonds'."
    )


@mcp.tool(name="as")
def as_tool(representation: str, selection: Optional[str] = "all") -> str:
    """
    Shows one representation while hiding all others for the specified selection
    """
    call_args = []
    if representation is not None: call_args.append(representation)
    if selection is not None: call_args.append(selection)
    
    res = send_request("as", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed as successfully."


@mcp.tool()
def putty_view(obj_name: str) -> str:
    """
    Visualizes protein flexibility using a putty (tube-width) representation.

    The cartoon tube radius scales linearly with crystallographic B-factor:
    thin/blue = rigid/ordered regions, thick/red = flexible/disordered regions.
    A 70%-transparent surface is shown, also colored by B-factor.
    Organic ligands are shown as sticks with yellow carbons. Black background.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon putty, {obj_name}"])
    send_request("spectrum", args=["b", "blue_white_red", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["cartoon_putty_scale_min", "0.3", obj_name])
    send_request("set", args=["cartoon_putty_scale_max", "3.0", obj_name])
    send_request("set", args=["cartoon_putty_transform", "0", obj_name])

    # Transparent surface colored by B-factor
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    send_request("spectrum", args=["b", "blue_white_red", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["transparency", "0.7", obj_name])

    # Organic ligands as sticks with yellow carbons
    send_request("show", args=["sticks", f"({obj_name}) and organic"])
    send_request("do", args=[f"util.cbaw ({obj_name}) and organic"])
    send_request("color", args=["yellow", f"({obj_name}) and organic and elem C"])

    send_request("do", args=["bg_color black"])
    send_request("orient", args=[obj_name])
    return f"Putty view applied to {obj_name}. Tube width and color scale with B-factor (blue=rigid, red=flexible)."


@mcp.tool()
def hydrophobic_surface_view(obj_name: str) -> str:
    """
    Colors the molecular surface by amino acid hydrophobicity.

    Orange = hydrophobic (ILE, VAL, LEU, PHE, MET, ALA, TRP, PRO),
    white = polar (SER, THR, CYS, TYR, ASN, GLN, GLY),
    sky blue = positively charged (ARG, LYS, HIS),
    salmon = negatively charged (ASP, GLU).
    A white cartoon is shown beneath a semi-transparent surface.
    Organic ligands shown as sticks with yellow carbons.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])

    # Color surface by residue hydrophobicity
    send_request("color", args=["orange",  f"({obj_name}) and (resn ILE+VAL+LEU+PHE+MET+ALA+TRP+PRO)"])
    send_request("color", args=["white",   f"({obj_name}) and (resn SER+THR+CYS+TYR+ASN+GLN+GLY)"])
    send_request("color", args=["skyblue", f"({obj_name}) and (resn ARG+LYS+HIS)"])
    send_request("color", args=["salmon",  f"({obj_name}) and (resn ASP+GLU)"])

    # White cartoon visible beneath surface
    send_request("set", args=["cartoon_color", "white", obj_name])
    send_request("set", args=["transparency", "0.15", obj_name])

    # Organic ligands as sticks with yellow carbons
    send_request("show", args=["sticks", f"({obj_name}) and organic"])
    send_request("do", args=[f"util.cbaw ({obj_name}) and organic"])
    send_request("color", args=["yellow", f"({obj_name}) and organic and elem C"])

    send_request("do", args=["bg_color black"])
    send_request("orient", args=[obj_name])
    return f"Hydrophobic surface view applied to {obj_name}. Orange=hydrophobic, white=polar, skyblue=positive, salmon=negative."


@mcp.tool()
def electrostatic_view(obj_name: str) -> str:
    """
    Colors the molecular surface by approximate residue-based electrostatics.

    Charges are assigned by residue pKa: ARG (+1.0), LYS (+0.9), HIS (+0.3),
    ASP (-0.9), GLU (-0.8), all others (0.0). Surface is colored red→white→blue
    via a B-factor spectrum. A white cartoon is shown beneath a semi-transparent
    surface. Organic ligands shown as sticks with yellow carbons.

    For a more accurate Poisson-Boltzmann electrostatic surface, use
    poisson_boltzmann_view (requires APBS and PDB2PQR to be installed).

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])

    # Assign charge values to B-factor by residue pKa
    send_request("do", args=[f"alter ({obj_name}) and polymer.protein, b=0.0"])
    send_request("do", args=[f"alter ({obj_name}) and resn ARG, b=1.0"])
    send_request("do", args=[f"alter ({obj_name}) and resn LYS, b=0.9"])
    send_request("do", args=[f"alter ({obj_name}) and resn HIS, b=0.3"])
    send_request("do", args=[f"alter ({obj_name}) and resn ASP, b=-0.9"])
    send_request("do", args=[f"alter ({obj_name}) and resn GLU, b=-0.8"])
    send_request("rebuild")
    send_request("do", args=[f"spectrum b, red_white_blue, ({obj_name}) and polymer.protein, minimum=-1, maximum=1"])

    # White cartoon visible beneath surface
    send_request("set", args=["cartoon_color", "white", obj_name])
    send_request("set", args=["transparency", "0.15", obj_name])

    # Organic ligands as sticks with yellow carbons
    send_request("show", args=["sticks", f"({obj_name}) and organic"])
    send_request("do", args=[f"util.cbaw ({obj_name}) and organic"])
    send_request("color", args=["yellow", f"({obj_name}) and organic and elem C"])

    send_request("do", args=["bg_color black"])
    send_request("orient", args=[obj_name])
    return f"Electrostatic view applied to {obj_name}. Red=negative, white=neutral, blue=positive (pKa-based approximation)."


@mcp.tool()
def poisson_boltzmann_view(obj_name: str) -> str:
    """
    Colors the molecular surface by true Poisson-Boltzmann electrostatic potential.

    Runs PDB2PQR (AMBER force field, pH 7.0) then APBS to compute the full
    electrostatic potential map. Surface is colored red→white→blue over the
    range ±20 kT/e. A white cartoon is shown beneath a semi-transparent surface.
    Organic ligands shown as sticks with yellow carbons.

    Requires APBS and PDB2PQR to be installed on the system:
        brew install brewsci/bio/apbs
        pip install pdb2pqr

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    import subprocess
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, f"{obj_name}.pdb")
        pqr_path = os.path.join(tmpdir, f"{obj_name}.pqr")
        apbs_in  = os.path.join(tmpdir, f"{obj_name}.in")
        dx_path  = pqr_path + ".dx"

        # Save protein from PyMOL
        send_request("do", args=[f"save {pdb_path}, ({obj_name}) and polymer.protein"])

        # PDB2PQR: assign charges/radii
        result = subprocess.run(
            ["pdb2pqr", "--ff=AMBER", f"--apbs-input={apbs_in}", "--with-ph=7.0",
             pdb_path, pqr_path],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            return f"PDB2PQR failed: {result.stderr[-500:]}"

        # APBS: compute electrostatic potential
        result = subprocess.run(
            ["apbs", apbs_in],
            capture_output=True, text=True, cwd=tmpdir
        )
        if result.returncode != 0:
            return f"APBS failed: {result.stderr[-500:]}"

        if not os.path.exists(dx_path):
            return f"APBS did not produce a .dx map. Check APBS output."

        # Load map and apply to surface
        send_request("do", args=[f"load {dx_path}, {obj_name}_esp_map"])
        send_request("do", args=[f"ramp_new {obj_name}_esp_ramp, {obj_name}_esp_map, [-20, 0, 20], [red, white, blue]"])

        send_request("hide", args=["everything", obj_name])
        send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
        send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
        send_request("do", args=[f"cartoon automatic, {obj_name}"])
        send_request("set", args=["cartoon_color", "white", obj_name])
        send_request("set", args=["surface_color", f"{obj_name}_esp_ramp", obj_name])
        send_request("set", args=["transparency", "0.15", obj_name])

        # Organic ligands as sticks with yellow carbons
        send_request("show", args=["sticks", f"({obj_name}) and organic"])
        send_request("do", args=[f"util.cbaw ({obj_name}) and organic"])
        send_request("color", args=["yellow", f"({obj_name}) and organic and elem C"])

        send_request("do", args=["bg_color black"])
        send_request("orient", args=[obj_name])

    return f"Poisson-Boltzmann electrostatic surface applied to {obj_name}. Red=negative, white=neutral, blue=positive (±20 kT/e)."


@mcp.tool()
def crosslink_view(obj_name: str) -> str:
    """
    Highlights structural cross-links: disulfide bonds, metals, and their coordination.

    Protein backbone shown as a thin grey cartoon. Cysteine side chains (CA→CB→SG)
    shown as yellow sticks, labeled by residue. Disulfide bonds drawn as yellow
    dashes. Metal ions shown as orange spheres. Metal coordination bonds drawn
    as dashed lines to nearby protein atoms. Black background.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])
    send_request("color", args=["grey70", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["cartoon_tube_radius", "0.2", obj_name])

    # Cysteine side chains: CA→CB→SG as yellow sticks
    cys_sc = f"({obj_name}) and resn CYS and (name CA+CB+SG)"
    send_request("show", args=["sticks", cys_sc])
    send_request("color", args=["yellow", cys_sc])

    # Label each CYS at CA
    send_request("label", args=[f"({obj_name}) and resn CYS and name CA", '"%s%s" % (resn, resi)'])
    send_request("set", args=["label_color", "white"])
    send_request("set", args=["label_size", "14"])

    # Disulfide bonds: SG–SG distances ≤ 2.5 Å
    send_request("do", args=[f"delete {obj_name}_disulfides"])
    send_request("do", args=[
        f"distance {obj_name}_disulfides, ({obj_name}) and resn CYS and name SG, "
        f"({obj_name}) and resn CYS and name SG, 2.5"
    ])
    send_request("color", args=["yellow", f"{obj_name}_disulfides"])
    send_request("hide", args=["labels", f"{obj_name}_disulfides"])
    send_request("set", args=["dash_width", "4", f"{obj_name}_disulfides"])
    send_request("set", args=["dash_gap", "0.1", f"{obj_name}_disulfides"])

    # Metal ions as orange spheres
    send_request("show", args=["spheres", f"({obj_name}) and metals"])
    send_request("color", args=["orange", f"({obj_name}) and metals"])
    send_request("set", args=["sphere_scale", "0.5", f"({obj_name}) and metals"])

    # Metal coordination bonds
    send_request("do", args=[f"delete {obj_name}_metalcoord"])
    send_request("do", args=[
        f"distance {obj_name}_metalcoord, ({obj_name}) and metals, "
        f"({obj_name}) and polymer.protein and (name N+O+S), 2.8"
    ])
    send_request("color", args=["orange", f"{obj_name}_metalcoord"])
    send_request("hide", args=["labels", f"{obj_name}_metalcoord"])
    send_request("set", args=["dash_width", "3", f"{obj_name}_metalcoord"])
    send_request("set", args=["dash_gap", "0.2", f"{obj_name}_metalcoord"])

    send_request("do", args=["bg_color black"])
    send_request("orient", args=[obj_name])
    return f"Crosslink view applied to {obj_name}. Yellow=disulfide bonds (CYS), orange=metal coordination."


@mcp.tool()
def set(setting: str, value: str, selection: Optional[str] = None) -> str:
    """
    Sets a PyMOL setting to a specified value
    """
    call_args = []
    if setting is not None: call_args.append(setting)
    if value is not None: call_args.append(value)
    if selection is not None: call_args.append(selection)
    
    res = send_request("set", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed set successfully."


@mcp.tool()
def cartoon(item_type: str, selection: Optional[str] = "all") -> str:
    """
    Sets the cartoon type for the specified selection
    """
    call_args = []
    if item_type is not None: call_args.append(item_type)
    if selection is not None: call_args.append(selection)
    
    res = send_request("cartoon", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed cartoon successfully."


@mcp.tool()
def spectrum(expression: str, palette: Optional[str] = "rainbow", selection: Optional[str] = "all") -> str:
    """
    Colors selection in a spectrum
    """
    call_args = []
    if expression is not None: call_args.append(expression)
    if palette is not None: call_args.append(palette)
    if selection is not None: call_args.append(selection)
    
    res = send_request("spectrum", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed spectrum successfully."


@mcp.tool()
def label(selection: str, expression: Optional[str] = "name") -> str:
    """
    Adds labels to atoms in the selection
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    if expression is not None: call_args.append(expression)
    
    res = send_request("label", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed label successfully."


@mcp.tool()
def angle(name: Optional[str] = None, selection1: Optional[str] = "(pk1)", selection2: Optional[str] = "(pk2)", selection3: Optional[str] = "(pk3)") -> str:
    """
    Measures the angle between three selections
    """
    call_args = []
    if name is not None: call_args.append(name)
    if selection1 is not None: call_args.append(selection1)
    if selection2 is not None: call_args.append(selection2)
    if selection3 is not None: call_args.append(selection3)
    
    res = send_request("angle", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed angle successfully."


@mcp.tool()
def dihedral(name: Optional[str] = None, selection1: Optional[str] = "(pk1)", selection2: Optional[str] = "(pk2)", selection3: Optional[str] = "(pk3)", selection4: Optional[str] = "(pk4)") -> str:
    """
    Measures the dihedral angle between four selections
    """
    call_args = []
    if name is not None: call_args.append(name)
    if selection1 is not None: call_args.append(selection1)
    if selection2 is not None: call_args.append(selection2)
    if selection3 is not None: call_args.append(selection3)
    if selection4 is not None: call_args.append(selection4)
    
    res = send_request("dihedral", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed dihedral successfully."


@mcp.tool()
def center(selection: Optional[str] = "all") -> str:
    """
    Centers the view on a selection
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("center", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed center successfully."


@mcp.tool()
def orient(selection: Optional[str] = "all") -> str:
    """
    Orients the view to align with principal axes of the selection
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("orient", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed orient successfully."


@mcp.tool()
def zoom(selection: Optional[str] = "all", buffer: Optional[str] = "5") -> str:
    """
    Zooms the view on a selection
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    if buffer is not None: call_args.append(buffer)
    
    res = send_request("zoom", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed zoom successfully."


@mcp.tool()
def reset(obj: Optional[str] = None) -> str:
    """
    Resets the view, optionally resetting an object's matrix
    """
    call_args = []
    if obj is not None: call_args.append(obj)
    
    res = send_request("reset", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed reset successfully."


@mcp.tool()
def turn(axis: str, angle: Optional[str] = "90") -> str:
    """
    Rotates the camera around an axis
    """
    call_args = []
    if axis is not None: call_args.append(axis)
    if angle is not None: call_args.append(angle)
    
    res = send_request("turn", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed turn successfully."


@mcp.tool()
def move(axis: str, distance: Optional[str] = "1") -> str:
    """
    Moves the camera along an axis
    """
    call_args = []
    if axis is not None: call_args.append(axis)
    if distance is not None: call_args.append(distance)
    
    res = send_request("move", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed move successfully."


@mcp.tool()
def clip(mode: str, distance: Optional[str] = "1") -> str:
    """
    Adjusts the clipping planes
    """
    call_args = []
    if mode is not None: call_args.append(mode)
    if distance is not None: call_args.append(distance)
    
    res = send_request("clip", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed clip successfully."


@mcp.tool()
def save(filename: str, selection: Optional[str] = "all", state: Optional[str] = "-1") -> str:
    """
    Saves data to a file
    """
    call_args = []
    if filename is not None: call_args.append(filename)
    if selection is not None: call_args.append(selection)
    if state is not None: call_args.append(state)
    
    res = send_request("save", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed save successfully."


@mcp.tool()
def png(filename: str, options: Optional[str] = None) -> str:
    """
    Saves a PNG image
    """
    call_args = []
    if filename is not None: call_args.append(filename)
    if options is not None: call_args.append(options)
    
    res = send_request("png", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed png successfully."


@mcp.tool()
def deselect() -> str:
    """
    Clears the current selection
    """
    call_args = []

    
    res = send_request("deselect", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed deselect successfully."


@mcp.tool()
def create(name: str, selection: Optional[str] = "all", source_state: Optional[str] = "1") -> str:
    """
    Creates a new object from a selection
    """
    call_args = []
    if name is not None: call_args.append(name)
    if selection is not None: call_args.append(selection)
    if source_state is not None: call_args.append(source_state)
    
    res = send_request("create", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed create successfully."


@mcp.tool()
def extract(name: str, selection: Optional[str] = "all") -> str:
    """
    Extracts a selection to a new object
    """
    call_args = []
    if name is not None: call_args.append(name)
    if selection is not None: call_args.append(selection)
    
    res = send_request("extract", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed extract successfully."


@mcp.tool()
def delete(name: str) -> str:
    """
    Deletes objects or selections
    """
    call_args = []
    if name is not None: call_args.append(name)
    
    res = send_request("delete", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed delete successfully."


@mcp.tool()
def align(mobile: str, target: Optional[str] = "all", options: Optional[str] = None) -> str:
    """
    Aligns one selection to another
    """
    call_args = []
    if mobile is not None: call_args.append(mobile)
    if target is not None: call_args.append(target)
    if options is not None: call_args.append(options)
    
    res = send_request("align", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed align successfully."


@mcp.tool()
def super(mobile: str, target: Optional[str] = "all", options: Optional[str] = None) -> str:
    """
    Superimposes one selection onto another
    """
    call_args = []
    if mobile is not None: call_args.append(mobile)
    if target is not None: call_args.append(target)
    if options is not None: call_args.append(options)
    
    res = send_request("super", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed super successfully."


@mcp.tool()
def intra_fit(selection: str) -> str:
    """
    Fits all states within an object
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("intra_fit", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed intra_fit successfully."


@mcp.tool()
def intra_rms(selection: str) -> str:
    """
    Calculates RMSD between states within an object
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("intra_rms", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed intra_rms successfully."


@mcp.tool()
def alter(selection: str, expression: str) -> str:
    """
    Alters atomic properties in a selection
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    if expression is not None: call_args.append(expression)
    
    res = send_request("alter", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed alter successfully."


@mcp.tool()
def alter_state(state: str, selection: str, expression: str) -> str:
    """
    Alters atomic coordinates in a state
    """
    call_args = []
    if state is not None: call_args.append(state)
    if selection is not None: call_args.append(selection)
    if expression is not None: call_args.append(expression)
    
    res = send_request("alter_state", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed alter_state successfully."


@mcp.tool()
def h_add(selection: Optional[str] = "all") -> str:
    """
    Adds hydrogens to a selection
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("h_add", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed h_add successfully."


@mcp.tool()
def h_fill(selection: Optional[str] = "all") -> str:
    """
    Adds hydrogens and adjusts valences
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("h_fill", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed h_fill successfully."


@mcp.tool()
def bond(atom1: str, atom2: str, order: Optional[str] = "1") -> str:
    """
    Creates a bond between two atoms
    """
    call_args = []
    if atom1 is not None: call_args.append(atom1)
    if atom2 is not None: call_args.append(atom2)
    if order is not None: call_args.append(order)
    
    res = send_request("bond", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed bond successfully."


@mcp.tool()
def unbond(atom1: str, atom2: str) -> str:
    """
    Removes a bond between two atoms
    """
    call_args = []
    if atom1 is not None: call_args.append(atom1)
    if atom2 is not None: call_args.append(atom2)
    
    res = send_request("unbond", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed unbond successfully."


@mcp.tool()
def rebuild(selection: Optional[str] = "all") -> str:
    """
    Regenerates all displayed geometry
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("rebuild", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed rebuild successfully."


@mcp.tool()
def refresh() -> str:
    """
    Refreshes the display
    """
    call_args = []

    
    res = send_request("refresh", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed refresh successfully."


@mcp.tool()
def util_cbc(selection: Optional[str] = "all") -> str:
    """
    Colors by chain (Color By Chain)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbc", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbc successfully."


@mcp.tool()
def util_cbaw(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, white carbons (Color By Atom, White)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbaw", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbaw successfully."


@mcp.tool()
def util_cbag(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, green carbons (Color By Atom, Green)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbag", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbag successfully."


@mcp.tool()
def util_cbac(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, cyan carbons (Color By Atom, Cyan)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbac", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbac successfully."


@mcp.tool()
def util_cbam(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, magenta carbons (Color By Atom, Magenta)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbam", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbam successfully."


@mcp.tool()
def util_cbay(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, yellow carbons (Color By Atom, Yellow)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbay", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbay successfully."


@mcp.tool()
def util_cbas(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, salmon carbons (Color By Atom, Salmon)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbas", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbas successfully."


@mcp.tool()
def util_cbab(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, slate carbons (Color By Atom, slateBLue)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbab", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbab successfully."


@mcp.tool()
def util_cbao(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, orange carbons (Color By Atom, Orange)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbao", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbao successfully."


@mcp.tool()
def util_cbap(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, purple carbons (Color By Atom, Purple)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbap", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbap successfully."


@mcp.tool()
def util_cbak(selection: Optional[str] = "all") -> str:
    """
    Colors by atom, pink carbons (Color By Atom, pinK)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.cbak", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.cbak successfully."


@mcp.tool()
def util_chainbow(selection: Optional[str] = "all") -> str:
    """
    Colors chains in rainbow gradient (CHAINs in rainBOW)
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.chainbow", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.chainbow successfully."


@mcp.tool()
def util_rainbow(selection: Optional[str] = "all") -> str:
    """
    Colors residues in rainbow from N to C terminus
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.rainbow", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.rainbow successfully."


@mcp.tool()
def util_ss(selection: Optional[str] = "all") -> str:
    """
    Colors by secondary structure
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.ss", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.ss successfully."


@mcp.tool()
def util_color_by_element(selection: Optional[str] = "all") -> str:
    """
    Colors atoms by their element
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.color_by_element", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.color_by_element successfully."


@mcp.tool()
def util_color_secondary(selection: Optional[str] = "all") -> str:
    """
    Colors secondary structure elements
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("util.color_secondary", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed util.color_secondary successfully."


@mcp.tool()
def spheroid(selection: Optional[str] = "all") -> str:
    """
    Displays atoms as smooth spheres
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    
    res = send_request("spheroid", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed spheroid successfully."


@mcp.tool()
def isomesh(name: str, map_object: str, level: str, selection: Optional[str] = "all") -> str:
    """
    Creates a mesh isosurface
    """
    call_args = []
    if name is not None: call_args.append(name)
    if map_object is not None: call_args.append(map_object)
    if level is not None: call_args.append(level)
    if selection is not None: call_args.append(selection)
    
    res = send_request("isomesh", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed isomesh successfully."


@mcp.tool()
def isosurface(name: str, map_object: str, level: str, selection: Optional[str] = "all") -> str:
    """
    Creates a solid isosurface
    """
    call_args = []
    if name is not None: call_args.append(name)
    if map_object is not None: call_args.append(map_object)
    if level is not None: call_args.append(level)
    if selection is not None: call_args.append(selection)
    
    res = send_request("isosurface", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed isosurface successfully."


@mcp.tool()
def sculpt_activate(obj: str) -> str:
    """
    Activates sculpting mode for an object
    """
    call_args = []
    if obj is not None: call_args.append(obj)
    
    res = send_request("sculpt_activate", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed sculpt_activate successfully."


@mcp.tool()
def sculpt_deactivate(obj: str) -> str:
    """
    Deactivates sculpting mode for an object
    """
    call_args = []
    if obj is not None: call_args.append(obj)
    
    res = send_request("sculpt_deactivate", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed sculpt_deactivate successfully."


@mcp.tool()
def sculpt_iterate(iterations: str, obj: Optional[str] = "all") -> str:
    """
    Performs sculpting iterations
    """
    call_args = []
    if iterations is not None: call_args.append(iterations)
    if obj is not None: call_args.append(obj)
    
    res = send_request("sculpt_iterate", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed sculpt_iterate successfully."


@mcp.tool()
def scene(key: str, action: Optional[str] = "recall") -> str:
    """
    Manages scenes for later recall
    """
    call_args = []
    if key is not None: call_args.append(key)
    if action is not None: call_args.append(action)
    
    res = send_request("scene", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed scene successfully."


@mcp.tool()
def scene_order(scene_list: str) -> str:
    """
    Sets the order of scenes
    """
    call_args = []
    if scene_list is not None: call_args.append(scene_list)
    
    res = send_request("scene_order", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed scene_order successfully."


@mcp.tool()
def mset(specification: str) -> str:
    """
    Defines a sequence of states for movie playback
    """
    call_args = []
    if specification is not None: call_args.append(specification)
    
    res = send_request("mset", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed mset successfully."


@mcp.tool()
def mplay() -> str:
    """
    Starts playing the movie
    """
    call_args = []

    
    res = send_request("mplay", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed mplay successfully."


@mcp.tool()
def mstop() -> str:
    """
    Stops the movie
    """
    call_args = []

    
    res = send_request("mstop", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed mstop successfully."


@mcp.tool()
def frame(frame_number: Optional[str] = None) -> str:
    """
    Sets or queries the current frame
    """
    call_args = []
    if frame_number is not None: call_args.append(frame_number)
    
    res = send_request("frame", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed frame successfully."


@mcp.tool()
def forward() -> str:
    """
    Advances one frame
    """
    call_args = []

    
    res = send_request("forward", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed forward successfully."


@mcp.tool()
def backward() -> str:
    """
    Goes back one frame
    """
    call_args = []

    
    res = send_request("backward", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed backward successfully."


@mcp.tool()
def rock() -> str:
    """
    Toggles a rocking animation
    """
    call_args = []

    
    res = send_request("rock", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed rock successfully."


@mcp.tool()
def ray(width: Optional[str] = None, height: Optional[str] = None) -> str:
    """
    Performs ray-tracing
    """
    call_args = []
    if width is not None: call_args.append(width)
    if height is not None: call_args.append(height)
    
    res = send_request("ray", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed ray successfully."


@mcp.tool()
def draw(width: Optional[str] = None, height: Optional[str] = None) -> str:
    """
    Uses OpenGL renderer (faster but lower quality)
    """
    call_args = []
    if width is not None: call_args.append(width)
    if height is not None: call_args.append(height)
    
    res = send_request("draw", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed draw successfully."


@mcp.tool()
def mpng(prefix: str) -> str:
    """
    Saves a series of PNG images for movie frames
    """
    call_args = []
    if prefix is not None: call_args.append(prefix)
    
    res = send_request("mpng", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed mpng successfully."


@mcp.tool()
def symexp(prefix: str, selection: str, cutoff: Optional[str] = "20", segi: Optional[str] = None) -> str:
    """
    Generates symmetry-related copies
    """
    call_args = []
    if prefix is not None: call_args.append(prefix)
    if selection is not None: call_args.append(selection)
    if cutoff is not None: call_args.append(cutoff)
    if segi is not None: call_args.append(segi)
    
    res = send_request("symexp", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed symexp successfully."


@mcp.tool()
def set_symmetry(selection: str, a: str, b: str, c: str, alpha: str, beta: str, gamma: str) -> str:
    """
    Sets symmetry parameters for an object
    """
    call_args = []
    if selection is not None: call_args.append(selection)
    if a is not None: call_args.append(a)
    if b is not None: call_args.append(b)
    if c is not None: call_args.append(c)
    if alpha is not None: call_args.append(alpha)
    if beta is not None: call_args.append(beta)
    if gamma is not None: call_args.append(gamma)
    
    res = send_request("set_symmetry", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed set_symmetry successfully."


@mcp.tool()
def fab(sequence: str, options: Optional[str] = None) -> str:
    """
    Creates a peptide chain from a sequence
    """
    call_args = []
    if sequence is not None: call_args.append(sequence)
    if options is not None: call_args.append(options)
    
    res = send_request("fab", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed fab successfully."


@mcp.tool()
def fragment(name: str) -> str:
    """
    Loads a molecular fragment
    """
    call_args = []
    if name is not None: call_args.append(name)
    
    res = send_request("fragment", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed fragment successfully."


@mcp.tool()
def full_screen() -> str:
    """
    Toggles fullscreen mode
    """
    call_args = []

    
    res = send_request("full_screen", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed full_screen successfully."


@mcp.tool()
def viewport(width: str, height: str) -> str:
    """
    Sets the viewport size
    """
    call_args = []
    if width is not None: call_args.append(width)
    if height is not None: call_args.append(height)
    
    res = send_request("viewport", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed viewport successfully."


@mcp.tool()
def cd(path: str) -> str:
    """
    Changes the current directory
    """
    call_args = []
    if path is not None: call_args.append(path)
    
    res = send_request("cd", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed cd successfully."


@mcp.tool()
def pwd() -> str:
    """
    Prints the current directory
    """
    call_args = []

    
    res = send_request("pwd", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed pwd successfully."


@mcp.tool()
def ls(path: Optional[str] = None) -> str:
    """
    Lists files in the current directory
    """
    call_args = []
    if path is not None: call_args.append(path)
    
    res = send_request("ls", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed ls successfully."


@mcp.tool()
def system(command: str) -> str:
    """
    Executes a system command
    """
    call_args = []
    if command is not None: call_args.append(command)
    
    res = send_request("system", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed system successfully."


@mcp.tool()
def help(command: Optional[str] = None) -> str:
    """
    Shows help for a command
    """
    call_args = []
    if command is not None: call_args.append(command)
    
    res = send_request("help", args=call_args)
    if res.get("status") == "error": return res.get("error")
    return f"Executed help successfully."


def main():
    mcp.run()

if __name__ == "__main__":
    main()
