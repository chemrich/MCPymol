import hashlib
import json
import math
import os
import socket
import time
import urllib.request
import urllib.parse
import urllib.error
from typing import Optional
from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("MCPymol")

HOST = '127.0.0.1'
# Port can be overridden via environment variable, e.g.: MCPYMOL_PORT=9867 uv run mcpymol
PORT = int(os.environ.get("MCPYMOL_PORT", 9876))

# MMseqs2 server for evolutionary conservation (ColabFold public API by default)
MMSEQS_URL = os.environ.get("MCPYMOL_MMSEQS_URL", "https://api.colabfold.com")

# Ordered green shades used by the ghost heart style
_GHOST_HEART_GREENS = ["forest", "limegreen", "chartreuse", "palegreen", "lime", "tv_green"]

# Standard amino acid alphabet for Shannon entropy calculation
_AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

# In-memory cache: MD5(sequence) → list of per-residue entropy values
# Avoids repeat MMseqs2 API calls when only the scale is being changed
_conservation_cache: dict[str, list[float]] = {}


def _run_mmseqs2(sequence: str, server_url: str = None, use_env: bool = True) -> str:
    """Submit a protein sequence to the ColabFold MMseqs2 API and return an A3M MSA.

    Uses a submit → poll → download pattern against the ColabFold public API
    (or a user-specified local server).

    Args:
        sequence: Amino acid sequence string.
        server_url: MMseqs2 API base URL. Defaults to MMSEQS_URL env/global.
        use_env: If True, search environmental databases too (mode "env").

    Returns:
        A3M-formatted multiple sequence alignment as a string.
    """
    host = server_url or MMSEQS_URL
    mode = "env" if use_env else "all"

    # 1. Submit the search job
    data = urllib.parse.urlencode({
        "q": f">query\n{sequence}\n",
        "mode": mode,
    }).encode()
    req = urllib.request.Request(f"{host}/ticket/msa", data=data)
    req.add_header("User-Agent", "mcpymol/1.0.0 (github.com/chemrich/MCPymol)")

    max_retries = 3
    for attempt in range(max_retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                ticket = json.loads(resp.read().decode())
            break
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == max_retries - 1:
                raise RuntimeError(f"Failed to submit MMseqs2 job after {max_retries} attempts: {e}")
            time.sleep(2 ** attempt)

    ticket_id = ticket.get("id")
    if not ticket_id:
        raise RuntimeError(f"MMseqs2 server returned no ticket ID: {ticket}")

    # 2. Poll until complete
    poll_url = f"{host}/ticket/{ticket_id}"
    for _ in range(120):  # up to ~10 minutes with 5s intervals
        time.sleep(5)
        try:
            with urllib.request.urlopen(poll_url, timeout=10) as resp:
                status = json.loads(resp.read().decode())
        except (urllib.error.URLError, TimeoutError):
            continue

        if status.get("status") == "COMPLETE":
            break
        if status.get("status") == "ERROR":
            raise RuntimeError(f"MMseqs2 search failed: {status}")
    else:
        raise RuntimeError("MMseqs2 search timed out after 10 minutes")

    # 3. Download the result
    dl_url = f"{host}/result/download/{ticket_id}"
    with urllib.request.urlopen(dl_url, timeout=30) as resp:
        import tarfile
        import io
        tar_bytes = resp.read()

    # The download is a tar.gz containing .a3m files
    a3m_content = ""
    tar_buf = io.BytesIO(tar_bytes)
    with tarfile.open(fileobj=tar_buf, mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.name.endswith(".a3m"):
                f = tar.extractfile(member)
                if f:
                    a3m_content += f.read().decode("utf-8")

    if not a3m_content:
        raise RuntimeError("No A3M alignment found in MMseqs2 results")

    return a3m_content


def _parse_a3m(a3m_text: str) -> list[list[str]]:
    """Parse an A3M-format MSA into a list of aligned sequences.

    A3M uses lowercase letters for insertions (relative to the query).
    We strip insertions so every sequence aligns column-by-column with
    the query.

    Returns:
        List of sequences, each a list of single-character residues
        aligned to the query. The first entry is the query itself.
    """
    sequences = []
    current = []
    for line in a3m_text.splitlines():
        if line.startswith(">"):
            if current:
                # Strip lowercase (insertions) and join
                seq = [ch for ch in "".join(current) if not ch.islower()]
                sequences.append(seq)
            current = []
        else:
            current.append(line.strip())
    if current:
        seq = [ch for ch in "".join(current) if not ch.islower()]
        sequences.append(seq)
    return sequences


def _compute_shannon_entropy(msa: list[list[str]]) -> list[float]:
    """Compute per-position Shannon entropy from an MSA.

    Lower entropy = more conserved.  The result is normalized to [0, 1]
    where 0 = perfectly conserved, 1 = maximum variability.

    Args:
        msa: List of aligned sequences (each a list of characters).

    Returns:
        List of normalized entropy values, one per query residue position.
    """
    if not msa:
        return []

    n_positions = len(msa[0])
    n_seqs = len(msa)
    max_entropy = math.log2(20)  # theoretical max for 20 amino acids
    entropies = []

    for col in range(n_positions):
        counts: dict[str, int] = {}
        total = 0
        for seq in msa:
            if col < len(seq):
                aa = seq[col].upper()
                if aa in _AA_ALPHABET:
                    counts[aa] = counts.get(aa, 0) + 1
                    total += 1

        if total == 0:
            entropies.append(1.0)
            continue

        entropy = 0.0
        for count in counts.values():
            p = count / total
            if p > 0:
                entropy -= p * math.log2(p)

        # Normalize to [0, 1]
        entropies.append(entropy / max_entropy if max_entropy > 0 else 0.0)

    return entropies

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

def send_request(action: str, args: list = None, kwargs: dict = None, timeout: float = 10.0) -> dict:
    """Send a JSON request to the PyMOL plugin socket server.
    
    Args:
        action: The PyMOL command or custom action to execute.
        args: Positional arguments for the action.
        kwargs: Keyword arguments for the action.
        timeout: Socket timeout in seconds. Defaults to 10.0.
    """
    payload = {
        "action": action,
        "args": args or [],
        "kwargs": kwargs or {}
    }
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(timeout)
            s.connect((HOST, PORT))
            s.sendall(json.dumps(payload).encode('utf-8'))
            data = s.recv(8192).decode('utf-8')
            return json.loads(data)
    except Exception as e:
        return {"status": "error", "error": f"Socket connection failed: {e}. Is the PyMOL plugin running?"}

def _apply_multimer_heuristic(name: str, cutoff: float = 5.0):
    """BFS expansion to find all connected chains in a multimer."""
    # 1. Get initial chains
    res = send_request("get_chains", args=[name])
    if res.get("status") != "success" or not res.get("result"):
        return
    
    all_chains = res.get("result")
    kept_chains = {all_chains[0]}
    
    # 2. Expand until stable
    while True:
        chain_sel = "+".join(list(kept_chains))
        # Find chains nearby the current set
        nearby_res = send_request("get_chains", args=[f"({name} and not chain {chain_sel}) and bychain (({name} and chain {chain_sel}) around {cutoff})"])
        
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
def fetch_structure(pdb_code: str, obj_name: Optional[str] = None, multimer_cutoff: float = 8.0) -> str:
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
def load_structure(file_path: str, obj_name: str, multimer_cutoff: float = 8.0) -> str:
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

    # Ligand as thick sticks with yellow carbons (organic) or spheres (inorganic)
    send_request("show", args=["sticks", f"({lig_sel}) and organic"])
    send_request("set", args=["stick_radius", "0.25", f"({lig_sel}) and organic"])
    send_request("show", args=["spheres", f"({lig_sel}) and inorganic"])
    send_request("set", args=["sphere_scale", "0.3", f"({lig_sel}) and inorganic"])
    
    send_request("do", args=[f"util.cbaw {lig_sel}"])
    send_request("color", args=["yellow", f"({lig_sel}) and organic and elem C"])
    send_request("color", args=["atomic", f"({lig_sel}) and inorganic"])

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
def conservation_view(
    obj_name: str,
    selection: str = "all",
    server_url: Optional[str] = None,
    use_env: bool = True,
    chain: Optional[str] = None,
    scale: str = "relative",
    force_refresh: bool = False,
) -> str:
    """
    Colors the structure by evolutionary conservation using Shannon entropy.

    Runs a full pipeline: extracts the protein sequence from the loaded
    structure, submits it to an MMseqs2 server (ColabFold public API by
    default) for multiple sequence alignment, computes per-residue Shannon
    entropy, and maps the conservation scores onto the structure via the
    B-factor column and spectrum coloring.

    Entropy scores are cached in memory by sequence, so changing the scale
    or re-running on the same protein does not require a repeat API call.

    Magenta/blue = highly conserved (low entropy), white = moderate,
    cyan/green = highly variable (high entropy).

    NOTE: The first call makes an external API call and may take 30 seconds
    to several minutes depending on the server and sequence length.
    Subsequent calls for the same sequence are instant.

    Args:
        obj_name: PyMOL object name (e.g. "1ubq")
        selection: PyMOL selection to analyze (default "all")
        server_url: Override the MMseqs2 server URL (defaults to ColabFold
                    public API, or MCPYMOL_MMSEQS_URL env var)
        use_env: Search environmental databases in addition to UniRef
                 (default True, gives deeper MSAs)
        chain: Specific chain ID to analyze. If None, uses the first
               protein chain found.
        scale: Color scaling mode. "relative" (default) maps the color
               gradient to the actual min/max entropy range of this protein,
               maximizing visual contrast. "absolute" uses the full
               theoretical entropy range (0 to log2(20)), useful when
               comparing conservation across different proteins.
        force_refresh: If True, bypass the cache and re-fetch the MSA from
                       the MMseqs2 server even if scores are cached.
    """
    # 1. Determine which chain to use
    if chain is None:
        chains_res = send_request("get_chains", args=[f"({obj_name}) and polymer.protein"])
        if chains_res.get("status") != "success" or not chains_res.get("result"):
            return f"Error: could not get protein chains from {obj_name}"
        chain = chains_res["result"][0]

    chain_sel = f"({obj_name}) and chain {chain} and polymer.protein"

    # 2. Extract the sequence via FASTA
    fasta_res = send_request("get_fastastr", args=[chain_sel])
    if fasta_res.get("status") != "success" or not fasta_res.get("result"):
        return f"Error: could not extract sequence for chain {chain} of {obj_name}"

    fasta_str = fasta_res["result"]
    # Parse the FASTA: skip header lines, join sequence lines
    seq_lines = [ln for ln in fasta_str.strip().splitlines() if not ln.startswith(">")]
    sequence = "".join(seq_lines).strip()

    if len(sequence) < 10:
        return f"Error: sequence too short ({len(sequence)} residues) for conservation analysis"

    # 3. Check cache; run MMseqs2 only on a miss or force_refresh
    cache_key = hashlib.md5(sequence.encode()).hexdigest()
    cache_hit = not force_refresh and cache_key in _conservation_cache

    if cache_hit:
        entropies = _conservation_cache[cache_key]
        msa_note = "cached"
    else:
        try:
            a3m_text = _run_mmseqs2(sequence, server_url=server_url, use_env=use_env)
        except RuntimeError as e:
            return f"Error running MMseqs2: {e}"

        # 4. Parse the MSA and compute Shannon entropy
        msa = _parse_a3m(a3m_text)
        if len(msa) < 2:
            return f"Warning: MSA contains only {len(msa)} sequence(s). Not enough homologs found for meaningful conservation analysis."

        entropies = _compute_shannon_entropy(msa)
        _conservation_cache[cache_key] = entropies
        msa_note = f"{len(msa)} sequences"

    # 5. Map conservation scores onto B-factor column
    #    We invert entropy so that high B = conserved (for intuitive spectrum coloring)
    n_residues = len(entropies)
    min_entropy = min(entropies)
    max_entropy = max(entropies)

    # Build a PyMOL command to alter B-factors residue by residue
    # First, zero out all B-factors
    send_request("do", args=[f"alter {chain_sel}, b=0"])

    # Get the residue indices from PyMOL to map entropy values correctly
    resi_res = send_request("do", args=[
        f'stored.resi_list = []\n'
        f'iterate {chain_sel} and name CA, stored.resi_list.append(resi)'
    ])

    # Compute per-residue conservation scores based on scaling mode
    # Both modes produce scores in [0, 100] where 100 = most conserved
    entropy_range = max_entropy - min_entropy

    for i, entropy in enumerate(entropies):
        if scale == "relative" and entropy_range > 0:
            # Relative: min_entropy → 100 (most conserved), max_entropy → 0
            conservation_score = (1.0 - (entropy - min_entropy) / entropy_range) * 100.0
        else:
            # Absolute (or relative with zero range, i.e. all positions identical)
            conservation_score = (1.0 - entropy) * 100.0

        # resi is 1-indexed in PDB convention; we use rank-order from the FASTA
        resi_idx = i + 1
        send_request("do", args=[
            f"alter ({chain_sel}) and resi {resi_idx} and name CA, b={conservation_score:.1f}"
        ])
        # Propagate CA B-factor to all atoms in the residue
        send_request("do", args=[
            f"alter ({chain_sel}) and resi {resi_idx} and not name CA, "
            f"b=cmd.get_model('({chain_sel}) and resi {resi_idx} and name CA').atom[0].b "
            f"if cmd.get_model('({chain_sel}) and resi {resi_idx} and name CA').atom else 0"
        ])

    # 6. Apply visualization
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", obj_name])

    # Spectrum: magenta (conserved, high b) → white → cyan (variable, low b)
    send_request("do", args=[f"spectrum b, cyan_white_magenta, {chain_sel}, minimum=0, maximum=100"])

    # Show other chains as gray ghost for context
    other_chains_sel = f"({obj_name}) and polymer.protein and not chain {chain}"
    send_request("color", args=["gray50", other_chains_sel])
    send_request("set", args=["cartoon_transparency", "0.5", other_chains_sel])

    send_request("do", args=["bg_color black"])
    send_request("center", args=[obj_name])
    send_request("do", args=[f"origin {obj_name}"])

    # Rebuild to apply B-factor changes
    send_request("do", args=["rebuild"])

    scale_label = "relative" if scale == "relative" else "absolute"
    return (
        f"Showing conservation view for chain {chain} of {obj_name} ({scale_label} scale). "
        f"Magenta = conserved, white = moderate, cyan = variable. "
        f"MSA: {msa_note}, {n_residues} residue positions scored. "
        f"Entropy range: {min_entropy:.3f} – {max_entropy:.3f}."
    )


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
def electrostatic_view(obj_name: str, mode: str = "atomic") -> str:
    """
    Colors the molecular surface by approximate residue-based electrostatics.

    Surface is colored red→white→blue via a B-factor spectrum. A white cartoon
    is shown beneath a semi-transparent surface. Organic ligands shown as sticks
    with yellow carbons.

    For a more accurate Poisson-Boltzmann electrostatic surface, use
    poisson_boltzmann_view (requires APBS and PDB2PQR to be installed).

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
        mode: Charge assignment strategy.
            "atomic" (default) — charges assigned only to terminal charged atoms
            (e.g. ARG NH1/NH2/NE, LYS NZ, ASP OD1/OD2, GLU OE1/OE2, HIS ND1/NE2).
            Produces localized color at charge centers with natural falloff to white.
            "residue" — charges assigned uniformly to all atoms in each charged residue.
            Produces saturated patches; useful for quickly locating charged regions.
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])

    # Zero out all B-factors first
    send_request("do", args=[f"alter ({obj_name}) and polymer.protein, b=0.0"])

    if mode == "atomic":
        # Assign charges only to the terminal charged atoms (actual charge centers)
        send_request("do", args=[f"alter ({obj_name}) and resn ARG and name NH1+NH2+NE, b=1.0"])
        send_request("do", args=[f"alter ({obj_name}) and resn LYS and name NZ, b=1.0"])
        send_request("do", args=[f"alter ({obj_name}) and resn HIS and name ND1+NE2, b=0.3"])
        send_request("do", args=[f"alter ({obj_name}) and resn ASP and name OD1+OD2, b=-1.0"])
        send_request("do", args=[f"alter ({obj_name}) and resn GLU and name OE1+OE2, b=-1.0"])
    else:
        # Assign charges uniformly to all atoms in each charged residue
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
    return f"Electrostatic view applied to {obj_name} (mode={mode}). Red=negative, white=neutral, blue=positive (pKa-based approximation)."


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
def pocket_view(obj_name: str, resn: str) -> str:
    """
    Visualizes the binding pocket cavity around a ligand as a colored surface.

    The pocket (all residues within 5 Å of the ligand) is shown as a
    semi-transparent surface colored by chemical character: orange=hydrophobic,
    white=polar, skyblue=positive, salmon=negative. Pocket residue sidechains
    are shown as sticks. The ligand is shown as yellow sticks. H-bonds between
    the ligand and pocket are drawn as cyan dashes. The protein backbone is
    shown as a thin grey cartoon for context.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
        resn: Ligand residue name (e.g. "ATP", "LIG", "ANP")
    """
    lig = f"({obj_name}) and resn {resn}"
    pocket_sel = f"({obj_name}) and polymer.protein and byres ({lig} around 5)"

    # Thin grey cartoon for whole protein
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])
    send_request("color", args=["grey60", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["cartoon_tube_radius", "0.2", obj_name])

    # Pocket cavity surface colored by chemical character
    send_request("show", args=["surface", pocket_sel])
    _hydrophobic = "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO"
    _polar       = "GLY+SER+THR+TYR+CYS+ASN+GLN"
    _positive    = "LYS+ARG+HIS"
    _negative    = "ASP+GLU"
    send_request("color", args=["orange", f"({pocket_sel}) and resn {_hydrophobic}"])
    send_request("color", args=["white",   f"({pocket_sel}) and resn {_polar}"])
    send_request("color", args=["skyblue", f"({pocket_sel}) and resn {_positive}"])
    send_request("color", args=["salmon",  f"({pocket_sel}) and resn {_negative}"])
    send_request("set", args=["transparency", "0.25", obj_name])

    # Pocket residue sidechains as sticks (element coloring)
    send_request("show", args=["sticks", f"({pocket_sel}) and not name N+C+O"])
    send_request("do", args=[f"util.cbaw ({pocket_sel})"])
    # Re-apply surface colors after util.cbaw recolored atoms
    send_request("color", args=["orange", f"({pocket_sel}) and resn {_hydrophobic}"])
    send_request("color", args=["white",   f"({pocket_sel}) and resn {_polar}"])
    send_request("color", args=["skyblue", f"({pocket_sel}) and resn {_positive}"])
    send_request("color", args=["salmon",  f"({pocket_sel}) and resn {_negative}"])

    # Labels at CA
    send_request("label", args=[f"({pocket_sel}) and name CA", '"%s%s" % (resn, resi)'])
    send_request("set", args=["label_color", "white"])
    send_request("set", args=["label_size", "12"])

    # Ligand as yellow sticks
    send_request("show", args=["sticks", lig])
    send_request("do", args=[f"util.cbaw {lig}"])
    send_request("color", args=["yellow", f"{lig} and elem C"])

    # H-bonds between ligand and pocket
    send_request("do", args=[f"delete {obj_name}_pocket_hbonds"])
    send_request("do", args=[
        f"distance {obj_name}_pocket_hbonds, ({lig}) and (elem N or elem O), "
        f"({pocket_sel}) and (elem N or elem O), 3.5"
    ])
    send_request("color", args=["cyan", f"{obj_name}_pocket_hbonds"])
    send_request("hide", args=["labels", f"{obj_name}_pocket_hbonds"])
    send_request("set", args=["dash_width", "2.5", f"{obj_name}_pocket_hbonds"])

    send_request("do", args=["bg_color black"])
    send_request("zoom", args=[lig, "6"])
    return (
        f"Pocket view applied: {resn} binding site in {obj_name}. "
        f"Orange=hydrophobic, white=polar, skyblue=positive, salmon=negative. "
        f"Cyan dashes=H-bonds."
    )


@mcp.tool()
def pharmacophore_view(obj_name: str, resn: str) -> str:
    """
    Colors a ligand by pharmacophore feature type.

    The ligand is shown as sticks color-coded by pharmacophore property:
    violet=ring/aromatic carbon, yellow=aliphatic carbon,
    skyblue=nitrogen (H-bond donor/acceptor), salmon=oxygen (H-bond acceptor),
    gold=sulfur, palegreen=halogen (F/Cl/Br/I). H-bonds to protein are shown
    as cyan dashes. Interacting residue sidechains are shown as element-colored
    sticks with CA labels. The pocket is shown as a semi-transparent grey
    surface for cavity context. The protein backbone is shown as a thin grey
    cartoon.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
        resn: Ligand residue name (e.g. "ATP", "LIG", "ANP")
    """
    lig = f"({obj_name}) and resn {resn}"
    pocket_sel = f"({obj_name}) and polymer.protein and byres ({lig} around 5)"

    # Thin grey cartoon
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])
    send_request("color", args=["grey60", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["cartoon_tube_radius", "0.2", obj_name])

    # Pocket semi-transparent surface for cavity context
    send_request("show", args=["surface", pocket_sel])
    send_request("color", args=["grey50", pocket_sel])
    send_request("set", args=["transparency", "0.6", obj_name])

    # Pocket sidechain sticks (element coloring, grey surface kept separate)
    send_request("show", args=["sticks", f"({pocket_sel}) and not name N+C+O"])
    send_request("do", args=[f"util.cbaw ({pocket_sel})"])
    send_request("set", args=["stick_radius", "0.15", pocket_sel])
    # Override surface color to grey after util.cbaw recolored atoms by element
    send_request("set", args=["surface_color", "grey50", pocket_sel])

    # Labels at CA
    send_request("label", args=[f"({pocket_sel}) and name CA", '"%s%s" % (resn, resi)'])
    send_request("set", args=["label_color", "white"])
    send_request("set", args=["label_size", "12"])

    # Ligand sticks
    send_request("show", args=["sticks", lig])
    send_request("set", args=["stick_radius", "0.2", lig])

    # Color by pharmacophore feature type
    # inring catches all ring carbons (PyMOL's 'aromatic' misses some due to bond-order perception)
    send_request("color", args=["violet",    f"{lig} and elem C and inring"])        # ring/aromatic
    send_request("color", args=["yellow",    f"{lig} and elem C and not inring"])    # aliphatic
    send_request("color", args=["skyblue",   f"{lig} and elem N"])                   # H-bond donor/acceptor
    send_request("color", args=["salmon",    f"{lig} and elem O"])                   # H-bond acceptor
    send_request("color", args=["gold",      f"{lig} and elem S"])                   # sulfur
    send_request("color", args=["palegreen", f"{lig} and (elem F or elem Cl or elem Br or elem I)"])  # halogen

    # H-bonds to protein
    send_request("do", args=[f"delete {obj_name}_pharm_hbonds"])
    send_request("do", args=[
        f"distance {obj_name}_pharm_hbonds, ({lig}) and (elem N or elem O or elem F), "
        f"({pocket_sel}) and (elem N or elem O), 3.5"
    ])
    send_request("color", args=["cyan", f"{obj_name}_pharm_hbonds"])
    send_request("hide", args=["labels", f"{obj_name}_pharm_hbonds"])
    send_request("set", args=["dash_width", "2.5", f"{obj_name}_pharm_hbonds"])

    send_request("do", args=["bg_color black"])
    send_request("zoom", args=[lig, "6"])
    return (
        f"Pharmacophore view applied to {resn} in {obj_name}. "
        f"Violet=ring/aromatic, yellow=aliphatic, skyblue=N (donor/acceptor), "
        f"salmon=O (acceptor), gold=S, palegreen=halogen. Cyan dashes=H-bonds."
    )


@mcp.tool()
def mutation_view(obj_name: str, mutations: str) -> str:
    """
    Highlights mutated residues on the protein structure.

    Given a comma-separated list of mutations in standard notation (e.g.
    "A123G,V45L,T200S"), the mutated residues are shown as magenta sticks
    and labeled. Nearby residues (within 4 Å) are shown as thin grey sticks
    for packing context. The protein backbone is shown as a grey cartoon.
    Organic ligands are shown as yellow sticks.

    Mutation format: <wildtype_aa><resi><mutant_aa>, e.g. "A123G" (Ala→Gly
    at position 123). Chain can optionally be prefixed: "A:A123G".

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
        mutations: Comma-separated mutation list (e.g. "A123G,V45L,T200S")
    """
    import re

    mut_list = [m.strip() for m in mutations.split(",")]
    resi_list = []
    parsed = []
    for m in mut_list:
        match = re.search(r'(\d+)', m)
        if match:
            resi_list.append(match.group(1))
            parsed.append(m)

    if not resi_list:
        return f"No valid mutations parsed from: {mutations}. Expected format: A123G,V45L"

    resi_sel = "+".join(resi_list)
    mut_residues = f"({obj_name}) and polymer.protein and resi {resi_sel}"

    # Grey cartoon for whole protein
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])
    send_request("color", args=["grey70", f"({obj_name}) and polymer.protein"])

    # Mutated residues: magenta sticks (sidechain only)
    send_request("show", args=["sticks", f"({mut_residues}) and not name N+C+O"])
    send_request("color", args=["magenta", mut_residues])

    # Labels at CA
    send_request("label", args=[f"({mut_residues}) and name CA", '"%s%s" % (resn, resi)'])
    send_request("set", args=["label_color", "white"])
    send_request("set", args=["label_size", "14"])

    # Context: nearby residues as thin element-colored sticks
    context_sel = (
        f"({obj_name}) and polymer.protein and byres ({mut_residues} around 4) "
        f"and not resi {resi_sel}"
    )
    send_request("show", args=["sticks", f"({context_sel}) and not name N+C+O"])
    send_request("do", args=[f"util.cbaw ({context_sel})"])
    send_request("set", args=["stick_radius", "0.1", context_sel])

    # Organic ligands as yellow sticks
    send_request("show", args=["sticks", f"({obj_name}) and organic"])
    send_request("do", args=[f"util.cbaw ({obj_name}) and organic"])
    send_request("color", args=["yellow", f"({obj_name}) and organic and elem C"])

    send_request("do", args=["bg_color black"])
    send_request("zoom", args=[mut_residues, "8"])
    return f"Mutation view applied to {obj_name}. Magenta = {', '.join(parsed)}."


@mcp.tool()
def textbook_view(obj_name: str) -> str:
    """
    Configures PyMOL for a crisp, cel-shaded illustrative look ("Textbook Illustration").

    This view transforms the structure into a bold, 2D illustrative style with sharp
    black outlines, ideal for presentations or textbook-style diagrams. It hides
    the interior complexities, showing a solid white cartoon and surface with heavy
    black edge contours. Ligands are styled similarly as opaque white sticks with outlines.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    
    # White background for print/textbook style
    send_request("do", args=["bg_color white"])
    
    # Show main structure as white cartoon and surface
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    send_request("color", args=["white", f"({obj_name}) and polymer.protein"])
    
    # Ligands as thick white sticks
    org_sel = f"({obj_name}) and organic"
    send_request("show", args=["sticks", org_sel])
    send_request("color", args=["white", org_sel])
    send_request("set", args=["stick_radius", "0.3", org_sel])
    
    # The "cel shading" effect
    send_request("set", args=["ray_trace_mode", "3"])      # 3 = comic-book style coloring
    send_request("set", args=["ray_trace_depth_factor", "0.4"]) 
    send_request("set", args=["ray_trace_disco_factor", "1.0"])
    
    # Heavy contour lines
    send_request("set", args=["antialias", "2"])
    
    # Improve surface appearance for cel shading
    send_request("set", args=["transparency", "0.0", obj_name])
    send_request("set", args=["surface_quality", "1", obj_name])
    
    send_request("orient", args=[obj_name])
    return f"Textbook Illustration view applied to {obj_name}. Note: the full cel-shaded outline effect requires rendering (use the 'ray' command)."


@mcp.tool()
def cinematic_view(obj_name: str) -> str:
    """
    Configures PyMOL for a depth-cued, cinematic look with dramatic lighting.

    This view emphasizes volume and scale using deep shadows, fog, and depth-cueing.
    The core of the structure emerges from a dark background, making massive
    complexes (like ribosomes or viral capsids) look dramatic and imposing.
    Protein uses standard coloring but with altered material properties.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    # Restore basic representation if not present
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    
    # Dramatic deep black background
    send_request("do", args=["bg_color black"])
    
    # Enable fog and depth cueing
    send_request("set", args=["depth_cue", "1"])
    send_request("set", args=["fog", "1"])
    send_request("set", args=["fog_start", "0.45"])   # Fog starts mid-structure
    send_request("set", args=["fog_color", "black"])
    
    # Cinematic lighting and shadows
    send_request("set", args=["light_count", "2"])
    send_request("set", args=["spec_reflect", "0.3"]) # Slightly glossy
    send_request("set", args=["ray_shadows", "1"])
    send_request("set", args=["ray_shadow_decay_factor", "0.1"])
    send_request("set", args=["ray_shadow_decay_range", "3"])
    
    # Enhance the surface material
    send_request("set", args=["transparency", "0.0", obj_name])
    
    return f"Cinematic view applied to {obj_name}. Fog and depth-cueing enabled. Use the 'ray' command to see the full dramatic shadow effect."


@mcp.tool()
def pointillist_view(obj_name: str) -> str:
    """
    Renders the structure as an artistic, abstract pointillist/starfield cloud.

    The continuous surface is replaced by thousands of individual dots representing
    the solvent-accessible surface, resembling a galaxy or pointillist painting.
    The protein backbone is hidden to emphasize the scattered volume. Ligands
    are shown as bright yellow spheres (stars) embedded in the cloud.

    Args:
        obj_name: PyMOL object name (e.g. "1abc")
    """
    send_request("hide", args=["everything", obj_name])
    send_request("do", args=["bg_color black"])
    
    # The "Starfield" point cloud
    send_request("show", args=["dots", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"]) # Default color recovery
    
    # Increase dot density for the pointillist effect
    send_request("set", args=["dot_density", "4"])
    send_request("set", args=["dot_width", "2"])
    
    # Optional: Light outline of the backbone trace 
    send_request("show", args=["ribbon", f"({obj_name}) and polymer.protein"])
    send_request("set", args=["ribbon_width", "0.5"])
    send_request("color", args=["grey30", f"({obj_name}) and polymer.protein and ribbon"])
    
    # Ligands as bright stars
    org_sel = f"({obj_name}) and organic"
    send_request("show", args=["spheres", org_sel])
    send_request("color", args=["yellow", org_sel])
    send_request("set", args=["sphere_scale", "0.4", org_sel])
    
    send_request("orient", args=[obj_name])
    return f"Pointillist/Starfield view applied to {obj_name}."

@mcp.tool(name="set")
def set_setting(setting: str, value: str, selection: Optional[str] = None) -> str:
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


@mcp.tool(name="super")
def super_tool(mobile: str, target: Optional[str] = "all", options: Optional[str] = None) -> str:
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
