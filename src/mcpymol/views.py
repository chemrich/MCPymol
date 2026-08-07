"""High-level view presets.

Each of these is a whole scene in one call — representation, colouring,
transparency, H-bonds, labels and camera — so the model does not have to
compose a dozen primitives to get a usable picture.
"""

import os
from typing import Annotated

from pydantic import Field

from mcpymol.app import mcp
from mcpymol.bridge import send_request
from mcpymol.pdbtext import parse_atoms

# Wall-clock ceiling for the external APBS/PDB2PQR binaries.  Without one a
# wedged solver hangs the whole MCP server, since tool calls are synchronous.
_PB_SUBPROCESS_TIMEOUT = float(os.environ.get("MCPYMOL_PB_TIMEOUT", 600.0))


@mcp.tool()
def ligand_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
    ligand_resn: Annotated[
        str, Field(description='3-letter residue name of the ligand (e.g. "ATP", "HEM", "LIG")')
    ],
) -> str:
    """
    Shows a binding-site view focused on a ligand.

    Protein rendered as a semi-transparent cartoon. Pocket residues (within 5Å
    of the ligand) shown as sticks with element coloring and lightblue carbons.
    Ligand shown as thick sticks with yellow carbons. H-bonds drawn as yellow
    dashes. Pocket residues labeled. View zooms to the ligand.
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

    return (
        f"Showing ligand view for {ligand_resn} in {obj_name}. H-bonds stored as 'hbonds' object."
    )


@mcp.tool()
def bfactor_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Colors the structure by crystallographic B-factor (temperature factor).

    Blue = rigid/ordered (low B), white = intermediate, red = flexible/disordered
    (high B). Useful for identifying dynamic loops, disordered termini, and
    rigid structural cores. Shown as cartoon on black background.
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", obj_name])
    send_request("do", args=[f"spectrum b, blue_white_red, {obj_name}"])
    send_request("do", args=["bg_color black"])
    send_request("center", args=[obj_name])
    send_request("do", args=[f"origin {obj_name}"])

    return f"Showing B-factor view for {obj_name}: blue=rigid, red=flexible."


# AlphaFold's published pLDDT confidence bands and their official colours.
# Ordered low → high; each entry is (name, lower_bound, rgb, label).
_PLDDT_BANDS = [
    ("plddt_very_low", 0.0, (1.000, 0.490, 0.271), "very low (<50)"),
    ("plddt_low", 50.0, (1.000, 0.859, 0.075), "low (50-70)"),
    ("plddt_confident", 70.0, (0.396, 0.796, 0.953), "confident (70-90)"),
    ("plddt_very_high", 90.0, (0.051, 0.341, 0.827), "very high (>90)"),
]


def _read_ca_bfactors(obj_name: str) -> list[float] | None:
    """Per-residue B-factors (one per CA), read back via a PDB dump.

    PyMOL has no "send me a variable" command, so the same trick
    ``list_ligands`` uses applies: ask for the atoms as PDB text and parse the
    fixed-width columns.
    """
    res = send_request("get_pdbstr", args=[f"({obj_name}) and name CA"], timeout=60.0)
    if res.get("status") == "error":
        return None
    return [a.bfactor for a in parse_atoms(res.get("result") or "", ca_only=True)]


@mcp.tool()
def plddt_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "AF_P69905")')],
) -> str:
    """
    Colors an AlphaFold model by pLDDT confidence, using the official palette.

    Dark blue = very high (>90), light blue = confident (70-90), yellow = low
    (50-70), orange = very low (<50). Regions below 70 are usually intrinsically
    disordered or simply unreliable — read them as "the model does not know",
    not as a flexible loop.

    AlphaFold stores pLDDT in the B-factor column, which is why ``bfactor_view``
    and ``putty_view`` get these models backwards: they assume low = rigid,
    whereas low pLDDT = low confidence. Use this instead for predicted models.
    """
    scores = _read_ca_bfactors(obj_name)
    if scores is None:
        return f"Error: could not read B-factors from '{obj_name}'."
    if not scores:
        return f"Error: no residues found in '{obj_name}'."

    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", obj_name])

    for i, (color_name, lower, (r, g, b), _label) in enumerate(_PLDDT_BANDS):
        send_request("do", args=[f"set_color {color_name}, [{r}, {g}, {b}]"])
        # Paint low-to-high; each band's upper edge is the next band's floor,
        # and the last band is open-ended.
        sel = f"({obj_name}) and b > {lower}"
        if i + 1 < len(_PLDDT_BANDS):
            sel += f" and b < {_PLDDT_BANDS[i + 1][1]}"
        send_request("color", args=[color_name, sel])

    send_request("do", args=["bg_color black"])
    send_request("center", args=[obj_name])
    send_request("do", args=[f"origin {obj_name}"])

    # A confidence breakdown is the number people actually want from a
    # predicted model, and it is cheap now that the scores are already here.
    total = len(scores)
    counts = []
    for i, (_name, lower, _rgb, label) in enumerate(_PLDDT_BANDS):
        upper = _PLDDT_BANDS[i + 1][1] if i + 1 < len(_PLDDT_BANDS) else 101.0
        n = sum(1 for s in scores if lower <= s < upper)
        counts.append(f"{label}: {100.0 * n / total:.0f}%")

    mean = sum(scores) / total
    warning = ""
    if max(scores) > 100.0 or mean == 0.0:
        warning = (
            " NOTE: these B-factors do not look like pLDDT (expected 0–100); "
            "this may be an experimental structure, where bfactor_view is the "
            "right tool."
        )

    return (
        f"Showing pLDDT confidence for {obj_name} over {total} residues "
        f"(mean {mean:.1f}). " + ", ".join(counts) + "." + warning
    )


@mcp.tool()
def interface_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
    chain_a: Annotated[str, Field(description='First chain ID (e.g. "A")')],
    chain_b: Annotated[str, Field(description='Second chain ID (e.g. "B")')],
) -> str:
    """
    Highlights the protein-protein binding interface between two chains.

    Chain A shown in marine blue, chain B in salmon. Interface residues (within
    4Å of the partner chain) shown as a solid surface patch with sticks.
    H-bonds across the interface drawn as yellow dashes.
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


@mcp.tool()
def putty_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Visualizes protein flexibility using a putty (tube-width) representation.

    The cartoon tube radius scales linearly with crystallographic B-factor:
    thin/blue = rigid/ordered regions, thick/red = flexible/disordered regions.
    A 70%-transparent surface is shown, also colored by B-factor.
    Organic ligands are shown as sticks with yellow carbons. Black background.
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
def hydrophobic_surface_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Colors the molecular surface by amino acid hydrophobicity.

    Orange = hydrophobic (ILE, VAL, LEU, PHE, MET, ALA, TRP, PRO),
    white = polar (SER, THR, CYS, TYR, ASN, GLN, GLY),
    sky blue = positively charged (ARG, LYS, HIS),
    salmon = negatively charged (ASP, GLU).
    A white cartoon is shown beneath a semi-transparent surface.
    Organic ligands shown as sticks with yellow carbons.
    """
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])

    # Color surface by residue hydrophobicity
    send_request(
        "color", args=["orange", f"({obj_name}) and (resn ILE+VAL+LEU+PHE+MET+ALA+TRP+PRO)"]
    )
    send_request("color", args=["white", f"({obj_name}) and (resn SER+THR+CYS+TYR+ASN+GLN+GLY)"])
    send_request("color", args=["skyblue", f"({obj_name}) and (resn ARG+LYS+HIS)"])
    send_request("color", args=["salmon", f"({obj_name}) and (resn ASP+GLU)"])

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
def electrostatic_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
    mode: Annotated[
        str,
        Field(
            description='Charge assignment strategy. "atomic" (default) — charges assigned only to terminal charged atoms (e.g. ARG NH1/NH2/NE, LYS NZ, ASP OD1/OD2, GLU OE1/OE2, HIS ND1/NE2). Produces localized color at charge centers with natural falloff to white. "residue" — charges assigned uniformly to all atoms in each charged residue. Produces saturated patches; useful for quickly locating charged regions.'
        ),
    ] = "atomic",
) -> str:
    """
    Colors the molecular surface by approximate residue-based electrostatics.

    Surface is colored red→white→blue via a B-factor spectrum. A white cartoon
    is shown beneath a semi-transparent surface. Organic ligands shown as sticks
    with yellow carbons.

    For a more accurate Poisson-Boltzmann electrostatic surface, use
    poisson_boltzmann_view (requires APBS and PDB2PQR to be installed).
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
    send_request(
        "do",
        args=[
            f"spectrum b, red_white_blue, ({obj_name}) and polymer.protein, minimum=-1, maximum=1"
        ],
    )

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
def poisson_boltzmann_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Colors the molecular surface by true Poisson-Boltzmann electrostatic potential.

    Runs PDB2PQR (AMBER force field, pH 7.0) then APBS to compute the full
    electrostatic potential map. Surface is colored red→white→blue over the
    range ±20 kT/e. A white cartoon is shown beneath a semi-transparent surface.
    Organic ligands shown as sticks with yellow carbons.

    Requires APBS and PDB2PQR to be installed on the system:
        brew install brewsci/bio/apbs
        pip install pdb2pqr
    """
    import os
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, f"{obj_name}.pdb")
        pqr_path = os.path.join(tmpdir, f"{obj_name}.pqr")
        apbs_in = os.path.join(tmpdir, f"{obj_name}.in")
        dx_path = pqr_path + ".dx"

        # Save protein from PyMOL.  Writing a big assembly can outrun the
        # default socket timeout, so allow longer and check the result — an
        # unnoticed failure here surfaces as a baffling PDB2PQR error about a
        # file that was never written.
        save_res = send_request(
            "do", args=[f"save {pdb_path}, ({obj_name}) and polymer.protein"], timeout=120.0
        )
        if save_res.get("status") == "error":
            return f"Error saving {obj_name} for electrostatics: {save_res.get('error')}"
        if not os.path.exists(pdb_path):
            return (
                f"PyMOL did not write {pdb_path}. The bridge and PyMOL must run on the "
                f"same machine for poisson_boltzmann_view to exchange files."
            )

        # PDB2PQR: assign charges/radii
        try:
            result = subprocess.run(
                [
                    "pdb2pqr",
                    "--ff=AMBER",
                    f"--apbs-input={apbs_in}",
                    "--with-ph=7.0",
                    pdb_path,
                    pqr_path,
                ],
                capture_output=True,
                text=True,
                timeout=_PB_SUBPROCESS_TIMEOUT,
            )
        except FileNotFoundError:
            return "PDB2PQR not found on PATH. Install it with: pip install pdb2pqr"
        except subprocess.TimeoutExpired:
            return f"PDB2PQR timed out after {_PB_SUBPROCESS_TIMEOUT}s."
        if result.returncode != 0:
            return f"PDB2PQR failed: {result.stderr[-500:]}"

        # APBS: compute electrostatic potential
        try:
            result = subprocess.run(
                ["apbs", apbs_in],
                capture_output=True,
                text=True,
                cwd=tmpdir,
                timeout=_PB_SUBPROCESS_TIMEOUT,
            )
        except FileNotFoundError:
            return "APBS not found on PATH. Install it with: brew install brewsci/bio/apbs"
        except subprocess.TimeoutExpired:
            return (
                f"APBS timed out after {_PB_SUBPROCESS_TIMEOUT}s. Large structures can "
                f"exceed this; raise MCPYMOL_PB_TIMEOUT to allow longer."
            )
        if result.returncode != 0:
            return f"APBS failed: {result.stderr[-500:]}"

        if not os.path.exists(dx_path):
            return "APBS did not produce a .dx map. Check APBS output."

        # Load map and apply to surface
        send_request("do", args=[f"load {dx_path}, {obj_name}_esp_map"])
        send_request(
            "do",
            args=[
                f"ramp_new {obj_name}_esp_ramp, {obj_name}_esp_map, [-20, 0, 20], [red, white, blue]"
            ],
        )

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
def crosslink_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Highlights structural cross-links: disulfide bonds, metals, and their coordination.

    Protein backbone shown as a thin grey cartoon. Cysteine side chains (CA→CB→SG)
    shown as yellow sticks, labeled by residue. Disulfide bonds drawn as yellow
    dashes. Metal ions shown as orange spheres. Metal coordination bonds drawn
    as dashed lines to nearby protein atoms. Black background.
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
    send_request(
        "do",
        args=[
            f"distance {obj_name}_disulfides, ({obj_name}) and resn CYS and name SG, "
            f"({obj_name}) and resn CYS and name SG, 2.5"
        ],
    )
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
    send_request(
        "do",
        args=[
            f"distance {obj_name}_metalcoord, ({obj_name}) and metals, "
            f"({obj_name}) and polymer.protein and (name N+O+S), 2.8"
        ],
    )
    send_request("color", args=["orange", f"{obj_name}_metalcoord"])
    send_request("hide", args=["labels", f"{obj_name}_metalcoord"])
    send_request("set", args=["dash_width", "3", f"{obj_name}_metalcoord"])
    send_request("set", args=["dash_gap", "0.2", f"{obj_name}_metalcoord"])

    send_request("do", args=["bg_color black"])
    send_request("orient", args=[obj_name])
    return f"Crosslink view applied to {obj_name}. Yellow=disulfide bonds (CYS), orange=metal coordination."


@mcp.tool()
def pocket_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
    resn: Annotated[str, Field(description='Ligand residue name (e.g. "ATP", "LIG", "ANP")')],
) -> str:
    """
    Visualizes the binding pocket cavity around a ligand as a colored surface.

    The pocket (all residues within 5 Å of the ligand) is shown as a
    semi-transparent surface colored by chemical character: orange=hydrophobic,
    white=polar, skyblue=positive, salmon=negative. Pocket residue sidechains
    are shown as sticks. The ligand is shown as yellow sticks. H-bonds between
    the ligand and pocket are drawn as cyan dashes. The protein backbone is
    shown as a thin grey cartoon for context.
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
    _polar = "GLY+SER+THR+TYR+CYS+ASN+GLN"
    _positive = "LYS+ARG+HIS"
    _negative = "ASP+GLU"
    send_request("color", args=["orange", f"({pocket_sel}) and resn {_hydrophobic}"])
    send_request("color", args=["white", f"({pocket_sel}) and resn {_polar}"])
    send_request("color", args=["skyblue", f"({pocket_sel}) and resn {_positive}"])
    send_request("color", args=["salmon", f"({pocket_sel}) and resn {_negative}"])
    send_request("set", args=["transparency", "0.25", obj_name])

    # Pocket residue sidechains as sticks (element coloring)
    send_request("show", args=["sticks", f"({pocket_sel}) and not name N+C+O"])
    send_request("do", args=[f"util.cbaw ({pocket_sel})"])
    # Re-apply surface colors after util.cbaw recolored atoms
    send_request("color", args=["orange", f"({pocket_sel}) and resn {_hydrophobic}"])
    send_request("color", args=["white", f"({pocket_sel}) and resn {_polar}"])
    send_request("color", args=["skyblue", f"({pocket_sel}) and resn {_positive}"])
    send_request("color", args=["salmon", f"({pocket_sel}) and resn {_negative}"])

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
    send_request(
        "do",
        args=[
            f"distance {obj_name}_pocket_hbonds, ({lig}) and (elem N or elem O), "
            f"({pocket_sel}) and (elem N or elem O), 3.5"
        ],
    )
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
def pharmacophore_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
    resn: Annotated[str, Field(description='Ligand residue name (e.g. "ATP", "LIG", "ANP")')],
) -> str:
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
    send_request("color", args=["violet", f"{lig} and elem C and inring"])  # ring/aromatic
    send_request("color", args=["yellow", f"{lig} and elem C and not inring"])  # aliphatic
    send_request("color", args=["skyblue", f"{lig} and elem N"])  # H-bond donor/acceptor
    send_request("color", args=["salmon", f"{lig} and elem O"])  # H-bond acceptor
    send_request("color", args=["gold", f"{lig} and elem S"])  # sulfur
    send_request(
        "color", args=["palegreen", f"{lig} and (elem F or elem Cl or elem Br or elem I)"]
    )  # halogen

    # H-bonds to protein
    send_request("do", args=[f"delete {obj_name}_pharm_hbonds"])
    send_request(
        "do",
        args=[
            f"distance {obj_name}_pharm_hbonds, ({lig}) and (elem N or elem O or elem F), "
            f"({pocket_sel}) and (elem N or elem O), 3.5"
        ],
    )
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
def mutation_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
    mutations: Annotated[
        str, Field(description='Comma-separated mutation list (e.g. "A123G,V45L,T200S")')
    ],
) -> str:
    """
    Highlights mutated residues on the protein structure.

    Given a comma-separated list of mutations in standard notation (e.g.
    "A123G,V45L,T200S"), the mutated residues are shown as magenta sticks
    and labeled. Nearby residues (within 4 Å) are shown as thin grey sticks
    for packing context. The protein backbone is shown as a grey cartoon.
    Organic ligands are shown as yellow sticks.

    Mutation format: <wildtype_aa><resi><mutant_aa>, e.g. "A123G" (Ala→Gly
    at position 123). Chain can optionally be prefixed: "A:A123G".
    """
    import re

    mut_list = [m.strip() for m in mutations.split(",")]
    resi_list = []
    parsed = []
    for m in mut_list:
        match = re.search(r"(\d+)", m)
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
def textbook_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Configures PyMOL for a crisp, cel-shaded illustrative look ("Textbook Illustration").

    This view transforms the structure into a bold, 2D illustrative style with sharp
    black outlines, ideal for presentations or textbook-style diagrams. It hides
    the interior complexities, showing a solid white cartoon and surface with heavy
    black edge contours. Ligands are styled similarly as opaque white sticks with outlines.
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
    send_request("set", args=["ray_trace_mode", "3"])  # 3 = comic-book style coloring
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
def cinematic_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Configures PyMOL for a depth-cued, cinematic look with dramatic lighting.

    This view emphasizes volume and scale using deep shadows, fog, and depth-cueing.
    The core of the structure emerges from a dark background, making massive
    complexes (like ribosomes or viral capsids) look dramatic and imposing.
    Protein uses standard coloring but with altered material properties.
    """
    # Restore basic representation if not present
    send_request("show", args=["cartoon", f"({obj_name}) and polymer.protein"])
    send_request("show", args=["surface", f"({obj_name}) and polymer.protein"])

    # Dramatic deep black background
    send_request("do", args=["bg_color black"])

    # Enable fog and depth cueing
    send_request("set", args=["depth_cue", "1"])
    send_request("set", args=["fog", "1"])
    send_request("set", args=["fog_start", "0.45"])  # Fog starts mid-structure
    send_request("set", args=["fog_color", "black"])

    # Cinematic lighting and shadows
    send_request("set", args=["light_count", "2"])
    send_request("set", args=["spec_reflect", "0.3"])  # Slightly glossy
    send_request("set", args=["ray_shadows", "1"])
    send_request("set", args=["ray_shadow_decay_factor", "0.1"])
    send_request("set", args=["ray_shadow_decay_range", "3"])

    # Enhance the surface material
    send_request("set", args=["transparency", "0.0", obj_name])

    return f"Cinematic view applied to {obj_name}. Fog and depth-cueing enabled. Use the 'ray' command to see the full dramatic shadow effect."


@mcp.tool()
def pointillist_view(
    obj_name: Annotated[str, Field(description='PyMOL object name (e.g. "1abc")')],
) -> str:
    """
    Renders the structure as an artistic, abstract pointillist/starfield cloud.

    The continuous surface is replaced by thousands of individual dots representing
    the solvent-accessible surface, resembling a galaxy or pointillist painting.
    The protein backbone is hidden to emphasize the scattered volume. Ligands
    are shown as bright yellow spheres (stars) embedded in the cloud.
    """
    send_request("hide", args=["everything", obj_name])
    send_request("do", args=["bg_color black"])

    # The "Starfield" point cloud
    send_request("show", args=["dots", f"({obj_name}) and polymer.protein"])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])  # Default color recovery

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
