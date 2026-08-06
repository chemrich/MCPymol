"""Evolutionary conservation: MSA, Shannon entropy, and the view that uses it.

The pipeline is: pull the chain's sequence out of PyMOL, submit it to an
MMseqs2 server (the ColabFold public API by default), parse the A3M alignment,
score per-position entropy, then push the scores back into the B-factor column
for spectrum colouring.
"""

import hashlib
import json
import math
import os
import time
import urllib.error
import urllib.parse
import urllib.request

from mcpymol.app import mcp
from mcpymol.bridge import send_request

# MMseqs2 server for evolutionary conservation (ColabFold public API by default)
MMSEQS_URL = os.environ.get("MCPYMOL_MMSEQS_URL", "https://api.colabfold.com")

# Standard amino acid alphabet for Shannon entropy calculation
_AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

# In-memory cache: MD5(sequence) → list of per-residue entropy values
# Avoids repeat MMseqs2 API calls when only the scale is being changed
_conservation_cache: dict[str, list[float]] = {}


def _run_mmseqs2(sequence: str, server_url: str | None = None, use_env: bool = True) -> str:
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
    data = urllib.parse.urlencode(
        {
            "q": f">query\n{sequence}\n",
            "mode": mode,
        }
    ).encode()
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
                raise RuntimeError(
                    f"Failed to submit MMseqs2 job after {max_retries} attempts: {e}"
                ) from e
            time.sleep(2**attempt)

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
        import io
        import tarfile

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
    sequences: list[list[str]] = []
    current: list[str] = []
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


@mcp.tool()
def conservation_view(
    obj_name: str,
    selection: str = "all",
    server_url: str | None = None,
    use_env: bool = True,
    chain: str | None = None,
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

    # 5. Map conservation scores onto B-factor column.
    #
    # We invert entropy so that high B = conserved (intuitive spectrum coloring).
    # Both scale modes produce scores in [0, 100]: 100 = most conserved.
    #
    # IMPORTANT: structures often have non-contiguous residue numbering (gaps,
    # insertion codes, modified termini), so we can't just use i+1 → resi.  We
    # let PyMOL walk the CAs in selection order and map *its* resi → our score
    # via a stored dict, then apply the mapping in a single alter pass.  This
    # collapses ~2N socket round-trips into one.
    n_residues = len(entropies)
    min_entropy = min(entropies)
    max_entropy = max(entropies)
    entropy_range = max_entropy - min_entropy

    scores: list[float] = []
    for entropy in entropies:
        if scale == "relative" and entropy_range > 0:
            score = (1.0 - (entropy - min_entropy) / entropy_range) * 100.0
        else:
            score = (1.0 - entropy) * 100.0
        scores.append(round(score, 2))

    # One do block: zero out, push scores, build resi → score map by walking
    # CAs in selection order, then alter all atoms in one pass.
    apply_script = "\n".join(
        [
            f"alter {chain_sel}, b=0",
            f"stored.cons_scores = {json.dumps(scores)}",
            "stored.cons_map = {}",
            "stored.cons_idx = 0",
            (
                f"iterate {chain_sel} and name CA, "
                "stored.cons_map[resi] = stored.cons_scores[stored.cons_idx] "
                "if stored.cons_idx < len(stored.cons_scores) else 0.0; "
                "stored.cons_idx = stored.cons_idx + 1"
            ),
            f"alter {chain_sel}, b=stored.cons_map.get(resi, 0.0)",
            f"rebuild {obj_name}",
        ]
    )
    send_request("do", args=[apply_script])

    # 6. Apply visualization
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", obj_name])

    # Spectrum: magenta (conserved, high b) → white → cyan (variable, low b)
    send_request(
        "do", args=[f"spectrum b, cyan_white_magenta, {chain_sel}, minimum=0, maximum=100"]
    )

    # Show other chains as gray ghost for context
    other_chains_sel = f"({obj_name}) and polymer.protein and not chain {chain}"
    send_request("color", args=["gray50", other_chains_sel])
    send_request("set", args=["cartoon_transparency", "0.5", other_chains_sel])

    send_request("do", args=["bg_color black"])
    send_request("center", args=[obj_name])
    send_request("do", args=[f"origin {obj_name}"])

    scale_label = "relative" if scale == "relative" else "absolute"
    return (
        f"Showing conservation view for chain {chain} of {obj_name} ({scale_label} scale). "
        f"Magenta = conserved, white = moderate, cyan = variable. "
        f"MSA: {msa_note}, {n_residues} residue positions scored. "
        f"Entropy range: {min_entropy:.3f} – {max_entropy:.3f}."
    )
