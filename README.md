# MCPymol — talk to PyMOL

![Nucleosome core particle (1AOI) rendered in MCPymol's ghost-heart style](assets/nucleosome.png)

**MCPymol** is a [Model Context Protocol](https://modelcontextprotocol.io/) server that lets you drive PyMOL with natural language. Load structures, set up analytical views, measure things, and explore proteins by talking to Claude or Gemini. The image above was made by typing *"show me a nucleosome"* into Claude Code. That was the whole prompt.

PyMOL is great, but its syntax is famously obscure — and despite the name, it isn't quite Python. MCPymol is for people who'd rather just look at structures.

## What you get

- **A vocabulary the LLM understands.** Tools like `fetch_structure`, `ligand_view`, `interface_view`, `mutation_view`, `conservation_view` do high-level setup in one call: pick the biological assembly, hide solvent, color sensibly, label the right residues, draw the right H-bonds.
- **~60 PyMOL primitives** exposed as individual tools (`show`, `hide`, `color`, `select`, `distance`, `align`, `spectrum`, …) so the model can compose finer motions when the high-level tools don't quite fit.
- **Scene introspection.** `list_objects`, `list_chains`, `list_ligands` let the model check what's actually loaded before guessing.
- **Smart structure prep.** Fetching a PDB code grabs the biological assembly when one exists, then runs a BFS heuristic over chain–chain contacts (default radius 5 Å) so sprawling functional multimers like the CRP pentamer or ferritin cage stay whole while crystallographic copies get dropped. Waters and crystallization additives are hidden automatically.
- **Two-process bridge.** PyMOL's GUI has its own Python; MCPymol works by running a tiny TCP listener *inside* PyMOL plus a separate FastMCP server *outside* it.

## How it talks

```
┌──────────────┐   MCP / stdio    ┌────────────────┐   JSON over    ┌──────────────┐
│ Claude /     │ ───────────────▶ │ mcpymol server │ ◀── TCP :9876 ─▶ │ PyMOL GUI    │
│ Gemini CLI   │                  │  (FastMCP)     │                │  (plugin.py) │
└──────────────┘                  └────────────────┘                └──────────────┘
```

The plugin half runs inside PyMOL and dispatches to `pymol.cmd`. The bridge half is what the MCP client launches.

## Try it

> "Fetch ubiquitin (1ubq) and show it as a cartoon."
> "Color the alpha helices red and the beta sheets blue."
> "What ligands are in 1HSG?"
> "Show the binding pocket around MK1 in 1HSG."
> "Highlight mutations E6V, K16E, and V67F in hemoglobin (4HHB)."
> "Run a Poisson-Boltzmann electrostatics calculation on 1LYZ."

Tip: if the model isn't sure what's loaded, ask it to *list the objects* — it'll call `list_objects` and ground itself before guessing names.

## Installation

There are two halves to wire up: the **native plugin** (runs inside PyMOL) and the **MCP bridge** (runs outside, and is what your AI assistant launches).

### 1. Start the native plugin

Open PyMOL, then in the PyMOL command line:

```pymol
run /path/to/MCPymol/src/mcpymol/plugin.py
```

You should see `MCPymol Native Plugin listening on 127.0.0.1:9876`.

To auto-load it on every PyMOL launch, add this to `~/.pymolrc.py`:

```python
from pymol import cmd
cmd.do("run /absolute/path/to/MCPymol/src/mcpymol/plugin.py")
```

**Changing the port.** Set `MCPYMOL_PORT` before launching **both** PyMOL and the bridge:

```bash
MCPYMOL_PORT=9867 open -a PyMOL       # macOS
MCPYMOL_PORT=9867 uv run mcpymol      # bridge
```

### 2. Register the bridge with your AI assistant

Pick one. All paths are to the cloned repo.

```bash
git clone https://github.com/chemrich/MCPymol.git
cd MCPymol
```

#### Claude Code CLI (macOS, `uv`)

```bash
uv sync
claude mcp add mcpymol -- uv --directory /absolute/path/to/MCPymol run mcpymol
```

Start a new Claude Code session.

#### Claude Desktop (macOS, `uv`)

```bash
uv sync
```

Add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "mcpymol": {
      "command": "uv",
      "args": ["--directory", "/absolute/path/to/MCPymol", "run", "mcpymol"]
    }
  }
}
```

Restart Claude Desktop.

#### Gemini CLI (macOS, `uv`)

```bash
uv sync
gemini mcp add mcpymol uv --directory /absolute/path/to/MCPymol run mcpymol
gemini mcp refresh
```

#### Restricted environments — no `uv`

If `uv` is blocked by your org's security policy, use a standard venv. On Linux you may need `sudo apt-get install python3-venv pymol` first.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .

# Point the assistant directly at the venv binary
claude mcp add mcpymol /absolute/path/to/MCPymol/.venv/bin/mcpymol
# or
gemini mcp add mcpymol /absolute/path/to/MCPymol/.venv/bin/mcpymol
```

If your network blocks PyPI, you may also need to tweak the repository URLs inside `uv.lock` to match your internal mirror.

## Troubleshooting

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| Every tool returns *"Socket connection failed. Is the PyMOL plugin running?"* | Plugin not loaded in PyMOL | `run /path/to/plugin.py` inside PyMOL, or add it to `~/.pymolrc.py` |
| *"Address already in use"* on plugin start | A previous PyMOL session left the port open, or another app uses 9876 | Quit lingering PyMOL processes, or set `MCPYMOL_PORT=9867` on both sides |
| Long `get_fastastr` / `get_chains` calls fail with a JSON parse error | You're on a pre-2026-05 version of MCPymol that capped recv() at 8 KB | Pull main — the bridge now drains the response in full |
| `conservation_view` is slow | First call hits the ColabFold MMseqs2 API (30 s–few min); subsequent calls for the same sequence hit a local cache | If you have an internal MMseqs2 server, set `MCPYMOL_MMSEQS_URL` |
| `poisson_boltzmann_view` fails | `apbs` or `pdb2pqr` missing | `brew install brewsci/bio/apbs` and `pip install pdb2pqr` |

## Tests

```bash
PYTHONPATH=src uv run pytest tests/
```

The suite mocks the socket layer and PyMOL's `cmd` module, so it runs without a PyMOL install. 130+ tests cover the socket payloads, the bridge framing, the conservation pipeline (A3M parsing, Shannon entropy, MMseqs2 mocking, MSA→B-factor mapping), and every view.

## The view library

These are the high-level visualization tools. Each does its own setup — coloring, transparency, H-bonds, labels — in one prompt.

### `ligand_view` — binding site

Pocket residues (within 5 Å of the ligand) as element-colored sticks with CA labels, ligand as yellow sticks, H-bonds as yellow dashes, protein as a translucent cartoon.

> Show me the ATP binding site in 1ATP

![Ligand view of ATP in cAMP-dependent kinase (1ATP)](assets/ligand_view.png)

### `interface_view` — protein–protein interface

Chain A marine, chain B salmon. Interface residues (within 4 Å of the partner) as solid surface patches with sidechain sticks and CA labels. Cross-chain H-bonds as yellow dashes.

> Show the interface between chain A and chain D in 1BRS

![Interface view of barnase–barstar complex (1BRS)](assets/interface_view.png)

### `putty_view` — B-factor flexibility

Tube radius *and* color scale with B-factor: blue/thin = rigid, red/thick = flexible. A 70%-transparent surface adds shape context.

> Show the B-factor flexibility of 1UBQ as a putty view

![Putty view of ubiquitin (1UBQ)](assets/putty_view.png)

### `bfactor_view` — B-factor flexibility, plain cartoon

The same color story as `putty_view` (blue → red on B-factor) but on a plain cartoon. Cheaper, cleaner in figures where you don't want the putty distortion.

> Color 1UBQ by B-factor

### `hydrophobic_surface_view` — surface chemistry

Surface colored by amino-acid chemistry: orange = hydrophobic, white = polar, sky blue = positive, salmon = negative. Useful for spotting hydrophobic patches, membrane belts, and charge complementarity.

> Show the hydrophobic surface of 1TCA

![Hydrophobic surface view of *Candida antarctica* Lipase B (1TCA)](assets/hydrophobic_surface_view.png)

### `electrostatic_view` — approximate electrostatics

Red→white→blue surface coloring driven by per-residue pKa-weighted partial charges. Two modes: `atomic` (charges on the actual charge-center atoms — localized, natural falloff) and `residue` (uniform across charged residues — saturated patches).

> Show the electrostatic surface of 1LYZ

![Electrostatic view of lysozyme (1LYZ)](assets/electrostatic_view.png)

### `poisson_boltzmann_view` — true PB electrostatics

Full Poisson-Boltzmann potential via [APBS](https://github.com/Electrostatics/apbs) and [PDB2PQR](https://github.com/Electrostatics/pdb2pqr), mapped onto the surface at ±20 kT/e. Physically correct, accounts for solvent screening and ionic strength.

> **Prerequisites** (must be on `PATH`):
> ```bash
> brew install brewsci/bio/apbs
> pip install pdb2pqr
> ```

> Run a Poisson-Boltzmann electrostatics calculation on 1LYZ

![Poisson-Boltzmann electrostatic surface of lysozyme (1LYZ)](assets/poisson_boltzmann_view.png)

### `conservation_view` — evolutionary conservation

Pipeline: extract the chain's sequence → submit to MMseqs2 (ColabFold public API by default; override with `MCPYMOL_MMSEQS_URL`) → parse the A3M alignment → compute per-position Shannon entropy → map onto B-factor and color via `cyan_white_magenta` spectrum.

Magenta = conserved, white = moderate, cyan = variable. First call takes 30 s – few minutes depending on sequence length; results are cached in memory by sequence, so re-running on the same protein (or changing only the color scale) is instant.

> Color 1ubq by conservation

### `crosslink_view` — disulfides & metal coordination

CYS sidechains and disulfide bonds in yellow, metal coordination bonds in orange, the rest of the protein as a thin grey cartoon.

> Show the disulfide bonds in 1CEL

![Crosslink view of cellulase (1CEL)](assets/crosslink_view.png)

### `pocket_view` — binding pocket surface

The pocket cavity (residues within 5 Å of the ligand) as a semi-transparent surface colored by chemistry. Sticks for the pocket sidechains, yellow sticks for the ligand, cyan dashes for H-bonds.

> Show the binding pocket around MK1 in 1HSG

![Pocket view of MK1 binding site in HIV-1 protease (1HSG)](assets/pocket_view.png)

### `pharmacophore_view` — ligand pharmacophore features

The ligand colored by pharmacophore type: violet = aromatic ring carbon, yellow = aliphatic carbon, sky blue = N (donor/acceptor), salmon = O (acceptor), gold = S, pale green = halogen. Interacting residue sticks with CA labels, cyan dashes for H-bonds.

> Show the pharmacophore features of MK1 in 1HSG

![Pharmacophore view of MK1 in HIV-1 protease (1HSG)](assets/pharmacophore_view.png)

### `mutation_view` — mutation hotspots

Grey cartoon, mutated sidechains as magenta sticks with white CA labels, neighboring residues (within 4 Å) as thin element-colored sticks for packing context. Standard `A123G` notation; optional chain prefix `A:A123G`.

> Highlight mutations E6V, K16E, and V67F in hemoglobin (4HHB)

![Mutation view of hemoglobin (4HHB) showing E6V, K16E, V67F](assets/mutation_view.png)

### `textbook_view` — cel-shaded illustration

White cartoon + surface with heavy ray-trace contours. The cel-shaded look kicks in after you run `ray` (or ask the model to render).

> Make 4HHB look like a textbook illustration

### `cinematic_view` — fog and shadows

Depth-cueing + fog + soft shadows on a black background. Best on big assemblies — ribosomes, capsids, nucleosomes — where you want a sense of scale. Run `ray` for the full effect.

> Give me a cinematic view of the ribosome

### `pointillist_view` — starfield surface

Replaces the surface with a dense dot cloud; ligands become bright yellow stars. More art than analysis.

## The name

My best friend in high school once shared an apartment with MC Chris, who voiced MC Pee Pants in Aqua Teen Hunger Force. I'm not saying that was the inspiration for the name of this project, but I'm not denying it either.

## Provenance

Built on macOS using the open-source PyMOL available via Homebrew. Started with Antigravity, then Gemini Pro 3.1 until I ran out of tokens, then Claude Code (Sonnet 4.6 thinking). Tested with Claude Code and Gemini CLI. Conservation analysis uses the [ColabFold](https://github.com/sokrypton/ColabFold) public MMseqs2 API (please don't hammer it).

## License

[MIT](LICENSE).
