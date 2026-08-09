"""Reading a heterogeneity job directory: ordered volumes, latents, and a method.

Every method in this space — cryoDRGN, 3DVA, RECOVAR, 3DFlex, DynaMight — puts
essentially the same thing on disk: **a directory of maps plus a table of latent
coordinates.** RECOVAR's own documentation says so in those words. So one reader
covers them, with per-method detection in front of it.

**Why the method name is load-bearing rather than cosmetic.** SPEC invariant
**I2** forbids an unlabelled latent rendering: every latent view has to name the
method that generated it and carry that method's interpretive caveat, because
the caveats are not interchangeable. cryoDRGN's latent *density* can bear no
relation to the truth; 3DVA's frames are linear interpolations that no particle
ever occupied. A viewer that draws both the same way is asserting something for
one of them that only holds for the other. So :func:`load_ensemble` identifies
the method, and :mod:`mcpymol.wiggles.latent` **refuses to render** when it
could not.

**Detection is evidence-based and never a guess.** Each method is recognised by
markers its own documentation describes — cryoDRGN's ``z.pkl`` beside a
``config.yaml``, RELION's STAR files, cryoSPARC's ``.cs`` files. Anything else
loads its volumes perfectly well and stays ``UNKNOWN``, and the report says what
was looked for and what was found. That is not a failure mode; it is I2 working.
The caller who knows what they ran can pass ``method=`` and say so.

**Provenance.** Volumes in a heterogeneity job are decoder output or subspace
interpolations — nothing observed exactly that density. :func:`load_ensemble`
therefore declares :data:`~mcpymol.wiggles.provenance.Provenance.GENERATED`
rather than leaving them ``UNKNOWN``, and the report states that it did so and
on what evidence. This is not the inference I1 forbids: I1 exists to stop a
generated volume being asserted as *measured*, and this asserts the opposite,
conservative direction. An explicit ``provenance=`` still wins.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from mcpymol.wiggles.mapinfo import MapHeader, read_map_header
from mcpymol.wiggles.port import PortError, PymolPort, call
from mcpymol.wiggles.provenance import Provenance, declare

#: Volume extensions a heterogeneity job writes, in the order we prefer them.
MAP_SUFFIXES = (".mrc", ".map", ".ccp4")

#: Ceiling on how many volumes one ensemble may load. A cryoDRGN analysis
#: directory can hold hundreds of sampled volumes, and loading all of them into
#: a live session is not a viewer operation. The report says what was dropped —
#: a silent truncation would read as "this is the whole trajectory".
DEFAULT_MAX_VOLUMES = 50


class Method(Enum):
    """The heterogeneity method that produced an ensemble.

    Named rather than ranked. The Flatiron blind challenge found no method
    consistently outperforms the others, so this package treats them as
    interchangeable producers of an ensemble and stays out of adjudicating
    them — see the compendium entry `benchmarks`.
    """

    CRYODRGN = "cryodrgn"
    THREE_DVA = "3dva"
    RECOVAR = "recovar"
    THREE_DFLEX = "3dflex"
    DYNAMIGHT = "dynamight"
    CRYOSPIRE = "cryospire"
    UNKNOWN = "unknown"

    @property
    def label(self) -> str:
        return {
            Method.CRYODRGN: "cryoDRGN",
            Method.THREE_DVA: "cryoSPARC 3D Variability Analysis",
            Method.RECOVAR: "RECOVAR",
            Method.THREE_DFLEX: "cryoSPARC 3DFlex",
            Method.DYNAMIGHT: "DynaMight (RELION-5)",
            Method.CRYOSPIRE: "CryoSPIRE",
            Method.UNKNOWN: "unidentified",
        }[self]

    @property
    def caveat(self) -> str:
        """The interpretive warning I2 requires this method's views to carry.

        Every one of these is the method's own documented limitation, taken
        from the compendium entry of the same name — not this package's
        opinion of it.
        """
        return {
            Method.CRYODRGN: (
                "DENSITY IN CRYODRGN'S LATENT SPACE IS NOT RELIABLY MEANINGFUL. "
                "RECOVAR's authors report cases where it bears no relation to the "
                "true underlying density. Cluster size here is not population, and "
                "a neural decoder can produce plausible density in poorly-sampled "
                "regions of latent space."
            ),
            Method.THREE_DVA: (
                "3DVA IS A LINEAR SUBSPACE MODEL. Real domain rotations are not "
                "linear, so a traversal frame is an interpolation along a straight "
                "line in volume space — not a conformation any particle occupied. "
                "Curved motion appears as apparent density loss at the extremes."
            ),
            Method.RECOVAR: (
                "RECOVAR IS A LINEAR-SUBSPACE METHOD and shares that ceiling: "
                "motions genuinely outside the subspace are not represented. Its "
                "own contribution is a better-founded density estimate, which is "
                "an argument about other methods' densities, not a guarantee of "
                "this one's."
            ),
            Method.THREE_DFLEX: (
                "3DFLEX FRAMES ARE MODEL OUTPUT. The deformation field is fitted, "
                "and the authors document limits on intricate motions. The motion "
                "is the well-supported claim here; the frames are its rendering, "
                "not observations."
            ),
            Method.DYNAMIGHT: (
                "DYNAMIGHT VOLUMES ARE DECODED FROM A DEFORMED GAUSSIAN MODEL. "
                "The authors are explicit about the method's limits, and estimated "
                "deformations should be read with the half-set uncertainty they "
                "provide wherever it is available."
            ),
            Method.CRYOSPIRE: (
                "CRYOSPIRE IS A 2025 PREPRINT and its on-disk conventions are not "
                "stable. Treat both this reader and the output as provisional."
            ),
            Method.UNKNOWN: (
                "METHOD NOT IDENTIFIED. Nothing in this directory matched a "
                "documented marker, so no caveat can be attached — and a latent "
                "view without its method's caveat is exactly what SPEC invariant "
                "I2 forbids."
            ),
        }[self]


#: Filename and directory markers that identify a method, each taken from what
#: the method's own documentation says it writes. Order matters only in that
#: the first method with a matching marker wins; the markers are disjoint in
#: practice. Detection deliberately stops here rather than guessing from
#: anything looser — see the module docstring.
_MARKERS: tuple[tuple[Method, tuple[str, ...], str], ...] = (
    (
        Method.CRYODRGN,
        ("z.pkl", "config.yaml", "weights.pkl"),
        "cryoDRGN writes z.pkl (per-particle latents) beside a config.yaml and model weights",
    ),
    (
        Method.RECOVAR,
        ("recovar", "eigenvolumes", "reordered_test_result.pkl"),
        "RECOVAR names its outputs after itself and emits eigenvolumes",
    ),
    (
        Method.DYNAMIGHT,
        ("deformations.star", "inverse_deformations", "run_data.star"),
        "DynaMight writes RELION-5 STAR conventions under the job directory",
    ),
    (
        Method.CRYOSPIRE,
        ("cryospire", "part_segmentation"),
        "CryoSPIRE emits a part segmentation over the density",
    ),
    (
        Method.THREE_DFLEX,
        ("flex_", "3dflex", "tetra_mesh"),
        "3DFlex carries a tetrahedral mesh and a flex-prefixed job name",
    ),
    (
        Method.THREE_DVA,
        ("3dva", "component_", "particles.cs"),
        "cryoSPARC 3DVA-Display writes per-component volume series beside .cs files",
    ),
)

#: Latent-table extensions this package reads without help.
_TEXT_LATENT = (".star", ".csv", ".tsv", ".txt")


@dataclass(frozen=True)
class LatentTable:
    """Per-particle latent coordinates, or a statement of why they are missing.

    Optional by design. The traversal is driven by the *volumes*, which are
    ordered on disk; the latent table informs the readout and nothing renders
    without it. An unreadable table is reported, never fatal.
    """

    source: Path | None = None
    rows: tuple[tuple[float, ...], ...] = ()
    unread_reason: str | None = None

    @property
    def dimensions(self) -> int:
        return len(self.rows[0]) if self.rows else 0

    @property
    def available(self) -> bool:
        return bool(self.rows)


@dataclass(frozen=True)
class Ensemble:
    """An ordered set of volumes from one heterogeneity job."""

    name: str
    directory: Path
    method: Method
    volumes: tuple[Path, ...]
    objects: tuple[str, ...]
    """Parallel to ``volumes``: frame *i* is ``volumes[i]`` on disk and
    ``objects[i]`` in the session."""

    headers: tuple[MapHeader, ...]
    latent: LatentTable
    evidence: tuple[str, ...] = ()
    dropped: int = 0
    method_declared: bool = False
    """True when the caller named the method rather than the directory revealing it."""

    @property
    def n_frames(self) -> int:
        return len(self.volumes)


# Process-global, like the other registries in this package. Correct for a
# single-session MCP server, which is what MCPymol is today; wrong the moment
# one process serves two sessions. Flagged rather than solved — see MOVING.md.
_ENSEMBLES: dict[str, Ensemble] = {}


def loaded_ensemble(name: str) -> Ensemble | None:
    """The ensemble registered under ``name``, or None."""
    return _ENSEMBLES.get(name)


def forget_ensemble(name: str | None = None) -> None:
    """Drop the record for ``name``, or all of them."""
    if name is None:
        _ENSEMBLES.clear()
    else:
        _ENSEMBLES.pop(name, None)


def _natural_key(path: Path) -> tuple:
    """Sort ``vol_2`` before ``vol_10``.

    Frame order *is* the trajectory. Lexicographic sorting would put frame 10
    between 1 and 2 and silently reorder the motion — the animation would still
    play, and it would be wrong.
    """
    return tuple(
        int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", path.name)
    )


def find_volumes(directory: Path) -> list[Path]:
    """Every volume in ``directory``, in trajectory order.

    Searches one level of subdirectories as well, because cryoDRGN puts sampled
    volumes under ``analyze.N/kmeans20/`` rather than at the top.
    """
    found: list[Path] = []
    for suffix in MAP_SUFFIXES:
        found.extend(directory.glob(f"*{suffix}"))
        found.extend(directory.glob(f"*/*{suffix}"))
        found.extend(directory.glob(f"*/*/*{suffix}"))
    # A half-map or a mask sitting in the job directory is not a frame.
    frames = [
        p for p in dict.fromkeys(found) if not re.search(r"(half[-_]?map|_msk|_mask)", p.name, re.I)
    ]
    return sorted(frames, key=_natural_key)


def detect_method(directory: Path) -> tuple[Method, list[str]]:
    """Identify the method from documented on-disk markers.

    Returns the method and the evidence for it. ``UNKNOWN`` with a list of what
    was searched is a legitimate, common answer — see the module docstring.
    """
    names = {p.name.lower() for p in directory.rglob("*") if p.is_file()}
    names |= {p.name.lower() for p in directory.rglob("*") if p.is_dir()}
    names.add(directory.name.lower())
    haystack = " ".join(sorted(names))

    for method, markers, why in _MARKERS:
        hit = next((m for m in markers if m in haystack), None)
        if hit is not None:
            return method, [f"found {hit!r} — {why}"]

    searched = sorted({m for _, markers, _ in _MARKERS for m in markers})
    return Method.UNKNOWN, [
        f"no documented marker found. Searched for: {', '.join(searched)}.",
        "The volumes still load; only latent views are withheld (SPEC I2). "
        "Pass method= to declare it if you know what produced this directory.",
    ]


def _read_text_latent(path: Path) -> LatentTable:
    """Parse a STAR loop or a delimited numeric table.

    STAR is handled because RELION-5 writes it; the rest is whitespace- or
    comma-delimited numbers, which is what everything else exports to.
    """
    try:
        text = path.read_text(errors="replace")
    except OSError as exc:
        return LatentTable(source=path, unread_reason=f"could not read: {exc}")

    rows: list[tuple[float, ...]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "data_", "loop_", "_", ";")):
            continue
        parts = re.split(r"[,\t ]+", stripped)
        try:
            values = tuple(float(p) for p in parts if p)
        except ValueError:
            continue  # a label column, or a header row
        if values:
            rows.append(values)

    if not rows:
        return LatentTable(source=path, unread_reason="no numeric rows found")

    width = min(len(r) for r in rows)
    return LatentTable(source=path, rows=tuple(r[:width] for r in rows))


def _read_npy_latent(path: Path) -> LatentTable:
    """Parse a plain ``.npy`` array without numpy and without unpickling.

    The format is documented and simple: a magic string, a header dict, then
    raw little-endian data. Object arrays are refused rather than unpickled —
    reading one executes whatever is inside it.
    """
    import ast
    import struct

    try:
        raw = path.read_bytes()
    except OSError as exc:
        return LatentTable(source=path, unread_reason=f"could not read: {exc}")

    if raw[:6] != b"\x93NUMPY":
        return LatentTable(source=path, unread_reason="not a .npy file")

    major = raw[6]
    header_len_size = 2 if major == 1 else 4
    start = 8 + header_len_size
    (header_len,) = struct.unpack("<H" if major == 1 else "<I", raw[8:start])
    try:
        header = ast.literal_eval(raw[start : start + header_len].decode("latin1").strip())
    except (ValueError, SyntaxError) as exc:
        return LatentTable(source=path, unread_reason=f"unreadable header: {exc}")

    descr, shape = header.get("descr"), header.get("shape")
    if not isinstance(descr, str) or descr.startswith(("|O", "O")):
        return LatentTable(
            source=path,
            unread_reason=(
                "array holds Python objects, which can only be read by unpickling, "
                "and unpickling runs whatever code the file contains"
            ),
        )
    fmt = {"<f4": "f", "<f8": "d", "<i4": "i", "<i8": "q"}.get(descr)
    if fmt is None or header.get("fortran_order"):
        return LatentTable(source=path, unread_reason=f"unsupported layout {descr!r}/fortran_order")
    if not (isinstance(shape, tuple) and len(shape) == 2):
        return LatentTable(source=path, unread_reason=f"expected a 2-D array, got shape {shape}")

    body = raw[start + header_len :]
    n_rows, n_cols = shape
    values = struct.unpack(
        f"<{n_rows * n_cols}{fmt}", body[: n_rows * n_cols * struct.calcsize(fmt)]
    )
    rows = tuple(
        tuple(float(v) for v in values[i * n_cols : (i + 1) * n_cols]) for i in range(n_rows)
    )
    return LatentTable(source=path, rows=rows)


def read_latent_table(directory: Path, *, trust_pickle: bool = False) -> LatentTable:
    """Find and read the latent coordinates, if they are readable at all.

    ``.pkl`` is **not read by default**, and that is a security position rather
    than a limitation: unpickling executes arbitrary code, and a job directory
    can come from a collaborator, a shared scratch volume or a download. Pass
    ``trust_pickle=True`` for output you produced yourself.
    """
    candidates = [
        p
        for p in sorted(directory.rglob("*"), key=_natural_key)
        if p.is_file() and re.search(r"(^|[_./])(z|latent|embedding|deformations)", p.name, re.I)
    ]
    if not candidates:
        return LatentTable(unread_reason="no latent table found in this directory")

    # Every candidate is tried before giving up, and a readable table always
    # wins over an unreadable one. Returning on the first failure would let a
    # `z.pkl` sitting beside a perfectly readable `z_values.txt` mask it —
    # the refusal is about the pickle, not about the directory.
    first_failure: LatentTable | None = None
    for path in candidates:
        suffix = path.suffix.lower()
        if suffix in _TEXT_LATENT:
            table = _read_text_latent(path)
        elif suffix in (".npy", ".cs"):
            table = _read_npy_latent(path)
        elif suffix == ".pkl":
            table = (
                _read_pickle_latent(path)
                if trust_pickle
                else LatentTable(
                    source=path,
                    unread_reason=(
                        "a .pkl latent table was found but not read. Unpickling runs "
                        "arbitrary code from the file, so it is opt-in: pass "
                        "trust_pickle=True if this is your own output. Nothing else "
                        "depends on it — the trajectory is ordered by the volumes."
                    ),
                )
            )
        else:
            continue
        if table.available:
            return table
        if first_failure is None:
            first_failure = table

    return first_failure or LatentTable(
        unread_reason="no readable latent table found in this directory"
    )


def _read_pickle_latent(path: Path) -> LatentTable:
    """Unpickle a latent table. Only ever called with the caller's consent."""
    import pickle

    try:
        with open(path, "rb") as fh:
            obj = pickle.load(fh)
    except Exception as exc:
        return LatentTable(source=path, unread_reason=f"could not unpickle: {exc}")

    try:
        rows = tuple(tuple(float(v) for v in row) for row in obj)
    except (TypeError, ValueError) as exc:
        return LatentTable(source=path, unread_reason=f"not a table of numbers: {exc}")
    return LatentTable(source=path, rows=rows)


def load_ensemble(
    port: PymolPort,
    directory: str | Path,
    name: str | None = None,
    *,
    method: Method | str | None = None,
    provenance: Provenance | str = Provenance.GENERATED,
    max_volumes: int = DEFAULT_MAX_VOLUMES,
    trust_pickle: bool = False,
) -> str:
    """Load a heterogeneity job's volumes in trajectory order.

    Args:
        port: A live or fake PyMOL port.
        directory: The job directory.
        name: Prefix for the loaded objects. Defaults to the directory name.
        method: Declare the method rather than letting detection decide.
        provenance: Defaults to ``generated`` — see the module docstring for why
            that is a declaration rather than the inference I1 forbids.
        max_volumes: Ceiling on frames loaded. Excess is reported, not hidden.
        trust_pickle: Permit reading a ``.pkl`` latent table.

    Returns:
        A report: the method and the evidence for it, the frame count, geometry,
        and what could not be read.

    Raises:
        PortError: the directory holds no volumes, or a load produced no object.
        ValueError: ``method`` or ``provenance`` is not a recognised value.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise PortError(f"{directory} is not a directory")

    volumes = find_volumes(directory)
    if not volumes:
        raise PortError(
            f"no MRC/CCP4 volumes under {directory}. A heterogeneity job should "
            f"contain a series of maps; looked for {', '.join(MAP_SUFFIXES)} at the "
            f"top level and two levels down, ignoring half-maps and masks."
        )

    dropped = max(0, len(volumes) - max_volumes)
    volumes = volumes[:max_volumes]

    if method is None:
        detected, evidence = detect_method(directory)
        declared = False
    else:
        detected = Method(method) if isinstance(method, str) else method
        evidence = [f"method declared by the caller as {detected.value!r}, not detected"]
        declared = True

    if isinstance(provenance, str):
        provenance = Provenance(provenance)

    prefix = name or re.sub(r"\W+", "_", directory.name) or "ensemble"
    width = max(2, len(str(len(volumes))))

    headers: list[MapHeader] = []
    objects: list[str] = []
    for index, path in enumerate(volumes, start=1):
        headers.append(read_map_header(path))  # before PyMOL, as load_map does
        obj = f"{prefix}_f{index:0{width}d}"
        call(port, "load", str(path), obj)
        objects.append(obj)

    names = port.query("get_names", "objects")
    missing = [o for o in objects if not (isinstance(names, (list, tuple)) and o in names)]
    if missing:
        for obj in objects:  # scope cleanup to what we made, never the session
            call(port, "delete", obj)
        raise PortError(
            f"{len(missing)} of {len(objects)} volumes did not arrive in the session "
            f"(first missing: {missing[0]!r}). Everything this call created was "
            f"removed; nothing else was touched."
        )

    for obj in objects:
        declare(obj, provenance)

    latent = read_latent_table(directory, trust_pickle=trust_pickle)
    ensemble = Ensemble(
        name=prefix,
        directory=directory,
        method=detected,
        volumes=tuple(volumes),
        objects=tuple(objects),
        headers=tuple(headers),
        latent=latent,
        evidence=tuple(evidence),
        dropped=dropped,
        method_declared=declared,
    )
    _ENSEMBLES[prefix] = ensemble

    return _report(ensemble, provenance)


def _report(ensemble: Ensemble, provenance: Provenance) -> str:
    first = ensemble.headers[0]
    vx, _, _ = first.voxel_size
    lines = [
        f"load_ensemble({ensemble.directory.name} -> {ensemble.name})",
        "",
        f"  Method: {ensemble.method.label}"
        + ("  [declared]" if ensemble.method_declared else "  [detected]"),
    ]
    lines += [f"    {reason}" for reason in ensemble.evidence]
    lines += [
        "",
        f"  {ensemble.n_frames} frames, in trajectory order, as "
        f"`{ensemble.objects[0]}` … `{ensemble.objects[-1]}`.",
        f"  Grid: {first.nx}x{first.ny}x{first.nz}, voxel "
        f"{'unknown' if vx is None else format(vx, '.4g')} Å.",
    ]

    if ensemble.dropped:
        lines += [
            f"  {ensemble.dropped} further volume(s) were NOT loaded — the ceiling is "
            f"{ensemble.n_frames}. Raise max_volumes to include them; what is shown "
            f"is a truncation, not the whole job.",
        ]

    lines.append("")
    if ensemble.latent.available:
        lines.append(
            f"  Latent table: {ensemble.latent.dimensions}-D, "
            f"{len(ensemble.latent.rows)} particles, from "
            f"{ensemble.latent.source.name if ensemble.latent.source else '?'}."
        )
    else:
        lines += [
            "  Latent table: not read.",
            f"    {ensemble.latent.unread_reason}",
            "    Nothing here depends on it — the trajectory is ordered by the volumes.",
        ]

    lines += [
        "",
        f"  Provenance: {provenance.value.upper()} (declared for all {ensemble.n_frames} frames)",
        "  " + provenance.caveat,
    ]
    if provenance is Provenance.GENERATED:
        lines += [
            "  Declared rather than inferred: these are a decoder's output or a",
            "  subspace interpolation, so asserting anything weaker would overstate",
            "  them. Pass provenance= to override.",
        ]

    if ensemble.method is Method.UNKNOWN:
        lines += [
            "",
            "  The method could not be identified, so latent views will REFUSE to",
            "  render (SPEC invariant I2 — no unlabelled latent plot). The volumes",
            "  themselves are loaded and usable. Pass method= to declare it.",
        ]

    warnings = first.warnings()
    if warnings:
        lines += ["", f"  Geometry warnings, first frame ({len(warnings)})"]
        lines += [f"    ! {w}" for w in warnings]

    return "\n".join(lines)
