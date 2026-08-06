"""3D-print export: watertight, multi-colour STL out of a PyMOL scene.

PyMOL's open-source build cannot write STL, and its OBJ exporter dumps the
whole visible scene as non-manifold surface soup.  This module isolates each
colour group, exports it, and rebuilds a single watertight manifold solid per
group in the shared coordinate frame.
"""

import os
from typing import cast

from mcpymol.app import mcp
from mcpymol.bridge import send_request

# ── 3D printing export ───────────────────────────────────────────────────────
#
# PyMOL's open-source build cannot write STL, and its OBJ exporter writes the
# whole *visible* scene (it ignores the selection argument). The mesh it does
# produce is a non-manifold surface soup — unprintable as-is. This tool works
# around all three: it isolates each colour group's surface, exports OBJ, then
# rebuilds a watertight manifold per group while keeping every group in the
# shared PyMOL coordinate frame so they assemble correctly as multi-material
# parts in a slicer.

_PRINT_DEPS_HINT = (
    "3D-print export needs extra libraries. Install them with:\n"
    "  uv sync --extra print     (from a MCPymol checkout)\n"
    "  uv pip install 'mcpymol[print]'\n"
    "Provides: trimesh, pymeshlab, scipy, scikit-image."
)


def _parse_groups(groups: str) -> "list[tuple[str, str]]":
    """Parse ``"label=selection; label2=selection2"`` into ordered pairs."""
    pairs = []
    for chunk in groups.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"bad group spec {chunk!r}; expected 'label=selection'")
        label, sel = chunk.split("=", 1)
        label, sel = label.strip(), sel.strip()
        if not label or not sel:
            raise ValueError(f"bad group spec {chunk!r}")
        pairs.append((label, sel))
    if not pairs:
        raise ValueError("no groups given")
    return pairs


def _repair_to_stl(
    src_obj: str, dst_stl: str, method: str, voxel_pitch: float, poisson_depth: int
) -> dict:
    """Rebuild a PyMOL surface OBJ into a watertight, manifold STL.

    ``method``: ``poisson`` keeps surface detail (good for bulky chains);
    ``voxel`` is robust for thin tubular geometry (nucleic acids) and slightly
    thickens fragile features; ``auto`` uses the cheapest path that works —
    compact structures (e.g. a GFP barrel) often export already-watertight, so
    a light cleanup beats Poisson, which can degrade an already-closed surface.
    Coordinates are preserved so every group stays in one frame.
    """
    import trimesh

    def _light(mesh, dst):
        # Already-watertight export: drop tiny internal shells (buried cavities
        # print as useless floating geometry), keep the largest body, tidy up.
        bodies = sorted(mesh.split(only_watertight=False), key=lambda b: len(b.faces), reverse=True)
        shell = bodies[0] if bodies else mesh
        shell.merge_vertices()
        shell.update_faces(shell.unique_faces())
        shell.update_faces(shell.nondegenerate_faces())
        shell.remove_unreferenced_vertices()
        trimesh.repair.fix_normals(shell)
        trimesh.repair.fill_holes(shell)
        shell.export(dst)
        return shell

    def _voxel(src, dst, pitch):
        # force="mesh" collapses a scene into a single Trimesh, but the
        # declared return type is the Geometry base class, so narrow it.
        m = cast("trimesh.Trimesh", trimesh.load(src, force="mesh"))
        vox = m.voxelized(pitch=pitch).fill()
        out = vox.marching_cubes
        # marching_cubes is in voxel-index space; map back to world coords.
        out.apply_transform(vox.transform)
        # Fragmented/open input (e.g. a PyMOL cartoon triangle-soup of
        # separate arrow and tube segments) marching-cubes into many loose,
        # non-watertight shells. Consolidate to the single largest solid and
        # close it — same philosophy as _light — so the result is one
        # watertight, printable body.
        bodies = sorted(
            out.split(only_watertight=False),
            key=lambda b: len(b.faces),
            reverse=True,
        )
        out = bodies[0] if bodies else out
        out.merge_vertices()
        out.update_faces(out.unique_faces())
        out.update_faces(out.nondegenerate_faces())
        out.remove_unreferenced_vertices()
        trimesh.repair.fix_normals(out)
        trimesh.repair.fill_holes(out)
        out.export(dst)
        return out

    def _poisson(src, dst, depth):
        import pymeshlab

        def apply_first(ms, names, **kw):
            last = None
            for n in names:
                try:
                    ms.apply_filter(n, **kw)
                    return
                except Exception as e:
                    last = e
            raise RuntimeError(f"none of {names} worked: {last}")

        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(src)
        apply_first(ms, ["meshing_remove_duplicate_vertices", "remove_duplicate_vertices"])
        apply_first(ms, ["meshing_remove_null_faces", "remove_zero_area_faces"])
        apply_first(ms, ["meshing_remove_unreferenced_vertices", "remove_unreferenced_vertices"])
        apply_first(ms, ["compute_normal_per_vertex", "re_orient_all_faces_coherently"])
        apply_first(
            ms,
            [
                "generate_surface_reconstruction_screened_poisson",
                "surface_reconstruction_screened_poisson",
            ],
            depth=depth,
            preclean=True,
        )
        ms.save_current_mesh(dst, binary=True)
        return trimesh.load(dst, force="mesh")

    used = method
    if method == "auto":
        raw = cast("trimesh.Trimesh", trimesh.load(src_obj, force="mesh", process=True))
        if raw.is_watertight:
            mesh = _light(raw, dst_stl)
            used = "light (already watertight)"
        else:
            try:
                mesh = _poisson(src_obj, dst_stl, poisson_depth)
                used = "poisson"
            except Exception:
                mesh = _voxel(src_obj, dst_stl, voxel_pitch)
                used = "voxel (poisson fallback)"
    elif method == "poisson":
        mesh = _poisson(src_obj, dst_stl, poisson_depth)
        used = "poisson"
    elif method == "voxel":
        mesh = _voxel(src_obj, dst_stl, voxel_pitch)
    else:
        raise ValueError(f"unknown method {method!r}")

    return {
        "method": used,
        "faces": len(mesh.faces),
        "watertight": bool(mesh.is_watertight),
    }


@mcp.tool()
def print_ribbon_view(obj_name: str, spine_radius: float = 0.9) -> str:
    """
    Chunky β-arrow ribbons plus a continuous backbone "spine", tuned for
    rigid, gap-free 3D printing.

    Configures the look developed for FDM printing: thick β-strand arrows
    and a fat helix on the main object with loop cartoon hidden, plus a
    separate ``<obj>_spine`` object showing PyMOL's ``cartoon tube`` (which
    ignores secondary structure) running unbroken through the whole
    backbone. The spine threads through the strand bodies, so when the two
    objects are exported together the voxel step fuses them into ONE
    watertight solid with no strand→loop discontinuity. The spine also acts
    as internal rebar, reinforcing the thin junctions for print rigidity.

    After calling this, export the fused solid with::

        print_export(obj_name="<obj>",
                      groups="<obj>=(<obj> or <obj>_spine)",
                      representation="cartoon",
                      method="voxel", voxel_pitch=0.2)

    Args:
        obj_name: PyMOL object name (e.g. "1abc").
        spine_radius: Radius (A) of the continuous backbone tube. Larger
                      values give more internal reinforcement. Default 0.9.
    """
    spine = f"{obj_name}_spine"

    # Main object: chunky β-arrows + fat helix.
    send_request("hide", args=["everything", obj_name])
    send_request("show", args=["cartoon", obj_name])
    send_request("do", args=[f"cartoon automatic, {obj_name}"])
    chunky = {
        "cartoon_loop_radius": "1.3",
        "cartoon_rect_width": "2.0",
        "cartoon_rect_length": "2.8",
        "cartoon_oval_width": "1.3",
        "cartoon_oval_length": "2.0",
        "cartoon_helix_radius": "1.3",
        "cartoon_fancy_sheets": "1",  # keep 3D arrowheads (strand direction)
        "cartoon_flat_sheets": "1",  # well-defined arrow plane
        "cartoon_smooth_loops": "1",
    }
    for setting, value in chunky.items():
        send_request("set", args=[setting, value])

    # Loops are represented by the spine only - drop their fat cartoon.
    send_request("hide", args=["cartoon", f"({obj_name}) and not (ss H or ss S)"])

    # Spine: one continuous tube through the entire backbone (no SS gaps).
    send_request("delete", args=[spine])
    send_request("create", args=[spine, obj_name])
    send_request("hide", args=["everything", spine])
    send_request("show", args=["cartoon", spine])
    send_request("do", args=[f"cartoon tube, {spine}"])
    send_request("set", args=["cartoon_tube_radius", str(spine_radius), spine])

    send_request("do", args=[f"util.chainbow('{obj_name}')"])
    send_request("do", args=[f"util.chainbow('{spine}')"])
    send_request("orient", args=[obj_name])

    return (
        f"Print-ribbon view applied to {obj_name} (spine object '{spine}' "
        f"created). Chunky β-arrows and helix with a continuous backbone "
        f"tube; loops are the spine only. Export as one fused, gap-free "
        f"solid with:\n"
        f'  print_export(obj_name="{obj_name}", '
        f'groups="{obj_name}=({obj_name} or {spine})", '
        f'representation="cartoon", method="voxel", voxel_pitch=0.2)'
    )


@mcp.tool()
def print_export(
    obj_name: str,
    groups: str,
    out_dir: str = ".",
    method: str = "auto",
    voxel_pitch: float = 0.7,
    poisson_depth: int = 10,
    representation: str = "surface",
) -> str:
    """
    Exports a structure as watertight STL files ready for multi-colour 3D printing.

    Each colour group becomes one STL file. PyMOL's OBJ exporter writes the
    whole visible scene, so each group is isolated on its own before export,
    then rebuilt into a single watertight, manifold solid. All groups share the
    same coordinate frame, so a slicer can load them as aligned multi-material
    parts (e.g. add the second STL as a "part" of the first in Bambu Studio).

    Requires the optional ``print`` extra (trimesh, pymeshlab); see the install
    hint returned if the libraries are missing.

    Args:
        obj_name: PyMOL object to export (e.g. "1abc").
        groups: Semicolon-separated ``label=selection`` pairs, one per colour.
                Example: "protein=polymer.protein; nucleic=polymer.nucleic".
        out_dir: Directory for the STL files (default: current directory).
        method: "auto" (light cleanup if the export is already watertight,
                else poisson with voxel fallback), "poisson" (keeps detail,
                best for bulky chains), or "voxel" (robust for thin nucleic
                acids).
        voxel_pitch: Voxel size in Angstrom for the voxel method. Smaller keeps
                     more detail; 0.7 keeps a ~10 A helix intact.
        poisson_depth: Screened-Poisson octree depth; higher = more detail.
        representation: What geometry to export. "surface" (default,
                        unchanged) isolates each group as its own temp
                        object and exports its molecular surface. "cartoon"
                        exports the *currently displayed* cartoon geometry of
                        the real objects, preserving per-residue rep flags
                        and per-object cartoon type (e.g. a `cartoon tube`
                        spine from :func:`print_ribbon_view`). In cartoon
                        mode each group must name whole object(s) — one
                        colour per object (e.g. "1abc or 1abc_spine") — and
                        groups are isolated by toggling object visibility, no
                        temp objects.
    """
    try:
        import trimesh  # noqa: F401
    except ImportError:
        return _PRINT_DEPS_HINT

    try:
        pairs = _parse_groups(groups)
    except ValueError as e:
        return f"Error: {e}"

    representation = representation.lower()
    if representation not in ("surface", "cartoon"):
        return f"Error: representation must be 'surface' or 'cartoon', got {representation!r}"

    os.makedirs(out_dir, exist_ok=True)
    tmp_objs = []
    results = []
    try:
        if representation == "surface":
            # Build an isolated temp object per group up front.
            for label, sel in pairs:
                tmp = f"_print_{label}"
                send_request("delete", args=[tmp])
                send_request("create", args=[tmp, f"({obj_name}) and ({sel})"])
                tmp_objs.append((label, sel, tmp))

            for label, sel, tmp in tmp_objs:
                # Isolate: only this temp object visible, only as a surface.
                send_request("do", args=["disable all"])
                send_request("enable", args=[tmp])
                send_request("hide", args=["everything", tmp])
                send_request("show", args=["surface", tmp])

                obj_path = os.path.join(out_dir, f"{obj_name}_{label}.obj")
                stl_path = os.path.join(out_dir, f"{obj_name}_{label}.stl")
                res = send_request("save", args=[obj_path, tmp], timeout=180.0)
                if res.get("status") == "error":
                    return f"Error exporting group '{label}': {res.get('error')}"

                try:
                    info = _repair_to_stl(obj_path, stl_path, method, voxel_pitch, poisson_depth)
                except ImportError:
                    return _PRINT_DEPS_HINT
                results.append((label, sel, stl_path, info))
        else:
            # Cartoon mode: export the displayed cartoon of the real objects,
            # so per-residue rep flags (hidden loops) and per-object cartoon
            # type (a `cartoon tube` spine) are preserved exactly. PyMOL's
            # OBJ exporter dumps the whole visible scene, so isolate a group
            # by enabling only its objects. One colour per object.
            for label, sel in pairs:
                res = send_request("get_object_list", args=[sel])
                if res.get("status") == "error":
                    return f"Error resolving group '{label}': {res.get('error')}"
                objs = res.get("result") or []
                if not objs:
                    return (
                        f"Error: cartoon-mode group '{label}' selection "
                        f"'{sel}' matched no objects. Each group must name "
                        f"whole object(s), e.g. '1abc or 1abc_spine'."
                    )

                send_request("do", args=["disable all"])
                for obj in objs:
                    send_request("enable", args=[obj])

                obj_path = os.path.join(out_dir, f"{obj_name}_{label}.obj")
                stl_path = os.path.join(out_dir, f"{obj_name}_{label}.stl")
                res = send_request("save", args=[obj_path, sel], timeout=180.0)
                if res.get("status") == "error":
                    return f"Error exporting group '{label}': {res.get('error')}"

                try:
                    info = _repair_to_stl(obj_path, stl_path, method, voxel_pitch, poisson_depth)
                except ImportError:
                    return _PRINT_DEPS_HINT
                results.append((label, sel, stl_path, info))
    finally:
        for _, _, tmp in tmp_objs:
            send_request("delete", args=[tmp])
        send_request("do", args=["enable all"])

    lines = [f"Exported {len(results)} group(s) for 3D printing to {out_dir}:"]
    for label, sel, stl_path, info in results:
        flag = "OK" if info["watertight"] else "NOT watertight - check mesh"
        lines.append(
            f"  {label} ({sel}) -> {os.path.basename(stl_path)}  "
            f"[{info['method']}, {info['faces']:,} faces, {flag}]"
        )
    lines.append(
        "All groups share one coordinate frame. In your slicer, load the first "
        "STL then add the others as parts (don't re-centre) and assign a "
        "filament per part."
    )
    return "\n".join(lines)
