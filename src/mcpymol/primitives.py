"""Thin wrappers over individual ``pymol.cmd`` functions.

One tool per PyMOL command, so the model can compose finer motions when the
high-level views in ``mcpymol.views`` do not quite fit.  Each wrapper is a real
``def`` — its signature and docstring are what the MCP client sees — delegating
to :func:`mcpymol.bridge._call` for the shared forwarding logic.
"""

from mcpymol.app import mcp
from mcpymol.bridge import _SLOW_OP_TIMEOUT, _call, format_measurement, send_request

# ── Primitive tools ──────────────────────────────────────────────────────────
#
# A note on PyMOL selection syntax (what to put in ``selection``):
#   "all"                      every atom in every loaded object
#   "1abc"                     all atoms of the object named 1abc
#   "1abc and chain A"         chain A of object 1abc
#   "1abc and resi 10-20"      residue range
#   "1abc and resn ATP"        residues with name ATP (e.g. ligands)
#   "1abc and name CA"         atom names
#   "1abc and ss H"            secondary-structure helix
#   "polymer.protein"          protein-only macro
#   "organic"                  small-molecule cofactors
#   "byres (a around 5)"       residues with any atom within 5 Å of selection a


@mcp.tool()
def show(representation: str, selection: str = "all") -> str:
    """Shows a graphical representation for a given selection.

    Args:
        representation: One of ``cartoon``, ``sticks``, ``spheres``, ``surface``,
            ``lines``, ``ribbon``, ``dots``, ``mesh``, ``nb_spheres``, ``labels``.
        selection: PyMOL selection string. See module-level note for syntax.
    """
    res = send_request("show", args=[representation, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Showing {representation} for selection '{selection}'"


@mcp.tool()
def hide(representation: str, selection: str = "all") -> str:
    """Hides a graphical representation for a given selection.

    Args:
        representation: Same vocabulary as :func:`show`, plus ``everything`` to
            hide all reps for the selection.
        selection: PyMOL selection string.
    """
    res = send_request("hide", args=[representation, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Hiding {representation} for selection '{selection}'"


@mcp.tool()
def color(color_name: str, selection: str = "all") -> str:
    """Sets the color for a selection.

    Args:
        color_name: A PyMOL color. Common names: ``red``, ``blue``, ``green``,
            ``yellow``, ``magenta``, ``cyan``, ``orange``, ``salmon``, ``marine``,
            ``forest``, ``palegreen``, ``skyblue``, ``violet``, ``grey50``,
            ``white``, ``black``. Use ``atomic`` to color non-carbon atoms by
            element while leaving carbons untouched.
        selection: PyMOL selection string.
    """
    res = send_request("color", args=[color_name, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Colored selection '{selection}' with {color_name}"


@mcp.tool()
def select(name: str, selection: str) -> str:
    """Creates (or replaces) a named selection for later reuse.

    Args:
        name: Identifier you'll refer to later (e.g. ``active_site``).
        selection: PyMOL selection expression to assign to that name.
    """
    res = send_request("select", args=[name, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Created named selection '{name}' for '{selection}'"


@mcp.tool()
def remove(selection: str) -> str:
    """
    Permanently removes the atoms matching the selection. This deletes atoms; it does not just
    hide them. To hide instead, use :func:`hide`.

    Args:
        selection: Atoms to delete permanently. To hide them instead, use hide.
    """
    res = send_request("remove", args=[selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Removed selection '{selection}'"


@mcp.tool()
def distance(name: str, selection1: str, selection2: str) -> str:
    """Measures the distance between two selections and returns it in Angstrom.

    Also draws the measurement in the viewport as a named distance object.
    With multi-atom selections PyMOL reports the average over the pairs it
    found within its default cutoff.

    Args:
        name: Name for the distance object (used to delete/recolor later).
        selection1: First selection.
        selection2: Second selection.
    """
    res = send_request("distance", args=[name, selection1, selection2])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return format_measurement(
        f"Distance between '{selection1}' and '{selection2}'",
        res.get("result"),
        "A",
        name,
    )


@mcp.tool()
def sasa(selection: str = "all", state: str | None = "1") -> str:
    """Measures solvent-accessible surface area, in square Angstrom.

    Reports the SASA of ``selection`` *in the context of the object it belongs
    to* — so a chain measured inside a complex is already partly occluded by
    its partner. To get the free (unbound) area, copy the chain to its own
    object first with ``create``, or use ``interface_report``, which does the
    bound/free bookkeeping for you.

    Accuracy depends on PyMOL's ``dot_solvent`` (0 = molecular surface,
    1 = solvent-accessible) and ``dot_density`` settings.

    Args:
        selection: What to measure (e.g. "1brs and chain A").
        state: Object state to measure. "1" is the first/only state.
    """
    return _call(
        "get_area",
        _measures=("Solvent-accessible surface area", "A^2"),
        selection=selection,
        state=state,
    )


@mcp.tool()
def rms_cur(mobile: str, target: str) -> str:
    """Measures RMSD between two selections *without* moving anything.

    Use this when the structures are already superposed, or when you want to
    know how far apart they are as currently positioned. Compare with ``align``
    and ``super``, which move the mobile structure to minimise RMSD before
    reporting it, and with ``superposition_view``, which shows where the
    difference is rather than summarising it as one number.

    Requires the two selections to have matching atom counts.

    Args:
        mobile: First selection.
        target: Second selection, same number of atoms.
    """
    return _call("rms_cur", _measures=("RMSD (as positioned)", "A"), mobile=mobile, target=target)


@mcp.tool()
def count_atoms(selection: str = "all") -> str:
    """Counts the atoms matching a selection.

    Handy for checking a selection expression does what you think before
    building a scene on top of it — an empty count means the expression is
    wrong, which is otherwise invisible until the picture comes out blank.

    Args:
        selection: PyMOL selection string to count.
    """
    res = send_request("count_atoms", args=[selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    try:
        n = int(res.get("result"))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return f"count_atoms returned an unexpected value: {res.get('result')!r}"
    if n == 0:
        return f"Selection '{selection}' matches no atoms. Check the object name and syntax."
    return f"Selection '{selection}' matches {n:,} atoms."


@mcp.tool()
def execute_pymol_command(command: str) -> str:
    """
    Executes a raw PyMOL command string (PyMOL CLI syntax). PREFER the dedicated tools when one
    exists — show, color, select, distance, ligand_view, interface_view, etc. They have better
    defaults, do compound setup in one call, and produce cleaner results. Reach for this tool
    only when no other tool covers what you need (e.g. ``set ray_shadow, 0``, ``bg_color
    grey20``, multi-statement scripts). Note: this accepts the PyMOL ``cmd.do`` mini-language,
    not Python.

    Args:
        command: A PyMOL CLI command, e.g. `set ray_shadow, 0`. Not Python.
    """
    res = send_request("do", args=[command])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Executed command: {command}"


@mcp.tool(name="as")
def as_tool(representation: str, selection: str | None = "all") -> str:
    """
    Shows one representation while hiding all others for the specified selection

    Args:
        representation: The one representation to leave visible: ``cartoon``, ``sticks``,
            ``spheres``, ``surface``, ``lines``, ``ribbon``, ``dots``, ``mesh``.
        selection: PyMOL selection string, e.g. "1abc and chain A". See the module note for
            syntax.
    """
    return _call("as", representation=representation, selection=selection)


@mcp.tool(name="set")
def set_setting(setting: str, value: str, selection: str | None = None) -> str:
    """
    Sets a PyMOL setting to a specified value

    Args:
        setting: PyMOL setting name, e.g. "cartoon_transparency", "ray_shadow", "sphere_scale".
            Hundreds exist; see the PyMOL settings reference.
        value: The value as a string, e.g. "0.5", "1", "black".
        selection: Limit the setting to this object or selection. Omit to set it globally.
    """
    return _call("set", setting=setting, value=value, selection=selection)


@mcp.tool()
def cartoon(item_type: str, selection: str | None = "all") -> str:
    """
    Sets the cartoon type for the specified selection

    Args:
        item_type: Cartoon style: ``automatic``, ``tube``, ``loop``, ``rectangle``, ``oval``,
            ``arrow``, ``dumbbell``, ``putty``, ``skip``. ``tube`` ignores secondary structure
            and runs unbroken through the backbone.
        selection: PyMOL selection string, e.g. "1abc and chain A". See the module note for
            syntax.
    """
    return _call("cartoon", item_type=item_type, selection=selection)


@mcp.tool()
def spectrum(
    expression: str, palette: str | None = "rainbow", selection: str | None = "all"
) -> str:
    """
    Colors selection in a spectrum

    Args:
        expression: Atom property to colour by: ``b`` (B-factor or pLDDT), ``q`` (occupancy),
            ``count``, ``resi``, ``pc`` (partial charge).
        palette: Colour ramp, e.g. ``rainbow``, ``blue_white_red``, ``cyan_white_magenta``,
            ``green_yellow_red``.
        selection: PyMOL selection string, e.g. "1abc and chain A". See the module note for
            syntax.
    """
    return _call("spectrum", expression=expression, palette=palette, selection=selection)


@mcp.tool()
def label(selection: str, expression: str | None = "name") -> str:
    """
    Adds labels to atoms in the selection

    Args:
        selection: Atoms to label. Use ``name CA`` to get one label per residue rather than one
            per atom.
        expression: Python expression over atom properties, e.g. ``"%s%s" % (resn, resi)`` or
            ``resi``.
    """
    return _call("label", selection=selection, expression=expression)


@mcp.tool()
def angle(
    name: str | None = None,
    selection1: str | None = "(pk1)",
    selection2: str | None = "(pk2)",
    selection3: str | None = "(pk3)",
) -> str:
    """Measures the angle between three selections and returns it in degrees.

    Also draws the measurement as a named angle object.

    Args:
        name: Name for the angle object.
        selection1: First selection (one arm).
        selection2: Second selection (the vertex).
        selection3: Third selection (the other arm).
    """
    return _call(
        "angle",
        _measures=("Angle", "degrees"),
        name=name,
        selection1=selection1,
        selection2=selection2,
        selection3=selection3,
    )


@mcp.tool()
def dihedral(
    name: str | None = None,
    selection1: str | None = "(pk1)",
    selection2: str | None = "(pk2)",
    selection3: str | None = "(pk3)",
    selection4: str | None = "(pk4)",
) -> str:
    """Measures the dihedral (torsion) angle between four selections, in degrees.

    Also draws the measurement as a named dihedral object. Useful for backbone
    phi/psi angles and ligand torsions.

    Args:
        name: Name for the dihedral object.
        selection1: First atom.
        selection2: Second atom.
        selection3: Third atom.
        selection4: Fourth atom.
    """
    return _call(
        "dihedral",
        _measures=("Dihedral", "degrees", True),
        name=name,
        selection1=selection1,
        selection2=selection2,
        selection3=selection3,
        selection4=selection4,
    )


@mcp.tool()
def center(selection: str | None = "all") -> str:
    """
    Centers the view on a selection

    Args:
        selection: What to centre the camera on.
    """
    return _call("center", selection=selection)


@mcp.tool()
def orient(selection: str | None = "all") -> str:
    """
    Orients the view to align with principal axes of the selection

    Args:
        selection: What to orient on. PyMOL aligns the longest axis of this selection with the
            screen's x-axis.
    """
    return _call("orient", selection=selection)


@mcp.tool()
def zoom(selection: str | None = "all", buffer: str | None = "5") -> str:
    """
    Zooms the view on a selection

    Args:
        selection: What to fill the view with.
        buffer: Extra padding around the selection, in Angstrom.
    """
    return _call("zoom", selection=selection, buffer=buffer)


@mcp.tool()
def reset(obj: str | None = None) -> str:
    """
    Resets the view, optionally resetting an object's matrix

    Args:
        obj: Object whose matrix to reset. Omit to reset only the camera.
    """
    return _call("reset", obj=obj)


@mcp.tool()
def turn(axis: str, angle: str | None = "90") -> str:
    """
    Rotates the camera around an axis

    Args:
        axis: Camera axis to rotate about: ``x``, ``y`` or ``z``.
        angle: Rotation in degrees. Negative reverses direction.
    """
    return _call("turn", axis=axis, angle=angle)


@mcp.tool()
def move(axis: str, distance: str | None = "1") -> str:
    """
    Moves the camera along an axis

    Args:
        axis: Camera axis to translate along: ``x``, ``y`` or ``z``.
        distance: Distance in Angstrom.
    """
    return _call("move", axis=axis, distance=distance)


@mcp.tool()
def clip(mode: str, distance: str | None = "1") -> str:
    """
    Adjusts the clipping planes

    Args:
        mode: Which plane to move: ``near``, ``far``, ``move`` (both), ``slab``, ``atoms``.
        distance: How far to move the plane, in Angstrom. Negative moves the other way.
    """
    return _call("clip", mode=mode, distance=distance)


@mcp.tool()
def save(filename: str, selection: str | None = "all", state: str | None = "-1") -> str:
    """
    Saves data to a file

    Args:
        filename: Output path. The extension picks the format: ``.pdb``, ``.cif``, ``.sdf``,
            ``.mol2``, ``.obj``, or ``.pse`` for a full session.
        selection: What to write out.
        state: Which state to save. "-1" is the current state, "0" writes all states.
    """
    return _call(
        "save", _timeout=_SLOW_OP_TIMEOUT, filename=filename, selection=selection, state=state
    )


@mcp.tool()
def png(filename: str, options: str | None = None) -> str:
    """
    Saves a PNG image

    Args:
        filename: Output path for the PNG.
        options: Extra arguments such as width, height, dpi, ray.
    """
    return _call("png", _timeout=_SLOW_OP_TIMEOUT, filename=filename, options=options)


@mcp.tool()
def deselect() -> str:
    """
    Clears the current selection
    """
    return _call("deselect")


@mcp.tool()
def create(name: str, selection: str | None = "all", source_state: str | None = "1") -> str:
    """
    Creates a new object from a selection

    Args:
        name: Name for the new object.
        selection: Atoms to copy into it. The original is left in place.
        source_state: State to copy from.
    """
    return _call("create", name=name, selection=selection, source_state=source_state)


@mcp.tool()
def extract(name: str, selection: str | None = "all") -> str:
    """
    Extracts a selection to a new object

    Args:
        name: Name for the new object.
        selection: Atoms to move into it. Unlike create, these are *removed* from the original
            object.
    """
    return _call("extract", name=name, selection=selection)


@mcp.tool()
def delete(name: str) -> str:
    """
    Deletes objects or selections

    Args:
        name: Object or named selection to delete. Accepts wildcards, e.g. ``tmp*``; ``all``
            clears everything.
    """
    return _call("delete", name=name)


@mcp.tool()
def align(mobile: str, target: str | None = "all", options: str | None = None) -> str:
    """
    Aligns one selection to another

    Args:
        mobile: Selection to move.
        target: Selection to align onto; this one stays put.
        options: Extra arguments, e.g. ``cycles=0`` to disable outlier rejection.
    """
    return _call("align", mobile=mobile, target=target, options=options)


@mcp.tool(name="super")
def super_tool(mobile: str, target: str | None = "all", options: str | None = None) -> str:
    """
    Superimposes one selection onto another

    Args:
        mobile: Selection to move.
        target: Selection to superimpose onto; this one stays put.
        options: Extra PyMOL arguments.
    """
    return _call("super", mobile=mobile, target=target, options=options)


@mcp.tool()
def intra_fit(selection: str) -> str:
    """
    Fits all states within an object

    Args:
        selection: Object whose states to fit onto its first state.
    """
    return _call("intra_fit", selection=selection)


@mcp.tool()
def intra_rms(selection: str) -> str:
    """
    Calculates RMSD between states within an object

    Args:
        selection: Object whose states to compare.
    """
    return _call("intra_rms", selection=selection)


@mcp.tool()
def alter(selection: str, expression: str) -> str:
    """
    Alters atomic properties in a selection

    Args:
        selection: Atoms to modify.
        expression: Assignment over atom properties, e.g. ``b=0``, ``chain='B'``,
            ``resi=str(int(resi)+100)``.
    """
    return _call("alter", selection=selection, expression=expression)


@mcp.tool()
def alter_state(state: str, selection: str, expression: str) -> str:
    """
    Alters atomic coordinates in a state

    Args:
        state: State to modify.
        selection: Atoms to modify.
        expression: Assignment over coordinates, e.g. ``x=x+10``.
    """
    return _call("alter_state", state=state, selection=selection, expression=expression)


@mcp.tool()
def h_add(selection: str | None = "all") -> str:
    """
    Adds hydrogens to a selection

    Args:
        selection: Where to add hydrogens.
    """
    return _call("h_add", selection=selection)


@mcp.tool()
def h_fill(selection: str | None = "all") -> str:
    """
    Adds hydrogens and adjusts valences

    Args:
        selection: Where to fill valences.
    """
    return _call("h_fill", selection=selection)


@mcp.tool()
def bond(atom1: str, atom2: str, order: str | None = "1") -> str:
    """
    Creates a bond between two atoms

    Args:
        atom1: First atom, as a selection matching exactly one atom.
        atom2: Second atom, as a selection matching exactly one atom.
        order: Bond order: ``1`` single, ``2`` double, ``3`` triple, ``4`` aromatic.
    """
    return _call("bond", atom1=atom1, atom2=atom2, order=order)


@mcp.tool()
def unbond(atom1: str, atom2: str) -> str:
    """
    Removes a bond between two atoms

    Args:
        atom1: First atom of the bond to remove.
        atom2: Second atom of the bond to remove.
    """
    return _call("unbond", atom1=atom1, atom2=atom2)


@mcp.tool()
def rebuild(selection: str | None = "all") -> str:
    """
    Regenerates all displayed geometry

    Args:
        selection: What to regenerate. Needed after altering coordinates or B-factors for the
            change to show.
    """
    return _call("rebuild", selection=selection)


@mcp.tool()
def refresh() -> str:
    """
    Refreshes the display
    """
    return _call("refresh")


@mcp.tool()
def util_cbc(selection: str | None = "all") -> str:
    """
    Colors by chain (Color By Chain)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbc", selection=selection)


@mcp.tool()
def util_cbaw(selection: str | None = "all") -> str:
    """
    Colors by atom, white carbons (Color By Atom, White)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbaw", selection=selection)


@mcp.tool()
def util_cbag(selection: str | None = "all") -> str:
    """
    Colors by atom, green carbons (Color By Atom, Green)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbag", selection=selection)


@mcp.tool()
def util_cbac(selection: str | None = "all") -> str:
    """
    Colors by atom, cyan carbons (Color By Atom, Cyan)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbac", selection=selection)


@mcp.tool()
def util_cbam(selection: str | None = "all") -> str:
    """
    Colors by atom, magenta carbons (Color By Atom, Magenta)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbam", selection=selection)


@mcp.tool()
def util_cbay(selection: str | None = "all") -> str:
    """
    Colors by atom, yellow carbons (Color By Atom, Yellow)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbay", selection=selection)


@mcp.tool()
def util_cbas(selection: str | None = "all") -> str:
    """
    Colors by atom, salmon carbons (Color By Atom, Salmon)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbas", selection=selection)


@mcp.tool()
def util_cbab(selection: str | None = "all") -> str:
    """
    Colors by atom, slate carbons (Color By Atom, slateBLue)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbab", selection=selection)


@mcp.tool()
def util_cbao(selection: str | None = "all") -> str:
    """
    Colors by atom, orange carbons (Color By Atom, Orange)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbao", selection=selection)


@mcp.tool()
def util_cbap(selection: str | None = "all") -> str:
    """
    Colors by atom, purple carbons (Color By Atom, Purple)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbap", selection=selection)


@mcp.tool()
def util_cbak(selection: str | None = "all") -> str:
    """
    Colors by atom, pink carbons (Color By Atom, pinK)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.cbak", selection=selection)


@mcp.tool()
def util_chainbow(selection: str | None = "all") -> str:
    """
    Colors chains in rainbow gradient (CHAINs in rainBOW)

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.chainbow", selection=selection)


@mcp.tool()
def util_rainbow(selection: str | None = "all") -> str:
    """
    Colors residues in rainbow from N to C terminus

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.rainbow", selection=selection)


@mcp.tool()
def util_ss(selection: str | None = "all") -> str:
    """
    Colors by secondary structure

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.ss", selection=selection)


@mcp.tool()
def util_color_by_element(selection: str | None = "all") -> str:
    """
    Colors atoms by their element

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.color_by_element", selection=selection)


@mcp.tool()
def util_color_secondary(selection: str | None = "all") -> str:
    """
    Colors secondary structure elements

    Args:
        selection: What to recolour. Defaults to everything.
    """
    return _call("util.color_secondary", selection=selection)


@mcp.tool()
def spheroid(selection: str | None = "all") -> str:
    """
    Displays atoms as smooth spheres

    Args:
        selection: What to render as smoothed spheres.
    """
    return _call("spheroid", selection=selection)


@mcp.tool()
def isomesh(name: str, map_object: str, level: str, selection: str | None = "all") -> str:
    """
    Creates a mesh isosurface

    Args:
        name: Name for the mesh object.
        map_object: Name of a loaded map (e.g. a CCP4 or DX density map).
        level: Contour level in sigma, e.g. ``1.0`` for 2Fo-Fc density.
        selection: Restrict the mesh to the region around this selection.
    """
    return _call("isomesh", name=name, map_object=map_object, level=level, selection=selection)


@mcp.tool()
def isosurface(name: str, map_object: str, level: str, selection: str | None = "all") -> str:
    """
    Creates a solid isosurface

    Args:
        name: Name for the surface object.
        map_object: Name of a loaded map.
        level: Contour level in sigma.
        selection: Restrict the surface to the region around this selection.
    """
    return _call("isosurface", name=name, map_object=map_object, level=level, selection=selection)


@mcp.tool()
def sculpt_activate(obj: str) -> str:
    """
    Activates sculpting mode for an object

    Args:
        obj: Object to enable real-space sculpting on.
    """
    return _call("sculpt_activate", obj=obj)


@mcp.tool()
def sculpt_deactivate(obj: str) -> str:
    """
    Deactivates sculpting mode for an object

    Args:
        obj: Object to stop sculpting.
    """
    return _call("sculpt_deactivate", obj=obj)


@mcp.tool()
def sculpt_iterate(iterations: str, obj: str | None = "all") -> str:
    """
    Performs sculpting iterations

    Args:
        iterations: Number of relaxation cycles to run.
        obj: Object to relax.
    """
    return _call("sculpt_iterate", iterations=iterations, obj=obj)


@mcp.tool()
def scene(key: str, action: str | None = "recall") -> str:
    """
    Manages scenes for later recall

    Args:
        key: Scene name, e.g. "F1", or "auto" to use the next free slot.
        action: ``store`` to save the current view, ``recall`` to restore it, ``clear``,
            ``update``, ``rename``, ``delete``, ``next``, ``previous``.
    """
    return _call("scene", key=key, action=action)


@mcp.tool()
def scene_order(scene_list: str) -> str:
    """
    Sets the order of scenes

    Args:
        scene_list: Space-separated scene names in the order wanted, e.g. "F2 F1 F3".
    """
    return _call("scene_order", scene_list=scene_list)


@mcp.tool()
def mset(specification: str) -> str:
    """
    Defines a sequence of states for movie playback

    Args:
        specification: Frame-to-state mapping, e.g. "1 x30" to hold state 1 for 30 frames, or "1
            -30" to sweep states 1 through 30.
    """
    return _call("mset", specification=specification)


@mcp.tool()
def mplay() -> str:
    """
    Starts playing the movie
    """
    return _call("mplay")


@mcp.tool()
def mstop() -> str:
    """
    Stops the movie
    """
    return _call("mstop")


@mcp.tool()
def frame(frame_number: str | None = None) -> str:
    """
    Sets or queries the current frame

    Args:
        frame_number: Frame to jump to. Omit to query the current frame.
    """
    return _call("frame", frame_number=frame_number)


@mcp.tool()
def forward() -> str:
    """
    Advances one frame
    """
    return _call("forward")


@mcp.tool()
def backward() -> str:
    """
    Goes back one frame
    """
    return _call("backward")


@mcp.tool()
def rock() -> str:
    """
    Toggles a rocking animation
    """
    return _call("rock")


@mcp.tool()
def ray(width: str | None = None, height: str | None = None) -> str:
    """
    Performs ray-tracing

    Args:
        width: Image width in pixels. Omit to use the viewport size.
        height: Image height in pixels.
    """
    return _call("ray", _timeout=_SLOW_OP_TIMEOUT, width=width, height=height)


@mcp.tool()
def draw(width: str | None = None, height: str | None = None) -> str:
    """
    Uses OpenGL renderer (faster but lower quality)

    Args:
        width: Image width in pixels.
        height: Image height in pixels.
    """
    return _call("draw", _timeout=_SLOW_OP_TIMEOUT, width=width, height=height)


@mcp.tool()
def mpng(prefix: str) -> str:
    """
    Saves a series of PNG images for movie frames

    Args:
        prefix: Filename prefix; PyMOL appends a zero-padded frame number.
    """
    return _call("mpng", _timeout=_SLOW_OP_TIMEOUT, prefix=prefix)


@mcp.tool()
def symexp(prefix: str, selection: str, cutoff: str | None = "20", segi: str | None = None) -> str:
    """
    Generates symmetry-related copies

    Args:
        prefix: Prefix for the generated symmetry-mate objects.
        selection: Selection whose crystallographic neighbours to build.
        cutoff: Distance in Angstrom out to which to generate mates.
        segi: Optional segment identifier for the new objects.
    """
    return _call("symexp", prefix=prefix, selection=selection, cutoff=cutoff, segi=segi)


@mcp.tool()
def set_symmetry(selection: str, a: str, b: str, c: str, alpha: str, beta: str, gamma: str) -> str:
    """
    Sets symmetry parameters for an object

    Args:
        selection: Object to assign the unit cell to.
        a: Unit cell edge a, in Angstrom.
        b: Unit cell edge b, in Angstrom.
        c: Unit cell edge c, in Angstrom.
        alpha: Unit cell angle alpha, in degrees.
        beta: Unit cell angle beta, in degrees.
        gamma: Unit cell angle gamma, in degrees.
    """
    return _call(
        "set_symmetry", selection=selection, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma
    )


@mcp.tool()
def fab(sequence: str, options: str | None = None) -> str:
    """
    Creates a peptide chain from a sequence

    Args:
        sequence: One-letter amino acid sequence, e.g. "ACDEFGH".
        options: Extra arguments, e.g. ``ss=1`` to build it as a helix.
    """
    return _call("fab", sequence=sequence, options=options)


@mcp.tool()
def fragment(name: str) -> str:
    """
    Loads a molecular fragment

    Args:
        name: Built-in fragment name, e.g. "benzene", "ala", "formamide".
    """
    return _call("fragment", name=name)


@mcp.tool()
def full_screen() -> str:
    """
    Toggles fullscreen mode
    """
    return _call("full_screen")


@mcp.tool()
def viewport(width: str, height: str) -> str:
    """
    Sets the viewport size

    Args:
        width: Viewport width in pixels.
        height: Viewport height in pixels.
    """
    return _call("viewport", width=width, height=height)


@mcp.tool()
def cd(path: str) -> str:
    """
    Changes the current directory

    Args:
        path: Directory to change into. PyMOL resolves relative paths from here.
    """
    return _call("cd", path=path)


@mcp.tool()
def pwd() -> str:
    """
    Prints the current directory
    """
    return _call("pwd")


@mcp.tool()
def ls(path: str | None = None) -> str:
    """
    Lists files in the current directory

    Args:
        path: Directory or glob to list. Omit for the current directory.
    """
    return _call("ls", path=path)


@mcp.tool()
def system(command: str) -> str:
    """
    Executes a system command

    Args:
        command: Shell command to run on the machine PyMOL is running on.
    """
    return _call("system", command=command)


@mcp.tool()
def help(command: str | None = None) -> str:
    """
    Shows help for a command

    Args:
        command: PyMOL command to describe. Omit for general help.
    """
    return _call("help", command=command)
