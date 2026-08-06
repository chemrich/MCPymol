"""Thin wrappers over individual ``pymol.cmd`` functions.

One tool per PyMOL command, so the model can compose finer motions when the
high-level views in ``mcpymol.views`` do not quite fit.  Each wrapper is a real
``def`` — its signature and docstring are what the MCP client sees — delegating
to :func:`mcpymol.bridge._call` for the shared forwarding logic.
"""

from mcpymol.app import mcp
from mcpymol.bridge import _SLOW_OP_TIMEOUT, _call, send_request

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
    """Permanently removes the atoms matching the selection.

    This deletes atoms; it does not just hide them. To hide instead, use
    :func:`hide`.
    """
    res = send_request("remove", args=[selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Removed selection '{selection}'"


@mcp.tool()
def distance(name: str, selection1: str, selection2: str) -> str:
    """Measures and draws a distance object between two selections.

    Args:
        name: Name for the distance object (used to delete/recolor later).
        selection1: First selection.
        selection2: Second selection.
    """
    res = send_request("distance", args=[name, selection1, selection2])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Measured distance between '{selection1}' and '{selection2}' as '{name}'"


@mcp.tool()
def execute_pymol_command(command: str) -> str:
    """Executes a raw PyMOL command string (PyMOL CLI syntax).

    PREFER the dedicated tools when one exists — show, color, select, distance,
    ligand_view, interface_view, etc. They have better defaults, do compound
    setup in one call, and produce cleaner results.

    Reach for this tool only when no other tool covers what you need (e.g.
    ``set ray_shadow, 0``, ``bg_color grey20``, multi-statement scripts).
    Note: this accepts the PyMOL ``cmd.do`` mini-language, not Python.
    """
    res = send_request("do", args=[command])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Executed command: {command}"


@mcp.tool(name="as")
def as_tool(representation: str, selection: str | None = "all") -> str:
    """
    Shows one representation while hiding all others for the specified selection
    """
    return _call("as", representation=representation, selection=selection)


@mcp.tool(name="set")
def set_setting(setting: str, value: str, selection: str | None = None) -> str:
    """
    Sets a PyMOL setting to a specified value
    """
    return _call("set", setting=setting, value=value, selection=selection)


@mcp.tool()
def cartoon(item_type: str, selection: str | None = "all") -> str:
    """
    Sets the cartoon type for the specified selection
    """
    return _call("cartoon", item_type=item_type, selection=selection)


@mcp.tool()
def spectrum(
    expression: str, palette: str | None = "rainbow", selection: str | None = "all"
) -> str:
    """
    Colors selection in a spectrum
    """
    return _call("spectrum", expression=expression, palette=palette, selection=selection)


@mcp.tool()
def label(selection: str, expression: str | None = "name") -> str:
    """
    Adds labels to atoms in the selection
    """
    return _call("label", selection=selection, expression=expression)


@mcp.tool()
def angle(
    name: str | None = None,
    selection1: str | None = "(pk1)",
    selection2: str | None = "(pk2)",
    selection3: str | None = "(pk3)",
) -> str:
    """
    Measures the angle between three selections
    """
    return _call(
        "angle", name=name, selection1=selection1, selection2=selection2, selection3=selection3
    )


@mcp.tool()
def dihedral(
    name: str | None = None,
    selection1: str | None = "(pk1)",
    selection2: str | None = "(pk2)",
    selection3: str | None = "(pk3)",
    selection4: str | None = "(pk4)",
) -> str:
    """
    Measures the dihedral angle between four selections
    """
    return _call(
        "dihedral",
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
    """
    return _call("center", selection=selection)


@mcp.tool()
def orient(selection: str | None = "all") -> str:
    """
    Orients the view to align with principal axes of the selection
    """
    return _call("orient", selection=selection)


@mcp.tool()
def zoom(selection: str | None = "all", buffer: str | None = "5") -> str:
    """
    Zooms the view on a selection
    """
    return _call("zoom", selection=selection, buffer=buffer)


@mcp.tool()
def reset(obj: str | None = None) -> str:
    """
    Resets the view, optionally resetting an object's matrix
    """
    return _call("reset", obj=obj)


@mcp.tool()
def turn(axis: str, angle: str | None = "90") -> str:
    """
    Rotates the camera around an axis
    """
    return _call("turn", axis=axis, angle=angle)


@mcp.tool()
def move(axis: str, distance: str | None = "1") -> str:
    """
    Moves the camera along an axis
    """
    return _call("move", axis=axis, distance=distance)


@mcp.tool()
def clip(mode: str, distance: str | None = "1") -> str:
    """
    Adjusts the clipping planes
    """
    return _call("clip", mode=mode, distance=distance)


@mcp.tool()
def save(filename: str, selection: str | None = "all", state: str | None = "-1") -> str:
    """
    Saves data to a file
    """
    return _call(
        "save", _timeout=_SLOW_OP_TIMEOUT, filename=filename, selection=selection, state=state
    )


@mcp.tool()
def png(filename: str, options: str | None = None) -> str:
    """
    Saves a PNG image
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
    """
    return _call("create", name=name, selection=selection, source_state=source_state)


@mcp.tool()
def extract(name: str, selection: str | None = "all") -> str:
    """
    Extracts a selection to a new object
    """
    return _call("extract", name=name, selection=selection)


@mcp.tool()
def delete(name: str) -> str:
    """
    Deletes objects or selections
    """
    return _call("delete", name=name)


@mcp.tool()
def align(mobile: str, target: str | None = "all", options: str | None = None) -> str:
    """
    Aligns one selection to another
    """
    return _call("align", mobile=mobile, target=target, options=options)


@mcp.tool(name="super")
def super_tool(mobile: str, target: str | None = "all", options: str | None = None) -> str:
    """
    Superimposes one selection onto another
    """
    return _call("super", mobile=mobile, target=target, options=options)


@mcp.tool()
def intra_fit(selection: str) -> str:
    """
    Fits all states within an object
    """
    return _call("intra_fit", selection=selection)


@mcp.tool()
def intra_rms(selection: str) -> str:
    """
    Calculates RMSD between states within an object
    """
    return _call("intra_rms", selection=selection)


@mcp.tool()
def alter(selection: str, expression: str) -> str:
    """
    Alters atomic properties in a selection
    """
    return _call("alter", selection=selection, expression=expression)


@mcp.tool()
def alter_state(state: str, selection: str, expression: str) -> str:
    """
    Alters atomic coordinates in a state
    """
    return _call("alter_state", state=state, selection=selection, expression=expression)


@mcp.tool()
def h_add(selection: str | None = "all") -> str:
    """
    Adds hydrogens to a selection
    """
    return _call("h_add", selection=selection)


@mcp.tool()
def h_fill(selection: str | None = "all") -> str:
    """
    Adds hydrogens and adjusts valences
    """
    return _call("h_fill", selection=selection)


@mcp.tool()
def bond(atom1: str, atom2: str, order: str | None = "1") -> str:
    """
    Creates a bond between two atoms
    """
    return _call("bond", atom1=atom1, atom2=atom2, order=order)


@mcp.tool()
def unbond(atom1: str, atom2: str) -> str:
    """
    Removes a bond between two atoms
    """
    return _call("unbond", atom1=atom1, atom2=atom2)


@mcp.tool()
def rebuild(selection: str | None = "all") -> str:
    """
    Regenerates all displayed geometry
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
    """
    return _call("util.cbc", selection=selection)


@mcp.tool()
def util_cbaw(selection: str | None = "all") -> str:
    """
    Colors by atom, white carbons (Color By Atom, White)
    """
    return _call("util.cbaw", selection=selection)


@mcp.tool()
def util_cbag(selection: str | None = "all") -> str:
    """
    Colors by atom, green carbons (Color By Atom, Green)
    """
    return _call("util.cbag", selection=selection)


@mcp.tool()
def util_cbac(selection: str | None = "all") -> str:
    """
    Colors by atom, cyan carbons (Color By Atom, Cyan)
    """
    return _call("util.cbac", selection=selection)


@mcp.tool()
def util_cbam(selection: str | None = "all") -> str:
    """
    Colors by atom, magenta carbons (Color By Atom, Magenta)
    """
    return _call("util.cbam", selection=selection)


@mcp.tool()
def util_cbay(selection: str | None = "all") -> str:
    """
    Colors by atom, yellow carbons (Color By Atom, Yellow)
    """
    return _call("util.cbay", selection=selection)


@mcp.tool()
def util_cbas(selection: str | None = "all") -> str:
    """
    Colors by atom, salmon carbons (Color By Atom, Salmon)
    """
    return _call("util.cbas", selection=selection)


@mcp.tool()
def util_cbab(selection: str | None = "all") -> str:
    """
    Colors by atom, slate carbons (Color By Atom, slateBLue)
    """
    return _call("util.cbab", selection=selection)


@mcp.tool()
def util_cbao(selection: str | None = "all") -> str:
    """
    Colors by atom, orange carbons (Color By Atom, Orange)
    """
    return _call("util.cbao", selection=selection)


@mcp.tool()
def util_cbap(selection: str | None = "all") -> str:
    """
    Colors by atom, purple carbons (Color By Atom, Purple)
    """
    return _call("util.cbap", selection=selection)


@mcp.tool()
def util_cbak(selection: str | None = "all") -> str:
    """
    Colors by atom, pink carbons (Color By Atom, pinK)
    """
    return _call("util.cbak", selection=selection)


@mcp.tool()
def util_chainbow(selection: str | None = "all") -> str:
    """
    Colors chains in rainbow gradient (CHAINs in rainBOW)
    """
    return _call("util.chainbow", selection=selection)


@mcp.tool()
def util_rainbow(selection: str | None = "all") -> str:
    """
    Colors residues in rainbow from N to C terminus
    """
    return _call("util.rainbow", selection=selection)


@mcp.tool()
def util_ss(selection: str | None = "all") -> str:
    """
    Colors by secondary structure
    """
    return _call("util.ss", selection=selection)


@mcp.tool()
def util_color_by_element(selection: str | None = "all") -> str:
    """
    Colors atoms by their element
    """
    return _call("util.color_by_element", selection=selection)


@mcp.tool()
def util_color_secondary(selection: str | None = "all") -> str:
    """
    Colors secondary structure elements
    """
    return _call("util.color_secondary", selection=selection)


@mcp.tool()
def spheroid(selection: str | None = "all") -> str:
    """
    Displays atoms as smooth spheres
    """
    return _call("spheroid", selection=selection)


@mcp.tool()
def isomesh(name: str, map_object: str, level: str, selection: str | None = "all") -> str:
    """
    Creates a mesh isosurface
    """
    return _call("isomesh", name=name, map_object=map_object, level=level, selection=selection)


@mcp.tool()
def isosurface(name: str, map_object: str, level: str, selection: str | None = "all") -> str:
    """
    Creates a solid isosurface
    """
    return _call("isosurface", name=name, map_object=map_object, level=level, selection=selection)


@mcp.tool()
def sculpt_activate(obj: str) -> str:
    """
    Activates sculpting mode for an object
    """
    return _call("sculpt_activate", obj=obj)


@mcp.tool()
def sculpt_deactivate(obj: str) -> str:
    """
    Deactivates sculpting mode for an object
    """
    return _call("sculpt_deactivate", obj=obj)


@mcp.tool()
def sculpt_iterate(iterations: str, obj: str | None = "all") -> str:
    """
    Performs sculpting iterations
    """
    return _call("sculpt_iterate", iterations=iterations, obj=obj)


@mcp.tool()
def scene(key: str, action: str | None = "recall") -> str:
    """
    Manages scenes for later recall
    """
    return _call("scene", key=key, action=action)


@mcp.tool()
def scene_order(scene_list: str) -> str:
    """
    Sets the order of scenes
    """
    return _call("scene_order", scene_list=scene_list)


@mcp.tool()
def mset(specification: str) -> str:
    """
    Defines a sequence of states for movie playback
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
    """
    return _call("ray", _timeout=_SLOW_OP_TIMEOUT, width=width, height=height)


@mcp.tool()
def draw(width: str | None = None, height: str | None = None) -> str:
    """
    Uses OpenGL renderer (faster but lower quality)
    """
    return _call("draw", _timeout=_SLOW_OP_TIMEOUT, width=width, height=height)


@mcp.tool()
def mpng(prefix: str) -> str:
    """
    Saves a series of PNG images for movie frames
    """
    return _call("mpng", _timeout=_SLOW_OP_TIMEOUT, prefix=prefix)


@mcp.tool()
def symexp(prefix: str, selection: str, cutoff: str | None = "20", segi: str | None = None) -> str:
    """
    Generates symmetry-related copies
    """
    return _call("symexp", prefix=prefix, selection=selection, cutoff=cutoff, segi=segi)


@mcp.tool()
def set_symmetry(selection: str, a: str, b: str, c: str, alpha: str, beta: str, gamma: str) -> str:
    """
    Sets symmetry parameters for an object
    """
    return _call(
        "set_symmetry", selection=selection, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma
    )


@mcp.tool()
def fab(sequence: str, options: str | None = None) -> str:
    """
    Creates a peptide chain from a sequence
    """
    return _call("fab", sequence=sequence, options=options)


@mcp.tool()
def fragment(name: str) -> str:
    """
    Loads a molecular fragment
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
    """
    return _call("viewport", width=width, height=height)


@mcp.tool()
def cd(path: str) -> str:
    """
    Changes the current directory
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
    """
    return _call("ls", path=path)


@mcp.tool()
def system(command: str) -> str:
    """
    Executes a system command
    """
    return _call("system", command=command)


@mcp.tool()
def help(command: str | None = None) -> str:
    """
    Shows help for a command
    """
    return _call("help", command=command)
