"""Thin wrappers over individual ``pymol.cmd`` functions.

One tool per PyMOL command, so the model can compose finer motions when the
high-level views in ``mcpymol.views`` do not quite fit.  Each wrapper is a real
``def`` — its signature and docstring are what the MCP client sees — delegating
to :func:`mcpymol.bridge._call` for the shared forwarding logic.
"""

from typing import Annotated

from pydantic import Field

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
def show(
    representation: Annotated[
        str,
        Field(
            description="One of ``cartoon``, ``sticks``, ``spheres``, ``surface``, ``lines``, ``ribbon``, ``dots``, ``mesh``, ``nb_spheres``, ``labels``."
        ),
    ],
    selection: Annotated[
        str, Field(description="PyMOL selection string. See module-level note for syntax.")
    ] = "all",
) -> str:
    """
    Shows a graphical representation for a given selection.
    """
    res = send_request("show", args=[representation, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Showing {representation} for selection '{selection}'"


@mcp.tool()
def hide(
    representation: Annotated[
        str,
        Field(
            description="Same vocabulary as :func:`show`, plus ``everything`` to hide all reps for the selection."
        ),
    ],
    selection: Annotated[str, Field(description="PyMOL selection string.")] = "all",
) -> str:
    """
    Hides a graphical representation for a given selection.
    """
    res = send_request("hide", args=[representation, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Hiding {representation} for selection '{selection}'"


@mcp.tool()
def color(
    color_name: Annotated[
        str,
        Field(
            description="A PyMOL color. Common names: ``red``, ``blue``, ``green``, ``yellow``, ``magenta``, ``cyan``, ``orange``, ``salmon``, ``marine``, ``forest``, ``palegreen``, ``skyblue``, ``violet``, ``grey50``, ``white``, ``black``. Use ``atomic`` to color non-carbon atoms by element while leaving carbons untouched."
        ),
    ],
    selection: Annotated[str, Field(description="PyMOL selection string.")] = "all",
) -> str:
    """
    Sets the color for a selection.
    """
    res = send_request("color", args=[color_name, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Colored selection '{selection}' with {color_name}"


@mcp.tool()
def select(
    name: Annotated[
        str, Field(description="Identifier you'll refer to later (e.g. ``active_site``).")
    ],
    selection: Annotated[
        str, Field(description="PyMOL selection expression to assign to that name.")
    ],
) -> str:
    """
    Creates (or replaces) a named selection for later reuse.
    """
    res = send_request("select", args=[name, selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Created named selection '{name}' for '{selection}'"


@mcp.tool()
def remove(
    selection: Annotated[
        str, Field(description="Atoms to delete permanently. To hide them instead, use hide.")
    ],
) -> str:
    """
    Permanently removes the atoms matching the selection. This deletes atoms; it does not just
    hide them. To hide instead, use :func:`hide`.
    """
    res = send_request("remove", args=[selection])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Removed selection '{selection}'"


@mcp.tool()
def distance(
    name: Annotated[
        str, Field(description="Name for the distance object (used to delete/recolor later).")
    ],
    selection1: Annotated[str, Field(description="First selection.")],
    selection2: Annotated[str, Field(description="Second selection.")],
) -> str:
    """
    Measures the distance between two selections and returns it in Angstrom.

    Also draws the measurement in the viewport as a named distance object.
    With multi-atom selections PyMOL reports the average over the pairs it
    found within its default cutoff.
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
def sasa(
    selection: Annotated[
        str, Field(description='What to measure (e.g. "1brs and chain A").')
    ] = "all",
    state: Annotated[
        str | None, Field(description='Object state to measure. "1" is the first/only state.')
    ] = "1",
) -> str:
    """
    Measures solvent-accessible surface area, in square Angstrom.

    Reports the SASA of ``selection`` *in the context of the object it belongs
    to* — so a chain measured inside a complex is already partly occluded by
    its partner. To get the free (unbound) area, copy the chain to its own
    object first with ``create``, or use ``interface_report``, which does the
    bound/free bookkeeping for you.

    Accuracy depends on PyMOL's ``dot_solvent`` (0 = molecular surface,
    1 = solvent-accessible) and ``dot_density`` settings.
    """
    return _call(
        "get_area",
        _measures=("Solvent-accessible surface area", "A^2"),
        selection=selection,
        state=state,
    )


@mcp.tool()
def rms_cur(
    mobile: Annotated[str, Field(description="First selection.")],
    target: Annotated[str, Field(description="Second selection, same number of atoms.")],
) -> str:
    """
    Measures RMSD between two selections *without* moving anything.

    Use this when the structures are already superposed, or when you want to
    know how far apart they are as currently positioned. Compare with ``align``
    and ``super``, which move the mobile structure to minimise RMSD before
    reporting it, and with ``superposition_view``, which shows where the
    difference is rather than summarising it as one number.

    Requires the two selections to have matching atom counts.
    """
    return _call("rms_cur", _measures=("RMSD (as positioned)", "A"), mobile=mobile, target=target)


@mcp.tool()
def count_atoms(
    selection: Annotated[str, Field(description="PyMOL selection string to count.")] = "all",
) -> str:
    """
    Counts the atoms matching a selection.

    Handy for checking a selection expression does what you think before
    building a scene on top of it — an empty count means the expression is
    wrong, which is otherwise invisible until the picture comes out blank.
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
def execute_pymol_command(
    command: Annotated[
        str, Field(description="A PyMOL CLI command, e.g. `set ray_shadow, 0`. Not Python.")
    ],
) -> str:
    """
    Executes a raw PyMOL command string (PyMOL CLI syntax). PREFER the dedicated tools when one
    exists — show, color, select, distance, ligand_view, interface_view, etc. They have better
    defaults, do compound setup in one call, and produce cleaner results. Reach for this tool
    only when no other tool covers what you need (e.g. ``set ray_shadow, 0``, ``bg_color
    grey20``, multi-statement scripts). Note: this accepts the PyMOL ``cmd.do`` mini-language,
    not Python.
    """
    res = send_request("do", args=[command])
    if res.get("status") == "error":
        return res.get("error", "Unknown error")
    return f"Executed command: {command}"


@mcp.tool(name="as")
def as_tool(
    representation: Annotated[
        str,
        Field(
            description="The one representation to leave visible: ``cartoon``, ``sticks``, ``spheres``, ``surface``, ``lines``, ``ribbon``, ``dots``, ``mesh``."
        ),
    ],
    selection: Annotated[
        str | None,
        Field(
            description='PyMOL selection string, e.g. "1abc and chain A". See the module note for syntax.'
        ),
    ] = "all",
) -> str:
    """
    Shows one representation while hiding all others for the specified selection
    """
    return _call("as", representation=representation, selection=selection)


@mcp.tool(name="set")
def set_setting(
    setting: Annotated[
        str,
        Field(
            description='PyMOL setting name, e.g. "cartoon_transparency", "ray_shadow", "sphere_scale". Hundreds exist; see the PyMOL settings reference.'
        ),
    ],
    value: Annotated[str, Field(description='The value as a string, e.g. "0.5", "1", "black".')],
    selection: Annotated[
        str | None,
        Field(
            description="Limit the setting to this object or selection. Omit to set it globally."
        ),
    ] = None,
) -> str:
    """
    Sets a PyMOL setting to a specified value
    """
    return _call("set", setting=setting, value=value, selection=selection)


@mcp.tool()
def cartoon(
    item_type: Annotated[
        str,
        Field(
            description="Cartoon style: ``automatic``, ``tube``, ``loop``, ``rectangle``, ``oval``, ``arrow``, ``dumbbell``, ``putty``, ``skip``. ``tube`` ignores secondary structure and runs unbroken through the backbone."
        ),
    ],
    selection: Annotated[
        str | None,
        Field(
            description='PyMOL selection string, e.g. "1abc and chain A". See the module note for syntax.'
        ),
    ] = "all",
) -> str:
    """
    Sets the cartoon type for the specified selection
    """
    return _call("cartoon", item_type=item_type, selection=selection)


@mcp.tool()
def spectrum(
    expression: Annotated[
        str,
        Field(
            description="Atom property to colour by: ``b`` (B-factor or pLDDT), ``q`` (occupancy), ``count``, ``resi``, ``pc`` (partial charge)."
        ),
    ],
    palette: Annotated[
        str | None,
        Field(
            description="Colour ramp, e.g. ``rainbow``, ``blue_white_red``, ``cyan_white_magenta``, ``green_yellow_red``."
        ),
    ] = "rainbow",
    selection: Annotated[
        str | None,
        Field(
            description='PyMOL selection string, e.g. "1abc and chain A". See the module note for syntax.'
        ),
    ] = "all",
) -> str:
    """
    Colors selection in a spectrum
    """
    return _call("spectrum", expression=expression, palette=palette, selection=selection)


@mcp.tool()
def label(
    selection: Annotated[
        str,
        Field(
            description="Atoms to label. Use ``name CA`` to get one label per residue rather than one per atom."
        ),
    ],
    expression: Annotated[
        str | None,
        Field(
            description='Python expression over atom properties, e.g. ``"%s%s" % (resn, resi)`` or ``resi``.'
        ),
    ] = "name",
) -> str:
    """
    Adds labels to atoms in the selection
    """
    return _call("label", selection=selection, expression=expression)


@mcp.tool()
def angle(
    name: Annotated[str | None, Field(description="Name for the angle object.")] = None,
    selection1: Annotated[str | None, Field(description="First selection (one arm).")] = "(pk1)",
    selection2: Annotated[
        str | None, Field(description="Second selection (the vertex).")
    ] = "(pk2)",
    selection3: Annotated[
        str | None, Field(description="Third selection (the other arm).")
    ] = "(pk3)",
) -> str:
    """
    Measures the angle between three selections and returns it in degrees.

    Also draws the measurement as a named angle object.
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
    name: Annotated[str | None, Field(description="Name for the dihedral object.")] = None,
    selection1: Annotated[
        str | None,
        Field(
            description="First atom of the torsion. For a backbone phi angle this is the preceding C."
        ),
    ] = "(pk1)",
    selection2: Annotated[
        str | None,
        Field(description="Second atom; the first of the two forming the rotatable bond."),
    ] = "(pk2)",
    selection3: Annotated[
        str | None,
        Field(description="Third atom; the second of the two forming the rotatable bond."),
    ] = "(pk3)",
    selection4: Annotated[
        str | None, Field(description="Fourth atom, at the far end of the torsion.")
    ] = "(pk4)",
) -> str:
    """
    Measures the dihedral (torsion) angle between four selections, in degrees.

    Also draws the measurement as a named dihedral object. Useful for backbone
    phi/psi angles and ligand torsions.
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
def center(
    selection: Annotated[str | None, Field(description="What to centre the camera on.")] = "all",
) -> str:
    """
    Centers the view on a selection
    """
    return _call("center", selection=selection)


@mcp.tool()
def orient(
    selection: Annotated[
        str | None,
        Field(
            description="What to orient on. PyMOL aligns the longest axis of this selection with the screen's x-axis."
        ),
    ] = "all",
) -> str:
    """
    Orients the view to align with principal axes of the selection
    """
    return _call("orient", selection=selection)


@mcp.tool()
def zoom(
    selection: Annotated[str | None, Field(description="What to fill the view with.")] = "all",
    buffer: Annotated[
        str | None, Field(description="Extra padding around the selection, in Angstrom.")
    ] = "5",
) -> str:
    """
    Zooms the view on a selection
    """
    return _call("zoom", selection=selection, buffer=buffer)


@mcp.tool()
def reset(
    obj: Annotated[
        str | None,
        Field(description="Object whose matrix to reset. Omit to reset only the camera."),
    ] = None,
) -> str:
    """
    Resets the view, optionally resetting an object's matrix
    """
    return _call("reset", obj=obj)


@mcp.tool()
def turn(
    axis: Annotated[str, Field(description="Camera axis to rotate about: ``x``, ``y`` or ``z``.")],
    angle: Annotated[
        str | None, Field(description="Rotation in degrees. Negative reverses direction.")
    ] = "90",
) -> str:
    """
    Rotates the camera around an axis
    """
    return _call("turn", axis=axis, angle=angle)


@mcp.tool()
def move(
    axis: Annotated[
        str, Field(description="Camera axis to translate along: ``x``, ``y`` or ``z``.")
    ],
    distance: Annotated[str | None, Field(description="Distance in Angstrom.")] = "1",
) -> str:
    """
    Moves the camera along an axis
    """
    return _call("move", axis=axis, distance=distance)


@mcp.tool()
def clip(
    mode: Annotated[
        str,
        Field(
            description="Which plane to move: ``near``, ``far``, ``move`` (both), ``slab``, ``atoms``."
        ),
    ],
    distance: Annotated[
        str | None,
        Field(description="How far to move the plane, in Angstrom. Negative moves the other way."),
    ] = "1",
) -> str:
    """
    Adjusts the clipping planes
    """
    return _call("clip", mode=mode, distance=distance)


@mcp.tool()
def save(
    filename: Annotated[
        str,
        Field(
            description="Output path. The extension picks the format: ``.pdb``, ``.cif``, ``.sdf``, ``.mol2``, ``.obj``, or ``.pse`` for a full session."
        ),
    ],
    selection: Annotated[str | None, Field(description="What to write out.")] = "all",
    state: Annotated[
        str | None,
        Field(description='Which state to save. "-1" is the current state, "0" writes all states.'),
    ] = "-1",
) -> str:
    """
    Saves data to a file
    """
    return _call(
        "save", _timeout=_SLOW_OP_TIMEOUT, filename=filename, selection=selection, state=state
    )


@mcp.tool()
def png(
    filename: Annotated[str, Field(description="Output path for the PNG.")],
    options: Annotated[
        str | None, Field(description="Extra arguments such as width, height, dpi, ray.")
    ] = None,
) -> str:
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
def create(
    name: Annotated[str, Field(description="Name for the new object.")],
    selection: Annotated[
        str | None, Field(description="Atoms to copy into it. The original is left in place.")
    ] = "all",
    source_state: Annotated[str | None, Field(description="State to copy from.")] = "1",
) -> str:
    """
    Creates a new object from a selection
    """
    return _call("create", name=name, selection=selection, source_state=source_state)


@mcp.tool()
def extract(
    name: Annotated[str, Field(description="Name for the new object.")],
    selection: Annotated[
        str | None,
        Field(
            description="Atoms to move into it. Unlike create, these are *removed* from the original object."
        ),
    ] = "all",
) -> str:
    """
    Extracts a selection to a new object
    """
    return _call("extract", name=name, selection=selection)


@mcp.tool()
def delete(
    name: Annotated[
        str,
        Field(
            description="Object or named selection to delete. Accepts wildcards, e.g. ``tmp*``; ``all`` clears everything."
        ),
    ],
) -> str:
    """
    Deletes objects or selections
    """
    return _call("delete", name=name)


@mcp.tool()
def align(
    mobile: Annotated[str, Field(description="Selection to move.")],
    target: Annotated[
        str | None, Field(description="Selection to align onto; this one stays put.")
    ] = "all",
    options: Annotated[
        str | None,
        Field(description="Extra arguments, e.g. ``cycles=0`` to disable outlier rejection."),
    ] = None,
) -> str:
    """
    Aligns one selection to another
    """
    return _call("align", mobile=mobile, target=target, options=options)


@mcp.tool(name="super")
def super_tool(
    mobile: Annotated[str, Field(description="Selection to move.")],
    target: Annotated[
        str | None, Field(description="Selection to superimpose onto; this one stays put.")
    ] = "all",
    options: Annotated[str | None, Field(description="Extra PyMOL arguments.")] = None,
) -> str:
    """
    Superimposes one selection onto another
    """
    return _call("super", mobile=mobile, target=target, options=options)


@mcp.tool()
def intra_fit(
    selection: Annotated[
        str, Field(description="Object whose states to fit onto its first state.")
    ],
) -> str:
    """
    Fits all states within an object
    """
    return _call("intra_fit", selection=selection)


@mcp.tool()
def intra_rms(
    selection: Annotated[str, Field(description="Object whose states to compare.")],
) -> str:
    """
    Calculates RMSD between states within an object
    """
    return _call("intra_rms", selection=selection)


@mcp.tool()
def alter(
    selection: Annotated[str, Field(description="Atoms to modify.")],
    expression: Annotated[
        str,
        Field(
            description="Assignment over atom properties, e.g. ``b=0``, ``chain='B'``, ``resi=str(int(resi)+100)``."
        ),
    ],
) -> str:
    """
    Alters atomic properties in a selection
    """
    return _call("alter", selection=selection, expression=expression)


@mcp.tool()
def alter_state(
    state: Annotated[str, Field(description="State to modify.")],
    selection: Annotated[str, Field(description="Atoms to modify.")],
    expression: Annotated[str, Field(description="Assignment over coordinates, e.g. ``x=x+10``.")],
) -> str:
    """
    Alters atomic coordinates in a state
    """
    return _call("alter_state", state=state, selection=selection, expression=expression)


@mcp.tool()
def h_add(
    selection: Annotated[str | None, Field(description="Where to add hydrogens.")] = "all",
) -> str:
    """
    Adds hydrogens to a selection
    """
    return _call("h_add", selection=selection)


@mcp.tool()
def h_fill(
    selection: Annotated[str | None, Field(description="Where to fill valences.")] = "all",
) -> str:
    """
    Adds hydrogens and adjusts valences
    """
    return _call("h_fill", selection=selection)


@mcp.tool()
def bond(
    atom1: Annotated[
        str, Field(description="First atom, as a selection matching exactly one atom.")
    ],
    atom2: Annotated[
        str, Field(description="Second atom, as a selection matching exactly one atom.")
    ],
    order: Annotated[
        str | None,
        Field(description="Bond order: ``1`` single, ``2`` double, ``3`` triple, ``4`` aromatic."),
    ] = "1",
) -> str:
    """
    Creates a bond between two atoms
    """
    return _call("bond", atom1=atom1, atom2=atom2, order=order)


@mcp.tool()
def unbond(
    atom1: Annotated[str, Field(description="First atom of the bond to remove.")],
    atom2: Annotated[str, Field(description="Second atom of the bond to remove.")],
) -> str:
    """
    Removes a bond between two atoms
    """
    return _call("unbond", atom1=atom1, atom2=atom2)


@mcp.tool()
def rebuild(
    selection: Annotated[
        str | None,
        Field(
            description="What to regenerate. Needed after altering coordinates or B-factors for the change to show."
        ),
    ] = "all",
) -> str:
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
def util_cbc(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by chain (Color By Chain)
    """
    return _call("util.cbc", selection=selection)


@mcp.tool()
def util_cbaw(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, white carbons (Color By Atom, White)
    """
    return _call("util.cbaw", selection=selection)


@mcp.tool()
def util_cbag(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, green carbons (Color By Atom, Green)
    """
    return _call("util.cbag", selection=selection)


@mcp.tool()
def util_cbac(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, cyan carbons (Color By Atom, Cyan)
    """
    return _call("util.cbac", selection=selection)


@mcp.tool()
def util_cbam(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, magenta carbons (Color By Atom, Magenta)
    """
    return _call("util.cbam", selection=selection)


@mcp.tool()
def util_cbay(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, yellow carbons (Color By Atom, Yellow)
    """
    return _call("util.cbay", selection=selection)


@mcp.tool()
def util_cbas(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, salmon carbons (Color By Atom, Salmon)
    """
    return _call("util.cbas", selection=selection)


@mcp.tool()
def util_cbab(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, slate carbons (Color By Atom, slateBLue)
    """
    return _call("util.cbab", selection=selection)


@mcp.tool()
def util_cbao(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, orange carbons (Color By Atom, Orange)
    """
    return _call("util.cbao", selection=selection)


@mcp.tool()
def util_cbap(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, purple carbons (Color By Atom, Purple)
    """
    return _call("util.cbap", selection=selection)


@mcp.tool()
def util_cbak(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by atom, pink carbons (Color By Atom, pinK)
    """
    return _call("util.cbak", selection=selection)


@mcp.tool()
def util_chainbow(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors chains in rainbow gradient (CHAINs in rainBOW)
    """
    return _call("util.chainbow", selection=selection)


@mcp.tool()
def util_rainbow(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors residues in rainbow from N to C terminus
    """
    return _call("util.rainbow", selection=selection)


@mcp.tool()
def util_ss(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors by secondary structure
    """
    return _call("util.ss", selection=selection)


@mcp.tool()
def util_color_by_element(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors atoms by their element
    """
    return _call("util.color_by_element", selection=selection)


@mcp.tool()
def util_color_secondary(
    selection: Annotated[
        str | None, Field(description="What to recolour. Defaults to everything.")
    ] = "all",
) -> str:
    """
    Colors secondary structure elements
    """
    return _call("util.color_secondary", selection=selection)


@mcp.tool()
def spheroid(
    selection: Annotated[
        str | None, Field(description="What to render as smoothed spheres.")
    ] = "all",
) -> str:
    """
    Displays atoms as smooth spheres
    """
    return _call("spheroid", selection=selection)


@mcp.tool()
def isomesh(
    name: Annotated[str, Field(description="Name for the mesh object.")],
    map_object: Annotated[
        str, Field(description="Name of a loaded map (e.g. a CCP4 or DX density map).")
    ],
    level: Annotated[
        str, Field(description="Contour level in sigma, e.g. ``1.0`` for 2Fo-Fc density.")
    ],
    selection: Annotated[
        str | None, Field(description="Restrict the mesh to the region around this selection.")
    ] = "all",
) -> str:
    """
    Creates a mesh isosurface
    """
    return _call("isomesh", name=name, map_object=map_object, level=level, selection=selection)


@mcp.tool()
def isosurface(
    name: Annotated[str, Field(description="Name for the surface object.")],
    map_object: Annotated[str, Field(description="Name of a loaded map.")],
    level: Annotated[str, Field(description="Contour level in sigma.")],
    selection: Annotated[
        str | None, Field(description="Restrict the surface to the region around this selection.")
    ] = "all",
) -> str:
    """
    Creates a solid isosurface
    """
    return _call("isosurface", name=name, map_object=map_object, level=level, selection=selection)


@mcp.tool()
def sculpt_activate(
    obj: Annotated[str, Field(description="Object to enable real-space sculpting on.")],
) -> str:
    """
    Activates sculpting mode for an object
    """
    return _call("sculpt_activate", obj=obj)


@mcp.tool()
def sculpt_deactivate(obj: Annotated[str, Field(description="Object to stop sculpting.")]) -> str:
    """
    Deactivates sculpting mode for an object
    """
    return _call("sculpt_deactivate", obj=obj)


@mcp.tool()
def sculpt_iterate(
    iterations: Annotated[str, Field(description="Number of relaxation cycles to run.")],
    obj: Annotated[str | None, Field(description="Object to relax.")] = "all",
) -> str:
    """
    Performs sculpting iterations
    """
    return _call("sculpt_iterate", iterations=iterations, obj=obj)


@mcp.tool()
def scene(
    key: Annotated[
        str, Field(description='Scene name, e.g. "F1", or "auto" to use the next free slot.')
    ],
    action: Annotated[
        str | None,
        Field(
            description="``store`` to save the current view, ``recall`` to restore it, ``clear``, ``update``, ``rename``, ``delete``, ``next``, ``previous``."
        ),
    ] = "recall",
) -> str:
    """
    Manages scenes for later recall
    """
    return _call("scene", key=key, action=action)


@mcp.tool()
def scene_order(
    scene_list: Annotated[
        str, Field(description='Space-separated scene names in the order wanted, e.g. "F2 F1 F3".')
    ],
) -> str:
    """
    Sets the order of scenes
    """
    return _call("scene_order", scene_list=scene_list)


@mcp.tool()
def mset(
    specification: Annotated[
        str,
        Field(
            description='Frame-to-state mapping, e.g. "1 x30" to hold state 1 for 30 frames, or "1 -30" to sweep states 1 through 30.'
        ),
    ],
) -> str:
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
def frame(
    frame_number: Annotated[
        str | None, Field(description="Frame to jump to. Omit to query the current frame.")
    ] = None,
) -> str:
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
def ray(
    width: Annotated[
        str | None, Field(description="Image width in pixels. Omit to use the viewport size.")
    ] = None,
    height: Annotated[str | None, Field(description="Image height in pixels.")] = None,
) -> str:
    """
    Performs ray-tracing
    """
    return _call("ray", _timeout=_SLOW_OP_TIMEOUT, width=width, height=height)


@mcp.tool()
def draw(
    width: Annotated[str | None, Field(description="Image width in pixels.")] = None,
    height: Annotated[str | None, Field(description="Image height in pixels.")] = None,
) -> str:
    """
    Uses OpenGL renderer (faster but lower quality)
    """
    return _call("draw", _timeout=_SLOW_OP_TIMEOUT, width=width, height=height)


@mcp.tool()
def mpng(
    prefix: Annotated[
        str, Field(description="Filename prefix; PyMOL appends a zero-padded frame number.")
    ],
) -> str:
    """
    Saves a series of PNG images for movie frames
    """
    return _call("mpng", _timeout=_SLOW_OP_TIMEOUT, prefix=prefix)


@mcp.tool()
def symexp(
    prefix: Annotated[str, Field(description="Prefix for the generated symmetry-mate objects.")],
    selection: Annotated[
        str, Field(description="Selection whose crystallographic neighbours to build.")
    ],
    cutoff: Annotated[
        str | None, Field(description="Distance in Angstrom out to which to generate mates.")
    ] = "20",
    segi: Annotated[
        str | None, Field(description="Optional segment identifier for the new objects.")
    ] = None,
) -> str:
    """
    Generates symmetry-related copies
    """
    return _call("symexp", prefix=prefix, selection=selection, cutoff=cutoff, segi=segi)


@mcp.tool()
def set_symmetry(
    selection: Annotated[str, Field(description="Object to assign the unit cell to.")],
    a: Annotated[str, Field(description="Unit cell edge a, in Angstrom.")],
    b: Annotated[str, Field(description="Unit cell edge b, in Angstrom.")],
    c: Annotated[str, Field(description="Unit cell edge c, in Angstrom.")],
    alpha: Annotated[str, Field(description="Unit cell angle alpha, in degrees.")],
    beta: Annotated[str, Field(description="Unit cell angle beta, in degrees.")],
    gamma: Annotated[str, Field(description="Unit cell angle gamma, in degrees.")],
) -> str:
    """
    Sets symmetry parameters for an object
    """
    return _call(
        "set_symmetry", selection=selection, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma
    )


@mcp.tool()
def fab(
    sequence: Annotated[str, Field(description='One-letter amino acid sequence, e.g. "ACDEFGH".')],
    options: Annotated[
        str | None, Field(description="Extra arguments, e.g. ``ss=1`` to build it as a helix.")
    ] = None,
) -> str:
    """
    Creates a peptide chain from a sequence
    """
    return _call("fab", sequence=sequence, options=options)


@mcp.tool()
def fragment(
    name: Annotated[
        str, Field(description='Built-in fragment name, e.g. "benzene", "ala", "formamide".')
    ],
) -> str:
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
def viewport(
    width: Annotated[str, Field(description="Viewport width in pixels.")],
    height: Annotated[str, Field(description="Viewport height in pixels.")],
) -> str:
    """
    Sets the viewport size
    """
    return _call("viewport", width=width, height=height)


@mcp.tool()
def cd(
    path: Annotated[
        str, Field(description="Directory to change into. PyMOL resolves relative paths from here.")
    ],
) -> str:
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
def ls(
    path: Annotated[
        str | None, Field(description="Directory or glob to list. Omit for the current directory.")
    ] = None,
) -> str:
    """
    Lists files in the current directory
    """
    return _call("ls", path=path)


@mcp.tool()
def system(
    command: Annotated[
        str, Field(description="Shell command to run on the machine PyMOL is running on.")
    ],
) -> str:
    """
    Executes a system command
    """
    return _call("system", command=command)


@mcp.tool()
def help(
    command: Annotated[
        str | None, Field(description="PyMOL command to describe. Omit for general help.")
    ] = None,
) -> str:
    """
    Shows help for a command
    """
    return _call("help", command=command)
