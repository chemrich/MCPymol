"""Wiggles: cryo-EM heterogeneity, occupancy and resolution views.

Ten tools for the things a structure viewer usually throws away — partial
occupancy, alternate conformations, ensemble spread, published Q-scores, map
geometry, and the per-voxel resolution field. Each one exists because a
plausible-looking render was hiding something a user needed to know.

The design rule the whole package follows: **say what is being shown, and
refuse when the picture would be meaningful-looking and wrong.** `morph_states`
declines to interpolate states that do not share a topology.
`local_resolution_view` declines to colour one map by another that does not
share its grid. `load_map` records provenance and never guesses it, because a
measured reconstruction and a neural-network hallucination are the same
isosurface once they are drawn.

Layout:

======================  ====================================================
``tools``               the ``@mcp.tool`` boundary — thin wrappers, no logic
``port``                the entire coupling surface to PyMOL
``atoms``, ``bfactors`` shared read layer and the B-factor stash
``occupancy``           ``occupancy_view``, ``altloc_view``
``ensembles``           ``ensemble_spread_view``, ``morph_states``
``qscore``              validation-report parsing and ``qscore_view``
``mapinfo``             MRC/CCP4 header reading; no PyMOL, no network
``maps``, ``density``   loading volumes, and contouring them honestly
``localres``            ``local_resolution_view``
``provenance``          where a volume came from; the default is UNKNOWN
======================  ====================================================

Everything above ``port`` is pure: the whole suite runs against ``FakePort``
with no PyMOL installed. Importing this package registers its tools, which is
why ``mcpymol.server`` imports it.

The research this was distilled from — a cited compendium of cryo-EM
heterogeneity methods — lives in a separate repository. Docstrings here name
the entry they came from (`occupancy-two-senses`, `local-resolution`) rather
than a path that would not resolve.
"""

from mcpymol.wiggles.atoms import Atom, fetch_atoms
from mcpymol.wiggles.bfactors import clear_stash
from mcpymol.wiggles.density import to_absolute, to_sigma
from mcpymol.wiggles.localres import grid_differences
from mcpymol.wiggles.mapinfo import MapHeader, read_map_header
from mcpymol.wiggles.port import BridgePort, FakePort, PortError, PymolPort, SendRequestPort
from mcpymol.wiggles.provenance import Provenance, declare, provenance_of

# Load-bearing: importing the tool module is what registers the ten tools on
# the shared FastMCP app. Imported last so the names above are available to it.
from mcpymol.wiggles.tools import (
    altloc_view,
    density_view,
    ensemble_spread_view,
    load_map,
    local_resolution_view,
    map_info,
    morph_states,
    occupancy_view,
    qscore_view,
    restore_bfactors,
)

__all__ = [  # noqa: RUF022
    # the ten tools, grouped by tier rather than alphabetically
    "occupancy_view",
    "altloc_view",
    "ensemble_spread_view",
    "morph_states",
    "qscore_view",
    "restore_bfactors",
    "map_info",
    "load_map",
    "density_view",
    "local_resolution_view",
    # the pieces worth reaching for directly
    "MapHeader",
    "read_map_header",
    "to_sigma",
    "to_absolute",
    "grid_differences",
    "Provenance",
    "declare",
    "provenance_of",
    "Atom",
    "fetch_atoms",
    "clear_stash",
    "PymolPort",
    "FakePort",
    "BridgePort",
    "SendRequestPort",
    "PortError",
]
