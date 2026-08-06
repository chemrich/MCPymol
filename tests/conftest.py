"""Shared test setup.

``mcpymol.plugin`` does ``from pymol import cmd`` at import time and binds the
result as a module global, so the mock has to be installed before *any* test
module imports it — and it has to be the *same* mock for every test file, or a
test asserting on its own mock silently watches an object the plugin never
uses.  conftest is imported before the test modules, which makes this the one
place that can guarantee both.
"""

import sys
from unittest.mock import MagicMock

mock_pymol = MagicMock()
sys.modules.setdefault("pymol", mock_pymol)
sys.modules.setdefault("pymol.cmd", mock_pymol.cmd)
