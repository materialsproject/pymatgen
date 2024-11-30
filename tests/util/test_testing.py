from __future__ import annotations

from pymatgen.util.testing import PymatgenTest


class TestPymatgenTestTestCase(PymatgenTest):
    """Baseline inspector for migration side effects.

    This is not a functional test but a utility to verify behaviors specific to
    unittest.TestCase. It ensures we're aware the side effects from the migration.
    """
