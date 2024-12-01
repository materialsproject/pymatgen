"""This is not a functional test but a utility to verify behaviors specific to
`unittest.TestCase`. It ensures we're aware the side effects from the migration.

TODO: remove this test module after migration (2026-01-01), see PR 4209.
"""

# ruff: noqa: PT009, PT027, FBT003

from __future__ import annotations

import pytest

from pymatgen.util.testing import PymatgenTest


@pytest.mark.filterwarnings("ignore", message="will not be supported by pytest after migration", category=FutureWarning)
class TestPymatgenTestTestCase(PymatgenTest):
    """Baseline inspector for migration side effects."""

    def test_unittest_testcase_specific_funcs(self):
        """Make sure TestCase-specific methods are available until migration,
        and FutureWarning is emitted.
        """
        msg = "will not be supported by pytest after migration"

        # Testing setUp and tearDown methods
        with pytest.warns(FutureWarning, match=msg):
            self.setUp()
        with pytest.warns(FutureWarning, match=msg):
            self.tearDown()

        # Testing class-level setUp and tearDown methods
        with pytest.warns(FutureWarning, match=msg):
            self.setUpClass()
        with pytest.warns(FutureWarning, match=msg):
            self.tearDownClass()

        # Test the assertion methods
        with pytest.warns(FutureWarning, match=msg):
            self.assertEqual(1, 1)

        with pytest.warns(FutureWarning, match=msg):
            self.assertNotEqual(1, 2)

        with pytest.warns(FutureWarning, match=msg):
            self.assertTrue(True)

        with pytest.warns(FutureWarning, match=msg):
            self.assertFalse(False)

        with pytest.warns(FutureWarning, match=msg):
            self.assertIsNone(None)

        with pytest.warns(FutureWarning, match=msg):
            self.assertIsNotNone("hello")

        with pytest.warns(FutureWarning, match=msg), self.assertRaises(ValueError):
            raise ValueError("hi")


class TestPymatgenTestPytest:
    def test_unittest_testcase_specific_funcs(self):
        """Test unittest.TestCase-specific methods for migration to pytest."""
        # Testing setUp and tearDown methods
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.setUp()
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.tearDown()

        # Testing class-level setUp and tearDown methods
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.setUpClass()
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.tearDownClass()

        # Test the assertion methods
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.assertTrue(True)
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.assertFalse(False)
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.assertEqual(1, 1)
