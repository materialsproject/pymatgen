"""This is not a functional test but a utility to verify behaviors specific to
`unittest.TestCase`. It ensures we're aware the side effects from the migration.

TODO: remove this test module after migration (2026-01-01), see PR 4209.

`unittest.TestCase`-specific features and brief migration guide:
- Setup/teardown methods (`setUp`, `setUpClass`, `tearDown`, `tearDownClass`):
    1. Recommended approach in pytest: Use fixtures.
       Documentation: https://docs.pytest.org/en/stable/reference/fixtures.html#fixture
    OR
    2. Use pytest's xUnit-style setup/teardown functions:
       `[setup/teardown]_[class/method/function]`.
       Documentation: https://docs.pytest.org/en/stable/how-to/xunit_setup.html

- Assertion methods (`assertTrue`, `assertFalse`, `assertEqual`, etc.):
    Replace with direct Python `assert` statements.
"""

# ruff: noqa: PT009, PT027, FBT003

from __future__ import annotations

import pytest

from pymatgen.util.testing import PymatgenTest


class TestPymatgenTestTestCase(PymatgenTest):
    """Baseline inspector for migration side effects."""

    def test_unittest_testcase_specific_funcs(self):
        """Make sure TestCase-specific methods are available until migration.
        TODO: check warnings
        """

        # Testing setUp and tearDown methods
        self.setUp()
        self.tearDown()

        # Testing class-level setUp and tearDown methods
        self.setUpClass()
        self.tearDownClass()

        # Test the assertion methods
        self.assertTrue(True)
        self.assertFalse(False)
        self.assertEqual(1, 1)


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
