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

from __future__ import annotations

import pytest

from pymatgen.util.testing import PymatgenTest


@pytest.mark.filterwarnings("ignore:PymatgenTest is scheduled for migration to pytest")
class TestPymatgenTestTestCase(PymatgenTest):
    """Baseline inspector for migration side effects."""

    def test_pmg_test_migration_warning(self):
        """Test PymatgenTest migration warning."""
        with pytest.warns(FutureWarning, match="PymatgenTest is scheduled for migration to pytest after 2026-01-01"):
            self.setUpClass()  #  invoke the setup phase

    def test_unittest_testcase_specific_funcs(self):
        """Make sure TestCase-specific methods are available until migration."""
        # Testing setUp and tearDown methods
        self.setUp()
        self.tearDown()

        # Testing class-level setUp and tearDown methods
        self.setUpClass()
        self.tearDownClass()

        # Test the assertion methods
        self.assertTrue(True)  # noqa: PT009, FBT003
        self.assertFalse(False)  # noqa: PT009, FBT003
        self.assertEqual(1, 1)  # noqa: PT009


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
            self.assertTrue(True)  # noqa: PT009, FBT003
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.assertFalse(False)  # noqa: PT009, FBT003
        with pytest.raises(AttributeError, match="'TestPymatgenTestPytest' object has no attribute"):
            self.assertEqual(1, 1)  # noqa: PT009
