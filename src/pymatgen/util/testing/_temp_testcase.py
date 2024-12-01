"""Temporary TestCase for migration to `pytest` framework,
inserted FutureWarning for unittest.TestCase-specific methods.

TODO: remove entire module after migration
"""

# ruff: noqa: PT009, PT027

from __future__ import annotations

import warnings
from unittest import TestCase


class _TempTestCase4Migrate(TestCase):
    @staticmethod
    def _issue_warning(method_name):
        warnings.warn(
            f"unittest {method_name=} will not be supported by pytest after migration by 2026-01-01, see PR4209.",
            FutureWarning,
            stacklevel=2,
        )

    def setUp(self, *args, **kwargs):
        self._issue_warning("setUp")
        super().setUp(*args, **kwargs)

    def tearDown(self, *args, **kwargs):
        self._issue_warning("tearDown")
        super().tearDown(*args, **kwargs)

    @classmethod
    def setUpClass(cls, *args, **kwargs):
        cls._issue_warning("setUpClass")
        super().setUpClass(*args, **kwargs)

    @classmethod
    def tearDownClass(cls, *args, **kwargs):
        cls._issue_warning("tearDownClass")
        super().tearDownClass(*args, **kwargs)

    def assertEqual(self, *args, **kwargs):
        self._issue_warning("assertEqual")
        return super().assertEqual(*args, **kwargs)

    def assertNotEqual(self, *args, **kwargs):
        self._issue_warning("assertNotEqual")
        return super().assertNotEqual(*args, **kwargs)

    def assertTrue(self, *args, **kwargs):
        self._issue_warning("assertTrue")
        return super().assertTrue(*args, **kwargs)

    def assertFalse(self, *args, **kwargs):
        self._issue_warning("assertFalse")
        return super().assertFalse(*args, **kwargs)

    def assertIsNone(self, *args, **kwargs):
        self._issue_warning("assertIsNone")
        return super().assertIsNone(*args, **kwargs)

    def assertIsNotNone(self, *args, **kwargs):
        self._issue_warning("assertIsNotNone")
        return super().assertIsNotNone(*args, **kwargs)

    def assertIn(self, *args, **kwargs):  # codespell:ignore
        self._issue_warning("assertIn")  # codespell:ignore
        return super().assertIn(*args, **kwargs)  # codespell:ignore

    def assertNotIn(self, *args, **kwargs):
        self._issue_warning("assertNotIn")
        return super().assertNotIn(*args, **kwargs)

    def assertIsInstance(self, *args, **kwargs):
        self._issue_warning("assertIsInstance")
        return super().assertIsInstance(*args, **kwargs)

    def assertNotIsInstance(self, *args, **kwargs):
        self._issue_warning("assertNotIsInstance")
        return super().assertNotIsInstance(*args, **kwargs)

    def assertRaises(self, *args, **kwargs):
        self._issue_warning("assertRaises")
        return super().assertRaises(*args, **kwargs)
