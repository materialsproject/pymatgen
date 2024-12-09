"""Test utilities for JDFTx tests shared between inputs and outputs parts.

This module will inherit everything from inputs_test_utils.py and outputs_test_utils.py upon final implementation.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pymatgen.util.testing import TEST_FILES_DIR

dump_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "tmp"


def assert_same_value(testval, knownval):
    if type(testval) not in [tuple, list]:
        assert isinstance(testval, type(knownval))
        if isinstance(testval, float):
            assert testval == pytest.approx(knownval)
        elif isinstance(testval, dict):
            for k in knownval:
                assert k in testval
                assert_same_value(testval[k], knownval[k])
        elif testval is None:
            assert knownval is None
        else:
            assert testval == knownval
    else:
        assert len(testval) == len(knownval)
        for i in range(len(testval)):
            assert_same_value(testval[i], knownval[i])
