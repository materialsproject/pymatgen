from __future__ import annotations

import logging
import struct

import pytest
from monty.io import zopen

from pymatgen.io.qchem.utils import lower_and_check_unique, process_parsed_hess
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Ryan Kingsbury, Samuel Blau"
__copyright__ = "Copyright 2018-2022, The Materials Project"

logger = logging.getLogger(__name__)


TEST_DIR = f"{TEST_FILES_DIR}/io/qchem/new_qchem_files"


class TestUtil(PymatgenTest):
    """test utils."""

    def test_lower_and_check_unique(self):
        dct = {"sVp": {"RHOISO": 0.0009}, "jobType": "SP"}
        d2 = lower_and_check_unique(dct)
        assert d2 == {"svp": {"RHOISO": 0.0009}, "job_type": "sp"}
        d3 = lower_and_check_unique(d2["svp"])
        assert d3 == {"rhoiso": "0.0009"}
        d4 = {"jobType": "SP", "JOBtype": "SP"}
        # should not raise an exception
        assert lower_and_check_unique(d4) == {"job_type": "sp"}
        d4["jobType"] = "opt"
        with pytest.raises(ValueError, match="Multiple instances of key"):
            lower_and_check_unique(d4)

    def test_process_parsed_hess(self):
        with zopen(f"{TEST_DIR}/parse_hess/132.0", mode="rb") as file:
            binary = file.read()
            data_132 = [struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0] for ii in range(int(len(binary) / 8))]

        with zopen(f"{TEST_DIR}/parse_hess/HESS", mode="rt", encoding="ISO-8859-1") as file:
            data_hess = file.readlines()

        processed_data_hess = process_parsed_hess(data_hess)

        assert len(data_132) == len(processed_data_hess)
        for ii, val in enumerate(data_132):
            diff = abs(val - processed_data_hess[ii])
            assert diff < 1e-15
