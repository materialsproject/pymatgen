# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import logging
import os
import struct
import unittest

import pytest
from monty.io import zopen

from pymatgen.io.qchem.utils import lower_and_check_unique, process_parsed_HESS
from pymatgen.util.testing import PymatgenTest

__author__ = "Ryan Kingsbury, Samuel Blau"
__copyright__ = "Copyright 2018-2022, The Materials Project"

logger = logging.getLogger(__name__)


test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "new_qchem_files")


class UtilTest(PymatgenTest):
    """
    test utils
    """

    def test_lower_and_check_unique(self):
        d = {"sVp": {"RHOISO": 0.0009}, "jobType": "SP"}
        d2 = lower_and_check_unique(d)
        assert d2 == {"svp": {"RHOISO": 0.0009}, "job_type": "sp"}
        d3 = lower_and_check_unique(d2["svp"])
        assert d3 == {"rhoiso": "0.0009"}
        d4 = {"jobType": "SP", "JOBtype": "SP"}
        # should not raise an exception
        assert lower_and_check_unique(d4) == {"job_type": "sp"}
        d4.update({"jobType": "opt"})
        with pytest.raises(ValueError, match="Multiple instances of key"):
            lower_and_check_unique(d4)

    def test_process_parsed_HESS(self):
        data_132 = []
        with zopen(os.path.join(test_dir, "parse_hess", "132.0"), mode="rb") as file:
            binary = file.read()
            for ii in range(int(len(binary) / 8)):
                data_132.append(struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0])

        data_HESS = []
        with zopen(
            os.path.join(test_dir, "parse_hess", "HESS"),
            mode="rt",
            encoding="ISO-8859-1",
        ) as f:
            data_HESS = f.readlines()

        processed_data_HESS = process_parsed_HESS(data_HESS)

        assert len(data_132) == len(processed_data_HESS)
        for ii, val in enumerate(data_132):
            diff = abs(val - processed_data_HESS[ii])
            assert diff < 1e-15


if __name__ == "__main__":
    unittest.main()
