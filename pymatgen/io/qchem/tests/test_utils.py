# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import logging

import pytest

from pymatgen.io.qchem.utils import lower_and_check_unique
from pymatgen.util.testing import PymatgenTest

__author__ = "Ryan Kingsbury"
__copyright__ = "Copyright 2018-2022, The Materials Project"

logger = logging.getLogger(__name__)


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
