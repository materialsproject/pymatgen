from __future__ import annotations

import os
import unittest
import warnings

from pytest import approx

from pymatgen.analysis.topological.spillage import SOCSpillage
from pymatgen.util.testing import PymatgenTest


class SolarTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_spillage_from_vasprun(self):
        wf_noso = os.path.join(os.path.dirname(__file__), "WAVECAR-NonSOC")
        wf_so = os.path.join(os.path.dirname(__file__), "WAVECAR-SOC")
        # JVASP-1044
        gamma_max = SOCSpillage(wf_noso=wf_noso, wf_so=wf_so).overlap_so_spinpol()

        assert gamma_max == approx(1.3634111271008775, abs=1e-5)


if __name__ == "__main__":
    unittest.main()
