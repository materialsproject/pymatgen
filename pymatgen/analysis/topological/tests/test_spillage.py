from pymatgen.util.testing import PymatgenTest
import unittest
import os
import warnings
from pymatgen.analysis.topological.spillage import overlap_so_spinpol


class SolarTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_slme_from_vasprun(self):
        wf_noso = os.path.join(os.path.dirname(__file__), "WAVECAR-NonSOC")
        wf_so = os.path.join(os.path.dirname(__file__), "WAVECAR-SOC")
        gamma_max = overlap_so_spinpol(wf_noso, wf_so)
        # print('spillage = ',gamma_max)
        self.assertEqual(gamma_max, 1.3634111271008749)


if __name__ == "__main__":
    unittest.main()
