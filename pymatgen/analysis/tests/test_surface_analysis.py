import unittest
import os
import random
import glob

import numpy as np

from pymatgen import SETTINGS
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.surface_analysis import SurfaceEnergyAnalyzer
from pymatgen.util.testing import PymatgenTest

__author__ = "Richard Tran"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "Aug 24, 2017"


def get_path(path_str):
    cwd = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(cwd, "..", "..", "..", "test_files", "surface_tests",
                        path_str)
    return path


class SurfaceEnergyAnalyzerTest(PymatgenTest):

    def setUp(self):



        vasprun_dict = {}
        for v in glob.glob(os.path.join(get_path(""), "*")):
            if ".xml.Cu" in v:
                vasprun_dict[tuple([int(i) for i in v[-6:].strip(".gz")])] = [Vasprun(v)]
        self.vasprun_dict = vasprun_dict
        self.Cu_analyzer = SurfaceEnergyAnalyzer("mp-30", self.vasprun_dict, "Cu")

    def test_gamma_calculator(self):

        # make sure we've loaded all our files correctly
        self.assertEqual(len(self.vasprun_dict.keys()), 13)
        for hkl, vaspruns in self.vasprun_dict.items():
            se_range = self.Cu_analyzer.calculate_gamma(vaspruns[0])
            # For a stoichiometric system, we expect surface
            # energy to be independent of chemical potential
            self.assertEqual(se_range[0], se_range[1])
            print(hkl, se_range[0])

    def test_get_intersections(self):

        # Just test if its working, for Cu, there are no different
        #  terminations so everything should be a nonetype
        for hkl in self.vasprun_dict.keys():
            self.assertFalse(self.Cu_analyzer.get_intersections(hkl))



if __name__ == "__main__":
    unittest.main()
