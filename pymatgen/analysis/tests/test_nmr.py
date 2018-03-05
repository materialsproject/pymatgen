import copy
import json

import unittest
from monty.json import MontyDecoder

import numpy as np
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.nmr import ChemicalShift


class TestChemicalShiftNotation(PymatgenTest):
    
    def test_construction(self):
        cs = ChemicalShift(np.arange(9).reshape((3,3)))
        self.assertEqual(cs.shape,(3,3))

        cs = ChemicalShift([1,2,3])
        self.assertEqual(cs.shape,(3,3))
        self.assertArrayEqual(np.diag(cs),[1,2,3])


    def test_principal_axis_system(self):
        cs = ChemicalShift([1,2,3])
        self.assertArrayEqual(cs.principal_axis_system,cs)

        cs = ChemicalShift(np.arange(9).reshape((3,3)))
        self.assertArrayAlmostEqual(np.diag(cs.principal_axis_system),[-1.3484692e+00, -1.1543332e-15,  1.3348469e+01],decimal=5)

    def test_notations(self):
        cs = ChemicalShift.from_maryland_notation(
            195.0788, 68.1733, 0.8337)
        hae1 = cs.haeberlen_values
        self.assertAlmostEqual(hae1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(hae1.delta_sigma, -65.33899505250002, places=5)
        self.assertAlmostEqual(hae1.zeta, -43.559330035000016, places=5)
        self.assertAlmostEqual(hae1.eta, 0.13013537835511454, places=5)
        meh1 = cs.mehring_values
        self.assertAlmostEqual(meh1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(meh1.sigma_11, 151.51946996499998, places=5)
        self.assertAlmostEqual(meh1.sigma_22, 214.02416007, places=5)
        self.assertAlmostEqual(meh1.sigma_33, 219.69276996500002, places=5)
        mary1 = cs.maryland_values
        self.assertAlmostEqual(mary1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(mary1.omega, 68.1733, places=5)
        self.assertAlmostEqual(mary1.kappa, 0.8337, places=5)
 

if __name__ == '__main__':
    unittest.main()