from __future__ import annotations

import unittest

import numpy as np

from pymatgen.analysis.nmr import ChemicalShielding, ElectricFieldGradient
from pymatgen.util.testing import PymatgenTest


class TestChemicalShieldingNotation(PymatgenTest):
    def test_construction(self):
        cs = ChemicalShielding(np.arange(9).reshape((3, 3)))
        self.assertEqual(cs.shape, (3, 3))

        cs = ChemicalShielding([1, 2, 3])
        self.assertEqual(cs.shape, (3, 3))
        self.assertArrayEqual(np.diag(cs), [1, 2, 3])

    def test_principal_axis_system(self):
        cs = ChemicalShielding([1, 2, 3])
        self.assertArrayEqual(cs.principal_axis_system, cs)

        cs = ChemicalShielding(np.arange(9).reshape((3, 3)))
        self.assertArrayAlmostEqual(
            np.diag(cs.principal_axis_system),
            [-1.74596669e00, -1.53807726e-15, 1.37459667e01],
            decimal=5,
        )

    def test_notations(self):
        cs = ChemicalShielding.from_maryland_notation(195.0788, 68.1733, 0.8337)
        hae1 = cs.haeberlen_values
        self.assertAlmostEqual(hae1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(hae1.delta_sigma_iso, -65.33899505250002, places=5)
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


class TestElectricFieldGradient(PymatgenTest):
    def test_construction(self):
        efg = ElectricFieldGradient(np.arange(9).reshape((3, 3)))
        self.assertEqual(efg.shape, (3, 3))

        efg = ElectricFieldGradient([1, 2, 3])
        self.assertEqual(efg.shape, (3, 3))

    def test_principal_axis_system(self):
        efg = ElectricFieldGradient([1, 2, 3])
        self.assertArrayEqual(efg.principal_axis_system, efg)

        efg = ElectricFieldGradient(np.arange(9).reshape((3, 3)))
        self.assertArrayAlmostEqual(
            np.diag(efg.principal_axis_system),
            [-1.3484692e00, -1.1543332e-15, 1.3348469e01],
            decimal=5,
        )

    def test_Attributes(self):

        efg = ElectricFieldGradient([[11.11, 1.371, 2.652], [1.371, 3.635, -3.572], [2.652, -3.572, -14.746]])
        self.assertAlmostEqual(efg.V_yy, 11.516, places=3)
        self.assertAlmostEqual(efg.V_xx, 4.204, places=3)
        self.assertAlmostEqual(efg.V_zz, -15.721, places=3)
        self.assertAlmostEqual(efg.asymmetry, 0.465, places=3)
        self.assertAlmostEqual(efg.coupling_constant("Al"), 5.573, places=3)


if __name__ == "__main__":
    unittest.main()
