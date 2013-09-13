import unittest
import json
import os
import numpy.testing as nptest
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BoltztrapAnalyzerTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "mp-22630_bz.json"), "rb") as f:
            d = json.loads(f.read())
            self.bz = BoltztrapAnalyzer.from_dict(d)

    def test_properties(self):
        self.assertAlmostEqual(self.bz.gap, 1.66449389)
        array = self.bz.cond[300][102]
        self.assertAlmostEqual(array[0][0], 7.7196281e+19)
        self.assertAlmostEqual(array[0][2], 8.3152525)
        self.assertAlmostEqual(array[1][0], -30.463203)
        self.assertAlmostEqual(array[2][2], 1.7804694e+19)
        array = self.bz.seebeck[300][22]
        self.assertAlmostEqual(array[0][1], -4.9346025e-22)
        self.assertAlmostEqual(array[1][1], -0.00031754687)
        self.assertAlmostEqual(array[1][2], -1.4829561e-23)
        self.assertAlmostEqual(array[2][2], -0.00030791881)
        array = self.bz.kappa[500][300]
        self.assertAlmostEqual(array[0][1], -0.0019343186)
        self.assertAlmostEqual(array[1][1], 329789550000000.0)
        self.assertAlmostEqual(array[1][2], -3.6551152e-05)
        self.assertAlmostEqual(array[2][2], 199006400000000.0)
        self.assertAlmostEqual(self.bz.hall[400][800][3], 2.9645623e-27)
        self.assertAlmostEqual(self.bz.hall[400][68][7], 5.5418434e-10)
        self.assertAlmostEqual(self.bz.doping['p'][3], 1e18)
        self.assertAlmostEqual(self.bz.mu_doping['p'][300][2], 0.155921299836)
        self.assertAlmostEqual(self.bz.mu_doping['n'][300][-1], 1.64765003579)
        self.assertAlmostEqual(self.bz.cond_doping['n'][800][3][1][1], 1.5255797e+16)
        self.assertAlmostEqual(self.bz.seebeck_doping['p'][600][2][0][1], -3.9286981e-22)

    def test_get_average_eff_mass_tensor(self):
        arrayp = self.bz.get_average_eff_mass_tensor()['p']
        refp = [[9.68887240e-01, -7.42874488e-20, 1.81738232e-19],
                [-7.42874488e-20, 2.94463908, -5.23810545e-18],
                [1.81738232e-19, -5.23810545e-18, 7.44808138e-01]]
        nptest.assert_allclose(arrayp, refp)
        arrayn = self.bz.get_average_eff_mass_tensor()['n']
        refn = [[1.01769157e+00, 5.27483483e-18, -5.08383843e-18],
                [5.27483483e-18, 1.01031349e+00, 4.17440894e-18],
                [-5.08383843e-18, 4.17440894e-18, 2.62933988e+00]]
        nptest.assert_allclose(arrayn, refn)

    def test_get_eig_average_eff_mass_tensor(self):
        list_p = self.bz.get_eig_average_eff_mass_tensor()['p']
        ref_p = [0.74480813774834398, 0.96888723962925827, 2.9446390767857635]
        for i in range(0, 3):
            self.assertAlmostEqual(list_p[i], ref_p[i])
        list_n = self.bz.get_eig_average_eff_mass_tensor()['n']
        ref_n = [1.0103134923591464, 1.0176915690738393, 2.6293398806340513]
        for i in range(0, 3):
            self.assertAlmostEqual(list_n[i], ref_n[i])


if __name__ == '__main__':
    unittest.main()
