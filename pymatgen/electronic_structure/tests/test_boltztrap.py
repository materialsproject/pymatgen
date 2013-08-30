import unittest
import json
import os
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BoltztrapAnalyzer_test(unittest.TestCase):

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
        refp = [[9.67252848e-01, -7.41621352e-20, 1.81431663e-19],
                [-7.41621352e-20, 2.93967184e+00, -5.22926942e-18],
                [1.81431663e-19, -5.22926942e-18, 7.43551739e-01]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(arrayp[i][j], refp[i][j])
        arrayn = self.bz.get_average_eff_mass_tensor()['n']
        refn = [[1.01597485e+00, 5.26593684e-18, -5.07526263e-18],
                [5.26593684e-18, 1.00860922e+00, 4.16736724e-18],
                [-5.07526263e-18, 4.16736724e-18, 2.62490451e+00]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(arrayn[i][j], refn[i][j])

    def test_get_eig_average_eff_mass_tensor(self):
        list_p = self.bz.get_eig_average_eff_mass_tensor()['p']
        ref_p = [0.74355173938236407, 0.96725284778659004, 2.9396718381950864]
        for i in range(0, 3):
            self.assertAlmostEqual(list_p[i], ref_p[i])
        list_n = self.bz.get_eig_average_eff_mass_tensor()['n']
        ref_n = [1.0086092195986951, 1.0159748504188562, 2.6249045124335639]
        for i in range(0, 3):
            self.assertAlmostEqual(list_n[i], ref_n[i])


if __name__ == '__main__':
    unittest.main()
