import unittest
import json
import os
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BoltztrapAnalyzer_test(unittest.TestCase):

    def setUp(self):
        self.bz = BoltztrapAnalyzer.from_files(os.path.join(test_dir, "boltztrap"))

    def test_properties(self):
        self.assertAlmostEqual(self.bz.gap, 1.66449389)
        array = self.bz.cond[300][102]
        self.assertAlmostEqual(array[0][0], 7.5756518e+19)
        self.assertAlmostEqual(array[0][2], -11.14679)
        self.assertAlmostEqual(array[1][0], -88.203286)
        self.assertAlmostEqual(array[2][2], 1.7133249e+19)
        array = self.bz.seebeck[300][22]
        self.assertAlmostEqual(array[0][1], 6.4546074e-22)
        self.assertAlmostEqual(array[1][1], -0.00032073711)
        self.assertAlmostEqual(array[1][2], -2.9868424e-24)
        self.assertAlmostEqual(array[2][2], -0.0003126543)
        array = self.bz.kappa[500][300]
        self.assertAlmostEqual(array[0][1], 0.00014524309)
        self.assertAlmostEqual(array[1][1], 328834400000000.0)
        self.assertAlmostEqual(array[1][2], 3.7758069e-05)
        self.assertAlmostEqual(array[2][2], 193943750000000.0)
        self.assertAlmostEqual(self.bz.hall[400][800][3], 9.5623749e-28)
        self.assertAlmostEqual(self.bz.hall[400][68][7], 6.5106975e-10)
        self.assertAlmostEqual(self.bz.doping['p'][3], 1e18)
        self.assertAlmostEqual(self.bz.mu_doping['p'][300][2], 0.155377071914)
        self.assertAlmostEqual(self.bz.mu_doping['n'][300][-1], 1.64860243466)
        self.assertAlmostEqual(self.bz.cond_doping['n'][800][3][1][1], 1.5564085e+16)
        self.assertAlmostEqual(self.bz.seebeck_doping['p'][600][2][0][1], 3.2860613e-23)
        self.assertAlmostEqual(self.bz.carrier_conc[500][67], 38.22832002)
        self.assertAlmostEqual(self.bz.vol, 612.975564135)

    def test_get_average_eff_mass_tensor(self):
        arrayp = self.bz.get_average_eff_mass_tensor()['p']
        refp = [[9.59811382e-01, -8.43977792e-19, -4.66376242e-19],
                [-8.43977792e-19, 2.93705659e+00, 2.98230675e-18],
                [-4.66376242e-19, 2.98230675e-18, 7.60149480e-01]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(arrayp[i][j], refp[i][j])
        arrayn = self.bz.get_average_eff_mass_tensor()['n']
        refn = [[9.93175685e-01, -1.19952865e-18, 2.09984689e-18],
                [-1.19952865e-18, 9.81375550e-01, -2.24252485e-18],
                [2.09984689e-18, -2.24252485e-18, 2.64678926e+00]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(arrayn[i][j], refn[i][j])

    def test_get_eig_average_eff_mass_tensor(self):
        list_p = self.bz.get_eig_average_eff_mass_tensor()['p']
        ref_p = [0.76014948015275929, 0.95981138218632356, 2.9370565873468486]
        for i in range(0, 3):
            self.assertAlmostEqual(list_p[i], ref_p[i])
        list_n = self.bz.get_eig_average_eff_mass_tensor()['n']
        ref_n = [0.98137554988503861, 0.9931756850257033, 2.6467892600243781]
        for i in range(0, 3):
            self.assertAlmostEqual(list_n[i], ref_n[i])


if __name__ == '__main__':
    unittest.main()
