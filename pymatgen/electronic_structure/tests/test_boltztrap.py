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
        array = self.bz.cond[300.0][102]
        self.assertAlmostEqual(array[0][0], 5.7935315e+17)
        self.assertAlmostEqual(array[0][2], -1.465385)
        self.assertAlmostEqual(array[1][0], -1.6783699)
        self.assertAlmostEqual(array[2][2], 4.8133152e+17)
        array = self.bz.seebeck[300.0][22]
        self.assertAlmostEqual(array[0][1], -1.6064909e-21)
        self.assertAlmostEqual(array[1][1], -0.00060416451)
        self.assertAlmostEqual(array[1][2], -7.344964e-24)
        self.assertAlmostEqual(array[2][2], -0.0005929311)
        array = self.bz.kappa[500.0][3000]
        self.assertAlmostEqual(array[0][1], -4.4284407e-10)
        self.assertAlmostEqual(array[1][1], 317231920.0)
        self.assertAlmostEqual(array[1][2], -7.9160729e-11)
        self.assertAlmostEqual(array[2][2], 88990472.0)
        self.assertAlmostEqual(self.bz.hall[400.0][1502][3], -5.8145763e-28)
        self.assertAlmostEqual(self.bz.hall[400.0][1502][7], -5.2528111e-10)
        self.assertAlmostEqual(self.bz.doping['p'][3], 1e18)
        self.assertAlmostEqual(self.bz.mu_doping['p'][300.0][2], 0.15605735681702)
        self.assertAlmostEqual(self.bz.mu_doping['n'][300.0][-1], 1.64778609277326)
        self.assertAlmostEqual(self.bz.cond_doping['n'][800.0][3][1][1], 1.5259202e+16)
        self.assertAlmostEqual(self.bz.seebeck_doping['p'][600.0][2][0][1], -5.327734e-22)

    def test_get_average_eff_mass_tensor(self):
        arrayp = self.bz.get_average_eff_mass_tensor()['p']
        refp = [[9.66529786e-01, -1.51068648e-18, 2.64045712e-19],
                [-1.51068648e-18, 2.95277153e+00, -6.05073958e-18],
                [2.64045712e-19, -6.05073958e-18, 7.56064465e-01]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(arrayp[i][j], refp[i][j])
        arrayn = self.bz.get_average_eff_mass_tensor()['n']
        refn = [[1.00705046e+00, 6.36676319e-18, -4.75261854e-18],
                [6.36676319e-18, 9.99191860e-01, 3.46451475e-18],
                [-4.75261854e-18, 3.46451475e-18, 2.64871412e+00]]
        for i in range(0, 3):
            for j in range(0, 3):
                self.assertAlmostEqual(arrayn[i][j], refn[i][j])

    def test_get_eig_average_eff_mass_tensor(self):
        list_p = self.bz.get_eig_average_eff_mass_tensor()['p']
        ref_p = [0.75606446476510314, 0.96652978558578007, 2.9527715280817657]
        for i in range(0, 3):
            self.assertAlmostEqual(list_p[i], ref_p[i])
        list_n = self.bz.get_eig_average_eff_mass_tensor()['n']
        ref_n = [0.99919186025868112, 1.0070504649439096, 2.6487141155664782]
        for i in range(0, 3):
            self.assertAlmostEqual(list_n[i], ref_n[i])


if __name__ == '__main__':
    unittest.main()