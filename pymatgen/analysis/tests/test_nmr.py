import copy
import json

import unittest
from monty.json import MontyDecoder

from pymatgen.analysis.nmr import NMRChemicalShiftNotation


class TestChemicalShiftNotation(unittest.TestCase):
    def test_nmr_chemical_shift_notation(self):
        cs1 = NMRChemicalShiftNotation.from_maryland_notation(
            195.0788, 68.1733, 0.8337)
        hae1 = cs1.haeberlen_values
        self.assertAlmostEqual(hae1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(hae1.delta_sigma, -65.33899505250002, places=5)
        self.assertAlmostEqual(hae1.zeta, -43.559330035000016, places=5)
        self.assertAlmostEqual(hae1.eta, 0.13013537835511454, places=5)
        meh1 = cs1.mehring_values
        self.assertAlmostEqual(meh1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(meh1.sigma_11, 151.51946996499998, places=5)
        self.assertAlmostEqual(meh1.sigma_22, 214.02416007, places=5)
        self.assertAlmostEqual(meh1.sigma_33, 219.69276996500002, places=5)
        mary1 = cs1.maryland_values
        self.assertAlmostEqual(mary1.sigma_iso, 195.0788, places=5)
        self.assertAlmostEqual(mary1.omega, 68.1733, places=5)
        self.assertAlmostEqual(mary1.kappa, 0.8337, places=5)
        cs2 = copy.deepcopy(cs1)
        cs2.sigma_22 = 235.0
        hae2 = cs2.haeberlen_values
        self.assertAlmostEqual(hae2.sigma_iso, 202.07074664333334, places=5)
        self.assertAlmostEqual(hae2.delta_sigma, -75.82691501750003, places=5)
        self.assertAlmostEqual(hae2.zeta, -50.55127667833335, places=5)
        self.assertAlmostEqual(hae2.eta, 0.30280600295028287, places=5)
        cs2.sigma_11 = 200.0
        cs2.sigma_33 = 270.0
        hae3 = cs2.haeberlen_values
        self.assertAlmostEqual(hae3.sigma_iso, 235.0, places=5)
        self.assertAlmostEqual(hae3.delta_sigma, -52.5, places=5)
        self.assertAlmostEqual(hae3.zeta, -35.0, places=5)
        self.assertAlmostEqual(hae3.eta, 1.0, places=5)
        cs2.sigma_11 = 235.0
        hae4 = cs2.haeberlen_values
        self.assertAlmostEqual(hae4.sigma_iso, 246.66666666666666, places=5)
        self.assertAlmostEqual(hae4.delta_sigma, 35.0, places=5)
        self.assertAlmostEqual(hae4.zeta, 23.333333333333343, places=5)
        self.assertAlmostEqual(hae4.eta, 0.0, places=5)
        cs3 = copy.deepcopy(cs1)
        cs3.sigma_11 = 1.0
        cs3.sigma_22 = 1000.0
        cs3.sigma_33 = 1000.0
        mary2 = cs3.maryland_values
        self.assertAlmostEqual(mary2.sigma_iso, 667.0, places=5)
        self.assertAlmostEqual(mary2.omega, 999.0, places=5)
        self.assertAlmostEqual(mary2.kappa, 1.0, places=5)
        cs3.sigma_22 = 1.0
        mary3 = cs3.maryland_values
        self.assertAlmostEqual(mary3.sigma_iso, 334.0, places=5)
        self.assertAlmostEqual(mary3.omega, 999.0, places=5)
        self.assertAlmostEqual(mary3.kappa, -1.0, places=5)
        cs3.sigma_22 = 500.5
        mary4 = cs3.maryland_values
        self.assertAlmostEqual(mary4.sigma_iso, 500.5, places=5)
        self.assertAlmostEqual(mary4.omega, 999.0, places=5)
        self.assertAlmostEqual(mary4.kappa, 0.0, places=5)
        cs4 = NMRChemicalShiftNotation(68.1733, 0.8337, 195.0788)
        cs5 = NMRChemicalShiftNotation(0.8337, 195.0788, 68.1733)
        self.assertAlmostEqual(cs5.sigma_11, cs4.sigma_11)
        self.assertAlmostEqual(cs5.sigma_22, cs4.sigma_22)
        self.assertAlmostEqual(cs5.sigma_33, cs4.sigma_33)
        cs6 = NMRChemicalShiftNotation(195.0788, 68.1733, 0.8337)
        self.assertAlmostEqual(cs6.sigma_11, cs4.sigma_11)
        self.assertAlmostEqual(cs6.sigma_22, cs4.sigma_22)
        self.assertAlmostEqual(cs6.sigma_33, cs4.sigma_33)
        d1 = cs1.as_dict()
        cs7 = NMRChemicalShiftNotation.from_dict(d1)
        hae7 = cs7.haeberlen_values
        self.assertAlmostEqual(hae1, hae7, places=5)
        j1 = cs1.to_json()
        cs8 = json.loads(j1, cls=MontyDecoder)
        hae8 = cs8.haeberlen_values
        self.assertAlmostEqual(hae1, hae8, places=5)
