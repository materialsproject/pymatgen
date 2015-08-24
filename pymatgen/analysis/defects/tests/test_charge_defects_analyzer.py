# coding: utf-8

from __future__ import unicode_literals

import unittest
import sys

from pymatgen.analysis.defects.point_defects import *
from pymatgen.analysis.defects.charge_defects_analyzer import *
#from pymatgen.core.structure import Structure
#from pymatgen.core.periodic_table import Element
#from pymatgen.analysis.bond_valence import BVAnalyzer
from monty.os.path import which
from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
#from pymatgen.io.cifio import CifParser

sxdefectalign_present = which('sxdefectalign')
test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class ParsedChargeDefectTest(unittest.TestCase):
    """
    This class is just an aggregate for the data. No tests are necessary.
    """
    pass


@unittest.skipIf(not sxdefectalign_present, "sxdefectalign not present.")
class CorrectionTest(unittest.TestCase):
    """
    Tests for the correction obtained Freysoldt method
    """
    def setUp(self):
        pass

    def test_get_correction(self):
        pass


class ChargeDefectsAnalyzerTest(unittest.TestCase):
    def setUp(self):
        self.data_defects = loadfn(os.path.join(test_dir,
            'InAs_cdt_output.json'), cls=MontyDecoder)
        self.corrections = loadfn(os.path.join(test_dir,
            'InAs_cdt_correction.json'), cls=MontyDecoder)
        self.e_vbm = self.data_defects['vbm']
        self.mu = {}
        mu_key = self.data_defects['mu_range'].keys()[0]
        mu = self.data_defects['mu_range'][mu_key]['poor']
        for key in mu:
            self.mu[key.encode('ascii','ignore')] = mu[key]
        self.band_gap = self.data_defects['gap']['energy']
        self.entry_bulk = self.data_defects['bulk_entry']
        self.defect_analyzer = ChargeDefectsAnalyzer(
                self.entry_bulk, self.e_vbm, self.mu, self.band_gap)

    def test_as_dict(self):
        d = self.defect_analyzer.as_dict()
        self.assertIn('formation_energies',d)
        self.assertIn('structure',d['entry_bulk'])
        self.defect_analyzer.add_parsed_defect(self.data_defects['defects'][0])
        d = self.defect_analyzer.as_dict()
        self.assertIsNotNone(d['formation_energies'])
        self.assertIsNotNone(d['defects'])

    def test_from_dict(self):
        orig_da = self.defect_analyzer
        orig_da.add_parsed_defect(self.data_defects['defects'][0])
        d = orig_da.as_dict()
        dupl_da = ChargeDefectsAnalyzer.from_dict(d)
        self.assertIsNotNone(dupl_da._formation_energies)
        self.assertAlmostEqual(dupl_da._e_vbm, orig_da._e_vbm)
        self.assertAlmostEqual(dupl_da._band_gap,orig_da._band_gap)
        for key in dupl_da._mu_elts:
            self.assertIn(key,orig_da._mu_elts)

    def test_add_parsed_defect(self):
        self.defect_analyzer.add_parsed_defect(self.data_defects['defects'][0])
        self.assertAlmostEqual(self.defect_analyzer._formation_energies[0],
                6.8216, places=4)

    def test_change_charge_correction(self):
        added_defect = self.data_defects['defects'][0]
        self.defect_analyzer.add_parsed_defect(added_defect)
        orig_form_energy = self.defect_analyzer._formation_energies[0]
        for defect in self.corrections[added_defect._name]:
            if defect['charge'] == added_defect._charge:
                correction = defect['correction'][0]
        self.defect_analyzer.change_charge_correction(0,correction)
        final_form_energy = self.defect_analyzer._formation_energies[0]
        self.assertNotEqual(orig_form_energy, final_form_energy)

    def test_change_charge_correction_zerocharge(self):
        added_defect = self.data_defects['defects'][0]
        added_defect._charge = 0
        self.defect_analyzer.add_parsed_defect(added_defect)
        orig_form_energy = self.defect_analyzer._formation_energies[0]
        for defect in self.corrections[added_defect._name]:
            if defect['charge'] == added_defect._charge:
                correction = defect['correction'][0]
        self.defect_analyzer.change_charge_correction(0,correction)
        final_form_energy = self.defect_analyzer._formation_energies[0]
        self.assertAlmostEqual(orig_form_energy, final_form_energy)

    def test_correct_bg_simple(self):
        added_defect = self.data_defects['defects'][0]
        for defect in self.corrections[added_defect._name]:
            if defect['charge'] == added_defect._charge:
                correction = defect['correction'][0]
        added_defect.charge_correction = correction
        self.defect_analyzer.add_parsed_defect(added_defect)
        orig_form_energy = self.defect_analyzer._formation_energies[0]
        self.defect_analyzer.correct_bg_simple(0.1,0.1)
        final_form_energy = self.defect_analyzer._formation_energies[0]
        self.assertAlmostEqual(abs(final_form_energy-orig_form_energy),0.1)

    def test_correct_bg(self):
        added_defect = self.data_defects['defects'][0]
        self.defect_analyzer.add_parsed_defect(added_defect)
        orig_form_energy = self.defect_analyzer._formation_energies[0]
        self.defect_analyzer.correct_bg({
            'vac_1_In': {'type': 'cbm_like', 'q*': -3}}, 0.1, 0.2)
        final_form_energy = self.defect_analyzer._formation_energies[0]
        self.assertAlmostEqual(abs(final_form_energy-orig_form_energy),0.3)

    def test_get_defects_concentration(self):
        self.defect_analyzer.add_parsed_defect(self.data_defects['defects'][0])
        conc = self.defect_analyzer.get_defects_concentration()[0]['conc']
        self.assertAlmostEqual(conc, 4.2818e-87, places=90)

    @unittest.skip("Not validated at present")
    def test_get_eq_ef(self):
        pass

    @unittest.skip("Not validated at present")
    def test_get_non_eq_ef(self):
        pass


if __name__ == "__main__":
    unittest.main()
