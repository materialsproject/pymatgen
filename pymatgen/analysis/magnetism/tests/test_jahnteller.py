# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.io.cif import CifParser
from pymatgen.analysis.magnetism.jahnteller import *

import unittest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class JahnTellerTest(unittest.TestCase):

    def setUp(self):
        self.jt = JahnTellerAnalyzer()

    def test_jahn_teller_species_analysis(self):
        # 1 d-shell electron
        m = self.jt.get_magnitude_of_effect_from_species('Ti3+', '', 'oct')
        self.assertEqual(m, "weak")

        # 2 d-shell electrons
        m = self.jt.get_magnitude_of_effect_from_species('Ti2+', '', 'oct')
        self.assertEqual(m, "weak")
        m = self.jt.get_magnitude_of_effect_from_species('V3+', '', 'oct')
        self.assertEqual(m, "weak")

        # 3
        m = self.jt.get_magnitude_of_effect_from_species('V2+', '', 'oct')
        self.assertEqual(m, "none")
        m = self.jt.get_magnitude_of_effect_from_species('Cr3+', '', 'oct')
        self.assertEqual(m, "none")

        # 4
        m = self.jt.get_magnitude_of_effect_from_species('Cr2+', 'high', 'oct')
        self.assertEqual(m, "strong")
        m = self.jt.get_magnitude_of_effect_from_species('Cr2+', 'low', 'oct')
        self.assertEqual(m, "weak")
        m = self.jt.get_magnitude_of_effect_from_species('Mn3+', 'high', 'oct')
        self.assertEqual(m, "strong")
        m = self.jt.get_magnitude_of_effect_from_species('Mn3+', 'low', 'oct')
        self.assertEqual(m, "weak")

        # 5
        m = self.jt.get_magnitude_of_effect_from_species('Mn2+', 'high', 'oct')
        self.assertEqual(m, "none")
        m = self.jt.get_magnitude_of_effect_from_species('Mn2+', 'low', 'oct')
        self.assertEqual(m, "weak")
        m = self.jt.get_magnitude_of_effect_from_species('Fe3+', 'high', 'oct')
        self.assertEqual(m, "none")
        m = self.jt.get_magnitude_of_effect_from_species('Fe3+', 'low', 'oct')
        self.assertEqual(m, "weak")

        # 6
        m = self.jt.get_magnitude_of_effect_from_species('Fe2+', 'high', 'oct')
        self.assertEqual(m, "weak")
        m = self.jt.get_magnitude_of_effect_from_species('Fe2+', 'low', 'oct')
        self.assertEqual(m, "none")
        m = self.jt.get_magnitude_of_effect_from_species('Co3+', 'high', 'oct')
        self.assertEqual(m, "weak")
        m = self.jt.get_magnitude_of_effect_from_species('Co3+', 'low', 'oct')
        self.assertEqual(m, "none")

        # 7
        m = self.jt.get_magnitude_of_effect_from_species('Co2+', 'high', 'oct')
        self.assertEqual(m, "weak")
        m = self.jt.get_magnitude_of_effect_from_species('Co2+', 'low', 'oct')
        self.assertEqual(m, "strong")

        # 8
        m = self.jt.get_magnitude_of_effect_from_species('Ni2+', '', 'oct')
        self.assertEqual(m, "none")

        # 9
        m = self.jt.get_magnitude_of_effect_from_species('Cu2+', '', 'oct')
        self.assertEqual(m, "strong")

        # 10
        m = self.jt.get_magnitude_of_effect_from_species('Cu+', '', 'oct')
        self.assertEqual(m, "none")
        m = self.jt.get_magnitude_of_effect_from_species('Zn2+', '', 'oct')
        self.assertEqual(m, "none")

    def test_jahn_teller_structure_analysis(self):
        parser = CifParser(os.path.join(test_dir, 'LiFePO4.cif'))
        LiFePO4 = parser.get_structures()[0]

        parser = CifParser(os.path.join(test_dir, 'Fe3O4.cif'))
        Fe3O4 = parser.get_structures()[0]

        self.assertTrue(self.jt.is_jahn_teller_active(LiFePO4))
        self.assertTrue(self.jt.is_jahn_teller_active(Fe3O4))

        LiFePO4_analysis = {
            'active': True,
            'strength': 'weak',
            'sites': [
                {
                    'ligand': 'O2-',
                    'ligand_bond_length_spread': 0.2111,
                    'ligand_bond_lengths': set([2.2951, 2.2215, 2.2383, 2.1382, 2.084, 2.0863]),
                    'strength': 'weak',
                    'motif': 'oct',
                    'motif_order_parameter': 0.1441,
                    'site_indices': [4, 5, 6, 7],
                    'species': 'Fe2+',
                    'spin_state': 'unknown'
                }
            ]
        }
        jt_predicted = self.jt.get_analysis(LiFePO4)
        # order does not matter
        jt_predicted['sites'][0]['ligand_bond_lengths'] = set(jt_predicted['sites'][0]['ligand_bond_lengths'])
        self.assertDictEqual(LiFePO4_analysis, jt_predicted)

    def test_mu_so(self):
        SpeciesCo = Specie(symbol='Co', oxidation_state=4)
        self.assertAlmostEqual(np.sqrt(3), JahnTellerAnalyzer.mu_so(SpeciesCo, 'oct', 'low'))
        self.assertAlmostEqual(np.sqrt(35), JahnTellerAnalyzer.mu_so(SpeciesCo, 'oct', 'high'))
        SpeciesNa = Specie(symbol='Na', oxidation_state=1)
        self.assertEqual(None, JahnTellerAnalyzer.mu_so(SpeciesNa, 'oct', 'high'))


if __name__ == '__main__':
    unittest.main()
