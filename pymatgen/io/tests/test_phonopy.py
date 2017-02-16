from __future__ import division, unicode_literals

import os
import unittest
from pymatgen.core.periodic_table import Element
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.phonopy import *

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "phonopy")


class PhonopyParserTest(PymatgenTest):

    def test_get_ph_bs(self):
        ph_bs = get_ph_bs_symm_line(os.path.join(test_dir, 'NaCl_band.yaml'), has_nac=True)
        
        self.assertAlmostEqual(ph_bs.bands[1][10], 0.7753555184)
        self.assertAlmostEqual(ph_bs.bands[5][100], 5.2548379776)
        self.assertArrayEqual(ph_bs.bands.shape, (6, 204))
        self.assertArrayEqual(ph_bs.eigendisplacements.shape, (6, 204, 2, 3))
        self.assertArrayAlmostEqual(ph_bs.eigendisplacements[3][50][0],
                              [0. + 0.j, 0.14166569 + 0.04098339j,
                               -0.14166569 - 0.04098339j])
        self.assertTrue(ph_bs.has_eigendisplacements, True)
        self.assertArrayEqual(ph_bs.min_freq()[0].frac_coords, [0, 0, 0])
        self.assertAlmostEqual(ph_bs.min_freq()[1], -0.03700895020)
        self.assertTrue(ph_bs.has_imaginary_freq())
        self.assertFalse(ph_bs.has_imaginary_freq(tol=0.5))
        self.assertArrayAlmostEqual(ph_bs.asr_breaking(), [-0.0370089502, -0.0370089502, -0.0221388897])
        self.assertEqual(ph_bs.nb_bands, 6)
        self.assertEqual(ph_bs.nb_qpoints, 204)
        self.assertArrayAlmostEqual(ph_bs.qpoints[1].frac_coords, [0.01, 0, 0])
        self.assertTrue(ph_bs.has_nac)
        self.assertAlmostEqual(ph_bs.get_nac_frequencies_along_dir([1, 1, 0])[3], 4.6084532143)
        self.assertIsNone(ph_bs.get_nac_frequencies_along_dir([1, 1, 1]))
        self.assertArrayAlmostEqual(ph_bs.get_nac_eigendisplacements_along_dir([1, 1, 0])[3][1],
                                    [(0.1063906409128248 + 0j), 0j, 0j])
        self.assertIsNone(ph_bs.get_nac_eigendisplacements_along_dir([1, 1, 1]))

    def test_get_ph_dos(self):
        dos = get_ph_dos(os.path.join(test_dir, 'NaCl_total_dos.dat'))

        self.assertAlmostEqual(dos.densities[15], 0.0001665998)
        self.assertAlmostEqual(dos.frequencies[20], 0.0894965119)
        self.assertAlmostEqual(dos.get_interpolated_value(3.), 1.2915532670115628)
        self.assertEqual(len(dos.frequencies), 201)
        self.assertEqual(len(dos.densities), 201)

    def test_get_complete_dos(self):
        cdos = get_complete_ph_dos(os.path.join(test_dir, 'NaCl_partial_dos.dat'),
                                   os.path.join(test_dir, 'NaCl_phonopy.yaml'))
        site_Na = cdos.structure[0]
        site_Cl = cdos.structure[1]

        self.assertEqual(len(cdos.frequencies), 201)
        self.assertAlmostEqual(cdos.pdos[site_Na][30],  0.008058208)
        self.assertAlmostEqual(cdos.pdos[site_Cl][30],  0.0119040783)

        self.assertIn(Element.Na, cdos.get_element_dos())
        self.assertIn(Element.Cl, cdos.get_element_dos())


if __name__ == '__main__':
    unittest.main()
