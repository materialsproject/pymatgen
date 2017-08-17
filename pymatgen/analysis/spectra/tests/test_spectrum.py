import os
from pymatgen.analysis.spectra.spectrum import XanesFEFF
import unittest

test_dir = os.path.join(os.path.dirname(__file__),"..", "..", "..", "..",
                        "test_files/spectrum_test")



class XanesFEFFSetTest(unittest.TestCase):
    def setUp(self):
        self.xaneobj = XanesFEFF.from_file(test_dir+'/Pd2O.json')

    def test_structure(self):
        self.assertEqual(self.xaneobj.absorbing_atom_index, 2)
        self.assertEqual(self.xaneobj.absorbing_specie, 'Pd')

    def test_e0_interpolate(self):
        self.xaneobj.e0_interpolate()
        self.assertAlmostEqual(24377.535, self.xaneobj.e0)