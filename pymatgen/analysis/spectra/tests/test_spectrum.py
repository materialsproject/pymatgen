import os
from pymatgen.analysis.spectra.spectrum import XANES
import unittest
import json

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files/spectrum_test")


class XANESSetTest(unittest.TestCase):
    def setUp(self):
        json_file = os.path.join(test_dir, 'Pd2O.json')
        with open(json_file, 'r') as f:
            data = json.load(f)
        self.xaneobj = XANES.from_dict(data)

    def test_structure(self):
        self.assertEqual(self.xaneobj.structure.composition.reduced_formula, 'Pd2O')
        self.assertEqual(self.xaneobj.absorption_specie, 'Pd')

    def test_find_e0(self):
        self.xaneobj.find_e0()
        self.assertAlmostEqual(24374.508999999998, self.xaneobj.e0)
