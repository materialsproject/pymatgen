# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
from pymatgen.analysis.spectra.spectrum import XANES
import unittest
import json
from monty.json import MontyDecoder

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files/spectrum_test")

with open(os.path.join(test_dir,'Pd2O.json')) as fp:
    spect_data_dict = json.load(fp, cls=MontyDecoder)


class XANESSetTest(unittest.TestCase):
    def setUp(self):
        self.xaneobj = XANES.from_dict(spect_data_dict)

    def test_structure(self):
        self.assertEqual(self.xaneobj.structure.composition.reduced_formula, 'Pd2O')
        self.assertEqual(self.xaneobj.absorption_specie, 'Pd')

    def test_find_e0(self):
        self.xaneobj.find_e0()
        self.assertAlmostEqual(24374.508999999998, self.xaneobj.e0)
