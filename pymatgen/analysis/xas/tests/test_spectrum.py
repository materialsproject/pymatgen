# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
from pymatgen.analysis.xas.spectrum import XANES
import unittest
from pymatgen.util.testing import PymatgenTest
import json
from monty.json import MontyDecoder
import numpy as np

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files/spectrum_test")

with open(os.path.join(test_dir, 'Pd2O.json')) as fp:
    spect_data_dict = json.load(fp, cls=MontyDecoder)


class XANESSetTest(PymatgenTest):
    def setUp(self):
        self.xaneobj = XANES.from_dict(spect_data_dict)

    def test_find_e0(self):
        self.xaneobj.find_e0()
        self.assertAlmostEqual(24374.508999999998, self.xaneobj.e0)

    def test_normalization(self):
        self.xaneobj.normalize(mode="sum")
        self.assertAlmostEqual(1.0, np.sum(self.xaneobj.y))

    def test_add_mul(self):
        scaled_spect = self.xaneobj + self.xaneobj
        scaled_spect2 = self.xaneobj * 3
        self.assertTrue(np.allclose(scaled_spect.y, 2 * self.xaneobj.y))
        self.assertTrue(np.allclose(scaled_spect2.y, 3 * self.xaneobj.y))
        self.assertAlmostEqual(0.348883, self.xaneobj.get_interpolated_value(24374.509))

    def test_to_from_dict(self):
        s = XANES.from_dict(self.xaneobj.as_dict())
        self.assertArrayAlmostEqual(s.y, self.xaneobj.y)


if __name__ == '__main__':
    unittest.main()