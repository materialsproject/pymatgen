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


class XANESTest(PymatgenTest):
    def setUp(self):
        self.xanes = XANES.from_dict(spect_data_dict)

    def test_e0(self):
        self.assertAlmostEqual(24374.508999999998, self.xanes.e0)

    def test_normalization(self):
        self.xanes.normalize(mode="sum")
        self.assertAlmostEqual(1.0, np.sum(self.xanes.y))

    def test_add_mul(self):
        scaled_spect = self.xanes + self.xanes
        scaled_spect2 = self.xanes * 3
        self.assertTrue(np.allclose(scaled_spect.y, 2 * self.xanes.y))
        self.assertTrue(np.allclose(scaled_spect2.y, 3 * self.xanes.y))
        self.assertAlmostEqual(0.348883,
                               self.xanes.get_interpolated_value(24374.509))

    def test_to_from_dict(self):
        s = XANES.from_dict(self.xanes.as_dict())
        self.assertArrayAlmostEqual(s.y, self.xanes.y)

    def test_attributes(self):
        self.assertArrayEqual(self.xanes.energy, self.xanes.x)
        self.assertArrayEqual(self.xanes.intensity, self.xanes.y)

    def test_str(self):
        self.assertIsNotNone(str(self.xanes))


if __name__ == '__main__':
    unittest.main()
