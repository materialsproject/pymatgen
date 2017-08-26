# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
from pymatgen.core.spectrum import Spectrum
import unittest
import json
from monty.json import MontyDecoder
import numpy as np


class SpectrumTest(unittest.TestCase):

    def setUp(self):
        self.spec1 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100))
        self.spec2 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100))

    def test_normalize(self):
        self.spec1.normalize(mode="max")
        self.assertAlmostEqual(np.max(self.spec1.y), 1)
        self.spec1.normalize(mode="sum")
        self.assertAlmostEqual(np.sum(self.spec1.y), 1)

    def test_operators(self):
        scaled_spect = 3 * self.spec1 + self.spec2
        self.assertTrue(np.allclose(scaled_spect.y, 3 * self.spec1.y + self.spec2.y))
        self.assertAlmostEqual(self.spec1.get_interpolated_value(0.05),
                               (self.spec1.y[0] + self.spec1.y[1]) / 2)

        scaled_spect = self.spec1 - self.spec2
        self.assertTrue(np.allclose(scaled_spect.y, self.spec1.y - self.spec2.y))

    def test_smear(self):
        y = np.array(self.spec1.y)
        self.spec1.smear(0.2)
        self.assertFalse(np.allclose(y, self.spec1.y))
        self.assertAlmostEqual(sum(y), sum(self.spec1.y))

if __name__ == '__main__':
    unittest.main()