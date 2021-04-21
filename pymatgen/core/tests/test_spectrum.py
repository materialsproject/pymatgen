# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

import numpy as np

from pymatgen.core.spectrum import Spectrum
from pymatgen.util.testing import PymatgenTest


class SpectrumTest(PymatgenTest):
    def setUp(self):
        self.spec1 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100))
        self.spec2 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100))

        self.multi_spec1 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100, 2))
        self.multi_spec2 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100, 2))

    def test_normalize(self):
        self.spec1.normalize(mode="max")
        self.assertAlmostEqual(np.max(self.spec1.y), 1)
        self.spec1.normalize(mode="sum")
        self.assertAlmostEqual(np.sum(self.spec1.y), 1)

        self.multi_spec1.normalize(mode="sum")
        self.assertAlmostEqual(np.sum(self.multi_spec1.y[:, 0]), 1)
        self.assertAlmostEqual(np.sum(self.multi_spec1.y[:, 1]), 1)

        # XRD style mode.
        self.spec1.normalize(mode="max", value=100)
        self.assertAlmostEqual(np.max(self.spec1.y), 100)

    def test_operators(self):
        scaled_spect = 3 * self.spec1 + self.spec2
        self.assertTrue(np.allclose(scaled_spect.y, 3 * self.spec1.y + self.spec2.y))
        self.assertAlmostEqual(
            self.spec1.get_interpolated_value(0.05),
            (self.spec1.y[0] + self.spec1.y[1]) / 2,
        )

        scaled_spect = self.spec1 - self.spec2
        self.assertTrue(np.allclose(scaled_spect.y, self.spec1.y - self.spec2.y))

        scaled_spect = self.spec1 / 3
        self.assertTrue(np.allclose(scaled_spect.y, self.spec1.y / 3))

        scaled_spect = 3 * self.multi_spec1 + self.multi_spec2
        self.assertTrue(np.allclose(scaled_spect.y, 3 * self.multi_spec1.y + self.multi_spec2.y))
        self.assertArrayAlmostEqual(
            self.multi_spec1.get_interpolated_value(0.05),
            (self.multi_spec1.y[0, :] + self.multi_spec1.y[1, :]) / 2,
        )

    def test_smear(self):
        y = np.zeros(100)
        y[25] = 1
        y[50] = 1
        y[75] = 1
        spec = Spectrum(np.linspace(-10, 10, 100), y)
        spec.smear(0.3)
        self.assertFalse(np.allclose(y, spec.y))
        self.assertAlmostEqual(sum(y), sum(spec.y))
        spec = Spectrum(np.linspace(-10, 10, 100), y)
        spec.smear(0.3, func="lorentzian")
        self.assertFalse(np.allclose(y, spec.y))
        self.assertAlmostEqual(sum(y), sum(spec.y))

        y = np.array(self.multi_spec1.y)
        self.multi_spec1.smear(0.2)
        self.assertFalse(np.allclose(y, self.multi_spec1.y))
        self.assertArrayAlmostEqual(np.sum(y, axis=0), np.sum(self.multi_spec1.y, axis=0))

    def test_str(self):
        # Just make sure that these methods work.
        self.assertIsNotNone(str(self.spec1))
        self.assertIsNotNone(str(self.multi_spec1))

    def test_copy(self):
        spec1copy = self.spec1.copy()
        spec1copy.y[0] = spec1copy.y[0] + 1
        self.assertNotEqual(spec1copy.y[0], self.spec1.y[0])
        self.assertEqual(spec1copy.y[1], self.spec1.y[1])


if __name__ == "__main__":
    unittest.main()
