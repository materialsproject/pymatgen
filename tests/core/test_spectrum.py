from __future__ import annotations

import numpy as np
from pytest import approx
from scipy import stats

from pymatgen.core.spectrum import Spectrum
from pymatgen.util.testing import PymatgenTest


class TestSpectrum(PymatgenTest):
    def setUp(self):
        self.spec1 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100))
        self.spec2 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100))

        self.multi_spec1 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100, 2))
        self.multi_spec2 = Spectrum(np.arange(0, 10, 0.1), np.random.randn(100, 2))

    def test_normalize(self):
        self.spec1.normalize(mode="max")
        assert np.max(self.spec1.y) == approx(1)
        self.spec1.normalize(mode="sum")
        assert np.sum(self.spec1.y) == approx(1)

        self.multi_spec1.normalize(mode="sum")
        assert np.sum(self.multi_spec1.y[:, 0]) == approx(1)
        assert np.sum(self.multi_spec1.y[:, 1]) == approx(1)

        # XRD style mode.
        self.spec1.normalize(mode="max", value=100)
        assert np.max(self.spec1.y) == approx(100)

    def test_operators(self):
        scaled_spect = 3 * self.spec1 + self.spec2
        assert np.allclose(scaled_spect.y, 3 * self.spec1.y + self.spec2.y)
        assert self.spec1.get_interpolated_value(0.05) == approx((self.spec1.y[0] + self.spec1.y[1]) / 2)

        scaled_spect = self.spec1 - self.spec2
        assert np.allclose(scaled_spect.y, self.spec1.y - self.spec2.y)

        scaled_spect = self.spec1 / 3
        assert np.allclose(scaled_spect.y, self.spec1.y / 3)

        scaled_spect = 3 * self.multi_spec1 + self.multi_spec2
        assert np.allclose(scaled_spect.y, 3 * self.multi_spec1.y + self.multi_spec2.y)
        assert np.allclose(
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
        assert not np.allclose(y, spec.y)
        assert sum(y) == approx(sum(spec.y))

        # Test direct callable use of smearing.
        spec2 = Spectrum(np.linspace(-10, 10, 100), y)
        spec2.smear(0, func=lambda x: stats.norm.pdf(x, scale=0.3))
        assert np.allclose(spec.y, spec2.y)

        spec = Spectrum(np.linspace(-10, 10, 100), y)
        spec.smear(0.3, func="lorentzian")
        assert not np.allclose(y, spec.y)
        assert sum(y) == approx(sum(spec.y))

        y = np.array(self.multi_spec1.y)
        self.multi_spec1.smear(0.2)
        assert not np.allclose(y, self.multi_spec1.y)
        assert np.allclose(np.sum(y, axis=0), np.sum(self.multi_spec1.y, axis=0))

    def test_str(self):
        # Just make sure that these methods work.
        assert str(self.spec1) is not None
        assert str(self.multi_spec1) is not None

    def test_copy(self):
        spec1copy = self.spec1.copy()
        spec1copy.y[0] = spec1copy.y[0] + 1
        assert spec1copy.y[0] != self.spec1.y[0]
        assert spec1copy.y[1] == self.spec1.y[1]
