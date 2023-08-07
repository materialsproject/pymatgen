from __future__ import annotations

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.analysis.nmr import ChemicalShielding, ElectricFieldGradient
from pymatgen.util.testing import PymatgenTest


class TestChemicalShieldingNotation(PymatgenTest):
    def test_construction(self):
        cs = ChemicalShielding(np.arange(9).reshape((3, 3)))
        assert cs.shape == (3, 3)

        cs = ChemicalShielding([1, 2, 3])
        assert cs.shape == (3, 3)
        assert_array_equal(np.diag(cs), [1, 2, 3])

    def test_principal_axis_system(self):
        cs = ChemicalShielding([1, 2, 3])
        assert_array_equal(cs.principal_axis_system, cs)

        cs = ChemicalShielding(np.arange(9).reshape((3, 3)))
        assert np.allclose(
            np.diag(cs.principal_axis_system),
            [-1.74596669e00, -1.53807726e-15, 1.37459667e01],
            atol=1e-5,
        )

    def test_notations(self):
        cs = ChemicalShielding.from_maryland_notation(195.0788, 68.1733, 0.8337)
        hae1 = cs.haeberlen_values
        assert hae1.sigma_iso == approx(195.0788, abs=1e-5)
        assert hae1.delta_sigma_iso == approx(-65.33899505250002, abs=1e-5)
        assert hae1.zeta == approx(-43.559330035000016, abs=1e-5)
        assert hae1.eta == approx(0.13013537835511454, abs=1e-5)
        meh1 = cs.mehring_values
        assert meh1.sigma_iso == approx(195.0788, abs=1e-5)
        assert meh1.sigma_11 == approx(151.51946996499998, abs=1e-5)
        assert meh1.sigma_22 == approx(214.02416007, abs=1e-5)
        assert meh1.sigma_33 == approx(219.69276996500002, abs=1e-5)
        mary1 = cs.maryland_values
        assert mary1.sigma_iso == approx(195.0788, abs=1e-5)
        assert mary1.omega == approx(68.1733, abs=1e-5)
        assert mary1.kappa == approx(0.8337, abs=1e-5)


class TestElectricFieldGradient(PymatgenTest):
    def test_construction(self):
        efg = ElectricFieldGradient(np.arange(9).reshape((3, 3)))
        assert efg.shape == (3, 3)

        efg = ElectricFieldGradient([1, 2, 3])
        assert efg.shape == (3, 3)

    def test_principal_axis_system(self):
        efg = ElectricFieldGradient([1, 2, 3])
        assert_array_equal(efg.principal_axis_system, efg)

        efg = ElectricFieldGradient(np.arange(9).reshape((3, 3)))
        assert np.allclose(
            np.diag(efg.principal_axis_system),
            [-1.3484692e00, -1.1543332e-15, 1.3348469e01],
            atol=1e-5,
        )

    def test_Attributes(self):
        efg = ElectricFieldGradient([[11.11, 1.371, 2.652], [1.371, 3.635, -3.572], [2.652, -3.572, -14.746]])
        assert efg.V_yy == approx(11.516, abs=1e-3)
        assert efg.V_xx == approx(4.204, abs=1e-3)
        assert efg.V_zz == approx(-15.721, abs=1e-3)
        assert efg.asymmetry == approx(0.465, abs=1e-3)
        assert efg.coupling_constant("Al") == approx(5.573, abs=1e-3)
