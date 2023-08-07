from __future__ import annotations

import numpy as np
import pytest
import scipy.special

from pymatgen.io.vasp.optics import DielectricFunctionCalculator, delta_func, delta_methfessel_paxton, step_func
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Jimmy-Xuan Shen"
__copyright__ = "Copyright 2022, The Materials Project"
__email__ = "jmmshn@gmail.com"


class TestVasprun(PymatgenTest):
    def test_optics(self):
        eps_data_path = f"{TEST_FILES_DIR}/reproduce_eps"
        vrun = Vasprun(f"{eps_data_path}/vasprun.xml")
        dfc = DielectricFunctionCalculator.from_directory(eps_data_path)
        egrid, eps = dfc.get_epsilon(0, 0)

        assert egrid[0] == 0
        assert egrid[-1] == 59.3802
        assert len(egrid) == len(eps) == 3000

        _, eps_real_ref, eps_imag_ref = vrun.dielectric
        eps_real_ref = np.array(eps_real_ref)[:, 0]
        eps_imag_ref = np.array(eps_imag_ref)[:, 0]
        assert np.max(np.abs(eps_real_ref - np.real(eps))) / np.max(np.abs(np.real(eps))) < 0.01

        # masking with all zeros
        mask = np.zeros_like(dfc.cder, dtype=float)
        _, eps = dfc.get_epsilon(0, 0, mask=mask)
        assert np.max(np.imag(eps)) < 1e-10

        # masking with all ones
        mask = np.ones_like(dfc.cder, dtype=float)
        _, eps = dfc.get_epsilon(0, 0, mask=mask)
        assert np.max(np.abs(eps_real_ref - np.real(eps))) / np.max(np.abs(np.real(eps))) < 0.01

        _, eps = dfc.get_epsilon(0, 1)
        assert np.max(np.abs(eps)) < 0.1

        # plotting
        mask = np.ones_like(dfc.cder, dtype=float)
        x_val, y_val, text = dfc.plot_weighted_transition_data(0, 0, mask=mask, min_val=0.001)
        assert len(x_val) == len(y_val) == len(text)


def test_delta_func():
    x = np.array([0, 1, 2, 3, 4, 5])

    # ismear < -1
    with pytest.raises(ValueError, match="Delta function not implemented for ismear < -1"):
        delta_func(x, -2)

    # ismear == -1
    assert np.all(delta_func(x, -1) == step_func(x, -1) * (1 - step_func(x, -1)))

    # ismear == 0
    assert np.all(delta_func(x, 0) == np.exp(-(x * x)) / np.sqrt(np.pi))

    # ismear > 0
    for ismear in [1, 2, 3]:
        assert np.all(delta_func(x, ismear) == delta_methfessel_paxton(x, ismear))


def test_step_func():
    # array of positive values
    x = np.array([1, 2, 3, 4, 5])
    assert np.allclose(step_func(x, -1), 1 / (1.0 + np.exp(-x)))

    # array of negative values
    x = np.array([-1, -2, -3, -4, -5])
    assert np.allclose(step_func(x, -1), 1 / (1.0 + np.exp(-x)))

    # array that includes zero
    x = np.array([-1, 0, 1])
    assert np.allclose(step_func(x, -1), 1 / (1.0 + np.exp(-x)))

    # ismear == 0
    x = np.array([1, 2, 3, 4, 5])
    assert np.allclose(step_func(x, 0), 0.5 + 0.5 * scipy.special.erf(x))
