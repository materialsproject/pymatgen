from __future__ import annotations

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.util.testing import PymatgenTest


class TestStress(PymatgenTest):
    def setUp(self):
        self.rand_stress = Stress(np.random.randn(3, 3))
        self.symm_stress = Stress([[0.51, 2.29, 2.42], [2.29, 5.14, 5.07], [2.42, 5.07, 5.33]])
        self.non_symm = Stress([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.2, 0.5, 0.5]])

    def test_properties(self):
        # mean_stress
        assert self.rand_stress.mean_stress == approx(
            (self.rand_stress[0, 0] + self.rand_stress[1, 1] + self.rand_stress[2, 2]) / 3
        )
        assert self.symm_stress.mean_stress == approx(3.66)
        # deviator_stress
        assert np.allclose(
            self.symm_stress.deviator_stress,
            Stress([[-3.15, 2.29, 2.42], [2.29, 1.48, 5.07], [2.42, 5.07, 1.67]]),
        )
        assert np.allclose(
            self.non_symm.deviator_stress,
            [[-0.2666666667, 0.2, 0.3], [0.4, 0.133333333, 0.6], [0.2, 0.5, 0.133333333]],
        )
        # deviator_principal_invariants
        assert np.allclose(self.symm_stress.dev_principal_invariants, [0, 44.2563, 111.953628])
        # von_mises
        assert self.symm_stress.von_mises == approx(11.52253878275)
        # piola_kirchoff 1, 2
        f = Deformation.from_index_amount((0, 1), 0.03)
        assert np.allclose(
            self.symm_stress.piola_kirchoff_1(f),
            [[0.4413, 2.29, 2.42], [2.1358, 5.14, 5.07], [2.2679, 5.07, 5.33]],
        )
        assert np.allclose(
            self.symm_stress.piola_kirchoff_2(f),
            [[0.377226, 2.1358, 2.2679], [2.1358, 5.14, 5.07], [2.2679, 5.07, 5.33]],
        )
        # voigt
        assert list(self.symm_stress.voigt) == [0.51, 5.14, 5.33, 5.07, 2.42, 2.29]

        with pytest.warns(
            UserWarning, match="Tensor is not symmetric, information may be lost in voigt conversion"
        ) as warns:
            _ = self.non_symm.voigt
        assert len(warns) == 1
