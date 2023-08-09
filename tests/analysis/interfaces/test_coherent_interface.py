from __future__ import annotations

import numpy as np

from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    from_2d_to_3d,
    get_2d_transform,
    get_rot_3d_for_2d,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest


class TestInterfaceBuilder(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        si_struct = cls.get_structure("Si")
        sio2_struct = cls.get_structure("SiO2")
        cls.si_conventional = SpacegroupAnalyzer(si_struct).get_conventional_standard_structure()
        cls.sio2_conventional = SpacegroupAnalyzer(sio2_struct).get_conventional_standard_structure()

    def test_utils(self):
        assert np.allclose(from_2d_to_3d([[1, 2], [3, 4]]), [[1, 2, 0], [3, 4, 0], [0, 0, 1]])
        assert np.allclose(get_2d_transform([[1, 0], [0, 1]], [[1, 2], [3, 4]]), [[1, 2], [3, 4]])
        assert np.allclose(
            get_rot_3d_for_2d([[1, 0, 0], [0, 1, 0]], [[1, 1, 0], [0, 1, 1]]),
            [
                [0.78867513, -0.21132487, 0.57735027],
                [0.57735027, 0.57735027, -0.57735027],
                [-0.21132487, 0.78867513, 0.57735027],
            ],
        )

    def test_coherent_interface_builder(self):
        builder = CoherentInterfaceBuilder(
            film_structure=self.sio2_conventional,
            substrate_structure=self.si_conventional,
            film_miller=(1, 0, 0),
            substrate_miller=(1, 1, 1),
        )

        assert len(builder.terminations) == 2
        # SP: I am commenting out this test which is super fragile and the result fluctuates between 6 and 30 for
        # no apparent reason. The author should fix this.
        # assert len(list(builder.get_interfaces(termination=("O2_Pmmm_1", "Si_R-3m_1")))) == 30
