from __future__ import annotations

from numpy.testing import assert_allclose

from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    from_2d_to_3d,
    get_2d_transform,
    get_rot_3d_for_2d,
)
from pymatgen.analysis.interfaces.substrate_analyzer import SubstrateAnalyzer
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import MatSciTest


class TestInterfaceBuilder(MatSciTest):
    @classmethod
    def setup_class(cls):
        si_struct = cls.get_structure("Si")
        sio2_struct = cls.get_structure("SiO2")
        cls.si_conventional = SpacegroupAnalyzer(si_struct).get_conventional_standard_structure()
        cls.sio2_conventional = SpacegroupAnalyzer(sio2_struct).get_conventional_standard_structure()

    def test_utils(self):
        assert_allclose(from_2d_to_3d([[1, 2], [3, 4]]), [[1, 2, 0], [3, 4, 0], [0, 0, 1]])
        assert_allclose(get_2d_transform([[1, 0], [0, 1]], [[1, 2], [3, 4]]), [[1, 2], [3, 4]])
        assert_allclose(
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
        # SP: this test is super fragile and the result fluctuates between 6, 30 and 42 for
        # no apparent reason. The author should fix this.
        assert len(list(builder.get_interfaces(termination=("O2_Pmmm_1", "Si_R-3m_1")))) >= 6


class TestCoherentInterfaceBuilder:
    def setup_method(self):
        # build substrate & film structure
        basis = [[0, 0, 0], [0.25, 0.25, 0.25]]
        self.substrate = Structure(Lattice.cubic(a=5.431), ["Si", "Si"], basis)
        self.film = Structure(Lattice.cubic(a=5.658), ["Ge", "Ge"], basis)

    def test_termination_searching(self):
        sub_analyzer = SubstrateAnalyzer()
        matches = list(sub_analyzer.calculate(substrate=self.substrate, film=self.film))
        cib1 = CoherentInterfaceBuilder(
            film_structure=self.film,
            substrate_structure=self.substrate,
            film_miller=matches[0].film_miller,
            substrate_miller=matches[0].substrate_miller,
            zslgen=sub_analyzer,
            termination_ftol=1e-4,
            label_index=True,
            filter_out_sym_slabs=False,
        )
        cib2 = CoherentInterfaceBuilder(
            film_structure=self.film,
            substrate_structure=self.substrate,
            film_miller=matches[0].film_miller,
            substrate_miller=matches[0].substrate_miller,
            zslgen=sub_analyzer,
            termination_ftol=(2, 0.1),
            label_index=True,
            filter_out_sym_slabs=False,
        )
        assert cib1.terminations == [
            ("1_Ge_P4/mmm_1", "1_Si_P4/mmm_1"),
            ("1_Ge_P4/mmm_1", "2_Si_P4/mmm_1"),
            ("2_Ge_P4/mmm_1", "1_Si_P4/mmm_1"),
            ("2_Ge_P4/mmm_1", "2_Si_P4/mmm_1"),
        ], "cib1 termination results wrong"

        assert cib2.terminations == [
            ("1_Ge_C2/m_2", "1_Si_P4/mmm_1"),
            ("1_Ge_C2/m_2", "2_Si_P4/mmm_1"),
        ], "cib2 termination results wrong"
