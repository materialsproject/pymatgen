# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from pymatgen.core.structure import Molecule, Species, Structure
from pymatgen.io.cp2k.sets import Cp2kInputSet, DftSet, HybridStaticSet, RelaxSet
from pymatgen.util.testing import PymatgenTest

Si_structure = Structure(
    lattice=[[0, 2.734364, 2.734364], [2.734364, 0, 2.734364], [2.734364, 2.734364, 0]],
    species=["Si", "Si"],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
)

nonsense_Structure = Structure(
    lattice=[[-1, -10, -100], [0.1, 0.01, 0.001], [7, 11, 21]], species=["X"], coords=[[-1, -1, -1]]
)

molecule = Molecule(species=["C", "H"], coords=[[0, 0, 0], [1, 1, 1]])


property_structure = Structure(
    lattice=[[10, 0, 0], [0, 10, 0], [0, 0, 10]],
    species=[
        Species("Ni", oxidation_state=4, properties={"spin": 0}),
        Species("O", oxidation_state=-2, properties={"spin": 0}),
        "Ni",
        "O",
    ],
    coords=[[0, 0, 0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [1, 1, 1]],
)


# TODO More comprehensive testing
class SetTest(PymatgenTest):
    def setUp(self):
        pass

    def test_all_sets(self):
        for s in [Si_structure, molecule]:
            cis = Cp2kInputSet(s)
            self.assertMSONable(cis)
            cis = Cp2kInputSet.from_dict(cis.as_dict())
            Cp2kInputSet.from_string(cis.get_string())

            ss = DftSet(s)

            assert not ss.check("motion")
            ss.activate_motion()
            assert ss.check("motion")

            ss.activate_hybrid()
            assert ss.check("force_eval/dft/xc/hf")
            assert ss.check("force_eval/dft/auxiliary_density_matrix_method")

    def test_aux_basis(self):
        Si_aux_bases = ["FIT", "cFIT", "pFIT", "cpFIT"]
        for s in Si_aux_bases:
            ss = HybridStaticSet(Si_structure, aux_basis={"Si": s})
            assert s in ss["force_eval"]["subsys"]["Si_1"]["BASIS_SET"][1].values[-1]

    def test_prints(self):
        cis = RelaxSet(Si_structure, print_ldos=False, print_pdos=False, print_v_hartree=False, print_e_density=False)
        self.assertFalse(cis.check("FORCE_EVAL/DFT/PRINT/PDOS"))
        cis.print_pdos()
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/PDOS"))

        self.assertFalse(cis.check("FORCE_EVAL/DFT/PRINT/PDOS/LDOS 1"))
        cis.print_ldos()
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/PDOS/LDOS 1"))

        self.assertFalse(cis.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"))
        cis.print_v_hartree()
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"))


if __name__ == "__main__":
    unittest.main()
