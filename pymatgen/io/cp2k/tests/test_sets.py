# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

from pymatgen.core.structure import Molecule, Species, Structure
from pymatgen.io.cp2k.sets import (
    CellOptSet,
    Cp2kInputSet,
    DftSet,
    HybridCellOptSet,
    HybridRelaxSet,
    HybridStaticSet,
    RelaxSet,
    StaticSet,
)
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

            DftSet(s)
            StaticSet(s)
            HybridStaticSet(s)
            RelaxSet(s)
            HybridRelaxSet(s)
            CellOptSet(s)
            HybridCellOptSet(s)

    def test_aux_basis(self):
        Si_aux_bases = ["FIT", "cFIT", "pFIT", "cpFIT"]
        for s in Si_aux_bases:
            HybridStaticSet(Si_structure, aux_basis={"Si": s})

    def test_prints(self):
        cis = RelaxSet(Si_structure, print_ldos=False, print_pdos=False)
        self.assertFalse(cis.check("FORCE_EVAL/DFT/PRINT/PRINT/PDOS"))
        cis = RelaxSet(Si_structure, print_ldos=True, print_hartree_potential=True)
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/PDOS/LDOS 1"))
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"))


if __name__ == "__main__":
    unittest.main()
