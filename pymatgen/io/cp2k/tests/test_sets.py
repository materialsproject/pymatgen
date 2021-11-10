# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pathlib import Path

from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.inputs import Cp2kInput
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
    lattice=[[-1.0, -10.0, -100.0], [0.1, 0.01, 0.001], [7.0, 11.0, 21.0]],
    species=["X"],
    coords=[[-1, -1, -1]],
)

molecule = Molecule(species=["C", "H"], coords=[[0, 0, 0], [1, 1, 1]])


class SetTest(PymatgenTest):
    def setUp(self):
        pass

    def test_all_sets(self):
        for s in [Si_structure, molecule]:
            cis = Cp2kInputSet(s)
            self.assertMSONable(cis)
            cis = Cp2kInputSet.from_dict(cis.as_dict())
            cis = Cp2kInputSet.from_string(cis.get_string())

            cis = DftSet(s)
            self.assertMSONable(cis)
            cis = DftSet.from_dict(cis.as_dict())
            cis = DftSet.from_string(cis.get_string())

            cis = StaticSet(s)
            self.assertMSONable(cis)
            cis = StaticSet.from_dict(cis.as_dict())
            cis = StaticSet.from_string(cis.get_string())

            cis = HybridStaticSet(s)
            self.assertMSONable(cis)
            cis = HybridStaticSet.from_dict(cis.as_dict())
            cis = HybridStaticSet.from_string(cis.get_string())

            cis = RelaxSet(s)
            self.assertMSONable(cis)
            cis = RelaxSet.from_dict(cis.as_dict())
            cis = RelaxSet.from_string(cis.get_string())

            cis = HybridRelaxSet(s)
            self.assertMSONable(cis)
            cis = HybridRelaxSet.from_dict(cis.as_dict())
            cis = HybridRelaxSet.from_string(cis.get_string())

            cis = CellOptSet(s)
            self.assertMSONable(cis)
            cis = CellOptSet.from_dict(cis.as_dict())
            cis = CellOptSet.from_string(cis.get_string())

            cis = HybridCellOptSet(s)
            self.assertMSONable(cis)
            cis = HybridCellOptSet.from_dict(cis.as_dict())
            cis = HybridCellOptSet.from_string(cis.get_string())

    def test_aux_basis(self):
        Si_aux_bases = ["FIT", "cFIT", "pFIT", "cpFIT"]
        for s in Si_aux_bases:
            cis = HybridStaticSet(Si_structure, aux_basis={"Si": s})
            cis = Cp2kInput.from_string(cis.get_string())

    def test_prints(self):
        cis = RelaxSet(Si_structure, print_ldos=False, print_pdos=False)
        self.assertFalse(cis.check("FORCE_EVAL/DFT/PRINT/PRINT/PDOS"))
        cis = RelaxSet(Si_structure, print_ldos=True, print_hartree_potential=True)
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/PDOS/LDOS 1"))
        self.assertTrue(cis.check("FORCE_EVAL/DFT/PRINT/V_HARTREE_CUBE"))


if __name__ == "__main__":
    unittest.main()
