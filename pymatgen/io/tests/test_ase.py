# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os

from pymatgen import Composition
from pymatgen.io.vasp.inputs import Poscar
import pymatgen.io.ase as aio

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class AseAtomsAdaptorTest(unittest.TestCase):

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        structure = p.structure
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        ase_composition = Composition(atoms.get_chemical_formula())
        self.assertEqual(ase_composition, structure.composition)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        atoms = aio.AseAtomsAdaptor.get_atoms(p.structure)
        self.assertEqual(aio.AseAtomsAdaptor.get_structure(atoms).formula,
                         "Fe4 P4 O16")


if __name__ == "__main__":
    unittest.main()
