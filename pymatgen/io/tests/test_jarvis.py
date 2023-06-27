from __future__ import annotations

import os
import unittest

import pymatgen.io.jarvis as jio
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.util.testing import PymatgenTest


class JarvisAtomsAdaptorTest(unittest.TestCase):
    @unittest.skipIf(not jio.Atoms, "JARVIS-tools not loaded.")
    def test_get_atoms_from_structure(self):
        structure = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR")).structure
        atoms = jio.JarvisAtomsAdaptor.get_atoms(structure)
        jarvis_composition = atoms.composition.reduced_formula
        assert jarvis_composition == structure.composition.reduced_formula
        assert atoms.lattice_mat is not None
        assert atoms.lattice_mat.any()

    @unittest.skipIf(not jio.Atoms, "JARVIS-tools not loaded.")
    def test_get_structure(self):
        structure = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR")).structure
        atoms = jio.JarvisAtomsAdaptor.get_atoms(structure)
        assert jio.JarvisAtomsAdaptor.get_structure(atoms).composition.reduced_formula == "FePO4"
