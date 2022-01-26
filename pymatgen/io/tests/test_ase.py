# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
import numpy as np

import pymatgen.io.ase as aio
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Molecule
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.util.testing import PymatgenTest


class AseAtomsAdaptorTest(unittest.TestCase):
    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        ase_composition = Composition(atoms.get_chemical_formula())
        self.assertEqual(ase_composition, structure.composition)
        self.assertTrue(atoms.cell is not None and atoms.cell.any())
        self.assertTrue(atoms.get_pbc() is not None and atoms.get_pbc().all())
        self.assertEqual(atoms.get_chemical_symbols(), [s.species_string for s in structure])
        self.assertFalse(atoms.has("initial_magmoms"))
        self.assertTrue(atoms.calc is None)

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        prop = [3.14] * len(structure)
        structure.add_site_property("prop", prop)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertEqual(atoms.get_array("prop").tolist(), prop)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_mags(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        mags = [1.0] * len(structure)
        structure.add_site_property("magmom", mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertFalse(atoms.has("initial_magmoms"))
        self.assertEqual(atoms.get_magnetic_moments().tolist(), mags)

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        initial_mags = [0.5] * len(structure)
        structure.add_site_property("initial_magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertEqual(atoms.get_initial_magnetic_moments().tolist(), initial_mags)

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        mags = [1.0] * len(structure)
        structure.add_site_property("magmom", mags)
        initial_mags = [2.0] * len(structure)
        structure.add_site_property("initial_magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertTrue(atoms.get_initial_magnetic_moments().tolist(), initial_mags)
        self.assertTrue(atoms.get_magnetic_moments().tolist(), mags)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_charge(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        charges = [1.0] * len(structure)
        structure.add_site_property("charge", charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertFalse(atoms.has("initial_charges"))
        self.assertEqual(atoms.get_charges().tolist(), charges)

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        initial_charges = [0.5] * len(structure)
        structure.add_site_property("initial_charge", initial_charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertEqual(atoms.get_initial_charges().tolist(), initial_charges)

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        charges = [1.0] * len(structure)
        structure.add_site_property("charge", charges)
        initial_charges = [2.0] * len(structure)
        structure.add_site_property("initial_charge", initial_charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertTrue(atoms.get_initial_charges().tolist(), initial_charges)
        self.assertTrue(atoms.get_charges().tolist(), charges)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_oxistates(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        oxi_states = [1.0] * len(structure)
        structure.add_oxidation_state_by_site(oxi_states)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertEqual(atoms.get_array("oxi_states").tolist(), oxi_states)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_dyn(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        structure.add_site_property("selective_dynamics", [[False] * 3] * len(structure))
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        self.assertEqual(atoms.constraints[0].get_indices().tolist(), [atom.index for atom in atoms])

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_molecule(self):
        m = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms = aio.AseAtomsAdaptor.get_atoms(m)
        ase_composition = Composition(atoms.get_chemical_formula())
        self.assertEqual(ase_composition, m.composition)
        self.assertTrue(atoms.cell is None or not atoms.cell.any())
        self.assertTrue(atoms.get_pbc() is None or not atoms.get_pbc().any())
        self.assertEqual(atoms.get_chemical_symbols(), [s.species_string for s in m])
        self.assertFalse(atoms.has("initial_magmoms"))
        self.assertTrue(atoms.calc is None)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_molecule_mags(self):
        molecule = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        mags = [1.0] * len(molecule)
        molecule.add_site_property("magmom", mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        self.assertFalse(atoms.has("initial_magmoms"))
        self.assertEqual(atoms.get_magnetic_moments().tolist(), mags)

        molecule = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        initial_mags = [0.5] * len(molecule)
        molecule.add_site_property("initial_magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        self.assertEqual(atoms.get_initial_magnetic_moments().tolist(), initial_mags)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_molecule_dyn(self):

        molecule = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        molecule.add_site_property("selective_dynamics", [[False] * 3] * len(molecule))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        self.assertEqual(atoms.constraints[0].get_indices().tolist(), [atom.index for atom in atoms])

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure(self):
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        self.assertEqual(struct.formula, "Fe4 P4 O16")
        self.assertEqual([s.species_string for s in struct], atoms.get_chemical_symbols())

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        prop = np.array([3.14] * len(atoms))
        atoms.set_array("prop", prop)
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        self.assertEqual(struct.site_properties["prop"], prop.tolist())

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure_mag(self):
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        mags = [1.0] * len(atoms)
        atoms.set_initial_magnetic_moments(mags)
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        self.assertEqual(structure.site_properties["initial_magmom"], mags)

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "OUTCAR"))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        self.assertEqual(structure.site_properties["magmom"], atoms.get_magnetic_moments().tolist())

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure_dyn(self):
        from ase.io import read
        from ase.constraints import FixAtoms

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        self.assertEqual(structure.site_properties["selective_dynamics"][-1][0], False)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_molecule(self):
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        self.assertEqual(molecule.formula, "H2 C2")
        self.assertEqual([s.species_string for s in molecule], atoms.get_chemical_symbols())
        self.assertEqual(molecule.charge, 0)
        self.assertEqual(molecule.spin_multiplicity, 1)

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        initial_charges = [2.0] * len(atoms)
        initial_mags = [1.0] * len(atoms)
        atoms.set_initial_charges(initial_charges)
        atoms.set_initial_magnetic_moments(initial_mags)
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        self.assertEqual(molecule.charge, np.sum(initial_charges))
        self.assertEqual(molecule.spin_multiplicity, np.sum(initial_mags) + 1)
        self.assertEqual(molecule.site_properties.get("initial_charge", None), initial_charges)
        self.assertEqual(molecule.site_properties.get("initial_magmom", None), initial_mags)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_back_forth(self):
        from ase.io import read
        from ase.constraints import FixAtoms

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "OUTCAR"))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        atoms.set_initial_charges([1.0] * len(atoms))
        atoms.set_initial_magnetic_moments([2.0] * len(atoms))
        atoms.set_array("prop", np.array([3.0] * len(atoms)))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(structure)
        structure_back = aio.AseAtomsAdaptor.get_structure(atoms)
        self.assertEqual(structure, structure_back)

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        atoms.set_initial_charges([1.0] * len(atoms))
        atoms.set_initial_magnetic_moments([2.0] * len(atoms))
        atoms.set_array("prop", np.array([3.0] * len(atoms)))
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(molecule)
        self.assertEqual(atoms, atoms_back)
        molecule_back = aio.AseAtomsAdaptor.get_molecule(atoms)
        self.assertEqual(molecule, molecule_back)


if __name__ == "__main__":
    unittest.main()
