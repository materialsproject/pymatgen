# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import unittest

import numpy as np
import pytest

import pymatgen.io.ase as aio
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Molecule, StructureError
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.util.testing import PymatgenTest


class AseAtomsAdaptorTest(unittest.TestCase):
    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        ase_composition = Composition(atoms.get_chemical_formula())
        assert ase_composition == structure.composition
        assert atoms.cell is not None
        assert atoms.cell.any()
        assert atoms.get_pbc() is not None
        assert atoms.get_pbc().all()
        assert atoms.get_chemical_symbols() == [s.species_string for s in structure]
        assert not atoms.has("initial_magmoms")
        assert atoms.calc is None

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        prop = [3.14] * len(structure)
        structure.add_site_property("prop", prop)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_array("prop").tolist() == prop

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_mags(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        mags = [1.0] * len(structure)
        structure.add_site_property("final_magmom", mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert not atoms.has("initial_magmoms")
        assert atoms.get_magnetic_moments().tolist() == mags

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        initial_mags = [0.5] * len(structure)
        structure.add_site_property("magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_magnetic_moments().tolist() == initial_mags
        assert atoms.calc is None

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        mags = [1.0] * len(structure)
        structure.add_site_property("final_magmom", mags)
        initial_mags = [2.0] * len(structure)
        structure.add_site_property("magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_magnetic_moments().tolist(), initial_mags
        assert atoms.get_magnetic_moments().tolist(), mags

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_charge(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        charges = [1.0] * len(structure)
        structure.add_site_property("final_charge", charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert not atoms.has("initial_charges")
        assert atoms.get_charges().tolist() == charges

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        charges = [0.5] * len(structure)
        structure.add_site_property("charge", charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.calc is None
        assert atoms.get_initial_charges().tolist() == charges

        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        charges = [1.0] * len(structure)
        structure.add_site_property("final_charge", charges)
        initial_charges = [2.0] * len(structure)
        structure.add_site_property("charge", initial_charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_charges().tolist(), initial_charges
        assert atoms.get_charges().tolist(), charges

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_oxistates(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        oxi_states = [1.0] * len(structure)
        structure.add_oxidation_state_by_site(oxi_states)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_array("oxi_states").tolist() == oxi_states

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_structure_dyn(self):
        p = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structure = p.structure
        structure.add_site_property("selective_dynamics", [[False] * 3] * len(structure))
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_molecule(self):
        m = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms = aio.AseAtomsAdaptor.get_atoms(m)
        ase_composition = Composition(atoms.get_chemical_formula())
        assert ase_composition == m.composition
        assert atoms.cell is None or not atoms.cell.any()
        assert atoms.get_pbc() is None or not atoms.get_pbc().any()
        assert atoms.get_chemical_symbols() == [s.species_string for s in m]
        assert not atoms.has("initial_magmoms")
        assert atoms.calc is None

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_molecule_mags(self):
        molecule = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        mags = [1.0] * len(molecule)
        molecule.add_site_property("final_magmom", mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert not atoms.has("initial_magmoms")
        assert atoms.get_magnetic_moments().tolist() == mags

        molecule = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        initial_mags = [0.5] * len(molecule)
        molecule.add_site_property("magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert atoms.calc is None
        assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_atoms_from_molecule_dyn(self):
        molecule = Molecule.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        molecule.add_site_property("selective_dynamics", [[False] * 3] * len(molecule))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure(self):
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        assert struct.formula == "Fe4 P4 O16"
        assert [s.species_string for s in struct] == atoms.get_chemical_symbols()

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        prop = np.array([3.14] * len(atoms))
        atoms.set_array("prop", prop)
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        assert struct.site_properties["prop"] == prop.tolist()

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR_overlap"))
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        assert [s.species_string for s in struct] == atoms.get_chemical_symbols()
        with pytest.raises(StructureError):
            struct = aio.AseAtomsAdaptor.get_structure(atoms, validate_proximity=True)

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure_mag(self):
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        mags = [1.0] * len(atoms)
        atoms.set_initial_magnetic_moments(mags)
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure.site_properties["magmom"] == mags
        assert "final_magmom" not in structure.site_properties
        assert "initial_magmoms" not in structure.site_properties

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "OUTCAR"))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure.site_properties["final_magmom"] == atoms.get_magnetic_moments().tolist()
        assert "magmom" not in structure.site_properties
        assert "initial_magmoms" not in structure.site_properties

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_structure_dyn(self):
        from ase.constraints import FixAtoms
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure.site_properties["selective_dynamics"][-1][0] is False

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_get_molecule(self):
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        assert molecule.formula == "H2 C2"
        assert [s.species_string for s in molecule] == atoms.get_chemical_symbols()
        assert molecule.charge == 0
        assert molecule.spin_multiplicity == 1

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        initial_charges = [2.0] * len(atoms)
        initial_mags = [1.0] * len(atoms)
        atoms.set_initial_charges(initial_charges)
        atoms.set_initial_magnetic_moments(initial_mags)
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        assert molecule.charge == np.sum(initial_charges)
        assert molecule.spin_multiplicity == np.sum(initial_mags) + 1
        assert molecule.site_properties.get("charge", None) == initial_charges
        assert molecule.site_properties.get("magmom", None) == initial_mags

    @unittest.skipIf(not aio.ase_loaded, "ASE not loaded.")
    def test_back_forth(self):
        from ase.constraints import FixAtoms
        from ase.io import read

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "OUTCAR"))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        atoms.set_initial_charges([1.0] * len(atoms))
        atoms.set_initial_magnetic_moments([2.0] * len(atoms))
        atoms.set_array("prop", np.array([3.0] * len(atoms)))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(structure)
        structure_back = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure == structure_back

        atoms = read(os.path.join(PymatgenTest.TEST_FILES_DIR, "acetylene.xyz"))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        atoms.set_initial_charges([1.0] * len(atoms))
        atoms.set_initial_magnetic_moments([2.0] * len(atoms))
        atoms.set_array("prop", np.array([3.0] * len(atoms)))
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert atoms == atoms_back
        molecule_back = aio.AseAtomsAdaptor.get_molecule(atoms)
        assert molecule == molecule_back


if __name__ == "__main__":
    unittest.main()
