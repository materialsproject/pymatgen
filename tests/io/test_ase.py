from __future__ import annotations

import numpy as np
import pytest

import pymatgen.io.ase as aio
from pymatgen.core import Composition, Lattice, Molecule, Structure
from pymatgen.core.structure import StructureError
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.util.testing import TEST_FILES_DIR


@pytest.fixture()
def poscar():
    return Poscar.from_file(TEST_FILES_DIR / "POSCAR")


class TestAseAtomsAdaptor:
    def test_get_atoms_from_structure(self, poscar):
        p = poscar
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
        assert not atoms.has("initial_charges")
        assert atoms.calc is None

        p = poscar
        structure = p.structure
        prop = [3.14] * len(structure)
        structure.add_site_property("prop", prop)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_array("prop").tolist() == prop

    def test_get_atoms_from_structure_mags(self, poscar):
        p = poscar
        structure = p.structure
        mags = [1.0] * len(structure)
        structure.add_site_property("final_magmom", mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert not atoms.has("initial_magmoms")
        assert atoms.get_magnetic_moments().tolist() == mags

        p = poscar
        structure = p.structure
        initial_mags = [0.5] * len(structure)
        structure.add_site_property("magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

        p = poscar
        structure = p.structure
        mags = [1.0] * len(structure)
        structure.add_site_property("final_magmom", mags)
        initial_mags = [2.0] * len(structure)
        structure.add_site_property("magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_magnetic_moments().tolist(), initial_mags
        assert atoms.get_magnetic_moments().tolist(), mags

    def test_get_atoms_from_structure_charge(self, poscar):
        p = poscar
        structure = p.structure
        charges = [1.0] * len(structure)
        structure.add_site_property("final_charge", charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert not atoms.has("initial_charges")
        assert atoms.get_charges().tolist() == charges

        p = poscar
        structure = p.structure
        charges = [0.5] * len(structure)
        structure.add_site_property("charge", charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_charges().tolist() == charges

        p = poscar
        structure = p.structure
        charges = [1.0] * len(structure)
        structure.add_site_property("final_charge", charges)
        initial_charges = [2.0] * len(structure)
        structure.add_site_property("charge", initial_charges)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_initial_charges().tolist(), initial_charges
        assert atoms.get_charges().tolist(), charges

    def test_get_atoms_from_structure_oxi_states(self, poscar):
        p = poscar
        structure = p.structure
        oxi_states = [1.0] * len(structure)
        structure.add_oxidation_state_by_site(oxi_states)
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.get_array("oxi_states").tolist() == oxi_states

    def test_get_atoms_from_structure_dyn(self, poscar):
        p = poscar
        structure = p.structure
        structure.add_site_property("selective_dynamics", [[False] * 3] * len(structure))
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]

    def test_get_atoms_from_molecule(self):
        m = Molecule.from_file(TEST_FILES_DIR / "acetylene.xyz")
        atoms = aio.AseAtomsAdaptor.get_atoms(m)
        ase_composition = Composition(atoms.get_chemical_formula())
        assert ase_composition == m.composition
        assert atoms.cell is None or not atoms.cell.any()
        assert atoms.get_pbc() is None or not atoms.get_pbc().any()
        assert atoms.get_chemical_symbols() == [s.species_string for s in m]
        assert not atoms.has("initial_magmoms")

    def test_get_atoms_from_molecule_mags(self):
        molecule = Molecule.from_file(TEST_FILES_DIR / "acetylene.xyz")
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        mags = [1.0] * len(molecule)
        molecule.add_site_property("final_magmom", mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert not atoms.has("initial_magmoms")
        assert atoms.get_magnetic_moments().tolist() == mags

        molecule = Molecule.from_file(TEST_FILES_DIR / "acetylene.xyz")
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        initial_mags = [0.5] * len(molecule)
        molecule.add_site_property("magmom", initial_mags)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

        molecule = Molecule.from_file(TEST_FILES_DIR / "acetylene.xyz")
        molecule.set_charge_and_spin(-2, spin_multiplicity=3)
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert atoms.calc is None
        assert atoms.get_initial_magnetic_moments().tolist() == [0] * len(molecule)
        assert atoms.charge == -2
        assert atoms.spin_multiplicity == 3

    def test_get_atoms_from_molecule_dyn(self):
        molecule = Molecule.from_file(TEST_FILES_DIR / "acetylene.xyz")
        molecule.add_site_property("selective_dynamics", [[False] * 3] * len(molecule))
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]

    def test_get_structure(self):
        from ase.io import read

        atoms = read(TEST_FILES_DIR / "POSCAR")
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        assert struct.formula == "Fe4 P4 O16"
        assert [s.species_string for s in struct] == atoms.get_chemical_symbols()

        atoms = read(TEST_FILES_DIR / "POSCAR")
        prop = np.array([3.14] * len(atoms))
        atoms.set_array("prop", prop)
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        assert struct.site_properties["prop"] == prop.tolist()

        atoms = read(TEST_FILES_DIR / "POSCAR_overlap")
        struct = aio.AseAtomsAdaptor.get_structure(atoms)
        assert [s.species_string for s in struct] == atoms.get_chemical_symbols()
        with pytest.raises(StructureError, match="Structure contains sites that are less than 0.01 Angstrom apart"):
            struct = aio.AseAtomsAdaptor.get_structure(atoms, validate_proximity=True)

    def test_get_structure_mag(self):
        from ase.io import read

        atoms = read(TEST_FILES_DIR / "POSCAR")
        mags = [1.0] * len(atoms)
        atoms.set_initial_magnetic_moments(mags)
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure.site_properties["magmom"] == mags
        assert "final_magmom" not in structure.site_properties
        assert "initial_magmoms" not in structure.site_properties

        atoms = read(TEST_FILES_DIR / "OUTCAR")
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure.site_properties["final_magmom"] == atoms.get_magnetic_moments().tolist()
        assert "magmom" not in structure.site_properties
        assert "initial_magmoms" not in structure.site_properties

    def test_get_structure_dyn(self):
        from ase.constraints import FixAtoms
        from ase.io import read

        atoms = read(TEST_FILES_DIR / "POSCAR")
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        assert structure.site_properties["selective_dynamics"][-1][0] is False

        # https://github.com/materialsproject/pymatgen/issues/3011
        for select_dyn in (
            [True, True, True],
            [False, False, False],
            np.array([True, True, True]),
            np.array([False, False, False]),
        ):
            structure = Structure(
                lattice=Lattice.cubic(5),
                species=("Fe", "O"),
                coords=((0, 0, 0), (0.5, 0.5, 0.5)),
                site_properties={"selective_dynamics": select_dyn},
            )
            structure.sites[0].selective_dynamics = select_dyn

            # mostly testing that this call doesn't raise
            ase_atoms = AseAtomsAdaptor.get_atoms(structure)

            assert len(ase_atoms) == len(structure)

    def test_get_molecule(self):
        from ase.io import read

        atoms = read(TEST_FILES_DIR / "acetylene.xyz")
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        assert molecule.formula == "H2 C2"
        assert [s.species_string for s in molecule] == atoms.get_chemical_symbols()
        assert molecule.charge == 0
        assert molecule.spin_multiplicity == 1

        atoms = read(TEST_FILES_DIR / "acetylene.xyz")
        initial_charges = [2.0] * len(atoms)
        initial_mags = [1.0] * len(atoms)
        atoms.set_initial_charges(initial_charges)
        atoms.set_initial_magnetic_moments(initial_mags)
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        assert molecule.charge == np.sum(initial_charges)
        assert molecule.spin_multiplicity == np.sum(initial_mags) + 1
        assert molecule.site_properties.get("charge") == initial_charges
        assert molecule.site_properties.get("magmom") == initial_mags

        atoms = read(TEST_FILES_DIR / "acetylene.xyz")
        atoms.spin_multiplicity = 3
        atoms.charge = 2
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        assert molecule.charge == 2
        assert molecule.spin_multiplicity == 3

    def test_back_forth(self):
        from ase.constraints import FixAtoms
        from ase.io import read

        # Atoms --> Structure --> Atoms --> Structure
        atoms = read(TEST_FILES_DIR / "OUTCAR")
        atoms.info = {"test": "hi"}
        atoms.set_tags([1] * len(atoms))
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        atoms.set_initial_charges([1.0] * len(atoms))
        atoms.set_initial_magnetic_moments([2.0] * len(atoms))
        atoms.set_array("prop", np.array([3.0] * len(atoms)))
        structure = aio.AseAtomsAdaptor.get_structure(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(structure)
        structure_back = aio.AseAtomsAdaptor.get_structure(atoms_back)
        assert structure_back == structure
        for k, v in atoms.todict().items():
            assert str(atoms_back.todict()[k]) == str(v)

        # Structure --> Atoms --> Structure --> Atoms
        structure = Structure.from_file(TEST_FILES_DIR / "POSCAR")
        structure.add_site_property("final_magmom", [1.0] * len(structure))
        structure.add_site_property("magmom", [2.0] * len(structure))
        structure.add_site_property("final_charge", [3.0] * len(structure))
        structure.add_site_property("charge", [4.0] * len(structure))
        structure.add_site_property("prop", [5.0] * len(structure))
        structure.info = {"test": "hi"}
        atoms = aio.AseAtomsAdaptor.get_atoms(structure)
        structure_back = aio.AseAtomsAdaptor.get_structure(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(structure_back)
        assert structure_back == structure
        for k, v in atoms.todict().items():
            assert str(atoms_back.todict()[k]) == str(v)

        # Atoms --> Molecule --> Atoms --> Molecule
        atoms = read(TEST_FILES_DIR / "acetylene.xyz")
        atoms.info = {"test": "hi"}
        atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
        atoms.set_initial_charges([1.0] * len(atoms))
        atoms.set_initial_magnetic_moments([2.0] * len(atoms))
        atoms.set_array("prop", np.array([3.0] * len(atoms)))
        atoms.set_tags([1] * len(atoms))
        molecule = aio.AseAtomsAdaptor.get_molecule(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(molecule)
        molecule_back = aio.AseAtomsAdaptor.get_molecule(atoms_back)
        for k, v in atoms.todict().items():
            assert str(atoms_back.todict()[k]) == str(v)
        assert molecule_back == molecule

        # Molecule --> Atoms --> Molecule --> Atoms
        molecule = Molecule.from_file(TEST_FILES_DIR / "acetylene.xyz")
        molecule.set_charge_and_spin(-2, spin_multiplicity=3)
        molecule.info = {"test": "hi"}
        atoms = aio.AseAtomsAdaptor.get_atoms(molecule)
        molecule_back = aio.AseAtomsAdaptor.get_molecule(atoms)
        atoms_back = aio.AseAtomsAdaptor.get_atoms(molecule_back)
        for k, v in atoms.todict().items():
            assert str(atoms_back.todict()[k]) == str(v)
        assert molecule_back == molecule
