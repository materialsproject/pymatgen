from __future__ import annotations

import numpy as np
import pytest
from ase.constraints import FixAtoms
from ase.io import read
from monty.json import MontyDecoder, jsanitize

from pymatgen.core import Composition, Lattice, Molecule, Structure
from pymatgen.core.structure import StructureError
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.util.testing import TEST_FILES_DIR

ase = pytest.importorskip("ase")
structure = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR")


def test_get_atoms_from_structure():
    atoms = AseAtomsAdaptor.get_atoms(structure)
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

    prop = [3.14] * len(structure)
    structure.add_site_property("prop", prop)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.get_array("prop").tolist() == prop


def test_get_atoms_from_structure_mags():
    mags = [1.0] * len(structure)
    structure.add_site_property("final_magmom", mags)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert not atoms.has("initial_magmoms")
    assert atoms.get_magnetic_moments().tolist() == mags

    initial_mags = [0.5] * len(structure)
    structure.add_site_property("magmom", initial_mags)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

    mags = [1.0] * len(structure)
    structure.add_site_property("final_magmom", mags)
    initial_mags = [2.0] * len(structure)
    structure.add_site_property("magmom", initial_mags)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.get_initial_magnetic_moments().tolist(), initial_mags
    assert atoms.get_magnetic_moments().tolist(), mags


def test_get_atoms_from_structure_charge():
    charges = [1.0] * len(structure)
    structure.add_site_property("final_charge", charges)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert not atoms.has("initial_charges")
    assert atoms.get_charges().tolist() == charges

    charges = [0.5] * len(structure)
    structure.add_site_property("charge", charges)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.get_initial_charges().tolist() == charges

    charges = [1.0] * len(structure)
    structure.add_site_property("final_charge", charges)
    initial_charges = [2.0] * len(structure)
    structure.add_site_property("charge", initial_charges)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.get_initial_charges().tolist(), initial_charges
    assert atoms.get_charges().tolist(), charges


def test_get_atoms_from_structure_oxi_states():
    oxi_states = [1.0] * len(structure)
    structure.add_oxidation_state_by_site(oxi_states)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.get_array("oxi_states").tolist() == oxi_states


def test_get_atoms_from_structure_dyn():
    structure.add_site_property("selective_dynamics", [[False] * 3] * len(structure))
    atoms = AseAtomsAdaptor.get_atoms(structure)
    assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]


def test_get_atoms_from_molecule():
    mol = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")
    atoms = AseAtomsAdaptor.get_atoms(mol)
    ase_composition = Composition(atoms.get_chemical_formula())
    assert ase_composition == mol.composition
    assert atoms.cell is None or not atoms.cell.any()
    assert atoms.get_pbc() is None or not atoms.get_pbc().any()
    assert atoms.get_chemical_symbols() == [s.species_string for s in mol]
    assert not atoms.has("initial_magmoms")


def test_get_atoms_from_molecule_mags():
    molecule = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    mags = [1.0] * len(molecule)
    molecule.add_site_property("final_magmom", mags)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert not atoms.has("initial_magmoms")
    assert atoms.get_magnetic_moments().tolist() == mags

    molecule = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    initial_mags = [0.5] * len(molecule)
    molecule.add_site_property("magmom", initial_mags)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

    molecule = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")
    molecule.set_charge_and_spin(-2, spin_multiplicity=3)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert atoms.calc is None
    assert atoms.get_initial_magnetic_moments().tolist() == [0] * len(molecule)
    assert atoms.charge == -2
    assert atoms.spin_multiplicity == 3


def test_get_atoms_from_molecule_dyn():
    molecule = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")
    molecule.add_site_property("selective_dynamics", [[False] * 3] * len(molecule))
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]


def test_get_structure():
    atoms = read(f"{TEST_FILES_DIR}/POSCAR")
    struct = AseAtomsAdaptor.get_structure(atoms)
    assert struct.formula == "Fe4 P4 O16"
    assert [s.species_string for s in struct] == atoms.get_chemical_symbols()

    atoms = read(f"{TEST_FILES_DIR}/POSCAR")
    prop = np.array([3.14] * len(atoms))
    atoms.set_array("prop", prop)
    struct = AseAtomsAdaptor.get_structure(atoms)
    assert struct.site_properties["prop"] == prop.tolist()

    atoms = read(f"{TEST_FILES_DIR}/POSCAR_overlap")
    struct = AseAtomsAdaptor.get_structure(atoms)
    assert [s.species_string for s in struct] == atoms.get_chemical_symbols()
    with pytest.raises(StructureError, match="Structure contains sites that are less than 0.01 Angstrom apart"):
        struct = AseAtomsAdaptor.get_structure(atoms, validate_proximity=True)


def test_get_structure_mag():
    atoms = read(f"{TEST_FILES_DIR}/POSCAR")
    mags = [1.0] * len(atoms)
    atoms.set_initial_magnetic_moments(mags)
    structure = AseAtomsAdaptor.get_structure(atoms)
    assert structure.site_properties["magmom"] == mags
    assert "final_magmom" not in structure.site_properties
    assert "initial_magmoms" not in structure.site_properties

    atoms = read(f"{TEST_FILES_DIR}/OUTCAR")
    structure = AseAtomsAdaptor.get_structure(atoms)
    assert structure.site_properties["final_magmom"] == atoms.get_magnetic_moments().tolist()
    assert "magmom" not in structure.site_properties
    assert "initial_magmoms" not in structure.site_properties


@pytest.mark.parametrize(
    "select_dyn",
    [[True, True, True], [False, False, False], np.array([True, True, True]), np.array([False, False, False])],
)
def test_get_structure_dyn(select_dyn):
    atoms = read(f"{TEST_FILES_DIR}/POSCAR")
    atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
    structure = AseAtomsAdaptor.get_structure(atoms)
    assert structure.site_properties["selective_dynamics"][-1][0] is False

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


def test_get_molecule():
    atoms = read(f"{TEST_FILES_DIR}/acetylene.xyz")
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    assert molecule.formula == "H2 C2"
    assert [s.species_string for s in molecule] == atoms.get_chemical_symbols()
    assert molecule.charge == 0
    assert molecule.spin_multiplicity == 1

    atoms = read(f"{TEST_FILES_DIR}/acetylene.xyz")
    initial_charges = [2.0] * len(atoms)
    initial_mags = [1.0] * len(atoms)
    atoms.set_initial_charges(initial_charges)
    atoms.set_initial_magnetic_moments(initial_mags)
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    assert molecule.charge == np.sum(initial_charges)
    assert molecule.spin_multiplicity == np.sum(initial_mags) + 1
    assert molecule.site_properties.get("charge") == initial_charges
    assert molecule.site_properties.get("magmom") == initial_mags

    atoms = read(f"{TEST_FILES_DIR}/acetylene.xyz")
    atoms.spin_multiplicity = 3
    atoms.charge = 2
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    assert molecule.charge == 2
    assert molecule.spin_multiplicity == 3


@pytest.mark.parametrize("filename", ["OUTCAR", "V2O3.cif"])
def test_back_forth(filename):
    # Atoms --> Structure --> Atoms --> Structure
    atoms = read(f"{TEST_FILES_DIR}/{filename}")
    atoms.info = {"test": "hi"}
    atoms.set_tags([1] * len(atoms))
    atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
    atoms.set_initial_charges([1.0] * len(atoms))
    atoms.set_initial_magnetic_moments([2.0] * len(atoms))
    atoms.set_array("prop", np.array([3.0] * len(atoms)))
    structure = AseAtomsAdaptor.get_structure(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(structure)
    structure_back = AseAtomsAdaptor.get_structure(atoms_back)
    assert structure_back == structure
    for k, v in atoms.todict().items():
        assert str(atoms_back.todict()[k]) == str(v)


def test_back_forth_v2():
    # Structure --> Atoms --> Structure --> Atoms
    structure = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR")
    structure.add_site_property("final_magmom", [1.0] * len(structure))
    structure.add_site_property("magmom", [2.0] * len(structure))
    structure.add_site_property("final_charge", [3.0] * len(structure))
    structure.add_site_property("charge", [4.0] * len(structure))
    structure.add_site_property("prop", [5.0] * len(structure))
    structure.properties = {"test": "hi"}
    atoms = AseAtomsAdaptor.get_atoms(structure)
    structure_back = AseAtomsAdaptor.get_structure(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(structure_back)
    assert structure_back == structure
    for k, v in atoms.todict().items():
        assert str(atoms_back.todict()[k]) == str(v)

    # test document can be jsanitized and decoded
    d = jsanitize(structure, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)


def test_back_forth_v3():
    # Atoms --> Molecule --> Atoms --> Molecule
    atoms = read(f"{TEST_FILES_DIR}/acetylene.xyz")
    atoms.info = {"test": "hi"}
    atoms.set_constraint(FixAtoms(mask=[True] * len(atoms)))
    atoms.set_initial_charges([1.0] * len(atoms))
    atoms.set_initial_magnetic_moments([2.0] * len(atoms))
    atoms.set_array("prop", np.array([3.0] * len(atoms)))
    atoms.set_tags([1] * len(atoms))
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(molecule)
    molecule_back = AseAtomsAdaptor.get_molecule(atoms_back)
    for k, v in atoms.todict().items():
        assert str(atoms_back.todict()[k]) == str(v)
    assert molecule_back == molecule


def test_back_forth_v4():
    # Molecule --> Atoms --> Molecule --> Atoms
    molecule = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")
    molecule.set_charge_and_spin(-2, spin_multiplicity=3)
    molecule.properties = {"test": "hi"}
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    molecule_back = AseAtomsAdaptor.get_molecule(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(molecule_back)
    for k, v in atoms.todict().items():
        assert str(atoms_back.todict()[k]) == str(v)
    assert molecule_back == molecule

    # test document can be jsanitized and decoded
    d = jsanitize(molecule, strict=True, enum_values=True)
    MontyDecoder().process_decoded(d)
