from __future__ import annotations

import importlib
from importlib.metadata import PackageNotFoundError
from unittest.mock import patch

import numpy as np
import pytest
from monty.json import MontyDecoder, jsanitize

from pymatgen.core import Composition, Lattice, Molecule, Structure
from pymatgen.core.structure import StructureError
from pymatgen.io.ase import AseAtomsAdaptor, MSONAtoms
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR

ase = pytest.importorskip("ase", reason="ase not installed")

STRUCTURE = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
XYZ_STRUCTURE = f"{TEST_FILES_DIR}/io/xyz/acetylene.xyz"


def test_get_atoms_from_structure():
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    ase_composition = Composition(atoms.get_chemical_formula())
    assert ase_composition == STRUCTURE.composition
    assert atoms.cell is not None
    assert atoms.cell.any()
    assert atoms.get_pbc() is not None
    assert atoms.get_pbc().all()
    assert atoms.get_chemical_symbols() == [s.species_string for s in STRUCTURE]
    assert not atoms.has("initial_magmoms")
    assert not atoms.has("initial_charges")
    assert atoms.calc is None

    prop = [3.14] * len(STRUCTURE)
    STRUCTURE.add_site_property("prop", prop)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.get_array("prop").tolist() == prop


def test_get_atoms_from_structure_mags():
    mags = [1.0] * len(STRUCTURE)
    STRUCTURE.add_site_property("final_magmom", mags)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert not atoms.has("initial_magmoms")
    assert atoms.get_magnetic_moments().tolist() == mags

    initial_mags = [0.5] * len(STRUCTURE)
    STRUCTURE.add_site_property("magmom", initial_mags)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

    mags = [1.0] * len(STRUCTURE)
    STRUCTURE.add_site_property("final_magmom", mags)
    initial_mags = [2.0] * len(STRUCTURE)
    STRUCTURE.add_site_property("magmom", initial_mags)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.get_initial_magnetic_moments().tolist(), initial_mags
    assert atoms.get_magnetic_moments().tolist(), mags


def test_get_atoms_from_structure_charge():
    charges = [1.0] * len(STRUCTURE)
    STRUCTURE.add_site_property("final_charge", charges)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert not atoms.has("initial_charges")
    assert atoms.get_charges().tolist() == charges

    charges = [0.5] * len(STRUCTURE)
    STRUCTURE.add_site_property("charge", charges)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.get_initial_charges().tolist() == charges

    charges = [1.0] * len(STRUCTURE)
    STRUCTURE.add_site_property("final_charge", charges)
    initial_charges = [2.0] * len(STRUCTURE)
    STRUCTURE.add_site_property("charge", initial_charges)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.get_initial_charges().tolist(), initial_charges
    assert atoms.get_charges().tolist(), charges


def test_get_atoms_from_structure_oxi_states():
    oxi_states = [1.0] * len(STRUCTURE)
    STRUCTURE.add_oxidation_state_by_site(oxi_states)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.get_array("oxi_states").tolist() == oxi_states


def test_get_atoms_from_structure_dyn():
    STRUCTURE.add_site_property("selective_dynamics", [[False] * 3] * len(STRUCTURE))
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]
    STRUCTURE.add_site_property("selective_dynamics", [[True] * 3] * len(STRUCTURE))
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert len(atoms.constraints) == 0
    rng = np.random.default_rng(seed=1234)
    sel_dyn = [[rng.random() < 0.5, rng.random() < 0.5, rng.random() < 0.5] for _ in STRUCTURE]
    STRUCTURE.add_site_property("selective_dynamics", sel_dyn)
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    for c in atoms.constraints:
        assert isinstance(c, ase.constraints.FixAtoms | ase.constraints.FixCartesian)
        ase_mask = c.mask if isinstance(c, ase.constraints.FixCartesian) else [True, True, True]
        assert len(c.index) == len([mask for mask in sel_dyn if np.array_equal(mask, ~np.array(ase_mask))])


def test_get_atoms_from_structure_spacegroup():
    # Get structure with space group dictionary in properties
    space_group_info = STRUCTURE.get_space_group_info()
    STRUCTURE.properties["spacegroup"] = {"number": space_group_info[1], "setting": 1}

    # Convert to atoms and check that structure was not modified
    atoms = AseAtomsAdaptor.get_atoms(STRUCTURE)
    assert isinstance(STRUCTURE.properties["spacegroup"], dict)
    assert isinstance(atoms.info["spacegroup"], ase.spacegroup.Spacegroup)

    # Check that space group matches
    assert atoms.info["spacegroup"].no == STRUCTURE.properties["spacegroup"]["number"]
    assert atoms.info["spacegroup"].setting == STRUCTURE.properties["spacegroup"]["setting"]


def test_get_atoms_from_molecule():
    mol = Molecule.from_file(XYZ_STRUCTURE)
    atoms = AseAtomsAdaptor.get_atoms(mol)
    ase_composition = Composition(atoms.get_chemical_formula())
    assert ase_composition == mol.composition
    assert atoms.cell is None or not atoms.cell.any()
    assert atoms.get_pbc() is None or not atoms.get_pbc().any()
    assert atoms.get_chemical_symbols() == [s.species_string for s in mol]
    assert not atoms.has("initial_magmoms")


def test_get_atoms_from_molecule_mags():
    molecule = Molecule.from_file(XYZ_STRUCTURE)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    mags = [1.0] * len(molecule)
    molecule.add_site_property("final_magmom", mags)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert not atoms.has("initial_magmoms")
    assert atoms.get_magnetic_moments().tolist() == mags

    molecule = Molecule.from_file(XYZ_STRUCTURE)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    initial_mags = [0.5] * len(molecule)
    molecule.add_site_property("magmom", initial_mags)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert atoms.get_initial_magnetic_moments().tolist() == initial_mags

    molecule = Molecule.from_file(XYZ_STRUCTURE)
    molecule.set_charge_and_spin(-2, spin_multiplicity=3)
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert atoms.calc is None
    assert atoms.get_initial_magnetic_moments().tolist() == [0] * len(molecule)
    assert atoms.charge == -2
    assert atoms.spin_multiplicity == 3


def test_get_atoms_from_molecule_dyn():
    molecule = Molecule.from_file(XYZ_STRUCTURE)
    molecule.add_site_property("selective_dynamics", [[False] * 3] * len(molecule))
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    assert atoms.constraints[0].get_indices().tolist() == [atom.index for atom in atoms]


def test_get_structure():
    atoms = ase.io.read(f"{VASP_IN_DIR}/POSCAR")
    struct = AseAtomsAdaptor.get_structure(atoms)
    assert struct.formula == "Fe4 P4 O16"
    assert [s.species_string for s in struct] == atoms.get_chemical_symbols()

    atoms = ase.io.read(f"{VASP_IN_DIR}/POSCAR")
    prop = np.array([3.14] * len(atoms))
    atoms.set_array("prop", prop)
    struct = AseAtomsAdaptor.get_structure(atoms)
    assert struct.site_properties["prop"] == prop.tolist()

    atoms = ase.io.read(f"{VASP_IN_DIR}/POSCAR_overlap")
    struct = AseAtomsAdaptor.get_structure(atoms)
    assert [s.species_string for s in struct] == atoms.get_chemical_symbols()
    with pytest.raises(
        StructureError,
        match=f"sites are less than {struct.DISTANCE_TOLERANCE} Angstrom apart",
    ):
        struct = AseAtomsAdaptor.get_structure(atoms, validate_proximity=True)

    # Test invalid class type
    with pytest.raises(TypeError, match="Unsupported cls="):
        struct = AseAtomsAdaptor.get_structure(atoms, cls=Lattice)


def test_get_structure_mag():
    atoms = ase.io.read(f"{VASP_IN_DIR}/POSCAR")
    mags = [1.0] * len(atoms)
    atoms.set_initial_magnetic_moments(mags)
    structure = AseAtomsAdaptor.get_structure(atoms)
    assert structure.site_properties["magmom"] == mags
    assert "final_magmom" not in structure.site_properties
    assert "initial_magmoms" not in structure.site_properties

    atoms = ase.io.read(f"{VASP_OUT_DIR}/OUTCAR.gz")
    structure = AseAtomsAdaptor.get_structure(atoms)
    assert structure.site_properties["final_magmom"] == atoms.get_magnetic_moments().tolist()
    assert "magmom" not in structure.site_properties
    assert "initial_magmoms" not in structure.site_properties


@pytest.mark.parametrize(
    "select_dyn",
    [
        [True, True, True],
        [True, False, True],
        [False, False, False],
        np.array([True, True, True]),
        np.array([True, False, True]),
        np.array([False, False, False]),
    ],
)
def test_get_structure_dyn(select_dyn):
    atoms = ase.io.read(f"{VASP_IN_DIR}/POSCAR")
    atoms.set_constraint(ase.constraints.FixAtoms(mask=[True] * len(atoms)))
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


def test_get_structure_spacegroup():
    # set up Atoms object with spacegroup information
    a = 4.05
    atoms = ase.spacegroup.crystal("Al", [(0, 0, 0)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90])
    assert "spacegroup" in atoms.info
    assert isinstance(atoms.info["spacegroup"], ase.spacegroup.Spacegroup)

    # Test that get_structure does not mutate atoms
    structure = AseAtomsAdaptor.get_structure(atoms)
    assert isinstance(atoms.info["spacegroup"], ase.spacegroup.Spacegroup)
    assert isinstance(structure.properties["spacegroup"], dict)

    # Test that spacegroup info matches
    assert atoms.info["spacegroup"].no == structure.properties["spacegroup"]["number"]
    assert atoms.info["spacegroup"].setting == structure.properties["spacegroup"]["setting"]


def test_get_molecule():
    atoms = ase.io.read(XYZ_STRUCTURE)
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    assert molecule.formula == "H2 C2"
    assert [s.species_string for s in molecule] == atoms.get_chemical_symbols()
    assert molecule.charge == 0
    assert molecule.spin_multiplicity == 1

    atoms = ase.io.read(XYZ_STRUCTURE)
    initial_charges = [2.0] * len(atoms)
    initial_mags = [1.0] * len(atoms)
    atoms.set_initial_charges(initial_charges)
    atoms.set_initial_magnetic_moments(initial_mags)
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    assert molecule.charge == np.sum(initial_charges)
    assert molecule.spin_multiplicity == np.sum(initial_mags) + 1
    assert molecule.site_properties.get("charge") == initial_charges
    assert molecule.site_properties.get("magmom") == initial_mags

    atoms = ase.io.read(XYZ_STRUCTURE)
    atoms.spin_multiplicity = 3
    atoms.charge = 2
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    assert molecule.charge == 2
    assert molecule.spin_multiplicity == 3


@pytest.mark.parametrize("filename", ["io/vasp/outputs/OUTCAR.gz", "cif/V2O3.cif"])
def test_back_forth(filename):
    # Atoms --> Structure --> Atoms --> Structure
    atoms = ase.io.read(f"{TEST_FILES_DIR}/{filename}")
    atoms.info = {"test": "hi"}
    atoms.set_tags([1] * len(atoms))
    atoms.set_constraint(ase.constraints.FixAtoms(mask=[True] * len(atoms)))
    atoms.set_initial_charges([1.0] * len(atoms))
    atoms.set_initial_magnetic_moments([2.0] * len(atoms))
    atoms.set_array("prop", np.array([3.0] * len(atoms)))
    structure = AseAtomsAdaptor.get_structure(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(structure)
    structure_back = AseAtomsAdaptor.get_structure(atoms_back)
    assert structure_back == structure
    for key, val in atoms.todict().items():
        assert str(atoms_back.todict()[key]) == str(val)


@pytest.mark.parametrize("filename", ["io/vasp/outputs/OUTCAR.gz", "cif/V2O3.cif"])
def test_back_forth_constraints(filename):
    # Atoms --> Structure --> Atoms --> Structure
    atoms = ase.io.read(f"{TEST_FILES_DIR}/{filename}")
    mask = [True] * len(atoms)
    mask[-1] = False
    atoms.set_constraint(ase.constraints.FixAtoms(mask=mask))
    structure = AseAtomsAdaptor.get_structure(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(structure)
    structure_back = AseAtomsAdaptor.get_structure(atoms_back)
    assert structure_back == structure
    for key, val in atoms.todict().items():
        assert str(atoms_back.todict()[key]) == str(val)


def test_back_forth_v2():
    # Structure --> Atoms --> Structure --> Atoms
    structure = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
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
    for key, val in atoms.todict().items():
        assert str(atoms_back.todict()[key]) == str(val)

    # test document can be jsanitized and decoded
    dct = jsanitize(structure, strict=True, enum_values=True)
    MontyDecoder().process_decoded(dct)


def test_back_forth_v3():
    # Atoms --> Molecule --> Atoms --> Molecule
    atoms = ase.io.read(XYZ_STRUCTURE)
    atoms.info = {"test": "hi"}
    atoms.set_constraint(ase.constraints.FixAtoms(mask=[True] * len(atoms)))
    atoms.set_initial_charges([1.0] * len(atoms))
    atoms.set_initial_magnetic_moments([2.0] * len(atoms))
    atoms.set_array("prop", np.array([3.0] * len(atoms)))
    atoms.set_tags([1] * len(atoms))
    molecule = AseAtomsAdaptor.get_molecule(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(molecule)
    molecule_back = AseAtomsAdaptor.get_molecule(atoms_back)
    for key, val in atoms.todict().items():
        assert str(atoms_back.todict()[key]) == str(val)
    assert molecule_back == molecule


def test_back_forth_v4():
    # Molecule --> Atoms --> Molecule --> Atoms
    molecule = Molecule.from_file(XYZ_STRUCTURE)
    molecule.set_charge_and_spin(-2, spin_multiplicity=3)
    molecule.properties = {"test": "hi"}
    atoms = AseAtomsAdaptor.get_atoms(molecule)
    molecule_back = AseAtomsAdaptor.get_molecule(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(molecule_back)
    for key, val in atoms.todict().items():
        assert str(atoms_back.todict()[key]) == str(val)
    assert molecule_back == molecule

    # test document can be jsanitized and decoded
    dct = jsanitize(molecule, strict=True, enum_values=True)
    MontyDecoder().process_decoded(dct)


def test_back_forth_v5():
    # Structure --> Atoms --> Structure --> Atoms
    structure = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
    # 2D structure; test if PBC is preserved in all conversions
    structure.lattice.pbc = (True, True, False)
    atoms = AseAtomsAdaptor.get_atoms(structure)
    structure_back = AseAtomsAdaptor.get_structure(atoms)
    atoms_back = AseAtomsAdaptor.get_atoms(structure_back)
    assert structure_back == structure
    for key, val in atoms.todict().items():
        assert str(atoms_back.todict()[key]) == str(val)


def test_msonable_atoms():
    structure = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")

    atoms = ase.io.read(f"{VASP_OUT_DIR}/OUTCAR.gz")
    atoms_info = {"test": "hi", "structure": structure}
    atoms.info = atoms_info
    assert not isinstance(atoms, MSONAtoms)

    ref_atoms = atoms.copy()
    ref_atoms.info = {}

    msonable_atoms = MSONAtoms(atoms)
    assert atoms == msonable_atoms

    msonable_atoms_dict = msonable_atoms.as_dict()
    assert msonable_atoms_dict == {
        "@module": "pymatgen.io.ase",
        "@class": "MSONAtoms",
        "atoms_json": ase.io.jsonio.encode(ref_atoms),
        "atoms_info": jsanitize(atoms_info, strict=True),
    }

    atoms_back = MSONAtoms.from_dict(msonable_atoms_dict)
    assert atoms_back == atoms
    assert atoms_back.info == atoms.info

    atoms = AseAtomsAdaptor.get_atoms(structure, msonable=True)
    assert callable(atoms.as_dict)
    assert callable(atoms.from_dict)
    assert isinstance(atoms, MSONAtoms)

    atoms = AseAtomsAdaptor.get_atoms(structure, msonable=False)
    assert not hasattr(atoms, "as_dict")
    assert not hasattr(atoms, "from_dict")
    assert isinstance(atoms, ase.Atoms)


def test_no_ase_err():
    import pymatgen.io.ase

    with patch.dict("sys.modules", {"ase.atoms": None}):
        importlib.reload(pymatgen.io.ase)
        from pymatgen.io.ase import MSONAtoms

        with pytest.raises(PackageNotFoundError, match="AseAtomsAdaptor requires the ASE package."):
            MSONAtoms()
