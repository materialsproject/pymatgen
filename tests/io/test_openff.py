from __future__ import annotations

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
import pytest
from numpy.testing import assert_allclose

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.core import Molecule
from pymatgen.io.openff import (
    add_conformer,
    assign_partial_charges,
    create_openff_mol,
    get_atom_map,
    infer_openff_mol,
    mol_graph_from_openff_mol,
    mol_graph_to_openff_mol,
)
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/io/openff/classical_md_mols"

tk = pytest.importorskip("openff.toolkit")
pybel = pytest.importorskip("openbabel.pybel")


@pytest.fixture()
def mol_files():
    return {
        "CCO_xyz": f"{TEST_DIR}/CCO.xyz",
        "CCO_charges": f"{TEST_DIR}/CCO.npy",
        "FEC_r_xyz": f"{TEST_DIR}/FEC-r.xyz",
        "FEC_s_xyz": f"{TEST_DIR}/FEC-s.xyz",
        "FEC_charges": f"{TEST_DIR}/FEC.npy",
        "PF6_xyz": f"{TEST_DIR}/PF6.xyz",
        "PF6_charges": f"{TEST_DIR}/PF6.npy",
        "Li_charges": f"{TEST_DIR}/Li.npy",
        "Li_xyz": f"{TEST_DIR}/Li.xyz",
    }


def test_mol_graph_from_atom_bonds(mol_files):
    pytest.importorskip("openff")

    pf6_openff = tk.Molecule.from_smiles("F[P-](F)(F)(F)(F)F")

    pf6_graph = mol_graph_from_openff_mol(pf6_openff)

    assert len(pf6_graph.molecule) == 7
    assert pf6_graph.molecule.charge == -1

    em = iso.categorical_edge_match("weight", 1)

    pf6_openff2 = mol_graph_to_openff_mol(pf6_graph)
    pf6_graph2 = mol_graph_from_openff_mol(pf6_openff2)
    assert nx.is_isomorphic(pf6_graph.graph, pf6_graph2.graph, edge_match=em)


def test_mol_graph_from_openff_mol_cco():
    atom_coords = np.array(
        [
            [1.000000, 1.000000, 0.000000],
            [-0.515000, 1.000000, 0.000000],
            [-0.999000, 1.000000, 1.335000],
            [1.390000, 1.001000, -1.022000],
            [1.386000, 0.119000, 0.523000],
            [1.385000, 1.880000, 0.526000],
            [-0.907000, 0.118000, -0.516000],
            [-0.897000, 1.894000, -0.501000],
            [-0.661000, 0.198000, 1.768000],
        ]
    )

    atoms = ["C", "C", "O", "H", "H", "H", "H", "H", "H"]

    cco_openff = tk.Molecule.from_smiles("CCO")
    cco_openff.assign_partial_charges("mmff94")

    cco_mol_graph_1 = mol_graph_from_openff_mol(cco_openff)

    assert len(cco_mol_graph_1.molecule) == 9
    assert cco_mol_graph_1.molecule.charge == 0
    assert len(cco_mol_graph_1.graph.edges) == 8

    cco_pmg = Molecule(atoms, atom_coords)
    cco_mol_graph_2 = MoleculeGraph.with_local_env_strategy(cco_pmg, OpenBabelNN())

    em = iso.categorical_edge_match("weight", 1)

    assert nx.is_isomorphic(cco_mol_graph_1.graph, cco_mol_graph_2.graph, edge_match=em)


def test_mol_graph_to_openff_pf6(mol_files):
    """transform a water MoleculeGraph to a OpenFF water molecule"""
    pf6_mol = Molecule.from_file(mol_files["PF6_xyz"])
    pf6_mol.set_charge_and_spin(charge=-1)
    pf6_mol_graph = MoleculeGraph.with_edges(
        pf6_mol,
        {
            (0, 1): {"weight": 1},
            (0, 2): {"weight": 1},
            (0, 3): {"weight": 1},
            (0, 4): {"weight": 1},
            (0, 5): {"weight": 1},
            (0, 6): {"weight": 1},
        },
    )

    pf6_openff_1 = tk.Molecule.from_smiles("F[P-](F)(F)(F)(F)F")

    pf6_openff_2 = mol_graph_to_openff_mol(pf6_mol_graph)
    assert pf6_openff_1 == pf6_openff_2


def test_mol_graph_to_openff_cco(mol_files):
    cco_pmg = Molecule.from_file(mol_files["CCO_xyz"])
    cco_mol_graph = MoleculeGraph.with_local_env_strategy(cco_pmg, OpenBabelNN())

    cco_openff_1 = mol_graph_to_openff_mol(cco_mol_graph)

    cco_openff_2 = tk.Molecule.from_smiles("CCO")
    cco_openff_2.assign_partial_charges("mmff94")

    assert cco_openff_1 == cco_openff_2


def test_openff_back_and_forth():
    cco_openff = tk.Molecule.from_smiles("CC(=O)O")
    cco_openff.assign_partial_charges("mmff94")

    cco_mol_graph_1 = mol_graph_from_openff_mol(cco_openff)

    assert len(cco_mol_graph_1.molecule) == 8
    assert cco_mol_graph_1.molecule.charge == 0
    assert len(cco_mol_graph_1.graph.edges) == 7

    cco_openff_2 = mol_graph_to_openff_mol(cco_mol_graph_1)

    assert tk.Molecule.is_isomorphic_with(cco_openff, cco_openff_2, bond_order_matching=True)
    assert max(bond.bond_order for bond in cco_openff_2.bonds) == 2


@pytest.mark.parametrize(
    ("xyz_path", "smile", "map_values"),
    [
        ("CCO_xyz", "CCO", [0, 1, 2, 3, 4, 5, 6, 7, 8]),
        ("FEC_r_xyz", "O=C1OC[C@@H](F)O1", [0, 1, 2, 3, 4, 6, 7, 9, 8, 5]),
        ("FEC_s_xyz", "O=C1OC[C@H](F)O1", [0, 1, 2, 3, 4, 6, 7, 9, 8, 5]),
        ("PF6_xyz", "F[P-](F)(F)(F)(F)F", [1, 0, 2, 3, 4, 5, 6]),
    ],
)
def test_get_atom_map(xyz_path, smile, map_values, mol_files):
    mol = Molecule.from_file(mol_files[xyz_path])
    inferred_mol = infer_openff_mol(mol)
    openff_mol = tk.Molecule.from_smiles(smile)
    isomorphic, atom_map = get_atom_map(inferred_mol, openff_mol)
    assert isomorphic
    assert map_values == list(atom_map.values())


@pytest.mark.parametrize(
    ("xyz_path", "n_atoms", "n_bonds"),
    [
        ("CCO_xyz", 9, 8),
        ("FEC_r_xyz", 10, 10),
        ("FEC_s_xyz", 10, 10),
        ("PF6_xyz", 7, 6),
    ],
)
def test_infer_openff_mol(xyz_path, n_atoms, n_bonds, mol_files):
    mol = Molecule.from_file(mol_files[xyz_path])
    openff_mol = infer_openff_mol(mol)
    assert isinstance(openff_mol, tk.Molecule)
    assert openff_mol.n_atoms == n_atoms
    assert openff_mol.n_bonds == n_bonds


def test_add_conformer(mol_files):
    openff_mol = tk.Molecule.from_smiles("CCO")
    geometry = Molecule.from_file(mol_files["CCO_xyz"])
    openff_mol, atom_map = add_conformer(openff_mol, geometry)
    assert openff_mol.n_conformers == 1
    assert list(atom_map.values()) == list(range(openff_mol.n_atoms))


def test_assign_partial_charges(mol_files):
    openff_mol = tk.Molecule.from_smiles("CCO")
    geometry = Molecule.from_file(mol_files["CCO_xyz"])
    openff_mol, atom_map = add_conformer(openff_mol, geometry)
    partial_charges = np.load(mol_files["CCO_charges"])
    openff_mol = assign_partial_charges(openff_mol, atom_map, "am1bcc", partial_charges)
    assert_allclose(openff_mol.partial_charges.magnitude, partial_charges)


def test_create_openff_mol(mol_files):
    smile = "CCO"
    geometry = mol_files["CCO_xyz"]
    partial_charges = np.load(mol_files["CCO_charges"])
    openff_mol = create_openff_mol(smile, geometry, 1.0, partial_charges, "am1bcc")
    assert isinstance(openff_mol, tk.Molecule)
    assert openff_mol.n_atoms == 9
    assert openff_mol.n_bonds == 8
    assert_allclose(openff_mol.partial_charges.magnitude, partial_charges)


def test_add_conformer_no_geometry():
    openff_mol = tk.Molecule.from_smiles("CCO")
    openff_mol, atom_map = add_conformer(openff_mol, None)
    assert openff_mol.n_conformers == 1
    assert list(atom_map.values()) == list(range(openff_mol.n_atoms))


def test_assign_partial_charges_single_atom(mol_files):
    openff_mol = tk.Molecule.from_smiles("[Li+]")
    geometry = Molecule.from_file(mol_files["Li_xyz"])
    openff_mol, atom_map = add_conformer(openff_mol, geometry)
    openff_mol = assign_partial_charges(openff_mol, atom_map, "am1bcc", None)
    assert_allclose(openff_mol.partial_charges.magnitude, [1.0])


def test_create_openff_mol_no_geometry():
    smile = "CCO"
    openff_mol = create_openff_mol(smile)
    assert isinstance(openff_mol, tk.Molecule)
    assert openff_mol.n_atoms == 9
    assert openff_mol.n_bonds == 8
    assert openff_mol.n_conformers == 1


def test_create_openff_mol_geometry_path(mol_files):
    smile = "CCO"
    geometry = mol_files["CCO_xyz"]
    openff_mol = create_openff_mol(smile, geometry)
    assert isinstance(openff_mol, tk.Molecule)
    assert openff_mol.n_atoms == 9
    assert openff_mol.n_bonds == 8
    assert openff_mol.n_conformers == 1


def test_create_openff_mol_partial_charges_no_geometry():
    smile = "CCO"
    partial_charges = [-0.4, 0.2, 0.2]
    with pytest.raises(ValueError, match="geometries must be set if partial_charges is set"):
        create_openff_mol(smile, partial_charges=partial_charges)


def test_create_openff_mol_partial_charges_length_mismatch(mol_files):
    smile = "CCO"
    geometry = mol_files["CCO_xyz"]
    partial_charges = [-0.4, 0.2]
    with pytest.raises(ValueError, match="partial charges must have same length & order as geometry"):
        create_openff_mol(smile, geometry, partial_charges=partial_charges)
