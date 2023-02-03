# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import unittest

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.molecule_matcher import (
    BruteForceOrderMatcher,
    GeneticOrderMatcher,
    HungarianOrderMatcher,
    InchiMolAtomMapper,
    IsomorphismMolAtomMapper,
    KabschMatcher,
    MoleculeMatcher,
)
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Lattice, Molecule, Structure
from pymatgen.util.testing import PymatgenTest

try:
    import openbabel as ob

except (ImportError, RuntimeError):
    ob = None

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "molecule_matcher")


obalign_missing = ob is None or "OBAlign" not in dir(ob)


def rotate(mol, seed):
    """Performs a random rotation of the sites in a structure.

    Args:
        mol (Molecule): The Molecule object which will be transformed.
        seed (int): The seed value for the random generator.
    """
    rng = np.random.default_rng(seed=seed)

    op = SymmOp.from_origin_axis_angle([0, 0, 0], rng.random(3), 360 * rng.random())
    for site in mol:
        site.coords = op.operate(site.coords)


def perturb(mol, scale, seed):
    """Performs a random perturbation of the sites in a structure.

    Args:
        scale (float): Distance in angstroms by which to perturb each site.
        seed (int): The seed value for the random generator.
    """
    rng = np.random.default_rng(seed=seed)

    dV = rng.normal(scale=scale, size=(len(mol), 3))
    for site, dv in zip(mol.sites, dV):
        site.coords += dv


def permute(mol, seed):
    """Performs a random permutation of the sites in a structure.

    Args:
        seed (int): The seed value for the random generator.
    """
    rng = np.random.default_rng(seed=seed)

    inds = rng.permutation(len(mol))
    mol._sites = [mol[i] for i in inds]


def generate_Si_cluster():
    from pymatgen.io.xyz import XYZ

    coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
    lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

    struct = Structure(lattice, ["Si", "Si"], coords)
    struct.make_supercell([2, 2, 2])

    # Creating molecule for testing
    mol = Molecule.from_sites(struct)
    XYZ(mol).write_file(os.path.join(test_dir, "Si_cluster.xyz"))

    # Rorate the whole molecule
    mol_rotated = mol.copy()
    rotate(mol_rotated, seed=42)
    XYZ(mol_rotated).write_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))

    # Perturbing the atom positions
    mol_perturbed = mol.copy()
    perturb(mol_perturbed, 0.3, seed=42)
    XYZ(mol_perturbed).write_file(os.path.join(test_dir, "Si_cluster_perturbed.xyz"))

    # Permuting the order of the atoms
    mol_permuted = mol.copy()
    permute(mol_permuted, seed=42)
    XYZ(mol_permuted).write_file(os.path.join(test_dir, "Si_cluster_permuted.xyz"))

    # All-in-one
    mol2 = mol.copy()
    rotate(mol2, seed=42)
    perturb(mol2, 0.3, seed=42)
    permute(mol2, seed=42)
    XYZ(mol2).write_file(os.path.join(test_dir, "Si_cluster_2.xyz"))


def generate_Si2O_cluster():
    from pymatgen.io.xyz import XYZ

    coords = [
        [0.625, 0.625, 0.625],
        [0.625, 0.625, 0.125],
        [0.625, 0.125, 0.625],
        [0.125, 0.625, 0.625],
        [0.500, 0.500, 0.500],
        [0.750, 0.750, 0.750],
    ]

    lattice = Lattice.from_parameters(a=6.61657069, b=6.61657069, c=6.61657069, alpha=60, beta=60, gamma=60)
    struct = Structure(lattice, ["Si", "Si", "Si", "Si", "O", "O"], coords)
    # struct.make_supercell([2, 2, 2])

    # Creating molecule for testing
    mol = Molecule.from_sites(struct)
    XYZ(mol).write_file(os.path.join(test_dir, "Si2O_cluster.xyz"))

    # Rorate the whole molecule
    mol_rotated = mol.copy()
    rotate(mol_rotated, seed=42)
    XYZ(mol_rotated).write_file(os.path.join(test_dir, "Si2O_cluster_rotated.xyz"))

    # Perturbing the atom positions
    mol_perturbed = mol.copy()
    perturb(mol_perturbed, 0.3, seed=42)
    XYZ(mol_perturbed).write_file(os.path.join(test_dir, "Si2O_cluster_perturbed.xyz"))

    # Permuting the order of the atoms
    mol_permuted = mol.copy()
    permute(mol_permuted, seed=42)
    XYZ(mol_permuted).write_file(os.path.join(test_dir, "Si2O_cluster_permuted.xyz"))

    # All-in-one
    mol2 = mol.copy()
    rotate(mol2, seed=42)
    perturb(mol2, 0.3, seed=42)
    permute(mol2, seed=42)
    XYZ(mol2).write_file(os.path.join(test_dir, "Si2O_cluster_2.xyz"))


@unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
class MoleculeMatcherTest(unittest.TestCase):
    def test_fit(self):
        self.fit_with_mapper(IsomorphismMolAtomMapper())
        self.fit_with_mapper(InchiMolAtomMapper())

    def test_get_rmsd(self):
        mm = MoleculeMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        assert f"{mm.get_rmsd(mol1, mol2):7.3}" == "0.00488"

    def test_group_molecules(self):
        mm = MoleculeMatcher(tolerance=0.001)
        with open(os.path.join(test_dir, "mol_list.txt")) as f:
            filename_list = [line.strip() for line in f.readlines()]
        mol_list = [Molecule.from_file(os.path.join(test_dir, f)) for f in filename_list]
        mol_groups = mm.group_molecules(mol_list)
        filename_groups = [[filename_list[mol_list.index(m)] for m in g] for g in mol_groups]
        with open(os.path.join(test_dir, "grouped_mol_list.txt")) as f:
            grouped_text = f.read().strip()
        assert str(filename_groups) == grouped_text

    def test_to_and_from_dict(self):
        mm = MoleculeMatcher(tolerance=0.5, mapper=InchiMolAtomMapper(angle_tolerance=50.0))
        d = mm.as_dict()
        mm2 = MoleculeMatcher.from_dict(d)
        assert d == mm2.as_dict()

        mm = MoleculeMatcher(tolerance=0.5, mapper=IsomorphismMolAtomMapper())
        d = mm.as_dict()
        mm2 = MoleculeMatcher.from_dict(d)
        assert d == mm2.as_dict()

    def fit_with_mapper(self, mapper):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]
        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)
        mm = MoleculeMatcher(mapper=mapper)
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        assert not mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        assert mm.fit(mol1, mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        assert mm.fit(mol1, mol2)

        mm = MoleculeMatcher(tolerance=0.001, mapper=mapper)
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        assert not mm.fit(mol1, mol2)

    def test_strange_inchi(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "k1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "k2.sdf"))
        assert mm.fit(mol1, mol2)

    def test_thiane(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane2.sdf"))
        assert not mm.fit(mol1, mol2)

    def test_thiane_ethynyl(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl2.sdf"))
        assert not mm.fit(mol1, mol2)

    def test_cdi_23(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "cdi_23_1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cdi_23_2.xyz"))
        assert not mm.fit(mol1, mol2)


class KabschMatcherTest(unittest.TestCase):
    def test_get_rmsd(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = KabschMatcher(mol1)
        _, _, rmsd = mm.match(mol2)
        assert rmsd == approx(0.0028172956033732936, abs=1e-6)

    def test_to_and_from_dict(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))

        mm_source = KabschMatcher(mol1)
        d_source = mm_source.as_dict()

        mm_target = KabschMatcher.from_dict(d_source)
        assert d_source == mm_target.as_dict()

    def test_rotated_molecule(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]

        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)

        mm = KabschMatcher(mol1)
        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_mismatched_atom_composition(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))

        mm = KabschMatcher(mol1)

        with pytest.raises(ValueError):
            _, rmsd = mm.fit(mol2)

    def test_missmatched_atom_order(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))

        mm = KabschMatcher(mol1)

        with pytest.raises(ValueError):
            _, rmsd = mm.fit(mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))

        mm = KabschMatcher(mol1)

        with pytest.raises(ValueError):
            _, rmsd = mm.fit(mol2)

    def test_fit(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = KabschMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0.0028172956033732936, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        mm = KabschMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0, abs=6)


class HungarianOrderMatcherTest(unittest.TestCase):
    def test_get_rmsd(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0.002825344731118855, abs=1e-6)

    def test_to_and_from_dict(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))

        mm_source = HungarianOrderMatcher(mol1)
        d_source = mm_source.as_dict()

        mm_target = HungarianOrderMatcher.from_dict(d_source)
        assert d_source == mm_target.as_dict()

    def test_rotated_molecule(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]

        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)

        mm = HungarianOrderMatcher(mol1)
        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_mismatched_atom_composition(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        with pytest.raises(ValueError):
            _, rmsd = mm.fit(mol2)

    def test_fit(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))

        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(1.4171601659148593e-05, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(9.479012116064961e-05, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0.002825344731118855, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(9.28245597473488e-05, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0.00021150729609276233, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0.0001445787263551832, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0.00012447269440740117, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        assert rmsd == approx(0, abs=6)


class GeneticOrderMatcherTest(unittest.TestCase):
    def test_get_rmsd(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = GeneticOrderMatcher(mol1, threshold=0.3)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0.0028172956033734615, abs=1e-6)

    def test_to_and_from_dict(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))

        mm_source = GeneticOrderMatcher(mol1, threshold=0.3)
        d_source = mm_source.as_dict()

        mm_target = GeneticOrderMatcher.from_dict(d_source)
        assert d_source == mm_target.as_dict()

    def test_rotated_molecule(self):
        coords = [
            [0.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, 1.089000],
            [1.026719, 0.000000, -0.363000],
            [-0.513360, -0.889165, -0.363000],
            [-0.513360, 0.889165, -0.363000],
        ]

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]

        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)

        mm = GeneticOrderMatcher(mol1, threshold=0.3)
        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0, abs=6)

    def test_mismatched_atom_composition(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.3)

        with pytest.raises(ValueError):
            _, rmsd = mm.fit(mol2)[0]

    def test_fit(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))

        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(7.061017534055039e-05, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(9.459575146593829e-05, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0.0028172956033734615, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(9.28245597473488e-05, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0.00019757961816426042, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.1)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0.0001398867874149986, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0.00012190586696474853, abs=1e-6)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        mm = GeneticOrderMatcher(mol1, threshold=0.01)

        _, rmsd = mm.fit(mol2)[0]
        assert rmsd == approx(0, abs=6)


class KabschMatcherSiTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si_cluster.xyz"))
        cls.mm = KabschMatcher(cls.mol1)

    def test_to_and_from_dict(self):
        d = self.mm.as_dict()
        mm = KabschMatcher.from_dict(d)
        assert d == mm.as_dict()

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster.xyz"))
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_perturbed.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.2232223954240079, abs=1e-6)

    def test_permuted_atoms_order(self):
        # This test shows very poor rmsd result, because the `KabschMatcher`
        # is not capable to handle arbitrary atom's order
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_permuted.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(2.7962454578966454, abs=1e-6)


class BruteForceOrderMatcherSiTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si_cluster.xyz"))
        cls.mm = BruteForceOrderMatcher(cls.mol1)

    def test_to_and_from_dict(self):
        d = self.mm.as_dict()
        mm = BruteForceOrderMatcher.from_dict(d)
        assert d == mm.as_dict()

    def test_random_match(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_2.xyz"))

        # ValueError: The number of all possible permuataions (20922789888000) is not feasible to run this method!
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)


class HungarianOrderMatcherSiTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si_cluster.xyz"))
        cls.mm = HungarianOrderMatcher(cls.mol1)

    def test_to_and_from_dict(self):
        d = self.mm.as_dict()
        mm = HungarianOrderMatcher.from_dict(d)
        assert d == mm.as_dict()

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_rotated.xyz"))
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):
        # TODO: Checking the cause of the large deviation
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(1.025066171481399, abs=1e-6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_perturbed.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.2232223954240077, abs=1e-6)

    def test_permuted_atoms_order(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_permuted.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_random_match(self):
        # TODO: Checking the cause of the large deviation
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_2.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(1.0177241485450828, abs=1e-6)


class GeneticOrderMatcherSiTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si_cluster.xyz"))
        cls.mm = GeneticOrderMatcher(cls.mol1, threshold=0.3)

    def test_to_and_from_dict(self):
        d = self.mm.as_dict()
        mm = GeneticOrderMatcher.from_dict(d)
        assert d == mm.as_dict()

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster.xyz"))
        with pytest.raises(ValueError):
            self.mm.fit(mol2)

    def test_rotated_molecule(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))
        res = self.mm.fit(mol2)
        assert res[0][-1] == approx(0, abs=6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_perturbed.xyz"))
        res = self.mm.fit(mol2)
        assert res[0][-1] == approx(0.2232223954240079, abs=1e-6)

    def test_permuted_atoms_order(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_permuted.xyz"))
        res = self.mm.fit(mol2)
        assert res[0][-1] == approx(0, abs=6)

    def test_random_match(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_2.xyz"))
        res = self.mm.fit(mol2)
        assert res[0][-1] == approx(0.22163169511782, abs=1e-6)


class KabschMatcherSi2OTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster.xyz"))
        cls.mm = KabschMatcher(cls.mol1)

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_rotated.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_perturbed.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.24340452336622473, abs=1e-6)

    def test_permuted_atoms_order(self):
        # This task should fail, because `KabschMatcher` is not capable
        # to handle arbitrary atom's order
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_permuted.xyz"))
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)


class BruteForceOrderMatcherSi2OTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster.xyz"))
        cls.mm = BruteForceOrderMatcher(cls.mol1)

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_rotated.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_perturbed.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.2434045087608993, abs=1e-6)

    def test_permuted_atoms_order(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_permuted.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_random_match(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_2.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.23051587697194997, abs=1e-6)


class HungarianOrderMatcherSi2OTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster.xyz"))
        cls.mm = HungarianOrderMatcher(cls.mol1)

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster_rotated.xyz"))
        with pytest.raises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_rotated.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_perturbed.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.24474957657894614, abs=1e-6)

    def test_permuted_atoms_order(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_permuted.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0, abs=6)

    def test_random_match(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_2.xyz"))
        _, rmsd = self.mm.fit(mol2)
        assert rmsd == approx(0.23231038877573124, abs=1e-6)


class GeneticOrderMatcherSi2OTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mol1 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster.xyz"))
        cls.mm = GeneticOrderMatcher(cls.mol1, threshold=0.3)

    def test_missmatched_atoms(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si_cluster.xyz"))
        with pytest.raises(ValueError):
            self.mm.fit(mol2)

    def test_rotated_molecule(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_rotated.xyz"))
        res = self.mm.fit(mol2)
        assert res[0][1] == approx(0, abs=6)

    def test_perturbed_atom_position(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_perturbed.xyz"))
        res = self.mm.fit(mol2)
        assert len(res) == 3
        assert res[0][1] == approx(0.24340452336622473, abs=1e-6)

    def test_permuted_atoms_order(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_permuted.xyz"))
        res = self.mm.fit(mol2)
        assert len(res) == 3
        assert res[0][1] == approx(0, abs=6)

    def test_random_match(self):
        mol2 = Molecule.from_file(os.path.join(test_dir, "Si2O_cluster_2.xyz"))
        res = self.mm.match(mol2)
        assert len(res) == 3
        assert res[0][0] == [5, 0, 4, 1, 3, 2]
        assert res[0][-1] == approx(0.2305159973457393, abs=1e-6)


if __name__ == "__main__":
    # Run the following code to generate test cases:
    # generate_Si_cluster()
    # generate_Si2O_cluster()

    unittest.main()
