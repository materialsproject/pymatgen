# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import numpy as np

try:
    import openbabel as ob
    from pymatgen.analysis.molecule_matcher import MoleculeMatcher
    from pymatgen.analysis.molecule_matcher import IsomorphismMolAtomMapper
    from pymatgen.analysis.molecule_matcher import InchiMolAtomMapper
except (ImportError, RuntimeError):
    ob = None

from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Lattice, Structure, Molecule
from pymatgen.analysis.molecule_matcher import KabschMatcher, BruteForceOrderMatcher, HungarianOrderMatcher, GeneticOrderMatcher


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules", "molecule_matcher")

obalign_missing = (ob is None) or ('OBAlign' not in dir(ob))


def perturb(self, scale, seed):
    """
    Performs a random perturbation of the sites in a structure.
    Args:
        scale (float): Distance in angstroms by which to perturb each site.
        rng (np.random.Generator): Random generator object.
    """
    np.random.seed(seed)
    dV = np.random.normal(scale=scale, size=(len(self), 3))
    for site, dv in zip(self.sites, dV):
        site.coords += dv


def permute(self, seed):
    """
    Performs a random permutation of the sites in a structure.
    Args:
        rng (np.random.Generator): Random generator object.
    """
    np.random.seed(seed)
    inds = np.random.permutation(len(self))
    self._sites = [self[i] for i in inds]


@unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
class MoleculeMatcherTest(unittest.TestCase):

    def test_fit(self):
        self.fit_with_mapper(IsomorphismMolAtomMapper())
        self.fit_with_mapper(InchiMolAtomMapper())

    def test_get_rmsd(self):
        mm = MoleculeMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        self.assertEqual('{0:7.3}'.format(mm.get_rmsd(mol1, mol2)), "0.00488")

    def test_group_molecules(self):
        mm = MoleculeMatcher(tolerance=0.001)
        with open(os.path.join(test_dir, "mol_list.txt")) as f:
            filename_list = [line.strip() for line in f.readlines()]
        mol_list = [Molecule.from_file(os.path.join(test_dir, f))
                    for f in filename_list]
        mol_groups = mm.group_molecules(mol_list)
        filename_groups = [[filename_list[mol_list.index(m)] for m in g]
                           for g in mol_groups]
        with open(os.path.join(test_dir, "grouped_mol_list.txt")) as f:
            grouped_text = f.read().strip()
        self.assertEqual(str(filename_groups), grouped_text)

    def test_to_and_from_dict(self):
        mm = MoleculeMatcher(tolerance=0.5,
                             mapper=InchiMolAtomMapper(angle_tolerance=50.0))
        d = mm.as_dict()
        mm2 = MoleculeMatcher.from_dict(d)
        self.assertEqual(d, mm2.as_dict())

        mm = MoleculeMatcher(tolerance=0.5, mapper=IsomorphismMolAtomMapper())
        d = mm.as_dict()
        mm2 = MoleculeMatcher.from_dict(d)
        self.assertEqual(d, mm2.as_dict())

    def fit_with_mapper(self, mapper):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)
        mm = MoleculeMatcher(mapper=mapper)
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        self.assertFalse(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mm = MoleculeMatcher(tolerance=0.001, mapper=mapper)
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        self.assertFalse(mm.fit(mol1, mol2))

    def test_strange_inchi(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "k1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "k2.sdf"))
        self.assertTrue(mm.fit(mol1, mol2))

    def test_thiane(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane2.sdf"))
        self.assertFalse(mm.fit(mol1, mol2))

    def test_thiane_ethynyl(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl2.sdf"))
        self.assertFalse(mm.fit(mol1, mol2))

    def test_cdi_23(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "cdi_23_1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cdi_23_2.xyz"))
        self.assertFalse(mm.fit(mol1, mol2))


class KabschMatcherTest(unittest.TestCase):

    def test_get_rmsd(self):

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = KabschMatcher(mol1)
        _, _, rmsd = mm.match(mol2)
        self.assertAlmostEqual(rmsd, 0.0028172956033732936, places=8)

    def test_to_and_from_dict(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))

        mm_source = KabschMatcher(mol1)
        d_source = mm_source.as_dict()

        mm_target = KabschMatcher.from_dict(d_source)
        self.assertDictEqual(d_source, mm_target.as_dict())

    def test_rotated_molecule(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]

        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)

        mm = KabschMatcher(mol1)
        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0., places=8)

    def test_mismatched_atom_composition(self):

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))

        mm = KabschMatcher(mol1)

        with self.assertRaises(ValueError):

            _, rmsd = mm.fit(mol2)

    def test_missmatched_atom_order(self):

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))

        mm = KabschMatcher(mol1)

        with self.assertRaises(ValueError):
            _, rmsd = mm.fit(mol2)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))

        mm = KabschMatcher(mol1)

        with self.assertRaises(ValueError):
            _, rmsd = mm.fit(mol2)

    def test_fit(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = KabschMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0028172956033732936, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        mm = KabschMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = KabschMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0028172956033732936, places=8)


class HungarianOrderMatcherTest(unittest.TestCase):

    def test_get_rmsd(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.002825344731118855, places=8)

    def test_to_and_from_dict(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))

        mm_source = HungarianOrderMatcher(mol1)
        d_source = mm_source.as_dict()

        mm_target = HungarianOrderMatcher.from_dict(d_source)
        self.assertDictEqual(d_source, mm_target.as_dict())

    def test_rotated_molecule(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]

        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)

        mm = HungarianOrderMatcher(mol1)
        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0., places=8)

    def test_mismatched_atom_composition(self):

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        with self.assertRaises(ValueError):
            _, rmsd = mm.fit(mol2)

    def test_fit(self):

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))

        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 1.4171601659148593e-05, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 9.479012116064961e-05, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.002825344731118855, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 9.28245597473488e-05, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.00021150729609276233, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0001445787263551832, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.00012447269440740117, places=8)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        mm = HungarianOrderMatcher(mol1)

        _, rmsd = mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0., places=8)


class KabschMatcherSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = KabschMatcher(cls.mol1)

    def test_to_and_from_dict(self):

        d = self.mm.as_dict()

        mm = KabschMatcher.from_dict(d)
        self.assertDictEqual(d, mm.as_dict())

    def test_missmatched_atoms(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        mol2 = Molecule(["C", "H", "H", "H", "H"], coords)

        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)

        mol2 = self.mol1.copy()
        for site in mol2:
            site.coords = op.operate(site.coords)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0., places=8)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.2628450748567651, places=8)

    def test_perturbed_atoms_order(self):
        # This test shows very poor rmsd result, because the `KabschMatcher`
        # is not capable to handle arbitrary atom's order

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertNotAlmostEqual(rmsd, 0.0, places=8)


class BruteForceOrderMatcherSmallSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        # struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = BruteForceOrderMatcher(cls.mol1)

    def test_to_and_from_dict(self):

        d = self.mm.as_dict()

        mm = BruteForceOrderMatcher.from_dict(d)
        self.assertDictEqual(d, mm.as_dict())

    def test_missmatched_atoms(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        mol2 = Molecule(["C", "H", "H", "H", "H"], coords)

        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)

        mol2 = self.mol1.copy()
        for site in mol2:
            site.coords = op.operate(site.coords)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0., places=8)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.04525207199851348, places=8)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=8)


class BruteForceOrderMatcherSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = BruteForceOrderMatcher(cls.mol1)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        # ValueError: The number of all possible permuataions (20922789888000) is not feasible to run this method!
        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(mol2)


class HungarianOrderMatcherSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = HungarianOrderMatcher(cls.mol1)

    def test_to_and_from_dict(self):

        d = self.mm.as_dict()

        mm = HungarianOrderMatcher.from_dict(d)
        self.assertDictEqual(d, mm.as_dict())

    def test_missmatched_atoms(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        mol2 = Molecule(["C", "H", "H", "H", "H"], coords)

        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(mol2)

    def test_rotated_molecule(self):

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)

        mol2 = self.mol1.copy()
        for site in mol2:
            site.coords = op.operate(site.coords)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0., places=8)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.26284507485676495, places=8)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=8)


class KabschMatcherSiO2Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [
            [0.625, 0.625, 0.625],
            [0.625, 0.625, 0.125],
            [0.625, 0.125, 0.625],
            [0.125, 0.625, 0.625],
            [0.500, 0.500, 0.500],
            [0.750, 0.750, 0.750]
        ]

        lattice = Lattice.from_parameters(a=6.61657069, b=6.61657069, c=6.61657069, alpha=60, beta=60, gamma=60)
        struct = Structure(lattice, ['Si', 'Si', 'Si', 'Si', 'O', 'O'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = KabschMatcher(cls.mol1)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.28115706155171244, places=8)

    def test_perturbed_atoms_order(self):
        # This task should fail, because `KabschMatcher` is not capable
        # to handle arbitrary atom's order

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(mol2)
            self.assertAlmostEqual(rmsd, 0.0, places=8)


class BruteForceOrderMatcherSmallSiO2Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [
            [0.625, 0.625, 0.625],
            [0.625, 0.625, 0.125],
            [0.625, 0.125, 0.625],
            [0.125, 0.625, 0.625],
            [0.500, 0.500, 0.500],
            [0.750, 0.750, 0.750]
        ]

        lattice = Lattice.from_parameters(a=6.61657069, b=6.61657069, c=6.61657069, alpha=60, beta=60, gamma=60)
        struct = Structure(lattice, ['Si', 'Si', 'Si', 'Si', 'O', 'O'], coords)
        # struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = BruteForceOrderMatcher(cls.mol1)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.1906105291112988, places=8)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=8)


class HungarianOrderMatcherSiO2Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [
            [0.625, 0.625, 0.625],
            [0.625, 0.625, 0.125],
            [0.625, 0.125, 0.625],
            [0.125, 0.625, 0.625],
            [0.500, 0.500, 0.500],
            [0.750, 0.750, 0.750]
        ]

        lattice = Lattice.from_parameters(a=6.61657069, b=6.61657069, c=6.61657069, alpha=60, beta=60, gamma=60)
        struct = Structure(lattice, ['Si', 'Si', 'Si', 'Si', 'O', 'O'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = HungarianOrderMatcher(cls.mol1)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.2813160547386018, places=8)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        _, rmsd = self.mm.fit(mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=8)


class GeneticOrderMatcherSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = GeneticOrderMatcher(cls.mol1, threshold=0.3)

    def test_to_and_from_dict(self):

        d = self.mm.as_dict()

        mm = GeneticOrderMatcher.from_dict(d)
        self.assertDictEqual(d, mm.as_dict())

    def test_missmatched_atoms(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        mol2 = Molecule(["C", "H", "H", "H", "H"], coords)

        with self.assertRaises(ValueError):
            res = self.mm.fit(mol2)

    def test_rotated_molecule(self):

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)

        mol2 = self.mol1.copy()
        for site in mol2:
            site.coords = op.operate(site.coords)

        res = self.mm.fit(mol2)

        self.assertEqual(len(res), 1)
        self.assertAlmostEqual(res[0][-1], 0., places=8)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        res = self.mm.fit(mol2)
        self.assertEqual(len(res), 1)
        self.assertAlmostEqual(res[0][-1], 0.2628450748567651, places=8)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        res = self.mm.fit(mol2)
        self.assertEqual(len(res), 1)
        self.assertAlmostEqual(res[0][-1], 0.0, places=8)

    def test_all(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)
        permute(mol2, seed=42)

        res = self.mm.fit(mol2)
        self.assertEqual(len(res), 1)

        self.assertAlmostEqual(res[0][-1], 0.2628450748567651, places=8)

# @unittest.skip('because')


class GeneticOrderMatcherSiO2Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [
            [0.625, 0.625, 0.625],
            [0.625, 0.625, 0.125],
            [0.625, 0.125, 0.625],
            [0.125, 0.625, 0.625],
            [0.500, 0.500, 0.500],
            [0.750, 0.750, 0.750]
        ]

        lattice = Lattice.from_parameters(a=6.61657069, b=6.61657069, c=6.61657069, alpha=60, beta=60, gamma=60)
        struct = Structure(lattice, ['Si', 'Si', 'Si', 'Si', 'O', 'O'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = GeneticOrderMatcher(cls.mol1, threshold=0.3)

    def test_perturbed_atom_position(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)

        res = self.mm.fit(mol2)
        self.assertEqual(len(res), 3)
        self.assertAlmostEqual(res[0][1], 0.28115706155171244, places=8)

    def test_perturbed_atoms_order(self):

        mol2 = self.mol1.copy()
        permute(mol2, seed=42)

        res = self.mm.fit(mol2)

        self.assertEqual(len(res), 3)
        self.assertAlmostEqual(res[0][1], 0.0, places=8)

    def test_all(self):

        mol2 = self.mol1.copy()
        perturb(mol2, 0.3, seed=42)
        permute(mol2, seed=42)

        res = self.mm.match(mol2)

        self.assertEqual(len(res), 3)
        self.assertEqual(res[0][0], 
                         [24, 8, 30, 27, 33, 12, 11, 43, 4, 46, 9, 21, 2, 26, 0, 25, 10, 6, 19, 15, 38, 45, 28, 18, 20,
                          42, 16, 32, 40, 39, 7, 36, 22, 34, 14, 5, 29, 47, 35, 13, 1, 31, 37, 17, 44, 23, 3, 41])
        self.assertAlmostEqual(res[0][-1], 0.28115706155171244, places=8)


if __name__ == '__main__':
    unittest.main()
