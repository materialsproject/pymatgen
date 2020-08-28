# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import random
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
from pymatgen.analysis.molecule_matcher import KabschMatcher, PermInvMatcher


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules", "molecule_matcher")

obalign_missing = (ob is None) or ('OBAlign' not in dir(ob))


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

        mm = KabschMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        _, _, rmsd = mm.match(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0028172956033732936, places=5)

    def test_to_and_from_dict(self):

        mm_source = KabschMatcher()
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

        mm = KabschMatcher()
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0., places=5)

    def test_mismatched_atom_composition(self):

        mm = KabschMatcher()

        with self.assertRaises(ValueError):

            mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

    def test_missmatched_atom_order(self):

        mm = KabschMatcher()

        with self.assertRaises(ValueError):
            mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

        with self.assertRaises(ValueError):
            mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

        with self.assertRaises(ValueError):
            mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

        with self.assertRaises(ValueError):

            mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

        with self.assertRaises(ValueError):

            mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

        with self.assertRaises(ValueError):

            mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

    def test_fit(self):

        mm = KabschMatcher()

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0028172956033732936, places=7)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.)

        mm = KabschMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0028172956033732936, places=7)

    @unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
    def test_strange_inchi(self):

        mm = KabschMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "k1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "k2.sdf"))

        _, _, rmsd = mm.match(mol1, mol2)
        self.assertTrue(rmsd < 0.05)

    @unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
    def test_thiane(self):
        mm = KabschMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane2.sdf"))

        _, _, rmsd = mm.match(mol1, mol2)
        self.assertTrue(rmsd < 0.8)
        self.assertFalse(rmsd < 0.05)

    @unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
    def test_thiane_ethynyl(self):
        mm = KabschMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl2.sdf"))

        _, _, rmsd = mm.match(mol1, mol2)
        self.assertTrue(rmsd < 0.5)
        self.assertFalse(rmsd < 0.05)

    def test_cdi_23(self):
        mm = KabschMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "cdi_23_1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cdi_23_2.xyz"))

        _, _, rmsd = mm.match(mol1, mol2)
        self.assertTrue(rmsd < 0.3)
        self.assertFalse(rmsd < 0.05)


class PermInvMatcherTest(unittest.TestCase):

    def test_get_rmsd(self):

        mm = PermInvMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))

        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.002825344731118855, places=5)

    def test_to_and_from_dict(self):

        mm_source = PermInvMatcher()
        d_source = mm_source.as_dict()

        mm_target = PermInvMatcher.from_dict(d_source)
        self.assertDictEqual(d_source, mm_target.as_dict())

        # equal is dangerous when you have float numbers in the dictionary
        pass

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

        mm = PermInvMatcher()
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0., places=5)

    def test_mismatched_atom_composition(self):

        mm = PermInvMatcher()

        with self.assertRaises(ValueError):

            mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
            mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
            _, rmsd = mm.fit(mol1, mol2)

    def test_fit(self):

        mm = PermInvMatcher()

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 1.4171601659148593e-05, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 9.479012116064961e-05, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.002825344731118855, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 9.28245597473488e-05, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.00021150729609276233, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0001445787263551832, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.00012447269440740117, places=5)

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0., places=5)

        mm = PermInvMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.002825344731118855, places=5)

    @unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
    def test_strange_inchi(self):

        mm = PermInvMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "k1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "k2.sdf"))

        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=5)

    @unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
    def test_thiane(self):

        mm = PermInvMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane2.sdf"))

        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.9295516991456388, places=5)

    @unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
    def test_thiane_ethynyl(self):

        mm = PermInvMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl2.sdf"))

        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.45282817632344835, places=5)

    def test_cdi_23(self):

        mm = PermInvMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "cdi_23_1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cdi_23_2.xyz"))

        _, rmsd = mm.fit(mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.3011042416249653, places=5)


class KabschMatcherSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = KabschMatcher()

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
            _, rmsd = self.mm.fit(self.mol1, mol2)

    def test_rotated_molecule(self):

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)

        mol2 = self.mol1.copy()
        for site in mol2:
            site.coords = op.operate(site.coords)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0., places=5)

    def test_perturbed_atom_position(self):

        np.random.seed(42)

        mol2 = self.mol1.copy()
        mol2.perturb(0.3)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.1648, places=3)

    def test_perturbed_atoms_order(self):
        # This test shows very poor rmsd result, because the `KabschMatcher` 
        # is not capable to handle arbitrary atom's order    

        random.seed(42)

        mol2 = self.mol1.copy()
        random.shuffle(mol2)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertNotAlmostEqual(rmsd, 0.0, places=5)


class PermInvMatcherSiTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)

        struct = Structure(lattice, ['Si', 'Si'], coords)
        struct.make_supercell([2, 2, 2])

        # Creating molecule for testing
        cls.mol1 = Molecule.from_sites(struct)
        cls.mm = PermInvMatcher()

    def test_to_and_from_dict(self):

        d = self.mm.as_dict()

        mm = PermInvMatcher.from_dict(d)
        self.assertDictEqual(d, mm.as_dict())

    def test_missmatched_atoms(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        mol2 = Molecule(["C", "H", "H", "H", "H"], coords)

        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(self.mol1, mol2)

    def test_rotated_molecule(self):

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)

        mol2 = self.mol1.copy()
        for site in mol2:
            site.coords = op.operate(site.coords)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0., places=5)

    def test_perturbed_atom_position(self):

        np.random.seed(42)

        mol2 = self.mol1.copy()
        mol2.perturb(0.3)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        print(rmsd)
        self.assertAlmostEqual(rmsd, 0.1648, places=3)

    def test_perturbed_atoms_order(self):

        random.seed(42)

        mol2 = self.mol1.copy()
        random.shuffle(mol2)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=5)


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
        cls.mm = KabschMatcher()

    def test_perturbed_atom_position(self):

        np.random.seed(42)

        mol2 = self.mol1.copy()
        mol2.perturb(0.3)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.17232715121269107, places=5)

    def test_perturbed_atoms_order(self):
        # This task should fail, because `KabschMatcher` is not capable 
        # to handle arbitrary atom's order    
        random.seed(42)

        mol2 = self.mol1.copy()
        random.shuffle(mol2)

        with self.assertRaises(ValueError):
            _, rmsd = self.mm.fit(self.mol1, mol2)
            self.assertAlmostEqual(rmsd, 0.0, places=5)


class PermInvMatcherSiO2Test(unittest.TestCase):

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
        cls.mm = PermInvMatcher()

    def test_perturbed_atom_position(self):

        np.random.seed(42)

        mol2 = self.mol1.copy()
        mol2.perturb(0.3)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.17236663588463758, places=5)

    def test_perturbed_atoms_order(self):

        random.seed(42)

        mol2 = self.mol1.copy()
        random.shuffle(mol2)

        _, rmsd = self.mm.fit(self.mol1, mol2)
        self.assertAlmostEqual(rmsd, 0.0, places=5)


if __name__ == '__main__':
    unittest.main()
