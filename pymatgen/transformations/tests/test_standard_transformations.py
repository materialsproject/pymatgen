#!/usr/bin/env python
from __future__ import division

'''
Created on Sep 23, 2011
'''

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import os
import unittest
import random

from pymatgen.transformations.standard_transformations import *
from pymatgen.io.vaspio import Poscar

class TransformationsTest(unittest.TestCase):

    def setUp(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Si"] * 2 , coords)

    def test_identity_transformation(self):
        t = IdentityTransformation()
        self.assertEqual(self.struct, t.apply_transformation(self.struct))

    def test_to_dict(self):
        t = IdentityTransformation()
        self.assertIn("version", t.to_dict)
        self.assertIn("init_args", t.to_dict)

    def test_rotation_transformation(self):
        t = RotationTransformation([0, 1, 0], 30, False)
        s2 = t.apply_transformation(self.struct)
        s1 = t.inverse.apply_transformation(s2)
        self.assertTrue((abs(s1.lattice.matrix - self.struct.lattice.matrix) < 1e-8).all())


class RemoveSpeciesTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = RemoveSpeciesTransformation(["Li+"])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "O2")


class SubstitutionTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = SubstitutionTransformation({"Li+":"Na+", "O2-":"S2-"})
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "Na2 S2")

    def test_fractional_substitution(self):
        t = SubstitutionTransformation({"Li+":"Na+", "O2-":{"S2-":0.5, "Se2-":0.5}})
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "Na2 Se1 S1")

class SupercellTransformationTest(unittest.TestCase):

    def setUp(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = SupercellTransformation([[2, 1, 0], [0, 2, 0], [1, 0, 2]])
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.composition.formula, "Li16 O16")

    def test_from_scaling_factors(self):
        scale_factors = [random.randint(1, 5) for i in xrange(3)]
        t = SupercellTransformation.from_scaling_factors(*scale_factors)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.num_sites, 4 * reduce(lambda a, b: a * b, scale_factors))

class OxidationStateDecorationTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = OxidationStateDecorationTransformation({"Li":1, "O":-2})
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li", "Li", "O", "O"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(str(s[0].specie), "Li+")
        self.assertEqual(str(s[2].specie), "O2-")


class PartialRemoveSpecieTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = PartialRemoveSpecieTransformation("Li+", 1.0 / 3, True)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "O2-"], coords)
        self.assertEqual(len(t.apply_transformation(struct, True)), 2)

    def test_apply_transformation_fast(self):
        t = PartialRemoveSpecieTransformation("Li+", 0.5)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.1, 0.1, 0.1])
        coords.append([0.3, 0.75, 0.3])
        lattice = Lattice([[ 10, 0.00, 0.00], [0, 10, 0.00], [0.00, 0, 10]])
        struct = Structure(lattice, ["Li+"] * 6, coords)
        fast_opt_s = t.apply_transformation(struct)
        t = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        slow_opt_s = t.apply_transformation(struct)
        self.assertAlmostEqual(EwaldSummation(fast_opt_s).total_energy, EwaldSummation(slow_opt_s).total_energy, 4)
        self.assertEqual(fast_opt_s, slow_opt_s)

    def test_apply_transformations_complete_ranking(self):

        module_dir = os.path.dirname(os.path.abspath(__file__))
        p = Poscar.from_file(os.path.join(module_dir, 'POSCAR.LiFePO4'))
        t1 = OxidationStateDecorationTransformation({"Li":1, "Fe":2, "P":5, "O":-2})
        s = t1.apply_transformation(p.struct)
        t = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        self.assertEqual(len(t.apply_transformation(s, True)), 3)

    def test_apply_transformations_best_first(self):

        module_dir = os.path.dirname(os.path.abspath(__file__))
        p = Poscar.from_file(os.path.join(module_dir, 'POSCAR.LiFePO4'))
        t1 = OxidationStateDecorationTransformation({"Li":1, "Fe":2, "P":5, "O":-2})
        s = t1.apply_transformation(p.struct)
        t = PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_BEST_FIRST)
        self.assertEqual(len(t.apply_transformation(s)), 26)

    def test_to_from_dict(self):
        json_str = json.dumps(PartialRemoveSpecieTransformation("Li+", 0.5, PartialRemoveSpecieTransformation.ALGO_BEST_FIRST).to_dict)
        t = transformation_from_json(json_str)
        module_dir = os.path.dirname(os.path.abspath(__file__))
        p = Poscar.from_file(os.path.join(module_dir, 'POSCAR.LiFePO4'))
        t1 = OxidationStateDecorationTransformation({"Li":1, "Fe":2, "P":5, "O":-2})
        s = t1.apply_transformation(p.struct)
        self.assertEqual(len(t.apply_transformation(s)), 26)

class OrderDisorderedStructureTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = OrderDisorderedStructureTransformation(num_structures=50)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])

        struct = Structure(lattice, [{"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}, {"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}, {"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}, {"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}] , coords)
        output = t.apply_transformation(struct, True)
        self.assertEqual(len(output), 12)
        self.assertIsInstance(output[0]['structure'], Structure)

        struct = Structure(lattice, [{"Si4+":0.5}, {"Si4+":0.5}, {"P5+":0.5, "O2-": 0.5}, {"P5+":0.5, "O2-": 0.5}] , coords)
        output = t.apply_transformation(struct, return_ranked_list=True)
        self.assertIsInstance(output, list)
        self.assertEqual(len(output), 4)
        self.assertEqual(t.lowest_energy_structure, output[0]['structure'])

        struct = Structure(lattice, [{"Si4+":0.5}, {"Si4+":0.5}, {"O2-": 0.5}, {"O2-": 0.5}] , coords)
        allstructs = t.apply_transformation(struct, True)
        self.assertEqual(len(allstructs), 4)

        struct = Structure(lattice, [{"Si4+":0.333}, {"Si4+":0.333}, {"Si4+":0.333}, "O2-"] , coords)
        allstructs = t.apply_transformation(struct, True)
        self.assertEqual(len(allstructs), 3)

class PrimitiveCellTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = PrimitiveCellTransformation()
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([.5, .5, .5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(len(s), 4)

class ChargeBalanceTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = ChargeBalanceTransformation('Li+')
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([.5, .5, .5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)

        self.assertAlmostEqual(s.charge, 0, 5)


class TransformationJsonTest(unittest.TestCase):

    def test_from_json(self):
        self.assertIsInstance(transformation_from_json('{"name": "IdentityTransformation", "init_args": {}}'), IdentityTransformation)
        self.assertIsInstance(transformation_from_json('{"name": "RotationTransformation", "init_args": {"angle": 30, "angle_in_radians": false, "axis": [0, 1, 0]}}'), RotationTransformation)

if __name__ == "__main__":
    unittest.main()
