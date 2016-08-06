# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import os
import random
import unittest2 as unittest
import json
import six

from pymatgen.core.lattice import Lattice
from pymatgen.core import PeriodicSite
from monty.json import MontyDecoder
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.transformations.standard_transformations import *
from pymatgen.symmetry.structure import SymmetrizedStructure

'''
Created on Sep 23, 2011
'''

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Sep 23, 2011"


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class RotationTransformationsTest(unittest.TestCase):
    def setUp(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Si"] * 2, coords)

    def test_as_from_dict(self):
        t = RotationTransformation([0, 1, 0], 30, False)
        d = t.as_dict()
        self.assertEqual(type(RotationTransformation.from_dict(d)),
                         RotationTransformation)

    def test_rotation_transformation(self):
        t = RotationTransformation([0, 1, 0], 30, False)
        s2 = t.apply_transformation(self.struct)
        s1 = t.inverse.apply_transformation(s2)
        self.assertTrue((abs(s1.lattice.matrix - self.struct.lattice.matrix)
                         < 1e-8).all())


class RemoveSpeciesTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = RemoveSpeciesTransformation(["Li+"])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "O2")

        d = t.as_dict()
        self.assertEqual(type(RemoveSpeciesTransformation.from_dict(d)),
                         RemoveSpeciesTransformation)


class SubstitutionTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = SubstitutionTransformation({"Li+": "Na+", "O2-": "S2-"})
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "Na2 S2")

    def test_fractional_substitution(self):
        t = SubstitutionTransformation({"Li+": "Na+",
                                        "O2-": {"S2-": 0.5, "Se2-": 0.5}})
        # test the to and from dict on the nested dictionary
        t = SubstitutionTransformation.from_dict(t.as_dict())
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
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
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = SupercellTransformation([[2, 1, 0], [0, 2, 0], [1, 0, 2]])
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.composition.formula, "Li16 O16")

    def test_from_scaling_factors(self):
        scale_factors = [random.randint(1, 5) for i in range(3)]
        t = SupercellTransformation.from_scaling_factors(*scale_factors)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.num_sites,
                         4 * six.moves.reduce(lambda a, b: a * b,
                                              scale_factors))


class OxidationStateDecorationTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = OxidationStateDecorationTransformation({"Li": 1, "O": -2})
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li", "Li", "O", "O"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s[0].species_string, "Li+")
        self.assertEqual(s[2].species_string, "O2-")
        d = t.as_dict()
        self.assertEqual(
            type(OxidationStateDecorationTransformation.from_dict(d)),
            OxidationStateDecorationTransformation)


class AutoOxiStateDecorationTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        t = AutoOxiStateDecorationTransformation()
        s = t.apply_transformation(p.structure)
        expected_oxi = {"Li": 1, "P": 5, "O": -2, "Fe": 2}
        for site in s:
            self.assertEqual(site.specie.oxi_state,
                             expected_oxi[site.specie.symbol])

    def test_as_from_dict(self):
        t = AutoOxiStateDecorationTransformation()
        d = t.as_dict()
        t = AutoOxiStateDecorationTransformation.from_dict(d)
        self.assertEqual(t.analyzer.dist_scale_factor, 1.015)


class OxidationStateRemovalTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = OxidationStateRemovalTransformation()
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "O2-", "O2-"], coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s[0].species_string, "Li")
        self.assertEqual(s[2].species_string, "O")

        d = t.as_dict()
        self.assertEqual(type(OxidationStateRemovalTransformation.from_dict(d)),
                         OxidationStateRemovalTransformation)


class PartialRemoveSpecieTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = PartialRemoveSpecieTransformation("Li+", 1.0 / 3, 3)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "O2-"], coords)
        self.assertEqual(len(t.apply_transformation(struct, 100)), 2)

        d = t.as_dict()
        self.assertEqual(type(PartialRemoveSpecieTransformation.from_dict(d)),
                         PartialRemoveSpecieTransformation)

    def test_apply_transformation_fast(self):
        t = PartialRemoveSpecieTransformation("Li+", 0.5)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.1, 0.1, 0.1])
        coords.append([0.3, 0.75, 0.3])
        lattice = Lattice([[10, 0.00, 0.00], [0, 10, 0.00], [0.00, 0, 10]])
        struct = Structure(lattice, ["Li+"] * 6, coords)
        fast_opt_s = t.apply_transformation(struct)
        t = PartialRemoveSpecieTransformation("Li+", 0.5,
                                              PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        slow_opt_s = t.apply_transformation(struct)
        self.assertAlmostEqual(EwaldSummation(fast_opt_s).total_energy,
                               EwaldSummation(slow_opt_s).total_energy, 4)
        self.assertEqual(fast_opt_s, slow_opt_s)

    def test_apply_transformations_complete_ranking(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        t1 = OxidationStateDecorationTransformation({"Li": 1, "Fe": 2, "P": 5,
                                                     "O": -2})
        s = t1.apply_transformation(p.structure)
        t = PartialRemoveSpecieTransformation("Li+", 0.5,
                                              PartialRemoveSpecieTransformation.ALGO_COMPLETE)
        self.assertEqual(len(t.apply_transformation(s, 10)), 6)

    def test_apply_transformations_best_first(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.LiFePO4'),
                             check_for_POTCAR=False)
        t1 = OxidationStateDecorationTransformation({"Li": 1, "Fe": 2, "P": 5,
                                                     "O": -2})
        s = t1.apply_transformation(p.structure)
        t = PartialRemoveSpecieTransformation("Li+", 0.5,
                                              PartialRemoveSpecieTransformation.ALGO_BEST_FIRST)
        self.assertEqual(len(t.apply_transformation(s)), 26)


class OrderDisorderedStructureTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = OrderDisorderedStructureTransformation()
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])

        struct = Structure(lattice, [{"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                                     {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                                     {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                                     {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}],
                           coords)
        output = t.apply_transformation(struct, return_ranked_list=50)
        self.assertEqual(len(output), 12)
        self.assertIsInstance(output[0]['structure'], Structure)

        struct = Structure(lattice, [{"Si4+": 0.5}, {"Si4+": 0.5},
                                     {"P5+": 0.5, "O2-": 0.5},
                                     {"P5+": 0.5, "O2-": 0.5}],
                           coords)
        output = t.apply_transformation(struct, return_ranked_list=50)
        self.assertIsInstance(output, list)
        self.assertEqual(len(output), 4)
        self.assertEqual(t.lowest_energy_structure, output[0]['structure'])

        struct = Structure(lattice, [{"Si4+": 0.5}, {"Si4+": 0.5}, {"O2-": 0.5},
                                     {"O2-": 0.5}], coords)
        allstructs = t.apply_transformation(struct, 50)
        self.assertEqual(len(allstructs), 4)

        struct = Structure(lattice, [{"Si4+": 0.333}, {"Si4+": 0.333},
                                     {"Si4+": 0.333}, "O2-"], coords)
        allstructs = t.apply_transformation(struct, 50)
        self.assertEqual(len(allstructs), 3)

        d = t.as_dict()
        self.assertEqual(
            type(OrderDisorderedStructureTransformation.from_dict(d)),
            OrderDisorderedStructureTransformation)

    def test_symmetrized_structure(self):
        t = OrderDisorderedStructureTransformation(symmetrized_structures=True)
        c = []
        sp = []
        c.append([0.5, 0.5, 0.5])
        sp.append('Si4+')
        c.append([0.45, 0.45, 0.45])
        sp.append({"Si4+": 0.5})
        c.append([0.56, 0.56, 0.56])
        sp.append({"Si4+": 0.5})
        c.append([0.25, 0.75, 0.75])
        sp.append({"Si4+": 0.5})
        c.append([0.75, 0.25, 0.25])
        sp.append({"Si4+": 0.5})
        l = Lattice.cubic(5)
        s = Structure(l, sp, c)
        test_site = PeriodicSite("Si4+", c[2], l)
        s = SymmetrizedStructure(s, 'not_real', [0, 1, 1, 2, 2])
        output = t.apply_transformation(s)
        self.assertTrue(test_site in output.sites)

    def test_too_small_cell(self):
        t = OrderDisorderedStructureTransformation()
        coords = list()
        coords.append([0.5, 0.5, 0.5])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, [{"X4+": 0.33, "O2-": 0.33, "P5+": 0.33}],
                           coords)
        self.assertRaises(ValueError, t.apply_transformation, struct)

    def test_best_first(self):
        t = OrderDisorderedStructureTransformation(algo=2)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.75, 0.75])
        coords.append([0.5, 0.5, 0.5])
        coords.append([0.25, 0.25, 0.25])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])

        struct = Structure(lattice, [{"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                                     {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                                     {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                                     {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}],
                           coords)
        output = t.apply_transformation(struct, return_ranked_list=3)
        self.assertAlmostEqual(output[0]['energy'], -234.57813667648315, 4)


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

        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+",
                                     "O2-", "O2-", "O2-", "O2-"],
                           coords)
        s = t.apply_transformation(struct)
        self.assertEqual(len(s), 4)

        with open(os.path.join(test_dir, "TiO2_super.json")) as f:
            s = json.load(f, cls=MontyDecoder)
            prim = t.apply_transformation(s)
            self.assertEqual(prim.formula, "Ti4 O8")

        d = t.as_dict()
        self.assertEqual(type(PrimitiveCellTransformation.from_dict(d)),
                         PrimitiveCellTransformation)


class PerturbStructureTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = PerturbStructureTransformation(0.05)
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([.5, .5, .5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = [[3.8401979337, 0.00, 0.00],
                   [1.9200989668, 3.3257101909, 0.00],
                   [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+",
                                     "O2-", "O2-", "O2-", "O2-"], coords)
        transformed_s = t.apply_transformation(struct)
        for i, site in enumerate(transformed_s):
            self.assertAlmostEqual(site.distance(struct[i]), 0.05)

        d = t.as_dict()
        self.assertEqual(type(PerturbStructureTransformation.from_dict(d)),
                         PerturbStructureTransformation)


class DeformStructureTransformationTest(unittest.TestCase):
    def test_apply_transformation(self):
        t = DeformStructureTransformation([[1., 0., 0.],
                                           [0., 1., 0.05],
                                           [0., 0., 1.]])
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.375, 0.375, 0.375])
        coords.append([.5, .5, .5])
        coords.append([0.875, 0.875, 0.875])
        coords.append([0.125, 0.125, 0.125])
        coords.append([0.25, 0.25, 0.25])
        coords.append([0.625, 0.625, 0.625])
        coords.append([0.75, 0.75, 0.75])

        lattice = [[3.8401979337, 0.00, 0.00],
                   [1.9200989668, 3.3257101909, 0.00],
                   [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+",
                                     "O2-", "O2-", "O2-", "O2-"], coords)
        transformed_s = t.apply_transformation(struct)
        self.assertAlmostEqual(transformed_s.lattice.a, 3.84019793)
        self.assertAlmostEqual(transformed_s.lattice.b, 3.84379750)
        self.assertAlmostEqual(transformed_s.lattice.c, 3.75022981)

        d = t.as_dict()
        self.assertEqual(type(DeformStructureTransformation.from_dict(d)),
                         DeformStructureTransformation)


if __name__ == "__main__":
    unittest.main()
