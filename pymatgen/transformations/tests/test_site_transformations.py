# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Mar 15, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 15, 2012"

import unittest

import numpy as np

from pymatgen import Lattice, Structure
from pymatgen.transformations.site_transformations import \
    InsertSitesTransformation, TranslateSitesTransformation, \
    ReplaceSiteSpeciesTransformation, RemoveSitesTransformation, \
    PartialRemoveSitesTransformation

from monty.os.path import which


enumlib_present = which('multienum.x') and which('makestr.x')


class TranslateSitesTransformationTest(unittest.TestCase):

    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-",
                                          "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = TranslateSitesTransformation([0], [0.1, 0.2, 0.3])
        s = t.apply_transformation(self.struct)
        self.assertTrue(np.allclose(s[0].frac_coords, [0.1, 0.2, 0.3]))
        inv_t = t.inverse
        s = inv_t.apply_transformation(s)
        self.assertTrue(np.allclose(s[0].frac_coords, [0, 0, 0]))
        str(t)

    def test_to_from_dict(self):
        d = TranslateSitesTransformation([0], [0.1, 0.2, 0.3]).as_dict()
        t = TranslateSitesTransformation.from_dict(d)
        s = t.apply_transformation(self.struct)
        self.assertTrue(np.allclose(s[0].frac_coords, [0.1, 0.2, 0.3]))
        str(t)


class ReplaceSiteSpeciesTransformationTest(unittest.TestCase):

    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-",
                                          "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = ReplaceSiteSpeciesTransformation({0: "Na"})
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Na1 Li3 O4")
        str(t)

    def test_to_from_dict(self):
        d = ReplaceSiteSpeciesTransformation({0: "Na"}).as_dict()
        t = ReplaceSiteSpeciesTransformation.from_dict(d)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Na1 Li3 O4")


class RemoveSitesTransformationTest(unittest.TestCase):

    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-",
                                          "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = RemoveSitesTransformation(range(2))
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O4")
        str(t)

    def test_to_from_dict(self):
        d = RemoveSitesTransformation(range(2)).as_dict()
        t = RemoveSitesTransformation.from_dict(d)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O4")


class InsertSitesTransformationTest(unittest.TestCase):

    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-",
                                          "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = InsertSitesTransformation(["Fe", "Mn"], [[0.1, 0, 0],
                                                     [0.1, 0.2, 0.2]])
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li4 Mn1 Fe1 O4")
        t = InsertSitesTransformation(["Fe", "Mn"], [[0.001, 0, 0],
                                                     [0.1, 0.2, 0.2]])
        #Test validate proximity
        self.assertRaises(ValueError, t.apply_transformation, self.struct)

    def test_to_from_dict(self):
        d = InsertSitesTransformation(["Fe", "Mn"],
                                      [[0.1, 0, 0], [0.1, 0.2, 0.2]]).as_dict()
        t = InsertSitesTransformation.from_dict(d)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li4 Mn1 Fe1 O4")


class PartialRemoveSitesTransformationTest(unittest.TestCase):

    def setUp(self):
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
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-",
                                          "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation_complete(self):
        t = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_COMPLETE
        )
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")
        s = t.apply_transformation(self.struct, 12)
        self.assertEqual(len(s), 12)

    @unittest.skipIf(not enumlib_present, "enum_lib not present.")
    def test_apply_transformation_enumerate(self):
        t = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_ENUMERATE
        )
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")
        s = t.apply_transformation(self.struct, 12)
        self.assertEqual(len(s), 12)

    def test_apply_transformation_best_first(self):
        t = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_BEST_FIRST
        )
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")

    def test_apply_transformation_fast(self):
        t = PartialRemoveSitesTransformation(
            [tuple(range(4)), tuple(range(4, 8))],
            [0.5, 0.5],
            PartialRemoveSitesTransformation.ALGO_FAST
        )
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")
        t = PartialRemoveSitesTransformation(
            [tuple(range(8))], [0.5],
            PartialRemoveSitesTransformation.ALGO_FAST
        )
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")

    def test_to_from_dict(self):
        d = PartialRemoveSitesTransformation([tuple(range(4))], [0.5]).as_dict()
        t = PartialRemoveSitesTransformation.from_dict(d)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O4")

    def test_str(self):
        d = PartialRemoveSitesTransformation([tuple(range(4))], [0.5]).as_dict()
        self.assertIsNotNone(str(d))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
