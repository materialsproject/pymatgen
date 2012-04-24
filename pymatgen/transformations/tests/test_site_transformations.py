#!/usr/bin/env python

'''
Created on Mar 15, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 15, 2012"

import unittest

import numpy as np
import json

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.transformations.site_transformations import TranslateSitesTransformation, ReplaceSiteSpeciesTransformation, RemoveSitesTransformation, PartialRemoveSitesTransformation
from pymatgen.transformations.standard_transformations import transformation_from_json

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

        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)


    def test_apply_transformation(self):
        t = TranslateSitesTransformation([0], [0.1, 0.2, 0.3])
        s = t.apply_transformation(self.struct)
        self.assertTrue(np.allclose(s[0].frac_coords, [0.1, 0.2, 0.3]))
        inv_t = t.inverse
        s = inv_t.apply_transformation(s)
        self.assertTrue(np.allclose(s[0].frac_coords, [0, 0, 0]))


    def test_to_from_dict(self):
        json_str = json.dumps(TranslateSitesTransformation([0], [0.1, 0.2, 0.3]).to_dict)
        t = transformation_from_json(json_str)
        s = t.apply_transformation(self.struct)
        self.assertTrue(np.allclose(s[0].frac_coords, [0.1, 0.2, 0.3]))

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

        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = ReplaceSiteSpeciesTransformation({0:"Na"})
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Na1 Li3 O4")

    def test_to_from_dict(self):
        json_str = json.dumps(ReplaceSiteSpeciesTransformation({0:"Na"}).to_dict)
        t = transformation_from_json(json_str)
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

        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = RemoveSitesTransformation(range(2))
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O4")

    def test_to_from_dict(self):
        json_str = json.dumps(RemoveSitesTransformation(range(2)).to_dict)
        t = transformation_from_json(json_str)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O4")


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

        lattice = Lattice([[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(lattice, ["Li+", "Li+", "Li+", "Li+", "O2-", "O2-", "O2-", "O2-"], coords)

    def test_apply_transformation(self):
        t = PartialRemoveSitesTransformation([tuple(range(4)), tuple(range(4, 8))], [0.5, 0.5], PartialRemoveSitesTransformation.ALGO_COMPLETE)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")
        s = t.apply_transformation(self.struct, 12)
        self.assertEqual(len(s), 12)

    def test_apply_transformation_best_first(self):
        t = PartialRemoveSitesTransformation([tuple(range(4)), tuple(range(4, 8))], [0.5, 0.5], PartialRemoveSitesTransformation.ALGO_BEST_FIRST)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")

    def test_apply_transformation_fast(self):
        t = PartialRemoveSitesTransformation([tuple(range(4)), tuple(range(4, 8))], [0.5, 0.5], PartialRemoveSitesTransformation.ALGO_FAST)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")
        t = PartialRemoveSitesTransformation([tuple(range(8))], [0.5], PartialRemoveSitesTransformation.ALGO_FAST)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O2")

    def test_to_from_dict(self):
        json_str = json.dumps(PartialRemoveSitesTransformation([tuple(range(4))], [0.5]).to_dict)
        t = transformation_from_json(json_str)
        s = t.apply_transformation(self.struct)
        self.assertEqual(s.formula, "Li2 O4")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
