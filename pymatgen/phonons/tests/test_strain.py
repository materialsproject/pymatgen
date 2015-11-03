#!/usr/bin/python

import unittest
import os
import random

import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.phonons.strain import Strain, Deformation, DeformedStructureSet
from pymatgen.util.testing import PymatgenTest
from numpy.testing import *
import warnings
import math

class DeformationTest(PymatgenTest):
    def setUp(self):
        self.norm_defo = Deformation.from_index_amount((0,0),0.02)
        self.ind_defo = Deformation.from_index_amount((0,1),0.02)
        self.non_ind_defo = Deformation([[1.0, 0.02, 0.02],
                                         [0.0, 1.0, 0.0],
                                         [0.0, 0.0, 1.0]])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], [[0, 0, 0], 
                                                           [0.75, 0.5, 0.75]])

    def test_properties(self):
        # green_lagrange_strain
        self.assertArrayAlmostEqual(self.ind_defo.green_lagrange_strain,
                                    [[0., 0.01, 0.],
                                     [0.01, 0.0002, 0.],
                                     [0., 0., 0.]])
        self.assertArrayAlmostEqual(self.non_ind_defo.green_lagrange_strain,
                                    [[0., 0.01, 0.01],
                                     [0.01, 0.0002, 0.0002],
                                     [0.01, 0.0002, 0.0002]])

    def test_check_independent(self):
        self.assertRaises(ValueError,self.non_ind_defo.check_independent)
        self.assertEqual(self.ind_defo.check_independent(),(0,1))

    def test_apply_to_structure(self):
        strained_norm = self.norm_defo.apply_to_structure(self.structure)
        strained_ind = self.ind_defo.apply_to_structure(self.structure)
        strained_non = self.non_ind_defo.apply_to_structure(self.structure)
        # Check lattices
        self.assertArrayAlmostEqual(strained_norm.lattice.matrix,
                                    [[3.9170018886, 0, 0],
                                     [1.958500946136, 3.32571019, 0],
                                     [0, -2.21713849,3.13550906]])
        self.assertArrayAlmostEqual(strained_ind.lattice.matrix,
                                    [[3.84019793, 0.07680396, 0],
                                     [1.92009897, 3.36411217, 0],
                                     [0, -2.21713849,3.13550906]])
        self.assertArrayAlmostEqual(strained_non.lattice.matrix,
                                    [[3.84019793, 0.07680396, 0.07680396],
                                     [1.92009897, 3.36411217, 0.0384019794],
                                     [0, -2.21713849,3.13550906]])
        # Check coordinates
        self.assertArrayAlmostEqual(strained_norm.sites[1].coords,
                                    [3.91700189, 1.224e-06, 2.3516318])
        self.assertArrayAlmostEqual(strained_ind.sites[1].coords,
                                    [3.84019793, 0.07680518, 2.3516318])
        self.assertArrayAlmostEqual(strained_non.sites[1].coords,
                                    [3.84019793, 0.07680518, 2.42843575])

class StrainTest(PymatgenTest):
    def setUp(self):
        self.norm_str = Strain.from_deformation([[1.02, 0, 0],
                                                 [0, 1, 0],
                                                 [0, 0, 1]])
        self.ind_str = Strain.from_deformation([[1, 0.02, 0],
                                                [0, 1, 0],
                                                [0, 0, 1]])
        self.non_ind_str = Strain.from_deformation([[1, 0.02, 0.02],
                                                    [0, 1, 0],
                                                    [0, 0, 1]])

    def test_new(self):
        # test warning message for constructing Strain without defo. matrix
        with warnings.catch_warnings(record=True) as w:
            self.no_dfm = Strain([[0., 0.01, 0.],
                                  [0.01, 0.0002, 0.],
                                  [0., 0., 0.]])
            assert len(w) == 1

    def test_from_deformation(self):
        self.assertArrayAlmostEqual(self.norm_str,
                                    [[0.0202, 0, 0],
                                     [0, 0, 0],
                                     [0, 0, 0]])
        self.assertArrayAlmostEqual(self.ind_str,
                                    [[0., 0.01, 0.],
                                     [0.01, 0.0002, 0.],
                                     [0., 0., 0.]])
        self.assertArrayAlmostEqual(self.non_ind_str,
                                    [[0., 0.01, 0.01],
                                     [0.01, 0.0002, 0.0002],
                                     [0.01, 0.0002, 0.0002]])

    def test_properties(self):
        pass

class DeformedStructureSetTest(PymatgenTest):
    def setUp(self):
        pass

    def test_init(self):
        pass

    def test_properties(self):
        pass

if __name__ == '__main__':
    unittest.main()
