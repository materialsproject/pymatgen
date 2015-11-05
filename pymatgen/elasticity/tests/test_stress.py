#!/usr/bin/python

import unittest
import os
import random

import numpy as np
from pymatgen.elasticity.stress import Stress
from pymatgen.util.testing import PymatgenTest
from numpy.testing import *
import math

class StressTest(PymatgenTest):
    def setUp(self):
        self.rand_stress = Stress(np.random.randn(3,3))
        self.symm_stress = Stress([[0.1,0.3,0.4],
                                   [0.3,0.5,0.2],
                                   [0.4,0.2,0.6]])
        self.non_symm = Stress([[0.1,0.2,0.3],
                                [0.4,0.5,0.6],
                                [0.2,0.5,0.5]])

    def test_new(self):
        non_sq_matrix = [[0.1,0.2,0.1],
                         [0.1,0.2,0.3],
                         [0.1,0.2,0.3],
                         [0.1,0.1,0.1]]
        bad_matrix = [[0.1,0.2],
                      [0.2,0.3,0.4],
                      [0.2,0.3,0.5]]
        self.assertRaises(ValueError, Stress, non_sq_matrix)
        self.assertRaises(ValueError, Stress, bad_matrix)

    def test_properties(self):
        # mean_stress
        self.assertEqual(self.rand_stress.mean_stress,
                         1./3.*(self.rand_stress[0,0] +
                                self.rand_stress[1,1] +
                                self.rand_stress[2,2]))
        self.assertAlmostEqual(self.symm_stress.mean_stress,0.4)
        # deviator_stress
        self.assertArrayAlmostEqual(self.symm_stress.deviator_stress,
                                    Stress([[-0.3,0.3,0.4],
                                            [0.3,0.1,0.2],
                                            [0.4,0.2,0.2]]))
        
if __name__ == '__main__':
    unittest.main()
