#!/usr/bin/python

import unittest
import os
import random

import numpy as np
from pymatgen.phonons.tensors import SQTensor
from pymatgen.util.testing import PymatgenTest

class SQTensorTest(PymatgenTest):

    def setUp(self):
        self.rand_sqtensor = SQTensor(np.random.randn(3,3))
        self.symm_sqtensor = SQTensor([[0.1,0.3,0.4],
                                      [0.3,0.5,0.2],
                                      [0.4,0.2,0.6]])
        self.non_invertible = SQTensor([[0.1,0,0],
                                        [0.2,0,0],
                                        [0,0,0]])

    def test_properties(self):

    def test_to
