# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Test for the piezo tensor class
"""

__author__ = "Handong Ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__status__ = "Development"
__date__ = "4/23/19"

import os
import unittest
import numpy as np
from pymatgen.analysis.piezo import PiezoTensor
from pymatgen.analysis.piezo_sensitivity import *
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry import site_symmetries as ss
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy.file_IO import write_disp_yaml
except ImportError:
    Phonopy = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', 'piezo_sensitivity')

class PiezoSensitivityTest(PymatgenTest):
	def setUp(self):
		self.piezo_struc = self.get_structure(os.path.join('Pb2TiZrO6'))
		self.IST = np.load(os.path.join(test_dir, "pztist.npy"))
		self.BEC = np.load(os.path.join(test_dir, "pztborn.npy"))
		self.FCM = np.load(os.path.join(test_dir, "pztfcm.npy"))
	def test_BornEffectiveChargeTensor(self):
		bt = BornEffectiveChargeTensor(self.BEC)
		self.assertArrayAlmostEqual(bt, self.BEC)
		bad_dim_array = np.zeros((3, 3))
		self.assertRaises(ValueError, BornEffectiveChargeTensor, bad_dim_array)
	# def test_get_rand_BEC(self):
	# 	btrand = BornEffectiveChargeTensor.get_rand_BEC()
	def test_InternalStrainTensor(self):
		ist = InternalStrainTensor(self.IST)
		self.assertArrayAlmostEqual(ist, self.IST)
		bad_dim_array = np.zeros((3, 3))
		self.assertRaises(ValueError, InternalStrainTensor, bad_dim_array)
	def test_get_rand_IST(self):
		def test_ForceConstantMatrix(self):
			fcmt = ForceConstantMatrix(self.FCM)
			self.assertArrayAlmostEqual(fcmt, self.FCM)
			bad_dim_array = np.zeros((3, 3))
			self.assertRaises(ValueError, ForceConstantMatrix, bad_dim_array)
	# def test_get_rand_FCM(self):








