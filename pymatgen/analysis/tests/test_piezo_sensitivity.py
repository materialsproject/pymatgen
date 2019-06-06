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

import pymatgen
from pymatgen.analysis.piezo import PiezoTensor
from pymatgen.analysis.piezo_sensitivity import *
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry import site_symmetries as ss
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

try:
	from phonopy import Phonopy
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
		self.pointops = np.load(os.path.join(test_dir, "pointops.npy"))
		self.sharedops = np.load(os.path.join(test_dir, "sharedops.npy"))
		self.BEC_IST_operations= np.load(os.path.join(test_dir, "becistops.npy"))
		self.FCM_operations= np.load(os.path.join(test_dir, "fcmops.npy"))
		self.piezo = np.array([[[  5.32649351e-03,  -1.33404642e-14,  -6.87580224e+02],
					        [ -1.33404642e-14,   4.95526253e-03,  -5.60353712e-13],
					        [ -6.87580224e+02,  -5.60353712e-13,   1.33209787e-02]],

					       [[  4.86622567e-03,   3.14840965e-13,  -7.41719319e-13],
					        [  3.14840965e-13,   5.23745666e-03,  -6.68536818e+02],
					        [ -7.41719319e-13,  -6.68536818e+02,   1.35025755e-02]],

					       [[ -1.01086427e+02,   3.20177004e-14,  -4.34079030e-14],
					        [  3.20177004e-14,  -1.01086427e+02,   1.22012318e-14],
					        [ -4.34079030e-14,   1.22012318e-14,  -5.32241086e+02]]])

	def test_BornEffectiveChargeTensor(self):
		bec = BornEffectiveCharge(self.piezo_struc, self.BEC, self.pointops)
		self.assertArrayAlmostEqual(self.BEC, bec.bec)

	def test_InternalStrainTensor(self):
		ist = InternalStrainTensor(self.piezo_struc, self.IST, self.pointops)
		self.assertArrayAlmostEqual(ist.ist, self.IST)

	def test_ForceConstantMatrix(self):
		fcmt = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		self.assertArrayAlmostEqual(fcmt.fcm, self.FCM)


	def test_get_BEC_IST_operations(self):
		bec = BornEffectiveCharge(self.piezo_struc, self.BEC, self.pointops)
		bec.get_BEC_IST_operations() 
		self.assertTrue(np.all(self.BEC_IST_operations == bec.BEC_operations))

		ist = InternalStrainTensor(self.piezo_struc, self.IST, self.pointops)
		ist.get_BEC_IST_operations(self.BEC) 
		self.assertTrue(np.all(self.BEC_IST_operations == ist.IST_operations))

	def test_get_rand_BEC(self):
		bec = BornEffectiveCharge(self.piezo_struc, self.BEC, self.pointops)
		bec.get_BEC_IST_operations() 
		rand_BEC = bec.get_rand_BEC()
		for i in range(len(self.BEC_IST_operations)):
			for j in range(len(self.BEC_IST_operations[i][2])):
				self.assertTrue(np.allclose(rand_BEC[self.BEC_IST_operations[i][0]],
					self.BEC_IST_operations[i][2][j].transform_tensor(rand_BEC[self.BEC_IST_operations[i][1]]), atol = 1e-03))

	def test_get_rand_IST(self):
		ist = InternalStrainTensor(self.piezo_struc, self.IST, self.pointops)
		ist.get_BEC_IST_operations(self.BEC) 
		rand_IST = ist.get_rand_IST()
		for i in range(len(self.BEC_IST_operations)):
			for j in range(len(self.BEC_IST_operations[i][2])):
				self.assertTrue(np.allclose(rand_IST[self.BEC_IST_operations[i][0]],
					self.BEC_IST_operations[i][2][j].transform_tensor(rand_IST[self.BEC_IST_operations[i][1]]), atol = 1e-03))
	def test_get_FCM_operations(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations() 
		self.assertTrue(np.all(fcm.FCM_operations == self.FCM_operations))
	def test_get_unstable_FCM(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations()
		rand_FCM = fcm.get_unstable_FCM()
		rand_FCM = np.reshape(rand_FCM, (10,3,10,3)).swapaxes(1,2)
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(rand_FCM[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
                    rand_FCM[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))
	def test_get_FCM_symmetry(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations()
		
		fcm = fcm.get_symmetrized_FCM(np.random.rand(30,30))
		fcm = np.reshape(fcm, (10,3,10,3)).swapaxes(1,2)
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(fcm[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
                    fcm[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))
	def test_get_asum_FCM(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations()
		rand_FCM = fcm.get_unstable_FCM()
		rand_FCM = fcm.get_asum_FCM(rand_FCM)
		rand_FCM = np.reshape(rand_FCM, (10,3,10,3)).swapaxes(1,2)
		
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(rand_FCM[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
                    rand_FCM[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))

		for i in range(len(rand_FCM)):
			asum1 = np.zeros([3,3])
			asum2 = np.zeros([3,3])
			for j in range(len(rand_FCM[i])):
				asum1 += rand_FCM[i][j]
				asum2 += rand_FCM[j][i]
			self.assertTrue(np.allclose(asum1, np.zeros([3,3]), atol = 1e-05))
			self.assertTrue(np.allclose(asum2, np.zeros([3,3]), atol = 1e-05))

	def test_get_stable_FCM(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations()
		rand_FCM = fcm.get_unstable_FCM()
		rand_FCM1 = fcm.get_stable_FCM(rand_FCM)

		eigs, vecs = np.linalg.eig(rand_FCM1)
		eigsort = np.argsort(np.abs(eigs))
		for i in range(3,len(eigs)):
			self.assertTrue(eigs[eigsort[i]] < 1e-06)

		rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)
		
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(rand_FCM1[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
				rand_FCM1[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))

		for i in range(len(rand_FCM1)):
			asum1 = np.zeros([3,3])
			asum2 = np.zeros([3,3])
			for j in range(len(rand_FCM1[i])):
				asum1 += rand_FCM1[i][j]
				asum2 += rand_FCM1[j][i]
			self.assertTrue(np.allclose(asum1, np.zeros([3,3]), atol = 1e-05))
			self.assertTrue(np.allclose(asum2, np.zeros([3,3]), atol = 1e-05))


	def test_get_stable_FCM(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations()
		rand_FCM = fcm.get_unstable_FCM()
		rand_FCM1 = fcm.get_stable_FCM(rand_FCM)

		eigs, vecs = np.linalg.eig(rand_FCM1)
		eigsort = np.argsort(np.abs(eigs))
		for i in range(3,len(eigs)):
			self.assertTrue(eigs[eigsort[i]] < 1e-06)

		rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)
		
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(rand_FCM1[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
				rand_FCM1[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))

		for i in range(len(rand_FCM1)):
			asum1 = np.zeros([3,3])
			asum2 = np.zeros([3,3])
			for j in range(len(rand_FCM1[i])):
				asum1 += rand_FCM1[i][j]
				asum2 += rand_FCM1[j][i]
			self.assertTrue(np.allclose(asum1, np.zeros([3,3]), atol = 1e-05))
			self.assertTrue(np.allclose(asum2, np.zeros([3,3]), atol = 1e-05))

	def test_rand_FCM(self):
		fcm = ForceConstantMatrix(self.piezo_struc, self.FCM, self.pointops, self.sharedops)
		fcm.get_FCM_operations()
		rand_FCM = fcm.get_rand_FCM()
		structure = pymatgen.io.phonopy.get_phonopy_structure(self.piezo_struc)
		pnstruc = Phonopy(structure, np.eye(3), np.eye(3))

		pnstruc.set_force_constants(rand_FCM)
		dyn = pnstruc.get_dynamical_matrix_at_q([0,0,0])
		dyn = np.reshape(dyn, (10,3,10,3)).swapaxes(1,2)
		dyn = np.real(dyn)
		numsites = len(self.piezo_struc)
		masses = []
		for j in range(numsites):
			masses.append(self.piezo_struc.sites[j].specie.atomic_mass)
		dynmass = np.zeros([numsites,numsites,3,3])
		for m in range(numsites):
			for n in range(numsites):
				dynmass[m][n] = dyn[m][n]/np.sqrt(masses[m])/np.sqrt(masses[n])

		dynmass = np.reshape(np.swapaxes(dynmass,1,2),(10*3,10*3))
		eigs, vecs = np.linalg.eig(dynmass)
		eigsort = np.argsort(np.abs(eigs))
		for i in range(3,len(eigs)):
			self.assertTrue(eigs[eigsort[i]] < 1e-06)
		# rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)
		
		dynmass = np.reshape(dynmass, (10,3,10,3)).swapaxes(1,2)
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(dynmass[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
				dynmass[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))

		for i in range(len(dynmass)):
			asum1 = np.zeros([3,3])
			asum2 = np.zeros([3,3])
			for j in range(len(dynmass[i])):
				asum1 += dynmass[i][j]
				asum2 += dynmass[j][i]
			self.assertTrue(np.allclose(asum1, np.zeros([3,3]), atol = 1e-05))
			self.assertTrue(np.allclose(asum2, np.zeros([3,3]), atol = 1e-05))

	def test_get_piezo(self):
		piezo = get_piezo(self.BEC, self.IST, self.FCM)
		self.assertTrue(np.allclose(piezo, self.piezo))

	def test_rand_piezo(self):
		rand_BEC, rand_IST, rand_FCM, piezo = rand_piezo(self.piezo_struc, self.pointops, self.sharedops, self.BEC, self.IST, self.FCM)

		for i in range(len(self.BEC_IST_operations)):
			for j in range(len(self.BEC_IST_operations[i][2])):
				self.assertTrue(np.allclose(rand_BEC[self.BEC_IST_operations[i][0]],
					self.BEC_IST_operations[i][2][j].transform_tensor(rand_BEC[self.BEC_IST_operations[i][1]]), atol = 1e-03))

		for i in range(len(self.BEC_IST_operations)):
			for j in range(len(self.BEC_IST_operations[i][2])):
				self.assertTrue(np.allclose(rand_IST[self.BEC_IST_operations[i][0]],
					self.BEC_IST_operations[i][2][j].transform_tensor(rand_IST[self.BEC_IST_operations[i][1]]), atol = 1e-03))

		structure = pymatgen.io.phonopy.get_phonopy_structure(self.piezo_struc)
		pnstruc = Phonopy(structure, np.eye(3), np.eye(3))

		pnstruc.set_force_constants(rand_FCM)
		dyn = pnstruc.get_dynamical_matrix_at_q([0,0,0])
		dyn = np.reshape(dyn, (10,3,10,3)).swapaxes(1,2)
		dyn = np.real(dyn)
		numsites = len(self.piezo_struc)
		masses = []
		for j in range(numsites):
			masses.append(self.piezo_struc.sites[j].specie.atomic_mass)
		dynmass = np.zeros([numsites,numsites,3,3])
		for m in range(numsites):
			for n in range(numsites):
				dynmass[m][n] = dyn[m][n]/np.sqrt(masses[m])/np.sqrt(masses[n])

		dynmass = np.reshape(np.swapaxes(dynmass,1,2),(10*3,10*3))
		eigs, vecs = np.linalg.eig(dynmass)
		eigsort = np.argsort(np.abs(eigs))
		for i in range(3,len(eigs)):
			self.assertTrue(eigs[eigsort[i]] < 1e-06)
		# rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)
		
		dynmass = np.reshape(dynmass, (10,3,10,3)).swapaxes(1,2)
		for i in range(len(self.FCM_operations)):
			for j in range(len(self.FCM_operations[i][4])):
				self.assertTrue(np.allclose(self.FCM_operations[i][4][j].transform_tensor(dynmass[self.FCM_operations[i][2]][self.FCM_operations[i][3]]),
				dynmass[self.FCM_operations[i][0]][self.FCM_operations[i][1]], atol = 1e-04))

		for i in range(len(dynmass)):
			asum1 = np.zeros([3,3])
			asum2 = np.zeros([3,3])
			for j in range(len(dynmass[i])):
				asum1 += dynmass[i][j]
				asum2 += dynmass[j][i]
			self.assertTrue(np.allclose(asum1, np.zeros([3,3]), atol = 1e-05))
			self.assertTrue(np.allclose(asum2, np.zeros([3,3]), atol = 1e-05))
if __name__ == '__main__':
    unittest.main()



