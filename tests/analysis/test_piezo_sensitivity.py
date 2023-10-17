"""Test for the piezo tensor class."""

from __future__ import annotations

import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pymatgen
from pymatgen.analysis.piezo_sensitivity import (
    BornEffectiveCharge,
    ForceConstantMatrix,
    InternalStrainTensor,
    get_piezo,
    rand_piezo,
)
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    from phonopy import Phonopy
except ImportError:
    Phonopy = None

__author__ = "Handong Ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__status__ = "Development"
__date__ = "4/23/19"

test_dir = f"{TEST_FILES_DIR}/piezo_sensitivity"


class TestPiezoSensitivity(PymatgenTest):
    def setUp(self):
        self.piezo_struct = self.get_structure("Pb2TiZrO6")
        self.IST = np.load(f"{test_dir}/pztist.npy", allow_pickle=True)
        self.BEC = np.load(f"{test_dir}/pztborn.npy", allow_pickle=True)
        self.FCM = np.load(f"{test_dir}/pztfcm.npy", allow_pickle=True)
        self.pointops = np.load(f"{test_dir}/pointops.npy", allow_pickle=True)
        self.sharedops = np.load(f"{test_dir}/sharedops.npy", allow_pickle=True)
        self.IST_operations = np.load(f"{test_dir}/istops.npy", allow_pickle=True)
        with open(f"{test_dir}/becops.pkl", "rb") as file:
            self.BEC_operations = pickle.load(file)
        with open(f"{test_dir}/fcmops.pkl", "rb") as file:
            self.FCM_operations = pickle.load(file)
        self.piezo = np.array(
            [
                [
                    [5.32649351e-3, -1.33404642e-14, -6.86958142e02],
                    [-1.33404642e-14, 4.95526253e-3, -5.60353712e-13],
                    [-6.86958142e02, -5.60353712e-13, 1.33209787e-2],
                ],
                [
                    [4.86622567e-3, 3.14840965e-13, -7.41608150e-13],
                    [3.14840965e-13, 5.23745666e-3, -6.68536818e02],
                    [-7.41608150e-13, -6.68536818e02, 1.35025755e-2],
                ],
                [
                    [-1.01086427e02, 3.20177004e-14, -3.68487214e-14],
                    [3.20177004e-14, -1.01086427e02, 1.22012318e-14],
                    [-3.68487214e-14, 1.22012318e-14, -5.32241086e02],
                ],
            ]
        )

    def test_born_effective_charge_tensor(self):
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.pointops)
        assert_allclose(self.BEC, bec.bec)

    def test_internal_strain_tensor(self):
        ist = InternalStrainTensor(self.piezo_struct, self.IST, self.pointops)
        assert_allclose(ist.ist, self.IST)

    def test_force_constant_matrix(self):
        fcmt = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        assert_allclose(fcmt.fcm, self.FCM)

    def test_get_bec_operations(self):
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.pointops)
        # update test file
        # with open(f"{test_dir}/becops.pkl", "wb") as file:
        #     pickle.dump(bec.get_BEC_operations(), file)
        bec.get_BEC_operations()
        assert np.all(self.BEC_operations == bec.BEC_operations)

    def test_get_rand_bec(self):
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.pointops)
        bec.get_BEC_operations()
        rand_BEC = bec.get_rand_BEC()
        for i in range(len(self.BEC_operations)):
            for j in range(len(self.BEC_operations[i][2])):
                assert_allclose(
                    rand_BEC[self.BEC_operations[i][0]],
                    self.BEC_operations[i][2][j].transform_tensor(rand_BEC[self.BEC_operations[i][1]]),
                    atol=1e-3,
                )

    def test_get_rand_ist(self):
        ist = InternalStrainTensor(self.piezo_struct, self.IST, self.pointops)
        ist.get_IST_operations()
        rand_IST = ist.get_rand_IST()
        for i in range(len(self.IST_operations)):
            for j in range(len(self.IST_operations[i])):
                assert_allclose(
                    rand_IST[i],
                    self.IST_operations[i][j][1].transform_tensor(rand_IST[self.IST_operations[i][j][0]]),
                    atol=1e-3,
                )

    def test_get_fcm_operations(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        # update test file
        # with open(f"{test_dir}/fcmops.pkl", "wb") as file:
        #     pickle.dump(fcm.get_FCM_operations(), file)
        fcm.get_FCM_operations()
        assert np.all(fcm.FCM_operations == self.FCM_operations)

    def test_get_unstable_fcm(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_unstable_FCM()
        rand_FCM = np.reshape(rand_FCM, (10, 3, 10, 3)).swapaxes(1, 2)
        for i in range(len(self.FCM_operations)):
            for j in range(len(self.FCM_operations[i][4])):
                assert_allclose(
                    self.FCM_operations[i][4][j].transform_tensor(
                        rand_FCM[self.FCM_operations[i][2]][self.FCM_operations[i][3]]
                    ),
                    rand_FCM[self.FCM_operations[i][0]][self.FCM_operations[i][1]],
                    atol=1e-4,
                )

    def test_get_fcm_symmetry(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        fcm.get_FCM_operations()

        fcm = fcm.get_symmetrized_FCM(np.random.rand(30, 30))
        fcm = np.reshape(fcm, (10, 3, 10, 3)).swapaxes(1, 2)
        for i in range(len(self.FCM_operations)):
            for j in range(len(self.FCM_operations[i][4])):
                assert_allclose(
                    self.FCM_operations[i][4][j].transform_tensor(
                        fcm[self.FCM_operations[i][2]][self.FCM_operations[i][3]]
                    ),
                    fcm[self.FCM_operations[i][0]][self.FCM_operations[i][1]],
                    atol=1e-4,
                )

    def test_get_asum_fcm(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_unstable_FCM()
        rand_FCM = fcm.get_asum_FCM(rand_FCM)
        rand_FCM = np.reshape(rand_FCM, (10, 3, 10, 3)).swapaxes(1, 2)

        for i in range(len(self.FCM_operations)):
            for j in range(len(self.FCM_operations[i][4])):
                assert_allclose(
                    self.FCM_operations[i][4][j].transform_tensor(
                        rand_FCM[self.FCM_operations[i][2]][self.FCM_operations[i][3]]
                    ),
                    rand_FCM[self.FCM_operations[i][0]][self.FCM_operations[i][1]],
                    atol=1e-4,
                )

        for i in range(len(rand_FCM)):
            asum1 = np.zeros([3, 3])
            asum2 = np.zeros([3, 3])
            for j in range(len(rand_FCM[i])):
                asum1 += rand_FCM[i][j]
                asum2 += rand_FCM[j][i]
            assert_allclose(asum1, np.zeros([3, 3]), atol=1e-5)
            assert_allclose(asum2, np.zeros([3, 3]), atol=1e-5)

    def test_get_stable_fcm(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_unstable_FCM()
        rand_FCM1 = fcm.get_stable_FCM(rand_FCM)

        eigs, vecs = np.linalg.eig(rand_FCM1)
        eigsort = np.argsort(np.abs(eigs))
        for i in range(3, len(eigs)):
            assert eigs[eigsort[i]] < 1e-6

        rand_FCM1 = np.reshape(rand_FCM1, (10, 3, 10, 3)).swapaxes(1, 2)

        for i in range(len(self.FCM_operations)):
            for j in range(len(self.FCM_operations[i][4])):
                assert_allclose(
                    self.FCM_operations[i][4][j].transform_tensor(
                        rand_FCM1[self.FCM_operations[i][2]][self.FCM_operations[i][3]]
                    ),
                    rand_FCM1[self.FCM_operations[i][0]][self.FCM_operations[i][1]],
                    atol=1e-4,
                )

        for i in range(len(rand_FCM1)):
            asum1 = np.zeros([3, 3])
            asum2 = np.zeros([3, 3])
            for j in range(len(rand_FCM1[i])):
                asum1 += rand_FCM1[i][j]
                asum2 += rand_FCM1[j][i]
            assert_allclose(asum1, np.zeros([3, 3]), atol=1e-5)
            assert_allclose(asum2, np.zeros([3, 3]), atol=1e-5)

    def test_rand_fcm(self):
        pytest.importorskip("phonopy")
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.pointops, self.sharedops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_rand_FCM()
        structure = pymatgen.io.phonopy.get_phonopy_structure(self.piezo_struct)
        pn_struct = Phonopy(structure, np.eye(3), np.eye(3))

        pn_struct.set_force_constants(rand_FCM)
        dyn = pn_struct.get_dynamical_matrix_at_q([0, 0, 0])
        dyn = np.reshape(dyn, (10, 3, 10, 3)).swapaxes(1, 2)
        dyn = np.real(dyn)
        n_sites = len(self.piezo_struct)
        masses = [site.specie.atomic_mass for site in self.piezo_struct]
        dyn_mass = np.zeros([n_sites, n_sites, 3, 3])
        for m in range(n_sites):
            for n in range(n_sites):
                dyn_mass[m][n] = dyn[m][n] / np.sqrt(masses[m]) / np.sqrt(masses[n])

        dyn_mass = np.reshape(np.swapaxes(dyn_mass, 1, 2), (10 * 3, 10 * 3))
        eigs, vecs = np.linalg.eig(dyn_mass)
        eigsort = np.argsort(np.abs(eigs))
        for i in range(3, len(eigs)):
            assert eigs[eigsort[i]] < 1e-6
        # rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)

        dyn_mass = np.reshape(dyn_mass, (10, 3, 10, 3)).swapaxes(1, 2)
        for i in range(len(self.FCM_operations)):
            for j in range(len(self.FCM_operations[i][4])):
                assert_allclose(
                    self.FCM_operations[i][4][j].transform_tensor(
                        dyn_mass[self.FCM_operations[i][2]][self.FCM_operations[i][3]]
                    ),
                    dyn_mass[self.FCM_operations[i][0]][self.FCM_operations[i][1]],
                    atol=1e-4,
                )

        for i in range(len(dyn_mass)):
            asum1 = np.zeros([3, 3])
            asum2 = np.zeros([3, 3])
            for j in range(len(dyn_mass[i])):
                asum1 += dyn_mass[i][j]
                asum2 += dyn_mass[j][i]
            assert_allclose(asum1, np.zeros([3, 3]), atol=1e-5)
            assert_allclose(asum2, np.zeros([3, 3]), atol=1e-5)

    def test_get_piezo(self):
        piezo = get_piezo(self.BEC, self.IST, self.FCM)
        assert_allclose(piezo, self.piezo, atol=1e-5)

    def test_rand_piezo(self):
        pytest.importorskip("phonopy")
        rand_BEC, rand_IST, rand_FCM, piezo = rand_piezo(
            self.piezo_struct, self.pointops, self.sharedops, self.BEC, self.IST, self.FCM
        )

        for i in range(len(self.BEC_operations)):
            for j in range(len(self.BEC_operations[i][2])):
                assert_allclose(
                    rand_BEC[self.BEC_operations[i][0]],
                    self.BEC_operations[i][2][j].transform_tensor(rand_BEC[self.BEC_operations[i][1]]),
                    atol=1e-3,
                )

        for i in range(len(self.IST_operations)):
            for j in range(len(self.IST_operations[i])):
                assert_allclose(
                    rand_IST[i],
                    self.IST_operations[i][j][1].transform_tensor(rand_IST[self.IST_operations[i][j][0]]),
                    atol=1e-3,
                )

        structure = pymatgen.io.phonopy.get_phonopy_structure(self.piezo_struct)
        pn_struct = Phonopy(structure, np.eye(3), np.eye(3))

        pn_struct.set_force_constants(rand_FCM)
        dyn = pn_struct.get_dynamical_matrix_at_q([0, 0, 0])
        dyn = np.reshape(dyn, (10, 3, 10, 3)).swapaxes(1, 2)
        dyn = np.real(dyn)
        n_sites = len(self.piezo_struct)
        masses = [site.specie.atomic_mass for site in self.piezo_struct]
        dyn_mass = np.zeros([n_sites, n_sites, 3, 3])
        for m in range(n_sites):
            for n in range(n_sites):
                dyn_mass[m][n] = dyn[m][n] / np.sqrt(masses[m]) / np.sqrt(masses[n])

        dyn_mass = np.reshape(np.swapaxes(dyn_mass, 1, 2), (10 * 3, 10 * 3))
        eigs, _eig_vecs = np.linalg.eig(dyn_mass)
        eig_sorted = np.argsort(np.abs(eigs))
        for i in range(3, len(eigs)):
            assert eigs[eig_sorted[i]] < 1e-6
        # rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)

        dyn_mass = np.reshape(dyn_mass, (10, 3, 10, 3)).swapaxes(1, 2)
        for i in range(len(self.FCM_operations)):
            for j in range(len(self.FCM_operations[i][4])):
                assert_allclose(
                    self.FCM_operations[i][4][j].transform_tensor(
                        dyn_mass[self.FCM_operations[i][2]][self.FCM_operations[i][3]]
                    ),
                    dyn_mass[self.FCM_operations[i][0]][self.FCM_operations[i][1]],
                    atol=1e-4,
                )

        for i in range(len(dyn_mass)):
            asum1 = np.zeros([3, 3])
            asum2 = np.zeros([3, 3])
            for j in range(len(dyn_mass[i])):
                asum1 += dyn_mass[i][j]
                asum2 += dyn_mass[j][i]
            assert_allclose(asum1, np.zeros([3, 3]), atol=1e-5)
            assert_allclose(asum2, np.zeros([3, 3]), atol=1e-5)
