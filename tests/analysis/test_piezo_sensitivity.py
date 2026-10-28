"""Test for the piezo tensor class."""

from __future__ import annotations

import pickle
import platform
import warnings
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pymatgen.analysis import piezo_sensitivity
from pymatgen.analysis.piezo_sensitivity import (
    BornEffectiveCharge,
    ForceConstantMatrix,
    InternalStrainTensor,
    get_piezo,
    rand_piezo,
)
from pymatgen.core import Lattice, Structure
from pymatgen.io.phonopy import get_phonopy_structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

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

TEST_DIR = f"{TEST_FILES_DIR}/analysis/piezo_sensitivity"


class _SeededNumpy:
    """Stand-in for the numpy module inside piezo_sensitivity whose
    random.default_rng() always returns one seeded generator. Patching the
    module-level name (rather than numpy itself) keeps third-party code such as
    phonopy from consuming the seeded stream.
    """

    def __init__(self, seed: int):
        rng = np.random.default_rng(seed)
        self.random = SimpleNamespace(default_rng=lambda *_args, **_kwargs: rng)

    def __getattr__(self, name):
        return getattr(np, name)


def seed_piezo_rng(monkeypatch: pytest.MonkeyPatch, seed: int) -> None:
    """(Re)seed all random number generation in piezo_sensitivity."""
    monkeypatch.setattr(piezo_sensitivity, "np", _SeededNumpy(seed))


class TestPiezoSensitivity(MatSciTest):
    @pytest.fixture(autouse=True)
    def _seed_rng(self, monkeypatch):
        seed_piezo_rng(monkeypatch, 42)

    def setup_method(self):
        self.piezo_struct = self.get_structure("Pb2TiZrO6")
        self.IST = np.load(f"{TEST_DIR}/pztist.npy", allow_pickle=True)
        self.BEC = np.load(f"{TEST_DIR}/pztborn.npy", allow_pickle=True)
        self.FCM = np.load(f"{TEST_DIR}/pztfcm.npy", allow_pickle=True)
        self.point_ops = np.load(f"{TEST_DIR}/pointops.npy", allow_pickle=True)
        self.shared_ops = np.load(f"{TEST_DIR}/sharedops.npy", allow_pickle=True)
        self.IST_operations = np.load(f"{TEST_DIR}/istops.npy", allow_pickle=True)
        with open(f"{TEST_DIR}/becops.pkl", "rb") as file:
            self.BEC_operations = pickle.load(file)
        with open(f"{TEST_DIR}/fcmops.pkl", "rb") as file:
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
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.point_ops)
        assert_allclose(self.BEC, bec.bec)

    def test_internal_strain_tensor(self):
        ist = InternalStrainTensor(self.piezo_struct, self.IST, self.point_ops)
        assert_allclose(ist.ist, self.IST)

    def test_force_constant_matrix(self):
        fcmt = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        assert_allclose(fcmt.fcm, self.FCM)

    def test_get_bec_operations(self):
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.point_ops)
        # update test file
        # with open(f"{TEST_DIR}/becops.pkl", "wb") as file:
        #     pickle.dump(bec.get_BEC_operations(), file)
        bec.get_BEC_operations()
        assert np.all(self.BEC_operations == bec.BEC_operations)

    def test_get_rand_bec(self):
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.point_ops)
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
        ist = InternalStrainTensor(self.piezo_struct, self.IST, self.point_ops)
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
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        # update test file
        # with open(f"{TEST_DIR}/fcmops.pkl", "wb") as file:
        #     pickle.dump(fcm.get_FCM_operations(), file)
        fcm.get_FCM_operations()
        assert np.all(fcm.FCM_operations == self.FCM_operations)

    def test_get_unstable_fcm(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
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
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        fcm.get_FCM_operations()

        fcm = fcm.get_symmetrized_FCM(np.random.default_rng(0).random((30, 30)))
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
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_unstable_FCM()
        rand_FCM = fcm.get_asum_FCM(rand_FCM)
        rand_FCM = np.reshape(rand_FCM, (10, 3, 10, 3)).swapaxes(1, 2)

        for ii in range(len(self.FCM_operations)):
            for jj in range(len(self.FCM_operations[ii][4])):
                assert_allclose(
                    self.FCM_operations[ii][4][jj].transform_tensor(
                        rand_FCM[self.FCM_operations[ii][2]][self.FCM_operations[ii][3]]
                    ),
                    rand_FCM[self.FCM_operations[ii][0]][self.FCM_operations[ii][1]],
                    atol=1e-4,
                )

        for ii in range(len(rand_FCM)):
            sum1 = np.zeros([3, 3])
            sum2 = np.zeros([3, 3])
            for jj in range(len(rand_FCM[ii])):
                sum1 += rand_FCM[ii][jj]
                sum2 += rand_FCM[jj][ii]
            assert_allclose(sum1, np.zeros([3, 3]), atol=1e-5)
            assert_allclose(sum2, np.zeros([3, 3]), atol=1e-5)

    def test_get_stable_fcm(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_unstable_FCM()
        rand_FCM1 = fcm.get_stable_FCM(rand_FCM)

        eig_vals, _vecs = np.linalg.eig(rand_FCM1)
        eig_sorted = np.argsort(np.abs(eig_vals))
        for i in range(3, len(eig_vals)):
            assert eig_vals[eig_sorted[i]] < 1e-6

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

    @pytest.mark.xfail(
        platform.system() == "Windows" and int(np.__version__[0]) >= 2,
        reason="See https://github.com/conda-forge/phonopy-feedstock/pull/158#issuecomment-2227506701",
    )
    def test_rand_fcm(self):
        pytest.importorskip("phonopy")
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        fcm.get_FCM_operations()
        rand_FCM = fcm.get_rand_FCM()
        structure = get_phonopy_structure(self.piezo_struct)
        pn_struct = Phonopy(structure, np.eye(3), np.eye(3))

        pn_struct.force_constants = rand_FCM
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
        eig_vals, _eig_vecs = np.linalg.eig(dyn_mass)
        eig_sort = np.argsort(np.abs(eig_vals))
        for i in range(3, len(eig_vals)):
            assert eig_vals[eig_sort[i]] < 1e-6
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

    @pytest.mark.xfail(
        platform.system() == "Windows" and int(np.__version__[0]) >= 2,
        reason="See https://github.com/conda-forge/phonopy-feedstock/pull/158#issuecomment-2227506701",
    )
    def test_rand_piezo(self):
        pytest.importorskip("phonopy")
        rand_BEC, rand_IST, rand_FCM, _piezo = rand_piezo(
            self.piezo_struct,
            self.point_ops,
            self.shared_ops,
            self.BEC,
            self.IST,
            self.FCM,
        )

        for ii in range(len(self.BEC_operations)):
            for jj in range(len(self.BEC_operations[ii][2])):
                assert_allclose(
                    rand_BEC[self.BEC_operations[ii][0]],
                    self.BEC_operations[ii][2][jj].transform_tensor(rand_BEC[self.BEC_operations[ii][1]]),
                    atol=1e-3,
                )

        for ii in range(len(self.IST_operations)):
            for jj in range(len(self.IST_operations[ii])):
                assert_allclose(
                    rand_IST[ii],
                    self.IST_operations[ii][jj][1].transform_tensor(rand_IST[self.IST_operations[ii][jj][0]]),
                    atol=1e-3,
                )

        structure = get_phonopy_structure(self.piezo_struct)
        pn_struct = Phonopy(structure, np.eye(3), np.eye(3))

        pn_struct.force_constants = rand_FCM
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
        eig_vals, _eig_vecs = np.linalg.eig(dyn_mass)
        eig_sort = np.argsort(np.abs(eig_vals))
        for idx in range(3, len(eig_vals)):
            assert eig_vals[eig_sort[idx]] < 1e-6
        # rand_FCM1 = np.reshape(rand_FCM1, (10,3,10,3)).swapaxes(1,2)

        dyn_mass = np.reshape(dyn_mass, (10, 3, 10, 3)).swapaxes(1, 2)
        for ii in range(len(self.FCM_operations)):
            for jj in range(len(self.FCM_operations[ii][4])):
                assert_allclose(
                    self.FCM_operations[ii][4][jj].transform_tensor(
                        dyn_mass[self.FCM_operations[ii][2]][self.FCM_operations[ii][3]]
                    ),
                    dyn_mass[self.FCM_operations[ii][0]][self.FCM_operations[ii][1]],
                    atol=1e-4,
                )

        for ii in range(len(dyn_mass)):
            asum1 = np.zeros([3, 3])
            asum2 = np.zeros([3, 3])
            for jj in range(len(dyn_mass[ii])):
                asum1 += dyn_mass[ii][jj]
                asum2 += dyn_mass[jj][ii]
            assert_allclose(asum1, np.zeros([3, 3]), atol=1e-5)
            assert_allclose(asum2, np.zeros([3, 3]), atol=1e-5)

    @pytest.mark.parametrize("delta", [0.1, -0.1])
    def test_bec_charge_neutrality_warning(self, delta):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            BornEffectiveCharge(self.piezo_struct, self.BEC, self.point_ops)

        bec = self.BEC.copy()
        bec[0, 2, 2] += delta
        # A violation of either sign must be reported
        with pytest.warns(UserWarning, match="does not satisfy charge neutrality"):
            BornEffectiveCharge(self.piezo_struct, bec, self.point_ops)

    def test_ist_symmetry_warning(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            InternalStrainTensor(self.piezo_struct, self.IST, self.point_ops)

        ist = self.IST.copy()
        ist[0, 0, 0, 1] += 0.1
        with pytest.warns(UserWarning, match="does not satisfy standard symmetries"):
            InternalStrainTensor(self.piezo_struct, ist, self.point_ops)

    def test_get_rand_bec_sum_rule_and_site_symmetry(self, monkeypatch):
        bec = BornEffectiveCharge(self.piezo_struct, self.BEC, self.point_ops)
        bec.get_BEC_operations()
        rand_BEC = bec.get_rand_BEC()
        assert rand_BEC.shape == (10, 3, 3)
        assert np.abs(rand_BEC).max() > 1e-2
        # Acoustic sum rule (charge neutrality): sum_k Z*_k = 0
        assert_allclose(rand_BEC.sum(axis=0), np.zeros((3, 3)), atol=1e-10)
        # Each tensor is invariant under its own site symmetry
        for site_idx, ops in enumerate(self.point_ops):
            for op in ops:
                assert_allclose(op.transform_tensor(rand_BEC[site_idx]), rand_BEC[site_idx], atol=1e-10)

        # Same seed, scaled output
        seed_piezo_rng(monkeypatch, 42)
        assert_allclose(bec.get_rand_BEC(max_charge=2.5), 2.5 * rand_BEC)

    def test_get_rand_ist_symmetries(self, monkeypatch):
        ist = InternalStrainTensor(self.piezo_struct, self.IST, self.point_ops)
        ist_ops = ist.get_IST_operations()
        assert ist_ops == ist.IST_operations
        assert len(ist_ops) == 10
        assert ist_ops[0] == []  # the first site has nothing before it to map from
        rand_IST = ist.get_rand_IST()
        assert rand_IST.shape == (10, 3, 3, 3)
        assert np.abs(rand_IST).max() > 1e-2
        # Symmetric in the strain (last two) indices
        assert_allclose(rand_IST, rand_IST.transpose(0, 1, 3, 2), atol=1e-12)
        for site_idx, ops in enumerate(self.point_ops):
            for op in ops:
                assert_allclose(op.transform_tensor(rand_IST[site_idx]), rand_IST[site_idx], atol=1e-10)

        seed_piezo_rng(monkeypatch, 42)
        assert_allclose(ist.get_rand_IST(max_force=3), 3 * rand_IST)

    def test_get_unstable_and_symmetrized_fcm_are_symmetric(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        fcm.get_FCM_operations()
        unstable = fcm.get_unstable_FCM()
        assert unstable.shape == (30, 30)
        assert_allclose(unstable, unstable.T)

        sym = fcm.get_symmetrized_FCM(np.random.default_rng(0).random((30, 30)))
        assert_allclose(sym, sym.T)
        # Symmetrizing an already symmetrized matrix is a no-op
        assert_allclose(fcm.get_symmetrized_FCM(sym.copy()), sym, atol=1e-12)

    def test_get_asum_fcm_requires_operations(self):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        with pytest.raises(RuntimeError, match="Run get_FCM_operations first"):
            fcm.get_asum_FCM(np.zeros((30, 30)))

    def test_get_stable_fcm_gives_up_after_20_retries(self, monkeypatch):
        fcm = ForceConstantMatrix(self.piezo_struct, self.FCM, self.point_ops, self.shared_ops)
        fcm.get_FCM_operations()
        n_calls = 0

        def always_unstable(_matrix):
            nonlocal n_calls
            n_calls += 1
            return np.eye(30)  # all eigenvalues positive, i.e. unstable

        monkeypatch.setattr(fcm, "get_symmetrized_FCM", always_unstable)
        monkeypatch.setattr(fcm, "get_asum_FCM", lambda matrix: matrix)
        out = fcm.get_stable_FCM(-np.eye(30))
        assert n_calls == 21
        assert_allclose(out, np.eye(30))

    def test_get_piezo_symmetry_and_linearity(self):
        piezo = get_piezo(self.BEC, self.IST, self.FCM)
        assert piezo.shape == (3, 3, 3)
        # e_ijk is symmetric in the strain indices
        assert_allclose(piezo, piezo.transpose(0, 2, 1), atol=1e-8)
        # Linear in the Born charges and in the internal strain tensor
        assert_allclose(get_piezo(2 * self.BEC, self.IST, self.FCM), 2 * piezo, atol=1e-8)
        assert_allclose(get_piezo(self.BEC, -self.IST, self.FCM), -piezo, atol=1e-8)

    def test_rand_piezo_consistency(self, monkeypatch):
        pytest.importorskip("phonopy")
        rand_BEC, rand_IST, rand_FCM, piezo = rand_piezo(
            self.piezo_struct, self.point_ops, self.shared_ops, self.BEC, self.IST, self.FCM
        )
        assert rand_BEC.shape == (10, 3, 3)
        assert rand_IST.shape == (10, 3, 3, 3)
        assert rand_FCM.shape == (10, 10, 3, 3)
        assert piezo.shape == (3, 3, 3)
        assert_allclose(rand_BEC.sum(axis=0), np.zeros((3, 3)), atol=1e-10)
        assert_allclose(
            piezo,
            get_piezo(rand_BEC, rand_IST, rand_FCM) * 16.0216559424 / self.piezo_struct.volume,
        )
        assert_allclose(piezo, piezo.transpose(0, 2, 1), atol=1e-12)

        # Reproducible for a fixed seed
        seed_piezo_rng(monkeypatch, 42)
        *_, piezo2 = rand_piezo(self.piezo_struct, self.point_ops, self.shared_ops, self.BEC, self.IST, self.FCM)
        assert_allclose(piezo2, piezo)


def test_get_fcm_operations_transpose_relation():
    """A pair block that is only related to an earlier block by transposition
    (F_02 = F_01^T) must be mapped onto the transposed pair (1, 0).
    """
    struct = Structure(
        Lattice.from_parameters(3.1, 4.2, 5.3, 81, 86, 97),
        ["Na", "Cl", "K"],
        [[0, 0, 0], [0.31, 0.22, 0.13], [0.6, 0.7, 0.45]],
    )
    assert SpacegroupAnalyzer(struct).get_space_group_symbol() == "P1"
    asym = np.array([[1.0, 0.5, 0.0], [0.0, 4.0, 0.3], [0.0, 0.0, 6.0]])
    fcm = np.zeros((3, 3, 3, 3))
    fcm[0, 0] = np.diag([10.0, 20.0, 30.0])
    fcm[1, 1] = np.diag([11.0, 21.0, 31.0])
    fcm[2, 2] = np.diag([13.0, 23.0, 33.0])
    fcm[1, 2] = fcm[2, 1] = np.diag([12.0, 22.0, 32.0])
    fcm[0, 1] = fcm[2, 0] = asym
    fcm[1, 0] = fcm[0, 2] = asym.T

    fcm_obj = ForceConstantMatrix(struct, fcm, [[] for _ in struct], None)
    ops = fcm_obj.get_FCM_operations()
    assert [op[:4] for op in ops] == [
        [0, 0, 0, 0],
        [0, 1, 1, 0],
        [0, 2, 1, 0],  # same eigenvalues as (0, 1) -> remapped onto the transpose (1, 0)
        [1, 1, 1, 1],
        [1, 2, 2, 1],
        [2, 2, 2, 2],
    ]
    # Only the identity (P1) relates F_02 to F_10^T = F_01
    (identity,) = ops[2][4]
    assert_allclose(identity.rotation_matrix, np.eye(3), atol=1e-10)
    assert_allclose(identity.transform_tensor(fcm[0, 1]).T, fcm[0, 2], atol=1e-12)
