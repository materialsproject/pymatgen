from __future__ import annotations

import json
import random
import warnings
from copy import deepcopy

import numpy as np
import pytest
from pytest import approx
from scipy.misc import central_diff_weights

from pymatgen.analysis.elasticity.elastic import (
    ComplianceTensor,
    ElasticTensor,
    ElasticTensorExpansion,
    NthOrderElasticTensor,
    diff_fit,
    find_eq_stress,
    generate_pseudo,
    get_diff_coeff,
    get_strain_state_dict,
)
from pymatgen.analysis.elasticity.strain import Deformation, Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.tensors import Tensor
from pymatgen.core.units import FloatWithUnit
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestElasticTensor(PymatgenTest):
    def setUp(self):
        self.voigt_1 = [
            [59.33, 28.08, 28.08, 0, 0, 0],
            [28.08, 59.31, 28.07, 0, 0, 0],
            [28.08, 28.07, 59.32, 0, 0, 0],
            [0, 0, 0, 26.35, 0, 0],
            [0, 0, 0, 0, 26.35, 0],
            [0, 0, 0, 0, 0, 26.35],
        ]
        mat = np.random.randn(6, 6)
        mat = mat + np.transpose(mat)
        self.rand_elastic_tensor = ElasticTensor.from_voigt(mat)
        self.ft = np.array(
            [
                [
                    [[59.33, 0, 0], [0, 28.08, 0], [0, 0, 28.08]],
                    [[0, 26.35, 0], [26.35, 0, 0], [0, 0, 0]],
                    [[0, 0, 26.35], [0, 0, 0], [26.35, 0, 0]],
                ],
                [
                    [[0, 26.35, 0], [26.35, 0, 0], [0, 0, 0]],
                    [[28.08, 0, 0], [0, 59.31, 0], [0, 0, 28.07]],
                    [[0, 0, 0], [0, 0, 26.35], [0, 26.35, 0]],
                ],
                [
                    [[0, 0, 26.35], [0, 0, 0], [26.35, 0, 0]],
                    [[0, 0, 0], [0, 0, 26.35], [0, 26.35, 0]],
                    [[28.08, 0, 0], [0, 28.07, 0], [0, 0, 59.32]],
                ],
            ]
        )

        self.elastic_tensor_1 = ElasticTensor(self.ft)
        filepath = f"{TEST_FILES_DIR}/Sn_def_stress.json"
        with open(filepath) as f:
            self.def_stress_dict = json.load(f)
        with open(f"{TEST_FILES_DIR}/test_toec_data.json") as f:
            self.toec_dict = json.load(f)
        self.structure = self.get_structure("Sn")

        warnings.simplefilter("always")

    def test_properties(self):
        # compliance tensor
        ct = ComplianceTensor.from_voigt(np.linalg.inv(self.elastic_tensor_1.voigt))
        self.assert_all_close(ct, self.elastic_tensor_1.compliance_tensor)
        # KG average properties
        assert self.elastic_tensor_1.k_voigt == approx(38.49111111111)
        assert self.elastic_tensor_1.g_voigt == approx(22.05866666666)
        assert self.elastic_tensor_1.k_reuss == approx(38.49110945133)
        assert self.elastic_tensor_1.g_reuss == approx(20.67146635306)
        assert self.elastic_tensor_1.k_vrh == approx(38.49111028122)
        assert self.elastic_tensor_1.g_vrh == approx(21.36506650986)

        # universal anisotropy
        assert self.elastic_tensor_1.universal_anisotropy == approx(0.33553509658699)
        # homogeneous Poisson
        assert self.elastic_tensor_1.homogeneous_poisson == approx(0.26579965576472)
        # voigt notation tensor
        self.assert_all_close(self.elastic_tensor_1.voigt, self.voigt_1)
        # young's modulus
        assert self.elastic_tensor_1.y_mod == approx(54087787667.160583)

        # prop dict
        prop_dict = self.elastic_tensor_1.property_dict
        assert prop_dict["homogeneous_poisson"] == approx(0.26579965576)
        for k, v in prop_dict.items():
            assert getattr(self.elastic_tensor_1, k) == approx(v)

    def test_directional_elastic_mod(self):
        assert self.elastic_tensor_1.directional_elastic_mod([1, 0, 0]) == approx(self.elastic_tensor_1.voigt[0, 0])
        assert self.elastic_tensor_1.directional_elastic_mod([1, 1, 1]) == approx(73.624444444)

    def test_compliance_tensor(self):
        stress = self.elastic_tensor_1.calculate_stress([0.01] + [0] * 5)
        comp = self.elastic_tensor_1.compliance_tensor
        strain = Strain(comp.einsum_sequence([stress]))
        self.assert_all_close(strain.voigt, [0.01] + [0] * 5)

    def test_directional_poisson_ratio(self):
        v_12 = self.elastic_tensor_1.directional_poisson_ratio([1, 0, 0], [0, 1, 0])
        assert v_12 == approx(0.321388)

    def test_structure_based_methods(self):
        # trans_velocity
        assert self.elastic_tensor_1.trans_v(self.structure) == approx(1996.35019877)
        # long_velocity
        assert self.elastic_tensor_1.long_v(self.structure) == approx(3534.68123832)
        # Snyder properties
        assert self.elastic_tensor_1.snyder_ac(self.structure) == approx(18.06127074)
        assert self.elastic_tensor_1.snyder_opt(self.structure) == approx(0.18937465)
        assert self.elastic_tensor_1.snyder_total(self.structure) == approx(18.25064540)
        # Clarke
        assert self.elastic_tensor_1.clarke_thermalcond(self.structure) == approx(0.3450307)
        # Cahill
        assert self.elastic_tensor_1.cahill_thermalcond(self.structure) == approx(0.37896275)
        # Debye
        assert self.elastic_tensor_1.debye_temperature(self.structure) == approx(198.8037985019)

        # structure-property dict
        sprop_dict = self.elastic_tensor_1.get_structure_property_dict(self.structure)
        assert sprop_dict["long_v"] == approx(3534.68123832)
        for val in sprop_dict.values():
            assert not isinstance(val, FloatWithUnit)
        for k, v in sprop_dict.items():
            if k == "structure":
                assert v == self.structure
            else:
                f = getattr(self.elastic_tensor_1, k)
                if callable(f):
                    assert getattr(self.elastic_tensor_1, k)(self.structure) == approx(v)
                else:
                    assert getattr(self.elastic_tensor_1, k) == approx(v)

        # Test other sprop dict modes
        sprop_dict = self.elastic_tensor_1.get_structure_property_dict(self.structure, include_base_props=False)
        assert "k_vrh" not in sprop_dict

        # Test ValueError being raised for structure properties
        test_et = deepcopy(self.elastic_tensor_1)
        test_et[0][0][0][0] = -100000
        prop_dict = test_et.property_dict
        for attr_name in sprop_dict:
            if attr_name not in ([*prop_dict, "structure"]):
                with pytest.raises(
                    ValueError, match="Bulk or shear modulus is negative, property cannot be determined"
                ):
                    getattr(test_et, attr_name)(self.structure)
        with pytest.raises(ValueError, match="Bulk or shear modulus is negative, property cannot be determined"):
            test_et.get_structure_property_dict(self.structure)
        noval_sprop_dict = test_et.get_structure_property_dict(self.structure, ignore_errors=True)
        assert noval_sprop_dict["snyder_ac"] is None

    def test_new(self):
        self.assert_all_close(self.elastic_tensor_1, ElasticTensor(self.ft))
        non_symm = self.ft
        non_symm[0, 1, 2, 2] += 1.0
        with pytest.warns(
            UserWarning, match="Input elastic tensor does not satisfy standard Voigt symmetries"
        ) as warns:
            ElasticTensor(non_symm)
        assert len(warns) == 1
        bad_tensor1 = np.zeros((3, 3, 3))
        bad_tensor2 = np.zeros((3, 3, 3, 2))
        with pytest.raises(ValueError, match="ElasticTensor input must be rank 4"):
            ElasticTensor(bad_tensor1)
        with pytest.raises(
            ValueError,
            match="Pymatgen only supports 3-dimensional tensors, and default tensor constructor uses standard notation",
        ):
            ElasticTensor(bad_tensor2)

    def test_from_pseudoinverse(self):
        strain_list = [Strain.from_deformation(def_matrix) for def_matrix in self.def_stress_dict["deformations"]]
        stress_list = list(self.def_stress_dict["stresses"])
        with pytest.warns(
            UserWarning,
            match="Pseudo-inverse fitting of Strain/Stress lists may yield questionable results from "
            "vasp data, use with caution",
        ):
            et_fl = -0.1 * ElasticTensor.from_pseudoinverse(strain_list, stress_list).voigt
            self.assert_all_close(
                et_fl.round(2),
                [
                    [59.29, 24.36, 22.46, 0, 0, 0],
                    [28.06, 56.91, 22.46, 0, 0, 0],
                    [28.06, 25.98, 54.67, 0, 0, 0],
                    [0, 0, 0, 26.35, 0, 0],
                    [0, 0, 0, 0, 26.35, 0],
                    [0, 0, 0, 0, 0, 26.35],
                ],
            )

    def test_from_independent_strains(self):
        strains = self.toec_dict["strains"]
        stresses = self.toec_dict["stresses"]
        with pytest.warns(UserWarning, match="No eq state found, returning zero voigt stress") as warns:
            et = ElasticTensor.from_independent_strains(strains, stresses)
        assert len(warns) == 2
        self.assert_all_close(et.voigt, self.toec_dict["C2_raw"], decimal=-1)

    def test_energy_density(self):
        film_elac = ElasticTensor.from_voigt(
            [
                [324.32, 187.3, 170.92, 0, 0, 0],
                [187.3, 324.32, 170.92, 0, 0, 0],
                [170.92, 170.92, 408.41, 0, 0, 0],
                [0, 0, 0, 150.73, 0, 0],
                [0, 0, 0, 0, 150.73, 0],
                [0, 0, 0, 0, 0, 238.74],
            ]
        )

        dfm = Deformation(
            [
                [-9.86004855e-01, 2.27539582e-01, -4.64426035e-17],
                [-2.47802121e-01, -9.91208483e-01, -7.58675185e-17],
                [-6.12323400e-17, -6.12323400e-17, 1.00000000e00],
            ]
        )

        assert film_elac.energy_density(dfm.green_lagrange_strain) == approx(0.00125664672793)

        film_elac.energy_density(
            Strain.from_deformation(
                [
                    [0.99774738, 0.11520994, -0],
                    [-0.11520994, 0.99774738, 0],
                    [-0, -0, 1],
                ]
            )
        )


class TestElasticTensorExpansion(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/test_toec_data.json") as f:
            self.data_dict = json.load(f)
        self.strains = [Strain(sm) for sm in self.data_dict["strains"]]
        self.pk_stresses = [Stress(d) for d in self.data_dict["pk_stresses"]]
        self.c2 = self.data_dict["C2_raw"]
        self.c3 = self.data_dict["C3_raw"]
        self.exp = ElasticTensorExpansion.from_voigt([self.c2, self.c3])
        self.cu = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.623), ["Cu"], [[0] * 3])
        indices = [(0, 0), (0, 1), (3, 3)]
        values = [167.8, 113.5, 74.5]
        cu_c2 = ElasticTensor.from_values_indices(values, indices, structure=self.cu, populate=True)
        indices = [(0, 0, 0), (0, 0, 1), (0, 1, 2), (0, 3, 3), (0, 5, 5), (3, 4, 5)]
        values = [-1507, -965, -71, -7, -901, 45]
        cu_c3 = Tensor.from_values_indices(values, indices, structure=self.cu, populate=True)
        self.exp_cu = ElasticTensorExpansion([cu_c2, cu_c3])
        cu_c4 = Tensor.from_voigt(self.data_dict["Cu_fourth_order"])
        self.exp_cu_4 = ElasticTensorExpansion([cu_c2, cu_c3, cu_c4])
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        c_ijkl = Tensor.from_voigt(self.c2)
        c_ijklmn = Tensor.from_voigt(self.c3)
        exp = ElasticTensorExpansion([c_ijkl, c_ijklmn])
        ElasticTensorExpansion.from_voigt([self.c2, self.c3])
        assert exp.order == 3

    def test_from_diff_fit(self):
        ElasticTensorExpansion.from_diff_fit(self.strains, self.pk_stresses)

    def test_calculate_stress(self):
        calc_stress = self.exp.calculate_stress(self.strains[0])
        self.assert_all_close(self.pk_stresses[0], calc_stress, decimal=2)

    def test_energy_density(self):
        e_density = self.exp.energy_density(self.strains[0])
        assert e_density == approx(1.36363099e-4)

    def test_gruneisen(self):
        # Get GGT
        ggt = self.exp_cu.get_ggt([1, 0, 0], [0, 1, 0])
        self.assert_all_close(np.eye(3) * np.array([4.92080537, 4.2852349, -0.7147651]), ggt)
        # Get TGT
        tgt = self.exp_cu.get_tgt()
        self.assert_all_close(tgt, np.eye(3) * 2.59631832)

        # Get heat capacity
        c0 = self.exp_cu.get_heat_capacity(0, self.cu, [1, 0, 0], [0, 1, 0])
        assert c0 == 0.0
        c = self.exp_cu.get_heat_capacity(300, self.cu, [1, 0, 0], [0, 1, 0])
        assert c == approx(8.285611958)

        # Get Gruneisen parameter
        gp = self.exp_cu.get_gruneisen_parameter()
        assert gp == approx(2.59631832)
        _ = self.exp_cu.get_gruneisen_parameter(temperature=200, structure=self.cu)

    def test_thermal_expansion_coeff(self):
        # TODO get rid of duplicates
        alpha_dp = self.exp_cu.thermal_expansion_coeff(self.cu, 300, mode="dulong-petit")
        alpha_dp_ground_truth = 6.3471959e-07 * np.ones((3, 3))
        alpha_dp_ground_truth[np.diag_indices(3)] = 2.2875769e-7
        self.assert_all_close(alpha_dp_ground_truth, alpha_dp, decimal=4)

        alpha_debye = self.exp_cu.thermal_expansion_coeff(self.cu, 300, mode="debye")
        alpha_comp = 5.9435148e-7 * np.ones((3, 3))
        alpha_comp[np.diag_indices(3)] = 21.4533472e-06
        self.assert_all_close(alpha_comp, alpha_debye)

    def test_get_compliance_expansion(self):
        ce_exp = self.exp_cu.get_compliance_expansion()
        et_comp = ElasticTensorExpansion(ce_exp)
        strain_orig = Strain.from_voigt([0.01, 0, 0, 0, 0, 0])
        stress = self.exp_cu.calculate_stress(strain_orig)
        strain_revert = et_comp.calculate_stress(stress)
        self.assert_all_close(strain_orig, strain_revert, decimal=4)

    def test_get_effective_ecs(self):
        # Ensure zero strain is same as SOEC
        test_zero = self.exp_cu.get_effective_ecs(np.zeros((3, 3)))
        self.assert_all_close(test_zero, self.exp_cu[0])
        s = np.zeros((3, 3))
        s[0, 0] = 0.02
        test_2percent = self.exp_cu.get_effective_ecs(s)
        diff = test_2percent - test_zero
        self.assert_all_close(self.exp_cu[1].einsum_sequence([s]), diff)

    def test_get_strain_from_stress(self):
        strain = Strain.from_voigt([0.05, 0, 0, 0, 0, 0])
        stress3 = self.exp_cu.calculate_stress(strain)
        strain_revert3 = self.exp_cu.get_strain_from_stress(stress3)
        self.assert_all_close(strain, strain_revert3, decimal=2)
        # fourth order
        stress4 = self.exp_cu_4.calculate_stress(strain)
        strain_revert4 = self.exp_cu_4.get_strain_from_stress(stress4)
        self.assert_all_close(strain, strain_revert4, decimal=2)

    def test_get_yield_stress(self):
        self.exp_cu_4.get_yield_stress([1, 0, 0])


class TestNthOrderElasticTensor(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/test_toec_data.json") as f:
            self.data_dict = json.load(f)
        self.strains = [Strain(sm) for sm in self.data_dict["strains"]]
        self.pk_stresses = [Stress(d) for d in self.data_dict["pk_stresses"]]
        self.c2 = NthOrderElasticTensor.from_voigt(self.data_dict["C2_raw"])
        self.c3 = NthOrderElasticTensor.from_voigt(self.data_dict["C3_raw"])

    def test_init(self):
        c2 = NthOrderElasticTensor(self.c2.tolist())
        c3 = NthOrderElasticTensor(self.c3.tolist())
        c4 = NthOrderElasticTensor(np.zeros([3] * 8))
        for n, c in enumerate([c2, c3, c4]):
            assert c.order == n + 2
        with pytest.raises(ValueError, match="ElasticTensor must have even rank"):
            NthOrderElasticTensor(np.zeros([3] * 5))

    def test_from_diff_fit(self):
        c3 = NthOrderElasticTensor.from_diff_fit(
            self.strains,
            self.pk_stresses,
            eq_stress=self.data_dict["eq_stress"],
            order=3,
        )
        self.assert_all_close(c3.voigt, self.data_dict["C3_raw"], decimal=2)

    def test_calculate_stress(self):
        calc_stress = self.c2.calculate_stress(self.strains[0])
        self.assert_all_close(self.pk_stresses[0], calc_stress, decimal=0)
        # Test calculation from voigt strain
        self.c2.calculate_stress(self.strains[0].voigt)

    def test_energy_density(self):
        self.c3.energy_density(self.strains[0])


class TestDiffFit(PymatgenTest):
    """Tests various functions related to diff fitting."""

    def setUp(self):
        with open(f"{TEST_FILES_DIR}/test_toec_data.json") as f:
            self.data_dict = json.load(f)
        self.strains = [Strain(sm) for sm in self.data_dict["strains"]]
        self.pk_stresses = [Stress(d) for d in self.data_dict["pk_stresses"]]

    def test_get_strain_state_dict(self):
        strain_inds = [(0,), (1,), (2,), (1, 3), (1, 2, 3)]
        vecs = {}
        strain_states = []
        for strain_ind in strain_inds:
            ss = np.zeros(6)
            np.put(ss, strain_ind, 1)
            strain_states.append(tuple(ss))
            vec = np.zeros((4, 6))
            rand_values = np.random.uniform(0.1, 1, 4)
            for i in strain_ind:
                vec[:, i] = rand_values
            vecs[strain_ind] = vec
        all_strains = [Strain.from_voigt(v).zeroed() for vec in vecs.values() for v in vec]
        random.shuffle(all_strains)
        all_stresses = [Stress.from_voigt(np.random.random(6)).zeroed() for s in all_strains]
        strain_dict = {k.tobytes(): v for k, v in zip(all_strains, all_stresses)}
        ss_dict = get_strain_state_dict(all_strains, all_stresses, add_eq=False)
        # Check length of ss_dict
        assert len(strain_inds) == len(ss_dict)
        # Check sets of strain states are correct
        assert set(strain_states) == set(ss_dict)
        for data in ss_dict.values():
            # Check correspondence of strains/stresses
            for strain, stress in zip(data["strains"], data["stresses"]):
                self.assert_all_close(
                    Stress.from_voigt(stress),
                    strain_dict[Strain.from_voigt(strain).tobytes()],
                )
        # Add test to ensure zero strain state doesn't cause issue
        strains, stresses = [Strain.from_voigt([-0.01] + [0] * 5)], [Stress(np.eye(3))]
        ss_dict = get_strain_state_dict(strains, stresses)
        self.assert_all_close(list(ss_dict), [[1, 0, 0, 0, 0, 0]])

    def test_find_eq_stress(self):
        test_strains = deepcopy(self.strains)
        test_stresses = deepcopy(self.pk_stresses)
        no_eq = find_eq_stress(test_strains, test_stresses)
        self.assert_all_close(no_eq, np.zeros((3, 3)))
        test_strains[3] = Strain.from_voigt(np.zeros(6))
        eq_stress = find_eq_stress(test_strains, test_stresses)
        self.assert_all_close(test_stresses[3], eq_stress)

    def test_get_diff_coeff(self):
        forward_11 = get_diff_coeff([0, 1], 1)
        forward_13 = get_diff_coeff([0, 1, 2, 3], 1)
        backward_26 = get_diff_coeff(np.arange(-6, 1), 2)
        central_29 = get_diff_coeff(np.arange(-4, 5), 2)
        self.assert_all_close(forward_11, [-1, 1])
        self.assert_all_close(forward_13, [-11 / 6, 3, -3 / 2, 1 / 3])
        self.assert_all_close(
            backward_26,
            [137 / 180, -27 / 5, 33 / 2, -254 / 9, 117 / 4, -87 / 5, 203 / 45],
        )
        self.assert_all_close(central_29, central_diff_weights(9, 2))

    def test_generate_pseudo(self):
        strain_states = np.eye(6).tolist()
        m2, abs = generate_pseudo(strain_states, order=2)
        m3, abs = generate_pseudo(strain_states, order=3)
        m4, abs = generate_pseudo(strain_states, order=4)

    def test_fit(self):
        diff_fit(self.strains, self.pk_stresses, self.data_dict["eq_stress"])
        reduced = [(e, pk) for e, pk in zip(self.strains, self.pk_stresses) if not (abs(abs(e) - 0.05) < 1e-10).any()]
        # Get reduced dataset
        r_strains, r_pk_stresses = zip(*reduced)
        c2 = diff_fit(r_strains, r_pk_stresses, self.data_dict["eq_stress"], order=2)
        c2, c3, c4 = diff_fit(r_strains, r_pk_stresses, self.data_dict["eq_stress"], order=4)
        c2, c3 = diff_fit(self.strains, self.pk_stresses, self.data_dict["eq_stress"], order=3)
        c2_red, c3_red = diff_fit(r_strains, r_pk_stresses, self.data_dict["eq_stress"], order=3)
        self.assert_all_close(c2.voigt, self.data_dict["C2_raw"])
        self.assert_all_close(c3.voigt, self.data_dict["C3_raw"], decimal=5)
        self.assert_all_close(c2, c2_red, decimal=0)
        self.assert_all_close(c3, c3_red, decimal=-1)
