from __future__ import absolute_import
from __future__ import print_function

import unittest2 as unittest
import os

import numpy as np
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Strain, IndependentStrain, Deformation
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.util.testing import PymatgenTest
import warnings
import json
from six.moves import zip

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class ElasticTensorTest(PymatgenTest):
    def setUp(self):
        self.voigt_1 = [[59.33, 28.08, 28.08, 0, 0, 0],
                        [28.08, 59.31, 28.07, 0, 0, 0],
                        [28.08, 28.07, 59.32, 0, 0, 0],
                        [0, 0, 0, 26.35, 0, 0],
                        [0, 0, 0, 0, 26.35, 0],
                        [0, 0, 0, 0, 0, 26.35]]
        mat = np.random.randn(6, 6)
        mat = mat + np.transpose(mat)
        self.rand_elastic_tensor = ElasticTensor.from_voigt(mat)
        self.ft = np.array([[[[59.33, 0, 0],
                              [0, 28.08, 0],
                              [0, 0, 28.08]],
                             [[0, 26.35, 0],
                              [26.35, 0, 0],
                              [0, 0, 0]],
                             [[0, 0, 26.35],
                              [0, 0, 0],
                              [26.35, 0, 0]]],
                            [[[0, 26.35, 0],
                              [26.35, 0, 0],
                              [0, 0, 0]],
                             [[28.08, 0, 0],
                              [0, 59.31, 0],
                              [0, 0, 28.07]],
                             [[0, 0, 0],
                              [0, 0, 26.35],
                              [0, 26.35, 0]]],
                            [[[0, 0, 26.35],
                              [0, 0, 0],
                              [26.35, 0, 0]],
                             [[0, 0, 0],
                              [0, 0, 26.35],
                              [0, 26.35, 0]],
                             [[28.08, 0, 0],
                              [0, 28.07, 0],
                              [0, 0, 59.32]]]])

        self.elastic_tensor_1 = ElasticTensor(self.ft)
        filepath = os.path.join(test_dir, 'Sn_def_stress.json')
        with open(filepath) as f:
            self.def_stress_dict = json.load(f)
        self.structure = self.get_structure("Sn")

        warnings.simplefilter("always")

    def test_properties(self):
        # compliance tensor
        self.assertArrayAlmostEqual(np.linalg.inv(self.elastic_tensor_1.voigt),
                                    self.elastic_tensor_1.compliance_tensor)
        # KG average properties
        self.assertAlmostEqual(38.49111111111, self.elastic_tensor_1.k_voigt)
        self.assertAlmostEqual(22.05866666666, self.elastic_tensor_1.g_voigt)
        self.assertAlmostEqual(38.49110945133, self.elastic_tensor_1.k_reuss)
        self.assertAlmostEqual(20.67146635306, self.elastic_tensor_1.g_reuss)
        self.assertAlmostEqual(38.49111028122, self.elastic_tensor_1.k_vrh)
        self.assertAlmostEqual(21.36506650986, self.elastic_tensor_1.g_vrh)
        self.assertArrayAlmostEqual(self.elastic_tensor_1.kg_average,
                                    [38.49111111111,
                                     22.05866666666,
                                     38.49110945133,
                                     20.67146635306,
                                     38.49111028122,
                                     21.36506650986])
        # universal anisotropy
        self.assertAlmostEqual(0.33553509658699,
                               self.elastic_tensor_1.universal_anisotropy)
        # homogeneous poisson
        self.assertAlmostEqual(0.26579965576472,
                               self.elastic_tensor_1.homogeneous_poisson)
        # voigt notation tensor
        self.assertArrayAlmostEqual(self.elastic_tensor_1.voigt,
                                    self.voigt_1)

        # young's modulus
        self.assertAlmostEqual(54087787667.160583,
                               self.elastic_tensor_1.y_mod)

    def test_structure_based_methods(self):
        # trans_velocity
        self.assertAlmostEqual(1996.35019877,
                               self.elastic_tensor_1.trans_v(self.structure))

        # long_velocity
        self.assertAlmostEqual(3534.68123832,
                               self.elastic_tensor_1.long_v(self.structure))
        # Snyder properties
        self.assertAlmostEqual(18.06127074,
                               self.elastic_tensor_1.snyder_ac(self.structure))
        self.assertAlmostEqual(0.18937465,
                               self.elastic_tensor_1.snyder_opt(self.structure))
        self.assertAlmostEqual(18.25064540,
                               self.elastic_tensor_1.snyder_total(self.structure))
        # Clarke
        self.assertAlmostEqual(0.3450307,
                               self.elastic_tensor_1.clarke_thermalcond(self.structure))
        # Cahill
        self.assertAlmostEqual(0.37896275,
                               self.elastic_tensor_1.cahill_thermalcond(self.structure))
        # Debye
        self.assertAlmostEqual(247.3058931,
                               self.elastic_tensor_1.debye_temperature(self.structure))
        self.assertAlmostEqual(189.05670205,
                               self.elastic_tensor_1.debye_temperature_gibbs(self.structure))

    def test_new(self):
        self.assertArrayAlmostEqual(self.elastic_tensor_1,
                                    ElasticTensor(self.ft))
        nonsymm = self.ft
        nonsymm[0, 1, 2, 2] += 1.0
        with warnings.catch_warnings(record=True) as w:
            ElasticTensor(nonsymm)
            self.assertEqual(len(w), 1)
        badtensor1 = np.zeros((3, 3, 3))
        badtensor2 = np.zeros((3, 3, 3, 2))
        self.assertRaises(ValueError, ElasticTensor, badtensor1)
        self.assertRaises(ValueError, ElasticTensor, badtensor2)

    def test_from_strain_stress_list(self):
        strain_list = [Strain.from_deformation(def_matrix)
                       for def_matrix in self.def_stress_dict['deformations']]
        stress_list = [stress for stress in self.def_stress_dict['stresses']]
        with warnings.catch_warnings(record=True):
            et_fl = -0.1*ElasticTensor.from_strain_stress_list(strain_list, 
                                                               stress_list).voigt
            self.assertArrayAlmostEqual(et_fl.round(2),
                                        [[59.29, 24.36, 22.46, 0, 0, 0],
                                         [28.06, 56.91, 22.46, 0, 0, 0],
                                         [28.06, 25.98, 54.67, 0, 0, 0],
                                         [0, 0, 0, 26.35, 0, 0],
                                         [0, 0, 0, 0, 26.35, 0],
                                         [0, 0, 0, 0, 0, 26.35]])

    def test_from_stress_dict(self):
        stress_dict = dict(list(zip([IndependentStrain(def_matrix) for def_matrix
                                in self.def_stress_dict['deformations']],
                                [Stress(stress_matrix) for stress_matrix
                                in self.def_stress_dict['stresses']])))
        with warnings.catch_warnings(record = True):
            et_from_sd = ElasticTensor.from_stress_dict(stress_dict)
        self.assertArrayAlmostEqual(et_from_sd.voigt_symmetrized.round(2),
                                    self.elastic_tensor_1)

    def test_energy_density(self):

        film_elac = ElasticTensor.from_voigt([
            [324.32,  187.3,   170.92,    0.,      0.,      0.],
            [187.3,   324.32,  170.92,    0.,      0.,      0.],
            [170.92,  170.92,  408.41,    0.,      0.,      0.],
            [0.,      0.,      0.,    150.73,    0.,      0.],
            [0.,      0.,      0.,      0.,    150.73,    0.],
            [0.,      0.,      0.,      0.,      0.,    238.74]])

        dfm = Deformation([[ -9.86004855e-01,2.27539582e-01,-4.64426035e-17],
                           [ -2.47802121e-01,-9.91208483e-01,-7.58675185e-17],
                           [ -6.12323400e-17,-6.12323400e-17,1.00000000e+00]])

        self.assertAlmostEqual(film_elac.energy_density(dfm.green_lagrange_strain),
            0.000125664672793)

if __name__ == '__main__':
    unittest.main()
