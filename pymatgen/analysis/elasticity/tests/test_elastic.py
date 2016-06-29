from __future__ import absolute_import, print_function

import json
import os
import warnings

import numpy as np
import unittest2 as unittest

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import (Deformation,
                                                 IndependentStrain, Strain)
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.util.testing import PymatgenTest
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
        # full tensor
        self.assertArrayAlmostEqual(self.elastic_tensor_1.voigt,
                                    self.voigt_1)

        # homogenous youngs
        self.assertAlmostEqual(self.elastic_tensor_1.homogeneous_youngs,54.087787667160576)

    def test_new(self):
        self.assertArrayAlmostEqual(self.elastic_tensor_1,
                                    ElasticTensor(self.ft))

    def test_from_voigt(self):
        with self.assertRaises(ValueError):
            ElasticTensor.from_voigt([[59.33, 28.08, 28.08, 0],
                                      [28.08, 59.31, 28.07, 0],
                                      [28.08, 28.07, 59.32, 0, 0],
                                      [0, 0, 0, 26.35, 0],
                                      [0, 0, 0, 0, 26.35]])
        with warnings.catch_warnings(record=True) as w:
            ElasticTensor.from_voigt([[59.33, 28.08, 28.08, 0, 0, 0],
                                      [0.0, 59.31, 28.07, 0, 0, 0],
                                      [28.08, 28.07, 59.32, 0, 0, 0],
                                      [0, 0, 0, 26.35, 0, 0],
                                      [0, 0, 0, 0, 26.35, 0],
                                      [0, 0, 0, 0, 0, 26.35]])
            self.assertEqual(len(w), 1)

    def test_from_strain_stress_list(self):
        strain_list = [Strain.from_deformation(def_matrix)
                       for def_matrix in self.def_stress_dict['deformations']]
        stress_list = [stress for stress in self.def_stress_dict['stresses']]
        with warnings.catch_warnings(record=True):
            et_fl = -0.1 * ElasticTensor.from_strain_stress_list(strain_list,
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
        with warnings.catch_warnings(record=True):
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

        dfm = Deformation([[-9.86004855e-01, 2.27539582e-01, -4.64426035e-17],
                           [-2.47802121e-01, -9.91208483e-01, -7.58675185e-17],
                           [-6.12323400e-17, -6.12323400e-17, 1.00000000e+00]])

        self.assertAlmostEqual(film_elac.energy_density(dfm.green_lagrange_strain),
                               0.000125664672793)

    def test_christoffel_tensor(self):
        self.assertArrayAlmostEqual(self.elastic_tensor_1.ChristoffelTensor([1, 0, 0]),
                                    [[59.33,  0.,     0.],
                                     [0.,    26.35,   0.],
                                     [0.,     0.,    26.35]])
        self.assertArrayAlmostEqual(self.elastic_tensor_1.ChristoffelTensor([1, 0, 1]),
                                    [[42.84,   0.,  27.215],
                                     [0.,  26.35,   0.],
                                     [27.215,   0.,  42.835]])

        self.assertArrayAlmostEqual(self.elastic_tensor_1.WaveVelocities([1, 0, 0]),
                                    [7.702597,  5.1332251,  5.1332251])
        self.assertArrayAlmostEqual(self.elastic_tensor_1.WaveVelocities([1, 0, 1]),
                                    [8.3697372,  3.9525308,  5.1332251])

    def test_stability(self):
        self.assertTrue(self.elastic_tensor_1.elasticically_stable)
        elac = ElasticTensor.from_voigt([[-86.,  25.,  18.,   0.,   0.,   0.],
                              [25.,  73.,  24.,   0.,   0.,   0.],
                              [18.,  24., -86.,   0.,   0.,   0.],
                              [0.,   0.,   0.,  24.,   0.,   0.],
                              [0.,   0.,   0.,   0.,  18.,   0.],
                              [0.,   0.,   0.,   0.,   0.,  25.]])

        self.assertFalse(elac.elasticically_stable)
if __name__ == '__main__':
    unittest.main()
