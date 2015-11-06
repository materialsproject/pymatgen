from __future__ import absolute_import
from __future__ import print_function

import unittest
import os
import random

import numpy as np
from pymatgen.elasticity.elastic import ElasticTensor
from pymatgen.elasticity.strain import Strain, Deformation, IndependentStrain
from pymatgen.elasticity.stress import Stress
from pymatgen.util.testing import PymatgenTest
from numpy.testing import *
import math
import warnings
import json
from six.moves import zip

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class ElasticTensorTest(PymatgenTest):
    def setUp(self):
        self.elastic_tensor_1 = ElasticTensor([[59.33, 28.08, 28.08, 0, 0, 0],
                                               [28.08, 59.31, 28.07, 0, 0, 0],
                                               [28.08, 28.07, 59.32, 0, 0, 0],
                                               [0, 0, 0, 26.35, 0, 0],
                                               [0, 0, 0, 0, 26.35, 0],
                                               [0, 0, 0, 0, 0, 26.35]])
        mat = np.random.randn(6, 6)
        mat = mat + np.transpose(mat)
        self.rand_elastic_tensor = ElasticTensor(mat)
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

        filepath = os.path.join(test_dir, 'Sn_def_stress.json')
        self.def_stress_dict = json.load(open(filepath))
        warnings.simplefilter("always")
    
    def test_new(self):
        with self.assertRaises(ValueError):
            ElasticTensor([[59.33, 28.08, 28.08, 0],
                           [28.08, 59.31, 28.07, 0],
                           [28.08, 28.07, 59.32, 0, 0],
                           [0, 0, 0, 26.35, 0],
                           [0, 0, 0, 0, 26.35]])
        with warnings.catch_warnings(record=True) as w:
            ElasticTensor([[59.33, 28.08, 28.08, 0, 0, 0],
                           [0.0, 59.31, 28.07, 0, 0, 0],
                           [28.08, 28.07, 59.32, 0, 0, 0],
                           [0, 0, 0, 26.35, 0, 0],
                           [0, 0, 0, 0, 26.35, 0],
                           [0, 0, 0, 0, 0, 26.35]])
            self.assertEqual(len(w),1)

    def test_properties(self):
        # compliance tensor
        self.assertArrayAlmostEqual(np.linalg.inv(self.elastic_tensor_1),
                                    self.elastic_tensor_1.compliance_tensor)
        # KG average properties
        et = self.elastic_tensor_1
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
        self.assertArrayAlmostEqual(self.elastic_tensor_1.full_tensor,
                                    self.ft)

    def test_from_full_tensor(self):
        self.assertArrayAlmostEqual(self.elastic_tensor_1,
                                    ElasticTensor.from_full_tensor(self.ft))

    def test_from_strain_stress_list(self):
        strain_list = [Strain.from_deformation(def_matrix) 
                       for def_matrix in self.def_stress_dict['deformations']]
        stress_list = [stress for stress in self.def_stress_dict['stresses']]
        with warnings.catch_warnings(record=True):
            et_fl = -0.1*ElasticTensor.from_strain_stress_list(strain_list,stress_list)
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
        et_from_sd = ElasticTensor.from_stress_dict(stress_dict)
        self.assertArrayAlmostEqual(et_from_sd.round(2),
                                    self.elastic_tensor_1)
if __name__ == '__main__':
    unittest.main()
