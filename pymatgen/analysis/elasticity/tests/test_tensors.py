from __future__ import absolute_import

import unittest2 as unittest
import math
import json
import os

import numpy as np

from pymatgen.analysis.elasticity.tensors import TensorBase, SquareTensor
from pymatgen.core.operations import SymmOp
from pymatgen.util.testing import PymatgenTest
from pymatgen import Structure

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

class TensorBaseTest(PymatgenTest):
    def setUp(self):
        self.vec = TensorBase([1., 0., 0.])
        self.rand_rank2 = TensorBase(np.random.randn(3,3))
        self.rand_rank3 = TensorBase(np.random.randn(3,3,3))
        self.rand_rank4 = TensorBase(np.random.randn(3,3,3,3))
        a = 3.14 * 42.5 / 180
        self.non_symm = SquareTensor([[0.1, 0.2, 0.3],
                                      [0.4, 0.5, 0.6],
                                      [0.2, 0.5, 0.5]])
        self.rotation = SquareTensor([[math.cos(a), 0, math.sin(a)],
                                      [0, 1, 0],
                                      [-math.sin(a), 0, math.cos(a)]])
        self.low_val = TensorBase([[1e-6, 1 + 1e-5, 1e-6],
                                   [1 + 1e-6, 1e-6, 1e-6],
                                   [1e-7, 1e-7, 1 + 1e-5]])
        self.symm_rank2 = TensorBase([[1, 2, 3],
                                      [2, 4, 5],
                                      [3, 5, 6]])
        self.symm_rank3 = TensorBase([[[1, 2, 3],
                                       [2, 4, 5],
                                       [3, 5, 6]],
                                      [[2, 4, 5],
                                       [4, 7, 8],
                                       [5, 8, 9]],
                                      [[3, 5, 6],
                                       [5, 8, 9],
                                       [6, 9, 10]]])
        self.symm_rank4 = TensorBase([[[[1.2, 0.4, -0.92],
                                        [0.4, 0.05, 0.11],
                                        [-0.92, 0.11, -0.02]],
                                       [[0.4, 0.05, 0.11],
                                        [0.05, -0.47, 0.09],
                                        [0.11, 0.09, -0.]],
                                       [[-0.92, 0.11, -0.02],
                                        [0.11, 0.09, 0.],
                                        [-0.02, 0., -0.3]]],
                                      [[[0.4, 0.05, 0.11],
                                        [0.05, -0.47, 0.09],
                                        [0.11, 0.09, 0.]],
                                       [[0.05, -0.47, 0.09],
                                        [-0.47, 0.17, 0.62],
                                        [0.09, 0.62, 0.3]],
                                       [[0.11, 0.09, 0.],
                                        [0.09, 0.62, 0.3],
                                        [0., 0.3, -0.18]]],
                                      [[[-0.92, 0.11, -0.02],
                                        [0.11, 0.09, 0.],
                                        [-0.02, 0, -0.3]],
                                       [[0.11, 0.09, 0.],
                                        [0.09, 0.62, 0.3],
                                        [0., 0.3, -0.18]],
                                       [[-0.02, 0., -0.3],
                                        [0., 0.3, -0.18],
                                        [-0.3, -0.18, -0.51]]]])

        # Structural symmetries tested using BaNiO3 piezo/elastic tensors
        self.fit_r3 = TensorBase([[[0., 0., 0.03839],
                                 [0., 0., 0.],
                                 [0.03839, 0., 0.]],
                                [[0., 0., 0.],
                                 [0., 0., 0.03839],
                                 [0., 0.03839, 0.]],
                                [[6.89822, 0., 0.],
                                 [0., 6.89822, 0.],
                                 [0., 0., 27.4628]]])
        self.fit_r4 = TensorBase([[[[157.9, 0., 0.],
                                    [0., 63.1, 0.],
                                    [0., 0., 29.4]],
                                   [[0., 47.4, 0.],
                                    [47.4, 0., 0.],
                                    [0., 0., 0.]],
                                   [[0., 0., 4.3],
                                    [0., 0., 0.],
                                    [4.3, 0., 0.]]],
                                  [[[0., 47.4, 0.],
                                    [47.4, 0., 0.],
                                    [0., 0., 0.]],
                                   [[63.1, 0., 0.],
                                    [0., 157.9, 0.],
                                    [0., 0., 29.4]],
                                   [[0., 0., 0.],
                                    [0., 0., 4.3],
                                    [0., 4.3, 0.]]],
                                  [[[0., 0., 4.3],
                                    [0., 0., 0.],
                                    [4.3, 0., 0.]],
                                   [[0., 0., 0.],
                                    [0., 0., 4.3],
                                    [0., 4.3, 0.]],
                                   [[29.4, 0., 0.],
                                    [0., 29.4, 0.],
                                    [0., 0., 207.6]]]])
        
        self.unfit4 = TensorBase([[[[161.26, 0., 0.],
                                    [0., 62.76, 0.],
                                    [0., 0., 30.18]],
                                   [[0., 47.08, 0.],
                                    [47.08, 0., 0.],
                                    [0., 0., 0.]],
                                   [[0., 0., 4.23],
                                    [0., 0., 0.],
                                    [4.23, 0., 0.]]],
                                  [[[0., 47.08, 0.],
                                    [47.08, 0., 0.],
                                    [0., 0., 0.]],
                                   [[62.76, 0., 0.],
                                    [0., 155.28, -0.06],
                                    [0., -0.06, 28.53]],
                                   [[0., 0., 0.],
                                    [0., -0.06, 4.44],
                                    [0., 4.44, 0.]]],
                                  [[[0., 0., 4.23],
                                    [0., 0., 0.],
                                    [4.23, 0., 0.]],
                                   [[0., 0., 0.],
                                    [0., -0.06, 4.44],
                                    [0., 4.44, 0.]],
                                   [[30.18, 0., 0.],
                                    [0., 28.53, 0.],
                                    [0., 0., 207.57]]]])

        self.structure = self.get_structure('BaNiO3')
        ieee_file_path = os.path.join(test_dir, "ieee_conversion_data.json")
        with open(ieee_file_path) as f:
            self.ieee_data = json.load(f)

    def test_new(self):
        bad_2 = np.zeros((4, 4))
        bad_3 = np.zeros((4, 4, 4))
        self.assertRaises(ValueError, TensorBase, bad_2)
        self.assertRaises(ValueError, TensorBase, bad_3)
        self.assertEqual(self.rand_rank2.rank, 2)
        self.assertEqual(self.rand_rank3.rank, 3)
        self.assertEqual(self.rand_rank4.rank, 4)

    def test_zeroed(self):
        self.assertArrayEqual(self.low_val.zeroed(),
                              TensorBase([[0, 1 + 1e-5, 0],
                                          [1 + 1e-6, 0, 0],
                                          [0, 0, 1 + 1e-5]]))
        self.assertArrayEqual(self.low_val.zeroed(tol=1e-6),
                              TensorBase([[1e-6, 1 + 1e-5, 1e-6],
                                          [1 + 1e-6, 1e-6, 1e-6],
                                          [0, 0, 1 + 1e-5]]))
        self.assertArrayEqual(TensorBase([[1e-6, -30, 1],
                                          [1e-7, 1, 0],
                                          [1e-8, 0, 1]]).zeroed(),
                              TensorBase([[0, -30, 1],
                                          [0, 1, 0],
                                          [0, 0, 1]]))

    def test_transform(self):
        # Rank 3
        tensor = TensorBase(np.arange(0, 27).reshape(3, 3, 3))
        symm_op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30,
                                                         False, [0, 0, 1])
        new_tensor = tensor.transform(symm_op)

        self.assertArrayAlmostEqual(new_tensor,
                                    [[[-0.871, -2.884, -1.928],
                                      [-2.152, -6.665, -4.196],
                                      [-1.026, -2.830, -1.572]],
                                     [[0.044, 1.531, 1.804],
                                      [4.263, 21.008, 17.928],
                                      [5.170, 23.026, 18.722]],
                                     [[1.679, 7.268, 5.821],
                                      [9.268, 38.321, 29.919],
                                      [8.285, 33.651, 26.000]]], 3)

    
    def test_rotate(self):
        self.assertArrayEqual(self.vec.rotate([[0, -1, 0],
                                               [1, 0, 0],
                                               [0, 0, 1]]),
                              [0, 1, 0])
        self.assertArrayAlmostEqual(self.non_symm.rotate(self.rotation),
                                    SquareTensor([[0.531, 0.485, 0.271],
                                                  [0.700, 0.5, 0.172],
                                                  [0.171, 0.233, 0.068]]),
                                    decimal=3)
        self.assertRaises(ValueError, self.non_symm.rotate, 
                          self.symm_rank2)

    def test_symmetrized(self):
        self.assertTrue(self.rand_rank2.symmetrized.is_symmetric())
        self.assertTrue(self.rand_rank3.symmetrized.is_symmetric())
        self.assertTrue(self.rand_rank4.symmetrized.is_symmetric())
    
    
    def test_is_symmetric(self):
        self.assertTrue(self.symm_rank2.is_symmetric())
        self.assertTrue(self.symm_rank3.is_symmetric())
        self.assertTrue(self.symm_rank4.is_symmetric())
        tol_test = self.symm_rank4
        tol_test[0, 1, 2, 2] += 1e-6
        self.assertFalse(self.low_val.is_symmetric(tol=1e-8))

    def test_fit_to_structure(self):
        new_fit = self.unfit4.fit_to_structure(self.structure)
        self.assertArrayAlmostEqual(new_fit, self.fit_r4, 1)

    def test_is_fit_to_structure(self):
        self.assertFalse(self.unfit4.is_fit_to_structure(self.structure))
        self.assertTrue(self.fit_r3.is_fit_to_structure(self.structure))
        self.assertTrue(self.fit_r4.is_fit_to_structure(self.structure))

    def test_convert_to_ieee(self):
        for xtal in self.ieee_data.keys():
            orig = TensorBase(self.ieee_data[xtal]['original_tensor'])
            ieee = TensorBase(self.ieee_data[xtal]['ieee_tensor'])
            struct = Structure.from_dict(self.ieee_data[xtal]['structure'])
            diff = np.max(abs(ieee - orig.convert_to_ieee(struct)))
            err_msg = "{} IEEE conversion failed with max diff {}".format(
                xtal, diff)
            self.assertArrayAlmostEqual(ieee, orig.convert_to_ieee(struct),
                                        err_msg = err_msg, decimal=3)

class SquareTensorTest(PymatgenTest):
    def setUp(self):
        self.rand_sqtensor = SquareTensor(np.random.randn(3, 3))
        self.symm_sqtensor = SquareTensor([[0.1, 0.3, 0.4],
                                           [0.3, 0.5, 0.2],
                                           [0.4, 0.2, 0.6]])
        self.non_invertible = SquareTensor([[0.1, 0, 0],
                                            [0.2, 0, 0],
                                            [0, 0, 0]])
        self.non_symm = SquareTensor([[0.1, 0.2, 0.3],
                                      [0.4, 0.5, 0.6],
                                      [0.2, 0.5, 0.5]])
        self.low_val = SquareTensor([[1e-6, 1 + 1e-5, 1e-6],
                                     [1 + 1e-6, 1e-6, 1e-6],
                                     [1e-7, 1e-7, 1 + 1e-5]])
        self.low_val_2 = SquareTensor([[1e-6, -1 - 1e-6, 1e-6],
                                       [1 + 1e-7, 1e-6, 1e-6],
                                       [1e-7, 1e-7, 1 + 1e-6]])
        a = 3.14 * 42.5 / 180
        self.rotation = SquareTensor([[math.cos(a), 0, math.sin(a)],
                                      [0, 1, 0],
                                      [-math.sin(a), 0, math.cos(a)]])
        
    def test_new(self):
        non_sq_matrix = [[0.1, 0.2, 0.1],
                         [0.1, 0.2, 0.3],
                         [0.1, 0.2, 0.3],
                         [0.1, 0.1, 0.1]]
        bad_matrix = [[0.1, 0.2],
                      [0.2, 0.3, 0.4],
                      [0.2, 0.3, 0.5]]
        too_high_rank = np.zeros((3,3,3))
        self.assertRaises(ValueError, SquareTensor, non_sq_matrix)
        self.assertRaises(ValueError, SquareTensor, bad_matrix)
        self.assertRaises(ValueError, SquareTensor, too_high_rank)

    def test_properties(self):
        # transpose
        self.assertArrayEqual(self.non_symm.trans, 
                              SquareTensor([[0.1, 0.4, 0.2],
                                            [0.2, 0.5, 0.5],
                                            [0.3, 0.6, 0.5]]))
        self.assertArrayEqual(self.rand_sqtensor.trans,
                              np.transpose(self.rand_sqtensor))
        self.assertArrayEqual(self.symm_sqtensor,
                              self.symm_sqtensor.trans)
        # inverse
        self.assertArrayEqual(self.non_symm.inv,
                              np.linalg.inv(self.non_symm))
        with self.assertRaises(ValueError):
            self.non_invertible.inv

        # determinant
        self.assertEqual(self.rand_sqtensor.det,
                         np.linalg.det(self.rand_sqtensor))
        self.assertEqual(self.non_invertible.det,
                         0.0)
        self.assertEqual(self.non_symm.det, 0.009)

        # symmetrized
        self.assertArrayEqual(self.rand_sqtensor.symmetrized,
                              0.5 * (self.rand_sqtensor + self.rand_sqtensor.trans))
        self.assertArrayEqual(self.symm_sqtensor,
                              self.symm_sqtensor.symmetrized)
        self.assertArrayAlmostEqual(self.non_symm.symmetrized,
                                    SquareTensor([[0.1, 0.3, 0.25],
                                                  [0.3, 0.5, 0.55],
                                                  [0.25, 0.55, 0.5]]))

        # invariants
        i1 = np.trace(self.rand_sqtensor)
        i2 = self.rand_sqtensor[0, 0] * self.rand_sqtensor[1, 1] + \
             self.rand_sqtensor[1, 1] * self.rand_sqtensor[2, 2] + \
             self.rand_sqtensor[2, 2] * self.rand_sqtensor[0, 0] - \
             self.rand_sqtensor[0, 1] * self.rand_sqtensor[1, 0] - \
             self.rand_sqtensor[0, 2] * self.rand_sqtensor[2, 0] - \
             self.rand_sqtensor[2, 1] * self.rand_sqtensor[1, 2]
        i3 = np.linalg.det(self.rand_sqtensor)
        self.assertArrayAlmostEqual([i1, i2, i3],
                                    self.rand_sqtensor.principal_invariants)

    def test_is_rotation(self):
        self.assertTrue(self.rotation.is_rotation())
        self.assertFalse(self.symm_sqtensor.is_rotation())
        self.assertTrue(self.low_val_2.is_rotation())
        self.assertFalse(self.low_val_2.is_rotation(tol=1e-8))

    def test_get_scaled(self):
        self.assertArrayEqual(self.non_symm.get_scaled(10.),
                              SquareTensor([[1, 2, 3], [4, 5, 6], [2, 5, 5]]))

    def test_polar_decomposition(self):
        u, p = self.rand_sqtensor.polar_decomposition()
        self.assertArrayAlmostEqual(np.dot(u, p), self.rand_sqtensor)
        self.assertArrayAlmostEqual(np.eye(3),
                                    np.dot(u, np.conjugate(np.transpose(u))))

if __name__ == '__main__':
    unittest.main()
