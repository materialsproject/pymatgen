import unittest
import math

from monty.serialization import MontyDecoder
from pymatgen.core.tensors import *
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class TensorTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.vec = Tensor([1., 0., 0.])
        self.rand_rank2 = Tensor(np.random.randn(3, 3))
        self.rand_rank3 = Tensor(np.random.randn(3, 3, 3))
        self.rand_rank4 = Tensor(np.random.randn(3, 3, 3, 3))
        a = 3.14 * 42.5 / 180
        self.non_symm = SquareTensor([[0.1, 0.2, 0.3],
                                      [0.4, 0.5, 0.6],
                                      [0.2, 0.5, 0.5]])
        self.rotation = SquareTensor([[math.cos(a), 0, math.sin(a)],
                                      [0, 1, 0],
                                      [-math.sin(a), 0, math.cos(a)]])
        self.low_val = Tensor([[1e-6, 1 + 1e-5, 1e-6],
                               [1 + 1e-6, 1e-6, 1e-6],
                               [1e-7, 1e-7, 1 + 1e-5]])
        self.symm_rank2 = Tensor([[1, 2, 3],
                                  [2, 4, 5],
                                  [3, 5, 6]])
        self.symm_rank3 = Tensor([[[1, 2, 3],
                                   [2, 4, 5],
                                   [3, 5, 6]],
                                  [[2, 4, 5],
                                   [4, 7, 8],
                                   [5, 8, 9]],
                                  [[3, 5, 6],
                                   [5, 8, 9],
                                   [6, 9, 10]]])
        self.symm_rank4 = Tensor([[[[1.2, 0.4, -0.92],
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
        self.fit_r3 = Tensor([[[0., 0., 0.03839],
                               [0., 0., 0.],
                               [0.03839, 0., 0.]],
                              [[0., 0., 0.],
                               [0., 0., 0.03839],
                               [0., 0.03839, 0.]],
                              [[6.89822, 0., 0.],
                               [0., 6.89822, 0.],
                               [0., 0., 27.4628]]])
        self.fit_r4 = Tensor([[[[157.9, 0., 0.],
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

        self.unfit4 = Tensor([[[[161.26, 0., 0.],
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
        self.ones = Tensor(np.ones((3, 3)))
        self.ieee_data = loadfn(ieee_file_path)

    def test_new(self):
        bad_2 = np.zeros((4, 4))
        bad_3 = np.zeros((4, 4, 4))
        self.assertRaises(ValueError, Tensor, bad_2)
        self.assertRaises(ValueError, Tensor, bad_3)
        self.assertEqual(self.rand_rank2.rank, 2)
        self.assertEqual(self.rand_rank3.rank, 3)
        self.assertEqual(self.rand_rank4.rank, 4)

    def test_zeroed(self):
        self.assertArrayEqual(self.low_val.zeroed(),
                              Tensor([[0, 1 + 1e-5, 0],
                                      [1 + 1e-6, 0, 0],
                                      [0, 0, 1 + 1e-5]]))
        self.assertArrayEqual(self.low_val.zeroed(tol=1e-6),
                              Tensor([[1e-6, 1 + 1e-5, 1e-6],
                                      [1 + 1e-6, 1e-6, 1e-6],
                                      [0, 0, 1 + 1e-5]]))
        self.assertArrayEqual(Tensor([[1e-6, -30, 1],
                                      [1e-7, 1, 0],
                                      [1e-8, 0, 1]]).zeroed(),
                              Tensor([[0, -30, 1],
                                      [0, 1, 0],
                                      [0, 0, 1]]))

    def test_transform(self):
        # Rank 3
        tensor = Tensor(np.arange(0, 27).reshape(3, 3, 3))
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

    def test_einsum_sequence(self):
        x = [1, 0, 0]
        test = Tensor(np.arange(0, 3 ** 4).reshape((3, 3, 3, 3)))
        self.assertArrayAlmostEqual([0, 27, 54], test.einsum_sequence([x] * 3))
        self.assertEqual(360, test.einsum_sequence([np.eye(3)] * 2))
        self.assertRaises(ValueError, test.einsum_sequence, Tensor(np.zeros(3)))

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
        for entry in self.ieee_data:
            xtal = entry['xtal']
            struct = entry['structure']
            orig = Tensor(entry['original_tensor'])
            ieee = Tensor(entry['ieee_tensor'])
            diff = np.max(abs(ieee - orig.convert_to_ieee(struct)))
            err_msg = "{} IEEE conversion failed with max diff {}. " \
                      "Numpy version: {}".format(xtal, diff, np.__version__)
            converted = orig.convert_to_ieee(struct, refine_rotation=False)
            self.assertArrayAlmostEqual(ieee, converted, err_msg=err_msg, decimal=3)
            converted_refined = orig.convert_to_ieee(struct, refine_rotation=True)
            err_msg = "{} IEEE conversion with refinement failed with max diff {}. " \
                      "Numpy version: {}".format(xtal, diff, np.__version__)
            self.assertArrayAlmostEqual(ieee, converted_refined, err_msg=err_msg,
                                        decimal=2)

    def test_structure_transform(self):
        # Test trivial case
        trivial = self.fit_r4.structure_transform(self.structure,
                                                  self.structure.copy())
        self.assertArrayAlmostEqual(trivial, self.fit_r4)

        # Test simple rotation
        rot_symm_op = SymmOp.from_axis_angle_and_translation([1, 1, 1], 55.5)
        rot_struct = self.structure.copy()
        rot_struct.apply_operation(rot_symm_op)
        rot_tensor = self.fit_r4.rotate(rot_symm_op.rotation_matrix)
        trans_tensor = self.fit_r4.structure_transform(self.structure, rot_struct)
        self.assertArrayAlmostEqual(rot_tensor, trans_tensor)

        # Test supercell
        bigcell = self.structure.copy()
        bigcell.make_supercell([2, 2, 3])
        trans_tensor = self.fit_r4.structure_transform(self.structure, bigcell)
        self.assertArrayAlmostEqual(self.fit_r4, trans_tensor)

        # Test rotated primitive to conventional for fcc structure
        sn = self.get_structure("Sn")
        sn_prim = SpacegroupAnalyzer(sn).get_primitive_standard_structure()
        sn_prim.apply_operation(rot_symm_op)
        rotated = self.fit_r4.rotate(rot_symm_op.rotation_matrix)
        transformed = self.fit_r4.structure_transform(sn, sn_prim)
        self.assertArrayAlmostEqual(rotated, transformed)

    def test_from_voigt(self):
        with self.assertRaises(ValueError):
            Tensor.from_voigt([[59.33, 28.08, 28.08, 0],
                               [28.08, 59.31, 28.07, 0],
                               [28.08, 28.07, 59.32, 0, 0],
                               [0, 0, 0, 26.35, 0],
                               [0, 0, 0, 0, 26.35]])
        # Rank 4
        Tensor.from_voigt([[59.33, 28.08, 28.08, 0, 0, 0],
                           [28.08, 59.31, 28.07, 0, 0, 0],
                           [28.08, 28.07, 59.32, 0, 0, 0],
                           [0, 0, 0, 26.35, 0, 0],
                           [0, 0, 0, 0, 26.35, 0],
                           [0, 0, 0, 0, 0, 26.35]])
        # Rank 3
        Tensor.from_voigt(np.zeros((3, 6)))
        # Rank 2
        Tensor.from_voigt(np.zeros(6))
        # Addresses occasional cast issues for integers
        Tensor.from_voigt(np.arange(6))

    def test_symmetry_reduce(self):
        tbs = [Tensor.from_voigt(row) for row in np.eye(6) * 0.01]
        reduced = symmetry_reduce(tbs, self.get_structure("Sn"))
        self.assertEqual(len(reduced), 2)
        self.assertArrayEqual([len(i) for i in reduced.values()], [2, 2])
        reconstructed = []
        for k, v in reduced.items():
            reconstructed.extend([k.voigt] + [k.transform(op).voigt for op in v])
        reconstructed = sorted(reconstructed, key=lambda x: np.argmax(x))
        self.assertArrayAlmostEqual([tb for tb in reconstructed], np.eye(6) * 0.01)

    def test_tensor_mapping(self):
        # Test get
        tbs = [Tensor.from_voigt(row) for row in np.eye(6) * 0.01]
        reduced = symmetry_reduce(tbs, self.get_structure("Sn"))
        tkey = Tensor.from_values_indices([0.01], [(0, 0)])
        tval = reduced[tkey]
        for tens_1, tens_2 in zip(tval, reduced[tbs[0]]):
            self.assertAlmostEqual(tens_1, tens_2)
        # Test set
        reduced[tkey] = 'test_val'
        self.assertEqual(reduced[tkey], 'test_val')
        # Test empty initialization
        empty = TensorMapping()
        self.assertEqual(empty._tensor_list, [])

    def test_populate(self):
        test_data = loadfn(os.path.join(test_dir, 'test_toec_data.json'))

        sn = self.get_structure("Sn")
        vtens = np.zeros((6, 6))
        vtens[0, 0] = 259.31
        vtens[0, 1] = 160.71
        vtens[3, 3] = 73.48
        et = Tensor.from_voigt(vtens)
        populated = et.populate(sn, prec=1e-3).voigt.round(2)
        self.assertAlmostEqual(populated[1, 1], 259.31)
        self.assertAlmostEqual(populated[2, 2], 259.31)
        self.assertAlmostEqual(populated[0, 2], 160.71)
        self.assertAlmostEqual(populated[1, 2], 160.71)
        self.assertAlmostEqual(populated[4, 4], 73.48)
        self.assertAlmostEqual(populated[5, 5], 73.48)
        # test a rank 6 example
        vtens = np.zeros([6] * 3)
        indices = [(0, 0, 0), (0, 0, 1), (0, 1, 2),
                   (0, 3, 3), (0, 5, 5), (3, 4, 5)]
        values = [-1271., -814., -50., -3., -780., -95.]
        for v, idx in zip(values, indices):
            vtens[idx] = v
        toec = Tensor.from_voigt(vtens)
        toec = toec.populate(sn, prec=1e-3, verbose=True)
        self.assertAlmostEqual(toec.voigt[1, 1, 1], -1271)
        self.assertAlmostEqual(toec.voigt[0, 1, 1], -814)
        self.assertAlmostEqual(toec.voigt[0, 2, 2], -814)
        self.assertAlmostEqual(toec.voigt[1, 4, 4], -3)
        self.assertAlmostEqual(toec.voigt[2, 5, 5], -3)
        self.assertAlmostEqual(toec.voigt[1, 2, 0], -50)
        self.assertAlmostEqual(toec.voigt[4, 5, 3], -95)

        et = Tensor.from_voigt(test_data["C3_raw"]).fit_to_structure(sn)
        new = np.zeros(et.voigt.shape)
        for idx in indices:
            new[idx] = et.voigt[idx]
        new = Tensor.from_voigt(new).populate(sn)
        self.assertArrayAlmostEqual(new, et, decimal=2)

    def test_from_values_indices(self):
        sn = self.get_structure("Sn")
        indices = [(0, 0), (0, 1), (3, 3)]
        values = [259.31, 160.71, 73.48]
        et = Tensor.from_values_indices(values, indices, structure=sn,
                                        populate=True).voigt.round(4)
        self.assertAlmostEqual(et[1, 1], 259.31)
        self.assertAlmostEqual(et[2, 2], 259.31)
        self.assertAlmostEqual(et[0, 2], 160.71)
        self.assertAlmostEqual(et[1, 2], 160.71)
        self.assertAlmostEqual(et[4, 4], 73.48)
        self.assertAlmostEqual(et[5, 5], 73.48)

    def test_serialization(self):
        # Test base serialize-deserialize
        d = self.symm_rank2.as_dict()
        new = Tensor.from_dict(d)
        self.assertArrayAlmostEqual(new, self.symm_rank2)

        d = self.symm_rank3.as_dict(voigt=True)
        new = Tensor.from_dict(d)
        self.assertArrayAlmostEqual(new, self.symm_rank3)

    def test_projection_methods(self):
        self.assertAlmostEqual(self.rand_rank2.project([1, 0, 0]),
                               self.rand_rank2[0, 0])
        self.assertAlmostEqual(self.rand_rank2.project([1, 1, 1]),
                               np.sum(self.rand_rank2) / 3)
        # Test integration
        self.assertArrayAlmostEqual(self.ones.average_over_unit_sphere(), 1)

    def test_summary_methods(self):
        self.assertEqual(set(self.ones.get_grouped_indices()[0]),
                         set(itertools.product(range(3), range(3))))
        self.assertEqual(self.ones.get_grouped_indices(voigt=True)[0],
                         [(i,) for i in range(6)])
        self.assertEqual(self.ones.get_symbol_dict(),
                         {"T_1": 1})
        self.assertEqual(self.ones.get_symbol_dict(voigt=False),
                         {"T_11": 1})

    def test_round(self):
        test = self.non_symm + 0.01
        rounded = test.round(1)
        self.assertArrayAlmostEqual(rounded, self.non_symm)
        self.assertTrue(isinstance(rounded, Tensor))


class TensorCollectionTest(PymatgenTest):
    def setUp(self):
        self.seq_tc = [t for t in np.arange(4 * 3 ** 3).reshape((4, 3, 3, 3))]
        self.seq_tc = TensorCollection(self.seq_tc)
        self.rand_tc = TensorCollection([t for t in np.random.random((4, 3, 3))])
        self.diff_rank = TensorCollection([np.ones([3] * i) for i in range(2, 5)])
        self.struct = self.get_structure("Si")
        ieee_file_path = os.path.join(test_dir, "ieee_conversion_data.json")
        self.ieee_data = loadfn(ieee_file_path)

    def list_based_function_check(self, attribute, coll, *args, **kwargs):
        """
        This function allows for more efficient testing of list-based
        functions in a "collection"-style class like TensorCollection

        It ensures that the test function
        """
        tc_orig = TensorCollection(coll)
        tc_mod = getattr(tc_orig, attribute)
        if callable(tc_mod):
            tc_mod = tc_mod(*args, **kwargs)
        for t_orig, t_mod in zip(tc_orig, tc_mod):
            this_mod = getattr(t_orig, attribute)
            if callable(this_mod):
                this_mod = this_mod(*args, **kwargs)
            if isinstance(this_mod, np.ndarray):
                self.assertArrayAlmostEqual(this_mod, t_mod)

    def test_list_based_functions(self):
        # zeroed
        tc = TensorCollection([1e-4 * Tensor(np.eye(3))] * 4)
        for t in tc.zeroed():
            self.assertArrayEqual(t, np.zeros((3, 3)))
        for t in tc.zeroed(1e-5):
            self.assertArrayEqual(t, 1e-4 * np.eye(3))
        self.list_based_function_check("zeroed", tc)
        self.list_based_function_check("zeroed", tc, tol=1e-5)

        # transform
        symm_op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30,
                                                         False, [0, 0, 1])
        self.list_based_function_check("transform", self.seq_tc, symm_op=symm_op)

        # symmetrized
        self.list_based_function_check("symmetrized", self.seq_tc)

        # rotation
        a = 3.14 * 42.5 / 180
        rotation = SquareTensor([[math.cos(a), 0, math.sin(a)], [0, 1, 0],
                                 [-math.sin(a), 0, math.cos(a)]])
        self.list_based_function_check("rotate", self.diff_rank, matrix=rotation)

        # is_symmetric
        self.assertFalse(self.seq_tc.is_symmetric())
        self.assertTrue(self.diff_rank.is_symmetric())

        # fit_to_structure
        self.list_based_function_check("fit_to_structure", self.diff_rank, self.struct)
        self.list_based_function_check("fit_to_structure", self.seq_tc, self.struct)

        # fit_to_structure
        self.list_based_function_check("fit_to_structure", self.diff_rank, self.struct)
        self.list_based_function_check("fit_to_structure", self.seq_tc, self.struct)

        # voigt
        self.list_based_function_check("voigt", self.diff_rank)

        # is_voigt_symmetric
        self.assertTrue(self.diff_rank.is_voigt_symmetric())
        self.assertFalse(self.seq_tc.is_voigt_symmetric())

        # Convert to ieee
        for entry in self.ieee_data[:2]:
            xtal = entry['xtal']
            tc = TensorCollection([entry['original_tensor']] * 3)
            struct = entry['structure']
            self.list_based_function_check("convert_to_ieee", tc, struct)

        # from_voigt
        tc_input = [t for t in np.random.random((3, 6, 6))]
        tc = TensorCollection.from_voigt(tc_input)
        for t_input, t in zip(tc_input, tc):
            self.assertArrayAlmostEqual(Tensor.from_voigt(t_input), t)

    def test_serialization(self):
        # Test base serialize-deserialize
        d = self.seq_tc.as_dict()
        new = TensorCollection.from_dict(d)
        for t, t_new in zip(self.seq_tc, new):
            self.assertArrayAlmostEqual(t, t_new)

        # Suppress vsym warnings and test voigt
        with warnings.catch_warnings(record=True):
            vsym = self.rand_tc.voigt_symmetrized
            d = vsym.as_dict(voigt=True)
            new_vsym = TensorCollection.from_dict(d)
            for t, t_new in zip(vsym, new_vsym):
                self.assertArrayAlmostEqual(t, t_new)


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
        too_high_rank = np.zeros((3, 3, 3))
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
        i2 = (self.rand_sqtensor[0, 0] * self.rand_sqtensor[1, 1] +
              self.rand_sqtensor[1, 1] * self.rand_sqtensor[2, 2] +
              self.rand_sqtensor[2, 2] * self.rand_sqtensor[0, 0] -
              self.rand_sqtensor[0, 1] * self.rand_sqtensor[1, 0] -
              self.rand_sqtensor[0, 2] * self.rand_sqtensor[2, 0] -
              self.rand_sqtensor[2, 1] * self.rand_sqtensor[1, 2])
        i3 = np.linalg.det(self.rand_sqtensor)
        self.assertArrayAlmostEqual([i1, i2, i3],
                                    self.rand_sqtensor.principal_invariants)

    def test_is_rotation(self):
        self.assertTrue(self.rotation.is_rotation())
        self.assertFalse(self.symm_sqtensor.is_rotation())
        self.assertTrue(self.low_val_2.is_rotation())
        self.assertFalse(self.low_val_2.is_rotation(tol=1e-8))

    def test_refine_rotation(self):
        self.assertArrayAlmostEqual(self.rotation, self.rotation.refine_rotation())
        new = self.rotation.copy()
        new[2, 2] += 0.02
        self.assertFalse(new.is_rotation())
        self.assertArrayAlmostEqual(self.rotation, new.refine_rotation())
        new[1] *= 1.05
        self.assertArrayAlmostEqual(self.rotation, new.refine_rotation())

    def test_get_scaled(self):
        self.assertArrayEqual(self.non_symm.get_scaled(10.),
                              SquareTensor([[1, 2, 3], [4, 5, 6], [2, 5, 5]]))

    def test_polar_decomposition(self):
        u, p = self.rand_sqtensor.polar_decomposition()
        self.assertArrayAlmostEqual(np.dot(u, p), self.rand_sqtensor)
        self.assertArrayAlmostEqual(np.eye(3),
                                    np.dot(u, np.conjugate(np.transpose(u))))

    def test_serialization(self):
        # Test base serialize-deserialize
        d = self.rand_sqtensor.as_dict()
        new = SquareTensor.from_dict(d)
        self.assertArrayAlmostEqual(new, self.rand_sqtensor)
        self.assertIsInstance(new, SquareTensor)

        # Ensure proper object-independent deserialization
        obj = MontyDecoder().process_decoded(d)
        self.assertIsInstance(obj, SquareTensor)

        with warnings.catch_warnings(record=True):
            vsym = self.rand_sqtensor.voigt_symmetrized
            d_vsym = vsym.as_dict(voigt=True)
            new_voigt = Tensor.from_dict(d_vsym)
            self.assertArrayAlmostEqual(vsym, new_voigt)


if __name__ == '__main__':
    unittest.main()
