from __future__ import annotations

import math

import numpy as np
import pytest
from monty.serialization import MontyDecoder, loadfn
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.core.operations import SymmOp
from pymatgen.core.tensors import SquareTensor, Tensor, TensorCollection, TensorMapping, itertools, symmetry_reduce
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestTensor(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.vec = Tensor([1.0, 0.0, 0.0])
        self.rand_rank2 = Tensor(np.random.randn(3, 3))
        self.rand_rank3 = Tensor(np.random.randn(3, 3, 3))
        self.rand_rank4 = Tensor(np.random.randn(3, 3, 3, 3))
        a = 3.14 * 42.5 / 180
        self.non_symm = SquareTensor([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.2, 0.5, 0.5]])
        self.rotation = SquareTensor([[math.cos(a), 0, math.sin(a)], [0, 1, 0], [-math.sin(a), 0, math.cos(a)]])
        self.low_val = Tensor([[1e-6, 1 + 1e-5, 1e-6], [1 + 1e-6, 1e-6, 1e-6], [1e-7, 1e-7, 1 + 1e-5]])
        self.symm_rank2 = Tensor([[1, 2, 3], [2, 4, 5], [3, 5, 6]])
        self.symm_rank3 = Tensor(
            [
                [[1, 2, 3], [2, 4, 5], [3, 5, 6]],
                [[2, 4, 5], [4, 7, 8], [5, 8, 9]],
                [[3, 5, 6], [5, 8, 9], [6, 9, 10]],
            ]
        )
        self.symm_rank4 = Tensor(
            [
                [
                    [[1.2, 0.4, -0.92], [0.4, 0.05, 0.11], [-0.92, 0.11, -0.02]],
                    [[0.4, 0.05, 0.11], [0.05, -0.47, 0.09], [0.11, 0.09, -0.0]],
                    [[-0.92, 0.11, -0.02], [0.11, 0.09, 0.0], [-0.02, 0.0, -0.3]],
                ],
                [
                    [[0.4, 0.05, 0.11], [0.05, -0.47, 0.09], [0.11, 0.09, 0.0]],
                    [[0.05, -0.47, 0.09], [-0.47, 0.17, 0.62], [0.09, 0.62, 0.3]],
                    [[0.11, 0.09, 0.0], [0.09, 0.62, 0.3], [0.0, 0.3, -0.18]],
                ],
                [
                    [[-0.92, 0.11, -0.02], [0.11, 0.09, 0.0], [-0.02, 0, -0.3]],
                    [[0.11, 0.09, 0.0], [0.09, 0.62, 0.3], [0.0, 0.3, -0.18]],
                    [[-0.02, 0.0, -0.3], [0.0, 0.3, -0.18], [-0.3, -0.18, -0.51]],
                ],
            ]
        )

        # Structural symmetries tested using BaNiO3 piezo/elastic tensors
        self.fit_r3 = Tensor(
            [
                [[0.0, 0.0, 0.03839], [0.0, 0.0, 0.0], [0.03839, 0.0, 0.0]],
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.03839], [0.0, 0.03839, 0.0]],
                [[6.89822, 0.0, 0.0], [0.0, 6.89822, 0.0], [0.0, 0.0, 27.4628]],
            ]
        )
        self.fit_r4 = Tensor(
            [
                [
                    [[157.9, 0.0, 0.0], [0.0, 63.1, 0.0], [0.0, 0.0, 29.4]],
                    [[0.0, 47.4, 0.0], [47.4, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 4.3], [0.0, 0.0, 0.0], [4.3, 0.0, 0.0]],
                ],
                [
                    [[0.0, 47.4, 0.0], [47.4, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[63.1, 0.0, 0.0], [0.0, 157.9, 0.0], [0.0, 0.0, 29.4]],
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 4.3], [0.0, 4.3, 0.0]],
                ],
                [
                    [[0.0, 0.0, 4.3], [0.0, 0.0, 0.0], [4.3, 0.0, 0.0]],
                    [[0.0, 0.0, 0.0], [0.0, 0.0, 4.3], [0.0, 4.3, 0.0]],
                    [[29.4, 0.0, 0.0], [0.0, 29.4, 0.0], [0.0, 0.0, 207.6]],
                ],
            ]
        )

        self.unfit4 = Tensor(
            [
                [
                    [[161.26, 0.0, 0.0], [0.0, 62.76, 0.0], [0.0, 0.0, 30.18]],
                    [[0.0, 47.08, 0.0], [47.08, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 4.23], [0.0, 0.0, 0.0], [4.23, 0.0, 0.0]],
                ],
                [
                    [[0.0, 47.08, 0.0], [47.08, 0.0, 0.0], [0.0, 0.0, 0.0]],
                    [[62.76, 0.0, 0.0], [0.0, 155.28, -0.06], [0.0, -0.06, 28.53]],
                    [[0.0, 0.0, 0.0], [0.0, -0.06, 4.44], [0.0, 4.44, 0.0]],
                ],
                [
                    [[0.0, 0.0, 4.23], [0.0, 0.0, 0.0], [4.23, 0.0, 0.0]],
                    [[0.0, 0.0, 0.0], [0.0, -0.06, 4.44], [0.0, 4.44, 0.0]],
                    [[30.18, 0.0, 0.0], [0.0, 28.53, 0.0], [0.0, 0.0, 207.57]],
                ],
            ]
        )

        self.structure = self.get_structure("BaNiO3")
        ieee_file_path = f"{TEST_FILES_DIR}/ieee_conversion_data.json"
        self.ones = Tensor(np.ones((3, 3)))
        self.ieee_data = loadfn(ieee_file_path)

    def test_new(self):
        bad_2 = np.zeros((4, 4))
        bad_3 = np.zeros((4, 4, 4))
        expected_msg = (
            "Pymatgen only supports 3-dimensional tensors, and default tensor constructor uses standard notation"
        )
        with pytest.raises(ValueError, match=expected_msg):
            Tensor(bad_2)
        with pytest.raises(ValueError, match=expected_msg):
            Tensor(bad_3)
        assert self.rand_rank2.rank == 2
        assert self.rand_rank3.rank == 3
        assert self.rand_rank4.rank == 4

    def test_zeroed(self):
        assert self.low_val.zeroed() == approx(Tensor([[0, 1 + 1e-5, 0], [1 + 1e-6, 0, 0], [0, 0, 1 + 1e-5]]))
        assert self.low_val.zeroed(tol=1e-6) == approx(
            Tensor([[1e-6, 1 + 1e-5, 1e-6], [1 + 1e-6, 1e-6, 1e-6], [0, 0, 1 + 1e-5]])
        )
        assert Tensor([[1e-6, -30, 1], [1e-7, 1, 0], [1e-8, 0, 1]]).zeroed() == approx(
            Tensor([[0, -30, 1], [0, 1, 0], [0, 0, 1]])
        )

    def test_transform(self):
        # Rank 3
        tensor = Tensor(np.arange(0, 27).reshape(3, 3, 3))
        symm_op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False, [0, 0, 1])
        new_tensor = tensor.transform(symm_op)

        assert np.allclose(
            new_tensor,
            [
                [
                    [-0.871, -2.884, -1.928],
                    [-2.152, -6.665, -4.196],
                    [-1.026, -2.830, -1.572],
                ],
                [
                    [0.044, 1.531, 1.804],
                    [4.263, 21.008, 17.928],
                    [5.170, 23.026, 18.722],
                ],
                [
                    [1.679, 7.268, 5.821],
                    [9.268, 38.321, 29.919],
                    [8.285, 33.651, 26.000],
                ],
            ],
            3,
        )

    def test_rotate(self):
        assert_array_equal(self.vec.rotate([[0, -1, 0], [1, 0, 0], [0, 0, 1]]), [0, 1, 0])
        assert np.allclose(
            self.non_symm.rotate(self.rotation),
            SquareTensor([[0.531, 0.485, 0.271], [0.700, 0.5, 0.172], [0.171, 0.233, 0.068]]),
            atol=1e-3,
        )
        with pytest.raises(ValueError, match="Rotation matrix is not valid"):
            self.non_symm.rotate(self.symm_rank2)

    def test_einsum_sequence(self):
        x = [1, 0, 0]
        test = Tensor(np.arange(0, 3**4).reshape((3, 3, 3, 3)))
        assert np.allclose([0, 27, 54], test.einsum_sequence([x] * 3))
        assert test.einsum_sequence([np.eye(3)] * 2) == 360
        with pytest.raises(ValueError, match="other tensors must be list of tensors or tensor input"):
            test.einsum_sequence(Tensor(np.zeros(3)))

    def test_symmetrized(self):
        assert self.rand_rank2.symmetrized.is_symmetric()
        assert self.rand_rank3.symmetrized.is_symmetric()
        assert self.rand_rank4.symmetrized.is_symmetric()

    def test_is_symmetric(self):
        assert self.symm_rank2.is_symmetric()
        assert self.symm_rank3.is_symmetric()
        assert self.symm_rank4.is_symmetric()
        tol_test = self.symm_rank4
        tol_test[0, 1, 2, 2] += 1e-6
        assert not self.low_val.is_symmetric(tol=1e-8)

    def test_fit_to_structure(self):
        new_fit = self.unfit4.fit_to_structure(self.structure)
        assert np.allclose(new_fit, self.fit_r4, atol=1e-1)

    def test_is_fit_to_structure(self):
        assert not self.unfit4.is_fit_to_structure(self.structure)
        assert self.fit_r3.is_fit_to_structure(self.structure)
        assert self.fit_r4.is_fit_to_structure(self.structure)

    def test_convert_to_ieee(self):
        for entry in self.ieee_data:
            entry["xtal"]
            struct = entry["structure"]
            orig = Tensor(entry["original_tensor"])
            ieee = Tensor(entry["ieee_tensor"])
            np.max(abs(ieee - orig.convert_to_ieee(struct)))
            converted = orig.convert_to_ieee(struct, refine_rotation=False)
            assert np.allclose(ieee, converted, atol=1e-2)
            converted_refined = orig.convert_to_ieee(struct, refine_rotation=True)
            assert np.allclose(ieee, converted_refined, atol=1e-2)

    def test_structure_transform(self):
        # Test trivial case
        trivial = self.fit_r4.structure_transform(self.structure, self.structure.copy())
        assert np.allclose(trivial, self.fit_r4)

        # Test simple rotation
        rot_symm_op = SymmOp.from_axis_angle_and_translation([1, 1, 1], 55.5)
        rot_struct = self.structure.copy()
        rot_struct.apply_operation(rot_symm_op)
        rot_tensor = self.fit_r4.rotate(rot_symm_op.rotation_matrix)
        trans_tensor = self.fit_r4.structure_transform(self.structure, rot_struct)
        assert np.allclose(rot_tensor, trans_tensor)

        # Test supercell
        bigcell = self.structure.copy()
        bigcell.make_supercell([2, 2, 3])
        trans_tensor = self.fit_r4.structure_transform(self.structure, bigcell)
        assert np.allclose(self.fit_r4, trans_tensor)

        # Test rotated primitive to conventional for fcc structure
        sn = self.get_structure("Sn")
        sn_prim = SpacegroupAnalyzer(sn).get_primitive_standard_structure()
        sn_prim.apply_operation(rot_symm_op)
        rotated = self.fit_r4.rotate(rot_symm_op.rotation_matrix)
        transformed = self.fit_r4.structure_transform(sn, sn_prim)
        assert np.allclose(rotated, transformed)

    def test_from_voigt(self):
        with pytest.raises(ValueError, match="The requested array has an inhomogeneous shape after 1 dimensions."):
            Tensor.from_voigt(
                [
                    [59.33, 28.08, 28.08, 0],
                    [28.08, 59.31, 28.07, 0],
                    [28.08, 28.07, 59.32, 0, 0],
                    [0, 0, 0, 26.35, 0],
                    [0, 0, 0, 0, 26.35],
                ]
            )
        # Rank 4
        Tensor.from_voigt(
            [
                [59.33, 28.08, 28.08, 0, 0, 0],
                [28.08, 59.31, 28.07, 0, 0, 0],
                [28.08, 28.07, 59.32, 0, 0, 0],
                [0, 0, 0, 26.35, 0, 0],
                [0, 0, 0, 0, 26.35, 0],
                [0, 0, 0, 0, 0, 26.35],
            ]
        )
        # Rank 3
        Tensor.from_voigt(np.zeros((3, 6)))
        # Rank 2
        Tensor.from_voigt(np.zeros(6))
        # Addresses occasional cast issues for integers
        Tensor.from_voigt(np.arange(6))

    def test_symmetry_reduce(self):
        tbs = [Tensor.from_voigt(row) for row in np.eye(6) * 0.01]
        reduced = symmetry_reduce(tbs, self.get_structure("Sn"))
        assert len(reduced) == 2
        assert [len(i) for i in reduced.values()] == [2, 2]
        reconstructed = []
        for k, v in reduced.items():
            reconstructed.extend([k.voigt] + [k.transform(op).voigt for op in v])
        reconstructed = sorted(reconstructed, key=lambda x: np.argmax(x))
        assert np.allclose(list(reconstructed), np.eye(6) * 0.01)

    def test_tensor_mapping(self):
        # Test get
        tbs = [Tensor.from_voigt(row) for row in np.eye(6) * 0.01]
        reduced = symmetry_reduce(tbs, self.get_structure("Sn"))
        tkey = Tensor.from_values_indices([0.01], [(0, 0)])
        tval = reduced[tkey]
        for tens_1, tens_2 in zip(tval, reduced[tbs[0]]):
            assert approx(tens_1) == tens_2
        # Test set
        reduced[tkey] = "test_val"
        assert reduced[tkey] == "test_val"
        # Test empty initialization
        empty = TensorMapping()
        assert empty._tensor_list == []

        # test adding to empty tensor mapping
        empty[tkey] = 1
        assert empty[tkey] == 1

    def test_populate(self):
        test_data = loadfn(f"{TEST_FILES_DIR}/test_toec_data.json")

        sn = self.get_structure("Sn")
        vtens = np.zeros((6, 6))
        vtens[0, 0] = 259.31
        vtens[0, 1] = 160.71
        vtens[3, 3] = 73.48
        et = Tensor.from_voigt(vtens)
        populated = et.populate(sn, prec=1e-3).voigt.round(2)
        assert populated[1, 1] == approx(259.31)
        assert populated[2, 2] == approx(259.31)
        assert populated[0, 2] == approx(160.71)
        assert populated[1, 2] == approx(160.71)
        assert populated[4, 4] == approx(73.48)
        assert populated[5, 5] == approx(73.48)
        # test a rank 6 example
        vtens = np.zeros([6] * 3)
        indices = [(0, 0, 0), (0, 0, 1), (0, 1, 2), (0, 3, 3), (0, 5, 5), (3, 4, 5)]
        values = [-1271.0, -814.0, -50.0, -3.0, -780.0, -95.0]
        for v, idx in zip(values, indices):
            vtens[idx] = v
        toec = Tensor.from_voigt(vtens)
        toec = toec.populate(sn, prec=1e-3, verbose=True)
        assert toec.voigt[1, 1, 1] == approx(-1271)
        assert toec.voigt[0, 1, 1] == approx(-814)
        assert toec.voigt[0, 2, 2] == approx(-814)
        assert toec.voigt[1, 4, 4] == approx(-3)
        assert toec.voigt[2, 5, 5] == approx(-3)
        assert toec.voigt[1, 2, 0] == approx(-50)
        assert toec.voigt[4, 5, 3] == approx(-95)

        et = Tensor.from_voigt(test_data["C3_raw"]).fit_to_structure(sn)
        new = np.zeros(et.voigt.shape)
        for idx in indices:
            new[idx] = et.voigt[idx]
        new = Tensor.from_voigt(new).populate(sn)
        assert np.allclose(new, et, atol=1e-2)

    def test_from_values_indices(self):
        sn = self.get_structure("Sn")
        indices = [(0, 0), (0, 1), (3, 3)]
        values = [259.31, 160.71, 73.48]
        et = Tensor.from_values_indices(values, indices, structure=sn, populate=True).voigt.round(4)
        assert et[1, 1] == approx(259.31)
        assert et[2, 2] == approx(259.31)
        assert et[0, 2] == approx(160.71)
        assert et[1, 2] == approx(160.71)
        assert et[4, 4] == approx(73.48)
        assert et[5, 5] == approx(73.48)

    def test_serialization(self):
        # Test base serialize-deserialize
        d = self.symm_rank2.as_dict()
        new = Tensor.from_dict(d)
        assert np.allclose(new, self.symm_rank2)

        d = self.symm_rank3.as_dict(voigt=True)
        new = Tensor.from_dict(d)
        assert np.allclose(new, self.symm_rank3)

    def test_projection_methods(self):
        assert self.rand_rank2.project([1, 0, 0]) == approx(self.rand_rank2[0, 0])
        assert self.rand_rank2.project([1, 1, 1]) == approx(np.sum(self.rand_rank2) / 3)
        # Test integration
        assert np.allclose(self.ones.average_over_unit_sphere(), 1)

    def test_summary_methods(self):
        assert set(self.ones.get_grouped_indices()[0]) == set(itertools.product(range(3), range(3)))
        assert self.ones.get_grouped_indices(voigt=True)[0] == [(i,) for i in range(6)]
        assert self.ones.get_symbol_dict() == {"T_1": 1}
        assert self.ones.get_symbol_dict(voigt=False) == {"T_11": 1}

    def test_round(self):
        test = self.non_symm + 0.01
        rounded = test.round(1)
        assert np.allclose(rounded, self.non_symm)
        assert isinstance(rounded, Tensor)


class TestTensorCollection(PymatgenTest):
    def setUp(self):
        self.seq_tc = list(np.arange(4 * 3**3).reshape((4, 3, 3, 3)))
        self.seq_tc = TensorCollection(self.seq_tc)
        self.rand_tc = TensorCollection(list(np.random.random((4, 3, 3))))
        self.diff_rank = TensorCollection([np.ones([3] * i) for i in range(2, 5)])
        self.struct = self.get_structure("Si")
        ieee_file_path = f"{TEST_FILES_DIR}/ieee_conversion_data.json"
        self.ieee_data = loadfn(ieee_file_path)

    def list_based_function_check(self, attribute, coll, *args, **kwargs):
        """
        This function allows for more efficient testing of list-based
        functions in a "collection"-style class like TensorCollection.

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
                assert np.allclose(this_mod, t_mod)

    def test_list_based_functions(self):
        # zeroed
        tc = TensorCollection([1e-4 * Tensor(np.eye(3))] * 4)
        for t in tc.zeroed():
            assert t == approx(0)
        for t in tc.zeroed(1e-5):
            assert t == approx(1e-4 * np.eye(3))
        self.list_based_function_check("zeroed", tc)
        self.list_based_function_check("zeroed", tc, tol=1e-5)

        # transform
        symm_op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False, [0, 0, 1])
        self.list_based_function_check("transform", self.seq_tc, symm_op=symm_op)

        # symmetrized
        self.list_based_function_check("symmetrized", self.seq_tc)

        # rotation
        a = 3.14 * 42.5 / 180
        rotation = SquareTensor([[math.cos(a), 0, math.sin(a)], [0, 1, 0], [-math.sin(a), 0, math.cos(a)]])
        self.list_based_function_check("rotate", self.diff_rank, matrix=rotation)

        # is_symmetric
        assert not self.seq_tc.is_symmetric()
        assert self.diff_rank.is_symmetric()

        # fit_to_structure
        self.list_based_function_check("fit_to_structure", self.diff_rank, self.struct)
        self.list_based_function_check("fit_to_structure", self.seq_tc, self.struct)

        # fit_to_structure
        self.list_based_function_check("fit_to_structure", self.diff_rank, self.struct)
        self.list_based_function_check("fit_to_structure", self.seq_tc, self.struct)

        # voigt
        self.list_based_function_check("voigt", self.diff_rank)

        # is_voigt_symmetric
        assert self.diff_rank.is_voigt_symmetric()
        assert not self.seq_tc.is_voigt_symmetric()

        # Convert to ieee
        for entry in self.ieee_data[:2]:
            entry["xtal"]
            tc = TensorCollection([entry["original_tensor"]] * 3)
            struct = entry["structure"]
            self.list_based_function_check("convert_to_ieee", tc, struct)

        # from_voigt
        tc_input = list(np.random.random((3, 6, 6)))
        tc = TensorCollection.from_voigt(tc_input)
        for t_input, t in zip(tc_input, tc):
            assert np.allclose(Tensor.from_voigt(t_input), t)

    def test_serialization(self):
        # Test base serialize-deserialize
        dct = self.seq_tc.as_dict()
        new = TensorCollection.from_dict(dct)
        for t, t_new in zip(self.seq_tc, new):
            assert np.allclose(t, t_new)

        voigt_symmetrized = self.rand_tc.voigt_symmetrized
        dct = voigt_symmetrized.as_dict(voigt=True)
        new_vsym = TensorCollection.from_dict(dct)
        for t, t_new in zip(voigt_symmetrized, new_vsym):
            assert np.allclose(t, t_new)


class TestSquareTensor(PymatgenTest):
    def setUp(self):
        self.rand_sqtensor = SquareTensor(np.random.randn(3, 3))
        self.symm_sqtensor = SquareTensor([[0.1, 0.3, 0.4], [0.3, 0.5, 0.2], [0.4, 0.2, 0.6]])
        self.non_invertible = SquareTensor([[0.1, 0, 0], [0.2, 0, 0], [0, 0, 0]])
        self.non_symm = SquareTensor([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.2, 0.5, 0.5]])
        self.low_val = SquareTensor([[1e-6, 1 + 1e-5, 1e-6], [1 + 1e-6, 1e-6, 1e-6], [1e-7, 1e-7, 1 + 1e-5]])
        self.low_val_2 = SquareTensor([[1e-6, -1 - 1e-6, 1e-6], [1 + 1e-7, 1e-6, 1e-6], [1e-7, 1e-7, 1 + 1e-6]])
        a = 3.14 * 42.5 / 180
        self.rotation = SquareTensor([[math.cos(a), 0, math.sin(a)], [0, 1, 0], [-math.sin(a), 0, math.cos(a)]])

    def test_new(self):
        non_sq_matrix = [
            [0.1, 0.2, 0.1],
            [0.1, 0.2, 0.3],
            [0.1, 0.2, 0.3],
            [0.1, 0.1, 0.1],
        ]
        bad_matrix = [[0.1, 0.2], [0.2, 0.3, 0.4], [0.2, 0.3, 0.5]]
        too_high_rank = np.zeros((3, 3, 3))
        with pytest.raises(
            ValueError,
            match="Pymatgen only supports 3-dimensional tensors, and default tensor constructor uses standard notation",
        ):
            SquareTensor(non_sq_matrix)
        with pytest.raises(ValueError, match="The requested array has an inhomogeneous shape after 1 dimensions."):
            SquareTensor(bad_matrix)
        with pytest.raises(ValueError, match="SquareTensor input must be rank 2"):
            SquareTensor(too_high_rank)

    def test_properties(self):
        # transpose
        assert self.non_symm.trans == approx(SquareTensor([[0.1, 0.4, 0.2], [0.2, 0.5, 0.5], [0.3, 0.6, 0.5]]))
        assert self.rand_sqtensor.trans == approx(np.transpose(self.rand_sqtensor))
        assert self.symm_sqtensor == approx(self.symm_sqtensor.trans)
        # inverse
        assert self.non_symm.inv == approx(np.linalg.inv(self.non_symm))
        with pytest.raises(ValueError, match="SquareTensor is non-invertible"):
            _ = self.non_invertible.inv

        # determinant
        assert self.rand_sqtensor.det == np.linalg.det(self.rand_sqtensor)
        assert self.non_invertible.det == 0.0
        assert self.non_symm.det == 0.009

        # symmetrized
        assert self.rand_sqtensor.symmetrized == approx(0.5 * (self.rand_sqtensor + self.rand_sqtensor.trans))
        assert self.symm_sqtensor == approx(self.symm_sqtensor.symmetrized)
        assert np.allclose(
            self.non_symm.symmetrized,
            SquareTensor([[0.1, 0.3, 0.25], [0.3, 0.5, 0.55], [0.25, 0.55, 0.5]]),
        )

        # invariants
        i1 = np.trace(self.rand_sqtensor)
        i2 = (
            self.rand_sqtensor[0, 0] * self.rand_sqtensor[1, 1]
            + self.rand_sqtensor[1, 1] * self.rand_sqtensor[2, 2]
            + self.rand_sqtensor[2, 2] * self.rand_sqtensor[0, 0]
            - self.rand_sqtensor[0, 1] * self.rand_sqtensor[1, 0]
            - self.rand_sqtensor[0, 2] * self.rand_sqtensor[2, 0]
            - self.rand_sqtensor[2, 1] * self.rand_sqtensor[1, 2]
        )
        i3 = np.linalg.det(self.rand_sqtensor)
        assert np.allclose([i1, i2, i3], self.rand_sqtensor.principal_invariants)

    def test_is_rotation(self):
        assert self.rotation.is_rotation()
        assert not self.symm_sqtensor.is_rotation()
        assert self.low_val_2.is_rotation()
        assert not self.low_val_2.is_rotation(tol=1e-8)

    def test_refine_rotation(self):
        assert np.allclose(self.rotation, self.rotation.refine_rotation())
        new = self.rotation.copy()
        new[2, 2] += 0.02
        assert not new.is_rotation()
        assert np.allclose(self.rotation, new.refine_rotation())
        new[1] *= 1.05
        assert np.allclose(self.rotation, new.refine_rotation())

    def test_get_scaled(self):
        assert self.non_symm.get_scaled(10.0) == approx(SquareTensor([[1, 2, 3], [4, 5, 6], [2, 5, 5]]))

    def test_polar_decomposition(self):
        u, p = self.rand_sqtensor.polar_decomposition()
        assert np.allclose(np.dot(u, p), self.rand_sqtensor)
        assert np.allclose(np.eye(3), np.dot(u, np.conjugate(np.transpose(u))))

    def test_serialization(self):
        # Test base serialize-deserialize
        dct = self.rand_sqtensor.as_dict()
        new = SquareTensor.from_dict(dct)
        assert np.allclose(new, self.rand_sqtensor)
        assert isinstance(new, SquareTensor)

        # Ensure proper object-independent deserialization
        obj = MontyDecoder().process_decoded(dct)
        assert isinstance(obj, SquareTensor)

        vsym = self.rand_sqtensor.voigt_symmetrized
        d_vsym = vsym.as_dict(voigt=True)
        new_voigt = Tensor.from_dict(d_vsym)
        assert np.allclose(vsym, new_voigt)
