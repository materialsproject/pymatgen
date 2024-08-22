from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose

from pymatgen.core.operations import MagSymmOp, SymmOp
from pymatgen.electronic_structure.core import Magmom
from pymatgen.util.testing import PymatgenTest


class TestSymmOp(PymatgenTest):
    def setUp(self):
        self.op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, translation_vec=[0, 0, 1])

    def test_properties(self):
        rot = self.op.rotation_matrix
        vec = self.op.translation_vector
        assert_allclose(rot, [[0.8660254, -0.5, 0.0], [0.5, 0.8660254, 0.0], [0.0, 0.0, 1.0]], 2)
        assert_allclose(vec, [0, 0, 1], 2)

    def test_operate(self):
        point = np.array([1, 2, 3])
        new_coord = self.op.operate(point)
        assert_allclose(new_coord, [-0.1339746, 2.23205081, 4.0], 2)

    def test_operate_multi(self):
        point = np.array([1, 2, 3])
        new_coords = self.op.operate_multi([point, point])
        assert_allclose(new_coords, [[-0.1339746, 2.23205081, 4.0]] * 2, 2)
        new_coords = self.op.operate_multi([[point, point]] * 2)
        assert_allclose(new_coords, [[[-0.1339746, 2.23205081, 4.0]] * 2] * 2, 2)

    def test_inverse(self):
        point = np.random.default_rng().random(3)
        new_coord = self.op.operate(point)
        assert_allclose(self.op.inverse.operate(new_coord), point, 2)

    def test_reflection(self):
        rng = np.random.default_rng()
        normal = rng.random(3)
        origin = rng.random(3)
        refl = SymmOp.reflection(normal, origin)
        point = rng.random(3)
        new_coord = refl.operate(point)
        # Distance to the plane should be negatives of each other.
        assert_allclose(np.dot(new_coord - origin, normal), -np.dot(point - origin, normal))

    def test_apply_rotation_only(self):
        point = np.random.default_rng().random(3)
        new_coord = self.op.operate(point)
        rotate_only = self.op.apply_rotation_only(point)
        assert_allclose(rotate_only + self.op.translation_vector, new_coord, 2)

    def test_transform_tensor(self):
        # Rank 2
        tensor = np.arange(0, 9).reshape(3, 3)
        new_tensor = self.op.transform_tensor(tensor)
        assert_allclose(
            new_tensor,
            [
                [-0.73205, -1.73205, -0.76794],
                [0.26795, 4.73205, 5.33013],
                [1.69615, 9.06218, 8.0],
            ],
            5,
        )

        # Rank 3
        tensor = np.arange(0, 27).reshape(3, 3, 3)
        new_tensor = self.op.transform_tensor(tensor)
        assert_allclose(
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
        # Rank 4
        tensor = np.arange(0, 81).reshape(3, 3, 3, 3)
        new_tensor = self.op.transform_tensor(tensor)
        assert_allclose(
            new_tensor,
            [
                [
                    [
                        [-0.981, -3.526, -2.514],
                        [-3.258, -11.660, -8.286],
                        [-2.184, -7.786, -5.517],
                    ],
                    [
                        [-2.454, -8.660, -6.090],
                        [-7.660, -26.722, -18.629],
                        [-4.858, -16.763, -11.588],
                    ],
                    [
                        [-1.194, -4.090, -2.811],
                        [-3.358, -11.165, -7.490],
                        [-1.909, -6.124, -3.983],
                    ],
                ],
                [
                    [
                        [-0.043, 0.340, 0.499],
                        [1.340, 6.866, 5.959],
                        [1.731, 7.825, 6.412],
                    ],
                    [
                        [4.340, 18.062, 14.155],
                        [21.794, 88.301, 68.123],
                        [18.754, 75.087, 57.517],
                    ],
                    [
                        [5.427, 21.620, 16.510],
                        [24.352, 95.979, 72.811],
                        [19.876, 77.909, 58.899],
                    ],
                ],
                [
                    [
                        [1.777, 6.999, 5.306],
                        [7.731, 30.218, 22.804],
                        [6.208, 24.170, 18.194],
                    ],
                    [
                        [9.927, 38.414, 28.804],
                        [41.146, 158.656, 118.694],
                        [32.170, 123.792, 92.488],
                    ],
                    [
                        [8.914, 34.268, 25.586],
                        [36.268, 139.086, 103.684],
                        [28.050, 107.416, 80.000],
                    ],
                ],
            ],
            3,
        )

    def test_are_symmetrically_related(self):
        point = np.random.default_rng().random(3)
        new_coord = self.op.operate(point)
        assert self.op.are_symmetrically_related(point, new_coord)
        assert self.op.are_symmetrically_related(new_coord, point)

    def test_are_symmetrically_related_vectors(self):
        tol = 0.001
        rng = np.random.default_rng()
        from_a = rng.random(3)
        to_a = rng.random(3)
        r_a = rng.integers(0, 10, 3)
        from_b = self.op.operate(from_a)
        to_b = self.op.operate(to_a)
        floored = np.floor([from_b, to_b])
        is_too_close = np.abs([from_b, to_b] - floored) > 1 - tol
        floored[is_too_close] += 1
        r_b = self.op.apply_rotation_only(r_a) - floored[0] + floored[1]
        from_b = from_b % 1
        to_b = to_b % 1
        assert self.op.are_symmetrically_related_vectors(from_a, to_a, r_a, from_b, to_b, r_b)[0]
        assert not self.op.are_symmetrically_related_vectors(from_a, to_a, r_a, from_b, to_b, r_b)[1]
        assert self.op.are_symmetrically_related_vectors(to_a, from_a, -r_a, from_b, to_b, r_b)[0]
        assert self.op.are_symmetrically_related_vectors(to_a, from_a, -r_a, from_b, to_b, r_b)[1]

    def test_as_from_dict(self):
        dct = self.op.as_dict()
        op = SymmOp.from_dict(dct)
        point = np.random.default_rng().random(3)
        new_coord = self.op.operate(point)
        assert op.are_symmetrically_related(point, new_coord)

    def test_inversion(self):
        rng = np.random.default_rng()
        origin = rng.random(3)
        op = SymmOp.inversion(origin)
        pt = rng.random(3)
        inv_pt = op.operate(pt)
        assert_allclose(pt - origin, origin - inv_pt)

    def test_xyz(self):
        op = SymmOp([[1, -1, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]])
        xyz_str = op.as_xyz_str()
        assert xyz_str == "x-y, -y, -z"
        assert op == SymmOp.from_xyz_str(xyz_str)

        op2 = SymmOp([[0, -1, 0, 0.5], [1, 0, 0, 0.5], [0, 0, 1, 0.5 + 1e-7], [0, 0, 0, 1]])
        s2 = op2.as_xyz_str()
        assert s2 == "-y+1/2, x+1/2, z+1/2"
        assert op2 == SymmOp.from_xyz_str(s2)

        op2 = SymmOp([[3, -2, -1, 0.5], [-1, 0, 0, 12.0 / 13], [0, 0, 1, 0.5 + 1e-7], [0, 0, 0, 1]])
        s2 = op2.as_xyz_str()
        assert s2 == "3x-2y-z+1/2, -x+12/13, z+1/2"
        assert op2 == SymmOp.from_xyz_str(s2)

        op3 = SymmOp.from_xyz_str("3x - 2y - z+1 /2 , -x+12/ 13, z+1/2")
        assert op2 == op3

        # Ensure strings can be read in any order
        op4 = SymmOp.from_xyz_str("1 /2 + 3X - 2y - z , 12/ 13-x, z+1/2")
        op5 = SymmOp.from_xyz_str("+1 /2 + 3x - 2y - z , 12/ 13-x, +1/2+z")
        assert op4 == op3
        assert op4 == op5
        assert op3 == op5

        self.assertWarns(UserWarning, self.op.as_xyz_str)

        symm_op = SymmOp.from_xyz_str("0.5+x, 0.25+y, 0.75+z")
        assert_allclose(symm_op.translation_vector, [0.5, 0.25, 0.75])
        symm_op = SymmOp.from_xyz_str("x + 0.5, y + 0.25, z + 0.75")
        assert_allclose(symm_op.translation_vector, [0.5, 0.25, 0.75])


class TestMagSymmOp(PymatgenTest):
    def test_xyzt_string(self):
        xyzt_strings = ["x, y, z, +1", "x, y, z, -1", "-y+1/2, x+1/2, x+1/2, +1"]

        for xyzt_string in xyzt_strings:
            op = MagSymmOp.from_xyzt_str(xyzt_string)
            xyzt_string_out = op.as_xyzt_str()
            assert xyzt_string == xyzt_string_out

        op = SymmOp([[3, -2, -1, 0.5], [-1, 0, 0, 12.0 / 13], [0, 0, 1, 0.5 + 1e-7], [0, 0, 0, 1]])

        magop = MagSymmOp.from_symmop(op, -1)
        magop_str = magop.as_xyzt_str()
        assert magop.time_reversal == -1
        assert magop_str == "3x-2y-z+1/2, -x+12/13, z+1/2, -1"

    def test_as_from_dict(self):
        op = SymmOp([[3, -2, -1, 0.5], [-1, 0, 0, 12.0 / 13], [0, 0, 1, 0.5 + 1e-7], [0, 0, 0, 1]])
        magop = MagSymmOp.from_symmop(op, -1)
        magop2 = MagSymmOp.from_dict(magop.as_dict())
        assert magop2.time_reversal == -1
        assert magop2.as_xyzt_str() == "3x-2y-z+1/2, -x+12/13, z+1/2, -1"

    def test_operate_magmom(self):
        # all test magmoms are the same
        magmoms = [
            Magmom([1, 2, 3]),  # as Magmom
            [1, 2, 3],  # as list
            Magmom([-3, 2, 1], saxis=[1, 0, 0]),
        ]  # as Magmom with non-default saxis

        xyzt_strings = ["x, y, z, +1", "x, y, z, -1", "x, -y, z, -1", "-x, -y, z, -1"]

        transformed_magmoms = [[1, 2, 3], [-1, -2, -3], [1, -2, 3], [1, 2, -3]]

        for xyzt_string, transformed_magmom in zip(xyzt_strings, transformed_magmoms, strict=True):
            for magmom in magmoms:
                op = MagSymmOp.from_xyzt_str(xyzt_string)
                assert_allclose(transformed_magmom, op.operate_magmom(magmom).global_moment)
