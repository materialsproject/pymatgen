# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import itertools

import numpy as np
import pytest

from pymatgen.core.lattice import Lattice, get_points_in_spheres
from pymatgen.core.operations import SymmOp
from pymatgen.util.testing import PymatgenTest


class LatticeTestCase(PymatgenTest):
    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.cubic = self.lattice
        self.tetragonal = Lattice.tetragonal(10, 20)
        self.orthorhombic = Lattice.orthorhombic(10, 20, 30)
        self.monoclinic = Lattice.monoclinic(10, 20, 30, 66)
        self.hexagonal = Lattice.hexagonal(10, 20)
        self.rhombohedral = Lattice.rhombohedral(10, 77)

        self.cubic_partial_pbc = Lattice.cubic(10.0, pbc=(True, True, False))

        family_names = [
            "cubic",
            "tetragonal",
            "orthorhombic",
            "monoclinic",
            "hexagonal",
            "rhombohedral",
        ]

        self.families = {}
        for name in family_names:
            self.families[name] = getattr(self, name)

    def test_equal(self):
        assert self.cubic == self.cubic
        assert self.cubic == self.lattice
        for name1, latt1 in self.families.items():
            for name2, latt2 in self.families.items():
                assert (name1 == name2) == (latt1 == latt2)

        # ensure partial periodic boundaries is unequal to all full periodic boundaries
        assert not any(self.cubic_partial_pbc == x for x in self.families.values())

    def test_format(self):
        assert "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]" == format(
            self.lattice, ".3fl"
        )
        assert """10.000 0.000 0.000
0.000 10.000 0.000
0.000 0.000 10.000""" == format(
            self.lattice, ".3f"
        )
        assert "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}" == format(self.lattice, ".1fp")

    def test_init(self):
        a = 9.026
        lattice = Lattice.cubic(a)
        assert lattice is not None, "Initialization from new_cubic failed"
        self.assertArrayEqual(lattice.pbc, (True, True, True))
        lattice2 = Lattice([[a, 0, 0], [0, a, 0], [0, 0, a]])
        for i in range(0, 3):
            for j in range(0, 3):
                assert (
                    round(abs(lattice.matrix[i][j] - lattice2.matrix[i][j]), 5) == 0
                ), "Inconsistent matrix from two inits!"
        self.assertArrayEqual(self.cubic_partial_pbc.pbc, (True, True, False))

    def test_copy(self):
        cubic_copy = self.cubic.copy()
        assert cubic_copy == self.cubic
        assert self.cubic_partial_pbc.copy() == self.cubic_partial_pbc
        assert cubic_copy._matrix is not self.cubic._matrix

    def test_get_cartesian_or_frac_coord(self):
        coord = self.lattice.get_cartesian_coords([0.15, 0.3, 0.4])
        self.assertArrayAlmostEqual(coord, [1.5, 3.0, 4.0])
        self.assertArrayAlmostEqual(
            self.tetragonal.get_fractional_coords([12.12312, 45.2134, 1.3434]),
            [1.212312, 4.52134, 0.06717],
        )

        # Random testing that get_cart and get_frac coords reverses each other.
        rand_coord = np.random.random_sample(3)
        coord = self.tetragonal.get_cartesian_coords(rand_coord)
        fcoord = self.tetragonal.get_fractional_coords(coord)
        self.assertArrayAlmostEqual(fcoord, rand_coord)

    def test_get_vector_along_lattice_directions(self):
        lattice_mat = np.array([[0.5, 0.0, 0.0], [0.5, np.sqrt(3) / 2.0, 0.0], [0.0, 0.0, 1.0]])
        lattice = Lattice(lattice_mat)
        cart_coord = np.array([0.5, np.sqrt(3) / 4.0, 0.5])
        latt_coord = np.array([0.25, 0.5, 0.5])
        from_direct = lattice.get_fractional_coords(cart_coord) * lattice.lengths
        self.assertArrayAlmostEqual(lattice.get_vector_along_lattice_directions(cart_coord), from_direct)
        self.assertArrayAlmostEqual(lattice.get_vector_along_lattice_directions(cart_coord), latt_coord)
        self.assertArrayEqual(
            lattice.get_vector_along_lattice_directions(cart_coord).shape,
            [
                3,
            ],
        )
        self.assertArrayEqual(
            lattice.get_vector_along_lattice_directions(cart_coord.reshape([1, 3])).shape,
            [1, 3],
        )

    def test_d_hkl(self):
        cubic_copy = self.cubic.copy()
        hkl = (1, 2, 3)
        dhkl = ((hkl[0] ** 2 + hkl[1] ** 2 + hkl[2] ** 2) / (cubic_copy.a**2)) ** (-1 / 2)
        assert dhkl == cubic_copy.d_hkl(hkl)

    def test_reciprocal_lattice(self):
        recip_latt = self.lattice.reciprocal_lattice
        self.assertArrayAlmostEqual(recip_latt.matrix, 0.628319 * np.eye(3), 5)
        self.assertArrayAlmostEqual(
            self.tetragonal.reciprocal_lattice.matrix,
            [[0.628319, 0.0, 0.0], [0.0, 0.628319, 0], [0.0, 0.0, 0.3141590]],
            5,
        )

        # Test the crystallographic version.
        recip_latt_xtal = self.lattice.reciprocal_lattice_crystallographic
        self.assertArrayAlmostEqual(recip_latt.matrix, recip_latt_xtal.matrix * 2 * np.pi, 5)

    def test_static_methods(self):
        lengths_c = [3.840198, 3.84019885, 3.8401976]
        angles_c = [119.99998575, 90, 60.00000728]
        mat_c = [
            [3.840198, 0.000000, 0.0000],
            [1.920099, 3.325710, 0.000000],
            [0.000000, -2.217138, 3.135509],
        ]
        # should give the lengths and angles above
        newlatt = Lattice(mat_c)
        lengths = newlatt.lengths
        angles = newlatt.angles
        for i in range(0, 3):
            assert round(abs(lengths[i] - lengths_c[i]), 5) == 0, "Lengths incorrect!"
            assert round(abs(angles[i] - angles_c[i]), 5) == 0, "Angles incorrect!"
        latt = Lattice.from_parameters(*lengths, *angles)
        lengths = latt.lengths
        angles = latt.angles
        for i in range(0, 3):
            assert round(abs(lengths[i] - lengths_c[i]), 5) == 0, "Lengths incorrect!"
            assert round(abs(angles[i] - angles_c[i]), 5) == 0, "Angles incorrect!"

    def test_attributes(self):
        """docstring for test_attributes"""
        lattice = Lattice.cubic(10.0)
        assert lattice.a == 10.0
        assert lattice.b == 10.0
        assert lattice.c == 10.0
        assert round(abs(lattice.volume - 1000.0), 7) == 0
        xyz = lattice.get_cartesian_coords([0.25, 0.35, 0.45])
        assert xyz[0] == 2.5
        assert xyz[1] == 3.5
        assert xyz[2] == 4.5

    def test_lattice_matrices(self):
        """
        If alpha == 90 and beta == 90, two matricies are identical.
        """

        def _identical(a, b, c, alpha, beta, gamma):
            mat1 = Lattice.from_parameters(a, b, c, alpha, beta, gamma, False).matrix
            mat2 = Lattice.from_parameters(a, b, c, alpha, beta, gamma, True).matrix
            # self.assertArrayAlmostEqual(mat1, mat2)
            return ((mat1 - mat2) ** 2).sum() < 1e-6

        assert _identical(2, 3, 4, 90, 90, 90)
        assert _identical(2, 3, 4, 90, 90, 80)
        assert _identical(2, 3, 4, 90, 90, 100)

        assert not _identical(2, 3, 4, 100, 90, 90)
        assert not _identical(2, 3, 4, 90, 100, 90)
        assert not _identical(2, 3, 4, 100, 100, 100)

    def test_get_lll_reduced_lattice(self):
        lattice = Lattice([1.0, 1, 1, -1.0, 0, 2, 3.0, 5, 6])
        reduced_latt = lattice.get_lll_reduced_lattice()

        expected_ans = Lattice(np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0, -2.0, 0.0, 1.0]).reshape((3, 3)))
        assert round(abs(np.linalg.det(np.linalg.solve(expected_ans.matrix, reduced_latt.matrix)) - 1), 7) == 0
        self.assertArrayAlmostEqual(sorted(reduced_latt.abc), sorted(expected_ans.abc))
        assert round(abs(reduced_latt.volume - lattice.volume), 7) == 0
        latt = [
            7.164750,
            2.481942,
            0.000000,
            -4.298850,
            2.481942,
            0.000000,
            0.000000,
            0.000000,
            14.253000,
        ]
        expected_ans = Lattice(
            np.array(
                [
                    -4.298850,
                    2.481942,
                    0.000000,
                    2.865900,
                    4.963884,
                    0.000000,
                    0.000000,
                    0.000000,
                    14.253000,
                ]
            )
        )
        reduced_latt = Lattice(latt).get_lll_reduced_lattice()
        assert round(abs(np.linalg.det(np.linalg.solve(expected_ans.matrix, reduced_latt.matrix)) - 1), 7) == 0
        self.assertArrayAlmostEqual(sorted(reduced_latt.abc), sorted(expected_ans.abc))

        expected_ans = Lattice([0.0, 10.0, 10.0, 10.0, 10.0, 0.0, 30.0, -30.0, 40.0])

        lattice = np.array([100.0, 0.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0, 10.0])
        lattice = lattice.reshape(3, 3)
        lattice = Lattice(lattice.T)
        reduced_latt = lattice.get_lll_reduced_lattice()
        assert round(abs(np.linalg.det(np.linalg.solve(expected_ans.matrix, reduced_latt.matrix)) - 1), 7) == 0
        self.assertArrayAlmostEqual(sorted(reduced_latt.abc), sorted(expected_ans.abc))

        random_latt = Lattice(np.random.random((3, 3)))
        if np.linalg.det(random_latt.matrix) > 1e-8:
            reduced_random_latt = random_latt.get_lll_reduced_lattice()
            assert round(abs(reduced_random_latt.volume - random_latt.volume), 7) == 0

    def test_get_niggli_reduced_lattice(self):
        latt = Lattice.from_parameters(3, 5.196, 2, 103 + 55 / 60, 109 + 28 / 60, 134 + 53 / 60)
        reduced_cell = latt.get_niggli_reduced_lattice()
        abc = reduced_cell.lengths
        angles = reduced_cell.angles
        assert round(abs(abc[0] - 2), 3) == 0
        assert round(abs(abc[1] - 3), 3) == 0
        assert round(abs(abc[2] - 3), 3) == 0
        assert round(abs(angles[0] - 116.382855225), 3) == 0
        assert round(abs(angles[1] - 94.769790287999996), 3) == 0
        assert round(abs(angles[2] - 109.466666667), 3) == 0

        mat = [[5.0, 0, 0], [0, 5.0, 0], [5.0, 0, 5.0]]
        latt = Lattice(np.dot([[1, 1, 1], [1, 1, 0], [0, 1, 1]], mat))
        reduced_cell = latt.get_niggli_reduced_lattice()
        abc = reduced_cell.lengths
        angles = reduced_cell.angles
        for l in abc:
            assert round(abs(l - 5), 3) == 0
        for a in angles:
            assert round(abs(a - 90), 3) == 0

        latt = Lattice(
            [
                1.432950,
                0.827314,
                4.751000,
                -1.432950,
                0.827314,
                4.751000,
                0.0,
                -1.654628,
                4.751000,
            ]
        )
        # ans = [[-1.432950, -2.481942, 0.0],
        #       [-2.8659, 0.0, 0.0],
        #       [-1.432950, -0.827314, -4.751000]]
        ans = [
            [-1.43295, -2.481942, 0.0],
            [-2.8659, 0.0, 0.0],
            [-1.43295, -0.827314, -4.751],
        ]
        self.assertArrayAlmostEqual(latt.get_niggli_reduced_lattice().matrix, ans)

        latt = Lattice.from_parameters(7.365450, 6.199506, 5.353878, 75.542191, 81.181757, 156.396627)
        ans = [
            [2.578932, 0.826965, 0.000000],
            [-0.831059, 2.067413, 1.547813],
            [-0.458407, -2.480895, 1.129126],
        ]
        self.assertArrayAlmostEqual(latt.get_niggli_reduced_lattice().matrix, np.array(ans), 5)

    def test_find_mapping(self):
        m = np.array([[0.1, 0.2, 0.3], [-0.1, 0.2, 0.7], [0.6, 0.9, 0.2]])
        latt = Lattice(m)

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [2, 3, 3], 35)
        rot = op.rotation_matrix
        scale = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])

        latt2 = Lattice(np.dot(rot, np.dot(scale, m).T).T)
        (aligned_out, rot_out, scale_out) = latt2.find_mapping(latt)
        assert round(abs(abs(np.linalg.det(rot)) - 1), 7) == 0

        rotated = SymmOp.from_rotation_and_translation(rot_out).operate_multi(latt.matrix)

        self.assertArrayAlmostEqual(rotated, aligned_out.matrix)
        self.assertArrayAlmostEqual(np.dot(scale_out, latt2.matrix), aligned_out.matrix)
        self.assertArrayAlmostEqual(aligned_out.parameters, latt.parameters)
        assert not np.allclose(aligned_out.parameters, latt2.parameters)

    def test_find_all_mappings(self):
        m = np.array([[0.1, 0.2, 0.3], [-0.1, 0.2, 0.7], [0.6, 0.9, 0.2]])
        latt = Lattice(m)

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [2, -1, 3], 40)
        rot = op.rotation_matrix
        scale = np.array([[0, 2, 0], [1, 1, 0], [0, 0, 1]])

        latt2 = Lattice(np.dot(rot, np.dot(scale, m).T).T)

        for (aligned_out, rot_out, scale_out) in latt.find_all_mappings(latt2):
            self.assertArrayAlmostEqual(np.inner(latt2.matrix, rot_out), aligned_out.matrix, 5)
            self.assertArrayAlmostEqual(np.dot(scale_out, latt.matrix), aligned_out.matrix)
            self.assertArrayAlmostEqual(aligned_out.parameters, latt2.parameters)
            assert not np.allclose(aligned_out.parameters, latt.parameters)

        latt = Lattice.orthorhombic(9, 9, 5)
        assert len(list(latt.find_all_mappings(latt))) == 16

        # catch the singular matrix error
        latt = Lattice.from_parameters(1, 1, 1, 10, 10, 10)
        for l, _, _ in latt.find_all_mappings(latt, ltol=0.05, atol=11):
            assert isinstance(l, Lattice)

    def test_mapping_symmetry(self):
        l = Lattice.cubic(1)
        l2 = Lattice.orthorhombic(1.1001, 1, 1)
        assert l.find_mapping(l2, ltol=0.1) is None
        assert l2.find_mapping(l, ltol=0.1) is None
        l2 = Lattice.orthorhombic(1.0999, 1, 1)
        assert l2.find_mapping(l, ltol=0.1) is not None
        assert l.find_mapping(l2, ltol=0.1) is not None

    def test_to_from_dict(self):
        d = self.tetragonal.as_dict()
        t = Lattice.from_dict(d)
        for i in range(3):
            assert t.abc[i] == self.tetragonal.abc[i]
            assert t.angles[i] == self.tetragonal.angles[i]
        # Make sure old style dicts work.
        d = self.tetragonal.as_dict(verbosity=1)
        del d["matrix"]
        t = Lattice.from_dict(d)
        for i in range(3):
            assert t.abc[i] == self.tetragonal.abc[i]
            assert t.angles[i] == self.tetragonal.angles[i]

    def test_scale(self):
        new_volume = 10
        for lattice in self.families.values():
            new_lattice = lattice.scale(new_volume)
            assert round(abs(new_lattice.volume - new_volume), 7) == 0
            self.assertArrayAlmostEqual(new_lattice.angles, lattice.angles)

    def test_get_wigner_seitz_cell(self):
        ws_cell = Lattice([[10, 0, 0], [0, 5, 0], [0, 0, 1]]).get_wigner_seitz_cell()
        assert 6 == len(ws_cell)
        for l in ws_cell[3]:
            assert [abs(i) for i in l] == [5.0, 2.5, 0.5]

    def test_dot_and_norm(self):
        frac_basis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        for lattice in self.families.values():
            # print(family_name)
            self.assertArrayAlmostEqual(lattice.norm(lattice.matrix, frac_coords=False), lattice.abc, 5)
            self.assertArrayAlmostEqual(lattice.norm(frac_basis), lattice.abc, 5)
            for (i, vec) in enumerate(frac_basis):
                length = lattice.norm(vec)
                self.assertArrayAlmostEqual(length[0], lattice.abc[i], 5)
                # We always get a ndarray.
                assert hasattr(length, "shape")

        # Passing complex arrays should raise TypeError
        with pytest.raises(TypeError):
            lattice.norm(np.zeros(3, dtype=np.complex128))

        # Cannot reshape the second argument.
        with pytest.raises(ValueError):
            lattice.dot(np.zeros(6), np.zeros(8))

        # Passing vectors of different length is invalid.
        with pytest.raises(ValueError):
            lattice.dot(np.zeros(3), np.zeros(6))

    def test_get_points_in_sphere(self):
        # This is a non-niggli representation of a cubic lattice
        latt = Lattice([[1, 5, 0], [0, 1, 0], [5, 0, 1]])
        # evenly spaced points array between 0 and 1
        pts = np.array(list(itertools.product(range(5), repeat=3))) / 5
        pts = latt.get_fractional_coords(pts)

        # Test getting neighbors within 1 neighbor distance of the origin
        fcoords, dists, inds, images = latt.get_points_in_sphere(pts, [0, 0, 0], 0.20001, zip_results=False)
        assert len(fcoords) == 7  # There are 7 neighbors
        assert np.isclose(dists, 0.2).sum() == 6  # 6 are at 0.2
        assert np.isclose(dists, 0).sum() == 1  # 1 is at 0
        assert len(set(inds)) == 7  # They have unique indices
        self.assertArrayEqual(images[np.isclose(dists, 0)], [[0, 0, 0]])

        # More complicated case, using the zip output
        result = latt.get_points_in_sphere(pts, [0.5, 0.5, 0.5], 1.0001)
        assert len(result) == 552
        assert len(result[0]) == 4  # coords, dists, ind, supercell

        # test pbc
        latt_pbc = Lattice([[1, 5, 0], [0, 1, 0], [5, 0, 1]], pbc=(True, True, False))
        fcoords, dists, inds, images = latt_pbc.get_points_in_sphere(pts, [0, 0, 0], 0.20001, zip_results=False)
        assert len(fcoords) == 6

    def test_get_all_distances(self):
        fcoords = np.array(
            [
                [0.3, 0.3, 0.5],
                [0.1, 0.1, 0.3],
                [0.9, 0.9, 0.8],
                [0.1, 0.0, 0.5],
                [0.9, 0.7, 0.0],
            ]
        )
        lattice = Lattice.from_parameters(8, 8, 4, 90, 76, 58)
        expected = np.array(
            [
                [0.000, 3.015, 4.072, 3.519, 3.245],
                [3.015, 0.000, 3.207, 1.131, 4.453],
                [4.072, 3.207, 0.000, 2.251, 1.788],
                [3.519, 1.131, 2.251, 0.000, 3.852],
                [3.245, 4.453, 1.788, 3.852, 0.000],
            ]
        )
        output = lattice.get_all_distances(fcoords, fcoords)
        self.assertArrayAlmostEqual(output, expected, 3)
        # test just one input point
        output2 = lattice.get_all_distances(fcoords[0], fcoords)
        self.assertArrayAlmostEqual(output2, [expected[0]], 2)
        # test distance when initial points are not in unit cell
        f1 = [0, 0, 17]
        f2 = [0, 0, 10]
        assert lattice.get_all_distances(f1, f2)[0, 0] == 0

        # test pbc
        lattice_pbc = Lattice.from_parameters(8, 8, 4, 90, 76, 58, pbc=(True, True, False))
        expected_pbc = np.array(
            [
                [0.000, 3.015, 4.072, 3.519, 4.089],
                [3.015, 0.000, 3.207, 1.131, 4.453],
                [4.072, 3.207, 0.000, 2.251, 3.578],
                [3.519, 1.131, 2.251, 0.000, 4.235],
            ]
        )
        output3 = lattice_pbc.get_all_distances(fcoords[:-1], fcoords)
        self.assertArrayAlmostEqual(output3, expected_pbc, 3)

    def test_monoclinic(self):
        a, b, c, alpha, beta, gamma = self.monoclinic.parameters
        assert round(abs(beta - 90), 7) != 0
        assert round(abs(alpha - 90), 7) == 0
        assert round(abs(gamma - 90), 7) == 0

    def test_is_hexagonal(self):
        assert not self.cubic.is_hexagonal()
        assert not self.tetragonal.is_hexagonal()
        assert not self.orthorhombic.is_hexagonal()
        assert not self.monoclinic.is_hexagonal()
        assert not self.rhombohedral.is_hexagonal()
        assert self.hexagonal.is_hexagonal()

    def test_get_distance_and_image(self):
        dist, image = self.cubic.get_distance_and_image([0, 0, 0.1], [0, 0.0, 0.9])
        assert round(abs(dist - 2), 7) == 0
        self.assertArrayAlmostEqual(image, [0, 0, -1])

    def test_get_distance_and_image_strict(self):
        for _ in range(10):
            lengths = [np.random.randint(1, 100) for i in range(3)]
            lattice = [np.random.rand(3) * lengths[i] for i in range(3)]
            lattice = Lattice(np.array(lattice))

            f1 = np.random.rand(3)
            f2 = np.random.rand(3)

            scope = list(range(-3, 4))
            min_image_dist = (float("inf"), None)
            for image in itertools.product(scope, scope, scope):
                cart = lattice.get_cartesian_coords(f1 - (f2 + image))
                dist = np.dot(cart, cart) ** 0.5
                if dist < min_image_dist[0]:
                    min_image_dist = (dist, image)

            pmg_result = lattice.get_distance_and_image(f1, f2)
            assert min_image_dist[0] + 1e-7 >= pmg_result[0]
            if abs(min_image_dist[0] - pmg_result[0]) < 1e-12:
                self.assertArrayAlmostEqual(min_image_dist[1], pmg_result[1])

    def test_lll_basis(self):
        a = np.array([1.0, 0.1, 0.0])
        b = np.array([0.0, 2.0, 0.0])
        c = np.array([0.0, 0.0, 3.0])

        l1 = Lattice([a, b, c])
        l2 = Lattice([a + b, b + c, c])

        ccoords = np.array([[1, 1, 2], [2, 2, 1.5]])
        l1_fcoords = l1.get_fractional_coords(ccoords)
        l2_fcoords = l2.get_fractional_coords(ccoords)

        self.assertArrayAlmostEqual(l1.matrix, l2.lll_matrix)
        self.assertArrayAlmostEqual(np.dot(l2.lll_mapping, l2.matrix), l1.matrix)

        self.assertArrayAlmostEqual(np.dot(l2_fcoords, l2.matrix), np.dot(l1_fcoords, l1.matrix))

        lll_fcoords = l2.get_lll_frac_coords(l2_fcoords)

        self.assertArrayAlmostEqual(lll_fcoords, l1_fcoords)
        self.assertArrayAlmostEqual(l1.get_cartesian_coords(lll_fcoords), np.dot(lll_fcoords, l2.lll_matrix))

        self.assertArrayAlmostEqual(l2.get_frac_coords_from_lll(lll_fcoords), l2_fcoords)

    def test_get_miller_index_from_sites(self):
        # test on a cubic system
        m = Lattice.cubic(1)
        s1 = np.array([0.5, -1.5, 3])
        s2 = np.array([0.5, 3.0, -1.5])
        s3 = np.array([2.5, 1.5, -4.0])
        assert m.get_miller_index_from_coords([s1, s2, s3]) == (2, 1, 1)

        # test on a hexagonal system
        m = Lattice([[2.319, -4.01662582, 0.0], [2.319, 4.01662582, 0.0], [0.0, 0.0, 7.252]])

        s1 = np.array([2.319, 1.33887527, 6.3455])
        s2 = np.array([1.1595, 0.66943764, 4.5325])
        s3 = np.array([1.1595, 0.66943764, 0.9065])
        hkl = m.get_miller_index_from_coords([s1, s2, s3])
        assert hkl == (2, -1, 0)

        # test for previous failing structure
        m = Lattice([10, 0, 0, 0, 10, 0, 0, 0, 10])
        sites = [[0.5, 0.8, 0.8], [0.5, 0.4, 0.2], [0.5, 0.3, 0.7]]

        hkl = m.get_miller_index_from_coords(sites, coords_are_cartesian=False)
        assert hkl == (1, 0, 0)

        # test for more than 3 sites
        sites = [[0.5, 0.8, 0.8], [0.5, 0.4, 0.2], [0.5, 0.3, 0.7], [0.5, 0.1, 0.2]]

        hkl = m.get_miller_index_from_coords(sites, coords_are_cartesian=False)
        assert hkl == (1, 0, 0)

    def test_points_in_spheres(self):
        points = [[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]]
        lattice = Lattice.cubic(3)
        center_points = [[1.5, 1.5, 1.5]]
        nns = get_points_in_spheres(
            all_coords=np.array(points),
            center_coords=np.array(center_points),
            r=3,
            pbc=np.array([0, 0, 0], dtype=int),
            lattice=lattice,
            numerical_tol=1e-8,
        )
        assert len(nns[0]) == 2  # two neighbors

        nns = get_points_in_spheres(
            all_coords=np.array(points),
            center_coords=np.array(center_points),
            r=3,
            pbc=[1, 1, 1],
            lattice=lattice,
            numerical_tol=1e-8,
            return_fcoords=True,
        )
        assert len(nns[0]) == 12

        nns = get_points_in_spheres(
            all_coords=np.array(points),
            center_coords=np.array(center_points),
            r=3,
            pbc=np.array([True, False, False], dtype=int),
            lattice=lattice,
        )
        assert len(nns[0]) == 4

    def test_selling_dist(self):
        # verification process described here
        # https://github.com/materialsproject/pymatgen/pull/1888#issuecomment-818072164
        np.testing.assert_(Lattice.selling_dist(Lattice.cubic(5), Lattice.cubic(5)) == 0)
        hex_lattice = Lattice.hexagonal(5, 8)
        triclinic_lattice = Lattice.from_parameters(4, 10, 11, 100, 110, 80)
        np.testing.assert_allclose(Lattice.selling_dist(hex_lattice, triclinic_lattice), 76, rtol=0.1)
        np.testing.assert_allclose(
            Lattice.selling_dist(Lattice.tetragonal(10, 12), Lattice.tetragonal(10.1, 11.9)),
            3.7,
            rtol=0.1,
        )
        np.testing.assert_allclose(
            Lattice.selling_dist(Lattice.cubic(5), Lattice.from_parameters(8, 10, 12, 80, 90, 95)),
            115.6,
            rtol=0.1,
        )

    def test_selling_vector(self):
        a1 = 10
        np.testing.assert_array_almost_equal(
            Lattice.cubic(a1).selling_vector.round(4),
            np.array([0, 0, 0, -(a1**2), -(a1**2), -(a1**2)]),
        )
        a2, c2 = 5, 8
        np.testing.assert_array_almost_equal(
            Lattice.tetragonal(a2, c2).selling_vector.round(4),
            np.array([0, 0, 0, -(a2**2), -(a2**2), -(c2**2)]),
        )
        a3, b3, c3 = 4, 6, 7
        np.testing.assert_array_almost_equal(
            Lattice.orthorhombic(a3, b3, c3).selling_vector.round(4),
            np.array([0, 0, 0, -(a3**2), -(b3**2), -(c3**2)]),
        )

    def test_is_3d_periodic(self):
        assert self.cubic.is_3d_periodic
        assert not self.cubic_partial_pbc.is_3d_periodic


if __name__ == "__main__":
    import unittest

    unittest.main()
