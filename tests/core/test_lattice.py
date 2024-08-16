from __future__ import annotations

import itertools

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.core.lattice import Lattice, get_points_in_spheres
from pymatgen.core.operations import SymmOp
from pymatgen.util.testing import PymatgenTest


class TestLattice(PymatgenTest):
    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.cubic = self.lattice
        self.tetragonal = Lattice.tetragonal(10, 20)
        self.orthorhombic = Lattice.orthorhombic(10, 20, 30)
        self.monoclinic = Lattice.monoclinic(10, 20, 30, 66)
        self.hexagonal = Lattice.hexagonal(10, 20)
        self.rhombohedral = Lattice.rhombohedral(10, 77)

        self.cubic_partial_pbc = Lattice.cubic(10.0, pbc=(True, True, False))

        self.families = {}
        for name in ("cubic", "tetragonal", "orthorhombic", "monoclinic", "hexagonal", "rhombohedral"):
            self.families[name] = getattr(self, name)

    def test_init(self):
        len_a = 9.026
        lattice = Lattice.cubic(len_a)
        assert lattice is not None, "Initialization from new_cubic failed"
        assert_array_equal(lattice.pbc, (True, True, True))
        assert hash(lattice) == -6887896986157825384

        lattice2 = Lattice(np.eye(3) * len_a)
        for ii in range(3):
            for jj in range(3):
                assert lattice.matrix[ii][jj] == lattice2.matrix[ii][jj], "Inconsistent matrix from two inits!"
        assert_array_equal(self.cubic_partial_pbc.pbc, (True, True, False))

        for bad_pbc in [(True, True), (True, True, True, True), (True, True, 2)]:
            with pytest.raises(ValueError, match="pbc must be a tuple of three True/False values, got"):
                Lattice(np.eye(3), pbc=bad_pbc)

    def test_equal(self):
        assert self.cubic == self.cubic
        assert self.cubic == self.lattice
        for name1, latt1 in self.families.items():
            for name2, latt2 in self.families.items():
                assert (name1 == name2) == (latt1 == latt2)

        # ensure partial periodic boundaries is unequal to all full periodic boundaries
        assert not any(self.cubic_partial_pbc == family for family in self.families.values())

    def test_format(self):
        assert (
            format(self.lattice, ".3fl") == "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]"
        )
        lattice_str = """10.000 0.000 0.000
0.000 10.000 0.000
0.000 0.000 10.000"""
        assert format(self.lattice, ".3f") == lattice_str
        assert format(self.lattice, ".1fp") == "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}"

    def test_copy(self):
        cubic_copy = self.cubic.copy()
        assert cubic_copy == self.cubic
        assert self.cubic_partial_pbc.copy() == self.cubic_partial_pbc
        assert cubic_copy._matrix is not self.cubic._matrix

    def test_get_cartesian_or_frac_coord(self):
        coord = self.lattice.get_cartesian_coords([0.15, 0.3, 0.4])
        assert_allclose(coord, [1.5, 3.0, 4.0])
        assert_allclose(
            self.tetragonal.get_fractional_coords([12.12312, 45.2134, 1.3434]),
            [1.212312, 4.52134, 0.06717],
        )

        # Random testing that get_cart and get_frac coords reverses each other.
        rand_coord = np.random.default_rng().random(3)
        coord = self.tetragonal.get_cartesian_coords(rand_coord)
        frac_coord = self.tetragonal.get_fractional_coords(coord)
        assert_allclose(frac_coord, rand_coord)

    def test_get_vector_along_lattice_directions(self):
        lattice_mat = np.array([[0.5, 0.0, 0.0], [0.5, np.sqrt(3) / 2.0, 0.0], [0.0, 0.0, 1.0]])
        lattice = Lattice(lattice_mat)
        cart_coord = np.array([0.5, np.sqrt(3) / 4.0, 0.5])
        latt_coord = np.array([0.25, 0.5, 0.5])
        from_direct = lattice.get_fractional_coords(cart_coord) * lattice.lengths
        assert_allclose(lattice.get_vector_along_lattice_directions(cart_coord), from_direct)
        assert_allclose(lattice.get_vector_along_lattice_directions(cart_coord), latt_coord)
        assert_array_equal(
            lattice.get_vector_along_lattice_directions(cart_coord).shape,
            [
                3,
            ],
        )
        assert_array_equal(
            lattice.get_vector_along_lattice_directions(cart_coord.reshape([1, 3])).shape,
            [1, 3],
        )

    def test_d_hkl(self):
        cubic_copy = self.cubic.copy()
        hkl = (1, 2, 3)
        d_hkl = ((hkl[0] ** 2 + hkl[1] ** 2 + hkl[2] ** 2) / (cubic_copy.a**2)) ** (-1 / 2)
        assert d_hkl == cubic_copy.d_hkl(hkl)

    def test_reciprocal_lattice(self):
        recip_lattice = self.lattice.reciprocal_lattice
        assert_allclose(recip_lattice.matrix, 0.628319 * np.eye(3), atol=1e-6)
        assert_allclose(
            self.tetragonal.reciprocal_lattice.matrix,
            [[0.628319, 0.0, 0.0], [0.0, 0.628319, 0], [0.0, 0.0, 0.3141590]],
            atol=1e-6,
        )

        # Test the crystallographic version.
        recip_lattice_crystallographic = self.lattice.reciprocal_lattice_crystallographic
        assert_allclose(recip_lattice.matrix, recip_lattice_crystallographic.matrix * 2 * np.pi, 5)

    def test_static_methods(self):
        expected_lengths = [3.840198, 3.84019885, 3.8401976]
        expected_angles = [119.99998575, 90, 60.00000728]
        matrix = [
            [3.840198, 0.000000, 0.0000],
            [1.920099, 3.325710, 0.000000],
            [0.000000, -2.217138, 3.135509],
        ]
        # should give the lengths and angles above
        lattice = Lattice(matrix)
        assert lattice.lengths == approx(expected_lengths)
        assert lattice.angles == approx(expected_angles)

        lattice = Lattice.from_parameters(*lattice.parameters)
        assert lattice.lengths == approx(expected_lengths)
        assert lattice.angles == approx(expected_angles)

    def test_attributes(self):
        """Docstring for test_attributes."""
        lattice = Lattice.cubic(10.0)
        assert lattice.a == 10.0
        assert lattice.b == 10.0
        assert lattice.c == 10.0
        assert lattice.volume == 1000.0
        xyz = lattice.get_cartesian_coords([0.25, 0.35, 0.45])
        assert xyz[0] == 2.5
        assert xyz[1] == 3.5
        assert xyz[2] == 4.5

    def test_lattice_matrices(self):
        """If alpha == 90 and beta == 90, two matrices are identical."""

        def _identical(a, b, c, alpha, beta, gamma):
            mat1 = Lattice.from_parameters(a, b, c, alpha, beta, gamma, vesta=False).matrix
            mat2 = Lattice.from_parameters(a, b, c, alpha, beta, gamma, vesta=True).matrix
            # self.assertArrayAlmostEqual(mat1, mat2)
            return ((mat1 - mat2) ** 2).sum() < 1e-6

        assert _identical(2, 3, 4, 90, 90, 90)
        assert _identical(2, 3, 4, 90, 90, 80)
        assert _identical(2, 3, 4, 90, 90, 100)

        assert not _identical(2, 3, 4, 100, 90, 90)
        assert not _identical(2, 3, 4, 90, 100, 90)
        assert not _identical(2, 3, 4, 100, 100, 100)

    def test_get_lll_reduced_lattice(self):
        lattice = Lattice([1, 1, 1, -1, 0, 2, 3, 5, 6])
        reduced_latt = lattice.get_lll_reduced_lattice()

        expected = Lattice([[0, 1, 0], [1, 0, 1], [-2, 0, 1]])
        assert np.linalg.det(np.linalg.solve(expected.matrix, reduced_latt.matrix)) == approx(1)
        assert_allclose(sorted(reduced_latt.abc), sorted(expected.abc))
        assert reduced_latt.volume == approx(lattice.volume)
        lattice = [7.164750, 2.481942, 0.000000, -4.298850, 2.481942, 0.000000, 0.000000, 0.000000, 14.253000]
        expected = Lattice([-4.298850, 2.481942, 0.000000, 2.865900, 4.963884, 0.000000, 0.000000, 0.000000, 14.253000])
        reduced_latt = Lattice(lattice).get_lll_reduced_lattice()
        assert np.linalg.det(np.linalg.solve(expected.matrix, reduced_latt.matrix)) == approx(1)
        assert_allclose(sorted(reduced_latt.abc), sorted(expected.abc))

        expected = Lattice([0.0, 10.0, 10.0, 10.0, 10.0, 0.0, 30.0, -30.0, 40.0])

        lattice = np.array([100.0, 0.0, 10.0, 10.0, 10.0, 20.0, 10.0, 10.0, 10.0])
        lattice = lattice.reshape(3, 3)
        lattice = Lattice(lattice.T)
        reduced_latt = lattice.get_lll_reduced_lattice()
        assert np.linalg.det(np.linalg.solve(expected.matrix, reduced_latt.matrix)) == approx(1)
        assert_allclose(sorted(reduced_latt.abc), sorted(expected.abc))

        random_latt = Lattice(np.random.default_rng().random((3, 3)))
        if np.linalg.det(random_latt.matrix) > 1e-8:
            reduced_random_latt = random_latt.get_lll_reduced_lattice()
            assert reduced_random_latt.volume == approx(random_latt.volume)

    def test_get_niggli_reduced_lattice(self):
        lattice = Lattice.from_parameters(3, 5.196, 2, 103 + 55 / 60, 109 + 28 / 60, 134 + 53 / 60)
        reduced_cell = lattice.get_niggli_reduced_lattice()
        abc = reduced_cell.lengths
        angles = reduced_cell.angles
        assert abc == approx([2, 3, 3], abs=1e-3)
        assert angles == approx([116.382855225, 94.769790287999996, 109.466666667])

        lattice = Lattice(np.dot([[1, 1, 1], [1, 1, 0], [0, 1, 1]], 5 * np.eye(3)))
        reduced_cell = lattice.get_niggli_reduced_lattice()
        assert reduced_cell.lengths == approx([5, 5, 5])
        assert reduced_cell.angles == approx([90, 90, 90])

        lattice = Lattice([1.432950, 0.827314, 4.751000, -1.432950, 0.827314, 4.751000, 0.0, -1.654628, 4.751000])
        expected = [
            [-1.43295, -2.481942, 0.0],
            [-2.8659, 0.0, 0.0],
            [-1.43295, -0.827314, -4.751],
        ]
        assert_allclose(lattice.get_niggli_reduced_lattice().matrix, expected)

        lattice = Lattice.from_parameters(7.365450, 6.199506, 5.353878, 75.542191, 81.181757, 156.396627)
        expected = [
            [2.578932, 0.826965, 0.000000],
            [-0.831059, 2.067413, 1.547813],
            [-0.458407, -2.480895, 1.129126],
        ]
        assert_allclose(lattice.get_niggli_reduced_lattice().matrix, expected, atol=1e-5)

        lattice = Lattice([-0.2590, 1.1866, -0.1235, 2.2166, 1.0065, 0.7327, 1.1439, -0.4686, -0.0229])
        expected = [
            [-0.8849, -0.718, 0.1464],
            [0.1878, 0.7571, 0.902],
            [-0.4468, 0.4295, -1.0255],
        ]
        assert_allclose(lattice.get_niggli_reduced_lattice().matrix, expected, atol=1e-5)

    def test_find_mapping(self):
        matrix = [[0.1, 0.2, 0.3], [-0.1, 0.2, 0.7], [0.6, 0.9, 0.2]]
        lattice = Lattice(matrix)

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [2, 3, 3], 35)
        rot = op.rotation_matrix
        scale = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])

        latt2 = Lattice(np.dot(rot, np.dot(scale, matrix).T).T)
        mapping = latt2.find_mapping(lattice)
        assert isinstance(mapping, tuple)
        aligned_out, rot_out, scale_out = mapping
        assert abs(np.linalg.det(rot)) == approx(1)

        rotated = SymmOp.from_rotation_and_translation(rot_out).operate_multi(lattice.matrix)

        assert_allclose(rotated, aligned_out.matrix)
        assert_allclose(np.dot(scale_out, latt2.matrix), aligned_out.matrix)
        assert_allclose(aligned_out.parameters, lattice.parameters)
        assert not np.allclose(aligned_out.parameters, latt2.parameters)

    def test_find_all_mappings(self):
        matrix = [[0.1, 0.2, 0.3], [-0.1, 0.2, 0.7], [0.6, 0.9, 0.2]]
        lattice = Lattice(matrix)

        op = SymmOp.from_origin_axis_angle([0, 0, 0], [2, -1, 3], 40)
        rot = op.rotation_matrix
        scale = np.array([[0, 2, 0], [1, 1, 0], [0, 0, 1]])

        latt2 = Lattice(np.dot(rot, np.dot(scale, matrix).T).T)

        for aligned_out, rot_out, scale_out in lattice.find_all_mappings(latt2):
            assert_allclose(np.inner(latt2.matrix, rot_out), aligned_out.matrix, atol=1e-12)
            assert_allclose(np.dot(scale_out, lattice.matrix), aligned_out.matrix)
            assert_allclose(aligned_out.parameters, latt2.parameters)
            assert not np.allclose(aligned_out.parameters, lattice.parameters)

        lattice = Lattice.orthorhombic(9, 9, 5)
        assert len(list(lattice.find_all_mappings(lattice))) == 16

        # catch the singular matrix error
        lattice = Lattice.from_parameters(1, 1, 1, 10, 10, 10)
        for latt, _, _ in lattice.find_all_mappings(lattice, ltol=0.05, atol=11):
            assert isinstance(latt, Lattice)

    def test_mapping_symmetry(self):
        lattice = Lattice.cubic(1)
        l2 = Lattice.orthorhombic(1.1001, 1, 1)
        assert lattice.find_mapping(l2, ltol=0.1) == l2.find_mapping(lattice, ltol=0.1)
        assert l2.find_mapping(lattice, ltol=0.1) is None
        l2 = Lattice.orthorhombic(1.0999, 1, 1)
        mapping = l2.find_mapping(lattice, ltol=0.1)
        assert isinstance(mapping, tuple)
        assert len(mapping) == 3
        assert lattice.find_mapping(l2, ltol=0.1) is not None

    def test_as_from_dict(self):
        dct = self.tetragonal.as_dict()
        expected_keys = {"matrix", "@class", "pbc", "@module"}
        assert {*dct} == expected_keys
        tetragonal = Lattice.from_dict(dct)
        assert tetragonal.abc == self.tetragonal.abc
        assert tetragonal.angles == self.tetragonal.angles

        # Make sure old style dicts work.
        dct = self.tetragonal.as_dict(verbosity=1)
        assert {*dct} == expected_keys | {"a", "b", "c", "alpha", "beta", "gamma", "volume"}
        del dct["matrix"]
        tetragonal = Lattice.from_dict(dct)
        assert tetragonal.abc == self.tetragonal.abc
        assert tetragonal.angles == self.tetragonal.angles

    def test_parameters(self):
        params_dict = self.tetragonal.params_dict
        assert params_dict == {"a": 10, "b": 10, "c": 20, "alpha": 90, "beta": 90, "gamma": 90}
        params = self.tetragonal.parameters
        assert params == (10, 10, 20, 90, 90, 90)
        assert tuple(params_dict.values()) == params

    def test_scale(self):
        new_volume = 10
        for lattice in self.families.values():
            new_lattice = lattice.scale(new_volume)
            assert new_lattice.volume == approx(new_volume)
            assert_allclose(new_lattice.angles, lattice.angles)

    def test_get_wigner_seitz_cell(self):
        ws_cell = Lattice(np.diag([10, 5, 1])).get_wigner_seitz_cell()
        assert len(ws_cell) == 6
        for vec in ws_cell[3]:
            assert [abs(i) for i in vec] == [5.0, 2.5, 0.5]

    def test_dot_and_norm(self):
        frac_basis = np.eye(3)

        for lattice in self.families.values():
            assert_allclose(lattice.norm(lattice.matrix, frac_coords=False), lattice.abc, 5)
            assert_allclose(lattice.norm(frac_basis), lattice.abc, 5)
            for i, vec in enumerate(frac_basis):
                length = lattice.norm(vec)
                assert_allclose(length[0], lattice.abc[i], 5)
                # We always get a ndarray.
                assert length.shape == (1,)

        # Passing complex arrays should raise TypeError
        with pytest.raises(TypeError, match="Complex array"):
            lattice.norm(np.zeros(3, dtype=np.complex128))

        # Cannot reshape the second argument.
        with pytest.raises(ValueError, match="cannot reshape array"):
            lattice.dot(np.zeros(6), np.zeros(8))

        # Passing vectors of different length is invalid.
        with pytest.raises(ValueError, match="Coordinates must have same length"):
            lattice.dot(np.zeros(3), np.zeros(6))

    def test_get_points_in_sphere(self):
        # This is a non-niggli representation of a cubic lattice
        lattice = Lattice([[1, 5, 0], [0, 1, 0], [5, 0, 1]])
        # evenly spaced points array between 0 and 1
        points = np.array(list(itertools.product(range(5), repeat=3))) / 5
        points = lattice.get_fractional_coords(points)

        # Test getting neighbors within 1 neighbor distance of the origin
        frac_coords, dists, indices, images = lattice.get_points_in_sphere(
            points, [0, 0, 0], 0.20001, zip_results=False
        )
        assert len(frac_coords) == 7  # There are 7 neighbors
        assert np.isclose(dists, 0.2).sum() == 6  # 6 are at 0.2
        assert np.isclose(dists, 0).sum() == 1  # 1 is at 0
        assert len(set(indices)) == 7  # They have unique indices
        assert images[np.isclose(dists, 0)].tolist() == [[0, 0, 0]]

        # More complicated case, using the zip output
        result = lattice.get_points_in_sphere(points, [0.5, 0.5, 0.5], 1.0001)
        assert isinstance(result, tuple)
        assert len(result) == 552
        assert len(result[0]) == 4  # coords, dists, indices, supercell

        # test pbc
        latt_pbc = Lattice([[1, 5, 0], [0, 1, 0], [5, 0, 1]], pbc=(True, True, False))
        frac_coords, dists, indices, images = latt_pbc.get_points_in_sphere(
            points, [0, 0, 0], 0.20001, zip_results=False
        )
        assert len(frac_coords) == 6

        # ensure consistent return type if zip_results=False and no points in sphere are found
        # https://github.com/materialsproject/pymatgen/issues/3794
        result = lattice.get_points_in_sphere(points, [0.5, 0.5, 0.5], 0.0001, zip_results=False)
        assert isinstance(result, tuple)
        assert len(result) == 4
        assert all(len(arr) == 0 for arr in result)
        types = {*map(type, result)}
        assert types == {np.ndarray}, f"Expected only np.ndarray, got {types}"

    def test_get_all_distances(self):
        frac_coords = np.array(
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
        output = lattice.get_all_distances(frac_coords, frac_coords)
        assert_allclose(output, expected, 3)
        # test just one input point
        output2 = lattice.get_all_distances(frac_coords[0], frac_coords)
        assert_allclose(output2, [expected[0]], 2)
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
        output3 = lattice_pbc.get_all_distances(frac_coords[:-1], frac_coords)
        assert_allclose(output3, expected_pbc, 3)

    def test_monoclinic(self):
        assert self.monoclinic.angles == approx([90, 66, 90])
        assert self.monoclinic.lengths == approx([10, 20, 30])

    def test_is_hexagonal(self):
        assert not self.cubic.is_hexagonal()
        assert not self.tetragonal.is_hexagonal()
        assert not self.orthorhombic.is_hexagonal()
        assert not self.monoclinic.is_hexagonal()
        assert not self.rhombohedral.is_hexagonal()
        assert self.hexagonal.is_hexagonal()

    def test_get_distance_and_image(self):
        dist, image = self.cubic.get_distance_and_image([0, 0, 0.1], [0, 0.0, 0.9])
        assert dist == approx(2)
        assert_allclose(image, [0, 0, -1])

    def test_get_distance_and_image_strict(self):
        rng = np.random.default_rng()
        for _ in range(10):
            lengths = rng.integers(1, 100, 3)
            lattice = rng.random((3, 3)) * lengths
            lattice = Lattice(lattice)

            f1 = rng.random(3)
            f2 = rng.random(3)

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
                assert_allclose(min_image_dist[1], pmg_result[1])

    def test_lll_basis(self):
        a = np.array([1.0, 0.1, 0.0])
        b = np.array([0.0, 2.0, 0.0])
        c = np.array([0.0, 0.0, 3.0])

        l1 = Lattice([a, b, c])
        l2 = Lattice([a + b, b + c, c])

        cart_coords = np.array([[1, 1, 2], [2, 2, 1.5]])
        l1_frac_coords = l1.get_fractional_coords(cart_coords)
        l2_frac_coords = l2.get_fractional_coords(cart_coords)

        assert_allclose(l1.matrix, l2.lll_matrix)
        assert_allclose(np.dot(l2.lll_mapping, l2.matrix), l1.matrix)

        assert_allclose(np.dot(l2_frac_coords, l2.matrix), np.dot(l1_frac_coords, l1.matrix))

        lll_frac_coords = l2.get_lll_frac_coords(l2_frac_coords)

        assert_allclose(lll_frac_coords, l1_frac_coords)
        assert_allclose(l1.get_cartesian_coords(lll_frac_coords), np.dot(lll_frac_coords, l2.lll_matrix))

        assert_allclose(l2.get_frac_coords_from_lll(lll_frac_coords), l2_frac_coords)

    def test_get_miller_index_from_sites(self):
        # test on a cubic system
        cubic = Lattice.cubic(1)
        s1 = np.array([0.5, -1.5, 3])
        s2 = np.array([0.5, 3.0, -1.5])
        s3 = np.array([2.5, 1.5, -4.0])
        assert cubic.get_miller_index_from_coords([s1, s2, s3]) == (2, 1, 1)

        # test on a hexagonal system
        hexagonal = Lattice([[2.319, -4.01662582, 0.0], [2.319, 4.01662582, 0.0], [0.0, 0.0, 7.252]])

        s1 = np.array([2.319, 1.33887527, 6.3455])
        s2 = np.array([1.1595, 0.66943764, 4.5325])
        s3 = np.array([1.1595, 0.66943764, 0.9065])
        hkl = hexagonal.get_miller_index_from_coords([s1, s2, s3])
        assert hkl == (2, -1, 0)

        # test for previous failing structure
        cubic_from_flat = Lattice([10, 0, 0, 0, 10, 0, 0, 0, 10])
        sites = [[0.5, 0.8, 0.8], [0.5, 0.4, 0.2], [0.5, 0.3, 0.7]]

        hkl = cubic_from_flat.get_miller_index_from_coords(sites, coords_are_cartesian=False)
        assert hkl == (1, 0, 0)

        # test for more than 3 sites
        sites = [[0.5, 0.8, 0.8], [0.5, 0.4, 0.2], [0.5, 0.3, 0.7], [0.5, 0.1, 0.2]]

        hkl = cubic_from_flat.get_miller_index_from_coords(sites, coords_are_cartesian=False)
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
        assert Lattice.selling_dist(Lattice.cubic(5), Lattice.cubic(5)) == 0
        hex_lattice = Lattice.hexagonal(5, 8)
        triclinic_lattice = Lattice.from_parameters(4, 10, 11, 100, 110, 80)
        assert Lattice.selling_dist(hex_lattice, triclinic_lattice) == approx(76, abs=0.1)
        assert Lattice.selling_dist(Lattice.tetragonal(10, 12), Lattice.tetragonal(10.1, 11.9)) == approx(3.7, abs=0.1)
        assert Lattice.selling_dist(Lattice.cubic(5), Lattice.from_parameters(8, 10, 12, 80, 90, 95)) == approx(
            125.99, abs=0.1
        )

    def test_selling_vector(self):
        a1 = 10
        assert_allclose(
            Lattice.cubic(a1).selling_vector.round(4),
            [0, 0, 0, -(a1**2), -(a1**2), -(a1**2)],
        )
        a2, c2 = 5, 8
        assert_allclose(
            Lattice.tetragonal(a2, c2).selling_vector.round(4),
            [0, 0, 0, -(a2**2), -(a2**2), -(c2**2)],
        )
        a3, b3, c3 = 4, 6, 7
        assert_allclose(
            Lattice.orthorhombic(a3, b3, c3).selling_vector.round(4),
            [0, 0, 0, -(a3**2), -(b3**2), -(c3**2)],
        )

    def test_is_3d_periodic(self):
        assert self.cubic.is_3d_periodic
        assert not self.cubic_partial_pbc.is_3d_periodic
