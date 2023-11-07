from __future__ import annotations

import random
import unittest

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.core.lattice import Lattice
from pymatgen.util import coord


class TestCoordUtils:
    def test_get_linear_interpolated_value(self):
        x_vals = [0, 1, 2, 3, 4, 5]
        y_vals = [3, 6, 7, 8, 10, 12]
        assert coord.get_linear_interpolated_value(x_vals, y_vals, 3.6) == 9.2
        with pytest.raises(ValueError, match="x is out of range of provided x_values"):
            coord.get_linear_interpolated_value(x_vals, y_vals, 6)

    def test_in_coord_list(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        assert not coord.in_coord_list(coords, test_coord)
        assert coord.in_coord_list(coords, test_coord, atol=0.15)
        assert not coord.in_coord_list([0.99, 0.99, 0.99], test_coord, atol=0.15)

    def test_is_coord_subset(self):
        c1 = [0, 0, 0]
        c2 = [0, 1.2, -1]
        c3 = [3, 2, 1]
        c4 = [3 - 9e-9, 2 - 9e-9, 1 - 9e-9]
        assert coord.is_coord_subset([c1, c2, c3], [c1, c4, c2])
        assert coord.is_coord_subset([c1], [c2, c1])
        assert coord.is_coord_subset([c1, c2], [c2, c1])
        assert not coord.is_coord_subset([c1, c2], [c2, c3])
        assert not coord.is_coord_subset([c1, c2], [c2])

    def test_coord_list_mapping(self):
        c1 = [0, 0.124, 0]
        c2 = [0, 1.2, -1]
        c3 = [3, 2, 1]
        a = np.array([c1, c2])
        b = np.array([c3, c2, c1])
        inds = coord.coord_list_mapping(a, b)
        assert_allclose(a, b[inds])
        with pytest.raises(ValueError, match="not a subset of superset"):
            coord.coord_list_mapping([c1, c2], [c2, c3])
        with pytest.raises(ValueError, match="Something wrong with the inputs, likely duplicates in superset"):
            coord.coord_list_mapping([c2], [c2, c2])

    def test_coord_list_mapping_pbc(self):
        c1 = [0.1, 0.2, 0.3]
        c2 = [0.2, 0.3, 0.3]
        c3 = [0.5, 0.3, 0.6]
        c4 = [1.5, -0.7, -1.4]

        a = np.array([c1, c3, c2])
        b = np.array([c4, c2, c1])

        inds = coord.coord_list_mapping_pbc(a, b)
        diff = a - b[inds]
        diff -= np.round(diff)
        assert_allclose(diff, 0)
        with pytest.raises(ValueError, match="not a subset of superset"):
            coord.coord_list_mapping_pbc([c1, c2], [c2, c3])
        with pytest.raises(ValueError, match="Something wrong with the inputs, likely duplicates in superset"):
            coord.coord_list_mapping_pbc([c2], [c2, c2])
        coord.coord_list_mapping_pbc([c1, c2], [c2, c1], pbc=(False, False, False))
        with pytest.raises(ValueError, match="not a subset of superset"):
            coord.coord_list_mapping_pbc(a, b, pbc=(True, True, False))

    def test_find_in_coord_list(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        assert coord.find_in_coord_list(coords, test_coord).size == 0
        assert coord.find_in_coord_list(coords, test_coord, atol=0.15)[0] == 0
        assert coord.find_in_coord_list([0.99, 0.99, 0.99], test_coord, atol=0.15).size == 0
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.1, 0.1, 0.1]]
        assert_array_equal(coord.find_in_coord_list(coords, test_coord, atol=0.15), [0, 2])

    def test_all_distances(self):
        coords1 = [[0, 0, 0], [0.5, 0.5, 0.5]]
        coords2 = [[1, 2, -1], [1, 0, 0], [1, 0, 0]]
        result = [[2.44948974, 1, 1], [2.17944947, 0.8660254, 0.8660254]]
        assert_allclose(coord.all_distances(coords1, coords2), result, 4)

    def test_pbc_diff(self):
        assert_allclose(coord.pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]), [-0.2, -0.4, 0.2])
        assert_allclose(coord.pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9]), [-0.4, -0.4, 0.11])
        assert_allclose(coord.pbc_diff([0.1, 0.6, 1.01], [0.6, 0.1, 0.9]), [-0.5, 0.5, 0.11])
        assert_allclose(coord.pbc_diff([100.1, 0.2, 0.3], [0123123.4, 0.5, 502312.6]), [-0.3, -0.3, -0.3])
        assert_allclose(coord.pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9], pbc=(True, True, False)), [-0.2, -0.4, -0.8])
        assert_allclose(coord.pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9], pbc=(True, True, False)), [-0.4, -0.4, 0.11])

    def test_in_coord_list_pbc(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        assert not coord.in_coord_list_pbc(coords, test_coord)
        assert coord.in_coord_list_pbc(coords, test_coord, atol=0.15)
        test_coord = [0.99, 0.99, 0.99]
        assert not coord.in_coord_list_pbc(coords, test_coord, atol=0.01)

    def test_find_in_coord_list_pbc(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        assert coord.find_in_coord_list_pbc(coords, test_coord).size == 0
        assert coord.find_in_coord_list_pbc(coords, test_coord, atol=0.15)[0] == 0
        test_coord = [0.99, 0.99, 0.99]
        assert coord.find_in_coord_list_pbc(coords, test_coord, atol=0.02)[0] == 0
        test_coord = [-0.499, -0.499, -0.499]
        assert coord.find_in_coord_list_pbc(coords, test_coord, atol=0.01)[0] == 1
        test_coord = [-0.5, -0.5, -0.5]
        pbc = (True, True, False)
        assert coord.find_in_coord_list_pbc(coords, test_coord, pbc=pbc).size == 0
        test_coord = [-0.5, -0.5, 0.5]
        assert coord.find_in_coord_list_pbc(coords, test_coord, pbc=pbc)[0] == 1

    def test_is_coord_subset_pbc(self):
        c1 = [0, 0, 0]
        c2 = [0, 1.2, -1]
        c3 = [2.3, 0, 1]
        c4 = [1.3 - 9e-9, -1 - 9e-9, 1 - 9e-9]
        assert coord.is_coord_subset_pbc([c1, c2, c3], [c1, c4, c2])
        assert coord.is_coord_subset_pbc([c1], [c2, c1])
        assert coord.is_coord_subset_pbc([c1, c2], [c2, c1])
        assert not coord.is_coord_subset_pbc([c1, c2], [c2, c3])
        assert not coord.is_coord_subset_pbc([c1, c2], [c2])

        # test tolerances
        c5 = [0.1, 0.1, 0.2]
        atol1 = [0.25, 0.15, 0.15]
        atol2 = [0.15, 0.15, 0.25]
        assert not coord.is_coord_subset_pbc([c1], [c5], atol1)
        assert coord.is_coord_subset_pbc([c1], [c5], atol2)

        # test mask
        mask1 = [[True]]
        assert not coord.is_coord_subset_pbc([c1], [c5], atol2, mask1)
        mask2 = [[True, False]]
        assert coord.is_coord_subset_pbc([c1], [c2, c1], mask=mask2)
        assert not coord.is_coord_subset_pbc([c1], [c1, c2], mask=mask2)
        mask3 = [[False, True]]
        assert not coord.is_coord_subset_pbc([c1], [c2, c1], mask=mask3)
        assert coord.is_coord_subset_pbc([c1], [c1, c2], mask=mask3)

        # test pbc
        c5 = [1.3 - 9e-9, 9e-9, 1 - 9e-9]
        pbc = (True, False, True)
        assert coord.is_coord_subset_pbc([c1], [c2, c1], pbc=pbc)
        assert not coord.is_coord_subset_pbc([c1, c2, c3], [c1, c4, c2], pbc=pbc)
        assert coord.is_coord_subset_pbc([c1, c2, c3], [c1, c5, c2], pbc=pbc)

    def test_lattice_points_in_supercell(self):
        supercell = np.array([[1, 3, 5], [-3, 2, 3], [-5, 3, 1]])
        points = coord.lattice_points_in_supercell(supercell)
        assert len(points) == approx(abs(np.linalg.det(supercell)))
        assert np.min(points) >= -1e-10
        assert np.max(points) <= 1 - 1e-10

        supercell = np.array([[-5, -5, -3], [0, -4, -2], [0, -5, -2]])
        points = coord.lattice_points_in_supercell(supercell)
        assert len(points) == approx(abs(np.linalg.det(supercell)))
        assert np.min(points) >= -1e-10
        assert np.max(points) <= 1 - 1e-10

    def test_barycentric(self):
        # 2d test
        simplex1 = np.array([[0.3, 0.1], [0.2, -1.2], [1.3, 2.3]])
        pts1 = np.array([[0.6, 0.1], [1.3, 2.3], [0.5, 0.5], [0.7, 1]])
        output1 = coord.barycentric_coords(pts1, simplex1)
        # do back conversion to cartesian
        o_dot_s = np.sum(output1[:, :, None] * simplex1[None, :, :], axis=1)
        assert_allclose(pts1, o_dot_s)

        # do 3d tests
        simplex2 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, 0]])
        pts2 = np.array([[0, 0, 1], [0, 0.5, 0.5], [1.0 / 3, 1.0 / 3, 1.0 / 3]])
        output2 = coord.barycentric_coords(pts2, simplex2)
        assert_allclose(output2[1], [0.5, 0.5, 0, 0])
        # do back conversion to cartesian
        o_dot_s = np.sum(output2[:, :, None] * simplex2[None, :, :], axis=1)
        assert_allclose(pts2, o_dot_s)
        # test single point
        assert_allclose(output2[2], coord.barycentric_coords(pts2[2], simplex2).squeeze())

    def test_pbc_shortest_vectors(self):
        frac_coords = [
            [0.3, 0.3, 0.5],
            [0.1, 0.1, 0.3],
            [0.9, 0.9, 0.8],
            [0.1, 0.0, 0.5],
            [0.9, 0.7, 0.0],
        ]

        lattice = Lattice.from_parameters(8, 8, 4, 90, 76, 58)
        expected = [
            [0.000, 3.015, 4.072, 3.519, 3.245],
            [3.015, 0.000, 3.207, 1.131, 4.453],
            [4.072, 3.207, 0.000, 2.251, 1.788],
            [3.519, 1.131, 2.251, 0.000, 3.852],
        ]

        vectors = coord.pbc_shortest_vectors(lattice, frac_coords[:-1], frac_coords)
        dists = np.sum(vectors**2, axis=-1) ** 0.5
        assert_allclose(dists, expected, 3)

        prev_threshold = coord.LOOP_THRESHOLD
        coord.LOOP_THRESHOLD = 0

        vectors = coord.pbc_shortest_vectors(lattice, frac_coords[:-1], frac_coords)
        dists = np.sum(vectors**2, axis=-1) ** 0.5
        assert_allclose(dists, expected, 3)

        coord.LOOP_THRESHOLD = prev_threshold

        lattice_pbc = Lattice.from_parameters(8, 8, 4, 90, 76, 58, pbc=(True, True, False))
        expected_pbc = np.array(
            [
                [0.000, 3.015, 4.072, 3.519, 4.089],
                [3.015, 0.000, 3.207, 1.131, 4.453],
                [4.072, 3.207, 0.000, 2.251, 3.578],
                [3.519, 1.131, 2.251, 0.000, 4.235],
            ]
        )
        vectors = coord.pbc_shortest_vectors(lattice_pbc, frac_coords[:-1], frac_coords)
        dists = np.sum(vectors**2, axis=-1) ** 0.5
        assert_allclose(dists, expected_pbc, 3)

    def test_get_angle(self):
        v1 = (1, 0, 0)
        v2 = (1, 1, 1)
        assert coord.get_angle(v1, v2) == approx(54.7356103172)
        assert coord.get_angle(v1, v2, units="radians") == approx(0.9553166181245092)


class TestSimplex(unittest.TestCase):
    def setUp(self):
        coords = [[0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]]
        self.simplex = coord.Simplex(coords)

    def test_equal(self):
        c2 = list(self.simplex.coords)
        random.shuffle(c2)
        assert coord.Simplex(c2) == self.simplex

    def test_in_simplex(self):
        assert self.simplex.in_simplex([0.1, 0.1, 0.1])
        assert not self.simplex.in_simplex([0.6, 0.6, 0.6])
        for _ in range(10):
            coord = np.random.random_sample(size=3) / 3
            assert self.simplex.in_simplex(coord)

    def test_2d_triangle(self):
        simplex = coord.Simplex([[0, 1], [1, 1], [1, 0]])
        assert_allclose(simplex.bary_coords([0.5, 0.5]), [0.5, 0, 0.5])
        assert_allclose(simplex.bary_coords([0.5, 1]), [0.5, 0.5, 0])
        assert_allclose(simplex.bary_coords([0.5, 0.75]), [0.5, 0.25, 0.25])
        assert_allclose(simplex.bary_coords([0.75, 0.75]), [0.25, 0.5, 0.25])

        simplex = coord.Simplex([[1, 1], [1, 0]])
        with pytest.raises(ValueError, match="Simplex is not full-dimensional"):
            simplex.bary_coords([0.5, 0.5])

    def test_volume(self):
        # Should be value of a right tetrahedron.
        assert self.simplex.volume == approx(1 / 6)

    def test_str(self):
        assert str(self.simplex).startswith("3-simplex in 4D space")
        assert repr(self.simplex).startswith("3-simplex in 4D space")

    def test_bary_coords(self):
        s = coord.Simplex([[0, 2], [3, 1], [1, 0]])
        point = [0.7, 0.5]
        bc = s.bary_coords(point)
        assert_allclose(bc, [0.26, -0.02, 0.76])
        new_point = s.point_from_bary_coords(bc)
        assert_allclose(point, new_point)

    def test_intersection(self):
        # simple test, with 2 intersections at faces
        s = coord.Simplex([[0, 2], [3, 1], [1, 0]])
        point1 = [0.7, 0.5]
        point2 = [0.5, 0.7]
        intersections = s.line_intersection(point1, point2)
        expected = np.array([[1.13333333, 0.06666667], [0.8, 0.4]])
        assert_allclose(intersections, expected)

        # intersection through point and face
        point1 = [0, 2]  # simplex point
        point2 = [1, 1]  # inside simplex
        expected = np.array([[1.66666667, 0.33333333], [0, 2]])
        intersections = s.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

        # intersection through point only
        point1 = [0, 2]  # simplex point
        point2 = [0.5, 0.7]
        expected = np.array([[0, 2]])
        intersections = s.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

        # 3d intersection through edge and face
        point1 = [0.5, 0, 0]  # edge point
        point2 = [0.5, 0.5, 0.5]  # in simplex
        expected = np.array([[0.5, 0.25, 0.25], [0.5, 0.0, 0.0]])
        intersections = self.simplex.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

        # 3d intersection through edge only
        point1 = [0.5, 0, 0]  # edge point
        point2 = [0.5, 0.5, -0.5]  # outside simplex
        expected = np.array([[0.5, 0.0, 0.0]])
        intersections = self.simplex.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

        # coplanar to face (no intersection)
        point1 = [-1, 2]
        point2 = [0, 0]
        expected = np.array([])
        intersections = s.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

        # coplanar to face (with intersection line)
        point1 = [0, 2]  # simplex point
        point2 = [1, 0]
        expected = np.array([[1, 0], [0, 2]])
        intersections = s.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

        # coplanar to face (with intersection points)
        point1 = [0.1, 2]
        point2 = [1.1, 0]
        expected = np.array([[1.08, 0.04], [0.12, 1.96]])
        intersections = s.line_intersection(point1, point2)
        assert_allclose(intersections, expected)

    def test_to_json(self):
        assert isinstance(self.simplex.to_json(), str)
