from __future__ import annotations

import itertools

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import (
    Plane,
    anticlockwise_sort,
    anticlockwise_sort_indices,
    changebasis,
    collinear,
    diamond_functions,
    function_comparison,
    get_lower_and_upper_f,
    is_anion_cation_bond,
    matrixTimesVector,
    quarter_ellipsis_functions,
    rectangle_surface_intersection,
    rotateCoords,
    rotateCoordsOpt,
    separation_in_list,
    solid_angle,
    sort_separation,
    sort_separation_tuple,
    spline_functions,
    vectorsToMatrix,
)
from pymatgen.util.testing import MatSciTest

__author__ = "David Waroquiers"


class TestPlanesUtils(MatSciTest):
    def setup_method(self):
        # Test of plane 4x + 2y - 4z + 3 = 0 (used in most test cases)
        self.expected_coefficients = np.array([4, 2, -4, 3], float)
        self.p1 = np.array([0, 0, 0.75])
        self.p2 = np.array([-0.75, 0, 0])
        self.p3 = np.array([0, -1.5, 0])
        self.plane = Plane.from_3points(self.p1, self.p2, self.p3)

    def test_factors_abcd_normal_vector(self):
        factors = self.plane.coefficients / self.expected_coefficients
        assert_allclose([factors[0]] * 4, list(factors))
        assert_allclose([2 / 3, 1 / 3, -2 / 3], self.plane.normal_vector)

    def test_from_npoints_plane(self):
        best_fits = ["least_square_distance", "maximum_distance"]
        delta = 0.0001
        for best_fit in best_fits:
            plane = Plane.from_npoints([self.p1, self.p2, self.p3], best_fit=best_fit)
            assert self.plane.is_same_plane_as(plane)
            points = [
                np.array([5.1, 0.3, -2.3]),
                np.array([-2, 4.3, -6.3]),
                np.array([3.1, 2.3, -21.3]),
                np.array([-2, -0.5, 0.05]),
                np.array([11, 12, -13]),
                np.array([10, 8.3, -6.32]),
            ]
            plane = Plane.from_npoints(points, best_fit=best_fit)
            fit_error_plane = plane.fit_error(points, fit=best_fit)
            coeffs = [
                [1 + delta, 1, 1, 1],
                [1, 1 + delta, 1, 1],
                [1, 1, 1 + delta, 1],
                [1, 1, 1, 1 + delta],
                [1 - delta, 1 + delta, 1 - delta, 1 + delta],
            ]
            for coeff in coeffs:
                plane_changed = Plane.from_coefficients(
                    coeff[0] * plane.a,
                    coeff[1] * plane.b,
                    coeff[2] * plane.c,
                    coeff[3] * plane.d,
                )
                fit_error = plane_changed.fit_error(points, fit=best_fit)
                assert fit_error > fit_error_plane
            coeff = [-2.1, -2.1, -2.1, -2.1]
            plane_not_changed = Plane.from_coefficients(
                coeff[0] * plane.a,
                coeff[1] * plane.b,
                coeff[2] * plane.c,
                coeff[3] * plane.d,
            )
            fit_error = plane_not_changed.fit_error(points, fit=best_fit)
            assert fit_error == approx(fit_error_plane)

    def test_is_in_plane(self):
        assert self.plane.is_in_plane(self.p1, 0.001)
        assert self.plane.is_in_plane(self.p2, 0.001)
        assert self.plane.is_in_plane(self.p3, 0.001)
        assert not self.plane.is_in_plane(np.zeros(3), 0.001)
        assert self.plane.is_in_plane(np.array([1, 1, 2.25]), 0.001)
        assert self.plane.is_in_plane(np.array([1, 1, 2.22]), 0.1)
        assert not self.plane.is_in_plane(np.array([1, 1, 2.22]), 0.001)
        assert self.plane.is_in_plane(self.p1 + self.plane.normal_vector * 1, 1.000001)
        assert not self.plane.is_in_plane(self.p1 + self.plane.normal_vector * 1.00001, 1)
        assert not self.plane.is_in_plane(self.p1 + self.plane.normal_vector * 1, 0.999999)
        assert self.plane.is_in_plane(self.plane.p1, 0.00001)
        assert self.plane.is_in_plane(self.plane.p2, 0.00001)
        assert self.plane.is_in_plane(self.plane.p3, 0.00001)

    def test_normal_vector_is_normed(self):
        assert np.isclose(np.linalg.norm(self.plane.normal_vector), 1)

    def test_orthonormal_vectors(self):
        ortho = self.plane.orthonormal_vectors()
        assert np.isclose(np.dot(ortho[0], self.plane.normal_vector), 0)
        assert np.isclose(np.dot(ortho[1], self.plane.normal_vector), 0)
        assert np.isclose(np.dot(ortho[2], self.plane.normal_vector), 1)
        assert np.isclose(np.dot(ortho[0], ortho[1]), 0)
        assert np.isclose(np.dot(ortho[1], ortho[2]), 0)
        assert np.isclose(np.dot(ortho[2], ortho[0]), 0)
        assert_allclose(np.cross(ortho[0], ortho[1]), ortho[2])
        assert_allclose(np.cross(ortho[0], ortho[1]), self.plane.normal_vector)
        assert_allclose(np.cross(ortho[1], ortho[2]), ortho[0])
        assert_allclose(np.cross(ortho[2], ortho[0]), ortho[1])
        assert not np.allclose(np.cross(ortho[1], ortho[0]), ortho[2])
        assert_allclose(np.cross(ortho[0], ortho[1]), self.plane.normal_vector)

    def test_plane_comparison(self):
        plane_test_1 = Plane.from_coefficients(4, 2, -4, 3)
        assert self.plane.is_same_plane_as(plane_test_1)
        plane_test_2 = Plane.from_coefficients(-4, -2, 4, -3)
        assert self.plane.is_same_plane_as(plane_test_2)
        plane_test_3 = Plane.from_coefficients(-12, -6, 12, -9)
        assert self.plane.is_same_plane_as(plane_test_3)
        plane_test_4 = Plane.from_coefficients(3, 0, 2, 4)
        assert not self.plane.is_same_plane_as(plane_test_4)

    def test_plane_is_in_list_of_planes(self):
        plane_test_1 = Plane.from_coefficients(-8.1, 2, -4, 3)
        plane_test_2 = Plane.from_coefficients(0, -2, 4, 0)
        plane_test_3 = Plane.from_coefficients(-12, -6, 12, -9)
        plane_test_4 = Plane.from_coefficients(3, 0, 0, 4)
        plane_list = [plane_test_1, plane_test_2, plane_test_3, plane_test_4]
        assert self.plane.is_in_list(plane_list)
        plane_list = [plane_test_1, plane_test_2, plane_test_4]
        assert not self.plane.is_in_list(plane_list)

    def test_plane_3_coefficients(self):
        plane_1 = Plane.from_coefficients(0, 2, -1, 3)
        assert plane_1.is_in_plane(plane_1.p1, 0.000001)
        assert plane_1.is_in_plane(plane_1.p2, 0.000001)
        assert plane_1.is_in_plane(plane_1.p3, 0.000001)
        plane_2 = Plane.from_coefficients(12, 0, 2, -4)
        assert plane_2.is_in_plane(plane_2.p1, 0.000001)
        assert plane_2.is_in_plane(plane_2.p2, 0.000001)
        assert plane_2.is_in_plane(plane_2.p3, 0.000001)
        plane_3 = Plane.from_coefficients(-8, 8, 0, 0)
        assert plane_3.is_in_plane(plane_3.p1, 0.000001)
        assert plane_3.is_in_plane(plane_3.p2, 0.000001)
        assert plane_3.is_in_plane(plane_3.p3, 0.000001)

    def test_plane_2_coefficients(self):
        plane_1 = Plane.from_coefficients(-21, 0, 0, 3)
        assert plane_1.is_in_plane(plane_1.p1, 0.000001)
        assert plane_1.is_in_plane(plane_1.p2, 0.000001)
        assert plane_1.is_in_plane(plane_1.p3, 0.000001)
        plane_2 = Plane.from_coefficients(0, 4, 0, -4)
        assert plane_2.is_in_plane(plane_2.p1, 0.000001)
        assert plane_2.is_in_plane(plane_2.p2, 0.000001)
        assert plane_2.is_in_plane(plane_2.p3, 0.000001)
        plane_3 = Plane.from_coefficients(0, 0, 3, 1)
        assert plane_3.is_in_plane(plane_3.p1, 0.000001)
        assert plane_3.is_in_plane(plane_3.p2, 0.000001)
        assert plane_3.is_in_plane(plane_3.p3, 0.000001)

    def test_indices_separate(self):
        # Test with the common test plane
        point_1 = np.array([0, 0, 0], float)
        point_2 = np.array([0, 0, 0.75], float)
        point_3 = np.array([-0.75, 0, 0], float)
        point_4 = np.array([1, 0, 0], float)
        point_5 = np.array([0, -1.5, 0], float)
        point_6 = np.array([10, 2, -20], float)
        point_7 = np.array([10, 10, 10], float)
        point_8 = np.array([100, 0, 0], float)
        plist = [point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8]
        sep = self.plane.indices_separate(plist, 0.000001)
        assert len(sep[0]) == 0
        assert len(sep[1]) == 3
        assert len(sep[2]) == 5
        assert_allclose(sep[1], [1, 2, 4])
        assert_allclose(sep[2], [0, 3, 5, 6, 7])
        sep = self.plane.indices_separate(plist, 10)
        assert len(sep[0]) == 0
        assert len(sep[1]) == 6
        assert len(sep[2]) == 2
        assert_allclose(sep[1], [0, 1, 2, 3, 4, 6])
        assert_allclose(sep[2], [5, 7])
        sep = self.plane.indices_separate(plist, 100000)
        assert len(sep[0]) == 0
        assert len(sep[1]) == 8
        assert len(sep[2]) == 0
        # Test with 2 coeff facets (Normal vector = [1, 0, 0] or [0, 1, 0] or [0, 0, 1])
        # Plane x-2=0 (perpendicular to x)
        plane = Plane.from_coefficients(-4, 0, 0, 8)
        sep = plane.indices_separate(plist, 0.00001)
        assert sep[0] == [0, 1, 2, 3, 4]
        assert sep[1] == []
        assert sep[2] == [5, 6, 7]
        sep = plane.indices_separate(plist, 1)
        assert sep[0] == [0, 1, 2, 4]
        assert sep[1] == [3]
        assert sep[2] == [5, 6, 7]
        # Plane 2y+1=0 (perpendicular to y)
        plane = Plane.from_coefficients(0, 2, 0, 1)
        sep = plane.indices_separate(plist, 0.00001)
        assert sep[0] == [4]
        assert sep[1] == []
        assert sep[2] == [0, 1, 2, 3, 5, 6, 7]
        sep = plane.indices_separate(plist, 1)
        assert sep[0] == []
        assert sep[1] == [0, 1, 2, 3, 4, 7]
        assert sep[2] == [5, 6]
        # Plane 4z-3=0 (perpendicular to z)
        plane = Plane.from_coefficients(0, 0, -4, 3)
        sep = plane.indices_separate(plist, 0.00001)
        assert sep[0] == [0, 2, 3, 4, 5, 7]
        assert sep[1] == [1]
        assert sep[2] == [6]
        sep = plane.indices_separate(plist, 0.75)
        assert sep[0] == [5]
        assert sep[1] == [0, 1, 2, 3, 4, 7]
        assert sep[2] == [6]
        # Test with 3 coeff facets (Normal vector = [0, a, b] or [a, 0, b] or [a, b, 0])
        # Plane 2y-z+4=0
        plane = Plane.from_coefficients(0, 2, -1, 0)
        sep = plane.indices_separate(plist, 0.00001)
        assert sep[0] == [1, 4]
        assert sep[1] == [0, 2, 3, 7]
        assert sep[2] == [5, 6]
        sep = plane.indices_separate(plist, 0.75)
        assert sep[0] == [4]
        assert sep[1] == [0, 1, 2, 3, 7]
        assert sep[2] == [5, 6]
        # Plane 2y-z+4=0
        plane = Plane.from_coefficients(4, 0, -2, -20)
        sep = plane.indices_separate(plist, 0.00001)
        assert sep[0] == [0, 1, 2, 3, 4]
        assert sep[1] == [6]
        assert sep[2] == [5, 7]
        # Plane 2y-z+4=0
        plane = Plane.from_coefficients(-2, 9, 0, 2)
        sep = plane.indices_separate(plist, 0.00001)
        assert sep[0] == [0, 1, 2, 6]
        assert sep[1] == [3, 5]
        assert sep[2] == [4, 7]

    def test_distances(self):
        # Test with the common test plane
        point_1 = np.array([0, 0, 0], float)
        point_2 = np.array([0, 0, 0.75], float)
        point_3 = np.array([-0.75, 0, 0], float)
        point_4 = np.array([1, 0, 0], float)
        point_5 = np.array([0, -1.5, 0], float)
        point_6 = np.array([10, 2, -20], float)
        point_7 = np.array([10, 10, 10], float)
        point_8 = np.array([100, 0, 0], float)
        plist = [point_1, point_2, point_4, point_6, point_7, point_8]
        distances, indices_sorted = self.plane.distances_indices_sorted(plist)
        assert_allclose(
            distances,
            [
                0.5,
                0,
                1.16666666666666,
                21.1666666666666,
                3.8333333333333,
                67.1666666666666,
            ],
            atol=1e-9,
        )
        assert indices_sorted == [1, 0, 2, 4, 3, 5]
        # Plane 2y+1=0 (perpendicular to y)
        plane = Plane.from_coefficients(0, 2, 0, 1)
        plist = [point_1, point_5, point_6, point_7]
        distances, indices_sorted = plane.distances_indices_sorted(plist)
        assert_allclose(distances, [0.5, -1, 2.5, 10.5])
        assert indices_sorted == [0, 1, 2, 3]
        plist = [point_1, point_5, point_6, point_7]
        distances, indices_sorted = plane.distances_indices_sorted(plist)
        assert_allclose(distances, [0.5, -1, 2.5, 10.5])
        assert indices_sorted == [0, 1, 2, 3]
        distances, indices_sorted = plane.distances_indices_sorted(plist, sign=True)
        assert_allclose(distances, [0.5, -1, 2.5, 10.5])
        assert indices_sorted == [(0, 1), (1, -1), (2, 1), (3, 1)]
        plist = [point_5, point_1, point_6, point_7]
        distances, indices_sorted, groups = plane.distances_indices_groups(plist)
        assert_allclose(distances, [-1, 0.5, 2.5, 10.5])
        assert indices_sorted == [1, 0, 2, 3]
        assert groups == [[1, 0], [2], [3]]
        plist = [point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8]
        distances, indices_sorted = plane.distances_indices_sorted(plist)
        assert_allclose(distances, [0.5, 0.5, 0.5, 0.5, -1, 2.5, 10.5, 0.5])
        assert set(indices_sorted[:5]) == {0, 1, 2, 3, 7}
        # Plane z=0 (plane xy)
        plane = Plane.from_coefficients(0, 0, 1, 0)
        zzs = [0.1, -0.2, 0.7, -2.1, -1.85, 0, -0.71, -0.82, -6.5, 1.8]
        plist = []
        rng = np.random.default_rng()
        for zz in zzs:
            plist.append([rng.uniform(-20, 20), rng.uniform(-20, 20), zz])
        distances, indices_sorted, groups = plane.distances_indices_groups(points=plist, delta=0.25)
        assert indices_sorted == [5, 0, 1, 2, 6, 7, 9, 4, 3, 8]
        assert groups == [[5, 0, 1], [2, 6, 7], [9, 4, 3], [8]]
        assert_allclose(distances, zzs)
        distances, indices_sorted, groups = plane.distances_indices_groups(points=plist, delta_factor=0.1)
        assert indices_sorted == [5, 0, 1, 2, 6, 7, 9, 4, 3, 8]
        assert groups == [[5, 0, 1, 2, 6, 7], [9, 4, 3], [8]]
        assert_allclose(distances, zzs)
        distances, indices_sorted, groups = plane.distances_indices_groups(points=plist, delta_factor=0.1, sign=True)
        assert indices_sorted == [
            (5, 0),
            (0, 1),
            (1, -1),
            (2, 1),
            (6, -1),
            (7, -1),
            (9, 1),
            (4, -1),
            (3, -1),
            (8, -1),
        ]
        assert groups == [
            [(5, 0), (0, 1), (1, -1), (2, 1), (6, -1), (7, -1)],
            [(9, 1), (4, -1), (3, -1)],
            [(8, -1)],
        ]
        assert_allclose(distances, zzs)

    def test_projections(self):
        # Projections of points that are already on the plane
        expected_projected_points = [
            self.p1,
            self.p2,
            self.p3,
            self.plane.p1,
            self.plane.p2,
            self.plane.p3,
        ]
        projected_points = self.plane.projectionpoints(expected_projected_points)
        projected_2d = self.plane.project_and_to2dim(expected_projected_points, "mean")
        for ii, pp in enumerate(expected_projected_points):
            assert_allclose(pp, projected_points[ii], atol=1e-9)
        for i1, i2 in itertools.combinations(list(range(len(expected_projected_points))), 2):
            assert np.isclose(
                np.linalg.norm(expected_projected_points[i1] - expected_projected_points[i2]),
                np.linalg.norm(projected_2d[i1] - projected_2d[i2]),
            )
        # Projections of random points (check on distances between the 3D points and the 2D points)
        points_to_project = [
            np.array([5.1, 0.3, -2.3]),
            np.array([-2, 4.3, -6.3]),
            np.array([3.1, 2.3, -21.3]),
            np.array([-2, -0.5, 0.05]),
            np.array([11, 12, -13]),
            np.array([10, 8.3, -6.32]),
        ]
        projected_points = self.plane.projectionpoints(points_to_project)
        mean_point = np.array([np.mean([pp[ii] for pp in points_to_project]) for ii in range(3)])
        projected_2d = self.plane.project_and_to2dim(points_to_project, "mean")
        projected_2d_bis = self.plane.project_and_to2dim(points_to_project, mean_point)
        for ii, pp in enumerate(projected_2d):
            assert_allclose(pp, projected_2d_bis[ii])
        for i1, i2 in itertools.combinations(list(range(len(projected_points))), 2):
            assert np.isclose(
                np.linalg.norm(projected_points[i1] - projected_points[i2]),
                np.linalg.norm(projected_2d[i1] - projected_2d[i2]),
            )
        for pp in points_to_project:
            projected_2d = self.plane.project_and_to2dim([pp], pp)
            assert_allclose(projected_2d[0], 0)
        # Check some specific projections
        points = [
            np.zeros(3, float),
            np.array([10, 10, 10], float),
            np.array([1.2, 2.3, 3.4], float),
            np.array([-1, -2, -3], float),
            np.array([-1, 1, -1], float),
        ]
        projected_points = self.plane.projectionpoints(points)
        expected_projected_points = [
            np.array([-0.33333333, -0.16666667, 0.33333333]),
            np.array([7.44444444, 8.72222222, 12.55555556]),
            np.array([1.33333333, 2.36666667, 3.26666667]),
            np.array([-1.77777778, -2.38888889, -2.22222222]),
            np.array([-1.55555556, 0.72222222, -0.44444444]),
        ]
        for ii, pp in enumerate(projected_points):
            assert_allclose(pp, expected_projected_points[ii])
            assert self.plane.is_in_plane(pp, 0.0000001)
        projected_2d_points_000 = self.plane.project_and_to2dim(points, [0, 0, 0])
        projected_2d_points_mean = self.plane.project_and_to2dim(points, "mean")
        for i1, i2 in itertools.combinations(list(range(len(projected_2d_points_000))), 2):
            norm_000 = np.linalg.norm(projected_2d_points_000[i1] - projected_2d_points_000[i2])
            norm_mean = np.linalg.norm(projected_2d_points_mean[i1] - projected_2d_points_mean[i2])
            norm_xyz_projected = np.linalg.norm(projected_points[i1] - projected_points[i2])
            assert np.isclose(norm_000, norm_mean)
            assert np.isclose(norm_000, norm_xyz_projected)


def test_plane_zero_normal_vector_raises():
    with pytest.raises(ValueError, match="Normal vector is equal to 0.0"):
        Plane.from_coefficients(0, 0, 0, 1)


def test_plane_properties_and_str():
    plane = Plane.from_coefficients(0, 0, 2, -2)  # z = 1
    assert plane.abcd == approx((0, 0, 1, -1))
    assert (plane.a, plane.b, plane.c, plane.d) == approx((0, 0, 1, -1))
    assert plane.distance_to_origin == approx(-1)
    assert not plane.crosses_origin
    assert Plane.from_coefficients(1, 1, 0, 0).crosses_origin
    assert str(plane) == (
        "Plane object\n  => Normal vector : [0. 0. 1.]\n  => Equation of the plane ax + by + cz + d = 0\n"
        "     with a = 0.0\n          b = 0.0\n          c = 1.0\n          d = -1.0"
    )


def test_plane_distances():
    plane = Plane.from_coefficients(0, 0, 1, -1)  # z = 1
    points = [np.array(pt, float) for pt in [(0, 0, 3), (1, 2, -1), (5, 5, 1), (3, -2, 1.5)]]
    assert plane.distance_to_point(points[1]) == approx(2)
    assert plane.distances(points) == approx([2, -2, 0, 0.5])
    distances, indices = plane.distances_indices_sorted(points, sign=True)
    assert distances == approx([2, -2, 0, 0.5])
    assert indices == [(2, 0), (3, 1), (0, 1), (1, -1)]

    points = [np.array([0, 0, z]) for z in (1, 1.01, 2, 2.05, -1)]
    for delta in (0.1, None):  # None -> 0.05 * max |distance| = 0.1
        _, indices, groups = plane.distances_indices_groups(points, delta=delta)
        assert indices == [0, 1, 2, 3, 4]
        assert groups == [[0, 1], [2, 3], [4]]
    _, _, groups = plane.distances_indices_groups(points, delta=0.1, sign=True)
    assert groups == [[(0, 0), (1, 1)], [(2, 1), (3, 1)], [(4, -1)]]


def test_plane_fit_error_unknown():
    plane = Plane.from_coefficients(0, 0, 1, 0)
    assert plane.fit_error([np.array([0, 0, 1.0])], fit="unknown") is None
    assert plane.fit_error([np.array([0, 0, 1.0]), np.array([0, 0, -2.0])]) == approx(5)
    assert plane.fit_error([np.array([0, 0, 1.0]), np.array([0, 0, -2.0])], fit="maximum_distance") == approx(2)


def test_plane_constructors():
    plane = Plane.from_npoints([np.array([1.0, 0, 0]), np.array([0, 1.0, 0])])
    assert plane.coefficients == approx([0, 0, 1, 0])
    assert plane.crosses_origin
    assert plane.is_same_plane_as(Plane.from_2points_and_origin(np.array([2.0, 1, 0]), np.array([-1.0, 3, 0])))

    bisector = Plane.perpendicular_bisector(np.array([0.0, 0, 0]), np.array([2.0, 0, 0]))
    assert bisector.coefficients == approx([1, 0, 0, -1])

    four_points = [np.array(pt, float) for pt in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0.1)]]
    with pytest.raises(ValueError, match="Cannot initialize Plane"):
        Plane.from_npoints(four_points, best_fit="unknown")


def test_plane_project_and_to2dim():
    plane = Plane.from_coefficients(0, 0, 1, 0)  # z = 0, in-plane basis is (x, y)
    points = [np.array(pt, float) for pt in [(1, 0, 5), (0, 1, -3), (-1, 0, 0), (0, -1, 2)]]
    expected = np.array([(1, 0), (0, 1), (-1, 0), (0, -1)], float)
    assert_allclose(plane.project_and_to2dim(points, None), expected, atol=1e-12)
    assert_allclose(plane.project_and_to2dim(points, "mean"), expected, atol=1e-12)
    assert_allclose(plane.project_and_to2dim(points, np.array([1.0, 1, 7])), expected - 1, atol=1e-12)
    # Angles: 0, pi/2, pi, -pi/2
    assert list(plane.project_and_to2dim_ordered_indices(points)) == [3, 0, 1, 2]


def test_function_comparison():
    def f_zero(x):
        return np.zeros_like(x)

    def f_relu(x):
        return np.maximum(x, 0.0)

    assert function_comparison(f_zero, lambda x: x + 3, 0, 1) == "<"
    assert function_comparison(lambda x: x + 3, f_zero, 0, 1) == ">"
    assert function_comparison(f_zero, f_zero, 0, 1) == "="
    assert function_comparison(f_zero, f_relu, -1, 1) == "<="
    assert function_comparison(f_relu, f_zero, -1, 1) == ">="
    with pytest.raises(RuntimeError, match="Error in comparing functions"):
        function_comparison(f_zero, lambda x: x, -1, 1)


@pytest.mark.parametrize(("xx", "yy"), [((1, 1), (3, 3)), ((3, 3), (1, 1)), ((1, 3), (3, 1)), ((3, 1), (1, 3))])
def test_quarter_ellipsis_functions(xx, yy):
    funcs = quarter_ellipsis_functions(xx, yy)
    p_left, p_right = sorted([xx, yy])
    for func in (funcs["lower"], funcs["upper"]):
        assert func(p_left[0]) == approx(p_left[1])
        assert func(p_right[0]) == approx(p_right[1])
    # Both quarter ellipses have semi-axes of length 2 around the corners of the [1, 3] x [1, 3] square
    assert funcs["lower"](2.0) == approx(3 - np.sqrt(3))
    assert funcs["upper"](2.0) == approx(1 + np.sqrt(3))


def test_quarter_ellipsis_and_diamond_invalid():
    with pytest.raises(RuntimeError, match="Invalid points for quarter_ellipsis_functions"):
        quarter_ellipsis_functions((1, 2), (1, 3))
    with pytest.raises(RuntimeError, match="Invalid points for diamond_functions"):
        diamond_functions((1, 2), (3, 2), y_x0=0.5, x_y0=0.5)


@pytest.mark.parametrize(
    ("xx", "yy", "lower_expected", "upper_expected"),
    [
        # Increasing diamond from (0, 0) to (2, 2): lower = x-b-q, upper = a-p-y
        ((0, 0), (2, 2), [0, 0, 0.5, 1], [1.5, 2, 2, 2]),
        ((2, 2), (0, 0), [0, 0, 0.5, 1], [1.5, 2, 2, 2]),
        # Decreasing diamond from (0, 2) to (2, 0): lower = a-p-y, upper = x-b-q
        ((0, 2), (2, 0), [0.5, 0, 0, 0], [2, 2, 1.5, 1]),
        ((2, 0), (0, 2), [0.5, 0, 0, 0], [2, 2, 1.5, 1]),
    ],
)
def test_diamond_functions(xx, yy, lower_expected, upper_expected):
    funcs = diamond_functions(xx, yy, y_x0=1.0, x_y0=1.0)
    x_vals = np.array([0.5, 1.0, 1.5, 2.0])
    assert_allclose(funcs["lower"](x_vals), lower_expected)
    assert_allclose(funcs["upper"](x_vals), upper_expected)


def test_get_lower_and_upper_f():
    options = {
        "distance_bounds": {"lower": 1.0, "upper": 2.0, "delta": 0.1},
        "angle_bounds": {"lower": 0.2, "upper": 0.8, "delta": 0.1},
    }
    elliptic = get_lower_and_upper_f({**options, "type": "standard_elliptic"})
    reference = quarter_ellipsis_functions((1.0, 0.8), (2.0, 0.2))
    for key in ("lower", "upper"):
        assert elliptic[key](1.5) == approx(reference[key](1.5))
        assert elliptic[key](1.0) == approx(0.8)
        assert elliptic[key](2.0) == approx(0.2)

    diamond = get_lower_and_upper_f({**options, "type": "standard_diamond"})
    reference = diamond_functions((1.0, 0.8), (2.0, 0.2), y_x0=0.1, x_y0=0.1)
    for key in ("lower", "upper"):
        assert diamond[key](1.5) == approx(reference[key](1.5))

    spline = get_lower_and_upper_f(
        {
            **options,
            "type": "standard_spline",
            "lower_points": [(1.0, 0.0), (2.0, 1.0), (3.0, 2.0)],
            "upper_points": [(1.0, 1.0), (2.0, 3.0), (3.0, 5.0)],
            "degree": 1,
        }
    )
    assert spline["lower"](1.5) == approx(0.5)
    assert spline["upper"](2.5) == approx(4.0)

    with pytest.raises(ValueError, match='Surface calculation of type "unknown" is not implemented'):
        get_lower_and_upper_f({**options, "type": "unknown"})


def test_spline_functions_cubic():
    points = [(x, x**3) for x in np.linspace(0, 2, 5)]
    funcs = spline_functions(points, [(x, y + 1) for x, y in points], degree=3)
    assert funcs["lower"](1.25) == approx(1.25**3)
    assert funcs["upper"](1.25) == approx(1.25**3 + 1)


def test_rectangle_surface_intersection():
    def f_lower(x):
        return np.zeros_like(x)

    def f_upper(x):
        return np.asarray(x, float)

    bounds = (0.0, 2.0)
    # Integral of min(x, 1) between 0 and 2
    area, _err = rectangle_surface_intersection(
        ((0.0, 2.0), (0.0, 1.0)), f_lower, f_upper, bounds_lower=bounds, bounds_upper=bounds
    )
    assert area == approx(1.5)
    # Rectangle partially outside of the bounds of the functions
    area, _err = rectangle_surface_intersection(
        ((1.0, 5.0), (0.5, 3.0)), f_lower, f_upper, bounds_lower=bounds, bounds_upper=bounds
    )
    # Integral of (x - 0.5) between 1 and 2
    assert area == approx(1.0)
    # Rectangle outside of the bounds of the functions
    assert rectangle_surface_intersection(
        ((3.0, 4.0), (0.0, 1.0)), f_lower, f_upper, bounds_lower=bounds, bounds_upper=bounds
    ) == (0.0, 0.0)
    # Without the check
    area, _err = rectangle_surface_intersection(((0.0, 1.0), (0.0, 1.0)), f_lower, f_upper, bounds, check=False)
    assert area == approx(0.5)


def test_rectangle_surface_intersection_errors():
    def f_lower(x):
        return np.zeros_like(x)

    def f_upper(x):
        return np.asarray(x, float)

    rectangle = ((0.0, 1.0), (0.0, 1.0))
    with pytest.raises(ValueError, match="Bounds should be identical"):
        rectangle_surface_intersection(rectangle, f_lower, f_upper, bounds_lower=(0, 1), bounds_upper=(0, 2))
    with pytest.raises(ValueError, match="Bounds are given for f_lower but not for f_upper"):
        rectangle_surface_intersection(rectangle, f_lower, f_upper, bounds_lower=(0, 1))
    with pytest.raises(ValueError, match="Bounds are given for f_upper but not for f_lower"):
        rectangle_surface_intersection(rectangle, f_lower, f_upper, bounds_upper=(0, 1))
    with pytest.raises(RuntimeError, match="within the domain defined by the functions bounds"):
        rectangle_surface_intersection(rectangle, f_upper, f_lower, bounds_lower=(0.5, 1), bounds_upper=(0.5, 1))
    with pytest.raises(RuntimeError, match="within the domain defined by x1 and x2"):
        rectangle_surface_intersection(((0.5, 1.0), (0.0, 1.0)), f_upper, f_lower)
    with pytest.raises(NotImplementedError, match="Bounds should be given right now"):
        rectangle_surface_intersection(rectangle, f_lower, f_upper)


@pytest.mark.parametrize("center", [(0.0, 0.0, 0.0), (1.0, -2.0, 3.0)])
def test_solid_angle(center):
    center = np.array(center)
    # One octant of the sphere: 4 pi / 8
    octant = [center + np.array(v, float) for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    assert solid_angle(center, octant) == approx(np.pi / 2)
    # Square face of a cube seen from its center: 4 pi / 6
    face = [center + np.array(v, float) for v in [(1, 1, 1), (1, -1, 1), (-1, -1, 1), (-1, 1, 1)]]
    assert solid_angle(center, face) == approx(2 * np.pi / 3)


def test_matrix_vector_helpers():
    aa = np.array([1.0, 2, 3])
    bb = np.array([-1.0, 0.5, 4])
    assert_allclose(vectorsToMatrix(aa, bb), np.outer(aa, bb))
    mat = np.arange(9, dtype=float).reshape(3, 3)
    assert_allclose(matrixTimesVector(mat, aa), mat @ aa)

    rot_z_90 = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
    coords = [np.array([1.0, 0, 0]), np.array([0.0, 1, 2])]
    expected = [[0, 1, 0], [-1, 0, 2]]
    assert_allclose(rotateCoords(coords, rot_z_90), expected, atol=1e-12)
    assert_allclose(rotateCoordsOpt(coords, rot_z_90), expected, atol=1e-12)

    # Basis (y, z, x): point (1, 2, 3) has coordinates (2, 3, 1)
    new_coords = changebasis(np.array([0, 1, 0]), np.array([0, 0, 1]), np.array([1, 0, 0]), [np.array([1.0, 2, 3])])
    assert_allclose(new_coords, [[2, 3, 1]])


@pytest.mark.parametrize(
    ("p1", "p2", "p3", "expected"),
    [
        ((1, 0, 0), (2, 0, 0), None, True),
        ((1, 0, 0), (0, 1, 0), None, False),
        ((1, 0, 0), (2, 0.01, 0), None, True),
        ((1, 0, 0), (2, 0, 0), (3, 0, 0), True),
        ((1, 0, 0), (0, 1, 0), (0, 0, 0), False),
        ((1, 1, 1), (2, 2, 2), (3, 3, 3.01), True),
    ],
)
def test_collinear(p1, p2, p3, expected):
    p3 = None if p3 is None else np.array(p3, float)
    assert collinear(np.array(p1, float), np.array(p2, float), p3) == expected


def test_anticlockwise_sort():
    pps = [np.array(pt, float) for pt in [(0, 1), (1, 0), (-1, 0), (0, -1)]]
    assert list(anticlockwise_sort_indices(pps)) == [3, 1, 0, 2]
    assert_allclose(anticlockwise_sort(pps), [(0, -1), (1, 0), (0, 1), (-1, 0)])


def test_sort_separation():
    assert sort_separation([[3, 1], [5], [0]]) == [[0], [5], [1, 3]]
    assert sort_separation([[2], [4, 3], [1, 0]]) == [[2], [3, 4], [0, 1]]
    assert sort_separation_tuple([[3, 1], [5], [0]]) == ((0,), (5,), (1, 3))
    assert sort_separation_tuple([[2], [4, 3], [1, 0]]) == ((2,), (3, 4), (0, 1))

    assert separation_in_list([[0], [2, 1], [3]], [[[4], [1, 2], [5]]])
    assert not separation_in_list([[0], [2, 1], [3]], [[[0], [1, 3], [2]]])
    assert not separation_in_list([[0], [2, 1], [3]], [[[0], [1], [2, 3]]])
    assert not separation_in_list([[0], [2, 1], [3]], [])


@pytest.mark.parametrize(
    ("valences", "ii", "jj", "expected"),
    [
        ("undefined", 0, 1, True),
        ([2, -2], 0, 1, True),
        ([-2, 2], 0, 1, True),
        ([2, 1], 0, 1, False),
        ([-1, -2], 0, 1, False),
        ([0, 2], 0, 1, True),
        ([2, 0], 0, 1, True),
    ],
)
def test_is_anion_cation_bond(valences, ii, jj, expected):
    assert is_anion_cation_bond(valences, ii, jj) is expected
