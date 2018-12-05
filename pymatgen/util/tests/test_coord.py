# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Created on Apr 25, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 25, 2012"

import random
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord import *
from pymatgen.util.testing import PymatgenTest


class CoordUtilsTest(PymatgenTest):

    def test_get_linear_interpolated_value(self):
        xvals = [0, 1, 2, 3, 4, 5]
        yvals = [3, 6, 7, 8, 10, 12]
        self.assertEqual(get_linear_interpolated_value(xvals, yvals, 3.6), 9.2)
        self.assertRaises(ValueError, get_linear_interpolated_value, xvals,
                          yvals, 6)

    def test_in_coord_list(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(in_coord_list(coords, test_coord))
        self.assertTrue(in_coord_list(coords, test_coord, atol=0.15))
        self.assertFalse(in_coord_list([0.99, 0.99, 0.99], test_coord,
                                       atol=0.15))

    def test_is_coord_subset(self):
        c1 = [0,0,0]
        c2 = [0,1.2,-1]
        c3 = [3,2,1]
        c4 = [3-9e-9, 2-9e-9, 1-9e-9]
        self.assertTrue(is_coord_subset([c1, c2, c3], [c1, c4, c2]))
        self.assertTrue(is_coord_subset([c1], [c2, c1]))
        self.assertTrue(is_coord_subset([c1, c2], [c2, c1]))
        self.assertFalse(is_coord_subset([c1, c2], [c2, c3]))
        self.assertFalse(is_coord_subset([c1, c2], [c2]))

    def test_coord_list_mapping(self):
        c1 = [0,.124,0]
        c2 = [0,1.2,-1]
        c3 = [3,2,1]
        a = np.array([c1, c2])
        b = np.array([c3, c2, c1])
        inds = coord_list_mapping(a, b)
        self.assertTrue(np.allclose(a, b[inds]))
        self.assertRaises(Exception, coord_list_mapping, [c1,c2], [c2,c3])
        self.assertRaises(Exception, coord_list_mapping, [c2], [c2,c2])

    def test_coord_list_mapping_pbc(self):
        c1 = [0.1, 0.2, 0.3]
        c2 = [0.2, 0.3, 0.3]
        c3 = [0.5, 0.3, 0.6]
        c4 = [1.5, -0.7, -1.4]

        a = np.array([c1, c3, c2])
        b = np.array([c4, c2, c1])

        inds = coord_list_mapping_pbc(a, b)
        diff = a - b[inds]
        diff -= np.round(diff)
        self.assertTrue(np.allclose(diff, 0))
        self.assertRaises(Exception, coord_list_mapping_pbc, [c1,c2], [c2,c3])
        self.assertRaises(Exception, coord_list_mapping_pbc, [c2], [c2,c2])

    def test_find_in_coord_list(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(find_in_coord_list(coords, test_coord))
        self.assertEqual(find_in_coord_list(coords, test_coord, atol=0.15)[0],
                         0)
        self.assertFalse(find_in_coord_list([0.99, 0.99, 0.99], test_coord,
                                            atol=0.15))
        coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.1, 0.1, 0.1]]
        self.assertArrayEqual(find_in_coord_list(coords, test_coord,
                                                 atol=0.15), [0, 2])

    def test_all_distances(self):
        coords1 = [[0, 0, 0], [0.5, 0.5, 0.5]]
        coords2 = [[1, 2, -1], [1, 0, 0], [1, 0, 0]]
        result = [[2.44948974, 1, 1], [2.17944947, 0.8660254, 0.8660254]]
        self.assertArrayAlmostEqual(all_distances(coords1, coords2), result, 4)

    def test_pbc_diff(self):
        self.assertArrayAlmostEqual(pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]),
                                    [-0.2, -0.4, 0.2])
        self.assertArrayAlmostEqual(pbc_diff([0.9, 0.1, 1.01],
                                             [0.3, 0.5, 0.9]),
                                    [-0.4, -0.4, 0.11])
        self.assertArrayAlmostEqual(pbc_diff([0.1, 0.6, 1.01],
                                             [0.6, 0.1, 0.9]),
                                    [-0.5, 0.5, 0.11])
        self.assertArrayAlmostEqual(pbc_diff([100.1, 0.2, 0.3],
                                             [0123123.4, 0.5, 502312.6]),
                                    [-0.3, -0.3, -0.3])

    def test_in_coord_list_pbc(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(in_coord_list_pbc(coords, test_coord))
        self.assertTrue(in_coord_list_pbc(coords, test_coord, atol=0.15))
        test_coord = [0.99, 0.99, 0.99]
        self.assertFalse(in_coord_list_pbc(coords, test_coord, atol=0.01))

    def test_find_in_coord_list_pbc(self):
        coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
        test_coord = [0.1, 0.1, 0.1]
        self.assertFalse(find_in_coord_list_pbc(coords, test_coord))
        self.assertEqual(find_in_coord_list_pbc(coords, test_coord,
                                                atol=0.15)[0], 0)
        test_coord = [0.99, 0.99, 0.99]
        self.assertEqual(
            find_in_coord_list_pbc(coords, test_coord, atol=0.02)[0], 0)
        test_coord = [-0.499, -0.499, -0.499]
        self.assertEqual(
            find_in_coord_list_pbc(coords, test_coord, atol=0.01)[0], 1)

    def test_is_coord_subset_pbc(self):
        c1 = [0, 0, 0]
        c2 = [0, 1.2, -1]
        c3 = [2.3, 0, 1]
        c4 = [1.3-9e-9, -1-9e-9, 1-9e-9]
        self.assertTrue(is_coord_subset_pbc([c1, c2, c3], [c1, c4, c2]))
        self.assertTrue(is_coord_subset_pbc([c1], [c2, c1]))
        self.assertTrue(is_coord_subset_pbc([c1, c2], [c2, c1]))
        self.assertFalse(is_coord_subset_pbc([c1, c2], [c2, c3]))
        self.assertFalse(is_coord_subset_pbc([c1, c2], [c2]))

        # test tolerances
        c5 = [0.1, 0.1, 0.2]
        atol1 = [0.25, 0.15, 0.15]
        atol2 = [0.15, 0.15, 0.25]
        self.assertFalse(is_coord_subset_pbc([c1], [c5], atol1))
        self.assertTrue(is_coord_subset_pbc([c1], [c5], atol2))

        # test mask
        mask1 = [[True]]
        self.assertFalse(is_coord_subset_pbc([c1], [c5], atol2, mask1))
        mask2 = [[True, False]]
        self.assertTrue(is_coord_subset_pbc([c1], [c2, c1], mask=mask2))
        self.assertFalse(is_coord_subset_pbc([c1], [c1, c2], mask=mask2))
        mask3 = [[False, True]]
        self.assertFalse(is_coord_subset_pbc([c1], [c2, c1], mask=mask3))
        self.assertTrue(is_coord_subset_pbc([c1], [c1, c2], mask=mask3))


    def test_lattice_points_in_supercell(self):
        supercell = np.array([[1, 3, 5], [-3, 2, 3], [-5, 3, 1]])
        points = lattice_points_in_supercell(supercell)
        self.assertAlmostEqual(len(points), abs(np.linalg.det(supercell)))
        self.assertGreaterEqual(np.min(points), -1e-10)
        self.assertLessEqual(np.max(points), 1-1e-10)

        supercell = np.array([[-5, -5, -3], [0, -4, -2], [0, -5, -2]])
        points = lattice_points_in_supercell(supercell)
        self.assertAlmostEqual(len(points), abs(np.linalg.det(supercell)))
        self.assertGreaterEqual(np.min(points), -1e-10)
        self.assertLessEqual(np.max(points), 1-1e-10)

    def test_barycentric(self):
        #2d test
        simplex1 = np.array([[0.3, 0.1], [0.2, -1.2], [1.3, 2.3]])
        pts1 = np.array([[0.6, 0.1], [1.3, 2.3], [0.5, 0.5], [.7, 1]])
        output1 = barycentric_coords(pts1, simplex1)
        #do back conversion to cartesian
        o_dot_s = np.sum(output1[:, :, None] * simplex1[None, :, :], axis=1)
        self.assertTrue(np.allclose(pts1, o_dot_s))

        #do 3d tests
        simplex2 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, 0]])
        pts2 = np.array([[0, 0, 1], [0, 0.5, 0.5], [1./3, 1./3, 1./3]])
        output2 = barycentric_coords(pts2, simplex2)
        self.assertTrue(np.allclose(output2[1], [0.5, 0.5, 0, 0]))
        #do back conversion to cartesian
        o_dot_s = np.sum(output2[:, :, None] * simplex2[None, :, :], axis=1)
        self.assertTrue(np.allclose(pts2, o_dot_s))
        #test single point
        self.assertTrue(np.allclose(output2[2],
                                    barycentric_coords(pts2[2], simplex2)))

    def test_pbc_shortest_vectors(self):
        fcoords = np.array([[0.3, 0.3, 0.5],
                            [0.1, 0.1, 0.3],
                            [0.9, 0.9, 0.8],
                            [0.1, 0.0, 0.5],
                            [0.9, 0.7, 0.0]])
        lattice = Lattice.from_lengths_and_angles([8, 8, 4],
                                                  [90, 76, 58])
        expected = np.array([[0.000, 3.015, 4.072, 3.519, 3.245],
                             [3.015, 0.000, 3.207, 1.131, 4.453],
                             [4.072, 3.207, 0.000, 2.251, 1.788],
                             [3.519, 1.131, 2.251, 0.000, 3.852]])

        vectors = pbc_shortest_vectors(lattice, fcoords[:-1], fcoords)
        dists = np.sum(vectors**2, axis = -1)**0.5
        self.assertArrayAlmostEqual(dists, expected, 3)

        #now try with small loop threshold
        from pymatgen.util import coord
        prev_threshold = coord.LOOP_THRESHOLD
        coord.LOOP_THRESHOLD = 0

        vectors = pbc_shortest_vectors(lattice, fcoords[:-1], fcoords)
        dists = np.sum(vectors**2, axis = -1)**0.5
        self.assertArrayAlmostEqual(dists, expected, 3)

        coord.LOOP_THRESHOLD = prev_threshold

    def test_get_angle(self):
        v1 = (1, 0, 0)
        v2 = (1, 1, 1)
        self.assertAlmostEqual(get_angle(v1, v2), 54.7356103172)
        self.assertAlmostEqual(get_angle(v1, v2, units="radians"),
                               0.9553166181245092)

class SimplexTest(PymatgenTest):

    def setUp(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0, 1, 0])
        coords.append([0, 0, 1])
        coords.append([1, 0, 0])
        self.simplex = Simplex(coords)

    def test_equal(self):
        c2 = list(self.simplex.coords)
        random.shuffle(c2)
        self.assertEqual(Simplex(c2), self.simplex)

    def test_in_simplex(self):
        self.assertTrue(self.simplex.in_simplex([0.1, 0.1, 0.1]))
        self.assertFalse(self.simplex.in_simplex([0.6, 0.6, 0.6]))
        for i in range(10):
            coord = np.random.random_sample(size=3) / 3
            self.assertTrue(self.simplex.in_simplex(coord))

    def test_2dtriangle(self):
        s = Simplex([[0, 1], [1, 1], [1, 0]])
        self.assertArrayAlmostEqual(s.bary_coords([0.5, 0.5]),
                                    [0.5, 0, 0.5])
        self.assertArrayAlmostEqual(s.bary_coords([0.5, 1]), [0.5, 0.5, 0])
        self.assertArrayAlmostEqual(s.bary_coords([0.5, 0.75]), [0.5, 0.25, 0.25])
        self.assertArrayAlmostEqual(s.bary_coords([0.75, 0.75]), [0.25, 0.5, 0.25])

        s = Simplex([[1, 1], [1, 0]])
        self.assertRaises(ValueError, s.bary_coords, [0.5, 0.5])

    def test_volume(self):
        # Should be value of a right tetrahedron.
        self.assertAlmostEqual(self.simplex.volume, 1/6)

    def test_str(self):
        self.assertTrue(str(self.simplex).startswith("3-simplex in 4D space"))
        self.assertTrue(repr(self.simplex).startswith("3-simplex in 4D space"))

    def test_bary_coords(self):
        s = Simplex([[0, 2], [3, 1], [1, 0]])
        point = [0.7, 0.5]
        bc = s.bary_coords(point)
        self.assertArrayAlmostEqual(bc, [0.26, -0.02, 0.76])
        new_point = s.point_from_bary_coords(bc)
        self.assertArrayAlmostEqual(point, new_point)

    def test_intersection(self):
        # simple test, with 2 intersections at faces
        s = Simplex([[0, 2], [3, 1], [1, 0]])
        point1 = [0.7, 0.5]
        point2 = [0.5, 0.7]
        intersections = s.line_intersection(point1, point2)
        expected = np.array([[1.13333333, 0.06666667],
                             [ 0.8,  0.4]])
        self.assertArrayAlmostEqual(intersections, expected)

        # intersection through point and face
        point1 = [0, 2]  # simplex point
        point2 = [1, 1]  # inside simplex
        expected = np.array([[1.66666667, 0.33333333],
                             [0, 2]])
        intersections = s.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)

        # intersection through point only
        point1 = [0, 2]  # simplex point
        point2 = [0.5, 0.7]
        expected = np.array([[0, 2]])
        intersections = s.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)

        # 3d intersection through edge and face
        point1 = [0.5, 0, 0] # edge point
        point2 = [0.5, 0.5, 0.5] # in simplex
        expected = np.array([[ 0.5, 0.25, 0.25],
                             [ 0.5, 0. , 0. ]])
        intersections = self.simplex.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)

        # 3d intersection through edge only
        point1 = [0.5, 0, 0] # edge point
        point2 = [0.5, 0.5, -0.5] # outside simplex
        expected = np.array([[0.5, 0., 0.]])
        intersections = self.simplex.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)

        # coplanar to face (no intersection)
        point1 = [-1, 2]
        point2 = [0, 0]
        expected = np.array([])
        intersections = s.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)

        # coplanar to face (with intersection line)
        point1 = [0, 2]  # simplex point
        point2 = [1, 0]
        expected = np.array([[1, 0],
                             [0, 2]])
        intersections = s.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)

        # coplanar to face (with intersection points)
        point1 = [0.1, 2]
        point2 = [1.1, 0]
        expected = np.array([[1.08, 0.04],
                             [0.12, 1.96]])
        intersections = s.line_intersection(point1, point2)
        self.assertArrayAlmostEqual(intersections, expected)


if __name__ == "__main__":
    import unittest
    unittest.main()
