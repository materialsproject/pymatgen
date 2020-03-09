# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest

from pymatgen.core import Lattice
from pymatgen.electronic_structure.core import Orbital, Spin, Magmom
import numpy as np


class SpinTest(unittest.TestCase):

    def test_init(self):
        self.assertEqual(int(Spin.up), 1)
        self.assertEqual(int(Spin.down), -1)

    def test_from_int(self):
        self.assertEqual(Spin(1), Spin.up)
        self.assertEqual(Spin(-1), Spin.down)
        self.assertRaises(ValueError, Spin, 0)

    def test_cached(self):
        self.assertEqual(id(Spin(1)), id(Spin.up))


class OrbitalTest(unittest.TestCase):

    def test_init(self):
        for orb in Orbital:
            self.assertEqual(Orbital(orb.value), orb)
        self.assertRaises(ValueError, Orbital, 100)

    def test_cached(self):
        self.assertEqual(id(Orbital(0)), id(Orbital.s))


class MagmomTest(unittest.TestCase):

    def test_init(self):
        # backwards compatibility for scalar-like magmoms
        magmom = Magmom(2.0)
        self.assertEqual(float(magmom), 2.0)
        # backwards compatibility for list-like magmoms
        magmom2 = Magmom([1, 2, 3])
        self.assertEqual(list(magmom2), [1, 2, 3])
        self.assertEqual(magmom2.global_moment.tolist(), [1, 2, 3])
        # non-default saxis, normalized internally
        magmom3 = Magmom([1, 2, 3], saxis=[1, 1, 1])
        self.assertTrue(np.allclose(magmom3.saxis, [np.sqrt(1 / 3.)] * 3))
        # test construction from known global moment and desired, non-default saxis
        magmom4 = Magmom.from_global_moment_and_saxis([1, 2, 3], saxis=[1, 0, 0])
        self.assertTrue(np.allclose(magmom4.moment, [-3, 2, 1]))
        # test global moments with non-default saxis
        magmom5 = Magmom([-3, 2, 1], saxis=[1, 0, 0])
        self.assertTrue(np.allclose(magmom5.global_moment, [1, 2, 3]))

    def test_get_moments(self):

        # simple cases
        magmom_along_x = Magmom([1, 0, 0])
        self.assertTrue(np.allclose(magmom_along_x.get_moment(saxis=[1, 0, 0]), [0, 0, 1]))

        magmom_along_y = Magmom([0, 1, 0])
        self.assertTrue(np.allclose(magmom_along_y.get_moment(saxis=[0, 1, 0]), [0, 0, 1]))

        # test transformations
        magmoms = [[0, 0, 0],
                   [0, 0, 1],
                   [0, 0, -1],
                   [1, 2, 3],
                   [-1, 2, 3],
                   [-1, -2, -3]]

        for magmom in magmoms:
            magmom1 = Magmom(magmom)
            # transform to non-default saxis
            magmom2 = magmom1.get_00t_magmom_with_xyz_saxis()
            # and back to default saxis
            magmom3 = magmom2.get_xyz_magmom_with_001_saxis()
            self.assertTrue(np.allclose(magmom1.moment, magmom))
            self.assertTrue(np.allclose(magmom1.saxis, [0, 0, 1]))
            self.assertTrue(np.allclose(magmom1.get_moment(saxis=magmom1.saxis), magmom1.moment))
            self.assertTrue(np.allclose(magmom1.get_moment(saxis=magmom2.saxis), magmom2.moment))
            self.assertTrue(np.allclose(magmom2.get_moment(saxis=[0, 0, 1]), magmom1.moment))
            self.assertTrue(np.allclose(magmom2.get_moment(saxis=magmom2.saxis), magmom2.moment))
            self.assertTrue(np.allclose(magmom3.moment, magmom1.moment))

    def test_is_collinear(self):
        magmoms_list = [[0, 0, 0],
                        [1, 1, 1],
                        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                        [[0, 0, 1], [0, 0, 1], [0, 0, 1]],
                        [[0, 0, -1], [0, 0, 1], [0, 0, 1]],
                        [[2, 2, 2], [-2, -2, -2], [2, 2, 2]]]
        for magmoms in magmoms_list:
            self.assertEqual(Magmom.are_collinear(magmoms), True)
        ncl_magmoms = [[[0, 0, 1], [0, 0, 1], [1, 2, 3]]]
        self.assertEqual(Magmom.are_collinear(ncl_magmoms), False)

    def test_have_consistent_saxis(self):
        magmom1 = Magmom([1, 2, 3])
        magmom2 = Magmom([1, 2, 3])
        magmom3 = Magmom([1, 2, 3], saxis=[0, 0, -1])
        magmom4 = Magmom([1, 2, 3], saxis=[1, 2, 3])
        self.assertTrue(Magmom.have_consistent_saxis([magmom1, magmom2]))
        self.assertFalse(Magmom.have_consistent_saxis([magmom1, magmom3]))
        self.assertFalse(Magmom.have_consistent_saxis([magmom1, magmom4]))

    def test_get_consistent_set_and_saxis(self):
        magmoms = [1, 1, 2, 2, 0, 0, 2]
        magmoms, saxis = Magmom.get_consistent_set_and_saxis(magmoms)
        self.assertTrue(np.allclose(saxis, [0, 0, 1]))

        magmoms = [[0, 0, 0],
                   [1, 1, 1],
                   [2, 2, 2]]
        magmoms, saxis = Magmom.get_consistent_set_and_saxis(magmoms)
        self.assertTrue(np.allclose(saxis, [np.sqrt(1 / 3.)] * 3))

    def test_relative_to_crystal_axes(self):
        lattice = Lattice.from_parameters(5, 10, 5, 90, 110, 90)
        moment = [1, 0, 2]
        magmom = Magmom.from_moment_relative_to_crystal_axes(moment, lattice)
        self.assertTrue(np.allclose(magmom.moment, [0.93969262, 0.0, 1.65797986]))
        self.assertTrue(np.allclose(magmom.get_moment_relative_to_crystal_axes(lattice), moment))

    def test_equality(self):
        self.assertTrue(Magmom([1, 1, 1]) == Magmom([1, 1, 1]))
        self.assertFalse(Magmom([1, 1, 2]) == Magmom([1, 1, 1]))
        self.assertTrue(Magmom([0, 0, 10]) == 10)

    def test_negative(self):
        self.assertEqual(-Magmom([1, 2, 3]), Magmom([-1, -2, -3]))


if __name__ == '__main__':
    unittest.main()
