# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np
import unittest
import os

from pymatgen.analysis.local_env import VoronoiNN, JMolNN, \
        MinimumDistanceNN, MinimumOKeeffeNN, MinimumVIRENN
from pymatgen import Element, Structure, Lattice
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class VoronoiNNTest(PymatgenTest):
    def setUp(self):
        self.s = self.get_structure('LiFePO4')
        self.nn = VoronoiNN(targets=[Element("O")])

    def test_get_voronoi_polyhedra(self):
        self.assertEqual(len(self.nn.get_voronoi_polyhedra(self.s, 0).items()), 8)

    def test_get_cn(self):
        self.assertAlmostEqual(self.nn.get_cn(
                self.s, 0, use_weights=True), 5.809265748999465, 7)

    def test_get_coordinated_sites(self):
        self.assertEqual(len(self.nn.get_nn(self.s, 0)), 8)


class JMolNNTest(PymatgenTest):

    def test_get_nn(self):
        s = self.get_structure('LiFePO4')

        # Test the default near-neighbor finder.
        nn = JMolNN()
        nsites_checked = 0

        for site_idx, site in enumerate(s):
            if site.specie == Element("Li"):
                self.assertEqual(nn.get_cn(s, site_idx), 0)
                nsites_checked += 1
            elif site.specie == Element("Fe"):
                self.assertEqual(nn.get_cn(s, site_idx), 6)
                nsites_checked += 1
            elif site.specie == Element("P"):
                self.assertEqual(nn.get_cn(s, site_idx), 4)
                nsites_checked += 1
        self.assertEqual(nsites_checked, 12)

        # Test a user override that would cause Li to show up as 6-coordinated
        nn = JMolNN(el_radius_updates={"Li": 1})
        self.assertEqual(nn.get_cn(s, 0), 6)

        # Verify get_nn function works
        self.assertEqual(len(nn.get_nn(s, 0)), 6)


class MiniDistNNTest(PymatgenTest):

    def test_all(self):
        diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264],
                     [0, 0, 2.528]]), ["C0+", "C0+"], [[2.554, 1.806, 4.423],
                                                       [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012],
                     [0, 0, 4.025]]), ["Na1+", "Cl1-"], [[0, 0, 0],
                                                         [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"], [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.assertAlmostEqual(MinimumDistanceNN().get_cn(
            diamond, 0), 4)
        self.assertAlmostEqual(MinimumDistanceNN().get_cn(
            nacl, 0), 6)
        self.assertAlmostEqual(MinimumDistanceNN(tol=0.01).get_cn(
            cscl, 0), 8)

        self.assertAlmostEqual(MinimumOKeeffeNN(tol=0.01).get_cn(
            diamond, 0), 4)
        self.assertAlmostEqual(MinimumOKeeffeNN(tol=0.01).get_cn(
            nacl, 0), 6)
        self.assertAlmostEqual(MinimumOKeeffeNN(tol=0.01).get_cn(
            cscl, 0), 8)

        self.assertAlmostEqual(MinimumVIRENN(tol=0.01).get_cn(
            diamond, 0), 4)
        self.assertAlmostEqual(MinimumVIRENN(tol=0.01).get_cn(
            nacl, 0), 6)
        self.assertAlmostEqual(MinimumVIRENN(tol=0.01).get_cn(
            cscl, 0), 8)


if __name__ == '__main__':
    unittest.main()
