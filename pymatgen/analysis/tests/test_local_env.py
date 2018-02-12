# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np
import unittest
import os

from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator, \
        VoronoiNN, JMolNN, \
        MinimumDistanceNN, MinimumOKeeffeNN, MinimumVIRENN, \
        get_neighbors_of_site_with_index, site_is_of_motif_type, \
        NearNeighbors, LocalStructOrderParas
from pymatgen import Element, Structure, Lattice
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.cif import CifParser

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class ValenceIonicRadiusEvaluatorTest(PymatgenTest):
    def setUp(self):
        """
        Setup MgO rocksalt structure for testing Vacancy
        """
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
                         [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, True,
                                 True)
        self._mgo_valrad_evaluator = ValenceIonicRadiusEvaluator(self._mgo_uc)

    def test_valences_ionic_structure(self):
        valence_dict = self._mgo_valrad_evaluator.valences
        for val in list(valence_dict.values()):
            self.assertTrue(val in {2, -2})

    def test_radii_ionic_structure(self):
        radii_dict = self._mgo_valrad_evaluator.radii
        for rad in list(radii_dict.values()):
            self.assertTrue(rad in {0.86, 1.26})

    def tearDown(self):
        del self._mgo_uc
        del self._mgo_valrad_evaluator

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

    def tearDown(self):
        del self.s
        del self.nn


class JMolNNTest(PymatgenTest):

    def setUp(self):
        self.jmol = JMolNN()
        self.jmol_update = JMolNN(el_radius_updates={"Li": 1})

    def test_get_nn(self):
        s = self.get_structure('LiFePO4')

        # Test the default near-neighbor finder.
        nsites_checked = 0

        for site_idx, site in enumerate(s):
            if site.specie == Element("Li"):
                self.assertEqual(self.jmol.get_cn(s, site_idx), 0)
                nsites_checked += 1
            elif site.specie == Element("Fe"):
                self.assertEqual(self.jmol.get_cn(s, site_idx), 6)
                nsites_checked += 1
            elif site.specie == Element("P"):
                self.assertEqual(self.jmol.get_cn(s, site_idx), 4)
                nsites_checked += 1
        self.assertEqual(nsites_checked, 12)

        # Test a user override that would cause Li to show up as 6-coordinated
        self.assertEqual(self.jmol_update.get_cn(s, 0), 6)

        # Verify get_nn function works
        self.assertEqual(len(self.jmol_update.get_nn(s, 0)), 6)

    def tearDown(self):
        del self.jmol
        del self.jmol_update


class MiniDistNNTest(PymatgenTest):

    def setUp(self):
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264],
                     [0, 0, 2.528]]), ["C0+", "C0+"], [[2.554, 1.806, 4.423],
                                                       [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012],
                     [0, 0, 4.025]]), ["Na1+", "Cl1-"], [[0, 0, 0],
                                                         [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"], [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.mos2 = Structure(
            Lattice([[3.19, 0, 0], [-1.595, 2.763, 0], [0, 0, 17.44]]),
            ['Mo', 'S', 'S'], [[-1e-06, 1.842, 3.72], [1.595, 0.92, 5.29], \
            [1.595, 0.92, 2.155]], coords_are_cartesian=True)

    def test_all_nn_classes(self):
        self.assertAlmostEqual(MinimumDistanceNN().get_cn(
            self.diamond, 0), 4)
        self.assertAlmostEqual(MinimumDistanceNN().get_cn(
            self.nacl, 0), 6)
        self.assertAlmostEqual(MinimumDistanceNN(tol=0.01).get_cn(
            self.cscl, 0), 8)
        self.assertAlmostEqual(MinimumDistanceNN(tol=0.1).get_cn(
            self.mos2, 0), 6)
        for image in MinimumDistanceNN(tol=0.1).get_nn_images(self.mos2, 0):
            self.assertTrue(image in [[0, 0, 0], [0, 1, 0], [-1, 0, 0], \
                    [0, 0, 0], [0, 1, 0], [-1, 0, 0]])

        self.assertAlmostEqual(MinimumOKeeffeNN(tol=0.01).get_cn(
            self.diamond, 0), 4)
        self.assertAlmostEqual(MinimumOKeeffeNN(tol=0.01).get_cn(
            self.nacl, 0), 6)
        self.assertAlmostEqual(MinimumOKeeffeNN(tol=0.01).get_cn(
            self.cscl, 0), 8)

        self.assertAlmostEqual(MinimumVIRENN(tol=0.01).get_cn(
            self.diamond, 0), 4)
        self.assertAlmostEqual(MinimumVIRENN(tol=0.01).get_cn(
            self.nacl, 0), 6)
        self.assertAlmostEqual(MinimumVIRENN(tol=0.01).get_cn(
            self.cscl, 0), 8)

    def tearDown(self):
        del self.diamond
        del self.nacl
        del self.cscl
        del self.mos2


class MotifIdentificationTest(PymatgenTest):

    def setUp(self):
        self.silicon = Structure(
                Lattice.from_lengths_and_angles(
                        [5.47, 5.47, 5.47],
                        [90.0, 90.0, 90.0]),
                ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
                [[0.000000, 0.000000, 0.500000],
                [0.750000, 0.750000, 0.750000],
                [0.000000, 0.500000, 1.000000],
                [0.750000, 0.250000, 0.250000],
                [0.500000, 0.000000, 1.000000],
                [0.250000, 0.750000, 0.250000],
                [0.500000, 0.500000, 0.500000],
                [0.250000, 0.250000, 0.750000]],
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=False, site_properties=None)
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264],
                     [0, 0, 2.528]]), ["C0+", "C0+"], [[2.554, 1.806, 4.423],
                                                       [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012],
                     [0, 0, 4.025]]), ["Na1+", "Cl1-"], [[0, 0, 0],
                                                         [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"], [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.square_pyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["C", "C", "C", "C", "C", "C"], [
            [0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], \
            [0, 0, 1]], validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.trigonal_bipyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["P", "Cl", "Cl", "Cl", "Cl", "Cl"], [
            [0, 0, 0], [0, 0, 2.14], [0, 2.02, 0], [1.74937, -1.01, 0], \
            [-1.74937, -1.01, 0], [0, 0, -2.14]], validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)

    def test_site_is_of_motif_type(self):
        for i in range(self.diamond.num_sites):
            self.assertEqual(site_is_of_motif_type(
                    self.diamond, i), "tetrahedral")
        for i in range(self.nacl.num_sites):
            self.assertEqual(site_is_of_motif_type(
                    self.nacl, i), "octahedral")
        for i in range(self.cscl.num_sites):
            self.assertEqual(site_is_of_motif_type(
                    self.cscl, i), "bcc")
        self.assertEqual(site_is_of_motif_type(
                self.square_pyramid, 0), "square pyramidal")
        for i in range(1, self.square_pyramid.num_sites):
            self.assertEqual(site_is_of_motif_type(
                    self.square_pyramid, i), "unrecognized")
        self.assertEqual(site_is_of_motif_type(
                self.trigonal_bipyramid, 0), "trigonal bipyramidal")
        for i in range(1, self.trigonal_bipyramid.num_sites):
            self.assertEqual(site_is_of_motif_type(
                    self.trigonal_bipyramid, i), "unrecognized")

    def test_get_neighbors_of_site_with_index(self):
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.diamond, 0)), 4)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.nacl, 0)), 6)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.cscl, 0)), 8)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.diamond, 0, delta=0.01)), 4)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.diamond, 0, cutoff=6)), 4)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.diamond, 0, approach="voronoi")), 4)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.diamond, 0, approach="min_OKeeffe")), 4)
        self.assertEqual(len(get_neighbors_of_site_with_index(
                self.diamond, 0, approach="min_VIRE")), 4)


    def tearDown(self):
        del self.silicon
        del self.diamond
        del self.nacl
        del self.cscl

class NearNeighborTest(PymatgenTest):

    def setUp(self):
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264],
                     [0, 0, 2.528]]), ["C0+", "C0+"], [[2.554, 1.806, 4.423],
                                                       [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)

    def set_nn_info(self):

        # check conformance
        # implicitly assumes that all NearNeighbors subclasses
        # will correctly identify bonds in diamond, if it
        # can't there are probably bigger problems
        subclasses = NearNeighbors.__subclasses__()
        for subclass in subclasses:
            nn_info = subclass().get_nn_info(self.diamond, 0)
            self.assertEqual(nn_info[0]['site_index'], 1)
            self.assertEqual(nn_info[0]['image'][0], 1)

    def tearDown(self):
        del self.diamond

class LocalStructOrderParasTest(PymatgenTest):
    def setUp(self):
        self.single_bond = Structure(
            Lattice.from_lengths_and_angles(
            [10, 10, 10], [90, 90, 90]),
            ["H", "H", "H"], [[1, 0, 0], [0, 0, 0], [6, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.linear = Structure(
            Lattice.from_lengths_and_angles(
            [10, 10, 10], [90, 90, 90]),
            ["H", "H", "H"], [[1, 0, 0], [0, 0, 0], [2, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.bent45 = Structure(
            Lattice.from_lengths_and_angles(
            [10, 10, 10], [90, 90, 90]), ["H", "H", "H"],
            [[0, 0, 0], [0.707, 0.707, 0], [0.707, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=True,
            site_properties=None)
        self.cubic = Structure(
            Lattice.from_lengths_and_angles(
            [1, 1, 1], [90, 90, 90]),
            ["H"], [[0, 0, 0]], validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=False,
            site_properties=None)
        self.bcc = Structure(
            Lattice.from_lengths_and_angles(
            [1, 1, 1], [90, 90, 90]),
            ["H", "H"], [[0, 0, 0], [0.5, 0.5, 0.5]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.fcc = Structure(
            Lattice.from_lengths_and_angles(
            [1, 1, 1], [90, 90, 90]), ["H", "H", "H", "H"],
            [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.hcp = Structure(
            Lattice.from_lengths_and_angles(
            [1, 1, 1.633], [90, 90, 120]), ["H", "H"],
            [[0.3333, 0.6667, 0.25], [0.6667, 0.3333, 0.75]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.diamond = Structure(
            Lattice.from_lengths_and_angles(
            [1, 1, 1], [90, 90, 90]), ["H", "H", "H", "H", "H", "H", "H", "H"],
            [[0, 0, 0.5], [0.75, 0.75, 0.75], [0, 0.5, 0], [0.75, 0.25, 0.25],
            [0.5, 0, 0], [0.25, 0.75, 0.25], [0.5, 0.5, 0.5],
            [0.25, 0.25, 0.75]], validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.trigonal_off_plane = Structure(
            Lattice.from_lengths_and_angles(
            [100, 100, 100], [90, 90, 90]),
            ["H", "H", "H", "H"],
            [[0.50, 0.50, 0.50], [0.25, 0.75, 0.25], \
            [0.25, 0.25, 0.75], [0.75, 0.25, 0.25]], \
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.regular_triangle = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["H", "H", "H", "H"],
            [[15, 15.28867, 15.65], [14.5, 15, 15], [15.5, 15, 15], \
            [15, 15.866, 15]], validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.trigonal_planar = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["H", "H", "H", "H"],
            [[15, 15.28867, 15], [14.5, 15, 15], [15.5, 15, 15], \
            [15, 15.866, 15]], validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.square_planar = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["H", "H", "H", "H", "H"],
            [[15, 15, 15], [14.75, 14.75, 15], [14.75, 15.25, 15], \
            [15.25, 14.75, 15], [15.25, 15.25, 15]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.square = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["H", "H", "H", "H", "H"],
            [[15, 15, 15.707], [14.75, 14.75, 15], [14.75, 15.25, 15], \
            [15.25, 14.75, 15], [15.25, 15.25, 15]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.T_shape = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["H", "H", "H", "H"],
            [[15, 15, 15], [15, 15, 15.5], [15, 15.5, 15],
            [15, 14.5, 15]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.square_pyramid = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["H", "H", "H", "H", "H", "H"],
            [[15, 15, 15], [15, 15, 15.3535], [14.75, 14.75, 15],
            [14.75, 15.25, 15], [15.25, 14.75, 15], [15.25, 15.25, 15]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.pentagonal_planar = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["Xe", "F", "F", "F", "F", "F"],
            [[0, -1.6237, 0], [1.17969, 0, 0], [-1.17969, 0, 0], \
            [1.90877, -2.24389, 0], [-1.90877, -2.24389, 0], [0, -3.6307, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.pentagonal_pyramid = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["Xe", "F", "F", "F", "F", "F", "F"],
            [[0, -1.6237, 0], [0, -1.6237, 1.17969], [1.17969, 0, 0], \
            [-1.17969, 0, 0], [1.90877, -2.24389, 0], \
            [-1.90877, -2.24389, 0], [0, -3.6307, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.pentagonal_bipyramid = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]),
            ["Xe", "F", "F", "F", "F", "F", "F", "F"],
            [[0, -1.6237, 0], [0, -1.6237, -1.17969], \
            [0, -1.6237, 1.17969], [1.17969, 0, 0], \
            [-1.17969, 0, 0], [1.90877, -2.24389, 0], \
            [-1.90877, -2.24389, 0], [0, -3.6307, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.hexagonal_planar = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]),
            ["H", "C", "C", "C", "C", "C", "C"],
            [[0, 0, 0], [0.71, 1.2298, 0],
            [-0.71, 1.2298, 0], [0.71, -1.2298, 0], [-0.71, -1.2298, 0],
            [1.4199, 0, 0], [-1.4199, 0, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.hexagonal_pyramid = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), \
            ["H", "Li", "C", "C", "C", "C", "C", "C"],
            [[0, 0, 0], [0, 0, 1.675], [0.71, 1.2298, 0], \
            [-0.71, 1.2298, 0], [0.71, -1.2298, 0], [-0.71, -1.2298, 0], \
            [1.4199, 0, 0], [-1.4199, 0, 0]], \
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.hexagonal_bipyramid = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), \
            ["H", "Li", "Li", "C", "C", "C", "C", "C", "C"],
            [[0, 0, 0], [0, 0, 1.675], [0, 0, -1.675], \
            [0.71, 1.2298, 0], [-0.71, 1.2298, 0], \
            [0.71, -1.2298, 0], [-0.71, -1.2298, 0], \
            [1.4199, 0, 0], [-1.4199, 0, 0]], \
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.trigonal_pyramid = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["P", "Cl", "Cl", "Cl", "Cl"],
            [[0, 0, 0], [0, 0, 2.14], [0, 2.02, 0],
            [1.74937, -1.01, 0], [-1.74937, -1.01, 0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.trigonal_bipyramidal = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]), ["P", "Cl", "Cl", "Cl", "Cl", "Cl"],
            [[0, 0, 0], [0, 0, 2.14], [0, 2.02, 0],
            [1.74937, -1.01, 0], [-1.74937, -1.01, 0], [0, 0, -2.14]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.cuboctahedron = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]),
            ["H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
            [[15, 15, 15], [15, 14.5, 14.5], [15, 14.5, 15.5],
            [15, 15.5, 14.5], [15, 15.5, 15.5],
            [14.5, 15, 14.5], [14.5, 15, 15.5], [15.5, 15, 14.5], [15.5, 15, 15.5],
            [14.5, 14.5, 15], [14.5, 15.5, 15], [15.5, 14.5, 15], [15.5, 15.5, 15]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)
        self.see_saw_rect = Structure(
            Lattice.from_lengths_and_angles(
            [30, 30, 30], [90, 90, 90]),
            ["H", "H", "H", "H", "H"],
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0 , 0.0],
            [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=True, site_properties=None)

    def test_init(self):
        self.assertIsNotNone(
            LocalStructOrderParas(["cn"], parameters=None, cutoff=0.99))

    def test_get_order_parameters(self):
        # Set up everything.
        op_types = ["cn", "bent", "bent", "tet", "oct", "bcc", "q2", "q4", \
            "q6", "reg_tri", "sq", "sq_pyr_legacy", "tri_bipyr", "sgl_bd", \
            "tri_plan", "sq_plan", "pent_plan", "sq_pyr", "tri_pyr", \
            "pent_pyr", "hex_pyr", "pent_bipyr", "hex_bipyr", "T", "cuboct", \
            "see_saw_rect", "hex_plan_max", "tet_max", "oct_max", "tri_plan_max", "sq_plan_max", \
            "pent_plan_max", "cuboct_max"]
        op_paras = [None for i in range(len(op_types))]
        op_paras[1] = {'TA': 1, 'IGW_TA': 1./0.0667}
        op_paras[2] = {'TA': 45./180, 'IGW_TA': 1./0.0667}
        ops_044 = LocalStructOrderParas(op_types, parameters=op_paras, cutoff=0.44)
        ops_071 = LocalStructOrderParas(op_types, parameters=op_paras, cutoff=0.71)
        ops_087 = LocalStructOrderParas(op_types, parameters=op_paras, cutoff=0.87)
        ops_099 = LocalStructOrderParas(op_types, parameters=op_paras, cutoff=0.99)
        ops_101 = LocalStructOrderParas(op_types, parameters=op_paras, cutoff=1.01)
        ops_501 = LocalStructOrderParas(op_types, parameters=op_paras, cutoff=5.01)
        ops_voro = LocalStructOrderParas(op_types, parameters=op_paras)

        # Single bond.
        op_vals = ops_101.get_order_parameters(self.single_bond, 0)
        self.assertAlmostEqual(int(op_vals[13] * 1000), 1000)
        op_vals = ops_501.get_order_parameters(self.single_bond, 0)
        self.assertAlmostEqual(int(op_vals[13] * 1000), 799)
        op_vals = ops_101.get_order_parameters(self.linear, 0)
        self.assertAlmostEqual(int(op_vals[13] * 1000), 0)

        # Linear motif.
        op_vals = ops_101.get_order_parameters(self.linear, 0)
        self.assertAlmostEqual(int(op_vals[1] * 1000), 1000)

        # 45 degrees-bent motif.
        op_vals = ops_101.get_order_parameters(self.bent45, 0)
        self.assertAlmostEqual(int(op_vals[2] * 1000), 1000)

        # T-shape motif.
        op_vals = ops_101.get_order_parameters(
            self.T_shape, 0, indices_neighs=[1,2,3])
        self.assertAlmostEqual(int(op_vals[23] * 1000), 1000)

        # Cubic structure.
        op_vals = ops_099.get_order_parameters(self.cubic, 0)
        self.assertAlmostEqual(op_vals[0], 0.0)
        self.assertIsNone(op_vals[3])
        self.assertIsNone(op_vals[4])
        self.assertIsNone(op_vals[5])
        self.assertIsNone(op_vals[6])
        self.assertIsNone(op_vals[7])
        self.assertIsNone(op_vals[8])
        op_vals = ops_101.get_order_parameters(self.cubic, 0)
        self.assertAlmostEqual(op_vals[0], 6.0)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 23)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 1000)
        self.assertAlmostEqual(int(op_vals[5] * 1000), 333)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[7] * 1000), 763)
        self.assertAlmostEqual(int(op_vals[8] * 1000), 353)
        self.assertAlmostEqual(int(op_vals[28] * 1000), 1000)

        # Bcc structure.
        op_vals = ops_087.get_order_parameters(self.bcc, 0)
        self.assertAlmostEqual(op_vals[0], 8.0)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 200)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 145)
        self.assertAlmostEqual(int(op_vals[5] * 1000 + 0.5), 1000)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[7] * 1000), 509)
        self.assertAlmostEqual(int(op_vals[8] * 1000), 628)

        # Fcc structure.
        op_vals = ops_071.get_order_parameters(self.fcc, 0)
        self.assertAlmostEqual(op_vals[0], 12.0)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 36)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 78)
        self.assertAlmostEqual(int(op_vals[5] * 1000), -2)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[7] * 1000), 190)
        self.assertAlmostEqual(int(op_vals[8] * 1000), 574)

        # Hcp structure.
        op_vals = ops_101.get_order_parameters(self.hcp, 0)
        self.assertAlmostEqual(op_vals[0], 12.0)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 33)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 82)
        self.assertAlmostEqual(int(op_vals[5] * 1000), -26)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[7] * 1000), 97)
        self.assertAlmostEqual(int(op_vals[8] * 1000), 484)

        # Diamond structure.
        op_vals = ops_044.get_order_parameters(self.diamond, 0)
        self.assertAlmostEqual(op_vals[0], 4.0)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 1000)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 37)
        self.assertAlmostEqual(op_vals[5], 0.75)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[7] * 1000), 509)
        self.assertAlmostEqual(int(op_vals[8] * 1000), 628)
        self.assertAlmostEqual(int(op_vals[27] * 1000), 1000)

        # Trigonal off-plane molecule.
        op_vals = ops_044.get_order_parameters(self.trigonal_off_plane, 0)
        self.assertAlmostEqual(op_vals[0], 3.0)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 1000)

        # Trigonal-planar motif.
        op_vals = ops_101.get_order_parameters(self.trigonal_planar, 0)
        self.assertEqual(int(op_vals[0] + 0.5), 3)
        self.assertAlmostEqual(int(op_vals[14] * 1000 + 0.5), 1000)
        self.assertAlmostEqual(int(op_vals[29] * 1000 + 0.5), 1000)

        # Regular triangle motif.
        op_vals = ops_101.get_order_parameters(self.regular_triangle, 0)
        self.assertAlmostEqual(int(op_vals[9] * 1000), 999)

        # Square-planar motif.
        op_vals = ops_101.get_order_parameters(self.square_planar, 0)
        self.assertAlmostEqual(int(op_vals[15] * 1000 + 0.5), 1000)
        self.assertAlmostEqual(int(op_vals[30] * 1000 + 0.5), 1000)

        # Square motif.
        op_vals = ops_101.get_order_parameters(self.square, 0)
        self.assertAlmostEqual(int(op_vals[10] * 1000), 1000)

        # Pentagonal planar.
        op_vals = ops_101.get_order_parameters(
                self.pentagonal_planar.sites, 0, indices_neighs=[1,2,3,4,5])
        self.assertAlmostEqual(int(op_vals[12] * 1000 + 0.5), 126)
        self.assertAlmostEqual(int(op_vals[16] * 1000 + 0.5), 1000)
        self.assertAlmostEqual(int(op_vals[31] * 1000 + 0.5), 1000)

        # Trigonal pyramid motif.
        op_vals = ops_101.get_order_parameters(
            self.trigonal_pyramid, 0, indices_neighs=[1,2,3,4])
        self.assertAlmostEqual(int(op_vals[18] * 1000 + 0.5), 1000)

        # Square pyramid motif.
        op_vals = ops_101.get_order_parameters(self.square_pyramid, 0)
        self.assertAlmostEqual(int(op_vals[11] * 1000 + 0.5), 1000)
        self.assertAlmostEqual(int(op_vals[12] * 1000 + 0.5), 667)
        self.assertAlmostEqual(int(op_vals[17] * 1000 + 0.5), 1000)

        # Pentagonal pyramid motif.
        op_vals = ops_101.get_order_parameters(
            self.pentagonal_pyramid, 0, indices_neighs=[1,2,3,4,5,6])
        self.assertAlmostEqual(int(op_vals[19] * 1000 + 0.5), 1000)

        # Hexagonal pyramid motif.
        op_vals = ops_101.get_order_parameters(
            self.hexagonal_pyramid, 0, indices_neighs=[1,2,3,4,5,6,7])
        self.assertAlmostEqual(int(op_vals[20] * 1000 + 0.5), 1000)

        # Trigonal bipyramidal.
        op_vals = ops_101.get_order_parameters(
            self.trigonal_bipyramidal.sites, 0, indices_neighs=[1,2,3,4,5])
        self.assertAlmostEqual(int(op_vals[12] * 1000 + 0.5), 1000)

        # Pentagonal bipyramidal.
        op_vals = ops_101.get_order_parameters(
            self.pentagonal_bipyramid.sites, 0,
            indices_neighs=[1,2,3,4,5,6,7])
        self.assertAlmostEqual(int(op_vals[21] * 1000 + 0.5), 1000)

        # Hexagonal bipyramid motif.
        op_vals = ops_101.get_order_parameters(
            self.hexagonal_bipyramid, 0, indices_neighs=[1,2,3,4,5,6,7,8])
        self.assertAlmostEqual(int(op_vals[22] * 1000 + 0.5), 1000)

        # Cuboctahedral motif.
        op_vals = ops_101.get_order_parameters(
            self.cuboctahedron, 0, indices_neighs=[i for i in range(1, 13)])
        self.assertAlmostEqual(int(op_vals[24] * 1000 + 0.5), 1000)
        self.assertAlmostEqual(int(op_vals[32] * 1000 + 0.5), 1000)

        # See-saw motif.
        op_vals = ops_101.get_order_parameters(
            self.see_saw_rect, 0, indices_neighs=[i for i in range(1, 5)])
        self.assertAlmostEqual(int(op_vals[25] * 1000 + 0.5), 1000)

        # Hexagonal planar motif.
        op_vals = ops_101.get_order_parameters(
            self.hexagonal_planar, 0, indices_neighs=[1,2,3,4,5,6])
        self.assertAlmostEqual(int(op_vals[26] * 1000 + 0.5), 1000)

        # Test providing explicit neighbor lists.
        op_vals = ops_101.get_order_parameters(self.bcc, 0, indices_neighs=[1])
        self.assertIsNotNone(op_vals[0])
        self.assertIsNone(op_vals[3])
        with self.assertRaises(ValueError):
            ops_101.get_order_parameters(self.bcc, 0, indices_neighs=[2])


    def tearDown(self):
        del self.single_bond
        del self.linear
        del self.bent45
        del self.cubic
        del self.fcc
        del self.bcc
        del self.hcp
        del self.diamond
        del self.regular_triangle
        del self.square
        del self.square_pyramid
        del self.trigonal_off_plane
        del self.trigonal_pyramid
        del self.trigonal_planar
        del self.square_planar
        del self.pentagonal_pyramid
        del self.hexagonal_pyramid
        del self.pentagonal_bipyramid
        del self.T_shape
        del self.cuboctahedron
        del self.see_saw_rect


if __name__ == '__main__':
    unittest.main()
