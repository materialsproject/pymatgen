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
        NearNeighbors
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

if __name__ == '__main__':
    unittest.main()
