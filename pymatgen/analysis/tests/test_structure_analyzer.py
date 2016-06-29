# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np
import unittest2 as unittest
import os

from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder, \
    solid_angle, contains_peroxide, RelaxationAnalyzer, VoronoiConnectivity, \
    oxide_type, sulfide_type, OrderParameters, average_coordination_number, \
    VoronoiAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen import Element, Structure, Lattice
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class VoronoiCoordFinderTest(PymatgenTest):
    def setUp(self):
        s = self.get_structure('LiFePO4')
        self.finder = VoronoiCoordFinder(s, [Element("O")])

    def test_get_voronoi_polyhedra(self):
        self.assertEqual(len(self.finder.get_voronoi_polyhedra(0).items()), 8)

    def test_get_coordination_number(self):
        self.assertAlmostEqual(self.finder.get_coordination_number(0),
                               5.809265748999465, 7)

    def test_get_coordinated_sites(self):
        self.assertEqual(len(self.finder.get_coordinated_sites(0)), 8)


class VoronoiAnalyzerTest(PymatgenTest):

    def setUp(self):
        self.ss = Xdatcar(os.path.join(test_dir, 'XDATCAR.MD')).structures
        self.s = self.ss[1]
        self.va = VoronoiAnalyzer(cutoff=4.0)

    def test_analyze(self):
        # Check for the Voronoi index of site i in Structure
        single_structure = self.va.analyze(self.s, n=5)
        self.assertIn(single_structure.view(),
                      np.array([4, 3, 3, 4, 2, 2, 1, 0]).view(),
                      "Cannot find the right polyhedron.")
        # Check for the presence of a Voronoi index and its frequency in
        # a ensemble (list) of Structures
        ensemble = self.va.analyze_structures(self.ss, step_freq=2,
                                              most_frequent_polyhedra=10)
        self.assertIn(('[1 3 4 7 1 0 0 0]', 3),
                      ensemble, "Cannot find the right polyhedron in ensemble.")


class RelaxationAnalyzerTest(unittest.TestCase):
    def setUp(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR.Li2O'),
                             check_for_POTCAR=False)
        s1 = p.structure
        p = Poscar.from_file(os.path.join(test_dir, 'CONTCAR.Li2O'),
                             check_for_POTCAR=False)
        s2 = p.structure
        self.analyzer = RelaxationAnalyzer(s1, s2)

    def test_vol_and_para_changes(self):
        for k, v in self.analyzer.get_percentage_lattice_parameter_changes().items():
            self.assertAlmostEqual(-0.0092040921155279731, v)
            latt_change = v
        vol_change = self.analyzer.get_percentage_volume_change()
        self.assertAlmostEqual(-0.0273589101391,
                               vol_change)
        # This is a simple cubic cell, so the latt and vol change are simply
        # Related. So let's test that.
        self.assertAlmostEqual((1 + latt_change) ** 3 - 1, vol_change)

    def test_get_percentage_bond_dist_changes(self):
        for k, v in self.analyzer.get_percentage_bond_dist_changes().items():
            for k2, v2 in v.items():
                self.assertAlmostEqual(-0.009204092115527862, v2)


class VoronoiConnectivityTest(PymatgenTest):
    def test_connectivity_array(self):
        vc = VoronoiConnectivity(self.get_structure("LiFePO4"))
        ca = vc.connectivity_array
        np.set_printoptions(threshold=np.NAN, linewidth=np.NAN, suppress=np.NAN)

        expected = np.array([0, 1.96338392, 0, 0.04594495])
        self.assertTrue(np.allclose(ca[15, :4, ca.shape[2] // 2], expected))

        expected = np.array([0, 0, 0])
        self.assertTrue(np.allclose(ca[1, -3:, 51], expected))

        site = vc.get_sitej(27, 51)
        self.assertEqual(site.specie, Element('O'))
        expected = np.array([-0.29158, 0.74889, 0.95684])
        self.assertTrue(np.allclose(site.frac_coords, expected))


class MiscFunctionTest(PymatgenTest):
    def test_average_coordination_number(self):
        xdatcar = Xdatcar(os.path.join(test_dir, 'XDATCAR.MD'))
        coordination_numbers = average_coordination_number(xdatcar.structures,
                                                           freq=1)
        self.assertAlmostEqual(coordination_numbers['Fe'], 4.771903318390836, 5,
                               "Coordination number not calculated properly.")

    def test_solid_angle(self):
        center = [2.294508207929496, 4.4078057081404, 2.299997773791287]
        coords = [[1.627286218099362, 3.081185538926995, 3.278749383217061],
                  [1.776793751092763, 2.93741167455471, 3.058701096568852],
                  [3.318412187495734, 2.997331084033472, 2.022167590167672],
                  [3.874524708023352, 4.425301459451914, 2.771990305592935],
                  [2.055778446743566, 4.437449313863041, 4.061046832034642]]
        self.assertAlmostEqual(solid_angle(center, coords), 1.83570965938, 7,
                               "Wrong result returned by solid_angle")

    def test_contains_peroxide(self):

        for f in ['LiFePO4', 'NaFePO4', 'Li3V2(PO4)3', 'Li2O']:
            self.assertFalse(contains_peroxide(self.get_structure(f)))

        for f in ['Li2O2', "K2O2"]:
            self.assertTrue(contains_peroxide(self.get_structure(f)))

    def test_oxide_type(self):
        el_li = Element("Li")
        el_o = Element("O")
        latt = Lattice([[3.985034, 0.0, 0.0],
                        [0.0, 4.881506, 0.0],
                        [0.0, 0.0, 2.959824]])
        elts = [el_li, el_li, el_o, el_o, el_o, el_o]
        coords = list()
        coords.append([0.500000, 0.500000, 0.500000])
        coords.append([0.0, 0.0, 0.0])
        coords.append([0.632568, 0.085090, 0.500000])
        coords.append([0.367432, 0.914910, 0.500000])
        coords.append([0.132568, 0.414910, 0.000000])
        coords.append([0.867432, 0.585090, 0.000000])
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "superoxide")

        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_o, el_o, el_o]
        latt = Lattice.from_parameters(3.999911, 3.999911, 3.999911, 133.847504,
                                       102.228244, 95.477342)
        coords = [[0.513004, 0.513004, 1.000000],
                  [0.017616, 0.017616, 0.000000],
                  [0.649993, 0.874790, 0.775203],
                  [0.099587, 0.874790, 0.224797]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "ozonide")

        latt = Lattice.from_parameters(3.159597, 3.159572, 7.685205, 89.999884,
                                       89.999674, 60.000510)
        el_li = Element("Li")
        el_o = Element("O")
        elts = [el_li, el_li, el_li, el_li, el_o, el_o, el_o, el_o]
        coords = [[0.666656, 0.666705, 0.750001],
                  [0.333342, 0.333378, 0.250001],
                  [0.000001, 0.000041, 0.500001],
                  [0.000001, 0.000021, 0.000001],
                  [0.333347, 0.333332, 0.649191],
                  [0.333322, 0.333353, 0.850803],
                  [0.666666, 0.666686, 0.350813],
                  [0.666665, 0.666684, 0.149189]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "peroxide")

        el_li = Element("Li")
        el_o = Element("O")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000,
                                       90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_o, el_o]
        coords = [[0.000000, 0.500000, 0.413969],
                  [0.500000, 0.000000, 0.586031],
                  [0.000000, 0.000000, 0.000000],
                  [0.500000, 0.500000, 0.000000],
                  [0.000000, 0.500000, 0.192672],
                  [0.500000, 0.000000, 0.807328]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "hydroxide")

        el_li = Element("Li")
        el_n = Element("N")
        el_h = Element("H")
        latt = Lattice.from_parameters(3.565276, 3.565276, 4.384277, 90.000000,
                                       90.000000, 90.000000)
        elts = [el_h, el_h, el_li, el_li, el_n, el_n]
        coords = [[0.000000, 0.500000, 0.413969],
                  [0.500000, 0.000000, 0.586031],
                  [0.000000, 0.000000, 0.000000],
                  [0.500000, 0.500000, 0.000000],
                  [0.000000, 0.500000, 0.192672],
                  [0.500000, 0.000000, 0.807328]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "None")

        el_o = Element("O")
        latt = Lattice.from_parameters(4.389828, 5.369789, 5.369789, 70.786622,
                                       69.244828, 69.244828)
        elts = [el_o, el_o, el_o, el_o, el_o, el_o, el_o, el_o]
        coords = [[0.844609, 0.273459, 0.786089],
                  [0.155391, 0.213911, 0.726541],
                  [0.155391, 0.726541, 0.213911],
                  [0.844609, 0.786089, 0.273459],
                  [0.821680, 0.207748, 0.207748],
                  [0.178320, 0.792252, 0.792252],
                  [0.132641, 0.148222, 0.148222],
                  [0.867359, 0.851778, 0.851778]]
        struct = Structure(latt, elts, coords)
        self.assertEqual(oxide_type(struct, 1.1), "None")

    def test_sulfide_type(self):
        # NaS2 -> polysulfide
        latt = Lattice.tetragonal(9.59650, 11.78850)
        species = ["Na"] * 2 + ["S"] * 2
        coords = [[0.00000, 0.00000, 0.17000],
                  [0.27600, 0.25000, 0.12500],
                  [0.03400, 0.25000, 0.29600],
                  [0.14700, 0.11600, 0.40000]]
        struct = Structure.from_spacegroup(122, latt, species, coords)
        self.assertEqual(sulfide_type(struct), "polysulfide")

        # NaCl type NaS -> sulfide
        latt = Lattice.cubic(5.75)
        species = ["Na", "S"]
        coords = [[0.00000, 0.00000, 0.00000],
                  [0.50000, 0.50000, 0.50000]]
        struct = Structure.from_spacegroup(225, latt, species, coords)
        self.assertEqual(sulfide_type(struct), "sulfide")

        # Na2S2O3 -> None (sulfate)
        latt = Lattice.monoclinic(6.40100, 8.10000, 8.47400, 96.8800)
        species = ["Na"] * 2 + ["S"] * 2 + ["O"] * 3
        coords = [[0.29706, 0.62396, 0.08575],
                  [0.37673, 0.30411, 0.45416],
                  [0.52324, 0.10651, 0.21126],
                  [0.29660, -0.04671, 0.26607],
                  [0.17577, 0.03720, 0.38049],
                  [0.38604, -0.20144, 0.33624],
                  [0.16248, -0.08546, 0.11608]]
        struct = Structure.from_spacegroup(14, latt, species, coords)
        self.assertEqual(sulfide_type(struct), None)

        # Na3PS3O -> sulfide
        latt = Lattice.orthorhombic(9.51050, 11.54630, 5.93230)
        species = ["Na"] * 2 + ["S"] * 2 + ["P", "O"]
        coords = [[0.19920, 0.11580, 0.24950],
                  [0.00000, 0.36840, 0.29380],
                  [0.32210, 0.36730, 0.22530],
                  [0.50000, 0.11910, 0.27210],
                  [0.50000, 0.29400, 0.35500],
                  [0.50000, 0.30300, 0.61140]]
        struct = Structure.from_spacegroup(36, latt, species, coords)
        self.assertEqual(sulfide_type(struct), "sulfide")


class OrderParametersTest(PymatgenTest):
    def setUp(self):
        self.cubic = Structure(
            Lattice.from_lengths_and_angles(
                [1.0, 1.0, 1.0], [90.0, 90.0, 90.0]),
            ["H"], [[0.0, 0.0, 0.0]], validate_proximity=False,
            to_unit_cell=False, coords_are_cartesian=False,
            site_properties=None)
        self.bcc = Structure(
            Lattice.from_lengths_and_angles(
                [1.0, 1.0, 1.0], [90.0, 90.0, 90.0]),
            ["H", "H"], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.fcc = Structure(
            Lattice.from_lengths_and_angles(
                [1.0, 1.0, 1.0], [90.0, 90.0, 90.0]),
            ["H", "H", "H", "H"], [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                                   [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.hcp = Structure(
            Lattice.from_lengths_and_angles(
                [1.0, 1.0, 1.633], [90.0, 90.0, 120.0]),
            ["H", "H"],
            [[0.3333, 0.6667, 0.25], [0.6667, 0.3333, 0.75]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)
        self.diamond = Structure(
            Lattice.from_lengths_and_angles(
                [1.0, 1.0, 1.0], [90.0, 90.0, 90.0]),
            ["H", "H", "H", "H", "H", "H", "H", "H"],
            [[0.0, 0.0, 0.5], [0.75, 0.75, 0.75],
             [0.0, 0.5, 0.0], [0.75, 0.25, 0.25],
             [0.5, 0.0, 0.0], [0.25, 0.75, 0.25],
             [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]],
            validate_proximity=False, to_unit_cell=False,
            coords_are_cartesian=False, site_properties=None)

    def test_init(self):
        self.assertIsNotNone(OrderParameters(["cn"], [[]], 0.99))

    def test_get_order_parameters(self):
        # Set up everything.
        op_types = ["cn", "tet", "oct", "bcc", "q2", "q4", "q6"]
        op_paras = [[], [], [], [], [], [], []]
        ops_044 = OrderParameters(op_types, op_paras, 0.44)
        ops_071 = OrderParameters(op_types, op_paras, 0.71)
        ops_087 = OrderParameters(op_types, op_paras, 0.87)
        ops_099 = OrderParameters(op_types, op_paras, 0.99)
        ops_101 = OrderParameters(op_types, op_paras, 1.01)
        ops_voro = OrderParameters(op_types, op_paras)

        # Cubic structure.
        op_vals = ops_099.get_order_parameters(self.cubic, 0)
        self.assertAlmostEqual(op_vals[0], 0.0)
        self.assertIsNone(op_vals[1])
        self.assertIsNone(op_vals[2])
        self.assertIsNone(op_vals[3])
        self.assertIsNone(op_vals[4])
        self.assertIsNone(op_vals[5])
        self.assertIsNone(op_vals[6])
        op_vals = ops_101.get_order_parameters(self.cubic, 0)
        self.assertAlmostEqual(op_vals[0], 6.0)
        self.assertAlmostEqual(int(op_vals[1] * 1000), -72)
        self.assertAlmostEqual(int(op_vals[2] * 1000), 1000)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 125)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[5] * 1000), 763)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 353)

        # Bcc structure.
        op_vals = ops_087.get_order_parameters(self.bcc, 0)
        self.assertAlmostEqual(op_vals[0], 8.0)
        self.assertAlmostEqual(int(op_vals[1] * 1000),
                               1968)  # 1.9688949183589557
        self.assertAlmostEqual(int(op_vals[2] * 1000),
                               125)  # 0.12540815310925768
        self.assertAlmostEqual(int(op_vals[3] * 1000),
                               975)  # 0.9753713330608598
        self.assertAlmostEqual(int(op_vals[4] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[5] * 1000),
                               509)  # 0.5091750772173156
        self.assertAlmostEqual(int(op_vals[6] * 1000),
                               628)  # 0.6285393610547088

        # Fcc structure.
        op_vals = ops_071.get_order_parameters(self.fcc, 0)
        self.assertAlmostEqual(op_vals[0], 12.0)
        self.assertAlmostEqual(int(op_vals[1] * 1000),
                               -998)  # -0.9989621462333275
        self.assertAlmostEqual(int(op_vals[2] * 1000),
                               -1012)  # -1.0125484381377454
        self.assertAlmostEqual(int(op_vals[3] * 1000),
                               0)  # -0.0007417813723164877
        self.assertAlmostEqual(int(op_vals[4] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[5] * 1000),
                               190)  # 0.1909406539564932
        self.assertAlmostEqual(int(op_vals[6] * 1000), 574)  # 0.57452425971407

        # Hcp structure.
        op_vals = ops_101.get_order_parameters(self.hcp, 0)
        self.assertAlmostEqual(op_vals[0], 12.0)
        self.assertAlmostEqual(int(op_vals[1] * 1000), -455)
        self.assertAlmostEqual(int(op_vals[2] * 1000), -735)
        self.assertAlmostEqual(int(op_vals[3] * 1000), -155)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[5] * 1000), 97)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 484)

        # Diamond structure.
        op_vals = ops_044.get_order_parameters(self.diamond, 0)
        self.assertAlmostEqual(op_vals[0], 4.0)
        self.assertAlmostEqual(int(op_vals[1] * 1000), 1000)
        self.assertAlmostEqual(int(op_vals[2] * 1000), -8)
        self.assertAlmostEqual(int(op_vals[3] * 1000), 80)
        self.assertAlmostEqual(int(op_vals[4] * 1000), 0)
        self.assertAlmostEqual(int(op_vals[5] * 1000), 509)
        self.assertAlmostEqual(int(op_vals[6] * 1000), 628)

        # Test providing explicit neighbor lists.
        op_vals = ops_101.get_order_parameters(self.bcc, 0, indeces_neighs=[1])
        self.assertIsNotNone(op_vals[0])
        self.assertIsNone(op_vals[1])
        with self.assertRaises(ValueError):
            ops_101.get_order_parameters(self.bcc, 0, indeces_neighs=[2])

    def tearDown(self):
        del self.cubic
        del self.fcc
        del self.bcc
        del self.hcp
        del self.diamond


if __name__ == '__main__':
    unittest.main()
