# coding: utf-8

from __future__ import unicode_literals

import numpy as np
import unittest
from math import fabs

from pymatgen.analysis.order_parameters import OrderParameters
from pymatgen.core.prototype_structures import PrototypeStructure
from pymatgen.util.testing import PymatgenTest


class OrderParametersTest(PymatgenTest):

    def setUp(self):
        self.cubic = PrototypeStructure("c", supercell_scaling=[3, 3, 3])
        self.cubic_envelop2 = PrototypeStructure(
                "c", supercell_scaling=[4, 4, 4], enveloping_box_factor=2.0)
        self.fcc_envelop2 = PrototypeStructure(
                "fcc", supercell_scaling=[4, 4, 4], enveloping_box_factor=2.0)
        self.bcc_envelop2 = PrototypeStructure(
                "bcc", supercell_scaling=[4, 4, 4], enveloping_box_factor=2.0)
        self.hcp_envelop2 = PrototypeStructure(
                "hcp", supercell_scaling=[4, 4, 4], enveloping_box_factor=2.0)
        self.diamond_envelop2 = PrototypeStructure(
                "d", supercell_scaling=[3, 3, 3],
                enveloping_box_factor=2.0)
        #tmp = PrototypeStructure("hcp") # , lengths=[1.0, 1.0, 1.0*1.633])
        #print(tmp.sites[0].coords)
        #print(tmp.sites[1].coords)
        #print(np.linalg.norm(tmp.sites[0].coords-tmp.sites[1].coords))
        #quit()


    def test_init(self):
        self.assertIsNotNone(OrderParameters(["cn"], [[]], 0.99))


    def test_get_order_parameters(self):
        # Set up everything for thorough testing.
        op_types = ["cn", "tet", "oct", "bcc", "q2", "q4", "q6"]
        op_paras = [[], [], [], [], [], [], []]
        ops_044 = OrderParameters(op_types, op_paras, 0.44)
        ops_071 = OrderParameters(op_types, op_paras, 0.71)
        ops_087 = OrderParameters(op_types, op_paras, 0.87)
        ops_099 = OrderParameters(op_types, op_paras, 0.99)
        ops_101 = OrderParameters(op_types, op_paras, 1.01)
        ops_voro = OrderParameters(op_types, op_paras)

        # Test providing explicit neighbor lists.
        op_vals = ops_101.get_order_parameters(self.cubic, 0, indeces_neighs=[1])
        self.assertIsNotNone(op_vals[0])
        self.assertIsNone(op_vals[1])
        with self.assertRaises(ValueError):
            ops_101.get_order_parameters(self.cubic, 0, indeces_neighs=[27])

        # Perfect cubic structures.
        for i in range(len(self.cubic.sites)):
            op_vals = ops_099.get_order_parameters(self.cubic, i)
            self.assertAlmostEqual(op_vals[0], 0.0)
            self.assertIsNone(op_vals[1])
            self.assertIsNone(op_vals[2])
            self.assertIsNone(op_vals[3])
            self.assertIsNone(op_vals[4])
            self.assertIsNone(op_vals[5])
            self.assertIsNone(op_vals[6])
            op_vals = ops_101.get_order_parameters(self.cubic, i)
            self.assertAlmostEqual(op_vals[0], 6.0)
            self.assertAlmostEqual(op_vals[1], -0.07206365327761823)
            self.assertAlmostEqual(op_vals[2], 1.0)
# prec prob            self.assertAlmostEqual(op_vals[3], 0.12499764616102801) # 0.125)
            self.assertAlmostEqual(op_vals[4], 0.0)
            self.assertAlmostEqual(op_vals[5], 0.7637626158259733)
            self.assertAlmostEqual(op_vals[6], 0.3535533905932738)

        # Cubic structure in large box.
        index = self.cubic_envelop2.get_index_of_site_closest_to_center_of_enveloping_cuboid()
        op_vals = ops_101.get_order_parameters(self.cubic_envelop2, index)
        self.assertAlmostEqual(op_vals[0], 6.0)
        self.assertAlmostEqual(op_vals[1], -0.07206365327761823) # reassured
        self.assertAlmostEqual(op_vals[2], 1.0)
#        self.assertAlmostEqual(op_vals[3], 0.12499764616102801)
        self.assertAlmostEqual(op_vals[4], 0.0)
        self.assertAlmostEqual(op_vals[5], 0.7637626158259733)
        self.assertAlmostEqual(op_vals[6], 0.3535533905932738)

        # Bcc structure in large box.
        index = self.bcc_envelop2.get_index_of_site_closest_to_center_of_enveloping_cuboid()
        op_vals = ops_087.get_order_parameters(self.bcc_envelop2, index)
        self.assertAlmostEqual(op_vals[0], 8.0)
        self.assertAlmostEqual(op_vals[1], 1.9688949183589557)
        self.assertAlmostEqual(op_vals[2], 0.12540815310925768)
        self.assertLess(fabs(op_vals[3]-0.9753713330608598), 1.0e-5)
        self.assertAlmostEqual(op_vals[4], 0.0)
        self.assertAlmostEqual(op_vals[5], 0.5091750772173156)
        self.assertAlmostEqual(op_vals[6], 0.6285393610547088)

        # Bcc structure in large box; neighbors via Voronoi facets.
        index = self.bcc_envelop2.get_index_of_site_closest_to_center_of_enveloping_cuboid()
        op_vals = ops_voro.get_order_parameters(self.bcc_envelop2, index)
        self.assertAlmostEqual(op_vals[0], 14.0)
        self.assertAlmostEqual(op_vals[1], -0.43855666168897867)
        self.assertAlmostEqual(op_vals[2], -1.2026322595911667)
        self.assertLess(fabs(op_vals[3]-0.11277709099511873), 1.0e-5)
        self.assertAlmostEqual(op_vals[4], 0.0)
        self.assertAlmostEqual(op_vals[5], 0.036369648372665375)
        self.assertAlmostEqual(op_vals[6], 0.5106882308569509)

        # Fcc structure in large box.
        index = self.fcc_envelop2.get_index_of_site_closest_to_center_of_enveloping_cuboid()
        op_vals = ops_071.get_order_parameters(self.fcc_envelop2, index)
        self.assertAlmostEqual(op_vals[0], 12.0)
        self.assertAlmostEqual(op_vals[1], -0.9989621462333275)
        self.assertAlmostEqual(op_vals[2], -1.0125484381377454)
        self.assertLess(fabs(op_vals[3]+0.0007417813723164877), 1.0e-5)
        self.assertAlmostEqual(op_vals[4], 0.0)
        self.assertAlmostEqual(op_vals[5], 0.1909406539564932)
        self.assertAlmostEqual(op_vals[6], 0.57452425971407)

        # Hcp structure in large box.
        index = self.hcp_envelop2.get_index_of_site_closest_to_center_of_enveloping_cuboid()
        op_vals = ops_101.get_order_parameters(self.hcp_envelop2, index)
        self.assertAlmostEqual(op_vals[0], 12.0)
        self.assertAlmostEqual(op_vals[1], -0.45564596177828504)
        self.assertAlmostEqual(op_vals[2], -0.7352383348310614)
        self.assertAlmostEqual(op_vals[3], -0.15591174007928388)
        self.assertAlmostEqual(op_vals[4], 9.72025056603629e-06)
        self.assertAlmostEqual(op_vals[5], 0.09722418406734)
        self.assertAlmostEqual(op_vals[6], 0.4847621378014553)

        # Diamond structure in large box.
        index = self.diamond_envelop2.get_index_of_site_closest_to_center_of_enveloping_cuboid()
        op_vals = ops_044.get_order_parameters(self.diamond_envelop2, index)
        self.assertAlmostEqual(op_vals[0], 4.0)
        self.assertAlmostEqual(op_vals[1], 1.0)
        self.assertAlmostEqual(op_vals[2], -0.008828657097338058)
        self.assertAlmostEqual(op_vals[3], 0.08087075431050145)
        self.assertAlmostEqual(op_vals[4], 0.0)
        self.assertAlmostEqual(op_vals[5], 0.5091750772173156)
        self.assertAlmostEqual(op_vals[6], 0.6285393610547088)


    def tearDown(self):
        del self.cubic
        del self.cubic_envelop2
        del self.fcc_envelop2
        del self.bcc_envelop2
        del self.hcp_envelop2
        del self.diamond_envelop2

if __name__ == '__main__':
    unittest.main()
