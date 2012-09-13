#!/usr/bin/env python

'''
Created on Jul 17, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 17, 2012"

import unittest
import numpy as np

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.core.lattice import Lattice

class SiteTest(unittest.TestCase):

    def setUp(self):
        self.ordered_site = Site(Element("Fe"), [0.25, 0.35, 0.45])
        self.disordered_site = Site({Element("Fe"):0.5, Element("Mn"):0.5}, [0.25, 0.35, 0.45])
        self.propertied_site = Site(Specie("Fe", 2), [0.25, 0.35, 0.45], {'magmom':5.1, 'charge':4.2})

    def test_init(self):
        self.assertRaises(ValueError, Site, Specie("Fe", 2), [0.25, 0.35, 0.45], {'mag':5.1})

    def test_properties(self):
        self.assertRaises(AttributeError, getattr, self.disordered_site, 'specie')
        self.assertIsInstance(self.ordered_site.specie, Element)
        self.assertEqual(self.propertied_site.magmom, 5.1)
        self.assertEqual(self.propertied_site.charge, 4.2)

    def test_to_from_dict(self):
        d = self.disordered_site.to_dict
        site = Site.from_dict(d)
        self.assertEqual(site, self.disordered_site)
        self.assertNotEqual(site, self.ordered_site)
        d = self.propertied_site.to_dict
        site = Site.from_dict(d)
        self.assertEqual(site.magmom, 5.1)
        self.assertEqual(site.charge, 4.2)


class PeriodicSiteTest(unittest.TestCase):

    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.si = Element("Si")
        self.site = PeriodicSite("Fe", np.array([0.25, 0.35, 0.45]),
                                 self.lattice)
        self.site2 = PeriodicSite({"Si":0.5}, np.array([0, 0, 0]), self.lattice)
        self.assertEquals(self.site2.species_and_occu, {Element('Si'): 0.5},
                          "Inconsistent site created!")
        self.propertied_site = PeriodicSite(Specie("Fe", 2), [0.25, 0.35, 0.45],
                                            self.lattice,
                                            properties={'magmom':5.1,
                                                        'charge':4.2})

    def test_properties(self):
        """
        Test the properties for a site
        """
        self.assertEquals(self.site.a, 0.25)
        self.assertEquals(self.site.b, 0.35)
        self.assertEquals(self.site.c, 0.45)
        self.assertEquals(self.site.x, 2.5)
        self.assertEquals(self.site.y, 3.5)
        self.assertEquals(self.site.z, 4.5)
        self.assertTrue(self.site.is_ordered)
        self.assertFalse(self.site2.is_ordered)
        self.assertEqual(self.propertied_site.magmom, 5.1)
        self.assertEqual(self.propertied_site.charge, 4.2)

    def test_distance(self):
        other_site = PeriodicSite("Fe", np.array([0, 0, 0]), self.lattice)
        self.assertAlmostEquals(self.site.distance(other_site), 6.22494979899, 5)

    def test_distance_from_point(self):
        self.assertNotAlmostEqual(self.site.distance_from_point(np.array([0.1, 0.1, 0.1])), 6.22494979899, 5)
        self.assertAlmostEqual(self.site.distance_from_point(np.array([0.1, 0.1, 0.1])), 6.0564015718906887, 5)

    def test_distance_and_image(self):
        other_site = PeriodicSite("Fe", np.array([1, 1, 1]), self.lattice)
        (distance, image) = self.site.distance_and_image(other_site)
        self.assertAlmostEquals(distance, 6.22494979899, 5)
        self.assertTrue(([-1, -1, -1] == image).all())
        (distance, image) = self.site.distance_and_image(other_site, [1, 0, 0])
        self.assertAlmostEquals(distance, 19.461500456028563, 5)
        # Test that old and new distance algo give the same ans for "standard lattices"
        lattice = Lattice(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        site1 = PeriodicSite("Fe", np.array([0.01, 0.02, 0.03]), lattice)
        site2 = PeriodicSite("Fe", np.array([0.99, 0.98, 0.97]), lattice)
        self.assertAlmostEqual(site1.distance_and_image_old(site2)[0], site1.distance_and_image(site2)[0])
        lattice = Lattice.from_parameters(1, 0.01, 1, 10, 10, 10)
        site1 = PeriodicSite("Fe", np.array([0.01, 0.02, 0.03]), lattice)
        site2 = PeriodicSite("Fe", np.array([0.99, 0.98, 0.97]), lattice)
        self.assertTrue(site1.distance_and_image_old(site2)[0] > site1.distance_and_image(site2)[0])
        site2 = PeriodicSite("Fe", np.random.rand(3), lattice)
        (dist_old, jimage_old) = site1.distance_and_image_old(site2)
        (dist_new, jimage_new) = site1.distance_and_image(site2)
        self.assertTrue(dist_old - dist_new > -1e-8, "New distance algo should always give smaller answers!")
        self.assertFalse((abs(dist_old - dist_new) < 1e-8) ^ (jimage_old == jimage_new).all(), "If old dist == new dist, the images returned must be the same!")

    def test_is_periodic_image(self):
        other = PeriodicSite("Fe", np.array([1.25, 2.35, 4.45]), self.lattice)
        self.assertTrue(self.site.is_periodic_image(other),
                        "This other site should be a periodic image.")
        other = PeriodicSite("Fe", np.array([1.25, 2.35, 4.46]), self.lattice)
        self.assertFalse(self.site.is_periodic_image(other),
                         "This other site should not be a periodic image.")
        other = PeriodicSite("Fe", np.array([1.25, 2.35, 4.45]),
                             Lattice.rhombohedral(2, 60))
        self.assertFalse(self.site.is_periodic_image(other),
                         "Different lattices should result in different periodic sites.")

    def test_equality(self):
        other_site = PeriodicSite("Fe", np.array([1, 1, 1]), self.lattice)
        self.assertTrue(self.site.__eq__(self.site))
        self.assertFalse(other_site.__eq__(self.site))
        self.assertFalse(self.site.__ne__(self.site))
        self.assertTrue(other_site.__ne__(self.site))

    def test_to_from_dict(self):
        d = self.site2.to_dict
        site = PeriodicSite.from_dict(d)
        self.assertEqual(site, self.site2)
        self.assertNotEqual(site, self.site)
        d = self.propertied_site.to_dict
        site = Site.from_dict(d)
        self.assertEqual(site.magmom, 5.1)
        self.assertEqual(site.charge, 4.2)
        site3 = PeriodicSite({"Si":0.5, "Fe":0.5}, np.array([0, 0, 0]), self.lattice)
        d = site3.to_dict
        site = PeriodicSite.from_dict(d)
        self.assertEqual(site.species_and_occu, site3.species_and_occu)

    def test_to_unit_cell(self):
        site = PeriodicSite("Fe", np.array([1.25, 2.35, 4.46]), self.lattice)
        site = site.to_unit_cell
        val = [0.25, 0.35, 0.46]
        self.assertTrue(np.allclose(site.frac_coords, val))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
