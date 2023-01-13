# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import pickle

import numpy as np
import pytest

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.electronic_structure.core import Magmom
from pymatgen.util.testing import PymatgenTest


class SiteTest(PymatgenTest):
    def setUp(self):
        self.ordered_site = Site("Fe", [0.25, 0.35, 0.45])
        self.disordered_site = Site({"Fe": 0.5, "Mn": 0.5}, [0.25, 0.35, 0.45])
        self.propertied_site = Site("Fe2+", [0.25, 0.35, 0.45], {"magmom": 5.1, "charge": 4.2})
        self.propertied_magmomvector_site = Site(
            "Fe2+",
            [0.25, 0.35, 0.45],
            {"magmom": Magmom([2.6, 2.6, 3.5]), "charge": 4.2},
        )
        self.dummy_site = Site("X", [0, 0, 0])

    def test_properties(self):
        with pytest.raises(AttributeError):
            self.disordered_site.specie
        assert isinstance(self.ordered_site.specie, Element)
        assert self.propertied_site.properties["magmom"] == 5.1
        assert self.propertied_site.properties["charge"] == 4.2

    def test_to_from_dict(self):
        d = self.disordered_site.as_dict()
        site = Site.from_dict(d)
        assert site == self.disordered_site
        assert site != self.ordered_site
        d = self.propertied_site.as_dict()
        site = Site.from_dict(d)
        assert site.properties["magmom"] == 5.1
        assert site.properties["charge"] == 4.2
        d = self.propertied_magmomvector_site.as_dict()
        site = Site.from_dict(d)
        assert site.properties["magmom"] == Magmom([2.6, 2.6, 3.5])
        assert site.properties["charge"] == 4.2
        d = self.dummy_site.as_dict()
        site = Site.from_dict(d)
        assert site.species == self.dummy_site.species

    def test_hash(self):
        assert hash(self.ordered_site) == 26
        assert hash(self.disordered_site) == 51

    def test_cmp(self):
        assert self.ordered_site > self.disordered_site

    def test_distance(self):
        osite = self.ordered_site
        assert round(abs(np.linalg.norm([0.25, 0.35, 0.45]) - osite.distance_from_point([0, 0, 0])), 7) == 0
        assert round(abs(osite.distance(self.disordered_site) - 0), 7) == 0

    def test_pickle(self):
        o = pickle.dumps(self.propertied_site)
        assert pickle.loads(o) == self.propertied_site

    def test_setters(self):
        self.disordered_site.species = "Cu"
        assert self.disordered_site.species == Composition("Cu")
        self.disordered_site.x = 1.25
        self.disordered_site.y = 1.35
        assert self.disordered_site.coords[0] == 1.25
        assert self.disordered_site.coords[1] == 1.35

        def set_bad_species():
            self.disordered_site.species = {"Cu": 0.5, "Gd": 0.6}

        with pytest.raises(ValueError):
            set_bad_species()


class PeriodicSiteTest(PymatgenTest):
    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.si = Element("Si")
        self.site = PeriodicSite("Fe", [0.25, 0.35, 0.45], self.lattice)
        self.site2 = PeriodicSite({"Si": 0.5}, [0, 0, 0], self.lattice)
        assert self.site2.species == Composition({Element("Si"): 0.5}), "Inconsistent site created!"
        self.propertied_site = PeriodicSite(
            Species("Fe", 2),
            [0.25, 0.35, 0.45],
            self.lattice,
            properties={"magmom": 5.1, "charge": 4.2},
        )
        self.dummy_site = PeriodicSite("X", [0, 0, 0], self.lattice)

    def test_properties(self):
        """
        Test the properties for a site
        """
        assert self.site.a == 0.25
        assert self.site.b == 0.35
        assert self.site.c == 0.45
        assert self.site.x == 2.5
        assert self.site.y == 3.5
        assert self.site.z == 4.5
        assert self.site.is_ordered
        assert not self.site2.is_ordered
        assert self.propertied_site.properties["magmom"] == 5.1
        assert self.propertied_site.properties["charge"] == 4.2

    def test_distance(self):
        other_site = PeriodicSite("Fe", np.array([0, 0, 0]), self.lattice)
        assert round(abs(self.site.distance(other_site) - 6.22494979899), 5) == 0

    def test_distance_from_point(self):
        assert round(abs(self.site.distance_from_point([0.1, 0.1, 0.1]) - 6.22494979899), 5) != 0
        assert round(abs(self.site.distance_from_point([0.1, 0.1, 0.1]) - 6.0564015718906887), 5) == 0

    def test_distance_and_image(self):
        other_site = PeriodicSite("Fe", np.array([1, 1, 1]), self.lattice)
        distance, image = self.site.distance_and_image(other_site)
        assert round(abs(distance - 6.22494979899), 5) == 0
        assert ([-1, -1, -1] == image).all()
        distance, image = self.site.distance_and_image(other_site, [1, 0, 0])
        assert round(abs(distance - 19.461500456028563), 5) == 0
        # Test that old and new distance algo give the same ans for
        # "standard lattices"
        lattice = Lattice(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        site1 = PeriodicSite("Fe", np.array([0.01, 0.02, 0.03]), lattice)
        site2 = PeriodicSite("Fe", np.array([0.99, 0.98, 0.97]), lattice)
        assert round(abs(get_distance_and_image_old(site1, site2)[0] - site1.distance_and_image(site2)[0]), 7) == 0
        lattice = Lattice.from_parameters(1, 0.01, 1, 10, 10, 10)
        site1 = PeriodicSite("Fe", np.array([0.01, 0.02, 0.03]), lattice)
        site2 = PeriodicSite("Fe", np.array([0.99, 0.98, 0.97]), lattice)
        assert get_distance_and_image_old(site1, site2)[0] > site1.distance_and_image(site2)[0]
        site2 = PeriodicSite("Fe", np.random.rand(3), lattice)
        dist_old, jimage_old = get_distance_and_image_old(site1, site2)
        dist_new, jimage_new = site1.distance_and_image(site2)
        assert dist_old - dist_new > -1e-8, "New distance algo should give smaller answers!"
        assert (
            not (abs(dist_old - dist_new) < 1e-8) ^ (jimage_old == jimage_new).all()
        ), "If old dist == new dist, images must be the same!"
        latt = Lattice.from_parameters(3.0, 3.1, 10.0, 2.96, 2.0, 1.0)
        site = PeriodicSite("Fe", [0.1, 0.1, 0.1], latt)
        site2 = PeriodicSite("Fe", [0.99, 0.99, 0.99], latt)
        dist, img = site.distance_and_image(site2)
        assert round(abs(dist - 0.15495358379511573), 7) == 0
        assert list(img) == [-11, 6, 0]

    def test_is_periodic_image(self):
        other = PeriodicSite("Fe", np.array([1.25, 2.35, 4.45]), self.lattice)
        assert self.site.is_periodic_image(other), "This other site should be a periodic image."
        other = PeriodicSite("Fe", np.array([1.25, 2.35, 4.46]), self.lattice)
        assert not self.site.is_periodic_image(other), "This other site should not be a periodic image."
        other = PeriodicSite("Fe", np.array([1.25, 2.35, 4.45]), Lattice.rhombohedral(2, 60))
        assert not self.site.is_periodic_image(other), "Different lattices should not be periodic images."

    def test_equality(self):
        assert self.site == self.site
        other_site = PeriodicSite("Fe", np.array([1, 1, 1]), self.lattice)
        assert other_site != self.site

    def test_as_from_dict(self):
        d = self.site2.as_dict()
        site = PeriodicSite.from_dict(d)
        assert site == self.site2
        assert site != self.site
        d = self.propertied_site.as_dict()
        site3 = PeriodicSite({"Si": 0.5, "Fe": 0.5}, [0, 0, 0], self.lattice)
        d = site3.as_dict()
        site = PeriodicSite.from_dict(d)
        assert site.species == site3.species

        d = self.dummy_site.as_dict()
        site = PeriodicSite.from_dict(d)
        assert site.species == self.dummy_site.species

    def test_to_unit_cell(self):
        site = PeriodicSite("Fe", np.array([1.25, 2.35, 4.46]), self.lattice)
        site.to_unit_cell(in_place=True)
        val = [0.25, 0.35, 0.46]
        self.assertArrayAlmostEqual(site.frac_coords, val)

        lattice_pbc = Lattice(self.lattice.matrix, pbc=(True, True, False))
        site = PeriodicSite("Fe", np.array([1.25, 2.35, 4.46]), lattice_pbc)
        site.to_unit_cell(in_place=True)
        val = [0.25, 0.35, 4.46]
        self.assertArrayAlmostEqual(site.frac_coords, val)

    def test_setters(self):
        site = self.propertied_site
        site.species = "Cu"
        assert site.species == Composition("Cu")
        site.x = 1.25
        site.y = 1.35
        assert site.coords[0] == 1.25
        assert site.coords[1] == 1.35
        assert site.a == 0.125
        assert site.b == 0.135
        site.lattice = Lattice.cubic(100)
        assert site.x == 12.5

        def set_bad_species():
            site.species = {"Cu": 0.5, "Gd": 0.6}

        with pytest.raises(ValueError):
            set_bad_species()

        site.frac_coords = [0, 0, 0.1]
        self.assertArrayAlmostEqual(site.coords, [0, 0, 10])
        site.coords = [1.5, 3.25, 5]
        self.assertArrayAlmostEqual(site.frac_coords, [0.015, 0.0325, 0.05])

    def test_repr(self):
        assert repr(self.propertied_site) == "PeriodicSite: Fe2+ (2.5000, 3.5000, 4.5000) [0.2500, 0.3500, 0.4500]"


def get_distance_and_image_old(site1, site2, jimage=None):
    """
    Gets distance between two sites assuming periodic boundary conditions.
    If the index jimage of two sites atom j is not specified it selects the
    j image nearest to the i atom and returns the distance and jimage
    indices in terms of lattice vector translations. If the index jimage of
    atom j is specified it returns the distance between the i atom and the
    specified jimage atom, the given jimage is also returned.

    Args:
        other:
            other site to get distance from.
        jimage:
            specific periodic image in terms of lattice translations,
            e.g., [1,0,0] implies to take periodic image that is one
            a-lattice vector away. If jimage is None, the image that is
            nearest to the site is found.

    Returns:
        (distance, jimage):
            distance and periodic lattice translations of the other site
            for which the distance applies.

    .. note::
        Assumes the primitive cell vectors are sufficiently not skewed such
        that the condition \\|a\\|cos(ab_angle) < \\|b\\| for all possible cell
        vector pairs. ** this method does not check this condition **
    """
    if jimage is None:
        # Old algorithm
        jimage = -np.array(np.around(site2.frac_coords - site1.frac_coords), int)
    mapped_vec = site1.lattice.get_cartesian_coords(jimage + site2.frac_coords - site1.frac_coords)
    dist = np.linalg.norm(mapped_vec)
    return dist, jimage


if __name__ == "__main__":
    import unittest

    unittest.main()
