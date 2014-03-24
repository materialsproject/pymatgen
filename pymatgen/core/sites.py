#!/usr/bin/env python

"""
This module defines classes representing non-periodic and periodic sites.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 17, 2012"

import collections
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, Specie, DummySpecie,\
    get_el_sp
from pymatgen.serializers.json_coders import MSONable
from pymatgen.util.coord_utils import pbc_diff
from pymatgen.core.composition import Composition


class Site(collections.Mapping, collections.Hashable, MSONable):
    """
    A generalized *non-periodic* site. This is essentially a composition
    at a point in space, with some optional properties associated with it. A
    Composition is used to represent the atoms and occupancy, which allows for
    disordered site representation. Coords are given in standard cartesian
    coordinates.
    """

    position_atol = 1e-5

    def __init__(self, atoms_n_occu, coords, properties=None):
        """
        Create a *non-periodic* site.

        Args:
            atoms_n_occu: Species on the site. Can be:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.
                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords: Cartesian coordinates of site.
            properties: Properties associated with the site as a dict, e.g.
                {"magmom": 5}. Defaults to None.
        """
        if issubclass(atoms_n_occu.__class__, collections.Mapping):
            self._species = Composition({get_el_sp(k): v
                                         for k, v in atoms_n_occu.items()})
            totaloccu = self._species.num_atoms
            if totaloccu > 1 + Composition.amount_tolerance:
                raise ValueError("Species occupancies sum to more than 1!")
            self._is_ordered = totaloccu == 1 and len(self._species) == 1
        else:
            self._species = Composition({get_el_sp(atoms_n_occu): 1})
            self._is_ordered = True

        self._coords = coords
        self._properties = properties if properties else {}

    @property
    def properties(self):
        """
        Returns a view of properties as a dict.
        """
        return {k: v for k, v in self._properties.items()}

    def __getattr__(self, a):
        #overriding getattr doens't play nice with pickle, so we
        #can't use self._properties
        p = object.__getattribute__(self, '_properties')
        if a in p:
            return p[a]
        raise AttributeError(a)

    def distance(self, other):
        """
        Get distance between two sites.

        Args:
            other: Other site.

        Returns:
            Distance (float)
        """
        return np.linalg.norm(other.coords - self.coords)

    def distance_from_point(self, pt):
        """
        Returns distance between the site and a point in space.

        Args:
            pt: Cartesian coordinates of point.

        Returns:
            Distance (float)
        """
        return np.linalg.norm(np.array(pt) - self._coords)

    @property
    def species_string(self):
        """
        String representation of species on the site.
        """
        if self._is_ordered:
            return str(self._species.keys()[0])
        else:
            sorted_species = sorted(self._species.keys())
            return ", ".join(["{}:{:.3f}".format(sp, self._species[sp])
                              for sp in sorted_species])

    @property
    def species_and_occu(self):
        """
        The species at the site, i.e., a Composition mapping type of
        element/species to occupancy.
        """
        return self._species

    @property
    def specie(self):
        """
        The Specie/Element at the site. Only works for ordered sites. Otherwise
        an AttributeError is raised. Use this property sparingly.  Robust
        design should make use of the property species_and_occu instead.

        Raises:
            AttributeError if Site is not ordered.
        """
        if not self._is_ordered:
            raise AttributeError("specie property only works for ordered "
                                 "sites!")
        return self._species.keys()[0]

    @property
    def coords(self):
        """
        A copy of the cartesian coordinates of the site as a numpy array.
        """
        return np.copy(self._coords)

    @property
    def is_ordered(self):
        """
        True if site is an ordered site, i.e., with a single species with
        occupancy 1.
        """
        return self._is_ordered

    @property
    def x(self):
        """
        Cartesian x coordinate
        """
        return self._coords[0]

    @property
    def y(self):
        """
        Cartesian y coordinate
        """
        return self._coords[1]

    @property
    def z(self):
        """
        Cartesian z coordinate
        """
        return self._coords[2]

    def __getitem__(self, el):
        """
        Get the occupancy for element
        """
        return self._species[el]

    def __eq__(self, other):
        """
        Site is equal to another site if the species and occupancies are the
        same, and the coordinates are the same to some tolerance.  numpy
        function `allclose` is used to determine if coordinates are close.
        """
        if other is None:
            return False
        return self._species == other._species and \
            np.allclose(self._coords, other._coords,
                        atol=Site.position_atol) and \
            self._properties == other._properties

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        hashcode = sum((el.Z * occu for el, occu in self._species.items()))
        return hashcode

    def __contains__(self, el):
        return el in self._species

    def __len__(self):
        return len(self._species)

    def __iter__(self):
        return self._species.__iter__()

    def __repr__(self):
        return "Site: {} ({:.4f}, {:.4f}, {:.4f})".format(
            self.species_string, *self._coords)

    def __cmp__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas.  For example, FeO4PLi is
        automatically sorted in LiFePO4.
        """

        if self._species.average_electroneg < \
                other._species.average_electroneg:
            return -1
        if self._species.average_electroneg > \
                other._species.average_electroneg:
            return 1
        if self.species_string < other.species_string:
            return -1
        if self.species_string > other.species_string:
            return 1
        return 0

    def __str__(self):
        return "{} {}".format(self._coords, self.species_string)

    @property
    def to_dict(self):
        """
        Json-serializable dict representation for Site.
        """
        species_list = []
        for spec, occu in self._species.items():
            d = spec.to_dict
            del d["@module"]
            del d["@class"]
            d["occu"] = occu
            species_list.append(d)
        return {"name": self.species_string, "species": species_list,
                "xyz": [float(c) for c in self._coords],
                "properties": self._properties,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        """
        Create Site from dict representation
        """
        atoms_n_occu = {}
        for sp_occu in d["species"]:
            if "oxidation_state" in sp_occu and Element.is_valid_symbol(
                    sp_occu["element"]):
                sp = Specie.from_dict(sp_occu)
            elif "oxidation_state" in sp_occu:
                sp = DummySpecie.from_dict(sp_occu)
            else:
                sp = Element(sp_occu["element"])
            atoms_n_occu[sp] = sp_occu["occu"]
        props = d.get("properties", None)
        return cls(atoms_n_occu, d["xyz"], properties=props)


class PeriodicSite(Site, MSONable):
    """
    Extension of generic Site object to periodic systems.
    PeriodicSite includes a lattice system.
    """

    def __init__(self, atoms_n_occu, coords, lattice, to_unit_cell=False,
                 coords_are_cartesian=False, properties=None):
        """
        Create a periodic site.

        Args:
            atoms_n_occu: Species on the site. Can be:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.
                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (3x1 array or sequence): Coordinates of site as fractional
                or cartesian coordinates.
            lattice: Lattice associated with the site
            to_unit_cell (bool): Translates fractional coordinate to the
                basic unit cell, i.e. all fractional coordinates satisfy 0
                <= a < 1. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                cartesian coordinates. Defaults to False.
            properties (dict): Properties associated with the PeriodicSite,
                e.g., {"magmom":5}. Defaults to None.
        """
        self._lattice = lattice
        if coords_are_cartesian:
            self._fcoords = self._lattice.get_fractional_coords(coords)
            c_coords = coords
        else:
            self._fcoords = coords
            c_coords = lattice.get_cartesian_coords(coords)

        if to_unit_cell:
            self._fcoords = np.mod(self._fcoords, 1)
            c_coords = lattice.get_cartesian_coords(self._fcoords)
        Site.__init__(self, atoms_n_occu, c_coords, properties)

    @property
    def lattice(self):
        """
        The lattice associated with the site.
        """
        return self._lattice

    @property
    def frac_coords(self):
        """
        A copy of the fractional coordinates of the site.
        """
        return np.copy(self._fcoords)

    @property
    def a(self):
        """
        Fractional a coordinate
        """
        return self._fcoords[0]

    @property
    def b(self):
        """
        Fractional b coordinate
        """
        return self._fcoords[1]

    @property
    def c(self):
        """
        Fractional c coordinate
        """
        return self._fcoords[2]

    @property
    def to_unit_cell(self):
        """
        Copy of PeriodicSite translated to the unit cell.
        """
        return PeriodicSite(self._species, np.mod(self._fcoords, 1),
                            self._lattice, properties=self._properties)

    def is_periodic_image(self, other, tolerance=1e-8, check_lattice=True):
        """
        Returns True if sites are periodic images of each other.

        Args:
            other (PeriodicSite): Other site
            tolerance (float): Tolerance to compare fractional coordinates
            check_lattice (bool): Whether to check if the two sites have the
                same lattice.

        Returns:
            bool: True if sites are periodic images of each other.
        """
        if check_lattice and self._lattice != other._lattice:
            return False
        if self._species != other._species:
            return False

        frac_diff = pbc_diff(self._fcoords, other._fcoords)
        return np.allclose(frac_diff, [0, 0, 0], atol=tolerance)

    def __eq__(self, other):
        return self._species == other._species and \
            self._lattice == other._lattice and \
            np.allclose(self._coords, other._coords,
                        atol=Site.position_atol) and \
            self._properties == other._properties

    def __ne__(self, other):
        return not self.__eq__(other)

    def distance_and_image_from_frac_coords(self, fcoords, jimage=None):
        """
        Gets distance between site and a fractional coordinate assuming
        periodic boundary conditions. If the index jimage of two sites atom j
        is not specified it selects the j image nearest to the i atom and
        returns the distance and jimage indices in terms of lattice vector
        translations. If the index jimage of atom j is specified it returns the
        distance between the i atom and the specified jimage atom, the given
        jimage is also returned.

        Args:
            fcoords (3x1 array): fcoords to get distance from.
            jimage (3x1 array): Specific periodic image in terms of
                lattice translations, e.g., [1,0,0] implies to take periodic
                image that is one a-lattice vector away. If jimage == None,
                the image that is nearest to the site is found.

        Returns:
            (distance, jimage): distance and periodic lattice translations
            of the other site for which the distance applies.
        """
        if jimage is None:
            #The following code is heavily vectorized to maximize speed.
            #Get the image adjustment necessary to bring coords to unit_cell.
            adj1 = np.floor(self._fcoords)
            adj2 = np.floor(fcoords)
            #Shift coords to unitcell
            coord1 = self._fcoords - adj1
            coord2 = fcoords - adj2
            # Generate set of images required for testing.
            # This is a cheat to create an 8x3 array of all length 3
            # combinations of 0,1
            test_set = np.unpackbits(np.array([5, 57, 119],
                                              dtype=np.uint8)).reshape(8, 3)
            images = np.copysign(test_set, coord1 - coord2)
            # Create tiled cartesian coords for computing distances.
            vec = np.tile(coord2 - coord1, (8, 1)) + images
            vec = self._lattice.get_cartesian_coords(vec)
            # Compute distances manually.
            dist = np.sqrt(np.sum(vec ** 2, 1)).tolist()
            # Return the minimum distance and the adjusted image corresponding
            # to the min distance.
            mindist = min(dist)
            ind = dist.index(mindist)
            return mindist, adj1 - adj2 + images[ind]

        mapped_vec = self._lattice.get_cartesian_coords(jimage + fcoords
                                                        - self._fcoords)
        return np.linalg.norm(mapped_vec), jimage

    def distance_and_image(self, other, jimage=None):
        """
        Gets distance and instance between two sites assuming periodic boundary
        conditions. If the index jimage of two sites atom j is not specified it
        selects the j image nearest to the i atom and returns the distance and
        jimage indices in terms of lattice vector translations. If the index
        jimage of atom j is specified it returns the distance between the ith
        atom and the specified jimage atom, the given jimage is also returned.

        Args:
            other (PeriodicSite): Other site to get distance from.
            jimage (3x1 array): Specific periodic image in terms of lattice
                translations, e.g., [1,0,0] implies to take periodic image
                that is one a-lattice vector away. If jimage == None,
                the image that is nearest to the site is found.

        Returns:
            (distance, jimage): distance and periodic lattice translations
            of the other site for which the distance applies.
        """
        return self.distance_and_image_from_frac_coords(other._fcoords, jimage)

    def distance(self, other, jimage=None):
        """
        Get distance between two sites assuming periodic boundary conditions.

        Args:
            other (PeriodicSite): Other site to get distance from.
            jimage (3x1 array): Specific periodic image in terms of lattice
                translations, e.g., [1,0,0] implies to take periodic image
                that is one a-lattice vector away. If jimage == None,
                the image that is nearest to the site is found.

        Returns:
            distance (float): Distance between the two sites
        """
        return self.distance_and_image(other, jimage)[0]

    def __repr__(self):
        return "PeriodicSite: {} ({:.4f}, {:.4f}, {:.4f}) [{:.4f}, {:.4f}, " \
               "{:.4f}]".format(self.species_string, self._coords[0],
                                self._coords[1], self._coords[2],
                                self._fcoords[0], self._fcoords[1],
                                self._fcoords[2])

    @property
    def to_dict(self):
        """
        Json-serializable dict representation of PeriodicSite.
        """
        species_list = []
        for spec, occu in self._species.items():
            d = spec.to_dict
            del d["@module"]
            del d["@class"]
            d["occu"] = occu
            species_list.append(d)
        return {"label": self.species_string, "species": species_list,
                "xyz": [float(c) for c in self._coords],
                "abc": [float(c) for c in self._fcoords],
                "lattice": self._lattice.to_dict,
                "properties": self._properties,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d, lattice=None):
        """
        Create PeriodicSite from dict representation.

        Args:
            d (dict): dict representation of PeriodicSite
            lattice: Optional lattice to override lattice specified in d.
                Useful for ensuring all sites in a structure share the same
                lattice.

        Returns:
            PeriodicSite
        """
        atoms_n_occu = {}
        for sp_occu in d["species"]:
            if "oxidation_state" in sp_occu and Element.is_valid_symbol(
                    sp_occu["element"]):
                sp = Specie.from_dict(sp_occu)
            elif "oxidation_state" in sp_occu:
                sp = DummySpecie.from_dict(sp_occu)
            else:
                sp = Element(sp_occu["element"])
            atoms_n_occu[sp] = sp_occu["occu"]
        props = d.get("properties", None)
        lattice = lattice if lattice else Lattice.from_dict(d["lattice"])
        return cls(atoms_n_occu, d["abc"], lattice, properties=props)
