# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import collections
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, Specie, DummySpecie,\
    get_el_sp
from monty.json import MSONable
from pymatgen.util.coord import pbc_diff
from pymatgen.core.composition import Composition
from monty.dev import deprecated

"""
This module defines classes representing non-periodic and periodic sites.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 17, 2012"


class Site(collections.abc.Hashable, MSONable):
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
                i.  A Composition object (preferred)
                ii. An  element / specie specified either as a string
                    symbols, e.g. "Li", "Fe2+", "P" or atomic numbers,
                    e.g., 3, 56, or actual Element or Specie objects.
                iii.Dict of elements/species and occupancies, e.g.,
                    {"Fe" : 0.5, "Mn":0.5}. This allows the setup of
                    disordered structures.
            coords: Cartesian coordinates of site.
            properties: Properties associated with the site as a dict, e.g.
                {"magmom": 5}. Defaults to None.
        """
        if isinstance(atoms_n_occu, Composition):
            # Compositions are immutable, so don't need to copy (much faster)
            species = atoms_n_occu
        else:
            try:
                species = Composition({get_el_sp(atoms_n_occu): 1})
            except TypeError:
                species = Composition(atoms_n_occu)
        totaloccu = species.num_atoms
        if totaloccu > 1 + Composition.amount_tolerance:
            raise ValueError("Species occupancies sum to more than 1!")
        self._species = species
        self.coords = np.array(coords)
        self.properties = properties or {}

    def __getattr__(self, a):
        # overriding getattr doens't play nice with pickle, so we
        # can't use self._properties
        p = object.__getattribute__(self, 'properties')
        if a in p:
            return p[a]
        raise AttributeError(a)

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, species):
        if not isinstance(species, Composition):
            try:
                species = Composition({get_el_sp(species): 1})
            except TypeError:
                species = Composition(species)
        totaloccu = species.num_atoms
        if totaloccu > 1 + Composition.amount_tolerance:
            raise ValueError("Species occupancies sum to more than 1!")
        self._species = species

    @property
    def x(self):
        """
        Cartesian x coordinate
        """
        return self.coords[0]

    @x.setter
    def x(self, x):
        self.coords[0] = x

    @property
    def y(self):
        """
        Cartesian y coordinate
        """
        return self.coords[1]

    @y.setter
    def y(self, y):
        self.coords[1] = y

    @property
    def z(self):
        """
        Cartesian z coordinate
        """
        return self.coords[2]

    @z.setter
    def z(self, z):
        self.coords[2] = z

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
        return np.linalg.norm(np.array(pt) - self.coords)

    @property
    def species_string(self):
        """
        String representation of species on the site.
        """
        if self.is_ordered:
            return list(self.species.keys())[0].__str__()
        else:
            sorted_species = sorted(self.species.keys())
            return ", ".join(["{}:{:.3f}".format(sp, self.species[sp])
                              for sp in sorted_species])

    @property
    @deprecated(message="Use site.species instead. This will be deprecated with effect from pymatgen 2020.")
    def species_and_occu(self):
        """
        The species at the site, i.e., a Composition mapping type of
        element/species to occupancy.
        """
        return self.species

    @property
    def specie(self):
        """
        The Specie/Element at the site. Only works for ordered sites. Otherwise
        an AttributeError is raised. Use this property sparingly.  Robust
        design should make use of the property species_and_occu instead.

        Raises:
            AttributeError if Site is not ordered.
        """
        if not self.is_ordered:
            raise AttributeError("specie property only works for ordered "
                                 "sites!")
        return list(self.species.keys())[0]

    @property
    def is_ordered(self):
        """
        True if site is an ordered site, i.e., with a single species with
        occupancy 1.
        """
        totaloccu = self.species.num_atoms
        return totaloccu == 1 and len(self.species) == 1

    def __getitem__(self, el):
        """
        Get the occupancy for element
        """
        return self.species[el]

    def __eq__(self, other):
        """
        Site is equal to another site if the species and occupancies are the
        same, and the coordinates are the same to some tolerance.  numpy
        function `allclose` is used to determine if coordinates are close.
        """
        if other is None:
            return False
        return (self.species == other.species and
                np.allclose(self.coords, other.coords,
                            atol=Site.position_atol) and
                self.properties == other.properties)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        return sum([el.Z for el in self.species.keys()])

    def __contains__(self, el):
        return el in self.species

    def __repr__(self):
        return "Site: {} ({:.4f}, {:.4f}, {:.4f})".format(
            self.species_string, *self.coords)

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas.  For example, FeO4PLi is
        automatically sorted in LiFePO4.
        """
        if self.species.average_electroneg < other.species.average_electroneg:
            return True
        if self.species.average_electroneg > other.species.average_electroneg:
            return False
        if self.species_string < other.species_string:
            return True
        if self.species_string > other.species_string:
            return False
        return False

    def __str__(self):
        return "{} {}".format(self.coords, self.species_string)

    def as_dict(self):
        """
        Json-serializable dict representation for Site.
        """
        species_list = []
        for spec, occu in self.species.items():
            d = spec.as_dict()
            del d["@module"]
            del d["@class"]
            d["occu"] = occu
            species_list.append(d)
        d = {"name": self.species_string, "species": species_list,
             "xyz": [float(c) for c in self.coords],
             "properties": self.properties,
             "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        if self.properties:
            d["properties"] = self.properties
        return d

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
        if coords_are_cartesian:
            frac_coords = lattice.get_fractional_coords(coords)
            cart_coords = coords
        else:
            frac_coords = np.array(coords)
            cart_coords = lattice.get_cartesian_coords(coords)

        if to_unit_cell:
            frac_coords = np.mod(frac_coords, 1)
            cart_coords = lattice.get_cartesian_coords(frac_coords)

        self._lattice = lattice
        self._frac_coords = frac_coords
        super(PeriodicSite, self).__init__(atoms_n_occu, cart_coords, properties)

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        return sum([el.Z for el in self.species.keys()])

    @property
    def lattice(self):
        """
        Lattice associated with PeriodicSite
        """
        return self._lattice

    @lattice.setter
    def lattice(self, lattice):
        """
        Sets Lattice associated with PeriodicSite
        """
        self._lattice = lattice
        self.coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def frac_coords(self):
        """
        Fractional coordinates
        """
        return self._frac_coords

    @frac_coords.setter
    def frac_coords(self, frac_coords):
        """
        Fractional a coordinate
        """
        self._frac_coords = np.array(frac_coords)
        self.coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def a(self):
        """
        Fractional a coordinate
        """
        return self._frac_coords[0]

    @a.setter
    def a(self, a):
        self._frac_coords[0] = a
        self.coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def b(self):
        """
        Fractional b coordinate
        """
        return self._frac_coords[1]

    @b.setter
    def b(self, b):
        self._frac_coords[1] = b
        self.coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def c(self):
        """
        Fractional c coordinate
        """
        return self._frac_coords[2]

    @c.setter
    def c(self, c):
        self._frac_coords[2] = c
        self.coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def x(self):
        """
        Cartesian x coordinate
        """
        return self.coords[0]

    @x.setter
    def x(self, x):
        self.coords[0] = x
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    @property
    def y(self):
        """
        Cartesian y coordinate
        """
        return self.coords[1]

    @y.setter
    def y(self, y):
        self.coords[1] = y
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    @property
    def z(self):
        """
        Cartesian z coordinate
        """
        return self.coords[2]

    @z.setter
    def z(self, z):
        self.coords[2] = z
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    def to_unit_cell(self, in_place=False):
        """
        Move frac coords to within the unit cell cell.
        """
        frac_coords = np.mod(self.frac_coords, 1)
        if in_place:
            self.frac_coords = frac_coords
        else:
            return PeriodicSite(self.species, frac_coords, self.lattice,
                                properties=self.properties)

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
        if check_lattice and self.lattice != other.lattice:
            return False
        if self.species != other.species:
            return False

        frac_diff = pbc_diff(self.frac_coords, other.frac_coords)
        return np.allclose(frac_diff, [0, 0, 0], atol=tolerance)

    def __eq__(self, other):
        return self.species == other.species and \
               self.lattice == other.lattice and \
               np.allclose(self.coords, other.coords,
                           atol=Site.position_atol) and \
               self.properties == other.properties

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
                image that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            (distance, jimage): distance and periodic lattice translations
            of the other site for which the distance applies.
        """
        return self.lattice.get_distance_and_image(self.frac_coords, fcoords,
                                                   jimage=jimage)

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
                that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            (distance, jimage): distance and periodic lattice translations
            of the other site for which the distance applies.
        """
        return self.distance_and_image_from_frac_coords(other.frac_coords, jimage)

    def distance(self, other, jimage=None):
        """
        Get distance between two sites assuming periodic boundary conditions.

        Args:
            other (PeriodicSite): Other site to get distance from.
            jimage (3x1 array): Specific periodic image in terms of lattice
                translations, e.g., [1,0,0] implies to take periodic image
                that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            distance (float): Distance between the two sites
        """
        return self.distance_and_image(other, jimage)[0]

    def __repr__(self):
        return "PeriodicSite: {} ({:.4f}, {:.4f}, {:.4f}) [{:.4f}, {:.4f}, " \
               "{:.4f}]".format(self.species_string, self.coords[0],
                                self.coords[1], self.coords[2],
                                *self._frac_coords)

    def as_dict(self, verbosity=0):
        """
        Json-serializable dict representation of PeriodicSite.

        Args:
            verbosity (int): Verbosity level. Default of 0 only includes the
                matrix representation. Set to 1 for more details such as
                cartesian coordinates, etc.
        """
        species_list = []
        for spec, occu in self.species.items():
            d = spec.as_dict()
            del d["@module"]
            del d["@class"]
            d["occu"] = occu
            species_list.append(d)

        d = {"species": species_list,
             "abc": [float(c) for c in self.frac_coords],
             "lattice": self.lattice.as_dict(verbosity=verbosity),
             "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}

        if verbosity > 0:
            d["xyz"] = [float(c) for c in self.coords]
            d["label"] = self.species_string

        d["properties"] = self.properties
        return d

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
