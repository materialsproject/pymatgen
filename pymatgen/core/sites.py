# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines classes representing non-periodic and periodic sites.
"""

from __future__ import annotations

import collections
import json
from typing import cast

import numpy as np
from monty.json import MontyDecoder, MontyEncoder, MSONable

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import DummySpecies, Element, Species, get_el_sp
from pymatgen.util.coord import pbc_diff
from pymatgen.util.typing import ArrayLike, CompositionLike, SpeciesLike


class Site(collections.abc.Hashable, MSONable):
    """
    A generalized *non-periodic* site. This is essentially a composition
    at a point in space, with some optional properties associated with it. A
    Composition is used to represent the atoms and occupancy, which allows for
    disordered site representation. Coords are given in standard Cartesian
    coordinates.
    """

    position_atol = 1e-5

    def __init__(
        self,
        species: SpeciesLike | CompositionLike,
        coords: ArrayLike,
        properties: dict = None,
        skip_checks: bool = False,
    ):
        """
        Creates a non-periodic Site.

        :param species: Species on the site. Can be:
            i.  A Composition-type object (preferred)
            ii. An  element / species specified either as a string
                symbols, e.g. "Li", "Fe2+", "P" or atomic numbers,
                e.g., 3, 56, or actual Element or Species objects.
            iii.Dict of elements/species and occupancies, e.g.,
                {"Fe" : 0.5, "Mn":0.5}. This allows the setup of
                disordered structures.
        :param coords: Cartesian coordinates of site.
        :param properties: Properties associated with the site as a dict, e.g.
            {"magmom": 5}. Defaults to None.
        :param skip_checks: Whether to ignore all the usual checks and just
            create the site. Use this if the Site is created in a controlled
            manner and speed is desired.
        """
        if not skip_checks:
            if not isinstance(species, Composition):
                try:
                    species = Composition({get_el_sp(species): 1})
                except TypeError:
                    species = Composition(species)
            totaloccu = species.num_atoms
            if totaloccu > 1 + Composition.amount_tolerance:
                raise ValueError("Species occupancies sum to more than 1!")
            coords = np.array(coords)
        self._species: Composition = species  # type: ignore
        self.coords: np.ndarray = coords  # type: ignore
        self.properties: dict = properties or {}

    def __getattr__(self, a):
        # overriding getattr doesn't play nice with pickle, so we
        # can't use self._properties
        p = object.__getattribute__(self, "properties")
        if a in p:
            return p[a]
        raise AttributeError(a)

    @property
    def species(self) -> Composition:
        """
        :return: The species on the site as a composition, e.g., Fe0.5Mn0.5.
        """
        return self._species

    @species.setter
    def species(self, species: SpeciesLike | CompositionLike):
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
    def x(self) -> float:
        """
        Cartesian x coordinate
        """
        return self.coords[0]

    @x.setter
    def x(self, x: float):
        self.coords[0] = x

    @property
    def y(self) -> float:
        """
        Cartesian y coordinate
        """
        return self.coords[1]

    @y.setter
    def y(self, y: float):
        self.coords[1] = y

    @property
    def z(self) -> float:
        """
        Cartesian z coordinate
        """
        return self.coords[2]

    @z.setter
    def z(self, z: float):
        self.coords[2] = z

    def distance(self, other) -> float:
        """
        Get distance between two sites.

        Args:
            other: Other site.

        Returns:
            Distance (float)
        """
        return np.linalg.norm(other.coords - self.coords)

    def distance_from_point(self, pt) -> float:
        """
        Returns distance between the site and a point in space.

        Args:
            pt: Cartesian coordinates of point.

        Returns:
            Distance (float)
        """
        return np.linalg.norm(np.array(pt) - self.coords)

    @property
    def species_string(self) -> str:
        """
        String representation of species on the site.
        """
        if self.is_ordered:
            return str(list(self.species)[0])
        sorted_species = sorted(self.species)
        return ", ".join([f"{sp}:{self.species[sp]:.3f}" for sp in sorted_species])

    @property
    def specie(self) -> Element | Species | DummySpecies:
        """
        The Species/Element at the site. Only works for ordered sites. Otherwise
        an AttributeError is raised. Use this property sparingly.  Robust
        design should make use of the property species instead. Note that the
        singular of species is also species. So the choice of this variable
        name is governed by programmatic concerns as opposed to grammar.

        Raises:
            AttributeError if Site is not ordered.
        """
        if not self.is_ordered:
            raise AttributeError("specie property only works for ordered sites!")
        return list(self.species)[0]

    @property
    def is_ordered(self) -> bool:
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

    def __eq__(self, other: object) -> bool:
        """
        Site is equal to another site if the species and occupancies are the
        same, and the coordinates are the same to some tolerance.  numpy
        function `allclose` is used to determine if coordinates are close.
        """
        needed_attrs = ("species", "coords", "properties")
        if not all(hasattr(self, attr) for attr in needed_attrs):
            return NotImplemented
        other = cast(Site, other)
        return (
            self.species == other.species
            and np.allclose(self.coords, other.coords, atol=Site.position_atol)
            and self.properties == other.properties
        )

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        return sum(el.Z for el in self.species)

    def __contains__(self, el):
        return el in self.species

    def __repr__(self):
        return f"Site: {self.species_string} ({self.coords[0]:.4f}, {self.coords[1]:.4f}, {self.coords[2]:.4f})"

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas. For example, FeO4PLi is
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
        return f"{self.coords} {self.species_string}"

    def as_dict(self) -> dict:
        """
        JSON-serializable dict representation for Site.
        """
        species_list = []
        for spec, occu in self.species.items():
            d = spec.as_dict()
            del d["@module"]
            del d["@class"]
            d["occu"] = occu
            species_list.append(d)
        d = {
            "name": self.species_string,
            "species": species_list,
            "xyz": [float(c) for c in self.coords],
            "properties": self.properties,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }
        if self.properties:
            d["properties"] = self.properties
        return d

    @classmethod
    def from_dict(cls, d: dict) -> Site:
        """
        Create Site from dict representation
        """
        atoms_n_occu = {}
        for sp_occu in d["species"]:
            if "oxidation_state" in sp_occu and Element.is_valid_symbol(sp_occu["element"]):
                sp = Species.from_dict(sp_occu)
            elif "oxidation_state" in sp_occu:
                sp = DummySpecies.from_dict(sp_occu)
            else:
                sp = Element(sp_occu["element"])  # type: ignore
            atoms_n_occu[sp] = sp_occu["occu"]
        props = d.get("properties", None)
        if props is not None:
            for key in props:
                props[key] = json.loads(json.dumps(props[key], cls=MontyEncoder), cls=MontyDecoder)
        return cls(atoms_n_occu, d["xyz"], properties=props)


class PeriodicSite(Site, MSONable):
    """
    Extension of generic Site object to periodic systems.
    PeriodicSite includes a lattice system.
    """

    def __init__(
        self,
        species: SpeciesLike | CompositionLike,
        coords: ArrayLike,
        lattice: Lattice,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        properties: dict = None,
        skip_checks: bool = False,
    ):
        """
        Create a periodic site.

        :param species: Species on the site. Can be:
            i.  A Composition-type object (preferred)
            ii. An  element / species specified either as a string
                symbols, e.g. "Li", "Fe2+", "P" or atomic numbers,
                e.g., 3, 56, or actual Element or Species objects.
            iii.Dict of elements/species and occupancies, e.g.,
                {"Fe" : 0.5, "Mn":0.5}. This allows the setup of
                disordered structures.
        :param coords: Cartesian coordinates of site.
        :param lattice: Lattice associated with the site.
        :param to_unit_cell: Translates fractional coordinate to the
            basic unit cell, i.e. all fractional coordinates satisfy 0
            <= a < 1. Defaults to False.
        :param coords_are_cartesian: Set to True if you are providing
            Cartesian coordinates. Defaults to False.
        :param properties: Properties associated with the site as a dict, e.g.
            {"magmom": 5}. Defaults to None.
        :param skip_checks: Whether to ignore all the usual checks and just
            create the site. Use this if the PeriodicSite is created in a
            controlled manner and speed is desired.
        """

        if coords_are_cartesian:
            frac_coords = lattice.get_fractional_coords(coords)
        else:
            frac_coords = coords  # type: ignore

        if to_unit_cell:
            frac_coords = [np.mod(f, 1) if p else f for p, f in zip(lattice.pbc, frac_coords)]

        if not skip_checks:
            frac_coords = np.array(frac_coords)
            if not isinstance(species, Composition):
                try:
                    species = Composition({get_el_sp(species): 1})
                except TypeError:
                    species = Composition(species)

            totaloccu = species.num_atoms
            if totaloccu > 1 + Composition.amount_tolerance:
                raise ValueError("Species occupancies sum to more than 1!")

        self._lattice: Lattice = lattice
        self._frac_coords: ArrayLike = frac_coords
        self._species: Composition = species  # type: ignore
        self._coords: np.ndarray | None = None
        self.properties: dict = properties or {}

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        return sum(el.Z for el in self.species)

    @property
    def lattice(self) -> Lattice:
        """
        Lattice associated with PeriodicSite
        """
        return self._lattice

    @lattice.setter
    def lattice(self, lattice: Lattice):
        """
        Sets Lattice associated with PeriodicSite
        """
        self._lattice = lattice
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property  # type: ignore
    def coords(self) -> np.ndarray:  # type: ignore
        """
        Cartesian coordinates
        """
        if self._coords is None:
            self._coords = self._lattice.get_cartesian_coords(self._frac_coords)
        return self._coords

    @coords.setter
    def coords(self, coords):
        """
        Set Cartesian coordinates
        """
        self._coords = np.array(coords)
        self._frac_coords = self._lattice.get_fractional_coords(self._coords)

    @property
    def frac_coords(self) -> np.ndarray:
        """
        Fractional coordinates
        """
        return self._frac_coords  # type: ignore

    @frac_coords.setter
    def frac_coords(self, frac_coords):
        """
        Set fractional coordinates
        """
        self._frac_coords = np.array(frac_coords)
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def a(self) -> float:
        """
        Fractional a coordinate
        """
        return self._frac_coords[0]  # type: ignore

    @a.setter
    def a(self, a: float):
        self._frac_coords[0] = a  # type: ignore
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def b(self) -> float:
        """
        Fractional b coordinate
        """
        return self._frac_coords[1]  # type: ignore

    @b.setter
    def b(self, b: float):
        self._frac_coords[1] = b  # type: ignore
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def c(self) -> float:
        """
        Fractional c coordinate
        """
        return self._frac_coords[2]  # type: ignore

    @c.setter
    def c(self, c: float):
        self._frac_coords[2] = c  # type: ignore
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def x(self) -> float:
        """
        Cartesian x coordinate
        """
        return self.coords[0]

    @x.setter
    def x(self, x: float):
        self.coords[0] = x
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    @property
    def y(self) -> float:
        """
        Cartesian y coordinate
        """
        return self.coords[1]

    @y.setter
    def y(self, y: float):
        self.coords[1] = y
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    @property
    def z(self) -> float:
        """
        Cartesian z coordinate
        """
        return self.coords[2]

    @z.setter
    def z(self, z: float):
        self.coords[2] = z
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    def to_unit_cell(self, in_place=False) -> PeriodicSite | None:
        """
        Move frac coords to within the unit cell cell.
        """
        frac_coords = [np.mod(f, 1) if p else f for p, f in zip(self.lattice.pbc, self.frac_coords)]
        if in_place:
            self.frac_coords = np.array(frac_coords)
            return None
        return PeriodicSite(self.species, frac_coords, self.lattice, properties=self.properties)

    def is_periodic_image(self, other: PeriodicSite, tolerance: float = 1e-8, check_lattice: bool = True) -> bool:
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

        frac_diff = pbc_diff(self.frac_coords, other.frac_coords, self.lattice.pbc)
        return np.allclose(frac_diff, [0, 0, 0], atol=tolerance)

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("species", "lattice", "properties", "coords")
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented

        other = cast(PeriodicSite, other)

        return (
            self.species == other.species
            and self.lattice == other.lattice
            and np.allclose(self.coords, other.coords, atol=Site.position_atol)
            and self.properties == other.properties
        )

    def distance_and_image_from_frac_coords(
        self, fcoords: ArrayLike, jimage: ArrayLike | None = None
    ) -> tuple[float, np.ndarray]:
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
        return self.lattice.get_distance_and_image(self.frac_coords, fcoords, jimage=jimage)

    def distance_and_image(self, other: PeriodicSite, jimage: ArrayLike | None = None) -> tuple[float, np.ndarray]:
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

    def distance(self, other: PeriodicSite, jimage: ArrayLike | None = None):
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
        return (
            f"PeriodicSite: {self.species_string} "
            f"({self.coords[0]:.4f}, {self.coords[1]:.4f}, {self.coords[2]:.4f}) "
            f"[{self._frac_coords[0]:.4f}, {self._frac_coords[1]:.4f}, {self._frac_coords[2]:.4f}]"
        )

    def as_dict(self, verbosity: int = 0) -> dict:
        """
        JSON-serializable dict representation of PeriodicSite.

        Args:
            verbosity (int): Verbosity level. Default of 0 only includes the
                matrix representation. Set to 1 for more details such as
                Cartesian coordinates, etc.
        """
        species_list = []
        for spec, occu in self._species.items():
            d = spec.as_dict()
            del d["@module"]
            del d["@class"]
            d["occu"] = occu
            species_list.append(d)

        d = {
            "species": species_list,
            "abc": [float(c) for c in self._frac_coords],  # type: ignore
            "lattice": self._lattice.as_dict(verbosity=verbosity),
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

        if verbosity > 0:
            d["xyz"] = [float(c) for c in self.coords]
            d["label"] = self.species_string

        d["properties"] = self.properties

        return d

    @classmethod
    def from_dict(cls, d, lattice=None) -> PeriodicSite:
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
        species = {}
        for sp_occu in d["species"]:
            if "oxidation_state" in sp_occu and Element.is_valid_symbol(sp_occu["element"]):
                sp = Species.from_dict(sp_occu)
            elif "oxidation_state" in sp_occu:
                sp = DummySpecies.from_dict(sp_occu)
            else:
                sp = Element(sp_occu["element"])  # type: ignore
            species[sp] = sp_occu["occu"]
        props = d.get("properties", None)
        if props is not None:
            for key in props:
                props[key] = json.loads(json.dumps(props[key], cls=MontyEncoder), cls=MontyDecoder)
        lattice = lattice if lattice else Lattice.from_dict(d["lattice"])
        return cls(species, d["abc"], lattice, properties=props)
