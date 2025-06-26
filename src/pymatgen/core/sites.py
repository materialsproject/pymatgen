"""This module defines classes representing non-periodic and periodic sites."""

from __future__ import annotations

import collections
import json
import warnings
from typing import TYPE_CHECKING, cast

import numpy as np
from monty.json import MontyDecoder, MontyEncoder, MSONable

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import DummySpecies, Element, Species, get_el_sp
from pymatgen.util.coord import pbc_diff
from pymatgen.util.misc import is_np_dict_equal

if TYPE_CHECKING:
    from typing import Any, Literal

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.util.typing import CompositionLike, SpeciesLike


class Site(collections.abc.Hashable, MSONable):
    """A generalized *non-periodic* site. This is essentially a composition
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
        properties: dict | None = None,
        label: str | None = None,
        skip_checks: bool = False,
    ) -> None:
        """Create a non-periodic Site.

        Args:
            species: Species on the site. Can be:
                i.  A Composition-type object (preferred)
                ii. An element / species specified either as a string
                    symbols, e.g. "Li", "Fe2+", "P" or atomic numbers,
                    e.g. 3, 56, or actual Element or Species objects.
                iii.Dict of elements/species and occupancies, e.g.
                    {"Fe" : 0.5, "Mn":0.5}. This allows the setup of
                    disordered structures.
            coords (NDArray): Cartesian coordinates of site.
            properties (dict): Properties associated with the site, e.g.
                {"magmom": 5}. Defaults to None.
            label (str): Label for the site. Defaults to None.
            skip_checks (bool): Whether to ignore all the usual checks and just
                create the site. Use this if the Site is created in a controlled
                manner and speed is desired.
        """
        if not skip_checks:
            if not isinstance(species, Composition):
                try:
                    species = Composition({get_el_sp(species): 1})  # type: ignore[arg-type]
                except TypeError:
                    species = Composition(species)
            total_occu = species.num_atoms
            if total_occu > 1 + Composition.amount_tolerance:
                raise ValueError("Species occupancies sum to more than 1!")
            coords = np.array(coords)

        self._species = species
        self.coords: NDArray[np.float64] = np.asarray(coords, dtype=np.float64)
        self.properties: dict = properties or {}
        self._label = label

    def __getattr__(self, attr: str) -> Any:
        # Override getattr doesn't play nicely with pickle,
        # so we can't use self._properties
        props = self.__getattribute__("properties")
        if attr in props:
            return props[attr]
        raise AttributeError(f"{attr=} not found on {type(self).__name__}")

    def __getitem__(self, el: Element) -> float:
        """Get the occupancy for element."""
        return self.species[el]

    def __eq__(self, other: object) -> bool:
        """Site is equal to another site if the species and occupancies are the
        same, and the coordinates are the same to some tolerance. `np.allclose`
        is used to determine if coordinates are close.
        """
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self.species == other.species
            and np.allclose(self.coords, other.coords, atol=type(self).position_atol)
            and is_np_dict_equal(self.properties, other.properties)
        )

    def __hash__(self) -> int:
        """Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        return sum(el.Z for el in self.species)

    def __contains__(self, el: Element) -> bool:
        return el in self.species

    def __repr__(self) -> str:
        name = self.species_string

        if self.label != name:
            name = f"{self.label} ({name})"

        return f"Site: {name} ({self.coords[0]:.4f}, {self.coords[1]:.4f}, {self.coords[2]:.4f})"

    def __lt__(self, other: Self) -> bool:
        """Set a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas. For example, FeO4PLi is
        automatically sorted in LiFePO4.
        """
        if self.species.average_electroneg < other.species.average_electroneg:
            return True
        if self.species.average_electroneg > other.species.average_electroneg:
            return False
        return self.species_string < other.species_string

    def __str__(self) -> str:
        return f"{self.coords} {self.species_string}"

    @property
    def species(self) -> Composition:
        """The species on the site as a composition, e.g. Fe0.5Mn0.5."""
        return cast("Composition", self._species)

    @species.setter
    def species(self, species: SpeciesLike | CompositionLike) -> None:
        if not isinstance(species, Composition):
            try:
                species = Composition({get_el_sp(species): 1})  # type: ignore[arg-type]
            except TypeError:
                species = Composition(species)

        total_occu = species.num_atoms
        if total_occu > 1 + Composition.amount_tolerance:
            raise ValueError("Species occupancies sum to more than 1!")

        self._species = cast("Composition", species)

    @property
    def label(self) -> str:
        """Site label."""
        return self._label if self._label is not None else self.species_string

    @label.setter
    def label(self, label: str | None) -> None:
        self._label = label

    @property
    def x(self) -> float:
        """Cartesian x coordinate."""
        return self.coords[0]

    @x.setter
    def x(self, x: float) -> None:
        self.coords[0] = x

    @property
    def y(self) -> float:
        """Cartesian y coordinate."""
        return self.coords[1]

    @y.setter
    def y(self, y: float) -> None:
        self.coords[1] = y

    @property
    def z(self) -> float:
        """Cartesian z coordinate."""
        return self.coords[2]

    @z.setter
    def z(self, z: float) -> None:
        self.coords[2] = z

    def distance(self, other: Site) -> float:
        """Get distance between two sites.

        Args:
            other: Other site.

        Returns:
            float: distance
        """
        return float(np.linalg.norm(other.coords - self.coords))

    def distance_from_point(self, pt: ArrayLike) -> float:
        """Get distance between the site and a point in space.

        Args:
            pt: Cartesian coordinates of point.

        Returns:
            float: distance
        """
        return float(np.linalg.norm(np.array(pt) - self.coords))

    @property
    def species_string(self) -> str:
        """String representation of species on the site."""
        if self.is_ordered:
            return str(next(iter(self.species)))
        return ", ".join(f"{sp}:{self.species[sp]:.3}" for sp in sorted(self.species))

    @property
    def specie(self) -> Element | Species | DummySpecies:
        """The Species/Element at the site. Only works for ordered sites. Otherwise
        an AttributeError is raised. Use this property sparingly. Robust
        design should make use of the property species instead. Note that the
        singular of species is also species. So the choice of this variable
        name is governed by programmatic concerns as opposed to grammar.

        Raises:
            AttributeError if Site is not ordered.
        """
        if not self.is_ordered:
            raise AttributeError("specie property only works for ordered sites!")
        return next(iter(self.species))

    @property
    def is_ordered(self) -> bool:
        """True if site is an ordered site, i.e., with a single species with
        occupancy 1.
        """
        total_occu = self.species.num_atoms
        return total_occu == len(self.species) == 1

    def as_dict(self) -> dict:
        """JSON-serializable dict representation for Site."""
        species = []
        for spec, occu in self.species.items():
            spec_dct = spec.as_dict()
            del spec_dct["@module"]
            del spec_dct["@class"]
            spec_dct["occu"] = occu
            species.append(spec_dct)

        dct = {
            "name": self.species_string,
            "species": species,
            "xyz": self.coords.astype(float).tolist(),
            "properties": self.properties,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "label": self.label,
        }
        if self.properties:
            dct["properties"] = self.properties

        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create Site from dict representation."""
        atoms_n_occu = {}
        for sp_occu in dct["species"]:
            if "oxidation_state" in sp_occu and Element.is_valid_symbol(sp_occu["element"]):
                sp: Species | DummySpecies | Element = Species.from_dict(sp_occu)
            elif "oxidation_state" in sp_occu:
                sp = DummySpecies.from_dict(sp_occu)
            else:
                sp = Element(sp_occu["element"])
            atoms_n_occu[sp] = sp_occu["occu"]
        props = dct.get("properties")
        if props is not None:
            for key in props:
                props[key] = json.loads(json.dumps(props[key], cls=MontyEncoder), cls=MontyDecoder)
        label = dct.get("label")
        return cls(atoms_n_occu, dct["xyz"], properties=props, label=label)


class PeriodicSite(Site, MSONable):
    """Extension of generic Site object to periodic systems.
    PeriodicSite includes a lattice system.
    """

    def __init__(
        self,
        species: SpeciesLike | CompositionLike,
        coords: ArrayLike,
        lattice: Lattice,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        properties: dict | None = None,
        label: str | None = None,
        skip_checks: bool = False,
    ) -> None:
        """Create a periodic site.

        Args:
            species: Species on the site. Can be:
                i.  A Composition-type object (preferred)
                ii. An element / species specified either as a string
                    symbols, e.g. "Li", "Fe2+", "P" or atomic numbers,
                    e.g. 3, 56, or actual Element or Species objects.
                iii.Dict of elements/species and occupancies, e.g.
                    {"Fe": 0.5, "Mn": 0.5}. This allows the setup of
                    disordered structures.
            coords (ArrayLike): Coordinates of site, fractional coordinates
                by default. See ``coords_are_cartesian`` for more details.
            lattice (Lattice): Lattice associated with the site.
            to_unit_cell (bool): Translates fractional coordinate to the
                basic unit cell, i.e. all fractional coordinates satisfy 0
                <= a < 1. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                Cartesian coordinates. Defaults to False.
            properties (dict): Properties associated with the site, e.g.
                {"magmom": 5}. Defaults to None.
            label (str): Label for the site. Defaults to None.
            skip_checks (bool): Whether to ignore all the usual checks and just
                create the site. Use this if the PeriodicSite is created in a
                controlled manner and speed is desired.
        """
        frac_coords = lattice.get_fractional_coords(coords) if coords_are_cartesian else coords

        if to_unit_cell:
            frac_coords = np.array([np.mod(f, 1) if p else f for p, f in zip(lattice.pbc, frac_coords, strict=True)])  # type: ignore[arg-type]

        if not skip_checks:
            frac_coords = np.array(frac_coords)
            if not isinstance(species, Composition):
                try:
                    species = Composition({get_el_sp(species): 1})  # type: ignore[arg-type]
                except TypeError:
                    species = Composition(species)

            total_occu = species.num_atoms
            if total_occu > 1 + Composition.amount_tolerance:
                raise ValueError("Species occupancies sum to more than 1!")

        self._lattice: Lattice = lattice
        self._frac_coords: NDArray[np.float64] = np.asarray(frac_coords, dtype=np.float64)
        self._species: Composition = cast("Composition", species)
        self._coords: NDArray[np.float64] | None = None
        self.properties: dict = properties or {}
        self._label = label

    def __hash__(self) -> int:
        """Minimally effective hash function that just distinguishes between Sites
        with different elements.
        """
        return sum(el.Z for el in self.species)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented

        return (
            self.species == other.species
            and self.lattice == other.lattice
            and np.allclose(self.coords, other.coords, atol=Site.position_atol)
            and is_np_dict_equal(self.properties, other.properties)
        )

    def __repr__(self) -> str:
        name = self.species_string

        if self.label != name:
            name = f"{self.label} ({name})"

        x, y, z = self.coords
        x_frac, y_frac, z_frac = map(float, self.frac_coords)
        cls_name = type(self).__name__
        return f"{cls_name}: {name} ({x:.4}, {y:.4}, {z:.4}) [{x_frac:.4}, {y_frac:.4}, {z_frac:.4}]"

    @property
    def lattice(self) -> Lattice:
        """Lattice associated with PeriodicSite."""
        return self._lattice

    @lattice.setter
    def lattice(self, lattice: Lattice) -> None:
        """Set Lattice associated with PeriodicSite."""
        self._lattice = lattice
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def coords(self) -> NDArray[np.float64]:
        """Cartesian coordinates."""
        if self._coords is None:
            self._coords = self._lattice.get_cartesian_coords(self._frac_coords)
        return self._coords

    @coords.setter
    def coords(self, coords: ArrayLike) -> None:
        """Set Cartesian coordinates."""
        self._coords = np.asarray(coords, dtype=np.float64)
        self._frac_coords = self._lattice.get_fractional_coords(self._coords)

    @property
    def frac_coords(self) -> NDArray[np.float64]:
        """Fractional coordinates."""
        return self._frac_coords

    @frac_coords.setter
    def frac_coords(self, frac_coords: ArrayLike) -> None:
        """Set fractional coordinates."""
        self._frac_coords = np.array(frac_coords, dtype=np.float64)
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def a(self) -> float:
        """Fractional a coordinate."""
        return self._frac_coords[0]

    @a.setter
    def a(self, a: float) -> None:
        self._frac_coords[0] = a
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def b(self) -> float:
        """Fractional b coordinate."""
        return self._frac_coords[1]

    @b.setter
    def b(self, b: float) -> None:
        self._frac_coords[1] = b
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def c(self) -> float:
        """Fractional c coordinate."""
        return self._frac_coords[2]

    @c.setter
    def c(self, c: float) -> None:
        self._frac_coords[2] = c
        self._coords = self._lattice.get_cartesian_coords(self._frac_coords)

    @property
    def x(self) -> float:
        """Cartesian x coordinate."""
        return self.coords[0]

    @x.setter
    def x(self, x: float) -> None:
        self.coords[0] = x
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    @property
    def y(self) -> float:
        """Cartesian y coordinate."""
        return self.coords[1]

    @y.setter
    def y(self, y: float) -> None:
        self.coords[1] = y
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    @property
    def z(self) -> float:
        """Cartesian z coordinate."""
        return self.coords[2]

    @z.setter
    def z(self, z: float) -> None:
        self.coords[2] = z
        self._frac_coords = self._lattice.get_fractional_coords(self.coords)

    def to_unit_cell(self, in_place: bool = False) -> Self | None:
        """Move frac coords to within the unit cell."""
        frac_coords = np.array(
            [np.mod(f, 1) if p else f for p, f in zip(self.lattice.pbc, self.frac_coords, strict=True)]
        )
        if in_place:
            self.frac_coords = frac_coords
            return None
        return type(self)(
            self.species,
            frac_coords,
            self.lattice,
            properties=self.properties,
            label=self.label,
        )

    def is_periodic_image(
        self,
        other: Self,
        tolerance: float = 1e-8,
        check_lattice: bool = True,
    ) -> bool:
        """Check if sites are periodic images of each other.

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
        return all(abs(diff) <= tolerance for diff in frac_diff)

    def distance_and_image_from_frac_coords(
        self,
        fcoords: ArrayLike,
        jimage: ArrayLike | None = None,
    ) -> tuple[float, NDArray[np.int_]]:
        """Get distance between site and a fractional coordinate assuming
        periodic boundary conditions. If the index jimage of two sites atom j
        is not specified it selects the j image nearest to the i atom and
        returns the distance and jimage indices in terms of lattice vector
        translations. If the index jimage of atom j is specified it returns the
        distance between the i atom and the specified jimage atom, the given
        jimage is also returned.

        Args:
            fcoords (3x1 array): fractional coordinates to get distance from.
            jimage (3x1 array): Specific periodic image in terms of
                lattice translations, e.g. [1, 0, 0] implies to take periodic
                image that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            tuple[float, np.ndarray]: distance and periodic lattice translations (jimage)
                of the other site for which the distance applies.
        """
        return self.lattice.get_distance_and_image(self.frac_coords, fcoords, jimage=jimage)

    def distance_and_image(
        self,
        other: Self,
        jimage: ArrayLike | None = None,
    ) -> tuple[float, NDArray[np.int_]]:
        """Get distance and instance between two sites assuming periodic boundary
        conditions. If the index jimage of two sites atom j is not specified it
        selects the j image nearest to the i atom and returns the distance and
        jimage indices in terms of lattice vector translations. If the index
        jimage of atom j is specified it returns the distance between the ith
        atom and the specified jimage atom, the given jimage is also returned.

        Args:
            other (PeriodicSite): Other site to get distance from.
            jimage (3x1 array): Specific periodic image in terms of lattice
                translations, e.g. [1, 0, 0] implies to take periodic image
                that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            tuple[float, np.ndarray]: distance and periodic lattice translations (jimage)
                of the other site for which the distance applies.
        """
        return self.distance_and_image_from_frac_coords(other.frac_coords, jimage)

    def distance(
        self,
        other: Self,
        jimage: ArrayLike | None = None,
    ) -> float:
        """Get distance between two sites assuming periodic boundary conditions.

        Args:
            other (PeriodicSite): Other site to get distance from.
            jimage (3x1 array): Specific periodic image in terms of lattice
                translations, e.g. [1, 0, 0] implies to take periodic image
                that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            float: distance between the two sites.
        """
        return self.distance_and_image(other, jimage)[0]

    def as_dict(self, verbosity: Literal[0, 1] = 0) -> dict:
        """JSON-serializable dict representation of PeriodicSite.

        Args:
            verbosity (0 | 1): Verbosity level. Default of 0 only includes the matrix
                representation. Set to 1 for more details such as Cartesian coordinates, etc.
        """
        species = []
        for spec, occu in self._species.items():
            dct = spec.as_dict()
            del dct["@module"]
            del dct["@class"]
            dct["occu"] = occu
            species.append(dct)

        dct = {
            "species": species,
            "abc": self._frac_coords.astype(float).tolist(),
            "lattice": self._lattice.as_dict(verbosity=verbosity),
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "properties": self.properties,
            "label": self.label,
        }

        if verbosity not in {0, 1}:
            warnings.warn(
                f"`verbosity={verbosity}` is deprecated and will be disallowed in a future version. "
                "Please use 0 (silent) or 1 (verbose) explicitly.",
                DeprecationWarning,
                stacklevel=2,
            )
        if verbosity > 0:  # TODO: explicitly check `verbosity == 1`
            dct["xyz"] = self.coords.astype(float).tolist()

        return dct

    @classmethod
    def from_dict(cls, dct: dict, lattice: Lattice | None = None) -> Self:
        """Create PeriodicSite from dict representation.

        Args:
            dct (dict): dict representation of PeriodicSite
            lattice: Optional lattice to override lattice specified in d.
                Useful for ensuring all sites in a structure share the same
                lattice.

        Returns:
            PeriodicSite
        """
        species = {}
        for sp_occu in dct["species"]:
            if "oxidation_state" in sp_occu and Element.is_valid_symbol(sp_occu["element"]):
                sp: Species | DummySpecies | Element = Species.from_dict(sp_occu)
            elif "oxidation_state" in sp_occu:
                sp = DummySpecies.from_dict(sp_occu)
            else:
                sp = Element(sp_occu["element"])
            species[sp] = sp_occu["occu"]

        props = dct.get("properties")
        if props is not None:
            for key in props:
                props[key] = json.loads(json.dumps(props[key], cls=MontyEncoder), cls=MontyDecoder)
        label = dct.get("label")
        lattice = lattice or Lattice.from_dict(dct["lattice"])
        return cls(species, dct["abc"], lattice, properties=props, label=label)
