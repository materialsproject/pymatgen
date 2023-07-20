"""
This module provides classes used to define a non-periodic molecule and a
periodic structure.
"""

from __future__ import annotations

import collections
import contextlib
import functools
import inspect
import io
import itertools
import json
import math
import os
import random
import re
import sys
import warnings
from abc import ABCMeta, abstractmethod
from fnmatch import fnmatch
from inspect import isclass
from io import StringIO
from typing import TYPE_CHECKING, Any, Callable, Iterable, Iterator, Literal, Sequence, SupportsIndex, cast

import numpy as np
from monty.dev import deprecated
from monty.io import zopen
from monty.json import MSONable
from ruamel.yaml import YAML
from tabulate import tabulate

from pymatgen.core.bonds import CovalentBond, get_bond_length
from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice, get_points_in_spheres
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import DummySpecies, Element, Species, get_el_sp
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.units import Length, Mass
from pymatgen.electronic_structure.core import Magmom
from pymatgen.symmetry.maggroups import MagneticSpaceGroup
from pymatgen.util.coord import all_distances, get_angle, lattice_points_in_supercell

if TYPE_CHECKING:
    from ase.calculators.calculator import Calculator
    from ase.io.trajectory import Trajectory
    from ase.optimize.optimize import Optimizer
    from matgl.ext.ase import TrajectoryObserver
    from numpy.typing import ArrayLike

    from pymatgen.util.typing import CompositionLike, SpeciesLike


class Neighbor(Site):
    """
    Simple Site subclass to contain a neighboring atom that skips all the unnecessary checks for speed. Can be
    used as a fixed-length tuple of size 3 to retain backwards compatibility with past use cases.

        (site, nn_distance, index).

    In future, usage should be to call attributes, e.g., Neighbor.index, Neighbor.distance, etc.
    """

    def __init__(
        self,
        species: Composition,
        coords: np.ndarray,
        properties: dict | None = None,
        nn_distance: float = 0.0,
        index: int = 0,
        label: str | None = None,
    ):
        """
        :param species: Same as Site
        :param coords: Same as Site, but must be fractional.
        :param properties: Same as Site
        :param nn_distance: Distance to some other Site.
        :param index: Index within structure.
        :param label: Label for the site. Defaults to None.
        """
        self.coords = coords
        self._species = species
        self.properties = properties or {}
        self.nn_distance = nn_distance
        self.index = index
        self.label = label if label is not None else self.species_string

    def __len__(self) -> Literal[3]:
        """Make neighbor Tuple-like to retain backwards compatibility."""
        return 3

    def __getitem__(self, idx: int):
        """Make neighbor Tuple-like to retain backwards compatibility."""
        return (self, self.nn_distance, self.index)[idx]

    def as_dict(self) -> dict:  # type: ignore
        """
        Note that method calls the super of Site, which is MSONable itself.

        Returns: dict
        """
        return super(Site, self).as_dict()  # pylint: disable=E1003

    @classmethod
    def from_dict(cls, d: dict) -> Neighbor:  # type: ignore
        """
        Returns a Neighbor from a dict.

        Args:
            d: MSONable dict format.

        Returns:
            Neighbor
        """
        return super(Site, cls).from_dict(d)  # pylint: disable=E1003


class PeriodicNeighbor(PeriodicSite):
    """
    Simple PeriodicSite subclass to contain a neighboring atom that skips all
    the unnecessary checks for speed. Can be used as a fixed-length tuple of
    size 4 to retain backwards compatibility with past use cases.

        (site, distance, index, image).

    In future, usage should be to call attributes, e.g., PeriodicNeighbor.index,
    PeriodicNeighbor.distance, etc.
    """

    def __init__(
        self,
        species: Composition,
        coords: np.ndarray,
        lattice: Lattice,
        properties: dict | None = None,
        nn_distance: float = 0.0,
        index: int = 0,
        image: tuple = (0, 0, 0),
        label: str | None = None,
    ):
        """
        Args:
            species (Composition): Same as PeriodicSite
            coords (np.ndarray): Same as PeriodicSite, but must be fractional.
            lattice (Lattice): Same as PeriodicSite
            properties (dict, optional): Same as PeriodicSite. Defaults to None.
            nn_distance (float, optional): Distance to some other Site.. Defaults to 0.0.
            index (int, optional): Index within structure.. Defaults to 0.
            image (tuple, optional): PeriodicImage. Defaults to (0, 0, 0).
            label (str, optional): Label for the site. Defaults to None.
        """
        self._lattice = lattice
        self._frac_coords = coords
        self._species = species
        self.properties = properties or {}
        self.nn_distance = nn_distance
        self.index = index
        self.image = image
        self.label = label if label is not None else self.species_string

    @property  # type: ignore
    def coords(self) -> np.ndarray:  # type: ignore
        """:return: Cartesian coords."""
        return self._lattice.get_cartesian_coords(self._frac_coords)

    def __len__(self):
        """Make neighbor Tuple-like to retain backwards compatibility."""
        return 4

    def __getitem__(self, idx: int | slice):
        """Make neighbor Tuple-like to retain backwards compatibility."""
        return (self, self.nn_distance, self.index, self.image)[idx]

    def as_dict(self) -> dict:  # type: ignore
        """
        Note that method calls the super of Site, which is MSONable itself.

        Returns: dict
        """
        return super(Site, self).as_dict()  # pylint: disable=E1003

    @classmethod
    def from_dict(cls, d: dict) -> PeriodicNeighbor:  # type: ignore
        """
        Returns a PeriodicNeighbor from a dict.

        Args:
            d: MSONable dict format.

        Returns:
            PeriodicNeighbor
        """
        return super(Site, cls).from_dict(d)  # pylint: disable=E1003


class SiteCollection(collections.abc.Sequence, metaclass=ABCMeta):
    """
    Basic SiteCollection. Essentially a sequence of Sites or PeriodicSites.
    This serves as a base class for Molecule (a collection of Site, i.e., no
    periodicity) and Structure (a collection of PeriodicSites, i.e.,
    periodicity). Not meant to be instantiated directly.
    """

    # Tolerance in Angstrom for determining if sites are too close.
    DISTANCE_TOLERANCE = 0.5

    @property
    @abstractmethod
    def sites(self) -> tuple[Site, ...]:
        """Returns a tuple of sites."""

    @abstractmethod
    def get_distance(self, i: int, j: int) -> float:
        """
        Returns distance between sites at index i and j.

        Args:
            i: Index of first site
            j: Index of second site

        Returns:
            Distance between sites at index i and index j.
        """

    @property
    def distance_matrix(self) -> np.ndarray:
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this is overwritten to return the nearest image
        distance.
        """
        return all_distances(self.cart_coords, self.cart_coords)

    @property
    def species(self) -> list[Element | Species]:
        """
        Only works for ordered structures.

        Raises:
            AttributeError: If structure is disordered.

        Returns:
            ([Species]) List of species at each site of the structure.
        """
        if not self.is_ordered:
            raise AttributeError("species property only supports ordered structures!")
        return [site.specie for site in self]

    @property
    def species_and_occu(self) -> list[Composition]:
        """List of species and occupancies at each site of the structure."""
        return [site.species for site in self]

    @property
    def ntypesp(self) -> int:
        """Number of types of atoms."""
        return len(self.types_of_species)

    @property
    def types_of_species(self) -> tuple[Element | Species | DummySpecies]:
        """List of types of specie."""
        # Cannot use set since we want a deterministic algorithm.
        types: list[Element | Species | DummySpecies] = []
        for site in self:
            for sp, amt in site.species.items():
                if amt != 0:
                    types.append(sp)
        return tuple(sorted(set(types)))  # type: ignore

    @property
    def types_of_specie(self) -> tuple[Element | Species | DummySpecies]:
        """Specie->Species rename. Maintained for backwards compatibility."""
        return self.types_of_species

    def group_by_types(self) -> Iterator[Site | PeriodicSite]:
        """Iterate over species grouped by type."""
        for sp_typ in self.types_of_species:
            for site in self:
                if site.specie == sp_typ:
                    yield site

    def indices_from_symbol(self, symbol: str) -> tuple[int, ...]:
        """
        Returns a tuple with the sequential indices of the sites
        that contain an element with the given chemical symbol.
        """
        return tuple((idx for idx, specie in enumerate(self.species) if specie.symbol == symbol))

    @property
    def symbol_set(self) -> tuple[str, ...]:
        """
        Tuple with the set of chemical symbols.
        Note that len(symbol_set) == len(types_of_specie).
        """
        return tuple(sorted(specie.symbol for specie in self.types_of_species))

    @property
    def atomic_numbers(self) -> tuple[int, ...]:
        """List of atomic numbers."""
        try:
            return tuple(site.specie.Z for site in self)
        except AttributeError:
            raise AttributeError("atomic_numbers available only for ordered Structures")

    @property
    def site_properties(self) -> dict[str, Sequence]:
        """Returns the site properties as a dict of sequences. E.g. {"magmom": (5,-5), "charge": (-4,4)}."""
        props: dict[str, Sequence] = {}
        prop_keys: set[str] = set()
        for site in self:
            prop_keys.update(site.properties)

        for key in prop_keys:
            props[key] = [site.properties.get(key) for site in self]
        return props

    @property
    def labels(self) -> list[str]:
        """Return site labels as a list."""
        return [site.label for site in self]

    def __contains__(self, site: object) -> bool:
        return site in self.sites

    def __iter__(self) -> Iterator[Site]:
        return self.sites.__iter__()

    def __getitem__(self, ind):
        return self.sites[ind]

    def __len__(self) -> int:
        return len(self.sites)

    def __hash__(self) -> int:
        # for now, just use the composition hash code.
        return hash(self.composition)

    @property
    def num_sites(self) -> int:
        """Number of sites."""
        return len(self)

    @property
    def cart_coords(self) -> np.ndarray:
        """
        Returns an np.array of the Cartesian coordinates of sites in the
        structure.
        """
        return np.array([site.coords for site in self])

    @property
    def formula(self) -> str:
        """(str) Returns the formula."""
        return self.composition.formula

    @property
    def composition(self) -> Composition:
        """(Composition) Returns the composition."""
        elem_map: dict[Species, float] = collections.defaultdict(float)
        for site in self:
            for species, occu in site.species.items():
                elem_map[species] += occu
        return Composition(elem_map)

    @property
    def charge(self) -> float:
        """
        Returns the net charge of the structure based on oxidation states. If
        Elements are found, a charge of 0 is assumed.
        """
        charge = 0
        for site in self:
            for specie, amt in site.species.items():
                charge += (getattr(specie, "oxi_state", 0) or 0) * amt

        return charge

    @property
    def is_ordered(self) -> bool:
        """
        Checks if structure is ordered, meaning no partial occupancies in any
        of the sites.
        """
        return all(site.is_ordered for site in self)

    def get_angle(self, i: int, j: int, k: int) -> float:
        """
        Returns angle specified by three sites.

        Args:
            i: Index of first site.
            j: Index of second site.
            k: Index of third site.

        Returns:
            Angle in degrees.
        """
        v1 = self[i].coords - self[j].coords
        v2 = self[k].coords - self[j].coords
        return get_angle(v1, v2, units="degrees")

    def get_dihedral(self, i: int, j: int, k: int, l: int) -> float:  # noqa: E741
        """
        Returns dihedral angle specified by four sites.

        Args:
            i: Index of first site
            j: Index of second site
            k: Index of third site
            l: Index of fourth site

        Returns:
            Dihedral angle in degrees.
        """
        v1 = self[k].coords - self[l].coords
        v2 = self[j].coords - self[k].coords
        v3 = self[i].coords - self[j].coords
        v23 = np.cross(v2, v3)
        v12 = np.cross(v1, v2)
        return math.degrees(math.atan2(np.linalg.norm(v2) * np.dot(v1, v23), np.dot(v12, v23)))

    def is_valid(self, tol: float = DISTANCE_TOLERANCE) -> bool:
        """
        True if SiteCollection does not contain atoms that are too close
        together. Note that the distance definition is based on type of
        SiteCollection. Cartesian distances are used for non-periodic
        Molecules, while PBC is taken into account for periodic structures.

        Args:
            tol (float): Distance tolerance. Default is 0.5A.

        Returns:
            (bool) True if SiteCollection does not contain atoms that are too
            close together.
        """
        if len(self) == 1:
            return True
        all_dists = self.distance_matrix[np.triu_indices(len(self), 1)]
        return bool(np.min(all_dists) > tol)

    @abstractmethod
    def to(self, filename: str = "", fmt: str = "") -> str | None:
        """
        Generates string representations (cif, json, poscar, ....) of SiteCollections (e.g.,
        molecules / structures). Should return str or None if written to a file.
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_str(cls, input_string: str, fmt: Any):
        """Reads in SiteCollection from a string."""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_file(cls, filename: str):
        """Reads in SiteCollection from a filename."""
        raise NotImplementedError

    def add_site_property(self, property_name: str, values: list):
        """
        Adds a property to a site. Note: This is the preferred method
        for adding magnetic moments, selective dynamics, and related
        site-specific properties to a structure/molecule object.

        Examples:
            structure.add_site_property("magmom", [1.0, 0.0])
            structure.add_site_property("selective_dynamics", [[True, True, True], [False, False, False]])

        Args:
            property_name (str): The name of the property to add.
            values (list): A sequence of values. Must be same length as
                number of sites.
        """
        if len(values) != len(self):
            raise ValueError(f"Values has length {len(values)} but there are {len(self)} sites! Must be same length.")
        for site, val in zip(self.sites, values):
            site.properties[property_name] = val

    def remove_site_property(self, property_name: str):
        """
        Removes a property to a site.

        Args:
            property_name (str): The name of the property to remove.
        """
        for site in self:
            del site.properties[property_name]

    def replace_species(self, species_mapping: dict[SpeciesLike, SpeciesLike | dict[SpeciesLike, float]]) -> None:
        """
        Swap species. Note that this method modifies the structure in place.

        Args:
            species_mapping (dict): dict of species to swap. Species can be elements too. E.g.,
                {Element("Li"): Element("Na")} performs a Li for Na substitution. The second species can
                be a sp_and_occu dict. For example, a site with 0.5 Si that is passed the mapping
                {Element('Si): {Element('Ge'): 0.75, Element('C'): 0.25} } will have .375 Ge and .125 C.
        """
        sp_mapping = {get_el_sp(k): v for k, v in species_mapping.items()}
        sp_to_replace = set(sp_mapping)
        sp_in_structure = set(self.composition)
        if not sp_in_structure.issuperset(sp_to_replace):
            warnings.warn(
                "Some species to be substituted are not present in structure. Pls check your input. Species to be "
                f"substituted = {sp_to_replace}; Species in structure = {sp_in_structure}"
            )

        for site in self:
            if sp_to_replace.intersection(site.species):
                comp = Composition()
                for sp, amt in site.species.items():
                    new_sp = sp_mapping.get(sp, sp)
                    try:
                        comp += Composition(new_sp) * amt
                    except Exception:
                        comp += {new_sp: amt}
                site.species = comp

    def add_oxidation_state_by_element(self, oxidation_states: dict[str, float]) -> None:
        """
        Add oxidation states.

        Args:
            oxidation_states (dict): Dict of oxidation states.
                E.g., {"Li":1, "Fe":2, "P":5, "O":-2}

        Raises:
            ValueError if oxidation states are not specified for all elements.
        """
        missing = {el.symbol for el in self.composition} - {*oxidation_states}
        if missing:
            raise ValueError(f"Oxidation states not specified for all elements, {missing=}")
        for site in self:
            new_sp = {}
            for el, occu in site.species.items():
                new_sp[Species(el.symbol, oxidation_states[el.symbol])] = occu
            site.species = Composition(new_sp)

    def add_oxidation_state_by_site(self, oxidation_states: list[float]) -> None:
        """
        Add oxidation states to a structure by site.

        Args:
            oxidation_states (list[float]): List of oxidation states.
                E.g. [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2]
        """
        if len(oxidation_states) != len(self):
            raise ValueError(
                f"Oxidation states of all sites must be specified, expected {len(self)} values, "
                f"got {len(oxidation_states)}"
            )
        for site, ox in zip(self.sites, oxidation_states):
            new_sp = {}
            for el, occu in site.species.items():
                sym = el.symbol
                new_sp[Species(sym, ox)] = occu
            site.species = Composition(new_sp)

    def remove_oxidation_states(self) -> None:
        """Removes oxidation states from a structure."""
        for site in self:
            new_sp: dict[Element, float] = collections.defaultdict(float)
            for el, occu in site.species.items():
                sym = el.symbol
                new_sp[Element(sym)] += occu
            site.species = Composition(new_sp)

    def add_oxidation_state_by_guess(self, **kwargs) -> None:
        """
        Decorates the structure with oxidation state, guessing
        using Composition.oxi_state_guesses().

        Args:
            **kwargs: parameters to pass into oxi_state_guesses()
        """
        oxid_guess = self.composition.oxi_state_guesses(**kwargs)
        oxid_guess = oxid_guess or [{e.symbol: 0 for e in self.composition}]
        self.add_oxidation_state_by_element(oxid_guess[0])

    def add_spin_by_element(self, spins: dict[str, float]) -> None:
        """
        Add spin states to structure.

        Args:
            spins (dict): Dict of spins associated with elements or species,
                e.g. {"Ni":+5} or {"Ni2+":5}
        """
        for site in self:
            new_species = {}
            for sp, occu in site.species.items():
                sym = sp.symbol
                oxi_state = getattr(sp, "oxi_state", None)
                props = {"spin": spins.get(str(sp), spins.get(sym))}
                species = Species(sym, oxidation_state=oxi_state, properties=props)
                new_species[species] = occu
            site.species = Composition(new_species)

    def add_spin_by_site(self, spins: list[float]) -> None:
        """
        Add spin states to structure by site.

        Args:
            spins (list): List of spins
                E.g., [+5, -5, 0, 0]
        """
        if len(spins) != len(self):
            raise ValueError(f"Spins for all sites must be specified, expected {len(self)} spins, got {len(spins)}")

        for site, spin in zip(self.sites, spins):
            new_species = {}
            for sp, occu in site.species.items():
                sym = sp.symbol
                oxi_state = getattr(sp, "oxi_state", None)
                new_species[Species(sym, oxidation_state=oxi_state, properties={"spin": spin})] = occu
            site.species = Composition(new_species)

    def remove_spin(self) -> None:
        """Remove spin states from structure."""
        for site in self:
            new_sp: dict[Element, float] = collections.defaultdict(float)
            for sp, occu in site.species.items():
                oxi_state = getattr(sp, "oxi_state", None)
                new_sp[Species(sp.symbol, oxidation_state=oxi_state)] += occu
            site.species = Composition(new_sp)

    def extract_cluster(self, target_sites: list[Site], **kwargs) -> list[Site]:
        """
        Extracts a cluster of atoms based on bond lengths.

        Args:
            target_sites (list[Site]): Initial sites from which to nucleate cluster.
            **kwargs: kwargs passed through to CovalentBond.is_bonded.

        Returns:
            list[Site/PeriodicSite] Cluster of atoms.
        """
        cluster = list(target_sites)
        others = [site for site in self if site not in cluster]
        size = 0
        while len(cluster) > size:
            size = len(cluster)
            new_others = []
            for site in others:
                for site2 in cluster:
                    if CovalentBond.is_bonded(site, site2, **kwargs):
                        cluster.append(site)
                        break
                else:
                    new_others.append(site)
            others = new_others
        return cluster

    def _calculate(self, calculator: str | Calculator, verbose: bool = False) -> Calculator:
        """
        Performs an ASE calculation.

        Args:
            calculator (str | Calculator): An ASE Calculator or a string from the following case-insensitive
                options: "m3gnet", "gfn2-xtb", "chgnet".
            verbose (bool): whether to print stdout. Defaults to False.
                Has no effect when calculator=='chgnet'.

        Returns:
            Structure | Molecule: Structure or Molecule with new calc attribute containing
                result of ASE calculation.
        """
        from pymatgen.io.ase import AseAtomsAdaptor

        if isinstance(self, Molecule) and isinstance(calculator, str) and calculator.lower() in ("chgnet", "m3gnet"):
            raise ValueError(f"Can't use {calculator=} for a Molecule.")
        calculator = self._prep_calculator(calculator)

        # Get Atoms object
        adaptor = AseAtomsAdaptor()
        atoms = adaptor.get_atoms(self)

        # Set calculator
        atoms.calc = calculator

        # Run calculation
        stream = sys.stdout if verbose else io.StringIO()
        with contextlib.redirect_stdout(stream):
            atoms.get_potential_energy()

        return calculator

    def _relax(
        self,
        calculator: str | Calculator,
        relax_cell: bool = True,
        optimizer: str | Optimizer = "FIRE",
        steps: int = 500,
        fmax: float = 0.1,
        stress_weight: float = 0.01,
        opt_kwargs: dict | None = None,
        return_trajectory: bool = False,
        verbose: bool = False,
    ) -> Structure | Molecule | tuple[Structure | Molecule, TrajectoryObserver | Trajectory]:
        """
        Performs a structure relaxation using an ASE calculator.

        Args:
            calculator (str | ase.Calculator): An ASE Calculator or a string from the following options: "M3GNet",
                "gfn2-xtb".
            relax_cell (bool): whether to relax the lattice cell. Defaults to True.
            optimizer (str): name of the ASE optimizer class to use
            steps (int): max number of steps for relaxation. Defaults to 500.
            fmax (float): total force tolerance for relaxation convergence.
                Here fmax is a sum of force and stress forces. Defaults to 0.1.
            stress_weight (float): the stress weight for relaxation with M3GNet.
                Defaults to 0.01.
            opt_kwargs (dict): kwargs for the ASE optimizer class.
            return_trajectory (bool): Whether to return the trajectory of relaxation.
                Defaults to False.
            verbose (bool): whether to print stdout. Defaults to False.

        Returns:
            Structure | Molecule: Relaxed structure or molecule
        """
        from ase import optimize
        from ase.constraints import ExpCellFilter
        from ase.io import read
        from ase.optimize.optimize import Optimizer

        from pymatgen.io.ase import AseAtomsAdaptor

        opt_kwargs = opt_kwargs or {}
        is_molecule = isinstance(self, Molecule)
        # UIP=universal interatomic potential
        run_uip = isinstance(calculator, str) and calculator.lower() in ("m3gnet", "chgnet")

        calc_params = dict(stress_weight=stress_weight) if not is_molecule else {}
        calculator = self._prep_calculator(calculator, **calc_params)

        # check str is valid optimizer key
        def is_ase_optimizer(key):
            return isclass(obj := getattr(optimize, key)) and issubclass(obj, Optimizer)

        valid_keys = [key for key in dir(optimize) if is_ase_optimizer(key)]
        if isinstance(optimizer, str):
            if optimizer not in valid_keys:
                raise ValueError(f"Unknown {optimizer=}, must be one of {valid_keys}")
            opt_class = getattr(optimize, optimizer)

        # Get Atoms object
        adaptor = AseAtomsAdaptor()
        atoms = adaptor.get_atoms(self)

        # Use a TrajectoryObserver if running M3GNet or CHGNet; otherwise, write a .traj file
        if return_trajectory:
            if run_uip:
                from matgl.ext.ase import TrajectoryObserver

                traj_observer = TrajectoryObserver(atoms)
            else:
                opt_kwargs.setdefault("trajectory", "opt.traj")

        # Attach calculator
        atoms.calc = calculator

        # Run relaxation
        stream = sys.stdout if verbose else io.StringIO()
        with contextlib.redirect_stdout(stream):
            if relax_cell:
                if is_molecule:
                    raise ValueError("Can't relax cell for a Molecule.")
                ecf = ExpCellFilter(atoms)
                dyn = opt_class(ecf, **opt_kwargs)
            else:
                dyn = opt_class(atoms, **opt_kwargs)
            dyn.run(fmax=fmax, steps=steps)

        # Get pymatgen Structure or Molecule
        system: Structure | Molecule = adaptor.get_molecule(atoms) if is_molecule else adaptor.get_structure(atoms)

        # Attach important ASE results
        system.calc = atoms.calc
        system.dynamics = dyn.todict()

        if return_trajectory:
            if run_uip:
                traj_observer()  # save properties of the Atoms during the relaxation
            else:
                traj_file = opt_kwargs["trajectory"]
                traj_observer = read(traj_file, index=":")
            return system, traj_observer

        return system

    def _prep_calculator(self, calculator: Literal["m3gnet", "gfn2-xtb"] | Calculator, **params) -> Calculator:
        """
        Convert string name of special ASE calculators into ASE calculator objects.

        Args:
            calculator: An ASE Calculator or a string from the following options: "m3gnet",
                "gfn2-xtb".
            **params: Parameters for the calculator.

        Returns:
            Calculator: ASE calculator object.
        """
        if inspect.isclass(calculator):
            return calculator(**params)
        if not isinstance(calculator, str):
            return calculator

        if calculator.lower() == "chgnet":
            try:
                from chgnet.model import CHGNetCalculator
            except ImportError:
                raise ImportError("chgnet not installed. Try `pip install chgnet`.")
            return CHGNetCalculator()

        if calculator.lower() == "m3gnet":
            try:
                import matgl
                from matgl.ext.ase import M3GNetCalculator
            except ImportError:
                raise ImportError("matgl not installed. Try `pip install matgl`.")
            potential = matgl.load_model("M3GNet-MP-2021.2.8-PES")
            return M3GNetCalculator(potential=potential, **params)

        if calculator.lower() == "gfn2-xtb":
            try:
                from tblite.ase import TBLite
            except ImportError:
                raise ImportError(
                    "Must install tblite[ase]. Try `pip install tblite[ase]` (Linux)"
                    "or `conda install -c conda-forge tblite-python` on (Mac/Windows)."
                )

            return TBLite(method="GFN2-xTB", **params)

        raise ValueError(f"Unknown {calculator=}.")


class IStructure(SiteCollection, MSONable):
    """
    Basic immutable Structure object with periodicity. Essentially a sequence
    of PeriodicSites having a common lattice. IStructure is made to be
    (somewhat) immutable so that they can function as keys in a dict. To make
    modifications, use the standard Structure object instead. Structure
    extends Sequence and Hashable, which means that in many cases,
    it can be used like any Python sequence. Iterating through a
    structure is equivalent to going through the sites in sequence.
    """

    def __init__(
        self,
        lattice: ArrayLike | Lattice,
        species: Sequence[CompositionLike],
        coords: Sequence[ArrayLike],
        charge: float | None = None,
        validate_proximity: bool = False,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        site_properties: dict | None = None,
        labels: list | None = None,
    ) -> None:
        """
        Create a periodic structure.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/Cartesian coordinates of
                each species.
            charge (int): overall charge of the structure. Defaults to behavior
                in SiteCollection where total charge is the sum of the oxidation
                states.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            to_unit_cell (bool): Whether to map all sites into the unit cell,
                i.e. fractional coords between 0 and 1. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g. {"magmom":[5, 5, 5, 5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            labels (list[str]): Labels associated with the sites as a
                list of strings, e.g. ['Li1', 'Li2']. Must have the same
                length as the species and fractional coords. Defaults to
                None for no labels.
        """
        if len(species) != len(coords):
            raise StructureError("atomic species and fractional coordinates must have same length")

        self._lattice = lattice if isinstance(lattice, Lattice) else Lattice(lattice)

        sites = []
        for idx, specie in enumerate(species):
            prop = None
            if site_properties:
                prop = {key: val[idx] for key, val in site_properties.items()}

            label = labels[idx] if labels else None

            site = PeriodicSite(
                specie,
                coords[idx],
                self._lattice,
                to_unit_cell,
                coords_are_cartesian=coords_are_cartesian,
                properties=prop,  # type: ignore
                label=label,
            )
            sites.append(site)
        self._sites: tuple[PeriodicSite, ...] = tuple(sites)
        if validate_proximity and not self.is_valid():
            raise StructureError("Structure contains sites that are less than 0.01 Angstrom apart!")
        self._charge = charge

    @classmethod
    def from_sites(
        cls,
        sites: list[PeriodicSite],
        charge: float | None = None,
        validate_proximity: bool = False,
        to_unit_cell: bool = False,
    ) -> IStructure:
        """
        Convenience constructor to make a Structure from a list of sites.

        Args:
            sites: Sequence of PeriodicSites. Sites must have the same
                lattice.
            charge: Charge of structure.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            to_unit_cell (bool): Whether to translate sites into the unit
                cell.

        Returns:
            (Structure) Note that missing properties are set as None.
        """
        if len(sites) < 1:
            raise ValueError(f"You need at least one site to construct a {cls}")
        prop_keys: list[str] = []
        props = {}
        labels = [site.label for site in sites]
        lattice = sites[0].lattice
        for idx, site in enumerate(sites):
            if site.lattice != lattice:
                raise ValueError("Sites must belong to the same lattice")
            for key, val in site.properties.items():
                if key not in prop_keys:
                    prop_keys.append(key)
                    props[key] = [None] * len(sites)
                props[key][idx] = val
        for key, val in props.items():
            if any(vv is None for vv in val):
                warnings.warn(f"Not all sites have property {key}. Missing values are set to None.")
        return cls(
            lattice,
            [site.species for site in sites],
            [site.frac_coords for site in sites],
            charge=charge,
            site_properties=props,
            validate_proximity=validate_proximity,
            to_unit_cell=to_unit_cell,
            labels=labels,
        )

    @classmethod
    def from_spacegroup(
        cls,
        sg: str | int,
        lattice: list | np.ndarray | Lattice,
        species: Sequence[str | Element | Species | DummySpecies | Composition],
        coords: Sequence[Sequence[float]],
        site_properties: dict[str, Sequence] | None = None,
        coords_are_cartesian: bool = False,
        tol: float = 1e-5,
    ) -> IStructure | Structure:
        """
        Generate a structure using a spacegroup. Note that only symmetrically
        distinct species and coords should be provided. All equivalent sites
        are generated from the spacegroup operations.

        Args:
            sg (str/int): The spacegroup. If a string, it will be interpreted
                as one of the notations supported by
                pymatgen.symmetry.groups.Spacegroup. E.g., "R-3c" or "Fm-3m".
                If an int, it will be interpreted as an international number.
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
                Note that no attempt is made to check that the lattice is
                compatible with the spacegroup specified. This may be
                introduced in a future version.
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            tol (float): A fractional tolerance to deal with numerical
               precision issues in determining if orbits are the same.
        """
        from pymatgen.symmetry.groups import SpaceGroup

        try:
            num = int(sg)
            spg = SpaceGroup.from_int_number(num)
        except ValueError:
            spg = SpaceGroup(sg)  # type: ignore

        latt = lattice if isinstance(lattice, Lattice) else Lattice(lattice)

        if not spg.is_compatible(latt):
            raise ValueError(
                f"Supplied lattice with parameters {latt.parameters} is incompatible with supplied spacegroup "
                f"{spg.symbol}!"
            )

        if len(species) != len(coords):
            raise ValueError(f"Supplied species and coords lengths ({len(species)} vs {len(coords)}) are different!")

        frac_coords = (
            np.array(coords, dtype=np.float_) if not coords_are_cartesian else latt.get_fractional_coords(coords)
        )

        props = {} if site_properties is None else site_properties

        all_sp: list[str | Element | Species | DummySpecies | Composition] = []
        all_coords: list[list[float]] = []
        all_site_properties: dict[str, list] = collections.defaultdict(list)
        for idx, (sp, c) in enumerate(zip(species, frac_coords)):
            cc = spg.get_orbit(c, tol=tol)
            all_sp.extend([sp] * len(cc))
            all_coords.extend(cc)  # type: ignore
            for k, v in props.items():
                all_site_properties[k].extend([v[idx]] * len(cc))

        return cls(latt, all_sp, all_coords, site_properties=all_site_properties)

    @classmethod
    def from_magnetic_spacegroup(
        cls,
        msg: str | MagneticSpaceGroup,
        lattice: list | np.ndarray | Lattice,
        species: Sequence[str | Element | Species | DummySpecies | Composition],
        coords: Sequence[Sequence[float]],
        site_properties: dict[str, Sequence],
        coords_are_cartesian: bool = False,
        tol: float = 1e-5,
    ) -> IStructure | Structure:
        """
        Generate a structure using a magnetic spacegroup. Note that only
        symmetrically distinct species, coords and magmoms should be provided.]
        All equivalent sites are generated from the spacegroup operations.

        Args:
            msg (str/list/:class:`pymatgen.symmetry.maggroups.MagneticSpaceGroup`):
                The magnetic spacegroup.
                If a string, it will be interpreted as one of the notations
                supported by MagneticSymmetryGroup, e.g., "R-3'c" or "Fm'-3'm".
                If a list of two ints, it will be interpreted as the number of
                the spacegroup in its Belov, Neronova and Smirnova (BNS) setting.
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
                Note that no attempt is made to check that the lattice is
                compatible with the spacegroup specified. This may be
                introduced in a future version.
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:
                i.  A sequence of element / species specified either as string
                symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                e.g., (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Unlike Structure.from_spacegroup(),
                this argument is mandatory, since magnetic moment information
                has to be included. Note that the *direction* of the supplied
                magnetic moment relative to the crystal is important, even if
                the resulting structure is used for collinear calculations.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            tol (float): A fractional tolerance to deal with numerical
                precision issues in determining if orbits are the same.

        Returns:
            Structure | IStructure
        """
        if "magmom" not in site_properties:
            raise ValueError("Magnetic moments have to be defined.")

        magmoms = [Magmom(m) for m in site_properties["magmom"]]

        if not isinstance(msg, MagneticSpaceGroup):
            msg = MagneticSpaceGroup(msg)

        latt = lattice if isinstance(lattice, Lattice) else Lattice(lattice)

        if not msg.is_compatible(latt):
            raise ValueError(
                f"Supplied lattice with parameters {latt.parameters} is incompatible with supplied spacegroup "
                f"{msg.sg_symbol}!"
            )

        if len(species) != len(coords):
            raise ValueError(f"Supplied species and coords lengths ({len(species)} vs {len(coords)}) are different!")

        if len(species) != len(magmoms):
            raise ValueError(f"Supplied species and magmom lengths ({len(species)} vs {len(magmoms)}) are different!")

        frac_coords = coords if not coords_are_cartesian else latt.get_fractional_coords(coords)

        all_sp: list[str | Element | Species | DummySpecies | Composition] = []
        all_coords: list[list[float]] = []
        all_magmoms: list[float] = []
        all_site_properties: dict[str, list] = collections.defaultdict(list)
        for i, (sp, c, m) in enumerate(zip(species, frac_coords, magmoms)):  # type: ignore
            cc, mm = msg.get_orbit(c, m, tol=tol)
            all_sp.extend([sp] * len(cc))
            all_coords.extend(cc)
            all_magmoms.extend(mm)
            for k, v in site_properties.items():
                if k != "magmom":
                    all_site_properties[k].extend([v[i]] * len(cc))

        all_site_properties["magmom"] = all_magmoms

        return cls(latt, all_sp, all_coords, site_properties=all_site_properties)

    def unset_charge(self):
        """Reset the charge to None, i.e., computed dynamically based on oxidation states."""
        self._charge = None

    @property
    def charge(self) -> float:
        """Overall charge of the structure."""
        formal_charge = super().charge
        if self._charge is None:
            return super().charge
        if formal_charge != self._charge:
            warnings.warn(
                f"Structure charge ({self._charge}) is set to be not equal to the sum of oxidation states"
                f" ({formal_charge}). Use `unset_charge` if this is not desired."
            )
        return self._charge

    @property
    def distance_matrix(self) -> np.ndarray:
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this should return the nearest image distance.
        """
        return self.lattice.get_all_distances(self.frac_coords, self.frac_coords)

    @property
    def sites(self) -> tuple[PeriodicSite, ...]:
        """Returns an iterator for the sites in the Structure."""
        return self._sites

    @property
    def lattice(self) -> Lattice:
        """Lattice of the structure."""
        return self._lattice

    @property
    def density(self) -> float:
        """Returns the density in units of g/cm^3."""
        m = Mass(self.composition.weight, "amu")
        return m.to("g") / (self.volume * Length(1, "ang").to("cm") ** 3)

    @property
    def pbc(self) -> tuple[bool, bool, bool]:
        """Returns the periodicity of the structure."""
        return self._lattice.pbc

    @property
    def is_3d_periodic(self) -> bool:
        """True if the Lattice is periodic in all directions."""
        return self._lattice.is_3d_periodic

    def get_space_group_info(self, symprec=1e-2, angle_tolerance=5.0) -> tuple[str, int]:
        """
        Convenience method to quickly get the spacegroup of a structure.

        Args:
            symprec (float): Same definition as in SpacegroupAnalyzer.
                Defaults to 1e-2.
            angle_tolerance (float): Same definition as in SpacegroupAnalyzer.
                Defaults to 5 degrees.

        Returns:
            spacegroup_symbol, international_number
        """
        # Import within method needed to avoid cyclic dependency.
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        spga = SpacegroupAnalyzer(self, symprec=symprec, angle_tolerance=angle_tolerance)
        return spga.get_space_group_symbol(), spga.get_space_group_number()

    def matches(self, other: IStructure | Structure, anonymous: bool = False, **kwargs) -> bool:
        """
        Check whether this structure is similar to another structure.
        Basically a convenience method to call structure matching.

        Args:
            other (IStructure/Structure): Another structure.
            anonymous (bool): Whether to use anonymous structure matching which allows distinct
                species in one structure to map to another.
            **kwargs: Same **kwargs as in
                :class:`pymatgen.analysis.structure_matcher.StructureMatcher`.

        Returns:
            bool: True if the structures are similar under some affine transformation.
        """
        from pymatgen.analysis.structure_matcher import StructureMatcher

        matcher = StructureMatcher(**kwargs)
        if anonymous:
            return matcher.fit_anonymous(self, other)
        return matcher.fit(self, other)

    def __eq__(self, other: object) -> bool:
        # check for valid operand following class Student example from official functools docs
        # https://docs.python.org/3/library/functools.html#functools.total_ordering
        if not isinstance(other, IStructure):
            return NotImplemented

        if other is self:
            return True
        if len(self) != len(other):
            return False
        if self.lattice != other.lattice:
            return False
        return all(site in other for site in self)

    def __hash__(self) -> int:
        # For now, just use the composition hash code.
        return hash(self.composition)

    def __mul__(self, scaling_matrix: int | Sequence[int] | Sequence[Sequence[int]]) -> Structure:
        """
        Makes a supercell. Allowing to have sites outside the unit cell.

        Args:
            scaling_matrix: A scaling matrix for transforming the lattice
                vectors. Has to be all integers. Several options are possible:

                a. A full 3x3 scaling matrix defining the linear combination
                   of the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
                   1]] generates a new structure with lattice vectors a' =
                   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
                   vectors of the original structure.
                b. A sequence of three scaling factors. E.g., [2, 1, 1]
                   specifies that the supercell should have dimensions 2a x b x
                   c.
                c. A number, which simply scales all lattice vectors by the
                   same factor.

        Returns:
            Supercell structure. Note that a Structure is always returned,
            even if the input structure is a subclass of Structure. This is
            to avoid different arguments signatures from causing problems. If
            you prefer a subclass to return its own type, you need to override
            this method in the subclass.
        """
        scale_matrix = np.array(scaling_matrix, int)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), int)
        new_lattice = Lattice(np.dot(scale_matrix, self._lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        for site in self:
            for v in c_lat:
                s = PeriodicSite(
                    site.species,
                    site.coords + v,
                    new_lattice,
                    properties=site.properties,
                    coords_are_cartesian=True,
                    to_unit_cell=False,
                    skip_checks=True,
                )
                new_sites.append(s)

        new_charge = self._charge * np.linalg.det(scale_matrix) if self._charge else None
        return Structure.from_sites(new_sites, charge=new_charge, to_unit_cell=True)

    def __rmul__(self, scaling_matrix):
        """Similar to __mul__ to preserve commutativeness."""
        return self.__mul__(scaling_matrix)

    @property
    def frac_coords(self):
        """Fractional coordinates as a Nx3 numpy array."""
        return np.array([site.frac_coords for site in self._sites])

    @property
    def volume(self) -> float:
        """Returns the volume of the structure in Angstrom^3."""
        return self._lattice.volume

    def get_distance(self, i: int, j: int, jimage=None) -> float:
        """
        Get distance between site i and j assuming periodic boundary
        conditions. If the index jimage of two sites atom j is not specified it
        selects the jimage nearest to the i atom and returns the distance and
        jimage indices in terms of lattice vector translations if the index
        jimage of atom j is specified it returns the distance between the i
        atom and the specified jimage atom.

        Args:
            i (int): Index of first site
            j (int): Index of second site
            jimage: Number of lattice translations in each lattice direction.
                Default is None for nearest image.

        Returns:
            distance
        """
        return self[i].distance(self[j], jimage)

    def get_sites_in_sphere(
        self,
        pt: ArrayLike,
        r: float,
        include_index: bool = False,
        include_image: bool = False,
    ) -> list[PeriodicNeighbor]:
        """
        Find all sites within a sphere from the point, including a site (if any)
        sitting on the point itself. This includes sites in other periodic
        images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelepiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            pt (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Radius of sphere.
            include_index (bool): Whether the non-supercell site index
                is included in the returned data
            include_image (bool): Whether to include the supercell image
                is included in the returned data

        Returns:
            PeriodicNeighbor
        """
        site_fcoords = np.mod(self.frac_coords, 1)
        neighbors: list[PeriodicNeighbor] = []
        for frac_coord, dist, i, img in self._lattice.get_points_in_sphere(site_fcoords, pt, r):
            nn_site = PeriodicNeighbor(
                self[i].species,
                frac_coord,
                self._lattice,
                properties=self[i].properties,
                nn_distance=dist,
                image=img,  # type: ignore
                index=i,
            )
            neighbors.append(nn_site)
        return neighbors

    def get_neighbors(
        self,
        site: PeriodicSite,
        r: float,
        include_index: bool = False,
        include_image: bool = False,
    ) -> list[PeriodicNeighbor]:
        """
        Get all neighbors to a site within a sphere of radius r. Excludes the
        site itself.

        Args:
            site (Site): Which is the center of the sphere.
            r (float): Radius of sphere.
            include_index (bool): Deprecated. Now, the non-supercell site index
                is always included in the returned data.
            include_image (bool): Deprecated. Now the supercell image
                is always included in the returned data.

        Returns:
            PeriodicNeighbor
        """
        return self.get_all_neighbors(r, include_index=include_index, include_image=include_image, sites=[site])[0]

    @deprecated(get_neighbors, "This is retained purely for checking purposes.")
    def get_neighbors_old(self, site, r, include_index=False, include_image=False):
        """
        Get all neighbors to a site within a sphere of radius r. Excludes the
        site itself.

        Args:
            site (Site): Which is the center of the sphere.
            r (float): Radius of sphere.
            include_index (bool): Whether the non-supercell site index
                is included in the returned data
            include_image (bool): Whether to include the supercell image
                is included in the returned data

        Returns:
            PeriodicNeighbor
        """
        nn = self.get_sites_in_sphere(site.coords, r, include_index=include_index, include_image=include_image)
        return [d for d in nn if site != d[0]]

    def _get_neighbor_list_py(
        self,
        r: float,
        sites: list[PeriodicSite] | None = None,
        numerical_tol: float = 1e-8,
        exclude_self: bool = True,
    ) -> tuple[np.ndarray, ...]:
        """
        A python version of getting neighbor_list. The returned values are a tuple of
        numpy arrays (center_indices, points_indices, offset_vectors, distances).
        Atom `center_indices[i]` has neighbor atom `points_indices[i]` that is
        translated by `offset_vectors[i]` lattice vectors, and the distance is
        `distances[i]`.

        Args:
            r (float): Radius of sphere
            sites (list of Sites or None): sites for getting all neighbors,
                default is None, which means neighbors will be obtained for all
                sites. This is useful in the situation where you are interested
                only in one subspecies type, and makes it a lot faster.
            numerical_tol (float): This is a numerical tolerance for distances.
                Sites which are < numerical_tol are determined to be coincident
                with the site. Sites which are r + numerical_tol away is deemed
                to be within r from the site. The default of 1e-8 should be
                ok in most instances.
            exclude_self (bool): whether to exclude atom neighboring with itself within
                numerical tolerance distance, default to True
        Returns: (center_indices, points_indices, offset_vectors, distances)
        """
        neighbors = self.get_all_neighbors_py(
            r=r, include_index=True, include_image=True, sites=sites, numerical_tol=1e-8
        )
        center_indices = []
        points_indices = []
        offsets = []
        distances = []
        for idx, nns in enumerate(neighbors):
            if len(nns) > 0:
                for nn in nns:
                    if exclude_self and idx == nn.index and nn.nn_distance <= numerical_tol:
                        continue
                    center_indices.append(idx)
                    points_indices.append(nn.index)
                    offsets.append(nn.image)
                    distances.append(nn.nn_distance)
        return tuple(map(np.array, (center_indices, points_indices, offsets, distances)))

    def get_neighbor_list(
        self,
        r: float,
        sites: Sequence[PeriodicSite] | None = None,
        numerical_tol: float = 1e-8,
        exclude_self: bool = True,
    ) -> tuple[np.ndarray, ...]:
        """
        Get neighbor lists using numpy array representations without constructing
        Neighbor objects. If the cython extension is installed, this method will
        be orders of magnitude faster than `get_all_neighbors_old` and 2-3x faster
        than `get_all_neighbors`.
        The returned values are a tuple of numpy arrays
        (center_indices, points_indices, offset_vectors, distances).
        Atom `center_indices[i]` has neighbor atom `points_indices[i]` that is
        translated by `offset_vectors[i]` lattice vectors, and the distance is
        `distances[i]`.

        Args:
            r (float): Radius of sphere
            sites (list of Sites or None): sites for getting all neighbors,
                default is None, which means neighbors will be obtained for all
                sites. This is useful in the situation where you are interested
                only in one subspecies type, and makes it a lot faster.
            numerical_tol (float): This is a numerical tolerance for distances.
                Sites which are < numerical_tol are determined to be coincident
                with the site. Sites which are r + numerical_tol away is deemed
                to be within r from the site. The default of 1e-8 should be
                ok in most instances.
            exclude_self (bool): whether to exclude atom neighboring with itself within
                numerical tolerance distance, default to True
        Returns: (center_indices, points_indices, offset_vectors, distances)
        """
        try:
            from pymatgen.optimization.neighbors import find_points_in_spheres
        except ImportError:
            return self._get_neighbor_list_py(r, sites, exclude_self=exclude_self)  # type: ignore
        else:
            if sites is None:
                sites = self.sites
            site_coords = np.ascontiguousarray([site.coords for site in sites], dtype=float)
            cart_coords = np.ascontiguousarray(self.cart_coords, dtype=float)
            lattice_matrix = np.ascontiguousarray(self.lattice.matrix, dtype=float)
            pbc = np.ascontiguousarray(self.pbc, dtype=int)
            center_indices, points_indices, images, distances = find_points_in_spheres(
                cart_coords,
                site_coords,
                r=float(r),
                pbc=pbc,
                lattice=lattice_matrix,
                tol=numerical_tol,
            )
            cond = np.array([True] * len(center_indices))
            if exclude_self:
                self_pair = (center_indices == points_indices) & (distances <= numerical_tol)
                cond = ~self_pair
            return (center_indices[cond], points_indices[cond], images[cond], distances[cond])

    def get_symmetric_neighbor_list(
        self,
        r: float,
        sg: str,
        unique: bool = False,
        numerical_tol: float = 1e-8,
        exclude_self: bool = True,
    ) -> tuple[np.ndarray, ...]:
        """
        Similar to 'get_neighbor_list' with sites=None, but the neighbors are
        grouped by symmetry. The returned values are a tuple of numpy arrays
        (center_indices, points_indices, offset_vectors, distances,
         symmetry_indices). Atom `center_indices[i]` has neighbor atom
        `points_indices[i]` that is translated by `offset_vectors[i]` lattice
        vectors, and the distance is `distances[i]`. Symmetry_idx groups the bonds
        that are related by a symmetry of the provided space group and symmetry_op
        is the operation that relates the first bond of the same symmetry_idx to
        the respective atom. The first bond maps onto itself via the Identity. The
        output is sorted w.r.t. to symmetry_indices. If unique is True only one of the
        two bonds connecting two points is given. Out of the two, the bond that does not
        reverse the sites is chosen.

        Args:
            r (float): Radius of sphere
            sg (str/int): The spacegroup the symmetry operations of which will be
                used to classify the neighbors. If a string, it will be interpreted
                as one of the notations supported by
                pymatgen.symmetry.groups.Spacegroup. E.g., "R-3c" or "Fm-3m".
                If an int, it will be interpreted as an international number.
                If None, 'get_space_group_info' will be used to determine the
                space group, default to None.
            unique (bool): Whether a bond is given for both, or only a single
                direction is given. The default is False.
            numerical_tol (float): This is a numerical tolerance for distances.
                Sites which are < numerical_tol are determined to be coincident
                with the site. Sites which are r + numerical_tol away is deemed
                to be within r from the site. The default of 1e-8 should be
                ok in most instances.
            exclude_self (bool): whether to exclude atom neighboring with itself within
                numerical tolerance distance, default to True
        Returns: (center_indices, points_indices, offset_vectors, distances,
                  symmetry_indices, symmetry_ops)
        """
        from pymatgen.symmetry.groups import SpaceGroup

        if sg is None:
            ops = SpaceGroup(self.get_space_group_info()[0]).symmetry_ops
        else:
            try:  # first assume sg is int
                sgp = SpaceGroup.from_int_number(int(sg))
            except ValueError:
                sgp = SpaceGroup(sg)
            ops = sgp.symmetry_ops

        latt = self.lattice

        if not sgp.is_compatible(latt):
            raise ValueError(
                f"Supplied lattice with parameters {latt.parameters} is incompatible with "
                f"supplied spacegroup {sgp.symbol}!"
            )

        # get a list of neighbors up to distance r
        bonds = self.get_neighbor_list(r)

        if unique:
            redundant = []
            # compare all neighbors pairwise to find the pairs that connect the same
            # two sites, but with an inverted vector (R=-R) that connects the two and add
            # one of each pair to the redundant list.
            for it, (i, j, R, d) in enumerate(zip(*bonds)):
                if it in redundant:
                    pass
                else:
                    for it2, (i2, j2, R2, d2) in enumerate(zip(*bonds)):
                        bool1 = i == j2
                        bool2 = j == i2
                        bool3 = (-R2 == R).all()
                        bool4 = np.isclose(d, d2, atol=numerical_tol)
                        if bool1 and bool2 and bool3 and bool4:
                            redundant.append(it2)

            # delete the redundant neighbors
            m = ~np.in1d(np.arange(len(bonds[0])), redundant)
            idcs_dist = np.argsort(bonds[3][m])
            bonds = (bonds[0][m][idcs_dist], bonds[1][m][idcs_dist], bonds[2][m][idcs_dist], bonds[3][m][idcs_dist])

        # expand the output tuple by symmetry_indices and symmetry_ops.
        nbonds = len(bonds[0])
        symmetry_indices = np.empty(nbonds)
        symmetry_indices[:] = np.NaN
        symmetry_ops = np.empty(len(symmetry_indices), dtype=object)
        symmetry_identity = SymmOp.from_rotation_and_translation(np.eye(3), np.zeros(3))
        symmetry_index = 0

        # Again, compare all neighbors pairwise. For each pair of neighbors, all the symmetry operations of the provided
        # space group are iterated over. If an operation is found that connects the two bonds, it is assigned the same
        # symmetry index it is compared to, and the symmetry operation that connects the two is saved. To compare two
        # neighbors 'SymmOp.are_symmetrically_related_vectors' is used. It is also checked whether applying the
        # connecting symmetry operation generates the neighbor-pair itself, or the equivalent version with the
        # sites exchanged and R reversed. The output is always reordered such that the former case is true.
        for it in range(nbonds):
            if np.isnan(symmetry_indices[it]):
                symmetry_indices[it] = symmetry_index
                symmetry_ops[it] = symmetry_identity
                for it2 in np.arange(nbonds)[np.isnan(symmetry_indices)]:
                    equal_distance = np.isclose(bonds[3][it], bonds[3][it2], atol=numerical_tol)
                    if equal_distance:
                        from_a = self[bonds[0][it]].frac_coords
                        to_a = self[bonds[1][it]].frac_coords
                        r_a = bonds[2][it]
                        from_b = self[bonds[0][it2]].frac_coords
                        to_b = self[bonds[1][it2]].frac_coords
                        r_b = bonds[2][it2]
                        for op in ops:
                            are_related, is_reversed = op.are_symmetrically_related_vectors(
                                from_a, to_a, r_a, from_b, to_b, r_b
                            )
                            if are_related and not is_reversed:
                                symmetry_indices[it2] = symmetry_index
                                symmetry_ops[it2] = op
                            elif are_related and is_reversed:
                                symmetry_indices[it2] = symmetry_index
                                symmetry_ops[it2] = op
                                bonds[0][it2], bonds[1][it2] = bonds[1][it2], bonds[0][it2]
                                bonds[2][it2] = -bonds[2][it2]

                symmetry_index += 1

        # the bonds are ordered by their symmetry index
        idcs_symid = np.argsort(symmetry_indices)
        bonds = (bonds[0][idcs_symid], bonds[1][idcs_symid], bonds[2][idcs_symid], bonds[3][idcs_symid])
        symmetry_indices = symmetry_indices[idcs_symid]
        symmetry_ops = symmetry_ops[idcs_symid]

        # the groups of neighbors with the same symmetry index are ordered such that neighbors
        # that are the first occurrence of a new symmetry index in the ordered output are the ones
        # that are assigned the Identity as a symmetry operation.
        idcs_symop = np.arange(nbonds)
        identity_idcs = np.where(symmetry_ops == symmetry_identity)[0]
        for symmetry_idx in np.unique(symmetry_indices):
            first_idx = np.argmax(symmetry_indices == symmetry_idx)
            for second_idx in identity_idcs:
                if symmetry_indices[second_idx] == symmetry_idx:
                    idcs_symop[first_idx], idcs_symop[second_idx] = idcs_symop[second_idx], idcs_symop[first_idx]

        return (
            bonds[0][idcs_symop],
            bonds[1][idcs_symop],
            bonds[2][idcs_symop],
            bonds[3][idcs_symop],
            symmetry_indices[idcs_symop],
            symmetry_ops[idcs_symop],
        )

    def get_all_neighbors(
        self,
        r: float,
        include_index: bool = False,
        include_image: bool = False,
        sites: Sequence[PeriodicSite] | None = None,
        numerical_tol: float = 1e-8,
    ) -> list[list[PeriodicNeighbor]]:
        """
        Get neighbors for each atom in the unit cell, out to a distance r
        Returns a list of list of neighbors for each site in structure.
        Use this method if you are planning on looping over all sites in the
        crystal. If you only want neighbors for a particular site, use the
        method get_neighbors as it may not have to build such a large supercell
        However if you are looping over all sites in the crystal, this method
        is more efficient since it only performs one pass over a large enough
        supercell to contain all possible atoms out to a distance r.
        The return type is a [(site, dist) ...] since most of the time,
        subsequent processing requires the distance.

        A note about periodic images: Before computing the neighbors, this
        operation translates all atoms to within the unit cell (having
        fractional coordinates within [0,1)). This means that the "image" of a
        site does not correspond to how much it has been translates from its
        current position, but which image of the unit cell it resides.

        Args:
            r (float): Radius of sphere.
            include_index (bool): Deprecated. Now, the non-supercell site index
                is always included in the returned data.
            include_image (bool): Deprecated. Now the supercell image
                is always included in the returned data.
            sites (list of Sites or None): sites for getting all neighbors,
                default is None, which means neighbors will be obtained for all
                sites. This is useful in the situation where you are interested
                only in one subspecies type, and makes it a lot faster.
            numerical_tol (float): This is a numerical tolerance for distances.
                Sites which are < numerical_tol are determined to be coincident
                with the site. Sites which are r + numerical_tol away is deemed
                to be within r from the site. The default of 1e-8 should be
                ok in most instances.

        Returns:
            [[:class:`pymatgen.core.structure.PeriodicNeighbor`], ..]
        """
        if sites is None:
            sites = self.sites
        center_indices, points_indices, images, distances = self.get_neighbor_list(
            r=r, sites=sites, numerical_tol=numerical_tol
        )
        if len(points_indices) < 1:
            return [[]] * len(sites)
        f_coords = self.frac_coords[points_indices] + images
        neighbor_dict: dict[int, list] = collections.defaultdict(list)
        lattice = self.lattice
        atol = Site.position_atol
        all_sites = self.sites
        for cindex, pindex, image, f_coord, d in zip(center_indices, points_indices, images, f_coords, distances):
            psite = all_sites[pindex]
            csite = sites[cindex]
            if (
                d > numerical_tol
                or
                # This simply compares the psite and csite. The reason why manual comparison is done is
                # for speed. This does not check the lattice since they are always equal. Also, the or construct
                # returns True immediately once one of the conditions are satisfied.
                psite.species != csite.species
                or (not np.allclose(psite.coords, csite.coords, atol=atol))
                or (psite.properties != csite.properties)
            ):
                neighbor_dict[cindex].append(
                    PeriodicNeighbor(
                        species=psite.species,
                        coords=f_coord,
                        lattice=lattice,
                        properties=psite.properties,
                        nn_distance=d,
                        index=pindex,
                        image=tuple(image),
                    )
                )

        neighbors: list[list[PeriodicNeighbor]] = []

        for i in range(len(sites)):
            neighbors.append(neighbor_dict[i])
        return neighbors

    def get_all_neighbors_py(
        self,
        r: float,
        include_index: bool = False,
        include_image: bool = False,
        sites: Sequence[PeriodicSite] | None = None,
        numerical_tol: float = 1e-8,
    ) -> list[list[PeriodicNeighbor]]:
        """
        Get neighbors for each atom in the unit cell, out to a distance r
        Returns a list of list of neighbors for each site in structure.
        Use this method if you are planning on looping over all sites in the
        crystal. If you only want neighbors for a particular site, use the
        method get_neighbors as it may not have to build such a large supercell
        However if you are looping over all sites in the crystal, this method
        is more efficient since it only performs one pass over a large enough
        supercell to contain all possible atoms out to a distance r.
        The return type is a [(site, dist) ...] since most of the time,
        subsequent processing requires the distance.

        A note about periodic images: Before computing the neighbors, this
        operation translates all atoms to within the unit cell (having
        fractional coordinates within [0,1)). This means that the "image" of a
        site does not correspond to how much it has been translates from its
        current position, but which image of the unit cell it resides.

        Args:
            r (float): Radius of sphere.
            include_index (bool): Deprecated. Now, the non-supercell site index
                is always included in the returned data.
            include_image (bool): Deprecated. Now the supercell image
                is always included in the returned data.
            sites (list of Sites or None): sites for getting all neighbors,
                default is None, which means neighbors will be obtained for all
                sites. This is useful in the situation where you are interested
                only in one subspecies type, and makes it a lot faster.
            numerical_tol (float): This is a numerical tolerance for distances.
                Sites which are < numerical_tol are determined to be coincident
                with the site. Sites which are r + numerical_tol away is deemed
                to be within r from the site. The default of 1e-8 should be
                ok in most instances.

        Returns:
            list[list[PeriodicNeighbor]]
        """
        if sites is None:
            sites = self.sites
        site_coords = np.array([site.coords for site in sites])
        point_neighbors = get_points_in_spheres(
            self.cart_coords,
            site_coords,
            r=r,
            pbc=self.pbc,
            numerical_tol=numerical_tol,
            lattice=self.lattice,
        )
        neighbors: list[list[PeriodicNeighbor]] = []
        for point_neighbor, site in zip(point_neighbors, sites):
            nns: list[PeriodicNeighbor] = []
            if len(point_neighbor) < 1:
                neighbors.append([])
                continue
            for n in point_neighbor:
                coord, d, index, image = n
                if (d > numerical_tol) or (self[index] != site):
                    neighbor = PeriodicNeighbor(
                        species=self[index].species,
                        coords=coord,
                        lattice=self.lattice,
                        properties=self[index].properties,
                        nn_distance=d,
                        index=index,
                        image=tuple(image),
                    )
                    nns.append(neighbor)
            neighbors.append(nns)
        return neighbors

    @deprecated(get_all_neighbors, "This is retained purely for checking purposes.")
    def get_all_neighbors_old(self, r, include_index=False, include_image=False, include_site=True):
        """
        Get neighbors for each atom in the unit cell, out to a distance r
        Returns a list of list of neighbors for each site in structure.
        Use this method if you are planning on looping over all sites in the
        crystal. If you only want neighbors for a particular site, use the
        method get_neighbors as it may not have to build such a large supercell
        However if you are looping over all sites in the crystal, this method
        is more efficient since it only performs one pass over a large enough
        supercell to contain all possible atoms out to a distance r.
        The return type is a [(site, dist) ...] since most of the time,
        subsequent processing requires the distance.

        A note about periodic images: Before computing the neighbors, this
        operation translates all atoms to within the unit cell (having
        fractional coordinates within [0,1)). This means that the "image" of a
        site does not correspond to how much it has been translates from its
        current position, but which image of the unit cell it resides.

        Args:
            r (float): Radius of sphere.
            include_index (bool): Whether to include the non-supercell site
                in the returned data
            include_image (bool): Whether to include the supercell image
                in the returned data
            include_site (bool): Whether to include the site in the returned
                data. Defaults to True.

        Returns:
            PeriodicNeighbor
        """
        # Use same algorithm as get_sites_in_sphere to determine supercell but
        # loop over all atoms in crystal
        recp_len = np.array(self.lattice.reciprocal_lattice.abc)
        maxr = np.ceil((r + 0.15) * recp_len / (2 * math.pi))
        nmin = np.floor(np.min(self.frac_coords, axis=0)) - maxr
        nmax = np.ceil(np.max(self.frac_coords, axis=0)) + maxr

        all_ranges = [np.arange(x, y) for x, y in zip(nmin, nmax)]
        latt = self._lattice
        matrix = latt.matrix
        neighbors = [[] for _ in range(len(self._sites))]
        all_fcoords = np.mod(self.frac_coords, 1)
        coords_in_cell = np.dot(all_fcoords, matrix)
        site_coords = self.cart_coords

        indices = np.arange(len(self))

        for image in itertools.product(*all_ranges):
            coords = np.dot(image, matrix) + coords_in_cell
            all_dists = all_distances(coords, site_coords)
            all_within_r = np.bitwise_and(all_dists <= r, all_dists > 1e-8)

            for j, d, within_r in zip(indices, all_dists, all_within_r):
                if include_site:
                    nnsite = PeriodicSite(
                        self[j].species,
                        coords[j],
                        latt,
                        properties=self[j].properties,
                        coords_are_cartesian=True,
                        skip_checks=True,
                    )

                for i in indices[within_r]:
                    item = []
                    if include_site:
                        item.append(nnsite)
                    item.append(d[i])
                    if include_index:
                        item.append(j)
                    # Add the image, if requested
                    if include_image:
                        item.append(image)
                    neighbors[i].append(item)
        return neighbors

    def get_neighbors_in_shell(
        self, origin: ArrayLike, r: float, dr: float, include_index: bool = False, include_image: bool = False
    ) -> list[PeriodicNeighbor]:
        """
        Returns all sites in a shell centered on origin (coords) between radii
        r-dr and r+dr.

        Args:
            origin (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Inner radius of shell.
            dr (float): Width of shell.
            include_index (bool): Deprecated. Now, the non-supercell site index
                is always included in the returned data.
            include_image (bool): Deprecated. Now the supercell image
                is always included in the returned data.

        Returns:
            [NearestNeighbor] where Nearest Neighbor is a named tuple containing
            (site, distance, index, image).
        """
        outer = self.get_sites_in_sphere(origin, r + dr, include_index=include_index, include_image=include_image)
        inner = r - dr
        return [t for t in outer if t.nn_distance > inner]

    def get_sorted_structure(self, key: Callable | None = None, reverse: bool = False) -> IStructure | Structure:
        """
        Get a sorted copy of the structure. The parameters have the same
        meaning as in list.sort. By default, sites are sorted by the
        electronegativity of the species.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        sites = sorted(self, key=key, reverse=reverse)
        return type(self).from_sites(sites, charge=self._charge)

    def get_reduced_structure(self, reduction_algo: Literal["niggli", "LLL"] = "niggli") -> IStructure | Structure:
        """
        Get a reduced structure.

        Args:
            reduction_algo ("niggli" | "LLL"): The lattice reduction algorithm to use.
                Defaults to "niggli".
        """
        if reduction_algo == "niggli":
            reduced_latt = self._lattice.get_niggli_reduced_lattice()
        elif reduction_algo == "LLL":
            reduced_latt = self._lattice.get_lll_reduced_lattice()
        else:
            raise ValueError(f"Invalid reduction algo : {reduction_algo}")

        if reduced_latt != self.lattice:
            return self.__class__(
                reduced_latt,
                self.species_and_occu,
                self.cart_coords,  # type: ignore
                coords_are_cartesian=True,
                to_unit_cell=True,
                site_properties=self.site_properties,
                labels=self.labels,
                charge=self._charge,
            )
        return self.copy()

    def copy(self, site_properties=None, sanitize=False):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Args:
            site_properties (dict): Properties to add or override. The
                properties are specified in the same way as the constructor,
                i.e., as a dict of the form {property: [values]}. The
                properties should be in the order of the *original* structure
                if you are performing sanitization.
            sanitize (bool): If True, this method will return a sanitized
                structure. Sanitization performs a few things: (i) The sites are
                sorted by electronegativity, (ii) a LLL lattice reduction is
                carried out to obtain a relatively orthogonalized cell,
                (iii) all fractional coords for sites are mapped into the
                unit cell.

        Returns:
            A copy of the Structure, with optionally new site_properties and
            optionally sanitized.
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)
        if not sanitize:
            return self.__class__(
                self._lattice,
                self.species_and_occu,
                self.frac_coords,
                charge=self._charge,
                site_properties=props,
                labels=self.labels,
            )
        reduced_latt = self._lattice.get_lll_reduced_lattice()
        new_sites = []
        for idx, site in enumerate(self):
            frac_coords = reduced_latt.get_fractional_coords(site.coords)
            site_props = {}
            for prop, val in props.items():
                site_props[prop] = val[idx]
            new_sites.append(
                PeriodicSite(
                    site.species,
                    frac_coords,
                    reduced_latt,
                    to_unit_cell=True,
                    properties=site_props,
                    label=site.label,
                    skip_checks=True,
                )
            )
        new_sites = sorted(new_sites)
        return type(self).from_sites(new_sites, charge=self._charge)

    def interpolate(
        self,
        end_structure: IStructure | Structure,
        nimages: int | Iterable = 10,
        interpolate_lattices: bool = False,
        pbc: bool = True,
        autosort_tol: float = 0,
    ) -> list[IStructure | Structure]:
        """
        Interpolate between this structure and end_structure. Useful for
        construction of NEB inputs.

        Args:
            end_structure (Structure): structure to interpolate between this
                structure and end.
            nimages (int,list): No. of interpolation images or a list of
                interpolation images. Defaults to 10 images.
            interpolate_lattices (bool): Whether to interpolate the lattices.
                Interpolates the lengths and angles (rather than the matrix)
                so orientation may be affected.
            pbc (bool): Whether to use periodic boundary conditions to find
                the shortest path between endpoints.
            autosort_tol (float): A distance tolerance in angstrom in
                which to automatically sort end_structure to match to the
                closest points in this particular structure. This is usually
                what you want in a NEB calculation. 0 implies no sorting.
                Otherwise, a 0.5 value usually works pretty well.

        Returns:
            List of interpolated structures. The starting and ending
            structures included as the first and last structures respectively.
            A total of (nimages + 1) structures are returned.
        """
        # Check length of structures
        if len(self) != len(end_structure):
            raise ValueError("Structures have different lengths!")

        if not (interpolate_lattices or self.lattice == end_structure.lattice):
            raise ValueError("Structures with different lattices!")

        images = np.arange(nimages + 1) / nimages if not isinstance(nimages, collections.abc.Iterable) else nimages

        # Check that both structures have the same species
        for i, site in enumerate(self):
            if site.species != end_structure[i].species:
                raise ValueError(f"Different species!\nStructure 1:\n{self}\nStructure 2\n{end_structure}")

        start_coords = np.array(self.frac_coords)
        end_coords = np.array(end_structure.frac_coords)

        if autosort_tol:
            dist_matrix = self.lattice.get_all_distances(start_coords, end_coords)
            site_mappings: dict[int, list[int]] = collections.defaultdict(list)
            unmapped_start_ind = []
            for i, row in enumerate(dist_matrix):
                ind = np.where(row < autosort_tol)[0]
                if len(ind) == 1:
                    site_mappings[i].append(ind[0])
                else:
                    unmapped_start_ind.append(i)

            if len(unmapped_start_ind) > 1:
                raise ValueError(f"Unable to reliably match structures with {autosort_tol = }, {unmapped_start_ind = }")

            sorted_end_coords = np.zeros_like(end_coords)
            matched = []
            for i, j in site_mappings.items():
                if len(j) > 1:
                    raise ValueError(
                        f"Unable to reliably match structures with auto_sort_tol = {autosort_tol}. "
                        "More than one site match!"
                    )
                sorted_end_coords[i] = end_coords[j[0]]
                matched.append(j[0])

            if len(unmapped_start_ind) == 1:
                i = unmapped_start_ind[0]
                j = list(set(range(len(start_coords))) - set(matched))[0]  # type: ignore
                sorted_end_coords[i] = end_coords[j]

            end_coords = sorted_end_coords

        vec = end_coords - start_coords
        if pbc:
            vec[:, self.pbc] -= np.round(vec[:, self.pbc])
        sp = self.species_and_occu
        structs = []

        if interpolate_lattices:
            # interpolate lattice matrices using polar decomposition
            from scipy.linalg import polar

            # u is a unitary rotation, p is stretch
            u, p = polar(np.dot(end_structure.lattice.matrix.T, np.linalg.inv(self.lattice.matrix.T)))
            lvec = p - np.identity(3)
            lstart = self.lattice.matrix.T

        for x in images:
            if interpolate_lattices:
                l_a = np.dot(np.identity(3) + x * lvec, lstart).T
                lat = Lattice(l_a)
            else:
                lat = self.lattice
            fcoords = start_coords + x * vec
            structs.append(self.__class__(lat, sp, fcoords, site_properties=self.site_properties))
        return structs

    def get_miller_index_from_site_indexes(self, site_ids, round_dp=4, verbose=True):
        """
        Get the Miller index of a plane from a set of sites indexes.

        A minimum of 3 sites are required. If more than 3 sites are given
        the best plane that minimises the distance to all points will be
        calculated.

        Args:
            site_ids (list of int): A list of site indexes to consider. A
                minimum of three site indexes are required. If more than three
                sites are provided, the best plane that minimises the distance
                to all sites will be calculated.
            round_dp (int, optional): The number of decimal places to round the
                miller index to.
            verbose (bool, optional): Whether to print warnings.

        Returns:
            (tuple): The Miller index.
        """
        return self.lattice.get_miller_index_from_coords(
            self.frac_coords[site_ids],
            coords_are_cartesian=False,
            round_dp=round_dp,
            verbose=verbose,
        )

    def get_primitive_structure(
        self, tolerance: float = 0.25, use_site_props: bool = False, constrain_latt: list | dict | None = None
    ):
        """
        This finds a smaller unit cell than the input. Sometimes it doesn"t
        find the smallest possible one, so this method is recursively called
        until it is unable to find a smaller cell.

        NOTE: if the tolerance is greater than 1/2 the minimum inter-site
        distance in the primitive cell, the algorithm will reject this lattice.

        Args:
            tolerance (float): Tolerance for each coordinate of a
                particular site in Angstroms. For example, [0.1, 0, 0.1] in cartesian
                coordinates will be considered to be on the same coordinates
                as [0, 0, 0] for a tolerance of 0.25. Defaults to 0.25.
            use_site_props (bool): Whether to account for site properties in
                differentiating sites.
            constrain_latt (list/dict): List of lattice parameters we want to
                preserve, e.g. ["alpha", "c"] or dict with the lattice
                parameter names as keys and values we want the parameters to
                be e.g. {"alpha": 90, "c": 2.5}.

        Returns:
            The most primitive structure found.
        """
        if constrain_latt is None:
            constrain_latt = []

        def site_label(site):
            if not use_site_props:
                return site.species_string
            parts = [site.species_string]
            for key in sorted(site.properties):
                parts.append(f"{key}={site.properties[key]}")
            return ", ".join(parts)

        # group sites by species string
        sites = sorted(self._sites, key=site_label)

        grouped_sites = [list(a[1]) for a in itertools.groupby(sites, key=site_label)]
        grouped_fcoords = [np.array([s.frac_coords for s in g]) for g in grouped_sites]

        # min_vecs are approximate periodicities of the cell. The exact
        # periodicities from the supercell matrices are checked against these
        # first
        min_fcoords = min(grouped_fcoords, key=lambda x: len(x))
        min_vecs = min_fcoords - min_fcoords[0]

        # fractional tolerance in the supercell
        super_ftol = np.divide(tolerance, self.lattice.abc)
        super_ftol_2 = super_ftol * 2

        def pbc_coord_intersection(fc1, fc2, tol):
            """
            Returns the fractional coords in fc1 that have coordinates
            within tolerance to some coordinate in fc2.
            """
            d = fc1[:, None, :] - fc2[None, :, :]
            d -= np.round(d)
            np.abs(d, d)
            return fc1[np.any(np.all(d < tol, axis=-1), axis=-1)]

        # here we reduce the number of min_vecs by enforcing that every
        # vector in min_vecs approximately maps each site onto a similar site.
        # The subsequent processing is O(fu^3 * min_vecs) = O(n^4) if we do no
        # reduction.
        # This reduction is O(n^3) so usually is an improvement. Using double
        # the tolerance because both vectors are approximate
        for g in sorted(grouped_fcoords, key=lambda x: len(x)):
            for f in g:
                min_vecs = pbc_coord_intersection(min_vecs, g - f, super_ftol_2)

        def get_hnf(fu):
            """
            Returns all possible distinct supercell matrices given a
            number of formula units in the supercell. Batches the matrices
            by the values in the diagonal (for less numpy overhead).
            Computational complexity is O(n^3), and difficult to improve.
            Might be able to do something smart with checking combinations of a
            and b first, though unlikely to reduce to O(n^2).
            """

            def factors(n):
                for i in range(1, n + 1):
                    if n % i == 0:
                        yield i

            for det in factors(fu):
                if det == 1:
                    continue
                for a in factors(det):
                    for e in factors(det // a):
                        g = det // a // e
                        yield det, np.array(
                            [
                                [[a, b, c], [0, e, f], [0, 0, g]]
                                for b, c, f in itertools.product(range(a), range(a), range(e))
                            ]
                        )

        # we can't let sites match to their neighbors in the supercell
        grouped_non_nbrs = []
        for gfcoords in grouped_fcoords:
            fdist = gfcoords[None, :, :] - gfcoords[:, None, :]
            fdist -= np.round(fdist)
            np.abs(fdist, fdist)
            non_nbrs = np.any(fdist > 2 * super_ftol[None, None, :], axis=-1)
            # since we want sites to match to themselves
            np.fill_diagonal(non_nbrs, True)
            grouped_non_nbrs.append(non_nbrs)

        num_fu = functools.reduce(math.gcd, map(len, grouped_sites))
        for size, ms in get_hnf(num_fu):
            inv_ms = np.linalg.inv(ms)

            # find sets of lattice vectors that are present in min_vecs
            dist = inv_ms[:, :, None, :] - min_vecs[None, None, :, :]
            dist -= np.round(dist)
            np.abs(dist, dist)
            is_close = np.all(dist < super_ftol, axis=-1)
            any_close = np.any(is_close, axis=-1)
            inds = np.all(any_close, axis=-1)

            for inv_m, m in zip(inv_ms[inds], ms[inds]):
                new_m = np.dot(inv_m, self.lattice.matrix)
                ftol = np.divide(tolerance, np.sqrt(np.sum(new_m**2, axis=1)))

                valid = True
                new_coords = []
                new_sp = []
                new_props = collections.defaultdict(list)
                new_labels = []
                for gsites, gfcoords, non_nbrs in zip(grouped_sites, grouped_fcoords, grouped_non_nbrs):
                    all_frac = np.dot(gfcoords, m)

                    # calculate grouping of equivalent sites, represented by
                    # adjacency matrix
                    fdist = all_frac[None, :, :] - all_frac[:, None, :]
                    fdist = np.abs(fdist - np.round(fdist))
                    close_in_prim = np.all(fdist < ftol[None, None, :], axis=-1)
                    groups = np.logical_and(close_in_prim, non_nbrs)

                    # check that groups are correct
                    if not np.all(np.sum(groups, axis=0) == size):
                        valid = False
                        break

                    # check that groups are all cliques
                    for g in groups:
                        if not np.all(groups[g][:, g]):
                            valid = False
                            break
                    if not valid:
                        break

                    # add the new sites, averaging positions
                    added = np.zeros(len(gsites))
                    new_fcoords = all_frac % 1
                    for i, group in enumerate(groups):
                        if not added[i]:
                            added[group] = True
                            inds = np.where(group)[0]
                            coords = new_fcoords[inds[0]]
                            for n, j in enumerate(inds[1:]):
                                offset = new_fcoords[j] - coords
                                coords += (offset - np.round(offset)) / (n + 2)
                            new_sp.append(gsites[inds[0]].species)
                            for k in gsites[inds[0]].properties:
                                new_props[k].append(gsites[inds[0]].properties[k])
                            new_labels.append(gsites[inds[0]].label)
                            new_coords.append(coords)

                if valid:
                    inv_m = np.linalg.inv(m)
                    new_latt = Lattice(np.dot(inv_m, self.lattice.matrix))
                    struct = Structure(
                        new_latt,
                        new_sp,
                        new_coords,
                        site_properties=new_props,
                        labels=new_labels,
                        coords_are_cartesian=False,
                    )

                    # Default behavior
                    p = struct.get_primitive_structure(
                        tolerance=tolerance,
                        use_site_props=use_site_props,
                        constrain_latt=constrain_latt,
                    ).get_reduced_structure()
                    if not constrain_latt:
                        return p

                    # Only return primitive structures that
                    # satisfy the restriction condition
                    p_latt, s_latt = p.lattice, self.lattice
                    if type(constrain_latt).__name__ == "list":
                        if all(getattr(p_latt, pp) == getattr(s_latt, pp) for pp in constrain_latt):
                            return p
                    elif type(constrain_latt).__name__ == "dict" and all(
                        getattr(p_latt, pp) == constrain_latt[pp] for pp in constrain_latt
                    ):
                        return p

        return self.copy()

    def __repr__(self):
        outs = ["Structure Summary", repr(self.lattice)]
        if self._charge:
            outs.append(f"Overall Charge: {self._charge:+}")
        for site in self:
            outs.append(repr(site))
        return "\n".join(outs)

    def __str__(self):
        outs = [
            f"Full Formula ({self.composition.formula})",
            f"Reduced Formula: {self.composition.reduced_formula}",
        ]

        def to_str(x):
            return f"{x:>10.6f}"

        outs.append("abc   : " + " ".join(to_str(i) for i in self.lattice.abc))
        outs.append("angles: " + " ".join(to_str(i) for i in self.lattice.angles))
        outs.append("pbc   : " + " ".join(str(p).rjust(10) for p in self.lattice.pbc))
        if self._charge:
            outs.append(f"Overall Charge: {self._charge:+}")
        outs.append(f"Sites ({len(self)})")
        data = []
        props = self.site_properties
        keys = sorted(props)
        for i, site in enumerate(self):
            row = [str(i), site.species_string]
            row.extend([to_str(j) for j in site.frac_coords])
            for k in keys:
                row.append(props[k][i])
            data.append(row)
        outs.append(
            tabulate(
                data,
                headers=["#", "SP", "a", "b", "c", *keys],
            )
        )
        return "\n".join(outs)

    def get_orderings(self, mode: Literal["enum", "sqs"] = "enum", **kwargs) -> list[Structure]:
        """
        Returns list of orderings for a disordered structure. If structure
        does not contain disorder, the default structure is returned.

        Args:
            mode ("enum" | "sqs"): Either "enum" or "sqs". If enum,
                the enumlib will be used to return all distinct
                orderings. If sqs, mcsqs will be used to return
                an sqs structure.
            kwargs: kwargs passed to either
                pymatgen.command_line..enumlib_caller.EnumlibAdaptor
                or pymatgen.command_line.mcsqs_caller.run_mcsqs.
                For run_mcsqs, a default cluster search of 2 cluster interactions
                with 1NN distance and 3 cluster interactions with 2NN distance
                is set.

        Returns:
            List[Structure]
        """
        if self.is_ordered:
            return [self]
        if mode.startswith("enum"):
            from pymatgen.command_line.enumlib_caller import EnumlibAdaptor

            adaptor = EnumlibAdaptor(self, **kwargs)
            adaptor.run()
            return adaptor.structures
        if mode == "sqs":
            from pymatgen.command_line.mcsqs_caller import run_mcsqs

            if "clusters" not in kwargs:
                disordered_sites = [site for site in self if not site.is_ordered]
                subset_structure = Structure.from_sites(disordered_sites)
                dist_matrix = subset_structure.distance_matrix
                dists = sorted(set(dist_matrix.ravel()))
                unique_dists = []
                for idx in range(1, len(dists)):
                    if dists[idx] - dists[idx - 1] > 0.1:
                        unique_dists.append(dists[idx])
                clusters = {(i + 2): d + 0.01 for i, d in enumerate(unique_dists) if i < 2}
                kwargs["clusters"] = clusters
            return [run_mcsqs(self, **kwargs).bestsqs]
        raise ValueError("Invalid mode!")

    def as_dict(self, verbosity=1, fmt=None, **kwargs) -> dict[str, Any]:
        """
        Dict representation of Structure.

        Args:
            verbosity (int): Verbosity level. Default of 1 includes both
                direct and Cartesian coordinates for all sites, lattice
                parameters, etc. Useful for reading and for insertion into a
                database. Set to 0 for an extremely lightweight version
                that only includes sufficient information to reconstruct the
                object.
            fmt (str): Specifies a format for the dict. Defaults to None,
                which is the default format used in pymatgen. Other options
                include "abivars".
            **kwargs: Allow passing of other kwargs needed for certain
            formats, e.g., "abivars".

        Returns:
            JSON-serializable dict representation.
        """
        if fmt == "abivars":
            """Returns a dictionary with the ABINIT variables."""
            from pymatgen.io.abinit.abiobjects import structure_to_abivars

            return structure_to_abivars(self, **kwargs)

        latt_dict = self._lattice.as_dict(verbosity=verbosity)
        del latt_dict["@module"]
        del latt_dict["@class"]
        sites = []
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "charge": self.charge,
            "lattice": latt_dict,
        }
        for site in self:
            site_dict = site.as_dict(verbosity=verbosity)  # type: ignore[call-arg]
            del site_dict["lattice"]
            del site_dict["@module"]
            del site_dict["@class"]
            sites.append(site_dict)
        dct["sites"] = sites
        return dct

    def as_dataframe(self):
        """
        Create a Pandas dataframe of the sites. Structure-level attributes are stored in DataFrame.attrs.

        Example:
            Species    a    b             c    x             y             z  magmom
            0    (Si)  0.0  0.0  0.000000e+00  0.0  0.000000e+00  0.000000e+00       5
            1    (Si)  0.0  0.0  1.000000e-07  0.0 -2.217138e-07  3.135509e-07      -5
        """
        data = []
        site_properties = self.site_properties
        prop_keys = list(site_properties)
        for site in self:
            row = [site.species, *site.frac_coords, *site.coords]
            for k in prop_keys:
                row.append(site.properties.get(k))
            data.append(row)
        import pandas as pd

        df = pd.DataFrame(data, columns=["Species", "a", "b", "c", "x", "y", "z", *prop_keys])
        df.attrs["Reduced Formula"] = self.composition.reduced_formula
        df.attrs["Lattice"] = self.lattice
        return df

    @classmethod
    def from_dict(cls, d: dict[str, Any], fmt: Literal["abivars"] | None = None) -> Structure:
        """
        Reconstitute a Structure object from a dict representation of Structure
        created using as_dict().

        Args:
            d (dict): Dict representation of structure.
            fmt ('abivars' | None): Use structure_from_abivars() to parse the dict. Defaults to None.

        Returns:
            Structure object
        """
        if fmt == "abivars":
            from pymatgen.io.abinit.abiobjects import structure_from_abivars

            return structure_from_abivars(cls=cls, **d)

        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        charge = d.get("charge")
        return cls.from_sites(sites, charge=charge)

    def to(self, filename: str = "", fmt: str = "", **kwargs) -> str | None:
        """
        Outputs the structure to a file or string.

        Args:
            filename (str): If provided, output will be written to a file. If
                fmt is not specified, the format is determined from the
                filename. Defaults is None, i.e. string output.
            fmt (str): Format to output to. Defaults to JSON unless filename
                is provided. If fmt is specifies, it overrides whatever the
                filename is. Options include "cif", "poscar", "cssr", "json",
                "xsf", "mcsqs", "prismatic", "yaml", "fleur-inpgen".
                Non-case sensitive.
            **kwargs: Kwargs passthru to relevant methods. E.g., This allows
                the passing of parameters like symprec to the
                CifWriter.__init__ method for generation of symmetric cifs.

        Returns:
            (str) if filename is None. None otherwise.
        """
        fmt = fmt.lower()

        if fmt == "cif" or fnmatch(filename.lower(), "*.cif*"):
            from pymatgen.io.cif import CifWriter

            writer = CifWriter(self, **kwargs)
        elif fmt == "mcif" or fnmatch(filename.lower(), "*.mcif*"):
            from pymatgen.io.cif import CifWriter

            writer = CifWriter(self, write_magmoms=True, **kwargs)
        elif fmt == "poscar" or fnmatch(filename, "*POSCAR*"):
            from pymatgen.io.vasp import Poscar

            writer = Poscar(self, **kwargs)
        elif fmt == "cssr" or fnmatch(filename.lower(), "*.cssr*"):
            from pymatgen.io.cssr import Cssr

            writer = Cssr(self)  # type: ignore
        elif fmt == "json" or fnmatch(filename.lower(), "*.json*"):
            dct = json.dumps(self.as_dict())
            if filename:
                with zopen(filename, "wt") as file:
                    file.write(dct)
            return dct
        elif fmt == "xsf" or fnmatch(filename.lower(), "*.xsf*"):
            from pymatgen.io.xcrysden import XSF

            string = XSF(self).to_str()
            if filename:
                with zopen(filename, "wt", encoding="utf8") as file:
                    file.write(string)
            return string
        elif (
            fmt == "mcsqs"
            or fnmatch(filename, "*rndstr.in*")
            or fnmatch(filename, "*lat.in*")
            or fnmatch(filename, "*bestsqs*")
        ):
            from pymatgen.io.atat import Mcsqs

            string = Mcsqs(self).to_str()
            if filename:
                with zopen(filename, "wt", encoding="ascii") as file:
                    file.write(string)
            return string
        elif fmt == "prismatic" or fnmatch(filename, "*prismatic*"):
            from pymatgen.io.prismatic import Prismatic

            return Prismatic(self).to_str()
        elif fmt == "yaml" or fnmatch(filename, "*.yaml*") or fnmatch(filename, "*.yml*"):
            yaml = YAML()
            if filename:
                with zopen(filename, "wt") as file:
                    yaml.dump(self.as_dict(), file)
                return None
            sio = StringIO()
            yaml.dump(self.as_dict(), sio)
            return sio.getvalue()
        # fleur support implemented in external namespace pkg https://github.com/JuDFTteam/pymatgen-io-fleur
        elif fmt == "fleur-inpgen" or fnmatch(filename, "*.in*"):
            from pymatgen.io.fleur import FleurInput

            writer = FleurInput(self, **kwargs)
        elif fmt == "res" or fnmatch(filename, "*.res"):
            from pymatgen.io.res import ResIO

            string = ResIO.structure_to_str(self)
            if filename:
                with zopen(filename, "wt", encoding="utf8") as file:
                    file.write(string)
                return None
            return string
        else:
            if fmt == "":
                raise ValueError(f"Format not specified and could not infer from {filename=}")
            raise ValueError(f"Invalid format={fmt!r}")

        if filename:
            writer.write_file(filename)
            return None
        return str(writer)

    @classmethod
    def from_str(
        cls,
        input_string: str,
        fmt: Literal["cif", "poscar", "cssr", "json", "yaml", "xsf", "mcsqs", "res"],
        primitive: bool = False,
        sort: bool = False,
        merge_tol: float = 0.0,
        **kwargs,
    ) -> Structure | IStructure:
        """
        Reads a structure from a string.

        Args:
            input_string (str): String to parse.
            fmt (str): A file format specification. One of "cif", "poscar", "cssr",
                "json", "yaml", "xsf", "mcsqs".
            primitive (bool): Whether to find a primitive cell. Defaults to
                False.
            sort (bool): Whether to sort the sites in accordance to the default
                ordering criteria, i.e., electronegativity.
            merge_tol (float): If this is some positive number, sites that
                are within merge_tol from each other will be merged. Usually
                0.01 should be enough to deal with common numerical issues.
            **kwargs: Passthrough to relevant parser.

        Returns:
            IStructure | Structure
        """
        fmt_low = fmt.lower()
        if fmt_low == "cif":
            from pymatgen.io.cif import CifParser

            parser = CifParser.from_str(input_string, **kwargs)
            struct = parser.get_structures(primitive=primitive)[0]
        elif fmt_low == "poscar":
            from pymatgen.io.vasp import Poscar

            struct = Poscar.from_str(input_string, False, read_velocities=False, **kwargs).structure
        elif fmt_low == "cssr":
            from pymatgen.io.cssr import Cssr

            cssr = Cssr.from_str(input_string, **kwargs)
            struct = cssr.structure
        elif fmt_low == "json":
            d = json.loads(input_string)
            struct = Structure.from_dict(d)
        elif fmt_low == "yaml":
            yaml = YAML()
            d = yaml.load(input_string)
            struct = Structure.from_dict(d)
        elif fmt_low == "xsf":
            from pymatgen.io.xcrysden import XSF

            struct = XSF.from_str(input_string, **kwargs).structure
        elif fmt_low == "mcsqs":
            from pymatgen.io.atat import Mcsqs

            struct = Mcsqs.structure_from_str(input_string, **kwargs)
        # fleur support implemented in external namespace pkg https://github.com/JuDFTteam/pymatgen-io-fleur
        elif fmt == "fleur-inpgen":
            from pymatgen.io.fleur import FleurInput

            struct = FleurInput.from_string(input_string, inpgen_input=True, **kwargs).structure
        elif fmt == "fleur":
            from pymatgen.io.fleur import FleurInput

            struct = FleurInput.from_string(input_string, inpgen_input=False).structure
        elif fmt == "res":
            from pymatgen.io.res import ResIO

            struct = ResIO.structure_from_str(input_string, **kwargs)
        else:
            raise ValueError(f"Unrecognized format `{fmt}`!")

        if sort:
            struct = struct.get_sorted_structure()
        if merge_tol:
            struct.merge_sites(merge_tol)
        return cls.from_sites(struct)

    @classmethod
    def from_file(cls, filename, primitive=False, sort=False, merge_tol=0.0, **kwargs) -> Structure | IStructure:
        """
        Reads a structure from a file. For example, anything ending in
        a "cif" is assumed to be a Crystallographic Information Format file.
        Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT,
        vasprun.xml, CSSR, Netcdf and pymatgen's JSON-serialized structures.

        Args:
            filename (str): The filename to read from.
            primitive (bool): Whether to convert to a primitive cell. Only available for cifs. Defaults to False.
            sort (bool): Whether to sort sites. Default to False.
            merge_tol (float): If this is some positive number, sites that are within merge_tol from each other will be
                merged. Usually 0.01 should be enough to deal with common numerical issues.
            kwargs: Passthrough to relevant reader. E.g. if the file has CIF format, the kwargs will be passed
                through to CifParser.

        Returns:
            Structure.
        """
        filename = str(filename)
        if filename.endswith(".nc"):
            # Read Structure from a netcdf file.
            from pymatgen.io.abinit.netcdf import structure_from_ncdata

            struct = structure_from_ncdata(filename, cls=cls)
            if sort:
                struct = struct.get_sorted_structure()
            return struct

        from pymatgen.io.exciting import ExcitingInput
        from pymatgen.io.lmto import LMTOCtrl
        from pymatgen.io.vasp import Chgcar, Vasprun

        fname = os.path.basename(filename)
        with zopen(filename, "rt", errors="replace") as f:
            contents = f.read()
        if fnmatch(fname.lower(), "*.cif*") or fnmatch(fname.lower(), "*.mcif*"):
            return cls.from_str(contents, fmt="cif", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)
        if fnmatch(fname, "*POSCAR*") or fnmatch(fname, "*CONTCAR*") or fnmatch(fname, "*.vasp"):
            struct = cls.from_str(contents, fmt="poscar", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)

        elif fnmatch(fname, "CHGCAR*") or fnmatch(fname, "LOCPOT*"):
            struct = Chgcar.from_file(filename, **kwargs).structure
        elif fnmatch(fname, "vasprun*.xml*"):
            struct = Vasprun(filename, **kwargs).final_structure
        elif fnmatch(fname.lower(), "*.cssr*"):
            return cls.from_str(contents, fmt="cssr", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)
        elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
            return cls.from_str(contents, fmt="json", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)
        elif fnmatch(fname, "*.yaml*") or fnmatch(fname, "*.yml*"):
            return cls.from_str(contents, fmt="yaml", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)
        elif fnmatch(fname, "*.xsf"):
            return cls.from_str(contents, fmt="xsf", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)
        elif fnmatch(fname, "input*.xml"):
            return ExcitingInput.from_file(fname, **kwargs).structure
        elif fnmatch(fname, "*rndstr.in*") or fnmatch(fname, "*lat.in*") or fnmatch(fname, "*bestsqs*"):
            return cls.from_str(contents, fmt="mcsqs", primitive=primitive, sort=sort, merge_tol=merge_tol, **kwargs)
        elif fnmatch(fname, "CTRL*"):
            return LMTOCtrl.from_file(filename=filename, **kwargs).structure
        elif fnmatch(fname, "inp*.xml") or fnmatch(fname, "*.in*") or fnmatch(fname, "inp_*"):
            from pymatgen.io.fleur import FleurInput

            struct = FleurInput.from_file(filename, **kwargs).structure
        elif fnmatch(fname, "*.res"):
            from pymatgen.io.res import ResIO

            struct = ResIO.structure_from_file(filename, **kwargs)
        else:
            raise ValueError("Unrecognized file extension!")
        if sort:
            struct = struct.get_sorted_structure()
        if merge_tol:
            struct.merge_sites(merge_tol)

        struct.__class__ = cls
        return struct


class IMolecule(SiteCollection, MSONable):
    """
    Basic immutable Molecule object without periodicity. Essentially a
    sequence of sites. IMolecule is made to be immutable so that they can
    function as keys in a dict. For a mutable molecule,
    use the :class:Molecule.

    Molecule extends Sequence and Hashable, which means that in many cases,
    it can be used like any Python sequence. Iterating through a molecule is
    equivalent to going through the sites in sequence.
    """

    def __init__(
        self,
        species: Sequence[CompositionLike],
        coords: Sequence[ArrayLike],
        charge: float = 0.0,
        spin_multiplicity: int | None = None,
        validate_proximity: bool = False,
        site_properties: dict | None = None,
        labels: Sequence[str | None] | None = None,
        charge_spin_check: bool = True,
    ) -> None:
        """
        Create a Molecule.

        Args:
            species: list of atomic species. Possible kinds of input include a
                list of dict of elements/species and occupancies, a List of
                elements/specie specified as actual Element/Species, Strings
                ("Fe", "Fe2+") or atomic numbers (1,56).
            coords (3x1 array): list of Cartesian coordinates of each species.
            charge (float): Charge for the molecule. Defaults to 0.
            spin_multiplicity (int): Spin multiplicity for molecule.
                Defaults to None, which means that the spin multiplicity is
                set to 1 if the molecule has no unpaired electrons and to 2
                if there are unpaired electrons.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 1 Ang apart. Defaults to False.
            site_properties (dict): Properties associated with the sites as
                a dict of sequences, e.g., {"magmom":[5,5,5,5]}. The
                sequences have to be the same length as the atomic species
                and fractional_coords. Defaults to None for no properties.
            labels (list[str]): Labels associated with the sites as a
                list of strings, e.g. ['Li1', 'Li2']. Must have the same
                length as the species and fractional coords. Defaults to
                None for no labels.
            charge_spin_check (bool): Whether to check that the charge and
                spin multiplicity are compatible with each other. Defaults
                to True.
        """
        if len(species) != len(coords):
            raise StructureError(
                f"len(species) != len(coords) ({len(species)} != {len(coords)}). List of atomic species must "
                "have same length as list of fractional coordinates."
            )

        self._charge_spin_check = charge_spin_check

        sites: list[Site] = []
        for idx in range(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[idx] for k, v in site_properties.items()}
            sites.append(Site(species[idx], coords[idx], properties=prop))

        self._sites = tuple(sites)
        if validate_proximity and not self.is_valid():
            raise StructureError("Molecule contains sites that are less than 0.01 Angstrom apart!")

        self._charge = charge
        nelectrons = 0.0
        for site in sites:
            for sp, amt in site.species.items():
                if not isinstance(sp, DummySpecies):
                    nelectrons += sp.Z * amt
        nelectrons -= charge
        self._nelectrons = nelectrons
        if spin_multiplicity:
            if charge_spin_check and (nelectrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    f"Charge of {self._charge} and spin multiplicity of {spin_multiplicity} "
                    "is not possible for this molecule!"
                )
            self._spin_multiplicity = spin_multiplicity
        else:
            self._spin_multiplicity = 1 if nelectrons % 2 == 0 else 2

    @property
    def charge(self) -> float:
        """Charge of molecule."""
        return self._charge

    @property
    def spin_multiplicity(self) -> float:
        """Spin multiplicity of molecule."""
        return self._spin_multiplicity

    @property
    def nelectrons(self) -> float:
        """Number of electrons in the molecule."""
        return self._nelectrons

    @property
    def center_of_mass(self) -> np.ndarray:
        """Center of mass of molecule."""
        center = np.zeros(3)
        total_weight: float = 0
        for site in self:
            wt = site.species.weight
            center += site.coords * wt
            total_weight += wt
        return center / total_weight

    @property
    def sites(self) -> tuple[Site, ...]:
        """Returns a tuple of sites in the Molecule."""
        return self._sites

    @classmethod
    def from_sites(
        cls,
        sites: Sequence[Site],
        charge: float = 0,
        spin_multiplicity: int | None = None,
        validate_proximity: bool = False,
        charge_spin_check: bool = True,
    ) -> IMolecule | Molecule:
        """
        Convenience constructor to make a Molecule from a list of sites.

        Args:
            sites ([Site]): Sequence of Sites.
            charge (int): Charge of molecule. Defaults to 0.
            spin_multiplicity (int): Spin multicipity. Defaults to None,
                in which it is determined automatically.
            validate_proximity (bool): Whether to check that atoms are too
                close.
            charge_spin_check (bool): Whether to check that the charge and
                spin multiplicity are compatible with each other. Defaults
                to True.
        """
        props = collections.defaultdict(list)
        for site in sites:
            for k, v in site.properties.items():
                props[k].append(v)
        labels = [site.label for site in sites]
        return cls(
            [site.species for site in sites],
            [site.coords for site in sites],
            charge=charge,
            spin_multiplicity=spin_multiplicity,
            validate_proximity=validate_proximity,
            site_properties=props,
            labels=labels,
            charge_spin_check=charge_spin_check,
        )

    def break_bond(self, ind1: int, ind2: int, tol: float = 0.2) -> tuple[IMolecule | Molecule, ...]:
        """
        Returns two molecules based on breaking the bond between atoms at index
        ind1 and ind2.

        Args:
            ind1 (int): Index of first site.
            ind2 (int): Index of second site.
            tol (float): Relative tolerance to test. Basically, the code
                checks if the distance between the sites is less than (1 +
                tol) * typical bond distances. Defaults to 0.2, i.e.,
                20% longer.

        Returns:
            Two Molecule objects representing the two clusters formed from
            breaking the bond.
        """
        clusters = [[self._sites[ind1]], [self._sites[ind2]]]

        sites = [site for i, site in enumerate(self._sites) if i not in (ind1, ind2)]

        def belongs_to_cluster(site, cluster):
            return any(CovalentBond.is_bonded(site, test_site, tol=tol) for test_site in cluster)

        while len(sites) > 0:
            unmatched = []
            for site in sites:
                for cluster in clusters:
                    if belongs_to_cluster(site, cluster):
                        cluster.append(site)
                        break
                else:
                    unmatched.append(site)

            if len(unmatched) == len(sites):
                raise ValueError("Not all sites are matched!")
            sites = unmatched

        return tuple(type(self).from_sites(cluster) for cluster in clusters)

    def get_covalent_bonds(self, tol: float = 0.2) -> list[CovalentBond]:
        """
        Determines the covalent bonds in a molecule.

        Args:
            tol (float): The tol to determine bonds in a structure. See
                CovalentBond.is_bonded.

        Returns:
            List of bonds
        """
        bonds = []
        for site1, site2 in itertools.combinations(self._sites, 2):
            if CovalentBond.is_bonded(site1, site2, tol):
                bonds.append(CovalentBond(site1, site2))
        return bonds

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("charge", "spin_multiplicity", "sites")

        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented

        other = cast(IMolecule, other)

        if len(self) != len(other):
            return False
        if self.charge != other.charge:
            return False
        if self.spin_multiplicity != other.spin_multiplicity:
            return False
        return all(site in other for site in self)

    def get_zmatrix(self):
        """Returns a z-matrix representation of the molecule."""
        # TODO: allow more z-matrix conventions for element/site description

        output = []
        output_var = []
        for idx, site in enumerate(self._sites):
            if idx == 0:
                output.append(f"{site.specie}")
            elif idx == 1:
                nn = self._find_nn_pos_before_site(idx)
                bond_length = self.get_distance(idx, nn[0])
                output.append(f"{self._sites[idx].specie} {nn[0] + 1} B{idx}")
                output_var.append(f"B{idx}={bond_length:.6f}")
            elif idx == 2:
                nn = self._find_nn_pos_before_site(idx)
                bond_length = self.get_distance(idx, nn[0])
                angle = self.get_angle(idx, nn[0], nn[1])
                output.append(f"{self._sites[idx].specie} {nn[0] + 1} B{idx} {nn[1] + 1} A{idx}")
                output_var.append(f"B{idx}={bond_length:.6f}")
                output_var.append(f"A{idx}={angle:.6f}")
            else:
                nn = self._find_nn_pos_before_site(idx)
                bond_length = self.get_distance(idx, nn[0])
                angle = self.get_angle(idx, nn[0], nn[1])
                dih = self.get_dihedral(idx, nn[0], nn[1], nn[2])
                output.append(f"{self._sites[idx].specie} {nn[0] + 1} B{idx} {nn[1] + 1} A{idx} {nn[2] + 1} D{idx}")
                output_var.append(f"B{idx}={bond_length:.6f}")
                output_var.append(f"A{idx}={angle:.6f}")
                output_var.append(f"D{idx}={dih:.6f}")
        return "\n".join(output) + "\n\n" + "\n".join(output_var)

    def _find_nn_pos_before_site(self, site_idx):
        """Returns index of nearest neighbor atoms."""
        all_dist = [(self.get_distance(site_idx, idx), idx) for idx in range(site_idx)]
        all_dist = sorted(all_dist, key=lambda x: x[0])
        return [d[1] for d in all_dist]

    def __hash__(self) -> int:
        # For now, just use the composition hash code.
        return hash(self.composition)

    def __repr__(self):
        return "Molecule Summary\n" + "\n".join(map(repr, self))

    def __str__(self):
        outs = [
            f"Full Formula ({self.composition.formula})",
            "Reduced Formula: " + self.composition.reduced_formula,
            f"Charge = {self._charge}, Spin Mult = {self._spin_multiplicity}",
            f"Sites ({len(self)})",
        ]
        for idx, site in enumerate(self):
            outs.append(f"{idx} {site.species_string} {' '.join([f'{coord:0.6f}'.rjust(12) for coord in site.coords])}")
        return "\n".join(outs)

    def as_dict(self):
        """JSON-serializable dict representation of Molecule."""
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
            "sites": [],
        }
        for site in self:
            site_dict = site.as_dict()
            del site_dict["@module"]
            del site_dict["@class"]
            d["sites"].append(site_dict)
        return d

    @classmethod
    def from_dict(cls, d) -> IMolecule | Molecule:
        """
        Reconstitute a Molecule object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of Molecule.

        Returns:
            Molecule object
        """
        sites = [Site.from_dict(sd) for sd in d["sites"]]
        charge = d.get("charge", 0)
        spin_multiplicity = d.get("spin_multiplicity")
        return cls.from_sites(sites, charge=charge, spin_multiplicity=spin_multiplicity)

    def get_distance(self, i: int, j: int) -> float:
        """
        Get distance between site i and j.

        Args:
            i (int): Index of first site
            j (int): Index of second site

        Returns:
            Distance between the two sites.
        """
        return self[i].distance(self[j])

    def get_sites_in_sphere(self, pt: ArrayLike, r: float) -> list[Neighbor]:
        """
        Find all sites within a sphere from a point.

        Args:
            pt (3x1 array): Cartesian coordinates of center of sphere
            r (float): Radius of sphere.

        Returns:
            Neighbor
        """
        neighbors = []
        for i, site in enumerate(self._sites):
            dist = site.distance_from_point(pt)
            if dist <= r:
                neighbors.append(Neighbor(site.species, site.coords, site.properties, dist, i))
        return neighbors

    def get_neighbors(self, site: Site, r: float) -> list[Neighbor]:
        """
        Get all neighbors to a site within a sphere of radius r. Excludes the
        site itself.

        Args:
            site (Site): Site at the center of the sphere.
            r (float): Radius of sphere.

        Returns:
            Neighbor
        """
        nns = self.get_sites_in_sphere(site.coords, r)
        return [nn for nn in nns if nn != site]

    def get_neighbors_in_shell(self, origin: ArrayLike, r: float, dr: float) -> list[Neighbor]:
        """
        Returns all sites in a shell centered on origin (coords) between radii
        r-dr and r+dr.

        Args:
            origin (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Inner radius of shell.
            dr (float): Width of shell.

        Returns:
            Neighbor
        """
        outer = self.get_sites_in_sphere(origin, r + dr)
        inner = r - dr
        return [nn for nn in outer if nn.nn_distance > inner]

    def get_boxed_structure(
        self,
        a: float,
        b: float,
        c: float,
        images: ArrayLike = (1, 1, 1),
        random_rotation: bool = False,
        min_dist: float = 1.0,
        cls=None,
        offset: ArrayLike | None = None,
        no_cross: bool = False,
        reorder: bool = True,
    ) -> IStructure | Structure:
        """
        Creates a Structure from a Molecule by putting the Molecule in the
        center of a orthorhombic box. Useful for creating Structure for
        calculating molecules using periodic codes.

        Args:
            a (float): a-lattice parameter.
            b (float): b-lattice parameter.
            c (float): c-lattice parameter.
            images: No. of boxed images in each direction. Defaults to
                (1, 1, 1), meaning single molecule with 1 lattice parameter
                in each direction.
            random_rotation (bool): Whether to apply a random rotation to
                each molecule. This jumbles all the molecules so that they
                are not exact images of each other.
            min_dist (float): The minimum distance that atoms should be from
                each other. This is only used if random_rotation is True.
                The randomized rotations are searched such that no two atoms
                are less than min_dist from each other.
            cls: The Structure class to instantiate (defaults to pymatgen
                structure)
            offset: Translation to offset molecule from center of mass coords
            no_cross: Whether to forbid molecule coords from extending beyond
                boundary of box.
            reorder: Whether to reorder the sites to be in electronegativity
                order.

        Returns:
            Structure containing molecule in a box.
        """
        if offset is None:
            offset = np.array([0, 0, 0])

        coords = np.array(self.cart_coords)
        x_range = max(coords[:, 0]) - min(coords[:, 0])
        y_range = max(coords[:, 1]) - min(coords[:, 1])
        z_range = max(coords[:, 2]) - min(coords[:, 2])

        if a <= x_range or b <= y_range or c <= z_range:
            raise ValueError("Box is not big enough to contain Molecule.")
        lattice = Lattice.from_parameters(a * images[0], b * images[1], c * images[2], 90, 90, 90)  # type: ignore
        nimages = images[0] * images[1] * images[2]  # type: ignore
        all_coords: list[ArrayLike] = []

        centered_coords = self.cart_coords - self.center_of_mass + offset

        for i, j, k in itertools.product(
            list(range(images[0])), list(range(images[1])), list(range(images[2]))  # type: ignore
        ):
            box_center = [(i + 0.5) * a, (j + 0.5) * b, (k + 0.5) * c]
            if random_rotation:
                while True:
                    op = SymmOp.from_origin_axis_angle(
                        (0, 0, 0),
                        axis=np.random.rand(3),
                        angle=random.uniform(-180, 180),
                    )
                    m = op.rotation_matrix
                    new_coords = np.dot(m, centered_coords.T).T + box_center
                    if no_cross:
                        x_max, x_min = max(new_coords[:, 0]), min(new_coords[:, 0])
                        y_max, y_min = max(new_coords[:, 1]), min(new_coords[:, 1])
                        z_max, z_min = max(new_coords[:, 2]), min(new_coords[:, 2])
                        if x_max > a or x_min < 0 or y_max > b or y_min < 0 or z_max > c or z_min < 0:
                            raise ValueError("Molecule crosses boundary of box")
                    if len(all_coords) == 0:
                        break
                    distances = lattice.get_all_distances(
                        lattice.get_fractional_coords(new_coords),
                        lattice.get_fractional_coords(all_coords),
                    )
                    if np.amin(distances) > min_dist:
                        break
            else:
                new_coords = centered_coords + box_center
                if no_cross:
                    x_max, x_min = max(new_coords[:, 0]), min(new_coords[:, 0])
                    y_max, y_min = max(new_coords[:, 1]), min(new_coords[:, 1])
                    z_max, z_min = max(new_coords[:, 2]), min(new_coords[:, 2])
                    if x_max > a or x_min < 0 or y_max > b or y_min < 0 or z_max > c or z_min < 0:
                        raise ValueError("Molecule crosses boundary of box")
            all_coords.extend(new_coords)
        sprops = {k: v * nimages for k, v in self.site_properties.items()}  # type: ignore

        if cls is None:
            cls = Structure

        if reorder:
            return cls(
                lattice,
                self.species * nimages,  # type: ignore
                all_coords,
                coords_are_cartesian=True,
                site_properties=sprops,
            ).get_sorted_structure()

        return cls(
            lattice,
            self.species * nimages,  # type: ignore
            coords,
            coords_are_cartesian=True,
            site_properties=sprops,
        )

    def get_centered_molecule(self) -> IMolecule | Molecule:
        """
        Returns a Molecule centered at the center of mass.

        Returns:
            Molecule centered with center of mass at origin.
        """
        center = self.center_of_mass
        new_coords = np.array(self.cart_coords) - center
        return self.__class__(
            self.species_and_occu,
            new_coords,
            charge=self._charge,
            spin_multiplicity=self._spin_multiplicity,
            site_properties=self.site_properties,
            charge_spin_check=self._charge_spin_check,
        )

    def to(self, filename: str = "", fmt: str = "") -> str | None:
        """
        Outputs the molecule to a file or string.

        Args:
            filename (str): If provided, output will be written to a file. If
                fmt is not specified, the format is determined from the
                filename. Defaults is None, i.e. string output.
            fmt (str): Format to output to. Defaults to JSON unless filename
                is provided. If fmt is specifies, it overrides whatever the
                filename is. Options include "xyz", "gjf", "g03", "json". If
                you have OpenBabel installed, any of the formats supported by
                OpenBabel. Non-case sensitive.

        Returns:
            (str) if filename is None. None otherwise.
        """
        from pymatgen.io.babel import BabelMolAdaptor
        from pymatgen.io.gaussian import GaussianInput
        from pymatgen.io.xyz import XYZ

        fmt = fmt.lower()
        writer: Any
        if fmt == "xyz" or fnmatch(filename.lower(), "*.xyz*"):
            writer = XYZ(self)
        elif any(fmt == ext or fnmatch(filename.lower(), f"*.{ext}*") for ext in ["gjf", "g03", "g09", "com", "inp"]):
            writer = GaussianInput(self)
        elif fmt == "json" or fnmatch(filename, "*.json*") or fnmatch(filename, "*.mson*"):
            if filename:
                with zopen(filename, "wt", encoding="utf8") as f:
                    json.dump(self.as_dict(), f)
                    return None
            else:
                return json.dumps(self.as_dict())
        elif fmt == "yaml" or fnmatch(filename, "*.yaml*"):
            yaml = YAML()
            if filename:
                with zopen(filename, "wt", encoding="utf8") as f:
                    return yaml.dump(self.as_dict(), f)
            else:
                sio = StringIO()
                yaml.dump(self.as_dict(), sio)
                return sio.getvalue()
        else:
            m = re.search(r"\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)", filename.lower())
            if (not fmt) and m:
                fmt = m.group(1)
            writer = BabelMolAdaptor(self)
            return writer.write_file(filename, file_format=fmt)

        if filename:
            writer.write_file(filename)
        return str(writer)

    @classmethod
    def from_str(
        cls, input_string: str, fmt: Literal["xyz", "gjf", "g03", "g09", "com", "inp", "json", "yaml"]
    ) -> IMolecule | Molecule:
        """
        Reads the molecule from a string.

        Args:
            input_string (str): String to parse.
            fmt (str): Format to output to. Defaults to JSON unless filename
                is provided. If fmt is specifies, it overrides whatever the
                filename is. Options include "xyz", "gjf", "g03", "json". If
                you have OpenBabel installed, any of the formats supported by
                OpenBabel. Non-case sensitive.

        Returns:
            IMolecule or Molecule.
        """
        from pymatgen.io.gaussian import GaussianInput
        from pymatgen.io.xyz import XYZ

        if fmt.lower() == "xyz":
            mol = XYZ.from_str(input_string).molecule
        elif fmt in ["gjf", "g03", "g09", "com", "inp"]:
            mol = GaussianInput.from_str(input_string).molecule
        elif fmt == "json":
            dct = json.loads(input_string)
            return cls.from_dict(dct)
        elif fmt == "yaml":
            yaml = YAML()
            dct = yaml.load(input_string)
            return cls.from_dict(dct)
        else:
            from pymatgen.io.babel import BabelMolAdaptor

            mol = BabelMolAdaptor.from_str(input_string, file_format=fmt).pymatgen_mol
        return cls.from_sites(mol)

    @classmethod
    def from_file(cls, filename):
        """
        Reads a molecule from a file. Supported formats include xyz,
        gaussian input (gjf|g03|g09|com|inp), Gaussian output (.out|and
        pymatgen's JSON-serialized molecules. Using openbabel,
        many more extensions are supported but requires openbabel to be
        installed.

        Args:
            filename (str): The filename to read from.

        Returns:
            Molecule
        """
        filename = str(filename)
        from pymatgen.io.gaussian import GaussianOutput

        with zopen(filename) as f:
            contents = f.read()
        fname = filename.lower()
        if fnmatch(fname, "*.xyz*"):
            return cls.from_str(contents, fmt="xyz")
        if any(fnmatch(fname.lower(), f"*.{r}*") for r in ["gjf", "g03", "g09", "com", "inp"]):
            return cls.from_str(contents, fmt="g09")
        if any(fnmatch(fname.lower(), f"*.{r}*") for r in ["out", "lis", "log"]):
            return GaussianOutput(filename).final_structure
        if fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
            return cls.from_str(contents, fmt="json")
        if fnmatch(fname, "*.yaml*"):
            return cls.from_str(contents, fmt="yaml")
        from pymatgen.io.babel import BabelMolAdaptor

        m = re.search(r"\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)", filename.lower())
        if m:
            new = BabelMolAdaptor.from_file(filename, m.group(1)).pymatgen_mol
            new.__class__ = cls
            return new
        raise ValueError("Cannot determine file type.")


class Structure(IStructure, collections.abc.MutableSequence):
    """Mutable version of structure."""

    __hash__ = None  # type: ignore

    def __init__(
        self,
        lattice: ArrayLike | Lattice,
        species: Sequence[CompositionLike],
        coords: Sequence[ArrayLike],
        charge: float | None = None,
        validate_proximity: bool = False,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        site_properties: dict | None = None,
        labels: list | None = None,
    ):
        """
        Create a periodic structure.

        Args:
            lattice: The lattice, either as a pymatgen.core.lattice.Lattice or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species: List of species on each site. Can take in flexible input,
                including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            charge (int): overall charge of the structure. Defaults to behavior
                in SiteCollection where total charge is the sum of the oxidation
                states.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            to_unit_cell (bool): Whether to map all sites into the unit cell,
                i.e., fractional coords between 0 and 1. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            labels (list[str]): Labels associated with the sites as a
                list of strings, e.g. ['Li1', 'Li2']. Must have the same
                length as the species and fractional coords. Defaults to
                None for no labels.
        """
        super().__init__(
            lattice,
            species,
            coords,
            charge=charge,
            validate_proximity=validate_proximity,
            to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
            labels=labels,
        )

        self._sites: list[PeriodicSite] = list(self._sites)  # type: ignore

    def __setitem__(  # type: ignore
        self, i: int | slice | Sequence[int] | SpeciesLike, site: SpeciesLike | PeriodicSite | Sequence
    ):
        """
        Modify a site in the structure.

        Args:
            i (int, [int], slice, Species-like): Indices to change. You can
                specify these as an int, a list of int, or a species-like
                string.
            site (PeriodicSite/Species/Sequence): Three options exist. You
                can provide a PeriodicSite directly (lattice will be
                checked). Or more conveniently, you can provide a
                specie-like object or a tuple of up to length 3.

        Examples:
            s[0] = "Fe"
            s[0] = Element("Fe")
            both replaces the species only.
            s[0] = "Fe", [0.5, 0.5, 0.5]
            Replaces site and *fractional* coordinates. Any properties
            are inherited from current site.
            s[0] = "Fe", [0.5, 0.5, 0.5], {"spin": 2}
            Replaces site and *fractional* coordinates and properties.

            s[(0, 2, 3)] = "Fe"
            Replaces sites 0, 2 and 3 with Fe.

            s[0::2] = "Fe"
            Replaces all even index sites with Fe.

            s["Mn"] = "Fe"
            Replaces all Mn in the structure with Fe. This is
            a short form for the more complex replace_species.

            s["Mn"] = "Fe0.5Co0.5"
            Replaces all Mn in the structure with Fe: 0.5, Co: 0.5, i.e.,
            creates a disordered structure!
        """
        if isinstance(i, int):
            indices = [i]
        elif isinstance(i, (str, Element, Species)):
            self.replace_species({i: site})  # type: ignore
            return
        elif isinstance(i, slice):
            to_mod = self[i]
            indices = [ii for ii, s in enumerate(self._sites) if s in to_mod]
        else:
            indices = list(i)

        for ii in indices:
            if isinstance(site, PeriodicSite):
                if site.lattice != self._lattice:
                    raise ValueError("PeriodicSite added must have same lattice as Structure!")
                if len(indices) != 1:
                    raise ValueError("Site assignments makes sense only for single int indices!")
                self._sites[ii] = site
            elif isinstance(site, str) or (not isinstance(site, collections.abc.Sequence)):
                self._sites[ii].species = site  # type: ignore
            else:
                self._sites[ii].species = site[0]  # type: ignore
                if len(site) > 1:
                    self._sites[ii].frac_coords = site[1]  # type: ignore
                if len(site) > 2:
                    self._sites[ii].properties = site[2]  # type: ignore

    def __delitem__(self, idx: SupportsIndex | slice) -> None:
        """Deletes a site from the Structure."""
        self._sites.__delitem__(idx)

    @property
    def lattice(self) -> Lattice:
        """:return: Lattice associated with structure."""
        return self._lattice

    @lattice.setter
    def lattice(self, lattice: ArrayLike | Lattice):
        if not isinstance(lattice, Lattice):
            lattice = Lattice(lattice)
        self._lattice = lattice
        for site in self._sites:
            site.lattice = lattice

    def append(  # type: ignore
        self,
        species: CompositionLike,
        coords: ArrayLike,
        coords_are_cartesian: bool = False,
        validate_proximity: bool = False,
        properties: dict | None = None,
    ):
        """
        Append a site to the structure.

        Args:
            species: Species of inserted site
            coords (3x1 array): Coordinates of inserted site
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.
            properties (dict): Properties of the site.

        Returns:
            New structure with inserted site.
        """
        return self.insert(
            len(self),
            species,
            coords,
            coords_are_cartesian=coords_are_cartesian,
            validate_proximity=validate_proximity,
            properties=properties,
        )

    def insert(  # type: ignore
        self,
        i: int,
        species: CompositionLike,
        coords: ArrayLike,
        coords_are_cartesian: bool = False,
        validate_proximity: bool = False,
        properties: dict | None = None,
    ):
        """
        Insert a site to the structure.

        Args:
            i (int): Index to insert site
            species (species-like): Species of inserted site
            coords (3x1 array): Coordinates of inserted site
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.
            properties (dict): Properties associated with the site.

        Returns:
            New structure with inserted site.
        """
        if not coords_are_cartesian:
            new_site = PeriodicSite(species, coords, self._lattice, properties=properties)
        else:
            frac_coords = self._lattice.get_fractional_coords(coords)
            new_site = PeriodicSite(species, frac_coords, self._lattice, properties=properties)

        if validate_proximity:
            for site in self:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing site!")

        self._sites.insert(i, new_site)

    def replace(
        self,
        i: int,
        species: CompositionLike,
        coords: ArrayLike | None = None,
        coords_are_cartesian: bool = False,
        properties: dict | None = None,
    ):
        """
        Replace a single site. Takes either a species or a dict of species and
        occupations.

        Args:
            i (int): Index of the site in the _sites list.
            species (species-like): Species of replacement site
            coords (3x1 array): Coordinates of replacement site. If None,
                the current coordinates are assumed.
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            properties (dict): Properties associated with the site.
        """
        if coords is None:
            frac_coords = self[i].frac_coords
        elif coords_are_cartesian:
            frac_coords = self._lattice.get_fractional_coords(coords)
        else:
            frac_coords = coords  # type: ignore

        new_site = PeriodicSite(species, frac_coords, self._lattice, properties=properties)
        self._sites[i] = new_site

    def substitute(self, index: int, func_group: IMolecule | Molecule | str, bond_order: int = 1) -> None:
        """
        Substitute atom at index with a functional group.

        Args:
            index (int): Index of atom to substitute.
            func_group: Substituent molecule. There are two options:

                1. Providing an actual Molecule as the input. The first atom
                   must be a DummySpecies X, indicating the position of
                   nearest neighbor. The second atom must be the next
                   nearest atom. For example, for a methyl group
                   substitution, func_group should be X-CH3, where X is the
                   first site and C is the second site. What the code will
                   do is to remove the index site, and connect the nearest
                   neighbor to the C atom in CH3. The X-C bond indicates the
                   directionality to connect the atoms.
                2. A string name. The molecule will be obtained from the
                   relevant template in func_groups.json.
            bond_order (int): A specified bond order to calculate the bond
                length between the attached functional group and the nearest
                neighbor site. Defaults to 1.
        """
        # Find the nearest neighbor that is not a terminal atom.
        all_non_terminal_nn = []
        for nn, dist, _, _ in self.get_neighbors(self[index], 3):
            # Check that the nn has neighbors within a sensible distance but
            # is not the site being substituted.
            for inn, dist2, _, _ in self.get_neighbors(nn, 3):
                if inn != self[index] and dist2 < 1.2 * get_bond_length(nn.specie, inn.specie):
                    all_non_terminal_nn.append((nn, dist))
                    break

        if len(all_non_terminal_nn) == 0:
            raise RuntimeError("Can't find a non-terminal neighbor to attach functional group to.")

        non_terminal_nn = min(all_non_terminal_nn, key=lambda d: d[1])[0]

        # Set the origin point to be the coordinates of the nearest
        # non-terminal neighbor.
        origin = non_terminal_nn.coords

        # Pass value of functional group--either from user-defined or from
        # functional.json
        if not isinstance(func_group, Molecule):
            # Check to see whether the functional group is in database.
            if func_group not in FunctionalGroups:
                raise RuntimeError("Can't find functional group in list. Provide explicit coordinate instead")
            fgroup = FunctionalGroups[func_group]
        else:
            fgroup = func_group

        # If a bond length can be found, modify func_grp so that the X-group
        # bond length is equal to the bond length.
        try:
            bl = get_bond_length(non_terminal_nn.specie, fgroup[1].specie, bond_order=bond_order)
        # Catches for case of incompatibility between Element(s) and Species(s)
        except TypeError:
            bl = None

        if bl is not None:
            fgroup = fgroup.copy()
            vec = fgroup[0].coords - fgroup[1].coords
            vec /= np.linalg.norm(vec)
            fgroup[0] = "X", fgroup[1].coords + float(bl) * vec

        # Align X to the origin.
        x = fgroup[0]
        fgroup.translate_sites(list(range(len(fgroup))), origin - x.coords)

        # Find angle between the attaching bond and the bond to be replaced.
        v1 = fgroup[1].coords - origin
        v2 = self[index].coords - origin
        angle = get_angle(v1, v2)

        if 1 < abs(angle % 180) < 179:
            # For angles which are not 0 or 180, we perform a rotation about
            # the origin along an axis perpendicular to both bonds to align
            # bonds.
            axis = np.cross(v1, v2)
            op = SymmOp.from_origin_axis_angle(origin, axis, angle)
            fgroup.apply_operation(op)
        elif abs(abs(angle) - 180) < 1:
            # We have a 180 degree angle. Simply do an inversion about the
            # origin
            for i, fg in enumerate(fgroup):
                fgroup[i] = (fg.species, origin - (fg.coords - origin))

        # Remove the atom to be replaced, and add the rest of the functional
        # group.
        del self[index]
        for site in fgroup[1:]:
            s_new = PeriodicSite(site.species, site.coords, self.lattice, coords_are_cartesian=True)
            self._sites.append(s_new)

    def remove_species(self, species: Sequence[SpeciesLike]) -> None:
        """
        Remove all occurrences of several species from a structure.

        Args:
            species: Sequence of species to remove, e.g., ["Li", "Na"].
        """
        new_sites = []
        species = [get_el_sp(s) for s in species]

        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species.items() if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(
                    PeriodicSite(
                        new_sp_occu,
                        site.frac_coords,
                        self._lattice,
                        properties=site.properties,
                    )
                )
        self._sites = new_sites

    def remove_sites(self, indices: Sequence[int]) -> None:
        """
        Delete sites with at indices.

        Args:
            indices: Sequence of indices of sites to delete.
        """
        self._sites = [s for i, s in enumerate(self._sites) if i not in indices]

    def apply_operation(self, symmop: SymmOp, fractional: bool = False) -> Structure:
        """
        Apply a symmetry operation to the structure in place and return the modified
        structure. The lattice is operated on by the rotation matrix only.
        Coords are operated in full and then transformed to the new lattice.

        Args:
            symmop (SymmOp): Symmetry operation to apply.
            fractional (bool): Whether the symmetry operation is applied in
                fractional space. Defaults to False, i.e., symmetry operation
                is applied in Cartesian coordinates.

        Returns:
            Structure: post-operation structure
        """
        if not fractional:
            self._lattice = Lattice([symmop.apply_rotation_only(row) for row in self._lattice.matrix])

            def operate_site(site):
                new_cart = symmop.operate(site.coords)
                new_frac = self._lattice.get_fractional_coords(new_cart)
                return PeriodicSite(
                    site.species,
                    new_frac,
                    self._lattice,
                    properties=site.properties,
                    skip_checks=True,
                )

        else:
            new_latt = np.dot(symmop.rotation_matrix, self._lattice.matrix)
            self._lattice = Lattice(new_latt)

            def operate_site(site):
                return PeriodicSite(
                    site.species,
                    symmop.operate(site.frac_coords),
                    self._lattice,
                    properties=site.properties,
                    skip_checks=True,
                )

        self._sites = [operate_site(s) for s in self._sites]

        return self

    def apply_strain(self, strain: ArrayLike) -> None:
        """
        Apply a strain to the lattice.

        Args:
            strain (float or list): Amount of strain to apply. Can be a float,
                or a sequence of 3 numbers. E.g., 0.01 means all lattice
                vectors are increased by 1%. This is equivalent to calling
                modify_lattice with a lattice with lattice parameters that
                are 1% larger.
        """
        strain_matrix = (1 + np.array(strain)) * np.eye(3)
        self.lattice = Lattice(np.dot(self._lattice.matrix.T, strain_matrix).T)

    def sort(self, key: Callable | None = None, reverse: bool = False) -> None:
        """
        Sort a structure in place. The parameters have the same meaning as in
        list.sort. By default, sites are sorted by the electronegativity of
        the species. The difference between this method and
        get_sorted_structure (which also works in IStructure) is that the
        latter returns a new Structure, while this just sorts the Structure
        in place.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        self._sites.sort(key=key, reverse=reverse)

    def translate_sites(
        self, indices: int | Sequence[int], vector: ArrayLike, frac_coords: bool = True, to_unit_cell: bool = True
    ) -> None:
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            indices: Integer or List of site indices on which to perform the
                translation.
            vector: Translation vector for sites.
            frac_coords (bool): Whether the vector corresponds to fractional or
                Cartesian coordinates.
            to_unit_cell (bool): Whether new sites are transformed to unit
                cell
        """
        if not isinstance(indices, collections.abc.Iterable):
            indices = [indices]

        for i in indices:
            site = self._sites[i]
            if frac_coords:
                fcoords = site.frac_coords + vector
            else:
                fcoords = self._lattice.get_fractional_coords(site.coords + vector)
            if to_unit_cell:
                fcoords = [np.mod(f, 1) if p else f for p, f in zip(self.lattice.pbc, fcoords)]
            self._sites[i].frac_coords = fcoords

    def rotate_sites(
        self,
        indices: list[int] | None = None,
        theta: float = 0.0,
        axis: ArrayLike | None = None,
        anchor: ArrayLike | None = None,
        to_unit_cell: bool = True,
    ) -> None:
        """
        Rotate specific sites by some angle around vector at anchor.

        Args:
            indices (list): List of site indices on which to perform the
                translation.
            theta (float): Angle in radians
            axis (3x1 array): Rotation axis vector.
            anchor (3x1 array): Point of rotation.
            to_unit_cell (bool): Whether new sites are transformed to unit
                cell
        """
        from numpy import cross, eye
        from numpy.linalg import norm
        from scipy.linalg import expm

        if indices is None:
            indices = list(range(len(self)))

        if axis is None:
            axis = [0, 0, 1]

        if anchor is None:
            anchor = [0, 0, 0]

        anchor = np.array(anchor)
        axis = np.array(axis)

        theta %= 2 * np.pi

        rm = expm(cross(eye(3), axis / norm(axis)) * theta)
        for idx in indices:
            site = self._sites[idx]
            coords = ((np.dot(rm, np.array(site.coords - anchor).T)).T + anchor).ravel()
            new_site = PeriodicSite(
                site.species,
                coords,
                self._lattice,
                to_unit_cell=to_unit_cell,
                coords_are_cartesian=True,
                properties=site.properties,
                skip_checks=True,
            )
            self._sites[idx] = new_site

    def perturb(self, distance: float, min_distance: float | None = None) -> None:
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance (float): Distance in angstroms by which to perturb each
                site.
            min_distance (None, int, or float): if None, all displacements will
                be equal amplitude. If int or float, perturb each site a
                distance drawn from the uniform distribution between
                'min_distance' and 'distance'.
        """

        def get_rand_vec():
            # deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            dist = distance
            if isinstance(min_distance, (float, int)):
                dist = np.random.uniform(min_distance, dist)
            return vector / vnorm * dist if vnorm != 0 else get_rand_vec()

        for idx in range(len(self._sites)):
            self.translate_sites([idx], get_rand_vec(), frac_coords=False)

    def make_supercell(self, scaling_matrix: ArrayLike, to_unit_cell: bool = True) -> None:
        """
        Create a supercell.

        Args:
            scaling_matrix: A scaling matrix for transforming the lattice
                vectors. Has to be all integers. Several options are possible:

                a. A full 3x3 scaling matrix defining the linear combination
                   the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
                   1]] generates a new structure with lattice vectors a' =
                   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
                   vectors of the original structure.
                b. An sequence of three scaling factors. E.g., [2, 1, 1]
                   specifies that the supercell should have dimensions 2a x b x
                   c.
                c. A number, which simply scales all lattice vectors by the
                   same factor.
            to_unit_cell: Whether or not to fall back sites into the unit cell
        """
        supercell = self * scaling_matrix
        if to_unit_cell:
            for site in supercell:
                site.to_unit_cell(in_place=True)
        self._sites = supercell.sites
        self._lattice = supercell.lattice

    def scale_lattice(self, volume: float) -> None:
        """
        Performs a scaling of the lattice vectors so that length proportions
        and angles are preserved.

        Args:
            volume (float): New volume of the unit cell in A^3.
        """
        self.lattice = self._lattice.scale(volume)

    def merge_sites(self, tol: float = 0.01, mode: Literal["sum", "delete", "average"] = "sum") -> None:
        """
        Merges sites (adding occupancies) within tol of each other.
        Removes site properties.

        Args:
            tol (float): Tolerance for distance to merge sites.
            mode ('sum' | 'delete' | 'average'): "delete" means duplicate sites are
                deleted. "sum" means the occupancies are summed for the sites.
                "average" means that the site is deleted but the properties are averaged
                Only first letter is considered.
        """
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform

        dist_mat = self.distance_matrix
        np.fill_diagonal(dist_mat, 0)
        clusters = fcluster(linkage(squareform((dist_mat + dist_mat.T) / 2)), tol, "distance")
        sites = []
        for c in np.unique(clusters):
            inds = np.where(clusters == c)[0]
            species = self[inds[0]].species
            coords = self[inds[0]].frac_coords
            props = self[inds[0]].properties
            for n, i in enumerate(inds[1:]):
                sp = self[i].species
                if mode.lower()[0] == "s":
                    species += sp
                offset = self[i].frac_coords - coords
                coords = coords + ((offset - np.round(offset)) / (n + 2)).astype(coords.dtype)
                for key in props:
                    if props[key] is not None and self[i].properties[key] != props[key]:
                        if mode.lower()[0] == "a" and isinstance(props[key], float):
                            # update a running total
                            props[key] = props[key] * (n + 1) / (n + 2) + self[i].properties[key] / (n + 2)
                        else:
                            props[key] = None
                            warnings.warn(
                                f"Sites with different site property {key} are merged. So property is set to none"
                            )
            sites.append(PeriodicSite(species, coords, self.lattice, properties=props))

        self._sites = sites

    def set_charge(self, new_charge: float = 0.0) -> None:
        """
        Sets the overall structure charge.

        Args:
            new_charge (float): new charge to set
        """
        self._charge = new_charge

    def relax(
        self,
        calculator: str | Calculator = "m3gnet",
        relax_cell: bool = True,
        optimizer: str | Optimizer = "FIRE",
        steps: int = 500,
        fmax: float = 0.1,
        stress_weight: float = 0.01,
        opt_kwargs: dict | None = None,
        return_trajectory: bool = False,
        verbose: bool = False,
    ) -> Structure | tuple[Structure, TrajectoryObserver | Trajectory]:
        """
        Performs a crystal structure relaxation using an ASE calculator.

        Args:
            calculator: An ASE Calculator or a string from the following options: "m3gnet".
                Defaults to 'm3gnet', i.e. the M3GNet universal potential.
            relax_cell (bool): whether to relax the lattice cell. Defaults to True.
            optimizer (str): name of the ASE optimizer class to use
            steps (int): max number of steps for relaxation. Defaults to 500.
            fmax (float): total force tolerance for relaxation convergence.
                Here fmax is a sum of force and stress forces. Defaults to 0.1.
            stress_weight (float): the stress weight for relaxation with M3GNet.
                Defaults to 0.01.
            opt_kwargs (dict): kwargs for the ASE optimizer class.
            return_trajectory (bool): Whether to return the trajectory of relaxation.
                Defaults to False.
            verbose (bool): whether to print out relaxation steps. Defaults to False.

        Returns:
            Structure | tuple[Structure, Trajectory]: Relaxed structure or if return_trajectory=True,
                2-tuple of Structure and matgl TrajectoryObserver.
        """
        return self._relax(
            calculator,
            relax_cell=relax_cell,
            optimizer=optimizer,
            steps=steps,
            fmax=fmax,
            stress_weight=stress_weight,
            opt_kwargs=opt_kwargs,
            return_trajectory=return_trajectory,
            verbose=verbose,
        )

    def calculate(self, calculator: str | Calculator = "m3gnet", verbose: bool = False) -> Calculator:
        """
        Performs an ASE calculation.

        Args:
            calculator: An ASE Calculator or a string from the following options: "m3gnet".
                Defaults to 'm3gnet', i.e. the M3GNet universal potential.
            verbose (bool): whether to print stdout. Defaults to False.

        Returns:
            Calculator: ASE Calculator instance with a results attribute containing the output.
        """
        return self._calculate(calculator, verbose=verbose)

    @classmethod
    def from_prototype(cls, prototype: str, species: Sequence, **kwargs) -> Structure:
        """
        Method to rapidly construct common prototype structures.

        Args:
            prototype: Name of prototype. E.g., cubic, rocksalt, perovksite etc.
            species: List of species corresponding to symmetrically distinct sites.
            **kwargs: Lattice parameters, e.g., a = 3.0, b = 4, c = 5. Only the required lattice parameters need to be
                specified. For example, if it is a cubic prototype, only a needs to be specified.

        Returns:
            Structure
        """
        prototype = prototype.lower()
        try:
            if prototype == "fcc":
                return Structure.from_spacegroup("Fm-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0]])
            if prototype == "bcc":
                return Structure.from_spacegroup("Im-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0]])
            if prototype == "hcp":
                return Structure.from_spacegroup(
                    "P6_3/mmc", Lattice.hexagonal(kwargs["a"], kwargs["c"]), species, [[1 / 3, 2 / 3, 1 / 4]]
                )
            if prototype == "diamond":
                return Structure.from_spacegroup("Fd-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0]])
            if prototype == "rocksalt":
                return Structure.from_spacegroup(
                    "Fm-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0], [0.5, 0.5, 0]]
                )
            if prototype == "perovskite":
                return Structure.from_spacegroup(
                    "Pm-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0]]
                )
            if prototype in ("cscl"):
                return Structure.from_spacegroup(
                    "Pm-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0], [0.5, 0.5, 0.5]]
                )
            if prototype in ("fluorite", "caf2"):
                return Structure.from_spacegroup(
                    "Fm-3m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0], [1 / 4, 1 / 4, 1 / 4]]
                )
            if prototype in ("antifluorite"):
                return Structure.from_spacegroup(
                    "Fm-3m", Lattice.cubic(kwargs["a"]), species, [[1 / 4, 1 / 4, 1 / 4], [0, 0, 0]]
                )
            if prototype in ("zincblende"):
                return Structure.from_spacegroup(
                    "F-43m", Lattice.cubic(kwargs["a"]), species, [[0, 0, 0], [1 / 4, 1 / 4, 3 / 4]]
                )
        except KeyError as exc:
            raise ValueError(f"Required parameter {exc} not specified as a kwargs!") from exc
        raise ValueError(f"Unsupported {prototype=}!")


class Molecule(IMolecule, collections.abc.MutableSequence):
    """
    Mutable Molecule. It has all the methods in IMolecule, but in addition,
    it allows a user to perform edits on the molecule.
    """

    __hash__ = None  # type: ignore

    def __init__(
        self,
        species: Sequence[SpeciesLike],
        coords: Sequence[ArrayLike],
        charge: float = 0.0,
        spin_multiplicity: int | None = None,
        validate_proximity: bool = False,
        site_properties: dict | None = None,
        labels: Sequence[str | None] | None = None,
        charge_spin_check: bool = True,
    ) -> None:
        """
        Creates a MutableMolecule.

        Args:
            species: list of atomic species. Possible kinds of input include a
                list of dict of elements/species and occupancies, a List of
                elements/specie specified as actual Element/Species, Strings
                ("Fe", "Fe2+") or atomic numbers (1,56).
            coords (3x1 array): list of Cartesian coordinates of each species.
            charge (float): Charge for the molecule. Defaults to 0.
            spin_multiplicity (int): Spin multiplicity for molecule.
                Defaults to None, which means that the spin multiplicity is
                set to 1 if the molecule has no unpaired electrons and to 2
                if there are unpaired electrons.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 1 Ang apart. Defaults to False.
            site_properties (dict): Properties associated with the sites as
                a dict of sequences, e.g., {"magmom":[5,5,5,5]}. The
                sequences have to be the same length as the atomic species
                and fractional_coords. Defaults to None for no properties.
            labels (list[str]): Labels associated with the sites as a
                list of strings, e.g. ['Li1', 'Li2']. Must have the same
                length as the species and fractional coords. Defaults to
                None for no labels.
            charge_spin_check (bool): Whether to check that the charge and
                spin multiplicity are compatible with each other. Defaults
                to True.
        """
        super().__init__(
            species,
            coords,
            charge=charge,
            spin_multiplicity=spin_multiplicity,
            validate_proximity=validate_proximity,
            site_properties=site_properties,
            labels=labels,
            charge_spin_check=charge_spin_check,
        )
        self._sites: list[Site] = list(self._sites)  # type: ignore

    def __setitem__(  # type: ignore
        self, idx: int | slice | Sequence[int] | SpeciesLike, site: SpeciesLike | Site | Sequence
    ) -> None:
        """
        Modify a site in the molecule.

        Args:
            idx (int, [int], slice, Species-like): Indices to change. You can
                specify these as an int, a list of int, or a species-like
                string.
            site (PeriodicSite/Species/Sequence): Three options exist. You can
                provide a Site directly, or for convenience, you can provide
                simply a Species-like string/object, or finally a (Species,
                coords) sequence, e.g., ("Fe", [0.5, 0.5, 0.5]).
        """
        if isinstance(idx, int):
            indices = [idx]
        elif isinstance(idx, (str, Element, Species)):
            self.replace_species({idx: site})  # type: ignore
            return
        elif isinstance(idx, slice):
            to_mod = self[idx]
            indices = [idx for idx, site in enumerate(self._sites) if site in to_mod]
        else:
            indices = list(idx)

        for ii in indices:
            if isinstance(site, Site):
                self._sites[ii] = site
            elif isinstance(site, str) or not isinstance(site, collections.abc.Sequence):
                self._sites[ii].species = site  # type: ignore
            else:
                self._sites[ii].species = site[0]  # type: ignore
                if len(site) > 1:
                    self._sites[ii].coords = site[1]  # type: ignore
                if len(site) > 2:
                    self._sites[ii].properties = site[2]  # type: ignore

    def __delitem__(self, idx: SupportsIndex | slice) -> None:
        """Deletes a site from the Structure."""
        self._sites.__delitem__(idx)

    def append(  # type: ignore
        self,
        species: CompositionLike,
        coords: ArrayLike,
        validate_proximity: bool = False,
        properties: dict | None = None,
    ) -> Molecule:
        """
        Appends a site to the molecule.

        Args:
            species: Species of inserted site
            coords: Coordinates of inserted site
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.
            properties (dict): A dict of properties for the Site.

        Returns:
            New molecule with inserted site.
        """
        return self.insert(
            len(self),
            species,
            coords,
            validate_proximity=validate_proximity,
            properties=properties,
        )

    def set_charge_and_spin(self, charge: float, spin_multiplicity: int | None = None):
        """
        Set the charge and spin multiplicity.

        Args:
            charge (int): Charge for the molecule. Defaults to 0.
            spin_multiplicity (int): Spin multiplicity for molecule.
                Defaults to None, which means that the spin multiplicity is
                set to 1 if the molecule has no unpaired electrons and to 2
                if there are unpaired electrons.
        """
        self._charge = charge
        n_electrons = 0.0
        for site in self._sites:
            for sp, amt in site.species.items():
                if not isinstance(sp, DummySpecies):
                    n_electrons += sp.Z * amt
        n_electrons -= charge
        self._nelectrons = n_electrons
        if spin_multiplicity:
            if self._charge_spin_check and (n_electrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    f"Charge of {self._charge} and spin multiplicity of {spin_multiplicity} is"
                    " not possible for this molecule"
                )
            self._spin_multiplicity = spin_multiplicity
        else:
            self._spin_multiplicity = 1 if n_electrons % 2 == 0 else 2

    def insert(  # type: ignore
        self,
        i: int,
        species: CompositionLike,
        coords: ArrayLike,
        validate_proximity: bool = False,
        properties: dict | None = None,
    ) -> Molecule:
        """
        Insert a site to the molecule.

        Args:
            i (int): Index to insert site
            species: species of inserted site
            coords (3x1 array): coordinates of inserted site
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to True.
            properties (dict): Dict of properties for the Site.

        Returns:
            New molecule with inserted site.
        """
        new_site = Site(species, coords, properties=properties)
        if validate_proximity:
            for site in self:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing site!")
        self._sites.insert(i, new_site)

        return self

    def remove_species(self, species: Sequence[SpeciesLike]):
        """
        Remove all occurrences of a species from a molecule.

        Args:
            species: Species to remove.
        """
        new_sites = []
        species = [get_el_sp(sp) for sp in species]
        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species.items() if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(Site(new_sp_occu, site.coords, properties=site.properties))
        self._sites = new_sites

    def remove_sites(self, indices: Sequence[int]):
        """
        Delete sites with at indices.

        Args:
            indices: Sequence of indices of sites to delete.
        """
        self._sites = [self._sites[i] for i in range(len(self._sites)) if i not in indices]

    def translate_sites(self, indices: Sequence[int] | None = None, vector: ArrayLike | None = None):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            indices (list): List of site indices on which to perform the
                translation.
            vector (3x1 array): Translation vector for sites.
        """
        if indices is None:
            indices = range(len(self))
        if vector is None:
            vector = [0, 0, 0]
        for i in indices:
            site = self._sites[i]
            new_site = Site(site.species, site.coords + vector, properties=site.properties)
            self._sites[i] = new_site

    def rotate_sites(
        self,
        indices: Sequence[int] | None = None,
        theta: float = 0.0,
        axis: ArrayLike | None = None,
        anchor: ArrayLike | None = None,
    ):
        """
        Rotate specific sites by some angle around vector at anchor.

        Args:
            indices (list): List of site indices on which to perform the
                translation.
            theta (float): Angle in radians
            axis (3x1 array): Rotation axis vector.
            anchor (3x1 array): Point of rotation.
        """
        from numpy import cross, eye
        from numpy.linalg import norm
        from scipy.linalg import expm

        if indices is None:
            indices = range(len(self))

        if axis is None:
            axis = [0, 0, 1]

        if anchor is None:
            anchor = [0, 0, 0]

        anchor = np.array(anchor)
        axis = np.array(axis)

        theta %= 2 * np.pi

        rm = expm(cross(eye(3), axis / norm(axis)) * theta)

        for idx in indices:
            site = self._sites[idx]
            s = ((np.dot(rm, (site.coords - anchor).T)).T + anchor).ravel()
            new_site = Site(site.species, s, properties=site.properties)
            self._sites[idx] = new_site

    def perturb(self, distance: float):
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance (float): Distance in angstroms by which to perturb each
                site.
        """

        def get_rand_vec():
            # deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            return vector / vnorm * distance if vnorm != 0 else get_rand_vec()

        for idx in range(len(self._sites)):
            self.translate_sites([idx], get_rand_vec())

    def apply_operation(self, symmop: SymmOp):
        """
        Apply a symmetry operation to the molecule.

        Args:
            symmop (SymmOp): Symmetry operation to apply.
        """

        def operate_site(site):
            new_cart = symmop.operate(site.coords)
            return Site(site.species, new_cart, properties=site.properties)

        self._sites = [operate_site(site) for site in self._sites]

    def copy(self):
        """
        Convenience method to get a copy of the molecule.

        Returns:
            A copy of the Molecule.
        """
        return type(self).from_sites(self)

    def substitute(self, index: int, func_group: IMolecule | Molecule | str, bond_order: int = 1):
        """
        Substitute atom at index with a functional group.

        Args:
            index (int): Index of atom to substitute.
            func_group: Substituent molecule. There are two options:

                1. Providing an actual molecule as the input. The first atom
                   must be a DummySpecies X, indicating the position of
                   nearest neighbor. The second atom must be the next
                   nearest atom. For example, for a methyl group
                   substitution, func_group should be X-CH3, where X is the
                   first site and C is the second site. What the code will
                   do is to remove the index site, and connect the nearest
                   neighbor to the C atom in CH3. The X-C bond indicates the
                   directionality to connect the atoms.
                2. A string name. The molecule will be obtained from the
                   relevant template in func_groups.json.
            bond_order (int): A specified bond order to calculate the bond
                length between the attached functional group and the nearest
                neighbor site. Defaults to 1.
        """
        # Find the nearest neighbor that is not a terminal atom.
        all_non_terminal_nn = []
        for nn in self.get_neighbors(self[index], 3):
            # Check that the nn has neighbors within a sensible distance but
            # is not the site being substituted.
            for nn2 in self.get_neighbors(nn, 3):
                if nn2 != self[index] and nn2.nn_distance < 1.2 * get_bond_length(nn.specie, nn2.specie):
                    all_non_terminal_nn.append(nn)
                    break

        if len(all_non_terminal_nn) == 0:
            raise RuntimeError("Can't find a non-terminal neighbor to attach functional group to.")

        non_terminal_nn = min(all_non_terminal_nn, key=lambda nn: nn.nn_distance)

        # Set the origin point to be the coordinates of the nearest
        # non-terminal neighbor.
        origin = non_terminal_nn.coords

        # Pass value of functional group--either from user-defined or from
        # functional.json
        if isinstance(func_group, Molecule):
            functional_group = func_group
        else:
            # Check to see whether the functional group is in database.
            if func_group not in FunctionalGroups:
                raise RuntimeError("Can't find functional group in list. Provide explicit coordinate instead")
            functional_group = FunctionalGroups[func_group]

        # If a bond length can be found, modify func_grp so that the X-group
        # bond length is equal to the bond length.
        bond_len = get_bond_length(non_terminal_nn.specie, functional_group[1].specie, bond_order=bond_order)
        if bond_len is not None:
            functional_group = functional_group.copy()
            vec = functional_group[0].coords - functional_group[1].coords
            vec /= np.linalg.norm(vec)
            functional_group[0] = "X", functional_group[1].coords + float(bond_len) * vec

        # Align X to the origin.
        x = functional_group[0]
        functional_group.translate_sites(list(range(len(functional_group))), origin - x.coords)

        # Find angle between the attaching bond and the bond to be replaced.
        v1 = functional_group[1].coords - origin
        v2 = self[index].coords - origin
        angle = get_angle(v1, v2)

        if 1 < abs(angle % 180) < 179:
            # For angles which are not 0 or 180, we perform a rotation about
            # the origin along an axis perpendicular to both bonds to align
            # bonds.
            axis = np.cross(v1, v2)
            op = SymmOp.from_origin_axis_angle(origin, axis, angle)
            functional_group.apply_operation(op)
        elif abs(abs(angle) - 180) < 1:
            # We have a 180 degree angle. Simply do an inversion about the
            # origin
            for i, fg in enumerate(functional_group):
                functional_group[i] = (fg.species, origin - (fg.coords - origin))

        # Remove the atom to be replaced, and add the rest of the functional
        # group.
        del self[index]
        self._sites += list(functional_group[1:])

    def relax(
        self,
        calculator: str | Calculator = "gfn2-xtb",
        optimizer: str | Optimizer = "FIRE",
        steps: int = 500,
        fmax: float = 0.1,
        opt_kwargs: dict | None = None,
        return_trajectory: bool = False,
        verbose: bool = False,
    ) -> Molecule | tuple[Molecule, TrajectoryObserver]:
        """
        Performs a molecule relaxation using an ASE calculator.

        Args:
            calculator: An ASE Calculator or a string from the following options: "gfn2-xtb".
                Defaults to 'gfn2-xtb'.
            optimizer (str): name of the ASE optimizer class to use
            steps (int): max number of steps for relaxation. Defaults to 500.
            fmax (float): total force tolerance for relaxation convergence.
                Defaults to 0.1 eV/A.
            opt_kwargs (dict): kwargs for the ASE optimizer class.
            return_trajectory (bool): Whether to return the trajectory of relaxation.
                Defaults to False.
            verbose (bool): whether to print out relaxation steps. Defaults to False.

        Returns:
            Molecule | tuple[Molecule, Trajectory]: Relaxed Molecule or if return_trajectory=True,
                2-tuple of Molecule and ASE TrajectoryObserver.
        """
        return self._relax(
            calculator,
            relax_cell=False,
            optimizer=optimizer,
            steps=steps,
            fmax=fmax,
            opt_kwargs=opt_kwargs,
            return_trajectory=return_trajectory,
            verbose=verbose,
        )

    def calculate(self, calculator: str | Calculator = "gfn2-xtb", verbose: bool = False) -> Calculator:
        """
        Performs an ASE calculation.

        Args:
            calculator: An ASE Calculator or a string from the following options: "gfn2-xtb".
                Defaults to 'gfn2-xtb'.
            verbose (bool): whether to print stdout. Defaults to False.

        Returns:
            Calculator: ASE Calculator instance with a results attribute containing the output.
        """
        return self._calculate(calculator, verbose=verbose)


class StructureError(Exception):
    """
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    """


with open(os.path.join(os.path.dirname(__file__), "func_groups.json")) as file:
    FunctionalGroups = {k: Molecule(v["species"], v["coords"]) for k, v in json.load(file).items()}
