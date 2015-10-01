# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides classes used to define a non-periodic molecule and a
periodic structure.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"


import math
import os
import json
import collections
import itertools
from abc import ABCMeta, abstractmethod, abstractproperty
import random
import warnings
from fnmatch import fnmatch
import re
from fractions import gcd

import six
from six.moves import map, zip

import numpy as np

import yaml
try:
    from yaml import CSafeDumper as Dumper, CLoader as Loader
except ImportError:
    from yaml import SafeDumper as Dumper, Loader

from pymatgen.core.operations import SymmOp
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, Specie, get_el_sp
from pymatgen.serializers.json_coders import PMGSONable
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.core.bonds import CovalentBond, get_bond_length
from pymatgen.core.composition import Composition
from pymatgen.util.coord_utils import get_angle, all_distances, \
    lattice_points_in_supercell
from monty.design_patterns import singleton
from pymatgen.core.units import Mass, Length, ArrayWithUnit
from pymatgen.symmetry.groups import SpaceGroup
from monty.io import zopen


class SiteCollection(six.with_metaclass(ABCMeta, collections.Sequence)):
    """
    Basic SiteCollection. Essentially a sequence of Sites or PeriodicSites.
    This serves as a base class for Molecule (a collection of Site, i.e., no
    periodicity) and Structure (a collection of PeriodicSites, i.e.,
    periodicity). Not meant to be instantiated directly.
    """

    #Tolerance in Angstrom for determining if sites are too close.
    DISTANCE_TOLERANCE = 0.01

    @abstractproperty
    def sites(self):
        """
        Returns a tuple of sites.
        """
        return

    @abstractmethod
    def get_distance(self, i, j):
        """
        Returns distance between sites at index i and j.

        Args:
            i (int): Index of first site
            j (int): Index of second site

        Returns:
            (float) Distance between sites at index i and index j.
        """
        return

    @property
    def distance_matrix(self):
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this is overwritten to return the nearest image
        distance.
        """
        return all_distances(self.cart_coords, self.cart_coords)

    @property
    def species(self):
        """
        Only works for ordered structures.
        Disordered structures will raise an AttributeError.

        Returns:
            ([Specie]) List of species at each site of the structure.
        """
        return [site.specie for site in self]

    @property
    def species_and_occu(self):
        """
        List of species and occupancies at each site of the structure.
        """
        return [site.species_and_occu for site in self]

    @property
    def ntypesp(self):
        """Number of types of atoms."""
        return len(self.types_of_specie)

    @property
    def types_of_specie(self):
        """
        List of types of specie. Only works for ordered structures.
        Disordered structures will raise an AttributeError.
        """
        # Cannot use set since we want a deterministic algorithm.
        types = []
        for site in self:
            if site.specie not in types:
                types.append(site.specie)
        return types

    def group_by_types(self):
        """Iterate over species grouped by type"""
        for t in self.types_of_specie:
            for site in self:
                if site.specie == t:
                    yield site

    def indices_from_symbol(self, symbol):
        """
        Returns a tuple with the sequential indices of the sites
        that contain an element with the given chemical symbol.
        """
        indices = []
        for i, specie in enumerate(self.species):
            if specie.symbol == symbol:
                indices.append(i)
        return tuple(indices)

    @property
    def symbol_set(self):
        """
        Tuple with the set of chemical symbols.
        Note that len(symbol_set) == len(types_of_specie)
        """
        return tuple([specie.symbol for specie in self.types_of_specie])

    @property
    def atomic_numbers(self):
        """List of atomic numbers."""
        return [site.specie.number for site in self]

    @property
    def site_properties(self):
        """
        Returns the site properties as a dict of sequences. E.g.,
        {"magmom": (5,-5), "charge": (-4,4)}.
        """
        props = {}
        prop_keys = set()
        for site in self:
            prop_keys.update(site.properties.keys())

        for k in prop_keys:
            props[k] = [site.properties.get(k, None) for site in self]
        return props

    def __contains__(self, site):
        return site in self.sites

    def __iter__(self):
        return self.sites.__iter__()

    def __getitem__(self, ind):
        return self.sites[ind]

    def __len__(self):
        return len(self.sites)

    def __hash__(self):
        # for now, just use the composition hash code.
        return self.composition.__hash__()

    @property
    def num_sites(self):
        """
        Number of sites.
        """
        return len(self)

    @property
    def cart_coords(self):
        """
        Returns a list of the cartesian coordinates of sites in the structure.
        """
        return np.array([site.coords for site in self])

    @property
    def formula(self):
        """
        (str) Returns the formula.
        """
        return self.composition.formula

    @property
    def composition(self):
        """
        (Composition) Returns the composition
        """
        elmap = collections.defaultdict(float)
        for site in self:
            for species, occu in site.species_and_occu.items():
                elmap[species] += occu
        return Composition(elmap)

    @property
    def charge(self):
        """
        Returns the net charge of the structure based on oxidation states. If
        Elements are found, a charge of 0 is assumed.
        """
        charge = 0
        for site in self:
            for specie, amt in site.species_and_occu.items():
                charge += getattr(specie, "oxi_state", 0) * amt
        return charge

    @property
    def is_ordered(self):
        """
        Checks if structure is ordered, meaning no partial occupancies in any
        of the sites.
        """
        return all((site.is_ordered for site in self))

    def get_angle(self, i, j, k):
        """
        Returns angle specified by three sites.

        Args:
            i (int): Index of first site.
            j (int): Index of second site.
            k (int): Index of third site.

        Returns:
            (float) Angle in degrees.
        """
        v1 = self[i].coords - self[j].coords
        v2 = self[k].coords - self[j].coords
        return get_angle(v1, v2, units="degrees")

    def get_dihedral(self, i, j, k, l):
        """
        Returns dihedral angle specified by four sites.

        Args:
            i (int): Index of first site
            j (int): Index of second site
            k (int): Index of third site
            l (int): Index of fourth site

        Returns:
            (float) Dihedral angle in degrees.
        """
        v1 = self[k].coords - self[l].coords
        v2 = self[j].coords - self[k].coords
        v3 = self[i].coords - self[j].coords
        v23 = np.cross(v2, v3)
        v12 = np.cross(v1, v2)
        return math.degrees(math.atan2(np.linalg.norm(v2) * np.dot(v1, v23),
                            np.dot(v12, v23)))

    def is_valid(self, tol=DISTANCE_TOLERANCE):
        """
        True if SiteCollection does not contain atoms that are too close
        together. Note that the distance definition is based on type of
        SiteCollection. Cartesian distances are used for non-periodic
        Molecules, while PBC is taken into account for periodic structures.

        Args:
            tol (float): Distance tolerance. Default is 0.01A.

        Returns:
            (bool) True if SiteCollection does not contain atoms that are too
            close together.
        """
        if len(self.sites) == 1:
            return True
        all_dists = self.distance_matrix[np.triu_indices(len(self), 1)]
        return bool(np.min(all_dists) > tol)

    @abstractmethod
    def to(self, fmt=None, filename=None):
        """
        Generates well-known string representations of SiteCollections (e.g.,
        molecules / structures). Should return a string type or write to a file.
        """
        pass

    @classmethod
    @abstractmethod
    def from_str(cls, input_string, fmt):
        """
        Reads in SiteCollection from a string.
        """
        pass

    @classmethod
    @abstractmethod
    def from_file(cls, filename):
        """
        Reads in SiteCollection from a filename.
        """
        pass


class IStructure(SiteCollection, PMGSONable):
    """
    Basic immutable Structure object with periodicity. Essentially a sequence
    of PeriodicSites having a common lattice. IStructure is made to be
    (somewhat) immutable so that they can function as keys in a dict. To make
    modifications, use the standard Structure object instead. Structure
    extends Sequence and Hashable, which means that in many cases,
    it can be used like any Python sequence. Iterating through a
    structure is equivalent to going through the sites in sequence.
    """

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):
        """
        Create a periodic structure.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Specie]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
        """
        if len(species) != len(coords):
            raise StructureError("The list of atomic species must be of the"
                                 " same length as the list of fractional"
                                 " coordinates.")

        if isinstance(lattice, Lattice):
            self._lattice = lattice
        else:
            self._lattice = Lattice(lattice)

        sites = []
        for i in range(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[i]
                        for k, v in site_properties.items()}

            sites.append(
                PeriodicSite(species[i], coords[i], self._lattice,
                             to_unit_cell,
                             coords_are_cartesian=coords_are_cartesian,
                             properties=prop))
        self._sites = tuple(sites)
        if validate_proximity and not self.is_valid():
            raise StructureError(("Structure contains sites that are ",
                                  "less than 0.01 Angstrom apart!"))

    @classmethod
    def from_sites(cls, sites, validate_proximity=False,
                   to_unit_cell=False):
        """
        Convenience constructor to make a Structure from a list of sites.

        Args:
            sites: Sequence of PeriodicSites. Sites must have the same
                lattice.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            to_unit_cell (bool): Whether to translate sites into the unit
                cell.

        Returns:
            (Structure) Note that missing properties are set as None.
        """
        prop_keys = []
        props = {}
        lattice = None
        for i, site in enumerate(sites):
            if not lattice:
                lattice = site.lattice
            elif site.lattice != lattice:
                raise ValueError("Sites must belong to the same lattice")
            for k, v in site.properties.items():
                if k not in prop_keys:
                    prop_keys.append(k)
                    props[k] = [None] * len(sites)
                props[k][i] = v
        for k, v in props.items():
            if any((vv is None for vv in v)):
                warnings.warn("Not all sites have property %s. Missing values "
                              "are set to None." % k)
        return cls(lattice, [site.species_and_occu for site in sites],
                   [site.frac_coords for site in sites],
                   site_properties=props,
                   validate_proximity=validate_proximity,
                   to_unit_cell=to_unit_cell)

    @classmethod
    def from_spacegroup(cls, sg, lattice, species, coords, site_properties=None,
                        coords_are_cartesian=False, tol=1e-5):
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
            species ([Specie]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            tol (float): A fractional tolerance to deal with numerical
               precision issues in determining if orbits are the same.
        """
        try:
            i = int(sg)
            sgp = SpaceGroup.from_int_number(i)
        except ValueError:
            sgp = SpaceGroup(sg)

        if isinstance(lattice, Lattice):
            latt = lattice
        else:
            latt = Lattice(lattice)

        if not sgp.is_compatible(latt):
            raise ValueError(
                "Supplied lattice with parameters %s is incompatible with "
                "supplied spacegroup %s!" % (latt.lengths_and_angles,
                                             sgp.symbol)
            )

        frac_coords = coords if not coords_are_cartesian else \
            lattice.get_fractional_coords(coords)

        props = {} if site_properties is None else site_properties

        all_sp = []
        all_coords = []
        all_site_properties = collections.defaultdict(list)
        for i, (sp, c) in enumerate(zip(species, frac_coords)):
            cc = sgp.get_orbit(c, tol=tol)
            all_sp.extend([sp] * len(cc))
            all_coords.extend(cc)
            for k, v in props.items():
                all_site_properties[k].extend([v[i]] * len(cc))

        return cls(latt, all_sp, all_coords,
                   site_properties=all_site_properties)

    @classmethod
    def from_xsf_string(cls, input_string):
        """
        Initialize a `Structure` object from a string with data in XSF format.
        See http://www.xcrysden.org/doc/XSF.html
        """
        # CRYSTAL                                        see (1)
        # these are primitive lattice vectors (in Angstroms) 
        # PRIMVEC
        #    0.0000000    2.7100000    2.7100000         see (2)
        #    2.7100000    0.0000000    2.7100000
        #    2.7100000    2.7100000    0.0000000

        # these are conventional lattice vectors (in Angstroms)
        # CONVVEC
        #    5.4200000    0.0000000    0.0000000         see (3)
        #    0.0000000    5.4200000    0.0000000
        #    0.0000000    0.0000000    5.4200000

        # these are atomic coordinates in a primitive unit cell  (in Angstroms)
        # PRIMCOORD
        # 2 1                                            see (4)
        # 16      0.0000000     0.0000000     0.0000000  see (5)
        # 30      1.3550000    -1.3550000    -1.3550000

        lattice, coords, species = [], [], []
        lines = input_string.splitlines()

        for i in range(len(lines)):
            if "PRIMVEC" in lines[i]:
                for j in range(i+1, i+4):
                    lattice.append([float(c) for c in lines[j].split()])

            if "PRIMCOORD" in lines[i]:
                num_sites = int(lines[i+1].split()[0])
     
                for j in range(i+2, i+2+num_sites):
                    tokens = lines[j].split()
                    species.append(int(tokens[0]))
                    coords.append([float(j) for j in tokens[1:]])
                break
        else:
            raise ValueError("Invalid XSF data")

        return cls(lattice, species, coords, coords_are_cartesian=True)

    @classmethod
    def from_abivars(cls, *args, **kwargs):
        """Build a :class:`Structure` object from a dictionary with ABINIT variables."""
        kwargs.update(dict(*args))
        d = kwargs

        lattice = Lattice.from_abivars(d)
        coords, coords_are_cartesian = d.get("xred", None), False

        if coords is None:
            coords = d.get("xcart", None)
            if coords is not None:
                coords = ArrayWithUnit(coords, "bohr").to("ang")
            else:
                coords = d.get("xangst", None)
            coords_are_cartesian = True

        if coords is None:
            raise ValueError("Cannot extract atomic coordinates from dict %s" % str(d))

        coords = np.reshape(coords, (-1,3))

        znucl_type, typat = d["znucl"], d["typat"]

        if not isinstance(znucl_type, collections.Iterable):
            znucl_type = [znucl_type]

        if not isinstance(typat, collections.Iterable):
            typat = [typat]

        assert len(typat) == len(coords)

        # Note Fortran --> C indexing
        #znucl_type = np.rint(znucl_type)
        species = [znucl_type[typ-1] for typ in typat]

        #print(lattice, species, coords)
        #print(species)

        return cls(lattice, species, coords, validate_proximity=False,
                   to_unit_cell=False, coords_are_cartesian=coords_are_cartesian)

    @property
    def distance_matrix(self):
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this should return the nearest image distance.
        """
        return self.lattice.get_all_distances(self.frac_coords,
                                              self.frac_coords)

    @property
    def sites(self):
        """
        Returns an iterator for the sites in the Structure.
        """
        return self._sites

    @property
    def lattice(self):
        """
        Lattice of the structure.
        """
        return self._lattice

    @property
    def reciprocal_lattice(self):
        """
        Reciprocal lattice of the structure.
        """
        return self._lattice.reciprocal_lattice

    def lattice_vectors(self, space="r"):
        """
        Returns the vectors of the unit cell in Angstrom.

        Args:
            space: "r" for real space vectors, "g" for reciprocal space basis
                vectors.
        """
        if space.lower() == "r":
            return self.lattice.matrix
        if space.lower() == "g":
            return self.lattice.reciprocal_lattice.matrix
        raise ValueError("Wrong value for space: %s " % space)

    @property
    def density(self):
        """
        Returns the density in units of g/cc
        """
        m = Mass(self.composition.weight, "amu")
        return m.to("g") / (self.volume * Length(1, "ang").to("cm") ** 3)

    def __eq__(self, other):
        if other is None:
            return False
        if len(self) != len(other):
            return False
        if self.lattice != other.lattice:
            return False
        for site in self:
            if site not in other:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # For now, just use the composition hash code.
        return self.composition.__hash__()

    @property
    def frac_coords(self):
        """
        Fractional coordinates as a Nx3 numpy array.
        """
        return np.array([site.frac_coords for site in self._sites])

    @property
    def volume(self):
        """
        Returns the volume of the structure.
        """
        return self._lattice.volume

    def get_distance(self, i, j, jimage=None):
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

    def get_sites_in_sphere(self, pt, r, include_index=False):
        """
        Find all sites within a sphere from the point. This includes sites
        in other periodic images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelpiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            pt (3x1 array): cartesian coordinates of center of sphere.
            r (float): Radius of sphere.
            include_index (bool): Whether the non-supercell site index
                is included in the returned data

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        site_fcoords = np.mod(self.frac_coords, 1)
        neighbors = []
        for fcoord, dist, i in self._lattice.get_points_in_sphere(
                site_fcoords, pt, r):
            nnsite = PeriodicSite(self[i].species_and_occu,
                                  fcoord, self._lattice,
                                  properties=self[i].properties)
            neighbors.append((nnsite, dist) if not include_index
                             else (nnsite, dist, i))
        return neighbors

    def get_neighbors(self, site, r, include_index=False):
        """
        Get all neighbors to a site within a sphere of radius r.  Excludes the
        site itself.

        Args:
            site:
                site, which is the center of the sphere.
            r:
                radius of sphere.
            include_index:
                boolean that determines whether the non-supercell site index
                is included in the returned data

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        nn = self.get_sites_in_sphere(site.coords, r,
                                      include_index=include_index)
        return [d for d in nn if site != d[0]]

    def get_all_neighbors(self, r, include_index=False):
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

        Args:
            r (float): Radius of sphere.
            include_index (bool): Whether to include the non-supercell site
                in the returned data

        Returns:
            A list of a list of nearest neighbors for each site, i.e.,
            [[(site, dist, index) ...], ..]
            Index only supplied if include_index = True.
            The index is the index of the site in the original (non-supercell)
            structure. This is needed for ewaldmatrix by keeping track of which
            sites contribute to the ewald sum.
        """
        # Use same algorithm as get_sites_in_sphere to determine supercell but
        # loop over all atoms in crystal
        recp_len = np.array(self.lattice.reciprocal_lattice.abc)
        maxr = np.ceil((r + 0.15) * recp_len / (2 * math.pi))
        nmin = np.floor(np.min(self.frac_coords, axis=0)) - maxr
        nmax = np.ceil(np.max(self.frac_coords, axis=0)) + maxr

        all_ranges = [np.arange(x, y) for x, y in zip(nmin, nmax)]

        latt = self._lattice
        neighbors = [list() for i in range(len(self._sites))]
        all_fcoords = np.mod(self.frac_coords, 1)
        coords_in_cell = latt.get_cartesian_coords(all_fcoords)
        site_coords = self.cart_coords

        indices = np.arange(len(self))
        for image in itertools.product(*all_ranges):
            coords = latt.get_cartesian_coords(image) + coords_in_cell
            all_dists = all_distances(coords, site_coords)
            all_within_r = np.bitwise_and(all_dists <= r, all_dists > 1e-8)

            for (j, d, within_r) in zip(indices, all_dists, all_within_r):
                nnsite = PeriodicSite(self[j].species_and_occu, coords[j],
                                      latt, properties=self[j].properties,
                                      coords_are_cartesian=True)
                for i in indices[within_r]:
                    item = (nnsite, d[i], j) if include_index else (
                        nnsite, d[i])
                    neighbors[i].append(item)
        return neighbors

    def get_neighbors_in_shell(self, origin, r, dr, include_index=False):
        """
        Returns all sites in a shell centered on origin (coords) between radii
        r-dr and r+dr.

        Args:
            origin (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Inner radius of shell.
            dr (float): Width of shell.
            include_index (bool): Whether to include the non-supercell site
                in the returned data

        Returns:
            [(site, dist, index) ...] since most of the time, subsequent
            processing
            requires the distance. Index only supplied if include_index = True.
            The index is the index of the site in the original (non-supercell)
            structure. This is needed for ewaldmatrix by keeping track of which
            sites contribute to the ewald sum.
        """
        outer = self.get_sites_in_sphere(origin, r + dr,
                                         include_index=include_index)
        inner = r - dr
        return [t for t in outer if t[1] > inner]

    def get_sorted_structure(self, key=None, reverse=False):
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
        return self.__class__.from_sites(sites)

    def get_reduced_structure(self, reduction_algo="niggli"):
        """
        Get a reduced structure.

        Args:
            reduction_algo (str): The lattice reduction algorithm to use.
                Currently supported options are "niggli" or "LLL".
        """
        if reduction_algo == "niggli":
            reduced_latt = self._lattice.get_niggli_reduced_lattice()
        elif reduction_algo == "LLL":
            reduced_latt = self._lattice.get_lll_reduced_lattice()
        else:
            raise ValueError("Invalid reduction algo : {}"
                             .format(reduction_algo))

        return self.__class__(reduced_latt, self.species_and_occu,
                              self.cart_coords,
                              coords_are_cartesian=True, to_unit_cell=True)

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
            return self.__class__(self._lattice,
                                  self.species_and_occu,
                                  self.frac_coords,
                                  site_properties=props)
        else:
            reduced_latt = self._lattice.get_lll_reduced_lattice()
            new_sites = []
            for i, site in enumerate(self):
                frac_coords = reduced_latt.get_fractional_coords(site.coords)
                site_props = {}
                for p in props:
                    site_props[p] = props[p][i]
                new_sites.append(PeriodicSite(site.species_and_occu,
                                              frac_coords, reduced_latt,
                                              to_unit_cell=True,
                                              properties=site_props))
            new_sites = sorted(new_sites)
            return self.__class__.from_sites(new_sites)

    def interpolate(self, end_structure, nimages=10,
                    interpolate_lattices=False, pbc=True, autosort_tol=0):
        """
        Interpolate between this structure and end_structure. Useful for
        construction of NEB inputs.

        Args:
            end_structure (Structure): structure to interpolate between this
                structure and end.
            nimages (int): No. of interpolation images. Defaults to 10 images.
            interpolate_lattices (bool): Whether to interpolate the lattices.
                Interpolates the lengths and angles (rather than the matrix)
                so orientation may be affected.
            pbc (bool): Whether to use periodic boundary conditions to find
                the shortest path between endpoints.
            auto_sort_tol (float): A distance tolerance in angstrom in
                which to automatically sort end_structure to match to the
                closest points in this particular structure. This is usually
                what you want in a NEB calculation. 0 implies no sorting.
                Otherwise, a 0.5 value usually works pretty well.

        Returns:
            List of interpolated structures. The starting and ending
            structures included as the first and last structures respectively.
            A total of (nimages + 1) structures are returned.
        """
        #Check length of structures
        if len(self) != len(end_structure):
            raise ValueError("Structures have different lengths!")

        if interpolate_lattices:
            #interpolate lattices
            lstart = np.array(self.lattice.lengths_and_angles)
            lend = np.array(end_structure.lattice.lengths_and_angles)
            lvec = lend - lstart

        #Check that both structures have the same lattice
        elif not self.lattice == end_structure.lattice:
            raise ValueError("Structures with different lattices!")

        #Check that both structures have the same species
        for i in range(len(self)):
            if self[i].species_and_occu != end_structure[i].species_and_occu:
                raise ValueError("Different species!\nStructure 1:\n" +
                                 str(self) + "\nStructure 2\n" +
                                 str(end_structure))

        start_coords = np.array(self.frac_coords)
        end_coords = np.array(end_structure.frac_coords)

        if autosort_tol:
            dist_matrix = self.lattice.get_all_distances(start_coords, end_coords)
            site_mappings = collections.defaultdict(list)
            unmapped_start_ind = []
            for i, row in enumerate(dist_matrix):
                ind = np.where(row < autosort_tol)[0]
                if len(ind) == 1:
                    site_mappings[i].append(ind[0])
                else:
                    unmapped_start_ind.append(i)

            if len(unmapped_start_ind) > 1:
                raise ValueError("Unable to reliably match structures "
                                 "with auto_sort_tol = %f. unmapped indices "
                                 "= %s" % (autosort_tol, unmapped_start_ind))

            sorted_end_coords = np.zeros_like(end_coords)
            matched = []
            for i, j in site_mappings.items():
                if len(j) > 1:
                    raise ValueError("Unable to reliably match structures "
                                     "with auto_sort_tol = %f. More than one "
                                     "site match!" % autosort_tol)
                sorted_end_coords[i] = end_coords[j[0]]
                matched.append(j[0])

            if len(unmapped_start_ind) == 1:
                i = unmapped_start_ind[0]
                j = list(set(range(len(start_coords))).difference(matched))[0]
                sorted_end_coords[i] = end_coords[j]

            end_coords = sorted_end_coords

        vec = end_coords - start_coords
        if pbc:
            vec -= np.round(vec)
        sp = self.species_and_occu
        structs = []
        for x in range(nimages + 1):
            if interpolate_lattices:
                l_a = lstart + x / nimages * lvec
                l = Lattice.from_lengths_and_angles(*l_a)
            else:
                l = self.lattice
            fcoords = start_coords + x / nimages * vec
            structs.append(self.__class__(l, sp, fcoords,
                           site_properties=self.site_properties))
        return structs

    def get_primitive_structure(self, tolerance=0.25):
        """
        This finds a smaller unit cell than the input. Sometimes it doesn"t
        find the smallest possible one, so this method is recursively called
        until it is unable to find a smaller cell.

        NOTE: if the tolerance is greater than 1/2 the minimum inter-site
        distance in the primitive cell, the algorithm will reject this lattice.

        Args:
            tolerance (float): Tolerance for each coordinate of a particular
                site. For example, [0.1, 0, 0.1] in cartesian coordinates
                will be considered to be on the same coordinates as
                [0, 0, 0] for a tolerance of 0.25. Defaults to 0.25.

        Returns:
            The most primitive structure found.
        """
        # group sites by species string
        k = lambda s: s.species_string
        sites = sorted(self._sites, key=k)
        grouped_sites = [list(a[1]) for a in itertools.groupby(sites, key=k)]
        grouped_fcoords = [np.array([s.frac_coords for s in g])
                           for g in grouped_sites]

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
            within tolerance to some coordinate in fc2
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
                for i in range(1, n+1):
                    if n % i == 0:
                        yield i

            for det in factors(fu):
                if det == 1:
                    continue
                for a in factors(det):
                    for e in factors(det // a):
                        g = det // a // e
                        yield det, np.array(
                            [[[a, b, c], [0, e, f], [0, 0, g]]
                            for b, c, f in itertools.product(range(a), range(a),
                                                             range(e))])

        # we cant let sites match to their neighbors in the supercell
        grouped_non_nbrs = []
        for gfcoords in grouped_fcoords:
            fdist = gfcoords[None, :, :] - gfcoords[:, None, :]
            fdist -= np.round(fdist)
            np.abs(fdist, fdist)
            non_nbrs = np.any(fdist > 2 * super_ftol[None, None, :], axis=-1)
            # since we want sites to match to themselves
            np.fill_diagonal(non_nbrs, True)
            grouped_non_nbrs.append(non_nbrs)

        num_fu = six.moves.reduce(gcd, map(len, grouped_sites))
        for size, ms in get_hnf(num_fu):
            inv_ms = np.linalg.inv(ms)

            # find sets of lattice vectors that are are present in min_vecs
            dist = inv_ms[:, :, None, :] - min_vecs[None, None, :, :]
            dist -= np.round(dist)
            np.abs(dist, dist)
            is_close = np.all(dist < super_ftol, axis=-1)
            any_close = np.any(is_close, axis=-1)
            inds = np.all(any_close, axis=-1)

            for inv_m, m in zip(inv_ms[inds], ms[inds]):
                new_m = np.dot(inv_m, self.lattice.matrix)
                ftol = np.divide(tolerance, np.sqrt(np.sum(new_m ** 2, axis=1)))

                valid = True
                new_coords = []
                new_sp = []
                for gsites, gfcoords, non_nbrs in zip(grouped_sites,
                                                      grouped_fcoords,
                                                      grouped_non_nbrs):
                    all_frac = np.dot(gfcoords, m)

                    # calculate grouping of equivalent sites, represented by
                    # adjacency matrix
                    fdist = all_frac[None, :, :] - all_frac[:, None, :]
                    fdist = np.abs(fdist - np.round(fdist))
                    close_in_prim = np.all(fdist < ftol[None, None, :], axis=-1)
                    groups = np.logical_and(close_in_prim, non_nbrs)

                    #check that groups are correct
                    if not np.all(np.sum(groups, axis=0) == size):
                        valid = False
                        break

                    #check that groups are all cliques
                    for g in groups:
                        if not np.all(groups[g][:, g]):
                            valid = False
                            break
                    if not valid:
                        break

                    #add the new sites, averaging positions
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
                            new_sp.append(gsites[inds[0]].species_and_occu)
                            new_coords.append(coords)

                if valid:
                    inv_m = np.linalg.inv(m)
                    new_l = Lattice(np.dot(inv_m, self.lattice.matrix))
                    s = Structure(new_l, new_sp, new_coords,
                                  coords_are_cartesian=False)

                    return s.get_primitive_structure(
                        tolerance).get_reduced_structure()

        return Structure.from_sites(self)

    def __repr__(self):
        outs = ["Structure Summary", repr(self.lattice)]
        for s in self:
            outs.append(repr(s))
        return "\n".join(outs)

    def __str__(self):
        outs = ["Full Formula ({s})".format(s=self.composition.formula),
                "Reduced Formula: {}"
                .format(self.composition.reduced_formula)]
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i), site.species_string,
                                  " ".join([to_s(j).rjust(12)
                                            for j in site.frac_coords])]))
        return "\n".join(outs)

    def as_dict(self):
        """
        Json-serializable dict representation of Structure
        """
        latt_dict = self._lattice.as_dict()
        del latt_dict["@module"]
        del latt_dict["@class"]

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice": latt_dict, "sites": []}
        for site in self:
            site_dict = site.as_dict()
            del site_dict["lattice"]
            del site_dict["@module"]
            del site_dict["@class"]
            d["sites"].append(site_dict)
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a Structure object from a dict representation of Structure
        created using as_dict().

        Args:
            d (dict): Dict representation of structure.

        Returns:
            Structure object
        """
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        return cls.from_sites(sites)

    def to_abivars(self, **kwargs):
        """Returns a dictionary with the ABINIT variables."""
        types_of_specie = self.types_of_specie
        natom = self.num_sites

        znucl_type = [specie.number for specie in types_of_specie]

        znucl_atoms = self.atomic_numbers

        typat = np.zeros(natom, np.int)
        for (atm_idx, site) in enumerate(self):
            typat[atm_idx] = types_of_specie.index(site.specie) + 1

        rprim = ArrayWithUnit(self.lattice.matrix, "ang").to("bohr")
        xred = np.reshape([site.frac_coords for site in self], (-1,3))

        # Set small values to zero. This usually happens when the CIF file
        # does not give structure parameters with enough digits.
        rprim = np.where(np.abs(rprim) > 1e-8, rprim, 0.0)
        xred = np.where(np.abs(xred) > 1e-8, xred, 0.0)

        # Info on atoms.
        d = dict(
            natom=natom,
            ntypat=len(types_of_specie),
            typat=typat,
            znucl=znucl_type,
            xred=xred,
        )

        # Add info on the lattice.
        # Should we use (rprim, acell) or (angdeg, acell) to specify the lattice?
        geomode = kwargs.pop("geomode", "rprim")
        #latt_dict = self.lattice.to_abivars(geomode=geomode)

        if geomode == "rprim":
            d.update(dict(
                acell=3 * [1.0],
                rprim=rprim))

        elif geomode == "angdeg":
            d.update(dict(
                acell=3 * [1.0],
                angdeg=angdeg))
        else:
            raise ValueError("Wrong value for geomode: %s" % geomode)

        return d

    def to_xsf_string(self):
        """
        Returns a string with the structure in XSF format
        See http://www.xcrysden.org/doc/XSF.html
        """
        lines = []
        app = lines.append

        app("CRYSTAL")
        app("# Primitive lattice vectors in Angstrom")
        app("PRIMVEC")
        cell = self.lattice_vectors(space="r") 
        for i in range(3):
            app(' %.14f %.14f %.14f' % tuple(cell[i]))

        cart_coords = self.cart_coords  
        app("# Cartesian coordinates in Angstrom.")
        app("PRIMCOORD")
        app(" %d 1" % len(cart_coords))

        for a in range(len(cart_coords)):
            sp = "%d" % self.atomic_numbers[a]
            app(sp + ' %20.14f %20.14f %20.14f' % tuple(cart_coords[a]))

        return "\n".join(lines)

    def dot(self, coords_a, coords_b, space="r", frac_coords=False):
        """
        Compute the scalar product of vector(s) either in real space or
        reciprocal space.

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            space (str): "r" for real space, "g" for reciprocal space.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        lattice = {"r": self.lattice,
                   "g": self.reciprocal_lattice}[space.lower()]
        return lattice.dot(coords_a, coords_b, frac_coords=frac_coords)

    def norm(self, coords, space="r", frac_coords=True):
        """
        Compute the norm of vector(s) either in real space or reciprocal space.

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            space (str): "r" for real space, "g" for reciprocal space.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        return np.sqrt(self.dot(coords, coords, space=space,
                                frac_coords=frac_coords))

    def to(self, fmt=None, filename=None):
        """
        Outputs the structure to a file or string.

        Args:
            fmt (str): Format to output to. Defaults to JSON unless filename
                is provided. If fmt is specifies, it overrides whatever the
                filename is. Options include "cif", "poscar", "cssr", "json".
                Non-case sensitive.
            filename (str): If provided, output will be written to a file. If
                fmt is not specified, the format is determined from the
                filename. Defaults is None, i.e. string output.

        Returns:
            (str) if filename is None. None otherwise.
        """
        from pymatgen.io.cif import CifWriter
        from pymatgen.io.vasp import Poscar
        from pymatgen.io.cssr import Cssr
        filename = filename or ""
        fmt = "" if fmt is None else fmt.lower()
        fname = os.path.basename(filename)

        if fmt == "cif" or fnmatch(fname, "*.cif*"):
            writer = CifWriter(self)
        elif fmt == "poscar" or fnmatch(fname, "POSCAR*"):
            writer = Poscar(self)
        elif fmt == "cssr" or fnmatch(fname.lower(), "*.cssr*"):
            writer = Cssr(self)
        elif fmt == "json" or fnmatch(fname.lower(), "*.json"):
            s = json.dumps(self.as_dict())
            if filename:
                with zopen(filename, "wt") as f:
                    # This complicated for handles unicode in both Py2 and 3.
                    f.write("%s" % s)
                return
            else:
                return s
        elif fmt == "xsf" or fnmatch(fname.lower(), "*.xsf*"):
            if filename:
                with zopen(fname, "wt", encoding='utf8') as f:
                    s = self.to_xsf_string()
                    f.write(s)
                    return s
            else:
                return self.to_xsf_string()
        else:
            if filename:
                with open(filename, "w") as f:
                    yaml.dump(self.as_dict(), f, Dumper=Dumper)
                return
            else:
                return yaml.dump(self.as_dict(), Dumper=Dumper)

        if filename:
            writer.write_file(filename)
        else:
            return writer.__str__()

    @classmethod
    def from_str(cls, input_string, fmt, primitive=False, sort=False):
        """
        Reads a structure from a string.

        Args:
            input_string (str): String to parse.
            fmt (str): A format specification.

        Returns:
            IStructure / Structure
        """
        from pymatgen.io.cif import CifParser
        from pymatgen.io.vasp import Poscar
        from pymatgen.io.cssr import Cssr
        fmt = fmt.lower()
        if fmt == "cif":
            parser = CifParser.from_string(input_string)
            s = parser.get_structures(primitive=primitive)[0]
        elif fmt == "poscar":
            s = Poscar.from_string(input_string, False).structure
        elif fmt == "cssr":
            cssr = Cssr.from_string(input_string)
            s = cssr.structure
        elif fmt == "json":
            d = json.loads(input_string)
            s = cls.from_dict(d)
        elif fmt == "yaml":
            d = yaml.load(input_string)
            s = cls.from_dict(d)
        elif fmt == "xsf":
            s = cls.from_xsf_string(input_string)
        else:
            raise ValueError("Unrecognized format `%s`!" % fmt)

        if sort:
            s = s.get_sorted_structure()
        return cls.from_sites(s)

    @classmethod
    def from_file(cls, filename, primitive=True, sort=False):
        """
        Reads a structure from a file. For example, anything ending in
        a "cif" is assumed to be a Crystallographic Information Format file.
        Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT,
        vasprun.xml, CSSR, Netcdf and pymatgen's JSON serialized structures.

        Args:
            filename (str): The filename to read from.
            primitive (bool): Whether to convert to a primitive cell
                Only available for cifs. Defaults to True.
            sort (bool): Whether to sort sites. Default to False.

        Returns:
            Structure.
        """
        if filename.endswith(".nc"):
            # Read Structure from a netcdf file.
            from pymatgen.io.abinit.netcdf import structure_from_ncdata
            s = structure_from_ncdata(filename, cls=cls)
            if sort:
                s = s.get_sorted_structure()
            return s

        from pymatgen.io.vasp import Vasprun, Chgcar
        from monty.io import zopen
        fname = os.path.basename(filename)
        with zopen(filename) as f:
            contents = f.read()
        if fnmatch(fname.lower(), "*.cif*"):
            return cls.from_str(contents, fmt="cif",
                                primitive=primitive, sort=sort)
        elif fnmatch(fname, "POSCAR*") or fnmatch(fname, "CONTCAR*"):
            return cls.from_str(contents, fmt="poscar",
                                primitive=primitive, sort=sort)

        elif fnmatch(fname, "CHGCAR*") or fnmatch(fname, "LOCPOT*"):
            s = Chgcar.from_file(filename).structure
        elif fnmatch(fname, "vasprun*.xml*"):
            s = Vasprun(filename).final_structure
        elif fnmatch(fname.lower(), "*.cssr*"):
            return cls.from_str(contents, fmt="cssr",
                                primitive=primitive, sort=sort)
        elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
            return cls.from_str(contents, fmt="json",
                                primitive=primitive, sort=sort)
        elif fnmatch(fname, "*.yaml*"):
            return cls.from_str(contents, fmt="yaml",
                                primitive=primitive, sort=sort)
        elif fnmatch(fname, "*.xsf"):
            return cls.from_str(contents, fmt="xsf",
                                primitive=primitive, sort=sort)
        else:
            raise ValueError("Unrecognized file extension!")
        if sort:
            s = s.get_sorted_structure()

        s.__class__ = cls
        return s


class IMolecule(SiteCollection, PMGSONable):
    """
    Basic immutable Molecule object without periodicity. Essentially a
    sequence of sites. IMolecule is made to be immutable so that they can
    function as keys in a dict. For a mutable molecule,
    use the :class:Molecule.

    Molecule extends Sequence and Hashable, which means that in many cases,
    it can be used like any Python sequence. Iterating through a molecule is
    equivalent to going through the sites in sequence.
    """

    def __init__(self, species, coords, charge=0,
                 spin_multiplicity=None, validate_proximity=False,
                 site_properties=None):
        """
        Creates a Molecule.

        Args:
            species: list of atomic species. Possible kinds of input include a
                list of dict of elements/species and occupancies, a List of
                elements/specie specified as actual Element/Specie, Strings
                ("Fe", "Fe2+") or atomic numbers (1,56).
            coords (3x1 array): list of cartesian coordinates of each species.
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
        """
        if len(species) != len(coords):
            raise StructureError(("The list of atomic species must be of the",
                                  " same length as the list of fractional ",
                                  "coordinates."))

        sites = []
        for i in range(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[i] for k, v in site_properties.items()}
            sites.append(Site(species[i], coords[i], properties=prop))

        self._sites = tuple(sites)
        if validate_proximity and not self.is_valid():
            raise StructureError(("Molecule contains sites that are ",
                                  "less than 0.01 Angstrom apart!"))

        self._charge = charge
        nelectrons = 0
        for site in sites:
            for sp, amt in site.species_and_occu.items():
                nelectrons += sp.Z * amt
        nelectrons -= charge
        self._nelectrons = nelectrons
        if spin_multiplicity:
            if (nelectrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    "Charge of {} and spin multiplicity of {} is"
                    " not possible for this molecule".format(
                    self._charge, spin_multiplicity))
            self._spin_multiplicity = spin_multiplicity
        else:
            self._spin_multiplicity = 1 if nelectrons % 2 == 0 else 2

    @property
    def charge(self):
        """
        Charge of molecule
        """
        return self._charge

    @property
    def spin_multiplicity(self):
        """
        Spin multiplicity of molecule.
        """
        return self._spin_multiplicity

    @property
    def nelectrons(self):
        """
        Number of electrons in the molecule.
        """
        return self._nelectrons

    @property
    def center_of_mass(self):
        """
        Center of mass of molecule.
        """
        center = np.zeros(3)
        total_weight = 0
        for site in self:
            wt = site.species_and_occu.weight
            center += site.coords * wt
            total_weight += wt
        return center / total_weight

    @property
    def sites(self):
        """
        Returns a tuple of sites in the Molecule.
        """
        return self._sites

    @classmethod
    def from_sites(cls, sites, charge=0, spin_multiplicity=None,
                   validate_proximity=False):
        """
        Convenience constructor to make a Molecule from a list of sites.

        Args:
            sites: Sequence of Sites.
        """
        props = collections.defaultdict(list)
        for site in sites:
            for k, v in site.properties.items():
                props[k].append(v)
        return cls([site.species_and_occu for site in sites],
                   [site.coords for site in sites],
                   charge=charge, spin_multiplicity=spin_multiplicity,
                   validate_proximity=validate_proximity,
                   site_properties=props)

    def break_bond(self, ind1, ind2, tol=0.2):
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
        sites = self._sites
        clusters = [[sites[ind1]], [sites[ind2]]]

        sites = [site for i, site in enumerate(sites) if i not in (ind1, ind2)]

        def belongs_to_cluster(site, cluster):
            for test_site in cluster:
                if CovalentBond.is_bonded(site, test_site, tol=tol):
                    return True
            return False

        while len(sites) > 0:
            unmatched = []
            for site in sites:
                found = False
                for cluster in clusters:
                    if belongs_to_cluster(site, cluster):
                        cluster.append(site)
                        found = True
                        break
                if not found:
                    unmatched.append(site)

            if len(unmatched) == len(sites):
                raise ValueError("Not all sites are matched!")
            sites = unmatched

        return (self.__class__.from_sites(cluster)
                for cluster in clusters)

    def get_covalent_bonds(self, tol=0.2):
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

    def __eq__(self, other):
        if other is None:
            return False
        if len(self) != len(other):
            return False
        if self._charge != other._charge:
            return False
        if self._spin_multiplicity != other._spin_multiplicity:
            return False
        for site in self:
            if site not in other:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # For now, just use the composition hash code.
        return self.composition.__hash__()

    def __repr__(self):
        outs = ["Molecule Summary"]
        for s in self:
            outs.append(s.__repr__())
        return "\n".join(outs)

    def __str__(self):
        outs = ["Full Formula ({s})".format(s=self.composition.formula),
                "Reduced Formula: " + self.composition.reduced_formula,
                "Charge = {}, Spin Mult = {}".format(
                    self._charge, self._spin_multiplicity)]
        to_s = lambda x: "%0.6f" % x
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i), site.species_string,
                                  " ".join([to_s(j).rjust(12) for j in
                                            site.coords])]))
        return "\n".join(outs)

    def as_dict(self):
        """
        Json-serializable dict representation of Molecule
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "charge": self._charge,
             "spin_multiplicity": self._spin_multiplicity,
             "sites": []}
        for site in self:
            site_dict = site.as_dict()
            del site_dict["@module"]
            del site_dict["@class"]
            d["sites"].append(site_dict)
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a Molecule object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of Molecule.

        Returns:
            Molecule object
        """
        species = []
        coords = []
        props = collections.defaultdict(list)

        for site_dict in d["sites"]:
            species.append({Specie(sp["element"], sp["oxidation_state"])
                            if "oxidation_state" in sp else
                            Element(sp["element"]): sp["occu"]
                            for sp in site_dict["species"]})
            coords.append(site_dict["xyz"])
            siteprops = site_dict.get("properties", {})
            for k, v in siteprops.items():
                props[k].append(v)

        return cls(species, coords, charge=d.get("charge", 0),
                   spin_multiplicity=d.get("spin_multiplicity"),
                   site_properties=props)

    def get_distance(self, i, j):
        """
        Get distance between site i and j.

        Args:
            i (int): Index of first site
            j (int): Index of second site

        Returns:
            Distance between the two sites.
        """
        return self[i].distance(self[j])

    def get_sites_in_sphere(self, pt, r):
        """
        Find all sites within a sphere from a point.

        Args:
            pt (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Radius of sphere.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        neighbors = []
        for site in self._sites:
            dist = site.distance_from_point(pt)
            if dist <= r:
                neighbors.append((site, dist))
        return neighbors

    def get_neighbors(self, site, r):
        """
        Get all neighbors to a site within a sphere of radius r.  Excludes the
        site itself.

        Args:
            site (Site): Site at the center of the sphere.
            r (float): Radius of sphere.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        nn = self.get_sites_in_sphere(site.coords, r)
        return [(s, dist) for (s, dist) in nn if site != s]

    def get_neighbors_in_shell(self, origin, r, dr):
        """
        Returns all sites in a shell centered on origin (coords) between radii
        r-dr and r+dr.

        Args:
            origin (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Inner radius of shell.
            dr (float): Width of shell.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        outer = self.get_sites_in_sphere(origin, r + dr)
        inner = r - dr
        return [(site, dist) for (site, dist) in outer if dist > inner]

    def get_boxed_structure(self, a, b, c, images=(1, 1, 1),
                            random_rotation=False, min_dist=1, cls=None):
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
            cls: The Structure class to instantiate (defaults to pymatgen structure)

        Returns:
            Structure containing molecule in a box.
        """
        coords = np.array(self.cart_coords)
        x_range = max(coords[:, 0]) - min(coords[:, 0])
        y_range = max(coords[:, 1]) - min(coords[:, 1])
        z_range = max(coords[:, 2]) - min(coords[:, 2])
        if a <= x_range or b <= y_range or c <= z_range:
            raise ValueError("Box is not big enough to contain Molecule.")
        lattice = Lattice.from_parameters(a * images[0], b * images[1],
                                          c * images[2],
                                          90, 90, 90)
        nimages = images[0] * images[1] * images[2]
        coords = []

        centered_coords = self.cart_coords - self.center_of_mass
        for i, j, k in itertools.product(list(range(images[0])), list(range(images[1])),
                                         list(range(images[2]))):
            box_center = [(i + 0.5) * a, (j + 0.5) * b, (k + 0.5) * c]
            if random_rotation:
                while True:
                    op = SymmOp.from_origin_axis_angle(
                        (0, 0, 0), axis=np.random.rand(3),
                        angle=random.uniform(-180, 180))
                    m = op.rotation_matrix
                    new_coords = np.dot(m, centered_coords.T).T + box_center
                    if len(coords) == 0:
                        break
                    distances = lattice.get_all_distances(
                        lattice.get_fractional_coords(new_coords),
                        lattice.get_fractional_coords(coords))
                    if np.amin(distances) > min_dist:
                        break
            else:
                new_coords = centered_coords + box_center
            coords.extend(new_coords)
        sprops = {k: v * nimages for k, v in self.site_properties.items()}

        if cls is None: cls = Structure
        return cls(lattice, self.species * nimages, coords,
                   coords_are_cartesian=True,
                   site_properties=sprops).get_sorted_structure()

    def get_centered_molecule(self):
        """
        Returns a Molecule centered at the center of mass.

        Returns:
            Molecule centered with center of mass at origin.
        """
        center = self.center_of_mass
        new_coords = np.array(self.cart_coords) - center
        return self.__class__(self.species_and_occu, new_coords,
                              charge=self._charge,
                              spin_multiplicity=self._spin_multiplicity,
                              site_properties=self.site_properties)

    def to(self, fmt=None, filename=None):
        """
        Outputs the molecule to a file or string.

        Args:
            fmt (str): Format to output to. Defaults to JSON unless filename
                is provided. If fmt is specifies, it overrides whatever the
                filename is. Options include "xyz", "gjf", "g03", "json". If
                you have OpenBabel installed, any of the formats supported by
                OpenBabel. Non-case sensitive.
            filename (str): If provided, output will be written to a file. If
                fmt is not specified, the format is determined from the
                filename. Defaults is None, i.e. string output.

        Returns:
            (str) if filename is None. None otherwise.
        """
        from pymatgen.io.xyz import XYZ
        from pymatgen.io.gaussian import GaussianInput
        from pymatgen.io.babel import BabelMolAdaptor

        fmt = "" if fmt is None else fmt.lower()
        fname = os.path.basename(filename or "")
        if fmt == "xyz" or fnmatch(fname.lower(), "*.xyz*"):
            writer = XYZ(self)
        elif any([fmt == r or fnmatch(fname.lower(), "*.{}*".format(r))
                  for r in ["gjf", "g03", "g09", "com", "inp"]]):
            writer = GaussianInput(self)
        elif fmt == "json" or fnmatch(fname, "*.json*") or fnmatch(fname,
                                                                  "*.mson*"):
            if filename:
                with zopen(filename, "wt", encoding='utf8') as f:
                    return json.dump(self.as_dict(), f)
            else:
                return json.dumps(self.as_dict())
        elif fmt == "yaml" or fnmatch(fname, "*.yaml*"):
            if filename:
                with zopen(fname, "wt", encoding='utf8') as f:
                    return yaml.dump(self.as_dict(), f, Dumper=Dumper)
            else:
                return yaml.dump(self.as_dict(), Dumper=Dumper)

        else:
            m = re.search("\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                          fname.lower())
            if (not fmt) and m:
                fmt = m.group(1)
            writer = BabelMolAdaptor(self)
            return writer.write_file(filename, file_format=fmt)

        if filename:
            writer.write_file(filename)
        else:
            return str(writer)

    @classmethod
    def from_str(cls, input_string, fmt):
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
        from pymatgen.io.xyz import XYZ
        from pymatgen.io.gaussian import GaussianInput
        if fmt.lower() == "xyz":
            m = XYZ.from_string(input_string).molecule
        elif fmt in ["gjf", "g03", "g09", "com", "inp"]:
            m = GaussianInput.from_string(input_string).molecule
        elif fmt == "json":
            d = json.loads(input_string)
            return cls.from_dict(d)
        elif fmt == "yaml":
            d = yaml.load(input_string, Loader=Loader)
            return cls.from_dict(d)
        else:
            from pymatgen.io.babel import BabelMolAdaptor
            m = BabelMolAdaptor.from_string(input_string,
                                            file_format=fmt).pymatgen_mol
        return cls.from_sites(m)

    @classmethod
    def from_file(cls, filename):
        """
        Reads a molecule from a file. Supported formats include xyz,
        gaussian input (gjf|g03|g09|com|inp), Gaussian output (.out|and
        pymatgen's JSON serialized molecules. Using openbabel,
        many more extensions are supported but requires openbabel to be
        installed.

        Args:
            filename (str): The filename to read from.

        Returns:
            Molecule
        """
        from pymatgen.io.gaussian import GaussianOutput
        with zopen(filename) as f:
            contents = f.read()
        fname = filename.lower()
        if fnmatch(fname, "*.xyz*"):
            return cls.from_str(contents, fmt="xyz")
        elif any([fnmatch(fname.lower(), "*.{}*".format(r))
                  for r in ["gjf", "g03", "g09", "com", "inp"]]):
            return cls.from_str(contents, fmt="g09")
        elif any([fnmatch(fname.lower(), "*.{}*".format(r))
                  for r in ["out", "lis", "log"]]):
            return GaussianOutput(filename).final_structure
        elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
            return cls.from_str(contents, fmt="json")
        elif fnmatch(fname, "*.yaml*"):
            return cls.from_str(contents, fmt="yaml")
        else:
            from pymatgen.io.babel import BabelMolAdaptor
            m = re.search("\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                          filename.lower())
            if m:
                new = BabelMolAdaptor.from_file(filename,
                                                m.group(1)).pymatgen_mol
                new.__class__ = cls
                return new

        raise ValueError("Unrecognized file extension!")


class Structure(IStructure, collections.MutableSequence):
    """
    Mutable version of structure.
    """
    __hash__ = None

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):
        """
        Create a periodic structure.

        Args:
            lattice: The lattice, either as a pymatgen.core.lattice.Lattice or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species: List of species on each site. Can take in flexible input,
                including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            fractional_coords: list of fractional coordinates of each species.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
        """
        super(Structure, self).__init__(lattice, species, coords,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties)

        self._sites = list(self._sites)

    def __setitem__(self, i, site):
        """
        Modify a site in the structure.

        Args:
            i (int): Index
            site (PeriodicSite/Specie/Sequence): Three options exist. You
                can provide a PeriodicSite directly (lattice will be
                checked). Or more conveniently, you can provide a
                specie-like object or a tuple of up to length 3. Examples:
                s[0] = "Fe"
                s[0] = Element("Fe")
                both replaces the species only.
                s[0] = "Fe", [0.5, 0.5, 0.5]
                Replaces site and *fractional* coordinates. Any properties
                are inherited from current site.
                s[0] = "Fe", [0.5, 0.5, 0.5], {"spin": 2}
                Replaces site and *fractional* coordinates and properties.
        """
        if isinstance(site, PeriodicSite):
            if site.lattice != self._lattice:
                raise ValueError("PeriodicSite added must have same lattice "
                                 "as Structure!")
            self._sites[i] = site
        else:
            if isinstance(site, six.string_types) or (not isinstance(site, \
                    collections.Sequence)):
                sp = site
                frac_coords = self._sites[i].frac_coords
                properties = self._sites[i].properties
            else:
                sp = site[0]
                frac_coords = site[1] if len(site) > 1 else self._sites[i]\
                    .frac_coords
                properties = site[2] if len(site) > 2 else self._sites[i]\
                    .properties

            self._sites[i] = PeriodicSite(sp, frac_coords, self._lattice,
                                          properties=properties)

    def __delitem__(self, i):
        """
        Deletes a site from the Structure.
        """
        self._sites.__delitem__(i)

    def append(self, species, coords, coords_are_cartesian=False,
               validate_proximity=False, properties=None):
        """
        Append a site to the structure.

        Args:
            species: Species of inserted site
            coords (3x1 array): Coordinates of inserted site
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.

        Returns:
            New structure with inserted site.
        """
        return self.insert(len(self), species, coords,
                           coords_are_cartesian=coords_are_cartesian,
                           validate_proximity=validate_proximity,
                           properties=properties)

    def insert(self, i, species, coords, coords_are_cartesian=False,
               validate_proximity=False, properties=None):
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

        Returns:
            New structure with inserted site.
        """
        if not coords_are_cartesian:
            new_site = PeriodicSite(species, coords, self._lattice,
                                    properties=properties)
        else:
            frac_coords = self._lattice.get_fractional_coords(coords)
            new_site = PeriodicSite(species, frac_coords, self._lattice,
                                    properties=properties)

        if validate_proximity:
            for site in self:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing "
                                     "site!")

        self._sites.insert(i, new_site)

    def add_site_property(self, property_name, values):
        """
        Adds a property to all sites.

        Args:
            property_name (str): The name of the property to add.
            values: A sequence of values. Must be same length as number of
                sites.
        """
        if len(values) != len(self._sites):
            raise ValueError("Values must be same length as sites.")
        for i in range(len(self._sites)):
            site = self._sites[i]
            props = site.properties
            if not props:
                props = {}
            props[property_name] = values[i]
            self._sites[i] = PeriodicSite(site.species_and_occu,
                                          site.frac_coords, self._lattice,
                                          properties=props)

    def replace_species(self, species_mapping):
        """
        Swap species in a structure.

        Args:
            species_mapping (dict): Dict of species to swap. Species can be
                elements too. e.g., {Element("Li"): Element("Na")} performs
                a Li for Na substitution. The second species can be a
                sp_and_occu dict. For example, a site with 0.5 Si that is
                passed the mapping {Element('Si): {Element('Ge'):0.75,
                Element('C'):0.25} } will have .375 Ge and .125 C.
        """
        latt = self._lattice
        species_mapping = {get_el_sp(k): v
                           for k, v in species_mapping.items()}

        def mod_site(site):
            c = Composition()
            for sp, amt in site.species_and_occu.items():
                new_sp = species_mapping.get(sp, sp)
                if isinstance(new_sp, collections.Mapping):
                    c += Composition(new_sp) * amt
                else:
                    c += {new_sp: amt}
            return PeriodicSite(c, site.frac_coords, latt,
                                properties=site.properties)

        self._sites = list(map(mod_site, self._sites))

    def replace(self, i, species, coords=None, coords_are_cartesian=False,
                properties=None):
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
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.
        """
        if coords is None:
            frac_coords = self[i].frac_coords
        elif coords_are_cartesian:
            frac_coords = self._lattice.get_fractional_coords(coords)
        else:
            frac_coords = coords

        new_site = PeriodicSite(species, frac_coords, self._lattice,
                                properties=properties)
        self._sites[i] = new_site

    def remove_species(self, species):
        """
        Remove all occurrences of several species from a structure.

        Args:
            species: Sequence of species to remove, e.g., ["Li", "Na"].
        """
        new_sites = []
        species = list(map(get_el_sp, species))

        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species_and_occu.items()
                           if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(PeriodicSite(
                    new_sp_occu, site.frac_coords, self._lattice,
                    properties=site.properties))
        self._sites = new_sites

    def remove_sites(self, indices):
        """
        Delete sites with at indices.

        Args:
            indices: Sequence of indices of sites to delete.
        """
        self._sites = [s for i, s in enumerate(self._sites)
                       if i not in indices]

    def apply_operation(self, symmop, fractional=False):
        """
        Apply a symmetry operation to the structure and return the new
        structure. The lattice is operated by the rotation matrix only.
        Coords are operated in full and then transformed to the new lattice.

        Args:
            symmop (SymmOp): Symmetry operation to apply.
            fractional (bool): Whether the symmetry operation is applied in
                fractional space. Defaults to False, i.e., symmetry operation
                is applied in cartesian coordinates.
        """
        if not fractional:
            self._lattice = Lattice([symmop.apply_rotation_only(row)
                                     for row in self._lattice.matrix])

            def operate_site(site):
                new_cart = symmop.operate(site.coords)
                new_frac = self._lattice.get_fractional_coords(new_cart)
                return PeriodicSite(site.species_and_occu, new_frac, self._lattice,
                                    properties=site.properties)

        else:
            new_latt = np.dot(symmop.rotation_matrix, self._lattice.matrix)
            self._lattice = Lattice(new_latt)

            def operate_site(site):
                return PeriodicSite(site.species_and_occu,
                                    symmop.operate(site.frac_coords),
                                    self._lattice,
                                    properties=site.properties)

        self._sites = [operate_site(s) for s in self._sites]



    def modify_lattice(self, new_lattice):
        """
        Modify the lattice of the structure.  Mainly used for changing the
        basis.

        Args:
            new_lattice (Lattice): New lattice
        """
        self._lattice = new_lattice
        new_sites = []
        for site in self._sites:
            new_sites.append(PeriodicSite(site.species_and_occu,
                                          site.frac_coords,
                                          self._lattice,
                                          properties=site.properties))
        self._sites = new_sites

    def apply_strain(self, strain):
        """
        Apply a strain to the lattice.

        Args:
            strain (float or list): Amount of strain to apply. Can be a float,
                or a sequence of 3 numbers. E.g., 0.01 means all lattice
                vectors are increased by 1%. This is equivalent to calling
                modify_lattice with a lattice with lattice parameters that
                are 1% larger.
        """
        s = (1 + np.array(strain)) * np.eye(3)
        self.modify_lattice(Lattice(np.dot(self._lattice.matrix.T, s).T))

    def sort(self, key=None, reverse=False):
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
        self._sites = sorted(self._sites, key=key, reverse=reverse)

    def translate_sites(self, indices, vector, frac_coords=True,
                        to_unit_cell=True):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            indices: Integer or List of site indices on which to perform the
                translation.
            vector: Translation vector for sites.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.
            to_unit_cell (bool): Whether new sites are transformed to unit
                cell
        """
        if not isinstance(indices, collections.Iterable):
            indices = [indices]

        for i in indices:
            site = self._sites[i]
            if frac_coords:
                fcoords = site.frac_coords + vector
            else:
                fcoords = self._lattice.get_fractional_coords(site.coords
                                                              + vector)
            new_site = PeriodicSite(site.species_and_occu, fcoords,
                                    self._lattice, to_unit_cell=to_unit_cell,
                                    coords_are_cartesian=False,
                                    properties=site.properties)
            self._sites[i] = new_site

    def perturb(self, distance):
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance (float): Distance in angstroms by which to perturb each
                site.
        """
        def get_rand_vec():
            #deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            return vector / vnorm * distance if vnorm != 0 else get_rand_vec()

        for i in range(len(self._sites)):
            self.translate_sites([i], get_rand_vec(), frac_coords=False)

    def add_oxidation_state_by_element(self, oxidation_states):
        """
        Add oxidation states to a structure.

        Args:
            oxidation_states (dict): Dict of oxidation states.
                E.g., {"Li":1, "Fe":2, "P":5, "O":-2}
        """
        try:
            for i, site in enumerate(self._sites):
                new_sp = {}
                for el, occu in site.species_and_occu.items():
                    sym = el.symbol
                    new_sp[Specie(sym, oxidation_states[sym])] = occu
                new_site = PeriodicSite(new_sp, site.frac_coords,
                                        self._lattice,
                                        coords_are_cartesian=False,
                                        properties=site.properties)
                self._sites[i] = new_site

        except KeyError:
            raise ValueError("Oxidation state of all elements must be "
                             "specified in the dictionary.")

    def add_oxidation_state_by_site(self, oxidation_states):
        """
        Add oxidation states to a structure by site.

        Args:
            oxidation_states (list): List of oxidation states.
                E.g., [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2]
        """
        try:
            for i, site in enumerate(self._sites):
                new_sp = {}
                for el, occu in site.species_and_occu.items():
                    sym = el.symbol
                    new_sp[Specie(sym, oxidation_states[i])] = occu
                new_site = PeriodicSite(new_sp, site.frac_coords,
                                        self._lattice,
                                        coords_are_cartesian=False,
                                        properties=site.properties)
                self._sites[i] = new_site

        except IndexError:
            raise ValueError("Oxidation state of all sites must be "
                             "specified in the dictionary.")

    def remove_oxidation_states(self):
        """
        Removes oxidation states from a structure.
        """
        for i, site in enumerate(self._sites):
            new_sp = collections.defaultdict(float)
            for el, occu in site.species_and_occu.items():
                sym = el.symbol
                new_sp[Element(sym)] += occu
            new_site = PeriodicSite(new_sp, site.frac_coords,
                                    self._lattice,
                                    coords_are_cartesian=False,
                                    properties=site.properties)
            self._sites[i] = new_site

    def make_supercell(self, scaling_matrix):
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
        """
        scale_matrix = np.array(scaling_matrix, np.int16)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
        new_lattice = Lattice(np.dot(scale_matrix, self._lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        for site in self:
            for v in c_lat:
                s = PeriodicSite(site.species_and_occu, site.coords + v,
                                 new_lattice, properties=site.properties,
                                 coords_are_cartesian=True, to_unit_cell=True)
                new_sites.append(s)
        self._sites = new_sites
        self._lattice = new_lattice

    def scale_lattice(self, volume):
        """
        Performs a scaling of the lattice vectors so that length proportions
        and angles are preserved.

        Args:
            volume (float): New volume of the unit cell in A^3.
        """
        self.modify_lattice(self._lattice.scale(volume))

    def merge_sites(self, tol=0.01):
        """
        Merges sites (adding occupancies) within tol of each other.
        Removes site properties
        """
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import fcluster, linkage

        d = self.distance_matrix
        np.fill_diagonal(d, 0)
        clusters = fcluster(linkage(squareform((d + d.T) / 2)),
                            tol, 'distance')

        sites = []
        for c in np.unique(clusters):
            inds = np.argwhere(clusters == c)
            species = self[inds[0]].species_and_occu
            coords = self[inds[0]].frac_coords
            for n, i in enumerate(inds[1:]):
                species += self[i].species_and_occu
                offset = self[i].frac_coords - coords
                coords += (offset - np.round(offset)) / (n + 2)
            sites.append(PeriodicSite(species, coords, self.lattice))

        self._sites = sites




class Molecule(IMolecule, collections.MutableSequence):
    """
    Mutable Molecule. It has all the methods in IMolecule, but in addition,
    it allows a user to perform edits on the molecule.
    """
    __hash__ = None

    def __init__(self, species, coords, charge=0,
                 spin_multiplicity=None, validate_proximity=False,
                 site_properties=None):
        """
        Creates a MutableMolecule.

        Args:
            species: list of atomic species. Possible kinds of input include a
                list of dict of elements/species and occupancies, a List of
                elements/specie specified as actual Element/Specie, Strings
                ("Fe", "Fe2+") or atomic numbers (1,56).
            coords (3x1 array): list of cartesian coordinates of each species.
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
        """
        super(Molecule, self).__init__(species, coords, charge=charge,
            spin_multiplicity=spin_multiplicity,
            validate_proximity=validate_proximity,
            site_properties=site_properties)
        self._sites = list(self._sites)

    def __setitem__(self, i, site):
        """
        Modify a site in the molecule.

        Args:
            i (int): Index
            site (PeriodicSite/Specie/Sequence): Three options exist. You can
                provide a Site directly, or for convenience, you can provide
                simply a Specie-like string/object, or finally a (Specie,
                coords) sequence, e.g., ("Fe", [0.5, 0.5, 0.5]).
        """
        if isinstance(site, Site):
            self._sites[i] = site
        else:
            if isinstance(site, six.string_types) or (not isinstance(site, \
                    collections.Sequence)):
                sp = site
                coords = self._sites[i].coords
                properties = self._sites[i].properties
            else:
                sp = site[0]
                coords = site[1] if len(site) > 1 else self._sites[i].coords
                properties = site[2] if len(site) > 2 else self._sites[i]\
                    .properties

            self._sites[i] = Site(sp, coords, properties=properties)

    def __delitem__(self, i):
        """
        Deletes a site from the Structure.
        """
        self._sites.__delitem__(i)

    def append(self, species, coords, validate_proximity=True,
               properties=None):
        """
        Appends a site to the molecule.

        Args:
            species: Species of inserted site
            coords: Coordinates of inserted site
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to True.
            properties (dict): A dict of properties for the Site.

        Returns:
            New molecule with inserted site.
        """
        return self.insert(len(self), species, coords,
                           validate_proximity=validate_proximity,
                           properties=properties)

    def set_charge_and_spin(self, charge, spin_multiplicity=None):
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
        nelectrons = 0
        for site in self._sites:
            for sp, amt in site.species_and_occu.items():
                nelectrons += sp.Z * amt
        nelectrons -= charge
        self._nelectrons = nelectrons
        if spin_multiplicity:
            if (nelectrons + spin_multiplicity) % 2 != 1:
                raise ValueError(
                    "Charge of {} and spin multiplicity of {} is"
                    " not possible for this molecule".format(
                    self._charge, spin_multiplicity))
            self._spin_multiplicity = spin_multiplicity
        else:
            self._spin_multiplicity = 1 if nelectrons % 2 == 0 else 2

    def insert(self, i, species, coords, validate_proximity=False,
               properties=None):
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
                    raise ValueError("New site is too close to an existing "
                                     "site!")
        self._sites.insert(i, new_site)

    def add_site_property(self, property_name, values):
        """
        Adds a property to a site.

        Args:
            property_name (str): The name of the property to add.
            values (list): A sequence of values. Must be same length as
                number of sites.
        """
        if len(values) != len(self._sites):
            raise ValueError("Values must be same length as sites.")
        for i in range(len(self._sites)):
            site = self._sites[i]
            props = site.properties
            if not props:
                props = {}
            props[property_name] = values[i]
            self._sites[i] = Site(site.species_and_occu, site.coords,
                                  properties=props)

    def replace_species(self, species_mapping):
        """
        Swap species in a molecule.

        Args:
            species_mapping (dict): dict of species to swap. Species can be
                elements too. E.g., {Element("Li"): Element("Na")} performs
                a Li for Na substitution. The second species can be a
                sp_and_occu dict. For example, a site with 0.5 Si that is
                passed the mapping {Element('Si): {Element('Ge'):0.75,
                Element('C'):0.25} } will have .375 Ge and .125 C.
        """
        species_mapping = {get_el_sp(k): v
                           for k, v in species_mapping.items()}

        def mod_site(site):
            new_atom_occu = dict()
            for sp, amt in site.species_and_occu.items():
                if sp in species_mapping:
                    if isinstance(species_mapping[sp], (Element, Specie)):
                        if species_mapping[sp] in new_atom_occu:
                            new_atom_occu[species_mapping[sp]] += amt
                        else:
                            new_atom_occu[species_mapping[sp]] = amt
                    elif isinstance(species_mapping[sp], collections.Mapping):
                        for new_sp, new_amt in species_mapping[sp].items():
                            if new_sp in new_atom_occu:
                                new_atom_occu[new_sp] += amt * new_amt
                            else:
                                new_atom_occu[new_sp] = amt * new_amt
                else:
                    if sp in new_atom_occu:
                        new_atom_occu[sp] += amt
                    else:
                        new_atom_occu[sp] = amt
            return Site(new_atom_occu, site.coords, properties=site.properties)
        self._sites = list(map(mod_site, self._sites))

    def remove_species(self, species):
        """
        Remove all occurrences of a species from a molecule.

        Args:
            species: Species to remove.
        """
        new_sites = []
        species = list(map(get_el_sp, species))
        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species_and_occu.items()
                           if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(Site(new_sp_occu, site.coords,
                                      properties=site.properties))
        self._sites = new_sites

    def remove_sites(self, indices):
        """
        Delete sites with at indices.

        Args:
            indices: Sequence of indices of sites to delete.
        """
        self._sites = [self._sites[i] for i in range(len(self._sites))
                       if i not in indices]

    def translate_sites(self, indices, vector):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            sites (list): List of site indices on which to perform the
                translation.
            vector (3x1 array): Translation vector for sites.
        """
        for i in indices:
            site = self._sites[i]
            new_site = Site(site.species_and_occu, site.coords + vector,
                            properties=site.properties)
            self._sites[i] = new_site

    def perturb(self, distance):
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance (float): Distance in angstroms by which to perturb each
                site.
        """
        def get_rand_vec():
            #deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            return vector / vnorm * distance if vnorm != 0 else get_rand_vec()

        for i in range(len(self._sites)):
            self.translate_sites([i], get_rand_vec())

    def apply_operation(self, symmop):
        """
        Apply a symmetry operation to the molecule.

        Args:
            symmop (SymmOp): Symmetry operation to apply.
        """
        def operate_site(site):
            new_cart = symmop.operate(site.coords)
            return Site(site.species_and_occu, new_cart,
                        properties=site.properties)
        self._sites = list(map(operate_site, self._sites))

    def copy(self):
        """
        Convenience method to get a copy of the molecule.

        Returns:
            A copy of the Molecule.
        """
        return self.__class__.from_sites(self)

    def substitute(self, index, func_grp, bond_order=1):
        """
        Substitute atom at index with a functional group.

        Args:
            index (int): Index of atom to substitute.
            func_grp: Substituent molecule. There are two options:

                1. Providing an actual molecule as the input. The first atom
                   must be a DummySpecie X, indicating the position of
                   nearest neighbor. The second atom must be the next
                   nearest atom. For example, for a methyl group
                   substitution, func_grp should be X-CH3, where X is the
                   first site and C is the second site. What the code will
                   do is to remove the index site, and connect the nearest
                   neighbor to the C atom in CH3. The X-C bond indicates the
                   directionality to connect the atoms.
                2. A string name. The molecule will be obtained from the
                   relevant template in functional_groups.json.
            bond_order: A specified bond order to calculate the bond length
                between the attached functional group and the nearest
                neighbor site. Defaults to 1.
        """

        # Find the nearest neighbor that is not a terminal atom.
        all_non_terminal_nn = []
        for nn, dist in self.get_neighbors(self[index], 3):
            # Check that the nn has neighbors within a sensible distance but
            # is not the site being substituted.
            for inn, dist2 in self.get_neighbors(nn, 3):
                if inn != self[index] and \
                        dist2 < 1.2 * get_bond_length(nn.specie, inn.specie):
                    all_non_terminal_nn.append((nn, dist))
                    break

        if len(all_non_terminal_nn) == 0:
            raise RuntimeError("Can't find a non-terminal neighbor to attach"
                               " functional group to.")

        non_terminal_nn = min(all_non_terminal_nn, key=lambda d: d[1])[0]

        # Set the origin point to be the coordinates of the nearest
        # non-terminal neighbor.
        origin = non_terminal_nn.coords

        # Pass value of functional group--either from user-defined or from
        # functional.json
        if isinstance(func_grp, Molecule):
            func_grp = func_grp
        else:
            # Check to see whether the functional group is in database.
            func_dic = FunctionalGroups()
            if func_grp not in func_dic:
                raise RuntimeError("Can't find functional group in list. "
                                   "Provide explicit coordinate instead")
            else:
                func_grp = func_dic[func_grp]

        # If a bond length can be found, modify func_grp so that the X-group
        # bond length is equal to the bond length.
        bl = get_bond_length(non_terminal_nn.specie, func_grp[1].specie,
                             bond_order=bond_order)
        if bl is not None:
            func_grp = func_grp.copy()
            vec = func_grp[0].coords - func_grp[1].coords
            func_grp[0] = "X", func_grp[1].coords + bl / np.linalg.norm(vec)\
                               * vec

        # Align X to the origin.
        x = func_grp[0]
        func_grp.translate_sites(list(range(len(func_grp))), origin - x.coords)

        #Find angle between the attaching bond and the bond to be replaced.
        v1 = func_grp[1].coords - origin
        v2 = self[index].coords - origin
        angle = get_angle(v1, v2)

        if 1 < abs(angle % 180) < 179:
            # For angles which are not 0 or 180, we perform a rotation about
            # the origin along an axis perpendicular to both bonds to align
            # bonds.
            axis = np.cross(v1, v2)
            op = SymmOp.from_origin_axis_angle(origin, axis, angle)
            func_grp.apply_operation(op)
        elif abs(abs(angle) - 180) < 1:
            # We have a 180 degree angle. Simply do an inversion about the
            # origin
            for i in range(len(func_grp)):
                func_grp[i] = (func_grp[i].species_and_occu,
                               origin - (func_grp[i].coords - origin))

        # Remove the atom to be replaced, and add the rest of the functional
        # group.
        del self[index]
        for site in func_grp[1:]:
            self._sites.append(site)


class StructureError(Exception):
    """
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    """
    pass


@singleton
class FunctionalGroups(dict):

    def __init__(self):
        """
        Loads functional group data from json file. Return list that can be
        easily converted into a Molecule object. The .json file, of course,
        has to be under the same directory of this function
        """
        dict.__init__(self)
        with open(os.path.join(os.path.dirname(__file__),
                               "func_groups.json"), "rt") as f:
            for k, v in json.load(f).items():
                self[k] = Molecule(v["species"], v["coords"])
