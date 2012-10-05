#!/usr/bin/env python

"""
This module provides classes used to define a non-periodic molecule and a
periodic structure.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import re
import abc
import math
import collections
import itertools
from fractions import gcd

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, Specie, \
    smart_element_or_specie
from pymatgen.util.string_utils import formula_double_format
from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.core.bonds import CovalentBond


class SiteCollection(collections.Sequence, collections.Hashable):
    __metaclass__ = abc.ABCMeta

    """
    Basic SiteCollection. Essentially a sequence of Sites or PeriodicSites.
    This serves as a base class for Molecule (a collection of Site, i.e., no
    periodicity) and Structure (a collection of PeriodicSites, i.e.,
    periodicity). Not meant to be instantiated directly.
    """

    """
    Tolerance in Angstrom for determining if sites are too close.
    """
    DISTANCE_TOLERANCE = 0.01

    @abc.abstractproperty
    def sites(self):
        """
        Returns a tuple of sites in the Structure.
        """
        return

    @abc.abstractmethod
    def get_distance(self, i, j):
        """
        Returns distance between sites at index i and j.
        """
        return

    @property
    def distance_matrix(self):
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this should return the nearest image distance.
        """
        nsites = len(self)
        distmatrix = np.zeros((nsites, nsites))
        for i, j in itertools.combinations(xrange(nsites), 2):
            dist = self.get_distance(i, j)
            distmatrix[i, j] = dist
            distmatrix[j, i] = dist
        return distmatrix

    @property
    def species(self):
        """
        List of species at each site of the structure.
        Only works for ordered structures.
        Disordered structures will raise an AttributeError.
        """
        return [site.specie for site in self]

    @property
    def species_and_occu(self):
        """
        List of species and occupancies at each site of the structure.
        """
        return [site.species_and_occu for site in self]

    @property
    def site_properties(self):
        """
        Returns the site properties as a dict of sequences. E.g.,
        {"magmom": (5,-5), "charge": (-4,4)}.
        """
        props = collections.defaultdict(list)
        for site in self:
            for k, v in site.properties.items():
                props[k].append(v)
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
        #for now, just use the composition hash code.
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
        return [site.coords for site in self]

    @property
    def formula(self):
        """
        Returns the formula.
        """
        return self.composition.formula

    @property
    def composition(self):
        """
        Returns the composition
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
        for site in self:
            if not site.is_ordered:
                return False
        return True

    def get_angle(self, i, j, k):
        """
        Returns angle specified by three sites.

        Args:
            i:
                Index of first site
            j:
                Index of second site
            k:
                Index of third site

        Returns:
            Angle in degrees.
        """
        v1 = self[i].coords - self[j].coords
        v2 = self[k].coords - self[j].coords
        ans = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
        """
        Corrects for stupid numerical error which may result in acos being
        operated on a number with absolute value larger than 1
        """
        if ans > 1:
            ans = 1
        elif ans < -1:
            ans = -1
        return math.acos(ans) * 180 / math.pi

    def get_dihedral(self, i, j, k, l):
        """
        Returns dihedral angle specified by four sites.

        Args:
            i:
                Index of first site
            j:
                Index of second site
            k:
                Index of third site
            l:
                Index of fourth site

        Returns:
            Dihedral angle in degrees.
        """
        v1 = self[k].coords - self[l].coords
        v2 = self[j].coords - self[k].coords
        v3 = self[i].coords - self[j].coords
        v23 = np.cross(v2, v3)
        v12 = np.cross(v1, v2)
        return math.atan2(np.linalg.norm(v2) * np.dot(v1, v23),
                          np.dot(v12, v23)) * 180 / math.pi


class Structure(SiteCollection, MSONable):
    """
    Basic Structure object with periodicity. Essentially a sequence of
    PeriodicSites having a common lattice. Structure is made to be immutable
    so that they can function as keys in a dict. Modifications should be done
    by making a new Structure using the structure_modifier module or your own
    methods. Structure extends Sequence and Hashable, which means that in many
    cases, it can be used like any Python sequence. Iterating through a
    structure is equivalent to going through the sites in sequence.
    """

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None):
        """
        Create a periodic structure.

        Args:
            lattice:
                The lattice, either as a pymatgen.core.lattice.Lattice or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species:
                List of species on each site. Can take in flexible input,
                including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            fractional_coords:
                list of fractional coordinates of each species.
            validate_proximity:
                Whether to check if there are sites that are less than 1 Ang
                apart. Defaults to False.
            coords_are_cartesian:
                Set to True if you are providing coordinates in cartesian
                coordinates. Defaults to False.
            site_properties:
                Properties associated with the sites as a dict of sequences,
                e.g., {"magmom":[5,5,5,5]}. The sequences have to be the same
                length as the atomic species and fractional_coords.
                Defaults to None for no properties.
        """
        if len(species) != len(coords):
            raise StructureError(("The list of atomic species must be of the",
                                  "same length as the list of fractional",
                                  " coordinates."))

        if isinstance(lattice, Lattice):
            self._lattice = lattice
        else:
            self._lattice = Lattice(lattice)

        self._sites = []
        for i in xrange(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[i] for k, v in site_properties.items()}
            self._sites.append(PeriodicSite(species[i], coords[i],
                                            self._lattice, to_unit_cell,
                                            coords_are_cartesian,
                                            properties=prop))

        if validate_proximity:
            for (s1, s2) in itertools.combinations(self._sites, 2):
                if s1.distance(s2) < SiteCollection.DISTANCE_TOLERANCE:
                    raise StructureError(("Structure contains sites that are ",
                                          "less than 0.01 Angstrom apart!"))
        self._sites = tuple(self._sites)

    @staticmethod
    def from_sites(sites):
        """
        Convenience constructor to make a Structure from a list of sites.

        Args:
            sites:
                Sequence of PeriodicSites. Sites must have the same lattice.
        """
        props = collections.defaultdict(list)
        lattice = None
        for site in sites:
            if not lattice:
                lattice = site.lattice
            elif site.lattice != lattice:
                raise ValueError("Sites must belong to the same lattice")
            for k, v in site.properties.items():
                props[k].append(v)
        return Structure(lattice,
                         [site.species_and_occu for site in sites],
                         [site.frac_coords for site in sites],
                         site_properties=props)

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
    def density(self):
        """
        Returns the density in units of g/cc
        """
        constant = 1.660468
        return self.composition.weight / self.volume * constant

    @property
    def site_properties(self):
        """
        Returns site properties as a dict of {property: [values]}.
        """
        props = collections.defaultdict(list)
        for site in self._sites:
            for k, v in site.properties.items():
                props[k].append(v)
        return props

    def __eq__(self, other):
        if other is None:
            return False
        if len(self) != len(other):
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
        Returns the fractional coordinates.
        """
        return [site.frac_coords for site in self._sites]

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
            i:
                Index of first site
            j:
                Index of second site
            jimage:
                Number of lattice translations in each lattice direction.

                Default is None for nearest image.

        Returns:
            distance
        """
        return self[i].distance(self[j], jimage)

    def get_sites_in_sphere(self, pt, r):
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
            pt:
                cartesian coordinates of center of sphere.
            r:
                radius of sphere.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        recp_len = self._lattice.reciprocal_lattice.abc
        sr = r + 0.15
        nmax = [sr * l / (2 * math.pi) for l in recp_len]
        pcoords = self._lattice.get_fractional_coords(pt)
        axis_ranges = []
        for i in range(3):
            rangemax = int(math.floor(pcoords[i] + nmax[i]))
            rangemin = int(math.floor(pcoords[i] - nmax[i]))
            axis_ranges.append(range(rangemin, rangemax + 1))
        neighbors = []
        n = len(self._sites)
        site_fcoords = np.array([site.to_unit_cell.frac_coords
                                 for site in self._sites])

        frac_2_cart = self._lattice.get_cartesian_coords
        pts = np.array([pt] * n)
        for image in itertools.product(*axis_ranges):
            shift = np.array([image] * n)
            fcoords = site_fcoords + shift
            coords = frac_2_cart(fcoords)
            dists = np.sqrt(np.sum((coords - pts) ** 2, axis=1))
            withindists = (dists <= r)
            for i in range(n):
                if withindists[i]:
                    nnsite = PeriodicSite(self._sites[i].species_and_occu,
                                          fcoords[i], self._lattice,
                                          properties=self._sites[i].properties)
                    neighbors.append((nnsite, dists[i]))
        return neighbors

    def get_neighbors(self, site, r):
        """
        Get all neighbors to a site within a sphere of radius r.  Excludes the
        site itself.

        Args:
            site:
                site, which is the center of the sphere.
            r:
                radius of sphere.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        nn = self.get_sites_in_sphere(site.coords, r)
        return [(s, dist) for (s, dist) in nn if site != s]

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
            r:
                radius of sphere.

            include_index:
                boolean that determines whether the non-supercell site index
                is included in the returned data

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
        recp_len = self.lattice.reciprocal_lattice.abc
        sr = r + 0.15
        nmax = [sr * l / (2 * math.pi) for l in recp_len]
        site_nminmax = []
        unit_cell_sites = [site.to_unit_cell for site in self._sites]

        for site in self._sites:
            pcoords = site.frac_coords
            inmax = [int(math.floor(pcoords[i] + nmax[i])) for i in xrange(3)]
            inmin = [int(math.floor(pcoords[i] - nmax[i])) for i in xrange(3)]
            site_nminmax.append(zip(inmin, inmax))

        nmin = [min([i[j][0] for i in site_nminmax]) for j in xrange(3)]
        nmax = [max([i[j][1] for i in site_nminmax]) for j in xrange(3)]

        all_ranges = [range(nmin[i], nmax[i] + 1) for i in xrange(3)]

        neighbors = [list() for i in range(len(self._sites))]

        site_coords = np.array(self.cart_coords)
        frac_2_cart = self._lattice.get_cartesian_coords
        n = len(self._sites)
        for image in itertools.product(*all_ranges):
            for (j, site) in enumerate(unit_cell_sites):
                fcoords = site.frac_coords + np.array(image)
                coords = frac_2_cart(fcoords)
                submat = [coords] * n
                dists = (site_coords - submat) ** 2
                dists = np.sqrt(dists.sum(axis=1))
                withindists = (dists <= r) * (dists > 1e-8)
                for i in range(n):
                    if withindists[i]:
                        nnsite = PeriodicSite(site.species_and_occu, fcoords,
                                              site.lattice,
                                              properties=site.properties)
                        item = (nnsite, dists[i], j) if include_index\
                            else (nnsite, dists[i])
                        neighbors[i].append(item)
        return neighbors

    def get_neighbors_in_shell(self, origin, r, dr):
        """
        Returns all sites in a shell centered on origin (coords) between radii
        r-dr and r+dr.

        Args:
            origin:
                cartesian coordinates of center of sphere.
            r:
                inner radius of shell.
            dr:
                width of shell.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        outer = self.get_sites_in_sphere(origin, r + dr)
        inner = r - dr
        return [(site, dist) for (site, dist) in outer if dist > inner]

    def get_sorted_structure(self):
        """
        Get a sorted copy of the structure.
        Sites are sorted by the electronegativity of the species.
        """
        sites = sorted(self)
        return Structure.from_sites(sites)

    def copy(self, site_properties=None, sanitize=False):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Args:
            site_properties:
                Properties to add or override. The properties are specified in
                the same way as the constructor, i.e., as a dict of the form
                {property: [values]}. The properties should be in the order of
                the *original* structure if you are performing sanitization.
            sanitize:
                If True, this method will return a sanitized structure.
                Sanitization performs a few things: (i) The sites are sorted
                by electronegativity, (ii) a LLL lattice reduction is carried
                out to obtain a relatively orthogonalized cell, (iii) all
                fractional coords for sites are mapped into the unit cell.

        Returns:
            A copy of the Structure, with optionally new site_properties and
            optionally sanitized.
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)
        if not sanitize:
            return Structure(self._lattice,
                             [site.species_and_occu for site in self],
                             [site.frac_coords for site in self],
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
            return Structure.from_sites(new_sites)

    def interpolate(self, end_structure, nimages=10):
        """
        Interpolate between this structure and end_structure. Useful for
        construction of NEB inputs.

        Args:
            end_structure:
                structure to interpolate between this structure and end.
            nimages:
                number of interpolation images. Defaults to 10 images.

        Returns:
            List of interpolated structures.
        """
        #Check length of structures
        if len(self) != len(end_structure):
            raise ValueError("Structures have different lengths!")

        #Check that both structures have the same lattice
        if not np.allclose(self.lattice.matrix, end_structure.lattice.matrix):
            raise ValueError("Structures with different lattices!")

        #Check that both structures have the same species
        for i in range(0, len(self)):
            if self[i].species_and_occu != end_structure[i].species_and_occu:
                raise ValueError("Different species!\nStructure 1:\n" +
                                 str(self) + "\nStructure 2\n" +
                                 str(end_structure))

        start_coords = np.array(self.frac_coords)
        end_coords = np.array(end_structure.frac_coords)

        vec = end_coords - start_coords
        structs = [Structure(self.lattice,
                             [site.species_and_occu for site in self._sites],
                             start_coords + float(x) / float(nimages) * vec,
                             site_properties=self.site_properties)
                   for x in range(0, nimages + 1)]
        return structs

    def __repr__(self):
        outs = []
        outs.append("Structure Summary")
        outs.append(repr(self.lattice))
        for s in self:
            outs.append(repr(s))
        return "\n".join(outs)

    def __str__(self):
        outs = ["Structure Summary ({s})".format(s=str(self.composition))]
        outs.append("Reduced Formula: {}"
                    .format(self.composition.reduced_formula))
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i + 1), site.species_string,
                                  " ".join([to_s(j).rjust(12)
                                            for j in site.frac_coords])]))
        return "\n".join(outs)

    @property
    def to_dict(self):
        """
        Json-serializable dict representation of Structure
        """
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["lattice"] = self._lattice.to_dict
        d["sites"] = []
        for site in self:
            site_dict = site.to_dict
            del site_dict["lattice"]
            del site_dict["@module"]
            del site_dict["@class"]
            d["sites"].append(site_dict)
        return d

    @staticmethod
    def from_dict(d):
        """
        Reconstitute a Structure object from a dict representation of Structure
        created using to_dict.

        Args:
            d:
                dict representation of structure.

        Returns:
            Structure object
        """
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        return Structure.from_sites(sites)


class Molecule(SiteCollection, MSONable):
    """
    Basic Molecule object without periodicity. Essentially a sequence of sites.
    Molecule is made to be immutable so that they can function as keys in a
    dict. Modifications should be done by making a new Molecule.
    Molecule extends Sequence and Hashable, which means that in many cases,
    it can be used like any Python sequence. Iterating through a molecule is
    equivalent to going through the sites in sequence.
    """

    def __init__(self, species, coords, validate_proximity=False,
                 site_properties=None):
        """
        Creates a Molecule.

        Args:
            species:
                list of atomic species. Possible kinds of input include a list
                of dict of elements/species and occupancies, a List of
                elements/specie specified as actual Element/Specie, Strings
                ("Fe", "Fe2+") or atomic numbers (1,56).
            coords:
                list of cartesian coordinates of each species.
            validate_proximity:
                Whether to check if there are sites that are less than 1 Ang
                apart. Defaults to False.
            site_properties:
                Properties associated with the sites as a dict of sequences,
                e.g., {"magmom":[5,5,5,5]}. The sequences have to be the same
                length as the atomic species and fractional_coords.
                Defaults to None for no properties.
        """
        if len(species) != len(coords):
            raise StructureError(("The list of atomic species must be of the",
                                  " same length as the list of fractional ",
                                  "coordinates."))

        sites = []
        for i in xrange(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[i] for k, v in site_properties.items()}
            sites.append(Site(species[i], coords[i], properties=prop))
        if validate_proximity:
            for (s1, s2) in itertools.combinations(sites, 2):
                if s1.distance(s2) < Structure.DISTANCE_TOLERANCE:
                    raise StructureError(("Molecule contains sites that are ",
                                          "less than 0.01 Angstrom apart!"))
        self._sites = tuple(sites)

    @property
    def sites(self):
        """
        Returns a tuple of sites in the Molecule.
        """
        return self._sites

    @staticmethod
    def from_sites(sites):
        """
        Convenience static constructor to make a Molecule from a list of sites.

        Args:
            sites:
                Sequence of Sites.
        """
        props = collections.defaultdict(list)
        for site in sites:
            for k, v in site.properties.items():
                props[k].append(v)
        return Molecule([site.species_and_occu for site in sites],
                        [site.coords for site in sites],
                        site_properties=props)

    def break_bond(self, ind1, ind2, tol=0.2):
        """
        Returns two molecules based on breaking the bond between atoms at index
        ind1 and ind2.

        Args:
            ind1:
                Index of first site.
            ind2:
                Index of second site.
            tol:
                Relative tolerance to test. Basically, the code checks if the
                distance between the sites is less than (1 + tol) * typical
                bond distances. Defaults to 0.2, i.e., 20% longer.

        Returns:
            Two Molecule objects representing the two clusters formed from
            breaking the bond.
        """
        sites = self._sites
        clusters = [[sites[ind1]], [sites[ind2]]]

        sites = [sites[i] for i in xrange(len(sites)) if i not in (ind1, ind2)]

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
                break
            sites = unmatched

        return (Molecule.from_sites(cluster) for cluster in clusters)

    def get_covalent_bonds(self, tol=0.2):
        """
        Determines the covalent bonds in a molecule.

        Args:
            tol:
                The tol to determine bonds in a structure. See
                CovalentBond.is_bonded.

        Returns:
            List of bonds
        """
        bonds = []
        for site1, site2 in itertools.combinations(self._sites, 2):
            if CovalentBond.is_bonded(site1, site2, tol):
                bonds.append(CovalentBond(site1, site2))
        return bonds

    def __repr__(self):
        outs = []
        outs.append("Molecule Summary")
        for s in self:
            outs.append(repr(s))
        return "\n".join(outs)

    def __str__(self):
        outs = ["Molecule Summary ({s})".format(s=str(self.composition))]
        outs.append("Reduced Formula: " + self.composition.reduced_formula)
        to_s = lambda x: "%0.6f" % x
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i + 1), site.species_string,
                        " ".join([to_s(j).rjust(12) for j in site.coords])]))
        return "\n".join(outs)

    @property
    def to_dict(self):
        """
        Json-serializable dict representation of Molecule
        """
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["sites"] = [site.to_dict for site in self]
        return d

    @staticmethod
    def from_dict(d):
        """
        Reconstitute a Molecule object from a dict representation created using
        to_dict.

        Args:
            d:
                dict representation of Molecule.

        Returns:
            Molecule object
        """
        species = []
        coords = []
        props = {}

        for site_dict in d["sites"]:
            sp = site_dict["species"]
            species.append({Specie(sp["element"], sp["oxidation_state"])
                            if "oxidation_state" in sp else
                            Element(sp["element"]): sp["occu"]
                            for sp in site_dict["species"]})
            coords.append(site_dict["xyz"])
            siteprops = site_dict.get("properties", {})
            for k, v in siteprops.items():
                if k not in props:
                    props[k] = [v]
                else:
                    props[k].append(v)
        return Molecule(species, coords, site_properties=props)

    def get_distance(self, i, j):
        """
        Get distance between site i and j.

        Args:
            i:
                Index of first site
            j:
                Index of second site

        Returns:
            Distance between the two sites.
        """
        return self[i].distance(self[j])

    def get_sites_in_sphere(self, pt, r):
        """
        Find all sites within a sphere from a point.

        Args:
            pt:
                cartesian coordinates of center of sphere.
            r:
                radius of sphere.

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
            site:
                site, which is the center of the sphere.
            r:
                radius of sphere.

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
            origin:
                cartesian coordinates of center of sphere.
            r:
                inner radius of shell.
            dr:
                width of shell.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        outer = self.get_sites_in_sphere(origin, r + dr)
        inner = r - dr
        return [(site, dist) for (site, dist) in outer if dist > inner]

    def get_boxed_structure(self, a, b, c):
        """
        Creates a Structure from a Molecule by putting the Molecule in a box.
        Useful for creating Structure for calculating molecules using periodic
        codes.

        Args:
            a:
                a-lattice parameter.
            b:
                b-lattice parameter.
            c:
                c-lattice parameter.

        Returns:
            Structure containing molecule in a box.
        """
        coords = np.array(self.cart_coords)
        x_range = max(coords[:, 0]) - min(coords[:, 0])
        y_range = max(coords[:, 1]) - min(coords[:, 1])
        z_range = max(coords[:, 2]) - min(coords[:, 2])
        if a <= x_range or b <= y_range or c <= z_range:
            raise ValueError("Box is not big enough to contain Molecule.")
        lattice = Lattice.from_parameters(a, b, c, 90, 90, 90)
        return Structure(lattice, self.species, self.cart_coords,
                         coords_are_cartesian=True,
                         site_properties=self.site_properties)


class StructureError(Exception):
    """
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return "Structure Error : " + self.msg


class Composition (collections.Mapping, collections.Hashable, MSONable):
    """
    Represents a Composition, which is essentially a {element:amount} dict.

    Note that the key can be either an Element or a Specie. Elements and Specie
    are treated differently. i.e., a Fe2+ is not the same as a Fe3+ Specie and
    would be put in separate keys. This differentiation is deliberate to
    support using Composition to determine the fraction of a particular Specie.

    Works almost completely like a standard python dictionary, except that
    __getitem__ is overridden to return 0 when an element is not found.
    (somewhat like a defaultdict, except it is immutable).

    Also adds more convenience methods relevant to compositions, e.g.,
    get_fraction.

    >>> comp = Composition("LiFePO4")
    >>> comp.get_atomic_fraction(Element("Li"))
    0.14285714285714285
    >>> comp.num_atoms
    7.0
    >>> comp.reduced_formula
    "LiFePO4"
    >>> comp.formula
    "Li1 Fe1 P1 O4"
    >>> comp.get_wt_fraction(Element("Li"))
    0.04399794666951898
    >>> comp.num_atoms
    7.0
    """

    """
    Tolerance in distinguishing different composition amounts.
    1e-8 is fairly tight, but should cut out most floating point arithmetic
    errors.
    """
    amount_tolerance = 1e-8

    """
    Special formula handling for peroxides and certain elements. This is so
    that formula output does not write LiO instead of Li2O2 for example.
    """
    special_formulas = {"LiO": "Li2O2", "NaO": "Na2O2", "KO": "K2O2",
                        "HO": "H2O2", "O": "O2", "F": "F2", "N": "N2",
                        "Cl": "Cl2", "H": "H2"}

    def __init__(self, *args, **kwargs):
        """
        Very flexible Composition construction, similar to the built-in Python
        dict(). Also extended to allow simple string init.

        Args:
            Any form supported by the Python built-in dict() function.

            1. A dict of either {Element/Specie: amount},

               {string symbol:amount}, or {atomic number:amount} or any mixture
               of these. E.g., {Element("Li"):2 ,Element("O"):1},
               {"Li":2, "O":1}, {3:2, 8:1} all result in a Li2O composition.
            2. Keyword arg initialization, similar to a dict, e.g.,

               Compostion(Li = 2, O = 1)

            In addition, the Composition constructor also allows a single
            string as an input formula. E.g., Composition("Li2O").
        """
        if len(args) == 1 and isinstance(args[0], basestring):
            elmap = self._parse_formula(args[0])
        else:
            elmap = dict(*args, **kwargs)
        if any([e < 0 for e in elmap.values()]):
            raise ValueError("Amounts in Composition cannot be negative!")
        self._elmap = {smart_element_or_specie(k): v for k, v in elmap.items()}
        self._natoms = sum(self._elmap.values())

    def __getitem__(self, el):
        """
        Get the amount for element.
        """
        return self._elmap.get(el, 0)

    def __eq__(self, other):
        for el in self.elements:
            if self[el] != other[el]:
                return False
        for el in other.elements:
            if self[el] != other[el]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        """
        Adds two compositions. For example, an Fe2O3 composition + an FeO
        composition gives a Fe3O4 composition.
        """
        new_el_map = {el: self[el] for el in self}
        for k in other.keys():
            el = smart_element_or_specie(k)
            if el in self:
                new_el_map[el] += other[k]
            else:
                new_el_map[el] = other[k]
        return Composition(new_el_map)

    def __sub__(self, other):
        """
        Subtracts two compositions. For example, an Fe2O3 composition - an FeO
        composition gives an FeO2 composition.

        Raises:
            ValueError if the subtracted composition is greater than the
            original composition in any of its elements.
        """
        new_el_map = {el: self[el] for el in self}
        for k in other.keys():
            el = smart_element_or_specie(k)
            if el in self and other[k] <= self[el]:
                new_el_map[el] -= other[k]
            else:
                raise ValueError(("All elements in subtracted composition "
                                  "must exist in original composition in "
                                  "equal or lesser amount!"))
        return Composition(new_el_map)

    def __mul__(self, other):
        """
        Multiply a Composition by an integer or a float.
        Fe2O3 * 4 -> Fe8O12
        """
        if not (isinstance(other, int) or isinstance(other, float)):
            raise ValueError("Multiplication can only be done for int/floats!")
        return Composition({el: self[el] * other for el in self})

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between
        Compositions with different elements.
        """
        hashcode = 0
        for el in self._elmap.keys():
            #Ignore elements with zero amounts.
            if self[el] > self.amount_tolerance:
                hashcode += el.Z
        return 7

    def __contains__(self, el):
        return el in self._elmap

    def __len__(self):
        return len(self._elmap)

    def __iter__(self):
        return self._elmap.__iter__()

    @property
    def is_element(self):
        """
        True if composition is for an element.
        """
        positive_amts = [amt for amt in self._elmap.values()
                         if amt > self.amount_tolerance]
        return len(positive_amts) == 1

    def copy(self):
        return Composition(self._elmap)

    @property
    def formula(self):
        """
        Returns a formula string, with elements sorted by electronegativity,
        e.g., Li4 Fe4 P4 O16.
        """
        sym_amt = self.to_dict
        syms = sorted(sym_amt.keys(),
                      key=lambda s: smart_element_or_specie(s).X)
        formula = []
        for s in syms:
            if sym_amt[s] != 0:
                formula.append(s + formula_double_format(sym_amt[s], False))
        return " ".join(formula)

    @property
    def alphabetical_formula(self):
        """
        Returns a formula string, with elements sorted by alphabetically
        e.g., Fe4 Li4 O16 P4.
        """
        sym_amt = self.to_dict
        syms = sorted(sym_amt.keys())
        formula = []
        for s in syms:
            if sym_amt[s] != 0:
                formula.append(s + formula_double_format(sym_amt[s], False))
        return " ".join(formula)

    @property
    def reduced_composition(self):
        """
        Returns the reduced composition,i.e. amounts normalized by greatest
        common denominator. e.g., Composition("FePO4") for
        Composition("Fe4P4O16").
        """
        return self.get_reduced_composition_and_factor()[0]

    def get_reduced_composition_and_factor(self):
        """
        Returns a normalized composition and a multiplicative factor,
        i.e., Li4Fe4P4O16 returns (LiFePO4, 4).
        """
        (formula, factor) = self.get_reduced_formula_and_factor()
        return (Composition.from_formula(formula), factor)

    def get_reduced_formula_and_factor(self):
        """
        Returns a pretty normalized formula and a multiplicative factor, i.e.,
        Li4Fe4P4O16 returns (LiFePO4, 4).
        """
        is_int = lambda x: x == int(x)
        all_int = all([is_int(x) for x in self._elmap.values()])
        if not all_int:
            return (re.sub("\s", "", self.formula), 1)

        sym_amt = self.to_dict
        syms = sorted(sym_amt.keys(),
                      key=lambda s: smart_element_or_specie(s).X)

        syms = filter(lambda s: sym_amt[s] != 0, syms)
        num_el = len(syms)
        contains_polyanion = False
        if num_el >= 3:
            contains_polyanion = (smart_element_or_specie(syms[num_el - 1]).X
                                  - smart_element_or_specie(syms[num_el - 2]).X
                                  < 1.65)

        factor = reduce(gcd, self._elmap.values())
        reduced_form = ""
        n = num_el
        if contains_polyanion:
            n -= 2

        for i in range(0, n):
            s = syms[i]
            normamt = sym_amt[s] * 1.0 / factor
            reduced_form += s + formula_double_format(normamt)

        if contains_polyanion:
            polyamounts = list()
            polyamounts.append(sym_amt[syms[num_el - 2]] / factor)
            polyamounts.append(sym_amt[syms[num_el - 1]] / factor)
            polyfactor = reduce(gcd, polyamounts)
            for i in range(n, num_el):
                s = syms[i]
                normamt = sym_amt[s] / factor / polyfactor
                if normamt != 1.0:
                    if normamt != int(normamt):
                        polyfactor = 1
                        break

            poly_form = ""

            for i in range(n, num_el):
                s = syms[i]
                normamt = sym_amt[s] / factor / polyfactor
                poly_form += s + formula_double_format(normamt)

            if polyfactor != 1:
                reduced_form += "({}){}".format(poly_form, int(polyfactor))
            else:
                reduced_form += poly_form

        if reduced_form in Composition.special_formulas:
            reduced_form = Composition.special_formulas[reduced_form]
            factor = factor / 2

        return (reduced_form, factor)

    def get_fractional_composition(self):
        """
        Returns the normalized composition which the number of species sum to 1.
        
        Returns:
            Normalized composition which the number of species sum to 1.
        """
        natoms = self._natoms
        frac_map = {k:v / natoms for k, v in self._elmap.items()}
        return Composition(frac_map)

    @property
    def reduced_formula(self):
        """
        Returns a pretty normalized formula, i.e., LiFePO4 instead of
        Li4Fe4P4O16.
        """
        return self.get_reduced_formula_and_factor()[0]

    @property
    def elements(self):
        """
        Returns view of elements in Composition.
        """
        return self._elmap.keys()

    def __str__(self):
        return self.formula

    @property
    def num_atoms(self):
        """
        Total number of atoms in Composition
        """
        return self._natoms

    @property
    def weight(self):
        """
        Total molecular weight of Composition
        """
        return sum([amount * el.atomic_mass
                    for el, amount in self._elmap.items()])

    def get_atomic_fraction(self, el):
        """
        Args:
            el:
                Element or Specie

        Returns:
            Atomic fraction for element el in Composition
        """
        return self[el] / self._natoms

    def get_wt_fraction(self, el):
        """
        Args:
            el:
                Element or Specie

        Returns:
            Weight fraction for element el in Composition
        """
        return el.atomic_mass * self[el] / self.weight

    def _parse_formula(self, formula):
        """
        Args:
            formula:
                A string formula, e.g. Fe2O3, Li3Fe2(PO4)3

        Returns:
            Composition with that formula.
        """
        def get_sym_dict(f, factor):
            sym_dict = {}
            for m in re.finditer(r"([A-Z][a-z]*)([\.\d]*)", f):
                el = m.group(1)
                amt = 1
                if m.group(2).strip() != "":
                    amt = float(m.group(2))
                if el in sym_dict:
                    sym_dict[el] += amt * factor
                else:
                    sym_dict[el] = amt * factor
                f = f.replace(m.group(), "", 1)
            if f.strip():
                raise ValueError("{} is an invalid formula!".format(f))
            return sym_dict
        m = re.search(r"\(([^\(\)]+)\)([\.\d]*)", formula)
        if m:
            factor = 1
            if m.group(2) != "":
                factor = float(m.group(2))
            unit_sym_dict = get_sym_dict(m.group(1), factor)
            expanded_sym = "".join(["{}{}".format(el, amt)
                                    for el, amt in unit_sym_dict.items()])
            expanded_formula = formula.replace(m.group(), expanded_sym)
            return self._parse_formula(expanded_formula)
        return get_sym_dict(formula, 1)

    @staticmethod
    def from_formula(formula):
        """
        .. deprecated:: 1.6.1

        Use Composition(formula) instead.
        """
        return Composition(formula)

    @property
    def anonymized_formula(self):
        """
        An anonymized formula. Unique species are arranged in ordering of
        increasing amounts and assigned ascending alphabets. Useful for
        prototyping formulas. For example, all stoichiometric perovskites have
        anonymized_formula ABC3.
        """
        reduced_comp = self.get_reduced_composition_and_factor()[0]
        els = sorted(reduced_comp.elements, key=lambda e: reduced_comp[e])
        ascii_code = 65
        anon_formula = []
        for e in els:
            amt = reduced_comp[e]
            if amt > 0:
                if amt == 1:
                    amt_str = ""
                elif abs(amt % 1) < 1e-8:
                    amt_str = str(int(amt))
                else:
                    amt_str = str(amt)
                anon_formula.append("{}{}".format(chr(ascii_code), amt_str))
                ascii_code += 1
        return "".join(anon_formula)

    def __repr__(self):
        return "Comp: " + self.formula

    @staticmethod
    def from_dict(d):
        """
        Creates a composition from a dict generated by to_dict. Strictly not
        necessary given that the standard constructor already takes in such an
        input, but this method preserves the standard pymatgen API of having
        from_dict methods to reconstitute objects generated by to_dict. Allows
        for easier introspection.

        Args:
            d:
                {symbol: amount} dict.
        """
        return Composition(d)

    @property
    def to_dict(self):
        """
        Returns:
            dict with element symbol and (unreduced) amount e.g.,
            {"Fe": 4.0, "O":6.0}
        """
        d = {}
        for e, a in self.items():
            if e.symbol in d:
                d[e.symbol] += a
            else:
                d[e.symbol] = a
        return d

    @property
    def to_reduced_dict(self):
        """
        Returns:
            dict with element symbol and reduced amount e.g.,
            {"Fe": 2.0, "O":3.0}
        """
        reduced_formula = self.reduced_formula
        c = Composition.from_formula(reduced_formula)
        return c.to_dict

    @property
    def to_data_dict(self):
        """
        Returns a dict with many composition-related properties.

        Returns:
            A dict with many keys and values relating to Composition/Formula
        """
        d = {}
        d["reduced_cell_composition"] = self.to_reduced_dict
        d["unit_cell_composition"] = self.to_dict
        d["reduced_cell_formula"] = self.reduced_formula
        d["elements"] = self.to_dict.keys()
        d["nelements"] = len(self.to_dict.keys())
        return d

    @staticmethod
    def ranked_compositions_from_indeterminate_formula(fuzzy_formula,
                                                       lock_if_strict=True):
        """
        Takes in a formula where capitilization might not be correctly entered,
        and suggests a ranked list of potential Composition matches.
        Author: Anubhav Jain

        Args:
            fuzzy_formula:
                A formula string, such as "co2o3" or "MN", that may or may not

                have multiple interpretations
            lock_if_strict:
                If true, a properly entered formula will only return the one

                correct interpretation. For example, "Co1" will only return

                "Co1" if true, but will return both "Co1" and "C1 O1" if false.

        Returns:
            A ranked list of potential Composition matches
        """

        #if we have an exact match and the user specifies lock_if_strict, just
        #return the exact match!
        if lock_if_strict:
            #the strict composition parsing might throw an error, we can ignore
            #it and just get on with fuzzy matching
            try:
                comp = Composition.from_formula(fuzzy_formula)
                return [comp]
            except:
                pass

        all_matches = Composition._comps_from_fuzzy_formula(fuzzy_formula)
        #remove duplicates
        all_matches = list(set(all_matches))
        #sort matches by rank descending
        all_matches = sorted(all_matches,
                             key=lambda match: match[1], reverse=True)
        all_matches = [m[0] for m in all_matches]
        return all_matches

    @staticmethod
    def _comps_from_fuzzy_formula(fuzzy_formula, m_dict={}, m_points=0,
                                  factor=1):
        """
        A recursive helper method for formula parsing that helps in
        interpreting and ranking indeterminate formulas.
        Author: Anubhav Jain

        Args:
            fuzzy_formula:
                A formula string, such as "co2o3" or "MN", that may or may not
                have multiple interpretations.
            m_dict:
                A symbol:amt dictionary from the previously parsed formula.
            m_points:
                Number of points gained from the previously parsed formula.
            factor:
                Coefficient for this parse, e.g. (PO4)2 will feed in PO4 as the
                fuzzy_formula with a coefficient of 2.

        Returns:
            A list of tuples, with the first element being a Composition and
            the second element being the number of points awarded that
            Composition intepretation.
        """

        def _parse_chomp_and_rank(m, f, m_dict, m_points):
            """
            A helper method for formula parsing that helps in interpreting and
            ranking indeterminate formulas
            Author: Anubhav Jain

            Args:
                m:
                    a regex match, with the first group being the element and
                    the second group being the amount
                f:
                    The formula part containing the match
                m_dict:
                    A symbol:amt dictionary from the previously parsed formula
                m_points:
                    Number of points gained from the previously parsed formula

            Returns:
                A tuple of (f, m_dict, points) where m_dict now contains data
                from the match and the match has been removed (chomped) from
                the formula f. The "goodness" of the match determines the
                number of points returned for chomping. Returns
                (None, None, None) if no element could be found...
            """

            points = 0
            # Points awarded if the first element of the element is correctly
            # specified as a capital
            points_first_capital = 100
            # Points awarded if the second letter of the element is correctly
            # specified as lowercase
            points_second_lowercase = 100

            #get element and amount from regex match
            el = m.group(1)
            if len(el) > 2 or len(el) < 1:
                raise ValueError("Invalid element symbol entered!")
            amt = float(m.group(2)) if m.group(2).strip() != "" else 1

            #convert the element string to proper [uppercase,lowercase] format
            #and award points if it is already in that format
            char1 = el[0]
            char2 = el[1] if len(el) > 1 else ""

            if char1 == char1.upper():
                points += points_first_capital
            if char2 and char2 == char2.lower():
                points += points_second_lowercase

            el = char1.upper() + char2.lower()

            #if it's a valid element, chomp and add to the points
            if Element.is_valid_symbol(el):
                if el in m_dict:
                    m_dict[el] += amt * factor
                else:
                    m_dict[el] = amt * factor
                return (f.replace(m.group(), "", 1), m_dict, m_points + points)

            #else return None
            return (None, None, None)

        fuzzy_formula = fuzzy_formula.strip()

        if len(fuzzy_formula) == 0:
            #The entire formula has been parsed into m_dict. Return the
            #corresponding Composition and number of points
            if m_dict:
                yield (Composition.from_dict(m_dict), m_points)
        else:
            #if there is a parenthesis, remove it and match the remaining stuff
            #with the appropriate factor
            for mp in re.finditer(r"\(([^\(\)]+)\)([\.\d]*)", fuzzy_formula):
                mp_points = m_points
                mp_form = fuzzy_formula.replace(mp.group(), " ", 1)
                mp_dict = dict(m_dict)
                mp_factor = 1 if mp.group(2) == "" else float(mp.group(2))
                #Match the stuff inside the parenthesis with the appropriate
                #factor
                for match in \
                    Composition._comps_from_fuzzy_formula(mp.group(1),
                                                          mp_dict,
                                                          mp_points,
                                                          factor=mp_factor):
                    only_me = True
                    # Match the stuff outside the parentheses and return the
                    # sum.

                    for match2 in \
                        Composition._comps_from_fuzzy_formula(mp_form,
                                                              mp_dict,
                                                              mp_points,
                                                              factor=1):
                        only_me = False
                        yield (match[0] + match2[0], match[1] + match2[1])
                    #if the stuff inside the parenthesis is nothing, then just
                    #return the stuff inside the parentheses
                    if only_me:
                        yield match
                return

            #try to match the single-letter elements
            m1 = re.match(r"([A-z])([\.\d]*)", fuzzy_formula)
            if m1:
                m_points1 = m_points
                m_form1 = fuzzy_formula
                m_dict1 = dict(m_dict)
                (m_form1, m_dict1, m_points1) = \
                    _parse_chomp_and_rank(m1, m_form1, m_dict1, m_points1)
                if m_dict1:
                    #there was a real match
                    for match in \
                        Composition._comps_from_fuzzy_formula(m_form1,
                                                              m_dict1,
                                                              m_points1,
                                                              factor):
                        yield match

            #try to match two-letter elements
            m2 = re.match(r"([A-z]{2})([\.\d]*)", fuzzy_formula)
            if m2:
                m_points2 = m_points
                m_form2 = fuzzy_formula
                m_dict2 = dict(m_dict)
                (m_form2, m_dict2, m_points2) = \
                    _parse_chomp_and_rank(m2, m_form2, m_dict2, m_points2)
                if m_dict2:
                    #there was a real match
                    for match in \
                        Composition._comps_from_fuzzy_formula(m_form2, m_dict2,
                                                              m_points2,
                                                              factor):
                        yield match

if __name__ == "__main__":
    import doctest
    doctest.testmod()
