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

import abc
import math
import collections
import itertools

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.core.bonds import CovalentBond
from pymatgen.core.physical_constants import AMU_TO_KG
from pymatgen.core.composition import Composition
from pymatgen.util.coord_utils import get_points_in_sphere_pbc, pbc_diff


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
        return all((site.is_ordered for site in self))

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

        #Correct for stupid numerical error which may result in acos being
        #operated on a number with absolute value larger than 1
        ans = min(ans, 1)
        ans = max(-1, ans)
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
            raise StructureError("The list of atomic species must be of the"
                                 "same length as the list of fractional"
                                 " coordinates.")

        if isinstance(lattice, Lattice):
            self._lattice = lattice
        else:
            self._lattice = Lattice(lattice)

        sites = []
        for i in xrange(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[i] for k, v in site_properties.items()}
            sites.append(PeriodicSite(species[i], coords[i],
                                      self._lattice, to_unit_cell,
                                      coords_are_cartesian,
                                      properties=prop))

        if validate_proximity:
            for (s1, s2) in itertools.combinations(sites, 2):
                if s1.distance(s2) < SiteCollection.DISTANCE_TOLERANCE:
                    raise StructureError(("Structure contains sites that are ",
                                          "less than 0.01 Angstrom apart!"))
        self._sites = tuple(sites)

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
        constant = AMU_TO_KG * 1000 / 1e-24
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
        if self._lattice != other._lattice:
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
            pt:
                cartesian coordinates of center of sphere.
            r:
                radius of sphere.
            include_index:
                boolean that determines whether the non-supercell site index
                is included in the returned data

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        site_fcoords = np.mod(self.frac_coords, 1)
        neighbors = []
        for fcoord, dist, i in get_points_in_sphere_pbc(self._lattice,
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
        floor = math.floor
        for site in self:
            pcoords = site.frac_coords
            inmax = [int(floor(pcoords[i] + nmax[i])) for i in xrange(3)]
            inmin = [int(floor(pcoords[i] - nmax[i])) for i in xrange(3)]
            site_nminmax.append(zip(inmin, inmax))

        nmin = [min([i[j][0] for i in site_nminmax]) for j in xrange(3)]
        nmax = [max([i[j][1] for i in site_nminmax]) for j in xrange(3)]

        all_ranges = [range(nmin[i], nmax[i] + 1) for i in xrange(3)]

        neighbors = [list() for i in xrange(len(self._sites))]
        all_fcoords = np.mod(self.frac_coords, 1)

        site_coords = np.array(self.cart_coords)
        latt = self._lattice
        frac_2_cart = latt.get_cartesian_coords
        n = len(self)
        indices = np.array(range(n))
        for image in itertools.product(*all_ranges):
            for (j, fcoord) in enumerate(all_fcoords):
                fcoords = fcoord + image
                coords = frac_2_cart(fcoords)
                submat = np.tile(coords, (n, 1))
                dists = (site_coords - submat) ** 2
                dists = np.sqrt(dists.sum(axis=1))
                withindists = (dists <= r) * (dists > 1e-8)
                sp = self[j].species_and_occu
                props = self[j].properties
                for i in indices[withindists]:
                    nnsite = PeriodicSite(sp, fcoords, latt,
                                          properties=props)
                    item = (nnsite, dists[i], j) if include_index else (
                        nnsite, dists[i])
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

    def get_reduced_structure(self, reduction_algo="niggli"):
        """
        Get a reduced structure.

        Args:
            reduction_algo:
                The lattice reduction algorithm to use. Currently supported
                options are "niggli" or "LLL".
        """
        if reduction_algo == "niggli":
            reduced_latt = self._lattice.get_niggli_reduced_lattice()
        elif reduction_algo == "LLL":
            reduced_latt = self._lattice.get_lll_reduced_lattice()
        else:
            raise ValueError("Invalid reduction algo : {}"
            .format(reduction_algo))

        return Structure(reduced_latt, self.species_and_occu, self.cart_coords,
                         coords_are_cartesian=True, to_unit_cell=True)

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


    def get_primitive_structure(self, tolerance=0.5):
        """
        This finds a smaller unit cell than the input. Sometimes it doesn"t
        find the smallest possible one, so this method is recursively called
        until it is unable to find a smaller cell.

        The method works by finding possible smaller translations
        and then using that translational symmetry instead of one of the
        lattice basis vectors if more than one vector is found (usually the
        case for large cells) the one with the smallest norm is used.

        Things are done in fractional coordinates because its easier to
        translate back to the unit cell.

        Args:
            tolerance:
                Tolerance for each coordinate of a particular site. For
                example, [0.5, 0, 0.5] in cartesian coordinates will be
                considered to be on the same coordinates as [0, 0, 0] for a
                tolerance of 0.5. Defaults to 0.5.

        Returns:
            The most primitive structure found. The returned structure is
            guanranteed to have len(new structure) <= len(structure).
        """
        original_volume = self.volume
        (reduced_formula, num_fu) =\
            self.composition.get_reduced_composition_and_factor()

        min_vol = original_volume * 0.5 / num_fu

        #get the possible symmetry vectors
        sites = sorted(self._sites, key=lambda site: site.species_string)
        grouped_sites = [list(a[1]) for a
                         in itertools.groupby(sites,
            key=lambda s: s.species_string)]
        min_site_list = min(grouped_sites, key=lambda group: len(group))

        min_site_list = [site.to_unit_cell for site in min_site_list]
        org = min_site_list[0].coords
        possible_vectors = [min_site_list[i].coords - org
                            for i in xrange(1, len(min_site_list))]

        #Let's try to use the shortest vector possible first. Allows for faster
        #convergence to primitive cell.
        possible_vectors = sorted(possible_vectors,
            key=lambda x: np.linalg.norm(x))

        # Pre-create a few varibles for faster lookup.
        all_coords = [site.coords for site in sites]
        all_sp = [site.species_and_occu for site in sites]
        new_structure = None
        for v, repl_pos in itertools.product(possible_vectors, xrange(3)):
            #Try combinations of new lattice vectors with existing lattice
            #vectors.
            latt = self._lattice.matrix
            latt[repl_pos] = v
            #Exclude coplanar lattices from consideration.
            if abs(np.linalg.det(latt)) > min_vol:
                latt = Lattice(latt)
                #Convert to fractional tol
                tol = [tolerance / l for l in latt.abc]
                new_frac = latt.get_fractional_coords(all_coords)
                grouped_sp = []
                grouped_frac = []

                for i, f in enumerate(new_frac):
                    found = False
                    for j, g in enumerate(grouped_frac):
                        if all_sp[i] != grouped_sp[j]:
                            continue
                        fdiff = np.abs(pbc_diff(g[0], f))
                        if np.all(fdiff < tol):
                            g.append(f)
                            found = True
                            break
                    if not found:
                        grouped_frac.append([f])
                        grouped_sp.append(all_sp[i])

                num_images = [len(c) for c in grouped_frac]
                nimages = num_images[0]
                if nimages > 1 and all([i == nimages for i in num_images]):
                    new_frac = [f[0] for f in grouped_frac]
                    new_structure = Structure(latt, grouped_sp, new_frac,
                        to_unit_cell=True)
                    break

        if new_structure and len(new_structure) != len(self):
            # If a more primitive structure has been found, try to find an
            # even more primitive structure again.
            return new_structure.get_primitive_structure(tolerance=tolerance)
        else:
            return self

    def __repr__(self):
        outs = ["Structure Summary", repr(self.lattice)]
        for s in self:
            outs.append(repr(s))
        return "\n".join(outs)

    def __str__(self):
        outs = ["Structure Summary ({s})".format(s=str(self.composition)),
                "Reduced Formula: {}"
                .format(self.composition.reduced_formula)]
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
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice": self._lattice.to_dict, "sites": []}
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
        outs = ["Molecule Summary"]
        for s in self:
            outs.append(repr(s))
        return "\n".join(outs)

    def __str__(self):
        outs = ["Molecule Summary ({s})".format(s=str(self.composition)),
                "Reduced Formula: " + self.composition.reduced_formula]
        to_s = lambda x: "%0.6f" % x
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i + 1), site.species_string,
                                  " ".join([to_s(j).rjust(12) for j in
                                            site.coords])]))
        return "\n".join(outs)

    @property
    def to_dict(self):
        """
        Json-serializable dict representation of Molecule
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "sites": [site.to_dict for site in self]}
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
        """
        Args:
            msg:
                The error message.
        """
        self.msg = msg

    def __str__(self):
        return self.msg
