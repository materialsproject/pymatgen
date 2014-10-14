# coding: utf-8

from __future__ import division, unicode_literals

"""
This module implements representations of slabs and surfaces, as well as
algorithms for generating them.
"""

__author__ = "Richard Tran, Zihan Xu, Shyue Ping Ong"
__copyright__ = "Copyright 2014, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "6/10/14"

from fractions import gcd
import math
import itertools
import logging

import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

from monty.fractions import lcm

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord_utils import in_coord_list


logger = logging.getLogger(__name__)


class Slab(Structure):
    """
    Subclass of Structure representing a Slab. Implements additional
    attributes pertaining to slabs, but the init method does not
    actually implement any algorithm that creates a slab. This is a
    DUMMY class that's only purpose is to hold information about a slab.

    .. attribute:: parent

        Parent structure from which Slab was derived.

    .. attribute:: min_slab_size

        Minimum size in angstroms of layers containing atoms

    .. attribute:: min_vac_size

        Minimize size in angstroms of layers containing vacuum

    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: normal

        Surface normal vector.
    """

    def __init__(self, structure, miller_index, shift, normal, scale_factor,
                 initial_structure, min_slab_size, min_vac_size):
        """
        Makes a Slab structure, a structure object
        with additional information pertaining to slabs.

        Args:
            normal (array): Surface normal vector.
            slab_scale_factor (array): scale_factor Final computed scale factor
                that brings the parent cell to thesurface cell.
            structure (Structure): From the initial structure, this structure
                has been refined into a slab
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            initial_structure (Structure): Initial input structure.
            min_slab_size (float): In Angstroms
            min_vac_size (float): In Angstroms
        """
        slab = structure
        self.parent = initial_structure
        self.min_slab_size = min_slab_size
        self.min_vac_size = min_vac_size
        self.miller_index = miller_index
        self.scale_factor = np.array(scale_factor)
        self.normal = normal
        self.shift = shift

        super(Slab, self).__init__(
            slab.lattice, slab.species_and_occu,
            slab.frac_coords, site_properties=slab.site_properties)

    @property
    def surface_area(self):
        m = self.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))

    def add_adsorbate_atom(self, indices, specie, distance):
        """
        Gets the structure of single atom adsorption.
        slab structure from the Slab class(in [0, 0, 1])

        Args:
            indices ([int]): Indices of sites on which to put the absorbate.
                Absorbed atom will be displaced relative to the center of
                these sites.
            specie (Specie/Element/str): adsorbed atom species
            distance (float): between centers of the adsorbed atom and the
                given site in Angstroms.
        """
        #Let's do the work in cartesian coords
        center = np.sum([self[i].coords for i in indices], axis=0) / len(
            indices)

        coords = center + self.normal * distance / np.linalg.norm(self.normal)

        self.append(specie, coords, coords_are_cartesian=True)

    def make_interface(self, sec_struct):

        # Creates a supercell of the secondary structure
        m = self.lattice.matrix
        A1 = np.linalg.norm(np.cross(m[0], m[1]))
        m = sec_struct.lattice.matrix
        A2 = np.linalg.norm(np.cross(m[0], m[1]))

        height = int(math.ceil((A1/self.lattice.a)
                               / (A2/sec_struct.lattice.a))*2)
        base = int(math.ceil(self.lattice.a/sec_struct.lattice.a)*2)
        c_basis = int(math.ceil(self.lattice.c/sec_struct.lattice.c)*2)
        interface_scale_factor = [base, height, c_basis]
        sec_struct.make_supercell(interface_scale_factor)

        # Carves the secondary supercell into a shape compatible for interfacing with the slab
        new_sites = []
        specimen = sec_struct.species
        scaled_species = []

        for i in xrange(len(sec_struct)):
            if 0 <= sec_struct[i].coords[0] <= self.lattice.a and \
               0 <= sec_struct[i].coords[1] <= self.lattice.b and \
               0 <= sec_struct[i].coords[2] <= self.lattice.c:
                new_sites.append(sec_struct[i])
                scaled_species.append(specimen[i])

        scaled_sites = []
        for site in new_sites:
            scaled_sites.append(site.coords)

        sec_struct = Structure(self.lattice, scaled_species,
                               scaled_sites, coords_are_cartesian = True)

        # Creates the interface between two structures
        interface_sites = []
        interface_species = []
        # c_list = []
        specimen = self.species

        for i in xrange(len(self)):
            interface_sites.append(self[i].coords)
            interface_species.append(specimen[i])

        specimen = sec_struct.species
        for i in xrange(len(sec_struct)):
            if self.true_vac_size/2 <= sec_struct[i].coords[2] \
                            <= self.lattice.c-self.true_vac_size/2:
                continue
            else:
                interface_sites.append(sec_struct[i].coords)
                interface_species.append(specimen[i])

        # For some reason self[i].coords[ii] has coordinates
        # a and b switched, figure out the source later
        for i in range(len(self)+1, len(interface_sites)):
            interface_sites[i][0], interface_sites[i][1] = \
                interface_sites[i][1], interface_sites[i][0]

        interface_system = Structure(self.lattice, interface_species,
                                     interface_sites,
                                     coords_are_cartesian=True)

        return interface_system

    def organize_along_c(self):
        # Organizes a slab's sites along the c direction
        self.sites = sorted(self.sites, key=lambda s: s.c)


class SurfaceGenerator(object):

    """
    Generates a supercell of the initial structure with its miller
    index of plane parallel to surface. Once the supercell is
    generated, three methods can be used to create a list of Slab
    objects with the supercell.

    .. attribute:: super

        A supercell of the parent structure with the miller
        index of plane parallel to surface

    .. attribute:: unit_cell

        A unit cell of the parent structure with the miller
        index of plane parallel to surface

    .. attribute:: parent

        Parent structure from which Slab was derived.

    .. attribute:: lll_reduce

        Whether or not the slabs will be orthogonalized

    .. attribute:: standardize

        Whether or not the slabs will be centered between
        the vacuum layer

    .. attribute:: scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: normal

        Surface normal vector.

    .. attribute:: length

        Length of the slab layer of the surface(along normal)

    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: min_slab_size

        Minimum size in angstroms of layers containing atoms

    .. attribute:: min_vac_size

        Minimize size in angstroms of layers containing vacuum

    """

    def __init__(self, initial_structure, miller_index, min_slab_size,
                 min_vacuum_size, lll_reduce=False, standardize=False):
        """
        Generates a supercell of the initial structure with its miller
        index of plane parallel to surface.

        Args:
            initial_structure (Structure): Initial input structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            min_slab_size (float): In Angstroms
            min_vac_size (float): In Angstroms
            lll_reduce (bool): Whether to perform an LLL reduction on the
                eventual structure.
            standardize (bool): Whether to center the slab in the cell with equal
                vacuum spacing from the top and bottom.

        """
        latt = initial_structure.lattice
        d = reduce(gcd, miller_index)
        miller_index = [int(i / d) for i in miller_index]
        #Calculate the surface normal using the reciprocal lattice vector.
        recp = latt.reciprocal_lattice_crystallographic
        normal = recp.get_cartesian_coords(miller_index)
        normal /= np.linalg.norm(normal)

        slab_scale_factor = []
        non_orth_ind = []
        eye = np.eye(3, dtype=np.int)
        dist = 0
        for i, j in enumerate(miller_index):
            if j == 0:
                # Lattice vector is perpendicular to surface normal, i.e.,
                # in plane of surface. We will simply choose this lattice
                # vector as one of the basis vectors.
                slab_scale_factor.append(eye[i])
            else:
                #Calculate projection of lattice vector onto surface normal.
                d = abs(np.dot(normal, latt.matrix[i])) / np.linalg.norm(
                    latt.matrix[i])
                non_orth_ind.append(i)
                if d > dist:
                    # We want the vector that has maximum magnitude in the
                    # direction of the surface normal as the c-direction.
                    # Results in a more "orthogonal" unit cell.
                    latt_index = i
                    dist = d

        if len(non_orth_ind) > 1:
            lcm_miller = lcm(*[miller_index[i] for i in non_orth_ind])
            for i, j in itertools.combinations(non_orth_ind, 2):
                l = [0, 0, 0]
                l[i] = -int(round(lcm_miller / miller_index[i]))
                l[j] = int(round(lcm_miller / miller_index[j]))
                slab_scale_factor.append(l)
                if len(slab_scale_factor) == 2:
                    break

        slab_scale_factor.append(eye[latt_index])
        single = initial_structure.copy()
        single.make_supercell(slab_scale_factor)

        self.oriented_unit_cell = Structure.from_sites(single,
                                                       to_unit_cell=True)
        self.parent = initial_structure
        self.lll_reduce = lll_reduce
        self.standardize = standardize
        self.slab_scale_factor = slab_scale_factor
        self.normal = normal
        self.miller_index = miller_index
        self.min_vac_size = min_vacuum_size
        self.min_slab_size = min_slab_size

    def get_slab(self, shift=0):
        """
        This method takes in shift value for the c lattice direction and
        generates a slab based on the given shift. You should rarely use this
        method. Instead, it is used by other generation algorithms to obtain
        all slabs.

        Arg:
            shift (float): A shift value in Angstrom that determines how much a
            slab should be shifted.

        Returns:
            (Slab) A Slab object with a particular shifted oriented unit cell.
        """
        dist = abs(np.dot(self.normal,
                          self.oriented_unit_cell.lattice.matrix[2]))

        nlayers_slab = int(math.ceil(self.min_slab_size / dist))
        nlayers_vac = int(math.ceil(self.min_vac_size / dist))
        nlayers = nlayers_slab + nlayers_vac

        slab = self.oriented_unit_cell.copy()
        slab.translate_sites(range(len(slab)), [0, 0, -shift])
        slab = Structure.from_sites(slab, to_unit_cell=True)

        slab.make_supercell([1, 1, nlayers])

        new_sites = []
        for site in slab:
            if 0 <= site.c < nlayers_slab / nlayers:
                new_sites.append(site)

        slab = Structure.from_sites(new_sites)
        scale_factor = self.slab_scale_factor
        if self.lll_reduce:
            lll_slab = slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(slab.lattice)
            scale_factor = np.dot(mapping[2], scale_factor)
            slab = lll_slab

        if self.standardize:
            frac_c_list = []
            for site in xrange(len(slab)):
                frac_c_list.append(slab.frac_coords[site][2])
            frac_c_list.sort()
            # move along c, distance = frac_difference between cell & layer center
            slab.translate_sites(range(0, len(frac_c_list)),
                                 [0, 0, 0.5 - (
                                     frac_c_list[0] + frac_c_list[-1]) / 2])

        return Slab(slab, self.miller_index, shift, self.normal,
                    scale_factor, self.parent, self.min_slab_size,
                    self.min_vac_size)

    def _calculate_possible_shifts(self, tol=1e-2):
        frac_coords = self.oriented_unit_cell.frac_coords
        n = len(frac_coords)

        # We cluster the sites according to the c coordinates. But we need to
        # take into account PBC. Let's compute a fractional c-coordinate
        # distance matrix that accounts for PBC.
        dist_matrix = np.zeros((n, n))
        for i in range(n - 1):
            for j in range(i + 1, n):
                cdist = frac_coords[i][2] - frac_coords[j][2]
                cdist = abs(cdist - round(cdist))
                dist_matrix[i, j] = cdist
                dist_matrix[j, i] = cdist

        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        clusters = fcluster(z, tol, criterion="distance")

        #Generate dict of cluster# to c val - doesn't matter what the c is.
        c_loc = {c: frac_coords[i][2] for i, c in enumerate(clusters)}

        #Put all c into the unit cell.
        possible_c = [c - math.floor(c) for c in sorted(c_loc.values())]

        # Calculate the shifts
        nshifts = len(possible_c)
        shifts = []
        for i in range(nshifts):
            if i == nshifts - 1:
                # There is an additional shift between the first and last c
                # coordinate. But this needs special handling because of PBC.
                shift = (possible_c[0] + 1 + possible_c[i]) * 0.5
                if shift > 1:
                    shift -= 1
            else:
                shift = (possible_c[i] + possible_c[i + 1]) * 0.5
            shifts.append(shift)
        shifts = sorted(shifts)
        return shifts

    def get_slabs(self, bonds=None, tol=1e-2):
        """
        A method that generates a list of shift values in order to create a list
        of slabs, each with a different surface. A surface is identified by
        measuring how close sites are to each other in order to be on the same
        surface. This is done using the scipy function, fclusterdata.

        Args:
            thresh (float): Threshold parameter in fcluster in order to check
                if two atoms are lying on the same plane. Default thresh set
                to 1e-2 in the c fractional coordinate.
            bonds ({(specie1, specie2): max_bond_dist}: bonds are
                specified as a dict of tuples: float of specie1, specie2
                and the max bonding distance. For example, PO4 groups may be
                defined as {("P", "O"): 3}.

        Returns:
            ([Slab]) List of all possible terminations of a particular surface.
        """
        forbidden_c_ranges = []
        if bonds is not None:
            #Convert to species first
            bonds = {(get_el_sp(s1), get_el_sp(s2)): dist for (s1, s2), dist in
                     bonds.items()}
            for site1, site2 in itertools.combinations(self.oriented_unit_cell, 2):
                all_sp = set(site1.species_and_occu.keys())
                all_sp.update(site2.species_and_occu.keys())
                for species, bond_dist in bonds.items():
                    if all_sp.issuperset(species):
                        dist, image = site1.distance_and_image(site2)
                        if dist < bond_dist:
                            min_c = site1.c
                            max_c = (site2.frac_coords + image)[2]
                            c_range = sorted([min_c, max_c])
                            if c_range[1] > 1:
                                forbidden_c_ranges.append((c_range[0], 1))
                                forbidden_c_ranges.append((0, c_range[1] -1))
                            elif c_range[0] < 0:
                                forbidden_c_ranges.append((0, c_range[1]))
                                forbidden_c_ranges.append((c_range[0] + 1, 1))
                            else:
                                forbidden_c_ranges.append(c_range)

        def shift_allowed(shift):
            for r in forbidden_c_ranges:
                if r[0] <= shift <= r[1]:
                    return False
            return True

        return [self.get_slab(shift)
                for shift in self._calculate_possible_shifts(tol=tol)
                if shift_allowed(shift)]


def generate_all_slabs(structure, max_index, min_slab_size, min_vacuum_size,
                       bonds=None, tol=1e-2, symprec=0.01):
    analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
    symm_ops = analyzer.get_symmetry_operations()
    processed = []
    def is_already_analyzed(miller_index):
        for op in symm_ops:
            if in_coord_list(processed, op.operate(miller_index)):
                return True
        return False

    r = range(-max_index, max_index + 1)
    all_slabs = []
    for miller in itertools.product(r, r, r):
        if any([i != 0 for i in miller]):
            d = reduce(gcd, miller)
            miller = tuple([int(i / d) for i in miller])
            if not is_already_analyzed(miller):
                gen = SurfaceGenerator(structure, miller, min_slab_size,
                                       min_vacuum_size)
                slabs = gen.get_slabs(bonds=bonds, tol=tol)
                if len(slabs) > 0:
                    logger.debug("%s has %d slabs... " % (miller, len(slabs)))
                    all_slabs.extend(slabs)
                processed.append(miller)
    return all_slabs
