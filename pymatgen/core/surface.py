# coding: utf-8

from __future__ import division, unicode_literals


"""
This module implements representations of slabs and surfaces, as well as
algorithms for generating them.
"""

__author__ = "Shyue Ping Ong, vivid0036, richardtran415"
__copyright__ = "Copyright 2014, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "6/10/14"

from fractions import gcd
import math
import itertools

from pymatgen import Structure

import numpy as np
from scipy.cluster.hierarchy import fclusterdata


from monty.fractions import lcm


def organize(struct):
    """Takes in a list of sites (ie. a structure object) and organizes the list
    based on the c coordinate from sites closest to c=0 to farthest. Then
    it returns this organized structure (list of sites) along with a list
    of organized coordinates along the c direction and an organized list
    of species corresponding the the site"""
    el = struct.species
    org_coords = struct.frac_coords.tolist()
    new_coord, b = [], []

    for i in struct.frac_coords:
        b.append(i[2])
        new_coord.append(b)
        b = []

    new_coord, org_coords, el = \
        zip(*sorted(zip(new_coord, org_coords, el)))

    return [Structure(struct.lattice, el, org_coords), new_coord, org_coords, el]


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
        return organize(self)[0]


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
        dist = float('inf')
        for i, j in enumerate(miller_index):
            if j == 0:
                # Lattice vector is perpendicular to surface normal, i.e.,
                # in plane of surface. We will simply choose this lattice
                # vector as one of the basis vectors.
                slab_scale_factor.append(eye[i])
            else:
                #Calculate projection of lattice vector onto surface normal.
                d = abs(np.dot(normal, latt.matrix[i]))
                non_orth_ind.append(i)
                if d < dist:
                    latt_index = i
                    dist = d

        if len(non_orth_ind) > 1:
            lcm_miller = lcm([miller_index[i] for i in non_orth_ind])
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

        self.oriented_unit_cell = single
        self.parent = initial_structure
        self.lll_reduce = lll_reduce
        self.standardize = standardize
        self.slab_scale_factor = slab_scale_factor
        self.normal = normal
        self.miller_index = miller_index
        self.min_vac_size = min_vacuum_size
        self.min_slab_size = min_slab_size
        self.dist = dist

    def get_slab(self, shift=0):
        """
        This method takes in a list of shift values created by
        the user and generates slabs based on the given shift values.

        Arg:
            shift_list (list of floats): Shift values in Angstrom
            determine how much a slab should be shifted.
        """
        slabs = []
        nlayers_slab = int(math.ceil(self.min_slab_size / self.dist))
        nlayers_vac = int(math.ceil(self.min_vac_size / self.dist))
        nlayers = nlayers_slab + nlayers_vac
        slab = self.oriented_unit_cell.copy()
        slab.make_supercell([1, 1, nlayers])

        new_sites = []
        for site in slab:
            if shift <= np.dot(site.coords, self.normal) < nlayers_slab * \
                    self.dist + shift:
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

    def get_all_slabs(self, thresh=0.00001, crit ='distance'):
        """
        A method that generates a list of shift values in order to create a list
        of slabs, each with a different surface. A surface is identified by
        measuring how close sites are to each other in order to be on the same
        surface. This is done using the scipy function, fclusterdata.

        Args:
            thresh (float): Threshold parameter in fclusterdata in order to determine
                the number of terminations to be found. Default thresh set to 0 to
                discern each atom with a unique c position as a unique termination.
            crit (str): The criterion to set for fclusterdata (see fcluster for
                description).
        """

        # Clusters sites together that belong in the same termination surface based
        # on their position in the c direction. How close the sites have to
        # be to belong in the same termination surface depends on the user input thresh.
        l = organize(self.oriented_unit_cell)
        tracker_index = fclusterdata(l[1], thresh, crit)
        new_coord, tracker_index, org_coords, el = \
            zip(*sorted(zip(l[1], tracker_index, l[2], l[3])))

        # Creates a list (term_index) that tells us which at
        # which site does a termination begin. For 1 unit cell.
        term_index = []
        for i in xrange(len(l[0])):
            if i == len(l[0]) - 1:
                term_index.append(i)
                break
            else:
                if tracker_index[i] != tracker_index[i+1]:
                    term_index.append(i)

        term_sites = []
        for i in xrange(len(term_index)):
            term_slab = self.oriented_unit_cell.copy()
            #term_slab.make_supercell(self.slab_scale_factor)
            term_slab = organize(term_slab)[0]
            term_sites.append(np.dot(term_slab[term_index[i]].coords, self.normal))
        print(len(term_sites))
        return [self.get_slab(shift) for shift in term_sites]

    def get_non_bond_breaking_slabs(self, specie1, specie2, max_bond=3,
                                    min_dist=0.01):
        """
        This method analysis whether we can generate a stable surface for a given miller
        index. Get the rotated unit cell from Class SurfaceGenerator().

            Args:
                specie1, specie2(str): elements of the bond. eg. "Li+"
                Search element2 from element1.
                max_bond(float): max length of bonds, In Angstroms. Default to 3.
                min_dist(float): min distance to be able to separate into two surfaces
        """

        # use SurfaceGenerate Class to get rotated unit cell, consider c factor of new cells
        # slab, the rotated unit cell
        slab = self.oriented_unit_cell
        latt = slab.lattice
        z2 = abs(latt.matrix[2][2])

        e1_fcoord = []
        e2_fcoord = []
        # get the coordinates lists of species1&2
        for site in slab:
            if site.species_string == specie1:
                e1_fcoord.append(site.frac_coords)
            else:
                if site.species_string == specie2:
                    e2_fcoord.append(site.frac_coords)
        print specie1,  len(e1_fcoord),\
            specie2, len(e2_fcoord)

        # the fcoord of atoms bonded are in a group, search specie2 from specie1
        # make projection of atoms(c factor)
        orig = []
        proj = []
        for i in xrange(len(e1_fcoord)):
            orig.append([e1_fcoord[i]])
            # make projection for specie1 atom
            e1_proj = e1_fcoord[i][2]
            proj.append([e1_proj])
            for j in xrange(len(e2_fcoord)):
                # periodic boundary condition, make images of species 2
                dist_img = latt.get_distance_and_image(e1_fcoord[i],
                                                       e2_fcoord[j])
                if dist_img[0] < max_bond:
                    # get the actual bonded species2(maybe a image)
                    e2_fc = e2_fcoord[j] + dist_img[1]
                    orig[i].append(e2_fc)
                    # make projection for specie2 atom
                    e2_proj = e2_fc[2]
                    # add fcoord of specie2 to the group of specie1
                    proj[-1].append(e2_proj)
        # sort projection of each group(bonded atoms listed by specie1) from low to high
        for group in proj:
            group.sort()
        # sort groups in proj according to the first item of each group(low boundary)
        proj.sort()

        for group in proj:
            print group[0], group[-1], len(group), group

        layer = []
        # each range of overlapped bonds is a group, layer[i][0]: bottom, layer[i][-1]: top
        layer.append([proj[0][0], proj[0][-1]])
        print self.miller_index, layer
        for i in xrange(len(proj)):
            # consider each group of proj,
            # if the min of group[i] is larger than the top of layer,
            # group[i] could be a new layer
            if layer[-1][-1] < proj[i][0]:
                layer.append([proj[i][0], proj[i][-1]])
                print i, '@-@'
            # layer[-1] and proj[i] are overlapped
            else:
                # change the top of layer[-1] to the top of proj[i], merge layer
                if layer[-1][-1] < proj[i][-1]:
                    layer[-1][-1] = proj[i][-1]
                    print i, '><'
        # consider the top from the *lowest* atoms, layer[0][0]
        top = layer[0][0] + 1
        if layer[-1][-1] < top-min_dist:
            layer.append([top-min_dist, top])
        else:
            layer[-1][-1] = top
        print len(layer), layer

        lyr = []
        for group in layer:
            group = list(np.array(group)*z2)
            lyr.append(group)
        stable_list = []
        if len(lyr) > 1:
            for i in xrange(len(lyr)-1):
                if lyr[i][1] + min_dist < lyr[i+1][0]:
                    stable_list.append(lyr[i][1]+0.5*min_dist)
            if len(stable_list):
                print 'stable surface:', self.miller_index

        return [self.get_slab(shift) for shift in stable_list]