# !/usr/bin/env python

"""
This module implements representations of slabs and surfaces.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, vivid0036, richardtran415"
__copyright__ = "Copyright 2014, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "6/10/14"

from fractions import gcd
import math
import numpy as np
import itertools

from pymatgen.core.structure import Structure
from pymatgen import Lattice, Structure

import numpy as np
import copy
from scipy.cluster.hierarchy import fclusterdata
import os


def lcm(numbers):
    """Return lowest common multiple."""
    def lcm(a, b):
        return (a * b) / gcd(a, b)
    return reduce(lcm, numbers, 1)

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



    def __init__(self, normal, slab_scale_factor, structure, miller_index,
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

        slab = structure.copy()

        self.parent = initial_structure
        self.min_slab_size = min_slab_size
        self.min_vac_size = min_vac_size
        self.miller_index = miller_index
        self.scale_factor = np.array(slab_scale_factor)
        self.normal = normal


        super(Slab, self).__init__(
            slab.lattice, slab.species_and_occu,
            slab.frac_coords, site_properties=slab.site_properties)


    @property
    def surface_area(self):
        m = self.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))

    @classmethod
    def adsorb_atom(cls, structure_a, site_a, atom, distance,
                    surface=[0, 0, 1], xyz=0):
        """
        Gets the structure of single atom adsorption.

        Args:
        structure_a: the slab structure for adsorption
        site_a:  given sites for adsorption.
             default(xyz=0): site_a = [a, b, c], within [0,1];
             xyz=1: site_a = [x, y, z], in Angstroms.
        atom: adsorbed atom species
        distance: between centers of the adsorbed atom and the given site.
             in Angstroms
        surface: direction of the surface where atoms are adsorbed to
             default: surface = [0, 0, 1]
        xyz: default is 0. 1 means site_a = [x, y, z], in Angstroms.

        """
        from pymatgen.transformations.site_transformations import \
            InsertSitesTransformation

        lattice_s = structure_a.lattice
        abc_s = lattice_s.abc
        b123_s = lattice_s.inv_matrix.T
        vector_o = np.dot(surface, b123_s)
        print vector_o
        lens_v = np.sqrt(np.dot(vector_o, vector_o.T))
        V_o = vector_o / lens_v * distance

        if xyz == 0:
            # site_a = [a, b, c]
            for i in xrange(3):
                if site_a[i]> 1 or site_a[i] < 0:
                    raise ValueError("site_a is outsite the cell.")
            site_abc = V_o / abc_s + site_a
        else:
            # site_a = [x, y, z]
            for i in xrange(3):
                if site_a[i] > abc_s[i]:
                    raise ValueError("sites_a is outside the cell.")
            site_a1 = np.array(site_a)

            # convert to a1, a2, a3
            #< site_a2 = np.dot(a123_s, site_a1.T) / abc_s
            #< site_abc = (V_o+site_a2) / abc_s
            site_a2 = np.dot(b123_s, site_a1.T)
            site_abc = V_o/abc_s+site_a2

        for i in xrange(3):
            if site_abc[i] < 0 or site_abc[i] > 1:
                raise ValueError("wrong distance, atom will be outside the cell.")

        print 'add_site:', site_abc, atom

        ist = InsertSitesTransformation(species=atom, coords=[site_abc])
        structure_ad = ist.apply_transformation(structure_a)

        return structure_ad

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

def surface_list_generator(shift_list, normal, slab_scale, length, miller_index, initial_structure,
                           lll_reduce, standardize, min_slab_size, min_vac_size):
    """

    This function creates different slabs using the slab class depending
    on how many different shift values in shift_list there are and stores
    the slabs inside a list.

    Args:
        shift_list (list of floats): List of shift values in Angstroms
        normal (array): Surface normal vector.
        slab_scale (array): scale_factor Final computed scale factor
            that brings the parent cell to the surface cell.
        length (float): length of the slab region of the structure in Angstrom
        miller_index ([h, k, l]): Miller index of plane parallel to
            surface. Note that this is referenced to the input structure. If
            you need this to be based on the conventional cell,
            you should supply the conventional structure.
        structure (Structure): From the initial structure, this structure
            has been refined into a slab
        lll_reduce (bool): Whether to perform an LLL reduction on the
            eventual structure.
        standardize (bool): Whether to center the slab in the cell with equal
            vacuum spacing from the top and bottom.
        initial_structure (Structure): Initial input structure.
        min_slab_size (float): In Angstroms
        min_vac_size (float): In Angstroms
        initial_structure (Structure): The initial structure used to create a slab

    """

    slab_list = []

    for i in range(0, len(shift_list)):
        b = slab_scale
        term_slab = initial_structure.copy()
        term_slab.make_supercell(b)
        # term_slab = organize(term_slab)[0]

        new_sites = []
        # Get the sites of the surface atoms
        # length is in the direction of surface normal(z), so is shift
        for site in term_slab:
            if shift_list[i] <= np.dot(site.coords, normal) < length +\
                    shift_list[i]:
                new_sites.append(site)

        #if len(new_sites) == 0:
        #    continue
            
        term_slab = Structure.from_sites(new_sites)

        if lll_reduce:
            lll_slab = term_slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(term_slab.lattice)
            slab_scale_factor = np.dot(mapping[2], b)
            term_slab = lll_slab
        else:
            slab_scale_factor = b

        slab = Slab(normal, slab_scale_factor, term_slab, miller_index,
                    initial_structure, min_slab_size, min_vac_size)

        if standardize:
            frac_c_list = []
            for site in xrange(len(slab)):
                frac_c_list.append(slab.frac_coords[site][2])
            frac_c_list.sort()
            # move along c, distance = frac_difference between cell & layer center
            slab.translate_sites(range(0, len(frac_c_list)),
                                 [0, 0, 0.5-(frac_c_list[0]+frac_c_list[-1])/2])

        slab_list.append(slab)

    return slab_list


class SurfaceGenerator():

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
                 min_vacuum_size, lll_reduce=True, standardize=True):
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

        # By generating the normal first, we can prevent
        # the loop that creates the slabs from iterating
        # through repetitive steps.
        nlayers_slab = int(math.ceil(min_slab_size / dist))
        nlayers_vac = int(math.ceil(min_vacuum_size / dist))
        nlayers = nlayers_slab + nlayers_vac
        slab = initial_structure.copy()
        single_slab = copy.deepcopy(slab_scale_factor)
        slab_scale_factor.append(eye[latt_index] * nlayers)
        single_slab.append(eye[latt_index]*1)
        single = slab.copy()
        single.make_supercell(single_slab)
        slab.make_supercell(slab_scale_factor)

        self.super = slab
        self.unit_cell = single
        self.parent = initial_structure
        self.lll_reduce = lll_reduce
        self.standardize = standardize
        self.slab_scale_factor = slab_scale_factor
        self.normal = normal
        self.length = nlayers_slab*dist
        self.miller_index = miller_index
        self.min_vac_size = min_vacuum_size
        self.min_slab_size = min_slab_size

    def get_all_slab(self, thresh=0.00001, crit ='distance'):

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
        l = organize(self.super)
        tracker_index = fclusterdata(l[1], thresh, crit)
        new_coord, tracker_index, org_coords, el = \
            zip(*sorted(zip(l[1], tracker_index, l[2], l[3])))

        # Creates a list (term_index) that tells us which at
        # which site does a termination begin. For 1 unit cell.
        gg = 0
        term_index = []
        for i in xrange(len(l[0])):
            gg+=1
            if i == len(l[0]) - 1:
                term_index.append(i)
                break
            else:
                if tracker_index[i] != tracker_index[i+1]:
                    term_index.append(i)
            if gg >= len(self.unit_cell) and tracker_index[i] != tracker_index[i+1]:
                term_index.append(i)
                break

        term_sites = []
        for i in xrange(len(term_index)):
            term_slab = self.super.copy()
            #term_slab.make_supercell(self.slab_scale_factor)
            term_slab = organize(term_slab)[0]
            term_sites.append(np.dot(term_slab[term_index[i]].coords, self.normal))

        return surface_list_generator(term_sites, self.normal, self.slab_scale_factor, self.length,
                                      self.miller_index, self.parent, self.lll_reduce, self.standardize,
                                      self.min_slab_size, self.min_vac_size)


    @classmethod
    def get_slab(self, shift_list):

        """

        This method takes in a list of shift values created by
        the user and generates slabs based on the given shift values.

        Arg:
            shift_list (list of floats): Shift values in Angstrom
            determine how much a slab should be shifted.

        """

        return surface_list_generator(shift_list, self.normal, self.slab_scale_factor, self.length,
                                      self.miller_index, self.parent, self.lll_reduce, self.standardize,
                                      self.min_slab_size, self.min_vac_size)


    @classmethod
    def get_non_bond_breaking_slab(self, specie1, specie2, max_bond=3,
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
        slab = self.unit_cell
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

        return surface_list_generator(stable_list, self.normal, self.slab_scale_factor, self.length,
                                      self.miller_index, self.parent, self.lll_reduce, self.standardize,
                                      self.min_slab_size, self.min_vac_size)

