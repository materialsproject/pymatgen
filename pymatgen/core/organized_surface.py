# Organization
# Info, authorship etc...
# Packages to import
# lcm() function
# organize() function

# Slab class (dummy class) return Slab (Structure)
    # surface_area() return np.linalg.norm(np.cross(m[0], m[1])) (float)
    # adsorb_atom() return structure_ad (Structure)
    # make_interface() return interface_system (Structure)

# surface_list_generator() function takes in a list of shifts to generate surfaces.
# Calls on Slab class to help make slabs then takes the slab that was spit out and
# further refines it with stuff like standardization, lll_reduction, etc...

# Surface_Generator class, initial parameters of a surface to be made
# 3 class methods
    # check_bounds() creates a list of shifts depending on stable surfaces, uses surface_list_generator() to make a list of surfaces
    # manual() creates a list of shifts manually selected by the user, uses surface_list_generator() to make a list of surfaces
    # distance() user specifies a criteria from fclusterdata to generate all surfaces
# return list_of_slabs (list of Slabs), list_of_terminations (a list of a list of sites on the surface of a slab)


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
    """Takes in a structure (ie. a list of sites) and organizes the list
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


    def __init__(self, lattice, species, normal, slab_scale_factor, new_sites, structure, miller_index):

        slab = structure.copy()

        self.parent = structure
        self.miller_index = miller_index
        self.normal = normal
        self.scale_factor = np.array(slab_scale_factor)


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
        # a123_s = lattice_s.matrix
        b123_s = lattice_s.inv_matrix.T
        # print surface
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

        for i in range(0, len(sec_struct)):
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

        for i in range(0, len(self)):
            interface_sites.append(self[i].coords)
            interface_species.append(specimen[i])

        specimen = sec_struct.species
        for i in range(0, len(sec_struct)):
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

def surface_list_generator(shift_list, normal, slab_scale, length, miller_index, structure, lll_reduce, standardize):

    slab_list = []

    for i in range(0, len(shift_list)):
        b = slab_scale
        term_slab = structure.copy()
        term_slab.make_supercell(b)
        term_slab = organize(term_slab)[0]

        new_sites = []

        # Get the sites of the surface atoms
        term_slab = organize(term_slab)[0]

        for site in term_slab:
            if shift_list[i] < np.dot(site.coords, normal) <= length +\
                    shift_list[i]:
                new_sites.append(site)

        term_slab = Structure.from_sites(new_sites)
        if lll_reduce:
            lll_slab = term_slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(term_slab.lattice)
            slab_scale_factor = np.dot(mapping[2], b)
            term_slab = lll_slab

        slab = Slab(term_slab.lattice, term_slab.species, normal, slab_scale_factor,
                    term_slab.sites, term_slab, miller_index)


        min_c = []
        index = []
        for ii in range(0, len(slab)):
            min_c.append(slab[ii].frac_coords[2])
            index.append(ii)
        min_length = min(min_c)
        max_length = max(min_c)
        slab.translate_sites(index, [0, 0, -min_length])

        if standardize:
            slab.translate_sites(index, [0, 0, (1-max_length+min_length)/2])

        slab_list.append(slab)

    return slab_list

class Surface_Generator():

    def __init__(self, structure, miller_index, min_slab_size,
                 min_vacuum_size, lll_reduce=True, standardize=True):

        latt = structure.lattice
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

        nlayers_slab = int(math.ceil(min_slab_size / dist))
        nlayers_vac = int(math.ceil(min_vacuum_size / dist))
        nlayers = nlayers_slab + nlayers_vac
        slab = structure.copy()
        single_slab = copy.deepcopy(slab_scale_factor)
        slab_scale_factor.append(eye[latt_index] * nlayers)
        single_slab.append(eye[latt_index]*1)
        single = slab.copy()
        single.make_supercell(single_slab)
        slab.make_supercell(slab_scale_factor)

        self.super = slab
        self.unit_cell = single
        self.struct = structure
        self.lll_reduce = lll_reduce
        self.standardize = standardize
        self.slab_scale_factor = slab_scale_factor
        self.normal = normal
        self.length = nlayers_slab*dist
        self.miller_index = miller_index

    def distance(self, thresh=0.00001, crit ='distance'):

# Clusters sites together that belong in the same termination surface based
# on their position in the c direction. Also organizes sites by ascending
# order from teh origin along the c direction. How close the sites have to
# be to belong in the same termination surface depends on the user input thresh.
        l = organize(self.super)
        tracker_index = fclusterdata(l[1], thresh, crit)
        new_coord, tracker_index, org_coords, el = \
            zip(*sorted(zip(l[1], tracker_index, l[2], l[3])))

# Creates a list (term_index) that tells us which at
# which site does a termination begin. For 1 unit cell.
        gg = 0
        term_index = []
        for i in range(0, len(l[0])):
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
        for i in range(0, len(term_index)):
            term_slab = self.super.copy()
            term_slab.make_supercell(self.slab_scale_factor)
            term_slab = organize(term_slab)[0]
            term_sites.append(np.dot(term_slab[term_index[i]].coords, self.normal))

        return surface_list_generator(term_sites, self.normal, self.slab_scale_factor, self.length,
                                      self.miller_index, self.struct, self.lll_reduce, self.standardize)


    def manual(self, shift_list):

        return surface_list_generator(shift_list, self.normal, self.slab_scale_factor, self.length,
                                      self.miller_index, self.struct, self.lll_reduce, self.standardize)


    def check_bonds(self, specie1, specie2, max_bond=3,
                    min_dist=0.01):
        """Args:
                structure (Structure): Initial input structure.
                miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
                specie1, specie2(str): elements of the bond. eg. "Li+"
                Search element2 from element1.
                max_bond(float): max length of bonds, In Angstroms. Default to 3.
                min_dist(float): min distance to be able to separate into two surfaces
        """

        # use Slab Class to get rotated unit cell, consider c factor of new cells
        # slab with 1 unit cell of layer and 0 unit cell of vacuum, rotated
        slab = self.unit_cell
        latt = slab.lattice
        print latt
        c2 = latt.c
        #c_new = latt.get_cartesian_coords((0, 0, 0.5)) / n1

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
        #print slab

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
                    #print e2_fc[2], dist_img[1]
                    orig[i].append(e2_fc)
                    # make projection for specie2 atom
                    e2_proj = e2_fc[2]
                    # add coords of specie2 to the group of specie1
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
        #consider the top from the *lowest* atoms, layer[0][0]
        top = layer[0][0] + 1
        if layer[-1][-1] < top-min_dist:
            layer.append([top-min_dist, top])
        else:
            layer[-1][-1] = top
        print len(layer), layer

        lyr = []
        for group in layer:
            group = list(np.array(group)*c2)
            lyr.append(group)
        # print c2, lyr
        #
        # from pymatgen import write_structure
        # slab1 = Slab(structure, miller_index, min_slab_size=1,
        #              min_vacuum_size=0)
        # write_structure(slab1, 'rucell'+str(miller_index)+str(slab.formula)+'.cif')

        stable_list = []
        if len(lyr) > 1:
            print 'stable surface'
            for i in range(0, len(lyr)):

                stable_list.append(lyr[i][1])
                for site in self.super:
                    if lyr[i][1] < np.dot(site.coords, self.normal) < lyr[i+1][1]:
                        stable_list.append(np.dot(site.coords, self.normal))
                stable_list.append(lyr[i+1][0])
                if i+1 == len(lyr)-1:
                    break

                #slab2 = Slab1(structure, miller_index, min_slab_size=1,
                #             min_vacuum_size=1, shift=lyr[0][0]-min_dist)
                #write_structure(slab2, 'rslab'+str(miller_index)+str(slab.formula)+'.cif')

        return surface_list_generator(stable_list, self.normal, self.slab_scale_factor, self.length,
                                      self.miller_index, self.struct, self.lll_reduce, self.standardize)


import unittest
from pymatgen.core.lattice import Lattice
from pymatgen.io.smartio import CifParser
from pymatgen import write_structure

# To run this test, it is assumed that the
# pymatgen folder is in your home directory
def get_path(path_str):
    file_name = "pymatgen/pymatgen/core/tests/surface_tests/" + path_str
    path = os.path.join(os.path.expanduser("~"), file_name)
    return path

class SlabTest(unittest.TestCase):

    def setUp(self):
        self.cu = Structure(Lattice.cubic(3), ["Cu", "Cu", "Cu", "Cu"],
                            [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5],
                             [0, 0.5, 0.5]])

        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])

        Li_Fe_P_O4 = CifParser(get_path("LiFePO4.cif"))
        self.lifepo4 = (Li_Fe_P_O4.get_structures(primitive = False)[0])


    def test_init(self):
        for hkl in itertools.product(xrange(4), xrange(4), xrange(4)):
            if any(hkl):
                ssize = 6
                vsize = 10
                gen = Surface_Generator(self.cu, hkl, ssize, vsize, standardize=False)
                s = gen.manual([0])[0]
                if hkl == [0, 1, 1]:
                    self.assertEqual(len(s), 13)
                    self.assertAlmostEqual(s.surface_area, 12.727922061357855)

                manual = self.cu.copy()
                manual.make_supercell(s.scale_factor)
                self.assertEqual(manual.lattice.lengths_and_angles,
                                 s.lattice.lengths_and_angles)

        # # For visual debugging
        # from pymatgen import write_structure
        # write_structure(s.parent, "cu.cif")
        # write_structure(s, "cu_slab_%s_%.3f_%.3f.cif" %
        #                  (str(hkl), ssize, vsize))

    def test_adsorb_atom(self):
        gen = Surface_Generator(self.cu,[0, 0, 1], 5, 5, standardize=False)
        s001 = gen.manual([0])[0]
        # print s001
        # O adsorb on 4Cu[0.5, 0.5, 0.25], abc = [3, 3, 12]
        # 1. test site_a = abc input
        s001_ad1 = Slab.adsorb_atom(structure_a=s001, site_a=[0.5, 0.5, 0.25], atom= ['O'],
                                    distance=2)
        self.assertEqual(len(s001_ad1), 9)
        for i in xrange(len(s001_ad1)):
            if str(s001_ad1[i].specie) == 'O':
                print s001_ad1[i].frac_coords
                self.assertAlmostEqual(s001_ad1[i].a, 0.5)
                self.assertAlmostEqual(s001_ad1[i].b, 0.5)
                self.assertAlmostEqual(s001_ad1[i].c, 0.4166667)
        self.assertEqual(s001_ad1.lattice.lengths_and_angles,
                                 s001.lattice.lengths_and_angles)
        # 2. test site_a = xyz input
        s001_ad2 = Slab.adsorb_atom(structure_a=s001, site_a=[1.5, 1.5, 3], atom= ['O'],
                                    distance=2, xyz=1)
        self.assertEqual(len(s001_ad2), 9)
        for i in xrange(len(s001_ad2)):
            if str(s001_ad2[i].specie) == 'O':
                print s001_ad2[i].frac_coords
                self.assertAlmostEqual(s001_ad2[i].a, 0.5)
                self.assertAlmostEqual(s001_ad2[i].b, 0.5)
                self.assertAlmostEqual(s001_ad2[i].c, 0.4166667)

    def test_make_terminations(self):
        # Need to make a more thorough test later
        hkl = [1, 0, 1]
        vsize = 10
        ssize = 15
        gen = Surface_Generator(self.lifepo4, hkl, ssize, vsize)
        LiFePO4 = gen.manual([0.478575453, 1.2545])[0]
        #
        # fileName = get_path("LiFePO4.cif")
        #
        #
        # Name = fileName[:-4] + "_%s_slab %s_vac %s_#" \
        #                           %(str(hkl),
        #                             ssize,
        #                             vsize)
        # fileType = ".cif"

        # For visual debugging
        # for ii in range(0, len(LiFePO4.slab_list)):
        #     name_num = str(ii)
        #     newFile = Name + name_num + fileType
        #     write_structure(LiFePO4.slab_list[ii], newFile)

        # Prints the coordinates and the species of each
        # atom along with the number of atoms on a
        # surface termination site and checks if the correct
        # number of sites are on a surface

        # for iii in range(0, len(LiFePO4.term_coords)):
        #     print(LiFePO4.term_coords[iii])
        #     print("%s atoms on this surface."
        #           %(len(LiFePO4.term_coords[iii])))
        # for i in range(0, 14):
        #     self.assertEqual(len(LiFePO4.term_coords[i]), 2)

    def test_make_interface(self):
        gen = Surface_Generator(self.lifepo4, [0, 0, 1], 10, 10)
        slab = gen.distance()[0]
        interface = slab.make_interface(self.libcc)

        # For visual debugging
        # write_structure(interface, "Li_LiFePO4_interface.cif")

        for i in range(0, len(slab)):
            self.assertEqual(slab[i].coords[0], interface[i].coords[0])
            self.assertEqual(slab[i].coords[1], interface[i].coords[1])
            self.assertEqual(slab[i].coords[2], interface[i].coords[2])

if __name__ == "__main__":
    unittest.main()