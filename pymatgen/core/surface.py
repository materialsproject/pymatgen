# !/usr/bin/env python

"""
This module implements representations of slabs and surfaces.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
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
    """
    Subclass of Structure representing a Slab. Implements additional
    attributes pertaining to slabs.

    .. attribute:: parent

        Parent structure from which Slab was derived.

    .. attribute:: min_slab_size

        Minimum size in angstroms of layers containing atoms

    .. attribute:: min_vac_size

        Minimize size in angstroms of layers containing vacuum

    .. attribute:: scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: normal

        Surface normal vector.

    .. attribute:: term_coords

        List of list of sites on a surface

    .. attribute:: slab_list

        List of structures with all possible different surfaces

    .. attribute:: unit_cell
        The unit cell of the structure relative to the transformation based
        on the given miller index
    """

    def __init__(self, structure, miller_index, min_slab_size,
                 min_vacuum_size, thresh=0.00001, crit ='distance',
                 lll_reduce=True, standardize=True, nth_term=0,
                 shift=0, shift_bool=False):
        """
        Makes a Slab structure. Note that the code will make a slab with
        whatever size that is specified, rounded upwards. The a and b lattice
        vectors are always in-plane, and the c lattice is always out of plane
        (though not necessarily orthogonal). Enumerates through the different
        terminations in a unit cell in a supercell. The user can retrieve any
        of the supercells after creating the slab by calling slab_list. For
        example, a = Slab(self, structure, miller_index, min_slab_size,
        min_vacuum_size, thresh). a.slab_list[4] returns the structure
        with the fourth termination as the surface of the slab at hte bottom.
        Similarly, a.term_coords[4] will return the sites of the atoms on the
        surface.

        Args:
            structure (Structure): Initial input structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            min_slab_size (float): In Angstroms
            min_vacuum_size (float): In Angstroms
            thresh (float): Threshold parameter in fclusterdata in order to determine
                the number of terminations to be found. Default thresh set to 0 to
                discern each atom with a unique c position as a unique termination.
            crit (str): The criterion to set for fclusterdata (see fcluster for
                description).
            lll_reduce (bool): Whether to perform an LLL reduction on the
                eventual structure.
            standardize (bool): Whether to center the slab portion of the structure in
                between the vacuum layer.
            nth_term (int): The selected structure from slab_list to perform the class
                methods on
            shift_bool (bool): Turns off the termination enumerating algorithm and
                allows the user to manually shift a slab
            shift (float): Determines how much to shift the slab if shift_bool is True

        """
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


    # Clusters sites together that belong in the same termination surface based
    # on their position in the c direction. Also organizes sites by ascending
    # order from teh origin along the c direction. How close the sites have to
    # be to belong in the same termination surface depends on the user input thresh.
        l = organize(slab)
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
            if gg >= len(single) and tracker_index[i] != tracker_index[i+1]:
                term_index.append(i)
                break

        slab_list, term_coords = [], []
        t, a = 0, 0

        for i in range(0, len(term_index)):
            b = slab_scale_factor
            y = []
            term_slab = structure.copy()
            term_slab.make_supercell(b)
            term_slab = organize(term_slab)[0]
            term_site = np.dot(term_slab[term_index[i]].coords, normal)
            new_sites = []

            # Get the sites of the surface atoms
            term_slab = organize(term_slab)[0]
            if i == len(term_index)-1:
                a = len(structure)
            elif len(term_slab) < a:
                a = len(term_slab)-1
            else:
                a = term_index[i]+1
            for iv in range(t, a):
                y.append(term_slab[iv])
            term_coords.append(y)
            if a == len(term_slab)-1:
                break
            t = a

            if shift_bool:
                term_site = shift

            for site in term_slab:
                if term_site < np.dot(site.coords, normal) <= nlayers_slab * dist +\
                        term_site:
                    new_sites.append(site)

            term_slab = Structure.from_sites(new_sites)
            if lll_reduce:
                lll_slab = term_slab.copy(sanitize=True)
                if i == len(term_index)-1:
                    mapping = lll_slab.lattice.find_mapping(term_slab.lattice)
                    slab_scale_factor = np.dot(mapping[2], b)
                term_slab = lll_slab

            min_c = []
            index = []
            for ii in range(0, len(term_slab)):
                min_c.append(term_slab[ii].frac_coords[2])
                index.append(ii)
            min_length = min(min_c)
            max_length = max(min_c)
            term_slab.translate_sites(index, [0, 0, -min_length])

            if standardize:
                term_slab.translate_sites(index, [0, 0, (1-max_length+min_length)/2])

            slab_list.append(term_slab)

            if shift_bool:
                break

        if 0> nth_term >=len(slab_list):
            raise IndexError("nth_term is out of range, set nth_term "
                             "between 0 and %s" %(len(slab_list)))

        slab = slab_list[nth_term].copy()

        self.min_slab_size = min_slab_size
        self.nlayers_slab = nlayers_slab
        self.min_vac_size = min_vacuum_size
        self.slab_list = slab_list
        self.parent = structure
        self.miller_index = miller_index
        self.term_coords = term_coords
        self.scale_factor = np.array(slab_scale_factor)
        self.normal = normal
        self.true_vac_size = nlayers_vac*dist
        self.unit_cell = single.copy()

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
                s = Slab(self.cu, hkl, ssize, vsize, shift_bool=True, standardize=False, lll_reduce=False)
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
        s001 = Slab(self.cu,[0, 0, 1], 5, 5)
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
        hkl = [0, 1, 1]
        vsize = 10
        ssize = 15
        LiFePO4 = Slab(self.lifepo4, hkl, ssize, vsize)

        fileName = get_path("LiFePO4.cif")


        Name = fileName[:-4] + "_%s_slab %s_vac %s_#" \
                                  %(str(hkl),
                                    ssize,
                                    vsize)
        fileType = ".cif"

        # For visual debugging
        # for ii in range(0, len(LiFePO4.slab_list)):
        #     name_num = str(ii)
        #     newFile = Name + name_num + fileType
        #     write_structure(LiFePO4.slab_list[ii], newFile)

        # Prints the coordinates and the species of each
        # atom along with the number of atoms on a
        # surface termination site and checks if the correct
        # number of sites are on a surface
        for iii in range(0, len(LiFePO4.term_coords)):
            print(LiFePO4.term_coords[iii])
            print("%s atoms on this surface."
                  %(len(LiFePO4.term_coords[iii])))
        for i in range(0, 14):
            self.assertEqual(len(LiFePO4.term_coords[i]), 2)

    def test_make_interface(self):
        slab = Slab(self.lifepo4, [0, 0, 1], 10, 10)
        interface = slab.make_interface(self.libcc)

        # For visual debugging
        # write_structure(interface, "Li_LiFePO4_interface.cif")

        for i in range(0, len(slab)):
            self.assertEqual(slab[i].coords[0], interface[i].coords[0])
            self.assertEqual(slab[i].coords[1], interface[i].coords[1])
            self.assertEqual(slab[i].coords[2], interface[i].coords[2])

if __name__ == "__main__":
    unittest.main()
