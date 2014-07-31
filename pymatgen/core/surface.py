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
import scipy.cluster.hierarchy




def lcm(numbers):
    """Return lowest common multiple."""
    def lcm(a, b):
        return (a * b) / gcd(a, b)
    return reduce(lcm, numbers, 1)


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
    """

    def __init__(self, structure, miller_index, min_slab_size,
                 min_vacuum_size, temp_thresh, lll_reduce=True, shift=0):
        """
        Makes a Slab structure. Note that the code will make a slab with
        whatever size that is specified, rounded upwards. The a and b lattice
        vectors are always in-plane, and the c lattice is always out of plane
        (though not necessarily orthogonal).

        Args:
            structure (Structure): Initial input structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            min_slab_size (float): In Angstroms
            min_vacuum_size (float): In Angstroms
            temp_thresh (float): A temporary parameter used for setting the
                t parameter in fclusterdata in order to determine the number of
                terminations to be found. An algorithm wil eventually be implemented
                so that the program can determine this value by itself depending on
                the nature of the structure object parameter.
            lll_reduce (bool): Whether to perform an LLL reduction on the
                eventual structure.
            shift (float): In Angstroms (shifting the origin)
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

        min_slab_size = float(min_slab_size)
        min_vacuum_size = float(min_vacuum_size)
        nlayers_slab = int(math.ceil(min_slab_size / dist))
        nlayers_vac = int(math.ceil(min_vacuum_size / dist))
        nlayers = nlayers_slab + nlayers_vac

        slab_scale_factor.append(eye[latt_index] * nlayers)
        slab = structure.copy()
        slab.make_supercell(slab_scale_factor)
        new_sites = []

        # This for loop needs to take into account fractional atoms in a unit cell when c = 0 or c = 1 .
        # The technique to do so is highly volatile, make sure to analyze it later if something goes wrong. Such as when
        # a termination is missing from one of the faces of the cell.
        for site in slab:
            if shift <= np.dot(site.coords, normal) <= min_slab_size + shift:
                new_sites.append(site)
            if shift == 0 and np.dot(site.coords, normal) <= 0:
                new_sites.append(site)
            if shift == 0 and (site.coords[2] == slab.lattice.c or site.coords[1] == slab.lattice.c or site.coords[0] == slab.lattice.c):
                new_sites.append(site)

        slab = Structure.from_sites(new_sites)

# After this line is where all the necessary algorithms need to enumerate the different terminations is inserted
        el = slab.species
        org_coords = slab.frac_coords.tolist()
        new_coord, b = [], []

        for i in slab.frac_coords:
            b.append(i[2])
            new_coord.append(b)
            b = []

        tracker_index = scipy.cluster.hierarchy.fclusterdata(new_coord, temp_thresh, criterion= 'distance')
        new_coord, tracker_index, org_coords, el = zip(*sorted(zip(new_coord, tracker_index, org_coords, el)))

# Creates a list that keeps track of repetitive values of c in different coordinates
        org_index = []
        for i in range(0, len(slab)):
            if i == len(slab) - 1:
                org_index.append(i)
            else:
                if tracker_index[i] != tracker_index[i+1]:
                    org_index.append(i)

        term_enum, term_coords = [], []
        a, i = 0, 0
        new_slab = Structure(slab.lattice, el, org_coords)

        for iii in org_index:
            y = []
            alt_slab = new_slab.copy()
            n = 0

            for ii in range(0, len(alt_slab)):
                index = []
                term_scale = 0
                index.append(ii)
                if alt_slab.frac_coords[ii][2] > alt_slab.frac_coords[iii][2]:
                    term_scale = min_vacuum_size/(dist*nlayers)
                if float(alt_slab.frac_coords[ii][2]) + term_scale > 1: # Out of bounds atom
                    print("Warning, the shift in the origin may cause some sites be pushed outside the slab")
                    continue
                alt_slab.translate_sites(index, [0, 0, term_scale])
                n += 1

            term_enum.append(alt_slab)

            for iv in range(a, iii + 1):
                    y.append(term_enum[i][iv])
            i += 1
            term_coords.append(y)
            a = iii+1

        if lll_reduce:
            lll_slab = slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(slab.lattice)
            slab_scale_factor = np.dot(mapping[2], slab_scale_factor)
            slab = lll_slab


        self.min_slab_size = min_slab_size
        self.min_vac_size = min_vacuum_size
        self.enum = term_enum # Holds a list of Structure objects of slabs with different terminations
        self.parent = structure
        self.miller_index = miller_index
        self.org_in = org_index
        self.term_coords = term_coords # Holds the corresponding list of sites on the surface terminations
        self.temp_thresh = temp_thresh

        super(Slab, self).__init__(
            slab.lattice, slab.species_and_occu, slab.frac_coords,
            site_properties=slab.site_properties)

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

from pymatgen.io.smartio import CifParser
from pymatgen import write_structure
import unittest

class SlabTest(unittest.TestCase): # 4 different unit cells will be tested, ZnO, BaTiO4, TeI, and LifePO4

    def setUp(self):

        Zn_O = CifParser("pymatgen/pymatgen/core/surface tests/001_terminations/ZnO-wz.cif")
        self.zno = (Zn_O.get_structures(primitive = False)[0])
        Te_I  = CifParser("pymatgen/pymatgen/core/surface tests/001_terminations/icsd_TeI.cif")
        self.tei= (Te_I.get_structures(primitive = False)[0])
        Ba_Ti_O3  = CifParser("pymatgen/pymatgen/core/surface tests/001_terminations/icsd_batio3.cif")
        self.batio3= (Ba_Ti_O3.get_structures(primitive = False)[0])
        Li_Fe_PO4  = CifParser("pymatgen/pymatgen/core/surface tests/001_terminations/LiFePO4.cif")
        self.lifepo4= (Li_Fe_PO4.get_structures(primitive = False)[0])

    def test_init(self):

        hkl001 = [0, 0, 1]
        z001 = Slab(self.zno, hkl001, 10, 3, 0.025, shift = 2)
        t001 = Slab(self.tei, hkl001, 10, 3, 0.025, shift = 5)
        b001 = Slab(self.batio3, hkl001, 10, 3, 0.01, shift = 2)
        l001 = Slab(self.lifepo4, hkl001, 10, 3, 0.01, shift = 2)

        hkl100 = [1, 0, 0]
        z100 = Slab(self.zno, hkl100, 10, 5, 0.01)
        t100 = Slab(self.tei, hkl100, 10, 5, 0.01)
        b100 = Slab(self.batio3, hkl100, 10, 5, 0.01, shift = 3)
        l100 = Slab(self.lifepo4, hkl100, 10, 3, 0.01, shift = 2)

        m = [z001, t001, b001, l001, z100, t100, b100, l100]

        fileName = ["pymatgen/pymatgen/core/surface tests/001_terminations/ZnO-wz.cif",
                    "pymatgen/pymatgen/core/surface tests/001_terminations/icsd_TeI.cif",
                    "pymatgen/pymatgen/core/surface tests/001_terminations/icsd_batio3.cif",
                    "pymatgen/pymatgen/core/surface tests/001_terminations/LiFePO4.cif",
                    "pymatgen/pymatgen/core/surface tests/100_terminations/ZnO-wz.cif",
                    "pymatgen/pymatgen/core/surface tests/100_terminations/icsd_TeI.cif",
                    "pymatgen/pymatgen/core/surface tests/100_terminations/icsd_batio3.cif",
                    "pymatgen/pymatgen/core/surface tests/100_terminations/LiFePO4.cif"]

        for i in range(0, len(fileName)):
            Name = fileName[i][:-4] + "_%s_%s_%s_%s_" \
                                      %(str(m[i].miller_index),
                                        m[i].min_slab_size,
                                        m[i].min_vac_size,
                                        m[i].temp_thresh)
            fileType = ".cif"

            # For visual debugging
            for iii in range(0, len(m[i].enum)):
                name_num = str(iii)
                newFile = Name + name_num + fileType
                write_structure(m[i].enum[iii], newFile)

                # Compares the newly created structure to one that was already made
                test_comp = CifParser("pymatgen/pymatgen/core/surface tests/tests/" + newFile[+68:])
                # Line 347 might change depending on where you cloned pymatgen
                self.compare_structure = (test_comp.get_structures(primitive = False)[0])
                test_new = CifParser(newFile)
                self.new_structure = (test_new.get_structures(primitive = False)[0])
                self.assertEqual(self.compare_structure, self.new_structure) #This may take a while to process

            # Prints the coordinates and the species of each atom along with the number of atoms on a
            # surface termination site
            for iiii in range(0, len(m[i].term_coords)):
                print(m[i].term_coords[iiii])
                print("%s atoms in this termination surface." %(len(m[i].term_coords[iiii])))

            print(" ")
        # Checks to see if the program generates the number of terminations we would expect it to
        self.assertEqual(len(z001.enum), 7)
        self.assertEqual(len(b001.enum), 5)
        self.assertEqual(len(t001.enum), 8)
        self.assertEqual(len(z100.enum), 7)
        self.assertEqual(len(b100.enum), 5)
        self.assertEqual(len(t100.enum), 8)


if __name__ == "__main__":
    unittest.main()


class SlabTest(unittest.TestCase):

    def setUp(self):
        self.cu = Structure(Lattice.cubic(3), ["Cu", "Cu", "Cu", "Cu"],
                            [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5],
                             [0, 0.5, 0.5]])

    def test_init(self):
        for hkl in itertools.product(xrange(4), xrange(4), xrange(4)):
            if any(hkl):
                ssize = 6
                vsize = 10
                s = Slab(self.cu, hkl, ssize, vsize)
                if hkl == [0, 1, 1]:
                    self.assertEqual(len(s), 13)
                    self.assertAlmostEqual(s.surface_area, 12.727922061357855)
                manual = self.cu.copy()
                manual.make_supercell(s.scale_factor)
                self.assertEqual(manual.lattice.lengths_and_angles,
                                 s.lattice.lengths_and_angles)

         # For visual debugging
        from pymatgen import write_structure
        write_structure(s.parent, "cu.cif")
        write_structure(s, "cu_slab_%s_%.3f_%.3f.cif" %
                          (str(hkl), ssize, vsize))

    def test_adsorb_atom(self):
        s001 = Slab(self.cu,[0, 0, 1], 5, 5)
        print s001
        #O adsorb on 4Cu[0.5, 0.5, 0.25], abc = [3, 3, 12]
         #1. test site_a = abc input
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



if __name__ == "__main__":
    unittest.main()

