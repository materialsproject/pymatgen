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
                 min_vacuum_size, lll_reduce=True, shift=0):
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

        nlayers_slab = int(math.ceil(min_slab_size / dist))
        nlayers_vac = int(math.ceil(min_vacuum_size / dist))
        nlayers = nlayers_slab + nlayers_vac
        slab_scale_factor.append(eye[latt_index] * nlayers)

        slab = structure.copy()

        slab.make_supercell(slab_scale_factor)
        new_sites = []
        for site in slab:
            if shift <= np.dot(site.coords, normal) < nlayers_slab * dist + \
                    shift:
                new_sites.append(site)
        slab = Structure.from_sites(new_sites)

        if lll_reduce:
            lll_slab = slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(slab.lattice)
            slab_scale_factor = np.dot(mapping[2], slab_scale_factor)
            slab = lll_slab

<<<<<<< HEAD
        self.parent = structure
=======
        n = 0
        term_slab = slab.copy()
        c = term_slab.lattice.c
        # For loop moves all sites down to compensate for the space opened up by the shift
        for site in term_slab:
            index = []
            index.append(n)
            term_slab.translate_sites(index, [0, 0, -shift/c])
            n+=1

        # Rescales the  lattice
        new_latt = Lattice.from_parameters(term_slab.lattice.a, term_slab.lattice.b, min_vacuum_size+nlayers_slab*dist,
                                           term_slab.lattice.alpha, term_slab.lattice.beta, term_slab.lattice.gamma)
        term_slab.modify_lattice(new_latt)

        el = term_slab.species
        org_coords = term_slab.frac_coords.tolist()
        new_coord, b = [], []

        for i in org_coords:
            i[2] = (i[2]*c)/term_slab.lattice.c

        for i in term_slab.frac_coords:
            b.append(i[2])
            new_coord.append(b)
            b = []

        # Clusters sites together that belong in the same termination surface based on their position in the
        # c direction. Also organizes sites by ascending order from teh origin along the c direction. How close
        # the sites have to be to belong in the same termination surface depends on the user input thresh.
        tracker_index = scipy.cluster.hierarchy.fclusterdata(new_coord, thresh, criterion= crit)
        new_coord, tracker_index, org_coords, el = zip(*sorted(zip(new_coord, tracker_index, org_coords, el)))

# Creates a list (term_index) that tells us which at which site does a termination begin. For 1 unit cell.
        term_index = []
        gg = 0
        for i in range(0, len(term_slab)):
            gg +=1
            if i == len(term_slab) - 1:
                term_index.append(i)
            else:
                if tracker_index[i] != tracker_index[i+1]:
                    term_index.append(i)
            if gg == len(structure):
                    break

        slab_list, term_coords = [], []
        a, i = 0, 0
        new_slab = Structure(term_slab.lattice, el, org_coords)

        for iii in range(0, len(term_index)):
            y = []
            alt_slab = new_slab.copy()

            for ii in range(0, len(alt_slab)):
                index = []
                term_scale = 0
                index.append(ii)

                if alt_slab.frac_coords[ii][2] > alt_slab.frac_coords[term_index[iii]][2]:
                    term_scale = min_vacuum_size

                alt_slab.translate_sites(index, [0, 0, term_scale/(new_slab.lattice.c)])

            if standardize:
                index = []
                for f in range(0, len(alt_slab)):
                    index.append(f)
                if alt_slab.frac_coords[f][2] > alt_slab.frac_coords[term_index[iii]][2]:
                    standard_shift = -(alt_slab.frac_coords[term_index[iii]][2] + (0.5*min_vacuum_size)/alt_slab.lattice.c)
                else:
                    standard_shift = 1 - alt_slab.frac_coords[term_index[iii]][2] - (0.5*min_vacuum_size)/alt_slab.lattice.c
                alt_slab.translate_sites(index, [0, 0, standard_shift])

            slab_list.append(alt_slab)

            for iv in range(a, term_index[iii] + 1):
                y.append(slab_list[i][iv])

            i += 1
            term_coords.append(y)
            a = term_index[iii]+1

>>>>>>> parent of 16269eb... Latest commit
        self.min_slab_size = min_slab_size
        self.min_vac_size = min_vacuum_size
        self.scale_factor = np.array(slab_scale_factor)
        self.normal = normal
        self.miller_index = miller_index

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


import unittest
from pymatgen.core.lattice import Lattice

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



if __name__ == "__main__":
    unittest.main()


