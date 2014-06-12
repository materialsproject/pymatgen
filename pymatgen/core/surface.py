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

        self.parent = structure
        self.min_slab_size = min_slab_size
        self.min_vac_size = min_vacuum_size
        self.scale_factor = np.array(slab_scale_factor)
        self.normal = normal
        self.miller_index = miller_index

        super(Structure, self).__init__(
            slab.lattice, slab.species_and_occu, slab.frac_coords,
            site_properties=slab.site_properties)

    @property
    def surface_area(self):
        m = self.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))


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


if __name__ == "__main__":
    unittest.main()


