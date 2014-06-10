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

    def __init__(self, structure, miller_index, min_slab_size,
                 min_vacuum_size, shift=0):
        """
        Makes a Slab structure. Note that the code will make a slab with
        whatever size that is specified, rounded upwards.

        Args:
            structure (Structure): Initial input structure.
            miller_index ([h, k, l]): Miller index of plane.
            min_slab_size (float): In Angstroms
            min_vacuum_size (float): In Angstroms
            shift (float): In Angstroms (shifting the origin)
        """
        self.parent = structure
        self.min_slab_size = min_slab_size
        self.min_vac_size = min_vacuum_size
        lattice = structure.lattice

        #Calculate the surface normal using the reciprocal lattice vector.
        recp = lattice.reciprocal_lattice_crystallographic
        normal = recp.get_cartesian_coords(miller_index)
        normal /= np.linalg.norm(normal)

        surf_scale = []
        non_orth_ind = []
        eye = np.eye(3, dtype=np.int)
        dist = float('inf')
        for i in xrange(3):
            d = np.dot(normal, lattice.matrix[i])
            if d < 1e-8:
                surf_scale.append(eye[i])
            else:
                non_orth_ind.append(i)
                if d < dist:
                    latt_index = i
                    dist = d

        lcm_index = lcm([i for i in miller_index if i != 0])
        if len(non_orth_ind) > 1:
            for i, j in itertools.combinations(non_orth_ind, 2):
                l = [0, 0, 0]
                l[i] = -int(round(lcm_index / miller_index[i]))
                l[j] = int(round(lcm_index / miller_index[j]))
                surf_scale.append(l)
                if len(surf_scale) == 2:
                    break

        nlayers_slab = int(math.ceil(min_slab_size / dist))
        nlayers_vac = int(math.ceil(min_vacuum_size / dist))
        nlayers = nlayers_slab + nlayers_vac
        surf_scale.append(eye[latt_index] * nlayers)

        slab = Structure.from_sites(structure)

        slab.make_supercell(surf_scale)
        new_sites = []
        for site in slab:
            if shift <= np.dot(site.coords, normal) < nlayers_slab * dist + \
                    shift + 1e-8:
                new_sites.append(site)
        slab = Structure.from_sites(new_sites)

        slab = slab.get_primitive_structure()
        slab = slab.copy(sanitize=True)
        super(Structure, self).__init__(slab.lattice,
                                        slab.species_and_occu,
                                        slab.frac_coords)
        self.normal = normal


import unittest
from pymatgen.core.lattice import Lattice

class SlabTest(unittest.TestCase):

    def setUp(self):
        self.cu = Structure(Lattice.cubic(3), ["Cu", "Cu", "Cu", "Cu"],
                            [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5],
                             [0, 0.5, 0.5]])

    def test_init(self):
        hkl = [1, 1, 1]
        ssize = 6
        vsize = 10
        s = Slab(self.cu, hkl, ssize, vsize)
        # For visual debugging
        from pymatgen import write_structure
        write_structure(s.parent, "cu.cif")
        write_structure(s, "cu_slab_%s_%.3f_%.3f.cif" %
                         (str(hkl), ssize, vsize))


if __name__ == "__main__":
    unittest.main()


