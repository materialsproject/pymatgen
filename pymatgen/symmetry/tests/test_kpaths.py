# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Aug 23, 2017
"""


__author__ = "Katie Latimer"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Katie Latimer"
__email__ = "klatimer@berkeley.edu"
__date__ = "Aug 23, 2017"

import unittest
import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.bandstructure import HighSymmKpath

test_dir_structs = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                            'test_files', 'space_group_structs')

class HighSymmKpathTest(PymatgenTest):

    def test_kpath_generation(self):
        triclinic = [1, 2]
        monoclinic = range(3, 16)
        orthorhombic = range(16, 75)
        tetragonal = range(75, 143)
        rhombohedral = range(143, 168)
        hexagonal = range(168, 195)
        cubic = range(195, 231)
        
        species = ['K', 'La', 'Ti']
        coords = [[.345, 5, .77298], [.1345, 5.1, .77298], [.7, .8, .9]]
        for i in range(230):
            sg_num = i + 1
            if sg_num in triclinic:
                lattice = Lattice([[3.0233057319441246,0,0], [0,7.9850357844548681,0], [0,0,8.1136762279561818]])
            elif sg_num in monoclinic:
                lattice = Lattice.monoclinic(2, 9, 1, 99)
            elif sg_num in orthorhombic:
                lattice = Lattice.orthorhombic(2, 9, 1)
            elif sg_num in tetragonal:
                lattice = Lattice.tetragonal(2, 9)
            elif sg_num in rhombohedral:
                lattice = Lattice.hexagonal(2, 95)
            elif sg_num in hexagonal:
                lattice = Lattice.hexagonal(2, 9)
            elif sg_num in cubic:
                lattice = Lattice.cubic(2)
        
            struct = Structure.from_spacegroup(sg_num, lattice, species, coords)
            kpath = HighSymmKpath(struct) #Throws error if something doesn't work, causing test to fail.

        struct_file_path = os.path.join(test_dir_structs, 'ICSD_170.cif')
        struct = Structure.from_file(struct_file_path)
        hkp = HighSymmKpath(struct)
        self.assertEqual(hkp.name, 'MCLC5')

    def test_kpath_acentered(self):
        species = ['K', 'La', 'Ti']
        coords = [[.345, 5, .77298], [.1345, 5.1, .77298], [.7, .8, .9]]
        lattice = Lattice.orthorhombic(2, 9, 1)
        struct = Structure.from_spacegroup(38, lattice, species, coords)
        kpath = HighSymmKpath(struct)
        
        kpoints = kpath._kpath['kpoints']
        labels = list(kpoints.keys())

        self.assertEqual(sorted(labels), sorted(['\\Gamma', 'A', 'A_1', 'R', 'S', 'T', 'X', 'X_1', 'Y', 'Z']))

        self.assertEqual(kpoints['\\Gamma'][0], 0.00000000)
        self.assertAlmostEqual(kpoints['\\Gamma'][1], 0.00000000)
        self.assertAlmostEqual(kpoints['\\Gamma'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['A'][0], 0.25308642)
        self.assertAlmostEqual(kpoints['A'][1], 0.25308642)
        self.assertAlmostEqual(kpoints['A'][2], 0.50000000)

        self.assertAlmostEqual(kpoints['A_1'][0], -0.25308642)
        self.assertAlmostEqual(kpoints['A_1'][1], 0.74691358)
        self.assertAlmostEqual(kpoints['A_1'][2], 0.50000000)

        self.assertAlmostEqual(kpoints['R'][0], 0.00000000)
        self.assertAlmostEqual(kpoints['R'][1], 0.50000000)
        self.assertAlmostEqual(kpoints['R'][2], 0.50000000)

        self.assertAlmostEqual(kpoints['S'][0], 0.00000000)
        self.assertAlmostEqual(kpoints['S'][1], 0.50000000)
        self.assertAlmostEqual(kpoints['S'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['T'][0], -0.50000000)
        self.assertAlmostEqual(kpoints['T'][1], 0.50000000)
        self.assertAlmostEqual(kpoints['T'][2], 0.50000000)

        self.assertAlmostEqual(kpoints['X'][0], 0.25308642)
        self.assertAlmostEqual(kpoints['X'][1], 0.25308642)
        self.assertAlmostEqual(kpoints['X'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['X_1'][0], -0.25308642)
        self.assertAlmostEqual(kpoints['X_1'][1], 0.74691358)
        self.assertAlmostEqual(kpoints['X_1'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['Y'][0], -0.50000000)
        self.assertAlmostEqual(kpoints['Y'][1], 0.50000000)
        self.assertAlmostEqual(kpoints['Y'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['Z'][0], 0.00000000)
        self.assertAlmostEqual(kpoints['Z'][1], 0.00000000)
        self.assertAlmostEqual(kpoints['Z'][2], 0.50000000)

if __name__ == "__main__":
    unittest.main()
