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
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath

test_dir_structs = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                            'test_files', 'space_group_structs')

class HighSymmKpathTest(PymatgenTest):

    def test_kpath_generation(self):
        for sg_num in range(230):
            struct = Structure.from_file(os.path.join(test_dir_structs, str(sg_num+1) + '.json'))
            kpath = HighSymmKpath(struct) #Throws error if something doesn't work, causing test to fail.

    def test_kpath_acentered(self):
        struct = Structure.from_file(os.path.join(test_dir_structs, '38.json'))
        kpath = HighSymmKpath(struct)
        
        kpoints = kpath._kpath['kpoints']
        labels = list(kpoints.keys())

        self.assertEqual(labels, ['\\Gamma', 'A', 'A_1', 'R', 'S', 'T', 'X', 'X_1', 'Y', 'Z'])

        self.assertEqual(kpoints['\\Gamma'][0], 0.00000000)
        self.assertAlmostEqual(kpoints['\\Gamma'][1], 0.00000000)
        self.assertAlmostEqual(kpoints['\\Gamma'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['A'][0], 0.26217672)
        self.assertAlmostEqual(kpoints['A'][1], 0.26217672)
        self.assertAlmostEqual(kpoints['A'][2], 0.50000000)

        self.assertAlmostEqual(kpoints['A_1'][0], -0.26217672)
        self.assertAlmostEqual(kpoints['A_1'][1], 0.73782328)
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

        self.assertAlmostEqual(kpoints['X'][0], 0.26217672)
        self.assertAlmostEqual(kpoints['X'][1], 0.26217672)
        self.assertAlmostEqual(kpoints['X'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['X_1'][0], -0.26217672)
        self.assertAlmostEqual(kpoints['X_1'][1], 0.73782328)
        self.assertAlmostEqual(kpoints['X_1'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['Y'][0], -0.50000000)
        self.assertAlmostEqual(kpoints['Y'][1], 0.50000000)
        self.assertAlmostEqual(kpoints['Y'][2], 0.00000000)

        self.assertAlmostEqual(kpoints['Z'][0], 0.00000000)
        self.assertAlmostEqual(kpoints['Z'][1], 0.00000000)
        self.assertAlmostEqual(kpoints['Z'][2], 0.50000000)

if __name__ == "__main__":
    unittest.main()
