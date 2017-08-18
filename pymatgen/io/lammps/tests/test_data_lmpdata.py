# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import subprocess
import os

import numpy as np
from monty.os.path import which
from pymatgen import Lattice, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from pymatgen.io.lammps.data import LMPData


class LMPDataTest(unittest.TestCase):

    def setUp(self):
        latt = Lattice(np.eye(3) * 10 + np.random.randn(3, 3))
        natoms = np.random.randint(1, 10)
        structure = Structure(latt, ['H'] * natoms, np.random.rand(natoms, 3))
        structure.add_oxidation_state_by_site(np.random.rand(natoms))
        ld = LMPData(structure)
        self.ld = ld
        self.structure = structure

    def test_init(self):
        self.assertEqual(len(self.ld._elements), 1)
        aligned_structure = self.ld.structure
        sm = StructureMatcher(scale=False)
        self.assertTrue(sm.fit(self.structure, aligned_structure))
        a_axis = aligned_structure.lattice.matrix[0]
        self.assertEqual(np.dot(a_axis, np.array([0, 1, 0])), 0)
        self.assertEqual(np.dot(a_axis, np.array([0, 0, 1])), 0)

    @unittest.skipIf(not which('lmp_serial'), 'No LAMMPS cmd found.')
    def test_write_file(self):
        self.ld.write_file('test')
        test_cmds = ['units metal', 'atom_style charge',
                     'box tilt large', 'read_data data.test']
        with open('in.test', 'w') as f:
            f.write('\n'.join(test_cmds))
        p = subprocess.Popen(['lmp_serial', '-in', 'in.test'],
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                             close_fds=True)
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def tearDown(self):
        all_tmp_files = ['data.test', 'in.test', 'log.lammps']
        tmp_files = [f for f in all_tmp_files if os.path.exists(f)]
        for t in tmp_files:
            os.remove(t)


if __name__ == '__main__':
    unittest.main()

