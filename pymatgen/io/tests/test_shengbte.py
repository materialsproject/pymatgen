# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest

from pymatgen.io.shengbte import Control
from pymatgen.util.testing import PymatgenTest



test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "test_files", "shengbte")

this_dir = os.path.dirname(os.path.abspath(__file__))


class TestShengBTE(PymatgenTest):

    def setUp(self):
        self.filename = os.path.join(test_dir, "CONTROL-CSLD_Si")
        self.test_dict = {
            'allocations':
                {'nelements': 1,
                 'natoms': 2,
                 'ngrid': [25, 25, 25],
                 'norientations': 0},
            'crystal':
                {'lfactor': 0.1,
                 'lattvec': [[0.0, 2.734363999, 2.734363999],
                             [2.734363999, 0.0, 2.734363999],
                             [2.734363999, 2.734363999, 0.0]],
                 'elements': 'Si',
                 'types': [1, 1],
                 'positions': [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]],
                 'scell': [5, 5, 5]},
            'parameters':
                {'T': 500,
                 'scalebroad': 0.5},
            'flags':
                {'isotopes': False,
                 'onlyharmonic': False,
                 'nonanalytic': False,
                 'nanowires': False}
        }

    def test_from_file(self):
        io = Control.from_file(self.filename)
        self.assertIsInstance(io.alloc_dict, dict)
        self.assertIsInstance(io.crystal_dict, dict)
        self.assertIsInstance(io.params_dict, dict)
        self.assertIsInstance(io.flags_dict, dict)
        self.assertEqual(io.alloc_dict['nelements'], 1)
        self.assertEqual(io.alloc_dict['natoms'], 2)
        self.assertArrayEqual(io.alloc_dict['ngrid'], [25, 25, 25])
        self.assertEqual(io.alloc_dict['norientations'], 0)
        self.assertEqual(io.crystal_dict['lfactor'], 0.1)
        self.assertEqual(io.crystal_dict['lattvec'][0], [0.0, 2.734363999, 2.734363999])
        self.assertEqual(io.crystal_dict['lattvec'][1], [2.734363999, 0.0, 2.734363999])
        self.assertEqual(io.crystal_dict['lattvec'][2], [2.734363999, 2.734363999, 0.0])
        self.assertIsInstance(io.crystal_dict['elements'], (list, str))
        if isinstance(io.crystal_dict['elements'], list):
            all_strings = all(isinstance(item, str) for item in io.crystal_dict['elements'])
            self.assertTrue(all_strings)
        self.assertIsInstance(io.crystal_dict['types'], (list, int))
        if isinstance(io.crystal_dict['types'], list):
            all_ints = all(isinstance(item, int) for item in io.crystal_dict['types'])
            self.assertTrue(all_ints)
        self.assertArrayEqual(io.crystal_dict['positions'], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
        self.assertArrayEqual(io.crystal_dict['scell'], [5, 5, 5])
        self.assertEqual(io.params_dict['T'], 500)
        self.assertEqual(io.params_dict['scalebroad'], 0.5)
        self.assertFalse(io.flags_dict['isotopes'])
        self.assertFalse(io.flags_dict['onlyharmonic'])
        self.assertFalse(io.flags_dict['nonanalytic'])
        self.assertFalse(io.flags_dict['nanowires'])

        if os.path.exists(os.path.join(test_dir,'test_control')):
            os.remove(os.path.join(test_dir,'test_control'))
        io.to_file(filename=os.path.join(test_dir,'test_control'))

        with open(os.path.join(test_dir,'test_control'), 'r') as file:
            test_string = file.read()
        with open(os.path.join(test_dir, "CONTROL-CSLD_Si"), 'r') as reference_file:
            reference_string = reference_file.read()
        self.assertMultiLineEqual(test_string, reference_string)
        os.remove(os.path.join(test_dir, 'test_control'))

    def test_from_dict(self):
        io = Control.from_dict(self.test_dict)
        if os.path.exists(os.path.join(test_dir,'test_control')):
            os.remove(os.path.join(test_dir,'test_control'))
        io.to_file(filename=os.path.join(test_dir,'test_control'))
        with open(os.path.join(test_dir,'test_control'), 'r') as file:
            test_string = file.read()
        with open(os.path.join(test_dir, "CONTROL-CSLD_Si"), 'r') as reference_file:
            reference_string = reference_file.read()
        self.assertMultiLineEqual(test_string, reference_string)
        os.remove(os.path.join(test_dir, 'test_control'))


