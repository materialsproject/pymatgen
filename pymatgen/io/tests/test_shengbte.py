# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest

from pymatgen.io.shengbte import ShengBTE_CONTROL_IO
from pymatgen.util.testing import PymatgenTest



test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "test_files", "shengbte")

this_dir = os.path.dirname(os.path.abspath(__file__))


class TestShengBTE(PymatgenTest):

    def setUp(self):
        self.sbte_io = ShengBTE_CONTROL_IO()
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

    def test_read_CONTROL(self):
        filename = os.path.join(test_dir, "CONTROL-CSLD_Si")
        sbte_dict = self.sbte_io.read_CONTROL(filename)

        self.assertIsInstance(sbte_dict, dict)
        self.assertEqual(sbte_dict['allocations']['nelements'], 1)
        self.assertEqual(sbte_dict['allocations']['natoms'], 2)
        self.assertArrayEqual(sbte_dict['allocations']['ngrid'], [25, 25, 25])
        self.assertEqual(sbte_dict['allocations']['norientations'], 0)
        self.assertEqual(sbte_dict['crystal']['lfactor'], 0.1)
        self.assertEqual(sbte_dict['crystal']['lattvec'][0], [0.0, 2.734363999, 2.734363999])
        self.assertEqual(sbte_dict['crystal']['lattvec'][1], [2.734363999, 0.0, 2.734363999])
        self.assertEqual(sbte_dict['crystal']['lattvec'][2], [2.734363999, 2.734363999, 0.0])
        self.assertIsInstance(sbte_dict['crystal']['elements'], (list, str))
        if isinstance(sbte_dict['crystal']['elements'], list):
            all_strings = all(isinstance(item, str) for item in sbte_dict['crystal']['elements'])
            self.assertTrue(all_strings)
        self.assertIsInstance(sbte_dict['crystal']['types'], (list, int))
        if isinstance(sbte_dict['crystal']['types'], list):
            all_ints = all(isinstance(item, int) for item in sbte_dict['crystal']['types'])
            self.assertTrue(all_ints)
        self.assertArrayEqual(sbte_dict['crystal']['positions'], [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
        self.assertArrayEqual(sbte_dict['crystal']['scell'], [5, 5, 5])
        self.assertEqual(sbte_dict['parameters']['T'], 500)
        self.assertEqual(sbte_dict['parameters']['scalebroad'], 0.5)
        self.assertFalse(sbte_dict['flags']['isotopes'])
        self.assertFalse(sbte_dict['flags']['onlyharmonic'])
        self.assertFalse(sbte_dict['flags']['nonanalytic'])
        self.assertFalse(sbte_dict['flags']['nanowires'])

    def test_file_writer_helper_func(self):

        if os.path.exists(os.path.join(test_dir,'test_control')):
            os.remove(os.path.join(test_dir,'test_control'))
        self.sbte_io.file_writer_helper_func(self.test_dict, filename=os.path.join(test_dir,'test_control'))

        with open(os.path.join(test_dir,'test_control'), 'r') as file:
            test_string = file.read()
        with open(os.path.join(test_dir, "CONTROL-CSLD_Si"), 'r') as reference_file:
            reference_string = reference_file.read()
        self.assertMultiLineEqual(test_string, reference_string)
        os.remove(os.path.join(test_dir, 'test_control'))

    def test_write_CONTROL_from_dict(self):
        if os.path.exists(os.path.join(test_dir,'test_control')):
            os.remove(os.path.join(test_dir,'test_control'))
        self.sbte_io.write_CONTROL_from_dict(self.test_dict, filename=os.path.join(test_dir,'test_control'))
        with open(os.path.join(test_dir,'test_control'), 'r') as file:
            test_string = file.read()
        with open(os.path.join(test_dir, "CONTROL-CSLD_Si"), 'r') as reference_file:
            reference_string = reference_file.read()
        self.assertMultiLineEqual(test_string, reference_string)
        os.remove(os.path.join(test_dir, 'test_control'))


