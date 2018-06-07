# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest

from pymatgen.core.structure import Structure
from pymatgen.command_line.critic2_caller import *
from monty.os.path import which

__author__ = "Matthew Horton"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "July 2017"


@unittest.skipIf(not which('critic2'), "critic2 executable not present")
class Critic2CallerTest(unittest.TestCase):

    def test_from_path(self):
        # uses chgcars
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files/bader')

        c2c = Critic2Caller.from_path(test_dir)

        # check we have some results!
        self.assertGreaterEqual(len(c2c._stdout), 500)

    def test_from_structure(self):
        # uses promolecular density
        structure = Structure.from_file(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                     'test_files/critic2/MoS2.cif'))

        c2c = Critic2Caller(structure)

        # check we have some results!
        self.assertGreaterEqual(len(c2c._stdout), 500)


class Critic2OutputTest(unittest.TestCase):

    def setUp(self):
        stdout_file = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                   'test_files/critic2/MoS2_critic2_stdout.txt')
        with open(stdout_file, 'r') as f:
            reference_stdout = f.read()

        structure = Structure.from_file(os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                        'test_files/critic2/MoS2.cif'))

        self.c2o = Critic2Output(structure, reference_stdout)

    def test_properties_to_from_dict(self):

        self.assertEqual(len(self.c2o.critical_points), 6)
        self.assertEqual(len(self.c2o.nodes), 14)
        self.assertEqual(len(self.c2o.edges), 10)

        # reference dictionary for c2o.critical_points[0].as_dict()
        # {'@class': 'CriticalPoint',
        #  '@module': 'pymatgen.command_line.critic2_caller',
        #  'coords': None,
        #  'field': 93848.0413,
        #  'field_gradient': 0.0,
        #  'field_hessian': [[-2593274446000.0, -3.873587547e-19, -1.704530713e-08],
        #                    [-3.873587547e-19, -2593274446000.0, 1.386877485e-18],
        #                    [-1.704530713e-08, 1.386877485e-18, -2593274446000.0]],
        #  'frac_coords': [0.333333, 0.666667, 0.213295],
        #  'index': 0,
        #  'multiplicity': 1.0,
        #  'point_group': 'D3h',
        #  'type': < CriticalPointType.nucleus: 'nucleus' >}

        self.assertEqual(str(self.c2o.critical_points[0].type), "CriticalPointType.nucleus")

        # test connectivity
        self.assertDictEqual(self.c2o.edges[3], {'from_idx': 1, 'from_lvec': (0, 0, 0),
                                                 'to_idx': 0, 'to_lvec': (1, 0, 0)})
        # test as/from dict
        d = self.c2o.as_dict()
        self.assertEqual(set(d.keys()), {'@module', '@class',
                                         'structure', 'critic2_stdout'})
        self.c2o.from_dict(d)

if __name__ == '__main__':
    unittest.main()
