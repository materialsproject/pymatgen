# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest

from pymatgen.command_line.aconvasp_caller import get_num_division_kpoints, \
    get_minkowski_red, get_vasp_kpoint_file_sym
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.periodic_table import Element
from monty.os.path import which


aconvasp_present = which('aconvasp')
aconvasp_present = False  # disable aconvasp testing for now.


@unittest.skipIf(not aconvasp_present, "aconvasp not present.")
class AconvaspCallerTest(unittest.TestCase):

    def setUp(self):
        self.si = Element("Si")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        self.lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        self.struct = Structure(self.lattice, [self.si, self.si], coords)

    def test_get_num_division_kpoints(self):
        self.assertListEqual(get_num_division_kpoints(self.struct, 500),
                             [6, 7, 6])

    def test_get_minkowski_red(self):
        new_struct = get_minkowski_red(self.struct)
        self.assertAlmostEqual(new_struct.lattice.a, 3.840198)
        self.assertAlmostEqual(new_struct.lattice.alpha, 60.0)
        self.assertEqual(new_struct.species_and_occu[0], Composition("Si"))
        self.assertEqual(new_struct.frac_coords[1][0], 0.25)

    def test_get_vasp_kpoint_file_sym(self):
        self.assertEqual(get_vasp_kpoint_file_sym(self.struct).split("\n")[0],
                         "FCC (face-centered cubic) G-X-W-K-G-L-U-W-L-K U-X")


if __name__ == '__main__':
    unittest.main()
