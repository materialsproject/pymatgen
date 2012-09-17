import unittest

from nose.exc import SkipTest

from pymatgen.command_line.aconvasp_caller import get_num_division_kpoints, get_minkowski_red, get_vasp_kpoint_file_sym
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.periodic_table import Element
import numpy as np
from pymatgen.util.io_utils import which


aconvasp_present = which('aconvasp')


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
        if not aconvasp_present:
            raise SkipTest("aconvasp not present. Skipping...")
        self.assertListEqual(get_num_division_kpoints(self.struct, 500), [6, 7, 6])

    def test_get_minkowski_red(self):
        if not aconvasp_present:
            raise SkipTest("aconvasp not present. Skipping...")
        new_struct = get_minkowski_red(self.struct)
        self.assertAlmostEqual(new_struct.lattice.a, 3.840198)
        self.assertAlmostEqual(new_struct.lattice.alpha, 60.0)
        self.assertEqual(new_struct.species_and_occu[0], {Element("Si"): 1})
        self.assertEqual(new_struct.frac_coords[1][0], 0.25)

    def test_get_vasp_kpoint_file_sym(self):
        if not aconvasp_present:
            raise SkipTest("aconvasp not present. Skipping...")
        self.assertEqual(get_vasp_kpoint_file_sym(self.struct).split("\n")[0], "FCC (face-centered cubic) G-X-W-K-G-L-U-W-L-K U-X")


if __name__ == '__main__':
    unittest.main()
