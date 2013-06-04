#!/usr/bin/env python

'''
'''

from __future__ import division


import unittest
import os

from pymatgen.io.zeoio import ZeoCssr
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.core.structure import Structure

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class ZeoCssrTest(unittest.TestCase):

    def setUp(self):

        filepath = os.path.join(test_dir, 'POSCAR')
        p = Poscar.from_file(filepath)
        self.zeocssr = ZeoCssr(p.structure)

    def test_str(self):
        expected_string = """10.4118 6.0672 4.7595
90.00 90.00 90.00 SPGR =  1 P 1    OPT = 1
24 0
0 Fe4 P4 O16
1 Fe 0.2187 0.7500 0.4749 0 0 0 0 0 0 0 0 0.0000
2 Fe 0.2813 0.2500 0.9749 0 0 0 0 0 0 0 0 0.0000
3 Fe 0.7187 0.7500 0.0251 0 0 0 0 0 0 0 0 0.0000
4 Fe 0.7813 0.2500 0.5251 0 0 0 0 0 0 0 0 0.0000
5 P 0.0946 0.2500 0.4182 0 0 0 0 0 0 0 0 0.0000
6 P 0.4054 0.7500 0.9182 0 0 0 0 0 0 0 0 0.0000
7 P 0.5946 0.2500 0.0818 0 0 0 0 0 0 0 0 0.0000
8 P 0.9054 0.7500 0.5818 0 0 0 0 0 0 0 0 0.0000
9 O 0.0434 0.7500 0.7071 0 0 0 0 0 0 0 0 0.0000
10 O 0.0966 0.2500 0.7413 0 0 0 0 0 0 0 0 0.0000
11 O 0.1657 0.0461 0.2854 0 0 0 0 0 0 0 0 0.0000
12 O 0.1657 0.4539 0.2854 0 0 0 0 0 0 0 0 0.0000
13 O 0.3343 0.5461 0.7854 0 0 0 0 0 0 0 0 0.0000
14 O 0.3343 0.9539 0.7854 0 0 0 0 0 0 0 0 0.0000
15 O 0.4034 0.7500 0.2413 0 0 0 0 0 0 0 0 0.0000
16 O 0.4566 0.2500 0.2071 0 0 0 0 0 0 0 0 0.0000
17 O 0.5434 0.7500 0.7929 0 0 0 0 0 0 0 0 0.0000
18 O 0.5966 0.2500 0.7587 0 0 0 0 0 0 0 0 0.0000
19 O 0.6657 0.0461 0.2146 0 0 0 0 0 0 0 0 0.0000
20 O 0.6657 0.4539 0.2146 0 0 0 0 0 0 0 0 0.0000
21 O 0.8343 0.5461 0.7146 0 0 0 0 0 0 0 0 0.0000
22 O 0.8343 0.9539 0.7146 0 0 0 0 0 0 0 0 0.0000
23 O 0.9034 0.7500 0.2587 0 0 0 0 0 0 0 0 0.0000
24 O 0.9566 0.2500 0.2929 0 0 0 0 0 0 0 0 0.0000"""
        self.assertEqual(str(self.zeocssr), expected_string)

    def test_from_file(self):
        filename = os.path.join(test_dir, "EDI.cssr")
        zeocssr = ZeoCssr.from_file(filename)
        self.assertIsInstance(zeocssr.structure, Structure)


if __name__ == "__main__":
    unittest.main()
