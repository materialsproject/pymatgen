# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


__author__ = "Nils Edvin Richard Zimmermann"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Nils Edvin Richard Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__date__ = "June 23, 2016"

import os
import unittest

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.xr import Xr
from pymatgen.util.testing import PymatgenTest


class XrTest(unittest.TestCase):
    def setUp(self):
        filepath = os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR")
        p = Poscar.from_file(filepath)
        self.xr = Xr(p.structure)

    def test_str(self):
        expected_string = """pymatgen   10.4118 6.0672 4.7595
90.000 90.000 90.000
24 0
0 Fe4 P4 O16
1 Fe 2.2773 4.5504 2.2601
2 Fe 2.9285 1.5168 4.6399
3 Fe 7.4832 4.5504 0.1196
4 Fe 8.1344 1.5168 2.4994
5 P 0.9851 1.5168 1.9906
6 P 4.2208 4.5504 4.3704
7 P 6.1910 1.5168 0.3891
8 P 9.4267 4.5504 2.7689
9 O 0.4516 4.5504 3.3656
10 O 1.0062 1.5168 3.5283
11 O 1.7253 0.2795 1.3583
12 O 1.7253 2.7541 1.3583
13 O 3.4806 3.3131 3.7380
14 O 3.4806 5.7876 3.7380
15 O 4.1997 4.5504 1.1486
16 O 4.7543 1.5168 0.9859
17 O 5.6575 4.5504 3.7736
18 O 6.2121 1.5168 3.6109
19 O 6.9312 0.2795 1.0215
20 O 6.9312 2.7541 1.0215
21 O 8.6864 3.3131 3.4012
22 O 8.6864 5.7876 3.4012
23 O 9.4055 4.5504 1.2312
24 O 9.9602 1.5168 1.3939
10.4118 0.0000 0.0000
0.0000 6.0672 0.0000
0.0000 0.0000 4.7595
10.4118 0.0000 0.0000
0.0000 6.0672 0.0000
0.0000 0.0000 4.7595"""
        self.assertEqual(str(self.xr), expected_string)

    def test_from_file(self):
        filename = os.path.join(PymatgenTest.TEST_FILES_DIR, "EDI.xr")
        xr = Xr.from_file(filename)
        self.assertIsInstance(xr.structure, Structure)
        xr2 = Xr.from_file(filename, use_cores=False)
        self.assertIsInstance(xr2.structure, Structure)


if __name__ == "__main__":
    unittest.main()
