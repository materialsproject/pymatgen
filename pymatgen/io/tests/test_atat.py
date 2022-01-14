# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os

from pymatgen.core.structure import Structure
from pymatgen.io.atat import Mcsqs
from pymatgen.util.testing import PymatgenTest
from pymatgen.util.testing import PymatgenTest


test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "mcsqs")


class AtatTest(PymatgenTest):
    def test_mcsqs_import(self):
        test_string = """1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
0.000000 -1.000000 -2.000000
2.000000 -1.000000 0.000000
-1.000000 -1.000000 1.000000
0.000000 -2.000000 -1.000000 Mn
1.000000 -2.000000 -1.000000 Mn
0.000000 -1.000000 -1.000000 Mn
-0.000000 -2.000000 0.000000 Mn
1.000000 -2.000000 0.000000 Mn
0.000000 -1.000000 0.000000 Mn
1.000000 -1.000000 0.000000 Fe
1.000000 -3.000000 -1.000000 Mn
0.500000 -1.500000 -0.500000 Sr
1.500000 -1.500000 -0.500000 Ca
-0.500000 -1.500000 0.500000 Ca
0.500000 -1.500000 0.500000 Ca
1.500000 -2.500000 -1.500000 Ca
0.500000 -1.500000 -1.500000 Sr
0.500000 -2.500000 -0.500000 Sr
-0.500000 -1.500000 -0.500000 Ca
0.000000 -1.500000 -1.000000 O
1.000000 -1.500000 -1.000000 O
1.000000 -2.500000 0.000000 O
-0.000000 -1.500000 0.000000 O
1.000000 -1.500000 0.000000 O
0.000000 -0.500000 0.000000 O
0.000000 -2.500000 -1.000000 O
1.000000 -2.500000 -1.000000 O
0.500000 -2.000000 -1.000000 O
1.500000 -2.000000 -1.000000 O
0.500000 -1.000000 -1.000000 O
0.500000 -2.000000 0.000000 O
-0.500000 -1.000000 0.000000 O
0.500000 -1.000000 0.000000 O
1.500000 -1.000000 0.000000 O
-0.500000 -2.000000 -1.000000 O
0.000000 -2.000000 -0.500000 O
1.000000 -2.000000 -0.500000 O
0.000000 -1.000000 -0.500000 O
1.000000 -1.000000 -0.500000 O
1.000000 -2.000000 0.500000 O
0.000000 -1.000000 0.500000 O
1.000000 -2.000000 -1.500000 O
0.000000 -1.000000 -1.500000 O
"""

        s = Mcsqs.structure_from_string(test_string)

        self.assertEqual(s.composition.formula, "Sr3 Ca5 Mn7 Fe1 O24")
        self.assertAlmostEqual(s.lattice.a, 2.2360679775)
        self.assertAlmostEqual(s.lattice.b, 2.2360679775)
        self.assertAlmostEqual(s.lattice.c, 1.73205080757)

    def test_mcsqs_export(self):
        s = self.get_structure("SrTiO3")
        s.replace_species({"Sr2+": {"Sr2+": 0.5, "Ca2+": 0.5}})

        ref_string = """3.905000 0.000000 0.000000
-0.000000 3.905000 0.000000
0.000000 0.000000 3.905000
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
0.500000 0.500000 0.500000 Sr2+=0.5,Ca2+=0.5
0.000000 0.000000 0.000000 Ti4+=1.0
0.000000 0.000000 0.500000 O2-=1.0
0.000000 0.500000 0.000000 O2-=1.0
0.500000 0.000000 0.000000 O2-=1.0"""

        self.assertEqual(Mcsqs(s).to_string(), ref_string)

    def test_mcsqs_cif_nacl(self):

        # cif file from str2cif (utility distributed with atat)
        struc_from_cif = Structure.from_file(os.path.join(test_dir, "bestsqs_nacl.cif"))

        # output file directly from mcsqs
        struc_from_out = Structure.from_file(os.path.join(test_dir, "bestsqs_nacl.out"))

        self.assertTrue(struc_from_cif.matches(struc_from_out))
        self.assertArrayAlmostEqual(
            struc_from_out.lattice.parameters,
            struc_from_cif.lattice.parameters,
            decimal=4,
        )

    def test_mcsqs_cif_pzt(self):
        # cif file from str2cif (utility distributed with atat)
        struc_from_cif = Structure.from_file(os.path.join(test_dir, "bestsqs_pzt.cif"))

        # output file directly from mcsqs
        struc_from_out = Structure.from_file(os.path.join(test_dir, "bestsqs_pzt.out"))

        self.assertTrue(struc_from_cif.matches(struc_from_out))
        self.assertArrayAlmostEqual(
            struc_from_out.lattice.parameters,
            struc_from_cif.lattice.parameters,
            decimal=4,
        )
