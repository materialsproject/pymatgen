# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from pymatgen.util.testing import PymatgenTest
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.io.xcrysden import XSF


class XSFTest(PymatgenTest):
    def test_xsf(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        structure = Structure(lattice, ["Si", "Si"], coords)
        xsf = XSF(structure)
        self.assertTrue(structure, XSF.from_string(xsf.to_string()))

    def test_xsf_symbolparse(self):
        """
        Ensure that the same structure is parsed
        even if the atomic symbol / number convention
        is different.
        """

        test_string = """
CRYSTAL
PRIMVEC
       11.45191956     0.00000000     0.00000000
        5.72596044     9.91765288     0.00000000
      -14.31490370    -8.26471287    23.37613199
PRIMCOORD
1 1
H     -0.71644986    -0.41364333     1.19898200     0.00181803     0.00084718     0.00804832
"""
        structure = XSF.from_string(test_string).structure
        self.assertEqual(str(structure.species[0]), 'H')
        test_string2 = """
CRYSTAL
PRIMVEC
       11.45191956     0.00000000     0.00000000
        5.72596044     9.91765288     0.00000000
      -14.31490370    -8.26471287    23.37613199
PRIMCOORD
1 1
1     -0.71644986    -0.41364333     1.19898200     0.00181803     0.00084718     0.00804832
"""

        structure2 = XSF.from_string(test_string2).structure
        self.assertEqual(structure, structure2)


if __name__ == "__main__":
    import unittest

    unittest.main()
