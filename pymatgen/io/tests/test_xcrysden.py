# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import unicode_literals, division, print_function


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
        self.assertTrue(structure, xsf.from_string(xsf.to_string()))


if __name__ == "__main__":
    import unittest
    unittest.main()
