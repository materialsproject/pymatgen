#!/usr/bin/python

from __future__ import division

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.units import Unit, Energy


class UnitTest(PymatgenTest):

    def test_energy(self):
        a = Energy(1.1, "eV")
        b = a.to("Ha")
        self.assertAlmostEqual(b, 0.0404242579378)
        c = Energy(3.14, "J")
        self.assertAlmostEqual(c.to("eV"), 1.95983393276e+19)
        self.assertRaises(ValueError, Energy, 1, "m")


if __name__ == '__main__':
    import unittest
    unittest.main()
