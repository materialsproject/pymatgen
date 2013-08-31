#!/usr/bin/python

from __future__ import division

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.units import Energy, Time, Length, unitized


class UnitTest(PymatgenTest):

    def test_energy(self):
        a = Energy(1.1, "eV")
        b = a.to("Ha")
        self.assertAlmostEqual(b, 0.0404242579378)
        c = Energy(3.14, "J")
        self.assertAlmostEqual(c.to("eV"), 1.95983393276e+19)
        self.assertRaises(ValueError, Energy, 1, "m")

        d = Energy(1, "Ha")
        self.assertAlmostEqual(a + d, 28.31138386)
        self.assertAlmostEqual(a - d, -26.11138386)
        self.assertEqual(a + 1, 2.1)

    def test_time(self):
        a = Time(20, "h")
        self.assertAlmostEqual(a.to("s"), 3600 * 20)
        #Test left and right multiplication.
        self.assertEqual(str(a * 3), "60.0 h")
        self.assertEqual(str(3 * a), "60.0 h")

    def test_length(self):
        x = Length(4.2, "ang")
        self.assertEqual(x.to("cm"), 4.2e-08)
        self.assertEqual(x.to("pm"), 420)
        self.assertEqual(str(x / 2), "2.1 ang")

    def test_unitized(self):

        @unitized("energy", "eV")
        def f():
            return [1, 2, 3]

        self.assertEqual(str(f()[0]), "1.0 eV")
        self.assertIsInstance(f(), list)

        @unitized("energy", "eV")
        def f():
            return 2, 3, 4

        self.assertEqual(str(f()[0]), "2.0 eV")
        self.assertIsInstance(f(), tuple)

if __name__ == '__main__':
    import unittest
    unittest.main()
