#!/usr/bin/python
from __future__ import division

import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.units import (Energy, Time, Length, unitized, Angle,
                                 EnergyArray, TimeArray, LengthArray,
                                 Unit, SUPPORTED_UNITS,
                                 )

import collections
import math

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

    def test_angle(self):
        a = Angle(90, "deg")
        self.assertEqual(a.to("rad"), math.pi / 2)

    def test_unitized(self):

        @unitized("eV")
        def f():
            return [1, 2, 3]

        self.assertEqual(str(f()[0]), "1.0 eV")
        self.assertIsInstance(f(), list)

        @unitized("eV")
        def g():
            return 2, 3, 4

        self.assertEqual(str(g()[0]), "2.0 eV")
        self.assertIsInstance(g(), tuple)

        @unitized("pm")
        def h():
            d = collections.OrderedDict()
            for i in range(3):
                d[i] = i * 20
            return d

        self.assertEqual(str(h()[1]), "20.0 pm")
        self.assertIsInstance(h(), collections.OrderedDict)

    def test_conversion_factors(self):
        for unit_type in SUPPORTED_UNITS:
            for unit1 in SUPPORTED_UNITS[unit_type]:
                a = Unit(1.0, unit1, unit_type)
                for unit2 in SUPPORTED_UNITS[unit_type]:
                    if unit1 == unit2:
                        self.assertEqual(a.to(unit2).to(unit1), 1.0)
                    else:
                        self.assert_almost_equal(a.to(unit2).to(unit1), 1.0, decimal=12)


class ArrayWithUnitTest(PymatgenTest):

    def test_energy(self):
        """
        Similar to UnitTest.test_energy.
        Check whether EnergyArray and Unit have same behavior.

        # TODO
        One can merge the two tests easily:

        for obj in [Energy, EnergyArray]:
            a = obj(...)
            self.assert(...)

        """
        a = EnergyArray(1.1, "eV")
        b = a.to("Ha")
        self.assertAlmostEqual(b, 0.0404242579378)
        c = EnergyArray(3.14, "J")
        self.assertAlmostEqual(c.to("eV"), 1.95983393276e+19)
        self.assertRaises(ValueError, Energy, 1, "m")

        d = EnergyArray(1, "Ha")
        self.assertAlmostEqual(a + d, 28.31138386)
        self.assertAlmostEqual(a - d, -26.11138386)
        self.assertEqual(a + 1, 2.1)

    def test_time(self):
        """
        Similar to UnitTest.test_time.
        Check whether EnergyArray and Unit have same behavior.
        """
        # here there's a minor difference because we have a ndarray with dtype=np.int.
        a = TimeArray(20, "h")
        self.assertAlmostEqual(a.to("s"), 3600 * 20)
        #Test left and right multiplication.
        self.assertEqual(str(a * 3), "60 h")
        self.assertEqual(str(3 * a), "60 h")

    def test_length(self):
        """
        Similar to UnitTest.test_time.
        Check whether EnergyArray and Unit have same behavior.
        """
        x = LengthArray(4.2, "ang")
        self.assertEqual(x.to("cm"), 4.2e-08)
        self.assertEqual(x.to("pm"), 420)
        self.assertEqual(str(x / 2), "2.1 ang")

    def test_array_algebra(self):
        ene_ha = EnergyArray([1, 2], "Ha")
        ene_ev = EnergyArray([1, 2], "eV")
        time_s = TimeArray([1, 2], "s")

        e1 = ene_ha.copy()
        e1 += 1
        e2 = ene_ha.copy()
        e2 -= 1
        e3 = ene_ha.copy()
        e3 /= 2
        e4 = ene_ha.copy()
        e4 *= 2

        objects_with_unit = [
            ene_ha + ene_ev,
            ene_ha - ene_ev,
            3 * ene_ha,
            ene_ha * 3,
            ene_ha / 3,
            ene_ha.copy(),
            ene_ha[0:1],
            e1,
            e2,
            e3,
            e4,
        ]

        for obj in objects_with_unit:
            self.assertTrue(obj.unit == "Ha")

        objects_without_unit = [
            ene_ha * time_s,
            ene_ha / ene_ev,
            #3 / ene_ha,
            #ene_ha // ene_ev,
            # Here we could return a Unit object but I prefer this since Unit extends float while we could have an int.
            #ene_ha[0],
        ]

        for obj in objects_without_unit:
            print(obj, type(obj))
            self.assertTrue(type(obj) == np.ndarray)

        with self.assertRaises(ValueError):
            ene_ha + time_s

if __name__ == '__main__':
    import unittest
    unittest.main()
