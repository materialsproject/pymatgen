#!/usr/bin/python
from __future__ import division

import numpy as np
import tempfile
import cPickle as pickle
import collections

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.units import (Energy, Time, Length, unitized, Mass,
                                 EnergyArray, TimeArray, LengthArray, Unit,
                                 FloatWithUnit, ArrayWithUnit, UnitError)


class UnitTest(PymatgenTest):

    def test_init(self):
        u1 = Unit((("m", 1), ("s", -1)))
        self.assertEqual(str(u1), "m s^-1")
        u2 = Unit("kg m ^ 2 s ^ -2")
        self.assertEqual(str(u2), "J")
        self.assertEqual(str(u1 * u2), "J m s^-1")
        self.assertEqual(str(u2 / u1), "J s m^-1")
        self.assertEqual(str(u1 / Unit("m")), "s^-1")
        self.assertEqual(str(u1 * Unit("s")), "m")

        acc = u1 / Unit("s")
        newton = Unit("kg") * acc
        self.assertEqual(str(newton * Unit("m")), "N m")

class FloatWithUnitTest(PymatgenTest):

    def test_energy(self):
        a = Energy(1.1, "eV")
        b = a.to("Ha")
        self.assertAlmostEqual(b, 0.0404242579378)
        c = Energy(3.14, "J")
        self.assertAlmostEqual(c.to("eV"), 1.9598339337836966e+19)
        self.assertRaises(UnitError, Energy, 1, "m")

        d = Energy(1, "Ha")
        self.assertAlmostEqual(a + d, 28.31138386)
        self.assertAlmostEqual(a - d, -26.11138386)
        self.assertEqual(a + 1, 2.1)
        self.assertEqual(str(a / d), "1.1 eV Ha^-1")

    def test_time(self):
        a = Time(20, "h")
        self.assertAlmostEqual(a.to("s"), 3600 * 20)
        #Test left and right multiplication.
        self.assertEqual(str(a * 3), "60.0 h")
        self.assertEqual(str(3 * a), "60.0 h")

    def test_length(self):
        x = Length(4.2, "ang")
        self.assertAlmostEqual(x.to("cm"), 4.2e-08)
        self.assertEqual(x.to("pm"), 420)
        self.assertEqual(str(x / 2), "2.1 ang")
        self.assertEqual(str(x ** 3), "74.088 ang^3")

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

    def test_compound_operations(self):
        g = 10 * Length(1, "m") / (Time(1, "s") ** 2)
        e = Mass(1, "kg") * g * Length(1, "m")
        self.assertEqual(str(e), "10.0 N m")
        form_e = FloatWithUnit(10, unit="kJ mol^-1")
        self.assertEqual(str(form_e.to("eV atom^-1")), "0.103642691905 eV atom^-1")
        self.assertRaises(UnitError, form_e.to, "m s^-1")
        a = FloatWithUnit(1.0, "Ha^3")
        self.assertEqual(str(a.to("J^3")), "8.28672661615e-53 J^3")
        a = FloatWithUnit(1.0, "Ha bohr^-2")
        self.assertEqual(str(a.to("J m^-2")), "1556.89291457 J m^-2")


class ArrayWithFloatWithUnitTest(PymatgenTest):

    def test_energy(self):
        """
        Similar to FloatWithUnitTest.test_energy.
        Check whether EnergyArray and FloatWithUnit have same behavior.

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
        self.assertAlmostEqual(c.to("eV"), 1.9598339337836966e+19)
        # self.assertRaises(ValueError, Energy, 1, "m")

        d = EnergyArray(1, "Ha")
        self.assertAlmostEqual(a + d, 28.31138386)
        self.assertAlmostEqual(a - d, -26.11138386)
        self.assertEqual(a + 1, 2.1)

    def test_time(self):
        """
        Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        # here there's a minor difference because we have a ndarray with dtype=np.int.
        a = TimeArray(20, "h")
        self.assertAlmostEqual(a.to("s"), 3600 * 20)
        #Test left and right multiplication.
        self.assertEqual(str(a * 3), "60 h")
        self.assertEqual(str(3 * a), "60 h")

    def test_length(self):
        """
        Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        x = LengthArray(4.2, "ang")
        self.assertAlmostEqual(x.to("cm"), 4.2e-08)
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
            3 / ene_ha,
            ene_ha * time_s,
            ene_ha / ene_ev,
            ene_ha.copy(),
            ene_ha[0:1],
            e1,
            e2,
            e3,
            e4,
        ]

        for i, obj in enumerate(objects_with_unit):
            #print(i, obj.unit)
            self.assertTrue(hasattr(obj, "unit"))
            #self.assertTrue(str(obj.unit) == "Ha")

        objects_without_unit = [
            # Here we could return a FloatWithUnit object but I prefer this
            # a bare scalar since FloatWithUnit extends float while we could have an int.
            ene_ha[0],
        ]

        for obj in objects_without_unit:
            self.assertFalse(hasattr(obj, "unit"))

        with self.assertRaises(UnitError):
            ene_ha + time_s

    def test_factors(self):
        e = EnergyArray([27.21138386, 1], "eV").to("Ha")
        self.assertTrue(str(e) == "[ 1.          0.03674933] Ha")
        l = LengthArray([1.0], "ang").to("bohr")
        self.assertTrue(str(l) == "[ 1.88972613] bohr")
        v = ArrayWithUnit([1,2,3], "bohr^3").to("ang^3")
        self.assertTrue(str(v) == '[ 0.14818471  0.29636942  0.44455413] '
                                  'ang^3')


class DataPersistenceTest(PymatgenTest):
    def test_pickle(self):
        """Test whether FloatWithUnit and ArrayWithUnit support pickle"""

        for cls in [FloatWithUnit, ArrayWithUnit]:
            a = cls(1, "eV")
            b = cls(10, "N bohr")
            objects = [a, b]

            new_objects_from_protocol = self.serialize_with_pickle(objects)

            for new_objects in new_objects_from_protocol:
                for old_item, new_item in zip(objects, new_objects):
                    self.assertTrue(str(old_item) == str(new_item))


if __name__ == '__main__':
    import unittest
    unittest.main()
