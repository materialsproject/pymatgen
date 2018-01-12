# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import collections

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.units import (Energy, Time, Length, unitized, Mass, Memory,
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
        self.assertEqual(str(u1 / Unit("m")), "Hz")
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
        self.assertAlmostEqual(c.to("eV"), 1.959833865527343e+19)
        self.assertRaises(UnitError, Energy, 1, "m")

        d = Energy(1, "Ha")
        self.assertAlmostEqual(a + d, 28.31138602063284)
        self.assertAlmostEqual(a - d, -26.111386020632835)
        self.assertEqual(a + 1, 2.1)
        self.assertEqual(str(a / d), "1.1 eV Ha^-1")

    def test_time(self):
        a = Time(20, "h")
        self.assertAlmostEqual(float(a.to("s")), 3600 * 20)
        #Test left and right multiplication.
        b = a * 3
        self.assertAlmostEqual(float(b), 60.0)
        self.assertEqual(str(b.unit), "h")
        self.assertEqual(float(3 * a), 60.0)

    def test_length(self):
        x = Length(4.2, "ang")
        self.assertAlmostEqual(x.to("cm"), 4.2e-08)
        self.assertEqual(x.to("pm"), 420)
        self.assertEqual(str(x / 2), "2.1 ang")
        y = x ** 3
        self.assertAlmostEqual(y, 74.088)
        self.assertEqual(str(y.unit), "ang^3")

    def test_memory(self):
        mega = Memory(1, "Mb")
        self.assertEqual(mega.to("byte"), 1024**2)
        self.assertEqual(mega, Memory(1, "mb"))

        same_mega = Memory.from_string("1Mb")
        self.assertEqual(same_mega.unit_type, "memory")

        other_mega = Memory.from_string("+1.0 mb")
        self.assertEqual(mega, other_mega)

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

        @unitized("kg")
        def i():
            return FloatWithUnit(5, "g")

        self.assertEqual(i(), FloatWithUnit(0.005, "kg"))

        @unitized("kg")
        def j():
            return ArrayWithUnit([5, 10], "g")

        j_out = j()
        self.assertEqual(j_out.unit, Unit("kg"))
        self.assertEqual(j_out[0], 0.005)
        self.assertEqual(j_out[1], 0.01)

    def test_compound_operations(self):
        g = 10 * Length(1, "m") / (Time(1, "s") ** 2)
        e = Mass(1, "kg") * g * Length(1, "m")
        self.assertEqual(str(e), "10.0 N m")
        form_e = FloatWithUnit(10, unit="kJ mol^-1").to("eV atom^-1")
        self.assertAlmostEqual(float(form_e), 0.103642691905)
        self.assertEqual(str(form_e.unit), "eV atom^-1")
        self.assertRaises(UnitError, form_e.to, "m s^-1")
        a = FloatWithUnit(1.0, "Ha^3")
        b = a.to("J^3")
        self.assertAlmostEqual(b, 8.28672661615e-53)
        self.assertEqual(str(b.unit), "J^3")
        a = FloatWithUnit(1.0, "Ha bohr^-2")
        b = a.to("J m^-2")
        self.assertAlmostEqual(b, 1556.893078472351)
        self.assertEqual(str(b.unit), "J m^-2")

    def test_as_base_units(self):
        x = FloatWithUnit(5, "MPa")
        self.assertEqual(FloatWithUnit(5000000, "Pa"), x.as_base_units)


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
        self.assertAlmostEqual(float(b), 0.0404242579378)
        c = EnergyArray(3.14, "J")
        self.assertAlmostEqual(float(c.to("eV")), 1.959833865527343e+19, 5)
        # self.assertRaises(ValueError, Energy, 1, "m")

        d = EnergyArray(1, "Ha")
        self.assertAlmostEqual(float(a + d), 28.31138602063284)
        self.assertAlmostEqual(float(a - d), -26.111386020632835)
        self.assertEqual(float(a + 1), 2.1)

    def test_time(self):
        """
        Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        # here there's a minor difference because we have a ndarray with
        # dtype=np.int.
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
        self.assertAlmostEqual(float(x.to("cm")), 4.2e-08)
        self.assertEqual(float(x.to("pm")), 420)
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
        #e3 /= 2
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
            self.assertTrue(hasattr(obj, "unit"))

        objects_without_unit = [
            # Here we could return a FloatWithUnit object but I prefer this
            # a bare scalar since FloatWithUnit extends float while we could
            # have an int.
            ene_ha[0],
        ]

        for obj in objects_without_unit:
            self.assertFalse(hasattr(obj, "unit"))

        with self.assertRaises(UnitError):
            ene_ha + time_s

    def test_factors(self):
        e = EnergyArray([27.21138386, 1], "eV").to("Ha")
        self.assertTrue(str(e).endswith("Ha"))
        l = LengthArray([1.0], "ang").to("bohr")
        self.assertTrue(str(l).endswith(" bohr"))
        v = ArrayWithUnit([1, 2, 3], "bohr^3").to("ang^3")
        self.assertTrue(str(v).endswith(' ang^3'))

    def test_as_base_units(self):
        x = ArrayWithUnit([5, 10], "MPa")
        self.assertArrayEqual(ArrayWithUnit([5000000, 10000000], "Pa"), x.as_base_units)


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
