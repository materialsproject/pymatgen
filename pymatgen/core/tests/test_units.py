# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import pytest

from pymatgen.core.units import (
    ArrayWithUnit,
    Energy,
    EnergyArray,
    FloatWithUnit,
    Length,
    LengthArray,
    Mass,
    Memory,
    Time,
    TimeArray,
    Unit,
    UnitError,
    unitized,
)
from pymatgen.util.testing import PymatgenTest


class UnitTest(PymatgenTest):
    def test_init(self):
        u1 = Unit((("m", 1), ("s", -1)))
        assert str(u1) == "m s^-1"
        u2 = Unit("kg m ^ 2 s ^ -2")
        assert str(u2) == "J"
        assert str(u1 * u2) == "J m s^-1"
        assert str(u2 / u1) == "J s m^-1"
        assert str(u1 / Unit("m")) == "Hz"
        assert str(u1 * Unit("s")) == "m"

        acc = u1 / Unit("s")
        newton = Unit("kg") * acc
        assert str(newton * Unit("m")) == "N m"


class FloatWithUnitTest(PymatgenTest):
    def test_energy(self):
        a = Energy(1.1, "eV")
        b = a.to("Ha")
        assert round(abs(b - 0.0404242579378), 7) == 0
        c = Energy(3.14, "J")
        assert round(abs(c.to("eV") - 1.9598338493806797e19), 7) == 0
        with pytest.raises(UnitError):
            Energy(1, "m")

        d = Energy(1, "Ha")
        assert round(abs(a + d - 28.311386245987997), 7) == 0
        assert round(abs(a - d - -26.111386245987994), 7) == 0
        assert a + 1 == 2.1
        assert str(a / d) == "1.1 eV Ha^-1"

        e = Energy(1, "kJ")
        f = e.to("kCal")
        assert round(abs(f - 0.2390057361376673), 7) == 0
        assert str(e + f) == "2.0 kJ"
        assert str(f + e) == "0.4780114722753346 kCal"

    def test_time(self):
        a = Time(20, "h")
        assert round(abs(float(a.to("s")) - 3600 * 20), 7) == 0
        # Test left and right multiplication.
        b = a * 3
        assert round(abs(float(b) - 60.0), 7) == 0
        assert str(b.unit) == "h"
        assert float(3 * a) == 60.0
        a = Time(0.5, "d")
        assert round(abs(float(a.to("s")) - 3600 * 24 * 0.5), 7) == 0

    def test_length(self):
        x = Length(4.2, "ang")
        assert round(abs(x.to("cm") - 4.2e-08), 7) == 0
        assert x.to("pm") == 420
        assert str(x / 2) == "2.1 ang"
        y = x**3
        assert round(abs(y - 74.088), 7) == 0
        assert str(y.unit) == "ang^3"

    def test_memory(self):
        mega = Memory(1, "Mb")
        assert mega.to("byte") == 1024**2
        assert mega == Memory(1, "mb")

        same_mega = Memory.from_string("1Mb")
        assert same_mega.unit_type == "memory"

        other_mega = Memory.from_string("+1.0 mb")
        assert mega == other_mega

    def test_unitized(self):
        @unitized("eV")
        def f():
            return [1, 2, 3]

        assert str(f()[0]) == "1.0 eV"
        assert isinstance(f(), list)

        @unitized("eV")
        def g():
            return 2, 3, 4

        assert str(g()[0]) == "2.0 eV"
        assert isinstance(g(), tuple)

        @unitized("pm")
        def h():
            d = {}
            for i in range(3):
                d[i] = i * 20
            return d

        assert str(h()[1]) == "20.0 pm"
        assert isinstance(h(), dict)

        @unitized("kg")
        def i():
            return FloatWithUnit(5, "g")

        assert i() == FloatWithUnit(0.005, "kg")

        @unitized("kg")
        def j():
            return ArrayWithUnit([5, 10], "g")

        j_out = j()
        assert j_out.unit == Unit("kg")
        assert j_out[0] == 0.005
        assert j_out[1] == 0.01

    def test_compound_operations(self):
        g = 10 * Length(1, "m") / (Time(1, "s") ** 2)
        e = Mass(1, "kg") * g * Length(1, "m")
        assert str(e) == "10.0 N m"
        form_e = FloatWithUnit(10, unit="kJ mol^-1").to("eV atom^-1")
        assert round(abs(float(form_e) - 0.103642691905), 7) == 0
        assert str(form_e.unit) == "eV atom^-1"
        with pytest.raises(UnitError):
            form_e.to("m s^-1")
        a = FloatWithUnit(1.0, "Ha^3")
        b = a.to("J^3")
        assert round(abs(b - 8.28672661615e-53), 7) == 0
        assert str(b.unit) == "J^3"
        a = FloatWithUnit(1.0, "Ha bohr^-2")
        b = a.to("J m^-2")
        assert round(abs(b - 1556.8931028218924), 7) == 0
        assert str(b.unit) == "J m^-2"

    def test_as_base_units(self):
        x = FloatWithUnit(5, "MPa")
        assert FloatWithUnit(5000000, "Pa") == x.as_base_units


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
        assert round(abs(float(b) - 0.0404242579378), 7) == 0
        c = EnergyArray(3.14, "J")
        assert round(abs(float(c.to("eV")) - 1.9598338493806797e19), 5) == 0
        # self.assertRaises(ValueError, Energy, 1, "m")

        d = EnergyArray(1, "Ha")
        assert round(abs(float(a + d) - 28.311386245987997), 7) == 0
        assert round(abs(float(a - d) - -26.111386245987994), 7) == 0
        assert float(a + 1) == 2.1

    def test_time(self):
        """
        Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        # here there's a minor difference because we have a ndarray with
        # dtype=int.
        a = TimeArray(20, "h")
        assert pytest.approx(a.to("s")) == 3600 * 20
        # Test left and right multiplication.
        assert str(a * 3) == "60 h"
        assert str(3 * a) == "60 h"

    def test_length(self):
        """
        Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        x = LengthArray(4.2, "ang")
        assert round(abs(float(x.to("cm")) - 4.2e-08), 7) == 0
        assert float(x.to("pm")) == 420
        assert str(x / 2) == "2.1 ang"

    def test_array_algebra(self):
        ene_ha = EnergyArray([1, 2], "Ha")
        ene_ev = EnergyArray([1, 2], "eV")
        time_s = TimeArray([1, 2], "s")

        e1 = ene_ha.copy()
        e1 += 1
        e2 = ene_ha.copy()
        e2 -= 1
        e3 = ene_ha.copy()
        # e3 /= 2
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

        for obj in objects_with_unit:
            assert hasattr(obj, "unit")

        objects_without_unit = [
            # Here we could return a FloatWithUnit object but I prefer this
            # a bare scalar since FloatWithUnit extends float while we could
            # have an int.
            ene_ha[0],
        ]

        for obj in objects_without_unit:
            assert not hasattr(obj, "unit")

        with pytest.raises(UnitError):
            ene_ha + time_s

    def test_factors(self):
        e = EnergyArray([27.21138386, 1], "eV").to("Ha")
        assert str(e).endswith("Ha")
        l = LengthArray([1.0], "ang").to("bohr")
        assert str(l).endswith(" bohr")
        v = ArrayWithUnit([1, 2, 3], "bohr^3").to("ang^3")
        assert str(v).endswith(" ang^3")

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
                    assert str(old_item) == str(new_item)


if __name__ == "__main__":
    import unittest

    unittest.main()
