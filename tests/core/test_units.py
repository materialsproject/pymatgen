from __future__ import annotations

import warnings

import pytest
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.core.units import (
    ArrayWithUnit,
    Energy,
    EnergyArray,
    FloatWithUnit,
    Ha_to_eV,
    Length,
    LengthArray,
    Mass,
    Memory,
    Ry_to_eV,
    Time,
    TimeArray,
    Unit,
    UnitError,
    amu_to_kg,
    bohr_to_angstrom,
    eV_to_Ha,
    unitized,
)
from pymatgen.util.testing import PymatgenTest


def test_unit_conversions():
    assert Ha_to_eV == approx(27.211386245988)
    assert eV_to_Ha == 1 / Ha_to_eV
    assert Ry_to_eV == approx(Ha_to_eV / 2)
    assert bohr_to_angstrom == approx(0.529177210903)
    assert amu_to_kg == approx(1.66053906660e-27)


class TestUnit(PymatgenTest):
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


class TestFloatWithUnit(PymatgenTest):
    def test_energy(self):
        a = Energy(1.1, "eV")
        b = a.to("Ha")
        assert b == approx(0.0404242579378)
        c = Energy(3.14, "J")
        assert c.to("eV") == approx(1.9598338493806797e19)
        with pytest.raises(UnitError, match="m is not a supported unit for energy"):
            Energy(1, "m")

        d = Energy(1, "Ha")
        assert a + d == approx(28.311386245987997)
        assert a - d == approx(-26.111386245987994)
        assert a + 1 == 2.1
        assert str(a / d) == "1.1 eV Ha^-1"

        e_kj = Energy(1, "kJ")
        e_kcal = e_kj.to("kCal")
        assert e_kcal == approx(0.2390057361376673)
        assert str(e_kj + e_kcal) == "2.0 kJ"
        assert str(e_kcal + e_kj) == "0.4780114722753346 kCal"

    def test_time(self):
        a = Time(20, "h")
        assert float(a.to("s")) == approx(3600 * 20)
        # Test left and right multiplication.
        b = a * 3
        assert float(b) == approx(60)
        assert str(b.unit) == "h"
        assert float(3 * a) == 60.0
        a = Time(0.5, "d")
        assert float(a.to("s")) == approx(3600 * 24 * 0.5)

    def test_length(self):
        x = Length(4.2, "ang")
        assert x.to("cm") == approx(4.2e-08)
        assert x.to("pm") == 420
        assert str(x / 2) == "2.1 ang"

        y = x**3
        assert y == approx(74.088)
        assert str(y.unit) == "ang^3"

    def test_memory(self):
        mega_0 = Memory(1, "MB")
        assert mega_0.to("byte") == 1024**2

        mega_1 = Memory.from_str("1 MB")
        assert mega_1.unit_type == "memory"

        mega_2 = Memory.from_str("1MB")
        assert mega_2.unit_type == "memory"

        mega_3 = Memory.from_str("+1.0 MB")
        assert mega_0 == mega_1 == mega_2 == mega_3

    def test_deprecated_memory(self):
        # TODO: remove after 2025-01-01
        for unit in ("Kb", "kb", "Mb", "mb", "Gb", "gb", "Tb", "tb"):
            with pytest.warns(DeprecationWarning, match=f"Unit {unit} is deprecated"):
                Memory(1, unit)

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            for unit in ("KB", "MB", "GB", "TB"):
                Memory(1, unit)

    def test_unitized(self):
        @unitized("eV")
        def func1():
            return [1, 2, 3]

        assert str(func1()[0]) == "1.0 eV"
        assert isinstance(func1(), list)

        @unitized("eV")
        def func2():
            return 2, 3, 4

        assert str(func2()[0]) == "2.0 eV"
        assert isinstance(func2(), tuple)

        @unitized("pm")
        def func3():
            return {idx: idx * 20 for idx in range(3)}

        assert str(func3()[1]) == "20.0 pm"
        assert isinstance(func3(), dict)

        @unitized("kg")
        def func4():
            return FloatWithUnit(5, "g")

        assert func4() == FloatWithUnit(0.005, "kg")

        @unitized("kg")
        def func5():
            return ArrayWithUnit([5, 10], "g")

        j_out = func5()
        assert j_out.unit == Unit("kg")
        assert j_out[0] == 0.005
        assert j_out[1] == 0.01

    def test_compound_operations(self):
        earth_acc = 9.81 * Length(1, "m") / (Time(1, "s") ** 2)
        e_pot = Mass(1, "kg") * earth_acc * Length(1, "m")
        assert str(e_pot) == "9.81 N m"
        form_e = FloatWithUnit(10, unit="kJ mol^-1").to("eV atom^-1")
        assert form_e == approx(0.103642691905)
        assert str(form_e.unit) == "eV atom^-1"
        with pytest.raises(UnitError) as exc:
            form_e.to("m s^-1")
        assert "Units ('mol', -1) and ('m', 1) are not compatible" in str(exc.value)
        a = FloatWithUnit(1.0, "Ha^3")
        b = a.to("J^3")
        assert b == approx(8.28672661615e-53)
        assert str(b.unit) == "J^3"
        a = FloatWithUnit(1.0, "Ha bohr^-2")
        b = a.to("J m^-2")
        assert b == approx(1556.8931028218924)
        assert str(b.unit) == "J m^-2"

    def test_as_base_units(self):
        x = FloatWithUnit(5, "MPa")
        assert FloatWithUnit(5000000, "Pa") == x.as_base_units

    def test_neg(self):
        x = FloatWithUnit(5, "MPa")
        assert FloatWithUnit(-5, "MPa") == -x


class TestArrayWithUnit(PymatgenTest):
    def test_energy(self):
        """Similar to TestFloatWithUnit.test_energy.
        Check whether EnergyArray and FloatWithUnit have same behavior.

        # TODO
        One can merge the two tests easily:

        for obj in [Energy, EnergyArray]:
            a = obj(...)
            self.assert(...)
        """
        e_in_ev = EnergyArray(1.1, "eV")
        e_in_ha = e_in_ev.to("Ha")
        assert e_in_ha == approx(0.0404242579378)
        e_in_j = EnergyArray(3.14, "J")
        assert (e_in_j.to("eV")) == approx(1.9598338493806797e19)
        # self.assertRaises(ValueError, Energy, 1, "m")

        e2_in_ha = EnergyArray(1, "Ha")
        assert (e_in_ev + e2_in_ha) == approx(28.311386245987997)
        assert (e_in_ev - e2_in_ha) == approx(-26.111386245987994)
        assert float(e_in_ev + 1) == 2.1

    def test_time(self):
        """Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        # here there's a minor difference because we have a ndarray with dtype=int
        a = TimeArray(20, "h")
        assert a.to("s") == 3600 * 20
        # Test left and right multiplication.
        assert str(a * 3) == "60 h"
        assert str(3 * a) == "60 h"

    def test_length(self):
        """Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        x = LengthArray(4.2, "ang")
        assert float(x.to("cm")) == approx(4.2e-08)
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
            ene_ha[:1],
            e1,
            e2,
            e3,
            e4,
        ]

        for obj in objects_with_unit:
            assert obj.unit in (Unit("Ha"), Unit("eV"), Unit("s^-1"), Unit("Ha s"), Unit("Ha eV^-1"))

        # Here we could return a FloatWithUnit object but I prefer this
        # a bare scalar since FloatWithUnit extends float while we could
        # have an int.
        objects_without_unit = [ene_ha[0]]

        for obj in objects_without_unit:
            assert not hasattr(obj, "unit")

        with pytest.raises(UnitError, match="Adding different types of units is not allowed"):
            _ = ene_ha + time_s

    def test_factors(self):
        e_arr = EnergyArray([27.21138386, 1], "eV").to("Ha")
        assert str(e_arr).endswith("Ha")
        len_arr = LengthArray([1.0], "ang").to("bohr")
        assert str(len_arr).endswith(" bohr")
        v = ArrayWithUnit([1, 2, 3], "bohr^3").to("ang^3")
        assert str(v).endswith(" ang^3")

    def test_as_base_units(self):
        pressure_arr = ArrayWithUnit([5, 10], "MPa")
        assert_array_equal(ArrayWithUnit([5000000, 10000000], "Pa"), pressure_arr.as_base_units)


class TestDataPersistence(PymatgenTest):
    def test_pickle(self):
        """Test whether FloatWithUnit and ArrayWithUnit support pickle."""
        for cls in [FloatWithUnit, ArrayWithUnit]:
            a = cls(1, "eV")
            b = cls(10, "N bohr")
            objects = [a, b]

            self.serialize_with_pickle(objects)
