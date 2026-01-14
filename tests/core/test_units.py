from __future__ import annotations

import pickle

import numpy as np
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
    obj_with_unit,
    unitized,
)
from pymatgen.util.testing import MatSciTest


def test_unit_conversions():
    assert Ha_to_eV == approx(27.211386245988)
    assert eV_to_Ha == 1 / Ha_to_eV
    assert Ry_to_eV == approx(Ha_to_eV / 2)
    assert bohr_to_angstrom == approx(0.529177210903)
    assert amu_to_kg == approx(1.66053906660e-27)


class TestUnit(MatSciTest):
    def test_init(self):
        u1 = Unit((("m", 1), ("s", -1)))
        assert str(u1) == "m s^-1"
        assert len(u1) == 2

        u2 = Unit("kg m ^ 2 s ^ -2")
        assert str(u2) == "J"
        assert str(u1 * u2) == "J m s^-1"
        assert str(u2 / u1) == "J s m^-1"
        assert str(u1 / Unit("m")) == "Hz"
        assert str(u1 * Unit("s")) == "m"

        acc = u1 / Unit("s")
        newton = Unit("kg") * acc
        assert str(newton * Unit("m")) == "N m"

    def test_as_base_units(self):
        # length
        u = Unit({"m": 1})  # meter
        base, factor = u.as_base_units
        assert base == {"m": 1}
        assert factor == 1.0

        u = Unit({"cm": 1})
        base, factor = u.as_base_units
        assert base == {"m": 1}
        assert pytest.approx(factor) == 0.01

        # force
        u = Unit({"N": 1})
        base, factor = u.as_base_units
        assert base == {"kg": 1, "m": 1, "s": -2}
        assert factor == 1.0

        # energy
        u = Unit({"J": 1})
        base, factor = u.as_base_units
        assert base == {"kg": 1, "m": 2, "s": -2}
        assert factor == 1.0

    def test_get_conversion_factor(self):
        # same unit
        u = Unit("m")
        factor = u.get_conversion_factor("m")
        assert factor == pytest.approx(1.0)

        # length
        u = Unit("cm")
        factor = u.get_conversion_factor("m")
        assert factor == pytest.approx(0.01)

        # mass
        u = Unit("g")
        factor = u.get_conversion_factor("kg")
        assert factor == pytest.approx(1e-3)

        # area
        u = Unit("cm^2")
        factor = u.get_conversion_factor("m^2")
        assert factor == pytest.approx(1e-4)

        # energy
        u = Unit("J")
        factor = u.get_conversion_factor("J")
        assert factor == pytest.approx(1.0)

        # incompatible conversion
        # TODO: `UnitError` in implementation doesn't work
        u = Unit("m")
        with pytest.raises(KeyError):
            u.get_conversion_factor("s")


class TestFloatWithUnit(MatSciTest):
    def test_rtruediv(self):
        # float / FloatWithUnit
        pressure = FloatWithUnit(5, "MPa")
        pressure_inv = 5 / pressure

        assert float(pressure_inv) == approx(1)
        assert pressure_inv.unit == Unit("MPa") ** -1
        assert isinstance(pressure_inv, FloatWithUnit)

        # FloatWithUnit / FloatWithUnit
        pres_pres = pressure / pressure
        assert float(pres_pres) == approx(1)
        assert isinstance(pres_pres, FloatWithUnit)

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
        assert a + 1 == approx(2.1)
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
        assert float(3 * a) == approx(60.0)
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

        with pytest.raises(ValueError, match="Unit is missing in string 1"):
            Memory.from_str("1")

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
        assert j_out[0] == approx(0.005)
        assert j_out[1] == approx(0.01)

        # Test None
        @unitized("eV")
        def func_none():
            return None

        assert func_none() is None

        # Test unsupported type
        @unitized("eV")
        def func_invalid():
            return object()

        with pytest.raises(TypeError, match="Don't know how to assign units"):
            func_invalid()

    def test_compound_operations(self):
        earth_acc = 9.81 * Length(1, "m") / (Time(1, "s") ** 2)
        e_pot = Mass(1, "kg") * earth_acc * Length(1, "m")
        assert str(e_pot) == "9.81 N m"

        form_e = FloatWithUnit(10, unit="kJ mol^-1").to("eV atom^-1")
        assert form_e == approx(0.103642691905)
        assert str(form_e.unit) == "eV atom^-1"
        with pytest.raises(UnitError, match="Units .* are not compatible"):
            form_e.to("m s^-1")
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

    def test_supported_units(self):
        energy_supported_units = Energy(1.1, "eV").supported_units
        assert energy_supported_units == ("eV", "meV", "Ha", "Ry", "J", "kJ", "kCal")

    def test_neg(self):
        x = FloatWithUnit(5, "MPa")
        assert FloatWithUnit(-5, "MPa") == -x

    def test_mul(self):
        # __mul__ with number
        energy = Energy(2, "eV")
        result = energy * 3
        assert isinstance(result, FloatWithUnit)
        assert pytest.approx(float(result)) == 6.0
        assert str(result.unit) == "eV"

        # __rmul__ with number
        energy = Energy(2, "eV")
        result = 3 * energy
        assert isinstance(result, FloatWithUnit)
        assert pytest.approx(float(result)) == 6.0
        assert str(result.unit) == "eV"

        # __mul__ with FloatWithUnit
        force = Energy(2, "eV") * Length(1, "ang") ** -1
        length = Length(3, "ang")
        result = force * length
        assert isinstance(result, FloatWithUnit)
        assert pytest.approx(float(result)) == 6.0
        assert str(result.unit) == "eV"

        # symmetry check for __rmul__ with FloatWithUnit
        result1 = force * length
        result2 = length * force
        assert pytest.approx(float(result1)) == float(result2)
        assert result1.unit == result2.unit

    def test_different_unit_type(self):
        energy_unit = Energy(1.1, "eV")
        length_unit = Length(4.2, "ang")

        # __add__
        with pytest.raises(UnitError, match="Adding different types of units"):
            energy_unit + length_unit

        # __sub__
        with pytest.raises(UnitError, match="Subtracting different units"):
            energy_unit - length_unit


class TestArrayWithUnit(MatSciTest):
    def test_neg(self):
        arr = LengthArray([1.0, 2.0, 3.0], "m")
        neg_arr = -arr

        assert isinstance(neg_arr, ArrayWithUnit)

        assert_array_equal(neg_arr, np.array([-1.0, -2.0, -3.0]))

        assert neg_arr.unit == arr.unit
        assert neg_arr.unit_type == arr.unit_type

    def test_mul(self):
        # __mul__ with scalar
        arr = LengthArray([1.0, 2.0, 3.0], "m")
        result = arr * 2
        assert isinstance(result, ArrayWithUnit)
        np.testing.assert_array_equal(result, np.array([2.0, 4.0, 6.0]))
        assert result.unit == arr.unit
        assert result.unit_type == arr.unit_type

        # __rmul__ with scalar
        result = 2 * arr
        assert isinstance(result, ArrayWithUnit)
        np.testing.assert_array_equal(result, np.array([2.0, 4.0, 6.0]))
        assert result.unit == arr.unit
        assert result.unit_type == arr.unit_type

        # __rmul__ with ArrayWithUnit
        length = LengthArray([1.0, 2.0], "m")
        time = TimeArray([3.0, 4.0], "s")

        result = length * time
        assert isinstance(result, ArrayWithUnit)
        np.testing.assert_array_equal(result, np.array([3.0, 8.0]))
        assert result.unit == length.unit * time.unit
        assert result.unit_type is None

        # symmetry check for __rmul__
        result2 = time * length
        np.testing.assert_array_equal(result, result2)
        assert result.unit == result2.unit

    def test_conversions(self):
        arr = LengthArray([1.0, 2.0, 3.0], "m")
        assert "[0.001 0.002 0.003] km" in arr.conversions()

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

        e2_in_ha = EnergyArray(1, "Ha")
        assert (e_in_ev + e2_in_ha) == approx(28.311386245987997)
        assert (e_in_ev - e2_in_ha) == approx(-26.111386245987994)
        assert float(e_in_ev + 1) == approx(2.1)

    def test_time(self):
        """Similar to FloatWithUnitTest.test_time.
        Check whether EnergyArray and FloatWithUnit have same behavior.
        """
        # here there's a minor difference because we have a ndarray with dtype=np.int64
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
            assert obj.unit in (
                Unit("Ha"),
                Unit("eV"),
                Unit("s^-1"),
                Unit("Ha s"),
                Unit("Ha eV^-1"),
            )

        # Here we could return a FloatWithUnit object but I prefer this
        # a bare scalar since FloatWithUnit extends float while we could
        # have an int.
        objects_without_unit = [ene_ha[0]]

        for obj in objects_without_unit:
            assert not hasattr(obj, "unit")

        with pytest.raises(UnitError, match="Adding different types of units is not allowed"):
            _ = ene_ha + time_s

        with pytest.raises(UnitError, match="Subtracting different units"):
            _ = ene_ha - time_s

    def test_factors(self):
        e_arr = EnergyArray([27.21138386, 1], "eV").to("Ha")
        assert str(e_arr).endswith("Ha")
        len_arr = LengthArray([1.0], "ang").to("bohr")
        assert str(len_arr).endswith(" bohr")
        v = ArrayWithUnit([1, 2, 3], "bohr^3").to("ang^3")
        assert str(v).endswith(" ang^3")

    def test_as_base_units(self):
        pressure_arr = ArrayWithUnit([5, 10], "MPa")
        assert_array_equal(ArrayWithUnit([5e6, 1e7], "Pa"), pressure_arr.as_base_units)

    def test_supported_units(self):
        pressure_arr = ArrayWithUnit([5, 10], "MPa", unit_type="pressure")
        assert pressure_arr.supported_units.keys() == {"Pa", "KPa", "MPa", "GPa"}

        unknown_unit_type_array = ArrayWithUnit([5, 10], "MPa")
        with pytest.raises(RuntimeError, match="Cannot get supported unit"):
            _ = unknown_unit_type_array.supported_units


class TestDataPersistence(MatSciTest):
    @pytest.mark.parametrize("cls", [FloatWithUnit, ArrayWithUnit])
    def test_pickle(self, cls):
        """Test whether FloatWithUnit and ArrayWithUnit support pickle."""
        a = cls(1, "eV")
        b = cls(10, "N bohr")
        objects = [a, b]

        dumped = pickle.dumps(objects)
        loaded = pickle.loads(dumped)

        # Type check
        assert all(isinstance(x, cls) for x in loaded)

        # Value and unit check
        for orig, new in zip(objects, loaded, strict=True):
            if isinstance(orig, ArrayWithUnit):
                assert_array_equal(orig, new)
            else:
                assert float(orig) == float(new)
            assert orig.unit == new.unit
            assert orig.unit_type == new.unit_type


class TestObjWithUnit:
    def test_scalar(self):
        e = obj_with_unit(2.5, "eV")
        assert isinstance(e, FloatWithUnit)
        assert float(e) == 2.5
        assert str(e.unit) == "eV"

    def test_array(self):
        arr = obj_with_unit([1.0, 2.0, 3.0], "m")
        assert isinstance(arr, ArrayWithUnit)
        np.testing.assert_array_equal(arr, np.array([1.0, 2.0, 3.0]))
        assert str(arr.unit) == "m"

    def test_mapping(self):
        d = obj_with_unit({"a": 1.0, "b": [2.0, 3.0]}, "Ha")
        assert isinstance(d, dict)
        assert isinstance(d["a"], FloatWithUnit)
        assert isinstance(d["b"], ArrayWithUnit)
        assert str(d["a"].unit) == "Ha"
        assert str(d["b"].unit) == "Ha"

    def test_nested_dict(self):
        d = obj_with_unit({"outer": {"inner": 1.0}}, "eV")
        assert isinstance(d["outer"], dict)
        assert isinstance(d["outer"]["inner"], FloatWithUnit)
        assert str(d["outer"]["inner"].unit) == "eV"
