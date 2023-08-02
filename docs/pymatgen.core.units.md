---
layout: default
title: pymatgen.core.units.md
nav_exclude: true
---

# pymatgen.core.units module

This module implements a FloatWithUnit, which is a subclass of float. It
also defines supported units for some commonly used units for energy, length,
temperature, time and charge. FloatWithUnit also support conversion to one
another, and additions and subtractions perform automatic conversion if
units are detected. An ArrayWithUnit is also implemented, which is a subclass
of numpy’s ndarray with similar unit features.


### _class_ pymatgen.core.units.ArrayWithUnit(input_array, unit, unit_type=None)
Bases: `ndarray`

Subclasses numpy.ndarray to attach a unit type. Typically, you should
use the pre-defined unit type subclasses such as EnergyArray,
LengthArray, etc. instead of using ArrayWithFloatWithUnit directly.

Supports conversion, addition and subtraction of the same unit type. E.g.,
1 m + 20 cm will be automatically converted to 1.2 m (units follow the
leftmost quantity).

```python
>>> a = EnergyArray([1, 2], "Ha")
>>> b = EnergyArray([1, 2], "eV")
>>> c = a + b
>>> print(c)
[ 1.03674933  2.07349865] Ha
>>> c.to("eV")
array([ 28.21138386,  56.42276772]) eV
```

Override __new__.


#### Error()
alias of `UnitError`


#### _property_ as_base_units()
Returns this ArrayWithUnit in base SI units, including derived units.


* **Returns**

    An ArrayWithUnit object in base SI units



#### conversions()
Returns a string showing the available conversions.
Useful tool in interactive mode.


#### _property_ supported_units()
Supported units for specific unit type.


#### to(new_unit)
Conversion to a new_unit.


* **Parameters**

    **new_unit** – New unit type.



* **Returns**

    A ArrayWithFloatWithUnit object in the new units.


Example usage:
>>> e = EnergyArray([1, 1.1], “Ha”)
>>> e.to(“eV”)
array([ 27.21138386,  29.93252225]) eV


#### _property_ unit(_: st_ )
The unit, e.g., “eV”.


* **Type**

    return



#### _property_ unit_type(_: st_ )
The type of unit. Energy, Charge, etc.


* **Type**

    return



### pymatgen.core.units.Charge(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='charge'_ )
A float with a charge unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., C, e (electron charge). Must be valid unit or UnitError
    is raised.



### pymatgen.core.units.Energy(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='energy'_ )
A float with an energy unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., eV, kJ, etc. Must be valid unit or UnitError is raised.



### _class_ pymatgen.core.units.FloatWithUnit(val, unit, unit_type=None)
Bases: `float`

Subclasses float to attach a unit type. Typically, you should use the
pre-defined unit type subclasses such as Energy, Length, etc. instead of
using FloatWithUnit directly.

Supports conversion, addition and subtraction of the same unit type. E.g.,
1 m + 20 cm will be automatically converted to 1.2 m (units follow the
leftmost quantity). Note that FloatWithUnit does not override the eq
method for float, i.e., units are not checked when testing for equality.
The reason is to allow this class to be used transparently wherever floats
are expected.

```python
>>> e = Energy(1.1, "Ha")
>>> a = Energy(1.1, "Ha")
>>> b = Energy(3, "eV")
>>> c = a + b
>>> print(c)
1.2102479761938871 Ha
>>> c.to("eV")
32.932522246000005 eV
```

Initializes a float with unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – A unit. E.g., “C”.


    * **unit_type** (*str*) – A type of unit. E.g., “charge”



#### Error()
alias of `UnitError`


#### _property_ as_base_units()
Returns this FloatWithUnit in base SI units, including derived units.


* **Returns**

    A FloatWithUnit object in base SI units



#### _classmethod_ from_str(s)
Parse string to FloatWithUnit.

Example: Memory.from_str(“1. Mb”)


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead

Use from_str instead.


#### _property_ supported_units()
Supported units for specific unit type.


#### to(new_unit)
Conversion to a new_unit. Right now, only supports 1 to 1 mapping of
units of each type.


* **Parameters**

    **new_unit** – New unit type.



* **Returns**

    A FloatWithUnit object in the new units.


Example usage:
>>> e = Energy(1.1, “eV”)
>>> e = Energy(1.1, “Ha”)
>>> e.to(“eV”)
29.932522246 eV


#### _property_ unit(_: st_ )
The unit, e.g., “eV”.


* **Type**

    return



#### _property_ unit_type(_: st_ )
The type of unit. Energy, Charge, etc.


* **Type**

    return



### pymatgen.core.units.Length(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='length'_ )
A float with a length unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., m, ang, bohr, etc. Must be valid unit or UnitError is
    raised.



### pymatgen.core.units.Mass(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='mass'_ )
A float with a mass unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., amu, kg, etc. Must be valid unit or UnitError is
    raised.



### pymatgen.core.units.Memory(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='memory'_ )
A float with a memory unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., Kb, Mb, Gb, Tb. Must be valid unit or UnitError
    is raised.



### pymatgen.core.units.Temp(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='temperature'_ )
A float with a temperature unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., K. Only K (kelvin) is supported.



### pymatgen.core.units.Time(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='time'_ )
A float with a time unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., s, min, h. Must be valid unit or UnitError is
    raised.



### _class_ pymatgen.core.units.Unit(unit_def)
Bases: `Mapping`

Represents a unit, e.g., “m” for meters, etc. Supports compound units.
Only integer powers are supported for units.

Constructs a unit.


* **Parameters**

    **unit_def** – A definition for the unit. Either a mapping of unit to
    powers, e.g., {“m”: 2, “s”: -1} represents “m^2 s^-1”,
    or simply as a string “kg m^2 s^-1”. Note that the supported
    format uses “^” as the power operator and all units must be
    space-separated.



#### Error()
alias of `UnitError`


#### _property_ as_base_units()
Converts all units to base SI units, including derived units.


* **Returns**

    (base_units_dict, scaling factor). base_units_dict will not
    contain any constants, which are gathered in the scaling factor.



#### get_conversion_factor(new_unit)
Returns a conversion factor between this unit and a new unit.
Compound units are supported, but must have the same powers in each
unit type.


* **Parameters**

    **new_unit** – The new unit.



### _exception_ pymatgen.core.units.UnitError()
Bases: `BaseException`

Exception class for unit errors.


### pymatgen.core.units.kb(_ = 8.617333262e-0_ )
Definitions of supported units. Values below are essentially scaling and
conversion factors. What matters is the relative values, not the absolute.
The SI units must have factor 1.


### pymatgen.core.units.obj_with_unit(obj: Any, unit: str)
Returns a FloatWithUnit instance if obj is scalar, a dictionary of
objects with units if obj is a dict, else an instance of
ArrayWithFloatWithUnit.


* **Parameters**


    * **obj** (*Any*) – Object to be given a unit.


    * **unit** (*str*) – Specific units (eV, Ha, m, ang, etc.).



### pymatgen.core.units.unitized(unit)
Useful decorator to assign units to the output of a function. You can also
use it to standardize the output units of a function that already returns
a FloatWithUnit or ArrayWithUnit. For sequences, all values in the sequences
are assigned the same unit. It works with Python sequences only. The creation
of numpy arrays loses all unit information. For mapping types, the values
are assigned units.


* **Parameters**

    **unit** – Specific unit (eV, Ha, m, ang, etc.).


Example usage:

```default
@unitized(unit="kg")
def get_mass():
    return 123.45
```