#!/usr/bin/env python

"""
This module implements a Unit, which is a subclass of float. It also defines
supported units for some commonly used units for energy, length, temperature,
time and charge. Units also support conversion to one another,
and additions and subtractions perform automatic conversion if units are
detected.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Aug 30, 2013"


import collections
import numbers
from functools import partial


"""
Definitions of supported units. Values below are essentially scaling and
conversion factors. What matters is the relative values, not the absolute.
"""
SUPPORTED_UNITS = {
    "energy": {
        "eV": 1,
        "Ha": 27.21138386,
        "Ry": 13.605698066,
        "J": 6.24150934e18
    },
    "length": {
        "ang": 1,
        "m": 1e10,
        "cm": 1e8,
        "pm": 1e-2,
        "bohr": 0.5291772083
    },
    "mass": {
        "kg": 1,
        "g": 1e-3,
        "amu": 1.660538921e-27
    },
    "temperature": {
        "K": 1
    },
    "time": {
        "s": 1,
        "min": 60,
        "h": 3600
    },
    "charge": {
        "C": 1,
        "e": 1.602176565e-19
    }
}


class Unit(float):
    """
    Subclasses float to attach a unit type. Typically, you should use the
    pre-defined unit type subclasses such as Energy, Length, etc. instead of
    using Unit directly.

    Supports conversion, addition and subtraction of the same unit type. E.g.,
    1 m + 20 cm will be automatically converted to 1.2 m (units follow the
    leftmost quantity).

    >>> e = Energy(1.1, "Ha")
    >>> a = Energy(1.1, "Ha")
    >>> b = Energy(3, "eV")
    >>> c = a + b
    >>> print c
    1.21024797619 Ha
    >>> c.to("eV")
    32.932522246 eV
    """

    def __new__(cls, val, unit, unit_type):
        return float.__new__(cls, val)

    def __init__(self, val, unit, unit_type):
        self.unit_type = unit_type
        if unit not in SUPPORTED_UNITS[unit_type]:
            raise ValueError(
                "{} is not a supported unit for {}".format(unit, unit_type))
        self.unit = unit
        super(Unit, self).__init__(val)

    def __repr__(self):
        s = super(Unit, self).__repr__()
        return "{} {}".format(s, self.unit)

    def __str__(self):
        return self.__repr__()

    def __add__(self, other):
        if not hasattr(other, "unit_type"):
            return super(Unit, self).__add__(other)
        if other.unit_type != self.unit_type:
            raise ValueError("Adding different types of units is not allowed")
        val = other
        if other.unit != self.unit:
            val = other.to(self.unit)
        return Unit(float(self) + val, unit_type=self.unit_type,
                    unit=self.unit)

    def __sub__(self, other):
        if not hasattr(other, "unit_type"):
            return super(Unit, self).__sub__(other)
        if other.unit_type != self.unit_type:
            raise ValueError("Subtracting different units is not allowed")
        val = other
        if other.unit != self.unit:
            val = other.to(self.unit)
        return Unit(float(self) - val, unit_type=self.unit_type,
                    unit=self.unit)

    def __mul__(self, other):
        val = super(Unit, self).__mul__(other)
        if isinstance(other, (float, int)):
            return Unit(val, unit_type=self.unit_type, unit=self.unit)
        return val

    def __rmul__(self, other):
        val = super(Unit, self).__rmul__(other)
        if isinstance(other, (float, int)):
            return Unit(val, unit_type=self.unit_type, unit=self.unit)
        return val

    def __div__(self, other):
        val = super(Unit, self).__div__(other)
        if isinstance(other, (float, int)):
            return Unit(val, unit_type=self.unit_type, unit=self.unit)
        return val

    def __truediv__(self, other):
        val = super(Unit, self).__truediv__(other)
        if isinstance(other, (float, int)):
            return Unit(val, unit_type=self.unit_type, unit=self.unit)
        return val

    def __neg__(self):
        return Unit(super(Unit, self).__neg__(), unit_type=self.unit_type,
                    unit=self.unit)

    def to(self, new_unit):
        """
        Conversion to a new_unit.

        Args:
            new_unit:
                New unit type.

        Returns:
            A Unit object in the new units.

        Example usage:
        >>> e = Energy(1.1, "eV")
        >>> e = Energy(1.1, "Ha")
        >>> e.to("eV")
        29.932522246 eV
        """
        if new_unit not in SUPPORTED_UNITS[self.unit_type]:
            raise ValueError(
                "{} is not a supported unit for {}".format(new_unit,
                                                           self.unit_type))
        conversion = SUPPORTED_UNITS[self.unit_type]
        return Unit(
            self / conversion[new_unit] * conversion[self.unit],
            unit_type=self.unit_type, unit=new_unit)

    @property
    def supported_units(self):
        """
        Supported units for specific unit type.
        """
        return SUPPORTED_UNITS[self.unit_type]


Energy = partial(Unit, unit_type="energy")

Length = partial(Unit, unit_type="length")

Mass = partial(Unit, unit_type="mass")

Temp = partial(Unit, unit_type="temperature")

Time = partial(Unit, unit_type="time")

Charge = partial(Unit, unit_type="charge")


def unitized(unit_type, unit):
    """
    Useful decorator to assign units to the output of a function. For
    sequences, all values in the sequences are assigned the same unit.

    Args:
        unit_type:
            Type of units (energy, length, mass, etc.).
        units:
            Specific units (eV, Ha, m, ang, etc.).
    """
    def wrap(f):
        def wrapped_f(*args, **kwargs):
            val = f(*args, **kwargs)
            if isinstance(val, collections.Sequence):
                # This complicated way is to ensure the sequence type is
                # preserved (list or tuple).
                return val.__class__([Unit(i, unit_type=unit_type,
                                           unit=unit) for i in val])
            elif isinstance(val, numbers.Number):
                return Unit(val, unit_type=unit_type, unit=unit)
            return val
        return wrapped_f
    return wrap


if __name__ == "__main__":
    import doctest
    doctest.testmod()