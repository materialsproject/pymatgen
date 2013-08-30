"""
Unit conversion tools. These functions are defined here instead of
core.physical_constants because we have to import and inspect the attributes
of the module.
"""
from __future__ import division

import numpy as np
from functools import partial


# The values below are essentially scaling and conversion factors. What matters
# is the relative values, not the absolute.

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
        "bohr": 0.5291772083
    },
    "mass": {
        "kg": 1,
        "amu": 1.660538921e-27
    },
    "temperature": {
        "K": 1
    },
    "time": {
        "s": 1,
        "min": 1 / 60
    },
    "charge": {
        "C": 1,
        "e": 1.602176565e-19
    }
}


class Unit(float):
    """
    Subclasses float to attach a unit type.
    """

    def __new__(cls, val, unit, unit_type):
        return float.__new__(cls, val)

    def __init__(self, val, unit, unit_type):
        self.val = val
        self.unit_type = unit_type
        if unit not in SUPPORTED_UNITS[unit_type]:
            raise ValueError(
                "{} is not a supported unit for {}".format(unit, unit_type))
        self.unit = unit
        super(Unit, self).__init__(val)

    def __repr__(self):
        return "{} {}".format(self.val, self.unit)

    def __str__(self):
        return "{} {}".format(self.val, self.unit)

    def __add__(self, other):
        if other.unit_type != self.unit_type:
            raise ValueError("Adding different types of units is not allowed")
        val = other.val
        if other.unit != self.unit:
            val = other.to(self.unit)
        return Unit(self.val + val, unit_type=self.unit_type,
                    unit=self.unit)

    def __sub__(self, other):
        if other.unit_type != self.unit_type:
            raise ValueError("Adding different types of units is not allowed")
        val = other.val
        if other.unit != self.unit:
            val = other.to(self.unit)
        return Unit(self.val - val, unit_type=self.unit_type,
                    unit=self.unit)

    def __neg__(self):
        return Unit(-self.val, unit_type=self.unit_type,
                    unit=self.unit)

    def to(self, new_unit):
        if new_unit not in SUPPORTED_UNITS[self.unit_type]:
            raise ValueError(
                "{} is not a supported unit for {}".format(new_unit,
                                                           self.unit_type))
        conversion = SUPPORTED_UNITS[self.unit_type]
        return Unit(
            self.val / conversion[new_unit] * conversion[self.unit],
            unit_type=self.unit_type, unit=new_unit)


Energy = partial(Unit, unit_type="energy")

Length = partial(Unit, unit_type="length")

Mass = partial(Unit, unit_type="mass")

Temp = partial(Unit, unit_type="temperature")

Time = partial(Unit, unit_type="time")
