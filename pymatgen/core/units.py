"""
Unit conversion tools. These functions are defined here instead of
core.physical_constants because we have to import and inspect the attributes
of the module.
"""
from __future__ import division

import numpy as np
import pymatgen.core.physical_constants as c
from functools import partial

#The values

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
        "u": 1.660538921e-27
    },
    "temperature": {
        "K": 1
    },
    "time": {
        "s": 1,
        "min": 1 / 60
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


def _nop(values):
    return np.asanyarray(values)

_any2Ha = {"Ha": _nop}
for var in vars(c):
    attr = c.__getattribute__(var)
    if var.endswith("2Ha") and callable(attr):
        _any2Ha[var[:-3]] = attr


def any2Ha(units):
    """
    Returns a callable that converts energie(s) from units to Hartree

    >>> any2Ha("Ha")([27.21138386, 1])
    array([ 27.21138386,   1.        ])

    >>> any2Ha("eV")([27.21138386, 1])
    array([ 1.        ,  0.03674933])
    """
    return _any2Ha[units]

_any2eV = {"eV": _nop}
for var in vars(c):
    attr = c.__getattribute__(var)
    if var.endswith("2eV") and callable(attr):
        _any2eV[var[:-3]] = attr


def any2eV(units):
    """
    Returns a callable that converts energie(s) from units to Hartree

    >>> any2eV("Ha")(any2Ha("eV")(1))
    1.0
    """
    return _any2eV[units]


_any2Bohr = {"Bohr": _nop}
for var in vars(c):
    attr = c.__getattribute__(var)
    if var.endswith("2Bohr") and callable(attr):
        _any2Bohr[var[:-5]] = attr


def any2Bohr(units):
    """
    Returns a callable that converts length(s) from units to Bohr.

    >>> any2Bohr("Bohr")(5)
    array(5)

    >>> any2Bohr("Ang")(1.0)
    1.8897261328856432
    """
    return _any2Bohr[units]


_any2Ang = {"Ang": _nop}
for var in vars(c):
    attr = c.__getattribute__(var)
    if var.endswith("2Ang") and callable(attr):
        _any2Ang[var[:-4]] = attr


def any2Ang(units):
    """
    Returns a callable that converts length(s) from units to Angstrom

    >>> any2Ang("Bohr")(any2Bohr("Ang")(1))
    1.0
    """
    return _any2Ang[units]


def any2Ang3(len_units):
    """
    Returns a callable that converts volume(s) given in unit length len_units to Angstrom**3.

    >>> any2Ang3("Ang")([1.,2.,3.])
    array([ 1.,  2.,  3.])

    >>> any2Ang3("Bohr")([1,2,3])
    array([ 0.14818471,  0.29636942,  0.44455413])
    """
    def func(values):
        func._len_units = len_units[:]
        values = np.asarray(values)
        if func._len_units == "Ang": 
            return values
        else:
            f = any2Ang(func._len_units)
            return f(values**(1/3.)) ** 3
    return func

