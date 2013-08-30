import numpy as np
import pymatgen.core.physical_constants as c


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