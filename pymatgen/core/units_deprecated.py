import numpy as np
import pymatgen.core.physical_constants as c


#def _nop(values):
#    return np.asanyarray(values)
#
#_any2Ha = {"Ha": _nop}
#for var in vars(c):
#    attr = c.__getattribute__(var)
#    if var.endswith("2Ha") and callable(attr):
#        _any2Ha[var[:-3]] = attr


#def any2Ha(units):
#    """
#    Returns a callable that converts energie(s) from units to Hartree
#
#    >>> any2Ha("Ha")([27.21138386, 1])
#    array([ 27.21138386,   1.        ])
#
#    >>> any2Ha("eV")([27.21138386, 1])
#    array([ 1.        ,  0.03674933])
#    """
#    return _any2Ha[units]
#
#_any2eV = {"eV": _nop}
#for var in vars(c):
#    attr = c.__getattribute__(var)
#    if var.endswith("2eV") and callable(attr):
#        _any2eV[var[:-3]] = attr
#
#
#def any2eV(units):
#    """
#    Returns a callable that converts energie(s) from units to Hartree
#
#    >>> any2eV("Ha")(any2Ha("eV")(1))
#    1.0
#    """
#    return _any2eV[units]


#_any2Bohr = {"Bohr": _nop}
#for var in vars(c):
#    attr = c.__getattribute__(var)
#    if var.endswith("2Bohr") and callable(attr):
#        _any2Bohr[var[:-5]] = attr
#
#
#    1.0
#    """
#    return _any2Ang[units]

