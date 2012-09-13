#!/usr/bin/env python

"""
This module provides core classes needed by all define electronic structure,
such as the Spin, Orbital, etc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"


from pymatgen.util.decorators import cached_class


class Spin(object):
    """
    Enum type for Spin.  Only up and down.
    Usage: Spin.up, Spin.down.
    """

    @cached_class
    class _SpinImpl(object):
        """
        Internal representation for Spin. Not to be instantiated directly.
        Use Spin enum types. Class is implemented as a cached class for
        memory efficiency.
        """
        def __init__(self, name):
            self._name = name

        def __int__(self):
            return 1 if self._name == "up" else -1

        def __repr__(self):
            return self._name

        def __eq__(self, other):
            if other == None:
                return False
            return self._name == other._name

        def __hash__(self):
            return self.__int__()

        def __str__(self):
            return self._name

    @staticmethod
    def from_int(i):
        """
        Provides the spin from an int. +1 == Spin.up, -1 == Spin.down.

        Args:
            i:
                integer representing direction of spin.
        """
        if i == 1:
            return Spin.up
        elif i == -1:
            return Spin.down
        else:
            raise ValueError("Spin integers must be 1 or -1")

    up = _SpinImpl("up")
    down = _SpinImpl("down")
    all_spins = (up, down)


class Orbital(object):
    """
    Enum type for OrbitalType. Indices are basically the azimutal quantum
    number, l.
    """

    @cached_class
    class _OrbitalImpl(object):
        """
        Internal representation of an orbital.  Do not use directly.
        Use the Orbital class enum types.  Class is implemented as a cached
        class for memory efficiency.
        """

        def __init__(self, name, vasp_index):
            self._name = name
            self._vasp_index = vasp_index

        def __int__(self):
            return self._vasp_index

        def __repr__(self):
            return self._name

        def __eq__(self, other):
            if other == None:
                return False
            return self._name == other._name

        def __hash__(self):
            return self.__int__()

        @property
        def orbital_type(self):
            return self._name[0].upper()

        def __str__(self):
            return self._name

    s = _OrbitalImpl("s", 0)
    py = _OrbitalImpl("py", 1)
    pz = _OrbitalImpl("pz", 2)
    px = _OrbitalImpl("px", 3)
    dxy = _OrbitalImpl("dxy", 4)
    dyz = _OrbitalImpl("dyz", 5)
    dz2 = _OrbitalImpl("dz2", 6)
    dxz = _OrbitalImpl("dxz", 7)
    dx2 = _OrbitalImpl("dx2", 8)
    f_3 = _OrbitalImpl("f_3", 9)
    f_2 = _OrbitalImpl("f_2", 10)
    f_1 = _OrbitalImpl("f_1", 11)
    f0 = _OrbitalImpl("f0", 12)
    f1 = _OrbitalImpl("f1", 13)
    f2 = _OrbitalImpl("f2", 14)
    f3 = _OrbitalImpl("f3", 15)

    all_orbitals = (s,
                    py, pz, px,
                    dxy, dyz, dz2, dxz, dx2,
                    f_3, f_2, f_1, f0, f1, f2, f3)

    @staticmethod
    def from_vasp_index(i):
        for orb in Orbital.all_orbitals:
            if int(orb) == i:
                return orb
        raise IndexError("Illegal exceeds supported orbital set")

    @staticmethod
    def from_string(orb_str):
        for orb in Orbital.all_orbitals:
            if str(orb) == orb_str:
                return orb
        raise ValueError("Illegal orbital definition!")
