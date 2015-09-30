# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides core classes needed by all define electronic structure,
such as the Spin, Orbital, etc.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import collections


class Spin(object):
    """
    Enum type for Spin.  Only up and down.
    Usage: Spin.up, Spin.down.
    """
    up, down = (1, -1)

    @staticmethod
    def from_int(i):
        """
        Provides the spin from an int. +1 == Spin.up, -1 == Spin.down.

        Args:
            i (1/-1): integer representing direction of spin.
        """
        if i == 1:
            return Spin.up
        elif i == -1:
            return Spin.down
        else:
            raise ValueError("Spin integers must be 1 or -1")

    all_spins = (up, down)


class _OrbitalImpl(collections.namedtuple('_Orbital', 'name vasp_index')):
    """
    Internal representation of an orbital.  Do not use directly.
    Use the Orbital class enum types.
    """
    __slots__ = ()

    def __int__(self):
        return self.vasp_index

    def __eq__(self, other):
        if other is None:
            return False
        return self.vasp_index == other.vasp_index

    def __hash__(self):
        return self.vasp_index

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    @property
    def orbital_type(self):
        """
        String indicating the type of orbital. Is always uppercase. E.g.,
        S, P, D, F, etc.
        """
        return self.name[0].upper()


class Orbital(object):
    """
    Enum type for OrbitalType. Indices are basically the azimutal quantum
    number, l.
    """

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
        """
        Returns an orbital based on the index of the orbital in VASP runs.
        """
        return Orbital.all_orbitals[i]

    @staticmethod
    def from_string(orb_str):
        """
        Returns an orbital from a string representation, e.g., "s", "px".
        """
        for orb in Orbital.all_orbitals:
            if str(orb) == orb_str:
                return orb
        raise ValueError("Illegal orbital definition!")
