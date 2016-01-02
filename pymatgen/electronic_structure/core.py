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

from enum import Enum


class Spin(Enum):
    """
    Enum type for Spin.  Only up and down.
    Usage: Spin.up, Spin.down.
    """
    up, down = (1, -1)

    def __int__(self):
        return self.value

    def __str__(self):
        return str(self.value)

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


class Orbital(Enum):
    """
    Enum type for OrbitalType. Indices are basically the azimutal quantum
    number, l.
    """

    s = "s", 0
    py = "py", 1
    pz = "pz", 2
    px = "px", 3
    dxy = "dxy", 4
    dyz = "dyz", 5
    dz2 = "dz2", 6
    dxz = "dxz", 7
    dx2 = "dx2", 8
    f_3 = "f_3", 9
    f_2 = "f_2", 10
    f_1 = "f_1", 11
    f0 = "f0", 12
    f1 = "f1", 13
    f2 = "f2", 14
    f3 = "f3", 15

    @property
    def vasp_index(self):
        return self.value[1]

    def __int__(self):
        return self.value[1]

    def __eq__(self, other):
        if other is None:
            return False
        return self.vasp_index == other.vasp_index

    def __hash__(self):
        return self.vasp_index

    def __repr__(self):
        return self.value[0]

    def __str__(self):
        return self.value[0]

    @property
    def orbital_type(self):
        """
        String indicating the type of orbital. Is always uppercase. E.g.,
        S, P, D, F, etc.
        """
        return self.value[0][0].upper()

    @staticmethod
    def from_vasp_index(i):
        """
        Returns an orbital based on the index of the orbital in VASP runs.
        """
        for orb in Orbital:
            if orb.value[1] == i:
                return orb
        raise ValueError("Invalid vasp index %s" % i)

    @staticmethod
    def from_string(orb_str):
        """
        Returns an orbital from a string representation, e.g., "s", "px".
        """
        for orb in Orbital:
            if orb.value[0] == orb_str:
                return orb
        raise ValueError("Illegal orbital definition!")
