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

from enum import Enum, unique

from monty.dev import deprecated


@unique
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
    @deprecated(message="from_int has been deprecated. Please use the Enum's "
                        "default API of calling Spin(1). Will be removed in "
                        "pymatgen 4.0.")
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


@unique
class OrbitalType(Enum):
    """
    Enum type for orbital type. Indices are basically the azimutal quantum
    number, l.
    """

    s = 0
    p = 1
    d = 2
    f = 3

    def __str__(self):
        return self.name


@unique
class Orbital(Enum):
    """
    Enum type for specific orbitals. The indices are basically the order in
    which the orbitals are reported in VASP and has no special meaning.
    """

    s = 0
    py = 1
    pz = 2
    px = 3
    dxy = 4
    dyz = 5
    dz2 = 6
    dxz = 7
    dx2 = 8
    f_3 = 9
    f_2 = 10
    f_1 = 11
    f0 = 12
    f1 = 13
    f2 = 14
    f3 = 15

    def __int__(self):
        return self.value

    def __str__(self):
        return self.name

    @property
    def orbital_type(self):
        """
        Returns OrbitalType of an orbital.
        """
        return OrbitalType[self.name[0]]

    @staticmethod
    @deprecated(message="from_vasp_index has been deprecated. Please use the "
                        "Enum's default API of calling Orbital(<vasp_index>). "
                        "Will be removed in pymatgen 4.0.")
    def from_vasp_index(i):
        """
        Returns an orbital based on the index of the orbital in VASP runs.
        """
        return Orbital(i)

