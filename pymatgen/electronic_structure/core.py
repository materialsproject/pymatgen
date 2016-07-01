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
