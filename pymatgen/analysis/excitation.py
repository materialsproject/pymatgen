# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.core.spectrum import Spectrum

"""
This module defines an excitation spectrum class.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "May 10 2018"


class ExcitationSpectrum(Spectrum):
    """
    Basic excitation spectrum object.

    Args:
        x: A sequence of x-ray energies in eV
        y: A sequence of intensity values

    .. attribute: x
        The sequence of energies

    .. attribute: y
        The sequence of mu(E)

    """
    XLABEL = 'Energy (eV)'
    YLABEL = 'Intensity'

    def __init__(self, x, y):
        super().__init__(x, y)
