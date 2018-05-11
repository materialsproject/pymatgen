# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.core.spectrum import Spectrum
import numpy as np

"""
This module defines classes to represent all xas
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
    XLABEL = 'Energy'
    YLABEL = 'Intensity'

    def __init__(self, x, y):
        super(ExcitationSpectrum, self).__init__(x, y)
