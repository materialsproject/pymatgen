# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from monty.json import MSONable
from abc import ABC
import six
from pymatgen.core.units import EnergyArray, Energy
from scipy.interpolate import interp1d
from pymatgen import Structure
import numpy as np

"""
This module defines classes to represent all spectra
"""

__author__ = "Chen Zheng"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Chen Zheng"
__email__ = "chz022@ucsd.edu"
__date__ = "Aug 9, 2017"


class Spectrum(ABC, MSONable):
    """
    Base class for Xane spectrum, Exaf spectrum, NMR.
    Not meant to be instantiated directly
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y

    @property
    def x_value(self):
        return self.x

    @property
    def y_value(self):
        return self.y

    @property
    def xlabel(self):
        return self.xlabel

    @property
    def ylabel(self):
        return self.ylabel

    def intensity_sum_norm(self):
        """
        Normalize the spectrum with respect to the sum of intensity
        :return:
        """
        self.y = self.y / np.sum(self.y)

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__


class XANES(Spectrum):
    """
    Basic XANES object.

    Args:
        energies: A sequence of x-ray energies in eV
        mu: A sequence of mu(E)

    .. attribute: energies
        The sequence of energies
    .. attribute: mu
        The sequence of mu(E)
    .. attribute: absorption_specie
        The absorption_species of the spectrum
    """

    def __init__(self, x, y, structure, absorption_specie, edge):
        super(XANES, self).__init__(x, y)
        self.x = EnergyArray(x, 'eV')
        self.structure = structure
        self.absorption_specie = absorption_specie
        self.edge = edge

    def find_e0(self):
        """
        Use the maximum gradient to find e0
        """
        self.e0 = self.x[np.argmax(np.gradient(self.y) / np.gradient(self.x))]

    @classmethod
    def from_dict(cls, d):
        """
        Return XANES object from dict representation of XANES
        """
        return XANES(d['energies'], d['mu'], d['structure'],
                     d['absorption_specie'], d['edge'])

    def as_dict(self):
        """
        Json-serializable dict representation of XANES
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "energeis": list(self.energies),
                "mu": list(self.mu),
                "absorption_specie": list(self.absorption_specie),
                "edge": self.edge
                }
