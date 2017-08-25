# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from monty.json import MSONable
from pymatgen.core.units import EnergyArray, Energy
import numpy as np
from pymatgen.util.coord_utils import get_linear_interpolated_value

"""
This module defines classes to represent all spectra
"""

__author__ = "Chen Zheng"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Chen Zheng"
__email__ = "chz022@ucsd.edu"
__date__ = "Aug 9, 2017"


class Spectrum(MSONable):
    """
    Base class for Xane spectrum, Exaf spectrum, NMR.
    Not meant to be instantiated directly
    """

    def __init__(self, x, y):
        self.x = np.array(x)
        self.y = np.array(y)

    def intensity_sum_norm(self):
        """
        Normalize the spectrum with respect to the sum of intensity
        :return:
        """
        self.y = self.y / np.sum(self.y)

    def get_smeared_density(self, sigma):
        """
        Return Gaussian smeared spectrum
        :param sigma: Std dev of Gaussian smearing function
        :return: Gaussian-smeared densities
        """

        from scipy.ndimage.filters import gaussian_filter1d
        diff = [self.x[i + 1] - self.x[i] for i in range(len(self.x) - 1)]
        avg_eV_per_step = np.sum(diff) / len(diff)
        rv = gaussian_filter1d(self.y, sigma / avg_eV_per_step).tolist()
        return rv

    def get_interpolated_value(self, x_value):
        """
        Returns an interpolated y value for a particular x value
        :param x_value: x value to return the y value for
        """
        y_value = get_linear_interpolated_value(self.x, self.y, x_value)
        return y_value

    def __add__(self, other):
        """
        Add two Spectrum object together. Checks that x scales are the same.
        Otherwise, a ValueError is thrown
        :param other: Another Spectrum object
        :return: Sum of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis Values of both spectra are not compatible!")
        y_value = self.y + other.y
        return Spectrum(self.x, y_value)

    def __mul__(self, other):
        """
        Scale the Spectrum's y values
        :param other: The scale amount
        :return: Spectrum object with y values scaled
        """
        return Spectrum(self.x, other * self.y)

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__


class XANES(Spectrum):
    """
    Basic XANES object.

    Args:
        x: A sequence of x-ray energies in eV
        y: A sequence of mu(E)

    .. attribute: x
        The sequence of energies
    .. attribute: y
        The sequence of mu(E)
    .. attribute: absorption_specie
        The absorption_species of the spectrum
    .. attribute: edge
        The edge of XANES spectrum
    """

    def __init__(self, x, y, structure, absorption_specie, edge):
        super(XANES, self).__init__(x, y)
        self.x = EnergyArray(x, 'eV')
        self.structure = structure
        self.absorption_specie = absorption_specie
        self.edge = edge
        self.xlabel = 'eV'
        self.ylabel = 'mu'

    def find_e0(self):
        """
        Use the maximum gradient to find e0
        """
        self.e0 = Energy(self.x[np.argmax(np.gradient(self.y) / np.gradient(self.x))], "eV")

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
