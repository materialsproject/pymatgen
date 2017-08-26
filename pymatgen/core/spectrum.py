# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from monty.json import MSONable
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

from pymatgen.util.coord_utils import get_linear_interpolated_value

"""
This module defines classes to represent any type of spectrum, essentially any
x y value pairs. 
"""

__author__ = "Chen Zheng"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Chen Zheng"
__email__ = "chz022@ucsd.edu"
__date__ = "Aug 9, 2017"


class Spectrum(MSONable):
    """
    Base class for any type of xas, essentially just x, y values. Examples
    include XRD patterns, XANES, EXAFS, NMR, DOS, etc.

    Implements basic tools like application of smearing, normalization, addition
    multiplication, etc.

    Subclasses should extend this object and ensure that super is called with
    ALL args and kwargs. That ensures subsequent things like add and mult work
    properly.
    """
    XLABEL = "x"
    YLABEL = "y"

    def __init__(self, x, y, *args, **kwargs):
        if len(x) != len(y):
            raise ValueError("x and y values have different lengths!")
        self.x = np.array(x)
        self.y = np.array(y)
        self._args = args
        self._kwargs = kwargs

    def normalize(self, mode="max"):
        """
        Normalize the spectrum with respect to the sum of intensity

        Args:
            mode (max): Normalization mode. Support modes are "max" (set the
                max y value to 1, e.g., in XRD patterns), "sum" (set the sum of
                y to 1, i.e., like a probability density)
        """
        factor = np.sum(self.y) if mode == "sum" else np.max(self.y)
        self.y = self.y / factor

    def smear(self, sigma):
        """
        Apply Gaussian smearing.
        """
        diff = [self.x[i + 1] - self.x[i] for i in range(len(self.x) - 1)]
        avg_x_per_step = np.sum(diff) / len(diff)
        self.y = gaussian_filter1d(self.y, sigma / avg_x_per_step)

    def get_interpolated_value(self, x_value):
        """
        Returns an interpolated y value for a particular x value
        :param x_value: x value to return the y value for
        """
        return get_linear_interpolated_value(self.x, self.y, x_value)

    def __add__(self, other):
        """
        Add two Spectrum object together. Checks that x scales are the same.
        Otherwise, a ValueError is thrown
        :param other: Another Spectrum object
        :return: Sum of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis values are not compatible!")
        return self.__class__(self.x, self.y + other.y, *self._args,
                              **self._kwargs)

    def __sub__(self, other):
        """
        Add two Spectrum object together. Checks that x scales are the same.
        Otherwise, a ValueError is thrown
        :param other: Another Spectrum object
        :return: Sum of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis values are not compatible!")
        return self.__class__(self.x, self.y - other.y, *self._args,
                              **self._kwargs)

    def __mul__(self, other):
        """
        Scale the Spectrum's y values
        :param other: The scale amount
        :return: Spectrum object with y values scaled
        """
        return self.__class__(self.x, other * self.y, *self._args,
                              **self._kwargs)
    __rmul__ = __mul__

    def __truediv__(self, other):
        """
        Scale the Spectrum's y values
        :param other: The scale amount
        :return: Spectrum object with y values scaled
        """
        return self.__class__(self.x, self.y / other, *self._args,
                              **self._kwargs)

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__
