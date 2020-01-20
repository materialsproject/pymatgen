# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines classes to represent any type of spectrum, essentially any
x y value pairs.
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

from monty.json import MSONable

from pymatgen.util.coord import get_linear_interpolated_value

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
        r"""
        Args:
            x (ndarray): A ndarray of N values.
            y (ndarray): A ndarray of N x k values. The first dimension must be
                the same as that of x. Each of the k values are interpreted as
            *args: All subclasses should provide args other than x and y
                when calling super, e.g., super().__init__(
                x, y, arg1, arg2, kwarg1=val1, ..). This guarantees the +, -, *,
                etc. operators work properly.
            **kwargs: Same as that for *args.
        """
        self.x = np.array(x)
        self.y = np.array(y)
        self.ydim = self.y.shape
        if self.x.shape[0] != self.ydim[0]:
            raise ValueError("x and y values have different first dimension!")
        self._args = args
        self._kwargs = kwargs

    def __getattr__(self, item):
        if item == self.XLABEL.lower():
            return self.x
        if item == self.YLABEL.lower():
            return self.y
        raise AttributeError("Invalid attribute name %s" % str(item))

    def __len__(self):
        return self.ydim[0]

    def normalize(self, mode="max", value=1):
        """
        Normalize the spectrum with respect to the sum of intensity

        Args:
            mode (str): Normalization mode. Supported modes are "max" (set the
                max y value to value, e.g., in XRD patterns), "sum" (set the
                sum of y to a value, i.e., like a probability density).
            value (float): Value to normalize to. Defaults to 1.
        """
        if mode.lower() == "sum":
            factor = np.sum(self.y, axis=0)
        elif mode.lower() == "max":
            factor = np.max(self.y, axis=0)
        else:
            raise ValueError("Unsupported normalization mode %s!" % mode)

        self.y /= factor / value

    def smear(self, sigma):
        """
        Apply Gaussian smearing to spectrum y value.

        Args:
            sigma: Std dev for Gaussian smear function
        """
        diff = [self.x[i + 1] - self.x[i] for i in range(len(self.x) - 1)]
        avg_x_per_step = np.sum(diff) / len(diff)
        if len(self.ydim) == 1:
            self.y = gaussian_filter1d(self.y, sigma / avg_x_per_step)
        else:
            self.y = np.array([
                gaussian_filter1d(self.y[:, k], sigma / avg_x_per_step)
                for k in range(self.ydim[1])]).T

    def get_interpolated_value(self, x):
        """
        Returns an interpolated y value for a particular x value.

        Args:
             x: x value to return the y value for

        Returns:
            Value of y at x
        """
        if len(self.ydim) == 1:
            return get_linear_interpolated_value(self.x, self.y, x)
        return [get_linear_interpolated_value(self.x, self.y[:, k], x)
                for k in range(self.ydim[1])]

    def copy(self):
        """
        Returns:
            Copy of Spectrum object.
        """
        return self.__class__(self.x, self.y, *self._args, **self._kwargs)

    def __add__(self, other):
        """
        Add two Spectrum object together. Checks that x scales are the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another Spectrum object

        Returns:
            Sum of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis values are not compatible!")
        return self.__class__(self.x, self.y + other.y, *self._args,
                              **self._kwargs)

    def __sub__(self, other):
        """
        Substract one Spectrum object from another. Checks that x scales are
        the same.
        Otherwise, a ValueError is thrown

        Args:
            other: Another Spectrum object

        Returns:
            Substraction of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis values are not compatible!")
        return self.__class__(self.x, self.y - other.y, *self._args,
                              **self._kwargs)

    def __mul__(self, other):
        """
        Scale the Spectrum's y values

        Args:
            other: scalar, The scale amount
        Returns:
            Spectrum object with y values scaled
        """
        return self.__class__(self.x, other * self.y, *self._args,
                              **self._kwargs)
    __rmul__ = __mul__

    def __truediv__(self, other):
        """
        True division of y
        Args:
            other: The divisor

        Returns:
            Spectrum object with y values divided
        """
        return self.__class__(self.x, self.y.__truediv__(other), *self._args,
                              **self._kwargs)

    def __floordiv__(self, other):
        """
        True division of y
        Args:
            other: The divisor

        Returns:
            Spectrum object with y values divided
        """
        return self.__class__(self.x, self.y.__floordiv__(other), *self._args,
                              **self._kwargs)

    __div__ = __truediv__

    def __str__(self):
        """
        Returns a string containing values and labels of spectrum object for
        plotting.
        """
        return "\n".join([self.__class__.__name__,
                          "%s: %s" % (self.XLABEL, self.x),
                          "%s: %s" % (self.YLABEL, self.y)])

    def __repr__(self):
        """
        Returns a printable representation of the class
        """
        return self.__str__()
