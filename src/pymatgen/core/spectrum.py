"""This module defines classes to represent any type of spectrum, essentially any
x y value pairs.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable
from scipy import stats
from scipy.ndimage import convolve1d

from pymatgen.util.coord import get_linear_interpolated_value

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Literal

    from numpy.typing import NDArray
    from typing_extensions import Self


def lorentzian(x: NDArray, x_0: float = 0, sigma: float = 1.0) -> NDArray:
    """The Lorentzian smearing function.

    Args:
        x: x values
        x_0: Center
        sigma: FWHM.

    Returns:
        Value of lorentzian at x.
    """
    return 1 / np.pi * 0.5 * sigma / ((x - x_0) ** 2 + (0.5 * sigma) ** 2)


class Spectrum(MSONable):
    """Base class for any type of XAS, essentially just x, y values. Examples
    include XRD patterns, XANES, EXAFS, NMR, DOS, etc.

    Implements basic tools like application of smearing, normalization, addition
    multiplication, etc.

    Subclasses should extend this object and ensure that super is called with
    ALL args and kwargs. That ensures subsequent things like add and mult work
    properly.
    """

    XLABEL = "x"
    YLABEL = "y"

    def __init__(self, x: NDArray, y: NDArray, *args, **kwargs) -> None:
        """
        Args:
            x (ndarray): A ndarray of N values.
            y (ndarray): A ndarray of N x k values. The first dimension must be
                the same as that of x. Each of the k values are interpreted as separate.
            *args: All subclasses should provide args other than x and y
                when calling super, e.g. super().__init__(
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

    def __getattr__(self, name: str) -> NDArray:
        if name == self.XLABEL.lower():
            return self.x
        if name == self.YLABEL.lower():
            return self.y
        raise AttributeError(f"Invalid attribute {name=}")

    def __len__(self) -> int:
        return self.ydim[0]

    def __add__(self, other: Self) -> Self:
        """Add two Spectrum object together. Checks that x scales are the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another Spectrum object

        Returns:
            Sum of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis values are not compatible!")

        return type(self)(self.x, self.y + other.y, *self._args, **self._kwargs)

    def __sub__(self, other: Self) -> Self:
        """Subtract one Spectrum object from another. Checks that x scales are
        the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another Spectrum object

        Returns:
            Subtraction of the two Spectrum objects
        """
        if not all(np.equal(self.x, other.x)):
            raise ValueError("X axis values are not compatible!")

        return type(self)(self.x, self.y - other.y, *self._args, **self._kwargs)

    def __mul__(self, other: Self) -> Self:
        """Scale the Spectrum's y values.

        Args:
            other: scalar, The scale amount

        Returns:
            Spectrum object with y values scaled
        """
        return type(self)(self.x, other * self.y, *self._args, **self._kwargs)

    __rmul__ = __mul__

    def __truediv__(self, other: Self) -> Self:
        """True division of y.

        Args:
            other: The divisor.

        Returns:
            Spectrum object with y values divided
        """
        return type(self)(self.x, self.y.__truediv__(other), *self._args, **self._kwargs)

    def __floordiv__(self, other: Self) -> Self:
        """True division of y.

        Args:
            other: The divisor.

        Returns:
            Spectrum object with y values divided
        """
        return type(self)(self.x, self.y.__floordiv__(other), *self._args, **self._kwargs)

    __div__ = __truediv__

    def __str__(self) -> str:
        """String containing values and labels of spectrum object for
        plotting.
        """
        return f"{type(self).__name__}\n{self.XLABEL}: {self.x}\n{self.YLABEL}: {self.y}"

    def __repr__(self) -> str:
        """A printable representation of the class."""
        return str(self)

    def normalize(
        self,
        mode: Literal["max", "sum"] = "max",
        value: float = 1.0,
    ) -> None:
        """Normalize the spectrum with respect to the sum of intensity.

        Args:
            mode ("max" | "sum"): Normalization mode. "max" sets the max y value to value,
                e.g. in XRD patterns. "sum" sets the sum of y to a value, i.e., like a
                probability density.
            value (float): Value to normalize to. Defaults to 1.
        """
        if mode.lower() == "sum":
            factor = np.sum(self.y, axis=0)
        elif mode.lower() == "max":
            factor = np.max(self.y, axis=0)
        else:
            raise ValueError(f"Unsupported normalization {mode=}!")

        self.y /= factor / value

    def smear(
        self,
        sigma: float = 0.0,
        func: Literal["gaussian", "lorentzian"] | Callable = "gaussian",
    ) -> None:
        """Apply Gaussian/Lorentzian smearing to spectrum y value.

        Args:
            sigma: Std dev for Gaussian smear function
            func: "gaussian" or "lorentzian" or a callable. If this is a callable, the sigma value is ignored. The
                callable should only take a single argument (a numpy array) and return a set of weights.
        """
        points = np.linspace(np.min(self.x) - np.mean(self.x), np.max(self.x) - np.mean(self.x), len(self.x))
        if callable(func):
            weights = func(points)
        elif func.lower() == "gaussian":
            weights = stats.norm.pdf(points, scale=sigma)
        elif func.lower() == "lorentzian":
            weights = lorentzian(points, sigma=sigma)
        else:
            raise ValueError(f"Invalid {func=}")
        weights /= np.sum(weights)
        if len(self.ydim) == 1:
            total = np.sum(self.y)
            self.y = convolve1d(self.y, weights)
            self.y *= total / np.sum(self.y)  # renormalize to maintain the same integrated sum as before.
        else:
            total = np.sum(self.y, axis=0)
            self.y = np.array([convolve1d(self.y[:, k], weights) for k in range(self.ydim[1])]).T
            self.y *= total / np.sum(self.y, axis=0)  # renormalize to maintain the same integrated sum as before.

    def get_interpolated_value(self, x: float) -> float | list[float]:
        """Get an interpolated y value for a particular x value.

        Args:
            x: x value to return the y value for

        Returns:
            Value of y at x
        """
        if len(self.ydim) == 1:
            return get_linear_interpolated_value(self.x, self.y, x)
        return [get_linear_interpolated_value(self.x, self.y[:, k], x) for k in range(self.ydim[1])]

    def copy(self) -> Self:
        """
        Returns:
            Copy of Spectrum object.
        """
        return type(self)(self.x, self.y, *self._args, **self._kwargs)
