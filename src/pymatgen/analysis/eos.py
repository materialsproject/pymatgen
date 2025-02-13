"""This module implements various equation of states.

Note: Most of the code were initially adapted from ASE and deltafactor by
@gmatteo but has since undergone major refactoring.
"""

from __future__ import annotations

import logging
import warnings
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import leastsq, minimize

try:
    from numpy.exceptions import RankWarning  # NPY2
except ImportError:
    from numpy import RankWarning  # NPY1

from pymatgen.core.units import FloatWithUnit
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig, pretty_plot

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, ClassVar

    import matplotlib.pyplot as plt

__author__ = "Kiran Mathew, gmatteo"
__credits__ = "Cormac Toher"

logger = logging.getLogger(__name__)


class EOSBase(ABC):
    """
    Abstract class that must be subclassed by all equation of state
    implementations.
    """

    def __init__(
        self,
        volumes: Sequence[float],
        energies: Sequence[float],
    ) -> None:
        """
        Args:
            volumes (Sequence[float]): in Ang^3.
            energies (Sequence[float]): in eV.
        """
        self.volumes = np.array(volumes)
        self.energies = np.array(energies)
        # minimum energy(e0), buk modulus(b0),
        # derivative of bulk modulus w.r.t. pressure(b1), minimum volume(v0)
        self._params: Sequence | None = None
        # the eos function parameters. It is the same as _params except for
        # equation of states that uses polynomial fits(delta_factor and
        # numerical_eos)
        self.eos_params: Sequence | None = None

    def __call__(self, volume: float) -> float:
        """
        Args:
            volume (float | list[float]): volume(s) in Ang^3.

        Returns:
            Compute EOS with this volume.
        """
        return self.func(volume)

    def _initial_guess(self) -> tuple[float, float, float, float]:
        """
        Quadratic fit to get an initial guess for the parameters.

        Returns:
            tuple[float, float, float, float]: e0, b0, b1, v0
        """
        a, b, c = np.polyfit(self.volumes, self.energies, 2)
        self.eos_params = [a, b, c]

        v0 = -b / (2 * a)
        e0 = a * (v0**2) + b * v0 + c
        b0 = 2 * a * v0
        b1 = 4  # b1 is usually a small number like 4

        vol_min, vol_max = min(self.volumes), max(self.volumes)

        if not vol_min < v0 and v0 < vol_max:
            raise EOSError("The minimum volume of a fitted parabola is not in the input volumes.")

        return e0, b0, b1, v0

    def fit(self) -> None:
        """
        Do the fitting. Does least square fitting. If you want to use custom
        fitting, must override this.
        """
        # the objective function that will be minimized in the least square fitting
        self._params = self._initial_guess()
        self.eos_params, ierr = leastsq(
            lambda pars, x, y: y - self._func(x, pars),
            self._params,
            args=(self.volumes, self.energies),
        )
        # e0, b0, b1, v0
        self._params = self.eos_params
        if ierr not in (1, 2, 3, 4):
            raise EOSError("Optimal parameters not found")

    @abstractmethod
    def _func(self, volume, params):
        """
        The equation of state function. This must be implemented by all classes
        that derive from this abstract class.

        Args:
            volume (float | list[float])
            params (list | tuple): values for the parameters other than the
                volume used by the eos.
        """

    def func(self, volume):
        """
        The equation of state function with the parameters other than volume set
        to the ones obtained from fitting.

        Args:
            volume (float | list[float]): volumes in Ang^3

        Returns:
            numpy.array
        """
        return self._func(np.array(volume), self.eos_params)

    @property
    def e0(self) -> float:
        """The min energy."""
        if self._params is None:
            raise RuntimeError("params have not be initialized.")

        return self._params[0]

    @property
    def b0(self) -> float:
        """The bulk modulus in units of energy/unit of volume^3."""
        if self._params is None:
            raise RuntimeError("params have not be initialized.")

        return self._params[1]

    @property
    def b0_GPa(self) -> FloatWithUnit:
        """The bulk modulus in GPa. This assumes the energy and volumes are in eV and Ang^3."""
        return FloatWithUnit(self.b0, "eV ang^-3").to("GPa")

    @property
    def b1(self):
        """The derivative of bulk modulus w.r.t. pressure(dimensionless)."""
        return self._params[2]

    @property
    def v0(self):
        """The minimum or the reference volume in Ang^3."""
        return self._params[3]

    @property
    def results(self) -> dict[str, Any]:
        """A summary dict."""
        return {"e0": self.e0, "b0": self.b0, "b1": self.b1, "v0": self.v0}

    def plot(
        self,
        width: float = 8,
        height: float | None = None,
        ax: plt.Axes = None,
        dpi: float | None = None,
        **kwargs,
    ) -> plt.Axes:
        """
        Plot the equation of state.

        Args:
            width (float): Width of plot in inches. Defaults to 8in.
            height (float): Height of plot in inches. Defaults to width *
                golden ratio.
            ax (plt.Axes): If supplied, changes will be made to the existing Axes.
                Otherwise, new Axes will be created.
            dpi (float): DPI.
            kwargs (dict): additional args fed to pyplot.plot.
                supported keys: style, color, text, label

        Returns:
            plt.Axes: The matplotlib axes.
        """
        ax = pretty_plot(width=width, height=height, ax=ax, dpi=dpi)

        color = kwargs.get("color", "r")
        label = kwargs.get("label", f"{type(self).__name__} fit")
        lines = [
            f"Equation of State: {type(self).__name__}",
            f"Minimum energy = {self.e0:1.2f} eV",
            f"Minimum or reference volume = {self.v0:1.2f} Ang^3",
            f"Bulk modulus = {self.b0:1.2f} eV/Ang^3 = {self.b0_GPa:1.2f} GPa",
            f"Derivative of bulk modulus w.r.t. pressure = {self.b1:1.2f}",
        ]
        text = "\n".join(lines)
        text = kwargs.get("text", text)

        # Plot input data.
        ax.plot(self.volumes, self.energies, linestyle="None", marker="o", color=color)

        # Plot eos fit.
        vmin, vmax = min(self.volumes), max(self.volumes)
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        vfit = np.linspace(vmin, vmax, 100)

        ax.plot(vfit, self.func(vfit), linestyle="dashed", color=color, label=label)

        ax.grid(visible=True)
        ax.set_xlabel("Volume $\\AA^3$")
        ax.set_ylabel("Energy (eV)")
        ax.legend(loc="best", shadow=True)
        # Add text with fit parameters.
        ax.text(0.4, 0.5, text, transform=ax.transAxes)

        return ax

    @add_fig_kwargs
    def plot_ax(
        self,
        ax: plt.Axes | None = None,
        fontsize: float = 12,
        **kwargs,
    ) -> plt.Figure:
        """
        Plot the equation of state on axis `ax`.

        Args:
            ax: matplotlib Axes or None if a new figure should be created.
            fontsize: Legend fontsize.

        Returns:
            plt.Figure: matplotlib figure.
        """
        ax, fig = get_ax_fig(ax=ax)

        color = kwargs.get("color", "r")
        label = kwargs.get("label", f"{type(self).__name__} fit")
        lines = [
            f"Equation of State: {type(self).__name__}",
            f"Minimum energy = {self.e0:1.2f} eV",
            f"Minimum or reference volume = {self.v0:1.2f} Ang^3",
            f"Bulk modulus = {self.b0:1.2f} eV/Ang^3 = {self.b0_GPa:1.2f} GPa",
            f"Derivative of bulk modulus w.r.t. pressure = {self.b1:1.2f}",
        ]
        text = "\n".join(lines)
        text = kwargs.get("text", text)

        # Plot input data.
        ax.plot(self.volumes, self.energies, linestyle="None", marker="o", color=color)

        # Plot eos fit.
        vmin, vmax = min(self.volumes), max(self.volumes)
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        vfit = np.linspace(vmin, vmax, 100)

        ax.plot(vfit, self.func(vfit), linestyle="dashed", color=color, label=label)

        ax.grid(visible=True)
        ax.set_xlabel("Volume $\\AA^3$")
        ax.set_ylabel("Energy (eV)")
        ax.legend(loc="best", shadow=True)
        # Add text with fit parameters.
        ax.text(
            0.5,
            0.5,
            text,
            fontsize=fontsize,
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
        )

        return fig


class Murnaghan(EOSBase):
    """Murnaghan EOS."""

    def _func(self, volume, params: tuple[float, float, float, float]):
        """From PRB 28,5480 (1983)."""
        e0, b0, b1, v0 = tuple(params)
        return e0 + b0 * volume / b1 * (((v0 / volume) ** b1) / (b1 - 1.0) + 1.0) - v0 * b0 / (b1 - 1.0)


class Birch(EOSBase):
    """Birch EOS."""

    def _func(self, volume, params: tuple[float, float, float, float]):
        """From Intermetallic compounds: Principles and Practice, Vol. I:
        Principles Chapter 9 pages 195-210 by M. Mehl. B. Klein,
        D. Papaconstantopoulos.
        case where n=0.
        """
        e0, b0, b1, v0 = tuple(params)
        return (
            e0
            + 9 / 8 * b0 * v0 * ((v0 / volume) ** (2 / 3.0) - 1.0) ** 2
            + 9 / 16 * b0 * v0 * (b1 - 4.0) * ((v0 / volume) ** (2 / 3.0) - 1.0) ** 3
        )


class BirchMurnaghan(EOSBase):
    """BirchMurnaghan EOS."""

    def _func(self, volume, params: tuple[float, float, float, float]):
        """BirchMurnaghan equation from PRB 70, 224107."""
        e0, b0, b1, v0 = tuple(params)
        eta = (v0 / volume) ** (1 / 3)
        return e0 + 9 * b0 * v0 / 16 * (eta**2 - 1) ** 2 * (6 + b1 * (eta**2 - 1.0) - 4 * eta**2)


class PourierTarantola(EOSBase):
    """Pourier-Tarantola EOS."""

    def _func(self, volume, params: tuple[float, float, float, float]):
        """Pourier-Tarantola equation from PRB 70, 224107."""
        e0, b0, b1, v0 = tuple(params)
        eta = (volume / v0) ** (1 / 3)
        squiggle = -3 * np.log(eta)
        return e0 + b0 * v0 * squiggle**2 / 6 * (3 + squiggle * (b1 - 2))


class Vinet(EOSBase):
    """Vinet EOS."""

    def _func(self, volume, params: tuple[float, float, float, float]):
        """Vinet equation from PRB 70, 224107."""
        e0, b0, b1, v0 = tuple(params)
        eta = (volume / v0) ** (1 / 3)
        return e0 + 2 * b0 * v0 / (b1 - 1.0) ** 2 * (
            2 - (5 + 3 * b1 * (eta - 1.0) - 3 * eta) * np.exp(-3 * (b1 - 1.0) * (eta - 1.0) / 2.0)
        )


class PolynomialEOS(EOSBase):
    """
    Derives from EOSBase. Polynomial based equations of states must subclass
    this.
    """

    def _func(self, volume, params):
        return np.poly1d(list(params))(volume)

    def fit(self, order: int) -> None:
        """
        Do polynomial fitting and set the parameters. Uses numpy polyfit.

        Args:
            order (int): order of the fit polynomial
        """
        self.eos_params = np.polyfit(self.volumes, self.energies, order)
        self._set_params()

    def _set_params(self) -> None:
        """
        Use the fit polynomial to compute the parameter e0, b0, b1 and v0
        and set to the _params attribute.
        """
        fit_poly = np.poly1d(self.eos_params)
        # the volume at min energy, used as the initial guess for the optimization w.r.t. volume.
        v_e_min = self.volumes[np.argmin(self.energies)]
        # evaluate e0, v0, b0 and b1
        min_wrt_v = minimize(fit_poly, v_e_min)
        e0, v0 = min_wrt_v.fun, min_wrt_v.x[0]
        pderiv2 = np.polyder(fit_poly, 2)
        pderiv3 = np.polyder(fit_poly, 3)
        b0 = v0 * np.poly1d(pderiv2)(v0)
        db0dv = np.poly1d(pderiv2)(v0) + v0 * np.poly1d(pderiv3)(v0)
        # db/dp
        b1 = -v0 * db0dv / b0
        self._params = [e0, b0, b1, v0]


class DeltaFactor(PolynomialEOS):
    """Fitting a polynomial EOS using delta factor."""

    def _func(self, volume, params):
        x = volume ** (-2 / 3.0)
        return np.poly1d(list(params))(x)

    def fit(self, order: int = 3) -> None:
        """Overridden since this eos works with volume**(2/3) instead of volume."""
        x = self.volumes ** (-2 / 3.0)
        self.eos_params = np.polyfit(x, self.energies, order)
        self._set_params()

    def _set_params(self):
        """
        Overridden to account for the fact the fit with volume**(2/3) instead
        of volume.
        """
        deriv0 = np.poly1d(self.eos_params)
        deriv1 = np.polyder(deriv0, 1)
        deriv2 = np.polyder(deriv1, 1)
        deriv3 = np.polyder(deriv2, 1)

        for x in np.roots(deriv1):
            if x > 0 and deriv2(x) > 0:
                v0 = x ** (-3 / 2.0)
                break
        else:
            raise EOSError("No minimum could be found")

        derivV2 = 4 / 9 * x**5 * deriv2(x)
        derivV3 = -20 / 9 * x ** (13 / 2.0) * deriv2(x) - 8 / 27 * x ** (15 / 2.0) * deriv3(x)
        b0 = derivV2 / x ** (3 / 2.0)
        b1 = -1 - x ** (-3 / 2.0) * derivV3 / derivV2

        # e0, b0, b1, v0
        self._params = [deriv0(v0 ** (-2 / 3.0)), b0, b1, v0]


class NumericalEOS(PolynomialEOS):
    """A numerical EOS."""

    def fit(
        self,
        min_ndata_factor: int = 3,
        max_poly_order_factor: int = 5,
        min_poly_order: int = 2,
    ) -> None:
        """Fit the input data to the 'numerical eos', the equation of state employed
        in the quasiharmonic Debye model described in the paper:
        10.1103/PhysRevB.90.174107.

        credits: Cormac Toher

        Args:
            min_ndata_factor (int): parameter that controls the minimum number
                of data points that will be used for fitting.
                minimum number of data points = total data points-2*min_ndata_factor
            max_poly_order_factor (int): parameter that limits the max order
                of the polynomial used for fitting.
                max_poly_order = number of data points used for fitting -
                max_poly_order_factor
            min_poly_order (int): minimum order of the polynomial to be
                considered for fitting.
        """
        warnings.simplefilter("ignore", RankWarning)

        def get_rms(x, y):
            return np.sqrt(np.sum((np.array(x) - np.array(y)) ** 2) / len(x))

        # list of (energy, volume) tuples
        e_v = list(zip(self.energies, self.volumes, strict=True))
        n_data = len(e_v)
        # minimum number of data points used for fitting
        n_data_min = max(n_data - 2 * min_ndata_factor, min_poly_order + 1)
        rms_min = np.inf
        # number of data points available for fit in each iteration
        n_data_fit = n_data
        # store the fit polynomial coefficients and the rms in a dict,
        # where the key=(polynomial order, number of data points used for
        # fitting)
        all_coeffs = {}

        # sort by energy
        e_v = sorted(e_v, key=lambda x: x[0])
        # minimum energy tuple
        e_min = e_v[0]
        # sort by volume
        e_v = sorted(e_v, key=lambda x: x[1])
        # index of minimum energy tuple in the volume sorted list
        e_min_idx = e_v.index(e_min)
        # the volume lower than the volume corresponding to minimum energy
        v_before = e_v[e_min_idx - 1][1]
        # the volume higher than the volume corresponding to minimum energy
        v_after = e_v[e_min_idx + 1][1]
        e_v_work = deepcopy(e_v)

        # loop over the data points.
        while (n_data_fit >= n_data_min) and (e_min in e_v_work):
            max_poly_order = n_data_fit - max_poly_order_factor
            energies = [ei[0] for ei in e_v_work]
            volumes = [ei[1] for ei in e_v_work]
            # loop over polynomial order
            for idx in range(min_poly_order, max_poly_order + 1):
                coeffs = np.polyfit(volumes, energies, idx)
                polyder = np.polyder(coeffs)
                a = np.poly1d(polyder)(v_before)
                b = np.poly1d(polyder)(v_after)
                if a * b < 0:
                    rms = get_rms(energies, np.poly1d(coeffs)(volumes))
                    rms_min = min(rms_min, rms * idx / n_data_fit)
                    all_coeffs[idx, n_data_fit] = [coeffs.tolist(), rms]
                    # store the fit coefficients small to large,
                    # i.e a0, a1, .. an
                    all_coeffs[idx, n_data_fit][0].reverse()
            # remove 1 data point from each end.
            e_v_work.pop()
            e_v_work.pop(0)
            n_data_fit = len(e_v_work)

        logger.info(f"total number of polynomials: {len(all_coeffs)}")

        norm = 0.0
        fit_poly_order = n_data
        # weight average polynomial coefficients.
        weighted_avg_coeffs = np.zeros((fit_poly_order,))

        # combine all the filtered polynomial candidates to get the final fit.
        for key, val in all_coeffs.items():
            # weighted rms = rms * polynomial order / rms_min / ndata_fit
            weighted_rms = val[1] * key[0] / rms_min / key[1]
            weight = np.exp(-(weighted_rms**2))
            norm += weight
            coeffs = np.array(val[0])
            # pad the coefficient array with zeros
            coeffs = np.pad(coeffs, (0, max(fit_poly_order - len(coeffs), 0)), "constant")
            weighted_avg_coeffs += weight * coeffs

        # normalization
        weighted_avg_coeffs /= norm
        weighted_avg_coeffs = weighted_avg_coeffs.tolist()
        # large to small(an, an-1, ..., a1, a0) as expected by np.poly1d
        weighted_avg_coeffs.reverse()

        self.eos_params = weighted_avg_coeffs
        self._set_params()


class EOS:
    """
    Convenient wrapper. Retained in its original state to ensure backward
    compatibility.

    Fit equation of state for bulk systems.

    The following equations are supported:

        murnaghan: PRB 28, 5480 (1983)

        birch: Intermetallic compounds: Principles and Practice, Vol I:
            Principles. pages 195-210

        birch_murnaghan: PRB 70, 224107

        pourier_tarantola: PRB 70, 224107

        vinet: PRB 70, 224107

        deltafactor

        numerical_eos: 10.1103/PhysRevB.90.174107.

    Usage:

       eos = EOS(eos_name='murnaghan')
       eos_fit = eos.fit(volumes, energies)
       eos_fit.plot()
    """

    MODELS: ClassVar[dict[str, Any]] = {
        "murnaghan": Murnaghan,
        "birch": Birch,
        "birch_murnaghan": BirchMurnaghan,
        "pourier_tarantola": PourierTarantola,
        "vinet": Vinet,
        "deltafactor": DeltaFactor,
        "numerical_eos": NumericalEOS,
    }

    def __init__(self, eos_name: str = "murnaghan") -> None:
        """
        Args:
            eos_name (str): Type of EOS to fit.
        """
        if eos_name not in self.MODELS:
            raise EOSError(
                f"The equation of state {eos_name!r} is not supported. "
                f"Please choose one from the following list: {list(self.MODELS)}"
            )
        self._eos_name = eos_name
        self.model = self.MODELS[eos_name]

    def fit(self, volumes: Sequence[float], energies: Sequence[float]) -> EOSBase:
        """Fit energies as function of volumes.

        Args:
            volumes (Sequence[float]): in Ang^3
            energies (Sequence[float]): in eV

        Returns:
            EOSBase: EOSBase object
        """
        eos_fit = self.model(np.array(volumes), np.array(energies))
        eos_fit.fit()
        return eos_fit


class EOSError(Exception):
    """Error class for EOS fitting."""
