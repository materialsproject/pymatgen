"""
This module implements various equation of states.

Note: Most of the code were initially adapted from ASE and deltafactor by
@gmatteo but has since undergone major refactoring.
"""

from __future__ import annotations

import logging
import warnings
from abc import ABCMeta, abstractmethod
from copy import deepcopy

import numpy as np
from scipy.optimize import leastsq, minimize

from pymatgen.core.units import FloatWithUnit
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt, pretty_plot

__author__ = "Kiran Mathew, gmatteo"
__credits__ = "Cormac Toher"

logger = logging.getLogger(__file__)


class EOSBase(metaclass=ABCMeta):
    """
    Abstract class that must be subclassed by all equation of state
    implementations.
    """

    def __init__(self, volumes, energies):
        """
        Args:
            volumes (list/numpy.array): volumes in Ang^3
            energies (list/numpy.array): energy in eV.
        """
        self.volumes = np.array(volumes)
        self.energies = np.array(energies)
        # minimum energy(e0), buk modulus(b0),
        # derivative of bulk modulus wrt pressure(b1), minimum volume(v0)
        self._params = None
        # the eos function parameters. It is the same as _params except for
        # equation of states that uses polynomial fits(deltafactor and
        # numerical_eos)
        self.eos_params = None

    def _initial_guess(self):
        """
        Quadratic fit to get an initial guess for the parameters.

        Returns:
            tuple: (e0, b0, b1, v0)
        """
        a, b, c = np.polyfit(self.volumes, self.energies, 2)
        self.eos_params = [a, b, c]

        v0 = -b / (2 * a)
        e0 = a * (v0**2) + b * v0 + c
        b0 = 2 * a * v0
        b1 = 4  # b1 is usually a small number like 4

        vmin, vmax = min(self.volumes), max(self.volumes)

        if not vmin < v0 and v0 < vmax:
            raise EOSError("The minimum volume of a fitted parabola is not in the input volumes\n.")

        return e0, b0, b1, v0

    def fit(self):
        """
        Do the fitting. Does least square fitting. If you want to use custom
        fitting, must override this.
        """
        # the objective function that will be minimized in the least square
        # fitting
        self._params = self._initial_guess()
        self.eos_params, ierr = leastsq(
            lambda pars, x, y: y - self._func(x, pars),
            self._params,
            args=(self.volumes, self.energies),
        )
        # e0, b0, b1, v0
        self._params = self.eos_params
        if ierr not in [1, 2, 3, 4]:
            raise EOSError("Optimal parameters not found")

    @abstractmethod
    def _func(self, volume, params):
        """
        The equation of state function. This must be implemented by all classes
        that derive from this abstract class.

        Args:
            volume (float/numpy.array)
             params (list/tuple): values for the parameters other than the
                volume used by the eos.
        """

    def func(self, volume):
        """
        The equation of state function with the parameters other than volume set
        to the ones obtained from fitting.

        Args:
             volume (list/numpy.array)

        Returns:
            numpy.array
        """
        return self._func(np.array(volume), self.eos_params)

    def __call__(self, volume):
        """
        Args:
            volume (): Volume.

        Returns:
            Compute EOS with this volume.
        """
        return self.func(volume)

    @property
    def e0(self):
        """Returns the min energy."""
        return self._params[0]

    @property
    def b0(self):
        """
        Returns the bulk modulus.
        Note: the units for the bulk modulus: unit of energy/unit of volume^3.
        """
        return self._params[1]

    @property
    def b0_GPa(self):
        """
        Returns the bulk modulus in GPa.
        Note: This assumes that the energy and volumes are in eV and Ang^3
            respectively.
        """
        return FloatWithUnit(self.b0, "eV ang^-3").to("GPa")

    @property
    def b1(self):
        """Returns the derivative of bulk modulus wrt pressure(dimensionless)."""
        return self._params[2]

    @property
    def v0(self):
        """Returns the minimum or the reference volume in Ang^3."""
        return self._params[3]

    @property
    def results(self):
        """
        Returns a summary dict.

        Returns:
            dict
        """
        return {"e0": self.e0, "b0": self.b0, "b1": self.b1, "v0": self.v0}

    def plot(self, width=8, height=None, plt=None, dpi=None, **kwargs):
        """
        Plot the equation of state.

        Args:
            width (float): Width of plot in inches. Defaults to 8in.
            height (float): Height of plot in inches. Defaults to width *
                golden ratio.
            plt (matplotlib.pyplot): If plt is supplied, changes will be made
                to an existing plot. Otherwise, a new plot will be created.
            dpi:
            kwargs (dict): additional args fed to pyplot.plot.
                supported keys: style, color, text, label

        Returns:
            Matplotlib plot object.
        """
        # pylint: disable=E1307
        plt = pretty_plot(width=width, height=height, plt=plt, dpi=dpi)

        color = kwargs.get("color", "r")
        label = kwargs.get("label", f"{type(self).__name__} fit")
        lines = [
            f"Equation of State: {type(self).__name__}",
            f"Minimum energy = {self.e0:1.2f} eV",
            f"Minimum or reference volume = {self.v0:1.2f} Ang^3",
            f"Bulk modulus = {self.b0:1.2f} eV/Ang^3 = {self.b0_GPa:1.2f} GPa",
            f"Derivative of bulk modulus wrt pressure = {self.b1:1.2f}",
        ]
        text = "\n".join(lines)
        text = kwargs.get("text", text)

        # Plot input data.
        plt.plot(self.volumes, self.energies, linestyle="None", marker="o", color=color)

        # Plot eos fit.
        vmin, vmax = min(self.volumes), max(self.volumes)
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        vfit = np.linspace(vmin, vmax, 100)

        plt.plot(vfit, self.func(vfit), linestyle="dashed", color=color, label=label)

        plt.grid(True)
        plt.xlabel("Volume $\\AA^3$")
        plt.ylabel("Energy (eV)")
        plt.legend(loc="best", shadow=True)
        # Add text with fit parameters.
        plt.text(0.4, 0.5, text, transform=plt.gca().transAxes)

        return plt

    @add_fig_kwargs
    def plot_ax(self, ax=None, fontsize=12, **kwargs):
        """
        Plot the equation of state on axis `ax`.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            fontsize: Legend fontsize.
            color (str): plot color.
            label (str): Plot label
            text (str): Legend text (options)

        Returns:
            Matplotlib figure object.
        """
        # pylint: disable=E1307
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        color = kwargs.get("color", "r")
        label = kwargs.get("label", f"{type(self).__name__} fit")
        lines = [
            f"Equation of State: {type(self).__name__}",
            f"Minimum energy = {self.e0:1.2f} eV",
            f"Minimum or reference volume = {self.v0:1.2f} Ang^3",
            f"Bulk modulus = {self.b0:1.2f} eV/Ang^3 = {self.b0_GPa:1.2f} GPa",
            f"Derivative of bulk modulus wrt pressure = {self.b1:1.2f}",
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

        ax.grid(True)
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

    def _func(self, volume, params):
        """From PRB 28,5480 (1983)."""
        e0, b0, b1, v0 = tuple(params)
        return e0 + b0 * volume / b1 * (((v0 / volume) ** b1) / (b1 - 1.0) + 1.0) - v0 * b0 / (b1 - 1.0)


class Birch(EOSBase):
    """Birch EOS."""

    def _func(self, volume, params):
        """
        From Intermetallic compounds: Principles and Practice, Vol. I:
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

    def _func(self, volume, params):
        """BirchMurnaghan equation from PRB 70, 224107."""
        e0, b0, b1, v0 = tuple(params)
        eta = (v0 / volume) ** (1 / 3)
        return e0 + 9 * b0 * v0 / 16 * (eta**2 - 1) ** 2 * (6 + b1 * (eta**2 - 1.0) - 4 * eta**2)


class PourierTarantola(EOSBase):
    """PourierTarantola EOS."""

    def _func(self, volume, params):
        """Pourier-Tarantola equation from PRB 70, 224107."""
        e0, b0, b1, v0 = tuple(params)
        eta = (volume / v0) ** (1 / 3)
        squiggle = -3 * np.log(eta)
        return e0 + b0 * v0 * squiggle**2 / 6 * (3 + squiggle * (b1 - 2))


class Vinet(EOSBase):
    """Vinet EOS."""

    def _func(self, volume, params):
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

    def fit(self, order):
        """
        Do polynomial fitting and set the parameters. Uses numpy polyfit.

        Args:
             order (int): order of the fit polynomial
        """
        self.eos_params = np.polyfit(self.volumes, self.energies, order)
        self._set_params()

    def _set_params(self):
        """
        Use the fit polynomial to compute the parameter e0, b0, b1 and v0
        and set to the _params attribute.
        """
        fit_poly = np.poly1d(self.eos_params)
        # the volume at min energy, used as the initial guess for the
        # optimization wrt volume.
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

    def fit(self, order=3):
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

    def fit(self, min_ndata_factor=3, max_poly_order_factor=5, min_poly_order=2):
        """
        Fit the input data to the 'numerical eos', the equation of state employed
        in the quasiharmonic Debye model described in the paper:
        10.1103/PhysRevB.90.174107.

        credits: Cormac Toher

        Args:
            min_ndata_factor (int): parameter that controls the minimum number
                of data points that will be used for fitting.
                minimum number of data points =
                    total data points-2*min_ndata_factor
            max_poly_order_factor (int): parameter that limits the max order
                of the polynomial used for fitting.
                max_poly_order = number of data points used for fitting -
                                 max_poly_order_factor
            min_poly_order (int): minimum order of the polynomial to be
                considered for fitting.
        """
        warnings.simplefilter("ignore", np.RankWarning)

        def get_rms(x, y):
            return np.sqrt(np.sum((np.array(x) - np.array(y)) ** 2) / len(x))

        # list of (energy, volume) tuples
        e_v = list(zip(self.energies, self.volumes))
        ndata = len(e_v)
        # minimum number of data points used for fitting
        ndata_min = max(ndata - 2 * min_ndata_factor, min_poly_order + 1)
        rms_min = np.inf
        # number of data points available for fit in each iteration
        ndata_fit = ndata
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
        while (ndata_fit >= ndata_min) and (e_min in e_v_work):
            max_poly_order = ndata_fit - max_poly_order_factor
            e = [ei[0] for ei in e_v_work]
            v = [ei[1] for ei in e_v_work]
            # loop over polynomial order
            for idx in range(min_poly_order, max_poly_order + 1):
                coeffs = np.polyfit(v, e, idx)
                pder = np.polyder(coeffs)
                a = np.poly1d(pder)(v_before)
                b = np.poly1d(pder)(v_after)
                if a * b < 0:
                    rms = get_rms(e, np.poly1d(coeffs)(v))
                    rms_min = min(rms_min, rms * idx / ndata_fit)
                    all_coeffs[(idx, ndata_fit)] = [coeffs.tolist(), rms]
                    # store the fit coefficients small to large,
                    # i.e a0, a1, .. an
                    all_coeffs[(idx, ndata_fit)][0].reverse()
            # remove 1 data point from each end.
            e_v_work.pop()
            e_v_work.pop(0)
            ndata_fit = len(e_v_work)

        logger.info(f"total number of polynomials: {len(all_coeffs)}")

        norm = 0.0
        fit_poly_order = ndata
        # weight average polynomial coefficients.
        weighted_avg_coeffs = np.zeros((fit_poly_order,))

        # combine all the filtered polynomial candidates to get the final fit.
        for k, v in all_coeffs.items():
            # weighted rms = rms * polynomial order / rms_min / ndata_fit
            weighted_rms = v[1] * k[0] / rms_min / k[1]
            weight = np.exp(-(weighted_rms**2))
            norm += weight
            coeffs = np.array(v[0])
            # pad the coefficient array with zeros
            coeffs = np.lib.pad(coeffs, (0, max(fit_poly_order - len(coeffs), 0)), "constant")
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

    The following equations are supported::

        murnaghan: PRB 28, 5480 (1983)

        birch: Intermetallic compounds: Principles and Practice, Vol I:
            Principles. pages 195-210

        birch_murnaghan: PRB 70, 224107

        pourier_tarantola: PRB 70, 224107

        vinet: PRB 70, 224107

        deltafactor

        numerical_eos: 10.1103/PhysRevB.90.174107.

    Usage::

       eos = EOS(eos_name='murnaghan')
       eos_fit = eos.fit(volumes, energies)
       eos_fit.plot()
    """

    MODELS = dict(
        murnaghan=Murnaghan,
        birch=Birch,
        birch_murnaghan=BirchMurnaghan,
        pourier_tarantola=PourierTarantola,
        vinet=Vinet,
        deltafactor=DeltaFactor,
        numerical_eos=NumericalEOS,
    )

    def __init__(self, eos_name="murnaghan"):
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

    def fit(self, volumes, energies):
        """
        Fit energies as function of volumes.

        Args:
            volumes (list/np.array)
            energies (list/np.array)

        Returns:
            EOSBase: EOSBase object
        """
        eos_fit = self.model(np.array(volumes), np.array(energies))
        eos_fit.fit()
        return eos_fit


class EOSError(Exception):
    """Error class for EOS fitting."""
