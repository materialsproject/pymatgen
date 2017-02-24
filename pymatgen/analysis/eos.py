# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

"""
This module implements of various equation of states.

Note: Except the numerical_eos, the rest were initially adapted from ASE and
deltafactor by @gmatteo.
"""

import logging
import warnings
import numpy as np
from scipy.optimize import leastsq, minimize
from monty.functools import return_none_if_raise
from copy import deepcopy

from pymatgen.core.units import FloatWithUnit, ArrayWithUnit, EnergyArray
from pymatgen.util.plotting_utils import get_publication_quality_plot

__author__ = "Guido Matteo, Kiran Mathew"
__credits__ = "Cormac Toher"

logger = logging.getLogger(__file__)


class EOSBase(object):

    def __init__(self, volumes, energies):
        self.volumes = volumes
        self.energies = energies

    def initial_guess(self):
        """
        Quadratic fit to get an initial guess for the parameters.

        Returns:
            tuple: (minimum energy(e0),
                buk modulus(b0),
                derivative of bulk modulus wrt pressure(b1),
                minimum volume(v0))
        """
        a, b, c = np.polyfit(self.volumes, self.energies, 2)
        self.eos_params = [a, b, c]

        v0 = -b / (2 * a)
        e0 = a * v0 ** 2 + b * v0 + c
        b0 = 2 * a * v0
        b1 = 4  # b1 is usually a small number like 4

        vmin, vmax = min(self.volumes), max(self.volumes)

        if not vmin < v0 and v0 < vmax:
            raise EOSError('The minimum volume of a fitted parabola is '
                           'not in the input volumes\n.')

        return e0, b0, b1, v0

    def fit(self):
        guess = self.initial_guess(self.volumes, self.energies)
        self.eos_params, ierr = leastsq(self.objective, guess, args=(self.volumes, self.energies))
        self._params = self.eos_params
        if ierr not in [1, 2, 3, 4]:
            raise EOSError("Optimal parameters not found")

    # the objective function that will be minimized
    def objective(self, pars, x, y):
        return y - self._func(x, *pars)

    def _func(self, x, E0, B0, B1, V0):
        pass

    def func(self, x):
        self._func(x, *self.eos_params)

    @property
    def e0(self):
        """
        Returns the min energy in eV.
        """
        return self._params[0]

    @property
    def v0(self):
        """
        Returns the minimum or the reference volume in Ang^3.
        """
        return self._params[3]

    @property
    def b0(self):
        """
        Returns the bulk modulus in eV/Ang^3
        """
        return self._params[1]

    @property
    def b0_GPa(self):
        """
        Returns the bulk modulus in GPa
        """
        return FloatWithUnit(self.b0, "eV ang^-3").to("GPa")

    @property
    def b1(self):
        """
        Returns the derivative of bulk modulus wrt pressure(dimensionless)
        """
        return self._params[2]

    @property
    @return_none_if_raise(AttributeError)
    def results(self):
        """
        Dictionary with the results. None if results are not available.

        Returns:
            dict
        """
        return dict(e0=self.e0, b0=self.b0, b1=self.b1, v0=self.v0)

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
        plt = get_publication_quality_plot(width=width, height=height, plt=plt,
                                           dpi=dpi)

        color = kwargs.get("color", "r")
        label = kwargs.get("label", "{} fit".format(self.eos_name))
        lines = ["Equation of State: %s" % self.eos_name,
                 "Minimum energy = %1.2f eV" % self.e0,
                 "Minimum or reference volume = %1.2f Ang^3" % self.v0,
                 "Bulk modulus = %1.2f eV/Ang^3 = %1.2f GPa" %
                 (self.b0, self.b0_GPa),
                 "Derivative of bulk modulus wrt pressure = %1.2f" % self.b1
                 ]
        text = "\n".join(lines)
        text = kwargs.get("text", text)

        # Plot input data.
        plt.plot(self.volumes, self.energies, linestyle="None", marker="o",
                 color=color)

        # Plot EOS.
        vmin, vmax = self.volumes.min(), self.volumes.max()
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        vfit = np.linspace(vmin, vmax, 100)

        plt.plot(vfit, self.func(vfit), linestyle="dashed", color=color, label=label)

        plt.grid(True)
        plt.xlabel("Volume $\AA^3$")
        plt.ylabel("Energy (eV)")
        plt.legend(loc="best", shadow=True)
        # Add text with fit parameters.
        plt.text(0.4, 0.5, text, transform=plt.gca().transAxes)

        return plt


class Quadratic(EOSBase):
    
    def _func(self, V, a, b, c):
        """
        Quadratic fit
        """
        return a * V**2 + b * V + c


class Murnaghan(EOSBase):

    def _func(self, V, E0, B0, B1, V0):
        """
        From PRB 28,5480 (1983)
        """
        return E0 + B0 * V / B1 * (((V0/V)**B1)/(B1-1)+1) - V0 * B0 / (B1-1)


class Birch(EOSBase):

    def _func(self, V, E0, B0, B1, V0):
        """
        From Intermetallic compounds: Principles and Practice, Vol. I: Principles
        Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos.
        case where n=0
        """
        return (E0
                + 9.0/8.0 * B0 * V0 * ((V0/V)**(2.0/3.0) - 1.0)**2
                + 9.0/16.0 * B0 * V0 * (B1-4.) * ((V0/V)**(2.0/3.0) - 1.0)**3)


class BirchMurnaghan(EOSBase):

    def _func(self, V, E0, B0, B1, V0):
        """
        BirchMurnaghan equation from PRB 70, 224107
        """
        eta = (V/V0)**(1./3.)
        return E0 + 9. * B0 * V0 / 16. * (eta**2-1)**2 * (6 + B1*(eta**2-1.) - 4. * eta**2)


class PourierTarantola(EOSBase):

    def _func(self, V, E0, B0, B1, V0):
        """
        Pourier-Tarantola equation from PRB 70, 224107
        """
        eta = (V/V0)**(1./3.)
        squiggle = -3.*np.log(eta)
        return E0 + B0 * V0 * squiggle**2 / 6. * (3. + squiggle * (B1 - 2))


class Vinet(EOSBase):

    def _func(self, V, E0, B0, B1, V0):
        """
        Vinet equation from PRB 70, 224107
        """
        eta = (V/V0)**(1./3.)
        return (E0 + 2. * B0 * V0 / (B1-1.)**2
                * (2. - (5. + 3. * B1 * (eta-1.) - 3.*eta)
                   * np.exp(-3. * (B1 - 1.) * (eta - 1.) / 2.)))


class PolynomialEOS(EOSBase):

    def func(self, volume):
         return np.poly1d(*self.eos_params)(volume)


class DeltaFactor(PolynomialEOS):

    def func(self, volume):
        return np.poly1d(*self.eos_params)(volume** (-2. / 3.))

    def fit(self):
        """
        This is the routine used to compute V0, B0, B1 in the deltafactor code.

        Args:
            volumes (list): list of volumes in Ang^3
            energies (list): list of energies in eV

        Returns:
            float, float, float, float, list: e0(min energy), b0(bulk modulus),
                b1(first derivative of bulk modulus), v0(min volume),
                eos_params(fit coefficients)
        """
        fitdata = np.polyfit(self.volumes**(-2./3.), self.energies, 3, full=True)

        deriv0 = np.poly1d(fitdata[0])
        deriv1 = np.polyder(deriv0, 1)
        deriv2 = np.polyder(deriv1, 1)
        deriv3 = np.polyder(deriv2, 1)

        for x in np.roots(deriv1):
            if x > 0 and deriv2(x) > 0:
                v0 = x**(-3./2.)
                break
        else:
            raise EOSError("No minimum could be found")

        derivV2 = 4./9. * x**5. * deriv2(x)
        derivV3 = (-20./9. * x**(13./2.) * deriv2(x) - 8./27. * x**(15./2.) * deriv3(x))
        b0 = derivV2 / x**(3./2.)
        b1 = -1 - x**(-3./2.) * derivV3 / derivV2

        # e0, b0, b1, v0, eos_params
        self._params = [np.poly1d(fitdata[0])(v0**(-2./3.)), b0, b1, v0]
        self.eos_params = fitdata[0]


class NumericalEOS(PolynomialEOS):

    def fit(self, min_ndata_factor=3, max_poly_order_factor=5,
                      min_poly_order=2):
        """
        Fit the input data to the 'numerical eos', the equation of state employed
        in the quasiharmonic Debye model described in the paper:
        10.1103/PhysRevB.90.174107.

        credits: Cormac Toher

        Args:
            volumes (list): list of volumes in Ang^3
            energies (list): list of energies in eV
            min_ndata_factor (int): parameter that controls the minimum number of
                data points that will be used for fitting.
                minimum number of data points = total data points-2*min_ndata_factor
            max_poly_order_factor (int): parameter that limits the max order of the
                polynomial used for fitting.
                max_poly_order = number of data points used for fitting - max_poly_order_factor
            min_poly_order (int): minimum order of the polynomial to be considered
                for fitting.

        Returns:
            (float, float, float, float, list): (
                min energy (eV),
                bulk modulus (eV/Ang^3),
                first derivative of bulk modulus wrt pressure (no unit),
                min volume (Ang^3),
                final fit coefficients)
        """
        warnings.simplefilter('ignore', np.RankWarning)

        get_rms = lambda x, y: np.sqrt(np.sum((np.array(x)-np.array(y))**2)/len(x))

        # list of (energy, volume) tuples
        e_v = [(i, j) for i, j in zip(self.energies, self.volumes)]
        ndata = len(e_v)
        # minimum number of data points used for fitting
        ndata_min = max(ndata - 2 * min_ndata_factor, min_poly_order + 1)
        rms_min = np.inf
        # number of data points available for fit in each iteration
        ndata_fit = ndata
        # store the fit polynomial coefficients and the rms in a dict,
        # where the key=(polynomial order, number of data points used for fitting)
        all_coeffs = {}

        # sort by energy
        e_v = sorted(e_v, key=lambda x: x[0])
        # minimum energy tuple
        e_min = e_v[0]
        # sort by volume
        e_v = sorted(e_v, key=lambda x: x[1])
        # index of minimum energy tuple in the volume sorted list
        emin_idx = e_v.index(e_min)
        # the volume lower than the volume corresponding to minimum energy
        v_before = e_v[emin_idx - 1][1]
        # the volume higher than the volume corresponding to minimum energy
        v_after = e_v[emin_idx + 1][1]
        e_v_work = deepcopy(e_v)

        # loop over the data points.
        while (ndata_fit >= ndata_min) and (e_min in e_v_work):
            max_poly_order = ndata_fit - max_poly_order_factor
            e = [ei[0] for ei in e_v_work]
            v = [ei[1] for ei in e_v_work]
            # loop over polynomial order
            for i in range(min_poly_order, max_poly_order + 1):
                coeffs = np.polyfit(v, e, i)
                pder = np.polyder(coeffs)
                a = np.poly1d(pder)(v_before)
                b = np.poly1d(pder)(v_after)
                if a * b < 0:
                    rms = get_rms(e, np.poly1d(coeffs)(v))
                    rms_min = min(rms_min, rms * i / ndata_fit)
                    all_coeffs[(i, ndata_fit)] = [coeffs.tolist(), rms]
                    # store the fit coefficients small to large, i.e a0, a1, .. an
                    all_coeffs[(i, ndata_fit)][0].reverse()
            # remove 1 data point from each end.
            e_v_work.pop()
            e_v_work.pop(0)
            ndata_fit = len(e_v_work)

        logger.info("total number of polynomials: {}".format(len(all_coeffs)))

        norm = 0.
        fit_poly_order = ndata
        # weight average polynomial coefficients.
        weighted_avg_coeffs = np.zeros((fit_poly_order,))

        # combine all the filtered polynomial candidates to get the final fit.
        for k, v in all_coeffs.items():
            # weighted rms = rms * polynomial order / rms_min / ndata_fit
            weighted_rms = v[1] * k[0] / rms_min / k[1]
            weight = np.exp(-(weighted_rms ** 2))
            norm += weight
            coeffs = np.array(v[0])
            # pad the coefficient array with zeros
            coeffs = np.lib.pad(coeffs, (0, max(fit_poly_order-len(coeffs), 0)),
                                'constant')
            weighted_avg_coeffs += weight * coeffs

        # normalization
        weighted_avg_coeffs /= norm
        weighted_avg_coeffs = weighted_avg_coeffs.tolist()
        # large to small(an, an-1, ..., a1, a0) as expected by np.poly1d
        weighted_avg_coeffs.reverse()
        fit_poly = np.poly1d(weighted_avg_coeffs)

        # evaluate e0, v0, b0 and b1
        min_wrt_v = minimize(fit_poly, e_min[1])
        e0, v0 = min_wrt_v.fun, min_wrt_v.x[0]
        pderiv2 = np.polyder(fit_poly, 2)
        pderiv3 = np.polyder(fit_poly, 3)
        b0 = v0 * np.poly1d(pderiv2)(v0)
        db0dv = np.poly1d(pderiv2)(v0) + v0 * np.poly1d(pderiv3)(v0)
        # db/dp
        b1 = - v0 * db0dv / b0

        self._params = [e0, b0, b1, v0]
        self.eos_params = weighted_avg_coeffs


class EOS(object):
    """
    Fit equation of state for bulk systems.

    The following equations are supported::

        quadratic:
            second order polynomial.

        murnaghan:
            PRB 28, 5480 (1983)

        birch:
            Intermetallic compounds: Principles and Practice,
            Vol I: Principles. pages 195-210

        birch_murnaghan:
            PRB 70, 224107

        pourier_tarantola:
            PRB 70, 224107

        vinet:
            PRB 70, 224107

        deltafactor

        numerical_eos:
            10.1103/PhysRevB.90.174107.

    Usage::

       eos = EOS(eos_name='murnaghan')
       fit = eos.fit(volumes, energies)
       print(fit)
       fit.plot()

    """

    MODELS = {
        "quadratic": Quadratic,
        "murnaghan": Murnaghan,
        "birch": Birch,
        "birch_murnaghan": BirchMurnaghan,
        "pourier_tarantola": PourierTarantola,
        "vinet": Vinet,
        "deltafactor": DeltaFactor,
        "numerical_eos": NumericalEOS
    }

    def __init__(self, eos_name='murnaghan'):
        if eos_name not in self.MODELS:
            raise EOSError(
                "The equation of state '{}' is not supported. "
                "Please choose one from the following list: {}".format(
                    eos_name, list(self.MODELS.keys())))
        self._eos_name = eos_name
        self.model = self.MODELS[eos_name]

    def fit(self, volumes, energies, vol_unit="ang^3", energy_unit="eV"):
        """
        Fit energies (in eV) as function of volumes (in Angstrom**3).

        Args:
            volumes (list/np.array)
            energies (list/np.array)
            vol_unit (str): volume units
            energy_unit (str): energy units

        Returns:
            EOSFit: EOSFit object that gives access to the optimal volume,
                the minumum energy, and the bulk modulus.

            Note: the units for the bulk modulus is eV/Angstrom^3.
        """
        # Convert volumes to Ang**3 and energies to eV (if needed).
        volumes = ArrayWithUnit(volumes, vol_unit).to("ang^3")
        energies = EnergyArray(energies, energy_unit).to("eV")
        eos = self.model(volumes, energies)
        eos.fit()
        return eos


class EOSError(Exception):
    """
    Exceptions raised by EOS.
    """
