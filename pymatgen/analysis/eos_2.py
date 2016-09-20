# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

"""
Tools to compute equations of states with different models.

Adapted from ASE and deltafactor.
"""

import logging
import numpy as np

from monty.functools import return_none_if_raise

import pymatgen.core.units as units
from pymatgen.core.units import FloatWithUnit
from pymatgen.util.plotting_utils import get_publication_quality_plot

logger = logging.getLogger(__file__)


def quadratic(V, a, b, c):
    """
    Quadratic fit
    """
    return a * V**2 + b * V + c


def murnaghan(V, E0, B0, B1, V0):
    """
    From PRB 28,5480 (1983)
    """
    return E0 + B0 * V / B1 * (((V0/V)**B1)/(B1-1)+1) - V0 * B0 / (B1-1)


def birch(V, E0, B0, B1, V0):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Principles
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos.
    case where n=0
    """
    return (E0
            + 9.0/8.0 * B0 * V0 * ((V0/V)**(2.0/3.0) - 1.0)**2
            + 9.0/16.0 * B0 * V0 * (B1-4.) * ((V0/V)**(2.0/3.0) - 1.0)**3)


def birch_murnaghan(V, E0, B0, B1, V0):
    """
    BirchMurnaghan equation from PRB 70, 224107
    """
    eta = (V/V0)**(1./3.)
    return E0 + 9. * B0 * V0 / 16. * (eta**2-1)**2 * (6 + B1*(eta**2-1.) - 4. * eta**2)


def pourier_tarantola(V, E0, B0, B1, V0):
    """
    Pourier-Tarantola equation from PRB 70, 224107
    """
    eta = (V/V0)**(1./3.)
    squiggle = -3.*np.log(eta)
    return E0 + B0 * V0 * squiggle**2 / 6. * (3. + squiggle * (B1 - 2))


def vinet(V, E0, B0, B1, V0):
    """
    Vinet equation from PRB 70, 224107
    """
    eta = (V/V0)**(1./3.)
    return (E0 + 2. * B0 * V0 / (B1-1.)**2
            * (2. - (5. + 3. * B1 * (eta-1.) - 3.*eta)
               * np.exp(-3. * (B1 - 1.) * (eta - 1.) / 2.)))


def deltafactor_polyfit(volumes, energies):
    """
    This is the routine used to compute V0, B0, B1 in the deltafactor code.
    """
    fitdata = np.polyfit(volumes**(-2./3.), energies, 3, full=True)

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
    return np.poly1d(fitdata[0])(v0**(-2./3.)), b0, b1, v0, fitdata[0]


class EOS(object):
    """
    Fit equation of state for bulk systems.

    The following equations are supported::

       murnaghan:
           PRB 28, 5480 (1983)

       birch:
           Intermetallic compounds: Principles and Practice, Vol I: Principles.
           pages 195-210

       birch_murnaghan:
           PRB 70, 224107

       pourier_tarantola:
           PRB 70, 224107

       vinet:
           PRB 70, 224107

    Use::

       eos = EOS(eos_name='murnaghan')
       fit = eos.fit(volumes, energies)
       print(fit)
       fit.plot()

    """

    MODELS = {
        "quadratic": quadratic,
        "murnaghan": murnaghan,
        "birch": birch,
        "birch_murnaghan": birch_murnaghan,
        "pourier_tarantola": pourier_tarantola,
        "vinet": vinet,
        "deltafactor": deltafactor_polyfit
    }

    def __init__(self, eos_name='murnaghan'):
        if eos_name not in self.MODELS:
            raise KeyError("The equation of state '{}' is not supported. "
                           "Please choose one from the following list: {}".
                             format(eos_name, list(self.MODELS.keys())))
        self._eos_name = eos_name
        self._func = self.MODELS[eos_name]

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
        volumes = units.ArrayWithUnit(volumes, vol_unit).to("ang^3")
        energies = units.EnergyArray(energies, energy_unit).to("eV")

        return EOSFit(volumes, energies, self._func, self._eos_name)


class EOSFit(object):
    """
    Performs the fit of E(V) and provides method to access the results of
    the fit.
    """

    def __init__(self, volumes, energies, func, eos_name):
        """
        Args:
            energies (list): energies in eV
            volumes (list): volumes in Angstrom^3
            func : callable function

        Attributes:
            eos_name: EOS name
            func: the equation of state function
            params: tuple(
                minimum energy(e0),
                buk modulus(b0),
                derivative of bulk modulus wrt pressure(b1),
                minimum or the reference volume(v0) )
            eos_params: the results returned by the fit.
        """
        assert len(volumes) == len(energies)

        self.volumes = np.array(volumes)
        self.energies = np.array(energies)

        self.func = func
        self.eos_name = eos_name
        self.params = None

        # initial guess: parabola
        guess = self.initial_guess()

        if self.eos_name in ["quadratic"]:
            self.params = guess
            logger.info("The initial guess is quadratic")
        else:
            self.params = self.fit(guess)

    def initial_guess(self):
        """
        Quadratic fit to get an initial guess for the parameters.

        Returns:
            tuple: (minimum energy(e0), buk modulus(b0),
                derivative of bulk modulus wrt pressure(b1),
                minimum volume(v0))
        """
        a, b, c = np.polyfit(self.volumes, self.energies, 2)
        self.eos_params = [a, b, c]

        v0 = -b / (2 * a)
        e0 = a * v0 ** 2 + b * v0 + c
        b0 = 2 * a * v0
        b1 = 4  # b1 is usually a small number like 4

        vmin, vmax = self.volumes.min(), self.volumes.max()

        if not vmin < v0 and v0 < vmax:
            raise EOSError('The minimum volume of a fitted parabola is '
                           'not in the input volumes\n.')

        return e0, b0, b1, v0

    def fit(self, guess):
        """
        Do the fitting

        Args:
            guess (tuple): initial guess for e0, b0, b1, v0(in that order)

        Returns:
            tuple: (minimum energy(e0), buk modulus(b0),
                derivative of bulk modulus wrt pressure(b1),
                reference volume(v0))
        """
        if self.eos_name == "deltafactor":
            try:
                results = deltafactor_polyfit(self.volumes, self.energies)
                self.eos_params = results[-1]
            except:
                raise EOSError()

        else:
            # the objective function that will be minimized
            def objective(pars, x, y):
                return y - self.func(x, *pars)

            from scipy.optimize import leastsq
            results, ierr = \
                leastsq(objective, guess, args=(self.volumes, self.energies))
            self.eos_params = results

            if ierr not in [1, 2, 3, 4]:
                raise EOSError("Optimal parameters not found")

        # e0, b0, b1, v0
        return results[0], results[1], results[2], results[3]

    def __str__(self):
        lines = ["Equation of State: %s" % self.name,
                 "Minimum energy = %1.2f eV" % self.params[0],
                 "Minimum or reference volume = %1.2f Ang^3" % self.params[3],
                 "Bulk modulus = %1.2f eV/Ang^3 = %1.2f GPa" %
                 (self.params[1], self.b0_GPa),
                 "Derivative of bulk modulus wrt pressure = %1.2f" %
                 self.params[2]
                 ]
        return "\n".join(lines)

    @property
    def name(self):
        """
        Returns:
            str: EOS name
        """
        return self.func.__name__

    @property
    def b0_GPa(self):
        """
        Returns the bulk modulus in GPa
        """
        return FloatWithUnit(self.params[1], "eV ang^-3").to("GPa")

    @property
    @return_none_if_raise(AttributeError)
    def results(self):
        """
        Dictionary with the results. None if results are not available.

        Returns:
            dict
        """
        return dict(e0=self.params[0], b0=self.params[1], b1=self.params[2],
                    v0=self.params[3])

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
        label = kwargs.get("label", "{} fit".format(self.name))
        text = kwargs.get("text", None)

        # Plot input data.
        plt.plot(self.volumes, self.energies, linestyle="None", marker="o",
                 color=color)

        # Plot EOS.
        vmin, vmax = self.volumes.min(), self.volumes.max()
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        vfit = np.linspace(vmin, vmax, 100)

        if self.eos_name == "deltafactor":
            plt.plot(vfit, np.polyval(self.eos_params, vfit**(-2./3.)),
                     linestyle="dashed", color=color, label=label)
        else:
            plt.plot(vfit, self.func(vfit, *self.eos_params),
                     linestyle="dashed", color=color, label=label)

        plt.grid(True)
        plt.xlabel("Volume $\AA^3$")
        plt.ylabel("Energy (eV)")
        plt.legend(loc="best", shadow=True)
        # Add text with fit parameters.
        if not text:
            plt.text(0.4, 0.5, str(self), transform=plt.gca().transAxes)

        return plt


class EOSError(Exception):
    """
    Exceptions raised by EOS.
    """
