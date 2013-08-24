"""Tools to compute equations of states with different models."""
from __future__ import division, print_function

import collections
import numpy as np

import pymatgen.core.physical_constants as const
import pymatgen.core.units as units

__all__ = [
    "EOS",
]

__version__ = "1.0"

##########################################################################################

def murnaghan(V, E0, B0, B1, V0):
    """From PRB 28,5480 (1983)"""

    E = E0 + B0*V/B1*(((V0/V)**B1)/(B1-1)+1) - V0*B0/(B1-1)
    return E

def birch(V, E0, B0, B1, V0):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Principles
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos paper downloaded from Web

    case where n=0
    """

    E = (E0
         + 9.0/8.0*B0*V0*((V0/V)**(2.0/3.0) - 1.0)**2
         + 9.0/16.0*B0*V0*(B1-4.)*((V0/V)**(2.0/3.0) - 1.0)**3)
    return E

def birch_murnaghan(V, E0, B0, B1, V0):
    """BirchMurnaghan equation from PRB 70, 224107"""

    eta = (V/V0)**(1./3.)
    E = E0 + 9.*B0*V0/16.*(eta**2-1)**2*(6 + B1*(eta**2-1.) - 4.*eta**2)
    return E

def pourier_tarantola(V, E0, B0, B1, V0):
    """Pourier-Tarantola equation from PRB 70, 224107"""

    eta = (V/V0)**(1./3.)
    squiggle = -3.*np.log(eta)

    E = E0 + B0*V0*squiggle**2/6.*(3. + squiggle*(B1 - 2))
    return E

def vinet(V, E0, B0, B1, V0):
    'Vinet equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)

    E = (E0 + 2.*B0*V0/(B1-1.)**2
         * (2. - (5. +3.*B1*(eta-1.)-3.*eta)*np.exp(-3.*(B1-1.)*(eta-1.)/2.)))
    return E


def deltafactor_polyfit(volumes, energies):
    """
    This is the routine used to compute V0, B0, B1 in the deltafactor code.

    Taken from deltafactor/eosfit.py
    """
    fitdata = np.polyfit(volumes**(-2./3.), energies, 3, full=True)
    ssr = fitdata[1]
    sst = np.sum((energies - np.average(energies))**2.)
    residuals0 = ssr/sst
    deriv0 = np.poly1d(fitdata[0])
    deriv1 = np.polyder(deriv0, 1)
    deriv2 = np.polyder(deriv1, 1)
    deriv3 = np.polyder(deriv2, 1)

    v0 = 0
    x = 0
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

    n = collections.namedtuple("DeltaFitResults", "v0 b0 b1 poly1d")
    return n(v0, b0, b1, fitdata[0])

##########################################################################################


class EOSError(Exception):
    """Exceptions raised by EOS."""


class EOS(object):
    """
    Fit equation of state for bulk systems.

    The following equation is used::

       murnaghan
           PRB 28, 5480 (1983)

       birch
           Intermetallic compounds: Principles and Practice, Vol I: Principles. pages 195-210

       birchmurnaghan
           PRB 70, 224107

       pouriertarantola
           PRB 70, 224107

       vinet
           PRB 70, 224107

    Use::

       eos = EOS(eos_name='murnaghan')
       fit = eos.fit(volumes, energies)
       print(fit)
       fit.plot()

    """
    Error = EOSError

    #: Models available.
    MODELS = {
        "murnaghan"        : murnaghan,
        "birch"            : birch,
        "birch_murnaghan"  : birch_murnaghan,
        "pourier_tarantola": pourier_tarantola,
        "vinet"            : vinet,
        "deltafactor"      : deltafactor_polyfit,
    }

    def __init__(self, eos_name='murnaghan'):
        self._eos_name = eos_name
        self._func = self.MODELS[eos_name]

    @staticmethod
    def Murnaghan():
        return EOS(eos_name='murnaghan')

    @staticmethod
    def Birch():
        return EOS(eos_name='birch')

    @staticmethod
    def Birch_Murnaghan():
        return EOS(eos_name='birch_murnaghan')

    @staticmethod
    def Pourier_Tarantola():
        return EOS(eos_name='pourier_tarantola')

    @staticmethod
    def Vinet():
        return EOS(eos_name='vinet')

    @staticmethod
    def DeltaFactor():
        return EOS(eos_name='deltafactor')

    def fit(self, volumes, energies, len_units="Ang", ene_units="eV"):
        """
        Fit energies [eV] as function of volumes [Angstrom^3].

        Returns EosFit instance that gives access to the optimal volume, 
        the minumum energy, and the bulk modulus.  
        Notice that the units for the bulk modulus is eV/Angstrom^3.
        """
        # Convert volumes to Ang**3 and energies to eV (if needed).
        volumes = units.any2Ang3(len_units)(volumes)
        energies = units.any2eV(ene_units)(energies)

        return EOS_Fit(volumes, energies, self._func, self._eos_name)

##########################################################################################


class EOS_Fit(object):
    """Performs the fit of E(V) and provides method to access the results of the fit."""

    def __init__(self, volumes, energies, func, eos_name):
        """ 
        args:
            energies: list of energies in eV 
            volumes: list of volumes in Angstrom^3
            func: callable function 
        """
        self.volumes  = np.array(volumes) 
        self.energies = np.array(energies)
        assert len(self.volumes) == len(self.energies)

        self.func     = func
        self.eos_name = eos_name
        self.exceptions = []

        if eos_name == "deltafactor":

            self.ierr = 0
            try:
                results = deltafactor_polyfit(self.volumes, self.energies)

                self.e0 = None
                self.v0 = results.v0
                self.b0 = results.b0
                self.b1 = results.b1
                self.p0 = results.poly1d
                self.eos_params = results.poly1d

            except EOSError as exc:
                self.ierr = 1
                print(str(exc))
                self.exceptions.append(exc)
                raise 

        else:
            # Objective function that will be minimized
            def objective(pars, x, y):
                return y - self.func(x, *pars)

            # Quadratic fit to get an initial guess for the parameters 
            a,b,c = np.polyfit(self.volumes, self.energies, 2) 
                                               
            v0 = -b/(2*a)
            e0 = a*v0**2 + b*v0 + c
            b0 = 2*a*v0
            b1 = 4  # b1 is usually a small number like 4

            vmin, vmax = self.volumes.min(), self.volumes.max()

            if not (vmin < v0 and v0 < vmax):
                exc = EOSError('The minimum volume of a fitted parabola is not in the input volumes\n.')
                print(str(exc))
                self.exceptions.append(exc)

            # Initial guesses for the parameters
            self.p0 = [e0, b0, b1, v0] 

            from scipy.optimize import leastsq
            self.eos_params, self.ierr = leastsq(objective, self.p0, args=(self.volumes, self.energies)) 

            if self.ierr not in [1,2,3,4]:
                exc = EOSError("Optimal parameters not found")
                print(str(exc))
                self.exceptions.append(exc)
                raise exc

            self.e0 = self.eos_params[0]
            self.b0 = self.eos_params[1]
            self.b1 = self.eos_params[2]
            self.v0 = self.eos_params[3]

    def __str__(self):
        lines = []
        app = lines.append
        app("Equation of State: %s" % self.name)
        app("Minimum volume = %1.2f Ang^3" % self.v0)
        app("Bulk modulus = %1.2f eV/Ang^3 = %1.2f GPa, b1 = %1.2f" % (self.b0, self.b0_GPa, self.b1))
        return "\n".join(lines)

    @property
    def name(self):
        return self.func.__name__

    @property
    def b0_GPa(self):
        return self.b0 * const.eVA3_GPa

    def plot(self, **kwargs):
        """
        Uses Matplotlib to plot the energy curve.  

        Args:
            show: 
                True to show the figure 
            savefig:
                'abc.png' or 'abc.eps' to save the figure to a file.

        Returns:
            Matplotlib figure.
        """
        import matplotlib.pyplot as plt
                                                                                                                                         
        vmin, vmax = self.volumes.min(), self.volumes.max()
        emin, emax = self.energies.min(), self.energies.max()
                                                                        
        vmin, vmax = (vmin - 0.01 * abs(vmin), vmax + 0.01 * abs(vmax))
        emin, emax = (emin - 0.01 * abs(emin), emax + 0.01 * abs(emax))
                                                                                                                                         
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        lines, legends = [], []

        # Plot input data.
        line, = ax.plot(self.volumes, self.energies, "ro")
        lines.append(line)
        legends.append("Input Data")

        # Plot EOS.
        vfit = np.linspace(vmin, vmax, 100)

        if self.eos_name == "deltafactor":
            xx = vfit**(-2./3.)
            line, = ax.plot(vfit, np.polyval(self.eos_params, xx) ,"b-")
        else:
            line, = ax.plot(vfit, self.func(vfit, *self.eos_params) ,"b-")

        lines.append(line)
        legends.append(self.name + ' fit')
                                                                                                                                         
        # Set xticks and labels.
        ax.grid(True)
        ax.set_xlabel("Volume $\AA^3$")
        ax.set_ylabel("Energy (eV)")
        ax.legend(lines, legends, 'upper right', shadow=True)

        # Add text with fit parameters.
        text = []; app = text.append
        app("Min Volume = %1.2f $\AA^3$" % self.v0)
        app("Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa" % (self.b0, self.b0_GPa))
        app("B1 = %1.2f" % self.b1)
        fig.text(0.4, 0.5, "\n".join(text), transform=ax.transAxes)

        if kwargs.pop("show", True):
            plt.show()

        savefig = kwargs.pop("savefig", None)
        if savefig is not None:
            fig.savefig(savefig)

        return fig

##########################################################################################
