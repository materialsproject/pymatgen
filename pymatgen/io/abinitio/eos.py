"Tools to compute equations of states with different models."
from __future__ import division, print_function

import numpy as np

__all__ = [
"EOS",
]

##########################################################################################

def murnaghan(V, E0, B0, BP, V0):
    'From PRB 28,5480 (1983)'

    E = E0 + B0*V/BP*(((V0/V)**BP)/(BP-1)+1) - V0*B0/(BP-1)
    return E

def birch(V, E0, B0, BP, V0):
    """
    From Intermetallic compounds: Principles and Practice, Vol. I: Principles
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos paper downloaded from Web

    case where n=0
    """

    E = (E0
         + 9.0/8.0*B0*V0*((V0/V)**(2.0/3.0) - 1.0)**2
         + 9.0/16.0*B0*V0*(BP-4.)*((V0/V)**(2.0/3.0) - 1.0)**3)
    return E

def birch_murnaghan(V, E0, B0, BP, V0):
    'BirchMurnaghan equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)
    E = E0 + 9.*B0*V0/16.*(eta**2-1)**2*(6 + BP*(eta**2-1.) - 4.*eta**2)
    return E

def pourier_tarantola(V, E0, B0, BP, V0):
    'Pourier-Tarantola equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)
    squiggle = -3.*np.log(eta)

    E = E0 + B0*V0*squiggle**2/6.*(3. + squiggle*(BP - 2))
    return E

def vinet(V, E0, B0, BP, V0):
    'Vinet equation from PRB 70, 224107'

    eta = (V/V0)**(1./3.)

    E = (E0 + 2.*B0*V0/(BP-1.)**2
         * (2. - (5. +3.*BP*(eta-1.)-3.*eta)*np.exp(-3.*(BP-1.)*(eta-1.)/2.)))
    return E

##########################################################################################
class EOSError(Exception):
    "Exceptions raised by EOS"

class EOS(object):
    """
    Fit equation of state for bulk systems.

    The following equation is used::

       murnaghan
           PRB 28, 5480 (1983)

       birch
           Intermetallic compounds: Principles and Practice,
           Vol I: Principles. pages 195-210

       birchmurnaghan
           PRB 70, 224107

       pouriertarantola
           PRB 70, 224107

       vinet
           PRB 70, 224107

    Use::

       eos = EOS(eos='murnaghan')
       fit = eos.fit(volumes, energies)
       print fit
       fit.plot()

    """
    Error = EOSError

    #: Models available.
    functions = {
        "murnaghan"        : murnaghan,
        "birch"            : birch,
        "birch_murnaghan"  : birch_murnaghan,
        "pourier_tarantola": pourier_tarantola,
        "vinet"            : vinet,
    }

    def __init__(self, eos_name='murnaghan'):
        self._func = self.functions[eos_name]

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

    def fit(self, volumes, energies):
        """
        Fit energies [eV] as function of volumes [Angstrom^3].

        Returns EosFit instance that gives access to the optimal volume, 
        the minumum energy, and the bulk modulus.  
        Notice that the units for the bulk modulus is eV/Angstrom^3.
        """
        return EOS_Fit(volumes, energies, self._func)

##########################################################################################


class EOS_Fit(object):
    "Performs the fit of E(V) and provides method to access the results of the fit."

    def __init__(self, volumes, energies, func):
        """ 
        args:
            energies: list of energies in eV 
            volumes: list of volumes in Angstrom^3
            func: callable function 
        """
        self.volumes  = np.array(volumes) 
        self.energies = np.array(energies)
        self.func     = func
        self.exceptions = []

        # objective function that will be minimized
        def objective(pars, x, y):
            return y - self.func(x, *pars)

        # quadratic fit to get an initial guess for the parameters 
        a,b,c = np.polyfit(volumes, energies, 2) 
                                           
        v0 = -b/(2*a)
        e0 = a*v0**2 + b*v0 + c
        b0 = 2*a*v0
        bP = 4  # Bp is usually a small number like 4

        vmin, vmax = self.volumes.min(), self.volumes.max()

        if not (vmin < v0 and v0 < vmax):
            exc = EOSError('The minimum volume of a fitted parabola is not in the input volumes\n.')
            self.exceptions.append(exc)
            print(str(exc))

        # Initial guesses for the parameters
        self.p0 = [e0, b0, bP, v0] 

        from scipy.optimize import leastsq
        self.eos_params, self.ierr = leastsq(objective, self.p0, args=(volumes, energies)) 

        if self.ierr not in [1,2,3,4]:
            exc = EOSError("Optimal parameters not found")
            self.exceptions.append(exc)
            raise exc

        self.e0 = self.eos_params[0]
        self.b  = self.eos_params[1]
        self.bp = self.eos_params[2]
        self.v0 = self.eos_params[3]

        #TODO: Add support for rich comparison so that we can define an order
        # based on the accuracy of the fit.

    def __str__(self):
        lines = []
        app = lines.append
        app("Equation of State: %s" % self.name)
        app("Minimum volume = %1.2f Ang^3" % self.v0)
        app("Bulk modulus = %1.2f eV/Ang^3 = %1.2f GPa, bp = %1.2f" % (self.b, self.b*160.21773, self.bp))
        return "\n".join(lines)

    @property
    def name(self):
        return self.func.__name__

    def plot(self, show=True, savefig=None):
        """
        Uses Matplotlib to plot the energy curve.  

        Args:
            show: 
                True to show the figure 
            savefig:
                'abc.png' or 'abc.eps' to save the figure to a file.
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
        line, = ax.plot(vfit, self.func(vfit, *self.eos_params) ,"b-")

        lines.append(line)
        legends.append(self.name + ' fit')
                                                                                                                                         
        # Set xticks and labels.
        ax.grid(True)
        ax.set_xlabel("Volume $\AA^3$")
        ax.set_ylabel("Energy (eV)")
        ax.legend(lines, legends, 'upper right', shadow=True)
                                                                                                                                         
        fig.text(0.4,0.5,'Min volume = %1.2f $\AA^3$' % self.v0, transform = ax.transAxes)
        fig.text(0.4,0.4,'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (self.b, self.b*160.21773), transform = ax.transAxes)
        fig.text(0.4,0.3,'Bp = %1.2f' % self.bp, transform = ax.transAxes)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

##########################################################################################
