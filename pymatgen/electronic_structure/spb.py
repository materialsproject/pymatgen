import scipy.integrate as integ
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

echarge = 1.602E-19
kB = 8.314 / 6.022E23
m0 = 9.10938188E-31
h = 6.626E-34
hbar = h / (2. * np.pi)
eps00 = 8.8541878176E-12


class SPBModel(object):

    def __init__(self):
        self.mass = None
        self.kay = 1.0
        self.mu0 = None
        
    def __getattribute__(self, attr):
        if attr.startswith('_'):
            return object.__getattribute__(self,attr)
        if callable(object.__getattribute__(self, attr)):
            return np.vectorize(object.__getattribute__(self, attr))
        else:
            return object.__getattribute__(self,attr)

    def n(self, eta, mass=None, T=300.):
        mass = mass or self.mass
        return 4 * np.pi * (2 * mass * m0 * kB * T / h ** 2) ** 1.5 * \
                fermiint(eta, 0.5) * 10**-6

    def rH(self, eta, kay=None):
        kay = kay or self.kay
        rH_tau = 1.5 / 2 * fermiint(eta, 0.5) * \
                fermiint(eta, -0.5) / (fermiint(eta, 0.)) ** 2
        rH_K = 3 * kay * (kay + 2) / (2 * kay + 1) ** 2  # anisotropy
        return rH_tau * rH_K

    def nH(self, eta, mass=None, T=300., kay=None, en=0.0, arH=0.0):
        """ cm-3 units """
        kay = kay or self.kay
        mass = mass or self.mass
        if not en:
            en = self.n(eta, mass, T)
        if not arH:
            arH = self.rH(eta, kay)
        return en / arH * 10**-6

    # mu0 = echarge / mass_conductivity/ m0 * (np.pi * hbar ** 4 * \
    # Cl * Nvs/(np.sqrt(2) * Edef ** 2. * (mass * m0 * kB * temp) ** 1.5)) \
    #* 100. ** 2; #cm**2/Vs  #Cl = vl**2 *density

    def muH(self, eta, mu0=None):
        mu0 = mu0 or self.mu0
        return mu0 * 0.5 * fermiint(eta, -0.5) / fermiint(eta, 0)

    def mu(self, eta, mu0=None):
        mu0 = mu0 or self.mu0
        return mu0 * 2. / 3 * fermiint(eta, 0.) / fermiint(eta, 0.5)

    def seebeck(self, eta):
        """ SPB Equation for Seebeck (units of uV /K)"""
        return kB / echarge * (2 * fermiint(eta, 1) / fermiint(eta, 0)\
                - eta) * 10 ** 6

    def lorenz(self, eta):
        return kB ** 2 / echarge ** 2 * \
                (3 * fermiint(eta, 0) * fermiint(eta, 2) - \
                2 ** 2 * (fermiint(eta, 1)) ** 2) \
                / (1 ** 2 * (fermiint(eta, 0)) ** 2)

    def etafromseebeck(self, seebeck):
        """ takes a value of seebeck and adjusts seebeck until it's equal
        Returns: eta (reduced chemical potential)"""
        return opt.fsolve(lambda x: (self.seebeck(x) - seebeck) ** 2, 1.0)[0]

    def etafromn(self, n, mass=None, T=300.):
        return opt.fsolve(lambda x: (np.log(self.n(x, mass, T)) -\
                                      np.log(n)) ** 2, 1.0)[0]

    def etafromnH(self, nH, mass=None, T=300., kay=None, en=0.0, arH=0.0):
        return opt.fsolve(lambda x: (np.log(self.nH(x, mass, T, kay, en, arH))\
                                      - np.log(nH)) ** 2, 1.0)[0]

    def massfromn(self, eta, n, T=300.):
        """  eta in kB*T units, n in cm-3, T in K
            returns mass in m0 units """
        return (n * 10 ** 6 / (4 * np.pi * fermiint(eta, 0.5))) ** (2. / 3)\
                / (2 * m0 * kB * T / h ** 2)
    def massfromseebn(self, seeb, n, T=300.):
        eta = self.etafromseebeck(seeb)
        mass = self.massfromn(eta, n, T)
        return mass

class Bands():
    def __init__(self,numbands):
        self.mass = [None] * numbands
        self.mu0 = [None] * numbands
        self.Edef = [None] * numbands
        self.kay = [1.0] * numbands
        self.results = {}
        self.bands = [None] * numbands
        ##Maybe should build in Edef at some point?
        ##Definitely  need to add in Nv
    def calczT(self, etas, T):
        self.results['seebeck'] = [None] * numbands
        self.results['n'] = [None] * numbands
        self.results['L'] = [None] * numbands
        self.results['mu'] = [None] * numbands
        self.results['rH'] = [None] * numbands
        self.results['nH'] = [None] * numbands
        self.results['muH'] = [None] * numbands
        self.results['RH'] = [None] * numbands
        self.results['kappaE'] = [None] * numbands
        self.results['all_seebeck'] = []
        self.results['all_sigma'] = []
        self.results['all_RH'] = []
        self.results['all_muH'] = []
        self.results['all_nH'] = []
        self.results['all_n'] = []
        self.results['all_kappaE_bipolar'] = []
        self.results['all_kappaE'] = []
        self.results['all_kappaE'] = []
        self.results['all_L'] = []
        for i in range(len(self.bands)):
            self.bands[i].mass = self.mass[i]
            self.bands[i].mu0 = self.mu0[i]
            self.bands[i].kay = self.kay[i]
            self.results['seebeck'][i] = self.bands[i].seebeck(etas)
            self.results['n'][i] = self.bands[i].n(etas, mass=None, T=T)#None is passed for mass, it just takes its own value that is set in self.bands[i].mass
            self.results['L'][i] = self.bands[i].lorenz(etas)
            self.results['mu'][i] = self.bands[i].L(etas, mu0=None)
            self.results['rH'][i] = self.bands[i].rH(etas)
            self.results['nH'][i] = self.results['n'][i] / self.results['rH'][i]
            self.results['muH'][i] = self.results['mu'][i] * self.results['rH'][i]
            self.results['sigma'][i] = self.results['mu'][i] * self.results['n'][i] * echarge/1000. #1/(mohm-cm)#cm^2/Vs * c * cm-3
            self.results['kappaE'][i] = self.results['sigma'][i]*1000*100 * self.results['L'][i] * T#w/m-k
            self.results['RH'][i] = 1./self.results['nH'][i] / echarge#nH = 1/RHe #cm^3/e
        self.results['all_seebeck'] = np.sum(np.array(self.results['seebeck'])*np.array(self.results['sigma']))
            
            
        
            
        

def fermi(eps, eta):
    """ eps is dimensionless energy (E/kbT), eta is the dimensionless
     chemical potential (mu/kBT)"""
    return (1. + np.exp(eps - eta)) ** -1.


def dfermi(eps, eta):
    return (1. / np.sqrt(np.exp(eps - eta)) + \
             np.sqrt(np.exp(eps - eta))) ** -2.


def fermiint(eta, j):
    """ Calculates the fermi energy for a vector of eta's, 
    uses the quad function, which seems to work the best"""
    try:
        a = np.array([])
        maxs = np.inf
        for etas in eta:
            calc = lambda t: t ** j * fermi(t - etas)
            a = np.append(a, integinf(calc))

    except TypeError:
        calc = lambda t: t ** j * fermi(t, eta)
        a = integinf(calc)
    return a

def integinf(func, epsabs=1E-12):
    """ Numerically integrates func from 0 to infinity """
    return integ.quad(func, 0., np.inf, epsrel=1.E-10,\
                       epsabs=epsabs, limlst=10)[0]


def makeintoarray(inputvar):
    """ Makes the input into array if it's just a single value, 
    necessary for 1D, doesn't work for 0D array (i think)"""
    try:
            inputvar.__iter__
            inputvar = np.array(inputvar) #this will mess up if you have a 0-d array..
    except AttributeError:
            inputvar = np.array([inputvar])
    return inputvar

if __name__ == "__main__":
    print "running it!"
    model = SPBModel()
    eta = np.linspace(-3, 10, 20)
