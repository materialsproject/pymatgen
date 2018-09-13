import numpy as np
import matplotlib.pyplot as plt
from monty.serialization import dumpfn
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.bandstructure import \
    BandStructureSymmLine, Kpoint, Spin
from pymatgen.io.vasp import Vasprun
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.electronic_structure.dos import Dos, CompleteDos, Orbital
from pymatgen.electronic_structure.boltztrap import BoltztrapError
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

try:
    from BoltzTraP2 import sphere
    from BoltzTraP2 import fite
    from BoltzTraP2 import bandlib as BL
    from BoltzTraP2 import units
except ImportError:
    raise BoltztrapError("BoltzTraP2 has to be installed and working")

"""
BoltzTraT2 is a python software interpolating band structures and
computing materials properties from dft band structure using Boltzmann
semi-classical transport theory.
This module provides a pymatgen interface to BoltzTraT2. 
Some of the code is written following the examples provided in BoltzTraP2

BoltzTraT2 has been developed by Georg Madsen, Jesús Carrete, Matthieu J. Verstraete.

https://gitlab.com/sousaw/BoltzTraP2
https://www.sciencedirect.com/science/article/pii/S0010465518301632

References are:

    Georg K.H.Madsen, Jesús Carrete, Matthieu J.Verstraete;
    BoltzTraP2, a program for interpolating band structures and 
    calculating semi-classical transport coefficients
    Computer Physics Communications 231, 140-145, 2018
    
    Madsen, G. K. H., and Singh, D. J. (2006).
    BoltzTraP. A code for calculating band-structure dependent quantities.
    Computer Physics Communications, 175, 67-71
    
TODO:
- spin polarized bands
- read first derivative of the eigenvalues from vasprun.xml (mommat)
- handle magnetic moments (magmom)
"""

__author__ = "Francesco Ricci"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Fracesco Ricci"
__email__ = "frankyricci@gmail.com"
__status__ = "Development"
__date__ = "August 2018"


class BandstructureLoader:
    """Loader for Bandsstrcture object"""

    def __init__(self, bs_obj, structure=None, nelect=None):
        """Structure and nelect is needed to be provide"""
        self.kpoints = np.array([kp.frac_coords for kp in bs_obj.kpoints])

        if structure is None:
            try:
                self.structure = bs_obj.structure
            except:
                BaseException('No structure found in the bs obj.')
        else:
            self.structure = structure

        self.atoms = AseAtomsAdaptor.get_atoms(self.structure)

        if len(bs_obj.bands) == 1:
            e = list(bs_obj.bands.values())[0]
            self.ebands = e * units.eV
            self.dosweight = 2.0
        elif len(bs_obj.bands) == 2:
            raise BaseException("spin bs case not implemented")

        self.lattvec = self.atoms.get_cell().T * units.Angstrom
        self.mommat = None
        self.magmom = None
        self.fermi = bs_obj.efermi * units.eV

        self.nelect = nelect
        self.UCvol = self.structure.volume * units.Angstrom ** 3

        self.vbm_idx = list(bs_obj.get_vbm()['band_index'].values())[0][-1]
        self.cbm_idx = list(bs_obj.get_cbm()['band_index'].values())[0][0]

    def get_lattvec(self):
        try:
            self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec

    def bandana(self, emin=-np.inf, emax=np.inf):
        """Cut out bands outside the range (emin,emax)"""
        bandmin = np.min(self.ebands, axis=1)
        bandmax = np.max(self.ebands, axis=1)
        ii = np.nonzero(bandmin < emax)
        nemax = ii[0][-1]
        ii = np.nonzero(bandmax > emin)
        nemin = ii[0][0]
        # BoltzTraP2.misc.info("BANDANA output")
        # for iband in range(len(self.ebands)):
        # BoltzTraP2.misc.info(iband, bandmin[iband], bandmax[iband], (
        # (bandmin[iband] < emax) & (bandmax[iband] > emin)))
        self.ebands = self.ebands[nemin:nemax]
        if self.mommat is not None:
            self.mommat = self.mommat[:, nemin:nemax, :]
        # Removing bands may change the number of valence electrons
        if self.nelect is not None:
            self.nelect -= self.dosweight * nemin
        return nemin, nemax

    def get_volume(self):
        try:
            self.UCvol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol


class VasprunLoader:
    """Loader for Vasprun object"""

    def __init__(self, vrun_obj=None):
        if vrun_obj:
            self.kpoints = np.array(vrun_obj.actual_kpoints)
            self.structure = vrun_obj.final_structure
            self.atoms = AseAtomsAdaptor.get_atoms(self.structure)
            self.proj = []
            if len(vrun_obj.eigenvalues) == 1:
                e = list(vrun_obj.eigenvalues.values())[0]
                self.ebands = e[:, :, 0].transpose() * units.eV
                self.dosweight = 2.0
                if vrun_obj.projected_eigenvalues:
                    self.proj = list(vrun_obj.projected_eigenvalues.values())[0]

            elif len(vrun_obj.eigenvalues) == 2:
                raise BoltztrapError("spin bs case not implemented")

            self.lattvec = self.atoms.get_cell().T * units.Angstrom

            # TODO: read mommat from vasprun
            self.mommat = None
            self.magmom = None
            self.fermi = vrun_obj.efermi * units.eV
            self.nelect = vrun_obj.parameters['NELECT']
            self.UCvol = self.structure.volume * units.Angstrom ** 3

    def from_file(self, vasprun_file):
        """Get a vasprun.xml file and return a VasprunLoader"""
        vrun_obj = Vasprun(vasprun_file, parse_projected_eigen=True)
        return VasprunLoader(vrun_obj)

    def get_lattvec(self):
        try:
            self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec

    def bandana(self, emin=-np.inf, emax=np.inf):
        """Cut out bands outside the range (emin,emax)"""
        bandmin = np.min(self.ebands, axis=1)
        bandmax = np.max(self.ebands, axis=1)
        ii = np.nonzero(bandmin < emax)
        nemax = ii[0][-1]
        ii = np.nonzero(bandmax > emin)
        nemin = ii[0][0]
        # BoltzTraP2.misc.info("BANDANA output")
        # for iband in range(len(self.ebands)):
        # BoltzTraP2.misc.info(iband, bandmin[iband], bandmax[iband], (
        # (bandmin[iband] < emax) & (bandmax[iband] > emin)))
        self.ebands = self.ebands[nemin:nemax]

        if isinstance(self.proj, np.ndarray):
            self.proj = self.proj[:,nemin:nemax,:,:]
            
        if self.mommat is not None:
            self.mommat = self.mommat[:, nemin:nemax, :]
        # Removing bands may change the number of valence electrons
        if self.nelect is not None:
            self.nelect -= self.dosweight * nemin
        return nemin, nemax

    def get_volume(self):
        try:
            self.UCvol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol


class BztInterpolator(object):
    """
        Interpolate the dft band structures
        
        Args:
            data: A loader
            lpfac: the number of interpolation points in the real space. By
                default 10 gives 10 time more points in the real space than
                the number of kpoints given in reciprocal space.
            energy_range: usually the interpolation is not needed on the entire energy
                range but on a specific range around the fermi level.
                This energy in eV fix the range around the fermi level (E_fermi-energy_range,E_fermi+energy_range) of
                bands that will be interpolated
                and taken into account to calculate the transport properties.
            curvature: boolean value to enable/disable the calculation of second
                derivative related trasport properties (Hall coefficient).
        Example:
            data = VasprunLoader().from_file('vasprun.xml')
            bztInterp = BztInterpolator(data)
    """

    def __init__(self, data, lpfac=10, energy_range=1.5, curvature=True):

        self.data = data
        num_kpts = self.data.kpoints.shape[0]
        self.efermi = self.data.fermi
        self.nemin, self.nemax = self.data.bandana(emin=self.efermi - (energy_range * units.eV),
                                                   emax=self.efermi + (energy_range * units.eV))
        self.equivalences = sphere.get_equivalences(self.data.atoms, self.data.magmom,
                                                    num_kpts * lpfac)
        self.coeffs = fite.fitde3D(self.data, self.equivalences)
        self.eband, self.vvband, self.cband = fite.getBTPbands(self.equivalences,
                                                               self.coeffs, self.data.lattvec,
                                                               curvature=curvature)

    def get_band_structure(self):
        """Return a BandStructureSymmLine object interpolating bands along a
        High symmetry path calculated from the structure using HighSymmKpath function"""

        kpath = HighSymmKpath(self.data.structure)
        kpt_line = [Kpoint(k, self.data.structure.lattice) for k
                    in
                    kpath.get_kpoints(coords_are_cartesian=False)[
                        0]]
        kpoints = np.array(
            [kp.frac_coords for kp in kpt_line])

        labels_dict = {l: k for k, l in zip(
            *kpath.get_kpoints(coords_are_cartesian=False)) if l}

        lattvec = self.data.get_lattvec()
        egrid, vgrid = fite.getBands(kpoints, self.equivalences, lattvec, self.coeffs)

        bands_dict = {Spin.up: (egrid / units.eV)}

        sbs = BandStructureSymmLine(kpoints, bands_dict,
                                    self.data.structure.lattice.reciprocal_lattice,
                                    self.efermi / units.eV,
                                    labels_dict=labels_dict)
        return sbs

    def get_dos(self, partial_dos=False, npts_mu=10000, T=None):
        """
            Return a Dos object interpolating bands
            
            Args:
                partial_dos: if True, projections will be interpolated as well
                    and partial doses will be return. Projections must be available
                    in the loader.
                npts_mu: number of energy points of the Dos
                T: parameter used to smooth the Dos
        """
        energies, densities, vvdos, cdos = BL.BTPDOS(self.eband, self.vvband, npts=npts_mu)
        if T is not None:
            densities = BL.smoothen_DOS(energies, densities, T)

        tdos = Dos(self.efermi / units.eV, energies / units.eV,
                   {Spin(1): densities})

        if partial_dos:
            tdos = self.get_partial_doses(tdos=tdos, npts_mu=npts_mu, T=T)

        return tdos

    def get_partial_doses(self, tdos, npts_mu, T):
        """
            Return a CompleteDos object interpolating the projections
            
            tdos: total dos previously calculated
            npts_mu: number of energy points of the Dos
            T: parameter used to smooth the Dos
        """
        if self.data.proj == []:
            raise BoltztrapError("No projections loaded.")

        bkp_data_ebands = np.copy(self.data.ebands)

        pdoss = {}
        # for spin in self.data.proj:
        for isite, site in enumerate(self.data.structure.sites):
            if site not in pdoss:
                pdoss[site] = {}
            for iorb, orb in enumerate(Orbital):
                if iorb == self.data.proj.shape[-1]: break

                if orb not in pdoss[site]:
                    pdoss[site][orb] = {}

                self.data.ebands = self.data.proj[:, :, isite, iorb].T
                coeffs = fite.fitde3D(self.data, self.equivalences)
                proj, vvproj, cproj = fite.getBTPbands(self.equivalences,
                                                       coeffs, self.data.lattvec)

                edos, pdos = BL.DOS(self.eband, npts=npts_mu, weights=np.abs(proj.real))

                if T is not None:
                    pdos = BL.smoothen_DOS(edos, pdos, T)

                pdoss[site][orb][Spin(1)] = pdos

        self.data.ebands = bkp_data_ebands

        return CompleteDos(self.data.structure, total_dos=tdos, pdoss=pdoss)

    def save(self):
        pass


class BztTransportProperties(object):
    """
        Compute Seebeck, Conductivity, Electrical part of thermal conductivity
        and Hall coefficient, conductivity effective mass, Power Factor tensors
        w.r.t. the chemical potential and temperatures, from dft band structure via
        interpolation.
        
        Args:
            BztInterpolator: a BztInterpolator previously generated
            temp_r: numpy array of temperatures at which to calculate trasport properties
            doping: doping levels at which to calculate trasport properties
            npts_mu: number of energy points at which to calculate trasport properties
            CRTA: constant value of the relaxation time
            
        Upon creation, it contains properties tensors w.r.t. the chemical potential
        of size (len(temp_r),npts_mu,3,3):
            Conductivity_mu (S/m), Seebeck_mu (microV/K), Kappa_mu (W/(m*K)),
            Power_Factor_mu (milliW/K);
            cond_Effective_mass_mu (m_e) calculated as Ref.
        Also:
            Carrier_conc_mu: carrier concentration of size (len(temp_r),npts_mu)
            Hall_carrier_conc_trace_mu: trace of Hall carrier concentration of size
                (len(temp_r),npts_mu)
            mu_r_eV: array of energies in eV and with E_fermi at 0.0
                where all the properties are calculated.
        
        Example:
            bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(100,1400,100))
    """

    def __init__(self, BztInterpolator, temp_r=np.arange(100, 1400, 100),
                 doping=10. ** np.arange(16, 23), npts_mu=4000, CRTA=1e-14, margin=None):

        self.CRTA = CRTA
        self.temp_r = temp_r
        self.doping = doping
        self.dosweight = BztInterpolator.data.dosweight
        lattvec = BztInterpolator.data.get_lattvec()

        self.epsilon, self.dos, self.vvdos, self.cdos = BL.BTPDOS(BztInterpolator.eband, BztInterpolator.vvband,
                                                                  npts=npts_mu, cband=BztInterpolator.cband)

        if margin is None:
            margin = 9. * units.BOLTZMANN * temp_r.max()

        mur_indices = np.logical_and(self.epsilon > self.epsilon.min() + margin,
                                     self.epsilon < self.epsilon.max() - margin)

        self.mu_r = self.epsilon[mur_indices]

        N, L0, L1, L2, Lm11 = BL.fermiintegrals(
            self.epsilon, self.dos, self.vvdos, mur=self.mu_r, Tr=temp_r, dosweight=self.dosweight, cdos=self.cdos)

        self.efermi = BztInterpolator.data.fermi / units.eV
        self.mu_r_eV = self.mu_r / units.eV - self.efermi
        self.nelect = BztInterpolator.data.nelect
        self.volume = BztInterpolator.data.get_volume()

        # Compute the Onsager coefficients from those Fermi integrals
        self.Conductivity_mu, self.Seebeck_mu, self.Kappa_mu, Hall_mu = BL.calc_Onsager_coefficients(L0, L1, L2,
                                                                                                     self.mu_r, temp_r,
                                                                                                     self.volume,
                                                                                                     Lm11=Lm11)

        # Common properties rescaling
        self.Conductivity_mu *= CRTA  # S / m
        self.Seebeck_mu *= 1e6  # microvolt / K
        self.Kappa_mu *= CRTA  # W / (m K)
        self.Hall_carrier_conc_trace_mu = units.Coulomb * 1e-6 / (np.abs(Hall_mu[:, :, 0, 1, 2] +
                                                                         Hall_mu[:, :, 2, 0, 1] +
                                                                         Hall_mu[:, :, 1, 2, 0]) / 3)
        self.Carrier_conc_mu = (N + self.nelect) / (self.volume / (units.Meter / 100.) ** 3)

        # Derived properties
        cond_eff_mass = np.zeros((len(self.temp_r), len(self.mu_r), 3, 3))
        for t in range(len(self.temp_r)):
            for i in range(len(self.mu_r)):
                try:
                    cond_eff_mass[t, i] = np.linalg.inv(self.Conductivity_mu[t, i]) * self.Carrier_conc_mu[
                        t, i] * units.qe_SI ** 2 / units.me_SI * 1e6
                except np.linalg.LinAlgError:
                    pass

        self.Effective_mass_mu = cond_eff_mass * CRTA

        self.Power_Factor_mu = (self.Seebeck_mu @ self.Seebeck_mu) @ self.Conductivity_mu
        self.Power_Factor_mu *= 1e-9  # milliWatt / m / K**2

        # self.props_as_dict()

    def compute_properties_doping(self, doping, temp_r=None):
        """
        Calculate all the properties w.r.t. the doping levels in input.

        Args:
            doping: numpy array specifing the doping levels

        When executed, it add the following variable at the BztTransportProperties
        object:
            Conductivity_doping, Seebeck_doping, Kappa_doping, Power_Factor_doping,
            cond_Effective_mass_doping are dictionaries with 'n' and 'p' keys and
            arrays of dim (len(temp_r),len(doping),3,3) as values
        doping_carriers: number of carriers for each doping level
        mu_doping_eV: the chemical potential corrispondent to each doping level
        """

        if temp_r is None:
            temp_r = self.temp_r

        self.Conductivity_doping, self.Seebeck_doping, self.Kappa_doping = {}, {}, {}
        # self.Hall_doping = {}

        self.Power_Factor_doping, self.Effective_mass_doping = {}, {}

        mu_doping = {}
        doping_carriers = [dop * (self.volume / (units.Meter / 100.) ** 3) for dop in doping]

        for dop_type in ['n', 'p']:
            sbk = np.zeros((len(temp_r), len(doping), 3, 3))
            cond = np.zeros((len(temp_r), len(doping), 3, 3))
            kappa = np.zeros((len(temp_r), len(doping), 3, 3))
            hall = np.zeros((len(temp_r), len(doping), 3, 3, 3))
            if dop_type == 'p':
                doping_carriers = [-dop for dop in doping_carriers]

            mu_doping[dop_type] = np.zeros((len(temp_r), len(doping)))
            for t, temp in enumerate(temp_r):
                for i, dop_car in enumerate(doping_carriers):
                    mu_doping[dop_type][t, i] = self.find_mu_doping(self.epsilon, self.dos, self.nelect + dop_car, temp,
                                                                    self.dosweight)

                N, L0, L1, L2, Lm11 = BL.fermiintegrals(self.epsilon, self.dos, self.vvdos, mur=mu_doping[dop_type][t],
                                                        Tr=np.array([temp]), dosweight=self.dosweight)
                cond[t], sbk[t], kappa[t], hall[t] = BL.calc_Onsager_coefficients(L0, L1, L2, mu_doping[dop_type][t],
                                                                                  np.array([temp]), self.volume, Lm11)

            self.Conductivity_doping[dop_type] = cond * self.CRTA  # S / m
            self.Seebeck_doping[dop_type] = sbk * 1e6  # microVolt / K
            self.Kappa_doping[dop_type] = kappa * self.CRTA  # W / (m K)
            # self.Hall_doping[dop_type] = hall

            self.Power_Factor_doping[dop_type] = (sbk @ sbk) @ cond * self.CRTA * 1e3

            cond_eff_mass = np.zeros((len(temp_r), len(doping), 3, 3))
            for t in range(len(temp_r)):
                for i, dop in enumerate(doping):
                    try:
                        cond_eff_mass[t, i] = np.linalg.inv(cond[t, i]) * dop * units.qe_SI ** 2 / units.me_SI * 1e6
                    except np.linalg.LinAlgError:
                        pass

            self.Effective_mass_doping[dop_type] = cond_eff_mass

        self.doping_carriers = doping_carriers
        self.doping = doping
        self.mu_doping = mu_doping
        self.mu_doping_eV = {k: v / units.eV - self.efermi for k, v in mu_doping.items()}

    def find_mu_doping(self, epsilon, dos, N0, T, dosweight=2.):
        delta = np.empty_like(epsilon)
        for i, e in enumerate(epsilon):
            delta[i] = BL.calc_N(epsilon, dos, e, T, dosweight) + N0
        delta = np.abs(delta)
        # Find the position optimizing this distance
        pos = np.abs(delta).argmin()
        return epsilon[pos]

    def props_as_dict(self):
        props = ("Conductivity", "Seebeck", "Kappa")  # ,"Hall"
        props_unit = (r"$\mathrm{kS\,m^{-1}}$", r"$\mu$V/K", r"")

        p_dict = {'Temps': self.temp_r, 'mu': self.mu_r / units.eV - self.efermi}
        for prop, unit in zip(props, props_unit):
            p_array = eval("self." + prop)
            if prop is not None:
                p_dict[prop] = {'units': unit}
            else:
                continue
            for it, temp in enumerate(self.temp_r):
                p_dict[prop][str(temp)] = {}
                p_dict[prop][str(temp)]['tensor'] = p_array[it]
                p_dict[prop][str(temp)]['eigs'] = np.linalg.eigh(p_array[it])[0]
                p_dict[prop][str(temp)]['avg_eigs'] = p_dict[prop][str(temp)]['eigs'].mean(axis=1)

        self.props_dict = p_dict

    def save(self, fname="Transport_Properties.json"):
        dumpfn(self.props_dict, fname)


class BztPlotter(object):
    """
        Plotter to plot transport properties, interpolated bands along some high
        symmetry k-path, and fermisurface
        
        Example:
            bztPlotter = BztPlotter(bztTransp,bztInterp)
    """

    def __init__(self, bzt_transP=None, bzt_interp=None):
        self.bzt_transP = bzt_transP
        self.bzt_interp = bzt_interp

    def plot_props(self, prop_y, prop_x, prop_z,
                   output='avg_eigs', dop_type='n', doping=None,
                   temps=None, xlim=(-2, 2), ax=None):

        """
            Function to plot the transport properties.
            
            Args:
                prop_y: property to plot among ("Conductivity","Seebeck","Kappa","Carrier_conc","Hall_carrier_conc_trace"). Abbreviations are possible, like "S" for "Seebeck"
                prop_x: independent variable in the x-axis among ('mu','doping','temp')
                prop_z: third variable to plot multiple curves ('doping','temp')
                output: 'avg_eigs' to plot the average of the eigenvalues of the properties
                    tensors; 'eigs' to plot the three eigenvalues of the properties
                    tensors.
                dop_type: 'n' or 'p' to specify the doping type in plots that use doping
                    levels as prop_x or prop_z
                doping: list of doping level to plot, useful to reduce the number of curves
                    when prop_z='doping'
                temps: list of temperatures to plot, useful to reduce the number of curves
                    when prop_z='temp'
                xlim: chemical potential range, useful when prop_x='mu'
                ax: figure.axes where to plot. If None, a new figure is produced.
                
            Example:
            bztPlotter.plot_props('S','mu','temp',temps=[600,900,1200]).show()
            more example are provided in the notebook
            "How to use Boltztra2 interface.ipynb".
        """

        props = ("Conductivity", "Seebeck", "Kappa", "Effective_mass",
                 "Power_Factor", "Carrier_conc", "Hall_carrier_conc_trace")
        props_lbl = ("Conductivity", "Seebeck", "$K_{el}$", "Effective mass",
                     "Power Factor", "Carrier concentration", "Hall carrier conc.")
        props_unit = (r"$(\mathrm{kS\,m^{-1}})$", r"($\mu$V/K)", r"$(W / (m \cdot K))$",
                      r"$(m_e)$", r"$( mW / (m\cdot K^2)$", r"$(cm^{-3})$", r"$(cm^{-3})$")

        props_short = [p[:len(prop_y)] for p in props]

        if prop_y not in props_short:
            raise BoltztrapError("prop_y not valid")

        if prop_x not in ('mu', 'doping', 'temp'):
            raise BoltztrapError("prop_x not valid")

        if prop_z not in ('doping', 'temp'):
            raise BoltztrapError("prop_z not valid")

        idx_prop = props_short.index(prop_y)

        leg_title = ""

        mu = self.bzt_transP.mu_r_eV

        if prop_z == 'doping' and prop_x == 'temp':
            p_array = eval("self.bzt_transP." + props[idx_prop] + '_' + prop_z)
        else:
            p_array = eval("self.bzt_transP." + props[idx_prop] + '_' + prop_x)

        if ax is None:
            fig = plt.figure(figsize=(10, 8))

        temps_all = self.bzt_transP.temp_r.tolist()
        if temps is None:
            temps = self.bzt_transP.temp_r.tolist()

        doping_all = self.bzt_transP.doping.tolist()
        if doping is None:
            doping = self.bzt_transP.doping.tolist()

        # special case of carrier and hall carrier concentration 2d arrays (temp,mu)
        if idx_prop in [5, 6]:
            if prop_z == 'temp' and prop_x == 'mu':
                for temp in temps:
                    ti = temps_all.index(temp)
                    prop_out = p_array[ti]
                    plt.semilogy(mu, prop_out, label=str(temp) + ' K')

                plt.xlabel(r"$\mu$ (eV)", fontsize=30)
                plt.xlim(xlim)
            else:
                raise BoltztrapError("only prop_x=mu and prop_z=temp are available for c.c. and Hall c.c.!")

        elif prop_z == 'temp' and prop_x == 'mu':
            for temp in temps:
                ti = temps_all.index(temp)
                prop_out = np.linalg.eigh(p_array[ti])[0]
                if output == 'avg_eigs':
                    plt.plot(mu, prop_out.mean(axis=1), label=str(temp) + ' K')
                elif output == 'eigs':
                    for i in range(3):
                        plt.plot(mu, prop_out[:, i],
                                 label='eig ' + str(i) + ' ' + str(temp) + ' K')

            plt.xlabel(r"$\mu$ (eV)", fontsize=30)
            plt.xlim(xlim)

        elif prop_z == 'temp' and prop_x == 'doping':
            for temp in temps:
                ti = temps_all.index(temp)
                prop_out = np.linalg.eigh(p_array[dop_type][ti])[0]
                if output == 'avg_eigs':
                    plt.semilogx(doping_all, prop_out.mean(axis=1), 's-',
                                 label=str(temp) + ' K')
                elif output == 'eigs':
                    for i in range(3):
                        plt.plot(doping_all, prop_out[:, i], 's-',
                                 label='eig ' + str(i) + ' ' + str(temp) + ' K')
            plt.xlabel(r"Carrier conc. $cm^{-3}$", fontsize=30)
            leg_title = dop_type + "-type"

        elif prop_z == 'doping' and prop_x == 'temp':
            for dop in doping:
                di = doping_all.index(dop)
                prop_out = np.linalg.eigh(p_array[dop_type][:, di])[0]
                if output == 'avg_eigs':
                    plt.plot(temps_all, prop_out.mean(axis=1),
                             's-', label=str(dop) + ' $cm^{-3}$')
                elif output == 'eigs':
                    for i in range(3):
                        plt.plot(temps_all, prop_out[:, i], 's-',
                                 label='eig ' + str(i) + ' ' + str(dop) + ' $cm^{-3}$')

            plt.xlabel(r"Temperature (K)", fontsize=30)
            leg_title = dop_type + "-type"

        plt.ylabel(props_lbl[idx_prop] + ' ' + props_unit[idx_prop], fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend(title=leg_title if leg_title != "" else "", fontsize=15)
        plt.tight_layout()
        plt.grid()
        return plt

    def plot_bands(self):
        """
            Plot a band structure on symmetry line using BSPlotter()
        """
        if self.bzt_interp is None:
            raise BoltztrapError("BztInterpolator not present")

        sbs = self.bzt_interp.get_band_structure()

        return BSPlotter(sbs).get_plot()

    def plot_dos(self, T=None, npoints=10000):
        """
            Plot the total Dos using DosPlotter()
        """
        if self.bzt_interp is None:
            raise BoltztrapError("BztInterpolator not present")
        
        tdos = self.bzt_interp.get_dos(T=T,npts_mu=npoints)
        #print(npoints)
        dosPlotter = DosPlotter()
        dosPlotter.add_dos('Total',tdos)
        
        return dosPlotter
