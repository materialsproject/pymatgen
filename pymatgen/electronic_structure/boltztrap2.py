import os
import os.path
import pickle
import logging
import spglib
from scipy.optimize import minimize
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
#from monty.serialization import loadfn,dumpfn
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.bandstructure import \
    BandStructureSymmLine, Kpoint, Spin
from pymatgen.io.vasp import Vasprun
from pymatgen.core.units import Energy, Length
from pymatgen.io.ase import AseAtomsAdaptor, Atoms
from pymatgen.electronic_structure.dos import Dos, Spin, CompleteDos, Orbital
from pymatgen.electronic_structure.boltztrap import BoltztrapError
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter

try:
    from BoltzTraP2 import dft as BTP
    from BoltzTraP2 import sphere
    from BoltzTraP2 import fite
    from BoltzTraP2 import bandlib as BL
    from BoltzTraP2 import serialization
    from BoltzTraP2 import units
    #from BoltzTraP2 import fermisurface
except ImportError:
    raise BoltztrapError("BoltzTraP2 has to be installed and working")

#print a progress bar of a loop
import sys
def progress(counter,lenght,step=0):
    if step == 0:
        step = np.ceil(lenght/10).astype('int')
    counter +=1
    perc = int(counter * 100. / lenght)
    l = range(step,100+step,step)
    if perc in l:
        s='|' + '#' * int(perc / step) + " " * int((100-perc) / step)+'| '
        s += str(counter) + "/" + str(lenght)
        print('\r'+str(perc)+'%\t'+s,end='')
        sys.stdout.flush()



class PMG_BS_Loader:
    def __init__(self, pmg_bs_obj,structure=None,nelect=None):
        self.kpoints = np.array([kp.frac_coords for kp in pmg_bs_obj.kpoints])
        
        if structure is None:
            try:
                self.structure = pmg_bs_obj.structure
            except:
                BaseException('No structure found in the bs obj.')
        
        self.atoms = AseAtomsAdaptor.get_atoms(self.structure)
        
        if len(pmg_bs_obj.bands) == 1:
            e = list(pmg_bs_obj.bands.values())[0]
            self.ebands = e * units.eV
            self.dosweight = 2.0
        elif len(pmg_bs_obj.bands) == 2:
            raise BaseException("spin bs case not implemented")
        
        self.lattvec = self.structure.lattice.matrix * units.Angstrom
        self.mommat = None
        self.fermi = pmg_bs_obj.efermi * units.eV
        
        self.nelect = nelect
        self.UCvol = self.structure.volume * units.Angstrom**3
        
    def get_lattvec(self):
        try:
            self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec
    
    def bandana(self, emin=-np.inf, emax=np.inf):
        bandmin = np.min(self.ebands, axis=1)
        bandmax = np.max(self.ebands, axis=1)
        II = np.nonzero(bandmin < emax)
        nemax = II[0][-1]
        II = np.nonzero(bandmax > emin)
        nemin = II[0][0]
        #BoltzTraP2.misc.info("BANDANA output")
        #for iband in range(len(self.ebands)):
            #BoltzTraP2.misc.info(iband, bandmin[iband], bandmax[iband], (
                #(bandmin[iband] < emax) & (bandmax[iband] > emin)))
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

class PMG_Vasprun_Loader:
    def __init__(self, vasprun_file):
        vrun = Vasprun(vasprun_file,parse_projected_eigen=True)
        self.kpoints = np.array(vrun.actual_kpoints)
        self.structure = vrun.final_structure
        self.atoms = AseAtomsAdaptor.get_atoms(self.structure)
        self.proj = []
        if len(vrun.eigenvalues) == 1:
            e = list(vrun.eigenvalues.values())[0]
            self.ebands = e[:,:,0].transpose() * units.eV
            self.dosweight = 2.0
            self.proj = list(vrun.projected_eigenvalues.values())[0]
            
        elif len(vrun.eigenvalues) == 2:
            raise BoltztrapError("spin bs case not implemented")
        
        self.lattvec = self.structure.lattice.matrix * units.Angstrom
        
        #TODO: read mommat from vasprun
        self.mommat = None
        self.fermi = vrun.efermi * units.eV
        self.nelect = vrun.parameters['NELECT']
        self.UCvol = self.structure.volume * units.Angstrom**3
        
    def get_lattvec(self):
        try:
            self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec
    
    def bandana(self, emin=-np.inf, emax=np.inf):
        bandmin = np.min(self.ebands, axis=1)
        bandmax = np.max(self.ebands, axis=1)
        II = np.nonzero(bandmin < emax)
        nemax = II[0][-1]
        II = np.nonzero(bandmax > emin)
        nemin = II[0][0]
        #BoltzTraP2.misc.info("BANDANA output")
        #for iband in range(len(self.ebands)):
            #BoltzTraP2.misc.info(iband, bandmin[iband], bandmax[iband], (
                #(bandmin[iband] < emax) & (bandmax[iband] > emin)))
        self.ebands = self.ebands[nemin:nemax]
        if self.proj != []:
            self.proj = self.proj[:,nemin:nemax,:,:]
            
        if self.mommat is not None:
            self.mommat = self.mommat[:, nemin:nemax, :]
        # Removing bands may change the number of valence electrons
        self.nelect -= self.dosweight * nemin
        return nemin, nemax

    def get_volume(self):
        try:
            self.UCvol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol


class MockDFTData:
    """Mock DFTData emulation class for the parabolic band example."""

    def __init__(self,nb = 1,de = [0.5], effm=[1], nelect=0,efermi=0.0):
        """Create a MockDFTData object based on global variables."""
        
        NK = 25
        
        self.mommat = None
        # Create a hypotetical monoatomic simple cubic structure with a lattice
        # parameter of 5 A.
        atoms = Atoms("Si", cell=5 * np.eye(3), pbc=True)
        self.structure = AseAtomsAdaptor.get_structure(atoms)
        lattvec = atoms.get_cell().T * units.Angstrom
        rlattvec = 2. * np.pi * np.linalg.inv(lattvec).T
        # Create a set of irreducible k points based on that structure
        fromspglib = spglib.get_ir_reciprocal_mesh([NK, NK, NK], atoms)
        indices = np.unique(fromspglib[0]).tolist()
        kpoints = fromspglib[1].T / float(NK)
        kpoints = kpoints[:, indices]
        cartesian = rlattvec @ kpoints
        k2 = (cartesian**2).sum(axis=0)
        rmax = rlattvec[0, 0] / 2.
        rmin = .8 * rmax
        bump = self.create_bump(rmin, rmax)
        weight = (bump(cartesian[0, :]) * bump(cartesian[1, :]) *
                bump(cartesian[2, :]))
        
        self.ebands = np.zeros((nb,kpoints.shape[1]))
        for inb in range(nb):
            # Compute the band energies in a purely parabolic model
            eband = np.abs(de[inb]) + k2 / 2. / effm[inb]
            bound = eband[k2 <= rmax * rmax].max()
            eband = np.minimum(eband, bound)
            eband = eband * weight + bound * (1. - weight)
            self.ebands[inb] = np.sign(de[inb])*eband.reshape((1, eband.size))
    
        self.atoms = atoms
        self.lattvec = lattvec
        self.kpoints = kpoints.T
        self.nelect = nelect
        self.proj = []
        self.fermi = efermi
        self.dosweight = 2
        
    # Weight the band by a bump function to make its derivative zero at the
    # boundary of the BZ.
    def create_bump(self,a, b):
        """Return a bump function f(x) equal to zero for |x| > b, equal to one for
        |x| < a, and providing a smooth transition in between.
        """
        if a <= 0. or b <= 0. or a >= b:
            raise ValueError("a and b must be positive numbers, with b > a")

        # Follow the prescription given by Loring W. Tu in
        # "An Introduction to Manifolds", 2nd Edition, Springer
        def f(t):
            """Auxiliary function used as a building block of the bump function."""
            if t <= 0.:
                return 0.
            else:
                return np.exp(-1. / t)

        f = np.vectorize(f)

        def nruter(x):
            """One-dimensional bump function."""
            arg = (x * x - a * a) / (b * b - a * a)
            return 1. - f(arg) / (f(arg) + f(1. - arg))

        return nruter



    def get_lattvec(self):
        """Return the matrix of lattice vectors."""
        return self.lattvec
    
    def bandana(self, emin=-np.inf, emax=np.inf):
        nemin = np.min(self.ebands, axis=1)
        nemax = np.max(self.ebands, axis=1)
        return nemin, nemax

    def get_volume(self):
        try:
            self.UCvol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol



class BZT_Interpolator(object):
    """
        Interpolate the dft band structures
    """
    #def __init__(self,vasp_dir, lpfac=5, bnd_energy_range=None):
        #self.data = BTP.DFTData(vasp_dir)
        #self.nemin, self.nemax = self.data.bandana(emin=self.data.fermi - .2, emax=self.data.fermi + .2)
        #self.equivalences = sphere.get_equivalences(self.data.atoms, len(self.data.kpoints) * lpfac)
        #self.coeffs = fite.fitde3D(self.data, self.equivalences)
        #self.efermi = self.data.fermi / units.eV
        
    def __init__(self, data, lpfac=10, energy_range=1.5,nelect=None):
        #self.data = PMG_Vasprun_Loader(vasprun_file)
        self.data = data
        num_kpts = self.data.kpoints.shape[0]
        self.efermi = self.data.fermi
        self.nemin, self.nemax = self.data.bandana(emin=self.efermi - (energy_range * units.eV), emax=self.efermi + (energy_range * units.eV))
        self.equivalences = sphere.get_equivalences(self.data.atoms, num_kpts * lpfac)
        self.coeffs = fite.fitde3D(self.data, self.equivalences)
        self.eband, self.vvband, self.cband = fite.getBTPbands(self.equivalences, 
                                                        self.coeffs, self.data.lattvec)

        
    def get_band_structure(self):
        
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
        egrid, vgrid, cgrid = fite.getBands(kpoints, self.equivalences, lattvec, self.coeffs)
        
        bands_dict = {Spin.up: (egrid / units.eV)}
        
        sbs = BandStructureSymmLine(kpoints, bands_dict,
                                    self.data.structure.lattice.reciprocal_lattice, 
                                    self.efermi / units.eV,
                                    labels_dict = labels_dict)
        return sbs
    
    
    def get_dos(self, partial_dos = False, npts_mu = 10000, T=None):

        energies, densities, vvdos, cdos = BL.BTPDOS(self.eband, self.vvband, npts=npts_mu)
        if T is not None:
            densities = BL.smoothen_DOS(energies, densities, T)
            
        tdos = Dos(self.efermi / units.eV, energies / units.eV,
                      {Spin(1): densities})
        
        if partial_dos:
            tdos = self.get_partial_doses(tdos=tdos, npts_mu = npts_mu, T=T)
        
        return tdos
    
    
    def get_partial_doses(self, tdos, npts_mu, T):

        bkp_data_ebands = np.copy(self.data.ebands)
        
        pdoss = {}
        #for spin in self.data.proj:
        cnt = 0
        tot_cnt = np.product(self.data.proj.shape[2:])
        for isite,site in enumerate(self.data.structure.sites):
            if site not in pdoss:
                pdoss[site] = {}
            for iorb, orb in enumerate(Orbital):
                if iorb == self.data.proj.shape[-1]: break
            
                if orb not in pdoss[site]: pdoss[site][orb] = {}
                
                cnt+=1
                progress(cnt,tot_cnt)
                
                self.data.ebands = self.data.proj[:,:,isite,iorb].T
                coeffs = fite.fitde3D(self.data, self.equivalences)
                proj, vvproj, cproj = fite.getBTPbands(self.equivalences, 
                                                        coeffs, self.data.lattvec)
                
                edos, pdos = BL.DOS(self.eband, npts=npts_mu, weights=np.abs(proj.real))

                if T is not None:
                    pdos = BL.smoothen_DOS(edos, pdos, T)
                
                pdoss[site][orb][Spin(1)] = pdos
            
        self.data.ebands = bkp_data_ebands
        
        return CompleteDos(self.data.structure, total_dos=tdos, pdoss=pdoss)

    
    def pockets_finder(self,bands,vbm_idx):
        
        
        def get_energy(kpt,bnd,equivalences,lattvec,coeffs):
            return fite.getBands(kpt[np.newaxis,:], equivalences, lattvec, coeffs)[0][bnd]
            
        R = []
        for i in it.combinations_with_replacement([0,0.05,-0.05],3):
            for j in it.permutations(i):
                if j not in R:
                    R.append(j)
        R=R[1:]
        out = []
        print('ok')
        
        for bnd in bands:
            out.append([])
            bnd -= self.nemin
            vbm_idx -= self.nemin
            cbm = True if bnd > vbm_idx else False

            cnt=0
            for i in np.linspace(0,1,7):
                for j in np.linspace(0,1,7):
                    for k in np.linspace(0,1,7):
                        sys.stdout.flush()
                        xkpt = np.array([i,j,k]) + (np.random.rand(3)*0.05-0.025)
                        suc=False
                        nit=0
                        while suc == False:
                            res=minimize(get_energy,xkpt[np.newaxis,:],
                                            args=(bnd,self.equivalences,
                                                  self.data.lattvec, self.coeffs),
                                            method='BFGS',tol=1e-06)
                            suc=res.success
                            xkpt=res.x
                            ene = res.fun
                            nit+=res.nit

                        found = True
                        
                        for r in R:
                            if ene > get_energy(xkpt,bnd,self.equivalences,
                                                  self.data.lattvec, self.coeffs):
                                found = False
                                break
                        
                        cnt+=1
                        if found:
                            out[-1].append( ' '.join(["step "+str(cnt), str(nit), str(np.array([i,j,k])), str(res.x), 
                            str(res.fun) if cbm else str(-res.fun)]))
                            #sys.stdout.flush()
                        break
                    break
                break
        return out
    
    
    def save(self):
        pass
    
class BZT_TransportProperties(object):
    """
        Compute Seebeck, Conductivity, Electrical part of thermal conductivity 
        and Hall tensor from dft band structure via interpolation.
    """
    def __init__(self, BZT_Interpolator, temp_r = np.arange(100,1300,100),
                 doping=10.**np.arange(16,23), npts_mu=4000, CRTA=1e-14, margin=None ):
        
        self.CRTA = CRTA
        self.temp_r = temp_r
        self.doping = doping
        self.dosweight = BZT_Interpolator.data.dosweight
        lattvec = BZT_Interpolator.data.get_lattvec()

        self.epsilon, self.dos, self.vvdos, self.cdos = BL.BTPDOS(BZT_Interpolator.eband, BZT_Interpolator.vvband, npts=npts_mu)
        
        if margin == None:
            margin = 9. * units.BOLTZMANN * temp_r.max()
            
        mur_indices = np.logical_and(self.epsilon > self.epsilon.min() + margin,
                                    self.epsilon < self.epsilon.max() - margin)

        self.mu_r = self.epsilon[mur_indices]
        
        N, L0, L1, L2, Lm11 = BL.fermiintegrals(
            self.epsilon, self.dos, self.vvdos, mur=self.mu_r, Tr=temp_r, dosweight=self.dosweight)

        self.efermi = BZT_Interpolator.data.fermi / units.eV
        self.mu_r_eV = self.mu_r /units.eV - self.efermi
        self.nelect = BZT_Interpolator.data.nelect
        self.volume = BZT_Interpolator.data.get_volume()
        
        # Compute the Onsager coefficients from those Fermi integrals
        self.Conductivity_mu, self.Seebeck_mu, self.Kappa_mu, self.Hall_mu = BL.calc_Onsager_coefficients(L0, L1, L2, self.mu_r, temp_r, self.volume, Lm11)
        
        self.Carrier_conc = (N+self.nelect)/(self.volume / (units.Meter / 100.)**3)
        
        cond_eff_mass = np.zeros((len(self.temp_r),len(self.mu_r),3,3))
        for t in range(len(self.temp_r)):
            for i in range(len(self.mu_r)):
                try:
                    cond_eff_mass[t,i] = np.linalg.inv(self.Conductivity_mu[t,i]) * self.Carrier_conc[t,i] * units.qe_SI**2 / units.me_SI * 1e6
                except np.linalg.LinAlgError:
                    pass
                
        self.cond_Effective_mass_mu = cond_eff_mass
        
        #TODO: try to use ArrayWithUnits ?
        self.Conductivity_mu *= CRTA # S / m
        self.Seebeck_mu *= 1e6  # microvolt / K
        
        self.Power_Factor_mu = (self.Seebeck_mu @ self.Seebeck_mu) @ self.Conductivity_mu
        self.Power_Factor_mu *= 1e6  # microwatt / m / K**2
        
        #self.props_as_dict()
    
    def compute_properties_doping(self, doping, temp_r=None):
        '''
            Calculate all the properties w.r.t. the doping levels in input.
        '''
        
        if temp_r == None:
            temp_r = self.temp_r
        
        self.Conductivity_doping, self.Seebeck_doping, self.Kappa_doping, self.Hall_doping = {},{},{},{}
        self.Power_Factor_doping, self.cond_Effective_mass_doping = {},{}
        
        mu_doping = {}
        doping_carriers = [dop*(self.volume / (units.Meter / 100.)**3) for dop in doping]
        
        for dop_type in ['n','p']:
            sbk = np.zeros((len(temp_r),len(doping),3,3))
            cond = np.zeros((len(temp_r),len(doping),3,3))
            kappa = np.zeros((len(temp_r),len(doping),3,3))
            hall = np.zeros((len(temp_r),len(doping),3,3,3))
            if dop_type == 'p':
                doping_carriers = [-dop for dop in doping_carriers]
            
            
            mu_doping[dop_type] = np.zeros((len(temp_r),len(doping)))
            for t,temp in enumerate(temp_r):
                for i, dop_car in enumerate(doping_carriers):
                    mu_doping[dop_type][t,i] = self.find_mu_doping(self.epsilon, self.dos, self.nelect + dop_car, temp, self.dosweight)

                N, L0, L1, L2, Lm11 = BL.fermiintegrals(self.epsilon, self.dos, self.vvdos, mur=mu_doping[dop_type][t], Tr=np.array([temp]), dosweight=self.dosweight)
                cond[t], sbk[t], kappa[t], hall[t] = BL.calc_Onsager_coefficients(L0, L1, L2, mu_doping[dop_type][t], np.array([temp]), self.volume, Lm11)
            
            self.Conductivity_doping[dop_type]  = cond * self.CRTA # S / m
            self.Seebeck_doping[dop_type] = sbk * 1e6  # microvolt / K
            self.Kappa_doping[dop_type] = kappa
            self.Hall_doping[dop_type] = hall
            
            self.Power_Factor_doping[dop_type] = (sbk @ sbk) @ cond
            
            cond_eff_mass = np.zeros((len(temp_r),len(doping),3,3))
            for t in range(len(temp_r)):
                for i,dop in enumerate(doping):
                    try:
                        cond_eff_mass[t,i] = np.linalg.inv(cond[t,i]) * dop * units.qe_SI**2 / units.me_SI * 1e6
                    except np.linalg.LinAlgError:
                        pass
                    
            self.cond_Effective_mass_doping[dop_type] = cond_eff_mass
            
        self.doping_carriers = doping_carriers
        self.doping = doping
        self.mu_doping = mu_doping
        self.mu_doping_eV = {k:v/units.eV - self.efermi for k,v in mu_doping.items()}


    def find_mu_doping(self,epsilon, dos, N0, T, dosweight=2.):
        delta = np.empty_like(epsilon)
        for i, e in enumerate(epsilon):
            delta[i] = BL.calc_N(epsilon, dos, e, T, dosweight) + N0
        delta = np.abs(delta)
        # Find the position optimizing this distance
        pos = np.abs(delta).argmin()
        return epsilon[pos]
    
    
    def props_as_dict(self):
        props = ("Conductivity","Seebeck","Kappa")#,"Hall"
        p_units = (r"$\mathrm{kS\,m^{-1}}$",r"$\mu$V/K",r"")
        
        
        p_dict = {'Temps':self.temp_r,'mu':self.mu_r/units.eV - self.efermi}
        for prop,unit in zip(props,p_units):
            p_array = eval("self." + prop)
            if prop != None:
                    p_dict[prop] = {'units':unit} 
            else:
                continue
            for it, temp in enumerate(self.temp_r):
                p_dict[prop][str(temp)] = {}
                p_dict[prop][str(temp)]['tensor'] = p_array[it]
                p_dict[prop][str(temp)]['eigs'] = np.linalg.eigh(p_array[it])[0]
                p_dict[prop][str(temp)]['avg_eigs'] = p_dict[prop][str(temp)]['eigs'].mean(axis=1)
        
        self.props_dict = p_dict
        
    def save(self,fname="Transport_Properties.json"):
        dumpfn(self.props_dict,fname)
    

class BZT_Plotter(object):
    """
        Plotter to plot transport properties, interpolated bands along some KPath,
        and fermisurface
    """
    def __init__(self, bzt_transP=None, bzt_interp=None):
        self.bzt_transP = bzt_transP
        self.bzt_interp = bzt_interp
        
    def plot_props(self, prop_y, prop_x, prop_z, 
                   output='avg_eigs', dop_type='n', doping=None, temps=None,xlim=(-2,2)):
        
        
        props = ("Conductivity","Seebeck","Kappa")#,"Hall"
        
        props_short = [p[:len(prop_y)] for p in props]
        
        if prop_y not in props_short:
            raise BoltztrapError("Property not valid")
        
        idx_prop = props_short.index(prop_y)
        
        p_units = (r"$\mathrm{kS\,m^{-1}}$",r"$\mu$V/K",r"")
        
        if prop_z == 'doping' and prop_x == 'temp':
            p_array = eval("self.bzt_transP." + props[idx_prop]+'_'+prop_z)
        else:
            p_array = eval("self.bzt_transP." + props[idx_prop]+'_'+prop_x)
        mu = self.bzt_transP.mu_r_eV
        
        fig = plt.figure(figsize=(10,8))
        temps_all = self.bzt_transP.temp_r.tolist()
        if temps == None:
            temps = self.bzt_transP.temp_r.tolist()
        
        doping_all = self.bzt_transP.doping.tolist()
        if doping == None:
            doping = self.bzt_transP.doping.tolist()
            
        if prop_z == 'temp' and prop_x == 'mu':
            for temp in temps:
                ti = temps_all.index(temp)
                prop_out = np.linalg.eigh(p_array[ti])[0]
                if output=='avg_eigs':
                    plt.plot(mu, prop_out.mean(axis=1),label=str(temp)+' K')
                elif output=='eigs':
                    for i in range(3):
                        plt.plot(mu, prop_out[:,i],
                                 label='eig_'+str(i)+' ' + str(temp)+' K')
                    
            plt.xlabel(r"$\mu$ (eV)", fontsize=30)
            plt.xlim(xlim)

        elif prop_z == 'temp' and prop_x == 'doping':
            for temp in temps:
                ti = temps_all.index(temp)
                prop_out = np.linalg.eigh(p_array[dop_type][ti])[0]
                if output=='avg_eigs':
                    plt.semilogx(doping_all, prop_out.mean(axis=1),'s-',
                                 label=str(temp)+' K')
                elif output=='eigs':
                    for i in range(3):
                        plt.plot(doping_all, prop_out[:,i],'s-',
                                 label='eig_'+str(i)+' '+str(temp)+' K')                    
            plt.xlabel(r"Carrier Conc $cm^{-3}$", fontsize=30)
            #plt.xlim()
        
        elif prop_z == 'doping' and prop_x == 'temp':
             for dop in doping:
                di = doping_all.index(dop)
                prop_out = np.linalg.eigh(p_array[dop_type][:,di])[0]
                if output=='avg_eigs':
                    plt.plot(temps_all, prop_out.mean(axis=1),
                             's-',label=str(dop)+' $cm^{-3}$')
                elif output=='eigs':
                    for i in range(3):
                        plt.plot(temps_all, prop_out[:,i],'s-',
                                 label='eig_'+str(i)+' '+str(dop)+' $cm^{-3}$')
                
             plt.xlabel(r"Temperature (K)", fontsize=30)
            
        plt.ylabel(props[idx_prop] + ' ' +p_units[idx_prop], fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend(fontsize=15)
        plt.tight_layout()
        plt.grid()
        return plt
    
    
    def plot_bands(self):
        if self.bzt_interp == None:
            raise BoltztrapError("BZT_Interpolator not present")
        
        sbs = self.bzt_interp.get_band_structure()
        
        return BSPlotter(sbs).get_plot()

    def plot_dos(self,T=None):
        if self.bzt_interp == None:
            raise BoltztrapError("BZT_Interpolator not present")
        
        tdos = self.bzt_interp.get_dos(T=T)
        
        dosPlotter = DosPlotter()
        dosPlotter.add_dos('Total',tdos)
        
        return dosPlotter
    
    
    def plot_fermisurface():
        return fermisurface.plot_fermisurface(bzt_interp.data,
                                              bzt_interp.equivalences,
                                              bzt_interp.eband,
                                              mu=bzt_interp.efermi*units/eV)
    
