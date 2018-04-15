#!/usr/bin/env python

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Shyam Dwaraknath"
__email__ = "dbroberg@berkeley.edu, shyamd@lbl.gov"
__status__ = "Development"
__date__ = "January 11, 2018"

import logging
import sys
import math
import numpy as np
norm = np.linalg.norm

from scipy import stats  #for statistical uncertainties of pot alignment
from monty.json import MSONable
from pymatgen.util.coord import pbc_shortest_vectors
# from pymatgen.analysis.defects.core import DefectCorrection
from core import DefectCorrection
# from pymatgen.entries import CompatibilityError

from utils import ang_to_bohr, hart_to_ev, eV_to_k, generate_reciprocal_vectors_squared, QModel, genrecip

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class FreysoldtCorrection(DefectCorrection):
    """
    A class for FreysoldtCorrection class. Largely adapated from PyCDT code
    Requires some parameters in the DefectEntry to properly function:
        axis_grid
            list of numpy arrays of x-axis values (in angstroms) corresponding to each avg esp supplied

        bulk_planar_averages

        defect_planar_averages


    axis_grid, pureavg, defavg, lattice, dielectricconst, q, defect_position, axis,
                          q_model=QModel(), madetol=0.0001, title=None, widthsample=1.0

    """

    def __init__(self, dielectric_const, q_model=None, energy_cutoff=520, madelung_energy_tolerance=0.0001):
        self.dielectric_const = dielectric_const
        self.q_model = QModel() if not q_model else q_model
        self.energy_cutoff = energy_cutoff
        self.madelung_energy_tolerance = madelung_energy_tolerance
        self.metadata = {}

        if isinstance(dielectric_const, int) or \
                isinstance(dielectric_const, float):
            self.dielectric = float(dielectric_const)
        else:
            self.dielectric = float(np.mean(np.diag(dielectric_const)))

    def get_correction(self, entry):
        """
        Gets the Freysoldt correction for a defect entry
        """

        list_axis_grid = np.array(entry.parameters["axis_grid"])
        list_bulk_plnr_avg_esp = np.array(entry.parameters["bulk_planar_averages"])
        list_defect_plnr_avg_esp = np.array(entry.parameters["defect_planar_averages"])
        list_axes = range(len(list_axis_grid))

        lattice = entry.defect.bulk_structure.lattice
        q = entry.defect.charge

        es_corr = self.perform_es_corr(
            lattice,
            self.dielectric,
            entry.charge,
            energy_cutoff=self.energy_cutoff,
            q_model=self.q_model,
            madetol=self.madelung_energy_tolerance)

        pot_corr_tracker = []
        for x, pureavg, defavg, axis in zip(list_axis_grid,
                                            list_bulk_plnr_avg_esp,
                                            list_defect_plnr_avg_esp,
                                            list_axes):
            tmp_pot_corr = self.perform_pot_corr(
                x,
                pureavg,
                defavg,
                lattice,
                self.dielectric,
                entry.charge,
                entry.site.coords,
                axis,
                q_model=self.q_model,
                madetol=self.madelung_energy_tolerance,
                title=None,
                widthsample=1.0)
            pot_corr_tracker.append(tmp_pot_corr)

        pot_corr = np.mean(pot_corr_tracker)

        return {"freysoldt_electrostatic": es_corr,
                "freysoldt_potential_alignment": pot_corr}


    def perform_es_corr(self, lattice, dielectric, q, energy_cutoff=520, q_model=QModel(), madetol=0.0001):
        """
        Peform Electrostatic Freysoldt Correction
        """
        logger = logging.getLogger(__name__)
        logger.info('Running Freysoldt 2011 PC calculation (should be '\
                     'equivalent to sxdefectalign)')
        logger.debug('defect lattice constants are (in angstroms)' \
                      + str(lattice.abc))

        [a1, a2, a3] = ang_to_bohr * np.array(lattice.get_cartesian_coords(1))
        logging.debug( 'In atomic units, lat consts are (in bohr):' \
                      + str([a1, a2, a3]))
        vol = np.dot(a1, np.cross(a2, a3))  #vol in bohr^3

        #compute isolated energy
        step = 1e-4
        encut1 = 20  #converge to some smaller encut first [eV]
        flag = 0
        converge = []
        while (flag != 1):
            eiso = 1.
            gcut = eV_to_k(encut1)  #gcut is in units of 1/A
            g = step  #initalize
            while g < (gcut + step):
                #simpson integration
                eiso += 4 * (q_model.rho_rec(g * g)**2)
                eiso += 2 * (q_model.rho_rec((g + step)**2)**2)
                g += 2 * step
            eiso -= q_model.rho_rec(gcut**2)**2
            eiso *= (q**2) * step / (3 * round(np.pi, 6))
            converge.append(eiso)
            if len(converge) > 2:
                if abs(converge[-1] - converge[-2]) < madetol:
                    flag = 1
                elif encut1 > energy_cutoff:
                    logger.error('Eiso did not converge before ' \
                                  + str(energy_cutoff) + ' eV')
                    raise
            encut1 += 20
        eiso = converge[-1]
        logger.debug('Eisolated : %f, converged at encut: %d', round(eiso, 5), encut1 - 20)

        #compute periodic energy;
        encut1 = 20  #converge to some smaller encut
        flag = 0
        converge = []
        while flag != 1:
            eper = 0.0
            for g2 in generate_reciprocal_vectors_squared(a1, a2, a3, encut1):
                eper += (q_model.rho_rec(g2)**2) / g2
            eper *= (q**2) * 2 * round(np.pi, 6) / vol
            eper += (q**2) *4* round(np.pi, 6) \
                    * q_model.rho_rec_limit0 / vol
            converge.append(eper)
            if len(converge) > 2:
                if abs(converge[-1] - converge[-2]) < madetol:
                    flag = 1
                elif encut1 > energy_cutoff:
                    logger.error('Eper did not converge before %d eV', energy_cutoff)
                    return
            encut1 += 20
        eper = converge[-1]

        logger.info('Eperiodic : %f hartree, converged at encut %d eV', round(eper, 5), encut1 - 20)
        logger.info('difference (periodic-iso) is %f hartree', round(eper - eiso, 6))
        logger.info('difference in (eV) is %f', round((eper - eiso) * hart_to_ev, 4))

        es_corr = round((eiso - eper) / dielectric * hart_to_ev, 6)
        logger.info('Defect Correction without alignment %f (eV): ', es_corr)
        return es_corr

    def perform_pot_corr(self,
                         axis_grid,
                         pureavg,
                         defavg,
                         lattice,
                         dielectricconst,
                         q,
                         defect_position,
                         axis,
                         q_model=QModel(),
                         madetol=0.0001,
                         title=None,
                         widthsample=1.0):
        """
        For performing planar averaging potential alignment

        title is for name of plot, if you dont want a plot then leave it as None
        widthsample is the width (in Angstroms) of the region in between defects where the potential alignment correction is averaged
        """
        logger = logging.getLogger(__name__)

        logging.debug('run Freysoldt potential alignment method for axis ' + str(axis))
        nx = len(axis_grid)

        #shift these planar averages to have defect at origin
        axfracval = lattice.get_fractional_coords(defect_position)[axis]
        axbulkval = axfracval * lattice.abc[axis]
        if axbulkval < 0:
            axbulkval += lattice.abc[axis]
        elif axbulkval > lattice.abc[axis]:
            axbulkval -= lattice.abc[axis]

        if axbulkval:
            for i in range(len(axis_grid)):
                if axbulkval < axis_grid[i]:
                    break
            rollind = len(axis_grid) - i
            pureavg = np.roll(pureavg, rollind)
            defavg = np.roll(defavg, rollind)

        #if not self._silence:
        logger.debug('calculating lr part along planar avg axis')
        reci_latt = lattice.reciprocal_lattice
        dg = reci_latt.abc[axis]
        dg /= ang_to_bohr  #convert to bohr to do calculation in atomic units

        v_G = np.empty(len(axis_grid), np.dtype('c16'))
        epsilon = dielectricconst
        # q needs to be that of the back ground
        v_G[0] = 4 * np.pi * -q / epsilon * q_model.rho_rec_limit0
        for i in range(1, nx):
            if (2 * i < nx):
                g = i * dg
            else:
                g = (i - nx) * dg
            g2 = g * g
            v_G[i] = 4 * np.pi / (epsilon * g2) * -q * q_model.rho_rec(g2)
        if not (nx % 2):
            v_G[nx // 2] = 0
        v_R = np.fft.fft(v_G)
        v_R_imag = np.imag(v_R)
        v_R /= (lattice.volume * ang_to_bohr**3)
        v_R = np.real(v_R) * hart_to_ev

        max_imag_vr = v_R_imag.max()
        if abs(max_imag_vr) > madetol:
            logging.error('imaginary part found to be %s', repr(max_imag_vr))
            sys.exit()

        #get correction
        short = (defavg - pureavg - v_R)
        checkdis = int((widthsample / 2) / (axis_grid[1] - axis_grid[0]))
        mid = int(len(short) / 2)

        tmppot = [short[i] for i in range(mid - checkdis, mid + checkdis + 1)]
        logger.debug('shifted defect position on axis (%s) to origin', repr(axbulkval))
        logger.debug('means sampling region is (%f,%f)', axis_grid[mid - checkdis], axis_grid[mid + checkdis])

        C = -np.mean(tmppot)
        logger.debug('C = %f', C)
        final_shift = [short[j] + C for j in range(len(v_R))]
        v_R = [elmnt - C for elmnt in v_R]

        logger.info('C value is averaged to be %f eV ', C)
        logger.info('Potentital alignment energy correction (-q*delta V):  %f (eV)', -q * C)
        self.pot_corr = -q * C

        #log plotting data:
        if 'pot_plot_data' not in self.metadata:   #x, v_R, dft_diff, final_shift, check,
            self.metadata['pot_plot_data'] = {axis: {'Vr':v_R, 'x': axis_grid, 'dft_diff': defavg - pureavg,
                                                     'final_shift':final_shift, 'check': [mid - checkdis, mid + checkdis + 1]}}
        else:
            self.metadata['pot_plot_data'][axis] = {'Vr':v_R, 'x': axis_grid, 'dft_diff': defavg - pureavg,
                                                    'final_shift':final_shift, 'check': [mid - checkdis, mid + checkdis + 1]}

        #log uncertainty:
        if 'pot_corr_uncertainty_md' not in self.metadata:
            self.metadata['pot_corr_uncertainty_md'] = {axis: {'stats': stats.describe(tmppot), 'potcorr': -q * C}}
        else:
            self.metadata['pot_corr_uncertainty_md'][axis] = {'stats': stats.describe(tmppot), 'potcorr': -q * C}

        return self.pot_corr


def freysoldt_plotter(x, v_R, dft_diff, final_shift, check, title = None, saved=False):
    """
    planar average electrostatic potential plotter for Freysoldt
    """
    plt.figure()
    plt.clf()
    plt.plot(x, v_R, c="green", zorder=1, label="long range from model")
    plt.plot(x, dft_diff, c="red", label="DFT locpot diff")
    plt.plot(x, final_shift, c="blue", label="short range (aligned)")

    tmpx = [x[i] for i in range(check[0], check[1])]
    plt.fill_between(tmpx, -100, 100, facecolor='red', alpha=0.15,
                     label='sampling region')

    plt.xlim(round(x[0]), round(x[-1]))
    ymin = min(min(v_R), min(dft_diff), min(final_shift))
    ymax = max(max(v_R), max(dft_diff), max(final_shift))
    plt.ylim(-0.2+ymin, 0.2+ymax)
    plt.xlabel('distance along axis ($\AA$)', fontsize=15)
    plt.ylabel('Potential (V)', fontsize=15)
    plt.legend(loc=9)
    plt.axhline(y=0, linewidth=0.2, color='black')
    plt.title(str(title) + ' defect potential', fontsize=18)
    plt.xlim(0, max(x))
    if saved:
        plt.savefig(str(title) + 'FreyplnravgPlot.pdf')
    else:
        return plt


"""
Below Here is for kumagai correction
"""

def find_optimal_gamma(structure, epsilon, tolerance = 0.0001, max_encut = 510):
    """
    Find optimal gamma by evaluating the brute force reciprocal
    summation and seeing when the values are on the order of 0.75,
    This calculation is the anisotropic Madelung potential at r = (0,0,0).
    """
    logger = logging.getLogger(__name__)
    angset = structure.lattice.get_cartesian_coords(1)
    [a1, a2, a3] = ang_to_bohr * angset  # convert to bohr
    vol = np.dot(a1, np.cross(a2, a3))

    #start with gamma s.t. gamma*L=5 (literature ref says this may be optimal)
    #optimizing gamma for the reciprocal sum to improve convergence
    maximum_gamma = 250 #choose a large gamma to stop trying this scheme at?
    gamma = 5.0/(vol ** (1/3.0))
    optimal_gamma_found = False

    def brute_force_recip_summation(tencut):
        recippart = 0.0
        cnt = 0
        for rec in genrecip(a1, a2, a3, tencut):
            Gdotdiel = np.dot(rec, np.dot(epsilon, rec))
            summand = math.exp(-Gdotdiel / (4 * (gamma ** 2))) / Gdotdiel
            recippart += summand
            cnt +=1
        recippart *= 4*np.pi/vol
        return recippart

    while not optimal_gamma_found:
        encut = 20 #start with small encut for expediency
        recippart1 = brute_force_recip_summation(encut)

        encut += 20
        recippart2 =  brute_force_recip_summation(encut)

        converge = [recippart1, recippart2]

        while abs(abs(converge[0]) - abs(converge[1])) * hart_to_ev > tolerance:
            encut += 20
            recippartnew = brute_force_recip_summation(encut)

            converge.reverse()
            converge[1] = recippartnew
            if encut > max_encut:
                msg = 'Optimal gamma not found before {} eV cutoff'.format(
                            max_encut)
                logger.error(msg)
                raise ValueError(msg)

        logger.debug('Reciprocal sum converged to %f eV',
                     converge[1] * hart_to_ev)
        logger.debug('Converging encut = %d eV', encut)

        if abs(converge[1]) * hart_to_ev < 0.75:
            logger.warning('Reciprocal summation value is less than 0.75 eV.')
            logger.warning('Might lead to errors')
            multiplier = 1.5 if gamma > 15. else 3.
            gamma *= multiplier
            logger.warning('Changing gamma to %f', gamma)
        else:
            logger.debug('Optimized gamma found to be %f', gamma)
            optimal_gamma_found = True

        if gamma > maximum_gamma:
            msg = 'Could not optimize gamma before gamma = %d'.format( maximum_gamma)
            raise ValueError(msg)

    return gamma


def generate_g_sum(structure, epsilon, dim, gamma):
    """
    Compute the reciprocal summation in the anisotropic Madelung
    potential.

    TODO: Get the input to fft cut by half by using rfft instead of fft
    """
    logger = logging.getLogger(__name__)
    logger.debug('Reciprocal summation in Madelung potential')
    latt = structure.lattice
    vol = latt.volume * ang_to_bohr ** 3 # in Bohr^3

    reci_latt = latt.reciprocal_lattice
    [b1, b2, b3] = reci_latt.get_cartesian_coords(1)
    b1 = np.array(b1) / ang_to_bohr # In 1/Bohr
    b2 = np.array(b2) / ang_to_bohr
    b3 = np.array(b3) / ang_to_bohr

    nx, ny, nz = dim
    logging.debug('nx: %d, ny: %d, nz: %d', nx, ny, nz)
    ind1 = np.arange(nx)
    for i in range(int(nx/2), nx):
        ind1[i] = i - nx
    ind2 = np.arange(ny)
    for i in range(int(ny/2), ny):
        ind2[i] = i - ny
    ind3 = np.arange(nz)
    for i in range(int(nz/2), nz):
        ind3[i] = i - nz

    g_array = np.zeros(dim, np.dtype('c16'))
    gamm2 = 4*(gamma**2)
    for i in ind1:
        for j in ind2:
            for k in ind3:
                g = i*b1 + j*b2 + k*b3
                g_eps_g = np.dot(g, np.dot(epsilon, g))
                if i == j == k == 0:
                    continue
                else:
                    g_array[i,j,k] = math.exp(-g_eps_g/gamm2) / g_eps_g

    r_array = np.fft.fftn(g_array)
    over_vol = 4*np.pi/vol # Multiply with q later
    r_array *= over_vol
    r_arr_real = np.real(r_array)
    r_arr_imag = np.imag(r_array)

    max_imag = r_arr_imag.max()
    logger.debug('Max imaginary part found to be %f', max_imag)

    return r_arr_real


class KumagaiCorrection(DefectCorrection):
    """
    A class for KumagaiCorrection class. Largely adapated from PyCDT code
    Requires some parameters in the DefectEntry to properly function:
        dim: dimensions for the NGX grid for the FFT to be done on

        bulk_atomic_site_averages:  list of bulk structure's atomic site averaged ESPs * charge,
            in same order as indices of bulk structure
            note this is list given by VASP's OUTCAR (so it is multiplied by a test charge of -1)

        defect_atomic_site_averages:  list of defect structure's atomic site averaged ESPs * charge,
            in same order as indices of defect structure
            note this is list given by VASP's OUTCAR (so it is multiplied by a test charge of -1)

        site_matching_indices (list):  list of paids of corresponding index values for bulk and defect site structures
            EXCLUDING the defect site itself  (ex. [[bulk structure site index, defect structure's corresponding site index], ... ]

        optional:
            sampling_radius
            madelung_energy_tolerance (tolerance in eV)
            gamma
            g_sum

    """

    def __init__(self, dielectric_tensor, sampling_radius=None, madelung_energy_tolerance=0.0001, gamma=None, g_sum=None):
        self.madelung_energy_tolerance = madelung_energy_tolerance
        self.metadata = {'gamma': gamma, 'g_sum': g_sum, 'sampling_radius': sampling_radius}

        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self.dielectric = np.identity(3) * dielectric_tensor
        else:
            self.dielectric = np.array(dielectric_tensor)

    def get_correction(self, entry):
        """
        Gets the Kumagai correction for a defect entry
        """
        dim = entry.parameters["dim"]
        bulk_atomic_site_averages = entry.parameters["bulk_atomic_site_averages"]
        defect_atomic_site_averages = entry.parameters["defect_atomic_site_averages"]
        site_matching_indices = entry.parameters["site_matching_indices"]
        bulk_structure = entry.defect.bulk_structure
        q = entry.defect.charge

        if not self.metadata['gamma']:
            if "gamma" in entry.parameters.keys():
                self.metadata['gamma'] = entry.parameters["gamma"]
            else:
                self.metadata['gamma'] = find_optimal_gamma(bulk_structure, self.dielectric,
                                                            self.madelung_energy_tolerance)

        if not len(self.metadata['g_sum']):
            if "g_sum" in entry.parameters.keys():
                self.metadata['g_sum'] = entry.parameters["g_sum"]
            else:
                self.metadata['g_sum'] = generate_g_sum(bulk_structure, self.dielectric,
                                                        dim, self.metadata['gamma'])

        es_corr = self.perform_es_corr( bulk_structure, q, self.metadata['g_sum'],
                                        self.metadata['gamma'], self.madelung_energy_tolerance)

        #create defective structure with entry
        defect_structure = entry.defect.generate_defect_structure()
        defect_position = entry.defect.site

        # if no sampling radius specified, then assuming Wigner-Seitz radius:
        if not self.metadata['sampling_radius']:
            wz = bulk_structure.lattice.get_wigner_seitz_cell()
            dist = []
            for facet in wz:
                midpt = np.mean(np.array(facet), axis=0)
                dist.append(norm(midpt))
            self.metadata['sampling_radius'] = min(dist)

        #assemble site_list as matching [[bulk_site object, defect_site object, Vqb for site], ... repeat for all non defective sites]
        site_list = []
        for bs_ind, ds_ind in site_matching_indices:
            Vqb = - (defect_atomic_site_averages[ds_ind] - bulk_atomic_site_averages[bs_ind])
            site_list.append([ bulk_structure[bs_ind], defect_structure[ds_ind], Vqb])

        pot_corr = self.perform_pot_corr( bulk_structure, defect_structure, defect_position, site_list,
                                          self.metadata['sampling_radius'], q,
                                          self.metadata['g_sum'], dim, self.metadata['gamma'],
                                          self.madelung_energy_tolerance)

        return {"kumagai_electrostatic": es_corr,
                "kumagai_potential_alignment": pot_corr}


    def perform_es_corr(self, structure, q, g_sum, gamma, madetol):
        """
        Peform Electrostatic Kumagai Correction
        Args:
            structure: Bulk pymatgen structure type
            dielectric: dielectric tensor
            q: Point charge (in units of e+)
            g_sum : comes from KumagaiBulkInit class
            gamma : convergence parameter optimized in KumagaiBulkInit class
            madetol: convergence for real sum term in eV
        """
        logger = logging.getLogger(__name__)
        angset = structure.lattice.get_cartesian_coords(1)
        [a1, a2, a3] = ang_to_bohr * angset  # convert to bohr
        vol = np.dot(a1, np.cross(a2, a3))
        determ = np.linalg.det(self.dielectric)
        invdiel = np.linalg.inv(self.dielectric)

        recip_part = q*g_sum[0,0,0]
        realpre = q / np.sqrt(determ)

        #Real space sum by converging with respect to real space vectors
        #create list of real space vectors that satisfy |i*a1+j*a2+k*a3|<=N
        Nmaxlength = 40  #max range to stopping considering for real space sum convergence (arbitrary, but "large")
        N = 2
        r_sums = []
        while N < Nmaxlength:
            real_part = 0.0
            for i in range(-N, N+1):
                for j in range(-N, N+1):
                    for k in range(-N, N+1):
                        if i == j == k == 0: #since this real sum involves r=0, skipping
                            continue
                        else:
                            r_vec = i*a1 + j*a2 + k*a3
                            loc_res = np.dot(r_vec, np.dot(invdiel, r_vec))
                            nmr = math.erfc(gamma * np.sqrt(loc_res))
                            dmr = np.sqrt(determ * loc_res)
                            real_part += nmr / dmr
            r_sums.append([N, realpre * real_part])

            if N == Nmaxlength-1:
                logging.getLogger(__name__).warning(
                    'Direct part could not converge with real space translation '
                    'tolerance of {} for gamma {}'.format(Nmaxlength-1, gamma))
                return
            elif len(r_sums) > 3:
                if abs(abs(r_sums[-1][1]) - abs(r_sums[-2][1])) * hart_to_ev < madetol: #converging in eV
                    real_part = r_sums[-1][1]
                    logging.debug("gamma is {}".format(gamma))
                    logging.getLogger(__name__).debug(
                        "convergence for real summation term occurs at step {} "
                        "where real sum is {}".format(N,  real_part * hart_to_ev))
                    break

            N += 1


        selfint = q*np.pi / (vol * (gamma**2)) #self interaction term
        surfterm = 2*gamma*q / np.sqrt(np.pi*determ) #surface term

        logger.debug('reciprocal part: {}'.format(recip_part * hart_to_ev))
        logger.debug('real part: {}'.format(real_part * hart_to_ev))
        logger.debug('self interaction part: {}'.format(selfint * hart_to_ev))
        logger.debug('surface term: {}'.format(surfterm * hart_to_ev))

        energy_pc = -(q*0.5*hart_to_ev) * (real_part + recip_part - selfint - surfterm)

        logger.debug('Final PC Energy term: %f eV (%f Hartree)',
                     energy_pc, energy_pc/hart_to_ev)

        return energy_pc


    def perform_pot_corr(self, bulk_structure, defect_structure, defect_position, site_list,
                         sampling_radius, q, g_sum, dim, gamma, madetol):
        """
        For performing potential alignment in manner described by Kumagai et al.
        Args:
            bulk_structure: Bulk pymatgen structure type
            defect_structure: Bulk pymatgen structure type
            defect_position: Defect Position, given as a Pymatgen Site object
            site_list: list of NON-defect site information in form:
                [[bulk_site object, defect_site object, Vqb for site], ... repeat for all non defective sites]
                where Vqb is the difference in the site averaged electrostatic potentials between defect and bulk
                    Note: In Vasp this is the negative value of what is given in OUTCAR:  Vqb = -(defect_cell_site - bulk_cell_site))
            sampling_radius: radius (in Angstrom) which sites must be outside of to be included in the correction
                Kumagai paper suggests Wigner-Seitz radius

            dielectric: dielectric tensor
            q: Point charge (in units of e+)
            g_sum : comes from KumagaiBulkInit class
            dim: dimlist for ngx
            gamma : convergence parameter optimized in KumagaiBulkInit class
            madetol: convergence for real sum term in eV


        """
        logger = logging.getLogger(__name__)
        logger.info('\nRunning potential alignment (atomic site averaging)')
        abc = bulk_structure.lattice.abc
        angset = bulk_structure.lattice.get_cartesian_coords(1)
        [a1, a2, a3] = ang_to_bohr * angset  # convert to bohr
        vol = np.dot(a1, np.cross(a2, a3))
        determ = np.linalg.det(self.dielectric)
        invdiel = np.linalg.inv(self.dielectric)
        selfint = q * np.pi / (vol * (gamma ** 2)) #self interaction term to include in Vpc

        #construct relevant information for all defects
        pot_dict = {}   #keys will be site index in the bulk structure
        for_correction = [] #region to sample for correction
        for bulk_cell_site, defect_cell_site, Vqb in site_list:
            bulk_index = bulk_structure.sites.index( bulk_cell_site)
            defect_cell_index = defect_structure.sites.index( defect_cell_site) #shoudl include a raise value error if the site_list contains sites not in the structure being analyzed...

            #calculate distance from defect's position to closest periodic image of current defect_cell site object
            relative_vector = pbc_shortest_vectors(bulk_structure.lattice, defect_position.frac_coords,
                                                   defect_cell_site.frac_coords)[0][0] #this is angstrom vector
            dist_to_defect = norm(relative_vector) #in angstroms

            pot_dict[bulk_index] = {'defect_cell_index': defect_cell_index, 'Vqb': Vqb,
                                    'dist_to_defect': dist_to_defect}

            #now compute Vpc term within Kumagai formalism
            #reciprocal part involves getting corresponding index from g_sum list
            relvec_in_fraccoord = bulk_structure.lattice.get_fractional_coords(relative_vector)
            recip_grid_indices = []
            for i in range(3):
                if relvec_in_fraccoord[i] < 0:
                    while relvec_in_fraccoord[i] < 0:
                        relvec_in_fraccoord[i] += 1
                elif relvec_in_fraccoord[i] >= 1:
                    while relvec_in_fraccoord[i] >= 1:
                        relvec_in_fraccoord[i] -= 1
                relvec_in_fraccoord[i] *= abc[i]
                num_pts = dim[i]
                x = [now_num / float(num_pts) * abc[i] for now_num in range(num_pts)]
                dx = x[1] - x[0]
                x_rprojection_delta_abs = np.absolute(x - relvec_in_fraccoord[i])
                ind = np.argmin(x_rprojection_delta_abs)
                if x_rprojection_delta_abs[ind] > dx*1.1: #to avoid numerical errors
                    logger = logging.getLogger(__name__)
                    logger.error("Input position not within the g_sum grid")
                    logger.error("%d, %d, %f", i, ind, relvec_in_fraccoord)
                    logger.error("%f", x_rprojection_delta_abs)
                    raise ValueError("Input position is not within the g_sum grid")
                recip_grid_indices.append(ind)

            i,j,k = recip_grid_indices
            recip_part = q * g_sum[i, j, k]

            #Real space sum by converging with respect to real space vectors
            #create list of real space vectors that satisfy |i*a1+j*a2+k*a3|<=N
            realpre = q / np.sqrt(determ)
            Nmaxlength = 40  #max range to stopping considering for real space sum convergence (arbitrary, but "large")
            N = 2
            r_sums = []
            while N < Nmaxlength:
                real_part = 0.0
                for i in range(-N, N+1):
                    for j in range(-N, N+1):
                        for k in range(-N, N+1):
                            r_vec = i*a1 + j*a2 + k*a3 - (relative_vector * ang_to_bohr) #in bohr
                            loc_res = np.dot(r_vec, np.dot(invdiel, r_vec))
                            nmr = math.erfc(gamma * np.sqrt(loc_res))
                            dmr = np.sqrt(determ * loc_res)
                            real_part += nmr / dmr
                r_sums.append([N, realpre * real_part])

                if N == Nmaxlength-1:
                    logging.getLogger(__name__).warning(
                        'Direct part could not converge with real space translation '
                        'tolerance of {} for gamma {}'.format(Nmaxlength-1, gamma))
                    return
                elif len(r_sums) > 3:
                    if abs(abs(r_sums[-1][1]) - abs(r_sums[-2][1])) * hart_to_ev < madetol: #converge in eV
                        real_part = r_sums[-1][1]
                        logging.debug("gamma is {}".format(gamma))
                        logging.getLogger(__name__).debug(
                            "convergence for real summatin term occurs at step {} "
                            "where real sum is {}".format(N,  real_part * hart_to_ev))
                        break

                N += 1

            #now add up total madelung potential part with self interaction term
            Vpc = hart_to_ev * (real_part + recip_part - selfint)
            pot_dict[bulk_index]['Vpc'] = Vpc

            logger.debug('Atom: %d, anisotropic madelung potential: %f',
                          bulk_index, Vpc)
            logger.debug('\treciprocal part: {}'.format(recip_part * hart_to_ev))
            logger.debug('\treal part: {}'.format(real_part * hart_to_ev))
            logger.debug('\tself interaction part: {}'.format(selfint * hart_to_ev))
            logger.debug('Atom: %d, bulk/defect potential difference = %f', bulk_index, Vqb)

            if dist_to_defect > sampling_radius:
                logger.debug('Atom: %d, has radius %f which is inside sampling radius %f',
                             dist_to_defect, sampling_radius)
                for_correction.append( Vqb - Vpc)


        pot_alignment = np.mean(for_correction)
        pot_corr = -q * pot_alignment

        #log uncertainty stats:
        self.metadata['pot_corr_uncertainty_md'] = {'stats': stats.describe(for_correction),
                                                    'number_sampled': len(for_correction), 'AllData': pot_dict}

        logger.info('Kumagai potential alignment (site averaging): %f', pot_alignment)
        logger.info('Kumagai potential alignment correction energy: %f eV', pot_corr)

        return pot_corr


def kumagai_plotter(r, Vqb, Vpc, eltnames, samplerad=None, title = None, saved=False):
    """
    atomic site  electrostatic potential plotter for Kumagai
    """
    plt.figure()
    plt.clf()
    collis = ['b', 'g', 'c', 'm', 'y', 'w', 'k']
    full = []
    for_corr = []
    for eltind, elt in enumerate(eltnames): #enumerate atom types for list
        if samplerad:
            for rval, site_vqb, site_vpc in zip( r[eltind], Vqb[eltind], Vpc[eltind]):
                full.append([ rval, site_vqb - site_vpc])
                if rval >= samplerad:
                    for_corr.append( site_vqb - site_vpc)
        plt.plot(r[eltind], Vqb[eltind], color=collis[eltind], marker='^', linestyle='None',
                 label=str(elt) + ': $V_{q/b}$')
        plt.plot(r[eltind], Vpc[eltind], color=collis[eltind], marker='o', linestyle='None',
                 label=str(elt) + ': $V_{pc}$')

    realfull = np.array(sorted(full, key=lambda x: x[0]))

    plt.plot( realfull[:,0], realfull[:,1], color=collis[-1], marker='x', linestyle='None',
             label='$V_{q/b}$ - $V_{pc}$')
    plt.xlabel('Distance from defect ($\AA$)',fontsize=15)
    plt.ylabel('Potential (V)',fontsize=15)

    x = np.arange(samplerad, max(realfull[:,0])*1.1, 0.01)
    plt.fill_between(x, min(realfull[:,1]) - .1, max(realfull[:,1]) + .1, facecolor='red',
                     alpha=0.15, label='sampling region')
    plt.axhline(y=np.mean(for_corr), linewidth=0.5, color='red',
                label='pot. align. / q')

    plt.legend(bbox_to_anchor=(1.05, 0.5))
    plt.axhline(y=0, linewidth=0.2, color='black')

    plt.title(str(title) + ' defect potential', fontsize=18)
    plt.xlim(0, max(x))
    if saved:
        plt.savefig(str(title) + 'KumagaisiteavgPlot.pdf')
    else:
        return plt


"""
Below Here is for other corrections
"""

class BandFillingCorrection(DefectCorrection):
    """
    A class for BandFillingCorrection class. Largely adapted from PyCDT code

    Requires some parameters in the DefectEntry to properly function:
        eigenvalues
            dictionary of defect eigenvalues, as stored in a Vasprun

        kpoint_weights
            kpoint weights corresponding to the dictionary of eigenvalues

        potalign
            potential alignment for the defect calculation
            Only applies to non-zero charge,
            When using potential alignment Correction (freysoldt or kumagai), need to divide by -q

        cbm
            CBM of bulk calculation (or band structure calculation of bulk);
            calculated on same level of theory as the eigenvalues list (ex. GGA defects -> need GGA cbm

        vbm
            VBM of bulk calculation (or band structure calculation of bulk);
            calculated on same level of theory as the eigenvalues list (ex. GGA defects -> need GGA vbm

    """
    def __init__(self):
        self.metadata = {'occupied_def_levels': [], 'unoccupied_def_levels': [],
                         'total_occupation_defect_levels': None,
                         'num_hole_vbm': None, 'num_elec_cbm': None, 'potalign': None }

    def get_correction(self, entry):
        """
        Gets the BandFilling correction for a defect entry
        """
        eigenvalues = entry.parameters["eigenvalues"]
        kpoint_weights = entry.parameters["kpoint_weights"]
        potalign = entry.parameters["potalign"]
        vbm = entry.parameters["vbm"]
        cbm = entry.parameters["cbm"]

        bf_corr = self.perform_bandfill_corr( eigenvalues, kpoint_weights, potalign, vbm, cbm)

        return {'bandfilling': bf_corr}


    def perform_bandfill_corr(self, eigenvalues, kpoint_weights, potalign, vbm, cbm):
        """
        This calculates the band filling correction based on excess of electrons/holes in CB/VB...

        Note that the total free holes and electrons may also be used for a "shallow donor/acceptor"
               correction with specified band shifts: +num_elec_cbm * Delta E_CBM (or -num_hole_vbm * Delta E_VBM)
               [this is done in the LevelShiftingCorrection class]
        """
        bf_corr = 0.

        self.metadata['potalign'] = potalign
        self.metadata['num_hole_vbm'] = 0.
        self.metadata['num_elec_cbm'] = 0.

        if len(eigenvalues.keys()) == 1: #needed because occupation of non-spin calcs is still 1... should be 2
            spinfctr = 1.
        elif len(eigenvalues.keys()) == 2:
            spinfctr = 2.
        else:
            raise ValueError("Eigenvalue keys greater than 2")

        #for tracking mid gap states...
        resolution = 0.01 #this is energy resolution to maintain for gap states
        occupied_midgap = {en:[] for en in np.arange( vbm, cbm+resolution, resolution)}
        occupation = {en:0. for en in np.arange( vbm, cbm+resolution, resolution)}
        unoccupied_midgap = {en:[] for en in np.arange( vbm, cbm+resolution, resolution)}

        for spinset in eigenvalues.values():
            for kptset, weight in zip(spinset, kpoint_weights):
                for eig in kptset: #eig[0] is eigenvalue and eig[1] is occupation
                    neweig = [eig[0] + potalign, eig[1]] #apply potential shift to eigenvalues
                    if (neweig[1] and (neweig[0] > cbm)): #donor MB correction
                        bf_corr += weight * spinfctr * neweig[1] * (neweig[0] - cbm) # "move the electrons down"
                        self.metadata['num_elec_cbm'] += weight  * spinfctr * neweig[1]
                    elif (neweig[1] !=1.) and (neweig[0] <= vbm): #acceptor MB correction
                        bf_corr += weight * spinfctr * (1.-neweig[1]) * (vbm - neweig[0]) #"move the holes up"
                        self.metadata['num_hole_vbm'] += weight * spinfctr * (1.-neweig[1])
                    elif (neweig[0] > vbm) and (neweig[0] < cbm):
                        for en in np.arange(vbm, cbm+resolution, resolution):
                            if (neweig[0]< en + resolution) and (neweig[0]>en):
                                if neweig[1]:
                                    occupied_midgap[en].append(neweig[0])
                                    occupation[en]+= neweig[1] * weight * spinfctr
                                else:
                                    unoccupied_midgap[en].append(neweig[0])
                                continue

        bf_corr *= -1  #need to take negative of this shift for energetic correction

        # summarize defect level results
        self.metadata['total_occupation_defect_levels'] = 0.
        self.metadata['occupied_def_levels'] = []
        self.metadata['unoccupied_def_levels'] = []
        for en in occupied_midgap.keys():
            if len(occupied_midgap[en]):
                self.metadata['occupied_def_levels'].append([np.mean(occupied_midgap[en]), occupation[en]])
                self.metadata['total_occupation_defect_levels'] += occupation[en]
            elif len(unoccupied_midgap[en]):
                self.metadata['unoccupied_def_levels'].append(np.mean(unoccupied_midgap[en]))

        return bf_corr


