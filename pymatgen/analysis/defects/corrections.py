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
import numpy as np
from scipy import stats #for statistical uncertainties of pot alignment

from pymatgen.io.vasp.outputs import Locpot
from utils import ang_to_bohr, hart_to_ev, eV_to_k, generate_reciprocal_vectors_squared
# from pymatgen.analysis.defects.utils import ang_to_bohr, hart_to_ev, eV_to_k, generate_reciprocal_vectors_squared

class OtherCorrection(object):
    #TODO: ask Shyam if this should be a formal pymatgen Correction class? Not using any ComputedEntries...
    """
    An abstract class for other types of Corrections (Shallow level etc...)

    """
    def __init__(self, **kwargs):
        self._correction = kwargs.get('correction', 0.)
        self.metadata = kwargs.get('metadata', dict())

    @property
    def correction(self):
        """
        Total correction to be applied to Formation energy
        """
        return self._correction



class ChargeCorrection(object):
    #TODO: ask Shyam if this should be a formal pymatgen Correction class? Not using any ComputedEntries...
    """
    An abstract class for Charge Corrections (Freysoldt, Kumagai etc. )

    All charge corrections have two parts: Electrostatic correction and Potential alignment correction
    """
    def __init__(self, **kwargs):
                 # es_corr=0., pot_corr=0., metadata={}):
        self.es_corr = kwargs.get('es_corr', 0.)
        self.pot_corr = kwargs.get('pot_corr', 0.)
        self.metadata = kwargs.get('metadata', dict())

    @property
    def correction(self):
        """
        Total correction to be applied to Formation energy
        """
        return self.es_corr + self.pot_corr

class QModel():
    """
    Model for the defect charge distribution.
    A combination of exponential tail and gaussian distribution is used
    (see Freysoldt (2011), DOI: 10.1002/pssb.201046289 )
    q_model(r) = q [x exp(-r/gamma) + (1-x) exp(-r^2/beta^2)]
            without normalization constants
    By default, gaussian distribution with 1 Bohr width is assumed.
    If defect charge is more delocalized, exponential tail is suggested.
    """
    def __init__(self, beta=1.0, expnorm=0.0, gamma=1.0):
        """
        Args:
            beta: Gaussian decay constant. Default value is 1 Bohr.
                  When delocalized (eg. diamond), 2 Bohr is more appropriate.
            expnorm: Weight for the exponential tail in the range of [0-1].
                     Default is 0.0 indicating no tail .
                     For delocalized charges ideal value is around 0.54-0.6.
            gamma: Exponential decay constant
        """
        self.beta2 = beta * beta
        self.x = expnorm
        self.gamma2 = gamma * gamma
        if expnorm and not gamma:
            raise ValueError("Please supply exponential decay constant.")

    def rho_rec(self, g2):
        """
        Reciprocal space model charge value
        for input squared reciprocal vector.
        Args:
            g2: Square of reciprocal vector

        Returns:
            Charge density at the reciprocal vector magnitude
        """
        return (self.x / np.sqrt(1+self.gamma2*g2)
                + (1-self.x) * np.exp(-0.25*self.beta2*g2))

    def rho_rec_limit0(self):
        """
        Reciprocal space model charge value
        close to reciprocal vector 0 .
        rho_rec(g->0) -> 1 + rho_rec_limit0 * g^2
        """
        return -2*self.gamma2*self.x - 0.25*self.beta2*(1-self.x)

class FreysoldtCorrection(ChargeCorrection):
    """
    A class for FreysoldtCorrection class. Largely adapated from PyCDT code
    """
    def __init__(self, **kwargs):
        super(self.__class__, self).__init__(**kwargs)

    @staticmethod
    def retrieve_correction_from_locpot_files(pure_locpot, defect_locpot, lattice, dielectricconst, q,
                                              defect_position, energy_cutoff=520, q_model=QModel(),
                                              madetol=0.0001, title=None, widthsample=1.0):
        """
        Args:
            pure_locpot:
                Bulk Locpot file path or locpot object
            defect_locpot:
                Defect Locpot file path or locpot object
            dielectricconst:
                Macroscopic dielectric tensor. Include ionic also if
                defect is relaxed, otherwise use ion clamped.
                Can be a matrix array or scalar.
            q (int or float):
                Charge associated with the defect (not of the homogeneous
                background). Typically integer
            energy_cutoff:
                Energy for plane wave cutoff (in eV).
                If not given, Materials Project default 520 eV is used.
            q_model (QModel object):
                User defined charge for correction.
                If not given, highly localized charge is assumed.
            madetol (float):
                Tolerance for convergence of energy terms (in eV)
            defect_position: Defect position as a pymatgen Site object in the bulk supercell structure
                    NOTE: this is optional but recommended, if not provided then analysis is done to find
                    the defect position; this analysis has been rigorously tested, but has broken in an example with
                    severe long range relaxation
                    (at which point you probably should not be including the defect in your analysis...)
        """
        logger = logging.getLogger(__name__)
        fc = FreysoldtCorrection()
        logger.info('This is Freysoldt Correction.')
        if not q:
            fc.es_corr = 0.
            fc.pot_corr = 0.
            return fc

        if isinstance(dielectricconst, int) or \
                isinstance(dielectricconst, float):
            dielectric = float(dielectricconst)
        else:
            dielectric = float(np.mean(np.diag(dielectricconst)))

        if not type(pure_locpot) is Locpot:
            logger.debug('Load bulk locpot')
            pure_locpot = Locpot.from_file(pure_locpot)

        if not type(defect_locpot) is Locpot:
            logger.debug('Load defect locpot')
            defect_locpot = Locpot.from_file(defect_locpot)

        list_axis_grid = []
        list_bulk_plnr_avg_esp = []
        list_defect_plnr_avg_esp = []
        list_axes = range(3)
        for axis in range(3):
            list_axis_grid.append( pure_locpot.get_axis_grid(axis))
            list_bulk_plnr_avg_esp.append( pure_locpot.get_average_along_axis(axis))
            list_defect_plnr_avg_esp.append( defect_locpot.get_average_along_axis(axis))

        fc = FreysoldtCorrection.retrieve_correction_from_plnr_avgs(list_axis_grid, list_bulk_plnr_avg_esp,
                                                                    list_defect_plnr_avg_esp, lattice, dielectricconst,
                                                                    q, defect_position, list_axes,
                                                                    energy_cutoff=energy_cutoff, q_model=q_model,
                                                                    madetol=madetol, title=title, widthsample=widthsample)

        return fc


    @staticmethod
    def retrieve_correction_from_plnr_avgs(list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp, lattice,
                                           dielectricconst, q, defect_position, list_axes, energy_cutoff=520,
                                           q_model=QModel(), madetol=0.0001, title=None, widthsample=1.0):
        """
        Args:
            list_axis_grid
                list of numpy arrays of x-axis values (in angstroms) corresponding to each avg esp supplied
            list_bulk_plnr_avg_esp:
                list of numpy arrays of planar averaged electrostatic potentials of bulk cell [units of eV]
            list_defect_plnr_avg_esp:
                list of numpy arrays of planar averaged electrostatic potential of defect cell [units of eV]
            lattice:
                Pymatgen Lattice object of structure's Supercell
            dielectricconst:
                Macroscopic dielectric tensor. Include ionic also if
                defect is relaxed, otherwise use ion clamped.
                Can be a matrix array or scalar.
            q (int or float):
                Charge associated with the defect (not of the homogeneous
                background). Typically integer
            defect_position:
                Defect position as a pymatgen Site object in the bulk supercell structure
            list_axes:
                list of axis numbers Freysoldt averaging is done over (zero-defined).
                required for calculating defect position alignment
            energy_cutoff:
                Energy for plane wave cutoff (in eV).
                If not given, Materials Project default 520 eV is used.
            madetol (float):
                Tolerance for convergence of energy terms (in eV)
            q_model (QModel object):
                User defined charge for correction.
                If not given, highly localized charge is assumed.
            title:
                set if you want to plot the planar averaged potential
            widthsample:
                is the width (in Angstroms) of the region in between defects where the potential alignment correction is averaged

        """
        logger = logging.getLogger(__name__)
        fc = FreysoldtCorrection()
        logger.info('This is Freysoldt Correction.')
        if not q:
            fc.es_corr = 0.
            fc.pot_corr = 0.
            return fc

        if isinstance(dielectricconst, int) or \
                isinstance(dielectricconst, float):
            dielectric = float(dielectricconst)
        else:
            dielectric = float(np.mean(np.diag(dielectricconst)))

        fc.metadata  = {'dielectric': dielectric, 'list_axis_grid':list_axis_grid,
                        'list_bulk_plnr_avg_esp': list_bulk_plnr_avg_esp, 'lattice': lattice,
                        'list_bulk_plnr_avg_esp': list_defect_plnr_avg_esp, 'list_axes': list_axes, 'charge': q,
                        'defect_position': defect_position, 'madetol': madetol, 'encut':energy_cutoff, 'q_model': q_model}

        fc.es_corr = fc._perform_es_corr(lattice, dielectric, q, energy_cutoff=energy_cutoff, q_model=q_model, madetol=madetol)
        logger.debug('PC calc done, correction = %f', round(fc.es_corr, 4))

        logger.debug('Now run potential alignment script') #note this add information about uncertainty of potalign to metadata
        pot_corr_tracker = []
        for x, pureavg, defavg, axis in zip(list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp, list_axes):
            tmp_pot_corr = fc._perform_pot_corr(x, pureavg, defavg, lattice, dielectric, q, defect_position, axis,
                                               q_model=q_model, madetol=madetol, title=title, widthsample=widthsample)
            pot_corr_tracker.append(tmp_pot_corr)

        fc.pot_corr = np.mean(pot_corr_tracker)

        logger.info('\n\nFreysoldt Correction details:')
        logger.info('PCenergy (E_lat) = %f', round(fc.es_corr, 5))
        logger.info('potential alignment (-q*delta V) = %f',
                     round(fc.pot_corr, 5))
        logger.info('TOTAL Freysoldt correction = %f',
                     round(fc.es_corr + fc.pot_corr, 5))

        return fc

    def _perform_es_corr(self, lattice, dielectricconst, q, energy_cutoff=520, q_model=QModel(), madetol=0.0001):
        """
        Peform Electrostatic Freysoldt Correction
        """
        logger = logging.getLogger(__name__)
        # ap = lattice.get_cartesian_coords(1)
        logger.info('Running Freysoldt 2011 PC calculation (should be '\
                     'equivalent to sxdefectalign)')
        logger.debug('defect lattice constants are (in angstroms)' \
                      + str(lattice.abc))
        [a1, a2, a3] = ang_to_bohr * np.array(lattice.abc)
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
                eiso += 4*(q_model.rho_rec(g*g) ** 2)
                eiso += 2*(q_model.rho_rec((g+step) ** 2) ** 2)
                g += 2 * step
            eiso -= q_model.rho_rec(gcut ** 2) ** 2
            eiso *= (q ** 2) * step / (3 * round(np.pi, 6))
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
        logger.debug('Eisolated : %f, converged at encut: %d',
                      round(eiso, 5), encut1-20)

        #compute periodic energy;
        encut1 = 20  #converge to some smaller encut
        flag = 0
        converge = []
        while flag != 1:
            eper = 0.0
            for g2 in generate_reciprocal_vectors_squared(a1, a2, a3, encut1):
                eper += (q_model.rho_rec(g2) ** 2) / g2
            eper *= (q**2) *2* round(np.pi, 6) / vol
            eper += (q**2) *4* round(np.pi, 6) \
                    * q_model.rho_rec_limit0() / vol
            converge.append(eper)
            if len(converge) > 2:
                if abs(converge[-1] - converge[-2]) < madetol:
                    flag = 1
                elif encut1 > energy_cutoff:
                    logger.error('Eper did not converge before %d eV',
                                  energy_cutoff)
                    return
            encut1 += 20
        eper = converge[-1]

        logger.info('Eperiodic : %f hartree, converged at encut %d eV',
                     round(eper, 5), encut1-20)
        logger.info('difference (periodic-iso) is %f hartree',
                     round(eper - eiso, 6))
        logger.info( 'difference in (eV) is %f',
                     round((eper-eiso) * hart_to_ev, 4))

        self.es_corr = round((eiso-eper) / dielectricconst * hart_to_ev, 6)
        logger.info('Defect Correction without alignment %f (eV): ', self.es_corr)
        return

    def _perform_pot_corr(self, axis_grid, pureavg, defavg, lattice, dielectricconst, q, defect_position, axis,
                          q_model=QModel(), madetol=0.0001, title=None, widthsample=1.0):
        """
        For performing planar averaging potential alignment

        title is for name of plot, if you dont want a plot then leave it as None
        widthsample is the width (in Angstroms) of the region in between defects where the potential alignment correction is averaged
        """
        logger = logging.getLogger(__name__)

        logging.debug('run Freysoldt potential alignment method for axis '+str(axis))
        nx = len(axis_grid)

        #shift these planar averages to have defect at origin
        axfracval=lattice.get_fractional_coords(defect_position)[axis]
        axbulkval=axfracval*lattice.abc[axis]
        if axbulkval<0:
            axbulkval += lattice.abc[axis]
        elif axbulkval > lattice.abc[axis]:
            axbulkval -= lattice.abc[axis]

        if axbulkval:
            for i in range(len(axis_grid)):
                if axbulkval < axis_grid[i]:
                    break
            rollind = len(axis_grid) - i
            pureavg = np.roll(pureavg,rollind)
            defavg = np.roll(defavg,rollind)

        #if not self._silence:
        logger.debug('calculating lr part along planar avg axis')
        reci_latt = lattice.reciprocal_lattice
        dg = reci_latt.abc[axis]
        dg /= ang_to_bohr #convert to bohr to do calculation in atomic units

        v_G = np.empty(len(axis_grid), np.dtype('c16'))
        epsilon = dielectricconst
        # q needs to be that of the back ground
        v_G[0] = 4*np.pi * -q / epsilon * q_model.rho_rec_limit0()
        for i in range(1,nx):
            if (2*i < nx):
                g = i * dg
            else:
                g = (i-nx) * dg
            g2 = g * g
            v_G[i] = 4*np.pi/(epsilon*g2) * -q * q_model.rho_rec(g2)
        if not (nx % 2):
            v_G[nx//2] = 0
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
        checkdis = int((widthsample/2) / (axis_grid[1]-axis_grid[0]))
        mid = int(len(short) / 2)

        tmppot = [short[i] for i in range(mid-checkdis, mid+checkdis+1)]
        logger.debug('shifted defect position on axis (%s) to origin',
                      repr(axbulkval))
        logger.debug('means sampling region is (%f,%f)',
                      axis_grid[mid-checkdis], axis_grid[mid+checkdis])

        C = -np.mean(tmppot)
        logger.debug('C = %f', C)
        final_shift = [short[j] + C for j in range(len(v_R))]
        v_R = [elmnt - C for elmnt in v_R]

        logger.info('C value is averaged to be %f eV ', C)
        logger.info('Potentital alignment energy correction (-q*delta V):  %f (eV)', -q*C)
        self.pot_corr = -q*C

        #log uncertainty:
        if 'pot_corr_uncertainty_md' not in self.metadata:
            self.metadata['pot_corr_uncertainty_md'] = {axis: {'stats':stats.describe(tmppot), 'potcorr': -q*C}}
        else:
            self.metadata['pot_corr_uncertainty_md'][axis] = {'stats':stats.describe(tmppot), 'potcorr': -q*C}

        #make plot, if desired
        if title:
            plotter = FreysoldtCorrPlotter(axis_grid, v_R, defavg-pureavg, final_shift,
                      np.array([mid-checkdis, mid+checkdis]))
            if title != 'written':
                plotter.plot(title=title)
            else:
                # TODO: Make this default fname more defect specific so it doesnt
                # over write previous defect data written
                fname = 'FreyAxisData' # Extension is npz
                plotter.to_datafile(fname)

        return

class FreysoldtCorrPlotter(object):
    def __init__(self, x, v_R, dft_diff, final_shift, check):
        self.x = x
        self.v_R = v_R
        self.dft_diff = dft_diff
        self.final_shift = final_shift
        self.check = check

    def plot(self, title='default'):
        """
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        plt.figure()
        plt.clf()
        plt.plot(self.x, self.v_R, c="green", zorder=1,
                 label="long range from model")
        plt.plot(self.x, self.dft_diff, c="red", label="DFT locpot diff")
        plt.plot(self.x, self.final_shift, c="blue",
                 label="short range (aligned)")
        tmpx = [self.x[i] for i in range(self.check[0], self.check[1])]
        plt.fill_between(tmpx, -100, 100, facecolor='red', alpha=0.15,
                         label='sampling region')

        plt.xlim(round(self.x[0]), round(self.x[-1]))
        ymin = min(min(self.v_R), min(self.dft_diff), min(self.final_shift))
        ymax = max(max(self.v_R), max(self.dft_diff), max(self.final_shift))
        plt.ylim(-0.2+ymin, 0.2+ymax)
        plt.xlabel('distance along axis ' + str(1) + ' ($\AA$)', fontsize=20)
        plt.ylabel('Potential (V)', fontsize=20)
        plt.legend(loc=9)
        plt.axhline(y=0, linewidth=0.2, color='black')
        plt.title(str(title) + ' defect potential')
        plt.xlim(0, max(self.x))

        plt.savefig(str(title) + 'FreyplnravgPlot.pdf')

    def to_datafile(self, file_name='FreyAxisData'):

        np.savez(file_name, x=self.x, v_R=self.v_R,
                 dft_diff=self.dft_diff, #defavg-pureavg,
                 final_shift=self.final_shift, #finalshift,
                 check_range=self.check) #np.array([mid-checkdis, mid+checkdis]))

    @classmethod
    def plot_from_datfile(cls, file_name='FreyAxisData.npz', title='default'):
        """
        Takes data file called 'name' and does plotting.
        Good for later plotting of locpot data after running run_correction()
        """
        with open(file_name) as f:
            plotvals = np.load(f)

            x = plotvals['x']
            v_R = plotvals['v_R']
            dft_diff = plotvals['dft_diff']
            final_shift = plotvals['final_shift']
            check = plotvals['check_range']

            plotter = cls(x, v_R, dft_diff, final_shift, check)
            plotter.plot(title)