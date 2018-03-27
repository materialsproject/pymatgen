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

from scipy import stats  #for statistical uncertainties of pot alignment
from monty.json import MSONable
from pymatgen.analysis.defects.core import DefectCorrection
from pymatgen.entries import CompatibilityError

from pymatgen.io.vasp.outputs import Locpot
from utils import ang_to_bohr, hart_to_ev, eV_to_k, generate_reciprocal_vectors_squared


class FreysoldtCorrection(DefectCorrection):
    """
    A class for FreysoldtCorrection class. Largely adapated from PyCDT code
    Requires some parameters in the DefectEntry to properly function:
        axis_grid
            list of numpy arrays of x-axis values (in angstroms) corresponding to each avg esp supplied

        bulk_planar_averages

        defect_planar_averages

        lattice:
            Pymatgen Lattice object of structure's Supercell


    axis_grid, pureavg, defavg, lattice, dielectricconst, q, defect_position, axis,
                          q_model=QModel(), madetol=0.0001, title=None, widthsample=1.0

    """

    def __init__(self, dielectric_const, q_model, energy_cutoff=520, madelung_energy_tolerance=0.0001):
        self.dielectric_const = dielectric_const
        self.q_model = q_model
        self.energy_cutoff = energy_cutoff
        self.madelung_energy_tolerance = madelung_energy_tolerance


        if isinstance(dielectric_const, int) or \
                isinstance(dielectric_const, float):
            self.dielectric = float(dielectric_const)
        else:
            self.dielectric = float(np.mean(np.diag(dielectric_const)))

    def get_correction(self, entry):
        """
        Gets the Freysoldt correction for a defect entry
        """

        list_axis_grid = [
            entry.parameters["bulk_planar_averages"][0][0],
            entry.parameters["bulk_planar_averages"][1][0],
            entry.parameters["bulk_planar_averages"][2][0],
            entry.parameters["defect_planar_averages"][0][0],
            entry.parameters["defect_planar_averages"][1][0],
            entry.parameters["defect_planar_averages"][2][0]
        ]

        list_axis_grid = []
        list_bulk_plnr_avg_esp = []
        list_defect_plnr_avg_esp = []
        list_axes = range(3)

        for axis in range(3):
            list_axis_grid.append(entry.parameters["bulk_planar_averages"][axis][0])
            list_bulk_plnr_avg_esp.append(entry.parameters["bulk_planar_averages"][axis][1])
            list_defect_plnr_avg_esp.append(entry.parameters["defect_planar_averages"][axis][1])

        lattice = entry.parameters["lattice"]
        dielectricconst = self.dielectric
        q = entry.defect.charge
        
        list_axes = list(range(3))

        self.dielectric_const = dielectric_const
        self.q_model = q_model
        self.energy_cutoff = energy_cutoff
        self.madelung_energy_tolerance = madelung_energy_tolerance

        es_corr = self.perform_es_corr(
            lattice,
            self.dielectric,
            self.defect.charge,
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
                self.defect.charge,
                entry.site.coords
                axis,
                q_model=self.q_model,
                madetol=self.madelung_energy_tolerance,
                title=None,
                widthsample=1.0)
            pot_corr_tracker.append(tmp_pot_corr)

        pot_corr = np.mean(pot_corr_tracker)

        return {"freysoldt_electrostatic": es_corr,
                "freysoldt_potential_alignment": pot_corr}


    def perform_es_corr(self, lattice, dielectricconst, q, energy_cutoff=520, q_model=QModel(), madetol=0.0001):
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
                    * q_model.rho_rec_limit0() / vol
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

        es_corr = round((eiso - eper) / dielectricconst * hart_to_ev, 6)
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
        v_G[0] = 4 * np.pi * -q / epsilon * q_model.rho_rec_limit0()
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

        #log uncertainty:
        if 'pot_corr_uncertainty_md' not in self.metadata:
            self.metadata['pot_corr_uncertainty_md'] = {axis: {'stats': stats.describe(tmppot), 'potcorr': -q * C}}
        else:
            self.metadata['pot_corr_uncertainty_md'][axis] = {'stats': stats.describe(tmppot), 'potcorr': -q * C}

        #make plot, if desired
        if title:
            plotter = FreysoldtCorrPlotter(axis_grid, v_R, defavg - pureavg, final_shift,
                                           np.array([mid - checkdis, mid + checkdis]))
            if title != 'written':
                plotter.plot(title=title)
            else:
                # TODO: Make this default fname more defect specific so it doesnt
                # over write previous defect data written
                fname = 'FreyAxisData'  # Extension is npz
                plotter.to_datafile(fname)

        return
