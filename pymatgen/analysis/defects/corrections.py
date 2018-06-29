# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
import sys
import numpy as np
import scipy
from scipy import stats
from pymatgen.analysis.defects.core import DefectCorrection
from pymatgen.analysis.defects.utils import ang_to_bohr, hart_to_ev, eV_to_k, \
    generate_reciprocal_vectors_squared, QModel

import matplotlib.pyplot as plt

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

logger = logging.getLogger(__name__)


class FreysoldtCorrection(DefectCorrection):
    """
    A class for FreysoldtCorrection class. Largely adapated from PyCDT code
    Requires some parameters in the DefectEntry to properly function:
        # TODO: @dbroberg please give better descriptions of what these are
        axis_grid
            list of numpy arrays of x-axis values (in angstroms) corresponding to each avg esp supplied
        bulk_planar_averages
        defect_planar_averages
        scaling_matrix
    """

    def __init__(self, dielectric_const, q_model=None, energy_cutoff=520, madelung_energy_tolerance=0.0001, axis=None):
        """
        Initializes the Freysoldt Correction
        Args:
            dielectric_const (float or 3x3 matrix): Dielectric constant for the structure
            q_mode (QModel): instantiated QModel object or None. Uses default parameters to instantiate QModel if None supplied
            energy_cutoff (int): Maximum energy in eV in recipripcol space to perform integration for potential correction
            madelung_energy_tolerance(float): Convergence criteria for the Madelung energy for potential correction
            axis (int): Axis to calculate correction. Averages over all three if not supplied
        """
        self.q_model = QModel() if not q_model else q_model
        self.energy_cutoff = energy_cutoff
        self.madelung_energy_tolerance = madelung_energy_tolerance

        if isinstance(dielectric_const, int) or \
                isinstance(dielectric_const, float):
            self.dielectric = float(dielectric_const)
        else:
            self.dielectric = float(np.mean(np.diag(dielectric_const)))

        self.axis = axis

        self.metadata = {"pot_plot_data": {}, "pot_corr_uncertainty_md": {}}

    def get_correction(self, entry):
        """
        Gets the Freysoldt correction for a defect entry
        Args:
            entry (DefectEntry): defect entry to compute Freysoldt correction on
        """
        if not self.axis:
            list_axis_grid = np.array(entry.parameters["axis_grid"])
            list_bulk_plnr_avg_esp = np.array(entry.parameters["bulk_planar_averages"])
            list_defect_plnr_avg_esp = np.array(entry.parameters["defect_planar_averages"])
            list_axes = range(len(list_axis_grid))
        else:
            list_axes = np.array(self.axis)
            list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp = [], [], []
            for ax in list_axes:
                list_axis_grid.append(np.array(entry.parameters["axis_grid"][ax]))
                list_bulk_plnr_avg_esp.append(np.array(entry.parameters["bulk_planar_averages"][ax]))
                list_defect_plnr_avg_esp.append(np.array(entry.parameters["defect_planar_averages"][ax]))

        bulk_struct = entry.defect.bulk_structure.copy()
        if "scaling_matrix" in entry.parameters.keys():
            bulk_struct.make_supercell(entry.parameters["scaling_matrix"])

        lattice = bulk_struct.lattice
        q = entry.defect.charge

        es_corr = self.perform_es_corr(
            lattice,
            self.dielectric,
            entry.charge,
            energy_cutoff=self.energy_cutoff,
            q_model=self.q_model,
            madetol=self.madelung_energy_tolerance)

        pot_corr_tracker = []

        for x, pureavg, defavg, axis in zip(list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp,
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

        entry.parameters["freysoldt_meta"] = dict(self.metadata)
        entry.parameters["potalign"] = pot_corr / (-q) if q else 0.

        return {"freysoldt_electrostatic": es_corr, "freysoldt_potential_alignment": pot_corr}

    def perform_es_corr(self, lattice, dielectric, q, energy_cutoff=520, q_model=QModel(), madetol=0.0001):
        """
        Peform Electrostatic Freysoldt Correction
        """
        logger = logging.getLogger(__name__)
        logger.info("Running Freysoldt 2011 PC calculation (should be " "equivalent to sxdefectalign)")
        logger.debug("defect lattice constants are (in angstroms)" + str(lattice.abc))

        [a1, a2, a3] = ang_to_bohr * np.array(lattice.get_cartesian_coords(1))
        logging.debug("In atomic units, lat consts are (in bohr):" + str([a1, a2, a3]))
        vol = np.dot(a1, np.cross(a2, a3))  # vol in bohr^3

        # compute isolated energy
        step = 1e-4
        encut1 = 20  # converge to some smaller encut first [eV]
        flag = 0
        converge = []
        while (flag != 1):
            eiso = 1.
            gcut = eV_to_k(encut1)  # gcut is in units of 1/A
            g = step  # initalize
            while g < (gcut + step):
                # simpson integration
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
                    raise Exception("Eiso did not converge before " + str(energy_cutoff) + " eV")
            encut1 += 20
        eiso = converge[-1]
        logger.debug("Eisolated : %f, converged at encut: %d", round(eiso, 5), encut1 - 20)

        # compute periodic energy;
        encut1 = 20  # converge to some smaller encut
        flag = 0
        converge = []
        while flag != 1:
            eper = 0.0
            for g2 in generate_reciprocal_vectors_squared(a1, a2, a3, encut1):
                eper += (q_model.rho_rec(g2)**2) / g2
            eper *= (q**2) * 2 * round(np.pi, 6) / vol
            eper += (q**2) * 4 * round(np.pi, 6) \
                * q_model.rho_rec_limit0 / vol
            converge.append(eper)
            if len(converge) > 2:
                if abs(converge[-1] - converge[-2]) < madetol:
                    flag = 1
                elif encut1 > energy_cutoff:
                    logger.error("Eper did not converge before %d eV", energy_cutoff)
                    return
            encut1 += 20
        eper = converge[-1]

        logger.info("Eperiodic : %f hartree, converged at encut %d eV", round(eper, 5), encut1 - 20)
        logger.info("difference (periodic-iso) is %f hartree", round(eper - eiso, 6))
        logger.info("difference in (eV) is %f", round((eper - eiso) * hart_to_ev, 4))

        es_corr = round((eiso - eper) / dielectric * hart_to_ev, 6)
        logger.info("Defect Correction without alignment %f (eV): ", es_corr)
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

        logging.debug("run Freysoldt potential alignment method for axis " + str(axis))
        nx = len(axis_grid)

        # shift these planar averages to have defect at origin
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

        # if not self._silence:
        logger.debug("calculating lr part along planar avg axis")
        reci_latt = lattice.reciprocal_lattice
        dg = reci_latt.abc[axis]
        dg /= ang_to_bohr  # convert to bohr to do calculation in atomic units

        v_G = np.empty(len(axis_grid), np.dtype("c16"))
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
            logging.error("imaginary part found to be %s", repr(max_imag_vr))
            sys.exit()

        # get correction
        short = (defavg - pureavg - v_R)
        checkdis = int((widthsample / 2) / (axis_grid[1] - axis_grid[0]))
        mid = int(len(short) / 2)

        tmppot = [short[i] for i in range(mid - checkdis, mid + checkdis + 1)]
        logger.debug("shifted defect position on axis (%s) to origin", repr(axbulkval))
        logger.debug("means sampling region is (%f,%f)", axis_grid[mid - checkdis], axis_grid[mid + checkdis])

        C = -np.mean(tmppot)
        logger.debug("C = %f", C)
        final_shift = [short[j] + C for j in range(len(v_R))]
        v_R = [elmnt - C for elmnt in v_R]

        logger.info("C value is averaged to be %f eV ", C)
        logger.info("Potentital alignment energy correction (-q*delta V):  %f (eV)", -q * C)
        self.pot_corr = -q * C

        # log plotting data:
        self.metadata["pot_plot_data"][axis] = {
            "Vr": v_R,
            "x": axis_grid,
            "dft_diff": defavg - pureavg,
            "final_shift": final_shift,
            "check": [mid - checkdis, mid + checkdis + 1]
        }

        # log uncertainty:
        self.metadata["pot_corr_uncertainty_md"][axis] = {"stats": stats.describe(tmppot)._asdict(), "potcorr": -q * C}

        return self.pot_corr


# TODO: Make this a part of above
def freysoldt_plotter(x, v_R, dft_diff, final_shift, check, title=None, saved=False):
    """
    planar average electrostatic potential plotter for Freysoldt
    """
    plt.figure()
    plt.clf()
    plt.plot(x, v_R, c="green", zorder=1, label="long range from model")
    plt.plot(x, dft_diff, c="red", label="DFT locpot diff")
    plt.plot(x, final_shift, c="blue", label="short range (aligned)")

    tmpx = [x[i] for i in range(check[0], check[1])]
    plt.fill_between(tmpx, -100, 100, facecolor="red", alpha=0.15, label="sampling region")

    plt.xlim(round(x[0]), round(x[-1]))
    ymin = min(min(v_R), min(dft_diff), min(final_shift))
    ymax = max(max(v_R), max(dft_diff), max(final_shift))
    plt.ylim(-0.2 + ymin, 0.2 + ymax)
    plt.xlabel("distance along axis ($\AA$)", fontsize=15)
    plt.ylabel("Potential (V)", fontsize=15)
    plt.legend(loc=9)
    plt.axhline(y=0, linewidth=0.2, color="black")
    plt.title(str(title) + " defect potential", fontsize=18)
    plt.xlim(0, max(x))
    if saved:
        plt.savefig(str(title) + "FreyplnravgPlot.pdf")
    else:
        return plt


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
        self.metadata = {
            "occupied_def_levels": [],
            "unoccupied_def_levels": [],
            "total_occupation_defect_levels": None,
            "num_hole_vbm": None,
            "num_elec_cbm": None,
            "potalign": None
        }

    def get_correction(self, entry):
        """
        Gets the BandFilling correction for a defect entry
        """
        eigenvalues = entry.parameters["eigenvalues"]
        kpoint_weights = entry.parameters["kpoint_weights"]
        potalign = entry.parameters["potalign"]
        vbm = entry.parameters["vbm"]
        cbm = entry.parameters["cbm"]

        bf_corr = self.perform_bandfill_corr(eigenvalues, kpoint_weights, potalign, vbm, cbm)

        entry.parameters["bandfilling_meta"] = dict(self.metadata)

        return {"bandfilling": bf_corr}

    def perform_bandfill_corr(self, eigenvalues, kpoint_weights, potalign, vbm, cbm):
        """
        This calculates the band filling correction based on excess of electrons/holes in CB/VB...

        Note that the total free holes and electrons may also be used for a "shallow donor/acceptor"
               correction with specified band shifts: +num_elec_cbm * Delta E_CBM (or -num_hole_vbm * Delta E_VBM)
               [this is done in the LevelShiftingCorrection class]
        """
        bf_corr = 0.

        self.metadata["potalign"] = potalign
        self.metadata["num_hole_vbm"] = 0.
        self.metadata["num_elec_cbm"] = 0.

        if len(eigenvalues.keys()) == 1:  # needed because occupation of non-spin calcs is still 1... should be 2
            spinfctr = 1.
        elif len(eigenvalues.keys()) == 2:
            spinfctr = 2.
        else:
            raise ValueError("Eigenvalue keys greater than 2")

        # for tracking mid gap states...
        resolution = 0.01  # this is energy resolution to maintain for gap states
        shifted_cbm = potalign + cbm  # shift cbm with potential alignment
        shifted_vbm = potalign + vbm  # shift vbm with potential alignment

        occupied_midgap = {en: [] for en in np.arange(shifted_vbm, shifted_cbm + resolution, resolution)}
        occupation = {en: 0. for en in np.arange(shifted_vbm, shifted_cbm + resolution, resolution)}
        unoccupied_midgap = {en: [] for en in np.arange(shifted_vbm, shifted_cbm + resolution, resolution)}

        for spinset in eigenvalues.values():
            for kptset, weight in zip(spinset, kpoint_weights):
                for eig in kptset:  # eig[0] is eigenvalue and eig[1] is occupation
                    if (eig[1] and (eig[0] > shifted_cbm)):  # donor MB correction
                        bf_corr += weight * spinfctr * eig[1] * (eig[0] - shifted_cbm)  # "move the electrons down"
                        self.metadata["num_elec_cbm"] += weight * spinfctr * eig[1]
                    elif (eig[1] != 1.) and (eig[0] <= shifted_vbm):  # acceptor MB correction
                        bf_corr += weight * spinfctr * (1. - eig[1]) * (shifted_vbm - eig[0])  # "move the holes up"
                        self.metadata["num_hole_vbm"] += weight * spinfctr * (1. - eig[1])
                    elif (eig[0] > shifted_vbm) and (eig[0] < shifted_cbm):
                        for en in np.arange(shifted_vbm, shifted_cbm + resolution, resolution):
                            if (eig[0] < en + resolution) and (eig[0] > en):
                                if eig[1]:
                                    occupied_midgap[en].append(eig[0])
                                    occupation[en] += eig[1] * weight * spinfctr
                                else:
                                    unoccupied_midgap[en].append(eig[0])
                                continue

        bf_corr *= -1  # need to take negative of this shift for energetic correction

        # summarize defect level results
        self.metadata["total_occupation_defect_levels"] = 0.
        self.metadata["occupied_def_levels"] = []
        self.metadata["unoccupied_def_levels"] = []
        for en in occupied_midgap.keys():
            if len(occupied_midgap[en]):
                self.metadata["occupied_def_levels"].append([np.mean(occupied_midgap[en]), occupation[en]])
                self.metadata["total_occupation_defect_levels"] += occupation[en]
            elif len(unoccupied_midgap[en]):
                self.metadata["unoccupied_def_levels"].append(np.mean(unoccupied_midgap[en]))

        return bf_corr


class BandEdgeShiftingCorrection(DefectCorrection):
    """
    A class for BandEdgeShiftingCorrection class. Largely adapted from PyCDT code

    Requires some parameters in the DefectEntry to properly function:
        hybrid_cbm
            CBM of HYBRID bulk calculation

        hybrid_vbm
            VBM of HYBRID bulk calculation

        cbm
            CBM of bulk calculation (or band structure calculation of bulk);
            calculated on same level of theory as the eigenvalues list (ex. GGA defects -> need GGA cbm

        vbm
            VBM of bulk calculation (or band structure calculation of bulk);
            calculated on same level of theory as the eigenvalues list (ex. GGA defects -> need GGA vbm

        num_hole_vbm
            number of free holes that were found in valence band for the defect calculation
            calculated in the metadata of the BandFilling Correction

        num_elec_cbm
            number of free electrons that were found in the conduction band for the defect calculation
            calculated in the metadata of the BandFilling Correction

    """

    def __init__(self):
        self.metadata = {
            "vbmshift": 0.,
            "cbmshift": 0.,
        }

    def get_correction(self, entry):
        """
        Gets the BandEdge correction for a defect entry
        """
        # TODO: in future will perform a defect level shift with projection method from kyle
        hybrid_cbm = entry.parameters["hybrid_cbm"]
        hybrid_vbm = entry.parameters["hybrid_vbm"]
        vbm = entry.parameters["vbm"]
        cbm = entry.parameters["cbm"]
        num_hole_vbm = entry.parameters["num_hole_vbm"]
        num_elec_cbm = entry.parameters["num_elec_cbm"]

        self.metadata["vbmshift"] = hybrid_vbm - vbm  # note vbmshift has UPWARD as positive convention
        self.metadata["cbmshift"] = hybrid_cbm - cbm  # note cbmshift has UPWARD as positive convention

        charge = entry.charge
        vbm_shift_correction = charge * self.metadata["vbmshift"]
        # negative sign has to do with fact that these are holes
        hole_vbm_shift_correction = -1. * num_hole_vbm * self.metadata["vbmshift"]
        elec_cbm_shift_correction = num_elec_cbm * self.metadata["cbmshift"]

        entry.parameters["bandshift_meta"] = dict(self.metadata)

        return {
            "vbm_shift_correction": vbm_shift_correction,
            "hole_vbm_shift_correction": hole_vbm_shift_correction,
            "elec_cbm_shift_correction": elec_cbm_shift_correction
        }
