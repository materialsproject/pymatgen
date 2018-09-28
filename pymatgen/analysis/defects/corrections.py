# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
import numpy as np
import scipy
from scipy import stats
from pymatgen.analysis.defects.core import DefectCorrection
from pymatgen.analysis.defects.utils import ang_to_bohr, hart_to_ev, eV_to_k, \
    generate_reciprocal_vectors_squared, QModel, converge

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
    """

    def __init__(self, dielectric_const, q_model=None, energy_cutoff=520, madetol=0.0001, axis=None):
        """
        Initializes the Freysoldt Correction
        Args:
            dielectric_const (float or 3x3 matrix): Dielectric constant for the structure
            q_mode (QModel): instantiated QModel object or None. Uses default parameters to instantiate QModel if None supplied
            energy_cutoff (int): Maximum energy in eV in recipripcol space to perform integration for potential correction
            madeltol(float): Convergence criteria for the Madelung energy for potential correction
            axis (int): Axis to calculate correction. Averages over all three if not supplied.
        """
        self.q_model = QModel() if not q_model else q_model
        self.energy_cutoff = energy_cutoff
        self.madetol = madetol
        self.dielectric_const = dielectric_const

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
            entry (DefectEntry): defect entry to compute Freysoldt correction on.
                Requires following parameters in the DefectEntry to exist:

                    axis_grid (3 x NGX where NGX is the length of the NGX grid
                    in the x,y and z axis directions. Same length as planar
                    average lists):
                        A list of 3 numpy arrays which contain the cartesian axis
                        values (in angstroms) that correspond to each planar avg
                        potential supplied.

                    bulk_planar_averages (3 x NGX where NGX is the length of
                    the NGX grid in the x,y and z axis directions.):
                        A list of 3 numpy arrays which contain the planar averaged
                        electrostatic potential for the bulk supercell.

                    defect_planar_averages (3 x NGX where NGX is the length of
                    the NGX grid in the x,y and z axis directions.):
                        A list of 3 numpy arrays which contain the planar averaged
                        electrostatic potential for the defective supercell.

                    scaling_matrix (3 x 1 matrix): scaling matrix required to convert the
                        entry.defect.bulk_structure object into the lattice which is used by
                        the bulk_planar_average and defect_planar_average

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

        es_corr = self.perform_es_corr(lattice, entry.charge)

        pot_corr_tracker = []

        for x, pureavg, defavg, axis in zip(list_axis_grid, list_bulk_plnr_avg_esp, list_defect_plnr_avg_esp,
                                            list_axes):
            tmp_pot_corr = self.perform_pot_corr(
                x, pureavg, defavg, lattice, entry.charge, entry.site.coords, axis, widthsample=1.0)
            pot_corr_tracker.append(tmp_pot_corr)

        pot_corr = np.mean(pot_corr_tracker)

        entry.parameters["freysoldt_meta"] = dict(self.metadata)
        entry.parameters["potalign"] = pot_corr / (-q) if q else 0.

        return {"freysoldt_electrostatic": es_corr, "freysoldt_potential_alignment": pot_corr}

    def perform_es_corr(self, lattice, q, step=1e-4):
        """
        Peform Electrostatic Freysoldt Correction
        """
        logger.info("Running Freysoldt 2011 PC calculation (should be " "equivalent to sxdefectalign)")
        logger.debug("defect lattice constants are (in angstroms)" + str(lattice.abc))

        [a1, a2, a3] = ang_to_bohr * np.array(lattice.get_cartesian_coords(1))
        logging.debug("In atomic units, lat consts are (in bohr):" + str([a1, a2, a3]))
        vol = np.dot(a1, np.cross(a2, a3))  # vol in bohr^3

        def e_iso(encut):
            gcut = eV_to_k(encut)  # gcut is in units of 1/A
            return scipy.integrate.quad(lambda g: self.q_model.rho_rec(g * g)**2, step, gcut)[0] * (q**2) / np.pi

        def e_per(encut):
            eper = 0
            for g2 in generate_reciprocal_vectors_squared(a1, a2, a3, encut):
                eper += (self.q_model.rho_rec(g2)**2) / g2
            eper *= (q**2) * 2 * round(np.pi, 6) / vol
            eper += (q**2) * 4 * round(np.pi, 6) \
                * self.q_model.rho_rec_limit0 / vol
            return eper

        eiso = converge(e_iso, 5, self.madetol, self.energy_cutoff)
        logger.debug("Eisolated : %f", round(eiso, 5))

        eper = converge(e_per, 5, self.madetol, self.energy_cutoff)

        logger.info("Eperiodic : %f hartree", round(eper, 5))
        logger.info("difference (periodic-iso) is %f hartree", round(eper - eiso, 6))
        logger.info("difference in (eV) is %f", round((eper - eiso) * hart_to_ev, 4))

        es_corr = round((eiso - eper) / self.dielectric * hart_to_ev, 6)
        logger.info("Defect Correction without alignment %f (eV): ", es_corr)
        return es_corr

    def perform_pot_corr(self,
                         axis_grid,
                         pureavg,
                         defavg,
                         lattice,
                         q,
                         defect_position,
                         axis,
                         madetol=0.0001,
                         widthsample=1.0):
        """
        For performing planar averaging potential alignment

        title is for name of plot, if you dont want a plot then leave it as None
        widthsample is the width (in Angstroms) of the region in between defects where the potential alignment correction is averaged
        """
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
            for i in range(nx):
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

        # Build background charge potential with defect at origin
        v_G = np.empty(len(axis_grid), np.dtype("c16"))
        v_G[0] = 4 * np.pi * -q / self.dielectric * self.q_model.rho_rec_limit0
        g = np.roll(np.arange(-nx / 2, nx / 2, 1, dtype=int), int(nx / 2)) * dg
        g2 = np.multiply(g, g)[1:]
        v_G[1:] = 4 * np.pi / (self.dielectric * g2) * -q * self.q_model.rho_rec(g2)
        v_G[nx // 2] = 0 if not (nx % 2) else v_G[nx // 2]

        # Get the real space potential by peforming a  fft and grabbing the imaginary portion
        v_R = np.fft.fft(v_G)

        if abs(np.imag(v_R).max()) > self.madetol:
            raise Exception("imaginary part found to be %s", repr(np.imag(v_R).max()))
        v_R /= (lattice.volume * ang_to_bohr**3)
        v_R = np.real(v_R) * hart_to_ev

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

    def plot(self, axis, title=None, saved=False):
        """
        Plots the planar average electrostatic potential against the Long range and short range models from Freysoldt

        """

        x = self.metadata['pot_plot_data'][axis]['x']
        v_R = self.metadata['pot_plot_data'][axis]['Vr']
        dft_diff = self.metadata['pot_plot_data'][axis]['dft_diff']
        final_shift = self.metadata['pot_plot_data'][axis]['final_shift']
        check = self.metadata['pot_plot_data'][axis]['check']

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

    def __init__(self, resolution=0.01):
        """
        Initializes the Bandfilling correction

        Args:
            resolution (float): energy resolution to maintain for gap states
        """
        self.resolution = resolution
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
            spinfctr = 2.
        elif len(eigenvalues.keys()) == 2:
            spinfctr = 1.
        else:
            raise ValueError("Eigenvalue keys greater than 2")

        # for tracking mid gap states...
        resolution = self.resolution
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
        # TODO: add smarter defect level shifting based on defect level projection onto host bands
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
