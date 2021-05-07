# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module implements DefectCompatibility analysis for consideration of
defects
"""

import logging

from monty.json import MSONable

from pymatgen.analysis.defects.core import Vacancy
from pymatgen.analysis.defects.corrections import (
    BandEdgeShiftingCorrection,
    BandFillingCorrection,
    FreysoldtCorrection,
    KumagaiCorrection,
)
from pymatgen.core import Structure

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

logger = logging.getLogger(__name__)


class DefectCompatibility(MSONable):
    """
    The DefectCompatibility class evaluates corrections and delocalization
    metrics on a DefectEntry. It can only parse based on the available
    parameters that already exist in the parameters dict of the DefectEntry.

    required settings in defect_entry.parameters for various types of analysis/correction:
        freysoldt: [ "dielectric", "axis_grid", "bulk_planar_averages", "defect_planar_averages",
                    "initial_defect_structure", "defect_frac_sc_coords"]
        kumagai: [ "dielectric", "bulk_atomic_site_averages", "defect_atomic_site_averages",
                   "site_matching_indices", "initial_defect_structure", "defect_frac_sc_coords"]
        bandfilling: ["eigenvalues", "kpoint_weights", "potalign", "vbm", "cbm"]
        bandshifting: ["hybrid_cbm", "hybrid_vbm", "vbm", "cbm"]
        defect relaxation/structure analysis: ["final_defect_structure", "initial_defect_structure",
                                              "sampling_radius", "defect_frac_sc_coords"]
    """

    def __init__(
        self,
        plnr_avg_var_tol=0.0001,
        plnr_avg_minmax_tol=0.1,
        atomic_site_var_tol=0.005,
        atomic_site_minmax_tol=0.1,
        tot_relax_tol=1.0,
        perc_relax_tol=50.0,
        defect_tot_relax_tol=2.0,
        preferred_cc="freysoldt",
        free_chg_cutoff=2.1,
        use_bandfilling=True,
        use_bandedgeshift=True,
    ):
        """
        Initializes the DefectCompatibility class

        Each argument helps decide whether a DefectEntry is flagged as compatible or not
        Args:
            plnr_avg_var_tol (float): compatibility tolerance for variance of the sampling region
                in the planar averaged electrostatic potential (FreysoldtCorrection)
            plnr_avg_minmax_tol (float): compatibility tolerance for max/min difference of the
                sampling region in the planar averaged electrostatic potential (FreysoldtCorrection)
            atomic_site_var_tol (float): compatibility tolerance for variance of the sampling
                region in the atomic site averaged electrostatic potential (KumagaiCorrection)
            atomic_site_minmax_tol (float): compatibility tolerance for max/min difference
                of the sampling region in the atomic site averaged electrostatic
                potential (KumagaiCorrection)
            tot_relax_tol (float): compatibility tolerance for total integrated relaxation
                amount outside of a given radius from the defect (in Angstrom).
                Radius is supplied as 'sampling_radius' within parameters of DefectEntry.
            perc_relax_tol (float): compatibility tolerance for percentage of total relaxation
                outside of a given radius from the defect (percentage amount),
                assuming a total integration relaxation greater than 1 Angstrom.
                Radius is supplied as 'sampling_radius' within parameters of DefectEntry.
            defect_tot_relax_tol (float): compatibility tolerance for displacement of defect site
                itself (in Angstrom).
            preferred_cc (str): Charge correction that is preferred to be used.
                If only one is available based on metadata, then that charge correction will be used.
                Options are: 'freysoldt' and 'kumagai'
            free_chg_cutoff (float): compatibility tolerance for total amount of host band occupation
                outside of band edges, given by eigenvalue data. Extra occupation in the CB would be
                free electrons, while lost occupation in VB would be free holes.
            use_bandfilling (bool): Whether to include BandFillingCorrection or not (assuming
                sufficient metadata is supplied to perform BandFillingCorrection).
            use_bandedgeshift (bool): Whether to perform a BandEdgeShiftingCorrection or not (assuming
                sufficient metadata is supplied to perform BandEdgeShiftingCorrection).
        """
        self.plnr_avg_var_tol = plnr_avg_var_tol
        self.plnr_avg_minmax_tol = plnr_avg_minmax_tol
        self.atomic_site_var_tol = atomic_site_var_tol
        self.atomic_site_minmax_tol = atomic_site_minmax_tol
        self.tot_relax_tol = tot_relax_tol
        self.perc_relax_tol = perc_relax_tol
        self.defect_tot_relax_tol = defect_tot_relax_tol

        self.preferred_cc = preferred_cc
        self.free_chg_cutoff = free_chg_cutoff
        self.use_bandfilling = use_bandfilling
        self.use_bandedgeshift = use_bandedgeshift

    def process_entry(self, defect_entry, perform_corrections=True):
        """
        Process a given DefectEntry with qualifiers given from initialization of class.
        Order of processing is:
            1) perform all possible defect corrections with information given
            2) consider delocalization analyses based on qualifier metrics
            given initialization of class. If delocalized, flag entry as delocalized
            3) update corrections to defect entry and flag as delocalized

        Corrections are applied based on:
            i) if free charges are more than free_chg_cutoff then will not apply charge correction,
                because it no longer is applicable
            ii) use charge correction set by preferred_cc
            iii) only use BandFilling correction if use_bandfilling is set to True
            iv) only use BandEdgeShift correction if use_bandedgeshift is set to True
        """
        for struct_key in [
            "bulk_sc_structure",
            "initial_defect_structure",
            "final_defect_structure",
        ]:
            if struct_key in defect_entry.parameters.keys() and isinstance(defect_entry.parameters[struct_key], dict):
                defect_entry.parameters[struct_key] = Structure.from_dict(defect_entry.parameters[struct_key])

        if perform_corrections:
            self.perform_all_corrections(defect_entry)

        self.delocalization_analysis(defect_entry)

        # apply corrections based on delocalization analysis
        corrections = {}
        skip_charge_corrections = False
        if "num_hole_vbm" in defect_entry.parameters.keys():
            if (self.free_chg_cutoff < defect_entry.parameters["num_hole_vbm"]) or (
                self.free_chg_cutoff < defect_entry.parameters["num_elec_cbm"]
            ):
                logger.info("Will not use charge correction because too many free charges")
                skip_charge_corrections = True

        if skip_charge_corrections:
            corrections.update({"charge_correction": 0.0})
        else:
            if ("freysoldt" in self.preferred_cc.lower()) and ("freysoldt_meta" in defect_entry.parameters.keys()):
                frey_meta = defect_entry.parameters["freysoldt_meta"]
                frey_corr = frey_meta["freysoldt_electrostatic"] + frey_meta["freysoldt_potential_alignment_correction"]
                corrections.update({"charge_correction": frey_corr})
            elif "kumagai_meta" in defect_entry.parameters.keys():
                kumagai_meta = defect_entry.parameters["kumagai_meta"]
                kumagai_corr = (
                    kumagai_meta["kumagai_electrostatic"] + kumagai_meta["kumagai_potential_alignment_correction"]
                )
                corrections.update({"charge_correction": kumagai_corr})
            else:
                logger.info("Could not use any charge correction because insufficient metadata was supplied.")

        if self.use_bandfilling:
            if "bandfilling_meta" in defect_entry.parameters.keys():
                bfc_corr = defect_entry.parameters["bandfilling_meta"]["bandfilling_correction"]
                corrections.update({"bandfilling_correction": bfc_corr})
            else:
                logger.info("Could not use band filling correction because insufficient metadata was supplied.")
        else:
            corrections.update({"bandfilling_correction": 0.0})

        if self.use_bandedgeshift and ("bandshift_meta" in defect_entry.parameters.keys()):
            corrections.update(
                {
                    "bandedgeshifting_correction": defect_entry.parameters["bandshift_meta"][
                        "bandedgeshifting_correction"
                    ]
                }
            )

            # also want to update relevant data for phase diagram
            defect_entry.parameters.update(
                {
                    "phasediagram_meta": {
                        "vbm": defect_entry.parameters["hybrid_vbm"],
                        "gap": defect_entry.parameters["hybrid_cbm"] - defect_entry.parameters["hybrid_vbm"],
                    }
                }
            )
        else:
            corrections.update({"bandedgeshifting_correction": 0.0})
            if isinstance(defect_entry.parameters["vbm"], float) and isinstance(defect_entry.parameters["cbm"], float):
                # still want to have vbm and gap ready for phase diagram
                defect_entry.parameters.update(
                    {
                        "phasediagram_meta": {
                            "vbm": defect_entry.parameters["vbm"],
                            "gap": defect_entry.parameters["cbm"] - defect_entry.parameters["vbm"],
                        }
                    }
                )

        defect_entry.corrections.update(corrections)

        return defect_entry

    def perform_all_corrections(self, defect_entry):
        """
        Perform all corrections for a defect.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        # consider running freysoldt correction
        required_frey_params = [
            "dielectric",
            "axis_grid",
            "bulk_planar_averages",
            "defect_planar_averages",
            "initial_defect_structure",
            "defect_frac_sc_coords",
        ]
        run_freysoldt = len(set(defect_entry.parameters.keys()).intersection(required_frey_params)) == len(
            required_frey_params
        )
        if not run_freysoldt:
            logger.info("Insufficient DefectEntry parameters exist for Freysoldt Correction.")
        else:
            defect_entry = self.perform_freysoldt(defect_entry)

        # consider running kumagai correction
        required_kumagai_params = [
            "dielectric",
            "bulk_atomic_site_averages",
            "defect_atomic_site_averages",
            "site_matching_indices",
            "initial_defect_structure",
            "defect_frac_sc_coords",
        ]
        run_kumagai = len(set(defect_entry.parameters.keys()).intersection(required_kumagai_params)) == len(
            required_kumagai_params
        )
        if not run_kumagai:
            logger.info("Insufficient DefectEntry parameters exist for Kumagai Correction.")
        else:
            try:
                defect_entry = self.perform_kumagai(defect_entry)
            except Exception:
                logger.info("Kumagai correction error occured! Wont perform correction.")

        # add potalign based on preferred correction setting if it does not already exist in defect entry
        if self.preferred_cc == "freysoldt":
            if "freysoldt_meta" in defect_entry.parameters.keys():
                potalign = defect_entry.parameters["freysoldt_meta"]["freysoldt_potalign"]
                defect_entry.parameters["potalign"] = potalign
            elif "kumagai_meta" in defect_entry.parameters.keys():
                logger.info(
                    "WARNING: was not able to use potalign from Freysoldt correction, "
                    "using Kumagai value for purposes of band filling correction."
                )
                potalign = defect_entry.parameters["kumagai_meta"]["kumagai_potalign"]
                defect_entry.parameters["potalign"] = potalign
        else:
            if "kumagai_meta" in defect_entry.parameters.keys():
                potalign = defect_entry.parameters["kumagai_meta"]["kumagai_potalign"]
                defect_entry.parameters["potalign"] = potalign
            elif "freysoldt_meta" in defect_entry.parameters.keys():
                logger.info(
                    "WARNING: was not able to use potalign from Kumagai correction, "
                    "using Freysoldt value for purposes of band filling correction."
                )
                potalign = defect_entry.parameters["freysoldt_meta"]["freysoldt_potalign"]
                defect_entry.parameters["potalign"] = potalign

        # consider running band filling correction
        required_bandfilling_params = [
            "eigenvalues",
            "kpoint_weights",
            "potalign",
            "vbm",
            "cbm",
        ]
        run_bandfilling = len(set(defect_entry.parameters.keys()).intersection(required_bandfilling_params)) == len(
            required_bandfilling_params
        )
        if run_bandfilling:
            if (
                (defect_entry.parameters["vbm"] is None)
                or (defect_entry.parameters["cbm"] is None)
                or (defect_entry.parameters["potalign"] is None)
            ):
                run_bandfilling = False

        if not run_bandfilling:
            logger.info("Insufficient DefectEntry parameters exist for BandFilling Correction.")
        else:
            defect_entry = self.perform_bandfilling(defect_entry)

        # consider running band edge shifting correction
        required_bandedge_shifting_params = ["hybrid_cbm", "hybrid_vbm", "vbm", "cbm"]
        run_bandedge_shifting = len(
            set(defect_entry.parameters.keys()).intersection(required_bandedge_shifting_params)
        ) == len(required_bandedge_shifting_params)
        if not run_bandedge_shifting:
            logger.info("Insufficient DefectEntry parameters exist for BandShifting Correction.")
        else:
            defect_entry = self.perform_band_edge_shifting(defect_entry)

        return defect_entry

    @staticmethod
    def perform_freysoldt(defect_entry):
        """
        Perform Freysoldt correction.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        FC = FreysoldtCorrection(defect_entry.parameters["dielectric"])
        freycorr = FC.get_correction(defect_entry)

        freysoldt_meta = FC.metadata.copy()
        freysoldt_meta["freysoldt_potalign"] = defect_entry.parameters["potalign"]
        freysoldt_meta["freysoldt_electrostatic"] = freycorr["freysoldt_electrostatic"]
        freysoldt_meta["freysoldt_potential_alignment_correction"] = freycorr["freysoldt_potential_alignment"]
        defect_entry.parameters.update({"freysoldt_meta": freysoldt_meta})
        return defect_entry

    @staticmethod
    def perform_kumagai(defect_entry):
        """
        Perform Kumagai correction.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        gamma = defect_entry.parameters["gamma"] if "gamma" in defect_entry.parameters.keys() else None
        sampling_radius = (
            defect_entry.parameters["sampling_radius"] if "sampling_radius" in defect_entry.parameters.keys() else None
        )

        KC = KumagaiCorrection(
            defect_entry.parameters["dielectric"],
            sampling_radius=sampling_radius,
            gamma=gamma,
        )
        kumagaicorr = KC.get_correction(defect_entry)

        kumagai_meta = dict(KC.metadata.items())
        kumagai_meta["kumagai_potalign"] = defect_entry.parameters["potalign"]
        kumagai_meta["kumagai_electrostatic"] = kumagaicorr["kumagai_electrostatic"]
        kumagai_meta["kumagai_potential_alignment_correction"] = kumagaicorr["kumagai_potential_alignment"]
        defect_entry.parameters.update({"kumagai_meta": kumagai_meta})
        return defect_entry

    @staticmethod
    def perform_bandfilling(defect_entry):
        """
        Perform bandfilling correction.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        BFC = BandFillingCorrection()
        bfc_dict = BFC.get_correction(defect_entry)

        bandfilling_meta = defect_entry.parameters["bandfilling_meta"].copy()
        bandfilling_meta.update({"bandfilling_correction": bfc_dict["bandfilling_correction"]})
        defect_entry.parameters.update(
            {
                "bandfilling_meta": bandfilling_meta,
                # also update free holes and electrons for shallow level shifting correction...
                "num_hole_vbm": bandfilling_meta["num_hole_vbm"],
                "num_elec_cbm": bandfilling_meta["num_elec_cbm"],
            }
        )
        return defect_entry

    @staticmethod
    def perform_band_edge_shifting(defect_entry):
        """
        Perform band edge shifting correction.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        BEC = BandEdgeShiftingCorrection()
        bec_dict = BEC.get_correction(defect_entry)

        bandshift_meta = defect_entry.parameters["bandshift_meta"].copy()
        bandshift_meta.update(bec_dict)

        defect_entry.parameters.update({"bandshift_meta": bandshift_meta})

        return defect_entry

    def delocalization_analysis(self, defect_entry):
        """
        Do delocalization analysis. To do this, one considers:
            i) sampling region of planar averaged electrostatic potential (freysoldt approach)
            ii) sampling region of atomic site averaged potentials (kumagai approach)
            iii) structural relaxation amount outside of radius considered in kumagai approach (default is wigner seitz
            radius)
            iv) if defect is not a vacancy type -> track to see how much the defect has moved

        calculations that fail delocalization get "is_compatibile" set to False in parameters
        also parameters recieves a "delocalization_meta" with following dict:
            plnr_avg = {'is_compatible': True/False, 'metadata': metadata used for determining this}
            atomic_site = {'is_compatible': True/False, 'metadata': metadata used for determining this}
            structure_relax = {'is_compatible': True/False, 'metadata': metadata used for determining this}
            defectsite_relax = {'is_compatible': True/False, 'metadata': metadata used for determing this}
        """
        defect_entry.parameters.update(
            {"is_compatible": True}
        )  # this will be switched to False if delocalization is detected

        if "freysoldt_meta" in defect_entry.parameters.keys():
            defect_entry = self.check_freysoldt_delocalized(defect_entry)
        else:
            logger.info(
                "Insufficient information provided for performing Freysoldt "
                "correction delocalization analysis.\n"
                "Cannot perform planar averaged electrostatic potential "
                "compatibility analysis."
            )

        if "kumagai_meta" in defect_entry.parameters.keys():
            defect_entry = self.check_kumagai_delocalized(defect_entry)
        else:
            logger.info(
                "Insufficient information provided for performing Kumagai "
                "correction delocalization analysis.\n"
                "Cannot perform atomic site averaged electrostatic "
                "potential compatibility analysis."
            )

        req_struct_delocal_params = [
            "final_defect_structure",
            "initial_defect_structure",
            "sampling_radius",
            "defect_frac_sc_coords",
        ]
        run_struct_delocal = len(set(defect_entry.parameters.keys()).intersection(req_struct_delocal_params)) == len(
            req_struct_delocal_params
        )
        if run_struct_delocal:
            defect_entry = self.check_final_relaxed_structure_delocalized(defect_entry)
        else:
            logger.info(
                "Insufficient information provided in defect_entry.parameters. "
                "Cannot perform full structure site relaxation compatibility analysis."
            )

        return defect_entry

    def check_freysoldt_delocalized(self, defect_entry):
        """
        Check for Freysoldt delocalization.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        plnr_avg_analyze_meta = {}
        plnr_avg_allows_compatible = True
        for ax in range(3):
            freystats = defect_entry.parameters["freysoldt_meta"]["pot_corr_uncertainty_md"][ax]["stats"]

            frey_variance_compatible = freystats["variance"] <= self.plnr_avg_var_tol
            frey_window = abs(freystats["minmax"][1] - freystats["minmax"][0])
            frey_minmax_compatible = frey_window <= self.plnr_avg_minmax_tol

            plnr_avg_analyze_meta.update(
                {
                    ax: {
                        "frey_variance_compatible": frey_variance_compatible,
                        "frey_variance": freystats["variance"],
                        "plnr_avg_var_tol": self.plnr_avg_var_tol,
                        "frey_minmax_compatible": frey_minmax_compatible,
                        "frey_minmax_window": frey_window,
                        "plnr_avg_minmax_tol": self.plnr_avg_minmax_tol,
                    }
                }
            )

            if (not frey_variance_compatible) or (not frey_minmax_compatible):
                plnr_avg_allows_compatible = False

        if "delocalization_meta" not in defect_entry.parameters.keys():
            defect_entry.parameters["delocalization_meta"] = {}
        defect_entry.parameters["delocalization_meta"].update(
            {
                "plnr_avg": {
                    "is_compatible": plnr_avg_allows_compatible,
                    "metadata": plnr_avg_analyze_meta,
                }
            }
        )

        if not plnr_avg_allows_compatible:
            defect_entry.parameters.update({"is_compatible": False})

        return defect_entry

    def check_kumagai_delocalized(self, defect_entry):
        """
        Check for Kumagai delocalization.

        Args:
            defect_entry (DefectEntry): Defect to correct.

        Returns:
            Corrected DefectEntry
        """
        atomic_site_analyze_meta = {}
        kumagaistats = defect_entry.parameters["kumagai_meta"]["pot_corr_uncertainty_md"]["stats"]

        kumagai_variance_compatible = kumagaistats["variance"] <= self.atomic_site_var_tol
        kumagai_window = abs(kumagaistats["minmax"][1] - kumagaistats["minmax"][0])
        kumagai_minmax_compatible = kumagai_window <= self.atomic_site_minmax_tol

        atomic_site_analyze_meta.update(
            {
                "kumagai_variance_compatible": kumagai_variance_compatible,
                "kumagai_variance": kumagaistats["variance"],
                "atomic_site_var_tol": self.atomic_site_var_tol,
                "kumagai_minmax_compatible": kumagai_minmax_compatible,
                "kumagai_minmax_window": kumagai_window,
                "plnr_avg_minmax_tol": self.atomic_site_minmax_tol,
            }
        )

        atomic_site_allows_compatible = kumagai_variance_compatible and kumagai_minmax_compatible

        if "delocalization_meta" not in defect_entry.parameters.keys():
            defect_entry.parameters["delocalization_meta"] = {}
        defect_entry.parameters["delocalization_meta"].update(
            {
                "atomic_site": {
                    "is_compatible": atomic_site_allows_compatible,
                    "metadata": atomic_site_analyze_meta,
                }
            }
        )

        if not atomic_site_allows_compatible:
            defect_entry.parameters.update({"is_compatible": False})

        return defect_entry

    def check_final_relaxed_structure_delocalized(self, defect_entry):
        """
        NOTE this assumes initial and final structures have sites indexed in same way

        :param defect_entry:
        :return:
        """
        structure_relax_analyze_meta = {}
        initial_defect_structure = defect_entry.parameters["initial_defect_structure"]
        final_defect_structure = defect_entry.parameters["final_defect_structure"]
        radius_to_sample = defect_entry.parameters["sampling_radius"]
        def_frac_coords = defect_entry.parameters["defect_frac_sc_coords"]

        initsites = [site.frac_coords for site in initial_defect_structure]
        finalsites = [site.frac_coords for site in final_defect_structure]
        distmatrix = initial_defect_structure.lattice.get_all_distances(finalsites, initsites)

        # calculate distance moved as a function of the distance from the defect
        distdata = []
        totpert = 0.0
        defindex = None
        for ind, site in enumerate(initial_defect_structure.sites):
            if site.distance_and_image_from_frac_coords(def_frac_coords)[0] < 0.01:
                defindex = ind
                continue

            totpert += distmatrix[ind, ind]
            # append [distance to defect, distance traveled, index in structure]
            distance_to_defect = initial_defect_structure.lattice.get_distance_and_image(
                def_frac_coords, initsites[ind]
            )[0]
            distdata.append([distance_to_defect, distmatrix[ind, ind], int(ind)])

        if defindex is None and not isinstance(defect_entry.defect, Vacancy):
            raise ValueError("fractional coordinate for defect could not be " "identified in initial_defect_structure")

        distdata.sort()
        tot_relax_outside_rad = 0.0
        perc_relax_outside_rad = 0.0
        for newind, d in enumerate(distdata):
            perc_relax = 100 * d[1] / totpert if totpert else 0.0
            d.append(perc_relax)  # percentage contribution to total relaxation
            if d[0] > radius_to_sample:
                tot_relax_outside_rad += d[1]
                perc_relax_outside_rad += d[3]

        structure_tot_relax_compatible = tot_relax_outside_rad <= self.tot_relax_tol
        structure_perc_relax_compatible = not (perc_relax_outside_rad > self.perc_relax_tol and totpert >= 1.0)
        structure_relax_analyze_meta.update(
            {
                "structure_tot_relax_compatible": structure_tot_relax_compatible,
                "tot_relax_outside_rad": tot_relax_outside_rad,
                "tot_relax_tol": self.tot_relax_tol,
                "structure_perc_relax_compatible": structure_perc_relax_compatible,
                "perc_relax_outside_rad": perc_relax_outside_rad,
                "perc_relax_tol": self.perc_relax_tol,
                "full_structure_relax_data": distdata,
                "defect_index": defindex,
            }
        )

        structure_relax_allows_compatible = structure_tot_relax_compatible and structure_perc_relax_compatible

        # NEXT: do single defect delocalization analysis (requires similar data, so might as well run in tandem
        # with structural delocalization)
        defectsite_relax_analyze_meta = {}
        if isinstance(defect_entry.defect, Vacancy):
            defectsite_relax_allows_compatible = True
            defectsite_relax_analyze_meta.update(
                {
                    "relax_amount": None,
                    "defect_tot_relax_tol": self.defect_tot_relax_tol,
                }
            )
        else:
            defect_relax_amount = distmatrix[defindex, defindex]
            defectsite_relax_allows_compatible = defect_relax_amount <= self.defect_tot_relax_tol
            defectsite_relax_analyze_meta.update(
                {
                    "relax_amount": defect_relax_amount,
                    "defect_tot_relax_tol": self.defect_tot_relax_tol,
                }
            )

        if "delocalization_meta" not in defect_entry.parameters.keys():
            defect_entry.parameters["delocalization_meta"] = {}
        defect_entry.parameters["delocalization_meta"].update(
            {
                "defectsite_relax": {
                    "is_compatible": defectsite_relax_allows_compatible,
                    "metadata": defectsite_relax_analyze_meta,
                }
            }
        )
        defect_entry.parameters["delocalization_meta"].update(
            {
                "structure_relax": {
                    "is_compatible": structure_relax_allows_compatible,
                    "metadata": structure_relax_analyze_meta,
                }
            }
        )

        if (not structure_relax_allows_compatible) or (not defectsite_relax_allows_compatible):
            defect_entry.parameters.update({"is_compatible": False})

        return defect_entry
