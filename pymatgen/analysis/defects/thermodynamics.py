# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Defect thermodynamics, such as defect phase diagrams, etc.
"""

import logging
import warnings
import copy
import json
import os
import numpy as np
from monty.json import MSONable
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import bisect
from scipy.constants.codata import value as _cd
from scipy.interpolate import interp1d
from itertools import chain
from typing import List, Union, Dict

from pymatgen.electronic_structure.dos import CompleteDos, FermiDos, f0
from pymatgen.analysis.defects.core import DefectEntry
from pymatgen.analysis.structure_matcher import PointDefectComparator
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import (
    GibbsComputedStructureEntry, CompositionEnergyAdjustment, ComputedStructureEntry
)
from pymatgen.core import Element

import matplotlib.pyplot as plt
import matplotlib.cm as cm

boltzmann_eV = _cd("Boltzmann constant in eV/K")

with open(os.path.join(os.path.dirname(__file__), '../../entries/data/g_els.json')) as f:
    G_ELEMS = json.load(f)
with open(os.path.join(os.path.dirname(__file__), '../../entries/data/nist_gas_gf.json')) as f:
    G_GASES = json.load(f)

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

logger = logging.getLogger(__name__)


class DefectPhaseDiagram(MSONable):
    """
    This is similar to a PhaseDiagram object in pymatgen,
    but has ability to do quick analysis of defect formation energies
    when fed DefectEntry objects.

    uses many of the capabilities from PyCDT's DefectsAnalyzer class...

    This class is able to get:
        a) stability of charge states for a given defect,
        b) list of all formation ens
        c) transition levels in the gap
    """

    def __init__(self, entries, vbm, band_gap, filter_compatible=True, metadata=None):
        """
        Args:
            dentries ([DefectEntry]): A list of DefectEntry objects
            vbm (float): Valence Band energy to use for all defect entries.
                NOTE if using band shifting-type correction then this VBM
                should still be that of the GGA calculation
                (the bandedgeshifting_correction accounts for shift's
                contribution to formation energy).
            band_gap (float): Band gap to use for all defect entries.
                NOTE if using band shifting-type correction then this gap
                should still be that of the Hybrid calculation you are shifting to.
            filter_compatible (bool): Whether to consider entries which were ruled
                incompatible by the DefectComaptibility class. Note this must be set to False
                if you desire a suggestion for larger supercell sizes.
                Default is True (to omit calculations which have "is_compatible"=False in
                    DefectEntry'sparameters)
            metadata (dict): Dictionary of metadata to store with the PhaseDiagram. Has
                no impact on calculations
        """
        self.vbm = vbm
        self.band_gap = band_gap
        self.filter_compatible = filter_compatible

        if filter_compatible:
            self.entries = [e for e in entries if e.parameters.get("is_compatible", True)]
        else:
            self.entries = entries

        for ent_ind, ent in enumerate(self.entries):
            if "vbm" not in ent.parameters.keys() or ent.parameters["vbm"] != vbm:
                logger.info(
                    "Entry {} did not have vbm equal to given DefectPhaseDiagram value."
                    " Manually overriding.".format(ent.name)
                )
                new_ent = ent.copy()
                new_ent.parameters["vbm"] = vbm
                self.entries[ent_ind] = new_ent

        self.metadata = metadata or {}
        self.find_stable_charges()

    def as_dict(self):
        """
        Returns:
            Json-serializable dict representation of DefectPhaseDiagram
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "entries": [entry.as_dict() for entry in self.entries],
            "vbm": self.vbm,
            "band_gap": self.band_gap,
            "filter_compatible": self.filter_compatible,
            "metadata": self.metadata,
        }
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a DefectPhaseDiagram object from a dict representation created using
        as_dict().

        Args:
            d (dict): dict representation of DefectPhaseDiagram.

        Returns:
            DefectPhaseDiagram object
        """
        entries = [DefectEntry.from_dict(entry_dict) for entry_dict in d.get("entries")]
        vbm = d["vbm"]
        band_gap = d["band_gap"]
        filter_compatible = d.get("filter_compatible", True)
        metadata = d.get("metadata", {})
        if "entry_id" in d.keys() and "entry_id" not in metadata:
            metadata["entry_id"] = d["entry_id"]

        return cls(
            entries,
            vbm,
            band_gap,
            filter_compatible=filter_compatible,
            metadata=metadata,
        )

    def find_stable_charges(self):
        """
        Sets the stable charges and transition states for a series of
        defect entries. This function uses scipy's HalfspaceInterection
        to oncstruct the polygons corresponding to defect stability as
        a function of the Fermi-level. The Halfspace Intersection
        constructs N-dimensional hyperplanes, in this case N=2,  based
        on the equation of defect formation energy with considering chemical
        potentials:
            E_form = E_0^{Corrected} + Q_{defect}*(E_{VBM} + E_{Fermi})

        Extra hyperplanes are constructed to bound this space so that
        the algorithm can actually find enclosed region.

        This code was modeled after the Halfspace Intersection code for
        the Pourbaix Diagram
        """

        def similar_defects(entryset):
            """
            Used for grouping similar defects of different charges
            Can distinguish identical defects even if they are not in same position
            """
            pdc = PointDefectComparator(check_charge=False, check_primitive_cell=True, check_lattice_scale=False)
            grp_def_sets = []
            grp_def_indices = []
            for ent_ind, ent in enumerate(entryset):
                # TODO: more pythonic way of grouping entry sets with PointDefectComparator.
                # this is currently most time intensive part of DefectPhaseDiagram
                matched_ind = None
                for grp_ind, defgrp in enumerate(grp_def_sets):
                    if pdc.are_equal(ent.defect, defgrp[0].defect):
                        matched_ind = grp_ind
                        break
                if matched_ind is not None:
                    grp_def_sets[matched_ind].append(ent.copy())
                    grp_def_indices[matched_ind].append(ent_ind)
                else:
                    grp_def_sets.append([ent.copy()])
                    grp_def_indices.append([ent_ind])

            return zip(grp_def_sets, grp_def_indices)

        # Limits for search
        # E_fermi = { -1 eV to band gap+1}
        # E_formation = { (min(Eform) - 30) to (max(Eform) + 30)}
        all_eform = [one_def.formation_energy(fermi_level=self.band_gap / 2.0) for one_def in self.entries]
        min_y_lim = min(all_eform) - 30
        max_y_lim = max(all_eform) + 30
        limits = [[-1, self.band_gap + 1], [min_y_lim, max_y_lim]]

        stable_entries = {}
        finished_charges = {}
        transition_level_map = {}

        # Grouping by defect types
        for defects, index_list in similar_defects(self.entries):
            defects = list(defects)

            # prepping coefficient matrix for half-space intersection
            # [-Q, 1, -1*(E_form+Q*VBM)] -> -Q*E_fermi+E+-1*(E_form+Q*VBM) <= 0  where E_fermi and E are the variables
            # in the hyperplanes
            hyperplanes = np.array(
                [
                    [
                        -1.0 * entry.charge,
                        1,
                        -1.0 * (entry.energy + entry.charge * self.vbm),
                    ]
                    for entry in defects
                ]
            )

            border_hyperplanes = [
                [-1, 0, limits[0][0]],
                [1, 0, -1 * limits[0][1]],
                [0, -1, limits[1][0]],
                [0, 1, -1 * limits[1][1]],
            ]
            hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])

            interior_point = [self.band_gap / 2, min(all_eform) - 1.0]

            hs_ints = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

            # Group the intersections and coresponding facets
            ints_and_facets = zip(hs_ints.intersections, hs_ints.dual_facets)
            # Only inlcude the facets corresponding to entries, not the boundaries
            total_entries = len(defects)
            ints_and_facets = filter(
                lambda int_and_facet: all(np.array(int_and_facet[1]) < total_entries),
                ints_and_facets,
            )
            # sort based on transition level
            ints_and_facets = list(sorted(ints_and_facets, key=lambda int_and_facet: int_and_facet[0][0]))

            # log a defect name for tracking (using full index list to avoid naming
            # in-equivalent defects with same name)
            str_index_list = [str(ind) for ind in sorted(index_list)]
            track_name = defects[0].name + "@" + str("-".join(str_index_list))

            if len(ints_and_facets):
                # Unpack into lists
                _, facets = zip(*ints_and_facets)
                # Map of transition level: charge states

                transition_level_map[track_name] = {
                    intersection[0]: [defects[i].charge for i in facet] for intersection, facet in ints_and_facets
                }

                stable_entries[track_name] = list({defects[i] for dual in facets for i in dual})

                finished_charges[track_name] = [defect.charge for defect in defects]
            else:
                # if ints_and_facets is empty, then there is likely only one defect...
                if len(defects) != 1:
                    # confirm formation energies dominant for one defect over other identical defects
                    name_set = [one_def.name + "_chg" + str(one_def.charge) for one_def in defects]
                    vb_list = [one_def.formation_energy(fermi_level=limits[0][0]) for one_def in defects]
                    cb_list = [one_def.formation_energy(fermi_level=limits[0][1]) for one_def in defects]

                    vbm_def_index = vb_list.index(min(vb_list))
                    name_stable_below_vbm = name_set[vbm_def_index]
                    cbm_def_index = cb_list.index(min(cb_list))
                    name_stable_above_cbm = name_set[cbm_def_index]

                    if name_stable_below_vbm != name_stable_above_cbm:
                        raise ValueError(
                            "HalfSpace identified only one stable charge out of list: {}\n"
                            "But {} is stable below vbm and {} is "
                            "stable above cbm.\nList of VBM formation energies: {}\n"
                            "List of CBM formation energies: {}"
                            "".format(
                                name_set,
                                name_stable_below_vbm,
                                name_stable_above_cbm,
                                vb_list,
                                cb_list,
                            )
                        )
                    logger.info("{} is only stable defect out of {}".format(name_stable_below_vbm, name_set))
                    transition_level_map[track_name] = {}
                    stable_entries[track_name] = list([defects[vbm_def_index]])
                    finished_charges[track_name] = [one_def.charge for one_def in defects]
                else:
                    transition_level_map[track_name] = {}

                    stable_entries[track_name] = list([defects[0]])

                    finished_charges[track_name] = [defects[0].charge]

        self.transition_level_map = transition_level_map
        self.transition_levels = {
            defect_name: list(defect_tls.keys()) for defect_name, defect_tls in transition_level_map.items()
        }
        self.stable_entries = stable_entries
        self.finished_charges = finished_charges
        self.stable_charges = {
            defect_name: [entry.charge for entry in entries] for defect_name, entries in stable_entries.items()
        }

    @property
    def defect_types(self):
        """
        List types of defects existing in the DefectPhaseDiagram
        """
        return list(self.finished_charges.keys())

    @property
    def all_stable_entries(self):
        """
        List all stable entries (defect+charge) in the DefectPhaseDiagram
        """
        return set(chain.from_iterable(self.stable_entries.values()))

    @property
    def all_unstable_entries(self):
        """
        List all unstable entries (defect+charge) in the DefectPhaseDiagram
        """
        all_stable_entries = self.all_stable_entries
        return [e for e in self.entries if e not in all_stable_entries]

    def defect_concentrations(self, chemical_potentials, temperature=300, fermi_level=0.0):
        """
        Give list of all concentrations at specified efermi in the DefectPhaseDiagram
        args:
            chemical_potentials = {Element: number} is dict of chemical potentials to provide formation energies for
            temperature = temperature to produce concentrations from
            fermi_level: (float) is fermi level relative to valence band maximum
                Default efermi = 0 = VBM energy
        returns:
            list of dictionaries of defect concentrations
        """
        concentrations = []
        for dfct in self.all_stable_entries:
            concentrations.append(
                {
                    "conc": dfct.defect_concentration(
                        chemical_potentials=chemical_potentials,
                        temperature=temperature,
                        fermi_level=fermi_level,
                    ),
                    "name": dfct.name,
                    "charge": dfct.charge,
                }
            )

        return concentrations

    def suggest_charges(self, tolerance=0.1):
        """
        Suggest possible charges for defects to compute based on proximity
        of known transitions from entires to VBM and CBM

        Args:
            tolerance (float): tolerance with respect to the VBM and CBM to
                    `          continue to compute new charges
        """
        recommendations = {}

        for def_type in self.defect_types:
            test_charges = np.arange(
                np.min(self.stable_charges[def_type]) - 1,
                np.max(self.stable_charges[def_type]) + 2,
            )
            test_charges = [charge for charge in test_charges if charge not in self.finished_charges[def_type]]

            if len(self.transition_level_map[def_type].keys()):
                # More positive charges will shift the minimum transition level down
                # Max charge is limited by this if its transition level is close to VBM
                min_tl = min(self.transition_level_map[def_type].keys())
                if min_tl < tolerance:
                    max_charge = max(self.transition_level_map[def_type][min_tl])
                    test_charges = [charge for charge in test_charges if charge < max_charge]

                # More negative charges will shift the maximum transition level up
                # Minimum charge is limited by this if transition level is near CBM
                max_tl = max(self.transition_level_map[def_type].keys())
                if max_tl > (self.band_gap - tolerance):
                    min_charge = min(self.transition_level_map[def_type][max_tl])
                    test_charges = [charge for charge in test_charges if charge > min_charge]
            else:
                test_charges = [charge for charge in test_charges if charge not in self.stable_charges[def_type]]

            recommendations[def_type] = test_charges

        return recommendations

    def suggest_larger_supercells(self, tolerance=0.1):
        """
        Suggest larger supercells for different defect+chg combinations based on use of
        compatibility analysis. Does this for any charged defects which have is_compatible = False,
        and the defect+chg formation energy is stable at fermi levels within the band gap.

        NOTE: Requires self.filter_compatible = False
        Args:
            tolerance (float): tolerance with respect to the VBM and CBM for considering
                               larger supercells for a given charge
        """
        if self.filter_compatible:
            raise ValueError("Cannot suggest larger supercells if filter_compatible is True.")

        recommendations = {}

        for def_type in self.defect_types:
            template_entry = self.stable_entries[def_type][0].copy()
            defect_indices = [int(def_ind) for def_ind in def_type.split("@")[-1].split("-")]

            for charge in self.finished_charges[def_type]:
                chg_defect = template_entry.defect.copy()
                chg_defect.set_charge(charge)

                for entry_index in defect_indices:
                    entry = self.entries[entry_index]
                    if entry.charge == charge:
                        break

                if entry.parameters.get("is_compatible", True):
                    continue

                # consider if transition level is within
                # tolerance of band edges
                suggest_bigger_supercell = True
                for tl, chgset in self.transition_level_map[def_type].items():
                    sorted_chgset = list(chgset)
                    sorted_chgset.sort(reverse=True)
                    if charge == sorted_chgset[0] and tl < tolerance:
                        suggest_bigger_supercell = False
                    elif charge == sorted_chgset[1] and tl > (self.band_gap - tolerance):
                        suggest_bigger_supercell = False

                if suggest_bigger_supercell:
                    if def_type not in recommendations:
                        recommendations[def_type] = []
                    recommendations[def_type].append(charge)

        return recommendations

    def solve_for_fermi_energy(self, temperature, chemical_potentials, bulk_dos):
        """
        Solve for the Fermi energy self-consistently as a function of T
        Observations are Defect concentrations, electron and hole conc
        Args:
            temperature: Temperature to equilibrate fermi energies for
            chemical_potentials: dict of chemical potentials to use for calculation fermi level
            bulk_dos: bulk system dos (pymatgen Dos object)
        Returns:
            Fermi energy dictated by charge neutrality
        """

        fdos = FermiDos(bulk_dos, bandgap=self.band_gap)
        _, fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):
            qd_tot = sum(
                [
                    d["charge"] * d["conc"]
                    for d in self.defect_concentrations(
                        chemical_potentials=chemical_potentials,
                        temperature=temperature,
                        fermi_level=ef,
                    )
                ]
            )
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)
            return qd_tot

        return bisect(_get_total_q, -1.0, self.band_gap + 1.0)

    def solve_for_non_equilibrium_fermi_energy(self, temperature, quench_temperature, chemical_potentials, bulk_dos):
        """
        Solve for the Fermi energy after quenching in the defect concentrations at a higher
        temperature (the quench temperature),
        as outlined in P. Canepa et al (2017) Chemistry of Materials (doi: 10.1021/acs.chemmater.7b02909)

        Args:
            temperature: Temperature to equilibrate fermi energy at after quenching in defects
            quench_temperature: Temperature to equilibrate defect concentrations at (higher temperature)
            chemical_potentials: dict of chemical potentials to use for calculation fermi level
            bulk_dos: bulk system dos (pymatgen Dos object)
        Returns:
            Fermi energy dictated by charge neutrality with respect to frozen in defect concentrations
        """

        high_temp_fermi_level = self.solve_for_fermi_energy(quench_temperature, chemical_potentials, bulk_dos)
        fixed_defect_charge = sum(
            [
                d["charge"] * d["conc"]
                for d in self.defect_concentrations(
                    chemical_potentials=chemical_potentials,
                    temperature=quench_temperature,
                    fermi_level=high_temp_fermi_level,
                )
            ]
        )

        fdos = FermiDos(bulk_dos, bandgap=self.band_gap)
        _, fdos_vbm = fdos.get_cbm_vbm()

        def _get_total_q(ef):
            qd_tot = fixed_defect_charge
            qd_tot += fdos.get_doping(fermi_level=ef + fdos_vbm, temperature=temperature)
            return qd_tot

        return bisect(_get_total_q, -1.0, self.band_gap + 1.0)

    def get_dopability_limits(self, chemical_potentials):
        """
        Find Dopability limits for a given chemical potential.
        This is defined by the defect formation energies which first cross zero
        in formation energies.
        This determine bounds on the fermi level.

        Does this by computing formation energy for every stable defect with non-zero charge.
        If the formation energy value changes sign on either side of the band gap, then
        compute the fermi level value where the formation energy is zero
        (formation energies are lines and basic algebra shows: x_crossing = x1 - (y1 / q)
        for fermi level, x1, producing formation energy y1)

        Args:
            chemical_potentials: dict of chemical potentials to use for calculation fermi level
        Returns:
             lower dopability limit, upper dopability limit
            (returns None if no limit exists for upper or lower i.e. no negative defect
            crossing before +/- 20 of band edges OR defect formation energies are entirely zero)
        """
        min_fl_range = -20.0
        max_fl_range = self.band_gap + 20.0

        lower_lim = None
        upper_lim = None
        for def_entry in self.all_stable_entries:
            min_fl_formen = def_entry.formation_energy(
                chemical_potentials=chemical_potentials, fermi_level=min_fl_range
            )
            max_fl_formen = def_entry.formation_energy(
                chemical_potentials=chemical_potentials, fermi_level=max_fl_range
            )

            if min_fl_formen < 0.0 and max_fl_formen < 0.0:
                logger.error(
                    "Formation energy is negative through entire gap for entry {} q={}."
                    " Cannot return dopability limits.".format(def_entry.name, def_entry.charge)
                )
                return None, None
            if np.sign(min_fl_formen) != np.sign(max_fl_formen):
                x_crossing = min_fl_range - (min_fl_formen / def_entry.charge)
                if min_fl_formen < 0.0:
                    if lower_lim is None or lower_lim < x_crossing:
                        lower_lim = x_crossing
                else:
                    if upper_lim is None or upper_lim > x_crossing:
                        upper_lim = x_crossing

        return lower_lim, upper_lim

    def plot(
        self,
        mu_elts=None,
        xlim=None,
        ylim=None,
        ax_fontsize=1.3,
        lg_fontsize=1.0,
        lg_position=None,
        fermi_level=None,
        title=None,
        saved=False,
    ):
        """
        Produce defect Formation energy vs Fermi energy plot
        Args:
            mu_elts:
                a dictionnary of {Element:value} giving the chemical
                potential of each element
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis
            ylim:
                Tuple (min,max) giving the range for the formation energy axis
            ax_fontsize:
                float  multiplier to change axis label fontsize
            lg_fontsize:
                float  multiplier to change legend label fontsize
            lg_position:
                Tuple (horizontal-position, vertical-position) giving the position
                to place the legend.
                Example: (0.5,-0.75) will likely put it below the x-axis.
            saved:


        Returns:
            a matplotlib object

        """
        if xlim is None:
            xlim = (-0.5, self.band_gap + 0.5)
        xy = {}
        lower_cap = -100.0
        upper_cap = 100.0
        y_range_vals = []  # for finding max/min values on y-axis based on x-limits
        for defnom, def_tl in self.transition_level_map.items():
            xy[defnom] = [[], []]
            if def_tl:
                org_x = list(def_tl.keys())  # list of transition levels
                org_x.sort()  # sorted with lowest first

                # establish lower x-bound
                first_charge = max(def_tl[org_x[0]])
                for chg_ent in self.stable_entries[defnom]:
                    if chg_ent.charge == first_charge:
                        form_en = chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=lower_cap)
                        fe_left = chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=xlim[0])

                xy[defnom][0].append(lower_cap)
                xy[defnom][1].append(form_en)
                y_range_vals.append(fe_left)

                # iterate over stable charge state transitions
                for fl in org_x:
                    charge = max(def_tl[fl])
                    for chg_ent in self.stable_entries[defnom]:
                        if chg_ent.charge == charge:
                            form_en = chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=fl)
                    xy[defnom][0].append(fl)
                    xy[defnom][1].append(form_en)
                    y_range_vals.append(form_en)

                # establish upper x-bound
                last_charge = min(def_tl[org_x[-1]])
                for chg_ent in self.stable_entries[defnom]:
                    if chg_ent.charge == last_charge:
                        form_en = chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=upper_cap)
                        fe_right = chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=xlim[1])
                xy[defnom][0].append(upper_cap)
                xy[defnom][1].append(form_en)
                y_range_vals.append(fe_right)
            else:
                # no transition - just one stable charge
                chg_ent = self.stable_entries[defnom][0]
                for x_extrem in [lower_cap, upper_cap]:
                    xy[defnom][0].append(x_extrem)
                    xy[defnom][1].append(chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=x_extrem))
                for x_window in xlim:
                    y_range_vals.append(chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=x_window))

        if ylim is None:
            window = max(y_range_vals) - min(y_range_vals)
            spacer = 0.1 * window
            ylim = (min(y_range_vals) - spacer, max(y_range_vals) + spacer)

        if len(xy) <= 8:
            colors = cm.Dark2(np.linspace(0, 1, len(xy)))  # pylint: disable=E1101
        else:
            colors = cm.gist_rainbow(np.linspace(0, 1, len(xy)))  # pylint: disable=E1101

        plt.figure()
        plt.clf()
        width = 12
        # plot formation energy lines
        for_legend = []
        for cnt, defnom in enumerate(xy.keys()):
            plt.plot(xy[defnom][0], xy[defnom][1], linewidth=3, color=colors[cnt])
            for_legend.append(self.stable_entries[defnom][0].copy())

        # plot transtition levels
        for cnt, defnom in enumerate(xy.keys()):
            x_trans, y_trans = [], []
            for x_val, chargeset in self.transition_level_map[defnom].items():
                x_trans.append(x_val)
                for chg_ent in self.stable_entries[defnom]:
                    if chg_ent.charge == chargeset[0]:
                        form_en = chg_ent.formation_energy(chemical_potentials=mu_elts, fermi_level=x_val)
                y_trans.append(form_en)
            if len(x_trans):
                plt.plot(
                    x_trans,
                    y_trans,
                    marker="*",
                    color=colors[cnt],
                    markersize=12,
                    fillstyle="full",
                )

        # get latex-like legend titles
        legends_txt = []
        for dfct in for_legend:
            flds = dfct.name.split("_")
            if flds[0] == "Vac":
                base = "$Vac"
                sub_str = "_{" + flds[1] + "}$"
            elif flds[0] == "Sub":
                flds = dfct.name.split("_")
                base = "$" + flds[1]
                sub_str = "_{" + flds[3] + "}$"
            elif flds[0] == "Int":
                base = "$" + flds[1]
                sub_str = "_{inter}$"
            else:
                base = dfct.name
                sub_str = ''

            legends_txt.append(base + sub_str)

        if not lg_position:
            plt.legend(legends_txt, fontsize=lg_fontsize * width, loc=0)
        else:
            plt.legend(
                legends_txt,
                fontsize=lg_fontsize * width,
                ncol=3,
                loc="lower center",
                bbox_to_anchor=lg_position,
            )

        plt.ylim(ylim)
        plt.xlim(xlim)

        plt.plot([xlim[0], xlim[1]], [0, 0], "k-")  # black dashed line for Eformation = 0
        plt.axvline(x=0.0, linestyle="--", color="k", linewidth=3)  # black dashed lines for gap edges
        plt.axvline(x=self.band_gap, linestyle="--", color="k", linewidth=3)

        if fermi_level is not None:
            plt.axvline(x=fermi_level, linestyle="-.", color="k", linewidth=2)  # smaller dashed lines for gap edges

        plt.xlabel("Fermi energy (eV)", size=ax_fontsize * width)
        plt.ylabel("Defect Formation\nEnergy (eV)", size=ax_fontsize * width)
        if title:
            plt.title("{}".format(title), size=ax_fontsize * width)

        if saved:
            plt.savefig(str(title) + "FreyplnravgPlot.pdf")
        else:
            return plt


class DefectPredominanceDiagram(MSONable):
    """
    Module for constructing defect predominance diagrams (DPDs). DPDs show the concentration of each
    defect type at a given chemical potential. Traditionally, this is for binary oxides, with the
    independent variable being the oxygen partial pressure. Defect predominance diagrams are related
    to Brouwer Diagrams, but where Brouwer diagrams are constructed by assuming a single dominant
    point defect in each region of interest, DPDs make no such assumption.
    """

    def __init__(self, defect_phase_diagram: DefectPhaseDiagram,
                 bulk_dos: CompleteDos, entries: List[ComputedStructureEntry]):
        """
        Initialize the DPD.

        Args:
            defect_phase_diagram (DefectPhaseDiagram): A defect phase diagram for the system of interest.
            bulk_dos (Dos): a DOS object for the bulk structure. If bulk dos has a different band gap
                than the band gap stored in defect_phase_diagram, then the bands will be scissored to
                reproduce the band gap in defect_phase_diagram. This allows one to use, for example, a GGA
                bulk dos with high kpoint density, but have the band gap come from hybrid calculations.
            entries ([ComputedStructureEntry]): Computed structure entries for the bulk compounds. These are
                used to construct the phase diagram at temperatures of interest and establish the bounds on
                the chemical potentials. Example: for defects in MgO, one should have a computed structure entry
                for MgO, O2, and Mg.

        """
        self.defect_phase_diagram = defect_phase_diagram
        self.bulk = self.defect_phase_diagram.entries[0].bulk_structure
        self.bulk_dos = bulk_dos
        self.entries = entries
        self.pd = PhaseDiagram(entries)
        self.chempots = {el: ent.energy_per_atom for el, ent in self.pd.el_refs.items()}
        self.limits = {e: self.pd.get_transition_chempots(e) for e in self.pd.elements}
        self.els = {e.defect.defect_composition for e in self.defect_phase_diagram.entries}

    def _get_pd_and_cp(self, temperature: Union[int, float]):
        """
        Helper function to get the finite temperature phase diagram and chemical potentials at
        a given temperature.
        """
        chempots = {}
        # Manually adjust the chempots of the elemental references
        for el in self.chempots:
            g_interp = interp1d([int(t) for t in G_ELEMS.keys()], [G_ELEMS[t][str(el)] for t in G_ELEMS])
            chempots[el] = self.chempots[el] + g_interp(temperature)

        # Now update the compounds (not affected by gf_sisso())
        cpt_entrs = copy.deepcopy(self.entries)
        for e in cpt_entrs:
            for el, amt in e.composition.items():
                e.energy_adjustments.append(
                    CompositionEnergyAdjustment(-chempots[el], n_atoms=amt, name="{} to 0eV adjustment".format(el))
                )

        gentries = GibbsComputedStructureEntry.from_entries(cpt_entrs, temp=temperature)

        return PhaseDiagram(gentries), chempots

    def solve_along_chempot_line(self, temperature: Union[int, float], element: Element):
        """
        Solve for the defect concentration with respect to the chemical potential of one element, e.g. oxygen,
        going from the minimum to maximum value of that element's chempot as determined from the phase
        diagram.

        Args:
            temperature (int or float): temperature of interest
            element (Element): the element who's chemical potential will be varied

        Returns:
            (mus, results) the chemical potential linspace used for the evaluation along with
            the an array of results (see solve())
        """
        if temperature >= 300:
            pd, chempots = self._get_pd_and_cp(temperature)
        elif temperature > 0:
            warnings.warn("Requested temperature is < 300K. Interpolation data does"
                          "not exist for low temps, and so 0K energies will be used")
            pd, chempots = self.pd, self.chempots
        else:
            raise ValueError("Illegal temperature. Only T > 0K are permitted.")

        bulk_energy = pd.get_hull_energy(self.bulk.composition.reduced_composition)

        def foo(x):
            mycopy = chempots.copy()
            mycopy.update({element: x})
            mycopy.update(
                {
                    c: bulk_energy - sum([v for k, v in mycopy.items() if k != c])
                    for c in mycopy if c != element
                }
            )
            return mycopy

        _min, _max = sorted(pd.get_transition_chempots(element))
        mus = np.linspace(_min, _max, 50)

        _results = [self.solve(temperature, foo(cp)) for cp in mus]

        results: Dict = {d['name']: [] for d in _results[0]}
        for r in _results:
            for d in r:
                results[d['name']].append(d['conc'])

        return mus, results

    def solve(self, temperature: Union[int, float], chempots: Dict):
        """
        Get the concentration of defects at a single point by solving for the fermi
        level self consistently.

        Args:
            temperature (int or float): temperature of interest
            chempots (dict): dictionary of chemical potentials for the various elements in
                the defect phase diagram.
        """
        ef_scf = self.defect_phase_diagram.solve_for_fermi_energy(temperature, chempots, self.dos)
        fd = FermiDos(self.dos, bandgap=self.defect_phase_diagram.band_gap)

        cb_integral = np.sum(
            fd.tdos[fd.idx_cbm:]
            * f0(fd.energies[fd.idx_cbm:], ef_scf+fd.get_cbm_vbm()[1], temperature)
            * fd.de[fd.idx_cbm:], axis=0)
        vb_integral = np.sum(
            fd.tdos[:fd.idx_vbm + 1]
            * f0(-fd.energies[:fd.idx_vbm + 1], -ef_scf-fd.get_cbm_vbm()[1], temperature)
            * fd.de[:fd.idx_vbm + 1], axis=0)

        p = vb_integral / (fd.volume * fd.A_to_cm ** 3)
        n = cb_integral / (fd.volume * fd.A_to_cm ** 3)

        conc = self.defect_phase_diagram.defect_concentrations(chempots, temperature, fermi_level=ef_scf)
        conc.extend([{'name': 'n', 'conc': n, 'charge': -1}, {'name': 'p', 'conc': p, 'charge': 1}])
        return conc

    def plot(self, temperature: Union[int, float],
             element: Element, partial_pressure: bool = False,
             concentration_threshold: Union[int, float] = 0):
        """
        Plot the defect predomenance diagram along the chempot line of the specified element.

        Args:
            temperature (int or float): temperature of interest
            element (Element): vary the chemical potential of this element
            partial_pressure (bool): whether or not to plot the abscissa in units of log(partial pressure)
                as opposed to chemical potential. Default=False.
            concentration_threshold (int or float): do not plot defects that only have concentrations
                below this value. Keeps the plot looking cleaner by ignoring very unstable defects.

        Returns:
            Matplotlib Axes
        """
        fig, ax1 = plt.subplots(figsize=(12, 8))

        ax1.minorticks_on()
        ax1.set_ylabel("log(C) $\\frac{1}{m^3}$", size=30)
        ax1.set_xlabel("$log \\left( p_{} \\right)$".format(element.symbol), size=30) if \
            partial_pressure else ax1.set_xlabel("$\\mu_{}$".format(element.symbol), size=30)
        ax1.tick_params(which='major', length=8, width=2, direction='in', top=True, right=True, labelsize=20)
        ax1.tick_params(which='minor', length=2, width=2, direction='in', top=True, right=True, labelsize=20)

        mx = 0
        mn = 0
        mus, results = self.solve_along_chempot_line(temperature=temperature, element=element)
        px2 = ((mus[:-1]) / (boltzmann_eV * temperature))*np.log(2.71828)

        for r in results:
            if all([_ < concentration_threshold for _ in results[r]]):
                continue
            cs = np.log10(np.multiply(1e6, results[r]))
            mx, mn = max((mx, max(cs))), min((mn, min(cs)))
            free_carrier = (r == 'n' or r == 'p')
            ax1.plot(
                px2 if partial_pressure else mus,
                cs[:-1] if partial_pressure else cs, '--' if free_carrier else '-',
                linewidth=3 if free_carrier else 5,
                alpha=.8, label=r
            )

        left = max(px2) if partial_pressure else max(mus)
        ax1.axvspan(left, left + .1, alpha=0.2, color='red')
        ax1.axvspan(min(px2), min(px2) - .1, alpha=0.2, color='red') if \
            partial_pressure else ax1.axvspan(min(mus), min(mus) - .1, alpha=0.2, color='red')

        ax1.set_ylim(bottom=0, top=mx + 5)

        plt.legend(loc='best', fontsize=12)
        return ax1
