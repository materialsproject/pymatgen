# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import numpy as np
from monty.json import MSONable
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import bisect
from itertools import chain

from pymatgen.electronic_structure.dos import FermiDos
from pymatgen.analysis.defects.core import DefectEntry
from pymatgen.analysis.structure_matcher import PointDefectComparator

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"


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
        d)

    Args:
        dentries ([DefectEntry]): A list of DefectEntry objects
    """

    def __init__(self, entries, vbm, band_gap, filter_compatible=True, metadata={}):
        self.vbm = vbm
        self.band_gap = band_gap
        self.filter_compatible = filter_compatible

        if filter_compatible:
            self.entries = [e for e in entries if e.parameters.get("is_compatible", True)]
        else:
            self.entries = entries

        self.metadata = metadata
        self.find_stable_charges()

    def as_dict(self):
        """
        Json-serializable dict representation of DefectPhaseDiagram
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "entries": [entry.as_dict() for entry in self.entries],
             "vbm": self.vbm,
             "band_gap": self.band_gap,
             "filter_compatible": self.filter_compatible,
             "metadata": self.metadata}
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
        if 'entry_id' in entry_dict.keys() and 'entry_id' not in metadata:
            metadata['entry_id'] = entry_dict['entry_id']

        return cls(entries, vbm, band_gap, filter_compatible=filter_compatible,
                   metadata=metadata)

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

        def similar_defects( entryset):
            """
            Used for grouping similar defects of different charges
            Can distinguish identical defects even if they are not in same position
            """
            pdc = PointDefectComparator( check_charge=False, check_primitive_cell=True,
                                         check_lattice_scale=False)
            grp_def_sets = []
            grp_def_indices = []
            for ent_ind, ent in enumerate( entryset):
                #TODO: more pythonic way of grouping entry sets with PointDefectComparator
                matched_ind = None
                for grp_ind, defgrp in enumerate(grp_def_sets):
                    if pdc.are_equal( ent.defect, defgrp[0].defect):
                        matched_ind = grp_ind
                        break
                if matched_ind is not None:
                    grp_def_sets[matched_ind].append( ent.copy())
                    grp_def_indices[matched_ind].append( ent_ind)
                else:
                    grp_def_sets.append( [ent.copy()])
                    grp_def_indices.append( [ent_ind])

            return zip(grp_def_sets, grp_def_indices)

        # Limits for search
        # E_fermi = { -1 eV to band gap+1}
        # E_formation = { (min(Eform) - 30) to (max(Eform) + 30)}
        all_eform = [one_def.formation_energy(fermi_level=self.band_gap/2.) for one_def in self.entries]
        min_y_lim = min(all_eform) - 30
        max_y_lim = max(all_eform) + 30
        limits = [[-1, self.band_gap + 1], [min_y_lim, max_y_lim]]

        stable_entries = {}
        finished_charges = {}
        transition_level_map = {}

        # Grouping by defect types
        for defects, index_list in similar_defects( self.entries):
            defects = list(defects)

            # prepping coefficient matrix for half-space intersection
            # [-Q, 1, -1*(E_form+Q*VBM)] -> -Q*E_fermi+E+-1*(E_form+Q*VBM) <= 0  where E_fermi and E are the variables in the hyperplanes
            hyperplanes = np.array(
                [[-1.0 * entry.charge, 1, -1.0 * (entry.energy + entry.charge * self.vbm)] for entry in defects])

            border_hyperplanes = [[-1, 0, limits[0][0]], [1, 0, -1 * limits[0][1]], [0, -1, limits[1][0]],
                                  [0, 1, -1 * limits[1][1]]]
            hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])

            interior_point = [self.band_gap / 2, min(all_eform) - 1.]

            hs_ints = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

            # Group the intersections and coresponding facets
            ints_and_facets = zip(hs_ints.intersections, hs_ints.dual_facets)
            # Only inlcude the facets corresponding to entries, not the boundaries
            total_entries = len(defects)
            ints_and_facets = filter(lambda int_and_facet: all(np.array(int_and_facet[1]) < total_entries),
                                     ints_and_facets)
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
                    intersection[0]: [defects[i].charge for i in facet]
                    for intersection, facet in ints_and_facets
                }

                stable_entries[track_name] = list(set([defects[i] for dual in facets for i in dual]))

                finished_charges[track_name] = [defect.charge for defect in defects]
            else:
                # if ints_and_facets is empty, then there is likely only one defect...
                if len(defects) != 1:
                    #confirm formation energies dominant for one defect over other identical defects
                    name_set = [one_def.name+'_chg'+str(one_def.charge) for one_def in defects]
                    vb_list = [one_def.formation_energy( fermi_level=limits[0][0]) for one_def in defects]
                    cb_list = [one_def.formation_energy( fermi_level=limits[0][1]) for one_def in defects]

                    vbm_def_index = vb_list.index( min(vb_list))
                    name_stable_below_vbm = name_set[vbm_def_index]
                    cbm_def_index = cb_list.index( min(cb_list))
                    name_stable_above_cbm = name_set[cbm_def_index]

                    if name_stable_below_vbm != name_stable_above_cbm:
                        raise ValueError("HalfSpace identified only one stable charge out of list: {}\n"
                                         "But {} is stable below vbm and {} is "
                                         "stable above cbm.\nList of VBM formation energies: {}\n"
                                         "List of CBM formation energies: {}"
                                         "".format(name_set, name_stable_below_vbm, name_stable_above_cbm,
                                                   vb_list, cb_list))
                    else:
                        print("{} is only stable defect out of {}".format( name_stable_below_vbm, name_set))
                        transition_level_map[track_name] = {}
                        stable_entries[track_name] = list([defects[vbm_def_index]])
                        finished_charges[track_name] = [one_def.charge for one_def in defects]
                else:
                    transition_level_map[track_name] = {}

                    stable_entries[track_name] = list([defects[0]])

                    finished_charges[track_name] = [defects[0].charge]

        self.transition_level_map = transition_level_map
        self.transition_levels = {
            defect_name: list(defect_tls.keys())
            for defect_name, defect_tls in transition_level_map.items()
        }
        self.stable_entries = stable_entries
        self.finished_charges = finished_charges
        self.stable_charges = {
            defect_name: [entry.charge for entry in entries]
            for defect_name, entries in stable_entries.items()
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

    def defect_concentrations(self, chemical_potentials, temperature=300, fermi_level=0.):
        """
        Give list of all concentrations at specified efermi in the DefectPhaseDiagram
        args:
            chemical_potentials = {Element: number} is dictionary of chemical potentials to provide formation energies for
            temperature = temperature to produce concentrations from
            fermi_level: (float) is fermi level relative to valence band maximum
                Default efermi = 0 = VBM energy
        returns:
            list of dictionaries of defect concentrations
        """
        concentrations = []
        for dfct in self.all_stable_entries:
            concentrations.append({
                'conc':
                dfct.defect_concentration(
                    chemical_potentials=chemical_potentials, temperature=temperature, fermi_level=fermi_level),
                'name':
                dfct.name,
                'charge':
                dfct.charge
            })

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
                np.max(self.stable_charges[def_type]) + 2)
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

    def solve_for_fermi_energy(self, temperature, chemical_potentials, bulk_dos):
        """
        Solve for the Fermi energy self-consistently as a function of T
        and p_O2
        Observations are Defect concentrations, electron and hole conc
        Args:
            bulk_dos: bulk system dos (pymatgen Dos object)
            gap: Can be used to specify experimental gap.
                Will be useful if the self consistent Fermi level
                is > DFT gap
        Returns:
            Fermi energy
        """

        fdos = FermiDos(bulk_dos, bandgap=self.band_gap)

        def _get_total_q(ef):

            qd_tot = sum([
                d['charge'] * d['conc']
                for d in self.defect_concentrations(
                    chemical_potentials=chemical_potentials, temperature=temperature, fermi_level=ef)
            ])
            qd_tot += fdos.get_doping(fermi=ef + self.vbm, T=temperature)
            return qd_tot

        return bisect(_get_total_q, -1., self.band_gap + 1.)
