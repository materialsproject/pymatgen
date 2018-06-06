#!/usr/bin/env python

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Shyam Dwaraknath"
__email__ = "dbroberg@berkeley.edu, shyamd@lbl.gov"
__status__ = "Development"
__date__ = "January 11, 2018"

from math import exp
import numpy as np
from monty.json import MSONable
from scipy import integrate
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import bisect
from itertools import groupby, chain

from pymatgen.analysis.defects.utils import kb


class DefectPhaseDiagram(MSONable):
    """
    This is similar to a PhaseDiagram object in pymatgen, but has ability to do quick analysis of defect formation energies
    when fed DefectEntry objects

    uses many of the capabilities from PyCDT's DefectsAnalyzer class...

    This class is able to get:
        a) stability of charge states for a given defect,
        b) list of all formation ens
        c) transition levels in the gap
        d)

    Args:
        dentries ([DefectEntry]): A list of DefectEntry objects
    """

    def __init__(self, entries, vbm, band_gap, filter_compatible=True):
        self.vbm = vbm
        self.band_gap = band_gap
        self.filter_compatible = filter_compatible

        if filter_compatible:
            self.entries = [e for e in entries if e.parameters.get("is_compatible", True)]
        else:
            self.entries = entries

        self.find_stable_charges()

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

        def similar_defect(a):
            """
            Used to filter out similar defects of different charges which
            are defined by the same type and location
            """
            return (a.name, a.site)

        # Limits for search
        # E_fermi = { -1 eV to band gap+1}
        # the 1 eV padding provides
        # E_formation. = { -21 eV to 20 eV}
        limits = [[-1, self.band_gap + 1], [-21, 20]]

        stable_entries = {}
        finished_charges = {}
        transition_level_map = {}

        # Grouping by defect types
        for _, defects in groupby(sorted(self.entries, key=similar_defect), similar_defect):
            defects = list(defects)
            for d in defects:
                print(d.name,d.charge, d.energy, self.vbm, d.energy + d.charge*self.vbm )

            # prepping coefficient matrix forx half-space intersection
            # [-Q, 1, -1*(E_form+Q*VBM)] -> -Q*E_fermi+E+-1*(E_form+Q*VBM) <= 0  where E_fermi and E are the variables in the hyperplanes
            hyperplanes = np.array(
                [[-1.0 * entry.charge, 1, -1.0 * (entry.energy + entry.charge * self.vbm)] for entry in defects])

            border_hyperplanes = [[-1, 0, limits[0][0]], [1, 0, -1 * limits[0][1]], [0, -1, limits[1][0]],
                                  [0, 1, -1 * limits[1][1]]]
            hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])

            interior_point = [self.band_gap / 2, -20]

            hs_ints = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

            # Group the intersections and coresponding facets
            ints_and_facets = zip(hs_ints.intersections, hs_ints.dual_facets)
            # Only inlcude the facets corresponding to entries, not the boundaries
            total_entries = len(defects)
            ints_and_facets = filter(lambda int_and_facet: all(np.array(int_and_facet[1]) < total_entries),
                                     ints_and_facets)
            # sort based on transition level
            ints_and_facets = list(sorted(ints_and_facets, key=lambda int_and_facet: int_and_facet[0][0]))
            # Unpack into lists
            intersections, facets = zip(*ints_and_facets)
            # Map of transition level: charge states

            transition_level_map[defects[0].name] = {
                intersection[0]: [defects[i].charge for i in facet]
                for intersection, facet in ints_and_facets
            }

            stable_entries[defects[0].name] = list(set([defects[i] for dual in facets for i in dual]))

            finished_charges[defects[0].name] = [defect.charge for defect in defects]

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

    def OLDfind_stable_charges(self):
        """
        Set the DefectPhaseDiagram attributes (which dont depend on chemical potential)

        THIS is a wrapper since find_stable_charges does not currently exist...will switch over when able
        """
        def similar_defect(a, b):
            """
            Used to filter out similar defects of different charges which
            are defined by the same type and location
            """
            return a.name == b.name and a.site == b.site

        self.stable_charges = {} #keys are defect names, items are list of charge states that are stable
        self.finished_charges = {} #keys are defect names, items are list of charge states that are included in the phase diagram
        self.transition_levels = {} #keys are defect names, items are list of [fermi level for transition, previous q, next q] sets
        xlim = (-0.1, self.band_gap+.1)
        nb_steps = 200 #SUPER bare bones because it is so slow...
        x = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nb_steps)
        list_elts = [elt  for dfct in self.entries  for elt in dfct.defect.generate_defect_structure().composition.elements]
        no_chempots = {elt: 0. for elt in set(list_elts)} #zerod chemical potentials for calculating stable defects

        #do a dumb grouping
        sim_defects_group = []
        for defectentry in self.entries:

            simdef_index = None
            for test_index, simdefset in enumerate(sim_defects_group):
                if (simdefset[0].name == defectentry.name) and (simdefset[0].site == defectentry.site):
                    simdef_index = test_index

            if simdef_index != None:
                sim_defects_group[simdef_index].append( defectentry)
            else:
                sim_defects_group.append( [defectentry])

        # Grouping by defect types
        for defects in sim_defects_group:
            defects = list(defects)
            t = defects[0].name
            print('Analyzing {}'.format(t))
            trans_level = []
            chg_type = []
            prev_min_q, cur_min_q = None, None
            for x_step in x:
                miny = 10000
                for dfct in defects:
                    val = dfct.formation_energy(chemical_potentials = no_chempots, fermi_level=x_step)
                    if val < miny:
                        miny = val
                        cur_min_q = dfct.charge

                if prev_min_q is not None:
                    if cur_min_q != prev_min_q:
                        trans_level.append((x_step, prev_min_q, cur_min_q))
                    if cur_min_q not in chg_type:
                        chg_type.append(cur_min_q)
                prev_min_q = cur_min_q

            self.stable_charges[t] = chg_type[:]
            self.finished_charges[t] = [e.charge for e in self.entries if e.name == t]
            self.transition_levels[t] = trans_level[:]

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

            # TODO: Should be this entry: concentration?
        return concentrations

    def suggest_charges(self, tolerance=0.1):
        """
        Suggest possible charges for defects to computee based on proximity
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
                for d in self.list_defect_concentrations(
                    chemical_potentials=chemical_potentials, temperature=temperature, fermi_level=ef)
            ])
            qd_tot += fdos.get_doping(fermi=ef, T=temperature)
            return qd_tot

        return bisect(_get_total_q, -1., self.band_gap + 1.)

    def solve_non_eq_fermilevel(self, lowT, highT, bulk_dos, show_carriers=True, hold_htemp_ef=True):
        """
        Solve for the Fermi energy at a low temperature value, but with defect concentrations frozen in
            at some higher temperature value This means that only the fermi level and free charge carriers concentrations
            will be solved self consistently (simulating a 'frozen-in' defect concentration approach)

        Args:
            lowT: the temperature at which you want to solve for fermi level and free carrier concentrations
            highT: the fixed temperature that sets the defect concentrations at a higher temperature
            bulk_dos: bulk system dos (pymatgen Dos object)
            show_carriers: Dictates the outputs of this function
                Set to True if you want the new carrier concentrations in addition to the fermi levels
            hold_htemp_ef: Specifies whether the fermi level will be fixed for the defect concentrations at high T
                Set to False if you want the fermi level of the high temperature defects to be the same as the fermi level
                for the low temperature charge carriers (see below note on implementation for explanation of this flag)
        Returns:
            Fermi energy

        A NOTE ON IMPLEMENTATION: while testing this function, it became apparent that negative defect formation energies
            can sometimes cause stochastic results for the fermi level as a result of the defect concentrations not
            being allowed to change within the self consistent result. To circumvent this, this function has a flag ('hold_htemp_ef')
            to allow for the fermi level of high temperature defects to vary within the self consistent approach.
            While this is no longer the fully physical 'frozen-in' approach, it helps remove the stocastic results
            that result from numerical issues

        """
        from scipy.optimize import bisect
        q_h_cont = IntrinsicCarrier(bulk_dos, exp_gap=self.band_gap)

        if hold_htemp_ef:
            self.solve_for_fermi_energy(bulk_dos, gap=self.band_gap)
            defectlist = self.get_defects_concentration(highT, self._ef, unitcell=False)

        def noneq_total_q(ef):
            qd_tot = 0
            if hold_htemp_ef:
                for d in defectlist:
                    qd_tot += d['charge'] * d['conc']
            else:
                for d in self.get_defects_concentration(highT, ef, unitcell=False):
                    qd_tot += d['charge'] * d['conc']
            ef_ref_cbm = ef - self.band_gap
            nden = q_h_cont.get_n_density(ef_ref_cbm, lowT, ref='CBM', unitcell=False)
            pden = q_h_cont.get_p_density(ef, lowT, unitcell=False)
            qd_tot += pden - nden
            return qd_tot

        non_eq_ef = bisect(lambda e: noneq_total_q(e), -1., self.band_gap + 1.)
        if show_carriers:
            ef_ref_cbm = non_eq_ef - self.band_gap
            nden = q_h_cont.get_n_density(ef_ref_cbm, lowT, ref='CBM', unitcell=False)
            pden = q_h_cont.get_p_density(non_eq_ef, lowT, unitcell=False)
            return non_eq_ef, nden, pden
        else:
            return non_eq_ef


