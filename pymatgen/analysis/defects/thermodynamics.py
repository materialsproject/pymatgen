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

    def __init__(self, entries, vbm, band_gap):
        self.entries = entries
        self.vbm = vbm
        self.band_gap = band_gap

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

        def similar_defect(a, b):
            """
            Used to filter out similar defects of different charges which
            are defined by the same type and location
            """
            return a.name == b.name and a.site == b.site

        # Limits for search
        # E_fermi = { -1 eV to band gap+1}
        # the 1 eV padding provides
        # E_formation. = { -1 eV to 20 eV}
        limits = [[-1, self.band_gap + 1], [-1, 20]]

        stable_entries = {}
        stable_charges = {}
        transition_level_map = {}

        # Grouping by defect types
        for _, defects in groupby(sorted(self.entries, key=similar_defect), similar_defect):
            defects = list(defects)

            # prepping coefficient matrix for half-space intersection
            hyperplanes = np.array([
                np.array([-1.0 * entry.charge, 1, -1.0 * (entry.energy + entry.charge * self.vbm)])
                for entry in defects
            ])

            border_hyperplanes = [[-1, 0, limits[0][0]],
                                  [1, 0, -limits[0][1]],
                                  [0, -1, limits[1][0]],
                                  [0, 1, -limits[1][1]]]
            hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])
            interior_point = np.array([self.band_gap / 2, 0])

            hs_ints = HalfspaceIntersection(halfspaces, np.array([0.0, 0.0]))

            # Group the intersections and coresponding facets
            ints_and_facets = zip(hs_int.intersections, hs_int.dual_facets)
            # Only inlcude the facets corresponding to entries, not the boundaries
            ints_and_facets = filter(ints_and_facets, lambda _, facet: all(facet < len(self.entries)))
            # sort based on transition level
            ints_and_facets = sorted(ints_and_facets, key=lambda intersection, _: intersection[0])
            # Unpack into lists
            intersections, facets = zip(*ints_and_facets)

            # Map of transition level: charge states
            transition_levels_map[defects[0].name] = {
                intersection[0]: [defects[i].charge for i in facet]
                for intersection, facet in ints_and_facets
            }

            stable_entries[defects[0].name] = [defects[i] for dual in facets for i in dual]

            finished_charges[defects[0].name] = [defect.charge for defect in defects]

        self.transition_levels_map = transition_levels_map
        self.transition_levels = {
            defect_name: list(defect_tls.keys())
            for defect_name, defect_tls in transition_levels_map.items()
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
                'conc': dfct.defect_concentration(chemical_potentials=chemical_potentials,
                                                  temperature=temperature,
                                                  fermi_level=fermi_level),
                'name': dfct.name,
                'charge': dfct.charge
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
            if min_tl < tolerance:
                min_tl = min(self.transition_levels_map.keys())
                max_charge = max(self.transition_levels_map[min_tl])
                test_charges = [charge for charge in test_charges if charge < max_charge]

            # More negative charges will shift the maximum transition level up
            # Minimum charge is limited by this if transition level is near CBM
            if max_tl > (self.band_gap - tolerance):
                max_tl = max(self.transition_levels_map.keys())
                min_charge = min(self.transition_levels_map[max_tl])
                test_charges = [charge for charge in test_charges if charge > min_charge]

            recommendations[def_type] = test_charges

        return recommendations


class GrandCanonicalDefectPhaseDiagram(DefectPhaseDiagram):
    """
    Adapts the DefectPhaseDiagram class to a grand canonical approach
    This allows for self consistently determining fermi levels,
    defect concentrations, and free carrier concentrations as a function
    of temperature.

    Args:
        mu_elts: {Element: float} is dictionary of chemical potentials, determined from a phase diagram
        T: Temperature (Kelvin)
        kwargs: Arguments to DefectsAnalyzer as keyword pair
    """

    def __init__(self, mu_elts, T=298.15, **kwargs):
        super(self.__class__, self).__init__(**kwargs)
        self.mu_elts = mu_elts
        self.T = T
        self._ef = None
        self._all_possible_mu_elts = None

    @staticmethod
    def generate_gcdpd_from_pda(dentries, pda, T=298.15):
        """
        A static method for instantiating a GrandCanonicalDefectPhaseDiagram from a DefectEntry list
        and a PhaseDiagramAnalyzer object from MaterialsProject
        """
        # use bulk_entry object from first dentry to find the composition
        #    that one requires chemical potentials from...
        CPA = ChemPotAnalyzer(bulk_ce=dentries[0].bulk)
        all_possible_mu_elts = CPA.get_chempots_from_pda(pda)
        print('Generated following possible regions for defining ' 'chemical potentials: ', all_possible_mu_elts.keys())
        mu_region, mu_elts = all_possible_mu_elts.items()[0]
        print('Proceeding with chemical potentials from ', mu_region, '\n', mu_elts)

        gcdpd = GrandCanonicalDefectPhaseDiagram(mu_elts, T=T, dentries=dentries)
        gcdpd.all_possible_mu_elts = all_possible_mu_elts

        return gcdpd

    @staticmethod
    def generate_gcdpd_from_MP(dentries, mpid=None, composition=None, T=298.15, mapi_key=None):
        """
        A static method for instantiating a GrandCanonicalDefectPhaseDiagram
        from either a mpid string or a pymatgen Composition object
            [If both are provided, defaults to use of mpid]
            [If neither are provided, uses bulk_entry of first DefectEntry to generate the composition objectt]

        Uses the MaterialsProject database to create a phase diagram
        """
        all_species = set([elt for dfct in dentries for elt in dfct._structure.composition.elements])
        bulk_species = set(dentries[0].bulk.composition.elements)
        sub_species = all_species - bulk_species

        if mpid:
            MPcpa = MPChemPotAnalyzer(sub_species=sub_species, mpid=mpid, mapi_key=mapi_key)
            all_possible_mu_elts = MPcpa.analyze_GGA_chempots()
        else:
            if not composition:
                composition = dentries[0].bulk.composition
            MPcpa = MPChemPotAnalyzer(sub_species=sub_species, mapi_key=mapi_key)
            all_possible_mu_elts = MPcpa.get_chempots_from_composition(composition)

        print('Generated following possible regions for defining ' 'chemical potentials: ', all_possible_mu_elts.keys())
        mu_region, mu_elts = all_possible_mu_elts.items()[0]
        print('Proceeding with chemical potentials from ', mu_region, '\n', mu_elts)

        gcdpd = GrandCanonicalDefectPhaseDiagram(mu_elts, T=T, dentries=dentries)
        gcdpd.all_possible_mu_elts = all_possible_mu_elts

        return gcdpd

    def set_T(self, T):
        self.T = T

    @property
    def fermi_energy(self):
        return self._ef

    def _get_total_q(self, ef, gap, bulk_dos):
        qd_tot = 0
        for d in self.list_defect_concentrations(self.mu_elts, temp=self.T, ef=ef):
            qd_tot += d['charge'] * d['conc']

        q_h_cont = IntrinsicCarrier(bulk_dos, exp_gap=gap)
        ef_ref_cbm = ef - gap
        nden = q_h_cont.get_n_density(ef_ref_cbm, self.T, ref='CBM', unitcell=False)
        pden = q_h_cont.get_p_density(ef, self.T, unitcell=False)
        qd_tot += pden - nden
        return qd_tot

    def solve_for_fermi_energy(self, bulk_dos):
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
        from scipy.optimize import bisect
        self._ef = bisect(lambda e: self._get_total_q(e, self.band_gap, bulk_dos), -1., self.band_gap + 1.)

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
