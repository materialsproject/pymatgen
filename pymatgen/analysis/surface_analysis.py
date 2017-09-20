# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
TODO:
    -Only works for monatomic adsorbates
    -Only works for slabs of three species or less, e.g.
        clean ternary, binary+adsorbate, clean binary,
        elemental+adsorbate or clean elemental slab
    -Will address these issues and make the analyzer more
        generalized in the future
"""

from __future__ import division, unicode_literals

import copy

import numpy as np
from scipy.stats import linregress
from matplotlib import cm
import itertools
import warnings
import random

from pymatgen.core.structure import Structure, Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import Slab
from pymatgen.analysis.wulff import WulffShape
from pymatgen import MPRester
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen import Element
from pymatgen.util.plotting import pretty_plot
from pymatgen.analysis.reaction_calculator import ComputedReaction

EV_PER_ANG2_TO_JOULES_PER_M2 = 16.0217656

__author__ = "Richard Tran"
__copyright__ = "Copyright 2014, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "8/24/17"


class SurfaceEnergyCalculator(object):
    """
    A class used for analyzing the surface energies of a material of a given
        material_id. By default, this will use entries calculated from the
        Materials Project to obtain chemical potential and bulk energy. As a
        result, the difference in VASP parameters between the user's entry
        and the parameters used by Materials Project, may lead to a rough
        estimate of the surface energy. For best results, it is recommend that
        the user calculates all decomposition components first, and insert the
        results into their own database as a pymatgen-db entry and use those
        entries instead (custom_entries). In addition, this code will only use
        one bulk entry to calculate surface energy. Ideally, to get the most
        accurate surface energy, the user should compare their slab energy to
        the energy of the oriented unit cell with both calculations containing
        consistent k-points to avoid convergence problems as the slab size is
        varied. See:
            Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
                Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.
        and
            Rogal, J., & Reuter, K. (2007). Ab Initio Atomistic Thermodynamics for
                Surfaces : A Primer. Experiment, Modeling and Simulation of Gas-Surface
                Interactions for Reactive Flows in Hypersonic Flights, 2–1 – 2–18.

    .. attribute:: ref_element

        All chemical potentials can be written in terms of the range of chemical
            potential of this element which will be used to calculate surface energy.

    .. attribute:: mprester

        Materials project rester for querying entries from the materials project.
            Requires user MAPIKEY.

    .. attribute:: ucell_entry

        Materials Project entry of the material of the slab.

    .. attribute:: x

        Reduced amount composition of decomposed compound A in the bulk.

    .. attribute:: y

        Reduced amount composition of ref_element in the bulk.

    .. attribute:: gbulk

        Gibbs free energy of the bulk per formula unit

    .. attribute:: e_of_element

        Energy per atom of ground state ref_element, eg. if ref_element=O,
            than e_of_element=1/2*E_O2.

    .. attribute:: adsorbate

        Composition of the adsorbate (if there is one).

    """

    def __init__(self, ucell_entry, comp1, ref_el_comp,
                 exclude_ids=[], custom_entries=[],
                 mapi_key=None, adsorbate_entry=None):
        """
        Analyzes surface energies and Wulff shape of a particular
            material using the chemical potential.
        Args:
            ucell_entry (material_id or computed_entry): Materials Project or entry
                of the bulk system the slab is based on (a string, e.g., mp-1234).
            comp1 (Composition): Composition to be considered as dependent
                variable.
            ref_el_comp (Composition): Composition to be considered as independent
                variable. E.g., if you want to show the stability
                ranges of all Li-Co-O phases wrt to uLi
            exclude_ids (list of material_ids): List of material_ids
                to exclude when obtaining the decomposition components
                to calculate the chemical potential
            custom_entries (list of pymatgen-db type entries): List of
                user specified pymatgen-db type entries to use in finding
                decomposition components for the chemical potential
            mapi_key (str): Materials Project API key for accessing the
                MP database via MPRester
            adsorbate_entry (Composition): Computed entry of adsorbate,
                defaults to None. Could be an isolated H2 or O2 molecule for gaseous
                    adsorption or just the bulk ground state structure for metals.
                    Either way, the extract term is going to be in energy per atom.
        """

        self.comp1 = comp1
        self.ref_el_comp = ref_el_comp
        self.mprester = MPRester(mapi_key) if mapi_key else MPRester()
        self.ucell_entry = self.mprester.get_entry_by_material_id( \
            ucell_entry, inc_structure=True) \
            if type(ucell_entry).__name__ == "str" else ucell_entry
        ucell_comp = self.ucell_entry.composition

        entries = [entry for entry in
                   self.mprester.get_entries_in_chemsys(list( \
                       ucell_comp.reduced_composition.as_dict().keys()),
                       property_data=["e_above_hull",
                                      "material_id"])
                   if entry.data["e_above_hull"] == 0 and
                   entry.data["material_id"] not in exclude_ids] \
            if not custom_entries else custom_entries

        # Get x and y, the number of compositions in a formula unit of the bulk
        self.reactants = [entry for entry in entries if
                          entry.composition.reduced_composition in [comp1,
                                                                    self.ref_el_comp]]
        rxn = ComputedReaction(self.reactants, [self.ucell_entry])
        if len(rxn.reactants) == len(rxn.products):
            x = y = 1
        else:
            y = abs(rxn.coeffs[rxn.all_comp.index(self.ref_el_comp)])
            x = abs(rxn.coeffs[rxn.all_comp.index(comp1)])

        # Calculate Gibbs free energy of the bulk per unit formula
        gbulk = self.ucell_entry.energy / \
                ucell_comp.get_integer_formula_and_factor()[1]

        self.x = x
        self.y = y
        self.gbulk = gbulk
        self.adsorbate_entry = adsorbate_entry
        self.adsorbate_as_str = self.adsorbate_entry.composition.elements[0].name \
            if self.adsorbate_entry else None
        self.decomp_entries = entries

    def chempot_range(self, full_chempot=True):
        """
        Calculates the chemical potential range allowed for
        non-stoichiometric clean surface energy calculations
        Args:
            full_chempot (bool): Whether or not to calculate
                the range based on chemical potential from the
                referance element to the compound of the ucell_entry
                (False), or to consider all decompositions.
        """

        pd = PhaseDiagram(self.decomp_entries)
        pda = PDAnalyzer(pd)

        chempot_ranges = pda.get_chempot_range_map(self.ref_el_comp.elements)
        # If no chemical potential is found, we return u=0, eg.
        # for a elemental system, the relative u of Cu for Cu is 0
        if full_chempot:
            if not chempot_ranges:
                chempot_range = [[-1, -1], [0, 0]]
            else:
                all_u = []
                for entry in chempot_ranges.keys():
                    all_u.extend(chempot_ranges[entry][0]._coords)
                chempot_range = [min(all_u), max(all_u)]
        else:
            ucell_comp = self.ucell_entry.composition
            chempot_range = [chempot_ranges[entry] for entry in chempot_ranges.keys()
                             if entry.composition ==
                             ucell_comp][0][0]._coords if \
                chempot_ranges else [[-1, -1], [0, 0]]

        chempot_range = list(chempot_range)
        return sorted([chempot_range[0][0], chempot_range[1][0]])

    def surface_energy_coefficients(self, clean_slab_entry, ads_slab_entry=None):
        """
        Calculates the coefficients of the surface energy for a single slab
            in the form of gamma=b1x1+b2x2+b3 where x1 is the chemical
            potential of ref_el_comp and x2 is the chemical potential of
            the adsorbate.
        Args:
            clean_slab_entry (entry): An entry object for the clean slab
            ads_slab_entry (entry): An optional entry object for the
                adsorbed slab, defaults to None.

        Returns (list): List of the coefficients for surface energy [b1, b2, b3]
        """

        reduced_comp = self.ucell_entry.composition.reduced_composition.as_dict()
        # Get the composition in the slab
        slab = clean_slab_entry.structure
        # Calculate surface area
        m = slab.lattice.matrix
        Aclean = np.linalg.norm(np.cross(m[0], m[1]))
        # For the adsorbed surface area, just set to 1 if no entry given
        if ads_slab_entry:
            adslab = ads_slab_entry.structure
            m = adslab.lattice.matrix
            Aads = np.linalg.norm(np.cross(m[0], m[1]))
        else:
            Aads = 1


        comp = slab.composition.as_dict()

        if len(reduced_comp.keys()) == 1:
            Nx = Ny = comp[self.ucell_entry.structure[0].species_string]
        else:
            Nx = slab.composition.as_dict()[str(self.comp1.elements[0])]
            Ny = slab.composition.as_dict()[str(self.ref_el_comp.elements[0])]


        # get number of adsorbates in slab
        Nads = self.Nads_in_slab(ads_slab_entry) if ads_slab_entry else 0
        # get number of adsorbed surfaces, set to 1 for clean to avoid division by 0
        Nsurfs = 1 if not ads_slab_entry else self.Nsurfs_ads_in_slab(ads_slab_entry)

        # return the coefficients of the surface energy
        # b1 is the coefficient for chempot of reference atom
        b1 = (-1 / (2 * Aclean))*(Ny - (self.y / self.x) * Nx)
        #b2 is the coefficient for chempot of adsorbate
        b2 = (-1*Nads) / (Nsurfs * Aads)
        gibbs_binding = self.gibbs_binding_energy(ads_slab_entry,
                                                  clean_slab_entry) if ads_slab_entry else 0
        # b3 is the intercept
        b3 = (1 / (2 * Aclean)) * (clean_slab_entry.energy - (Nx/self.x)*self.gbulk - \
                              (Ny - (self.y / self.x) * Nx) * [entry.energy_per_atom \
                                                               for entry in self.reactants
                                                               if entry.composition.reduced_composition
                                                               == self.ref_el_comp][0]) + \
             gibbs_binding*(Nads/(Nsurfs*Aads))

        return [b1, b2, b3]

    def calculate_gamma_at_u(self, clean_slab_entry, ads_slab_entry=None, u_ref=0, u_ads=0):
        """
        Quickly calculates the surface energy for the slab of the vasprun file
            at a specific chemical potential for the reference and adsorbate.
            Both chemical potentials are set to 0 by default (i.e. clean
            stoichiometric slab).
        args:
            slab_entry (entry): Entry containing the final energy and structure
                of the slab whose surface energy we want to calculate
            ads_slab_entry (entry): An optional entry object for the
                adsorbed slab, defaults to None.
            u_ref (float): The chemical potential of the reference element
            u_ads (float): The chemical potential of the adsorbate

        Returns (float): surface energy
        """

        coeffs = self.surface_energy_coefficients(clean_slab_entry,
                                                  ads_slab_entry=ads_slab_entry)
        return np.dot([u_ref, u_ads, 1], coeffs)

    def get_monolayer(self, ads_slab_entry, clean_slab_entry):

        m = ads_slab_entry.structure.lattice.matrix
        A_ads = np.linalg.norm(np.cross(m[0], m[1]))
        m = clean_slab_entry.structure.lattice.matrix
        A_clean = np.linalg.norm(np.cross(m[0], m[1]))
        n = (A_ads / A_clean)

        return n

    def gibbs_binding_energy(self, ads_slab_entry, clean_slab_entry):
        """
        Calculates the adsorption energy or Gibb's
        binding energy of an adsorbate on a surface
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
            clean_slab_entry (entry): The entry of the clean slab
        """

        n = self.get_monolayer(ads_slab_entry, clean_slab_entry)
        Nads = self.Nads_in_slab(ads_slab_entry)

        return (ads_slab_entry.energy - n * clean_slab_entry.energy) / Nads \
               - self.adsorbate_entry.energy_per_atom

    def Nads_in_slab(self, ads_slab_entry):
        """
        Counts the TOTAL number of adsorbates in the slab on BOTH sides
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
        """

        return ads_slab_entry.composition.as_dict()\
            [str(self.adsorbate_entry.composition.\
                 reduced_composition.elements[0])]

    def Nsurfs_ads_in_slab(self, ads_slab_entry):
        """
        Counts the TOTAL number of adsorbed surfaces in the slab
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
        """

        struct = ads_slab_entry.structure
        weights = [s.species_and_occu.weight for s in struct]
        center_of_mass = np.average(struct.frac_coords,
                                    weights=weights, axis=0)

        Nsurfs = 0
        if any([site.species_string == self.adsorbate_as_str for site in
                struct if site.frac_coords[2] > center_of_mass[2]]):
            Nsurfs += 1
        if any([site.species_string == self.adsorbate_as_str for site in
                struct if site.frac_coords[2] < center_of_mass[2]]):
            Nsurfs += 1

        return Nsurfs


class SurfaceEnergyPlotter(object):
    """
    A class used for plotting the surface energies of a material by taking in
        a SurfaceEnergyCalculator object. Produces stability maps of different
        slab configurations, phases diagrams of two parameters to determine
        stability of configurations, and Wulff shapes

    .. attribute:: surface_energy_calculator

        SurfaceEnergyCalculator object to calculate quantities such as surface
            energy, binding energy, adsorption energy, etc. Also contains a
            vasprun dict as an attribute.

    .. attribute:: chempot_range

        List of the min and max chemical potential of ref_element.

    .. attribute:: entry_dict

        Nested dictionary containing a list of entries for slab calculations as
            items and the corresponding Miller index of the slab as the key.
            To account for adsorption, each value is a sub-dictionary with the
            entry of a clean slab calculation as the sub-key and a list of
            entries for adsorption calculations as the sub-value. The sub-value
            can contain different adsorption configurations such as a different
            site or a different coverage, however, ordinarily only the most stable
            configuration for a particular coverage will be considered as the
            function of the adsorbed surface energy has an intercept dependent on
            the adsorption energy (ie an adsorption site with a higher adsorption
            energy will always provide a higher surface energy than a site with a
            lower adsorption energy). An example parameter is provided:
            {(h1,k1,l1): {clean_entry1: [ads_entry1, ads_entry2, ...],
                          clean_entry2: [...], ...}, (h2,k2,l2): {...}}
            where clean_entry1 can be a pristine surface and clean_entry2 can be a
            reconstructed surface while ads_entry1 can be adsorption at site 1 with
            a 2x2 coverage while ads_entry2 can have a 3x3 coverage. If adsorption
            entries are present (i.e. if entry_dict[(h,k,l)][clean_entry1]), we
            consider adsorption in all plots and analysis for this particular facet.

    .. attribute:: ref_element

        All chemical potentials can be written in terms of the range of chemical
            potential of this element which will be used to calculate surface energy.

    """

    def __init__(self, entry_dict, surface_energy_calculator):
        """
        Object for plotting surface energy in different ways for clean and
            adsorbed surfaces.
        Args:
            entry_dict (dict): Dictionary containing a list of entries
                for slab calculations. See attributes.
            surface_energy_calculator (SurfaceEnergyCalculator):
                Object for calculating thermodynamic quantities related to surfaces
        """

        self.se_calculator = surface_energy_calculator
        self.chempot_range = surface_energy_calculator.chempot_range()
        self.entry_dict = entry_dict
        self.ref_el_comp = str(self.se_calculator.ref_el_comp.elements[0])

    def wulff_shape_from_chempot(self, u_ref=0, u_ads=0, symprec=1e-5):
        """
        Method to get the Wulff shape at a specific chemical potential.
        Args:
            chempot (float): The chemical potential the Wulff Shape exist in.
        """

        # Check if the user provided chemical potential is within the
        # predetermine range of chemical potential. If not, raise a warning
        if not max(self.chempot_range) >= u_ref >= min(self.chempot_range):
            warnings.warn("The provided chemical potential is outside the range "
                          "of chemical potential (%s to %s). The resulting Wulff "
                          "shape might not be reasonable." % (min(self.chempot_range),
                                                              max(self.chempot_range)))

        latt = SpacegroupAnalyzer(self.se_calculator.ucell_entry.structure). \
            get_conventional_standard_structure().lattice

        miller_list = self.entry_dict.keys()
        e_surf_list = []
        for hkl in miller_list:
            # At each possible configuration, we calculate surface energy as a
            # function of u and take the lowest surface energy (corresponds to
            # the most stable slab termination at that particular u)
            e_list = []
            for entry in self.entry_dict[hkl]:
                e_list.append(self.se_calculator.\
                              calculate_gamma_at_u(entry, u_ref=u_ref))
                for ads_entry in self.entry_dict[hkl][entry]:
                    e_list.append(self.se_calculator. \
                                  calculate_gamma_at_u(entry, ads_slab_entry=ads_entry,
                                                       u_ref=u_ref, u_ads=u_ads))

            e_surf_list.append(min(e_list))

        return WulffShape(latt, miller_list, e_surf_list, symprec=symprec)

    def return_stable_slab_entry_at_u(self, miller_index, u_ref=0, u_ads=0):
        """
        Returns the vasprun corresponding to the most
        stable slab for a particular facet at a specific u
        Args:
            miller_index ((h,k,l)): The facet to find the most stable slab in
            u (float): The chemical potential to look for the most stable slab
        """

        all_entries, all_gamma = [], []
        for entry in self.entry_dict[miller_index].keys():
            gamma = self.se_calculator.calculate_gamma_at_u(entry, u_ref=u_ref)
            all_entries.append(entry)
            all_gamma.append(gamma)
            for ads_entry in self.entry_dict[miller_index][entry]:
                gamma = self.se_calculator.calculate_gamma_at_u(entry, u_ref=u_ref,
                                                                u_ads=u_ads)
                all_entries.append(ads_entry)
                all_gamma.append(gamma)

        return all_entries[all_gamma.index(min(all_gamma))], min(all_gamma)

    def get_intersects_hkl(self, miller_index, x_is_u_ads=False,
                           const_u=0, return_entries=False):
        """
            Returns intersections for a specific facet. A list of lists, each
                entry in the list represents an intersect with u and gamma(u).
                Useful for finding when the configuration of a particular facet
                changes. Finds intersections by solving a matrix of equations,
                2 eqns for clean surfaces (b2=0 for no adsorption) and 3 eqns
                for adsorbed surfaces

            Args:
                miller_index ((h, k, l)): Miller index of the facet we
                    are interested in. Optional parameter that looks for
                    intersections for a specfic facet only when given.
        """

        # First lets calculate the range of surface energies for
        # all terminations of a specific facet.
        all_coefficients = []
        for clean_entry in self.entry_dict[miller_index]:
            all_coefficients.append(self.se_calculator. \
                                    surface_energy_coefficients(clean_entry))
            gbind_list, ads_entries = [], []
            for ads_entry in self.entry_dict[miller_index][clean_entry]:
                ads_entries.append(ads_entry)
                all_coefficients.append(self.se_calculator.\
                                        surface_energy_coefficients(clean_entry,
                                                                    ads_slab_entry=ads_entry))

        if len(all_coefficients) == 1:
            print("only one coeffs")
            return []

        # Now get all possible intersection coordinates for each pair of lines
        intersections = []
        # Find the intersections of two lines if no
        # adsorption, solve for three lines if adsorption
        for pair_lines in itertools.combinations(all_coefficients, 2):
            # matA is a 2x2 coefficient matrix, matB is a 2x1 intercept
            # matrix. For matA, we use the coeff of the ref chempot if
            # x is the ref chempot, else use the adsorbate chempot
            u = (pair_lines[1][2]-pair_lines[0][2])/(pair_lines[0][0]-pair_lines[1][0]) if not \
                x_is_u_ads else (pair_lines[1][2]-pair_lines[0][2])/(pair_lines[0][1]-pair_lines[1][1])

            if x_is_u_ads:
                u_ref = const_u
                u_ads = u
            else:
                u_ref = u
                u_ads = const_u

            # if the intersection is beyond the chemical potential
            # range or if the lines are parallel, we ignore it
            if not x_is_u_ads and any([u_ref < min(self.chempot_range),
                                       u_ref > max(self.chempot_range)]):
                print('skip')
                continue

            # If the surface energies at this u for both lines facet is unstable, ignore it
            gamma = np.dot(pair_lines[0], [u_ref, u_ads, 1])
            e1, gamma1 = self.return_stable_slab_entry_at_u(miller_index,
                                                            u_ref, u_ads=u_ads)
            e2, gamma2 = self.return_stable_slab_entry_at_u(miller_index,
                                                            u_ref, u_ads=u_ads)

            # +10e-9 to handle floating point comparison for equal values
            if all([gamma > gamma1 + 10e-9, gamma > gamma2 + 10e-9]):
                continue
            # Each list is [u_ref, u_ads, gamma, hkl1, hkl2]
            intersect = [u, gamma]
            if return_entries:
                intersect.append([e1, e2])
            intersections.append(intersect)

        return sorted(intersections, key=lambda ints: ints[0]) if not x_is_u_ads \
            else sorted(intersections, key=lambda ints: ints[1])

    # def get_intersections(self, miller_index=(), clean_only=True):
    #     """
    #     Returns all intersections for a specific facet or for all facets.
    #         A list of lists, each entry in the list represents an intersect
    #         with u, gamma(u) and [hkl1, hkl2]. Useful for finding when the
    #         configuration of a particular facet changes or when one facet
    #         is more stable than the other. Finds intersections by solving a
    #         matrix of equations, 2 eqns for clean surfaces (b2=0 for no
    #         adsorption) and 3 eqns for adsorbed surfaces
    #
    #     Args:
    #         miller_index ((h, k, l)): Miller index of the facet we
    #             are interested in. Optional parameter that looks for
    #             intersections for a specfic facet only when given.
    #     """
    #
    #     # First lets calculate the range of surface energies for
    #     # all terminations of a specific facet or all facets.
    #     all_coefficients = []
    #     for hkl in self.entry_dict.keys():
    #         if miller_index and hkl != tuple(miller_index):
    #             continue
    #         for entry in self.entry_dict[hkl]:
    #             all_coefficients.append([self.se_calculator.\
    #                                     surface_energy_coefficients(entry),
    #                                      hkl])
    #             if not clean_only:
    #                 for ads_entry in self.entry_dict[hkl][entry]:
    #                     all_coefficients.append([\
    #                         self.se_calculator.\
    #                             surface_energy_coefficients(entry,
    #                                                         ads_slab_entry=ads_entry),
    #                         hkl])
    #
    #     if len(all_coefficients) == 1:
    #         return []
    #
    #     # Now get all possible intersection coordinates for each pair of lines
    #     intersections = []
    #     # Find the intersections of two lines if no
    #     # adsorption, solve for three lines if adsorption
    #     comb_type = 2 if clean_only else 3
    #     for pair_lines in itertools.combinations(all_coefficients, comb_type):
    #         # matA is the coefficient matrix, matB is the intercept matrix
    #         matA, matB = [], []
    #         for coeffs in pair_lines:
    #             coeffs = coeffs[0]
    #             matB.append(coeffs[2])
    #             if clean_only:
    #                 matA.append([coeffs[0], -1])
    #             else:
    #                 matA.append([coeffs[0], coeffs[1] -1])
    #
    #         # If the coefficient vectors in matrix A are equal, continue
    #         if all([tuple(coeff_vect[0])==tuple(coeff_vect[1])
    #                 for coeff_vect in itertools.combinations(matA, 2)]):
    #             continue
    #
    #         # Calculate the intersection
    #         if clean_only:
    #             u_ref, gamma = np.linalg.solve(matA, matB)
    #         else:
    #             u_ref, u_ads, gamma = np.linalg.solve(matA, matB)
    #
    #         # if the intersection is beyond the chemical potential
    #         # range or if the lines are parallel, we ignore it
    #         if u_ref < min(self.chempot_range) \
    #                 or u_ref > max(self.chempot_range):
    #             continue
    #
    #         # If the surface energies at this u for both lines facet is unstable, ignore it
    #         u_ads = 0 if clean_only else u_ads
    #         e1, gamma1 = self.return_stable_slab_entry_at_u(pair_lines[0][1], u_ref,
    #                                                         u_ads=u_ads, clean_only=clean_only)
    #         e2, gamma2 = self.return_stable_slab_entry_at_u(pair_lines[1][1], u_ref,
    #                                                         u_ads=u_ads, clean_only=clean_only)
    #
    #         # +10e-9 to handle floating point comparison for equal values
    #         if all([gamma > gamma1 + 10e-9, gamma > gamma2 + 10e-9]):
    #             continue
    #         # Each list is [u_ref, u_ads, gamma, hkl1, hkl2]
    #         intersections.append([u_ref, u_ads, gamma,
    #                               [pair_lines[0][1], pair_lines[1][1]]])
    #
    #     return sorted(intersections, key=lambda ints: ints[0]) if clean_only \
    #         else sorted(intersections, key=lambda ints: ints[1])

    # def get_stable_surf_regions(self, miller_index, clean_only=True):
    #     """
    #     For a specific facet, returns an energy stability range of u for each
    #         facet as a slope and intercept. e.g. if config 1 stable between
    #         -1<u<-0.5 and config 2 is stable between -0.5<u<0 we return:
    #         [[-1, -0.5, slope1, intercept1], [-0.5, 0, slope2, intercept2]].
    #     """
    #
    #     # First get the intercepts at the stable configs
    #     stable_slab_entry, gamma = self.return_stable_slab_entry_at_u(miller_index,
    #                                                                   min(self.chempot_range),
    #                                                                   clean_only=clean_only)
    #     ulist, gamma_list = [min(self.chempot_range)], [gamma]
    #     intersections = self.get_intersections(miller_index=miller_index, clean_only=clean_only)
    #     for int in intersections:
    #         stable_slab_entry, gamma = self.return_stable_slab_entry_at_u(miller_index, int[0],
    #                                                                       clean_only=clean_only)
    #         ulist.append(int[0])
    #         gamma_list.append(gamma)
    #
    #     # Next, build a stability map with the range in u, slope and intercept
    #     stability_map = []
    #     ent, final_gamma = self.return_stable_slab_entry_at_u(miller_index,
    #                                                           max(self.chempot_range),
    #                                                           clean_only=clean_only)
    #     for i, u in enumerate(ulist):
    #         high_u = max(self.chempot_range) if i == len(ulist) - 1 else ulist[i + 1]
    #         high_gamma = final_gamma if i == len(ulist) - 1 else gamma_list[i + 1]
    #         stability_map.append([[u, high_u],
    #                               [gamma_list[i], high_gamma]])
    #
    #     return stability_map

    # def wulff_shape_dict(self, symprec=1e-5, at_intersections=False):
    #     """
    #     As the surface energy is a function of chemical potential, so too is the
    #         Wulff shape. This methods generates a dictionary of Wulff shapes at
    #         certain chemical potentials where a facet goes through a transition.
    #         Returns a dict, eg. {chempot1: WulffShape1, chempot2: WulffShape2}
    #
    #     Args:
    #         symprec (float): for recp_operation, default is 1e-5.
    #         at_intersections (bool): Whether to generate a Wulff shape for each
    #             intersection of surface energy for a specific facet (eg. at the
    #             point where a (111) stoichiometric surface energy plot intersects
    #             with the (111) nonstoichiometric plot) or to just generate two
    #             Wulff shapes, one at the min and max chemical potential.
    #     """
    #
    #     # First lets get the Wulff shape at the
    #     # minimum and maximum chemical potential
    #     wulff_dict = {self.chempot_range[0]: \
    #                       self.wulff_shape_from_chempot(self.chempot_range[0],
    #                                                     symprec=symprec),
    #                   self.chempot_range[1]: \
    #                       self.wulff_shape_from_chempot(self.chempot_range[1],
    #                                                     symprec=symprec)}
    #
    #     # Now we get the Wulff shape each time a facet changes its configuration
    #     # (ie, adsorption coverage, stoichiometric to nonstoichiometric, etc)
    #     if at_intersections:
    #         # Get all values of chemical potential where an intersection occurs
    #         all_intersections = self.get_intersections()
    #         # Get a Wulff shape for each intersection. The change in the Wulff shape
    #         # will vary if the rate of change in surface energy for any facet changes
    #         for int in all_intersections:
    #             wulff = self.wulff_shape_from_chempot(int[0], symprec=symprec)
    #             if any([wulff.area_fraction_dict[hkl] != 0 for hkl in int[2]]):
    #                 wulff_dict[int[0]] = wulff
    #
    #     return wulff_dict

    def area_frac_vs_chempot_plot(self, cmap=cm.jet, at_intersections=False, xrange=None,
                                  increments=10, x_is_u_ads=False):
        """
        Plots the change in the area contribution of
        each facet as a function of chemical potential.
        Args:
            cmap (cm): A matplotlib colormap object, defaults to jet.
            at_intersections (bool): Whether to generate a Wulff shape for each
                intersection of surface energy for a specific facet (eg. at the
                point where a (111) stoichiometric surface energy plot intersects
                with the (111) nonstoichiometric plot) or to just generate two
                Wulff shapes, one at the min and max chemical potential.
            increments (bool): Number of data points between min/max or point
                of intersection. Defaults to 5 points.
        """

        # Choose unique colors for each facet
        f = [int(i) for i in np.linspace(0, 255, len(self.entry_dict.keys()))]

        # Get all points of min/max chempot and intersections
        # chempot_intersections = []
        # chempot_intersections.extend(self.chempot_range)
        # for hkl in self.entry_dict.keys():
        #     chempot_intersections.extend([ints[0] for ints in
        #                                   self.get_intersections(hkl)])
        # chempot_intersections = sorted(chempot_intersections)

        # # Get all chempots
        # if at_intersections:
        #     all_chempots = []
        #     for i, intersection in enumerate(chempot_intersections):
        #         if i < len(chempot_intersections) - 1:
        #             all_chempots.extend(np.linspace(intersection,
        #                                             chempot_intersections[i + 1],
        #                                             increments))
        # else:
        xrange = self.chempot_range if not xrange else xrange
        all_chempots = np.linspace(min(xrange),
                                   max(xrange), increments)

        # initialize a dictionary of lists of fractional areas for each hkl
        hkl_area_dict = {}
        for hkl in self.entry_dict.keys():
            hkl_area_dict[hkl] = []

        # Get plot points for each Miller index
        for u in all_chempots:
            wulffshape = self.wulff_shape_from_chempot(u_ads=u) if \
                x_is_u_ads else self.wulff_shape_from_chempot(u_ref=u)

            for hkl in wulffshape.area_fraction_dict.keys():
                hkl_area_dict[hkl].append(wulffshape.area_fraction_dict[hkl])

        # Plot the area fraction vs chemical potential for each facet
        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        for i, hkl in enumerate(self.entry_dict.keys()):
            # Ignore any facets that never show up on the
            # Wulff shape regardless of chemical potential
            if all([a == 0 for a in hkl_area_dict[hkl]]):
                continue
            else:
                plt.plot(all_chempots, hkl_area_dict[hkl],
                         '--', color=cmap(f[i]), label=str(hkl))

        # Make the figure look nice
        plt.ylabel(r"Fractional area $A^{Wulff}_{hkl}/A^{Wulff}$")
        self.chempot_plot_addons(plt, xrange, axes, pad=5,
                                 rect=[-0.0, 0, 0.95, 1],
                                 x_is_u_ads=x_is_u_ads)

        return plt

    def chempot_vs_gamma_plot_one(self, plt, clean_entry, u_ref_range, label='',
                                  ads_entry=None, color='r', JPERM2=False,
                                  x_is_u_ads=False, u_ads_range=[-5,0]):

        # use dashed lines for slabs that are not stoichiometric
        # wrt bulk. Label with formula if nonstoichiometric
        mark = '--' if clean_entry.composition.reduced_composition != \
                       self.se_calculator.ucell_entry. \
                           composition.reduced_composition else '-'

        # Get the rise and run of the plot
        gamma_range = [self.se_calculator.calculate_gamma_at_u(clean_entry,
                                                               ads_slab_entry=ads_entry,
                                                               u_ref=u_ref_range[0],
                                                               u_ads=u_ads_range[0]),
                       self.se_calculator.calculate_gamma_at_u(clean_entry,
                                                               ads_slab_entry=ads_entry,
                                                               u_ref=u_ref_range[1],
                                                               u_ads=u_ads_range[1])]
        print(gamma_range)
        se_range = np.array(gamma_range) * EV_PER_ANG2_TO_JOULES_PER_M2 \
            if JPERM2 else gamma_range
        xrange = u_ref_range if not x_is_u_ads else u_ads_range
        plt.plot(xrange, se_range, mark, color=color, label=label)

        return plt

    def chempot_vs_gamma_clean(self, miller_index=(), cmap=cm.jet, JPERM2=False):
        """
        Plots the surface energy of all facets as a function of chemical potential.
            Each facet will be associated with its own distinct colors. Dashed lines
            will represent stoichiometries different from that of the mpid's compound.

        Args:
            cmap (cm): A matplotlib colormap object, defaults to jet.
            highlight_stability (bool): For each facet, there may be various
                terminations or stoichiometries and the relative stability of
                these different slabs may change with chemical potential. This
                dict (if provided) will highlight the stability line for a facet.
                The key of the dict is the Miller index we want to highlight while
                the value is the color of the highlight e.g. {(1,1,1): 'r'}.
        """

        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        # Choose unique colors for each facet
        colors = self.color_palette(miller_index=miller_index, cmap=cmap)

        # Now we plot each individual slab surface energy
        already_labelled = []
        for hkl in self.entry_dict.keys():
            if miller_index and hkl != tuple(miller_index):
                continue

            # Plot the clean slabs
            for entry in self.entry_dict[hkl]:
                # Generate a label for the type of slab
                label = self.create_slab_label(entry, miller_index=\
                    () if miller_index else hkl)

                if label in already_labelled:
                    c = colors[already_labelled.index(label)]
                    label = None
                else:
                    c = colors[len(already_labelled)]
                    already_labelled.append(label)

                self.chempot_vs_gamma_plot_one(plt, entry, self.chempot_range, color=c,
                                               label=label, JPERM2=JPERM2, x_is_u_ads=False)

        # Make the figure look nice
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, self.ref_el_comp,
                                       self.chempot_range, axes)

        return plt

    def chempot_vs_gamma_facet(self, miller_index, cmap=cm.jet, x_is_u_ads=False,
                               JPERM2=False, show_stable=False, const_u=0):

        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        # Choose unique colors for each facet
        colors = self.color_palette(cmap=cmap,
                                    miller_index=miller_index)
        intersections = self.get_intersects_hkl(miller_index, x_is_u_ads=x_is_u_ads,
                                                const_u=const_u, return_entries=True)
        u_ref_range = [const_u, const_u] if x_is_u_ads else self.chempot_range
        u_ads_range = [const_u, const_u] if not x_is_u_ads else [-5,0]

        already_labelled = []
        for clean_entry in self.entry_dict[miller_index]:

            label = self.create_slab_label(clean_entry)
            if label in already_labelled:
                c = colors[already_labelled.index(label)]
                label = None
            else:
                c = colors[len(already_labelled)]
                already_labelled.append(label)

            xrange = []
            if show_stable:
                for intersect in intersections:
                    if clean_entry in intersect[2]:
                        xrange.append(intersect[0])
                if x_is_u_ads:
                    u_ads_range = xrange
                else:
                    u_ref_range = xrange
                if len(xrange) < 2:
                    continue
                print("xrange", xrange)

            self.chempot_vs_gamma_plot_one(plt, clean_entry, u_ref_range,
                                           color=c, label=label, JPERM2=JPERM2,
                                           x_is_u_ads=x_is_u_ads, u_ads_range=u_ads_range)

            for ads_entry in self.entry_dict[miller_index][clean_entry]:
                # Plot the adsorbed slabs
                # Generate a label for the type of slab
                label = self.create_slab_label(clean_entry, ads_entry=ads_entry)
                if label in already_labelled:
                    c = colors[already_labelled.index(label)]
                    label = None
                else:
                    c = colors[len(already_labelled)]
                    already_labelled.append(label)

                xrange = []
                if show_stable:
                    for intersect in intersections:
                        if ads_entry in intersect[2]:
                            xrange.append(intersect[0])
                    if x_is_u_ads:
                        u_ads_range = xrange
                    else:
                        u_ref_range = xrange

                self.chempot_vs_gamma_plot_one(plt, clean_entry, u_ref_range,
                                               ads_entry=ads_entry, color=c,
                                               label=label, JPERM2=JPERM2,
                                               x_is_u_ads=x_is_u_ads)

        # Make the figure look nice
        xrange = self.chempot_range if not x_is_u_ads else [-5, 0]
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, xrange, axes,
                                       x_is_u_ads=x_is_u_ads)
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        clean_se = self.se_calculator.calculate_gamma_at_u(clean_entry)
        plt.annotate(miller_index, xy=[np.mean(xlim), clean_se+max(ylim)*0.1],
                     xytext=[np.mean(xlim), clean_se+max(ylim)*0.1],
                     fontsize=20)

        return plt

    def color_palette(self, miller_index=(), cmap=cm.jet):

        total_surfaces = 0
        for hkl in self.entry_dict.keys():
            if miller_index and hkl != tuple(miller_index):
                continue
            total_surfaces += len(self.entry_dict[hkl].keys())
            for entry in self.entry_dict[hkl].keys():
                total_surfaces += len(self.entry_dict[hkl][entry])
        return [cmap(int(i)) for i in np.linspace(0, 255, total_surfaces)]

    def chempot_plot_addons(self, plt, xrange, axes,
                            pad=2.4, rect=[-0.047, 0, 0.84, 1], x_is_u_ads=False):

        # Make the figure look nice
        x_species = self.ref_el_comp if not \
            x_is_u_ads else self.se_calculator.adsorbate_as_str
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        axes.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" % (x_species))

        ylim = axes.get_ylim()
        plt.xticks(rotation=60)
        plt.ylim(ylim)
        xlim = axes.get_xlim()
        plt.xlim(xlim)
        plt.tight_layout(pad=pad, rect=rect)
        plt.plot([xrange[0], xrange[0]], ylim, '--k')
        plt.plot([xrange[1], xrange[1]], ylim, '--k')
        xy = [np.mean([xrange[1]]), np.mean(ylim)]
        plt.annotate("%s-rich" % (x_species), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)
        xy = [np.mean([xlim[0]]), np.mean(ylim)]
        plt.annotate("%s-poor" % (x_species), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)

        return plt

    def create_slab_label(self, entry, miller_index=(), ads_entry=None):

        label = str(miller_index) if miller_index else ""
        if entry.composition.reduced_composition != \
                self.se_calculator.ucell_entry.composition.reduced_composition:
            label += " %s" % (entry.composition.reduced_composition)
        if ads_entry:
            label += " +%s, ML=%.2f" %(self.se_calculator.adsorbate_as_str,
                                       self.se_calculator.get_monolayer(ads_entry,
                                                                        entry))
            print(label, ads_entry.composition, entry.composition)

        return label

    def surface_phase_diagram(self, y_param, x_param, miller_index):

        """
        Builds a 2D surface phase diagram of two parameters for a specific
            facet. Parameters can be chemical potentials (e.g. u_a, u_b,
            u_c), pressure, temperature, electric potential, pH, etc.

        Args:
        """

        return

    def wulff_shape_extrapolated_model(self):
        # the area fraction of all facets for a single configuration can be modelled as
        # a arctan function with a pseudo asymptote approaching 0. This model gives us a
        return

    def surface_pourbaix_diagram(self):

        return

    def surface_u_vs_u_phase_diagram(self):

        return

    def surface_p_vs_t_phase_diagram(self):

        return

    def broken_bond_vs_gamma(self):

        return

def vaspruns_to_entry_dict(vaspruns):
    """
    Helper function to generate the entry_dict parameter for
        SurfaceEnergyCalculator from a list of Vasprun objects.
    """

    return