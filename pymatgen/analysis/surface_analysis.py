# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
TODO:
    -Only works for monatomic adsorbates
    -Need a method to automatically get chempot range when
        dealing with non-stoichiometric slabs
    -Simplify or clean the stable_u_range_dict() method
        for the SurfaceEnergyPlotter class
    -Combine chempot_vs_gamma_clean(), chempot_vs_gamma_plot_one()
        and chempot_vs_gamma_facet() into one class method.
    -Annotations for which configuration is on a plot
"""

from __future__ import division, unicode_literals

import numpy as np
import itertools
import warnings
import random, copy
from sympy import Symbol, Number

from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.wulff import WulffShape
from pymatgen.util.plotting import pretty_plot
from pymatgen.analysis.reaction_calculator import Reaction

EV_PER_ANG2_TO_JOULES_PER_M2 = 16.0217656

__author__ = "Richard Tran"
__copyright__ = "Copyright 2017, The Materials Virtual Lab"
__version__ = "0.2"
__maintainer__ = "Richard Tran"
__credits__ = "Joseph Montoya, Xianguo Li"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "8/24/17"


"""
This module defines tools to analyze surface and adsorption related
quantities as well as related plots. If you use this module, please
consider citing the following works::

    R. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,
        S. P. Ong, "Surface Energies of Elemental Crystals", Scientific
        Data, 2016, 3:160080, doi: 10.1038/sdata.2016.80.

    and

    Montoya, J. H., & Persson, K. A. (2017). A high-throughput framework
        for determining adsorption energies on solid surfaces. Npj
        Computational Materials, 3(1), 14.
        https://doi.org/10.1038/s41524-017-0017-z
"""

class MatrixError(Exception):
    def __init__(self):
        print('Number of equations (slab_entries) is not \
        equal to number of free variables (chempot and surface energy)')


class SlabEntry(ComputedStructureEntry):
    """
    An object encompassing all data relevant to slab for surface calculations.

    Args:
        entry (ComputedEntry/ComputedStructureEntry): An
            entry object
    """

    def __init__(self, entry, miller_index, name=None,
                 coverage=None, adsorbate=None):
        self.entry = entry
        self.miller_index = miller_index
        self.name = name
        self.coverage = coverage
        self.adsorbate = adsorbate

        super(SlabEntry, self).__init__(
            entry.structure, entry.energy, correction=0.0,
            parameters=None, data=None, entry_id=None)

    def as_dict(self):
        """
        Returns dict which contains Slab Entry data.
        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["entry"] = self.entry.as_dict()
        d["miller_index"] = miller_index
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Returns a SlabEntry by reading in an slab
        """

        entry = SlabEntry.from_dict(d["entry"])
        miller_index = d["miller_index"]
        name = d["name"]
        coverage = d["coverage"]
        adsorbate = d["adsorbate"]

        return SlabEntry(entry, miller_index, name=name,
                         coverage=coverage, adsorbate=adsorbate)

    @property
    def surface_area(self):
        """
        Calculates the surface area of the slab
        """
        m = self.structure.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))


class SurfaceEnergyCalculator(object):
    """
    A class used for analyzing the surface energies, adsorbate binding energies
        and other quantities related to surface thermodynamics for a material.
        By default, this code will only use one bulk entry to calculate surface
        energy. Ideally, to get the most accurate surface energy, the user should
        compare their slab energy to the energy of the oriented unit cell with
        both calculations containing consistent k-points to avoid convergence
        problems as the slab size is varied. Note that this is for symmetric
        surfaces, i.e. for the case of adsorption, it is assumed both sides are
        adsorbed. See:

            Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
                Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.
        and
            Rogal, J., & Reuter, K. (2007). Ab Initio Atomistic Thermodynamics for
                Surfaces : A Primer. Experiment, Modeling and Simulation of Gas-Surface
                Interactions for Reactive Flows in Hypersonic Flights, 2–1 – 2–18.

    .. attribute:: ref_el_entry

        All chemical potentials can be written in terms of the range of chemical
            potential of this element which will be used to calculate surface energy.

    .. attribute:: ref_el_as_str

        Reference element as a string (see ref_el_entry).

    .. attribute:: ucell_entry

        Bulk entry of the base material of the slab.

    .. attribute:: x

        Reduced amount composition of decomposed compound A in the
            bulk. e.g. x=1 for FePO4 in LiFePO4

    .. attribute:: y

        Reduced amount composition of ref_element in the bulk. e.g.
            y=2 for Li in LiFePO4

    .. attribute:: gbulk

        Gibbs free energy of the bulk per formula unit

    .. attribute:: e_of_element

        Energy per atom of ground state ref_element, eg. if ref_element=O,
            than e_of_element=1/2*E_O2.

    .. attribute:: adsorbate_entry

        Entry of the adsorbate (if there is one).

    .. attribute:: adsorbate_as_str

        Adsorbate formula as a string (see adsorbate_entry).
    """

    def __init__(self, ucell_entry, ref_entries=[], adsorbate_entries=[]):
        """
        Analyzes surface energies and Wulff shape of a particular
            material using the chemical potential.
        Args:
            ucell_entry (material_id or computed_entry): Materials Project or entry
                of the bulk system the slab is based on (a string, e.g., mp-1234).
            ref_el_entry (ComputedStructureEntry): Entry to be considered as
                independent variable. E.g., if you want to show the stability
                ranges of all Li-Co-O phases wrt to uLi
            adsorbate_entries (ComputedStructureEntry): Computed entry of adsorbate,
                defaults to None. Could be an isolated H2 or O2 molecule for gaseous
                    adsorption or just the bulk ground state structure for metals.
                    Either way, the extracted term is going to be in energy per atom.
        """

        # Set up basic attributes
        self.ucell_entry = ucell_entry
        self.ref_entries = ref_entries
        self.ucell_comp = self.ucell_entry.composition
        self.ucell_reduced_comp = self.ucell_comp.reduced_composition
        self.adsorbate_entries_dict = {list(ads_entry.composition.as_dict().keys())[0]: \
                                           ads_entry for ads_entry in adsorbate_entries}

        # Calculate Gibbs free energy of the bulk per unit formula
        self.gbulk = self.ucell_entry.energy / \
                     self.ucell_comp.get_integer_formula_and_factor()[1]

    def surface_energy_coefficients(self, slab_entry):
        """
        Calculates the surface energy for a single slab.
        Args:
            clean_slab_entry (entry): An entry object for the clean slab
            bulk_entry (entry): An entry object for the bulk
            ref_entries (list: [entry]): A list of entries for each type
                of element to be used as a reservoir for nonstoichiometric
                systems. The length of this list MUST be n-1 where n is the
                number of different elements in the bulk entry. The chempot
                of the element ref_entry that is not in the list will be
                treated as a variable.
            ref_adsorbate_entry (entry): Entry for the adsorbate to be used
                as a reservoir during adsorption.
            chempot (dict): A dictionary with the key being an element
                existing in the slab system. All elements will have a
                chempot variable, however most of these variables will cancel
                out, so the chempot will have no effect. By default, if a
                chempot is not set for an element, we set it to 0. This only
                effects the surface energy if a reservoir (either for an
                adsorbate or an element in the bulk) is needed.

        Returns (float): Surface energy
        """

        # Add missing chemical potentials for elements
        # not present in chempot, defaults to 0.
        chempot_vars = {el: Symbol(el) for el in \
                        slab_entry.composition.as_dict().keys()}

        # Calculate Gibbs free energy of the bulk per unit formula
        ref_entries_dict = {list(entry.composition.as_dict().keys())[0]: \
                                entry for entry in self.ref_entries}
        bulk_energy = 0
        # First we get the contribution to the bulk energy
        # from each element with an existing ref_entry.
        for el in ref_entries_dict.keys():
            entry = ref_entries_dict[el]
            N = slab_entry.composition.as_dict()[el]
            bulk_energy += N * (chempot_vars[el] + entry.energy_per_atom)

        # Add the reservoir from the adsorbate to the bulk energy
        if slab_entry.adsorbate:
            N = slab_entry.composition.as_dict()[slab_entry.adsorbate]
            bulk_energy += N * (chempot_vars[slab_entry.adsorbate] - \
                                self.adsorbate_entries_dict[slab_entry.adsorbate] \
                                .energy_per_atom)

        # Next, we add the contribution to the bulk energy from
        # the variable element (the element without a ref_entry),
        # as a function of the other elements
        for el in self.ucell_entry.composition.as_dict().keys():
            if el not in ref_entries_dict.keys():
                ref_el = el

        N = slab_entry.composition.as_dict()[ref_el]
        n = self.ucell_reduced_comp.as_dict()[ref_el]
        bulk_energy += (N / n) * (self.gbulk - \
                                  sum([self.ucell_reduced_comp[el] * \
                                       (chempot_vars[el] + \
                                        ref_entries_dict[el].energy_per_atom) \
                                       for el in ref_entries_dict.keys()]))

        # Full equation of the surface energy (constant if stoichiometric)
        se = (slab_entry.energy - bulk_energy) / (2 * slab_entry.surface_area)

        # Return dict of a constant if the equation is constant
        if type(se).__name__ != "Add":
            return {Number(1): se}
        else:
            # Remove any variables with coefficient of 0
            # (Sympy doesn't handle cancellation very well)
            coefficients = {el: coeff for el, coeff in \
                            se.as_coefficients_dict().items() if abs(coeff) > 1e-6}
            return coefficients

    def calculate_gamma(self, slab_entry, u_dict={}, u_default=0):
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

        coeffs = self.surface_energy_coefficients(slab_entry)
        coeff_vals = [coeffs[el] for el in coeffs.keys()]
        u_vals = []
        for el in coeffs.keys():
            if el in u_dict.keys():
                u_vals.append(u_dict[el])
            elif el == Number(1):
                u_vals.append(1)
            else:
                u_vals.append(u_default)
        return np.dot(u_vals, coeff_vals)

    def get_unit_primitive_area(self, ads_slab_entry, clean_slab_entry):
        """
        Returns the surface area of the adsorbed system per
        unit area of the primitive slab system.
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
            clean_slab_entry (entry): The entry of the clean slab
        """

        A_ads = ads_slab_entry.surface_area
        A_clean = clean_slab_entry.surface_area
        n = (A_ads / A_clean)
        return n

    def get_monolayer(self, ads_slab_entry, clean_slab_entry):
        """
        Returns the primitive unit surface area density of the
            adsorbate.
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
            clean_slab_entry (entry): The entry of the clean slab
        """

        unit_a = self.get_unit_primitive_area(ads_slab_entry, clean_slab_entry)
        Nsurfs = self.Nsurfs_ads_in_slab(ads_slab_entry)
        Nads = self.Nads_in_slab(ads_slab_entry)
        return Nads / (unit_a * Nsurfs)

    def gibbs_binding_energy(self, ads_slab_entry,
                             clean_slab_entry, eads=False):
        """
        Returns the adsorption energy or Gibb's binding energy
            of an adsorbate on a surface
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
            clean_slab_entry (entry): The entry of the clean slab
            eads (bool): Whether to calculate the adsorption energy
                (True) or the binding energy (False) which is just
                adsorption energy normalized by number of adsorbates.
        """

        n = self.get_unit_primitive_area(ads_slab_entry, clean_slab_entry)
        Nads = self.Nads_in_slab(ads_slab_entry)

        BE = (ads_slab_entry.energy - n * clean_slab_entry.energy) / Nads \
             - self.adsorbate_entries_dict[ads_slab_entry.adsorbate].energy_per_atom
        return BE * Nads if eads else BE

    def Nads_in_slab(self, ads_slab_entry):
        """
        Returns the TOTAL number of adsorbates in the slab on BOTH sides
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
        """
        print(ads_slab_entry.adsorbate, ads_slab_entry.composition.as_dict())
        return ads_slab_entry.composition.as_dict()[ads_slab_entry.adsorbate]

    def Nsurfs_ads_in_slab(self, ads_slab_entry):
        """
        Returns the TOTAL number of adsorbed surfaces in the slab
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
        """

        struct = ads_slab_entry.structure
        weights = [s.species_and_occu.weight for s in struct]
        center_of_mass = np.average(struct.frac_coords,
                                    weights=weights, axis=0)

        Nsurfs = 0
        if any([site.species_string == ads_slab_entry.adsorbate for site in
                struct if site.frac_coords[2] > center_of_mass[2]]):
            Nsurfs += 1
        if any([site.species_string == ads_slab_entry.adsorbate for site in
                struct if site.frac_coords[2] < center_of_mass[2]]):
            Nsurfs += 1

        return Nsurfs

    def get_surface_equilibrium(self, slab_entries, u_dict={}):

        """
        Takes in a list of SlabEntries and calculates the chemical potentials
            at which all slabs in the list coexists simultaneously. Useful for
            building surface phase diagrams. Note that to solve for x equations
            (x slab_entries), there must be x free variables (chemical potentials).
            Adjust u_dict as need be to get the correct number of free variables.
        Args:
            slab_entries (array): The coefficients of the first equation
            u_dict (dict): Dictionary of chemical potentials to keep constant.
                The key is the element and the value is the chemical potential.

        Returns:
            (array): Array containing a solution to x equations with x
                variables (x-1 chemical potential and 1 surface energy)
        """

        # Generate all possible coefficients
        all_parameters = []
        all_coeffs = []
        for slab_entry in slab_entries:
            coeffs = self.surface_energy_coefficients(slab_entry)

            # remove the free chempots we wish to keep constant
            for el in u_dict.keys():
                if el in coeffs.keys():
                    coeffs[1] += coeffs[el] * u_dict[el]
                    del coeffs[el]
            all_coeffs.append(coeffs)

        # Find all possible free variables
        for coeffs in all_coeffs:
            for el in coeffs.keys():
                if el not in all_parameters and el != Number(1):
                    all_parameters.append(el)

        # if there are no variables in the parameter list, then
        # it means the equations are constant (i.e. parallel)
        if not all_parameters:
            print("lines are parallel", all_parameters)
            return None

        print(all_parameters, "all_parameters", slab_entries[0].composition)
        # Check if its even possible for the system
        # of equations to even have a solution
        if len(slab_entries) != len(all_parameters) + 1:
            raise MatrixError()

        # Set up the matrix
        coeffs_matrix, const_vector = [], []
        for coeffs in all_coeffs:
            coeff_vector = []
            for el in all_parameters:
                if el in coeffs.keys():
                    coeff_vector.append(float(coeffs[el]))
                else:
                    coeff_vector.append(float(0))
            coeff_vector.append(float(-1))
            coeffs_matrix.append(coeff_vector)
            const_vector.append(float(coeffs[1]))

        # Solve for the chempots and surface energy
        print(coeffs_matrix, const_vector, "A, B")
        solution = np.linalg.solve(coeffs_matrix, const_vector)

        soln = {str(el): solution[i] for i, el in enumerate(all_parameters) \
                if i < len(solution)-1}
        soln["gamma"] = solution[-1]
        print("solution", soln, solution, all_parameters)
        return soln


# class SurfaceEnergyPlotter(object):
#     """
#     A class used for generating plots to analyze the thermodynamics of surfaces
#         of a material by taking in a SurfaceEnergyCalculator object. Produces
#         stability maps of different slab configurations, phases diagrams of two
#         parameters to determine stability of configurations, and Wulff shapes.
#
#     .. attribute:: se_calculator
#
#         SurfaceEnergyCalculator object to calculate quantities such as surface
#             energy, binding energy, adsorption energy, etc.
#
#     .. attribute:: chempot_range
#
#         List of the min and max chemical potential of ref_element.
#
#     .. attribute:: entry_dict
#
#         Nested dictionary containing a list of entries for slab calculations as
#             items and the corresponding Miller index of the slab as the key.
#             To account for adsorption, each value is a sub-dictionary with the
#             entry of a clean slab calculation as the sub-key and a list of
#             entries for adsorption calculations as the sub-value. The sub-value
#             can contain different adsorption configurations such as a different
#             site or a different coverage, however, ordinarily only the most stable
#             configuration for a particular coverage will be considered as the
#             function of the adsorbed surface energy has an intercept dependent on
#             the adsorption energy (ie an adsorption site with a higher adsorption
#             energy will always provide a higher surface energy than a site with a
#             lower adsorption energy). An example parameter is provided:
#             {(h1,k1,l1): {clean_entry1: [ads_entry1, ads_entry2, ...],
#                           clean_entry2: [...], ...}, (h2,k2,l2): {...}}
#             where clean_entry1 can be a pristine surface and clean_entry2 can be a
#             reconstructed surface while ads_entry1 can be adsorption at site 1 with
#             a 2x2 coverage while ads_entry2 can have a 3x3 coverage. If adsorption
#             entries are present (i.e. if entry_dict[(h,k,l)][clean_entry1]), we
#             consider adsorption in all plots and analysis for this particular facet.
#
#     ..attribute:: color_dict
#
#         Dictionary of colors (r,g,b,a) when plotting surface energy stability. The
#             keys are individual surface entries where clean surfaces have a solid color
#             while the corresponding adsorbed surface will be transparent.
#
#     .. attribute:: ref_el_comp
#
#         Composition of the ref_element. All chemical potentials can be written in terms
#             of the range of chemical potential of this element which will be used to
#             calculate surface energy.
#
#     """
#
#     def __init__(self, entry_dict, surface_energy_calculator, chempot_range):
#         """
#         Object for plotting surface energy in different ways for clean and
#             adsorbed surfaces.
#         Args:
#             entry_dict (dict): Dictionary containing a list of entries
#                 for slab calculations. See attributes.
#             surface_energy_calculator (SurfaceEnergyCalculator):
#                 Object for calculating thermodynamic quantities related to surfaces
#             chempot_range (list): Max and min range for the chemical potential for
#                 the ref_element.
#         """
#
#         self.se_calculator = surface_energy_calculator
#         self.chempot_range = chempot_range
#         self.entry_dict = entry_dict
#         self.color_dict = self.color_palette_dict()
#         self.ref_el_comp = str(self.se_calculator.ref_el_entry.composition.elements[0].name) if \
#             self.se_calculator.ref_el_entry else \
#             self.se_calculator.ucell_entry.composition.get_integer_formula_and_factor()[0]
#
#     def chempot_range_adsorption(self, ads_slab_entry, clean_slab_entry,
#                                  const_u_ref, buffer=0.2):
#         """
#         Returns the chemical potential range as a list for the adsorbate. The min
#             chemical potential will be located where the clean and adsorbed surface
#             energy lines intersect while the max chemical potential will depend on
#             the formation energy of the material.
#         Args:
#             ads_slab_entry (entry): The entry of the adsorbed slab
#             clean_slab_entry (entry): The entry of the clean slab
#             const_u_ref (float): Set the chemical potential of the
#                 ref element to a constant value.
#             buffer (float): A buffer fo r the x axis (chemical
#                 potential range). For plotting.
#         """
#
#         c1 = self.se_calculator.surface_energy_coefficients(clean_slab_entry)
#         c2 = self.se_calculator.surface_energy_coefficients(clean_slab_entry,
#                                               ads_slab_entry=ads_slab_entry)
#         umin, gamma = self.se_calculator.solve_2_linear_eqns(c1, c2, x_is_u_ads=True,
#                                                const_u=const_u_ref)
#         umin = -1*10e-5 if not umin else umin
#         # Make sure upper limit of u doesn't lead to negative or 0 values.
#         # Substract a value approaching 0 to avoid surface energies of 0
#         umax = (1*10e-5-const_u_ref*c2[0]-c2[2])/c2[1] if \
#             0 > self.se_calculator.calculate_gamma(clean_slab_entry,
#                                           ads_slab_entry=ads_slab_entry,
#                                           u_ref=const_u_ref, u_ads=0) else 0
#         return [umin-abs(umin)*buffer, umax]
#
#
#     def max_adsorption_chempot_range(self, const_u_ref, buffer=0.1):
#         """
#         Returns the chemical potential range as a list for the adsorbate among all
#             facets. The min chemical potential will be located where the clean and
#             adsorbed surface energy lines intersect while the max chemical potential
#             will depend on the formation energy of the material.
#         Args:
#             const_u_ref (float): Set the chemical potential of the
#                 ref element to a constant value.
#             buffer (float): A buffer fo r the x axis (chemical
#                 potential range). For plotting.
#         """
#
#         all_ranges = []
#         for hkl in self.entry_dict.keys():
#             for clean_entry in self.entry_dict[hkl].keys():
#                 for ads_entry in self.entry_dict[hkl][clean_entry]:
#                     all_ranges.append(self.\
#                                       chempot_range_adsorption(ads_entry, clean_entry,
#                                                                const_u_ref=const_u_ref,
#                                                                buffer=buffer))
#
#         if not all_ranges:
#             # If there is no intersection, the range is [-1,0]
#             return [-1,0]
#         # ensure our lower limit is at an intersection with the clean slab
#         all_ranges = sorted(all_ranges, key=lambda r: r[0])
#         max_range = [all_ranges[0][0]]
#         # ensure our upper limit corresponds to gamm > 0 or if gamma > 0 when u = 0
#         all_ranges = sorted(all_ranges, key=lambda r: r[1])
#         max_range.append(all_ranges[0][1])
#         return max_range
#
#     def wulff_shape_from_chempot(self, u_ref=0, u_ads=0, symprec=1e-5):
#         """
#         Method to get the Wulff shape at a specific chemical potential.
#         Args:
#             u_ref (float): The chemical potential of the reference element
#             u_ads (float): The chemical potential of the adsorbate
#             symprec (float): See WulffShape.
#
#         Returns:
#             (WulffShape): The WulffShape at u_ref and u_ads.
#         """
#
#         # Check if the user provided chemical potential is within the
#         # predetermine range of chemical potential. If not, raise a warning
#         if not max(self.chempot_range) >= u_ref >= min(self.chempot_range):
#             warnings.warn("The provided chemical potential is outside the range "
#                           "of chemical potential (%s to %s). The resulting Wulff "
#                           "shape might not be reasonable." % (min(self.chempot_range),
#                                                               max(self.chempot_range)))
#
#         latt = SpacegroupAnalyzer(self.se_calculator.ucell_entry.structure). \
#             get_conventional_standard_structure().lattice
#
#         miller_list = self.entry_dict.keys()
#         e_surf_list = []
#         for hkl in miller_list:
#             # At each possible configuration, we calculate surface energy as a
#             # function of u and take the lowest surface energy (corresponds to
#             # the most stable slab termination at that particular u)
#             entry, gamma = self.return_stable_slab_entry_at_u(hkl, u_ref=u_ref,
#                                                               u_ads=u_ads)
#             e_surf_list.append(gamma)
#
#         return WulffShape(latt, miller_list, e_surf_list, symprec=symprec)
#
#     def return_stable_slab_entry_at_u(self, miller_index, u_ref=0, u_ads=0):
#         """
#         Returns the entry corresponding to the most stable
#         slab for a particular facet at a specific chempot.
#
#         Args:
#             miller_index ((h,k,l)): The facet to find the most stable slab in
#             u_ref (float): The chemical potential of the reference element
#             u_ads (float): The chemical potential of the adsorbate
#         """
#
#         all_entries, all_gamma = [], []
#         for entry in self.entry_dict[miller_index].keys():
#             gamma = self.se_calculator.calculate_gamma(entry, u_ref=u_ref)
#             all_entries.append(entry)
#             all_gamma.append(gamma)
#             for ads_entry in self.entry_dict[miller_index][entry]:
#                 gamma = self.se_calculator.calculate_gamma(entry, u_ref=u_ref, u_ads=u_ads,
#                                                                 ads_slab_entry=ads_entry)
#                 all_entries.append(ads_entry)
#                 all_gamma.append(gamma)
#
#         return all_entries[all_gamma.index(min(all_gamma))], min(all_gamma)
#
#     def area_frac_vs_chempot_plot(self, const_u=0, increments=10,
#                                   x_is_u_ads=False):
#         """
#         Plots the change in the area contribution of
#         each facet as a function of chemical potential.
#
#         Args:
#             const_u (float): A chemical potential to hold constant if there are
#                 more than one parameters.
#             increments (bool): Number of data points between min/max or point
#                 of intersection. Defaults to 5 points.
#             x_is_u_ads (bool): Whether or not to set the adsorption chempot as
#                 a free variable (False).
#
#         Returns:
#             (Pylab): Plot of area frac on the Wulff shape
#                 for each facet vs chemical potential.
#         """
#
#         xrange = self.chempot_range if not x_is_u_ads \
#             else self.max_adsorption_chempot_range(const_u)
#         all_chempots = np.linspace(min(xrange), max(xrange),
#                                    increments)
#
#         # initialize a dictionary of lists of fractional areas for each hkl
#         hkl_area_dict = {}
#         for hkl in self.entry_dict.keys():
#             hkl_area_dict[hkl] = []
#
#         # Get plot points for each Miller index
#         for u in all_chempots:
#             u_ads = u if x_is_u_ads else const_u
#             u_ref = u if not x_is_u_ads else const_u
#             wulffshape = self.wulff_shape_from_chempot(u_ads=u_ads, u_ref=u_ref)
#
#             for hkl in wulffshape.area_fraction_dict.keys():
#                 hkl_area_dict[hkl].append(wulffshape.area_fraction_dict[hkl])
#
#         # Plot the area fraction vs chemical potential for each facet
#         plt = pretty_plot(width=8, height=7)
#         axes = plt.gca()
#
#         for i, hkl in enumerate(self.entry_dict.keys()):
#             clean_entry = list(self.entry_dict[hkl].keys())[0]
#             # Ignore any facets that never show up on the
#             # Wulff shape regardless of chemical potential
#             if all([a == 0 for a in hkl_area_dict[hkl]]):
#                 continue
#             else:
#                 plt.plot(all_chempots, hkl_area_dict[hkl],
#                          '--', color=self.color_dict[clean_entry],
#                          label=str(hkl))
#
#         # Make the figure look nice
#         plt.ylabel(r"Fractional area $A^{Wulff}_{hkl}/A^{Wulff}$")
#         self.chempot_plot_addons(plt, xrange, axes, pad=5,
#                                  rect=[-0.0, 0, 0.95, 1],
#                                  x_is_u_ads=x_is_u_ads,
#                                  ylim=[0,1])
#
#         return plt
#
#     def chempot_vs_gamma_plot_one(self, plt, clean_entry, label='',
#                                   ads_entry=None, JPERM2=False,
#                                   x_is_u_ads=False, const_u=0,
#                                   urange=None, ylim=[]):
#         """
#         Helper function to  help plot the surface energy of a
#         single surface entry as a function of chemical potential.
#
#         Args:
#             plt (Plot): A plot.
#             clean_entry (entry): Entry containing the final energy and structure
#                 of the slab whose surface energy we want to calculate
#             ads_entry (entry): An optional entry object for the adsorbed slab,
#                 defaults to None.
#             JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
#                 eV/A^2 (False)
#             x_is_u_ads (bool): Whether or not to set the adsorption chempot as
#                 a free variable (False).
#             const_u (float): A chemical potential for the fixed chempot.
#             urange (list): Chemical potential range for the free variable.
#
#         Returns:
#             (Plot): Plot of surface energy vs chemical potential for one entry.
#         """
#
#
#         if x_is_u_ads:
#             u_ref_range = [const_u, const_u]
#             u_ads_range = urange if urange else \
#                 self.max_adsorption_chempot_range(const_u)
#         else:
#             u_ref_range = urange if urange else self.chempot_range
#             u_ads_range = [const_u, const_u]
#
#         # use dashed lines for slabs that are not stoichiometric
#         # wrt bulk. Label with formula if nonstoichiometric
#         mark = '--' if clean_entry.composition.reduced_composition != \
#                        self.se_calculator.ucell_entry. \
#                            composition.reduced_composition else '-'
#
#         gamma_range = [self.se_calculator.calculate_gamma(clean_entry,
#                                                                ads_slab_entry=ads_entry,
#                                                                u_ref=u_ref_range[0],
#                                                                u_ads=u_ads_range[0]),
#                        self.se_calculator.calculate_gamma(clean_entry,
#                                                                ads_slab_entry=ads_entry,
#                                                                u_ref=u_ref_range[1],
#                                                                u_ads=u_ads_range[1])]
#
#         se_range = np.array(gamma_range) * EV_PER_ANG2_TO_JOULES_PER_M2 \
#             if JPERM2 else gamma_range
#         if not urange:
#             urange = u_ref_range if not x_is_u_ads else u_ads_range
#         color = self.color_dict[ads_entry] if ads_entry else self.color_dict[clean_entry]
#         plt.plot(urange, se_range, mark, color=color, label=label)
#
#         return plt
#
#     def chempot_vs_gamma_clean(self, miller_index=(), JPERM2=False, plt=None, ylim=[]):
#         """
#         Plots the surface energy of all facets as a function of chemical potential.
#             Each facet will be associated with its own distinct colors. Dashed lines
#             will represent stoichiometries different from that of the mpid's compound.
#
#         Args:
#             cmap (cm): A matplotlib colormap object, defaults to jet.
#             highlight_stability (bool): For each facet, there may be various
#                 terminations or stoichiometries and the relative stability of
#                 these different slabs may change with chemical potential. This
#                 dict (if provided) will highlight the stability line for a facet.
#                 The key of the dict is the Miller index we want to highlight while
#                 the value is the color of the highlight e.g. {(1,1,1): 'r'}.
#         """
#
#         plt = plt if plt else pretty_plot(width=8, height=7)
#         axes = plt.gca()
#
#         # Now we plot each individual slab surface energy
#         already_labelled = []
#         for hkl in self.entry_dict.keys():
#             if miller_index and hkl != tuple(miller_index):
#                 continue
#
#             # Plot the clean slabs
#             for entry in self.entry_dict[hkl]:
#                 # Generate a label for the type of slab
#                 label = self.create_slab_label(entry, miller_index=\
#                     () if miller_index else hkl)
#
#                 if label in already_labelled:
#                     label = None
#                 else:
#                     already_labelled.append(label)
#
#                 self.chempot_vs_gamma_plot_one(plt, entry, label=label,
#                                                JPERM2=JPERM2, x_is_u_ads=False)
#
#         # Make the figure look nice
#         plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
#             else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
#         plt = self.chempot_plot_addons(plt, self.chempot_range, axes, ylim=ylim)
#
#         return plt
#
#     def stable_u_range_dict(self, clean_only=True, const_u=0,
#                             buffer=0.1, miller_index=()):
#
#         """
#         Creates a dictionary where each entry is a key pointing to a chemical potential
#         range where the surface of that entry is stable. Does so by enumerating through
#         all possible solutions (intersect) for surface energies of a specific facet.
#
#         TODO:
#             -Make this simpler/cleaner.
#
#         Args:
#             clean_only (bool): Only get the range for clean surface entries if True.
#             const_u (float): A chemical potential for the fixed chempot.
#             buffer (float): A buffer fo r the x axis (chemical
#                 potential range). For plotting.
#             miller_index (list): Miller index for a specific facet to get a
#                 dictionary for.
#         """
#
#         all_intesects_dict = {}
#         stable_urange_dict = {}
#         max_ads_range = self.max_adsorption_chempot_range(const_u,
#                                                           buffer=buffer)
#         standard_range = max_ads_range if not clean_only else self.chempot_range
#         standard_range = sorted(standard_range)
#         x_is_u_ads = False if clean_only else True
#
#         # Get all entries for a specific facet
#         for hkl in self.entry_dict.keys():
#
#             # bool to check if at least one surface exists for a facet
#             entry_exists = False
#
#             if miller_index and hkl != tuple(miller_index):
#                 continue
#             entries_in_hkl = [entry for entry in self.entry_dict[hkl]]
#             if not clean_only:
#                 ads_entries_in_hkl = []
#                 for entry in self.entry_dict[hkl]:
#                     ads_entries_in_hkl.extend([ads_entry for ads_entry in
#                                                self.entry_dict[hkl][entry]])
#                 entries_in_hkl.extend(ads_entries_in_hkl)
#
#             # if there is only one entry for this facet, then just give it the
#             # default urange, you can't make combinations with just 1 item
#             if len(entries_in_hkl) == 1:
#                 stable_urange_dict[entries_in_hkl[0]] = standard_range
#                 continue
#             for pair in itertools.combinations(entries_in_hkl, 2):
#                 # Check if entry is adsorbed entry or clean, this is
#                 # a hassle so figure out a cleaner way to do this
#                 clean_ads_coeffs = []
#                 for p in pair:
#                     p1 = self.get_clean_ads_entry_pair(p)
#                     c = self.se_calculator.surface_energy_coefficients(p1[0],
#                                                                        ads_slab_entry=p1[1])
#                     clean_ads_coeffs.append(c)
#
#                 u, gamma = self.se_calculator.solve_2_linear_eqns(clean_ads_coeffs[0],
#                                                                   clean_ads_coeffs[1],
#                                                                   x_is_u_ads=x_is_u_ads,
#                                                                   const_u=const_u)
#
#                 if u:
#                     # If the u in a list is beyond the standard
#                     # range, set it to one of the limits
#                     if u < standard_range[0]:
#                         u_new = standard_range[0]
#                     elif u > standard_range[1]:
#                         u_new = standard_range[1]
#                     else:
#                         u_new = u
#
#                     for entry in pair:
#                         if entry not in all_intesects_dict.keys():
#                             all_intesects_dict[entry] = []
#                         all_intesects_dict[entry].append(u_new)
#
#             # Now that we have every single intersection for a given
#             # slab, find the range of u where each slab is stable
#             for entry in all_intesects_dict.keys():
#                 if entry not in stable_urange_dict.keys():
#                     stable_urange_dict[entry] = []
#                 for u_int in all_intesects_dict[entry]:
#                     stable_entries = []
#                     for i in [-1, 1]:
#                         u_ads = u_int+i*(10e-6) if x_is_u_ads else const_u
#                         u_ref = u_int+i*(10e-6) if not x_is_u_ads else const_u
#                         # return_stable_slab_entry_at_u() will only return one
#                         # entry for one gamma, since u is at an intersection,
#                         # this entry is ambiguous, we need to get the entry
#                         # slightly above and below u and check if the current
#                         # entry is any of these entries. Another inconvenience
#                         # that needs to be fixed
#                         stable_entry, gamma = \
#                             self.return_stable_slab_entry_at_u(hkl, u_ads=u_ads,
#                                                                u_ref=u_ref)
#                         stable_entries.append(stable_entry)
#                     # If this entry is stable at this u, append u
#                     if entry in stable_entries:
#                         stable_urange_dict[entry].append(u_int)
#                         entry_exists = True
#
#             # Now check for entries with only one intersection
#             # is it stable below or above u_intersect
#             for entry in stable_urange_dict.keys():
#                 # First lets check if all the u values for
#                 # an entry are the same as the standard range,
#                 # if so, just set it to the standard range
#                 if stable_urange_dict[entry]:
#                     if all([u in standard_range for u in stable_urange_dict[entry]]):
#                         stable_urange_dict[entry] = standard_range
#
#                 # If only one u, its stable from u to +-inf.
#                 # If no u, this entry is never stable
#                 if len(stable_urange_dict[entry]) == 1:
#                     u = stable_urange_dict[entry][0]
#                     for i in [-1, 1]:
#                         u_ads = u+i*(10e-6) if x_is_u_ads else const_u
#                         u_ref = u+i*(10e-6) if not x_is_u_ads else const_u
#
#                         e, se = self.return_stable_slab_entry_at_u(hkl,
#                                                                    u_ads=u_ads,
#                                                                    u_ref=u_ref)
#                         if e == entry:
#                             # If the entry stable below u, assume it is
#                             # stable at -inf, otherwise its stable at +inf
#                             u2 = standard_range[0] if i == -1 else standard_range[1]
#                             stable_urange_dict[entry].append(u2)
#                             entry_exists = True
#
#                 # now sort the ranges for each entry
#                 stable_urange_dict[entry] = sorted(stable_urange_dict[entry])
#             # Now we make sure that each facet has at least
#             # one entry, if no entries exist, this means there
#             # is not intersection, get the most stable surface
#             if not entry_exists:
#                 e, se = self.return_stable_slab_entry_at_u(hkl,
#                                                            u_ads=const_u,
#                                                            u_ref=const_u)
#                 stable_urange_dict[e] = standard_range
#
#         return stable_urange_dict
#
#     def chempot_vs_gamma_facet(self, miller_index=(), const_u=0,
#                                JPERM2=False, show_unstable=False, ylim=[]):
#         """
#         Plots the surface energy as a function of chemical potential.
#             Each facet will be associated with its own distinct colors.
#             Dashed lines will represent stoichiometries different from that
#             of the mpid's compound. Transparent lines indicates adsorption.
#
#         Args:
#             miller_index (list): Miller index for a specific facet to get a
#                 dictionary for.
#             const_u (float): A chemical potential to hold constant if there are
#                 more than one parameters.
#             plt (Plot): A plot.
#             JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
#                 eV/A^2 (False)
#             show_unstable (bool): Whether or not to show parts of the surface
#                 energy plot outside the region of stability.
#
#         Returns:
#             (Plot): Plot of surface energy vs chemical potential for all entries.
#         """
#
#         plt = pretty_plot(width=8, height=7)
#         axes = plt.gca()
#
#         # Plot wrt to adsorption chempot if any adsorption entries exist
#         x_is_u_ads = False
#         for hkl in self.entry_dict.keys():
#             for clean_entry in self.entry_dict[hkl].keys():
#                 if self.entry_dict[hkl][clean_entry]:
#                     x_is_u_ads = True
#
#         for hkl in self.entry_dict.keys():
#             if miller_index and hkl != tuple(miller_index):
#                 continue
#             if not show_unstable:
#                 clean_only = False if x_is_u_ads else True
#                 stable_u_range_dict = self.stable_u_range_dict(clean_only=clean_only,
#                                                                const_u=const_u,
#                                                                miller_index=hkl)
#
#             already_labelled = []
#             for clean_entry in self.entry_dict[hkl]:
#
#                 urange = stable_u_range_dict[clean_entry] if not show_unstable else None
#                 if urange != []:
#
#                     label = self.create_slab_label(clean_entry, miller_index=hkl)
#                     if label in already_labelled:
#                         label = None
#                     else:
#                         already_labelled.append(label)
#
#                     self.chempot_vs_gamma_plot_one(plt, clean_entry, label=label,
#                                                    JPERM2=JPERM2, x_is_u_ads=x_is_u_ads,
#                                                    const_u=const_u, urange=urange)
#
#                 for ads_entry in self.entry_dict[hkl][clean_entry]:
#                     # Plot the adsorbed slabs
#                     # Generate a label for the type of slab
#                     urange = stable_u_range_dict[ads_entry] if not show_unstable else None
#                     if urange != []:
#                         self.chempot_vs_gamma_plot_one(plt, clean_entry, JPERM2=JPERM2,
#                                                        ads_entry=ads_entry,
#                                                        const_u=const_u, urange=urange,
#                                                        x_is_u_ads=x_is_u_ads)
#
#         const_species = self.se_calculator.adsorbate_as_str  if not \
#             x_is_u_ads else self.ref_el_comp
#         if miller_index:
#             plt.title(r"%s,  $\Delta\mu_{%s}=%.2f$" %(str(miller_index),
#                                                       const_species, const_u), fontsize=20)
#
#         # Make the figure look nice
#         xrange = self.chempot_range if not x_is_u_ads \
#             else self.max_adsorption_chempot_range(const_u)
#         plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
#             else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
#         plt = self.chempot_plot_addons(plt, xrange, axes,
#                                        x_is_u_ads=x_is_u_ads, ylim=ylim)
#
#         return plt
#     def monolayer_vs_BE(self, plot_eads=False):
#         """
#         Plots the binding energy energy as a function of monolayers (ML),
#             i.e. the fractional area adsorbate density for all facets.
#         Args:
#             plot_eads (bool): Option to plot the adsorption energy (binding energy
#                 multiplied by number of adsorbates) instead.
#
#         Returns:
#             (Plot): Plot of binding energy vs monolayer for all facets.
#         """
#
#         plt = pretty_plot(width=8, height=7)
#         for hkl in self.entry_dict.keys():
#             for clean_entry in self.entry_dict[hkl].keys():
#                 if self.entry_dict[hkl][clean_entry]:
#                     monolayer = [self.se_calculator.get_monolayer(ads_entry, clean_entry)\
#                                  for ads_entry in self.entry_dict[hkl][clean_entry]]
#                     BEs = [self.se_calculator.gibbs_binding_energy(ads_entry, clean_entry,
#                                                                    eads=plot_eads) \
#                            for ads_entry in self.entry_dict[hkl][clean_entry]]
#                     # sort the binding energies and monolayers
#                     # in order to properly draw a line plot
#                     monolayer, BEs = zip(*sorted(zip(monolayer, BEs)))
#                     plt.plot(monolayer, BEs, '-o', c=self.color_dict[clean_entry], label=hkl)
#
#         plt.xlabel("%s Coverage (ML)" %(self.se_calculator.adsorbate_as_str))
#         plt.ylabel("Adsorption Energy (eV)") if plot_eads else plt.ylabel("Binding Energy (eV)")
#         plt.legend()
#         plt.tight_layout()
#
#         return plt
#
#     def chempot_plot_addons(self, plt, xrange, axes,
#                             pad=2.4, rect=[-0.047, 0, 0.84, 1],
#                             x_is_u_ads=False, ylim=[]):
#         """
#         Helper function to a chempot plot look nicer.
#
#         Args:
#             plt (Plot) Plot to add things to.
#             xrange (list): xlim parameter
#             axes(axes) Axes object from matplotlib
#             pad (float) For tight layout
#             rect (list): For tight layout
#             x_is_u_ads (bool): Whether or not to set the adsorption
#                 chempot as a free variable (False).
#             ylim (ylim parameter):
#
#         return (Plot): Modified plot with addons.
#         """
#
#         # Make the figure look nice
#         x_species = self.ref_el_comp if not \
#             x_is_u_ads else self.se_calculator.adsorbate_as_str
#         plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
#         axes.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" % (x_species))
#
#         ylim = ylim if ylim else axes.get_ylim()
#         plt.xticks(rotation=60)
#         plt.ylim(ylim)
#         xlim = axes.get_xlim()
#         plt.xlim(xlim)
#         plt.tight_layout(pad=pad, rect=rect)
#         plt.plot([xrange[0], xrange[0]], ylim, '--k')
#         plt.plot([xrange[1], xrange[1]], ylim, '--k')
#         xy = [np.mean([xrange[1]]), np.mean(ylim)]
#         plt.annotate("%s-rich" % (x_species), xy=xy,
#                      xytext=xy, rotation=90, fontsize=17)
#         xy = [np.mean([xlim[0]]), np.mean(ylim)]
#         plt.annotate("%s-poor" % (x_species), xy=xy,
#                      xytext=xy, rotation=90, fontsize=17)
#
#         return plt
#
#     def create_slab_label(self, entry, miller_index=(), ads_entry=None):
#
#         label = str(miller_index) if miller_index else ""
#         if entry.composition.reduced_composition != \
#                 self.se_calculator.ucell_entry.composition.reduced_composition:
#             label += " %s" % (entry.composition.reduced_composition)
#         if ads_entry:
#             label += r"+%s" %(self.se_calculator.adsorbate_as_str)
#         return label
#
#     def color_palette_dict(self, alpha=0.35):
#         """
#         Helper function to assign each facet a unique color using a dictionary.
#
#         Args:
#             alpha (float): Degree of transparency
#
#         return (dict): Dictionary of colors (r,g,b,a) when plotting surface
#             energy stability. The keys are individual surface entries where
#             clean surfaces have a solid color while the corresponding adsorbed
#             surface will be transparent.
#         """
#
#         color_dict = {}
#         for hkl in self.entry_dict.keys():
#             rgb_indices = [0, 1, 2]
#             color = [0, 0, 0, 1]
#             random.shuffle(rgb_indices)
#             for i, ind in enumerate(rgb_indices):
#                 if i == 2:
#                     break
#                 color[ind] = np.random.uniform(0, 1)
#
#             # Get the clean (solid) colors first
#             clean_list = np.linspace(0, 1, len(self.entry_dict[hkl]))
#             for i, clean in enumerate(self.entry_dict[hkl].keys()):
#                 c = copy.copy(color)
#                 c[rgb_indices[2]] = clean_list[i]
#                 color_dict[clean] = c
#
#                 # Now get the adsorbed (transparent) colors
#                 for ads_entry in self.entry_dict[hkl][clean]:
#                     c_ads = copy.copy(c)
#                     c_ads[3] = alpha
#                     color_dict[ads_entry] = c_ads
#
#         return color_dict
#
#     def get_clean_ads_entry_pair(self, entry):
#         """
#         Helper function that returns a pair of entries, the first is
#         the clean entry, the second is either nonetype (if the
#         initial entry is clean) or the corresponding adsorbed entry.
#
#         Args:
#             entry (ComputedStructureEntry): entry for a slab.
#         """
#
#         for hkl in self.entry_dict.keys():
#             if entry in list(self.entry_dict[hkl].keys()):
#                 # its a clean entry
#                 return [entry, None]
#             else:
#                 for clean_entry in self.entry_dict[hkl].keys():
#                     if entry in self.entry_dict[hkl][clean_entry]:
#                         return [clean_entry, entry]
#
#     def BE_vs_SE(self, plot_eads=False, const_u=0,
#                  annotate_monolayer=True, JPERM2=False):
#         """
#         For each facet, plot the clean surface energy against the most
#             stable binding energy.
#         Args:
#             plot_eads (bool): Option to plot the adsorption energy (binding
#                 energy multiplied by number of adsorbates) instead.
#             const_u (float): Ref element chemical potential as a constant.
#             annotate_monolayer (bool): Whether or not to label each data point
#                 with its monolayer (adsorbate density per unit primiitve area)
#             JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
#                 eV/A^2 (False)
#
#         Returns:
#             (Plot): Plot of clean surface energy vs binding energy for
#                 all facets.
#         """
#
#         plt = pretty_plot(width=8, height=7)
#         for hkl in self.entry_dict.keys():
#             for clean_entry in self.entry_dict[hkl].keys():
#                 if self.entry_dict[hkl][clean_entry]:
#
#                     clean_se = self.se_calculator.calculate_gamma(clean_entry,
#                                                                        u_ref=const_u)
#                     for ads_entry in self.entry_dict[hkl][clean_entry]:
#                         ml = self.se_calculator.get_monolayer(ads_entry, clean_entry)
#                         be = self.se_calculator.gibbs_binding_energy(ads_entry,
#                                                                      clean_entry,
#                                                                      eads=plot_eads)
#
#                         # Now plot the surface energy vs binding energy
#                         plt.scatter(clean_se, be)
#                         if annotate_monolayer:
#                             plt.annotate("%.2f" %(ml), xy=[clean_se, be],
#                                          xytext=[clean_se, be])
#
#         plt.xlabel(r"Surface energy ($J/m^2$)") if JPERM2 \
#             else plt.xlabel(r"Surface energy ($eV/\AA^2$)")
#         plt.ylabel("Adsorption Energy (eV)") if plot_eads \
#             else plt.ylabel("Binding Energy (eV)")
#         plt.tight_layout()
#         plt.xticks(rotation=60)
#
#         return plt
#
#     def nanoscale_stability(self):
#
#         return
#
#     def H_wrt_bulk(self):
#
#         return
#
#     def surface_phase_diagram(self, y_param, x_param, miller_index):
#         return
#
#     def wulff_shape_extrapolated_model(self):
#         return
#
#     def surface_pourbaix_diagram(self):
#
#         return
#
#     def surface_u_vs_u_phase_diagram(self):
#
#         return
#
#     def surface_p_vs_t_phase_diagram(self):
#
#         return
#
#     def broken_bond_vs_gamma(self):
#
#         return