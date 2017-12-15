# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
TODO:
    -Still assumes individual elements have their own chempots
        in a molecular adsorbate instead of considering a single
        chempot for a single molecular adsorbate. E.g. for an OH
        adsorbate, the surface energy is a function of delu_O and
        delu_H instead of delu_OH
    -Need a method to automatically get chempot range when
        dealing with non-stoichiometric slabs
    -Simplify or clean the stable_u_range_dict() method
        for the SurfaceEnergyPlotter class
    -Annotations for which configuration is on a plot
    -Simplify the input for SurfaceEnergyPlotter such that the
        user does not need to generate a dict
    -u_dict and ref_u should not be parameters at the same time? 
        save the user from himself.
"""

from __future__ import division, unicode_literals

import numpy as np
import itertools
import warnings
import random, copy
from sympy import Symbol, Number
from sympy.solvers import linsolve

from pymatgen.core.composition import Composition
from pymatgen import Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.wulff import WulffShape
from pymatgen.util.plotting import pretty_plot

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


class SlabEntry(ComputedStructureEntry):
    """
    An object encompassing all data relevant to slab for surface calculations.

    Args:
        entry (ComputedEntry/ComputedStructureEntry): An
            entry object
    """

    def __init__(self, structure, energy, miller_index, correction=0.0,
                 parameters=None, data=None, entry_id=None, label=None,
                 coverage=None, adsorbates=[], clean_entry=None):

        self.miller_index = miller_index
        self.label = label
        self.coverage = coverage
        self.adsorbates = adsorbates
        self.clean_entry = clean_entry
        self.ads_entries_dict = {str(list(ads.composition.as_dict().keys())[0]): \
                                ads for ads in adsorbates}

        super(SlabEntry, self).__init__(
            structure, energy, correction=correction,
            parameters=parameters, data=data, entry_id=entry_id)

    def as_dict(self):
        """
        Returns dict which contains Slab Entry data.
        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d["structure"] = self.structure
        d["energy"] = self.energy
        d["miller_index"] = self.miller_index
        d["label"] = self.label
        d["coverage"] = self.coverage
        d["adsorbates"] = self.adsorbates
        d["clean_entry"] = self.clean_entry

        return d

    def gibbs_binding_energy(self, eads=False):
        """
        Returns the adsorption energy or Gibb's binding energy
            of an adsorbate on a surface
        Args:
            eads (bool): Whether to calculate the adsorption energy
                (True) or the binding energy (False) which is just
                adsorption energy normalized by number of adsorbates.
        """

        n = self.get_unit_primitive_area
        Nads = self.Nads_in_slab

        BE = (self.energy - n * self.clean_entry.energy) / Nads - \
             sum([ads.energy_per_atom for ads in self.adsorbates])
        return BE * Nads if eads else BE

    def surface_energy(self, ucell_entry, ref_entries=[]):
        """
        Calculates the surface energy of this SlabEntry.
        Args:
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

        Returns (Add (Sympy class)): Surface energy
        """

        # Set up
        gamma = (Symbol("E_surf")- Symbol("Ebulk"))/(2*Symbol("A"))
        ucell_comp = ucell_entry.composition
        ucell_reduced_comp = ucell_comp.reduced_composition
        ref_entries_dict = {str(list(ref.composition.as_dict().keys())[0]): \
                                ref for ref in ref_entries}
        ref_entries_dict.update(self.ads_entries_dict)

        # Calculate Gibbs free energy of the bulk per unit formula
        gbulk = ucell_entry.energy / \
                ucell_comp.get_integer_formula_and_factor()[1]

        # First we get the contribution to the bulk energy
        # from each element with an existing ref_entry.
        bulk_energy, gbulk_eqn = 0, 0
        for el, ref in ref_entries_dict.items():
            N, delu = self.composition.as_dict()[el], Symbol("delu_"+str(el))
            if el in ucell_comp.as_dict().keys():
                gbulk_eqn += ucell_reduced_comp[el] * (delu + ref.energy_per_atom)
            bulk_energy += N * (Symbol("delu_"+el) + ref.energy_per_atom)

        # Next, we add the contribution to the bulk energy from
        # the variable element (the element without a ref_entry),
        # as a function of the other elements
        for ref_el in ucell_comp.as_dict().keys():
            if str(ref_el) not in ref_entries_dict.keys():
                delu = Symbol("delu_" + str(ref_el))
                break
        refEperA = (gbulk-gbulk_eqn)/ucell_reduced_comp.as_dict()[ref_el] - delu
        bulk_energy += self.composition.as_dict()[ref_el] * (delu + refEperA)

        return gamma.subs({Symbol("E_surf"): self.energy, Symbol("Ebulk"): bulk_energy,
                           Symbol("A"): self.surface_area})

    @property
    def get_unit_primitive_area(self):
        """
        Returns the surface area of the adsorbed system per
        unit area of the primitive slab system.
        """

        A_ads = self.surface_area
        A_clean = self.clean_entry.surface_area
        n = (A_ads / A_clean)
        return n

    @property
    def get_monolayer(self):
        """
        Returns the primitive unit surface area density of the
            adsorbate.
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
            clean_slab_entry (entry): The entry of the clean slab
        """

        unit_a = self.get_unit_primitive_area
        Nsurfs = self.Nsurfs_ads_in_slab
        Nads = self.Nads_in_slab
        return Nads / (unit_a * Nsurfs)

    @property
    def Nads_in_slab(self):
        """
        Returns the TOTAL number of adsorbates in the slab on BOTH sides
        Args:
            ads_slab_entry (entry): The entry of the adsorbed slab
        """
        return sum([self.composition.as_dict()[a] for a \
                    in self.ads_entries_dict.keys()])

    @property
    def Nsurfs_ads_in_slab(self):
        """
        Returns the TOTAL number of adsorbed surfaces in the slab
        """

        struct = self.structure
        weights = [s.species_and_occu.weight for s in struct]
        center_of_mass = np.average(struct.frac_coords,
                                    weights=weights, axis=0)

        Nsurfs = 0
        # Are there adsorbates on top surface?
        if any([site.species_string in self.ads_entries_dict.keys() for \
                site in struct if site.frac_coords[2] > center_of_mass[2]]):
            Nsurfs += 1
        # Are there adsorbates on bottom surface?
        if any([site.species_string in self.ads_entries_dict.keys() for \
                site in struct if site.frac_coords[2] < center_of_mass[2]]):
            Nsurfs += 1

        return Nsurfs

    @classmethod
    def from_dict(cls, d):
        """
        Returns a SlabEntry by reading in an slab
        """

        entry = SlabEntry.from_dict(d["entry"])
        miller_index = d["miller_index"]
        label = d["label"]
        coverage = d["coverage"]
        adsorbates = d["adsorbates"]
        clean_entry = d["clean_entry"] = self.clean_entry

        return SlabEntry(entry, miller_index, label=label,
                         coverage=coverage, adsorbates=adsorbates,
                         clean_entry=clean_entry)

    @property
    def surface_area(self):
        """
        Calculates the surface area of the slab
        """
        m = self.structure.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))

    @property
    def create_slab_label(self):

        if "label" in self.data.keys():
            return self.data["label"]

        label = str(self.miller_index)
        ads_strs = list(self.ads_entries_dict.keys())

        cleaned = self.structure.copy()
        cleaned.remove_species(ads_strs)
        label += " %s" % (cleaned.composition.reduced_composition)

        if self.adsorbates:
            for ads in ads_strs:
                label += r"+%s" %(ads)
            label += r", %.3f ML" %(self.get_monolayer)
        return label


class SurfaceEnergyPlotter(object):
    """
    A class used for generating plots to analyze the thermodynamics of surfaces
        of a material by taking in a SurfaceEnergyCalculator object. Produces
        stability maps of different slab configurations, phases diagrams of two
        parameters to determine stability of configurations, and Wulff shapes.

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

    ..attribute:: color_dict

        Dictionary of colors (r,g,b,a) when plotting surface energy stability. The
            keys are individual surface entries where clean surfaces have a solid color
            while the corresponding adsorbed surface will be transparent.

    .. attribute:: ref_el_comp

        Composition of the ref_element. All chemical potentials can be written in terms
            of the range of chemical potential of this element which will be used to
            calculate surface energy.

    """

    def __init__(self, entry_dict, ucell_entry, ref_entries=[]):
        """
        Object for plotting surface energy in different ways for clean and
            adsorbed surfaces.
        Args:
            entry_dict (dict): Dictionary containing a list of entries
                for slab calculations. See attributes.
            surface_energy_calculator (SurfaceEnergyCalculator):
                Object for calculating thermodynamic quantities related to surfaces
            chempot_range (list): Max and min range for the chemical potential for
                the ref_element.
        """

        self.ucell_entry = ucell_entry
        self.ref_entries = ref_entries
        self.entry_dict = entry_dict
        self.color_dict = self.color_palette_dict()

    def set_all_variables(self, entry, u_dict, u_default):

        # Set up the variables
        all_u_dict = {}
        for el in entry.composition.as_dict().keys():
            if Symbol("delu_%s" %(el)) in u_dict.keys():
                all_u_dict[Symbol("delu_%s" %(el))] = u_dict[Symbol("delu_%s" %(el))]
            else:
                all_u_dict[Symbol("delu_%s" %(el))] = u_default

        return all_u_dict

    def get_stable_entry_at_u(self, miller_index, u_dict={}, u_default=0):
        """
        Returns the entry corresponding to the most stable slab for a particular
            facet at a specific chempot. We assume that surface energy is constant
            so all free variables must be set with u_dict, otherwise they are
            assumed to be equal to u_default.

        Args:
            miller_index ((h,k,l)): The facet to find the most stable slab in
            u_dict (dict): Dictionary of chemical potentials to keep constant.
                The key is the element and the value is the chemical potential.
            u_default (float): Chempot value if a chempot has no assigned value.

        Returns:
            SlabEntry, surface_energy (float)
        """

        all_entries, all_gamma = [], []
        for entry in self.entry_dict[miller_index].keys():
            all_u_dict = self.set_all_variables(entry, u_dict, u_default)
            gamma = entry.surface_energy(self.ucell_entry,
                                         ref_entries=self.ref_entries)
            all_entries.append(entry)
            all_gamma.append(gamma.subs(all_u_dict))
            for ads_entry in self.entry_dict[miller_index][entry]:
                all_u_dict = self.set_all_variables(ads_entry, u_dict, u_default)
                gamma = ads_entry.surface_energy(self.ucell_entry,
                                                 ref_entries=self.ref_entries)
                all_entries.append(ads_entry)
                all_gamma.append(gamma.subs(all_u_dict))

        return all_entries[all_gamma.index(min(all_gamma))], float(min(all_gamma))

    def wulff_shape_from_chempot(self, u_dict={}, u_default=0, symprec=1e-5):
        """
        Method to get the Wulff shape at a specific chemical potential.
        Args:
            u_dict (dict): Dictionary of chemical potentials to keep constant.
                The key is the element and the value is the chemical potential.
            symprec (float): See WulffShape.

        Returns:
            (WulffShape): The WulffShape at u_ref and u_ads.
        """

        latt = SpacegroupAnalyzer(self.ucell_entry.structure). \
            get_conventional_standard_structure().lattice

        miller_list = self.entry_dict.keys()
        e_surf_list = []
        for hkl in miller_list:
            # At each possible configuration, we calculate surface energy as a
            # function of u and take the lowest surface energy (corresponds to
            # the most stable slab termination at that particular u)
            entry, gamma = self.get_stable_entry_at_u(hkl, u_dict=u_dict,
                                                              u_default=u_default)
            e_surf_list.append(gamma)

        return WulffShape(latt, miller_list, e_surf_list, symprec=symprec)

    def area_frac_vs_chempot_plot(self, ref_delu, chempot_range, u_dict={},
                                  u_default=0, increments=10):
        """
        1D plot. Plots the change in the area contribution
        of each facet as a function of chemical potential.

        Args:
            ref_delu (str): The free variable chempot.
            chempot_range (list): Min/max range of chemical potential to plot along
            u_dict (dict): Dictionary of chemical potentials to keep constant.
                The key is the element and the value is the chemical potential.
            u_default (float): Chempot value if a chempot has no assigned value.
            increments (bool): Number of data points between min/max or point
                of intersection. Defaults to 10 points.

        Returns:
            (Pylab): Plot of area frac on the Wulff shape
                for each facet vs chemical potential.
        """

        all_chempots = np.linspace(min(chempot_range), max(chempot_range),
                                   increments)

        # initialize a dictionary of lists of fractional areas for each hkl
        hkl_area_dict = {}
        for hkl in self.entry_dict.keys():
            hkl_area_dict[hkl] = []

        # Get plot points for each Miller index
        for u in all_chempots:
            u_dict[ref_delu] = u
            wulffshape = self.wulff_shape_from_chempot(u_dict=u_dict,
                                                       u_default=u_default)

            for hkl in wulffshape.area_fraction_dict.keys():
                hkl_area_dict[hkl].append(wulffshape.area_fraction_dict[hkl])

        # Plot the area fraction vs chemical potential for each facet
        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        for i, hkl in enumerate(self.entry_dict.keys()):
            clean_entry = list(self.entry_dict[hkl].keys())[0]
            # Ignore any facets that never show up on the
            # Wulff shape regardless of chemical potential
            if all([a == 0 for a in hkl_area_dict[hkl]]):
                continue
            else:
                plt.plot(all_chempots, hkl_area_dict[hkl],
                         '--', color=self.color_dict[clean_entry],
                         label=str(hkl))

        # Make the figure look nice
        plt.ylabel(r"Fractional area $A^{Wulff}_{hkl}/A^{Wulff}$")
        self.chempot_plot_addons(plt, chempot_range, str(ref_delu).split("_")[1],
                                 axes, rect=[-0.0, 0, 0.95, 1], pad=5, ylim=[0,1])


        return plt

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
                The key is the element and the value is the chemical potential

        Returns:
            (array): Array containing a solution to x equations with x
                variables (x-1 chemical potential and 1 surface energy)
        """

        # Generate all possible coefficients
        all_parameters = []
        all_eqns = []
        for slab_entry in slab_entries:
            se = slab_entry.surface_energy(self.ucell_entry,
                                           ref_entries=self.ref_entries)

            # remove the free chempots we wish to keep constant and
            # set the equation to 0 (subtract gamma from both sides)
            all_eqns.append(se.subs(u_dict) - Symbol("gamma"))
            all_parameters.extend([p for p in list(se.free_symbols)
                                   if p not in all_parameters])
        all_parameters.append(Symbol("gamma"))
        # Now solve the system of linear eqns to find the chempot
        # where the slabs are at equilibrium with each other

        soln = linsolve(all_eqns, all_parameters)
        if not soln:
            warnings.warn("No solution")
            return soln
        return {p: list(soln)[0][i] for i, p in enumerate(all_parameters)}

    def stable_u_range_dict(self, chempot_range, ref_u, clean_only=True,
                            u_dict={}, miller_index=()):
        """
        Creates a dictionary where each entry is a key pointing to a
        chemical potential range where the surface of that entry is stable.
        Does so by enumerating through all possible solutions (intersect)
        for surface energies of a specific facet.

        TODO:
            -Make this simpler/cleaner.

        Args:
            clean_only (bool): Only get the range for clean slab entries if True.
            const_u (float): A chemical potential for the fixed chempot.
            buffer (float): A buffer fo r the x axis (chemical
                potential range). For plotting.
            miller_index (list): Miller index for a specific facet to get a
                dictionary for.
        """

        all_intesects_dict, stable_urange_dict = {}, {}

        # Get all entries for a specific facet
        for hkl in self.entry_dict.keys():

            # Skip this facet if this is not the facet we want
            if miller_index and hkl != tuple(miller_index):
                continue

            entries_in_hkl = [clean for clean in self.entry_dict[hkl]]
            if not clean_only:
                for entry in self.entry_dict[hkl]:
                    entries_in_hkl.extend([ads_entry for ads_entry in
                                           self.entry_dict[hkl][entry]])

            for entry in entries_in_hkl:
                stable_urange_dict[entry] = []

            # if there is only one entry for this facet, then just give it the
            # default urange, you can't make combinations with just 1 item
            if len(entries_in_hkl) == 1:
                stable_urange_dict[entries_in_hkl[0]] = chempot_range
                continue

            for pair in itertools.combinations(entries_in_hkl, 2):
                # I'm assuming ref_u was not set in u_dict,
                # so the solution should be for ref_u
                solution = self.get_surface_equilibrium(pair, u_dict=u_dict)

                # Check if this solution is stable
                if not solution:
                    continue
                new_u_dict = u_dict.copy()
                new_u_dict[ref_u] = solution[ref_u]
                stable_entry, gamma = self.get_stable_entry_at_u(hkl, new_u_dict)
                if stable_entry not in pair:
                    continue

                # Now check if the solution is within the chempot range
                if not (chempot_range[0] <= solution[ref_u] <= chempot_range[1]):
                    continue

                for entry in pair:
                    stable_urange_dict[entry].append(solution[ref_u])

            # Now check if all entries have 2 chempot values. If only
            # one, we need to set the other value as either the upper
            # limit or lower limit of the user provided chempot_range
            new_u_dict = u_dict.copy()
            for u in chempot_range:
                print(u)
                new_u_dict[ref_u] = u
                entry, gamma = self.get_stable_entry_at_u(hkl, u_dict=new_u_dict)
                stable_urange_dict[entry].append(u)

        # sort the chempot ranges for each facet
        for entry in stable_urange_dict.keys():
            stable_urange_dict[entry] = sorted(stable_urange_dict[entry])

        return stable_urange_dict

    def color_palette_dict(self, alpha=0.35):
        """
        Helper function to assign each facet a unique color using a dictionary.

        Args:
            alpha (float): Degree of transparency

        return (dict): Dictionary of colors (r,g,b,a) when plotting surface
            energy stability. The keys are individual surface entries where
            clean surfaces have a solid color while the corresponding adsorbed
            surface will be transparent.
        """

        color_dict = {}
        for hkl in self.entry_dict.keys():
            rgb_indices = [0, 1, 2]
            color = [0, 0, 0, 1]
            random.shuffle(rgb_indices)
            for i, ind in enumerate(rgb_indices):
                if i == 2:
                    break
                color[ind] = np.random.uniform(0, 1)

            # Get the clean (solid) colors first
            clean_list = np.linspace(0, 1, len(self.entry_dict[hkl]))
            for i, clean in enumerate(self.entry_dict[hkl].keys()):
                c = copy.copy(color)
                c[rgb_indices[2]] = clean_list[i]
                color_dict[clean] = c

                # Now get the adsorbed (transparent) colors
                for ads_entry in self.entry_dict[hkl][clean]:
                    c_ads = copy.copy(c)
                    c_ads[3] = alpha
                    color_dict[ads_entry] = c_ads

        return color_dict

    def chempot_vs_gamma_plot_one(self, plt, entry, ref_u, chempot_range,
                                  u_dict={}, u_default=0, label='', JPERM2=False):
        """
        Helper function to  help plot the surface energy of a
        single surface entry as a function of chemical potential.

        Args:
            plt (Plot): A plot.
            clean_entry (entry): Entry containing the final energy and structure
                of the slab whose surface energy we want to calculate
            ads_entry (entry): An optional entry object for the adsorbed slab,
                defaults to None.
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)
            x_is_u_ads (bool): Whether or not to set the adsorption chempot as
                a free variable (False).
            const_u (float): A chemical potential for the fixed chempot.
            urange (list): Chemical potential range for the free variable.

        Returns:
            (Plot): Plot of surface energy vs chemical potential for one entry.
        """

        # use dashed lines for slabs that are not stoichiometric
        # wrt bulk. Label with formula if nonstoichiometric
        ucell_comp = self.ucell_entry.composition.reduced_composition.as_dict()
        is_stoich = []
        for el in ucell_comp.keys():
            if entry.adsorbates:
                s = entry.structure
                s.remove_species(entry.ads_entries_dict.keys())
                is_stoich.append(ucell_comp[el] == \
                                 s.composition.reduced_composition.as_dict()[el])
            else:
                is_stoich.append(ucell_comp[el] == \
                                 entry.composition.reduced_composition.as_dict()[el])
        mark = '--' if not all(is_stoich) else '-'

        u_dict = self.set_all_variables(entry, u_dict, u_default)
        u_dict[ref_u] = chempot_range[0]
        gamma_min = entry.surface_energy(self.ucell_entry,
                                         ref_entries=self.ref_entries).subs(u_dict)
        u_dict[ref_u] = chempot_range[1]
        gamma_max = entry.surface_energy(self.ucell_entry,
                                         ref_entries=self.ref_entries).subs(u_dict)
        gamma_range = [gamma_min, gamma_max]

        se_range = np.array(gamma_range) * EV_PER_ANG2_TO_JOULES_PER_M2 \
            if JPERM2 else gamma_range

        plt.plot(chempot_range, se_range, mark,
                 color=self.color_dict[entry], label=label)

        return plt

    def chempot_vs_gamma(self, ref_u, chempot_range, miller_index=(),
                         u_dict={}, u_default=0, JPERM2=False,
                         show_unstable=False, ylim=[], clean_only=False):
        """
        Plots the surface energy as a function of chemical potential.
            Each facet will be associated with its own distinct colors.
            Dashed lines will represent stoichiometries different from that
            of the mpid's compound. Transparent lines indicates adsorption.

        Args:
            miller_index (list): Miller index for a specific facet to get a
                dictionary for.
            const_u (float): A chemical potential to hold constant if there are
                more than one parameters.
            plt (Plot): A plot.
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)
            show_unstable (bool): Whether or not to show parts of the surface
                energy plot outside the region of stability.

        Returns:
            (Plot): Plot of surface energy vs chemical potential for all entries.
        """

        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        for hkl in self.entry_dict.keys():
            if miller_index and hkl != tuple(miller_index):
                continue
            if not show_unstable:
                stable_u_range_dict = self.stable_u_range_dict(chempot_range, ref_u,
                                                               clean_only=clean_only,
                                                               u_dict=u_dict,
                                                               miller_index=hkl)

            already_labelled = []
            for clean_entry in self.entry_dict[hkl]:

                urange = stable_u_range_dict[clean_entry] if \
                    not show_unstable else chempot_range
                # Don't plot if the slab is unstable, plot if it is.
                if urange != []:

                    label = clean_entry.label
                    if label in already_labelled:
                        label = None
                    else:
                        already_labelled.append(label)

                    plt = self.chempot_vs_gamma_plot_one(plt, clean_entry, ref_u,
                                                         urange, u_dict=u_dict,
                                                         u_default=u_default,
                                                         label=label, JPERM2=JPERM2)
                if not clean_only:
                    for ads_entry in self.entry_dict[hkl][clean_entry]:
                        # Plot the adsorbed slabs
                        # Generate a label for the type of slab
                        urange = stable_u_range_dict[ads_entry] \
                            if not show_unstable else chempot_range
                        if urange != []:
                            plt = self.chempot_vs_gamma_plot_one(plt, ads_entry,
                                                                 ref_u, urange,
                                                                 u_dict=u_dict,
                                                                 u_default=u_default,
                                                                 label=label,
                                                                 JPERM2=JPERM2)

        # Make the figure look nice
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, chempot_range, str(ref_u).split("_")[1],
                                       axes, ylim=ylim)

        return plt

    def monolayer_vs_BE(self, plot_eads=False):
        """
        Plots the binding energy energy as a function of monolayers (ML),
            i.e. the fractional area adsorbate density for all facets.
            For each facet at a specific monlayer, only plot the lowest
            binding energy.
        Args:
            plot_eads (bool): Option to plot the adsorption energy (binding energy
                multiplied by number of adsorbates) instead.

        Returns:
            (Plot): Plot of binding energy vs monolayer for all facets.
        """

        plt = pretty_plot(width=8, height=7)
        for hkl in self.entry_dict.keys():
            ml_be_dict = {}
            for clean_entry in self.entry_dict[hkl].keys():
                if self.entry_dict[hkl][clean_entry]:
                    for ads_entry in self.entry_dict[hkl][clean_entry]:
                        if ads_entry.get_monolayer not in ml_be_dict.keys():
                            ml_be_dict[ads_entry.get_monolayer] = 1000
                        be = ads_entry.gibbs_binding_energy(eads=plot_eads)
                        if be < ml_be_dict[ads_entry.get_monolayer]:
                            ml_be_dict[ads_entry.get_monolayer] = be
            # sort the binding energies and monolayers
            # in order to properly draw a line plot
            vals = sorted(ml_be_dict.items())
            monolayers, BEs = zip(*vals)
            plt.plot(monolayers, BEs, '-o',
                     c=self.color_dict[clean_entry], label=hkl)

        adsorbates = tuple(ads_entry.ads_entries_dict.keys())
        plt.xlabel(" %s"*len(adsorbates) %adsorbates + " Coverage (ML)")
        plt.ylabel("Adsorption Energy (eV)") if plot_eads \
            else plt.ylabel("Binding Energy (eV)")
        plt.legend()
        plt.tight_layout()

        return plt

    def chempot_plot_addons(self, plt, xrange, ref_el, axes, pad=2.4,
                            rect=[-0.047, 0, 0.84, 1], ylim=[]):

        """
        Helper function to a chempot plot look nicer.

        Args:
            plt (Plot) Plot to add things to.
            xrange (list): xlim parameter
            axes(axes) Axes object from matplotlib
            pad (float) For tight layout
            rect (list): For tight layout
            x_is_u_ads (bool): Whether or not to set the adsorption
                chempot as a free variable (False).
            ylim (ylim parameter):

        return (Plot): Modified plot with addons.
        """

        # Make the figure look nice
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        axes.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" % (ref_el))

        ylim = ylim if ylim else axes.get_ylim()
        plt.xticks(rotation=60)
        plt.ylim(ylim)
        xlim = axes.get_xlim()
        plt.xlim(xlim)
        plt.tight_layout(pad=pad, rect=rect)
        plt.plot([xrange[0], xrange[0]], ylim, '--k')
        plt.plot([xrange[1], xrange[1]], ylim, '--k')
        xy = [np.mean([xrange[1]]), np.mean(ylim)]
        plt.annotate("%s-rich" % (ref_el), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)
        xy = [np.mean([xlim[0]]), np.mean(ylim)]
        plt.annotate("%s-poor" % (ref_el), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)

        return plt

    def BE_vs_clean_SE(self, u_dict, u_default=0, plot_eads=False,
                       annotate_monolayer=True, JPERM2=False):
        """
        For each facet, plot the clean surface energy against the most
            stable binding energy.
        Args:
            plot_eads (bool): Option to plot the adsorption energy (binding
                energy multiplied by number of adsorbates) instead.
            const_u (float): Ref element chemical potential as a constant.
            annotate_monolayer (bool): Whether or not to label each data point
                with its monolayer (adsorbate density per unit primiitve area)
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)

        Returns:
            (Plot): Plot of clean surface energy vs binding energy for
                all facets.
        """

        plt = pretty_plot(width=8, height=7)
        for hkl in self.entry_dict.keys():
            for clean_entry in self.entry_dict[hkl].keys():
                if self.entry_dict[hkl][clean_entry]:

                    clean_se = self.se_calculator.calculate_gamma(clean_entry,
                                                                       u_ref=const_u)
                    for ads_entry in self.entry_dict[hkl][clean_entry]:
                        ml = self.se_calculator.get_monolayer(ads_entry, clean_entry)
                        be = self.se_calculator.gibbs_binding_energy(ads_entry,
                                                                     clean_entry,
                                                                     eads=plot_eads)

                        # Now plot the surface energy vs binding energy
                        plt.scatter(clean_se, be)
                        if annotate_monolayer:
                            plt.annotate("%.2f" %(ml), xy=[clean_se, be],
                                         xytext=[clean_se, be])

        plt.xlabel(r"Surface energy ($J/m^2$)") if JPERM2 \
            else plt.xlabel(r"Surface energy ($eV/\AA^2$)")
        plt.ylabel("Adsorption Energy (eV)") if plot_eads \
            else plt.ylabel("Binding Energy (eV)")
        plt.tight_layout()
        plt.xticks(rotation=60)

        return plt

    # def nanoscale_stability(self):
    #
    #     return
    #
    # def H_wrt_bulk(self):
    #
    #     return
    #
    # def surface_phase_diagram(self, y_param, x_param, miller_index):
    #     return
    #
    # def wulff_shape_extrapolated_model(self):
    #     return
    #
    # def surface_pourbaix_diagram(self):
    #
    #     return
    #
    # def surface_u_vs_u_phase_diagram(self):
    #
    #     return
    #
    # def surface_p_vs_t_phase_diagram(self):
    #
    #     return
    #
    # def broken_bond_vs_gamma(self):
    #
    #     return


# class GetChempotRange(object):
#     def __init__(self, entry):
#         self.entry = entry
#
#
# class SlabEntryGenerator(object):
#     def __init__(self, entry):
#         self.entry = entry
