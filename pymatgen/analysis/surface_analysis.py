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
    -Simplify the input for SurfaceEnergyPlotter such that the
        user does not need to generate a dict
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

    Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale
        stabilization of sodium oxides: Implications for Na-O2 batteries.
        Nano Letters, 14(2), 1016–1020. https://doi.org/10.1021/nl404557w

    and

    Montoya, J. H., & Persson, K. A. (2017). A high-throughput framework
        for determining adsorption energies on solid surfaces. Npj
        Computational Materials, 3(1), 14.
        https://doi.org/10.1038/s41524-017-0017-z

"""


class SlabEntry(ComputedStructureEntry):
    """
    A ComputedStructureEntry object encompassing all data relevant to a 
        slab for analyzing surface thermodynamics.
        
    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: label

        Brief description for this slab.

    .. attribute:: adsorbates
    
        List of ComputedStructureEntry for the types of adsorbates
        
    ..attribute:: clean_entry

        SlabEntry for the corresponding clean slab for an adsorbed slab
        
    ..attribute:: ads_entries_dict
    
        Dictionary where the key is the reduced composition of the 
            adsorbate entry and value is the entry itself
    """

    def __init__(self, structure, energy, miller_index, correction=0.0,
                 parameters=None, data=None, entry_id=None, label=None,
                 adsorbates=[], clean_entry=None):
        
        """
        Make a SlabEntry containing all relevant surface thermodynamics data.
        
        Args:
            structure (Slab): The primary slab associated with this entry.
            energy (float): Energy from total energy calculation
            miller_index (tuple(h, k, l)): Miller index of plane parallel 
                to surface
            correction (float): See ComputedSlabEntry
            parameters (dict): See ComputedSlabEntry
            data (dict): See ComputedSlabEntry
            entry_id (obj): See ComputedSlabEntry
            data (dict): See ComputedSlabEntry
            entry_id (str): See ComputedSlabEntry
            label (str): Any particular label for this slab, e.g. "Tasker 2", 
                "non-stoichiometric", "reconstructed"
            adsorbates ([ComputedStructureEntry]): List of reference entries  
                for the adsorbates on the slab, can be an isolated molecule 
                (e.g. O2 for O or O2 adsorption), a bulk structure (eg. fcc 
                Cu for Cu adsorption) or anything.
            clean_entry (ComputedStructureEntry): If the SlabEntry is for an 
                adsorbed slab, this is the corresponding SlabEntry for the 
                clean slab
        """

        self.miller_index = miller_index
        self.label = label
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
            ucell_entry (entry): An entry object for the bulk
            ref_entries (list: [entry]): A list of entries for each type
                of element to be used as a reservoir for nonstoichiometric
                systems. The length of this list MUST be n-1 where n is the
                number of different elements in the bulk entry. The chempot
                of the element ref_entry that is not in the list will be
                treated as a variable.

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
        """

        unit_a = self.get_unit_primitive_area
        Nsurfs = self.Nsurfs_ads_in_slab
        Nads = self.Nads_in_slab
        return Nads / (unit_a * Nsurfs)

    @property
    def Nads_in_slab(self):
        """
        Returns the TOTAL number of adsorbates in the slab on BOTH sides
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
        Returns a SlabEntry by reading in an dictionary
        """

        structure = SlabEntry.from_dict(d["structure"])
        energy = SlabEntry.from_dict(d["energy"])
        miller_index = d["miller_index"]
        label = d["label"]
        coverage = d["coverage"]
        adsorbates = d["adsorbates"]
        clean_entry = d["clean_entry"] = self.clean_entry

        return SlabEntry(structure, energy, miller_index, label=label,
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
    def cleaned_up_slab(self):
        """
        Returns a slab with the adsorbates removed
        """
        ads_strs = list(self.ads_entries_dict.keys())
        cleaned = self.structure.copy()
        cleaned.remove_species(ads_strs)
        return cleaned

    @property
    def create_slab_label(self):
        """
        Returns a label (str) for this particular slab based 
            on composition, coverage and Miller index.
        """

        if "label" in self.data.keys():
            return self.data["label"]

        label = str(self.miller_index)
        ads_strs = list(self.ads_entries_dict.keys())

        cleaned = self.cleaned_up_slab
        label += " %s" % (cleaned.composition.reduced_composition)

        if self.adsorbates:
            for ads in ads_strs:
                label += r"+%s" %(ads)
            label += r", %.3f ML" %(self.get_monolayer)
        return label


class SurfaceEnergyPlotter(object):
    """
    A class used for generating plots to analyze the thermodynamics of surfaces
        of a material. Produces stability maps of different slab configurations, 
        phases diagrams of two parameters to determine stability of configurations 
        (future release), and Wulff shapes.

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

    .. attribute:: ucell_entry

        ComputedStructureEntry of the bulk reference for this particular material.
        
    .. attribute:: ref_entries
    
        List of ComputedStructureEntries to be used for calculating chemical potential.
        
    .. attribute:: color_dict
    
        Randomly generated dictionary of colors associated with each facet.
    """

    def __init__(self, entry_dict, ucell_entry, ref_entries=[]):
        """
        Object for plotting surface energy in different ways for clean and
            adsorbed surfaces.
        Args:
            entry_dict (dict): Dictionary containing a list of entries
                for slab calculations. See attributes.
            ucell_entry (ComputedStructureEntry): ComputedStructureEntry 
                of the bulk reference for this particular material.
            ref_entries ([ComputedStructureEntries]): A list of entries for 
                each type of element to be used as a reservoir for 
                nonstoichiometric systems. The length of this list MUST be 
                n-1 where n is the number of different elements in the bulk 
                entry. The chempot of the element ref_entry that is not in 
                the list will be treated as a variable. e.g. if your ucell_entry 
                is for LiFePO4 than your ref_entries should have an entry for Li, 
                Fe, and P if you want to use the chempot of  O as the variable.
        """

        self.ucell_entry = ucell_entry
        self.ref_entries = ref_entries
        self.entry_dict = entry_dict
        self.color_dict = self.color_palette_dict()

    def set_all_variables(self, entry, u_dict, u_default):
        """
        Sets all chemical potential values and returns a dictionary where 
            the key is a sympy Symbol and the value is a float (chempot).
            
        Args:
            entry (SlabEntry): Computed structure entry of the slab 
            u_dict (Dict): Dictionary of the chemical potentials to be set as 
                constant. Note the key should be a sympy Symbol object of the 
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials 
        
        Returns:
            Dictionary of set chemical potential values
        """

        # Set up the variables
        all_u_dict = {}
        for el in entry.composition.as_dict().keys():
            if Symbol("delu_%s" %(el)) in u_dict.keys():
                all_u_dict[Symbol("delu_%s" %(el))] = u_dict[Symbol("delu_%s" %(el))]
            else:
                all_u_dict[Symbol("delu_%s" %(el))] = u_default

        return all_u_dict

    def get_stable_entry_at_u(self, miller_index, u_dict={}, u_default=0,
                              no_doped=False, no_clean=False):
        """
        Returns the entry corresponding to the most stable slab for a particular
            facet at a specific chempot. We assume that surface energy is constant
            so all free variables must be set with u_dict, otherwise they are
            assumed to be equal to u_default.

        Args:
            miller_index ((h,k,l)): The facet to find the most stable slab in
            u_dict (Dict): Dictionary of the chemical potentials to be set as 
                constant. Note the key should be a sympy Symbol object of the 
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials
            no_doped (bool): Consider stability of clean slabs only.
            no_clean (bool): Consider stability of doped slabs only.

        Returns:
            SlabEntry, surface_energy (float)
        """

        all_entries, all_gamma = [], []
        for entry in self.entry_dict[miller_index].keys():
            all_u_dict = self.set_all_variables(entry, u_dict, u_default)
            gamma = entry.surface_energy(self.ucell_entry,
                                         ref_entries=self.ref_entries)
            if not no_clean:
                all_entries.append(entry)
                all_gamma.append(gamma.subs(all_u_dict))
            for ads_entry in self.entry_dict[miller_index][entry]:
                all_u_dict = self.set_all_variables(ads_entry, u_dict, u_default)
                gamma = ads_entry.surface_energy(self.ucell_entry,
                                                 ref_entries=self.ref_entries)
                if not no_doped:
                    all_entries.append(ads_entry)
                    all_gamma.append(gamma.subs(all_u_dict))

        return all_entries[all_gamma.index(min(all_gamma))], float(min(all_gamma))

    def wulff_from_chempot(self, u_dict={}, u_default=0, symprec=1e-5,
                           no_clean=False, no_doped=False):
        """
        Method to get the Wulff shape at a specific chemical potential.
        
        Args:
            u_dict (Dict): Dictionary of the chemical potentials to be set as 
                constant. Note the key should be a sympy Symbol object of the 
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials 
            symprec (float): See WulffShape.
            no_doped (bool): Consider stability of clean slabs only.
            no_clean (bool): Consider stability of doped slabs only.

        Returns:
            (WulffShape): The WulffShape at u_ref and u_ads.
        """

        latt = SpacegroupAnalyzer(self.ucell_entry.structure). \
            get_conventional_standard_structure().lattice

        miller_list = self.entry_dict.keys()
        e_surf_list = []
        for hkl in miller_list:
            # For all configurations, calculate surface energy as a
            # function of u. Use the lowest surface energy (corresponds
            # to the most stable slab termination at that particular u)
            entry, gamma = self.get_stable_entry_at_u(hkl, u_dict=u_dict,
                                                      u_default=u_default,
                                                      no_clean=no_clean,
                                                      no_doped=no_doped)
            e_surf_list.append(gamma)

        return WulffShape(latt, miller_list, e_surf_list, symprec=symprec)

    def area_frac_vs_chempot_plot(self, ref_delu, chempot_range, u_dict={},
                                  u_default=0, increments=10):
        """
        1D plot. Plots the change in the area contribution
        of each facet as a function of chemical potential.

        Args:
            ref_delu (sympy Symbol): The free variable chempot with the format: 
                Symbol("delu_el") where el is the name of the element.
            chempot_range (list): Min/max range of chemical potential to plot along
            u_dict (Dict): Dictionary of the chemical potentials to be set as 
                constant. Note the key should be a sympy Symbol object of the 
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials 
            increments (int): Number of data points between min/max or point
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
            wulffshape = self.wulff_from_chempot(u_dict=u_dict,
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
            u_dict (Dict): Dictionary of the chemical potentials to be set as 
                constant. Note the key should be a sympy Symbol object of the 
                format: Symbol("delu_el") where el is the name of the element.

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

    def stable_u_range_dict(self, chempot_range, ref_delu, no_doped=True,
                            u_dict={}, miller_index=()):
        """
        Creates a dictionary where each entry is a key pointing to a
        chemical potential range where the surface of that entry is stable.
        Does so by enumerating through all possible solutions (intersect)
        for surface energies of a specific facet.
        
        Args:
            chempot_range ([max_chempot, min_chempot]): Range to consider the
                stability of the slabs.
            ref_delu (sympy Symbol): The range stability of each slab is based
                on the chempot range of this chempot. Should be a sympy Symbol
                object of the format: Symbol("delu_el") where el is the name of
                the element
            no_doped (bool): Consider stability of clean slabs only.
            no_clean (bool): Consider stability of doped slabs only.
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the 
                format: Symbol("delu_el") where el is the name of the element.
            miller_index (list): Miller index for a specific facet to get a
                dictionary for.
        """

        all_intesects_dict, stable_urange_dict = {}, {}

        # Get all entries for a specific facet
        for hkl in self.entry_dict.keys():

            # Skip this facet if this is not the facet we want
            if miller_index and hkl != tuple(miller_index):
                continue

            entries_in_hkl = [clean for clean in self.entry_dict[hkl].keys()]
            if not no_doped:
                for entry in self.entry_dict[hkl].keys():
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
                # I'm assuming ref_delu was not set in u_dict,
                # so the solution should be for ref_delu
                solution = self.get_surface_equilibrium(pair, u_dict=u_dict)

                # Check if this solution is stable
                if not solution:
                    continue
                new_u_dict = u_dict.copy()
                new_u_dict[ref_delu] = solution[ref_delu]
                stable_entry, gamma = self.get_stable_entry_at_u(hkl, new_u_dict)
                if stable_entry not in pair:
                    continue

                # Now check if the solution is within the chempot range
                if not (chempot_range[0] <= solution[ref_delu] <= chempot_range[1]):
                    continue

                for entry in pair:
                    stable_urange_dict[entry].append(solution[ref_delu])

            # Now check if all entries have 2 chempot values. If only
            # one, we need to set the other value as either the upper
            # limit or lower limit of the user provided chempot_range
            new_u_dict = u_dict.copy()
            for u in chempot_range:
                new_u_dict[ref_delu] = u
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

    def chempot_vs_gamma_plot_one(self, plt, entry, ref_delu, chempot_range,
                                  u_dict={}, u_default=0, label='', JPERM2=False):
        """
        Helper function to  help plot the surface energy of a
        single SlabEntry as a function of chemical potential.

        Args:
            plt (Plot): A plot.
            entry (SlabEntry): Entry of the slab whose surface energy we want
                to plot
            ref_delu (sympy Symbol): The range stability of each slab is based
                on the chempot range of this chempot. Should be a sympy Symbol
                object of the format: Symbol("delu_el") where el is the name of
                the element
            chempot_range ([max_chempot, min_chempot]): Range to consider the
                stability of the slabs.
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials
            label (str): Label of the slab for the legend.
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)

        Returns:
            (Plot): Plot of surface energy vs chemical potential for one entry.
        """

        # use dashed lines for slabs that are not stoichiometric
        # wrt bulk. Label with formula if nonstoichiometric
        ucell_comp = self.ucell_entry.composition.reduced_composition
        if entry.adsorbates:
            s = entry.cleaned_up_slab
            clean_comp = s.composition.reduced_composition
        else:
            clean_comp = entry.composition.reduced_composition
        mark = '--' if ucell_comp != clean_comp else '-'

        u_dict = self.set_all_variables(entry, u_dict, u_default)
        u_dict[ref_delu] = chempot_range[0]
        gamma_min = entry.surface_energy(self.ucell_entry,
                                         ref_entries=self.ref_entries).subs(u_dict)
        u_dict[ref_delu] = chempot_range[1]
        gamma_max = entry.surface_energy(self.ucell_entry,
                                         ref_entries=self.ref_entries).subs(u_dict)
        gamma_range = [gamma_min, gamma_max]

        se_range = np.array(gamma_range) * EV_PER_ANG2_TO_JOULES_PER_M2 \
            if JPERM2 else gamma_range

        plt.plot(chempot_range, se_range, mark,
                 color=self.color_dict[entry], label=label)

        return plt

    def chempot_vs_gamma(self, ref_delu, chempot_range, miller_index=(),
                         u_dict={}, u_default=0, JPERM2=False,
                         show_unstable=False, ylim=[],
                         no_clean=False, no_doped=False):
        """
        Plots the surface energy as a function of chemical potential.
            Each facet will be associated with its own distinct colors.
            Dashed lines will represent stoichiometries different from that
            of the mpid's compound. Transparent lines indicates adsorption.

        Args:
            ref_delu (sympy Symbol): The range stability of each slab is based
                on the chempot range of this chempot. Should be a sympy Symbol
                object of the format: Symbol("delu_el") where el is the name of
                the element
            chempot_range ([max_chempot, min_chempot]): Range to consider the
                stability of the slabs.
            miller_index (list): Miller index for a specific facet to get a
                dictionary for.
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)
            show_unstable (bool): Whether or not to show parts of the surface
                energy plot outside the region of stability.
            ylim ([ymax, ymin]): Range of y axis
            no_doped (bool): Whether to plot for the clean slabs only.
            no_clean (bool): Whether to plot for the doped slabs only.

        Returns:
            (Plot): Plot of surface energy vs chempot for all entries.
        """

        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        for hkl in self.entry_dict.keys():
            if miller_index and hkl != tuple(miller_index):
                continue
            # Get the chempot range of each surface if we only
            # want to show the region where each slab is stable
            if not show_unstable:
                stable_u_range_dict = self.stable_u_range_dict(chempot_range, ref_delu,
                                                               no_doped=no_doped,
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
                    if not no_clean:
                        print("urange", urange)
                        plt = self.chempot_vs_gamma_plot_one(plt, clean_entry, ref_delu,
                                                             urange, u_dict=u_dict,
                                                             u_default=u_default,
                                                             label=label, JPERM2=JPERM2)
                if not no_doped:
                    for ads_entry in self.entry_dict[hkl][clean_entry]:
                        # Plot the adsorbed slabs
                        # Generate a label for the type of slab
                        urange = stable_u_range_dict[ads_entry] \
                            if not show_unstable else chempot_range
                        if urange != []:

                            plt = self.chempot_vs_gamma_plot_one(plt, ads_entry,
                                                                 ref_delu, urange,
                                                                 u_dict=u_dict,
                                                                 u_default=u_default,
                                                                 label=label,
                                                                 JPERM2=JPERM2)

        # Make the figure look nice
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, chempot_range, str(ref_delu).split("_")[1],
                                       axes, ylim=ylim)

        return plt

    def monolayer_vs_BE(self, plot_eads=False):
        """
        Plots the binding energy energy as a function of monolayers (ML), i.e.
            the fractional area adsorbate density for all facets. For each
            facet at a specific monlayer, only plot the lowest binding energy.

        Args:
            plot_eads (bool): Option to plot the adsorption energy (binding
                 energy multiplied by number of adsorbates) instead.

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
            ref_el (str): Element of the referenced chempot.
            axes(axes) Axes object from matplotlib
            pad (float) For tight layout
            rect (list): For tight layout
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
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials
            plot_eads (bool): Option to plot the adsorption energy (binding
                energy multiplied by number of adsorbates) instead.
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
                all_u_dict = self.set_all_variables(clean_entry, u_dict, u_default)
                if self.entry_dict[hkl][clean_entry]:

                    clean_se = clean_entry.surface_energy(self.ucell_entry,
                                                          ref_entries=self.ref_entries)
                    se = clean_se.subs(all_u_dict)
                    for ads_entry in self.entry_dict[hkl][clean_entry]:
                        ml = ads_entry.get_monolayer
                        be = ads_entry.gibbs_binding_energy(eads=plot_eads)

                        # Now plot the surface energy vs binding energy
                        plt.scatter(se, be)
                        if annotate_monolayer:
                            plt.annotate("%.2f" %(ml), xy=[se, be],
                                         xytext=[se, be])

        plt.xlabel(r"Surface energy ($J/m^2$)") if JPERM2 \
            else plt.xlabel(r"Surface energy ($eV/\AA^2$)")
        plt.ylabel("Adsorption Energy (eV)") if plot_eads \
            else plt.ylabel("Binding Energy (eV)")
        plt.tight_layout()
        plt.xticks(rotation=60)

        return plt

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


class NanoscaleStability(object):
    """
    A class for analyzing the stability of nanoparticles of different
        polymorphs with respect to size. The Wulff shape will be the
        model for the nanoparticle. Stability will be determined by
        an energetic competition between the weighted surface energy
        (surface energy of the Wulff shape) and the bulk energy. A
        future release will include a 2D phase diagram (e.g. wrt size
        vs chempot for adsorbed or nonstoichiometric surfaces). Based
        on the following work:

        Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale
            stabilization of sodium oxides: Implications for Na-O2
            batteries. Nano Letters, 14(2), 1016–1020.
            https://doi.org/10.1021/nl404557w

    .. attribute:: se_analyzers

        List of SurfaceEnergyPlotter objects. Each item corresponds to a
            different polymorph.

    .. attribute:: symprec

        See WulffShape.
    """

    def __init__(self, se_analyzers, symprec=1e-5):

        """
        Analyzes the nanoscale stability of different polymorphs.
        """

        self.se_analyzers = se_analyzers
        self.symprec = symprec

    def solve_equilibrium_point(self, analyzer1, analyzer2,
                                u_dict={}, u_default=0):
        """
        Gives the radial size of two particles where equilibrium is reached
            between both particles. NOTE: the solution here is not the same
            as the solution visualized in the plot because solving for r
            requires that both the total surface area and volume of the
            particles are functions of r.

        Args:
            analyzer1 (SurfaceEnergyPlotter): Analyzer associated with the
                first polymorph
            analyzer2 (SurfaceEnergyPlotter): Analyzer associated with the
                second polymorph
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials

        Returns:
            Particle radius in nm
        """

        # Set up
        s1 = analyzer1.ucell_entry.structure
        s2 = analyzer2.ucell_entry.structure
        E1 = analyzer1.ucell_entry.energy_per_atom
        E2 = analyzer2.ucell_entry.energy_per_atom
        wulff1 = analyzer1.wulff_from_chempot(u_dict=u_dict,
                                              u_default=u_default,
                                              symprec=self.symprec)
        wulff2 = analyzer2.wulff_from_chempot(u_dict=u_dict,
                                              u_default=u_default,
                                              symprec=self.symprec)

        # Now calculate r
        delta_gamma = wulff1.weighted_surface_energy - wulff2.weighted_surface_energy
        delta_E = (len(s1) / s1.lattice.volume) * E1 - (len(s2) / s2.lattice.volume) * E2

        return ((-3 * delta_gamma) / (delta_E))/10

    def wulff_gform_and_r(self, wulffshape, bulk_entry,
                          r, from_sphere_area=False):
        """
        Calculates the formation energy of the particle with arbitrary radius r.

        Args:
            wulffshape (WulffShape): Initial, unscaled WulffShape
            bulk_entry (ComputedStructureEntry): Entry of the corresponding bulk.
            r (float (Ang)): Arbitrary effective radius of the WulffShape
            from_sphere_area (bool): There are two ways to calculate the bulk
                formation energy. Either by treating the volume and thus surface
                area of the particle as a perfect sphere, or as a Wulff shape.

        Returns:
            particle formation energy (float in keV), effective radius (float in nm)
        """

        # Set up
        miller_se_dict = wulffshape.miller_energy_dict
        new_wulff = self.scaled_wulff(wulffshape, r)
        new_wulff_area = new_wulff.miller_area_dict

        # calculate surface energy of the particle
        if not from_sphere_area:
            # By approximating the particle as a Wulff shape
            tot_wulff_se = 0
            for hkl in new_wulff_area.keys():
                tot_wulff_se += miller_se_dict[hkl] * new_wulff_area[hkl]
            Ebulk = self.bulk_gform(bulk_entry)*new_wulff.volume
            new_r = new_wulff.effective_radius
        else:
            # By approximating the particle as a perfect sphere
            sphere_sa = 4 * np.pi * r ** 2
            tot_wulff_se = wulffshape.weighted_surface_energy * sphere_sa
            wulff_v = (4/3)*np.pi*r**3
            Ebulk = self.bulk_gform(bulk_entry)*wulff_v
            new_r = r

        return (Ebulk + tot_wulff_se)/1000, new_r/10

    def bulk_gform(self, bulk_entry):
        """
        Returns the formation energy of the bulk
        Args:
            bulk_entry (ComputedStructureEntry): Entry of the corresponding bulk.
        """

        ucell = bulk_entry.structure
        N_per_ucell = len(ucell) / ucell.lattice.volume
        Gform = N_per_ucell * bulk_entry.energy_per_atom

        return Gform

    def scaled_wulff(self, wulffshape, r):
        """
        Scales the Wulff shape with an effective radius r. Note that the resulting
            Wulff does not neccesarily have the same effective radius as the one
            provided. The Wulff shape is scaled by its surface energies where first
            the surface energies are scale by the minimum surface energy and then
            multiplied by the given effective radius.

        Args:
            wulffshape (WulffShape): Initial, unscaled WulffShape
            r (float): Arbitrary effective radius of the WulffShape

        Returns:
            WulffShape (scaled by r)
        """

        miller_list = wulffshape.miller_energy_dict.keys()
        # Normalize the magnitude of the facet normal vectors
        # of the Wulff shape by the minimum surface energy.
        se_list = np.array(list(wulffshape.miller_energy_dict.values()))
        # Scale the magnitudes by r
        scaled_se = (se_list / min(se_list)) * r

        return WulffShape(wulffshape.lattice, miller_list,
                          scaled_se, symprec=self.symprec)

    def plot_one_stability_map(self, analyzer, max_r, u_dict={}, label="",
                               increments=50, u_default=0, plt=None,
                               from_sphere_area=False):
        """
        Returns the plot of the formation energy of a particle against its
            effect radius

        Args:
            analyzer (SurfaceEnergyPlotter): Analyzer associated with the
                first polymorph
            max_r (float): The maximum radius of the particle to plot up to.
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            label (str): Label of the plot for legend
            increments (int): Number of plot points
            u_default (float): Default value for all unset chemical potentials
            plt (pylab): Plot
            from_sphere_area (bool): There are two ways to calculate the bulk
                formation energy. Either by treating the volume and thus surface
                area of the particle as a perfect sphere, or as a Wulff shape.
        """

        plt = plt if plt else pretty_plot(width=8, height=7)

        wulffshape = analyzer.wulff_from_chempot(u_dict=u_dict,
                                                 u_default=u_default,
                                                 symprec=self.symprec)

        gform_list, r_list = [], []
        for r in np.linspace(1e-6, max_r, increments):
            gform, r = self.wulff_gform_and_r(wulffshape,
                                              analyzer.ucell_entry, r,
                                              from_sphere_area=from_sphere_area)
            gform_list.append(gform)
            r_list.append(r)
        plt.xlabel("Particle radius (nm)")
        plt.ylabel(r"$\Delta \bar{G}_{form}$ (keV)")

        plt.plot(r_list, gform_list, label=label)

        return plt

    def plot_all_stability_map(self, max_r, increments=50, u_dict={},
                               u_default=0, plt=None, labels=[],
                               from_sphere_area=False):
        """
        Returns the plot of the formation energy of a particles
            of different polymorphs against its effect radius

        Args:
            max_r (float): The maximum radius of the particle to plot up to.
            increments (int): Number of plot points
            u_dict (Dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            u_default (float): Default value for all unset chemical potentials
            plt (pylab): Plot
            labels (list): List of labels for each plot, corresponds to the
                list of se_analyzers
            from_sphere_area (bool): There are two ways to calculate the bulk
                formation energy. Either by treating the volume and thus surface
                area of the particle as a perfect sphere, or as a Wulff shape.
        """

        plt = plt if plt else pretty_plot(width=8, height=7)

        for i, analyzer in enumerate(self.se_analyzers):
            label = labels[i] if labels else ""
            plt = self.plot_one_stability_map(analyzer, max_r, u_dict,
                                              label=label, plt=plt,
                                              increments=increments,
                                              u_default=u_default,
                                              from_sphere_area=from_sphere_area)

        return plt


# class GetChempotRange(object):
#     def __init__(self, entry):
#         self.entry = entry
#
#
# class SlabEntryGenerator(object):
#     def __init__(self, entry):
#         self.entry = entry
