# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

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
    Nano Letters, 14(2), 1016-1020. https://doi.org/10.1021/nl404557w

    and

    Montoya, J. H., & Persson, K. A. (2017). A high-throughput framework
        for determining adsorption energies on solid surfaces. Npj
        Computational Materials, 3(1), 14.
        https://doi.org/10.1038/s41524-017-0017-z

Todo:
- Still assumes individual elements have their own chempots
    in a molecular adsorbate instead of considering a single
    chempot for a single molecular adsorbate. E.g. for an OH
    adsorbate, the surface energy is a function of delu_O and
    delu_H instead of delu_OH
- Need a method to automatically get chempot range when
    dealing with non-stoichiometric slabs
- Simplify the input for SurfaceEnergyPlotter such that the
    user does not need to generate a dict
"""

from __future__ import annotations

import copy
import itertools
import random
import warnings

import numpy as np
from sympy import Symbol
from sympy.solvers import linsolve, solve

from pymatgen.analysis.wulff import WulffShape
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_slab_regions
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.vasp.outputs import Locpot, Outcar, Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.plotting import pretty_plot

EV_PER_ANG2_TO_JOULES_PER_M2 = 16.0217656

__author__ = "Richard Tran"
__credits__ = "Joseph Montoya, Xianguo Li"


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

    def __init__(
        self,
        structure,
        energy,
        miller_index,
        correction=0.0,
        parameters=None,
        data=None,
        entry_id=None,
        label=None,
        adsorbates=None,
        clean_entry=None,
        marker=None,
        color=None,
    ):
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
            marker (str): Custom marker for gamma plots ("--" and "-" are typical)
            color (str or rgba): Custom color for gamma plots
        """
        self.miller_index = miller_index
        self.label = label
        self.adsorbates = adsorbates if adsorbates else []
        self.clean_entry = clean_entry
        self.ads_entries_dict = {str(list(ads.composition.as_dict())[0]): ads for ads in self.adsorbates}
        self.mark = marker
        self.color = color

        super().__init__(
            structure,
            energy,
            correction=correction,
            parameters=parameters,
            data=data,
            entry_id=entry_id,
        )

    def as_dict(self):
        """
        Returns dict which contains Slab Entry data.
        """
        d = {"@module": type(self).__module__, "@class": type(self).__name__}
        d["structure"] = self.structure
        d["energy"] = self.energy
        d["miller_index"] = self.miller_index
        d["label"] = self.label
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

        BE = (self.energy - n * self.clean_entry.energy) / Nads - sum(ads.energy_per_atom for ads in self.adsorbates)
        return BE * Nads if eads else BE

    def surface_energy(self, ucell_entry, ref_entries=None):
        """
        Calculates the surface energy of this SlabEntry.

        Args:
            ucell_entry (entry): An entry object for the bulk
            ref_entries (list: [entry]): A list of entries for each type
                of element to be used as a reservoir for non-stoichiometric
                systems. The length of this list MUST be n-1 where n is the
                number of different elements in the bulk entry. The chempot
                of the element ref_entry that is not in the list will be
                treated as a variable.

        Returns (Add (Sympy class)): Surface energy
        """
        # Set up
        ref_entries = ref_entries if ref_entries else []

        # Check if appropriate ref_entries are present if the slab is non-stoichiometric
        # TODO: There should be a way to identify which specific species are
        # non-stoichiometric relative to the others in systems with more than 2 species
        slab_comp = self.composition.as_dict()
        ucell_entry_comp = ucell_entry.composition.reduced_composition.as_dict()
        slab_clean_comp = Composition({el: slab_comp[el] for el in ucell_entry_comp})
        if slab_clean_comp.reduced_composition != ucell_entry.composition.reduced_composition:
            list_els = [list(entry.composition.as_dict())[0] for entry in ref_entries]
            if not any(el in list_els for el in ucell_entry.composition.as_dict()):
                warnings.warn("Elemental references missing for the non-dopant species.")

        gamma = (Symbol("E_surf") - Symbol("Ebulk")) / (2 * Symbol("A"))
        ucell_comp = ucell_entry.composition
        ucell_reduced_comp = ucell_comp.reduced_composition
        ref_entries_dict = {str(list(ref.composition.as_dict())[0]): ref for ref in ref_entries}
        ref_entries_dict.update(self.ads_entries_dict)

        # Calculate Gibbs free energy of the bulk per unit formula
        gbulk = ucell_entry.energy / ucell_comp.get_integer_formula_and_factor()[1]

        # First we get the contribution to the bulk energy
        # from each element with an existing ref_entry.
        bulk_energy, gbulk_eqn = 0, 0
        for el, ref in ref_entries_dict.items():
            N, delu = self.composition.as_dict()[el], Symbol("delu_" + str(el))
            if el in ucell_comp.as_dict():
                gbulk_eqn += ucell_reduced_comp[el] * (delu + ref.energy_per_atom)
            bulk_energy += N * (Symbol("delu_" + el) + ref.energy_per_atom)

        # Next, we add the contribution to the bulk energy from
        # the variable element (the element without a ref_entry),
        # as a function of the other elements
        for ref_el in ucell_comp.as_dict():
            if str(ref_el) not in ref_entries_dict:
                break
        refEperA = (gbulk - gbulk_eqn) / ucell_reduced_comp.as_dict()[ref_el]
        bulk_energy += self.composition.as_dict()[ref_el] * refEperA
        se = gamma.subs(
            {
                Symbol("E_surf"): self.energy,
                Symbol("Ebulk"): bulk_energy,
                Symbol("A"): self.surface_area,
            }
        )

        return float(se) if type(se).__name__ == "Float" else se

    @property
    def get_unit_primitive_area(self):
        """
        Returns the surface area of the adsorbed system per
        unit area of the primitive slab system.
        """
        A_ads = self.surface_area
        A_clean = self.clean_entry.surface_area
        n = A_ads / A_clean
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
        return sum(self.composition.as_dict()[a] for a in self.ads_entries_dict)

    @property
    def Nsurfs_ads_in_slab(self):
        """
        Returns the TOTAL number of adsorbed surfaces in the slab
        """
        struct = self.structure
        weights = [s.species.weight for s in struct]
        center_of_mass = np.average(struct.frac_coords, weights=weights, axis=0)

        Nsurfs = 0
        # Are there adsorbates on top surface?
        if any(
            site.species_string in self.ads_entries_dict for site in struct if site.frac_coords[2] > center_of_mass[2]
        ):
            Nsurfs += 1
        # Are there adsorbates on bottom surface?
        if any(
            site.species_string in self.ads_entries_dict for site in struct if site.frac_coords[2] < center_of_mass[2]
        ):
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
        adsorbates = d["adsorbates"]
        clean_entry = d["clean_entry"]

        return cls(
            structure,
            energy,
            miller_index,
            label=label,
            adsorbates=adsorbates,
            clean_entry=clean_entry,
        )

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
        ads_strs = list(self.ads_entries_dict)
        cleaned = self.structure.copy()
        cleaned.remove_species(ads_strs)
        return cleaned

    @property
    def create_slab_label(self):
        """
        Returns a label (str) for this particular slab based on composition, coverage and Miller index.
        """
        if "label" in self.data:
            return self.data["label"]

        label = str(self.miller_index)
        ads_strs = list(self.ads_entries_dict)

        cleaned = self.cleaned_up_slab
        label += f" {cleaned.composition.reduced_composition}"

        if self.adsorbates:
            for ads in ads_strs:
                label += f"+{ads}"
            label += f", {self.get_monolayer:.3f} ML"
        return label

    @staticmethod
    def from_computed_structure_entry(entry, miller_index, label=None, adsorbates=None, clean_entry=None, **kwargs):
        """
        Returns SlabEntry from a ComputedStructureEntry
        """
        return SlabEntry(
            entry.structure,
            entry.energy,
            miller_index,
            label=label,
            adsorbates=adsorbates,
            clean_entry=clean_entry,
            **kwargs,
        )


class SurfaceEnergyPlotter:
    """
    A class used for generating plots to analyze the thermodynamics of surfaces
        of a material. Produces stability maps of different slab configurations,
        phases diagrams of two parameters to determine stability of configurations
        (future release), and Wulff shapes.

    .. attribute:: all_slab_entries

        Either a list of SlabEntry objects (note for a list, the SlabEntry must
            have the adsorbates and clean_entry parameter pulgged in) or a Nested
            dictionary containing a list of entries for slab calculations as
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
            entries are present (i.e. if all_slab_entries[(h,k,l)][clean_entry1]), we
            consider adsorption in all plots and analysis for this particular facet.

    ..attribute:: color_dict

        Dictionary of colors (r,g,b,a) when plotting surface energy stability. The
            keys are individual surface entries where clean surfaces have a solid
            color while the corresponding adsorbed surface will be transparent.

    .. attribute:: ucell_entry

        ComputedStructureEntry of the bulk reference for this particular material.

    .. attribute:: ref_entries

        List of ComputedStructureEntries to be used for calculating chemical potential.

    .. attribute:: color_dict

        Randomly generated dictionary of colors associated with each facet.
    """

    def __init__(self, all_slab_entries, ucell_entry, ref_entries=None):
        """
        Object for plotting surface energy in different ways for clean and
            adsorbed surfaces.

        Args:
            all_slab_entries (dict or list): Dictionary or list containing
                all entries for slab calculations. See attributes.
            ucell_entry (ComputedStructureEntry): ComputedStructureEntry
                of the bulk reference for this particular material.
            ref_entries ([ComputedStructureEntries]): A list of entries for
                each type of element to be used as a reservoir for
                non-stoichiometric systems. The length of this list MUST be
                n-1 where n is the number of different elements in the bulk
                entry. The bulk energy term in the grand surface potential can
                be defined by a summation of the chemical potentials for each
                element in the system. As the bulk energy is already provided,
                one can solve for one of the chemical potentials as a function
                of the other chemical potetinals and bulk energy. i.e. there
                are n-1 variables (chempots). e.g. if your ucell_entry is for
                LiFePO4 than your ref_entries should have an entry for Li, Fe,
                and P if you want to use the chempot of O as the variable.
        """
        self.ucell_entry = ucell_entry
        self.ref_entries = ref_entries
        self.all_slab_entries = (
            all_slab_entries if type(all_slab_entries).__name__ == "dict" else entry_dict_from_list(all_slab_entries)
        )
        self.color_dict = self.color_palette_dict()

        se_dict, as_coeffs_dict = {}, {}
        for hkl in self.all_slab_entries:
            for clean in self.all_slab_entries[hkl]:
                se = clean.surface_energy(self.ucell_entry, ref_entries=self.ref_entries)
                if type(se).__name__ == "float":
                    se_dict[clean] = se
                    as_coeffs_dict[clean] = {1: se}
                else:
                    se_dict[clean] = se
                    as_coeffs_dict[clean] = se.as_coefficients_dict()
                for dope in self.all_slab_entries[hkl][clean]:
                    se = dope.surface_energy(self.ucell_entry, ref_entries=self.ref_entries)
                    if type(se).__name__ == "float":
                        se_dict[dope] = se
                        as_coeffs_dict[dope] = {1: se}
                    else:
                        se_dict[dope] = se
                        as_coeffs_dict[dope] = se.as_coefficients_dict()
        self.surfe_dict = se_dict
        self.as_coeffs_dict = as_coeffs_dict

        list_of_chempots = []
        for v in self.as_coeffs_dict.values():
            if type(v).__name__ == "float":
                continue
            for du in v:
                if du not in list_of_chempots:
                    list_of_chempots.append(du)
        self.list_of_chempots = list_of_chempots

    def get_stable_entry_at_u(
        self,
        miller_index,
        delu_dict=None,
        delu_default=0,
        no_doped=False,
        no_clean=False,
    ):
        """
        Returns the entry corresponding to the most stable slab for a particular
            facet at a specific chempot. We assume that surface energy is constant
            so all free variables must be set with delu_dict, otherwise they are
            assumed to be equal to delu_default.

        Args:
            miller_index ((h,k,l)): The facet to find the most stable slab in
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            no_doped (bool): Consider stability of clean slabs only.
            no_clean (bool): Consider stability of doped slabs only.

        Returns:
            SlabEntry, surface_energy (float)
        """
        all_delu_dict = self.set_all_variables(delu_dict, delu_default)

        def get_coeffs(e):
            coeffs = []
            for du in all_delu_dict:
                if type(self.as_coeffs_dict[e]).__name__ == "float":
                    coeffs.append(self.as_coeffs_dict[e])
                elif du in self.as_coeffs_dict[e]:
                    coeffs.append(self.as_coeffs_dict[e][du])
                else:
                    coeffs.append(0)
            return np.array(coeffs)

        all_entries, all_coeffs = [], []
        for entry in self.all_slab_entries[miller_index]:
            if not no_clean:
                all_entries.append(entry)
                all_coeffs.append(get_coeffs(entry))
            if not no_doped:
                for ads_entry in self.all_slab_entries[miller_index][entry]:
                    all_entries.append(ads_entry)
                    all_coeffs.append(get_coeffs(ads_entry))

        du_vals = np.array(list(all_delu_dict.values()))
        all_gamma = list(np.dot(all_coeffs, du_vals.T))

        return all_entries[all_gamma.index(min(all_gamma))], float(min(all_gamma))

    def wulff_from_chempot(
        self,
        delu_dict=None,
        delu_default=0,
        symprec=1e-5,
        no_clean=False,
        no_doped=False,
    ):
        """
        Method to get the Wulff shape at a specific chemical potential.

        Args:
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            symprec (float): See WulffShape.
            no_doped (bool): Consider stability of clean slabs only.
            no_clean (bool): Consider stability of doped slabs only.

        Returns:
            (WulffShape): The WulffShape at u_ref and u_ads.
        """
        latt = SpacegroupAnalyzer(self.ucell_entry.structure).get_conventional_standard_structure().lattice

        miller_list = list(self.all_slab_entries)
        e_surf_list = []
        for hkl in miller_list:
            # For all configurations, calculate surface energy as a
            # function of u. Use the lowest surface energy (corresponds
            # to the most stable slab termination at that particular u)
            gamma = self.get_stable_entry_at_u(
                hkl,
                delu_dict=delu_dict,
                delu_default=delu_default,
                no_clean=no_clean,
                no_doped=no_doped,
            )[1]
            e_surf_list.append(gamma)

        return WulffShape(latt, miller_list, e_surf_list, symprec=symprec)

    def area_frac_vs_chempot_plot(
        self,
        ref_delu,
        chempot_range,
        delu_dict=None,
        delu_default=0,
        increments=10,
        no_clean=False,
        no_doped=False,
    ):
        """
        1D plot. Plots the change in the area contribution
        of each facet as a function of chemical potential.

        Args:
            ref_delu (sympy Symbol): The free variable chempot with the format:
                Symbol("delu_el") where el is the name of the element.
            chempot_range (list): Min/max range of chemical potential to plot along
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            increments (int): Number of data points between min/max or point
                of intersection. Defaults to 10 points.

        Returns:
            (Pylab): Plot of area frac on the Wulff shape
                for each facet vs chemical potential.
        """
        delu_dict = delu_dict or {}
        chempot_range = sorted(chempot_range)
        all_chempots = np.linspace(min(chempot_range), max(chempot_range), increments)

        # initialize a dictionary of lists of fractional areas for each hkl
        hkl_area_dict = {}
        for hkl in self.all_slab_entries:
            hkl_area_dict[hkl] = []

        # Get plot points for each Miller index
        for u in all_chempots:
            delu_dict[ref_delu] = u
            wulffshape = self.wulff_from_chempot(
                delu_dict=delu_dict,
                no_clean=no_clean,
                no_doped=no_doped,
                delu_default=delu_default,
            )

            for hkl in wulffshape.area_fraction_dict:
                hkl_area_dict[hkl].append(wulffshape.area_fraction_dict[hkl])

        # Plot the area fraction vs chemical potential for each facet
        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        for hkl in self.all_slab_entries:
            clean_entry = list(self.all_slab_entries[hkl])[0]
            # Ignore any facets that never show up on the
            # Wulff shape regardless of chemical potential
            if all(a == 0 for a in hkl_area_dict[hkl]):
                continue
            plt.plot(
                all_chempots,
                hkl_area_dict[hkl],
                "--",
                color=self.color_dict[clean_entry],
                label=str(hkl),
            )

        # Make the figure look nice
        plt.ylabel(r"Fractional area $A^{Wulff}_{hkl}/A^{Wulff}$")
        self.chempot_plot_addons(
            plt,
            chempot_range,
            str(ref_delu).split("_")[1],
            axes,
            rect=[-0.0, 0, 0.95, 1],
            pad=5,
            ylim=[0, 1],
        )

        return plt

    def get_surface_equilibrium(self, slab_entries, delu_dict=None):
        """
        Takes in a list of SlabEntries and calculates the chemical potentials
            at which all slabs in the list coexists simultaneously. Useful for
            building surface phase diagrams. Note that to solve for x equations
            (x slab_entries), there must be x free variables (chemical potentials).
            Adjust delu_dict as need be to get the correct number of free variables.

        Args:
            slab_entries (array): The coefficients of the first equation
            delu_dict (dict): Dictionary of the chemical potentials to be set as
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
            se = self.surfe_dict[slab_entry]

            # remove the free chempots we wish to keep constant and
            # set the equation to 0 (subtract gamma from both sides)
            if type(se).__name__ == "float":
                all_eqns.append(se - Symbol("gamma"))
            else:
                se = sub_chempots(se, delu_dict) if delu_dict else se
                all_eqns.append(se - Symbol("gamma"))
                all_parameters.extend([p for p in list(se.free_symbols) if p not in all_parameters])

        all_parameters.append(Symbol("gamma"))
        # Now solve the system of linear eqns to find the chempot
        # where the slabs are at equilibrium with each other

        soln = linsolve(all_eqns, all_parameters)
        if not soln:
            warnings.warn("No solution")
            return soln
        return {p: list(soln)[0][i] for i, p in enumerate(all_parameters)}

    def stable_u_range_dict(
        self,
        chempot_range,
        ref_delu,
        no_doped=True,
        no_clean=False,
        delu_dict=None,
        miller_index=(),
        dmu_at_0=False,
        return_se_dict=False,
    ):
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
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            miller_index (list): Miller index for a specific facet to get a
                dictionary for.
            dmu_at_0 (bool): If True, if the surface energies corresponding to
                the chemical potential range is between a negative and positive
                value, the value is a list of three chemical potentials with the
                one in the center corresponding a surface energy of 0. Uselful
                in identifying unphysical ranges of surface energies and their
                chemical potential range.
            return_se_dict (bool): Whether or not to return the corresponding
                dictionary of surface energies
        """
        if delu_dict is None:
            delu_dict = {}
        chempot_range = sorted(chempot_range)
        stable_urange_dict, se_dict = {}, {}

        # Get all entries for a specific facet
        for hkl in self.all_slab_entries:
            entries_in_hkl = []
            # Skip this facet if this is not the facet we want
            if miller_index and hkl != tuple(miller_index):
                continue
            if not no_clean:
                entries_in_hkl.extend(self.all_slab_entries[hkl])
            if not no_doped:
                for entry in self.all_slab_entries[hkl]:
                    entries_in_hkl.extend(self.all_slab_entries[hkl][entry])

            for entry in entries_in_hkl:
                stable_urange_dict[entry] = []
                se_dict[entry] = []
            # if there is only one entry for this facet, then just give it the
            # default urange, you can't make combinations with just 1 item
            if len(entries_in_hkl) == 1:
                stable_urange_dict[entries_in_hkl[0]] = chempot_range
                u1, u2 = delu_dict.copy(), delu_dict.copy()
                u1[ref_delu], u2[ref_delu] = chempot_range[0], chempot_range[1]
                se = self.as_coeffs_dict[entries_in_hkl[0]]
                se_dict[entries_in_hkl[0]] = [
                    sub_chempots(se, u1),
                    sub_chempots(se, u2),
                ]
                continue

            for pair in itertools.combinations(entries_in_hkl, 2):
                # I'm assuming ref_delu was not set in delu_dict,
                # so the solution should be for ref_delu
                solution = self.get_surface_equilibrium(pair, delu_dict=delu_dict)

                # Check if this solution is stable
                if not solution:
                    continue
                new_delu_dict = delu_dict.copy()
                new_delu_dict[ref_delu] = solution[ref_delu]
                stable_entry, gamma = self.get_stable_entry_at_u(
                    hkl, new_delu_dict, no_doped=no_doped, no_clean=no_clean
                )
                if stable_entry not in pair:
                    continue

                # Now check if the solution is within the chempot range
                if not chempot_range[0] <= solution[ref_delu] <= chempot_range[1]:
                    continue

                for entry in pair:
                    stable_urange_dict[entry].append(solution[ref_delu])
                    se_dict[entry].append(gamma)

            # Now check if all entries have 2 chempot values. If only
            # one, we need to set the other value as either the upper
            # limit or lower limit of the user provided chempot_range
            new_delu_dict = delu_dict.copy()
            for u in chempot_range:
                new_delu_dict[ref_delu] = u
                entry, gamma = self.get_stable_entry_at_u(
                    hkl, delu_dict=new_delu_dict, no_doped=no_doped, no_clean=no_clean
                )
                stable_urange_dict[entry].append(u)
                se_dict[entry].append(gamma)

        if dmu_at_0:
            for entry, v in se_dict.items():
                # if se are of opposite sign, determine chempot when se=0.
                # Useful for finding a chempot range where se is unphysical
                if not stable_urange_dict[entry]:
                    continue
                if v[0] * v[1] < 0:
                    # solve for gamma=0
                    se = self.as_coeffs_dict[entry]
                    v.append(0)
                    stable_urange_dict[entry].append(solve(sub_chempots(se, delu_dict), ref_delu)[0])

        # sort the chempot ranges for each facet
        for entry, v in stable_urange_dict.items():
            se_dict[entry] = [se for i, se in sorted(zip(v, se_dict[entry]))]
            stable_urange_dict[entry] = sorted(v)

        if return_se_dict:
            return stable_urange_dict, se_dict
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
        for hkl in self.all_slab_entries:
            rgb_indices = [0, 1, 2]
            color = [0, 0, 0, 1]
            random.shuffle(rgb_indices)
            for i, ind in enumerate(rgb_indices):
                if i == 2:
                    break
                color[ind] = np.random.uniform(0, 1)

            # Get the clean (solid) colors first
            clean_list = np.linspace(0, 1, len(self.all_slab_entries[hkl]))
            for i, clean in enumerate(self.all_slab_entries[hkl]):
                c = copy.copy(color)
                c[rgb_indices[2]] = clean_list[i]
                color_dict[clean] = c

                # Now get the adsorbed (transparent) colors
                for ads_entry in self.all_slab_entries[hkl][clean]:
                    c_ads = copy.copy(c)
                    c_ads[3] = alpha
                    color_dict[ads_entry] = c_ads

        return color_dict

    def chempot_vs_gamma_plot_one(
        self,
        plt,
        entry,
        ref_delu,
        chempot_range,
        delu_dict=None,
        delu_default=0,
        label="",
        JPERM2=False,
    ):
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
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            label (str): Label of the slab for the legend.
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)

        Returns:
            (Plot): Plot of surface energy vs chemical potential for one entry.
        """
        if delu_dict is None:
            delu_dict = {}
        chempot_range = sorted(chempot_range)

        # use dashed lines for slabs that are not stoichiometric
        # wrt bulk. Label with formula if non-stoichiometric
        ucell_comp = self.ucell_entry.composition.reduced_composition
        if entry.adsorbates:
            s = entry.cleaned_up_slab
            clean_comp = s.composition.reduced_composition
        else:
            clean_comp = entry.composition.reduced_composition

        mark = "--" if ucell_comp != clean_comp else "-"

        delu_dict = self.set_all_variables(delu_dict, delu_default)
        delu_dict[ref_delu] = chempot_range[0]
        gamma_min = self.as_coeffs_dict[entry]
        gamma_min = gamma_min if type(gamma_min).__name__ == "float" else sub_chempots(gamma_min, delu_dict)
        delu_dict[ref_delu] = chempot_range[1]
        gamma_max = self.as_coeffs_dict[entry]
        gamma_max = gamma_max if type(gamma_max).__name__ == "float" else sub_chempots(gamma_max, delu_dict)
        gamma_range = [gamma_min, gamma_max]

        se_range = np.array(gamma_range) * EV_PER_ANG2_TO_JOULES_PER_M2 if JPERM2 else gamma_range

        mark = entry.mark if entry.mark else mark
        c = entry.color if entry.color else self.color_dict[entry]
        plt.plot(chempot_range, se_range, mark, color=c, label=label)

        return plt

    def chempot_vs_gamma(
        self,
        ref_delu,
        chempot_range,
        miller_index=(),
        delu_dict=None,
        delu_default=0,
        JPERM2=False,
        show_unstable=False,
        ylim=None,
        plt=None,
        no_clean=False,
        no_doped=False,
        use_entry_labels=False,
        no_label=False,
    ):
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
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            JPERM2 (bool): Whether to plot surface energy in /m^2 (True) or
                eV/A^2 (False)
            show_unstable (bool): Whether or not to show parts of the surface
                energy plot outside the region of stability.
            ylim ([ymax, ymin]): Range of y axis
            no_doped (bool): Whether to plot for the clean slabs only.
            no_clean (bool): Whether to plot for the doped slabs only.
            use_entry_labels (bool): If True, will label each slab configuration
                according to their given label in the SlabEntry object.
            no_label (bool): Option to turn off labels.

        Returns:
            (Plot): Plot of surface energy vs chempot for all entries.
        """
        if delu_dict is None:
            delu_dict = {}
        chempot_range = sorted(chempot_range)

        plt = plt if plt else pretty_plot(width=8, height=7)
        axes = plt.gca()

        for hkl in self.all_slab_entries:
            if miller_index and hkl != tuple(miller_index):
                continue
            # Get the chempot range of each surface if we only
            # want to show the region where each slab is stable
            if not show_unstable:
                stable_u_range_dict = self.stable_u_range_dict(
                    chempot_range, ref_delu, no_doped=no_doped, delu_dict=delu_dict, miller_index=hkl
                )

            already_labelled = []
            label = ""
            for clean_entry in self.all_slab_entries[hkl]:
                urange = stable_u_range_dict[clean_entry] if not show_unstable else chempot_range
                # Don't plot if the slab is unstable, plot if it is.
                if urange != []:
                    label = clean_entry.label
                    if label in already_labelled:
                        label = None
                    else:
                        already_labelled.append(label)
                    if not no_clean:
                        if use_entry_labels:
                            label = clean_entry.label
                        if no_label:
                            label = ""
                        plt = self.chempot_vs_gamma_plot_one(
                            plt,
                            clean_entry,
                            ref_delu,
                            urange,
                            delu_dict=delu_dict,
                            delu_default=delu_default,
                            label=label,
                            JPERM2=JPERM2,
                        )
                if not no_doped:
                    for ads_entry in self.all_slab_entries[hkl][clean_entry]:
                        # Plot the adsorbed slabs
                        # Generate a label for the type of slab
                        urange = stable_u_range_dict[ads_entry] if not show_unstable else chempot_range
                        if urange != []:
                            if use_entry_labels:
                                label = ads_entry.label
                            if no_label:
                                label = ""
                            plt = self.chempot_vs_gamma_plot_one(
                                plt,
                                ads_entry,
                                ref_delu,
                                urange,
                                delu_dict=delu_dict,
                                delu_default=delu_default,
                                label=label,
                                JPERM2=JPERM2,
                            )

        # Make the figure look nice
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, chempot_range, str(ref_delu).split("_")[1], axes, ylim=ylim)

        return plt

    def monolayer_vs_BE(self, plot_eads=False):
        """
        Plots the binding energy as a function of monolayers (ML), i.e.
            the fractional area adsorbate density for all facets. For each
            facet at a specific monlayer, only plot the lowest binding energy.

        Args:
            plot_eads (bool): Option to plot the adsorption energy (binding
                 energy multiplied by number of adsorbates) instead.

        Returns:
            (Plot): Plot of binding energy vs monolayer for all facets.
        """
        plt = pretty_plot(width=8, height=7)
        for hkl in self.all_slab_entries:
            ml_be_dict = {}
            for clean_entry in self.all_slab_entries[hkl]:
                if self.all_slab_entries[hkl][clean_entry]:
                    for ads_entry in self.all_slab_entries[hkl][clean_entry]:
                        if ads_entry.get_monolayer not in ml_be_dict:
                            ml_be_dict[ads_entry.get_monolayer] = 1000
                        be = ads_entry.gibbs_binding_energy(eads=plot_eads)
                        if be < ml_be_dict[ads_entry.get_monolayer]:
                            ml_be_dict[ads_entry.get_monolayer] = be
            # sort the binding energies and monolayers
            # in order to properly draw a line plot
            vals = sorted(ml_be_dict.items())
            monolayers, BEs = zip(*vals)
            plt.plot(monolayers, BEs, "-o", c=self.color_dict[clean_entry], label=hkl)

        adsorbates = tuple(ads_entry.ads_entries_dict)
        plt.xlabel(f"{' '.join(adsorbates)} Coverage (ML)")
        plt.ylabel("Adsorption Energy (eV)") if plot_eads else plt.ylabel("Binding Energy (eV)")
        plt.legend()
        plt.tight_layout()

        return plt

    @staticmethod
    def chempot_plot_addons(plt, xrange, ref_el, axes, pad=2.4, rect=None, ylim=None):
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
        return (Plot): Modified plot with addons.
        """
        # Make the figure look nice
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
        axes.set_xlabel(rf"Chemical potential $\Delta\mu_{{{ref_el}}}$ (eV)")

        ylim = ylim or axes.get_ylim()
        plt.xticks(rotation=60)
        plt.ylim(ylim)
        xlim = axes.get_xlim()
        plt.xlim(xlim)
        plt.tight_layout(pad=pad, rect=rect or [-0.047, 0, 0.84, 1])
        plt.plot([xrange[0], xrange[0]], ylim, "--k")
        plt.plot([xrange[1], xrange[1]], ylim, "--k")
        xy = [np.mean([xrange[1]]), np.mean(ylim)]
        plt.annotate(f"{ref_el}-rich", xy=xy, xytext=xy, rotation=90, fontsize=17)
        xy = [np.mean([xlim[0]]), np.mean(ylim)]
        plt.annotate(f"{ref_el}-poor", xy=xy, xytext=xy, rotation=90, fontsize=17)

        return plt

    def BE_vs_clean_SE(
        self,
        delu_dict,
        delu_default=0,
        plot_eads=False,
        annotate_monolayer=True,
        JPERM2=False,
    ):
        """
        For each facet, plot the clean surface energy against the most
            stable binding energy.

        Args:
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
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
        for hkl in self.all_slab_entries:
            for clean_entry in self.all_slab_entries[hkl]:
                all_delu_dict = self.set_all_variables(delu_dict, delu_default)
                if self.all_slab_entries[hkl][clean_entry]:
                    clean_se = self.as_coeffs_dict[clean_entry]
                    se = sub_chempots(clean_se, all_delu_dict)
                    for ads_entry in self.all_slab_entries[hkl][clean_entry]:
                        ml = ads_entry.get_monolayer
                        be = ads_entry.gibbs_binding_energy(eads=plot_eads)

                        # Now plot the surface energy vs binding energy
                        plt.scatter(se, be)
                        if annotate_monolayer:
                            plt.annotate(f"{ml:.2f}", xy=[se, be], xytext=[se, be])

        plt.xlabel(r"Surface energy ($J/m^2$)") if JPERM2 else plt.xlabel(r"Surface energy ($eV/\AA^2$)")
        plt.ylabel("Adsorption Energy (eV)") if plot_eads else plt.ylabel("Binding Energy (eV)")
        plt.tight_layout()
        plt.xticks(rotation=60)

        return plt

    def surface_chempot_range_map(
        self,
        elements,
        miller_index,
        ranges,
        incr=50,
        no_doped=False,
        no_clean=False,
        delu_dict=None,
        plt=None,
        annotate=True,
        show_unphyiscal_only=False,
        fontsize=10,
    ):
        """
        Adapted from the get_chempot_range_map() method in the PhaseDiagram
            class. Plot the chemical potential range map based on surface
            energy stability. Currently works only for 2-component PDs. At
            the moment uses a brute force method by enumerating through the
            range of the first element chempot with a specified increment
            and determines the chempot rangeo fht e second element for each
            SlabEntry. Future implementation will determine the chempot range
            map first by solving systems of equations up to 3 instead of 2.

        Args:
            elements (list): Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges of
                all Li-Co-O phases wrt to duLi and duO, you will supply
                [Element("Li"), Element("O")]
            miller_index ([h, k, l]): Miller index of the surface we are interested in
            ranges ([[range1], [range2]]): List of chempot ranges (max and min values)
                for the first and second element.
            incr (int): Number of points to sample along the range of the first chempot
            no_doped (bool): Whether or not to include doped systems.
            no_clean (bool): Whether or not to include clean systems.
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            plt (Plot): Plot object to plot on. If None, will create a new plot.
            annotate (bool): Whether to annotate each "phase" with the label of
                the entry. If no label, uses the reduced formula
            show_unphyiscal_only (bool): Whether to only show the shaded region where
                surface energy is negative. Useful for drawing other chempot range maps.
            fontsize (int): Font size of the annotation
        """
        # Set up
        delu_dict = delu_dict or {}
        plt = plt if plt else pretty_plot(12, 8)
        el1, el2 = str(elements[0]), str(elements[1])
        delu1 = Symbol(f"delu_{str(elements[0])}")
        delu2 = Symbol(f"delu_{str(elements[1])}")
        range1 = ranges[0]
        range2 = ranges[1]

        # Find a range map for each entry (surface). This part is very slow, will
        # need to implement a more sophisticated method of getting the range map
        vertices_dict = {}
        for dmu1 in np.linspace(range1[0], range1[1], incr):
            # Get chemical potential range of dmu2 for each increment of dmu1
            new_delu_dict = delu_dict.copy()
            new_delu_dict[delu1] = dmu1
            range_dict, se_dict = self.stable_u_range_dict(
                range2,
                delu2,
                dmu_at_0=True,
                miller_index=miller_index,
                no_doped=no_doped,
                no_clean=no_clean,
                delu_dict=new_delu_dict,
                return_se_dict=True,
            )

            # Save the chempot range for dmu1 and dmu2
            for entry, v in range_dict.items():
                if not v:
                    continue
                if entry not in vertices_dict:
                    vertices_dict[entry] = []

                selist = se_dict[entry]
                vertices_dict[entry].append({delu1: dmu1, delu2: [v, selist]})

        # Plot the edges of the phases
        for entry, v in vertices_dict.items():
            xvals, yvals = [], []

            # Plot each edge of a phase within the borders
            for ii, pt1 in enumerate(v):
                # Determine if the surface energy at this lower range
                # of dmu2 is negative. If so, shade this region.
                if len(pt1[delu2][1]) == 3:
                    if pt1[delu2][1][0] < 0:
                        neg_dmu_range = [pt1[delu2][0][0], pt1[delu2][0][1]]
                    else:
                        neg_dmu_range = [pt1[delu2][0][1], pt1[delu2][0][2]]
                    # Shade the threshold and region at which se<=0
                    plt.plot([pt1[delu1], pt1[delu1]], neg_dmu_range, "k--")
                elif pt1[delu2][1][0] < 0 and pt1[delu2][1][1] < 0 and not show_unphyiscal_only:
                    # Any chempot at this point will result
                    # in se<0, shade the entire y range
                    plt.plot([pt1[delu1], pt1[delu1]], range2, "k--")

                if ii == len(v) - 1:
                    break
                pt2 = v[ii + 1]
                if not show_unphyiscal_only:
                    plt.plot(
                        [pt1[delu1], pt2[delu1]],
                        [pt1[delu2][0][0], pt2[delu2][0][0]],
                        "k",
                    )

                # Need these values to get a good position for labelling phases
                xvals.extend([pt1[delu1], pt2[delu1]])
                yvals.extend([pt1[delu2][0][0], pt2[delu2][0][0]])

            # Plot the edge along the max x value
            pt = v[-1]
            delu1, delu2 = pt
            xvals.extend([pt[delu1], pt[delu1]])
            yvals.extend(pt[delu2][0])
            if not show_unphyiscal_only:
                plt.plot([pt[delu1], pt[delu1]], [pt[delu2][0][0], pt[delu2][0][-1]], "k")

            if annotate:
                # Label the phases
                x = np.mean([max(xvals), min(xvals)])
                y = np.mean([max(yvals), min(yvals)])
                label = entry.label if entry.label else entry.composition.reduced_formula
                plt.annotate(label, xy=[x, y], xytext=[x, y], fontsize=fontsize)

        # Label plot
        plt.xlim(range1)
        plt.ylim(range2)
        plt.xlabel(rf"$\Delta\mu_{{{el1}}} (eV)$", fontsize=25)
        plt.ylabel(rf"$\Delta\mu_{{{el2}}} (eV)$", fontsize=25)
        plt.xticks(rotation=60)

        return plt

    def set_all_variables(self, delu_dict, delu_default):
        """
        Sets all chemical potential values and returns a dictionary where
            the key is a sympy Symbol and the value is a float (chempot).

        Args:
            entry (SlabEntry): Computed structure entry of the slab
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials

        Returns:
            Dictionary of set chemical potential values
        """
        # Set up the variables
        all_delu_dict = {}
        for du in self.list_of_chempots:
            if delu_dict and du in delu_dict:
                all_delu_dict[du] = delu_dict[du]
            elif du == 1:
                all_delu_dict[du] = du
            else:
                all_delu_dict[du] = delu_default

        return all_delu_dict

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
        # def surface_p_vs_t_phase_diagram(self):
        #
        #     return
        #
        # def broken_bond_vs_gamma(self):
        #
        #     return


def entry_dict_from_list(all_slab_entries):
    """
    Converts a list of SlabEntry to an appropriate dictionary. It is
    assumed that if there is no adsorbate, then it is a clean SlabEntry
    and that adsorbed SlabEntry has the clean_entry parameter set.

    Args:
        all_slab_entries (list): List of SlabEntry objects

    Returns:
        (dict): Dictionary of SlabEntry with the Miller index as the main
            key to a dictionary with a clean SlabEntry as the key to a
            list of adsorbed SlabEntry.
    """
    entry_dict = {}

    for entry in all_slab_entries:
        hkl = tuple(entry.miller_index)
        if hkl not in entry_dict:
            entry_dict[hkl] = {}
        clean = entry.clean_entry if entry.clean_entry else entry
        if clean not in entry_dict[hkl]:
            entry_dict[hkl][clean] = []
        if entry.adsorbates:
            entry_dict[hkl][clean].append(entry)

    return entry_dict


class WorkFunctionAnalyzer:
    """
    A class used for calculating the work function
        from a slab model and visualizing the behavior
        of the local potential along the slab.

    .. attribute:: efermi

        The Fermi energy

    .. attribute:: locpot_along_c

        Local potential in eV along points along the  axis

    .. attribute:: vacuum_locpot

        The maximum local potential along the c direction for
            the slab model, ie the potential at the vacuum

    .. attribute:: work_function

        The minimum energy needed to move an electron from the
            surface to infinity. Defined as the difference between
            the potential at the vacuum and the Fermi energy.

    .. attribute:: slab

        The slab structure model

    .. attribute:: along_c

        Points along the c direction with same
            increments as the locpot in the c axis

    .. attribute:: ave_locpot

        Mean of the minimum and maximmum (vacuum) locpot along c

    .. attribute:: sorted_sites

        List of sites from the slab sorted along the c direction

    .. attribute:: ave_bulk_p

        The average locpot of the slab region along the c direction
    """

    def __init__(self, structure: Structure, locpot_along_c, efermi, shift=0, blength=3.5):
        """
        Initializes the WorkFunctionAnalyzer class.

        Args:
            structure (Structure): Structure object modelling the surface
            locpot_along_c (list): Local potential along the c direction
            outcar (MSONable): Outcar vasp output object
            shift (float): Parameter to translate the slab (and
                therefore the vacuum) of the slab structure, thereby
                translating the plot along the x axis.
            blength (float (Ang)): The longest bond length in the material.
                Used to handle pbc for noncontiguous slab layers
        """
        # ensure shift between 0 and 1
        if shift < 0:
            shift += -1 * int(shift) + 1
        elif shift >= 1:
            shift -= int(shift)
        self.shift = shift

        # properties that can be shifted
        slab = structure.copy()
        slab.translate_sites([i for i, site in enumerate(slab)], [0, 0, self.shift])
        self.slab = slab
        self.sorted_sites = sorted(self.slab, key=lambda site: site.frac_coords[2])

        # Get the plot points between 0 and c
        # increments of the number of locpot points
        self.along_c = np.linspace(0, 1, num=len(locpot_along_c))

        # Get the plot points between 0 and c
        # increments of the number of locpot points
        locpot_along_c_mid, locpot_end, locpot_start = [], [], []
        for i, s in enumerate(self.along_c):
            j = s + self.shift
            if j > 1:
                locpot_start.append(locpot_along_c[i])
            elif j < 0:
                locpot_end.append(locpot_along_c[i])
            else:
                locpot_along_c_mid.append(locpot_along_c[i])
        self.locpot_along_c = locpot_start + locpot_along_c_mid + locpot_end

        # identify slab region
        self.slab_regions = get_slab_regions(self.slab, blength=blength)
        # get the average of the signal in the bulk-like region of the
        # slab, i.e. the average of the oscillating region. This gives
        # a rough appr. of the potential in the interior of the slab
        bulk_p = []
        for r in self.slab_regions:
            bulk_p.extend([p for i, p in enumerate(self.locpot_along_c) if r[1] >= self.along_c[i] > r[0]])
        if len(self.slab_regions) > 1:
            bulk_p.extend([p for i, p in enumerate(self.locpot_along_c) if self.slab_regions[1][1] <= self.along_c[i]])
            bulk_p.extend([p for i, p in enumerate(self.locpot_along_c) if self.slab_regions[0][0] >= self.along_c[i]])
        self.ave_bulk_p = np.mean(bulk_p)

        # shift independent quantities
        self.efermi = efermi
        self.vacuum_locpot = max(self.locpot_along_c)
        # get the work function
        self.work_function = self.vacuum_locpot - self.efermi
        # for setting ylim and annotating
        self.ave_locpot = (self.vacuum_locpot - min(self.locpot_along_c)) / 2

    def get_locpot_along_slab_plot(self, label_energies=True, plt=None, label_fontsize=10):
        """
        Returns a plot of the local potential (eV) vs the
            position along the c axis of the slab model (Ang)

        Args:
            label_energies (bool): Whether to label relevant energy
                quantities such as the work function, Fermi energy,
                vacuum locpot, bulk-like locpot
            plt (plt): Matplotlib pylab object
            label_fontsize (float): Fontsize of labels

        Returns plt of the locpot vs c axis
        """
        plt = plt if plt else pretty_plot(width=6, height=4)

        # plot the raw locpot signal along c
        plt.plot(self.along_c, self.locpot_along_c, "b--")

        # Get the local averaged signal of the locpot along c
        xg, yg = [], []
        for i, p in enumerate(self.locpot_along_c):
            # average signal is just the bulk-like potential when in the slab region
            in_slab = False
            for r in self.slab_regions:
                if r[0] <= self.along_c[i] <= r[1]:
                    in_slab = True
            if len(self.slab_regions) > 1:
                if self.along_c[i] >= self.slab_regions[1][1]:
                    in_slab = True
                if self.along_c[i] <= self.slab_regions[0][0]:
                    in_slab = True

            if in_slab or p < self.ave_bulk_p:
                yg.append(self.ave_bulk_p)
                xg.append(self.along_c[i])
            else:
                yg.append(p)
                xg.append(self.along_c[i])
        xg, yg = zip(*sorted(zip(xg, yg)))
        plt.plot(xg, yg, "r", linewidth=2.5, zorder=-1)

        # make it look nice
        if label_energies:
            plt = self.get_labels(plt, label_fontsize=label_fontsize)
        plt.xlim([0, 1])
        plt.ylim([min(self.locpot_along_c), self.vacuum_locpot + self.ave_locpot * 0.2])
        plt.xlabel(r"Fractional coordinates ($\hat{c}$)", fontsize=25)
        plt.xticks(fontsize=15, rotation=45)
        plt.ylabel(r"Potential (eV)", fontsize=25)
        plt.yticks(fontsize=15)

        return plt

    def get_labels(self, plt, label_fontsize=10):
        """
        Handles the optional labelling of the plot with relevant quantities
        Args:
            plt (plt): Plot of the locpot vs c axis
            label_fontsize (float): Fontsize of labels
        Returns Labelled plt
        """
        # center of vacuum and bulk region
        if len(self.slab_regions) > 1:
            label_in_vac = (self.slab_regions[0][1] + self.slab_regions[1][0]) / 2
            if abs(self.slab_regions[0][0] - self.slab_regions[0][1]) > abs(
                self.slab_regions[1][0] - self.slab_regions[1][1]
            ):
                label_in_bulk = self.slab_regions[0][1] / 2
            else:
                label_in_bulk = (self.slab_regions[1][1] + self.slab_regions[1][0]) / 2
        else:
            label_in_bulk = (self.slab_regions[0][0] + self.slab_regions[0][1]) / 2
            if self.slab_regions[0][0] > 1 - self.slab_regions[0][1]:
                label_in_vac = self.slab_regions[0][0] / 2
            else:
                label_in_vac = (1 + self.slab_regions[0][1]) / 2

        plt.plot([0, 1], [self.vacuum_locpot] * 2, "b--", zorder=-5, linewidth=1)
        xy = [label_in_bulk, self.vacuum_locpot + self.ave_locpot * 0.05]
        plt.annotate(
            f"$V_{{vac}}={self.vacuum_locpot:.2f}$",
            xy=xy,
            xytext=xy,
            color="b",
            fontsize=label_fontsize,
        )

        # label the fermi energy
        plt.plot([0, 1], [self.efermi] * 2, "g--", zorder=-5, linewidth=3)
        xy = [label_in_bulk, self.efermi + self.ave_locpot * 0.05]
        plt.annotate(
            f"$E_F={self.efermi:.2f}$",
            xytext=xy,
            xy=xy,
            fontsize=label_fontsize,
            color="g",
        )

        # label the bulk-like locpot
        plt.plot([0, 1], [self.ave_bulk_p] * 2, "r--", linewidth=1.0, zorder=-1)
        xy = [label_in_vac, self.ave_bulk_p + self.ave_locpot * 0.05]
        plt.annotate(
            f"$V^{{interior}}_{{slab}}={self.ave_bulk_p:.2f}$",
            xy=xy,
            xytext=xy,
            color="r",
            fontsize=label_fontsize,
        )

        # label the work function as a barrier
        plt.plot(
            [label_in_vac] * 2,
            [self.efermi, self.vacuum_locpot],
            "k--",
            zorder=-5,
            linewidth=2,
        )
        xy = [label_in_vac, self.efermi + self.ave_locpot * 0.05]
        plt.annotate(
            rf"$\Phi={self.work_function:.2f}$",
            xy=xy,
            xytext=xy,
            fontsize=label_fontsize,
        )

        return plt

    def is_converged(self, min_points_frac=0.015, tol: float = 0.0025):
        """
        A well converged work function should have a flat electrostatic
            potential within some distance (min_point) about where the peak
            electrostatic potential is found along the c direction of the
            slab. This is dependent on the size of the slab.

        Args:
            min_point (fractional coordinates): The number of data points
                +/- the point of where the electrostatic potential is at
                its peak along the c direction.
            tol (float): If the electrostatic potential stays the same
                within this tolerance, within the min_points, it is converged.

        Returns a bool (whether or not the work function is converged)
        """
        conv_within = tol * (max(self.locpot_along_c) - min(self.locpot_along_c))
        min_points = int(min_points_frac * len(self.locpot_along_c))
        peak_i = self.locpot_along_c.index(self.vacuum_locpot)
        all_flat = []
        for i in range(len(self.along_c)):
            if peak_i - min_points < i < peak_i + min_points:
                if abs(self.vacuum_locpot - self.locpot_along_c[i]) > conv_within:
                    all_flat.append(False)
                else:
                    all_flat.append(True)
        return all(all_flat)

    @staticmethod
    def from_files(poscar_filename, locpot_filename, outcar_filename, shift=0, blength=3.5):
        """
        :param poscar_filename: POSCAR file
        :param locpot_filename: LOCPOT file
        :param outcar_filename: OUTCAR file
        :param shift: shift
        :param blength: The longest bond length in the material.
            Used to handle pbc for noncontiguous slab layers
        :return: WorkFunctionAnalyzer
        """
        poscar = Poscar.from_file(poscar_filename)
        locpot = Locpot.from_file(locpot_filename)
        outcar = Outcar(outcar_filename)
        return WorkFunctionAnalyzer(
            poscar.structure,
            locpot.get_average_along_axis(2),
            outcar.efermi,
            shift=shift,
            blength=blength,
        )


class NanoscaleStability:
    """
    A class for analyzing the stability of nanoparticles of different
        polymorphs with respect to size. The Wulff shape will be the
        model for the nanoparticle. Stability will be determined by
        an energetic competition between the weighted surface energy
        (surface energy of the Wulff shape) and the bulk energy. A
        future release will include a 2D phase diagram (e.g. wrt size
        vs chempot for adsorbed or non-stoichiometric surfaces). Based
        on the following work:

        Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale
            stabilization of sodium oxides: Implications for Na-O2
            batteries. Nano Letters, 14(2), 1016-1020.
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

    def solve_equilibrium_point(self, analyzer1, analyzer2, delu_dict=None, delu_default=0, units="nanometers"):
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
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            units (str): Can be nanometers or Angstrom

        Returns:
            Particle radius in nm
        """
        # Set up
        wulff1 = analyzer1.wulff_from_chempot(
            delu_dict=delu_dict or {}, delu_default=delu_default, symprec=self.symprec
        )
        wulff2 = analyzer2.wulff_from_chempot(
            delu_dict=delu_dict or {}, delu_default=delu_default, symprec=self.symprec
        )

        # Now calculate r
        delta_gamma = wulff1.weighted_surface_energy - wulff2.weighted_surface_energy
        delta_E = self.bulk_gform(analyzer1.ucell_entry) - self.bulk_gform(analyzer2.ucell_entry)
        r = (-3 * delta_gamma) / (delta_E)

        return r / 10 if units == "nanometers" else r

    def wulff_gform_and_r(
        self,
        wulffshape,
        bulk_entry,
        r,
        from_sphere_area=False,
        r_units="nanometers",
        e_units="keV",
        normalize=False,
        scale_per_atom=False,
    ):
        """
        Calculates the formation energy of the particle with arbitrary radius r.

        Args:
            wulffshape (WulffShape): Initial, unscaled WulffShape
            bulk_entry (ComputedStructureEntry): Entry of the corresponding bulk.
            r (float (Ang)): Arbitrary effective radius of the WulffShape
            from_sphere_area (bool): There are two ways to calculate the bulk
                formation energy. Either by treating the volume and thus surface
                area of the particle as a perfect sphere, or as a Wulff shape.
            r_units (str): Can be nanometers or Angstrom
            e_units (str): Can be keV or eV
            normalize (bool): Whether or not to normalize energy by volume
            scale_per_atom (True): Whether or not to normalize by number of
                atoms in the particle

        Returns:
            particle formation energy (float in keV), effective radius
        """
        # Set up
        miller_se_dict = wulffshape.miller_energy_dict
        new_wulff = self.scaled_wulff(wulffshape, r)
        new_wulff_area = new_wulff.miller_area_dict

        # calculate surface energy of the particle
        if not from_sphere_area:
            # By approximating the particle as a Wulff shape
            w_vol = new_wulff.volume
            tot_wulff_se = 0
            for hkl, v in new_wulff_area.items():
                tot_wulff_se += miller_se_dict[hkl] * v
            Ebulk = self.bulk_gform(bulk_entry) * w_vol
            new_r = new_wulff.effective_radius

        else:
            # By approximating the particle as a perfect sphere
            w_vol = (4 / 3) * np.pi * r**3
            sphere_sa = 4 * np.pi * r**2
            tot_wulff_se = wulffshape.weighted_surface_energy * sphere_sa
            Ebulk = self.bulk_gform(bulk_entry) * w_vol
            new_r = r

        new_r = new_r / 10 if r_units == "nanometers" else new_r
        e = Ebulk + tot_wulff_se
        e = e / 1000 if e_units == "keV" else e
        e = e / ((4 / 3) * np.pi * new_r**3) if normalize else e
        bulk_struct = bulk_entry.structure
        density = len(bulk_struct) / bulk_struct.lattice.volume
        e = e / (density * w_vol) if scale_per_atom else e

        return e, new_r

    @staticmethod
    def bulk_gform(bulk_entry):
        """
        Returns the formation energy of the bulk
        Args:
            bulk_entry (ComputedStructureEntry): Entry of the corresponding bulk.
        """
        return bulk_entry.energy / bulk_entry.structure.lattice.volume

    def scaled_wulff(self, wulffshape, r):
        """
        Scales the Wulff shape with an effective radius r. Note that the resulting
            Wulff does not necessarily have the same effective radius as the one
            provided. The Wulff shape is scaled by its surface energies where first
            the surface energies are scale by the minimum surface energy and then
            multiplied by the given effective radius.

        Args:
            wulffshape (WulffShape): Initial, unscaled WulffShape
            r (float): Arbitrary effective radius of the WulffShape

        Returns:
            WulffShape (scaled by r)
        """
        # get the scaling ratio for the energies
        r_ratio = r / wulffshape.effective_radius
        miller_list = list(wulffshape.miller_energy_dict)
        # Normalize the magnitude of the facet normal vectors
        # of the Wulff shape by the minimum surface energy.
        se_list = np.array(list(wulffshape.miller_energy_dict.values()))
        # Scale the magnitudes by r_ratio
        scaled_se = se_list * r_ratio

        return WulffShape(wulffshape.lattice, miller_list, scaled_se, symprec=self.symprec)

    def plot_one_stability_map(
        self,
        analyzer,
        max_r,
        delu_dict=None,
        label="",
        increments=50,
        delu_default=0,
        plt=None,
        from_sphere_area=False,
        e_units="keV",
        r_units="nanometers",
        normalize=False,
        scale_per_atom=False,
    ):
        """
        Returns the plot of the formation energy of a particle against its
            effect radius

        Args:
            analyzer (SurfaceEnergyPlotter): Analyzer associated with the
                first polymorph
            max_r (float): The maximum radius of the particle to plot up to.
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            label (str): Label of the plot for legend
            increments (int): Number of plot points
            delu_default (float): Default value for all unset chemical potentials
            plt (pylab): Plot
            from_sphere_area (bool): There are two ways to calculate the bulk
                formation energy. Either by treating the volume and thus surface
                area of the particle as a perfect sphere, or as a Wulff shape.
            r_units (str): Can be nanometers or Angstrom
            e_units (str): Can be keV or eV
            normalize (str): Whether or not to normalize energy by volume
        """
        plt = plt or pretty_plot(width=8, height=7)

        wulffshape = analyzer.wulff_from_chempot(delu_dict=delu_dict, delu_default=delu_default, symprec=self.symprec)

        gform_list, r_list = [], []
        for r in np.linspace(1e-6, max_r, increments):
            gform, r = self.wulff_gform_and_r(
                wulffshape,
                analyzer.ucell_entry,
                r,
                from_sphere_area=from_sphere_area,
                r_units=r_units,
                e_units=e_units,
                normalize=normalize,
                scale_per_atom=scale_per_atom,
            )
            gform_list.append(gform)
            r_list.append(r)

        ru = "nm" if r_units == "nanometers" else r"\AA"
        plt.xlabel(rf"Particle radius (${ru}$)")
        eu = f"${e_units}/{ru}^3$"
        plt.ylabel(rf"$G_{{form}}$ ({eu})")

        plt.plot(r_list, gform_list, label=label)

        return plt

    def plot_all_stability_map(
        self,
        max_r,
        increments=50,
        delu_dict=None,
        delu_default=0,
        plt=None,
        labels=None,
        from_sphere_area=False,
        e_units="keV",
        r_units="nanometers",
        normalize=False,
        scale_per_atom=False,
    ):
        """
        Returns the plot of the formation energy of a particles
            of different polymorphs against its effect radius

        Args:
            max_r (float): The maximum radius of the particle to plot up to.
            increments (int): Number of plot points
            delu_dict (dict): Dictionary of the chemical potentials to be set as
                constant. Note the key should be a sympy Symbol object of the
                format: Symbol("delu_el") where el is the name of the element.
            delu_default (float): Default value for all unset chemical potentials
            plt (pylab): Plot
            labels (list): List of labels for each plot, corresponds to the
                list of se_analyzers
            from_sphere_area (bool): There are two ways to calculate the bulk
                formation energy. Either by treating the volume and thus surface
                area of the particle as a perfect sphere, or as a Wulff shape.
        """
        plt = plt or pretty_plot(width=8, height=7)

        for i, analyzer in enumerate(self.se_analyzers):
            label = labels[i] if labels else ""
            plt = self.plot_one_stability_map(
                analyzer,
                max_r,
                delu_dict,
                label=label,
                plt=plt,
                increments=increments,
                delu_default=delu_default,
                from_sphere_area=from_sphere_area,
                e_units=e_units,
                r_units=r_units,
                normalize=normalize,
                scale_per_atom=scale_per_atom,
            )

        return plt

        # class GetChempotRange:
        #     def __init__(self, entry):
        #         self.entry = entry
        #
        #
        # class SlabEntryGenerator:
        #     def __init__(self, entry):
        #         self.entry = entry


def sub_chempots(gamma_dict, chempots):
    """
    Uses dot product of numpy array to sub chemical potentials
        into the surface grand potential. This is much faster
        than using the subs function in sympy.

    Args:
        gamma_dict (dict): Surface grand potential equation
            as a coefficient dictionary
        chempots (dict): Dictionary assigning each chemical
            potential (key) in gamma a value
    Returns:
        Surface energy as a float
    """
    coeffs = [gamma_dict[k] for k in gamma_dict]
    chempot_vals = []
    for k in gamma_dict:
        if k not in chempots:
            chempot_vals.append(k)
        elif k == 1:
            chempot_vals.append(1)
        else:
            chempot_vals.append(chempots[k])

    return np.dot(coeffs, chempot_vals)
