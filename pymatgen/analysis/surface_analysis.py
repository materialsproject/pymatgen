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

import numpy as np
from scipy.stats import linregress
import itertools
import warnings
import random, copy

from pymatgen.core.structure import Structure, Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import Slab
from pymatgen.analysis.wulff import WulffShape
from pymatgen import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
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
        chempot_ranges = pd.get_chempot_range_map(self.ref_el_comp.elements)
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
            # For elemental system
            ucell_comp = self.ucell_entry.composition
            chempot_range = [chempot_ranges[entry] for entry in chempot_ranges.keys()
                             if entry.composition ==
                             ucell_comp][0][0]._coords if \
                chempot_ranges else [[-1, -1], [0, 0]]

        chempot_range = list(chempot_range)
        return sorted([chempot_range[0][0], chempot_range[1][0]])

    def chempot_range_adsorption(self, ads_slab_entry, clean_slab_entry,
                                 const_u_ref, buffer=0.2):

        c1 = self.surface_energy_coefficients(clean_slab_entry)
        c2 = self.surface_energy_coefficients(clean_slab_entry,
                                              ads_slab_entry=ads_slab_entry)
        umin, gamma = self.solve_2_linear_eqns(c1, c2, x_is_u_ads=True,
                                               const_u=const_u_ref)
        umin = -1*10e-5 if not umin else umin
        # Make sure upper limit of u doesn't lead to negative or 0 values.
        # Substract a value approaching 0 to avoid surface energies of 0
        umax = (1*10e-5-const_u_ref*c2[0]-c2[2])/c2[1] if \
            0 > self.calculate_gamma_at_u(clean_slab_entry,
                                          ads_slab_entry=ads_slab_entry,
                                          u_ref=const_u_ref, u_ads=0) else 0
        return [umin-abs(umin)*buffer, umax]

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
        # b3 is the intercept (clean surface energy for a clean, stoichiometric system)
        g_y = [entry.energy_per_atom for entry in self.reactants
               if entry.composition.reduced_composition == self.ref_el_comp][0]
        b3 = (1 / (2 * Aclean)) * (clean_slab_entry.energy - (Nx/self.x)*self.gbulk - \
                              (Ny - (self.y / self.x) * Nx) * g_y) + \
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

    def solve_2_linear_eqns(self, c1, c2, x_is_u_ads=False, const_u=0):

        """
        Helper method to solve the intersect for two linear equations.
        """

        # set one of the terms as a constant (either the adsorption
        # or clean term) to get 2 linear eqns of 1 variable
        i = 0 if x_is_u_ads else 1
        b11 = np.dot([c1[i], c1[2]], [const_u, 1])
        b12 = np.dot([c2[i], c2[2]], [const_u, 1])
        i = 1 if x_is_u_ads else 0

        # If the slopes are equal, i.e. lines are parallel
        if c1[i] == c2[i]:
            return [None, None]
        # Now solve the two eqns, return [del_u, gamma]
        return np.linalg.solve([[c1[i], -1], [c2[i], -1]],
                               [-1*b11, -1*b12])


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

    def __init__(self, entry_dict, surface_energy_calculator, custom_chempot_range=[]):
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
        self.chempot_range = custom_chempot_range if custom_chempot_range \
            else surface_energy_calculator.chempot_range()
        self.entry_dict = entry_dict
        self.color_dict = self.color_palette_dict()
        self.ref_el_comp = str(self.se_calculator.ref_el_comp.elements[0])

    def max_adsorption_chempot_range(self, const_u_ref, buffer=0.1):

        all_ranges = []
        for hkl in self.entry_dict.keys():
            for clean_entry in self.entry_dict[hkl].keys():
                for ads_entry in self.entry_dict[hkl][clean_entry]:
                    all_ranges.append(self.se_calculator.\
                                      chempot_range_adsorption(ads_entry, clean_entry,
                                                               const_u_ref=const_u_ref,
                                                               buffer=buffer))

        if not all_ranges:
            # If there is no intersection, the range is [-1,0]
            return [-1,0]
        # ensure our lower limit is at an intersection with the clean slab
        all_ranges = sorted(all_ranges, key=lambda r: r[0])
        max_range = [all_ranges[0][0]]
        # ensure our upper limit corresponds to gamm > 0 or if gamma > 0 when u = 0
        all_ranges = sorted(all_ranges, key=lambda r: r[1])
        max_range.append(all_ranges[0][1])
        return max_range

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
            entry, gamma = self.return_stable_slab_entry_at_u(hkl, u_ref=u_ref,
                                                              u_ads=u_ads)
            e_surf_list.append(gamma)

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
                gamma = self.se_calculator.calculate_gamma_at_u(entry, u_ref=u_ref, u_ads=u_ads,
                                                                ads_slab_entry=ads_entry)
                all_entries.append(ads_entry)
                all_gamma.append(gamma)

        return all_entries[all_gamma.index(min(all_gamma))], min(all_gamma)

    def area_frac_vs_chempot_plot(self, u_const=0, increments=10,
                                  x_is_u_ads=False):
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

        xrange = self.chempot_range if not x_is_u_ads \
            else self.max_adsorption_chempot_range(u_const)
        all_chempots = np.linspace(min(xrange), max(xrange),
                                   increments)

        # initialize a dictionary of lists of fractional areas for each hkl
        hkl_area_dict = {}
        for hkl in self.entry_dict.keys():
            hkl_area_dict[hkl] = []

        # Get plot points for each Miller index
        for u in all_chempots:
            u_ads = u if x_is_u_ads else u_const
            u_ref = u if not x_is_u_ads else u_const
            wulffshape = self.wulff_shape_from_chempot(u_ads=u_ads, u_ref=u_ref)

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
        self.chempot_plot_addons(plt, xrange, axes, pad=5,
                                 rect=[-0.0, 0, 0.95, 1],
                                 x_is_u_ads=x_is_u_ads,
                                 ylim=[0,1])

        return plt

    def chempot_vs_gamma_plot_one(self, plt, clean_entry, label='',
                                  ads_entry=None, JPERM2=False,
                                  x_is_u_ads=False, const_u=0,
                                  urange=None):

        if x_is_u_ads:
            u_ref_range = [const_u, const_u]
            u_ads_range = urange if urange else self.max_adsorption_chempot_range(const_u)
        else:
            u_ref_range = urange if urange else self.chempot_range
            u_ads_range = [const_u, const_u]

        # use dashed lines for slabs that are not stoichiometric
        # wrt bulk. Label with formula if nonstoichiometric
        mark = '--' if clean_entry.composition.reduced_composition != \
                       self.se_calculator.ucell_entry. \
                           composition.reduced_composition else '-'

        gamma_range = [self.se_calculator.calculate_gamma_at_u(clean_entry,
                                                               ads_slab_entry=ads_entry,
                                                               u_ref=u_ref_range[0],
                                                               u_ads=u_ads_range[0]),
                       self.se_calculator.calculate_gamma_at_u(clean_entry,
                                                               ads_slab_entry=ads_entry,
                                                               u_ref=u_ref_range[1],
                                                               u_ads=u_ads_range[1])]

        se_range = np.array(gamma_range) * EV_PER_ANG2_TO_JOULES_PER_M2 \
            if JPERM2 else gamma_range
        if not urange:
            urange = u_ref_range if not x_is_u_ads else u_ads_range
        color = self.color_dict[ads_entry] if ads_entry else self.color_dict[clean_entry]
        plt.plot(urange, se_range, mark, color=color, label=label)

        return plt

    def chempot_vs_gamma_clean(self, miller_index=(), JPERM2=False, plt=None):
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

        plt = plt if plt else pretty_plot(width=8, height=7)
        axes = plt.gca()

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
                    label = None
                else:
                    already_labelled.append(label)

                self.chempot_vs_gamma_plot_one(plt, entry, label=label,
                                               JPERM2=JPERM2, x_is_u_ads=False)

        # Make the figure look nice
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, self.chempot_range, axes)

        return plt

    def chempot_vs_gamma_facet(self, miller_index=(), const_u=0,
                               JPERM2=False, show_unstable=False):

        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        # Plot wrt to adsorption chempot if any adsorption entries exist
        x_is_u_ads = False
        for hkl in self.entry_dict.keys():
            for clean_entry in self.entry_dict[hkl].keys():
                if self.entry_dict[hkl][clean_entry]:
                    x_is_u_ads = True

        for hkl in self.entry_dict.keys():
            if miller_index and hkl != tuple(miller_index):
                continue
            if not show_unstable:
                clean_only = False if x_is_u_ads else True
                stable_u_range_dict = self.stable_u_range_dict(clean_only=clean_only,
                                                               const_u=const_u,
                                                               miller_index=hkl)

            already_labelled = []
            for clean_entry in self.entry_dict[hkl]:

                urange = stable_u_range_dict[clean_entry] if not show_unstable else None
                if urange != []:

                    label = self.create_slab_label(clean_entry, miller_index=hkl)
                    if label in already_labelled:
                        label = None
                    else:
                        already_labelled.append(label)

                    self.chempot_vs_gamma_plot_one(plt, clean_entry, label=label,
                                                   JPERM2=JPERM2, x_is_u_ads=x_is_u_ads,
                                                   const_u=const_u, urange=urange)

                for ads_entry in self.entry_dict[hkl][clean_entry]:
                    # Plot the adsorbed slabs
                    # Generate a label for the type of slab
                    urange = stable_u_range_dict[ads_entry] if not show_unstable else None
                    if urange != []:
                        self.chempot_vs_gamma_plot_one(plt, clean_entry, JPERM2=JPERM2,
                                                       ads_entry=ads_entry,
                                                       const_u=const_u, urange=urange,
                                                       x_is_u_ads=x_is_u_ads)

        const_species = self.se_calculator.adsorbate_as_str  if not \
            x_is_u_ads else self.ref_el_comp
        if miller_index:
            plt.title(r"%s,  $\Delta\mu_{%s}=%.2f$" %(str(miller_index),
                                                      const_species, const_u))

        # Make the figure look nice
        xrange = self.chempot_range if not x_is_u_ads \
            else self.max_adsorption_chempot_range(const_u)
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        plt = self.chempot_plot_addons(plt, xrange, axes,
                                       x_is_u_ads=x_is_u_ads)


        return plt

    def chempot_plot_addons(self, plt, xrange, axes,
                            pad=2.4, rect=[-0.047, 0, 0.84, 1],
                            x_is_u_ads=False, ylim=[]):

        # Make the figure look nice
        x_species = self.ref_el_comp if not \
            x_is_u_ads else self.se_calculator.adsorbate_as_str
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        axes.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" % (x_species))

        ylim = ylim if ylim else axes.get_ylim()
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
            label += r"+%s" %(self.se_calculator.adsorbate_as_str)
        return label

    def color_palette_dict(self, alpha=0.35):

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

    def stable_u_range_dict(self, clean_only=True, const_u=0,
                            buffer=0.1, miller_index=()):

        all_intesects_dict = {}
        stable_urange_dict = {}
        max_ads_range = self.max_adsorption_chempot_range(const_u,
                                                          buffer=buffer)
        standard_range = max_ads_range if not clean_only else self.chempot_range
        standard_range = sorted(standard_range)
        x_is_u_ads = False if clean_only else True

        # Get all entries for a specific facet
        for hkl in self.entry_dict.keys():

            # bool to check if at least one surface exists for a facet
            entry_exists = False

            if miller_index and hkl != tuple(miller_index):
                continue
            entries_in_hkl = [entry for entry in self.entry_dict[hkl]]
            if not clean_only:
                ads_entries_in_hkl = []
                for entry in self.entry_dict[hkl]:
                    ads_entries_in_hkl.extend([ads_entry for ads_entry in
                                               self.entry_dict[hkl][entry]])
                entries_in_hkl.extend(ads_entries_in_hkl)

            # if there is only one entry for this facet, then just give it the
            # default urange, you can't make combinations with just 1 item
            if len(entries_in_hkl) == 1:
                stable_urange_dict[entries_in_hkl[0]] = standard_range
                continue
            for pair in itertools.combinations(entries_in_hkl, 2):
                # Check if entry is adsorbed entry or clean, this is
                # a hassle so figure out a cleaner way to do this
                clean_ads_coeffs = []
                for p in pair:
                    p1 = self.get_clean_ads_entry_pair(p)
                    c = self.se_calculator.surface_energy_coefficients(p1[0],
                                                                       ads_slab_entry=p1[1])
                    clean_ads_coeffs.append(c)

                u, gamma = self.se_calculator.solve_2_linear_eqns(clean_ads_coeffs[0],
                                                                  clean_ads_coeffs[1],
                                                                  x_is_u_ads=x_is_u_ads,
                                                                  const_u=const_u)

                if u:
                    # If the u in a list is beyond the standard
                    # range, set it to one of the limits
                    if u < standard_range[0]:
                        u_new = standard_range[0]
                    elif u > standard_range[1]:
                        u_new = standard_range[1]
                    else:
                        u_new = u

                    for entry in pair:
                        if entry not in all_intesects_dict.keys():
                            all_intesects_dict[entry] = []
                        all_intesects_dict[entry].append(u_new)

            # Now that we have every single intersection for a given
            # slab, find the range of u where each slab is stable
            for entry in all_intesects_dict.keys():
                if entry not in stable_urange_dict.keys():
                    stable_urange_dict[entry] = []
                for u_int in all_intesects_dict[entry]:
                    stable_entries = []
                    for i in [-1, 1]:
                        u_ads = u_int+i*(10e-6) if x_is_u_ads else const_u
                        u_ref = u_int+i*(10e-6) if not x_is_u_ads else const_u
                        # return_stable_slab_entry_at_u() will only return one entry for one gamma,
                        # since u is at an intersection, this entry is ambiguous, we need to get
                        # the entry slightly above and below u and check if the current entry is
                        # any of these entries. Another inconvenience that needs to be fixed

                        stable_entry, gamma = self.return_stable_slab_entry_at_u(hkl,
                                                                                 u_ads=u_ads,
                                                                                 u_ref=u_ref)
                        stable_entries.append(stable_entry)
                    # If this entry is stable at this u, append u
                    if entry in stable_entries:
                        stable_urange_dict[entry].append(u_int)
                        entry_exists = True

            # Now check for entries with only one intersection
            # is it stable below or above u_intersect
            for entry in stable_urange_dict.keys():
                # First lets check if all the u values for
                # an entry are the same as the standard range,
                # if so, just set it to the standard range
                if stable_urange_dict[entry]:
                    if all([u in standard_range for u in stable_urange_dict[entry]]):
                        stable_urange_dict[entry] = standard_range

                # If only one u, its stable from u to +-inf.
                # If no u, this entry is never stable
                if len(stable_urange_dict[entry]) == 1:
                    u = stable_urange_dict[entry][0]
                    for i in [-1, 1]:
                        u_ads = u+i*(10e-6) if x_is_u_ads else const_u
                        u_ref = u+i*(10e-6) if not x_is_u_ads else const_u

                        e, se = self.return_stable_slab_entry_at_u(hkl,
                                                                   u_ads=u_ads,
                                                                   u_ref=u_ref)
                        if e == entry:
                            # If the entry stable below u, assume it is
                            # stable at -inf, otherwise its stable at +inf
                            u2 = standard_range[0] if i == -1 else standard_range[1]
                            stable_urange_dict[entry].append(u2)
                            entry_exists = True

                # now sort the ranges for each entry
                stable_urange_dict[entry] = sorted(stable_urange_dict[entry])
            # Now we make sure that each facet has at least
            # one entry, if no entries exist, this means there
            # is not intersection, get the most stable surface
            if not entry_exists:
                e, se = self.return_stable_slab_entry_at_u(hkl,
                                                           u_ads=const_u,
                                                           u_ref=const_u)
                stable_urange_dict[e] = standard_range

        return stable_urange_dict

    def get_clean_ads_entry_pair(self, entry):

        """
        Returns a pair of entries, the first is the clean entry, the second
        is either nonetype (if the initial entry is clean) or the
        corresponding adsorbed entry
        """

        for hkl in self.entry_dict.keys():
            if entry in list(self.entry_dict[hkl].keys()):
                # its a clean entry
                return [entry, None]
            else:
                for clean_entry in self.entry_dict[hkl].keys():
                    if entry in self.entry_dict[hkl][clean_entry]:
                        return [clean_entry, entry]

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