# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

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


class SurfaceEnergyAnalyzer(object):

    """
    A class used for analyzing the surface energies of a material of a given
        material_id. By default, this will use entries calculated from the
        Materials Project to obtain chemical potential and bulk energy. As a
        result, the difference in VASP parameters between the user's entry
        (vasprun_dict) and the parameters used by Materials Project, may lead
        to a rough estimate of the surface energy. For best results, it is
        recommend that the user calculates all decomposition components first,
        and insert the results into their own database as a pymatgen-db entry
        and use those entries instead (custom_entries). In addition, this code
        will only use one bulk entry to calculate surface energy. Ideally, to
        get the most accurate surface energy, the user should compare their
        slab energy to the energy of the oriented unit cell with both calculations
        containing consistent k-points to avoid converegence problems as the
        slab size is varied. See:
            Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
                Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.
        and
            Rogal, J., & Reuter, K. (2007). Ab Initio Atomistic Thermodynamics for
                Surfaces : A Primer. Experiment, Modeling and Simulation of Gas-Surface
                Interactions for Reactive Flows in Hypersonic Flights, 2–1 – 2–18.

    .. attribute:: ref_element

        All chemical potentials cna be written in terms of the range of chemical
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

    .. attribute:: chempot_range

        List of the min and max chemical potential of ref_element.

    .. attribute:: e_of_element

        Energy per atom of ground state ref_element, eg. if ref_element=O,
            than e_of_element=1/2*E_O2.

    .. attribute:: vasprun_dict

        Nested dictionary containing a list of Vaspruns for slab calculations as
            items and the corresponding Miller index of the slab as the key.
            To account for adsorption, each value is a sub-dictionary with the
            vasprun of a clean slab calculation as the sub-key and a list of
            vaspruns for adsorption calculations as the sub-value. The sub-value
            can contain different adsorption configurations such as a different
            site or a different coverage, however, ordinarily only the most stable
            configuration for a particular coverage will be considered as the
            function of the adsorbed surface energy has an intercept dependent on
            the adsorption energy (ie an adsorption site with a higher adsorption
            energy will always provide a higher surface energy than a site with a
            lower adsorption energy). An example parameter is provided:
            {(h1,k1,l1): {clean_vrun1: [ads_vrun1, ads_vrun2, ...],
                          clean_vrun2: [...], ...}, (h2,k2,l2): {...}}
            where clean_vrun1 can be a pristine surface and clean_vrun2 can be a
            reconstructed surface while ads_vrun1 can be adsorption at site 1 with
            a 2x2 coverage while ads_vrun2 can have a 3x3 coverage.

    .. attribute:: adsorbate

        Composition of the adsorbate (if there is one).

    """

    def __init__(self, ucell_entry, vasprun_dict, comp1, ref_el_comp,
                 exclude_ids=[], custom_entries=[], mapi_key=None,
                 full_chempot=False, adsorbate=None):
        """
        Analyzes surface energies and Wulff shape of a particular
            material using the chemical potential.
        Args:
            ucell_entry (material_id or computed_entry): Materials Project or entry
                of the bulk system the slab is based on (a string, e.g., mp-1234).
            vasprun_dict (dict): Dictionary containing a list of Vaspruns
                for slab calculations. See attributes.
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
            adsorbate (Composition): Composition of adsorbate, defaults to None
        """

        self.comp1 = comp1
        self.ref_el_comp = ref_el_comp
        self.mprester = MPRester(mapi_key) if mapi_key else MPRester()
        self.ucell_entry = self.mprester.get_entry_by_material_id(\
            ucell_entry, inc_structure=True) \
            if type(ucell_entry).__name__ == "str" else ucell_entry
        ucell_comp = self.ucell_entry.composition

        entries = [entry for entry in
                   self.mprester.get_entries_in_chemsys(list(\
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
        gbulk = self.ucell_entry.energy /\
                ucell_comp.get_integer_formula_and_factor()[1]

        pd = PhaseDiagram(entries)
        pda = PDAnalyzer(pd)
        chempot_ranges = pda.get_chempot_range_map(self.ref_el_comp.elements)
        # If no chemical potential is found, we return u=0, eg.
        # for a elemental system, the relative u of Cu for Cu is 0
        if full_chempot:
            if not chempot_ranges:
                chempot_range = [[0,0], [0,0]]
            else:
                all_u = []
                for entry in chempot_ranges.keys():
                    all_u.extend(chempot_ranges[entry][0]._coords)
                chempot_range = [min(all_u), max(all_u)]
        else:
            chempot_range = [chempot_ranges[entry] for entry in chempot_ranges.keys()
                             if entry.composition ==
                             ucell_comp][0][0]._coords if \
                chempot_ranges else [[0,0], [0,0]]

        self.x = x
        self.y = y
        self.gbulk = gbulk
        chempot_range = list(chempot_range)
        self.chempot_range = sorted([chempot_range[0][0], chempot_range[1][0]])
        self.vasprun_dict = vasprun_dict
        self.adsorbate = adsorbate

    def calculate_slope_and_intercept(self, vasprun):
        """
        Calculates the slope and intercept of the surface energy for a single slab.
        Args:
            vasprun (Vasprun): A Vasprun object

        Returns (list): The surface energy for the minimum/maximun
            chemical potential and the second list gives the range
            of the chemical potential
        """

        reduced_comp = self.ucell_entry.composition.reduced_composition.as_dict()
        # Get the composition in the slab
        slab = vasprun.final_structure
        comp = slab.composition.as_dict()

        if len(reduced_comp.keys()) == 1:
            Nx = Ny = comp[self.ucell_entry.structure[0].species_string]
        else:
            # Ny = abs(r.coeffs[r.all_comp.index(self.ref_el_comp)])
            # Nx = abs(r.coeffs[r.all_comp.index(self.comp1)])
            Nx = slab.composition.as_dict()[str(self.comp1.elements[0])]
            Ny = slab.composition.as_dict()[str(self.ref_el_comp.elements[0])]

        # Calculate surface area
        m = slab.lattice.matrix
        A = np.linalg.norm(np.cross(m[0], m[1]))

        # return the slope and intercept
        slope = (-1 / (2 * A)) * (Ny - (self.y / self.x) * Nx)
        intercept = (1/(2*A))*(vasprun.final_energy-(Nx/self.x)*self.gbulk-\
                               (Ny-(self.y/self.x)*Nx)*[entry.energy_per_atom for entry in
                                                        self.reactants if
                                                        entry.composition.reduced_composition
                                                        == self.ref_el_comp][0])
        return slope, intercept

    def calculate_Eads(self, vasprun_ads, vasprun_clean, adsorbate):
        """
        Calculates the adsorption energy or Gibb's
        binding energy of an adsorbate on a surface
        Args:
            vasprun_ads (Vasprun): The Vasprun of the adsorbed slab
            vasprun_clean (Vasprun): The Vasprun of the clean slab
            adsorbate (str): The adsorbate as a string
        """

        m = vasprun_ads.final_structure.lattice.matrix
        A_ads = np.linalg.norm(np.cross(m[0], m[1]))
        m = vasprun_clean.final_structure.lattice.matrix
        A_clean = np.linalg.norm(np.cross(m[0], m[1]))
        n = (A_ads/A_clean)
        Nads = vasprun_ads.get_computed_entry().composition.as_dict()[adsorbate]
        return vasprun_ads.final_energy - n*vasprun_clean.final_energy \
               - Nads*self.ucell_entry.energy_per_atom

    def calculate_gamma_ads_range(self, vasprun_ads, vasprun_clean, adsorbate, u):
        """
        Calculates the surface energy for an adsorbed slab using the
        adsorption energy and chemical potential of the adsorbate
        Args:
            vasprun_ads (Vasprun): The Vasprun of the adsorbed slab
            vasprun_clean (Vasprun): The Vasprun of the clean slab
            adsorbate (str): The adsorbate as a string
            u (float): The chemical potential of the
                adsorbate to calculate the surface energy at
        """

        struct = vasprun_ads.final_structure
        weights = [s.species_and_occu.weight for s in struct]
        center_of_mass = np.average(struct.frac_coords,
                                    weights=weights, axis=0)
        nsurfs = 0
        if any([site.species_string == adsorbate for site in
                struct if site.frac_coords[2] > center_of_mass]):
            nsurfs += 1
        if any([site.species_string == adsorbate for site in
                struct if site.frac_coords[2] < center_of_mass]):
            nsurfs += 1

        gamma_clean = self.calculate_gamma_at_u(vasprun_clean, u)
        Eads = self.calculate_Eads(vasprun_ads, vasprun_clean, adsorbate)
        Nads = vasprun_ads.get_computed_entry().composition.as_dict()[adsorbate]
        m = vasprun_ads.final_structure.lattice.matrix
        A = np.linalg.norm(np.cross(m[0], m[1]))

        return gamma_clean - (Nads*u - Eads)/(nsurfs*A)

    def calculate_gamma_at_u(self, vasprun, u):
        """
        Quickly calculates the surface energy for
        the slab of the vasprun file at a specific u
        args:
            vasprun (Vasprun): Vasprun containing the final energy and structure
                of the slab whose surface energy we want ot calculate
            u (float): The chemical potential at which
                we want to claculate the surface energy

        Returns (float): surface energy
        """

        slope, intercept = self.calculate_slope_and_intercept(vasprun)

        return slope*u+intercept

    def wulff_shape_from_chempot(self, chempot, symprec=1e-5):
        """
        Method to get the Wulff shape at a specific chemical potential.
        Args:
            chempot (float): The chemical potential the Wulff Shape exist in.
        """

        # Check if the user provided chemical potential is within the
        # predetermine range of chemical potential. If not, raise a warning
        if not max(self.chempot_range) >= chempot >= min(self.chempot_range):
            warnings.warn("The provided chemical potential is outside the range "
                          "of chemical potential (%s to %s). The resulting Wulff "
                          "shape might not be reasonable." %(min(self.chempot_range),
                                                             max(self.chempot_range)))

        latt = SpacegroupAnalyzer(self.ucell_entry.structure).\
            get_conventional_standard_structure().lattice

        miller_list = self.vasprun_dict.keys()
        e_surf_list = []
        for hkl in miller_list:
            # At each possible configuration, we calculate surface energy as a
            # function of u and take the lowest surface energy (corresponds to
            # the most stable slab termination at that particular u)
            e_list = []
            for vasprun in self.vasprun_dict[hkl]:
                slope, intercept = self.calculate_slope_and_intercept(vasprun)
                e_list.append(slope * chempot + intercept)
            e_surf_list.append(min(e_list))

        return WulffShape(latt, miller_list, e_surf_list, symprec=symprec)

    def wulff_shape_dict(self, symprec=1e-5, at_intersections=False):
        """
        As the surface energy is a function of chemical potential, so too is the
            Wulff shape. This methods generates a dictionary of Wulff shapes at
            certain chemical potentials where a facet goes through a transition.
            Returns a dict, eg. {chempot1: WulffShape1, chempot2: WulffShape2}

        Args:
            symprec (float): for recp_operation, default is 1e-5.
            at_intersections (bool): Whether to generate a Wulff shape for each
                intersection of surface energy for a specific facet (eg. at the
                point where a (111) stoichiometric surface energy plot intersects
                with the (111) nonstoichiometric plot) or to just generate two
                Wulff shapes, one at the min and max chemical potential.
        """

        # First lets get the Wulff shape at the
        # minimum and maximum chemical potential
        wulff_dict = {self.chempot_range[0]: \
                          self.wulff_shape_from_chempot(self.chempot_range[0],
                                                        symprec=symprec),
                      self.chempot_range[1]: \
                          self.wulff_shape_from_chempot(self.chempot_range[1],
                                                        symprec=symprec)}

        # Now we get the Wulff shape each time a facet changes its configuration
        # (ie, adsorption coverage, stoichiometric to nonstoichiometric, etc)
        if at_intersections:
            # Get all values of chemical potential where an intersection occurs
            all_intersections = self.get_intersections()
            # Get a Wulff shape for each intersection. The change in the Wulff shape
            # will vary if the rate of change in surface energy for any facet changes
            for int in all_intersections:
                wulff = self.wulff_shape_from_chempot(int[0], symprec=symprec)
                if any([wulff.area_fraction_dict[hkl] != 0 for hkl in int[2]]):
                    wulff_dict[int[0]] = wulff

        return wulff_dict

    def get_stable_surf_regions(self, miller_index):
        """
        For a specific facet, returns an energy stability range of u for each
            facet as a slope and intercept. e.g. if config 1 stable between
            -1<u<-0.5 and config 2 is stable between -0.5<u<0 we return:
            [[-1, -0.5, slope1, intercept1], [-0.5, 0, slope2, intercept2]].
        """

        # First get the intercepts at the stable configs
        vasprun, gamma = self.return_stable_slab_at_u(miller_index,
                                                      min(self.chempot_range))
        ulist, gamma_list = [min(self.chempot_range)], [gamma]
        intersections = self.get_intersections(miller_index=miller_index)
        for int in intersections:
            vasprun, gamma = self.return_stable_slab_at_u(miller_index, int[0])
            ulist.append(int[0])
            gamma_list.append(gamma)

        # Next, build a stability map with the range in u, slope and intercept
        stability_map = []
        v, final_gamma = self.return_stable_slab_at_u(miller_index,
                                                      max(self.chempot_range))
        for i, u in enumerate(ulist):
            high_u = max(self.chempot_range) if i == len(ulist)-1 else ulist[i+1]
            high_gamma = final_gamma if i == len(ulist)-1 else gamma_list[i+1]
            stability_map.append([[u, high_u],
                                  [gamma_list[i], high_gamma]])

        return stability_map

    def return_stable_slab_at_u(self, miller_index, u):
        """
        Returns the vasprun corresponding to the most
        stable slab for a particular facet at a specific u
        Args:
            miller_index ((h,k,l)): The facet to find the most stable slab in
            u (float): The chemical potential to look for the most stable slab
        """

        all_gamma = [self.calculate_gamma_at_u(v, u)
                     for v in self.vasprun_dict[miller_index]]
        return self.vasprun_dict[miller_index][all_gamma.index(min(all_gamma))], min(all_gamma)


    def get_intersections(self, miller_index=()):
        """
        Returns a all intersections for a specific facet or for all facets.
            A list of lists, each entry in the list represents an intersect
            with u, gamma(u) and [hkl1, hkl2]. Useful for finding when the
            configuration of a particular facet changes or when one facet
            is more stable than the other.

        Args:
            miller_index ((h, k, l)): Miller index of the facet we
                are interested in. Optional parameter that looks for
                intersections for a specfic facet only when given.
        """

        # First lets calculate the range of surface energies for
        # all terminations of a specific facet or all facets.
        if miller_index:
            all_slope_intercepts = [[self.calculate_slope_and_intercept(vasprun), miller_index]
                                    for vasprun in self.vasprun_dict[miller_index]]
        else:
            all_slope_intercepts = []
            for hkl in self.vasprun_dict.keys():
                slope_intercept = [[self.calculate_slope_and_intercept(vasprun), hkl]
                                   for vasprun in self.vasprun_dict[hkl]]
                all_slope_intercepts.extend(slope_intercept)

        if len(all_slope_intercepts) == 1:
            return []

        # Now get all possible intersection coordinates for each pair of lines
        intersections = []
        for pair_lines in itertools.combinations(all_slope_intercepts, 2):
            slope1, intercept1 = pair_lines[0][0]
            slope2, intercept2 = pair_lines[1][0]

            if slope1 - slope2 == 0:
                # i.e. lines are parallel
                continue

            # Calculate the intersection
            u, gamma = np.linalg.solve([[slope1, -1],[slope2, -1]],
                                       [-1*intercept1, -1*intercept2])

            # if the intersection is beyond the chemical potential
            # range or if the lines are parallel, we ignore it
            if u < min(self.chempot_range) \
                    or u > max(self.chempot_range):
                continue

            # If the surface energies at this u for both lines facet is unstable, ignore it
            v1, gamma1 = self.return_stable_slab_at_u(pair_lines[0][1], u)
            v2, gamma2 = self.return_stable_slab_at_u(pair_lines[1][1], u)
            # +10e-9 to handle floating point comparison for equal values
            if all([gamma > gamma1+10e-9, gamma > gamma2+10e-9]):
                continue
            # Each list is [u, gamma, hkl1, hkl2]
            intersections.append([u, gamma,
                                 [pair_lines[0][1], pair_lines[1][1]]])

        return sorted(intersections, key=lambda ints: ints[0])

    def area_frac_vs_chempot_plot(self, cmap=cm.jet, at_intersections=False,
                                  increments=10):
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
        f = [int(i) for i in np.linspace(0, 255, len(self.vasprun_dict.keys()))]

        # Get all points of min/max chempot and intersections
        chempot_intersections = []
        chempot_intersections.extend(self.chempot_range)
        for hkl in self.vasprun_dict.keys():
            chempot_intersections.extend([ints[0] for ints in
                                          self.get_intersections(hkl)])
        chempot_intersections = sorted(chempot_intersections)

        # Get all chempots
        if at_intersections:
            all_chempots = []
            for i, intersection in enumerate(chempot_intersections):
                if i < len(chempot_intersections)-1:
                    all_chempots.extend(np.linspace(intersection,
                                                    chempot_intersections[i+1],
                                                    increments))
        else:
            all_chempots = np.linspace(min(self.chempot_range),
                                       max(self.chempot_range), increments)

        # initialize a dictionary of lists of fractional areas for each hkl
        hkl_area_dict = {}
        for hkl in self.vasprun_dict.keys():
            hkl_area_dict[hkl] = []

        # Get plot points for each Miller index
        for u in all_chempots:
            wulffshape = self.wulff_shape_from_chempot(u)
            for hkl in wulffshape.area_fraction_dict.keys():
                hkl_area_dict[hkl].append(wulffshape.area_fraction_dict[hkl])

        # Plot the area fraction vs chemical potential for each facet
        plt = pretty_plot(width=8, height=7)
        axes = plt.gca()

        for i, hkl in enumerate(self.vasprun_dict.keys()):
            # Ignore any facets that never show up on the
            # Wulff shape regardless of chemical potential
            if all([a == 0 for a in hkl_area_dict[hkl]]):
                continue
            else:
                plt.plot(all_chempots, hkl_area_dict[hkl],
                         '--', color=cmap(f[i]), label=str(hkl))

        # Make the figure look nice
        # ax2 = ax1.twiny()
        # ax2.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" %(str(self.comp1.elements[0])))
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        axes.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)"
                        %(str(self.ref_el_comp.elements[0])))

        ylim = axes.get_ylim()
        plt.xticks(rotation=60)
        plt.ylim(ylim)
        xlim = axes.get_xlim()
        plt.xlim(xlim)
        plt.tight_layout(pad=5, rect=[-0.0, 0, 0.95, 1])
        plt.plot([self.chempot_range[0], self.chempot_range[0]], ylim, '--k')
        plt.plot([self.chempot_range[1], self.chempot_range[1]], ylim, '--k')
        xy = [np.mean([self.chempot_range[1]]), np.mean(ylim)]
        plt.annotate("%s-rich" %(str(self.ref_el_comp.elements[0])), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)
        xy = [np.mean([xlim[0]]), np.mean(ylim)]
        plt.annotate("%s-poor" %(str(self.ref_el_comp.elements[0])), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)
        plt.ylabel(r"Fractional area $A^{Wulff}_{hkl}/A^{Wulff}$")

        return plt

    def chempot_vs_gamma_plot(self, cmap=cm.jet, JPERM2=False, highlight_stability={}):
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
        # box = axes.get_position()
        # axes.set_position([box.x0, box.y0, box.width * 0.7, box.height])
        # Choose unique colors for each facet
        f = [int(i) for i in np.linspace(0, 255, sum([len(vaspruns) for vaspruns in
                                                      self.vasprun_dict.values()]))]
        i, already_labelled, colors = 0, [], []
        for hkl in self.vasprun_dict.keys():
            for vasprun in self.vasprun_dict[hkl]:
                slab = vasprun.final_structure
                # Generate a label for the type of slab
                label = str(hkl)
                # use dashed lines for slabs that are not stoichiometric
                # wrt bulk. Label with formula if nonstoichiometric
                if slab.composition.reduced_composition != \
                        self.ucell_entry.composition.reduced_composition:
                    mark = '--'
                    label += " %s" % (slab.composition.reduced_composition)
                else:
                    mark = '-'

                # label the chemical environment at the surface if different from the bulk.
                # First get the surface sites, then get the reduced composition at the surface
                # s = vasprun.final_structure
                # ucell = SpacegroupAnalyzer(self.ucell_entry.structure).\
                #     get_conventional_standard_structure()
                # slab = Slab(s.lattice, s.species, s.frac_coords, hkl, ucell, 0, None)
                # surf_comp = slab.surface_composition()
                #
                # if surf_comp.reduced_composition != ucell.composition.reduced_composition:
                #     label += " %s" %(surf_comp.reduced_composition)

                if label in already_labelled:
                    c = colors[already_labelled.index(label)]
                    label = None
                else:
                    already_labelled.append(label)
                    c = cmap(f[i])
                    colors.append(c)

                gamma_range = [self.calculate_gamma_at_u(vasprun, self.chempot_range[0]),
                               self.calculate_gamma_at_u(vasprun, self.chempot_range[1])]
                se_range = np.array(gamma_range)*EV_PER_ANG2_TO_JOULES_PER_M2 \
                    if JPERM2 else gamma_range
                plt.plot(self.chempot_range, se_range, mark, color=c, label=label)
                i += 1

        for hkl in highlight_stability.keys():
            stab_map = self.get_stable_surf_regions(hkl)
            for i, s in enumerate(stab_map):
                label = "%s stability" % (str(hkl)) if i == 1 else None
                plt.plot(s[0], s[1], c=highlight_stability[hkl], linewidth=5, alpha=0.2,
                         label=label)

        # Make the figure look nice
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        plt.ylabel(r"Surface energy (J/$m^{2}$)") if JPERM2 \
            else plt.ylabel(r"Surface energy (eV/$\AA^{2}$)")
        # ax2 = axes.twiny()
        # ax2.tick_params(labelsize=20)
        # ax2.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" %(str(self.comp1.elements[0])), fontsize=24)
        axes.set_xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" %(str(self.ref_el_comp.elements[0])))

        ylim = axes.get_ylim()
        plt.xticks(rotation=60)
        plt.ylim(ylim)
        xlim = axes.get_xlim()
        plt.xlim(xlim)
        plt.tight_layout(pad=2.4, rect=[-0.047, 0, 0.84, 1])
        plt.plot([self.chempot_range[0], self.chempot_range[0]], ylim, '--k')
        plt.plot([self.chempot_range[1], self.chempot_range[1]], ylim, '--k')
        xy = [np.mean([self.chempot_range[1]]), np.mean(ylim)]
        plt.annotate("%s-rich" %(str(self.ref_el_comp.elements[0])), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)
        xy = [np.mean([xlim[0]]), np.mean(ylim)]
        plt.annotate("%s-poor" %(str(self.ref_el_comp.elements[0])), xy=xy,
                     xytext=xy, rotation=90, fontsize=17)

        return plt

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
