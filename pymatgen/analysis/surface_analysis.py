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

from pymatgen.core.structure import Structure, Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import Slab
from pymatgen.analysis.wulff import WulffShape
from pymatgen import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen import Element
from pymatgen.util.plotting import pretty_plot


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

        Dictionary containing a list of Vaspruns for slab calculations as
            items and the corresponding Miller index of the slab as the key

    """

    def __init__(self, material_id, vasprun_dict, ref_element,
                 exclude_ids=[], custom_entries=[], mapi_key=None):
        """
        Analyzes surface energies and Wulff shape of a particular
            material using the chemical potential.
        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).
            vasprun_dict (dict): Dictionary containing a list of Vaspruns
                for slab calculations as items and the corresponding Miller
                index of the slab as the key.
                eg. vasprun_dict = {(1,1,1): [vasprun_111_1, vasprun_111_2,
                vasprun_111_3], (1,1,0): [vasprun_111_1, vasprun_111_2], ...}
            element: element to be considered as independent
                variables. E.g., if you want to show the stability
                ranges of all Li-Co-O phases wrt to uLi
            exclude_ids (list of material_ids): List of material_ids
                to exclude when obtaining the decomposition components
                to calculate the chemical potential
            custom_entries (list of pymatgen-db type entries): List of
                user specified pymatgen-db type entries to use in finding
                decomposition components for the chemical potential
            mapi_key (str): Materials Project API key for accessing the
                MP database via MPRester
        """

        self.ref_element = ref_element
        self.mprester = MPRester(mapi_key) if mapi_key else MPRester()
        self.ucell_entry = \
            self.mprester.get_entry_by_material_id(material_id,
                                                   inc_structure=True,
                                                   property_data=
                                                   ["formation_energy_per_atom"])
        ucell = self.ucell_entry.structure

        # Get x and y, the number of species in a formula unit of the bulk
        reduced_comp = ucell.composition.reduced_composition.as_dict()
        if len(reduced_comp.keys()) == 1:
            x = y = reduced_comp[ucell[0].species_string]
        else:
            for el in reduced_comp.keys():
                if self.ref_element == el:
                    y = reduced_comp[el]
                else:
                    x = reduced_comp[el]

        # Calculate Gibbs free energy of the bulk per unit formula
        gbulk = self.ucell_entry.energy /\
                (len([site for site in ucell
                      if site.species_string == self.ref_element]) / y)

        entries = [entry for entry in
                   self.mprester.get_entries_in_chemsys(list(reduced_comp.keys()),
                                                        property_data=["e_above_hull",
                                                                       "material_id"])
                   if entry.data["e_above_hull"] == 0 and
                   entry.data["material_id"] not in exclude_ids] \
            if not custom_entries else custom_entries

        pd = PhaseDiagram(entries)
        chempot_ranges = pd.get_chempot_range_map([Element(self.ref_element)])
        # If no chemical potential is found, we return u=0, eg.
        # for a elemental system, the relative u of Cu for Cu is 0
        chempot_range = [chempot_ranges[entry] for entry in chempot_ranges.keys()
                         if entry.composition ==
                         self.ucell_entry.composition][0][0]._coords if \
            chempot_ranges else [[0,0], [0,0]]

        e_of_element = [entry.energy_per_atom for entry in
                        entries if str(entry.composition.reduced_composition)
                        == self.ref_element + "1"][0]

        self.x = x
        self.y = y
        self.gbulk = gbulk
        chempot_range = list(chempot_range)
        self.chempot_range = sorted([chempot_range[0][0], chempot_range[1][0]])
        self.e_of_element = e_of_element
        self.vasprun_dict = vasprun_dict

    def calculate_gamma(self, vasprun):
        """
        Calculates the surface energy for a single slab.
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
            Ny = comp[self.ucell_entry.structure[0].species_string]
            Nx = Ny
        else:
            for el in reduced_comp.keys():
                if self.ref_element == el:
                    Ny = comp[el]
                else:
                    Nx = comp[el]

        # Calculate surface area
        m = slab.lattice.matrix
        A = np.linalg.norm(np.cross(m[0], m[1]))

        # calculate the surface energy for the max and min chemical potential
        return [(1 / (2 * A)) * (vasprun.final_energy - (Nx / self.x)
                                 * self.gbulk - (Ny - (self.y / self.x) * Nx)
                                 * (delu + self.e_of_element))
                for delu in self.chempot_range]

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
            surf_e_range_list = [self.calculate_gamma(vasprun)
                                 for vasprun in self.vasprun_dict[hkl]]
            e_list = []
            for e_range in surf_e_range_list:
                slope, intercept = self.get_slope_and_intercept(e_range)
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
            u_at_intersection = [self.get_intersections(hkl)[0] for hkl in
                                 self.vasprun_dict.keys()
                                 if self.get_intersections(hkl)]
            # Get a Wulff shape for each intersection. The change in the Wulff shape
            # will vary if the rate of change in surface energy for any facet changes
            for u in u_at_intersection:
                wulff_dict[u] = self.wulff_shape_from_chempot(u, symprec=symprec)

        return wulff_dict

    def get_slope_and_intercept(self, surf_e_pair):
        """
        Returns the slope and intercept of the surface
            energy vs chemical potential line
        Args:
            surf_e_pair ([e_at_min_u, e_at_max_u]): The surface energy at the
                minimum chemical potential and maximum chemical potential
        """

        slope, intercept, r_value, p_value, std_err = \
            linregress(self.chempot_range, surf_e_pair)
        slope = 0 if str(slope) == 'nan' else slope
        intercept = surf_e_pair[0] if str(intercept) == 'nan' else intercept
        return slope, intercept

    def get_intersections(self, miller_index):
        """
        Returns a all intersections for a specific facet. Useful for
            finding when the configuration of a particular facet changes.

        Args:
            miller_index ((h, k, l)): Miller index of the facet we
                are interested in
        """

        # First lets calculate the range of surface
        # energies for all terminations of a specific facet
        all_se_ranges = [self.calculate_gamma(vasprun) for vasprun
                         in self.vasprun_dict[miller_index]]

        if len(all_se_ranges) == 1:
            return []

        # Now get all possible intersection coordinates for each pair of lines
        intersections = []
        for pair_ranges in itertools.combinations(all_se_ranges, 2):
            slope1, intercept1 = self.get_slope_and_intercept(pair_ranges[0])
            slope2, intercept2 = self.get_slope_and_intercept(pair_ranges[1])
            # Calculate the intersection coordinates
            u = (intercept1-intercept2)/(slope2-slope1)
            # if the intersection is beyond the chemical potential
            # range or if the lines are parallel, we ignore it
            if slope1-slope2 == 0 or u < min(self.chempot_range) \
                    or u > max(self.chempot_range):
                continue
            intersections.append([u, slope1 * u + intercept1])

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
        plt = pretty_plot()
        for i, hkl in enumerate(self.vasprun_dict.keys()):
            # Ignore any facets that never show up on the
            # Wulff shape regardless of chemical potential
            if all([a == 0 for a in hkl_area_dict[hkl]]):
                continue
            else:
                plt.plot(all_chempots, hkl_area_dict[hkl],
                         '--', color=cmap(f[i]), label=str(hkl))

        # Make the figure look nice
        plt.ylim([0,1])
        plt.xlim(self.chempot_range)
        plt.ylabel(r"Fractional area $A^{Wulff}_{hkl}/A^{Wulff}$")
        plt.xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" %(self.ref_element))
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

        return plt

    def chempot_vs_gamma_plot(self, cmap=cm.jet, show_unstable_points=False):
        """
        Plots the surface energy of all facets as a function of chemical potential.
            Each facet will be associated with its own distinct colors. Dashed lines
            will represent stoichiometries different from that of the mpid's compound.

        Args:
            cmap (cm): A matplotlib colormap object, defaults to jet.
            show_unstable_points (bool): For each facet, there may be various
                terminations or stoichiometries and the relative stability of
                these different slabs may change with chemical potential. This
                option will only plot the most stable surface energy for a
                given chemical potential.
        """

        plt = pretty_plot()
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

                se_range = self.calculate_gamma(vasprun)
                plt.plot(self.chempot_range, se_range, mark, color=c, label=label)
                i += 1

        # Make the figure look nice
        axes = plt.gca()
        ylim = axes.get_ylim()
        plt.ylim(ylim)
        plt.xlim(self.chempot_range)
        plt.ylabel(r"Surface energy (eV/$\AA$)")
        plt.xlabel(r"Chemical potential $\Delta\mu_{%s}$ (eV)" %(self.ref_element))
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

        return plt

    def broken_bond_vs_gamma(self):

        return
