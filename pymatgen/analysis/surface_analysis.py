# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import copy

import numpy as np
from scipy.stats import linregress
from matplotlib import cm
import itertools

from pymatgen.core.structure import Structure, Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import Slab
from pymatgen.analysis.wulff import WulffShape
from pymatgen import MPRester
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
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
        material_id. By default, this will use entries calculated from the Materials
        Project to obtain chemical potential and bulk energy. As a result, the
        difference in VASP parameters between the user's entry (vasprun_dict) and the
        parameters used by Materials Project, may lead to a rough estimate of the
        surface energy. For best results, it is recommend that the user calculates all
        decomposition components first, and insert the results into their own database
        as a pymatgen-db entry and use those entries instead (custom_entries). In
        addition, this code will only use one bulk entry to calculate surface energy.
        Ideally, to get the most accurate surface energy, the user should compare their
        slab energy to the energy of the oriented unit cell with both calculations
        containing consistent k-points to avoid converegence problems as the slab size
        is varied. See:
            Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
                Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.
        and
            Rogal, J., & Reuter, K. (2007). Ab Initio Atomistic Thermodynamics for
                Surfaces : A Primer. Experiment, Modeling and Simulation of Gas-Surface
                Interactions for Reactive Flows in Hypersonic Flights, 2–1 – 2–18.

    .. attribute:: ref_element

        All chemical potentials cna be written in terms of the range of chemical potential
            of this element which will be used to calculate surface energy.

    .. attribute:: mprester

        Materials project rester for querying entries from the materials project. Requires
            user MAPIKEY.

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
            y = reduced_comp[ucell[0].species_string]
            x = y
        else:
            for el in reduced_comp.keys():
                if self.ref_element == el:
                    y = reduced_comp[el]
                else:
                    x = reduced_comp[el]

        gbulk = self.ucell_entry.energy / (len([site for site in ucell
                                                if site.species_string == self.ref_element]) / y)

        entries = [entry for entry in
                   self.mprester.get_entries_in_chemsys(list(reduced_comp.keys()),
                                                        property_data=["e_above_hull",
                                                                       "material_id"])
                   if entry.data["e_above_hull"] == 0 and
                   entry.data["material_id"] not in exclude_ids] \
            if not custom_entries else custom_entries

        pd = PhaseDiagram(entries)
        pda = PDAnalyzer(pd)
        chempot_ranges = pda.get_chempot_range_map([Element(self.ref_element)])
        # If no chemical potential is found, we return u=0, eg.
        # for a elemental system, the relative u of Cu for Cu is 0
        chempot_range = [chempot_ranges[entry] for entry in chempot_ranges.keys()
                         if entry.composition == self.ucell_entry.composition][0][0]._coords if \
            chempot_ranges else [[0,0], [0,0]]

        e_of_element = [entry.energy_per_atom for entry in
                        entries if str(entry.composition.reduced_composition)
                        == self.ref_element + "1"][0]

        self.x = x
        self.y = y
        self.gbulk = gbulk
        chempot_range = list(chempot_range)
        self.chempot_range = [chempot_range[0][0], chempot_range[1][0]]
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
        comp = vasprun.final_structure.composition.as_dict()
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
        m = vasprun.final_structure.lattice.matrix
        A = np.linalg.norm(np.cross(m[0], m[1]))

        # calculate the surface energy for the max and min chemical potential
        return [(1 / (2 * A)) * (vasprun.final_energy - (Nx / self.x)
                                 * self.gbulk - (Ny - (self.y / self.x) * Nx)
                                 * (delu + self.e_of_element))
                for delu in self.chempot_range]

    def get_wulff_shape_dict(self, symprec=1e-5, sample_intersections=False):
        """
        As the surface energy is a function of chemical potential, so too is the
            Wulff shape. This methods generates a dictionary of Wulff shapes with
            the keys being the chemical potential and the value is the
            corresponding Wulff shape.

        Args:
            symprec (float): for recp_operation, default is 1e-5.
            sample_intersections (bool): Whether to generate a Wulff shape for each
                intersection of surface energy for a specific facet (eg. at the point
                where a (111) stoichiometric surface energy plot intersects with the
                (111) nonstoichiometric plot) or to just generate two Wulff shapes,
                one at the min and max chemical potential.

        """

        miller_list = self.vasprun_dict.keys()
        surf_e_list = [self.calculate_min_gamma_hkl(hkl) for
                       hkl in miller_list]

        ucell = SpacegroupAnalyzer(self.ucell_entry.structure).\
            get_conventional_standard_structure()
        wulffshape = WulffShape(ucell.lattice, miller_list,
                                surf_e_list, symprec=symprec)

        return wulffshape

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
            slope1, intercept1, r_value, p_value, std_err = \
                linregress(self.chempot_range, pair_ranges[0])
            slope1 = 0 if str(slope1) == 'nan' else slope1
            intercept1 = 0 if str(intercept1) == 'nan' else intercept1
            slope2, intercept2, r_value, p_value, std_err = \
                linregress(self.chempot_range, pair_ranges[1])
            slope2 = 0 if str(slope2) == 'nan' else slope2
            intercept2 = 0 if str(intercept2) == 'nan' else intercept2
            # Calculate the intersection coordinates
            u = (intercept1-intercept2)/(slope2-slope1)
            # if the intersection is beyond the chemical potential
            # range or if the lines are parallel, we ignore it
            if slope1-slope2 == 0 or u < min(self.chempot_range) \
                    or u > max(self.chempot_range):
                continue
            intersections.append([u, slope1 * u + intercept1])

        return sorted(intersections, key=lambda ints: ints[0])

    def area_frac_vs_chempot_plot(self):

        return

    def chempot_vs_gamma_plot(self, cmap=cm.jet, show_unstable_points=False):
        """
        Plots the surface energy of all facets as a function of chemical potential.
            Each facet will be associated with its own distinct colors. Dashed lines
            will represent stoichiometries different from that of the mpid's compound.

        Args:
            cmap (cm): A matplotlib colormap object, defaults to jet.
            show_unstable_points (bool): For each facet, there may be various terminations
                or stoichiometries and the relative stability of these different slabs may
                change with chemical potential. This option will only plot the most stable
                surface energy for a given chemical potential.
        """

        plt = pretty_plot()
        # Choose unique colors for each facet
        f = [int(i) for i in np.linspace(0, 255, sum([len(vaspruns) for vaspruns in
                                                      self.vasprun_dict.values()]))]
        i, already_labelled, colors = 0, [], []
        for hkl in self.vasprun_dict.keys():
            for vasprun in self.vasprun_dict[hkl]:
                # Generate a label for the type of slab
                label = str(hkl)
                # use dashed lines for slabs that are not stoichiometric
                # wrt bulk. Label with formula if nonstoichiometric
                if vasprun.final_structure.composition.reduced_composition != \
                    self.ucell_entry.composition.reduced_composition:
                    mark = '--'
                    label += " %s" % (vasprun.final_structure.composition.reduced_composition)
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

# class SlabAnalyzer(object):
#     def __init__(self):
#
#         None
#
#     def bulk_coordination(self, slab, bondlength, bond):
#
#         # -If bond is list of 2 species, we find coordination for
#         #   the first specie relative to the second only [s1, s2]
#         # -If list of one species, find coordination for all neighbors [s]
#         # -If first item is a list, we are looking for
#         #   the coordination of a polyhedron [[s1a, s1b], s2]
#         # IMPORTANT NOTE, cannot specify the specific bondlength of your polyhedron
#         # if you are looking for the coordination of a polyhedron. The bondlength
#         # will be the same as that of the polyhedron vertex and the next species
#
#         center_ion = bond[0] if type(bond[0]).__name__ == 'str' else bond[0][0]
#         mean_cn = []
#         ucell = slab.oriented_unit_cell
#         for i, site in enumerate(ucell):
#             cn = 0
#             if site.species_string == center_ion:
#                 nn = ucell.get_neighbors(site, bondlength,
#                                          include_index=True)
#
#                 for n in nn:
#                     # If we're dealing with the coordination of a single atoms
#                     if type(bond[0]).__name__ == 'str':
#                         if len(bond) == 2:
#                             if n[0].species_string == bond[1]:
#                                 cn += 1
#                         else:
#                             cn += 1
#                     # If we're dealing with the coordination of a polyhedron
#                     else:
#                         # Check if the vertex of the polyhedron is the correct species
#                         if n[0].species_string == bond[0][1]:
#                             # Get coordinated sites with that vertex
#                             vert_n = ucell.get_neighbors(ucell[n[2]], bondlength,
#                                                          include_index=True)
#                             for nnn in vert_n:
#                                 # Check the coordinated site of vertex is
#                                 # not the center of the polyhedron
#                                 if nnn[2] == i:
#                                     continue
#                                 if len(bond) == 2:
#                                     if nnn[0].species_string == bond[1]:
#                                         cn += 1
#                                 else:
#                                     cn += 1
#
#                 mean_cn.append(cn)
#
#         return min(mean_cn)
#
#
#     def surface_coordination(self, slab, bonds, top=True):
#
#         """
#         A function that analyzes the coordination environment of bulk atoms, surface atoms
#             and broken bonds. Returns a dictionary describing each type of environment for
#             each type of bond: eg. {"bulk": {(species1, species2): 12}, "surface": {(species1,
#             species2): 6}, "broken": {(species1, species2): 3}}
#
#         Args:
#             slab (Slab): Initial input slab.
#             bonds ({(specie1, specie2): max_bond_dist}: bonds are
#                 specified as a dict of tuples: float of specie1, specie2
#                 and the max bonding distance. For example, PO4 groups may be
#                 defined as {("P", "O"): 3}.
#         """
#
#         bulk_bonds, surface_bonds, broken_bonds = {}, {}, {}
#         for bond in bonds.keys():
#
#             # First we check the cn of the bulk for each type of bond
#             cn = bulk_coordination(slab, bonds[bond], bond)
#
#             # Next we use the cn of the bulk as
#             # reference to find the number of broken bonds
#             center_ion = bond[0] if type(bond[0]).__name__ == 'str' else bond[0][0]
#             cnb, tot_surf_cn, bb_cn = 0, 0, 0
#             for i, site in enumerate(slab):
#                 if str(site.specie) == center_ion:
#                     nn = slab.get_neighbors(site, bonds[bond],
#                                             include_index=True)
#
#                     def count_cnb(cnb, tot_surf_cn, bb_cn):
#                         slab_cn = 0
#                         for n in nn:
#
#                             # If we're dealing with the coordination of a single atoms
#                             if type(bond[0]).__name__ == 'str':
#                                 if len(bond) == 2:
#                                     if n[0].species_string == bond[1]:
#                                         slab_cn += 1
#                                 else:
#                                     slab_cn += 1
#
#                             # If we're dealing with the coordination of a polyhedron
#                             else:
#                                 # Check if the vertex of the polyhedron is the correct species
#                                 if n[0].species_string == bond[0][1]:
#                                     # Get coordinated sites with that vertex
#                                     vert_n = slab.get_neighbors(slab[n[2]], bonds[bond],
#                                                                 include_index=True)
#                                     for nnn in vert_n:
#                                         # Check the coordinated site of vertex is
#                                         # not the center of the polyhedron
#                                         if nnn[2] == i:
#                                             continue
#                                         if len(bond) == 2:
#                                             if nnn[0].species_string == bond[1]:
#                                                 slab_cn += 1
#                                         else:
#                                             slab_cn += 1
#
#                         cnb += cn
#                         tot_surf_cn += slab_cn
#                         bb_cn += cn - slab_cn
#                         return cnb, tot_surf_cn, bb_cn
#
#                     if top and site.frac_coords[2] > slab.center_of_mass:
#                         cnb, tot_surf_cn, bb_cn = count_cnb(cnb, tot_surf_cn, bb_cn)
#                     if not top and site.frac_coords[2] < slab.center_of_mass:
#                         cnb, tot_surf_cn, bb_cn = count_cnb(cnb, tot_surf_cn, bb_cn)
#
#             bulk_bonds[bond] = cnb / slab.surface_area
#             surface_bonds[bond] = tot_surf_cn / slab.surface_area
#             broken_bonds[bond] = bb_cn / slab.surface_area
#
#         return {"bulk": bulk_bonds, "surface": surface_bonds, "broken": broken_bonds}



# def get_cn_dict_ouc(slab, decimal=5):
#     """
#     For each species, get all unique coordination numbers.
#     Returns a dictionary with the species_string as key
#     and the value as a list of cn
#
#         Args:
#             decimal (int): The decimal place to determine a unique
#                 coordination number
#     """
#
#     cn_dict = {}
#     ucell = slab.oriented_unit_cell
#     v = VoronoiCoordFinder(ucell)
#     for i, site in enumerate(ucell):
#         if ucell[i].species_string not in cn_dict.keys():
#             cn_dict[ucell[i].species_string] = []
#         # Since this will get the cn as a result of the weighted
#         # polyhedra, the slightest difference in cn will indicate
#         # a different coordination environment for a species, eg.
#         # bond distance of each neighbor or neighbor species.
#         cn = v.get_coordination_number(i)
#         cn = round(cn, decimal)
#         if cn not in cn_dict[ucell[i].species_string]:
#             cn_dict[ucell[i].species_string].append(cn)
#
#     return cn_dict
#
#
# def get_surface_sites(slab, decimal=5, top=True):
#     """
#     Returns a list of sites with undercoordination. Due to
#     pathological error resulting from some surface sites, we
#     assume any site that has this error is a surface site as well
#
#         Args:
#             decimal (int): The decimal place to determine a unique
#                 coordination number
#     """
#
#     cn_dict = get_cn_dict_ouc(slab)
#
#     v = VoronoiCoordFinder(slab)
#     surface_sites = []
#     for i, site in enumerate(slab):
#         if site.frac_coords[2] > 0.5 and not top:
#             continue
#         try:
#             cn = round(v.get_coordination_number(i), decimal)
#             if cn not in cn_dict[site.species_string]:
#                 surface_sites.append([site, i])
#         except RuntimeError:
#             surface_sites.append([site, i])
#
#     return surface_sites
#
#
# def are_surfaces_equal(slab):
#     # Check if we have same number of equivalent sites on both surfaces.
#     # This is an alternative to checking Laue symmetry if we want to
#     # ensure both surfaces in the slab are the same
#
#     # create a structure isolating the
#     # two surfaces from the rest of the slab
#     surfsites = get_surface_sites(slab)
#     species = [s[0].specie for s in surfsites]
#     coords = [s[0].frac_coords for s in surfsites]
#     surface = Structure(slab.lattice, species, coords)
#     a = SpacegroupAnalyzer(surface)
#     symm_structure = a.get_symmetrized_structure()
#
#     # ensure each site on one surface has a
#     # corresponding equivalent site on the other
#     equal_surf_sites, total = [], 0
#     for equ in symm_structure.equivalent_sites:
#         top, bottom = 0, 0
#         for s in equ:
#             total += 1
#             if s.frac_coords[2] > 0.5:
#                 top += 1
#             else:
#                 bottom += 1
#         equal_surf_sites.append(top == bottom)
#     equal_surf_sites.append(total == len(surface))
#
#     return all(equal_surf_sites)