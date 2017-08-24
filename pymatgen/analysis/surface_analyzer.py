from __future__ import division, unicode_literals

import numpy as np
from pymatgen import Structure
from pymatgen.core.surface import SlabGenerator, Slab
import copy
from pymatgen.io.vasp.sets import MVLSlabSet
from pymacy.surface_adsorption.surface_updater import UpdateRepositoriesAndDBs, MPRester

from matplotlib import pylab as plt
from scipy.optimize import curve_fit
from pymatgen.analysis.wulff import WulffShape

__author__ = "Richard Tran"
__copyright__ = "Copyright 2014, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "8/24/17"


class SurfaceEnergyAnalyzer(object):

    def __init__(self, material_id, mapi_key):
        """
        Analyzes surface energy
        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).
            mapi_key (str): Materials Project API key for accessing the
                MP database via MPRester
        """

        self.material_id = material_id
        self.mprester = MPRester(mapi_key)

    def calculate_gamma(self, vasprun, element, exclude_ids=[], custom_entries=[]):
        """
        Calculates the surface energy for a single slab.
        Args:
            vasprun (Vasprun): A Vasprun object
            element: element to be considered as independent
                variables. E.g., if you want to show the stability
                ranges of all Li-Co-O phases wrt to uLi
            exclude_ids (list of material_ids): List of material_ids
                to exclude when obtaining the decomposition components
                to calculate the chemical potential
            custom_entries (list of pymatgen-db type entries): List of
                user specified pymatgen-db type entries to use in finding
                decomposition components for the chemical potential
        """

        ucell_entry = self.mprester.get_entry_by_material_id(self.material_id, inc_structure=True,
                                                             property_data=["formation_energy_per_atom"])
        ucell = ucell_entry.structure

        # Calculate surface area
        m = vasprun.structure.lattice.matrix
        A = np.linalg.norm(np.cross(m[0], m[1]))

        # Get x and y, the number of species in a formula unit of the bulk
        reduced_comp = ucell.composition.reduced_composition.as_dict()
        for el in reduced_comp.keys():
            if element == el:
                y = reduced_comp[el]
            else:
                x = reduced_comp[el]

        gbulk = ucell_entry.energy / (len([site for site in ucell
                                           if site.species_string == element]) / y)

        entries = [entry for entry in
                   self.mprester.get_entries_in_chemsys(list(reduced_comp.keys()),
                                                        property_data=["e_above_hull",
                                                                       "material_id"])
                   if entry.data["e_above_hull"] == 0 and
                   entry.data["material_id"] not in exclude_ids] \
            if not custom_entries else custom_entries

        pd = PhaseDiagram(entries)
        pda = PDAnalyzer(pd)
        chempot_ranges = pda.get_chempot_range_map([Element(element)])

        chempot_range = [chempot_ranges[entry] for entry in chempot_ranges.keys()
                         if entry.composition == ucell_entry.composition][0]
        print("range: ", chempot_range)
        print(chempot_range[0])
        print("gbulk", gbulk)

        # Get the composition in the slab
        comp = vasprun.structure.composition.as_dict()
        for el in reduced_comp.keys():
            if element == el:
                Ny = comp[el]
            else:
                Nx = comp[el]

        e_of_element = [entry.energy_per_atom for entry in
                        entries if str(entry.composition.reduced_composition)
                        == element + "1"][0]
        print(e_of_element)

        return [(1 / (2 * A)) * (vasprun.final_energy - (Nx / x) * gbulk - \
                                 (Ny - (y / x) * Nx) * (delu + e_of_element))
                for delu in chempot_range[0]._coords]


    def build_wulff_shape(self):

        return

    def chemical_potential_gamma_plot(self):

        return

    def broken_bond_gamma_plot(self):

        return

    def broken_bond_phi_plot(self):

        return

    def calculate_weighted_phi(self):

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