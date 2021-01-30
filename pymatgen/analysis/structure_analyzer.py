# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to perform topological analyses of structures.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman"
__copyright__ = "Copyright 2011, The Materials Project"

import collections
import itertools
from math import acos, pi
from warnings import warn

import numpy as np
from monty.dev import deprecated
from scipy.spatial import Voronoi

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.local_env import JmolNN, VoronoiNN
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.num import abs_cap


def average_coordination_number(structures, freq=10):
    """
    Calculates the ensemble averaged Voronoi coordination numbers
    of a list of Structures using VoronoiNN.
    Typically used for analyzing the output of a Molecular Dynamics run.

    Args:
        structures (list): list of Structures.
        freq (int): sampling frequency of coordination number [every freq steps].

    Returns:
        Dictionary of elements as keys and average coordination numbers as values.
    """
    coordination_numbers = {}
    for spec in structures[0].composition.as_dict().keys():
        coordination_numbers[spec] = 0.0
    count = 0
    for i, s in enumerate(structures):
        if i % freq != 0:
            continue
        count += 1
        vnn = VoronoiNN()
        for j, atom in enumerate(s):
            cn = vnn.get_cn(s, j, use_weights=True)
            coordination_numbers[atom.species_string] += cn
    elements = structures[0].composition.as_dict()
    for el in coordination_numbers:
        coordination_numbers[el] = coordination_numbers[el] / elements[el] / count
    return coordination_numbers


class VoronoiAnalyzer:
    """
    Performs a statistical analysis of Voronoi polyhedra around each site.
    Each Voronoi polyhedron is described using Schaefli notation.
    That is a set of indices {c_i} where c_i is the number of faces with i
    number of vertices.  E.g. for a bcc crystal, there is only one polyhedron
    notation of which is [0,6,0,8,0,0,...].
    In perfect crystals, these also corresponds to the Wigner-Seitz cells.
    For distorted-crystals, liquids or amorphous structures, rather than one-type,
    there is a statistical distribution of polyhedra.
    See ref: Microstructure and its relaxation in Fe-B amorphous system
    simulated by molecular dynamics,
        Stepanyuk et al., J. Non-cryst. Solids (1993), 159, 80-87.
    """

    def __init__(self, cutoff=5.0, qhull_options="Qbb Qc Qz"):
        """
        Args:
            cutoff (float): cutoff distance to search for neighbors of a given atom
                (default = 5.0)
            qhull_options (str): options to pass to qhull (optional)
        """
        self.cutoff = cutoff
        self.qhull_options = qhull_options

    def analyze(self, structure, n=0):
        """
        Performs Voronoi analysis and returns the polyhedra around atom n
        in Schlaefli notation.

        Args:
            structure (Structure): structure to analyze
            n (int): index of the center atom in structure

        Returns:
            voronoi index of n: <c3,c4,c6,c6,c7,c8,c9,c10>
                where c_i denotes number of facets with i vertices.
        """
        center = structure[n]
        neighbors = structure.get_sites_in_sphere(center.coords, self.cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        qvoronoi_input = np.array([s.coords for s in neighbors])
        voro = Voronoi(qvoronoi_input, qhull_options=self.qhull_options)
        vor_index = np.array([0, 0, 0, 0, 0, 0, 0, 0])

        for key in voro.ridge_dict:
            if 0 in key:  # This means if the center atom is in key
                if -1 in key:  # This means if an infinity point is in key
                    raise ValueError("Cutoff too short.")
                try:
                    vor_index[len(voro.ridge_dict[key]) - 3] += 1
                except IndexError:
                    # If a facet has more than 10 edges, it's skipped here.
                    pass
        return vor_index

    def analyze_structures(self, structures, step_freq=10, most_frequent_polyhedra=15):
        """
        Perform Voronoi analysis on a list of Structures.
        Note that this might take a significant amount of time depending on the
        size and number of structures.

        Args:
            structures (list): list of Structures
            cutoff (float: cutoff distance around an atom to search for
                neighbors
            step_freq (int): perform analysis every step_freq steps
            qhull_options (str): options to pass to qhull
            most_frequent_polyhedra (int): this many unique polyhedra with
                highest frequences is stored.

        Returns:
            A list of tuples in the form (voronoi_index,frequency)
        """
        voro_dict = {}
        step = 0
        for structure in structures:
            step += 1
            if step % step_freq != 0:
                continue

            v = []
            for n in range(len(structure)):
                v.append(str(self.analyze(structure, n=n).view()))
            for voro in v:
                if voro in voro_dict:
                    voro_dict[voro] += 1
                else:
                    voro_dict[voro] = 1
        return sorted(voro_dict.items(), key=lambda x: (x[1], x[0]), reverse=True)[:most_frequent_polyhedra]

    @staticmethod
    def plot_vor_analysis(voronoi_ensemble):
        """
        Plot the Voronoi analysis.

        :param voronoi_ensemble:
        :return: matplotlib.pyplot
        """
        t = list(zip(*voronoi_ensemble))
        labels = t[0]
        val = list(t[1])
        tot = np.sum(val)
        val = [float(j) / tot for j in val]
        pos = np.arange(len(val)) + 0.5  # the bar centers on the y axis
        import matplotlib.pyplot as plt

        plt.figure()
        plt.barh(pos, val, align="center", alpha=0.5)
        plt.yticks(pos, labels)
        plt.xlabel("Count")
        plt.title("Voronoi Spectra")
        plt.grid(True)
        return plt


class RelaxationAnalyzer:
    """
    This class analyzes the relaxation in a calculation.
    """

    def __init__(self, initial_structure, final_structure):
        """
        Please note that the input and final structures should have the same
        ordering of sites. This is typically the case for most computational
        codes.

        Args:
            initial_structure (Structure): Initial input structure to
                calculation.
            final_structure (Structure): Final output structure from
                calculation.
        """
        if final_structure.formula != initial_structure.formula:
            raise ValueError("Initial and final structures have different " + "formulas!")
        self.initial = initial_structure
        self.final = final_structure

    def get_percentage_volume_change(self):
        """
        Returns the percentage volume change.

        Returns:
            Volume change in percentage, e.g., 0.055 implies a 5.5% increase.
        """
        initial_vol = self.initial.lattice.volume
        final_vol = self.final.lattice.volume
        return final_vol / initial_vol - 1

    def get_percentage_lattice_parameter_changes(self):
        """
        Returns the percentage lattice parameter changes.

        Returns:
            A dict of the percentage change in lattice parameter, e.g.,
            {'a': 0.012, 'b': 0.021, 'c': -0.031} implies a change of 1.2%,
            2.1% and -3.1% in the a, b and c lattice parameters respectively.
        """
        initial_latt = self.initial.lattice
        final_latt = self.final.lattice
        d = {l: getattr(final_latt, l) / getattr(initial_latt, l) - 1 for l in ["a", "b", "c"]}
        return d

    def get_percentage_bond_dist_changes(self, max_radius=3.0):
        """
        Returns the percentage bond distance changes for each site up to a
        maximum radius for nearest neighbors.

        Args:
            max_radius (float): Maximum radius to search for nearest
               neighbors. This radius is applied to the initial structure,
               not the final structure.

        Returns:
            Bond distance changes as a dict of dicts. E.g.,
            {index1: {index2: 0.011, ...}}. For economy of representation, the
            index1 is always less than index2, i.e., since bonding between
            site1 and siten is the same as bonding between siten and site1,
            there is no reason to duplicate the information or computation.
        """
        data = collections.defaultdict(dict)
        for inds in itertools.combinations(list(range(len(self.initial))), 2):
            (i, j) = sorted(inds)
            initial_dist = self.initial[i].distance(self.initial[j])
            if initial_dist < max_radius:
                final_dist = self.final[i].distance(self.final[j])
                data[i][j] = final_dist / initial_dist - 1
        return data


class VoronoiConnectivity:
    """
    Computes the solid angles swept out by the shared face of the voronoi
    polyhedron between two sites.
    """

    def __init__(self, structure, cutoff=10):
        """
        Args:
            structure (Structure): Input structure
            cutoff (float) Cutoff distance.
        """
        self.cutoff = cutoff
        self.s = structure
        recp_len = np.array(self.s.lattice.reciprocal_lattice.abc)
        i = np.ceil(cutoff * recp_len / (2 * pi))
        offsets = np.mgrid[-i[0] : i[0] + 1, -i[1] : i[1] + 1, -i[2] : i[2] + 1].T
        self.offsets = np.reshape(offsets, (-1, 3))
        # shape = [image, axis]
        self.cart_offsets = self.s.lattice.get_cartesian_coords(self.offsets)

    @property
    def connectivity_array(self):
        """
        Provides connectivity array.

        Returns:
            connectivity: An array of shape [atomi, atomj, imagej]. atomi is
            the index of the atom in the input structure. Since the second
            atom can be outside of the unit cell, it must be described
            by both an atom index and an image index. Array data is the
            solid angle of polygon between atomi and imagej of atomj
        """
        # shape = [site, axis]
        cart_coords = np.array(self.s.cart_coords)
        # shape = [site, image, axis]
        all_sites = cart_coords[:, None, :] + self.cart_offsets[None, :, :]
        vt = Voronoi(all_sites.reshape((-1, 3)))
        n_images = all_sites.shape[1]
        cs = (len(self.s), len(self.s), len(self.cart_offsets))
        connectivity = np.zeros(cs)
        vts = np.array(vt.vertices)
        for (ki, kj), v in vt.ridge_dict.items():
            atomi = ki // n_images
            atomj = kj // n_images

            imagei = ki % n_images
            imagej = kj % n_images

            if imagei != n_images // 2 and imagej != n_images // 2:
                continue

            if imagei == n_images // 2:
                # atomi is in original cell
                val = solid_angle(vt.points[ki], vts[v])
                connectivity[atomi, atomj, imagej] = val

            if imagej == n_images // 2:
                # atomj is in original cell
                val = solid_angle(vt.points[kj], vts[v])
                connectivity[atomj, atomi, imagei] = val

            if -10.101 in vts[v]:
                warn("Found connectivity with infinite vertex. " "Cutoff is too low, and results may be " "incorrect")
        return connectivity

    @property
    def max_connectivity(self):
        """
        returns the 2d array [sitei, sitej] that represents
        the maximum connectivity of site i to any periodic
        image of site j
        """
        return np.max(self.connectivity_array, axis=2)

    def get_connections(self):
        """
        Returns a list of site pairs that are Voronoi Neighbors, along
        with their real-space distances.
        """
        con = []
        maxconn = self.max_connectivity
        for ii in range(0, maxconn.shape[0]):
            for jj in range(0, maxconn.shape[1]):
                if maxconn[ii][jj] != 0:
                    dist = self.s.get_distance(ii, jj)
                    con.append([ii, jj, dist])
        return con

    def get_sitej(self, site_index, image_index):
        """
        Assuming there is some value in the connectivity array at indices
        (1, 3, 12). sitei can be obtained directly from the input structure
        (structure[1]). sitej can be obtained by passing 3, 12 to this function

        Args:
            site_index (int): index of the site (3 in the example)
            image_index (int): index of the image (12 in the example)
        """
        atoms_n_occu = self.s[site_index].species
        lattice = self.s.lattice
        coords = self.s[site_index].frac_coords + self.offsets[image_index]
        return PeriodicSite(atoms_n_occu, coords, lattice)


def solid_angle(center, coords):
    """
    Helper method to calculate the solid angle of a set of coords from the
    center.

    Args:
        center (3x1 array): Center to measure solid angle from.
        coords (Nx3 array): List of coords to determine solid angle.

    Returns:
        The solid angle.
    """
    o = np.array(center)
    r = [np.array(c) - o for c in coords]
    r.append(r[0])
    n = [np.cross(r[i + 1], r[i]) for i in range(len(r) - 1)]
    n.append(np.cross(r[1], r[0]))
    vals = []
    for i in range(len(n) - 1):
        v = -np.dot(n[i], n[i + 1]) / (np.linalg.norm(n[i]) * np.linalg.norm(n[i + 1]))
        vals.append(acos(abs_cap(v)))
    phi = sum(vals)
    return phi + (3 - len(r)) * pi


def get_max_bond_lengths(structure, el_radius_updates=None):
    """
    Provides max bond length estimates for a structure based on the JMol
    table and algorithms.

    Args:
        structure: (structure)
        el_radius_updates: (dict) symbol->float to update atomic radii

    Returns: (dict) - (Element1, Element2) -> float. The two elements are
        ordered by Z.
    """
    # jmc = JMolCoordFinder(el_radius_updates)
    jmnn = JmolNN(el_radius_updates=el_radius_updates)

    bonds_lens = {}
    els = sorted(structure.composition.elements, key=lambda x: x.Z)

    for i1, el1 in enumerate(els):
        for i2 in range(len(els) - i1):
            bonds_lens[el1, els[i1 + i2]] = jmnn.get_max_bond_distance(el1.symbol, els[i1 + i2].symbol)

    return bonds_lens


@deprecated(
    message=(
        "find_dimension has been moved to"
        "pymatgen.analysis.dimensionality.get_dimensionality_gorai"
        " this method will be removed in pymatgen v2019.1.1."
    )
)
def get_dimensionality(
    structure,
    max_hkl=2,
    el_radius_updates=None,
    min_slab_size=5,
    min_vacuum_size=5,
    standardize=True,
    bonds=None,
):
    """
    This method returns whether a structure is 3D, 2D (layered), or 1D (linear
    chains or molecules) according to the algorithm published in Gorai, P.,
    Toberer, E. & Stevanovic, V. Computational Identification of Promising
    Thermoelectric Materials Among Known Quasi-2D Binary Compounds. J. Mater.
    Chem. A 2, 4136 (2016).

    Note that a 1D structure detection might indicate problems in the bonding
    algorithm, particularly for ionic crystals (e.g., NaCl)

    Users can change the behavior of bonds detection by passing either
    el_radius_updates to update atomic radii for auto-detection of max bond
    distances, or bonds to explicitly specify max bond distances for atom pairs.
    Note that if you pass both, el_radius_updates are ignored.

    Args:
        structure: (Structure) structure to analyze dimensionality for
        max_hkl: (int) max index of planes to look for layers
        el_radius_updates: (dict) symbol->float to update atomic radii
        min_slab_size: (float) internal surface construction parameter
        min_vacuum_size: (float) internal surface construction parameter
        standardize (bool): whether to standardize the structure before
            analysis. Set to False only if you already have the structure in a
            convention where layers / chains will be along low <hkl> indexes.
        bonds ({(specie1, specie2): max_bond_dist}: bonds are
                specified as a dict of tuples: float of specie1, specie2
                and the max bonding distance. For example, PO4 groups may be
                defined as {("P", "O"): 3}.

    Returns: (int) the dimensionality of the structure - 1 (molecules/chains),
        2 (layered), or 3 (3D)

    """
    if standardize:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()

    if not bonds:
        bonds = get_max_bond_lengths(structure, el_radius_updates)

    num_surfaces = 0
    for h in range(max_hkl):
        for k in range(max_hkl):
            for l in range(max_hkl):
                if max([h, k, l]) > 0 and num_surfaces < 2:
                    sg = SlabGenerator(
                        structure,
                        (h, k, l),
                        min_slab_size=min_slab_size,
                        min_vacuum_size=min_vacuum_size,
                    )
                    slabs = sg.get_slabs(bonds)
                    for _ in slabs:
                        num_surfaces += 1

    return 3 - min(num_surfaces, 2)


def contains_peroxide(structure, relative_cutoff=1.1):
    """
    Determines if a structure contains peroxide anions.

    Args:
        structure (Structure): Input structure.
        relative_cutoff: The peroxide bond distance is 1.49 Angstrom.
            Relative_cutoff * 1.49 stipulates the maximum distance two O
            atoms must be to each other to be considered a peroxide.

    Returns:
        Boolean indicating if structure contains a peroxide anion.
    """
    return oxide_type(structure, relative_cutoff) == "peroxide"


class OxideType:
    """
    Separate class for determining oxide type.
    """

    def __init__(self, structure, relative_cutoff=1.1):
        """
        Args:
            structure: Input structure.
            relative_cutoff: Relative_cutoff * act. cutoff stipulates the max.
                distance two O atoms must be from each other. Default value is
                1.1. At most 1.1 is recommended, nothing larger, otherwise the
                script cannot distinguish between superoxides and peroxides.
        """
        self.structure = structure
        self.relative_cutoff = relative_cutoff
        self.oxide_type, self.nbonds = self.parse_oxide()

    def parse_oxide(self):
        """
        Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.

        Returns:
            oxide_type (str): Type of oxide
            ozonide/peroxide/superoxide/hydroxide/None.
            nbonds (int): Number of peroxide/superoxide/hydroxide bonds in
            structure.
        """
        structure = self.structure
        relative_cutoff = self.relative_cutoff
        o_sites_frac_coords = []
        h_sites_frac_coords = []
        lattice = structure.lattice

        if isinstance(structure.composition.elements[0], Element):
            comp = structure.composition
        elif isinstance(structure.composition.elements[0], Species):
            elmap = collections.defaultdict(float)
            for site in structure:
                for species, occu in site.species.items():
                    elmap[species.element] += occu
            comp = Composition(elmap)
        if Element("O") not in comp or comp.is_element:
            return "None", 0

        for site in structure:
            syms = [sp.symbol for sp in site.species.keys()]
            if "O" in syms:
                o_sites_frac_coords.append(site.frac_coords)
            if "H" in syms:
                h_sites_frac_coords.append(site.frac_coords)

        if h_sites_frac_coords:
            dist_matrix = lattice.get_all_distances(o_sites_frac_coords, h_sites_frac_coords)
            if np.any(dist_matrix < relative_cutoff * 0.93):
                return (
                    "hydroxide",
                    len(np.where(dist_matrix < relative_cutoff * 0.93)[0]) / 2.0,
                )
        dist_matrix = lattice.get_all_distances(o_sites_frac_coords, o_sites_frac_coords)
        np.fill_diagonal(dist_matrix, 1000)
        is_superoxide = False
        is_peroxide = False
        is_ozonide = False
        if np.any(dist_matrix < relative_cutoff * 1.35):
            bond_atoms = np.where(dist_matrix < relative_cutoff * 1.35)[0]
            is_superoxide = True
        elif np.any(dist_matrix < relative_cutoff * 1.49):
            is_peroxide = True
            bond_atoms = np.where(dist_matrix < relative_cutoff * 1.49)[0]
        if is_superoxide:
            if len(bond_atoms) > len(set(bond_atoms)):
                is_superoxide = False
                is_ozonide = True
        try:
            nbonds = len(set(bond_atoms))
        except UnboundLocalError:
            nbonds = 0.0
        if is_ozonide:
            str_oxide = "ozonide"
        elif is_superoxide:
            str_oxide = "superoxide"
        elif is_peroxide:
            str_oxide = "peroxide"
        else:
            str_oxide = "oxide"
        if str_oxide == "oxide":
            nbonds = comp["O"]
        return str_oxide, nbonds


def oxide_type(structure, relative_cutoff=1.1, return_nbonds=False):
    """
    Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide

    Args:
        structure (Structure): Input structure.
        relative_cutoff (float): Relative_cutoff * act. cutoff stipulates the
            max distance two O atoms must be from each other.
        return_nbonds (bool): Should number of bonds be requested?
    """

    ox_obj = OxideType(structure, relative_cutoff)
    if return_nbonds:
        return ox_obj.oxide_type, ox_obj.nbonds
    return ox_obj.oxide_type


def sulfide_type(structure):
    """
    Determines if a structure is a sulfide/polysulfide/sulfate

    Args:
        structure (Structure): Input structure.

    Returns:
        (str) sulfide/polysulfide or None if structure is a sulfate.
    """
    structure = structure.copy()
    structure.remove_oxidation_states()
    s = Element("S")
    comp = structure.composition
    if comp.is_element or s not in comp:
        return None

    finder = SpacegroupAnalyzer(structure, symprec=0.1)
    symm_structure = finder.get_symmetrized_structure()
    s_sites = [sites[0] for sites in symm_structure.equivalent_sites if sites[0].specie == s]

    def process_site(site):

        # in an exceptionally rare number of structures, the search
        # radius needs to be increased to find a neighbor atom
        search_radius = 4
        neighbors = []
        while len(neighbors) == 0:
            neighbors = structure.get_neighbors(site, search_radius)
            search_radius *= 2
            if search_radius > max(structure.lattice.abc) * 2:
                break

        neighbors = sorted(neighbors, key=lambda n: n.nn_distance)
        dist = neighbors[0].nn_distance
        coord_elements = [nn.specie for nn in neighbors if nn.nn_distance < dist + 0.4][:4]
        avg_electroneg = np.mean([e.X for e in coord_elements])
        if avg_electroneg > s.X:
            return "sulfate"
        if avg_electroneg == s.X and s in coord_elements:
            return "polysulfide"
        return "sulfide"

    types = {process_site(site) for site in s_sites}
    if "sulfate" in types:
        return None
    if "polysulfide" in types:
        return "polysulfide"
    return "sulfide"
