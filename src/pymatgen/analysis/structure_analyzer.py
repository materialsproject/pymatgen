"""This module provides classes to perform topological analyses of structures."""

from __future__ import annotations

import itertools
from collections import defaultdict
from math import acos, pi
from typing import TYPE_CHECKING
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi

from pymatgen.analysis.local_env import JmolNN, VoronoiNN
from pymatgen.core import Composition, Element, PeriodicSite, Species
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from pymatgen.core import Structure

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman"
__copyright__ = "Copyright 2011, The Materials Project"


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
    for spec in structures[0].composition.as_dict():
        coordination_numbers[spec] = 0.0
    count = 0
    for idx, site in enumerate(structures):
        if idx % freq != 0:
            continue
        count += 1
        vnn = VoronoiNN()
        for j, atom in enumerate(site):
            cn = vnn.get_cn(site, j, use_weights=True)
            coordination_numbers[atom.species_string] += cn
    elements = structures[0].composition.as_dict()
    return {el: v / elements[el] / count for el, v in coordination_numbers.items()}


class VoronoiAnalyzer:
    """
    Performs a statistical analysis of Voronoi polyhedra around each site.
    Each Voronoi polyhedron is described using Schaefli notation.
    That is a set of indices {c_i} where c_i is the number of faces with i
    number of vertices. E.g. for a bcc crystal, there is only one polyhedron
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
            qhull_options (str): options to pass to qhull (optional).
        """
        self.cutoff = cutoff
        self.qhull_options = qhull_options

    def analyze(self, structure: Structure, n=0):
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
                highest frequencies is stored.

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
    def plot_vor_analysis(voronoi_ensemble: list[tuple[str, float]]) -> plt.Axes:
        """Plot the Voronoi analysis.

        Args:
            voronoi_ensemble (list[tuple[str, float]]): List of tuples containing labels and
                values for Voronoi analysis.

        Returns:
            plt.Axes: Matplotlib Axes object with the plotted Voronoi analysis.
        """
        labels, val = zip(*voronoi_ensemble, strict=True)
        arr = np.array(val, dtype=float)
        arr /= np.sum(arr)
        pos = np.arange(len(arr)) + 0.5  # the bar centers on the y axis

        _fig, ax = plt.subplots()
        ax.barh(pos, arr, align="center", alpha=0.5)
        ax.set_yticks(pos)
        ax.set_yticklabels(labels)
        ax.set(title="Voronoi Spectra", xlabel="Count")
        ax.grid(visible=True)
        return ax


class RelaxationAnalyzer:
    """This class analyzes the relaxation in a calculation."""

    def __init__(self, initial_structure: Structure, final_structure: Structure) -> None:
        """Please note that the input and final structures should have the same
        ordering of sites. This is typically the case for most computational codes.

        Args:
            initial_structure (Structure): Initial input structure to
                calculation.
            final_structure (Structure): Final output structure from
                calculation.

        Raises:
            ValueError: If initial and final structures have different formulas.
        """
        if final_structure.formula != initial_structure.formula:
            raise ValueError("Initial and final structures have different formulas!")
        self.initial = initial_structure
        self.final = final_structure

    def get_percentage_volume_change(self) -> float:
        """Get the percentage volume change.

        Returns:
            float: Volume change in percent. 0.055 means a 5.5% increase.
        """
        return self.final.volume / self.initial.volume - 1

    def get_percentage_lattice_parameter_changes(self) -> dict[str, float]:
        """Get the percentage lattice parameter changes.

        Returns:
            dict[str, float]: Percent changes in lattice parameter, e.g.
                {'a': 0.012, 'b': 0.021, 'c': -0.031} implies a change of 1.2%,
                2.1% and -3.1% in the a, b and c lattice parameters respectively.
        """
        initial_latt = self.initial.lattice
        final_latt = self.final.lattice
        return {length: getattr(final_latt, length) / getattr(initial_latt, length) - 1 for length in ["a", "b", "c"]}

    def get_percentage_bond_dist_changes(self, max_radius: float = 3.0) -> dict[int, dict[int, float]]:
        """Get the percentage bond distance changes for each site up to a
        maximum radius for nearest neighbors.

        Args:
            max_radius (float): Maximum radius to search for nearest
               neighbors. This radius is applied to the initial structure,
               not the final structure.

        Returns:
            dict[int, dict[int, float]]: Bond distance changes in the form {index1: {index2: 0.011, ...}}.
                For economy of representation, the index1 is always less than index2, i.e., since bonding
                between site1 and site_n is the same as bonding between site_n and site1, there is no
                reason to duplicate the information or computation.
        """
        data: dict[int, dict[int, float]] = defaultdict(dict)
        for indices in itertools.combinations(list(range(len(self.initial))), 2):
            ii, jj = sorted(indices)
            initial_dist = self.initial[ii].distance(self.initial[jj])
            if initial_dist < max_radius:
                final_dist = self.final[ii].distance(self.final[jj])
                data[ii][jj] = final_dist / initial_dist - 1
        return data


class VoronoiConnectivity:
    """
    Computes the solid angles swept out by the shared face of the voronoi
    polyhedron between two sites.
    """

    def __init__(self, structure: Structure, cutoff=10):
        """
        Args:
            structure (Structure): Input structure
            cutoff (float) Cutoff distance.
        """
        self.cutoff = cutoff
        self.structure = structure
        recip_vec = np.array(self.structure.lattice.reciprocal_lattice.abc)
        cutoff_vec = np.ceil(cutoff * recip_vec / (2 * pi))
        offsets = np.mgrid[
            -cutoff_vec[0] : cutoff_vec[0] + 1,
            -cutoff_vec[1] : cutoff_vec[1] + 1,
            -cutoff_vec[2] : cutoff_vec[2] + 1,
        ].T
        self.offsets = np.reshape(offsets, (-1, 3))
        # shape = [image, axis]
        self.cart_offsets = self.structure.lattice.get_cartesian_coords(self.offsets)

    @property
    def connectivity_array(self):
        """The connectivity array of shape [atom_i, atom_j, image_j]. atom_i is the index of the
        atom in the input structure. Since the second atom can be outside of the unit cell, it
        must be described by both an atom index and an image index. Array data is the solid
        angle of polygon between atom_i and image_j of atom_j.
        """
        # shape = [site, axis]
        cart_coords = np.array(self.structure.cart_coords)
        # shape = [site, image, axis]
        all_sites = cart_coords[:, None, :] + self.cart_offsets[None, :, :]
        vt = Voronoi(all_sites.reshape((-1, 3)))
        n_images = all_sites.shape[1]
        cs = (len(self.structure), len(self.structure), len(self.cart_offsets))
        connectivity = np.zeros(cs)
        vts = np.array(vt.vertices)
        for (ki, kj), v in vt.ridge_dict.items():
            atom_i = ki // n_images
            atom_j = kj // n_images

            image_i = ki % n_images
            image_j = kj % n_images

            if image_i != n_images // 2 and image_j != n_images // 2:
                continue

            if image_i == n_images // 2:
                # atom_i is in original cell
                val = solid_angle(vt.points[ki], vts[v])
                connectivity[atom_i, atom_j, image_j] = val

            if image_j == n_images // 2:
                # atom_j is in original cell
                val = solid_angle(vt.points[kj], vts[v])
                connectivity[atom_j, atom_i, image_i] = val

            if -10.101 in vts[v]:
                warn("Found connectivity with infinite vertex. Cutoff is too low, and results may be incorrect")
        return connectivity

    @property
    def max_connectivity(self):
        """The 2d array [site_i, site_j] that represents the maximum connectivity of
        site i to any periodic image of site j.
        """
        return np.max(self.connectivity_array, axis=2)

    def get_connections(self):
        """Get a list of site pairs that are Voronoi Neighbors, along
        with their real-space distances.
        """
        con = []
        max_conn = self.max_connectivity
        for ii in range(max_conn.shape[0]):
            for jj in range(max_conn.shape[1]):
                if max_conn[ii][jj] != 0:
                    dist = self.structure.get_distance(ii, jj)
                    con.append([ii, jj, dist])
        return con

    def get_sitej(self, site_index, image_index):
        """
        Assuming there is some value in the connectivity array at indices
        (1, 3, 12). site_i can be obtained directly from the input structure
        (structure[1]). site_j can be obtained by passing 3, 12 to this function.

        Args:
            site_index (int): index of the site (3 in the example)
            image_index (int): index of the image (12 in the example)
        """
        atoms_n_occu = self.structure[site_index].species
        lattice = self.structure.lattice
        coords = self.structure[site_index].frac_coords + self.offsets[image_index]
        return PeriodicSite(atoms_n_occu, coords, lattice)


def solid_angle(center, coords):
    """
    Helper method to calculate the solid angle of a set of coords from the center.

    Args:
        center (3x1 array): Center to measure solid angle from.
        coords (Nx3 array): List of coords to determine solid angle.

    Returns:
        float: The solid angle.
    """
    origin = np.array(center)
    radii = [np.array(c) - origin for c in coords]
    radii.append(radii[0])
    cross_products = [np.cross(radii[i + 1], radii[i]) for i in range(len(radii) - 1)]
    cross_products.append(np.cross(radii[1], radii[0]))
    vals = []
    for i in range(len(cross_products) - 1):
        v = -np.dot(cross_products[i], cross_products[i + 1]) / (
            np.linalg.norm(cross_products[i]) * np.linalg.norm(cross_products[i + 1])
        )
        vals.append(acos(np.clip(v, -1, 1)))
    phi = sum(vals)
    return phi + (3 - len(radii)) * pi


def get_max_bond_lengths(structure, el_radius_updates=None):
    """
    Provides max bond length estimates for a structure based on the JMol
    table and algorithms.

    Args:
        structure: (structure)
        el_radius_updates: (dict) symbol->float to update atom_ic radii

    Returns:
        dict[(Element1, Element2)], float]: The two elements are ordered by Z.
    """
    # jmc = JMolCoordFinder(el_radius_updates)
    jm_nn = JmolNN(el_radius_updates=el_radius_updates)

    bonds_lens = {}
    els = sorted(structure.elements, key=lambda x: x.Z)

    for i1, el1 in enumerate(els):
        for i2 in range(len(els) - i1):
            bonds_lens[el1, els[i1 + i2]] = jm_nn.get_max_bond_distance(el1.symbol, els[i1 + i2].symbol)

    return bonds_lens


def contains_peroxide(structure, relative_cutoff=1.1):
    """
    Determines if a structure contains peroxide anions.

    Args:
        structure (Structure): Input structure.
        relative_cutoff: The peroxide bond distance is 1.49 Angstrom.
            Relative_cutoff * 1.49 stipulates the maximum distance two O
            atoms must be to each other to be considered a peroxide.

    Returns:
        bool: True if structure contains a peroxide anion.
    """
    return oxide_type(structure, relative_cutoff) == "peroxide"


class OxideType:
    """Separate class for determining oxide type."""

    def __init__(self, structure: Structure, relative_cutoff=1.1):
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

    def parse_oxide(self) -> tuple[str, int]:
        """
        Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.

        Returns:
            tuple[str, int]: Type of oxide (ozonide/peroxide/superoxide/hydroxide/None) and number of
                peroxide/superoxide/hydroxide bonds in structure.
        """
        structure = self.structure
        relative_cutoff = self.relative_cutoff
        o_sites_frac_coords = []
        h_sites_frac_coords = []
        lattice = structure.lattice

        if isinstance(structure.elements[0], Element):
            comp = structure.composition
        elif isinstance(structure.elements[0], Species):
            elem_map: dict[Element, float] = defaultdict(float)
            for site in structure:
                for species, occu in site.species.items():
                    elem_map[species.element] += occu
            comp = Composition(elem_map)
        else:
            raise TypeError("Invalid type for element.")

        if Element("O") not in comp or comp.is_element:
            return "None", 0

        for site in structure:
            syms = [sp.symbol for sp in site.species]
            if "O" in syms:
                o_sites_frac_coords.append(site.frac_coords)
            if "H" in syms:
                h_sites_frac_coords.append(site.frac_coords)

        if h_sites_frac_coords:
            dist_matrix = lattice.get_all_distances(o_sites_frac_coords, h_sites_frac_coords)
            if np.any(dist_matrix < relative_cutoff * 0.93):
                return "hydroxide", int(len(np.where(dist_matrix < relative_cutoff * 0.93)[0]) / 2)
        dist_matrix = lattice.get_all_distances(o_sites_frac_coords, o_sites_frac_coords)
        np.fill_diagonal(dist_matrix, 1000)
        is_superoxide = is_peroxide = is_ozonide = False
        bond_atoms = []
        if np.any(dist_matrix < relative_cutoff * 1.35):
            bond_atoms = np.where(dist_matrix < relative_cutoff * 1.35)[0]
            is_superoxide = True
        elif np.any(dist_matrix < relative_cutoff * 1.49):
            is_peroxide = True
            bond_atoms = np.where(dist_matrix < relative_cutoff * 1.49)[0]
        if is_superoxide and len(bond_atoms) > len(set(bond_atoms)):
            is_superoxide = False
            is_ozonide = True

        n_bonds = len(set(bond_atoms))

        if is_ozonide:
            str_oxide = "ozonide"
        elif is_superoxide:
            str_oxide = "superoxide"
        elif is_peroxide:
            str_oxide = "peroxide"
        else:
            str_oxide = "oxide"
        if str_oxide == "oxide":
            n_bonds = int(comp["O"])
        return str_oxide, n_bonds


def oxide_type(
    structure: Structure, relative_cutoff: float = 1.1, return_nbonds: bool = False
) -> str | tuple[str, int]:
    """
    Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.

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
    Determines if a structure is a sulfide/polysulfide/sulfate.

    Args:
        structure (Structure): Input structure.

    Returns:
        str: sulfide/polysulfide or None if structure is a sulfate.
    """
    structure = structure.copy().remove_oxidation_states()
    sulphur = Element("S")
    comp = structure.composition
    if comp.is_element or sulphur not in comp:
        return None

    try:
        finder = SpacegroupAnalyzer(structure, symprec=0.1)
        symm_structure = finder.get_symmetrized_structure()
        s_sites = [sites[0] for sites in symm_structure.equivalent_sites if sites[0].specie == sulphur]
    except Exception:
        # Sometimes the symmetry analyzer fails for some tolerance or other issues. This is a fall back that simply
        # analyzes all S sites.
        s_sites = [site for site in structure if site.specie == sulphur]

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
        avg_electroneg = np.mean([elem.X for elem in coord_elements])
        if avg_electroneg > sulphur.X:
            return "sulfate"
        if avg_electroneg == sulphur.X and sulphur in coord_elements:
            return "polysulfide"
        return "sulfide"

    types = {process_site(site) for site in s_sites}
    if "sulfate" in types:
        return None
    if "polysulfide" in types:
        return "polysulfide"
    return "sulfide"
