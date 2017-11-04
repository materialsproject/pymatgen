# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import warnings
import six
import ruamel.yaml as yaml
import os

"""
This module provides classes to perform topological analyses of structures.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import math
from math import pi, sin, asin, atan, sqrt, exp, cos, acos
import numpy as np
import itertools
import collections

from warnings import warn
from scipy.spatial import Voronoi
from pymatgen import PeriodicSite
from pymatgen import Element, Specie, Composition
from pymatgen.util.num import abs_cap
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.surface import Slab, SlabGenerator


default_op_paras = {}
with open(os.path.join(os.path.dirname(
        __file__),'op_paras.yaml'), "rt") as f:
    default_op_paras = yaml.safe_load(f)
    f.close()

class VoronoiCoordFinder(object):
    """
    Uses a Voronoi algorithm to determine the coordination for each site in a
    structure.

    Args:
        structure (Structure): Input structure
        target ([Element/Specie]): A list of target species to determine
            coordination for.
        cutoff (float): Radius in Angstrom cutoff to look for coordinating
            atoms. Defaults to 10.0.
        allow_pathological (bool): whether to allow infinite vertices in
            determination of Voronoi coordination
    """

    def __init__(self, structure, target=None, cutoff=10.0,
                 allow_pathological=False):
        warnings.warn("VoronoiCoordFinder will be moved to local_env.py"
                      " with pymatgen >= 2018")

        self._structure = structure
        self.cutoff = cutoff
        self.allow_pathological = allow_pathological
        if target is None:
            self._target = structure.composition.elements
        else:
            self._target = target

    def get_voronoi_polyhedra(self, n):
        """
        Gives a weighted polyhedra around a site. This uses the voronoi
        construction with solid angle weights.
        See ref: A Proposed Rigorous Definition of Coordination Number,
        M. O'Keeffe, Acta Cryst. (1979). A35, 772-775

        Args:
            n (int): Site index

        Returns:
            A dict of sites sharing a common Voronoi facet with the site
            n and their solid angle weights
        """
        localtarget = self._target
        center = self._structure[n]
        neighbors = self._structure.get_sites_in_sphere(
            center.coords, self.cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        qvoronoi_input = [s.coords for s in neighbors]
        voro = Voronoi(qvoronoi_input)
        all_vertices = voro.vertices

        results = {}
        for nn, vind in voro.ridge_dict.items():
            if 0 in nn:
                if -1 in vind:
                    if self.allow_pathological:
                        continue
                    else:
                        raise RuntimeError("This structure is pathological,"
                                           " infinite vertex in the voronoi "
                                           "construction")

                facets = [all_vertices[i] for i in vind]
                results[neighbors[sorted(nn)[1]]] = solid_angle(
                    center.coords, facets)

        maxangle = max(results.values())

        resultweighted = {}
        for nn, angle in results.items():
            # is nn site is ordered use "nn.specie" to get species, else use "nn.species_and_occu" to get species
            if nn.is_ordered:
                if nn.specie in localtarget:
                    resultweighted[nn] = angle / maxangle
            else:  # is nn site is disordered
                for disordered_sp in nn.species_and_occu.keys():
                    if disordered_sp in localtarget:
                        resultweighted[nn] = angle / maxangle

        return resultweighted

    def get_coordination_number(self, n):
        """
        Returns the coordination number of site with index n.

        Args:
            n (int): Site index
        """
        return sum(self.get_voronoi_polyhedra(n).values())

    def get_coordinated_sites(self, n, tol=0, target=None):
        """
        Returns the sites that are in the coordination radius of site with
        index n.

        Args:
            n (int): Site index.
            tol (float): Weight tolerance to determine if a particular pair is
                considered a neighbor.
            target (Element): Target element

        Returns:
            Sites coordinating input site.
        """
        coordinated_sites = []
        for site, weight in self.get_voronoi_polyhedra(n).items():
            if weight > tol and (target is None or site.specie == target):
                coordinated_sites.append(site)
        return coordinated_sites


class JMolCoordFinder:
    """
    Determine coordinated sites and coordination number using an emulation of 
    JMol's default autoBond() algorithm. This version of the algorithm does not 
    take into account any information regarding known charge states.
    """

    def __init__(self, el_radius_updates=None):
        """
        Initialize coordination finder parameters (atomic radii)
        
        Args:
            el_radius_updates: (dict) symbol->float to override default atomic 
                radii table values 
        """

        warnings.warn("JMolCoordFinder will be moved to local_env.py"
                " with pymatgen >= 2018")

        # load elemental radii table
        bonds_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "bonds_jmol_ob.yaml")
        with open(bonds_file, 'r') as f:
            self.el_radius = yaml.safe_load(f)

        # update any user preference elemental radii
        if el_radius_updates:
            self.el_radius.update(el_radius_updates)

    def get_max_bond_distance(self, el1_sym, el2_sym, constant=0.56):
        """
        Use JMol algorithm to determine bond length from atomic parameters
        Args:
            el1_sym: (str) symbol of atom 1
            el2_sym: (str) symbol of atom 2
            constant: (float) factor to tune model

        Returns: (float) max bond length

        """
        return sqrt(
            (self.el_radius[el1_sym] + self.el_radius[el2_sym] + constant) ** 2)

    def get_coordination_number(self, structure, n, tol=1E-3):
        """
        Get the coordination number of a site
        Args:
            structure: (Structure)
            n: (int) index of site in the structure to get CN for
            tol: (float) a numerical tolerance to extend search

        Returns: (int) the coordination number
        """

        return len(self.get_coordinated_sites(structure, n, tol))

    def get_coordinated_sites(self, structure, n, tol=1E-3):
        """
        Get the coordinated sites for a site
        Args:
            structure: (Structure)
            n: (int) index of site in the structure to analyze
            tol: (float) a numerical tolerance to extend search

        Returns: ([sites]) a list of coordinated sites
        """
        site = structure[n]

        # determine relevant bond lengths based on atomic radii table
        bonds = {}
        for el in structure.composition.elements:
            bonds[site.specie, el] = self.get_max_bond_distance(
                site.specie.symbol, el.symbol)

        # search for neighbors up to max bond length + tolerance
        max_rad = max(bonds.values()) + tol

        all_neighbors = []
        for neighb, dist in structure.get_neighbors(site, max_rad):
            # confirm neighbor based on bond length specific to atom pair
            if dist <= bonds[(site.specie, neighb.specie)] + tol:
                all_neighbors.append(neighb)

        return all_neighbors


def average_coordination_number(structures, freq=10):
    """
    Calculates the ensemble averaged Voronoi coordination numbers
    of a list of Structures using VoronoiCoordFinder.
    Typically used for analyzing the output of a Molecular Dynamics run.
    Args:
        structures (list): list of Structures.
        freq (int): sampling frequency of coordination number [every freq steps].
    Returns:
        Dictionary of elements as keys and average coordination numbers as values.
    """
    coordination_numbers = {}
    for el in structures[0].composition.elements:
        coordination_numbers[el.name] = 0.0
    count = 0
    for t in range(len(structures)):
        if t % freq != 0:
            continue
        count += 1
        vor = VoronoiCoordFinder(structures[t])
        for atom in range(len(structures[0])):
            cn = vor.get_coordination_number(atom)
            coordination_numbers[structures[t][atom].species_string] += cn
    elements = structures[0].composition.as_dict()
    for el in coordination_numbers:
        coordination_numbers[el] = coordination_numbers[el] / elements[
            el] / count
    return coordination_numbers


class VoronoiAnalyzer(object):
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

    Args:
        cutoff (float): cutoff distance to search for neighbors of a given atom
            (default = 5.0)
        qhull_options (str): options to pass to qhull (optional)
    """

    def __init__(self, cutoff=5.0, qhull_options="Qbb Qc Qz"):
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
                else:
                    try:
                        vor_index[len(voro.ridge_dict[key]) - 3] += 1
                    except IndexError:
                        # If a facet has more than 10 edges, it's skipped here.
                        pass
        return vor_index

    def analyze_structures(self, structures, step_freq=10,
                           most_frequent_polyhedra=15):
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
        return sorted(voro_dict.items(),
                      key=lambda x: (x[1], x[0]),
                      reverse=True)[:most_frequent_polyhedra]

    @staticmethod
    def plot_vor_analysis(voronoi_ensemble):
        t = zip(*voronoi_ensemble)
        labels = t[0]
        val = list(t[1])
        tot = np.sum(val)
        val = [float(j) / tot for j in val]
        pos = np.arange(len(val)) + .5  # the bar centers on the y axis
        import matplotlib.pyplot as plt
        plt.figure()
        plt.barh(pos, val, align='center', alpha=0.5)
        plt.yticks(pos, labels)
        plt.xlabel('Count')
        plt.title('Voronoi Spectra')
        plt.grid(True)
        return plt


class RelaxationAnalyzer(object):
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
            raise ValueError("Initial and final structures have different " +
                             "formulas!")
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
        d = {l: getattr(final_latt, l) / getattr(initial_latt, l) - 1
             for l in ["a", "b", "c"]}
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


class VoronoiConnectivity(object):
    """
    Computes the solid angles swept out by the shared face of the voronoi
    polyhedron between two sites.

    Args:
        structure (Structure): Input structure
        cutoff (float) Cutoff distance.
    """

    # Radius in Angstrom cutoff to look for coordinating atoms

    def __init__(self, structure, cutoff=10):
        self.cutoff = cutoff
        self.s = structure
        recp_len = np.array(self.s.lattice.reciprocal_lattice.abc)
        i = np.ceil(cutoff * recp_len / (2 * pi))
        offsets = np.mgrid[-i[0]:i[0] + 1, -i[1]:i[1] + 1, -i[2]:i[2] + 1].T
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
                warn('Found connectivity with infinite vertex. '
                     'Cutoff is too low, and results may be '
                     'incorrect')
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
        atoms_n_occu = self.s[site_index].species_and_occu
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
        v = -np.dot(n[i], n[i + 1]) \
            / (np.linalg.norm(n[i]) * np.linalg.norm(n[i + 1]))
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
    jmc = JMolCoordFinder(el_radius_updates)

    bonds_lens = {}
    els = sorted(structure.composition.elements, key=lambda x: x.Z)

    for i1 in range(len(els)):
        for i2 in range(len(els) - i1):
            bonds_lens[els[i1], els[i1 + i2]] = jmc.get_max_bond_distance(
                els[i1].symbol, els[i1 + i2].symbol)

    return bonds_lens


def get_dimensionality(structure, max_hkl=2, el_radius_updates=None,
                       min_slab_size=5, min_vacuum_size=5,
                       standardize=True, bonds=None):

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
        structure = SpacegroupAnalyzer(structure).\
            get_conventional_standard_structure()

    if not bonds:
        bonds = get_max_bond_lengths(structure, el_radius_updates)

    num_surfaces = 0
    for h in range(max_hkl):
        for k in range(max_hkl):
            for l in range(max_hkl):
                if max([h, k, l]) > 0 and num_surfaces < 2:
                    sg = SlabGenerator(structure, (h, k, l),
                                       min_slab_size=min_slab_size,
                                       min_vacuum_size=min_vacuum_size)
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
    ox_type = oxide_type(structure, relative_cutoff)
    if ox_type == "peroxide":
        return True
    else:
        return False


class OxideType(object):
    """
    Separate class for determining oxide type.

    Args:
        structure: Input structure.
        relative_cutoff: Relative_cutoff * act. cutoff stipulates the max.
            distance two O atoms must be from each other. Default value is
            1.1. At most 1.1 is recommended, nothing larger, otherwise the
            script cannot distinguish between superoxides and peroxides.
    """

    def __init__(self, structure, relative_cutoff=1.1):
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
        elif isinstance(structure.composition.elements[0], Specie):
            elmap = collections.defaultdict(float)
            for site in structure:
                for species, occu in site.species_and_occu.items():
                    elmap[species.element] += occu
            comp = Composition(elmap)
        if Element("O") not in comp or comp.is_element:
            return "None", 0

        for site in structure:
            syms = [sp.symbol for sp in site.species_and_occu.keys()]
            if "O" in syms:
                o_sites_frac_coords.append(site.frac_coords)
            if "H" in syms:
                h_sites_frac_coords.append(site.frac_coords)

        if h_sites_frac_coords:
            dist_matrix = lattice.get_all_distances(o_sites_frac_coords,
                                                    h_sites_frac_coords)
            if np.any(dist_matrix < relative_cutoff * 0.93):
                return "hydroxide", len(
                    np.where(dist_matrix < relative_cutoff * 0.93)[0]) / 2.0
        dist_matrix = lattice.get_all_distances(o_sites_frac_coords,
                                                o_sites_frac_coords)
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
    else:
        return ox_obj.oxide_type


def sulfide_type(structure):
    """
    Determines if a structure is a sulfide/polysulfide

    Args:
        structure (Structure): Input structure.

    Returns:
        (str) sulfide/polysulfide/sulfate
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
        neighbors = structure.get_neighbors(site, 4)
        neighbors = sorted(neighbors, key=lambda n: n[1])
        nn, dist = neighbors[0]
        coord_elements = [site.specie for site, d in neighbors
                          if d < dist + 0.4][:4]
        avg_electroneg = np.mean([e.X for e in coord_elements])
        if avg_electroneg > s.X:
            return "sulfate"
        elif avg_electroneg == s.X and s in coord_elements:
            return "polysulfide"
        else:
            return "sulfide"

    types = set([process_site(site) for site in s_sites])
    if "sulfate" in types:
        return None
    elif "polysulfide" in types:
        return "polysulfide"
    else:
        return "sulfide"


def gramschmidt(vin, uin):
    """
    Returns that part of the first input vector
    that is orthogonal to the second input vector.
    The output vector is not normalized.

    Args:
        vin (numpy array):
            first input vector
        uin (numpy array):
            second input vector
    """

    vin_uin = np.inner(vin, uin)
    uin_uin = np.inner(uin, uin)
    if uin_uin <= 0.0:
        raise ValueError("Zero or negative inner product!")
    return vin - (vin_uin / uin_uin) * uin


class OrderParameters(object):
    """
    This class permits the calculation of various types of local
    structure order parameters.
    """

    __supported_types = (
            "cn", "sgl_bd", "bent", "tri_plan", "reg_tri", "sq_plan", \
            "pent_plan", "sq", "tet", "tri_pyr", \
            "sq_pyr", "sq_pyr_legacy", "tri_bipyr", "sq_bipyr", "oct", \
            "oct_legacy", "pent_pyr", "hex_pyr", "pent_bipyr", "hex_bipyr", \
            "T", "cuboct", "see_saw", "bcc", "q2", "q4", "q6")

    def __init__(self, types, parameters=None, cutoff=-10.0):
        """
        Args:
            types ([string]): list of strings representing the types of
                order parameters to be calculated. Note that multiple
                mentions of the same type may occur. Currently available
                types recognize following environments:
                  "cn": simple coordination number---normalized
                        if desired;
                  "sgl_bd": single bonds;
                  "bent": bent (angular) coordinations
                          (Zimmermann & Jain, in progress, 2017);
                  "T": T-shape coordinations;
                  "see_saw": see saw-like coordinations;
                  "tet": tetrahedra
                         (Zimmermann et al., submitted, 2017);
                  "oct": octahedra
                         (Zimmermann et al., submitted, 2017);
                  "bcc": body-centered cubic environments (Peters,
                         J. Chem. Phys., 131, 244103, 2009);
                  "tri_plan": trigonal planar environments;
                  "sq_plan": square planar environments;
                  "pent_plan": pentagonal planar environments;
                  "tri_pyr": trigonal pyramids (coordinated atom is in
                             the center of the basal plane);
                  "sq_pyr": square pyramids;
                  "pent_pyr": pentagonal pyramids;
                  "hex_pyr": hexagonal pyramids;
                  "tri_bipyr": trigonal bipyramids;
                  "sq_bipyr": square bipyramids;
                  "pent_bipyr": pentagonal bipyramids;
                  "hex_bipyr": hexagonal bipyramids;
                  "cuboct": cuboctahedra;
                  "q2": motif-unspecific bond orientational order
                        parameter (BOOP) of weight l=2 (Steinhardt
                        et al., Phys. Rev. B, 28, 784-805, 1983);
                  "q4": BOOP of weight l=4;
                  "q6": BOOP of weight l=6.
                  "reg_tri": regular triangle with varying height
                             to basal plane;
                  "sq": square coordination (cf., "reg_tri");
                  "oct_legacy": original Peters-style OP recognizing
                                octahedral coordination environments
                                (Zimmermann et al., J. Am. Chem. Soc.,
                                137, 13352-13361, 2015) that can, however,
                                produce small negative values sometimes.
                  "sq_pyr_legacy": square pyramids (legacy);
            parameters ([dict]): list of dictionaries
                that store float-type parameters associated with the
                definitions of the different order parameters
                (length of list = number of OPs). If an entry
                is None, default values are used that are read from
                the op_paras.yaml file. With few exceptions, 9 different
                parameters are used across all OPs:
                  "norm": normalizing constant (used in "cn"
                      (default value: 1)).
                  "TA": target angle (TA) in fraction of 180 degrees
                      ("bent" (1), "tet" (0.6081734479693927),
                      "tri_plan" (0.66666666667), "pent_plan" (0.6),
                      "sq_pyr_legacy" (0.5)).
                  "IGW_TA": inverse Gaussian width (IGW) for penalizing
                      angles away from the target angle in inverse
                      fractions of 180 degrees to ("bent" and "tet" (15),
                      "tri_plan" (13.5), "pent_plan" (18),
                      "sq_pyr_legacy" (30)).
                  "IGW_EP": IGW for penalizing angles away from the
                      equatorial plane (EP) at 90 degrees ("T", "see_saw",
                      "oct", "sq_plan", "tri_pyr", "sq_pyr", "pent_pyr",
                      "hex_pyr", "tri_bipyr", "sq_bipyr", "pent_bipyr",
                      "hex_bipyr", and "oct_legacy" (18)).
                  "fac_AA": factor applied to azimuth angle (AA) in cosine
                      term ("T", "tri_plan", and "sq_plan" (1), "tet",
                      "tri_pyr", and "tri_bipyr" (1.5), "oct", "sq_pyr",
                      "sq_bipyr", and "oct_legacy" (2), "pent_pyr"
                      and "pent_bipyr" (2.5), "hex_pyr" and
                      "hex_bipyr" (3)).
                  "exp_cos_AA": exponent applied to cosine term of AA
                      ("T", "tet", "oct", "tri_plan", "sq_plan",
                      "tri_pyr", "sq_pyr", "pent_pyr", "hex_pyr",
                      "tri_bipyr", "sq_bipyr", "pent_bipyr", "hex_bipyr",
                      and "oct_legacy" (2)).
                  "min_SPP": smallest angle (in radians) to consider
                      a neighbor to be
                      at South pole position ("see_saw", "oct", "bcc",
                      "sq_plan", "tri_bipyr", "sq_bipyr", "pent_bipyr",
                      "hex_bipyr", "cuboct", and "oct_legacy"
                      (2.792526803190927)).
                  "IGW_SPP": IGW for penalizing angles away from South
                      pole position ("see_saw", "oct", "bcc", "sq_plan",
                      "tri_bipyr", "sq_bipyr", "pent_bipyr", "hex_bipyr",
                      "cuboct", and "oct_legacy" (15)).
                  "w_SPP": weight for South pole position relative to
                      equatorial positions ("see_saw" and "sq_plan" (1),
                      "cuboct" (1.8), "tri_bipyr" (2), "oct",
                      "sq_bipyr", and "oct_legacy" (3), "pent_bipyr" (4),
                      "hex_bipyr" (5), "bcc" (6)).
            cutoff (float): Cutoff radius to determine which nearest
                neighbors are supposed to contribute to the order
                parameters. If the value is negative the neighboring
                sites found by distance and cutoff radius are further
                pruned using the get_coordinated_sites method from the
                VoronoiCoordFinder class.
        """
        for t in types:
            if t not in OrderParameters.__supported_types:
                raise ValueError("Unknown order parameter type (" + \
                                 t + ")!")
        self._types = tuple(types)

        self._paras = []
        for i, t in enumerate(self._types):
            d = default_op_paras[t].copy() if default_op_paras[t] is not None \
                else None
            if parameters is None:
                self._paras.append(d)
            elif parameters[i] is None:
                self._paras.append(d)
            else:
                self._paras.append(parameters[i].copy())

        self._computerijs = self._computerjks = self._geomops = False
        self._geomops2 = self._boops = False
        self._max_trig_order = -1

        # Add here any additional flags to be used during calculation.
        if "sgl_bd" in self._types:
            self._computerijs = True
        if not set(self._types).isdisjoint(
            ["tet", "oct", "bcc", "sq_pyr", "sq_pyr_legacy", \
             "tri_bipyr", "sq_bipyr", "oct_legacy", "tri_plan", \
             "sq_plan", "pent_plan", "tri_pyr", "pent_pyr", "hex_pyr", \
             "pent_bipyr", "hex_bipyr", "T", "cuboct"]):
            self._computerijs = self._geomops = True
        if not set(self._types).isdisjoint(["reg_tri", "sq"]):
            self._computerijs = self._computerjks = self._geomops2 = True
        if not set(self._types).isdisjoint(["q2", "q4", "q6"]):
            self._computerijs = self._boops = True
        if "q2" in self._types:
            self._max_trig_order = 2
        if "q4" in self._types:
            self._max_trig_order = 4
        if "q6" in self._types:
            self._max_trig_order = 6

        # Finish parameter treatment.
        if cutoff < 0.0:
            self._cutoff = -cutoff
            self._voroneigh = True
        elif cutoff > 0.0:
            self._cutoff = cutoff
            self._voroneigh = False
        else:
            raise ValueError("Cutoff radius is zero!")

        # Further variable definitions.
        self._last_nneigh = -1
        self._pow_sin_t = {}
        self._pow_cos_t = {}
        self._sin_n_p = {}
        self._cos_n_p = {}

    @property
    def num_ops(self):

        """"
        Returns:
            int: the number of different order parameters that are targeted
                to be calculated.
        """

        return len(self._types)

    @property
    def last_nneigh(self):

        """"
        Returns:
            int: the number of neighbors encountered during the most
                recent order parameter calculation. A value of -1 indicates
                that no such calculation has yet been performed for this
                instance.
        """

        return len(self._last_nneigh)

    def compute_trigonometric_terms(self, thetas, phis):

        """"
        Computes trigonometric terms that are required to
        calculate bond orientational order parameters using
        internal variables.

        Args:
            thetas ([float]): polar angles of all neighbors in radians.
            phis ([float]): azimuth angles of all neighbors in radians.
                The list of
                azimuth angles of all neighbors in radians.  The list of
                azimuth angles is expected to have the same size as the
                list of polar angles; otherwise, a ValueError is raised.
                Also, the two lists of angles have to be coherent in
                order. That is, it is expected that the order in the list
                of azimuth angles corresponds to a distinct sequence of
                neighbors. And, this sequence has to equal the sequence
                of neighbors in the list of polar angles.

        """

        if len(thetas) != len(phis):
            raise ValueError("List of polar and azimuthal angles have to be"
                             " equal!")

        self._pow_sin_t.clear()
        self._pow_cos_t.clear()
        self._sin_n_p.clear()
        self._cos_n_p.clear()

        self._pow_sin_t[1] = [sin(float(t)) for t in thetas]
        self._pow_cos_t[1] = [cos(float(t)) for t in thetas]
        self._sin_n_p[1] = [sin(float(p)) for p in phis]
        self._cos_n_p[1] = [cos(float(p)) for p in phis]

        for i in range(2, self._max_trig_order + 1):
            self._pow_sin_t[i] = [e[0] * e[1] for e in zip(
                self._pow_sin_t[i - 1], self._pow_sin_t[1])]
            self._pow_cos_t[i] = [e[0] * e[1] for e in zip(
                self._pow_cos_t[i - 1], self._pow_cos_t[1])]
            self._sin_n_p[i] = [sin(float(i) * float(p)) \
                                for p in phis]
            self._cos_n_p[i] = [cos(float(i) * float(p)) \
                                for p in phis]

    def get_q2(self, thetas=None, phis=None):

        """
        Calculates the value of the bond orientational order parameter of
        weight l=2.  If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh.  Otherwise, it is expected that the
        compute_trigonometric_terms function has been just called.

        Args:
            thetas ([float]): polar angles of all neighbors in radians.
            phis ([float]): azimuth angles of all neighbors in radians.

        Returns:
            float: bond orientational order parameter of weight l=2
                corresponding to the input angles thetas and phis.
        """

        if thetas is not None and phis is not None:
            self.compute_trigonometric_terms(thetas, phis)
        nnn = len(self._pow_sin_t[1])
        nnn_range = range(nnn)

        sqrt_15_2pi = sqrt(15.0 / (2.0 * pi))
        sqrt_5_pi = sqrt(5.0 / pi)

        pre_y_2_2 = [0.25 * sqrt_15_2pi * val for val in self._pow_sin_t[2]]
        pre_y_2_1 = [0.5 * sqrt_15_2pi * val[0] * val[1]
                     for val in zip(self._pow_sin_t[1], self._pow_cos_t[1])]

        acc = 0.0

        # Y_2_-2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_2_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_2_-1
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_2_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_2_0
        real = imag = 0.0
        for i in nnn_range:
            real += 0.25 * sqrt_5_pi * (3.0 * self._pow_cos_t[2][i] - 1.0)
        acc += (real * real)

        # Y_2_1
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_2_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_2_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_2_2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_2[i] * self._cos_n_p[2][i]
            imag += pre_y_2_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        q2 = sqrt(4.0 * pi * acc / (5.0 * float(nnn * nnn)))
        return q2

    def get_q4(self, thetas=None, phis=None):

        """
        Calculates the value of the bond orientational order parameter of
        weight l=4.  If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh.  Otherwise, it is expected that the
        compute_trigonometric_terms function has been just called.

        Args:
            thetas ([float]): polar angles of all neighbors in radians.
            phis ([float]): azimuth angles of all neighbors in radians.

        Returns:
            float: bond orientational order parameter of weight l=4
                corresponding to the input angles thetas and phis.
        """

        if thetas is not None and phis is not None:
            self.compute_trigonometric_terms(thetas, phis)
        nnn = len(self._pow_sin_t[1])
        nnn_range = range(nnn)

        i16_3 = 3.0 / 16.0
        i8_3 = 3.0 / 8.0

        sqrt_35_pi = sqrt(35.0 / pi)
        sqrt_35_2pi = sqrt(35.0 / (2.0 * pi))
        sqrt_5_pi = sqrt(5.0 / pi)
        sqrt_5_2pi = sqrt(5.0 / (2.0 * pi))
        sqrt_1_pi = sqrt(1.0 / pi)

        pre_y_4_4 = [i16_3 * sqrt_35_2pi * val for val in self._pow_sin_t[4]]
        pre_y_4_3 = [i8_3 * sqrt_35_pi * val[0] * val[1] \
                     for val in zip(self._pow_sin_t[3], self._pow_cos_t[1])]
        pre_y_4_2 = [i8_3 * sqrt_5_2pi * val[0] * (7.0 * val[1] - 1.0) \
                     for val in zip(self._pow_sin_t[2], self._pow_cos_t[2])]
        pre_y_4_1 = [i8_3 * sqrt_5_pi * val[0] * (7.0 * val[1] - 3.0 * val[2]) \
                     for val in zip(self._pow_sin_t[1], self._pow_cos_t[3], \
                                    self._pow_cos_t[1])]

        acc = 0.0

        # Y_4_-4
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_4_4[i] * self._sin_n_p[4][i]
        acc += (real * real + imag * imag)

        # Y_4_-3
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_4_3[i] * self._sin_n_p[3][i]
        acc += (real * real + imag * imag)

        # Y_4_-2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_4_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_4_-1
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_4_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_4_0
        real = imag = 0.0
        for i in nnn_range:
            real += i16_3 * sqrt_1_pi * (35.0 * self._pow_cos_t[4][i] - \
                                         30.0 * self._pow_cos_t[2][i] + 3.0)
        acc += (real * real)

        # Y_4_1
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_4_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_4_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_4_2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_2[i] * self._cos_n_p[2][i]
            imag += pre_y_4_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_4_3
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_4_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_4_3[i] * self._sin_n_p[3][i]
        acc += (real * real + imag * imag)

        # Y_4_4
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_4[i] * self._cos_n_p[4][i]
            imag += pre_y_4_4[i] * self._sin_n_p[4][i]
        acc += (real * real + imag * imag)

        q4 = sqrt(4.0 * pi * acc / (9.0 * float(nnn * nnn)))
        return q4

    def get_q6(self, thetas=None, phis=None):

        """
        Calculates the value of the bond orientational order parameter of
        weight l=6.  If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh.  Otherwise, it is expected that the
        compute_trigonometric_terms function has been just called.

        Args:
            thetas ([float]): polar angles of all neighbors in radians.
            phis ([float]): azimuth angles of all neighbors in radians.

        Returns:
            float: bond orientational order parameter of weight l=6
                corresponding to the input angles thetas and phis.
        """

        if thetas is not None and phis is not None:
            self.compute_trigonometric_terms(thetas, phis)
        nnn = len(self._pow_sin_t[1])
        nnn_range = range(nnn)

        i64 = 1.0 / 64.0
        i32 = 1.0 / 32.0
        i32_3 = 3.0 / 32.0
        i16 = 1.0 / 16.0

        sqrt_3003_pi = sqrt(3003.0 / pi)
        sqrt_1001_pi = sqrt(1001.0 / pi)
        sqrt_91_2pi = sqrt(91.0 / (2.0 * pi))
        sqrt_1365_pi = sqrt(1365.0 / pi)
        sqrt_273_2pi = sqrt(273.0 / (2.0 * pi))
        sqrt_13_pi = sqrt(13.0 / pi)

        pre_y_6_6 = [i64 * sqrt_3003_pi * val for val in self._pow_sin_t[6]]
        pre_y_6_5 = [i32_3 * sqrt_1001_pi * val[0] * val[1] \
                     for val in zip(self._pow_sin_t[5], self._pow_cos_t[1])]
        pre_y_6_4 = [i32_3 * sqrt_91_2pi * val[0] * (11.0 * val[1] - 1.0) \
                     for val in zip(self._pow_sin_t[4], self._pow_cos_t[2])]
        pre_y_6_3 = [
            i32 * sqrt_1365_pi * val[0] * (11.0 * val[1] - 3.0 * val[2]) \
            for val in zip(self._pow_sin_t[3], self._pow_cos_t[3], \
                           self._pow_cos_t[1])]
        pre_y_6_2 = [i64 * sqrt_1365_pi * val[0] * (33.0 * val[1] - \
                                                    18.0 * val[2] + 1.0) for val
                     in zip(self._pow_sin_t[2], \
                            self._pow_cos_t[4], self._pow_cos_t[2])]
        pre_y_6_1 = [i16 * sqrt_273_2pi * val[0] * (33.0 * val[1] - \
                                                    30.0 * val[2] + 5.0 * val[
                                                        3]) for val in
                     zip(self._pow_sin_t[1], \
                         self._pow_cos_t[5], self._pow_cos_t[3],
                         self._pow_cos_t[1])]

        acc = 0.0

        # Y_6_-6
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i]  # cos(x) =  cos(-x)
            imag -= pre_y_6_6[i] * self._sin_n_p[6][i]  # sin(x) = -sin(-x)
        acc += (real * real + imag * imag)

        # Y_6_-5
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += (real * real + imag * imag)

        # Y_6_-4
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += (real * real + imag * imag)

        # Y_6_-3
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += (real * real + imag * imag)

        # Y_6_-2
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_6_-1
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_6_0
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += i32 * sqrt_13_pi * (231.0 * self._pow_cos_t[6][i] - \
                                        315.0 * self._pow_cos_t[4][i] + 105.0 *
                                        self._pow_cos_t[2][i] - 5.0)
        acc += (real * real)

        # Y_6_1
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_6_2
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag += pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_6_3
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += (real * real + imag * imag)

        # Y_6_4
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag += pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += (real * real + imag * imag)

        # Y_6_5
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += (real * real + imag * imag)

        # Y_6_6
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i]
            imag += pre_y_6_6[i] * self._sin_n_p[6][i]
        acc += (real * real + imag * imag)

        q6 = sqrt(4.0 * pi * acc / (13.0 * float(nnn * nnn)))
        return q6

    def get_type(self, index):

        """
        Return type of order parameter at the index provided and
        represented by a short string.

        Args:
            index (int): index of order parameter for which type is
                to be returned.
        Returns:
            str: OP type.
        """
        if index < 0 or index >= len(self._types):
            raise ValueError("Index for getting order parameter type"
                             " out-of-bounds!")
        return self._types[index]

    def get_parameters(self, index):

        """
        Returns list of floats that represents
        the parameters associated with calculation of the order
        parameter that was defined at the index provided.
        Attention: the parameters do not need to equal those originally
        inputted because of processing out of efficiency reasons.

        Args:
            index (int):
                index of order parameter for which associated parameters
                are to be returned.
        Returns:
            [float]: parameters of a given OP.
        """
        if index < 0 or index >= len(self._types):
            raise ValueError("Index for getting parameters associated with"
                             " order parameter calculation out-of-bounds!")
        return self._paras[index]

    def get_order_parameters(self, structure, n, indices_neighs=None, \
                             tol=0.0, target_spec=None):

        """
        Compute all order parameters of site n.

        Args:
            structure (Structure): input structure.
            n (int): index of site in input structure,
                for which OPs are to be
                calculated.  Note that we do not use the sites iterator
                here, but directly access sites via struct[index].
            indices_neighs ([int]): list of indices of those neighbors
                in Structure object
                structure that are to be considered for OP computation.
                This optional argument overwrites the way neighbors are
                to be determined as defined in the constructor (i.e.,
                Voronoi coordination finder via negative cutoff radius
                vs constant cutoff radius if cutoff was positive).
                We do not use information about the underlying
                structure lattice if the neighbor indices are explicitly
                provided.  This has two important consequences.  First,
                the input Structure object can, in fact, be a
                simple list of Site objects.  Second, no nearest images
                of neighbors are determined when providing an index list.
                Note furthermore that this neighbor
                determination type ignores the optional target_spec
                argument.
            tol (float): threshold of weight
                (= solid angle / maximal solid angle)
                to determine if a particular pair is
                considered neighbors; this is relevant only in the case
                when Voronoi polyhedra are used to determine coordination
            target_spec (Specie): target species to be considered
                when calculating the order
                parameters of site n; None includes all species of input
                structure.

        Returns:
            [floats]: representing order parameters.  Should it not be
            possible to compute a given OP for a conceptual reason, the
            corresponding entry is None instead of a float.  For Steinhardt
            et al.'s bond orientational OPs and the other geometric OPs
            ("tet", "oct", "bcc", etc.),
            this can happen if there is a single
            neighbor around site n in the structure because that
            does not permit calculation of angles between multiple
            neighbors.
        """

        # Do error-checking and initialization.
        if n < 0:
            raise ValueError("Site index smaller zero!")
        if n >= len(structure):
            raise ValueError("Site index beyond maximum!")
        if indices_neighs is not None:
            for index in indices_neighs:
                if index >= len(structure):
                    raise ValueError("Neighbor site index beyond maximum!")
        if tol < 0.0:
            raise ValueError("Negative tolerance for weighted solid angle!")

        left_of_unity = 1.0 - 1.0e-12
        # The following threshold has to be adapted to non-Angstrom units.
        very_small = 1.0e-12

        # Find central site and its neighbors.
        # Note that we adopt the same way of accessing sites here as in
        # VoronoiCoordFinder; that is, not via the sites iterator.
        centsite = structure[n]
        if indices_neighs is not None:
            neighsites = [structure[index] for index in indices_neighs]
        elif self._voroneigh:
            vorocf = VoronoiCoordFinder(structure)
            neighsites = vorocf.get_coordinated_sites(n, tol, target_spec)
        else:
            # Structure.get_sites_in_sphere --> also other periodic images
            neighsitestmp = [i[0] for i in structure.get_sites_in_sphere(
                centsite.coords, self._cutoff)]
            neighsites = []
            if centsite not in neighsitestmp:
                raise ValueError("Could not find center site!")
            else:
                neighsitestmp.remove(centsite)
            if target_spec is None:
                neighsites = list(neighsitestmp)
            else:
                neighsites[:] = [site for site in neighsitestmp \
                                 if site.specie.symbol == target_spec]
        nneigh = len(neighsites)
        self._last_nneigh = nneigh

        # Prepare angle calculations, if applicable.
        rij = []
        rjk = []
        rijnorm = []
        rjknorm = []
        dist = []
        distjk_unique = []
        distjk = []
        centvec = centsite.coords
        if self._computerijs:
            for j, neigh in enumerate(neighsites):
                rij.append((neigh.coords - centvec))
                dist.append(np.linalg.norm(rij[j]))
                rijnorm.append((rij[j] / dist[j]))
        if self._computerjks:
            for j, neigh in enumerate(neighsites):
                rjk.append([])
                rjknorm.append([])
                distjk.append([])
                kk = 0
                for k in range(len(neighsites)):
                    if j != k:
                        rjk[j].append(neighsites[k].coords - neigh.coords)
                        distjk[j].append(np.linalg.norm(rjk[j][kk]))
                        if k > j:
                            distjk_unique.append(distjk[j][kk])
                        rjknorm[j].append(rjk[j][kk] / distjk[j][kk])
                        kk = kk + 1
        # Initialize OP list and, then, calculate OPs.
        ops = [0.0 for t in self._types]
        norms = [[0.0 for i in range(nneigh)] for t in self._types]

        # First, coordination number and distance-based OPs.
        for i, t in enumerate(self._types):
            if t == "cn":
                ops[i] = nneigh / self._paras[i]['norm']
            elif t == "sgl_bd":
                dist_sorted = sorted(dist)
                if len(dist_sorted) == 1:
                    ops[i] = 1.0
                elif len(dist_sorted) > 1:
                    ops[i] = 1.0 - dist_sorted[0]/dist_sorted[1]

        # Then, bond orientational OPs based on spherical harmonics
        # according to Steinhardt et al., Phys. Rev. B, 28, 784-805, 1983.
        if self._boops:
            thetas = []
            phis = []
            for j, vec in enumerate(rijnorm):

                # z is North pole --> theta between vec and (0, 0, 1)^T.
                # Because vec is normalized, dot product is simply vec[2].
                thetas.append(acos(max(-1.0, min(vec[2], 1.0))))
                tmpphi = 0.0

                # Compute phi only if it is not (almost) perfectly
                # aligned with z-axis.
                if vec[2] < left_of_unity and vec[2] > - (left_of_unity):
                    # x is prime meridian --> phi between projection of vec
                    # into x-y plane and (1, 0, 0)^T
                    tmpphi = acos(max(-1.0, min(vec[0] / (sqrt(
                        vec[0] * vec[0] + vec[1] * vec[1])), 1.0)))
                    if vec[1] < 0.0:
                        tmpphi = -tmpphi
                phis.append(tmpphi)

            # Note that None flags that we have too few neighbors
            # for calculating BOOPS.
            for i, t in enumerate(self._types):
                if t == "q2":
                    ops[i] = self.get_q2(thetas, phis) if len(
                        thetas) > 0 else None
                elif t == "q4":
                    ops[i] = self.get_q4(thetas, phis) if len(
                        thetas) > 0 else None
                elif t == "q6":
                    ops[i] = self.get_q6(thetas, phis) if len(
                        thetas) > 0 else None

        # Then, deal with the Peters-style OPs that are tailor-made
        # to recognize common structural motifs
        # (Peters, J. Chem. Phys., 131, 244103, 2009;
        #  Zimmermann et al., J. Am. Chem. Soc., under revision, 2015).
        if self._geomops:
            gaussthetak = [0.0 for t in self._types]  # not used by all OPs
            qsptheta = [[] for t in self._types]  # not used by all OPs
            ipi = 1.0 / pi
            piover2 = pi / 2.0
            tetangoverpi = acos(-1.0 / 3.0) * ipi # xxx: delete
            itetangminuspihalfoverpi = 1.0 / (tetangoverpi - 0.5)
            onethird = 1.0 / 3.0
            twothird = 2.0 / 3.0
            for j in range(nneigh):  # Neighbor j is put to the North pole.
                zaxis = rijnorm[j]
                for i, t in enumerate(self._types):
                    qsptheta[i].append(0.0)
                for k in range(nneigh):  # From neighbor k, we construct
                    if j != k:  # the prime meridian.
                        tmp = max(
                            -1.0, min(np.inner(zaxis, rijnorm[k]), 1.0))
                        thetak = acos(tmp)
                        xaxistmp = gramschmidt(rijnorm[k], zaxis)
                        if np.linalg.norm(xaxistmp) < very_small:
                            flag_xaxis = True
                        else:
                            xaxis = xaxistmp / np.linalg.norm(xaxistmp)
                            flag_xaxis = False

                        # Contributions of j-i-k angles, where i represents the central atom
                        # and j and k two of the neighbors.
                        for i, t in enumerate(self._types):
                            if t in ["bent", "sq_pyr_legacy"]:
                                tmp = self._paras[i]['IGW_TA'] * (
                                    thetak * ipi - self._paras[i]['TA'])
                                qsptheta[i][j] += exp(-0.5 * tmp * tmp)
                                norms[i][j] += 1
                            elif t in ["tri_plan", "tet"]:
                                tmp = self._paras[i]['IGW_TA'] * (
                                    thetak * ipi - self._paras[i]['TA'])
                                gaussthetak[i] = exp(-0.5 * tmp * tmp)
                            elif t in ["sq_plan", "oct", "oct_legacy", \
                                    "see_saw", "tri_bipyr", "sq_bipyr", \
                                    "pent_bipyr", "hex_bipyr", "cuboct"]:
                                if thetak >= self._paras[i]['min_SPP']:
                                    tmp = self._paras[i]['IGW_SPP'] * (
                                        thetak * ipi - 1.0)
                                    if t in ["tri_bipyr", "sq_bipyr", \
                                            "pent_bipyr", "hex_bipyr"]:
                                        qsptheta[i][j] = (self._paras[i]['w_SPP'] * \
                                            exp(-0.5 * tmp * tmp))
                                    else:
                                        qsptheta[i][j] += (self._paras[i]['w_SPP'] * \
                                            exp(-0.5 * tmp * tmp))
                                    norms[i][j] += self._paras[i]['w_SPP']
                            elif t == "pent_plan":
                                if thetak <= self._paras[i]['TA'] * pi:
                                    tmp = self._paras[i]['IGW_TA'] * (
                                        thetak * ipi - 0.4)
                                    gaussthetak[i] = exp(-0.5 * tmp * tmp)
                            elif t == "bcc" and j < k:
                                if thetak >= self._paras[i]['min_SPP']:
                                    tmp = self._paras[i]['IGW_SPP'] * (
                                        thetak * ipi - 1.0)
                                    qsptheta[i][j] += (self._paras[i]['w_SPP'] * \
                                        exp(-0.5 * tmp * tmp))
                                    norms[i][j] += self._paras[i]['w_SPP']

                        for m in range(nneigh):
                            if (m != j) and (m != k) and (not flag_xaxis):
                                tmp = max(
                                    -1.0, min(np.inner(zaxis, rijnorm[m]), 1.0))
                                thetam = acos(tmp)
                                xtwoaxistmp = gramschmidt(rijnorm[m], zaxis)
                                l = np.linalg.norm(xtwoaxistmp)
                                if l < very_small:
                                    flag_xtwoaxis = True
                                else:
                                    xtwoaxis = xtwoaxistmp / l
                                    phi = acos(max(
                                        -1.0,
                                        min(np.inner(xtwoaxis, xaxis), 1.0)))
                                    flag_xtwoaxis = False

                                # Contributions of j-i-m angle and
                                # angles between plane j-i-k and i-m vector.
                                if not flag_xaxis and not flag_xtwoaxis:
                                    for i, t in enumerate(self._types):
                                        if t in ["tri_plan", "tet"]:
                                            tmp = self._paras[i]['IGW_TA'] * (
                                                thetam * ipi - \
                                                self._paras[i]['TA'])
                                            tmp2 = cos(self._paras[i]['fac_AA'] * \
                                                phi) ** self._paras[i]['exp_cos_AA']
                                            ops[i] += gaussthetak[i] * exp(
                                                -0.5 * tmp * tmp) * tmp2
                                            norms[i][j] += 1
                                        elif t == "pent_plan":
                                            if thetak <= self._paras[i]['TA'] * pi and \
                                                    thetam >= self._paras[i]['TA'] * pi:
                                                tmp = self._paras[i]['IGW_TA'] * (
                                                    thetam * ipi - 0.8)
                                                tmp2 = cos(phi)
                                                ops[i] += gaussthetak[i] * \
                                                    exp(-0.5 * tmp * \
                                                    tmp) * tmp2 * tmp2
                                                norms[i][j] += 1
                                        elif t in ["T", "tri_pyr", "sq_pyr", \
                                                "pent_pyr", "hex_pyr"]:
                                            tmp = cos(self._paras[i]['fac_AA'] * \
                                                phi) ** self._paras[i]['exp_cos_AA']
                                            tmp2 = self._paras[i]['IGW_EP'] * (
                                                thetak * ipi - 0.5)
                                            tmp3 = self._paras[i]['IGW_EP'] * (
                                                thetam * ipi - 0.5)
                                            qsptheta[i][j] += tmp * exp(
                                                -0.5 * tmp2 * tmp2) * exp(
                                                -0.5 * tmp3 * tmp3)
                                            norms[i][j] += 1
                                        elif t in ["sq_plan", "oct", \
                                                "oct_legacy", "tri_bipyr", \
                                                "sq_bipyr", "pent_bipyr", \
                                                "hex_bipyr"]:
                                            if thetak < self._paras[i]['min_SPP'] and \
                                                    thetam < self._paras[i]['min_SPP']:
                                                tmp = cos(self._paras[i]['fac_AA'] * \
                                                    phi) ** self._paras[i]['exp_cos_AA']
                                                tmp2 = self._paras[i]['IGW_EP'] * (
                                                    thetam * ipi - 0.5)
                                                qsptheta[i][j] += tmp * \
                                                    exp(-0.5 * tmp2 * tmp2)
                                                if t == "oct_legacy":
                                                    qsptheta[i][j] -= tmp * \
                                                        self._paras[i][6] * \
                                                        self._paras[i][7]
                                                norms[i][j] += 1
                                        elif t == "bcc" and j < k:
                                            if thetak < self._paras[i]['min_SPP']:
                                                if thetak > piover2:
                                                    fac = 1.0
                                                else:
                                                    fac = -1.0
                                                tmp = (thetam - piover2) / (
                                                    19.47 * pi / 180.0)
                                                qsptheta[i][j] += fac * cos(
                                                    3.0 * phi) * 1.6 * tmp * \
                                                    exp(-0.5 * tmp * tmp)
                                        elif t == "see_saw":
                                            if thetak < self._paras[i]['min_SPP'] and \
                                                    thetam < self._paras[i]['min_SPP']:
                                                tmp = self._paras[i]['IGW_EP'] * (
                                                    phi * ipi - 0.5)
                                                tmp2 = self._paras[i]['IGW_EP'] * (
                                                    thetam * ipi - 0.5)
                                                qsptheta[i][j] += exp(
                                                    -0.5 * tmp * tmp) * \
                                                    exp(-0.5 * tmp2 * tmp2)
                                                norms[i][j] += 1.0
                                        elif t == "cuboct":
                                            if thetam < self._paras[i]['min_SPP'] and \
                                                    thetak > self._paras[i][4] and \
                                                    thetak < self._paras[i][2]:
                                                if thetam > self._paras[i][4] and \
                                                        thetam < self._paras[i][2]:
                                                    tmp = cos(phi)
                                                    tmp2 = self._paras[i][5] * (
                                                        thetam * ipi - 0.5)
                                                    qsptheta[i][j] += tmp * tmp * \
                                                        exp(
                                                        -0.5 * tmp2 * tmp2)
                                                    norms[i][j] += 1.0
                                                elif thetam < self._paras[i][4]:
                                                    tmp = 0.0556 * (
                                                        cos(phi-0.5*pi) - 0.81649658)
                                                    tmp2 = self._paras[i][6] * (
                                                        thetam * ipi - onethird)
                                                    qsptheta[i][j] += exp(-0.5 * tmp * tmp) * \
                                                        exp(
                                                        -0.5 * tmp2 * tmp2)
                                                    norms[i][j] += 1.0
                                                elif thetam > self._paras[i][2]:
                                                    tmp = 0.0556 * (
                                                        cos(phi-0.5*pi) - 0.81649658)
                                                    tmp2 = self._paras[i][6] * (
                                                        thetam * ipi - twothird)
                                                    qsptheta[i][j] += exp(-0.5 * tmp * tmp) * \
                                                        exp(
                                                        -0.5 * tmp2 * tmp2)
                                                    norms[i][j] += 1.0


            # Normalize Peters-style OPs.
            for i, t in enumerate(self._types):
                if t in ["tri_plan", "tet", \
                         "pent_plan"]:
                    ops[i] = ops[i] / sum(norms[i]) \
                        if sum(norms[i]) > 1.0e-12 else None
                elif t in ["bent", "sq_plan", "oct", "oct_legacy", \
                        "cuboct"]:
                    ops[i] = sum(qsptheta[i]) / sum(norms[i]) \
                        if sum(norms[i]) > 1.0e-12 else None
                elif t in ["T", "tri_pyr", "see_saw", "sq_pyr", "tri_bipyr", \
                           "sq_bipyr", "pent_pyr", "hex_pyr", "pent_bipyr", \
                           "hex_bipyr"]:
                    for j in range(nneigh):
                        qsptheta[i][j] = qsptheta[i][j] / norms[i][j] \
                            if norms[i][j] > 1.0e-12 else 0.0
                    ops[i] = max(qsptheta[i]) if len(qsptheta[i]) > 0 else None
                elif t == "bcc":
                    ops[i] = sum(qsptheta[i]) / float(0.5 * float(
                            nneigh * (6 + (nneigh - 2) * (nneigh - 3)))) \
                            if nneigh > 3 else None
                elif t == "sq_pyr_legacy":
                    if nneigh > 1:
                        dmean = np.mean(dist)
                        acc = 0.0
                        for d in dist:
                            tmp = self._paras[i][2] * (d - dmean)
                            acc = acc + exp(-0.5 * tmp * tmp)
                        ops[i] = acc * max(qsptheta[i]) / float(
                                nneigh * (nneigh - 1))
                    else:
                        ops[i] = None

        # Then, deal with the new-style OPs that require vectors between
        # neighbors.
        if self._geomops2:
            # Compute all (unique) angles and sort the resulting list.
            aij = []
            for ir, r in enumerate(rijnorm):
                for j in range(ir+1, len(rijnorm)):
                    aij.append(acos(max(-1.0, min(np.inner(
                            r, rijnorm[j]), 1.0))))
            aijs = sorted(aij)

            # Compute height, side and diagonal length estimates.
            neighscent = np.array([0.0, 0.0, 0.0])
            for j, neigh in enumerate(neighsites):
                neighscent = neighscent + neigh.coords
            if nneigh > 0:
                neighscent = (neighscent / float(nneigh))
            h = np.linalg.norm(neighscent - centvec)
            b = min(distjk_unique) if len(distjk_unique) > 0 else 0
            dhalf = max(distjk_unique) / 2.0 if len(distjk_unique) > 0 else 0

            for i, t in enumerate(self._types):
                if t == "reg_tri" or t == "sq":
                    if nneigh < 3:
                        ops[i] = None
                    else:
                        ops[i] = 1.0
                        if t == "reg_tri":
                            a = 2.0 * asin(b / (2.0 * sqrt(h*h + (b / (
                                    2.0 * cos(3.0 * pi / 18.0)))**2.0)))
                            nmax = 3
                        elif t == "sq":
                            a = 2.0 * asin(b / (2.0 * sqrt(h*h + dhalf*dhalf)))
                            nmax = 4
                        for j in range(min([nneigh,nmax])):
                            ops[i] = ops[i] * exp(-0.5 * ((
                                    aijs[j] - a) * self._paras[i][0])**2)

        return ops
