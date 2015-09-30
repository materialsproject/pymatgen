# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

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
import numpy as np
import itertools
import collections

from warnings import warn
from pyhull.voronoi import VoronoiTess
from pymatgen import PeriodicSite
from pymatgen import Element, Specie, Composition
from pymatgen.util.num_utils import abs_cap


class VoronoiCoordFinder(object):
    """
    Uses a Voronoi algorithm to determine the coordination for each site in a
    structure.

    Args:
        structure (Structure): Input structure
        target ([Element/Specie]): A list of target species to determine
            coordination for.
    """

    # Radius in Angstrom cutoff to look for coordinating atoms
    default_cutoff = 10.0

    def __init__(self, structure, target=None):
        self._structure = structure
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
            center.coords, VoronoiCoordFinder.default_cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        qvoronoi_input = [s.coords for s in neighbors]
        voro = VoronoiTess(qvoronoi_input)
        all_vertices = voro.vertices

        results = {}
        for nn, vind in voro.ridges.items():
            if 0 in nn:
                if 0 in vind:
                    raise RuntimeError("This structure is pathological,"
                                       " infinite vertex in the voronoi "
                                       "construction")

                facets = [all_vertices[i] for i in vind]
                results[neighbors[nn[1]]] = solid_angle(center.coords, facets)

        maxangle = max(results.values())

        resultweighted = {}
        for nn, angle in results.items():
            if nn.specie in localtarget:
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
        i = np.ceil(cutoff * recp_len / (2 * math.pi))
        offsets = np.mgrid[-i[0]:i[0] + 1, -i[1]:i[1] + 1, -i[2]:i[2] + 1].T
        self.offsets = np.reshape(offsets, (-1, 3))
        #shape = [image, axis]
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
        #shape = [site, axis]
        cart_coords = np.array(self.s.cart_coords)
        #shape = [site, image, axis]
        all_sites = cart_coords[:, None, :] + self.cart_offsets[None, :, :]
        vt = VoronoiTess(all_sites.reshape((-1, 3)))
        n_images = all_sites.shape[1]
        cs = (len(self.s), len(self.s), len(self.cart_offsets))
        connectivity = np.zeros(cs)
        vts = np.array(vt.vertices)
        for (ki, kj), v in vt.ridges.items():
            atomi = ki // n_images
            atomj = kj // n_images

            imagei = ki % n_images
            imagej = kj % n_images

            if imagei != n_images // 2 and imagej != n_images // 2:
                continue

            if imagei == n_images // 2:
                #atomi is in original cell
                val = solid_angle(vt.points[ki], vts[v])
                connectivity[atomi, atomj, imagej] = val

            if imagej == n_images // 2:
                #atomj is in original cell
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
        v = -np.dot(n[i], n[i + 1])\
            / (np.linalg.norm(n[i]) * np.linalg.norm(n[i + 1]))
        vals.append(math.acos(abs_cap(v)))
    phi = sum(vals)
    return phi + (3 - len(r)) * math.pi


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
            syms = [sp. symbol for sp in site.species_and_occu.keys()]
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
