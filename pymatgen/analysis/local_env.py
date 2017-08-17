# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import six
import ruamel.yaml as yaml
import os

"""
This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures.
To do:
- Insert LocalStructOrderParas class here.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman,"+\
    " Nils E. R. Zimmermann"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Nils E. R. Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__status__ = "Production"
__date__ = "August 17, 2017"

from math import pow, pi, asin, atan, sqrt, exp, cos, acos
import numpy as np

from scipy.spatial import Voronoi
from pymatgen import Element
from pymatgen.util.num import abs_cap
from pymatgen.analysis.bond_valence import BV_PARAMS
from pymatgen.analysis.defects.point_defects import \
        ValenceIonicRadiusEvaluator


class NearNeighbors(object):
    """
    Base class to determine near neighbors that typically include nearest
    neighbors and others that are within some tolerable distance.
    """

    def __init__(self):
        pass

    def get_cn(self, structure, n, use_weights=False):
        """
        Get coordination number, CN, of site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine CN.
            use_weights (boolean): flag indicating whether (True)
                to use weights for computing the coordination number
                or not (False, default: each coordinated site has equal
                weight).
        Returns:
            cn (integer or float): coordination number.
        """

        siw = self.get_nn_info(structure, n)
        return sum([e['weight'] for e in siw]) if use_weights else len(siw)

    def get_nn(self, structure, n):
        """
        Get near neighbors of site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (integer): index of site in structure for which to determine
                    neighbors.
        Returns:
            sites (list of Site objects): near neighbors.
        """

        return [e['site'] for e in self.get_nn_info(structure, n)]

    def get_weights_of_nn_sites(self, n):
        """
        Get weight associated with each near neighbor of site with
        index n in structure.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine the weights.
        Returns:
            weights (list of floats): near-neighbor weights.
        """

        return [e['weight'] for e in self.get_nn_info(structure, n)]

    def get_nn_images(self, structure, n):
        """
        Get image location of all near neighbors of site with index n in
        structure.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine the image
                location of near neighbors.
        Returns:
            images (list of 3D integer array): image locations of
                near neighbors as lattice translational vectors.
        """

        return [e['image'] for e in self.get_nn_info(structure, n)]

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near-neighbor
                information.
 
        Returns:
            siw (list of dicts): each dictionary provides information
                about a single near neighbor, where key 'site' gives
                access to the corresponding Site object, 'image' gives
                the image location (lattice translation vector), and
                'weight' provides the weight that a given near-neighbor
                site contributes to the coordination number (1 or smaller).
        """

        raise NotImplementedError("get_nn_info(structure, n)"
                " is not defined!")

class VoronoiNN(NearNeighbors):
    """
    Uses a Voronoi algorithm to determine near neighbors for each site in a
    structure.

    Args:
        tol (float): tolerance parameter for near-neighbor finding.
        targets (Element or list of Elements): target element(s).
        cutoff (float): cutoff radius in Angstrom to look for near-neighbor
            atoms. Defaults to 10.0.
        allow_pathological (bool): whether to allow infinite vertices in
            determination of Voronoi coordination.
    """

    def __init__(self, tol=None, targets=None, cutoff=10.0,
                 allow_pathological=False):
        self.tol = tol
        self.cutoff = cutoff
        self.allow_pathological = allow_pathological
        self.targets = targets

    def get_voronoi_polyhedra(self, structure, n):
        """
        Gives a weighted polyhedra around a site. This uses the Voronoi
        construction with solid angle weights.
        See ref: A Proposed Rigorous Definition of Coordination Number,
        M. O'Keeffe, Acta Cryst. (1979). A35, 772-775

        Args:
            structure (Structure): structure for which to evaluate the
                coordination environment.
            n (integer): site index.

        Returns:
            A dict of sites sharing a common Voronoi facet with the site
            n and their solid angle weights
        """
        if self.targets is None:
            targets = structure.composition.elements
        else:
            targets = self.targets
        center = structure[n]
        neighbors = structure.get_sites_in_sphere(
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
                if nn.specie in targets:
                    resultweighted[nn] = angle / maxangle
            else:  # is nn site is disordered
                for disordered_sp in nn.species_and_occu.keys():
                    if disordered_sp in targets:
                        resultweighted[nn] = angle / maxangle

        return resultweighted


    def get_nn_info(self, structure, n):
        """"
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure
        using Voronoi decomposition.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near-neighbor
                sites.
 
        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """

        if self.targets is None:
            targets = structure.composition.elements
        else:
            targets = self.targets
        siw = []
        for site, weight in self.get_voronoi_polyhedra(
                structure, n).items():
            if weight > self.tol and site.specie in targets:
                dist, image = structure.sites[n].distance_and_image(site)
                siw.append({'site': site, 'image': image, 'weight': weight})
        return siw


class JMolNN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using an emulation
    of JMol's default autoBond() algorithm. This version of the algorithm
    does not take into account any information regarding known charge
    states.

    Args:
        tol (float): tolerance parameter for bond determination
            (default: 1E-3).
        el_radius_updates: (dict) symbol->float to override default atomic 
            radii table values 
    """

    def __init__(self, tol=None, el_radius_updates=None):

        self.tol = 1E-3 if tol is None else tol

        # Load elemental radii table
        bonds_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "bonds_jmol_ob.yaml")
        with open(bonds_file, 'r') as f:
            self.el_radius = yaml.safe_load(f)

        # Update any user preference elemental radii
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

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the bond identification
        algorithm underlying JMol.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near
                neighbors.
 
        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """

        site = structure[n]

        # Determine relevant bond lengths based on atomic radii table
        bonds = {}
        for el in structure.composition.elements:
            bonds[site.specie, el] = self.get_max_bond_distance(
                site.specie.symbol, el.symbol)

        # Search for neighbors up to max bond length + tolerance
        max_rad = max(bonds.values()) + self.tol
        min_rad = min(bonds.values())

        siw = []
        for neighb, dist in structure.get_neighbors(site, max_rad):
            # Confirm neighbor based on bond length specific to atom pair
            if dist <= bonds[(site.specie, neighb.specie)] + self.tol:
                weight = min_rad / dist
                d, image = site.distance_and_image(neighb)
                siw.append({'site': neighb, 'image': image, 'weight': weight})
        return siw


class MinimumDistanceNN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    nearest neighbor(s) at distance, d_min, plus all neighbors
    within a distance (1 + delta) * d_min, where delta is a
    (relative) distance tolerance parameter.

    Args:
        tol (float): tolerance parameter for neighbor identification
            (default: 0.1).
        cutoff (float): cutoff radius in Angstrom to look for trial
            near-neighbor sites (default: 10.0).
    """

    def __init__(self, tol=0.1, cutoff=10.0):

        self.tol = tol
        self.cutoff = cutoff

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest neighbor
        distance-based method.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near
                neighbors.
 
        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """

        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        min_dist = min([dist for neigh, dist in neighs_dists])

        siw = []
        for s, dist in neighs_dists:
            if dist < (1.0 + self.tol) * min_dist:
                w = min_dist / dist
                d, i = site.distance_and_image(s)
                siw.append({'site': s, 'image': i, 'weight': w})
        return siw


class MinimumOKeeffeNN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    neighbor(s) at closest relative distance, d_min_OKeffee, plus some
    relative tolerance, where bond valence parameters from O'Keeffe's
    bond valence method (J. Am. Chem. Soc. 1991, 3226-3229) are used
    to calculate relative distances.

    Args:
        tol (float): tolerance parameter for neighbor identification
            (default: 0.1).
        cutoff (float): cutoff radius in Angstrom to look for trial
            near-neighbor sites (default: 10.0).
    """

    def __init__(self, tol=0.1, cutoff=10.0):

        self.tol = tol
        self.cutoff = cutoff

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest relative
        neighbor distance-based method with O'Keeffe parameters.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near
                neighbors.
 
        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """

        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        try:
            eln = site.specie.element
        except:
            eln = site.species_string

        reldists_neighs = []
        for neigh, dist in neighs_dists:
            try:
                el2 = neigh.specie.element
            except:
                el2 = neigh.species_string
            reldists_neighs.append([dist / get_okeeffe_distance_prediction(
                    eln, el2), neigh])

        siw = []
        min_reldist = min([reldist for reldist, neigh in reldists_neighs])
        for reldist, s in reldists_neighs:
            if reldist < (1.0 + self.tol) * min_reldist:
                w = min_reldist / reldist
                d, i = site.distance_and_image(s)
                siw.append({'site': s, 'image': i, 'weight': w})

        return siw


class MinimumVIRENN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    neighbor(s) at closest relative distance, d_min_VIRE, plus some
    relative tolerance, where atom radii from the
    Pymatgen's ValenceIonicRadiusEvaluator (VIRE) are used
    to calculate relative distances.

    Args:
        tol (float): tolerance parameter for neighbor identification
            (default: 0.1).
        cutoff (float): cutoff radius in Angstrom to look for trial
            near-neighbor sites (default: 10.0).
    """

    def __init__(self, tol=0.1, cutoff=10.0):

        self.tol = tol
        self.cutoff = cutoff

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest relative
        neighbor distance-based method with VIRE atomic/ionic radii.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near
                neighbors.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """

        vire = ValenceIonicRadiusEvaluator(structure)
        site = vire.structure[n]
        neighs_dists = vire.structure.get_neighbors(site, self.cutoff)
        rn = vire.radii[vire.structure[n].species_string]

        reldists_neighs = []
        for neigh, dist in neighs_dists:
            reldists_neighs.append([dist / (
                    vire.radii[neigh.species_string] + rn), neigh])

        siw = []
        min_reldist = min([reldist for reldist, neigh in reldists_neighs])
        for reldist, s in reldists_neighs:
            if reldist < (1.0 + self.tol) * min_reldist:
                w = min_reldist / reldist
                d, i = site.distance_and_image(s)
                siw.append({'site': s, 'image': i, 'weight': w})

        return siw


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


def get_okeeffe_params(el_symbol):
    """
    Returns the elemental parameters related to atom size and
    electronegativity which are used for estimating bond-valence
    parameters (bond length) of pairs of atoms on the basis of data
    provided in 'Atoms Sizes and Bond Lengths in Molecules and Crystals'
    (O'Keeffe & Brese, 1991).

    Args:
        el_symbol (str): element symbol.
    Returns:
        (dict): atom-size ('r') and electronegativity-related ('c')
                parameter.
    """

    el = Element(el_symbol)
    if el not in list(BV_PARAMS.keys()):
        raise RuntimeError("Could not find O'Keeffe parameters for element"
                           " \"{}\" in \"BV_PARAMS\"dictonary"
                           " provided by pymatgen".format(el_symbol))

    return BV_PARAMS[el]


def get_okeeffe_distance_prediction(el1, el2):
    """
    Returns an estimate of the bond valence parameter (bond length) using
    the derived parameters from 'Atoms Sizes and Bond Lengths in Molecules
    and Crystals' (O'Keeffe & Brese, 1991). The estimate is based on two
    experimental parameters: r and c. The value for r  is based off radius,
    while c is (usually) the Allred-Rochow electronegativity. Values used
    are *not* generated from pymatgen, and are found in
    'okeeffe_params.json'.

    Args:
        el1, el2 (Element): two Element objects
    Returns:
        a float value of the predicted bond length
    """
    el1_okeeffe_params = get_okeeffe_params(el1)
    el2_okeeffe_params = get_okeeffe_params(el2)

    r1 = el1_okeeffe_params['r']
    r2 = el2_okeeffe_params['r']
    c1 = el1_okeeffe_params['c']
    c2 = el2_okeeffe_params['c']

    return r1 + r2 - r1 * r2 * pow(
            sqrt(c1) - sqrt(c2), 2) / (c1 * r1 + c2 * r2)

