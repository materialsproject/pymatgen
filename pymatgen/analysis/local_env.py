# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import six
import ruamel.yaml as yaml
import os
import json

"""
This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures.
To do:
- Insert LocalStructOrderParas class here.
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman,"+\
    " Nils E. R. Zimmermann, Bharat Medasani"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Nils E. R. Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__status__ = "Production"
__date__ = "August 17, 2017"

from math import pow, pi, asin, atan, sqrt, exp, cos, acos
import numpy as np

from bisect import bisect_left
from scipy.spatial import Voronoi
from pymatgen import Element
from pymatgen.core.structure import Structure
from pymatgen.util.num import abs_cap
from pymatgen.analysis.bond_valence import BV_PARAMS
from pymatgen.analysis.structure_analyzer import OrderParameters


file_dir = os.path.dirname(__file__)
rad_file = os.path.join(file_dir, 'ionic_radii.json')
with open(rad_file, 'r') as fp:
    _ion_radii = json.load(fp)


class ValenceIonicRadiusEvaluator(object):
    """
    Computes site valences and ionic radii for a structure using bond valence
    analyzer

    Args:
        structure: pymatgen.core.structure.Structure
    """

    def __init__(self, structure):
        self._structure = structure.copy()
        self._valences = self._get_valences()
        self._ionic_radii = self._get_ionic_radii()

    @property
    def radii(self):
        """
        List of ionic radii of elements in the order of sites.
        """
        el = [site.species_string for site in self._structure.sites]
        radii_dict = dict(zip(el, self._ionic_radii))
        #print radii_dict
        return radii_dict

    @property
    def valences(self):
        """
        List of oxidation states of elements in the order of sites.
        """
        el = [site.species_string for site in self._structure.sites]
        valence_dict = dict(zip(el, self._valences))
        return valence_dict

    @property
    def structure(self):
        """
        Returns oxidation state decorated structure.
        """
        return self._structure.copy()


    def _get_ionic_radii(self):
        """
        Computes ionic radii of elements for all sites in the structure.
        If valence is zero, atomic radius is used.
        """
        radii = []
        vnn = VoronoiNN() # self._structure)

        def nearest_key(sorted_vals, key):
            i = bisect_left(sorted_vals, key)
            if i == len(sorted_vals):
                return sorted_vals[-1]
            if i == 0:
                return sorted_vals[0]
            before = sorted_vals[i-1]
            after = sorted_vals[i]
            if after-key < key-before:
                return after
            else:
                return before

        for i in range(len(self._structure.sites)):
            site = self._structure.sites[i]
            if isinstance(site.specie,Element):
                radius = site.specie.atomic_radius
                # Handle elements with no atomic_radius
                # by using calculated values instead.
                if radius is None:
                    radius = site.specie.atomic_radius_calculated
                if radius is None:
                    raise ValueError(
                            "cannot assign radius to element {}".format(
                            site.specie))
                radii.append(radius)
                continue

            el = site.specie.symbol
            oxi_state = int(round(site.specie.oxi_state))
            coord_no = int(round(vnn.get_cn(self._structure, i)))
            try:
                tab_oxi_states = sorted(map(int, _ion_radii[el].keys()))
                oxi_state = nearest_key(tab_oxi_states, oxi_state)
                radius = _ion_radii[el][str(oxi_state)][str(coord_no)]
            except KeyError:
                if vnn.get_cn(self._structure, i)-coord_no > 0:
                    new_coord_no = coord_no + 1
                else:
                    new_coord_no = coord_no - 1
                try:
                    radius = _ion_radii[el][str(oxi_state)][str(new_coord_no)]
                    coord_no = new_coord_no
                except:
                    tab_coords = sorted(map(int, _ion_radii[el][str(oxi_state)].keys()))
                    new_coord_no = nearest_key(tab_coords, coord_no)
                    i = 0
                    for val in tab_coords:
                        if  val > coord_no:
                            break
                        i = i + 1
                    if i == len(tab_coords):
                        key = str(tab_coords[-1])
                        radius = _ion_radii[el][str(oxi_state)][key]
                    elif i == 0:
                        key = str(tab_coords[0])
                        radius = _ion_radii[el][str(oxi_state)][key]
                    else:
                        key = str(tab_coords[i-1])
                        radius1 = _ion_radii[el][str(oxi_state)][key]
                        key = str(tab_coords[i])
                        radius2 = _ion_radii[el][str(oxi_state)][key]
                        radius = (radius1+radius2)/2

            #implement complex checks later
            radii.append(radius)
        return radii

    def _get_valences(self):
        """
        Computes ionic valences of elements for all sites in the structure.
        """
        try:
            bv = BVAnalyzer()
            self._structure = bv.get_oxi_state_decorated_structure(self._structure)
            valences = bv.get_valences(self._structure)
        except:
            try:
                bv = BVAnalyzer(symm_tol=0.0)
                self._structure = bv.get_oxi_state_decorated_structure(self._structure)
                valences = bv.get_valences(self._structure)
            except:
                valences = []
                for site in self._structure.sites:
                    if len(site.specie.common_oxidation_states) > 0:
                        valences.append(site.specie.common_oxidation_states[0])
                    # Handle noble gas species
                    # which have no entries in common_oxidation_states.
                    else:
                        valences.append(0)
                if sum(valences):
                    valences = [0]*self._structure.num_sites
                else:
                    self._structure.add_oxidation_state_by_site(valences)
                #raise

        #el = [site.specie.symbol for site in self._structure.sites]
        #el = [site.species_string for site in self._structure.sites]
        #el = [site.specie for site in self._structure.sites]
        #valence_dict = dict(zip(el, valences))
        #print valence_dict
        return valences


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

    def get_weights_of_nn_sites(self, structure, n):
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
                near neighbors.
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
                the image location, and 'weight' provides the weight
                that a given near-neighbor site contributes
                to the coordination number (1 or smaller), 'site_index'
                gives index of the corresponding site in
                the original structure.
        """

        raise NotImplementedError("get_nn_info(structure, n)"
                " is not defined!")

    @staticmethod
    def _get_image(frac_coords):
        """Private convenience method for get_nn_info,
        gives lattice image from provided PeriodicSite."""
        return [int(f) if f >= 0 else int(f - 1)
                for f in frac_coords]

    @staticmethod
    def _get_original_site(structure, site):
        """Private convenience method for get_nn_info,
        gives original site index from ProvidedPeriodicSite."""
        is_periodic_image = [site.is_periodic_image(s) for s in structure]
        return is_periodic_image.index(True)

class VoronoiNN(NearNeighbors):
    """
    Uses a Voronoi algorithm to determine near neighbors for each site in a
    structure.

    Args:
        tol (float): tolerance parameter for near-neighbor finding
            (default: 0).
        targets (Element or list of Elements): target element(s).
        cutoff (float): cutoff radius in Angstrom to look for near-neighbor
            atoms. Defaults to 10.0.
        allow_pathological (bool): whether to allow infinite vertices in
            determination of Voronoi coordination.
    """

    def __init__(self, tol=0, targets=None, cutoff=10.0,
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
                siw.append({'site': site,
                            'image': self._get_image(site.frac_coords),
                            'weight': weight,
                            'site_index': self._get_original_site(structure, site)})
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

    def __init__(self, tol=1E-3, el_radius_updates=None):

        self.tol = tol

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
                siw.append({'site': neighb,
                            'image': self._get_image(neighb.frac_coords),
                            'weight': weight,
                            'site_index': self._get_original_site(structure, neighb)})
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
                siw.append({'site': s,
                            'image': self._get_image(s.frac_coords),
                            'weight': w,
                            'site_index': self._get_original_site(structure, s)})
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
                siw.append({'site': s,
                            'image': self._get_image(s.frac_coords),
                            'weight': w,
                            'site_index': self._get_original_site(structure, s)})

        return siw


class MinimumVIRENN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    neighbor(s) at closest relative distance, d_min_VIRE, plus some
    relative tolerance, where atom radii from the
    ValenceIonicRadiusEvaluator (VIRE) are used
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
                siw.append({'site': s,
                            'image': self._get_image(s.frac_coords),
                            'weight': w,
                            'site_index': self._get_original_site(structure, s)})

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


def get_neighbors_of_site_with_index(struct, n, approach="min_dist", delta=0.1, \
        cutoff=10.0):
    """
    Returns the neighbors of a given site using a specific neighbor-finding
    method.

    Args:
        struct (Structure): input structure.
        n (int): index of site in Structure object for which motif type
                is to be determined.
        approach (str): type of neighbor-finding approach, where
              "min_dist" will use the MinimumDistanceNN class,
              "voronoi" the VoronoiNN class, "min_OKeeffe" the
              MinimumOKeeffe class, and "min_VIRE" the MinimumVIRENN class.
        delta (float): tolerance involved in neighbor finding.
        cutoff (float): (large) radius to find tentative neighbors.

    Returns: neighbor sites.
    """

    if approach == "min_dist":
        return MinimumDistanceNN(tol=delta, cutoff=cutoff).get_nn(
                struct, n)
    elif approach == "voronoi":
        return VoronoiNN(tol=delta, cutoff=cutoff).get_nn(
                struct, n)
    elif approach == "min_OKeeffe":
        return MinimumOKeeffeNN(tol=delta, cutoff=cutoff).get_nn(
                struct, n)
    elif approach == "min_VIRE":
        return MinimumVIRENN(tol=delta, cutoff=cutoff).get_nn(
                struct, n)
    else:
        raise RuntimeError("unsupported neighbor-finding method ({}).".format(
                approach))


def site_is_of_motif_type(struct, n, approach="min_dist", delta=0.1, \
        cutoff=10.0, thresh=None):
    """
    Returns the motif type of the site with index n in structure struct;
    currently featuring "tetrahedral", "octahedral", "bcc", and "cp"
    (close-packed: fcc and hcp) as well as "square pyramidal" and
    "trigonal bipyramidal".  If the site is not recognized,
    "unrecognized" is returned.  If a site should be assigned to two
    different motifs, "multiple assignments" is returned.

    Args:
        struct (Structure): input structure.
        n (int): index of site in Structure object for which motif type
                is to be determined.
        approach (str): type of neighbor-finding approach, where
              "min_dist" will use the MinimumDistanceNN class,
              "voronoi" the VoronoiNN class, "min_OKeeffe" the
              MinimumOKeeffe class, and "min_VIRE" the MinimumVIRENN class.
        delta (float): tolerance involved in neighbor finding.
        cutoff (float): (large) radius to find tentative neighbors.
        thresh (dict): thresholds for motif criteria (currently, required
                keys and their default values are "qtet": 0.5,
                "qoct": 0.5, "qbcc": 0.5, "q6": 0.4).

    Returns: motif type (str).
    """

    if thresh is None:
        thresh = {
            "qtet": 0.5, "qoct": 0.5, "qbcc": 0.5, "q6": 0.4,
            "qtribipyr": 0.8, "qsqpyr": 0.8}

    ops = OrderParameters([
            "cn", "tet", "oct", "bcc", "q6", "sq_pyr", "tri_bipyr"])

    neighs_cent = get_neighbors_of_site_with_index(
            struct, n, approach=approach, delta=delta, cutoff=cutoff)
    neighs_cent.append(struct.sites[n])
    opvals = ops.get_order_parameters(
            neighs_cent, len(neighs_cent)-1, indices_neighs=[
            i for i in range(len(neighs_cent)-1)])
    cn = int(opvals[0] + 0.5)
    motif_type = "unrecognized"
    nmotif = 0

    if cn == 4 and opvals[1] > thresh["qtet"]:
        motif_type = "tetrahedral"
        nmotif += 1
    if cn == 5 and opvals[5] > thresh["qsqpyr"]:
       motif_type = "square pyramidal"
       nmotif += 1
    if cn == 5 and opvals[6] > thresh["qtribipyr"]:
       motif_type = "trigonal bipyramidal"
       nmotif += 1
    if cn == 6 and opvals[2] > thresh["qoct"]:
        motif_type = "octahedral"
        nmotif += 1
    if cn == 8 and (opvals[3] > thresh["qbcc"] and opvals[1] < thresh["qtet"]):
        motif_type = "bcc"
        nmotif += 1
    if cn == 12 and (opvals[4] > thresh["q6"] and opvals[1] < thresh["q6"] and \
                                 opvals[2] < thresh["q6"] and opvals[3] < thresh["q6"]):
        motif_type = "cp"
        nmotif += 1

    if nmotif > 1:
        motif_type = "multiple assignments"

    return motif_type
