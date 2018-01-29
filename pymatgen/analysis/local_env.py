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
"""

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman,"+\
    " Nils E. R. Zimmermann, Bharat Medasani"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Nils E. R. Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__status__ = "Production"
__date__ = "August 17, 2017"

from math import pow, pi, asin, atan, sqrt, exp, sin, cos, acos, fabs
import numpy as np

from bisect import bisect_left
from scipy.spatial import Voronoi
from pymatgen import Element
from pymatgen.core.structure import Structure
from pymatgen.util.num import abs_cap
from pymatgen.analysis.bond_valence import BV_PARAMS


default_op_paras = {}
with open(os.path.join(os.path.dirname(
        __file__), 'op_paras.yaml'), "rt") as f:
    default_op_paras = yaml.safe_load(f)
    f.close()

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
        vnn = VoronoiNN()

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

    ops = LocalStructOrderParas([
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


class LocalStructOrderParas(object):
    """
    This class permits the calculation of various types of local
    structure order parameters.
    """

    __supported_types = (
        "cn", "sgl_bd", "bent", "tri_plan", "tri_plan_max", "reg_tri", "sq_plan", \
        "sq_plan_max", "pent_plan", "pent_plan_max", "sq", "tet", "tet_max", "tri_pyr", \
        "sq_pyr", "sq_pyr_legacy", "tri_bipyr", "sq_bipyr", "oct", \
        "oct_legacy", "pent_pyr", "hex_pyr", "pent_bipyr", "hex_bipyr", \
        "T", "cuboct", "cuboct_max", "see_saw_rect", "bcc", "q2", "q4", "q6", "oct_max", "hex_plan_max")

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
                  "see_saw_rect": see saw-like coordinations;
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
                      equatorial plane (EP) at 90 degrees ("T", "see_saw_rect",
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
                      at South pole position ("see_saw_rect", "oct", "bcc",
                      "sq_plan", "tri_bipyr", "sq_bipyr", "pent_bipyr",
                      "hex_bipyr", "cuboct", and "oct_legacy"
                      (2.792526803190927)).
                  "IGW_SPP": IGW for penalizing angles away from South
                      pole position ("see_saw_rect", "oct", "bcc", "sq_plan",
                      "tri_bipyr", "sq_bipyr", "pent_bipyr", "hex_bipyr",
                      "cuboct", and "oct_legacy" (15)).
                  "w_SPP": weight for South pole position relative to
                      equatorial positions ("see_saw_rect" and "sq_plan" (1),
                      "cuboct" (1.8), "tri_bipyr" (2), "oct",
                      "sq_bipyr", and "oct_legacy" (3), "pent_bipyr" (4),
                      "hex_bipyr" (5), "bcc" (6)).
            cutoff (float): Cutoff radius to determine which nearest
                neighbors are supposed to contribute to the order
                parameters. If the value is negative the neighboring
                sites found by distance and cutoff radius are further
                pruned using the get_nn method from the
                VoronoiNN class.
        """
        for t in types:
            if t not in LocalStructOrderParas.__supported_types:
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
                ["tet", "oct", "bcc", "sq_pyr", "sq_pyr_legacy",
                 "tri_bipyr", "sq_bipyr", "oct_legacy", "tri_plan",
                 "sq_plan", "pent_plan",  "tri_pyr", "pent_pyr", "hex_pyr",
                 "pent_bipyr", "hex_bipyr", "T", "cuboct", "oct_max", "tet_max",
                 "tri_plan_max", "sq_plan_max", "pent_plan_max", "cuboct_max",
                 "bent", "see_saw_rect", "hex_plan_max"]):
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
        pre_y_6_5 = [i32_3 * sqrt_1001_pi * val[0] * val[1]
                     for val in zip(self._pow_sin_t[5], self._pow_cos_t[1])]
        pre_y_6_4 = [i32_3 * sqrt_91_2pi * val[0] * (11.0 * val[1] - 1.0)
                     for val in zip(self._pow_sin_t[4], self._pow_cos_t[2])]
        pre_y_6_3 = [
            i32 * sqrt_1365_pi * val[0] * (11.0 * val[1] - 3.0 * val[2])
            for val in zip(self._pow_sin_t[3], self._pow_cos_t[3],
                           self._pow_cos_t[1])]
        pre_y_6_2 = [i64 * sqrt_1365_pi * val[0] * (33.0 * val[1] -
                                                    18.0 * val[2] + 1.0) for val
                     in zip(self._pow_sin_t[2],
                            self._pow_cos_t[4], self._pow_cos_t[2])]
        pre_y_6_1 = [i16 * sqrt_273_2pi * val[0] * (33.0 * val[1] -
                                                    30.0 * val[2] + 5.0 * val[
                                                        3]) for val in
                     zip(self._pow_sin_t[1],
                         self._pow_cos_t[5], self._pow_cos_t[3],
                         self._pow_cos_t[1])]

        acc = 0.0

        # Y_6_-6
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i]  # cos(x) =  cos(-x)
            imag -= pre_y_6_6[i] * self._sin_n_p[6][i]  # sin(x) = -sin(-x)
        acc += (real * real + imag * imag)

        # Y_6_-5
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += (real * real + imag * imag)

        # Y_6_-4
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += (real * real + imag * imag)

        # Y_6_-3
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += (real * real + imag * imag)

        # Y_6_-2
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_6_-1
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_6_0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += i32 * sqrt_13_pi * (231.0 * self._pow_cos_t[6][i] -
                                        315.0 * self._pow_cos_t[4][i] + 105.0 *
                                        self._pow_cos_t[2][i] - 5.0)
        acc += (real * real)

        # Y_6_1
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += (real * real + imag * imag)

        # Y_6_2
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag += pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += (real * real + imag * imag)

        # Y_6_3
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += (real * real + imag * imag)

        # Y_6_4
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag += pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += (real * real + imag * imag)

        # Y_6_5
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += (real * real + imag * imag)

        # Y_6_6
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
        the parameters associated
        with calculation of the order
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
        fac_bcc = 1.0 / exp(-0.5)

        # Find central site and its neighbors.
        # Note that we adopt the same way of accessing sites here as in
        # VoronoiNN; that is, not via the sites iterator.
        centsite = structure[n]
        if indices_neighs is not None:
            neighsites = [structure[index] for index in indices_neighs]
        elif self._voroneigh:
            vnn = VoronoiNN(tol=tol, targets=target_spec)
            neighsites = vnn.get_nn(structure, n)
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
        #norms = [[[] for j in range(nneigh)] for t in self._types]

        # First, coordination number and distance-based OPs.
        for i, t in enumerate(self._types):
            if t == "cn":
                ops[i] = nneigh / self._paras[i]['norm']
            elif t == "sgl_bd":
                dist_sorted = sorted(dist)
                if len(dist_sorted) == 1:
                    ops[i] = 1.0
                elif len(dist_sorted) > 1:
                    ops[i] = 1.0 - dist_sorted[0] / dist_sorted[1]

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
                if -left_of_unity < vec[2] < left_of_unity:
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
            qsptheta = [[[] for j in range(nneigh)] for t in self._types]
            norms = [[[] for j in range(nneigh)] for t in self._types]
            ipi = 1.0 / pi
            piover2 = pi / 2.0
            tetangoverpi = acos(-1.0 / 3.0) * ipi  # xxx: delete
            itetangminuspihalfoverpi = 1.0 / (tetangoverpi - 0.5)
            onethird = 1.0 / 3.0
            twothird = 2.0 / 3.0
            for j in range(nneigh):  # Neighbor j is put to the North pole.
                zaxis = rijnorm[j]
                kc = 0
                for k in range(nneigh):  # From neighbor k, we construct
                    if j != k:  # the prime meridian.
                        for i in range(len(self._types)):
                            qsptheta[i][j].append(0.0)
                            norms[i][j].append(0)
                        tmp = max(
                            -1.0, min(np.inner(zaxis, rijnorm[k]), 1.0))
                        thetak = acos(tmp)
                        xaxistmp = gramschmidt(rijnorm[k], zaxis)
                        if np.linalg.norm(xaxistmp) < very_small:
                            flag_xaxis = True
                        else:
                            xaxis = xaxistmp / np.linalg.norm(xaxistmp)
                            flag_xaxis = False

                        # Contributions of j-i-k angles, where i represents the
                        # central atom and j and k two of the neighbors.
                        for i, t in enumerate(self._types):
                            if t in ["bent", "sq_pyr_legacy"]:
                                tmp = self._paras[i]['IGW_TA'] * (
                                        thetak * ipi - self._paras[i]['TA'])
                                qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                norms[i][j][kc] += 1
                            elif t in ["tri_plan", "tri_plan_max", "tet", "tet_max"]:
                                tmp = self._paras[i]['IGW_TA'] * (
                                        thetak * ipi - self._paras[i]['TA'])
                                gaussthetak[i] = exp(-0.5 * tmp * tmp)
                                if t in ["tri_plan_max", "tet_max"]:
                                    qsptheta[i][j][kc] += gaussthetak[i]
                                    norms[i][j][kc] += 1
                            elif t in ["T", "tri_pyr", "sq_pyr", "pent_pyr", "hex_pyr"]:
                                tmp = self._paras[i]['IGW_EP'] * (thetak * ipi - 0.5)
                                qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                norms[i][j][kc] += 1
                            elif t in ["sq_plan", "oct", "oct_legacy",
                                       "cuboct", "cuboct_max"]:
                                if thetak >= self._paras[i]['min_SPP']:
                                    tmp = self._paras[i]['IGW_SPP'] * (
                                            thetak * ipi - 1.0)
                                    qsptheta[i][j][kc] += (
                                            self._paras[i]['w_SPP'] *
                                            exp(-0.5 * tmp * tmp))
                                    norms[i][j][kc] += self._paras[i]['w_SPP']
                            elif t in ["see_saw_rect", "tri_bipyr", "sq_bipyr",
                                       "pent_bipyr", "hex_bipyr", "oct_max",
                                       "sq_plan_max", "hex_plan_max"]:
                                if thetak < self._paras[i]['min_SPP']:
                                    tmp = self._paras[i]['IGW_EP'] * (
                                            thetak * ipi - 0.5) if t != "hex_plan_max" else \
                                            self._paras[i]['IGW_TA'] * (
                                            fabs(thetak * ipi - 0.5) - self._paras[i]['TA'])
                                    qsptheta[i][j][kc] += exp(
                                        -0.5 * tmp * tmp)
                                    norms[i][j][kc] += 1
                            elif t in ["pent_plan", "pent_plan_max"]:
                                tmp = 0.4 if thetak <= self._paras[i]['TA'] * pi \
                                        else 0.8
                                tmp2 = self._paras[i]['IGW_TA'] * (
                                        thetak * ipi - tmp)
                                gaussthetak[i] = exp(-0.5 * tmp2 * tmp2)
                                if t == "pent_plan_max":
                                    qsptheta[i][j][kc] += gaussthetak[i]
                                    norms[i][j][kc] += 1
                            elif t == "bcc" and j < k:
                                if thetak >= self._paras[i]['min_SPP']:
                                    tmp = self._paras[i]['IGW_SPP'] * (
                                            thetak * ipi - 1.0)
                                    qsptheta[i][j][kc] += (self._paras[i]['w_SPP'] *
                                                       exp(-0.5 * tmp * tmp))
                                    norms[i][j][kc] += self._paras[i]['w_SPP']

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

                                # South pole contributions of m.
                                if t in ["tri_bipyr", "sq_bipyr", "pent_bipyr",
                                         "hex_bipyr", "oct_max", "sq_plan_max",
                                         "hex_plan_max", "see_saw_rect"]:
                                    if thetam >= self._paras[i]['min_SPP']:
                                        tmp = self._paras[i]['IGW_SPP'] * (
                                                thetam * ipi - 1.0)
                                        qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                        norms[i][j][kc] += 1

                                # Contributions of j-i-m angle and
                                # angles between plane j-i-k and i-m vector.
                                if not flag_xaxis and not flag_xtwoaxis:
                                    for i, t in enumerate(self._types):
                                        if t in ["tri_plan", "tri_plan_max", \
                                                 "tet", "tet_max"]:
                                            tmp = self._paras[i]['IGW_TA'] * (
                                                thetam * ipi -
                                                self._paras[i]['TA'])
                                            tmp2 = cos(
                                                self._paras[i]['fac_AA'] *
                                                phi) ** self._paras[i][
                                                'exp_cos_AA']
                                            tmp3 = 1 if t in ["tri_plan_max", "tet_max"] \
                                                else gaussthetak[i]
                                            qsptheta[i][j][kc] += tmp3 * exp(
                                                -0.5 * tmp * tmp) * tmp2
                                            norms[i][j][kc] += 1
                                        elif t in ["pent_plan", "pent_plan_max"]:
                                            tmp = 0.4 if thetam <= self._paras[i]['TA'] * pi \
                                                    else 0.8
                                            tmp2 = self._paras[i]['IGW_TA'] * (
                                                    thetam * ipi - tmp)
                                            tmp3 = cos(phi)
                                            tmp4 = 1 if t == "pent_plan_max" \
                                                else gaussthetak[i]
                                            qsptheta[i][j][kc] += tmp4 * exp(
                                                    -0.5 * tmp2 * tmp2) * tmp3 * tmp3
                                            norms[i][j][kc] += 1
                                        elif t in ["T", "tri_pyr", "sq_pyr",
                                                   "pent_pyr", "hex_pyr"]:
                                            tmp = cos(self._paras[i]['fac_AA'] *
                                                      phi) ** self._paras[i][
                                                      'exp_cos_AA']
                                            tmp3 = self._paras[i]['IGW_EP'] * (
                                                    thetam * ipi - 0.5)
                                            qsptheta[i][j][kc] += tmp * exp(
                                                -0.5 * tmp3 * tmp3)
                                            norms[i][j][kc] += 1
                                        elif t in ["sq_plan", "oct", "oct_legacy"]:
                                            if thetak < self._paras[i]['min_SPP'] and \
                                                    thetam < self._paras[i]['min_SPP']:
                                                tmp = cos(self._paras[i]['fac_AA'] *
                                                        phi) ** self._paras[i]['exp_cos_AA']
                                                tmp2 = self._paras[i]['IGW_EP'] * (
                                                        thetam * ipi - 0.5)
                                                qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp2 * tmp2)
                                                if t == "oct_legacy":
                                                    qsptheta[i][j][kc] -= tmp * self._paras[i][6] * self._paras[i][7]
                                                norms[i][j][kc] += 1
                                        elif t in ["tri_bipyr", "sq_bipyr", "pent_bipyr",
                                                   "hex_bipyr", "oct_max", "sq_plan_max",
                                                   "hex_plan_max"]:
                                            if thetam < self._paras[i]['min_SPP']:
                                                if thetak < self._paras[i]['min_SPP']:
                                                    tmp = cos(self._paras[i]['fac_AA'] *
                                                            phi) ** self._paras[i]['exp_cos_AA']
                                                    tmp2 = self._paras[i]['IGW_EP'] * (
                                                            thetam * ipi - 0.5) if t != "hex_plan_max" else \
                                                            self._paras[i]['IGW_TA'] * (
                                                            fabs(thetam * ipi - 0.5) - self._paras[i]['TA'])
                                                    qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1
                                        elif t == "bcc" and j < k:
                                            if thetak < self._paras[i]['min_SPP']:
                                                if thetak > piover2:
                                                    fac = 1.0
                                                else:
                                                    fac = -1.0
                                                tmp = (thetam - piover2) / asin(1/3)
                                                qsptheta[i][j][kc] += fac * cos(
                                                    3.0 * phi) * fac_bcc * \
                                                    tmp * exp(-0.5 * tmp * tmp)
                                                norms[i][j][kc] += 1
                                        elif t == "see_saw_rect":
                                            if thetam < self._paras[i]['min_SPP']:
                                                if thetak < self._paras[i]['min_SPP'] and phi < 0.75 * pi:
                                                    tmp = cos(self._paras[i]['fac_AA'] *
                                                            phi) ** self._paras[i]['exp_cos_AA']
                                                    tmp2 = self._paras[i]['IGW_EP'] * (
                                                            thetam * ipi - 0.5)
                                                    qsptheta[i][j][kc] += tmp * \
                                                            exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1.0
                                        elif t in ["cuboct", "cuboct_max"]:
                                            if thetam < self._paras[i]['min_SPP'] and \
                                                    thetak > self._paras[i][4] and \
                                                    thetak < self._paras[i][2]:
                                                if thetam > self._paras[i][4] and \
                                                        thetam < self._paras[i][2]:
                                                    tmp = cos(phi)
                                                    tmp2 = self._paras[i][5] * (thetam * ipi - 0.5)
                                                    qsptheta[i][j][kc] += tmp * tmp * exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1.0
                                                elif thetam < self._paras[i][4]:
                                                    tmp = 0.0556 * (cos(
                                                            phi - 0.5 * pi) - 0.81649658)
                                                    tmp2 = self._paras[i][6] * (
                                                            thetam * ipi - onethird)
                                                    qsptheta[i][j][kc] += exp(
                                                        -0.5 * tmp * tmp) * \
                                                        exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1.0
                                                elif thetam > self._paras[i][2]:
                                                    tmp = 0.0556 * (cos(phi - 0.5 * pi) - \
                                                            0.81649658)
                                                    tmp2 = self._paras[i][6] * (thetam * ipi - \
                                                            twothird)
                                                    qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp) * \
                                                            exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1.0
                        kc += 1

            # Normalize Peters-style OPs.
            for i, t in enumerate(self._types):
                #if t == "pent_plan":
                #    ops[i] = ops[i] / sum(norms[i]) \
                #        if sum(norms[i]) > 1.0e-12 else None
                if t in ["tri_plan", "tet", "bent", "sq_plan",
                           "oct", "oct_legacy", "cuboct", "pent_plan"]:
                    ops[i] = tmp_norm = 0.0
                    for j in range(nneigh):
                        ops[i] += sum(qsptheta[i][j])
                        tmp_norm += float(sum(norms[i][j]))
                    ops[i] = ops[i] / tmp_norm if tmp_norm > 1.0e-12 else None
                elif t in ["T", "tri_pyr", "see_saw_rect", "sq_pyr", "tri_bipyr",
                        "sq_bipyr", "pent_pyr", "hex_pyr", "pent_bipyr",
                        "hex_bipyr", "oct_max", "tri_plan_max", "tet_max",
                        "sq_plan_max", "pent_plan_max", "cuboct_max", "hex_plan_max"]:
                    ops[i] = None
                    if nneigh > 1:
                        for j in range(nneigh):
                            for k in range(len(qsptheta[i][j])):
                                qsptheta[i][j][k] = qsptheta[i][j][k] / norms[i][j][k] \
                                    if norms[i][j][k] > 1.0e-12 else 0.0
                            ops[i] = max(qsptheta[i][j]) if j == 0 \
                                    else max(ops[i], max(qsptheta[i][j]))
                    #ops[i] = max(qsptheta[i]) if len(qsptheta[i]) > 0 else None
                elif t == "bcc":
                    ops[i] = 0.0
                    for j in range(nneigh):
                        ops[i] += sum(qsptheta[i][j])
                    ops[i] = ops[i] / float(0.5 * float(
                        nneigh * (6 + (nneigh - 2) * (nneigh - 3)))) \
                        if nneigh > 3 else None
                elif t == "sq_pyr_legacy":
                    if nneigh > 1:
                        dmean = np.mean(dist)
                        acc = 0.0
                        for d in dist:
                            tmp = self._paras[i][2] * (d - dmean)
                            acc = acc + exp(-0.5 * tmp * tmp)
                        for j in range(nneigh):
                            ops[i] = max(qsptheta[i][j]) if j == 0 \
                                    else max(ops[i], max(qsptheta[i][j]))
                        ops[i] = acc * ops[i] / float(nneigh)
                            #nneigh * (nneigh - 1))
                    else:
                        ops[i] = None

        # Then, deal with the new-style OPs that require vectors between
        # neighbors.
        if self._geomops2:
            # Compute all (unique) angles and sort the resulting list.
            aij = []
            for ir, r in enumerate(rijnorm):
                for j in range(ir + 1, len(rijnorm)):
                    aij.append(acos(max(-1.0, min(np.inner(r, rijnorm[j]), 1.0))))
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
                            a = 2.0 * asin(b / (2.0 * sqrt(h * h + (b / (
                                    2.0 * cos(3.0 * pi / 18.0))) ** 2.0)))
                            nmax = 3
                        elif t == "sq":
                            a = 2.0 * asin(
                                b / (2.0 * sqrt(h * h + dhalf * dhalf)))
                            nmax = 4
                        for j in range(min([nneigh, nmax])):
                            ops[i] = ops[i] * exp(-0.5 * ((
                                                                  aijs[j] - a) *
                                                          self._paras[i][
                                                              0]) ** 2)

        return ops
