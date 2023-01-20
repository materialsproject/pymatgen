# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures.
"""

from __future__ import annotations

import json
import math
import os
import warnings
from bisect import bisect_left
from collections import defaultdict, namedtuple
from copy import deepcopy
from functools import lru_cache
from math import acos, asin, atan2, cos, exp, fabs, pi, pow, sin, sqrt
from typing import Any, Literal, get_args

import numpy as np
from monty.dev import requires
from monty.serialization import loadfn
from ruamel.yaml import YAML
from scipy.spatial import Voronoi

from pymatgen.analysis.bond_valence import BV_PARAMS, BVAnalyzer
from pymatgen.analysis.graphs import MoleculeGraph, StructureGraph
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.core.composition import SpeciesLike
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.structure import IStructure, PeriodicNeighbor, Structure

try:
    from openbabel import openbabel
except Exception:
    openbabel = None

__author__ = "Shyue Ping Ong, Geoffroy Hautier, Sai Jayaraman, "
__author__ += "Nils E. R. Zimmermann, Bharat Medasani, Evan Spotte-Smith"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Nils E. R. Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__status__ = "Production"
__date__ = "August 17, 2017"

_directory = os.path.join(os.path.dirname(__file__))
yaml = YAML()

with open(os.path.join(_directory, "op_params.yaml")) as f:
    default_op_params = yaml.load(f)

with open(os.path.join(_directory, "cn_opt_params.yaml")) as f:
    cn_opt_params = yaml.load(f)

with open(os.path.join(_directory, "ionic_radii.json")) as fp:
    _ion_radii = json.load(fp)


class ValenceIonicRadiusEvaluator:
    """
    Computes site valences and ionic radii for a structure using bond valence
    analyzer
    """

    def __init__(self, structure):
        """
        Args:
            structure: pymatgen.core.structure.Structure
        """
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
        # print radii_dict
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

        def nearest_key(sorted_vals, skey):
            n = bisect_left(sorted_vals, skey)
            if n == len(sorted_vals):
                return sorted_vals[-1]
            if n == 0:
                return sorted_vals[0]
            before = sorted_vals[n - 1]
            after = sorted_vals[n]
            if after - skey < skey - before:
                return after
            return before

        for i, site in enumerate(self._structure.sites):
            if isinstance(site.specie, Element):
                radius = site.specie.atomic_radius
                # Handle elements with no atomic_radius
                # by using calculated values instead.
                if radius is None:
                    radius = site.specie.atomic_radius_calculated
                if radius is None:
                    raise ValueError(f"cannot assign radius to element {site.specie}")
                radii.append(radius)
                continue

            el = site.specie.symbol
            oxi_state = int(round(site.specie.oxi_state))
            coord_no = int(round(vnn.get_cn(self._structure, i)))
            try:
                tab_oxi_states = sorted(map(int, _ion_radii[el]))
                oxi_state = nearest_key(tab_oxi_states, oxi_state)
                radius = _ion_radii[el][str(oxi_state)][str(coord_no)]
            except KeyError:
                if vnn.get_cn(self._structure, i) - coord_no > 0:
                    new_coord_no = coord_no + 1
                else:
                    new_coord_no = coord_no - 1
                try:
                    radius = _ion_radii[el][str(oxi_state)][str(new_coord_no)]
                    coord_no = new_coord_no
                except Exception:
                    tab_coords = sorted(map(int, _ion_radii[el][str(oxi_state)]))
                    new_coord_no = nearest_key(tab_coords, coord_no)
                    i = 0
                    for val in tab_coords:
                        if val > coord_no:
                            break
                        i = i + 1
                    if i == len(tab_coords):
                        key = str(tab_coords[-1])
                        radius = _ion_radii[el][str(oxi_state)][key]
                    elif i == 0:
                        key = str(tab_coords[0])
                        radius = _ion_radii[el][str(oxi_state)][key]
                    else:
                        key = str(tab_coords[i - 1])
                        radius1 = _ion_radii[el][str(oxi_state)][key]
                        key = str(tab_coords[i])
                        radius2 = _ion_radii[el][str(oxi_state)][key]
                        radius = (radius1 + radius2) / 2

            # implement complex checks later
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
        except Exception:
            try:
                bv = BVAnalyzer(symm_tol=0.0)
                self._structure = bv.get_oxi_state_decorated_structure(self._structure)
                valences = bv.get_valences(self._structure)
            except Exception:
                valences = []
                for site in self._structure.sites:
                    if len(site.specie.common_oxidation_states) > 0:
                        valences.append(site.specie.common_oxidation_states[0])
                    # Handle noble gas species
                    # which have no entries in common_oxidation_states.
                    else:
                        valences.append(0)
                if sum(valences):
                    valences = [0] * self._structure.num_sites
                else:
                    self._structure.add_oxidation_state_by_site(valences)
                # raise

        # el = [site.specie.symbol for site in self._structure.sites]
        # el = [site.species_string for site in self._structure.sites]
        # el = [site.specie for site in self._structure.sites]
        # valence_dict = dict(zip(el, valences))
        # print valence_dict
        return valences


on_disorder_options = Literal["take_majority_strict", "take_majority_drop", "take_max_species", "error"]


def _handle_disorder(structure: Structure, on_disorder: on_disorder_options):
    """What to do in bonding and coordination number analysis if a site is disordered."""

    if all(site.is_ordered for site in structure):
        return structure

    if on_disorder == "error":
        raise ValueError(
            f"""Generating StructureGraphs for disordered Structures is unsupported. Pass on_disorder='take
            majority' | 'take_max_species' | 'error'. 'take_majority_strict' considers only the majority species from
            each site in the bonding algorithm and raises ValueError in case there is no majority (e.g. as in {{Fe:
            0.4, O: 0.4, C: 0.2}}) whereas 'take_majority_drop' just ignores the site altogether when computing bonds as
            if it didn't exist. 'take_max_species' extracts the first max species on each site (Fe in prev.
            example since Fe and O have equal occupancy and Fe comes first). 'error' raises an error in case
            of disordered structure. Offending {structure = }
        """
        )
    elif on_disorder.startswith("take_"):
        # disordered structures raise AttributeError when passed to NearNeighbors.get_cn()
        # or NearNeighbors.get_bonded_structure() (and probably others too, see GH-2070).
        # As a workaround, we create a new structure with majority species on each site.
        structure = structure.copy()  # make a copy so we don't mutate the original structure
        for idx, site in enumerate(structure):
            max_specie = max(site.species, key=site.species.get)  # type: ignore
            max_val = site.species[max_specie]
            if max_val <= 0.5:
                if on_disorder == "take_majority_strict":
                    raise ValueError(
                        f"Site {idx} has no majority species, the max species is {max_specie} with occupancy {max_val}"
                    )
                elif on_disorder == "take_majority_drop":
                    continue

            # this is the take_max_species case
            site.species = max_specie  # set site species in copied structure to max specie
    else:
        raise ValueError(f"Unexpected {on_disorder = }, should be one of {get_args(on_disorder_options)}")

    return structure


class NearNeighbors:
    """
    Base class to determine near neighbors that typically include nearest
    neighbors and others that are within some tolerable distance.
    """

    def __eq__(self, other: object) -> bool:
        if isinstance(other, type(self)):
            return self.__dict__ == other.__dict__
        return False

    def __hash__(self) -> int:
        return len(self.__dict__.items())

    @property
    def structures_allowed(self) -> bool:
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        raise NotImplementedError("structures_allowed is not defined!")

    @property
    def molecules_allowed(self) -> bool:
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        raise NotImplementedError("molecules_allowed is not defined!")

    @property
    def extend_structure_molecules(self) -> bool:
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        raise NotImplementedError("extend_structures_molecule is not defined!")

    def get_cn(
        self,
        structure: Structure,
        n: int,
        use_weights: bool = False,
        on_disorder: on_disorder_options = "take_majority_strict",
    ) -> float:
        """
        Get coordination number, CN, of site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine CN.
            use_weights (bool): flag indicating whether (True) to use weights for computing the coordination
                number or not (False, default: each coordinated site has equal weight).
            on_disorder ('take_majority_strict' | 'take_majority_drop' | 'take_max_species' | 'error'):
                What to do when encountering a disordered structure. 'error' will raise ValueError.
                'take_majority_strict' will use the majority specie on each site and raise
                ValueError if no majority exists. 'take_max_species' will use the first max specie
                on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, 'error' and 'take_majority_strict'
                will raise ValueError, while 'take_majority_drop' ignores this site altogether and
                'take_max_species' will use Fe as the site specie.
        Returns:
            cn (int or float): coordination number.
        """
        structure = _handle_disorder(structure, on_disorder)
        siw = self.get_nn_info(structure, n)
        return sum(e["weight"] for e in siw) if use_weights else len(siw)

    def get_cn_dict(self, structure: Structure, n: int, use_weights: bool = False):
        """
        Get coordination number, CN, of each element bonded to site with index n in structure

        Args:
            structure (Structure): input structure
            n (int): index of site for which to determine CN.
            use_weights (bool): flag indicating whether (True)
                to use weights for computing the coordination number
                or not (False, default: each coordinated site has equal
                weight).

        Returns:
            cn (dict): dictionary of CN of each element bonded to site
        """
        siw = self.get_nn_info(structure, n)

        cn_dict = {}
        for idx in siw:
            site_element = idx["site"].species_string
            if site_element not in cn_dict:
                if use_weights:
                    cn_dict[site_element] = idx["weight"]
                else:
                    cn_dict[site_element] = 1
            else:
                if use_weights:
                    cn_dict[site_element] += idx["weight"]
                else:
                    cn_dict[site_element] += 1
        return cn_dict

    def get_nn(self, structure: Structure, n: int):
        """
        Get near neighbors of site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site in structure for which to determine
                    neighbors.
        Returns:
            sites (list of Site objects): near neighbors.
        """
        return [e["site"] for e in self.get_nn_info(structure, n)]

    def get_weights_of_nn_sites(self, structure: Structure, n: int):
        """
        Get weight associated with each near neighbor of site with
        index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine the weights.
        Returns:
            weights (list of floats): near-neighbor weights.
        """
        return [e["weight"] for e in self.get_nn_info(structure, n)]

    def get_nn_images(self, structure: Structure, n: int):
        """
        Get image location of all near neighbors of site with index n in
        structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine the image
                location of near neighbors.
        Returns:
            images (list of 3D integer array): image locations of
                near neighbors.
        """
        return [e["image"] for e in self.get_nn_info(structure, n)]

    def get_nn_info(self, structure: Structure, n: int) -> list[dict]:
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                information.

        Returns:
            siw (list[dict]): each dictionary provides information
                about a single near neighbor, where key 'site' gives access to the
                corresponding Site object, 'image' gives the image location, and
                'weight' provides the weight that a given near-neighbor site contributes
                to the coordination number (1 or smaller), 'site_index' gives index of
                the corresponding site in the original structure.
        """
        raise NotImplementedError("get_nn_info(structure, n) is not defined!")

    def get_all_nn_info(self, structure):
        """Get a listing of all neighbors for all sites in a structure

        Args:
            structure (Structure): Input structure
        Return:
            List of NN site information for each site in the structure. Each
                entry has the same format as `get_nn_info`
        """
        return [self.get_nn_info(structure, n) for n in range(len(structure))]

    def get_nn_shell_info(self, structure: Structure, site_idx, shell):
        """Get a certain nearest neighbor shell for a certain site.

        Determines all non-backtracking paths through the neighbor network
        computed by `get_nn_info`. The weight is determined by multiplying
        the weight of the neighbor at each hop through the network. For
        example, a 2nd-nearest-neighbor that has a weight of 1 from its
        1st-nearest-neighbor and weight 0.5 from the original site will
        be assigned a weight of 0.5.

        As this calculation may involve computing the nearest neighbors of
        atoms multiple times, the calculation starts by computing all of the
        neighbor info and then calling `_get_nn_shell_info`. If you are likely
        to call this method for more than one site, consider calling `get_all_nn`
        first and then calling this protected method yourself.

        Args:
            structure (Structure): Input structure
            site_idx (int): index of site for which to determine neighbor
                information.
            shell (int): Which neighbor shell to retrieve (1 == 1st NN shell)
        Returns:
            list of dictionaries. Each entry in the list is information about
                a certain neighbor in the structure, in the same format as
                `get_nn_info`.
        """
        all_nn_info = self.get_all_nn_info(structure)
        sites = self._get_nn_shell_info(structure, all_nn_info, site_idx, shell)

        # Update the site positions
        #   Did not do this during NN options because that can be slower
        output = []
        for info in sites:
            orig_site = structure[info["site_index"]]
            info["site"] = PeriodicSite(
                orig_site.species,
                np.add(orig_site.frac_coords, info["image"]),
                structure.lattice,
                properties=orig_site.properties,
            )
            output.append(info)
        return output

    def _get_nn_shell_info(
        self,
        structure,
        all_nn_info,
        site_idx,
        shell,
        _previous_steps=frozenset(),
        _cur_image=(0, 0, 0),
    ):
        """Private method for computing the neighbor shell information

        Args:
            structure (Structure) - Structure being assessed
            all_nn_info ([[dict]]) - Results from `get_all_nn_info`
            site_idx (int) - index of site for which to determine neighbor
                information.
            shell (int) - Which neighbor shell to retrieve (1 == 1st NN shell)
            _previous_steps ({(site_idx, image}) - Internal use only: Set of
                sites that have already been traversed.
            _cur_image (tuple) - Internal use only Image coordinates of current atom
        Returns:
            list of dictionaries. Each entry in the list is information about
                a certain neighbor in the structure, in the same format as
                `get_nn_info`. Does not update the site positions
        """
        if shell <= 0:
            raise ValueError("Shell must be positive")

        # Append this site to the list of previously-visited sites
        _previous_steps = _previous_steps | {(site_idx, _cur_image)}

        # Get all the neighbors of this site
        possible_steps = list(all_nn_info[site_idx])
        for i, step in enumerate(possible_steps):
            # Update the image information
            #  Note: We do not update the site position yet, as making a
            #    PeriodicSite for each intermediate step is too costly
            step = dict(step)
            step["image"] = tuple(np.add(step["image"], _cur_image).tolist())
            possible_steps[i] = step

        # Get only the non-backtracking steps
        allowed_steps = [x for x in possible_steps if (x["site_index"], x["image"]) not in _previous_steps]

        # If we are the last step (i.e., shell == 1), done!
        if shell == 1:
            # No further work needed, just package these results
            return allowed_steps

        # If not, Get the N-1 NNs of these allowed steps
        terminal_neighbors = [
            self._get_nn_shell_info(
                structure,
                all_nn_info,
                x["site_index"],
                shell - 1,
                _previous_steps,
                x["image"],
            )
            for x in allowed_steps
        ]

        # Each allowed step results in many terminal neighbors
        #  And, different first steps might results in the same neighbor
        #  Now, we condense those neighbors into a single entry per neighbor
        all_sites = {}
        for first_site, term_sites in zip(allowed_steps, terminal_neighbors):
            for term_site in term_sites:
                key = (term_site["site_index"], tuple(term_site["image"]))

                # The weight for this site is equal to the weight of the
                #  first step multiplied by the weight of the terminal neighbor
                term_site["weight"] *= first_site["weight"]

                # Check if this site is already known
                value = all_sites.get(key)
                if value is not None:
                    # If so, add to its weight
                    value["weight"] += term_site["weight"]
                else:
                    # If not, prepare to add it
                    value = term_site
                all_sites[key] = value
        return list(all_sites.values())

    @staticmethod
    def _get_image(structure, site):
        """Private convenience method for get_nn_info,
        gives lattice image from provided PeriodicSite and Structure.

        Image is defined as displacement from original site in structure to a given site.
        i.e. if structure has a site at (-0.1, 1.0, 0.3), then (0.9, 0, 2.3) -> jimage = (1, -1, 2).
        Note that this method takes O(number of sites) due to searching an original site.

        Args:
            structure: Structure Object
            site: PeriodicSite Object

        Returns:
            image: ((int)*3) Lattice image
        """
        original_site = structure[NearNeighbors._get_original_site(structure, site)]
        image = np.around(np.subtract(site.frac_coords, original_site.frac_coords))
        image = tuple(image.astype(int))
        return image

    @staticmethod
    def _get_original_site(structure, site):
        """Private convenience method for get_nn_info,
        gives original site index from ProvidedPeriodicSite."""
        if isinstance(structure, (IStructure, Structure)):
            for i, s in enumerate(structure):
                if site.is_periodic_image(s):
                    return i
        else:
            for i, s in enumerate(structure):
                if site == s:
                    return i
        raise Exception("Site not found!")

    def get_bonded_structure(
        self,
        structure: Structure,
        decorate: bool = False,
        weights: bool = True,
        edge_properties: bool = False,
        on_disorder: on_disorder_options = "take_majority_strict",
    ) -> StructureGraph | MoleculeGraph:
        """
        Obtain a StructureGraph object using this NearNeighbor
        class. Requires the optional dependency networkx
        (pip install networkx).

        Args:
            structure: Structure object.
            decorate (bool): whether to annotate site properties with order parameters using neighbors
                determined by this NearNeighbor class
            weights (bool): whether to include edge weights from NearNeighbor class in StructureGraph
            edge_properties (bool) whether to include further edge properties from NearNeighbor class in StructureGraph
            on_disorder ('take_majority_strict' | 'take_majority_drop' | 'take_max_species' | 'error'):
                What to do when encountering a disordered structure. 'error' will raise ValueError.
                'take_majority_strict' will use the majority specie on each site and raise
                ValueError if no majority exists. 'take_max_species' will use the first max specie
                on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, 'error' and 'take_majority_strict'
                will raise ValueError, while 'take_majority_drop' ignores this site altogether and
                'take_max_species' will use Fe as the site specie.

        Returns: a pymatgen.analysis.graphs.StructureGraph object
        """
        structure = _handle_disorder(structure, on_disorder)

        if decorate:
            # Decorate all sites in the underlying structure
            # with site properties that provides information on the
            # coordination number and coordination pattern based
            # on the (current) structure of this graph.
            order_parameters = [self.get_local_order_parameters(structure, n) for n in range(len(structure))]
            structure.add_site_property("order_parameters", order_parameters)

        sg = StructureGraph.with_local_env_strategy(structure, self, weights=weights, edge_properties=edge_properties)

        # sets the attributes
        sg.set_node_attributes()
        return sg

    def get_local_order_parameters(self, structure: Structure, n: int):
        """
        Calculate those local structure order parameters for
        the given site whose ideal CN corresponds to the
        underlying motif (e.g., CN=4, then calculate the
        square planar, tetrahedral, see-saw-like,
        rectangular see-saw-like order parameters).

        Args:
            structure: Structure object
            n (int): site index.

        Returns (dict[str, float]):
            A dict of order parameters (values) and the
            underlying motif type (keys; for example, tetrahedral).
        """
        # code from @nisse3000, moved here from graphs to avoid circular
        # import, also makes sense to have this as a general NN method
        cn = self.get_cn(structure, n)
        int_cn = [int(k_cn) for k_cn in cn_opt_params]
        if cn in int_cn:
            names = list(cn_opt_params[cn])
            types = []
            params = []
            for name in names:
                types.append(cn_opt_params[cn][name][0])
                tmp = cn_opt_params[cn][name][1] if len(cn_opt_params[cn][name]) > 1 else None
                params.append(tmp)
            lostops = LocalStructOrderParams(types, parameters=params)
            sites = [structure[n]] + self.get_nn(structure, n)
            lostop_vals = lostops.get_order_parameters(sites, 0, indices_neighs=list(range(1, cn + 1)))  # type: ignore
            d = {}
            for i, lostop in enumerate(lostop_vals):
                d[names[i]] = lostop
            return d
        return None


class VoronoiNN(NearNeighbors):
    """
    Uses a Voronoi algorithm to determine near neighbors for each site in a
    structure.
    """

    def __init__(
        self,
        tol=0,
        targets=None,
        cutoff=13.0,
        allow_pathological=False,
        weight="solid_angle",
        extra_nn_info=True,
        compute_adj_neighbors=True,
    ):
        """
        Args:
            tol (float): tolerance parameter for near-neighbor finding. Faces that are
                smaller than `tol` fraction of the largest face are not included in the
                tessellation. (default: 0).
            targets (Element or list of Elements): target element(s).
            cutoff (float): cutoff radius in Angstrom to look for near-neighbor
                atoms. Defaults to 13.0.
            allow_pathological (bool): whether to allow infinite vertices in
                determination of Voronoi coordination.
            weight (string) - Statistic used to weigh neighbors (see the statistics
                available in get_voronoi_polyhedra)
            extra_nn_info (bool) - Add all polyhedron info to `get_nn_info`
            compute_adj_neighbors (bool) - Whether to compute which neighbors are
                adjacent. Turn off for faster performance
        """
        super().__init__()
        self.tol = tol
        self.cutoff = cutoff
        self.allow_pathological = allow_pathological
        self.targets = targets
        self.weight = weight
        self.extra_nn_info = extra_nn_info
        self.compute_adj_neighbors = compute_adj_neighbors

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_voronoi_polyhedra(self, structure: Structure, n: int):
        """
        Gives a weighted polyhedra around a site.

        See ref: A Proposed Rigorous Definition of Coordination Number,
        M. O'Keeffe, Acta Cryst. (1979). A35, 772-775

        Args:
            structure (Structure): structure for which to evaluate the
                coordination environment.
            n (int): site index.

        Returns:
            A dict of sites sharing a common Voronoi facet with the site
            n mapped to a directory containing statistics about the facet:
                - solid_angle - Solid angle subtended by face
                - angle_normalized - Solid angle normalized such that the
                    faces with the largest
                - area - Area of the facet
                - face_dist - Distance between site n and the facet
                - volume - Volume of Voronoi cell for this face
                - n_verts - Number of vertices on the facet
        """
        # Assemble the list of neighbors used in the tessellation
        #   Gets all atoms within a certain radius
        if self.targets is None:
            targets = structure.composition.elements
        else:
            targets = self.targets
        center = structure[n]

        cutoff = self.cutoff

        # max cutoff is the longest diagonal of the cell + room for noise
        corners = [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        d_corners = [np.linalg.norm(structure.lattice.get_cartesian_coords(c)) for c in corners]
        max_cutoff = max(d_corners) + 0.01

        while True:
            try:
                neighbors = structure.get_sites_in_sphere(center.coords, cutoff)
                neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]

                # Run the Voronoi tessellation
                qvoronoi_input = [s.coords for s in neighbors]

                voro = Voronoi(qvoronoi_input)  # can give seg fault if cutoff is too small

                # Extract data about the site in question
                cell_info = self._extract_cell_info(0, neighbors, targets, voro, self.compute_adj_neighbors)
                break

            except RuntimeError as e:
                if cutoff >= max_cutoff:
                    if e.args and "vertex" in e.args[0]:
                        # pass through the error raised by _extract_cell_info
                        raise e
                    raise RuntimeError("Error in Voronoi neighbor finding; max cutoff exceeded")
                cutoff = min(cutoff * 2, max_cutoff + 0.001)
        return cell_info

    def get_all_voronoi_polyhedra(self, structure):
        """Get the Voronoi polyhedra for all site in a simulation cell

        Args:
            structure (Structure): Structure to be evaluated
        Returns:
            A dict of sites sharing a common Voronoi facet with the site
            n mapped to a directory containing statistics about the facet:
                - solid_angle - Solid angle subtended by face
                - angle_normalized - Solid angle normalized such that the
                    faces with the largest
                - area - Area of the facet
                - face_dist - Distance between site n and the facet
                - volume - Volume of Voronoi cell for this face
                - n_verts - Number of vertices on the facet
        """
        # Special case: For atoms with 1 site, the atom in the root image is not
        # included in the get_all_neighbors output. Rather than creating logic to add
        # that atom to the neighbor list, which requires detecting whether it will be
        # translated to reside within the unit cell before neighbor detection, it is
        # less complex to just call the one-by-one operation
        if len(structure) == 1:
            return [self.get_voronoi_polyhedra(structure, 0)]

        # Assemble the list of neighbors used in the tessellation
        if self.targets is None:
            targets = structure.composition.elements
        else:
            targets = self.targets

        # Initialize the list of sites with the atoms in the origin unit cell
        # The `get_all_neighbors` function returns neighbors for each site's image in
        # the original unit cell. We start off with these central atoms to ensure they
        # are included in the tessellation

        sites = [x.to_unit_cell() for x in structure]
        indices = [(i, 0, 0, 0) for i, _ in enumerate(structure)]

        # Get all neighbors within a certain cutoff
        #   Record both the list of these neighbors, and the site indices
        all_neighs = structure.get_all_neighbors(self.cutoff, include_index=True, include_image=True)
        for neighs in all_neighs:
            sites.extend([x[0] for x in neighs])
            indices.extend([(x[2],) + x[3] for x in neighs])

        # Get the non-duplicates (using the site indices for numerical stability)
        indices = np.array(indices, dtype=int)
        indices, uniq_inds = np.unique(indices, return_index=True, axis=0)
        sites = [sites[i] for i in uniq_inds]

        # Sort array such that atoms in the root image are first
        #   Exploit the fact that the array is sorted by the unique operation such that
        #   the images associated with atom 0 are first, followed by atom 1, etc.
        (root_images,) = np.nonzero(np.abs(indices[:, 1:]).max(axis=1) == 0)

        del indices  # Save memory (tessellations can be costly)

        # Run the tessellation
        qvoronoi_input = [s.coords for s in sites]
        voro = Voronoi(qvoronoi_input)

        # Get the information for each neighbor
        return [
            self._extract_cell_info(i, sites, targets, voro, self.compute_adj_neighbors) for i in root_images.tolist()
        ]

    def _extract_cell_info(self, site_idx, sites, targets, voro, compute_adj_neighbors=False):
        """Get the information about a certain atom from the results of a tessellation

        Args:
            site_idx (int) - Index of the atom in question
            sites ([Site]) - List of all sites in the tessellation
            targets ([Element]) - Target elements
            voro - Output of qvoronoi
            compute_adj_neighbors (boolean) - Whether to compute which neighbors are adjacent
        Returns:
            A dict of sites sharing a common Voronoi facet. Key is facet id
             (not useful) and values are dictionaries containing statistics
             about the facet:
                - site: Pymatgen site
                - solid_angle - Solid angle subtended by face
                - angle_normalized - Solid angle normalized such that the
                    faces with the largest
                - area - Area of the facet
                - face_dist - Distance between site n and the facet
                - volume - Volume of Voronoi cell for this face
                - n_verts - Number of vertices on the facet
                - adj_neighbors - Facet id's for the adjacent neighbors
        """
        # Get the coordinates of every vertex
        all_vertices = voro.vertices

        # Get the coordinates of the central site
        center_coords = sites[site_idx].coords

        # Iterate through all the faces in the tessellation
        results = {}
        for nn, vind in voro.ridge_dict.items():
            # Get only those that include the site in question
            if site_idx in nn:
                other_site = nn[0] if nn[1] == site_idx else nn[1]
                if -1 in vind:
                    # -1 indices correspond to the Voronoi cell
                    #  missing a face
                    if self.allow_pathological:
                        continue

                    raise RuntimeError("This structure is pathological, infinite vertex in the Voronoi construction")

                # Get the solid angle of the face
                facets = [all_vertices[i] for i in vind]
                angle = solid_angle(center_coords, facets)

                # Compute the volume of associated with this face
                volume = 0
                # qvoronoi returns vertices in CCW order, so I can break
                # the face up in to segments (0,1,2), (0,2,3), ... to compute
                # its area where each number is a vertex size
                for j, k in zip(vind[1:], vind[2:]):
                    volume += vol_tetra(
                        center_coords,
                        all_vertices[vind[0]],
                        all_vertices[j],
                        all_vertices[k],
                    )

                # Compute the distance of the site to the face
                face_dist = np.linalg.norm(center_coords - sites[other_site].coords) / 2

                # Compute the area of the face (knowing V=Ad/3)
                face_area = 3 * volume / face_dist

                # Compute the normal of the facet
                normal = np.subtract(sites[other_site].coords, center_coords)
                normal /= np.linalg.norm(normal)

                # Store by face index
                results[other_site] = {
                    "site": sites[other_site],
                    "normal": normal,
                    "solid_angle": angle,
                    "volume": volume,
                    "face_dist": face_dist,
                    "area": face_area,
                    "n_verts": len(vind),
                }

                # If we are computing which neighbors are adjacent, store the vertices
                if compute_adj_neighbors:
                    results[other_site]["verts"] = vind

        # all sites should have at least two connected ridges in periodic system
        if len(results) == 0:
            raise ValueError("No Voronoi neighbors found for site - try increasing cutoff")

        # Get only target elements
        result_weighted = {}
        for nn_index, nn_stats in results.items():
            # Check if this is a target site
            nn = nn_stats["site"]
            if nn.is_ordered:
                if nn.specie in targets:
                    result_weighted[nn_index] = nn_stats
            else:  # if nn site is disordered
                for disordered_sp in nn.species:
                    if disordered_sp in targets:
                        result_weighted[nn_index] = nn_stats

        # If desired, determine which neighbors are adjacent
        if compute_adj_neighbors:
            # Initialize storage for the adjacent neighbors
            adj_neighbors = {i: [] for i in result_weighted}

            # Find the neighbors that are adjacent by finding those
            #  that contain exactly two vertices
            for a_ind, a_nn_info in result_weighted.items():
                # Get the indices for this site
                a_verts = set(a_nn_info["verts"])

                # Loop over all neighbors that have an index lower that this one
                # The goal here is to exploit the fact that neighbor adjacency is
                # symmetric (if A is adj to B, B is adj to A)
                for b_ind, b_nninfo in result_weighted.items():
                    if b_ind > a_ind:
                        continue
                    if len(a_verts.intersection(b_nninfo["verts"])) == 2:
                        adj_neighbors[a_ind].append(b_ind)
                        adj_neighbors[b_ind].append(a_ind)

            # Store the results in the nn_info
            for key, neighbors in adj_neighbors.items():
                result_weighted[key]["adj_neighbors"] = neighbors

        return result_weighted

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure
        using Voronoi decomposition.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        # Run the tessellation
        nns = self.get_voronoi_polyhedra(structure, n)

        # Extract the NN info
        return self._extract_nn_info(structure, nns)

    def get_all_nn_info(self, structure):
        """
        Args:
            structure (Structure): input structure.

        Returns:
            All nn info for all sites.
        """
        all_voro_cells = self.get_all_voronoi_polyhedra(structure)
        return [self._extract_nn_info(structure, cell) for cell in all_voro_cells]

    def _extract_nn_info(self, structure: Structure, nns):
        """Given Voronoi NNs, extract the NN info in the form needed by NearestNeighbors

        Args:
            structure (Structure): Structure being evaluated
            nns ([dicts]): Nearest neighbor information for a structure
        Returns:
            (list of tuples (Site, array, float)): See nn_info
        """
        # Get the target information
        if self.targets is None:
            targets = structure.composition.elements
        else:
            targets = self.targets

        # Extract the NN info
        siw = []
        max_weight = max(nn[self.weight] for nn in nns.values())
        for nstats in nns.values():
            site = nstats["site"]
            if nstats[self.weight] > self.tol * max_weight and _is_in_targets(site, targets):
                nn_info = {
                    "site": site,
                    "image": self._get_image(structure, site),
                    "weight": nstats[self.weight] / max_weight,
                    "site_index": self._get_original_site(structure, site),
                }

                if self.extra_nn_info:
                    # Add all the information about the site
                    poly_info = nstats
                    del poly_info["site"]
                    nn_info["poly_info"] = poly_info
                siw.append(nn_info)
        return siw


class IsayevNN(VoronoiNN):
    """
    Uses the algorithm defined in 10.1038/ncomms15679

    Sites are considered neighbors if (i) they share a Voronoi facet and (ii) the
    bond distance is less than the sum of the Cordero covalent radii + 0.25 Å.
    """

    def __init__(
        self,
        tol: float = 0.25,
        targets: Element | list[Element] | None = None,
        cutoff: float = 13.0,
        allow_pathological: bool = False,
        extra_nn_info: bool = True,
        compute_adj_neighbors: bool = True,
    ):
        """
        Args:
            tol: Tolerance in Å for bond distances that are considered coordinated.
            targets: Target element(s).
            cutoff: Cutoff radius in Angstrom to look for near-neighbor atoms.
            allow_pathological: Whether to allow infinite vertices in Voronoi
                coordination.
            extra_nn_info: Add all polyhedron info to `get_nn_info`.
            compute_adj_neighbors: Whether to compute which neighbors are adjacent. Turn
                off for faster performance.
        """
        super().__init__()
        self.tol = tol
        self.cutoff = cutoff
        self.allow_pathological = allow_pathological
        self.targets = targets
        self.extra_nn_info = extra_nn_info
        self.compute_adj_neighbors = compute_adj_neighbors

    def get_nn_info(self, structure: Structure, n: int) -> list[dict[str, Any]]:
        """
        Get all near-neighbor site information.

        Gets the associated image locations and weights of the site with index n
        in structure using Voronoi decomposition and distance cutoff.

        Args:
            structure: Input structure.
            n: Index of site for which to determine near-neighbor sites.

        Returns:
            List of dicts containing the near-neighbor information. Each dict has the
            keys:

            - "site": The near-neighbor site.
            - "image": The periodic image of the near-neighbor site.
            - "weight": The face weight of the Voronoi decomposition.
            - "site_index": The index of the near-neighbor site in the original
              structure.
        """
        nns = self.get_voronoi_polyhedra(structure, n)
        return self._filter_nns(structure, n, nns)

    def get_all_nn_info(self, structure: Structure) -> list[list[dict[str, Any]]]:
        """
        Args:
            structure (Structure): input structure.

        Returns:
            List of near neighbor information for each site. See get_nn_info for the
            format of the data for each site.
        """
        all_nns = self.get_all_voronoi_polyhedra(structure)
        return [self._filter_nns(structure, n, nns) for n, nns in enumerate(all_nns)]

    def _filter_nns(self, structure: Structure, n: int, nns: dict[str, Any]) -> list[dict[str, Any]]:
        """Extract and filter the NN info into the format needed by NearestNeighbors.

        Args:
            structure: The structure.
            n: The central site index.
            nns: Nearest neighbor information for the structure.

        Returns:
            See get_nn_info for the format of the returned data.
        """
        # Get the target information
        if self.targets is None:
            targets = structure.composition.elements
        else:
            targets = self.targets

        site = structure[n]

        # Extract the NN info
        siw = []
        max_weight = max(nn["area"] for nn in nns.values())
        for nstats in nns.values():
            nn = nstats.pop("site")

            # use the Cordero radius if it is available, otherwise the atomic radius
            cov_distance = _get_default_radius(site) + _get_default_radius(nn)
            nn_distance = np.linalg.norm(site.coords - nn.coords)

            # by default VoronoiNN only returns neighbors which share a Voronoi facet
            # therefore we don't need do to additional filtering based on the weight
            if _is_in_targets(nn, targets) and nn_distance <= cov_distance + self.tol:
                nn_info = {
                    "site": nn,
                    "image": self._get_image(structure, nn),
                    "weight": nstats["area"] / max_weight,
                    "site_index": self._get_original_site(structure, nn),
                }

                if self.extra_nn_info:
                    nn_info["poly_info"] = nstats
                siw.append(nn_info)
        return siw


def _is_in_targets(site, targets):
    """
    Test whether a site contains elements in the target list

    Args:
        site (Site): Site to assess
        targets ([Element]) List of elements
    Returns:
         (boolean) Whether this site contains a certain list of elements
    """
    elems = _get_elements(site)
    for elem in elems:
        if elem not in targets:
            return False
    return True


def _get_elements(site):
    """
    Get the list of elements for a Site

    Args:
         site (Site): Site to assess
    Returns:
        [Element]: List of elements
    """
    try:
        if isinstance(site.specie, Element):
            return [site.specie]
        return [Element(site.specie)]
    except Exception:
        return site.species.elements


class JmolNN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using an emulation
    of Jmol's default autoBond() algorithm. This version of the algorithm
    does not take into account any information regarding known charge
    states.
    """

    def __init__(
        self,
        tol: float = 0.45,
        min_bond_distance: float = 0.4,
        el_radius_updates: dict[SpeciesLike, float] | None = None,
    ):
        """
        Args:
            tol (float): tolerance parameter for bond determination
                Defaults to 0.56.
            min_bond_distance (float): minimum bond distance for consideration
                Defaults to 0.4.
            el_radius_updates: (dict) symbol->float to override default atomic
                radii table values
        """
        self.tol = tol
        self.min_bond_distance = min_bond_distance

        # Load elemental radii table
        bonds_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bonds_jmol_ob.yaml")
        with open(bonds_file) as f:
            yaml = YAML()
            self.el_radius = yaml.load(f)

        # Update any user preference elemental radii
        if el_radius_updates:
            self.el_radius.update(el_radius_updates)

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return True

    def get_max_bond_distance(self, el1_sym, el2_sym):
        """
        Use Jmol algorithm to determine bond length from atomic parameters
        Args:
            el1_sym: (str) symbol of atom 1
            el2_sym: (str) symbol of atom 2

        Returns: (float) max bond length
        """
        return sqrt((self.el_radius[el1_sym] + self.el_radius[el2_sym] + self.tol) ** 2)

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the bond identification
        algorithm underlying Jmol.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near
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
            bonds[site.specie, el] = self.get_max_bond_distance(site.specie.symbol, el.symbol)

        # Search for neighbors up to max bond length + tolerance
        max_rad = max(bonds.values()) + self.tol
        min_rad = min(bonds.values())

        siw = []
        for nn in structure.get_neighbors(site, max_rad):
            dist = nn.nn_distance
            # Confirm neighbor based on bond length specific to atom pair
            if dist <= (bonds[(site.specie, nn.specie)]) and (nn.nn_distance > self.min_bond_distance):
                weight = min_rad / dist
                siw.append(
                    {
                        "site": nn,
                        "image": self._get_image(structure, nn),
                        "weight": weight,
                        "site_index": self._get_original_site(structure, nn),
                    }
                )
        return siw


class MinimumDistanceNN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    nearest neighbor(s) at distance, d_min, plus all neighbors
    within a distance (1 + tol) * d_min, where tol is a
    (relative) distance tolerance parameter.
    """

    def __init__(self, tol: float = 0.1, cutoff=10.0, get_all_sites=False):
        """
        Args:
            tol (float): tolerance parameter for neighbor identification
                (default: 0.1).
            cutoff (float): cutoff radius in Angstrom to look for trial
                near-neighbor sites (default: 10.0).
            get_all_sites (bool): If this is set to True then the neighbor
                sites are only determined by the cutoff radius, tol is ignored
        """
        self.tol = tol
        self.cutoff = cutoff
        self.get_all_sites = get_all_sites

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return True

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest neighbor
        distance-based method.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near
                neighbors.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """
        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        is_periodic = isinstance(structure, (Structure, IStructure))
        siw = []
        if self.get_all_sites:
            for nn in neighs_dists:
                w = nn.nn_distance
                siw.append(
                    {
                        "site": nn,
                        "image": self._get_image(structure, nn) if is_periodic else None,
                        "weight": w,
                        "site_index": self._get_original_site(structure, nn),
                    }
                )
        else:
            min_dist = min(nn.nn_distance for nn in neighs_dists)
            for nn in neighs_dists:
                dist = nn.nn_distance
                if dist < (1.0 + self.tol) * min_dist:
                    w = min_dist / dist
                    siw.append(
                        {
                            "site": nn,
                            "image": self._get_image(structure, nn) if is_periodic else None,
                            "weight": w,
                            "site_index": self._get_original_site(structure, nn),
                        }
                    )
        return siw


class OpenBabelNN(NearNeighbors):
    """
    Determine near-neighbor sites and bond orders using OpenBabel API.

    NOTE: This strategy is only appropriate for molecules, and not for
    structures.
    """

    @requires(
        openbabel,
        "BabelMolAdaptor requires openbabel to be installed with "
        "Python bindings. Please get it at http://openbabel.org "
        "(version >=3.0.0).",
    )
    def __init__(self, order=True):
        """
        Args:
            order (bool): True if bond order should be returned as a weight, False
            if bond length should be used as a weight.
        """
        self.order = order

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return False

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return False

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites and weights (orders) of bonds for a given
        atom.

        Args:
            structure: Molecule object.
            n: index of site for which to determine near neighbors.

        Returns:
            (dict): representing a neighboring site and the type of
            bond present between site n and the neighboring site.
        """
        from pymatgen.io.babel import BabelMolAdaptor

        obmol = BabelMolAdaptor(structure).openbabel_mol

        siw = []

        # Get only the atom of interest
        site_atom = [
            a
            for i, a in enumerate(openbabel.OBMolAtomDFSIter(obmol))
            if [a.GetX(), a.GetY(), a.GetZ()] == list(structure[n].coords)
        ][0]

        for neighbor in openbabel.OBAtomAtomIter(site_atom):
            coords = [neighbor.GetX(), neighbor.GetY(), neighbor.GetZ()]
            site = [a for a in structure if list(a.coords) == coords][0]
            index = structure.index(site)

            bond = site_atom.GetBond(neighbor)

            if self.order:
                obmol.PerceiveBondOrders()
                weight = bond.GetBondOrder()
            else:
                weight = bond.GetLength()

            siw.append(
                {
                    "site": site,
                    "image": (0, 0, 0),
                    "weight": weight,
                    "site_index": index,
                }
            )

        return siw

    def get_bonded_structure(self, structure: Structure, decorate: bool = False) -> StructureGraph:  # type: ignore
        """
        Obtain a MoleculeGraph object using this NearNeighbor
        class. Requires the optional dependency networkx
        (pip install networkx).

        Args:
            structure: Molecule object.
            decorate (bool): whether to annotate site properties
            with order parameters using neighbors determined by
            this NearNeighbor class

        Returns: a pymatgen.analysis.graphs.MoleculeGraph object
        """
        if decorate:
            # Decorate all sites in the underlying structure
            # with site properties that provides information on the
            # coordination number and coordination pattern based
            # on the (current) structure of this graph.
            order_parameters = [self.get_local_order_parameters(structure, n) for n in range(len(structure))]
            structure.add_site_property("order_parameters", order_parameters)

        mg = MoleculeGraph.with_local_env_strategy(structure, self)

        return mg

    def get_nn_shell_info(self, structure: Structure, site_idx, shell):
        """Get a certain nearest neighbor shell for a certain site.

        Determines all non-backtracking paths through the neighbor network
        computed by `get_nn_info`. The weight is determined by multiplying
        the weight of the neighbor at each hop through the network. For
        example, a 2nd-nearest-neighbor that has a weight of 1 from its
        1st-nearest-neighbor and weight 0.5 from the original site will
        be assigned a weight of 0.5.

        As this calculation may involve computing the nearest neighbors of
        atoms multiple times, the calculation starts by computing all of the
        neighbor info and then calling `_get_nn_shell_info`. If you are likely
        to call this method for more than one site, consider calling `get_all_nn`
        first and then calling this protected method yourself.

        Args:
            structure (Molecule): Input structure
            site_idx (int): index of site for which to determine neighbor
                information.
            shell (int): Which neighbor shell to retrieve (1 == 1st NN shell)
        Returns:
            list of dictionaries. Each entry in the list is information about
                a certain neighbor in the structure, in the same format as
                `get_nn_info`.
        """
        all_nn_info = self.get_all_nn_info(structure)
        sites = self._get_nn_shell_info(structure, all_nn_info, site_idx, shell)

        # Update the site positions
        #   Did not do this during NN options because that can be slower
        output = []
        for info in sites:
            orig_site = structure[info["site_index"]]
            info["site"] = Site(orig_site.species, orig_site._coords, properties=orig_site.properties)
            output.append(info)
        return output


class CovalentBondNN(NearNeighbors):
    """
    Determine near-neighbor sites and bond orders using built-in
    pymatgen.Molecule CovalentBond functionality.

    NOTE: This strategy is only appropriate for molecules, and not for
    structures.
    """

    def __init__(self, tol: float = 0.2, order=True):
        """
        Args:
            tol (float): Tolerance for covalent bond checking.
            order (bool): If True (default), this class will compute bond
                orders. If False, bond lengths will be computed
        """
        self.tol = tol
        self.order = order

        self.bonds = None

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return False

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return False

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites and weights (orders) of bonds for a given
        atom.

        :param structure: input Molecule.
        :param n: index of site for which to determine near neighbors.
        :return: [dict] representing a neighboring site and the type of
            bond present between site n and the neighboring site.
        """
        # This is unfortunately inefficient, but is the best way to fit the
        # current NearNeighbors scheme
        self.bonds = bonds = structure.get_covalent_bonds(tol=self.tol)

        siw = []

        for bond in bonds:
            capture_bond = False
            if bond.site1 == structure[n]:
                site = bond.site2
                capture_bond = True
            elif bond.site2 == structure[n]:
                site = bond.site1
                capture_bond = True

            if capture_bond:
                index = structure.index(site)
                if self.order:
                    weight = bond.get_bond_order()
                else:
                    weight = bond.length

                siw.append({"site": site, "image": (0, 0, 0), "weight": weight, "site_index": index})

        return siw

    def get_bonded_structure(self, structure: Structure, decorate: bool = False) -> MoleculeGraph:  # type: ignore
        """
        Obtain a MoleculeGraph object using this NearNeighbor
        class.

        Args:
            structure: Molecule object.
            decorate (bool): whether to annotate site properties
            with order parameters using neighbors determined by
            this NearNeighbor class

        Returns: a pymatgen.analysis.graphs.MoleculeGraph object
        """
        # requires optional dependency which is why it's not a top-level import
        from pymatgen.analysis.graphs import MoleculeGraph

        if decorate:
            # Decorate all sites in the underlying structure
            # with site properties that provides information on the
            # coordination number and coordination pattern based
            # on the (current) structure of this graph.
            order_parameters = [self.get_local_order_parameters(structure, n) for n in range(len(structure))]
            structure.add_site_property("order_parameters", order_parameters)

        mg = MoleculeGraph.with_local_env_strategy(structure, self)

        return mg

    def get_nn_shell_info(self, structure: Structure, site_idx, shell):
        """Get a certain nearest neighbor shell for a certain site.

        Determines all non-backtracking paths through the neighbor network
        computed by `get_nn_info`. The weight is determined by multiplying
        the weight of the neighbor at each hop through the network. For
        example, a 2nd-nearest-neighbor that has a weight of 1 from its
        1st-nearest-neighbor and weight 0.5 from the original site will
        be assigned a weight of 0.5.

        As this calculation may involve computing the nearest neighbors of
        atoms multiple times, the calculation starts by computing all of the
        neighbor info and then calling `_get_nn_shell_info`. If you are likely
        to call this method for more than one site, consider calling `get_all_nn`
        first and then calling this protected method yourself.

        Args:
            structure (Molecule): Input structure
            site_idx (int): index of site for which to determine neighbor
                information.
            shell (int): Which neighbor shell to retrieve (1 == 1st NN shell)
        Returns:
            list of dictionaries. Each entry in the list is information about
                a certain neighbor in the structure, in the same format as
                `get_nn_info`.
        """
        all_nn_info = self.get_all_nn_info(structure)
        sites = self._get_nn_shell_info(structure, all_nn_info, site_idx, shell)

        # Update the site positions
        #   Did not do this during NN options because that can be slower
        output = []
        for info in sites:
            orig_site = structure[info["site_index"]]
            info["site"] = Site(orig_site.species, orig_site._coords, properties=orig_site.properties)
            output.append(info)
        return output


class MinimumOKeeffeNN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    neighbor(s) at closest relative distance, d_min_OKeffee, plus some
    relative tolerance, where bond valence parameters from O'Keeffe's
    bond valence method (J. Am. Chem. Soc. 1991, 3226-3229) are used
    to calculate relative distances.
    """

    def __init__(self, tol: float = 0.1, cutoff=10.0):
        """
        Args:
            tol (float): tolerance parameter for neighbor identification
                (default: 0.1).
            cutoff (float): cutoff radius in Angstrom to look for trial
                near-neighbor sites (default: 10.0).
        """
        self.tol = tol
        self.cutoff = cutoff

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return True

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest relative
        neighbor distance-based method with O'Keeffe parameters.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near
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
        except Exception:
            eln = site.species_string

        reldists_neighs = []
        for nn in neighs_dists:
            neigh = nn
            dist = nn.nn_distance
            try:
                el2 = neigh.specie.element
            except Exception:
                el2 = neigh.species_string
            reldists_neighs.append([dist / get_okeeffe_distance_prediction(eln, el2), neigh])

        siw = []
        min_reldist = min(reldist for reldist, neigh in reldists_neighs)
        for reldist, s in reldists_neighs:
            if reldist < (1.0 + self.tol) * min_reldist:
                w = min_reldist / reldist
                siw.append(
                    {
                        "site": s,
                        "image": self._get_image(structure, s),
                        "weight": w,
                        "site_index": self._get_original_site(structure, s),
                    }
                )

        return siw


class MinimumVIRENN(NearNeighbors):
    """
    Determine near-neighbor sites and coordination number using the
    neighbor(s) at closest relative distance, d_min_VIRE, plus some
    relative tolerance, where atom radii from the
    ValenceIonicRadiusEvaluator (VIRE) are used
    to calculate relative distances.
    """

    def __init__(self, tol: float = 0.1, cutoff=10.0):
        """
        Args:
            tol (float): tolerance parameter for neighbor identification
                (default: 0.1).
            cutoff (float): cutoff radius in Angstrom to look for trial
                near-neighbor sites (default: 10.0).
        """
        self.tol = tol
        self.cutoff = cutoff

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest relative
        neighbor distance-based method with VIRE atomic/ionic radii.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near
                neighbors.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """
        vire = _get_vire(structure)
        site = vire.structure[n]
        neighs_dists = vire.structure.get_neighbors(site, self.cutoff)
        rn = vire.radii[vire.structure[n].species_string]

        reldists_neighs = []
        for nn in neighs_dists:
            reldists_neighs.append([nn.nn_distance / (vire.radii[nn.species_string] + rn), nn])

        siw = []
        min_reldist = min(reldist for reldist, neigh in reldists_neighs)
        for reldist, s in reldists_neighs:
            if reldist < (1.0 + self.tol) * min_reldist:
                w = min_reldist / reldist
                siw.append(
                    {
                        "site": s,
                        "image": self._get_image(vire.structure, s),
                        "weight": w,
                        "site_index": self._get_original_site(vire.structure, s),
                    }
                )

        return siw


def _get_vire(structure: Structure | IStructure):
    """Get the ValenceIonicRadiusEvaluator object for an structure taking
    advantage of caching.

    Args:
        structure: A structure.

    Returns:
        Output of `ValenceIonicRadiusEvaluator(structure)`
    """
    # pymatgen does not hash Structure objects, so we need
    # to cast from Structure to the immutable IStructure
    if isinstance(structure, Structure):
        structure = IStructure.from_sites(structure)

    return _get_vire_istructure(structure)


@lru_cache(maxsize=1)
def _get_vire_istructure(structure: IStructure):
    """Get the ValenceIonicRadiusEvaluator object for an immutable structure
    taking advantage of caching.

    Args:
        structure: A structure.

    Returns:
        Output of `ValenceIonicRadiusEvaluator(structure)`
    """
    return ValenceIonicRadiusEvaluator(structure)


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

    # Compute the displacement from the center
    r = [np.subtract(c, center) for c in coords]

    # Compute the magnitude of each vector
    r_norm = [np.linalg.norm(i) for i in r]

    # Compute the solid angle for each tetrahedron that makes up the facet
    #  Following: https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
    angle = 0
    for i in range(1, len(r) - 1):
        j = i + 1
        tp = np.abs(np.dot(r[0], np.cross(r[i], r[j])))
        de = (
            r_norm[0] * r_norm[i] * r_norm[j]
            + r_norm[j] * np.dot(r[0], r[i])
            + r_norm[i] * np.dot(r[0], r[j])
            + r_norm[0] * np.dot(r[i], r[j])
        )
        if de == 0:
            my_angle = 0.5 * pi if tp > 0 else -0.5 * pi
        else:
            my_angle = np.arctan(tp / de)
        angle += (my_angle if my_angle > 0 else my_angle + np.pi) * 2

    return angle


def vol_tetra(vt1, vt2, vt3, vt4):
    """
    Calculate the volume of a tetrahedron, given the four vertices of vt1,
    vt2, vt3 and vt4.
    Args:
        vt1 (array-like): coordinates of vertex 1.
        vt2 (array-like): coordinates of vertex 2.
        vt3 (array-like): coordinates of vertex 3.
        vt4 (array-like): coordinates of vertex 4.
    Returns:
        (float): volume of the tetrahedron.
    """
    vol_tetra = np.abs(np.dot((vt1 - vt4), np.cross((vt2 - vt4), (vt3 - vt4)))) / 6
    return vol_tetra


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
    if el not in list(BV_PARAMS):
        raise RuntimeError(
            "Could not find O'Keeffe parameters for element"
            f' {el_symbol!r} in "BV_PARAMS"dictionary'
            " provided by pymatgen"
        )

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

    r1 = el1_okeeffe_params["r"]
    r2 = el2_okeeffe_params["r"]
    c1 = el1_okeeffe_params["c"]
    c2 = el2_okeeffe_params["c"]

    return r1 + r2 - r1 * r2 * pow(sqrt(c1) - sqrt(c2), 2) / (c1 * r1 + c2 * r2)


def get_neighbors_of_site_with_index(struct, n, approach="min_dist", delta=0.1, cutoff=10.0):
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
        return MinimumDistanceNN(tol=delta, cutoff=cutoff).get_nn(struct, n)

    if approach == "voronoi":
        return VoronoiNN(tol=delta, cutoff=cutoff).get_nn(struct, n)

    if approach == "min_OKeeffe":
        return MinimumOKeeffeNN(tol=delta, cutoff=cutoff).get_nn(struct, n)

    if approach == "min_VIRE":
        return MinimumVIRENN(tol=delta, cutoff=cutoff).get_nn(struct, n)
    raise RuntimeError(f"unsupported neighbor-finding method ({approach}).")


def site_is_of_motif_type(struct, n, approach="min_dist", delta=0.1, cutoff=10.0, thresh=None):
    """
    Returns the motif type of the site with index n in structure struct;
    currently featuring "tetrahedral", "octahedral", "bcc", and "cp"
    (close-packed: fcc and hcp) as well as "square pyramidal" and
    "trigonal bipyramidal". If the site is not recognized,
    "unrecognized" is returned. If a site should be assigned to two
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
            "qtet": 0.5,
            "qoct": 0.5,
            "qbcc": 0.5,
            "q6": 0.4,
            "qtribipyr": 0.8,
            "qsqpyr": 0.8,
        }

    ops = LocalStructOrderParams(["cn", "tet", "oct", "bcc", "q6", "sq_pyr", "tri_bipyr"])

    neighs_cent = get_neighbors_of_site_with_index(struct, n, approach=approach, delta=delta, cutoff=cutoff)
    neighs_cent.append(struct.sites[n])
    opvals = ops.get_order_parameters(
        neighs_cent,
        len(neighs_cent) - 1,
        indices_neighs=list(range(len(neighs_cent) - 1)),
    )
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
    if cn == 12 and (
        opvals[4] > thresh["q6"] and opvals[1] < thresh["q6"] and opvals[2] < thresh["q6"] and opvals[3] < thresh["q6"]
    ):
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


class LocalStructOrderParams:
    """
    This class permits the calculation of various types of local
    structure order parameters.
    """

    __supported_types = (
        "cn",
        "sgl_bd",
        "bent",
        "tri_plan",
        "tri_plan_max",
        "reg_tri",
        "sq_plan",
        "sq_plan_max",
        "pent_plan",
        "pent_plan_max",
        "sq",
        "tet",
        "tet_max",
        "tri_pyr",
        "sq_pyr",
        "sq_pyr_legacy",
        "tri_bipyr",
        "sq_bipyr",
        "oct",
        "oct_legacy",
        "pent_pyr",
        "hex_pyr",
        "pent_bipyr",
        "hex_bipyr",
        "T",
        "cuboct",
        "cuboct_max",
        "see_saw_rect",
        "bcc",
        "q2",
        "q4",
        "q6",
        "oct_max",
        "hex_plan_max",
        "sq_face_cap_trig_pris",
    )

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
                the op_params.yaml file. With few exceptions, 9 different
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
            if t not in LocalStructOrderParams.__supported_types:
                raise ValueError("Unknown order parameter type (" + t + ")!")
        self._types = tuple(types)

        self._comp_azi = False
        self._params = []
        for i, t in enumerate(self._types):
            d = deepcopy(default_op_params[t]) if default_op_params[t] is not None else None
            if parameters is None:
                self._params.append(d)
            elif parameters[i] is None:
                self._params.append(d)
            else:
                self._params.append(deepcopy(parameters[i]))

        self._computerijs = self._computerjks = self._geomops = False
        self._geomops2 = self._boops = False
        self._max_trig_order = -1

        # Add here any additional flags to be used during calculation.
        if "sgl_bd" in self._types:
            self._computerijs = True
        if not set(self._types).isdisjoint(
            [
                "tet",
                "oct",
                "bcc",
                "sq_pyr",
                "sq_pyr_legacy",
                "tri_bipyr",
                "sq_bipyr",
                "oct_legacy",
                "tri_plan",
                "sq_plan",
                "pent_plan",
                "tri_pyr",
                "pent_pyr",
                "hex_pyr",
                "pent_bipyr",
                "hex_bipyr",
                "T",
                "cuboct",
                "oct_max",
                "tet_max",
                "tri_plan_max",
                "sq_plan_max",
                "pent_plan_max",
                "cuboct_max",
                "bent",
                "see_saw_rect",
                "hex_plan_max",
                "sq_face_cap_trig_pris",
            ]
        ):
            self._computerijs = self._geomops = True
        if "sq_face_cap_trig_pris" in self._types:
            self._comp_azi = True
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
        """
        Returns:
            int: the number of different order parameters that are targeted
                to be calculated.
        """
        return len(self._types)

    @property
    def last_nneigh(self):
        """
        Returns:
            int: the number of neighbors encountered during the most
                recent order parameter calculation. A value of -1 indicates
                that no such calculation has yet been performed for this
                instance.
        """
        return len(self._last_nneigh)

    def compute_trigonometric_terms(self, thetas, phis):
        """
        Computes trigonometric terms that are required to
        calculate bond orientational order parameters using
        internal variables.

        Args:
            thetas ([float]): polar angles of all neighbors in radians.
            phis ([float]): azimuth angles of all neighbors in radians.
                The list of
                azimuth angles of all neighbors in radians. The list of
                azimuth angles is expected to have the same size as the
                list of polar angles; otherwise, a ValueError is raised.
                Also, the two lists of angles have to be coherent in
                order. That is, it is expected that the order in the list
                of azimuth angles corresponds to a distinct sequence of
                neighbors. And, this sequence has to equal the sequence
                of neighbors in the list of polar angles.
        """
        if len(thetas) != len(phis):
            raise ValueError("List of polar and azimuthal angles have to be equal!")

        self._pow_sin_t.clear()
        self._pow_cos_t.clear()
        self._sin_n_p.clear()
        self._cos_n_p.clear()

        self._pow_sin_t[1] = [sin(float(t)) for t in thetas]
        self._pow_cos_t[1] = [cos(float(t)) for t in thetas]
        self._sin_n_p[1] = [sin(float(p)) for p in phis]
        self._cos_n_p[1] = [cos(float(p)) for p in phis]

        for i in range(2, self._max_trig_order + 1):
            self._pow_sin_t[i] = [e[0] * e[1] for e in zip(self._pow_sin_t[i - 1], self._pow_sin_t[1])]
            self._pow_cos_t[i] = [e[0] * e[1] for e in zip(self._pow_cos_t[i - 1], self._pow_cos_t[1])]
            self._sin_n_p[i] = [sin(float(i) * float(p)) for p in phis]
            self._cos_n_p[i] = [cos(float(i) * float(p)) for p in phis]

    def get_q2(self, thetas=None, phis=None):
        """
        Calculates the value of the bond orientational order parameter of
        weight l=2. If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh. Otherwise, it is expected that the
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

        sqrt_15_2pi = sqrt(15.0 / (2 * pi))
        sqrt_5_pi = sqrt(5.0 / pi)

        pre_y_2_2 = [0.25 * sqrt_15_2pi * val for val in self._pow_sin_t[2]]
        pre_y_2_1 = [0.5 * sqrt_15_2pi * val[0] * val[1] for val in zip(self._pow_sin_t[1], self._pow_cos_t[1])]

        acc = 0.0

        # Y_2_-2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_2_2[i] * self._sin_n_p[2][i]
        acc += real * real + imag * imag

        # Y_2_-1
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_2_1[i] * self._sin_n_p[1][i]
        acc += real * real + imag * imag

        # Y_2_0
        real = imag = 0.0
        for i in nnn_range:
            real += 0.25 * sqrt_5_pi * (3 * self._pow_cos_t[2][i] - 1.0)
        acc += real * real

        # Y_2_1
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_2_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_2_1[i] * self._sin_n_p[1][i]
        acc += real * real + imag * imag

        # Y_2_2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_2[i] * self._cos_n_p[2][i]
            imag += pre_y_2_2[i] * self._sin_n_p[2][i]
        acc += real * real + imag * imag

        q2 = sqrt(4 * pi * acc / (5 * float(nnn * nnn)))
        return q2

    def get_q4(self, thetas=None, phis=None):
        """
        Calculates the value of the bond orientational order parameter of
        weight l=4. If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh. Otherwise, it is expected that the
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
        sqrt_35_2pi = sqrt(35.0 / (2 * pi))
        sqrt_5_pi = sqrt(5.0 / pi)
        sqrt_5_2pi = sqrt(5.0 / (2 * pi))
        sqrt_1_pi = sqrt(1.0 / pi)

        pre_y_4_4 = [i16_3 * sqrt_35_2pi * val for val in self._pow_sin_t[4]]
        pre_y_4_3 = [i8_3 * sqrt_35_pi * val[0] * val[1] for val in zip(self._pow_sin_t[3], self._pow_cos_t[1])]
        pre_y_4_2 = [
            i8_3 * sqrt_5_2pi * val[0] * (7 * val[1] - 1.0) for val in zip(self._pow_sin_t[2], self._pow_cos_t[2])
        ]
        pre_y_4_1 = [
            i8_3 * sqrt_5_pi * val[0] * (7 * val[1] - 3 * val[2])
            for val in zip(self._pow_sin_t[1], self._pow_cos_t[3], self._pow_cos_t[1])
        ]

        acc = 0.0

        # Y_4_-4
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_4_4[i] * self._sin_n_p[4][i]
        acc += real * real + imag * imag

        # Y_4_-3
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_4_3[i] * self._sin_n_p[3][i]
        acc += real * real + imag * imag

        # Y_4_-2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_4_2[i] * self._sin_n_p[2][i]
        acc += real * real + imag * imag

        # Y_4_-1
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_4_1[i] * self._sin_n_p[1][i]
        acc += real * real + imag * imag

        # Y_4_0
        real = imag = 0.0
        for i in nnn_range:
            real += i16_3 * sqrt_1_pi * (35 * self._pow_cos_t[4][i] - 30 * self._pow_cos_t[2][i] + 3.0)
        acc += real * real

        # Y_4_1
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_4_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_4_1[i] * self._sin_n_p[1][i]
        acc += real * real + imag * imag

        # Y_4_2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_2[i] * self._cos_n_p[2][i]
            imag += pre_y_4_2[i] * self._sin_n_p[2][i]
        acc += real * real + imag * imag

        # Y_4_3
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_4_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_4_3[i] * self._sin_n_p[3][i]
        acc += real * real + imag * imag

        # Y_4_4
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_4[i] * self._cos_n_p[4][i]
            imag += pre_y_4_4[i] * self._sin_n_p[4][i]
        acc += real * real + imag * imag

        q4 = sqrt(4 * pi * acc / (9 * float(nnn * nnn)))
        return q4

    def get_q6(self, thetas=None, phis=None):
        """
        Calculates the value of the bond orientational order parameter of
        weight l=6. If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh. Otherwise, it is expected that the
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
        sqrt_91_2pi = sqrt(91.0 / (2 * pi))
        sqrt_1365_pi = sqrt(1365.0 / pi)
        sqrt_273_2pi = sqrt(273.0 / (2 * pi))
        sqrt_13_pi = sqrt(13.0 / pi)

        pre_y_6_6 = [i64 * sqrt_3003_pi * val for val in self._pow_sin_t[6]]
        pre_y_6_5 = [i32_3 * sqrt_1001_pi * val[0] * val[1] for val in zip(self._pow_sin_t[5], self._pow_cos_t[1])]
        pre_y_6_4 = [
            i32_3 * sqrt_91_2pi * val[0] * (11 * val[1] - 1.0) for val in zip(self._pow_sin_t[4], self._pow_cos_t[2])
        ]
        pre_y_6_3 = [
            i32 * sqrt_1365_pi * val[0] * (11 * val[1] - 3 * val[2])
            for val in zip(self._pow_sin_t[3], self._pow_cos_t[3], self._pow_cos_t[1])
        ]
        pre_y_6_2 = [
            i64 * sqrt_1365_pi * val[0] * (33 * val[1] - 18 * val[2] + 1.0)
            for val in zip(self._pow_sin_t[2], self._pow_cos_t[4], self._pow_cos_t[2])
        ]
        pre_y_6_1 = [
            i16 * sqrt_273_2pi * val[0] * (33 * val[1] - 30 * val[2] + 5 * val[3])
            for val in zip(
                self._pow_sin_t[1],
                self._pow_cos_t[5],
                self._pow_cos_t[3],
                self._pow_cos_t[1],
            )
        ]

        acc = 0.0

        # Y_6_-6
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i]  # cos(x) =  cos(-x)
            imag -= pre_y_6_6[i] * self._sin_n_p[6][i]  # sin(x) = -sin(-x)
        acc += real * real + imag * imag

        # Y_6_-5
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += real * real + imag * imag

        # Y_6_-4
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += real * real + imag * imag

        # Y_6_-3
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += real * real + imag * imag

        # Y_6_-2
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += real * real + imag * imag

        # Y_6_-1
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += real * real + imag * imag

        # Y_6_0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += (
                i32
                * sqrt_13_pi
                * (231 * self._pow_cos_t[6][i] - 315 * self._pow_cos_t[4][i] + 105 * self._pow_cos_t[2][i] - 5.0)
            )
        acc += real * real

        # Y_6_1
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += real * real + imag * imag

        # Y_6_2
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag += pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += real * real + imag * imag

        # Y_6_3
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += real * real + imag * imag

        # Y_6_4
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag += pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += real * real + imag * imag

        # Y_6_5
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += real * real + imag * imag

        # Y_6_6
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i]
            imag += pre_y_6_6[i] * self._sin_n_p[6][i]
        acc += real * real + imag * imag

        q6 = sqrt(4 * pi * acc / (13 * float(nnn * nnn)))
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
            raise ValueError("Index for getting order parameter type out-of-bounds!")
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
            raise ValueError("Index for getting parameters associated with order parameter calculation out-of-bounds!")
        return self._params[index]

    def get_order_parameters(
        self, structure: Structure, n: int, indices_neighs: list[int] | None = None, tol: float = 0.0, target_spec=None
    ):
        """
        Compute all order parameters of site n.

        Args:
            structure (Structure): input structure.
            n (int): index of site in input structure,
                for which OPs are to be
                calculated. Note that we do not use the sites iterator
                here, but directly access sites via struct[index].
            indices_neighs (list[int]): list of indices of those neighbors
                in Structure object
                structure that are to be considered for OP computation.
                This optional argument overwrites the way neighbors are
                to be determined as defined in the constructor (i.e.,
                Voronoi coordination finder via negative cutoff radius
                vs constant cutoff radius if cutoff was positive).
                We do not use information about the underlying
                structure lattice if the neighbor indices are explicitly
                provided. This has two important consequences. First,
                the input Structure object can, in fact, be a
                simple list of Site objects. Second, no nearest images
                of neighbors are determined when providing an index list.
                Note furthermore that this neighbor
                determination type ignores the optional target_spec
                argument.
            tol (float): threshold of weight
                (= solid angle / maximal solid angle)
                to determine if a particular pair is
                considered neighbors; this is relevant only in the case
                when Voronoi polyhedra are used to determine coordination
            target_spec (Species): target species to be considered
                when calculating the order
                parameters of site n; None includes all species of input
                structure.

        Returns:
            [floats]: representing order parameters. Should it not be
            possible to compute a given OP for a conceptual reason, the
            corresponding entry is None instead of a float. For Steinhardt
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
            neighsitestmp = [i[0] for i in structure.get_sites_in_sphere(centsite.coords, self._cutoff)]
            neighsites = []
            if centsite not in neighsitestmp:
                raise ValueError("Could not find center site!")
            neighsitestmp.remove(centsite)

            if target_spec is None:
                neighsites = list(neighsitestmp)
            else:
                neighsites[:] = [site for site in neighsitestmp if site.specie.symbol == target_spec]
        nneigh = len(neighsites)
        self._last_nneigh = nneigh

        # Prepare angle calculations, if applicable.
        rij: list[np.ndarray] = []
        rjk: list[list[np.ndarray]] = []
        rijnorm: list[list[float]] = []
        rjknorm: list[list[np.ndarray]] = []
        dist: list[float] = []
        distjk_unique: list[float] = []
        distjk: list[list[float]] = []
        centvec = centsite.coords
        if self._computerijs:
            for j, neigh in enumerate(neighsites):
                rij.append(neigh.coords - centvec)
                dist.append(float(np.linalg.norm(rij[j])))
                rijnorm.append(rij[j] / dist[j])  # type: ignore
        if self._computerjks:
            for j, neigh in enumerate(neighsites):
                rjk.append([])
                rjknorm.append([])
                distjk.append([])
                kk = 0
                for k, neigh_2 in enumerate(neighsites):
                    if j != k:
                        rjk[j].append(neigh_2.coords - neigh.coords)
                        distjk[j].append(float(np.linalg.norm(rjk[j][kk])))
                        if k > j:
                            distjk_unique.append(distjk[j][kk])
                        rjknorm[j].append(rjk[j][kk] / distjk[j][kk])
                        kk = kk + 1
        # Initialize OP list and, then, calculate OPs.
        ops = [0.0 for t in self._types]
        # norms = [[[] for j in range(nneigh)] for t in self._types]

        # First, coordination number and distance-based OPs.
        for i, t in enumerate(self._types):
            if t == "cn":
                ops[i] = nneigh / self._params[i]["norm"]
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
            for vec in rijnorm:

                # z is North pole --> theta between vec and (0, 0, 1)^T.
                # Because vec is normalized, dot product is simply vec[2].
                thetas.append(acos(max(-1.0, min(vec[2], 1.0))))
                tmpphi = 0.0

                # Compute phi only if it is not (almost) perfectly
                # aligned with z-axis.
                if -left_of_unity < vec[2] < left_of_unity:
                    # x is prime meridian --> phi between projection of vec
                    # into x-y plane and (1, 0, 0)^T
                    tmpphi = acos(
                        max(
                            -1.0,
                            min(vec[0] / (sqrt(vec[0] * vec[0] + vec[1] * vec[1])), 1.0),
                        )
                    )
                    if vec[1] < 0.0:
                        tmpphi = -tmpphi
                phis.append(tmpphi)

            # Note that None flags that we have too few neighbors
            # for calculating BOOPS.
            for i, t in enumerate(self._types):
                if t == "q2":
                    ops[i] = self.get_q2(thetas, phis) if len(thetas) > 0 else None
                elif t == "q4":
                    ops[i] = self.get_q4(thetas, phis) if len(thetas) > 0 else None
                elif t == "q6":
                    ops[i] = self.get_q6(thetas, phis) if len(thetas) > 0 else None

        # Then, deal with the Peters-style OPs that are tailor-made
        # to recognize common structural motifs
        # (Peters, J. Chem. Phys., 131, 244103, 2009;
        #  Zimmermann et al., J. Am. Chem. Soc., under revision, 2015).
        if self._geomops:
            gaussthetak = [0.0 for t in self._types]  # not used by all OPs
            qsptheta = [[[] for j in range(nneigh)] for t in self._types]  # type: ignore
            norms = [[[] for j in range(nneigh)] for t in self._types]  # type: ignore
            ipi = 1.0 / pi
            piover2 = pi / 2.0
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
                        tmp = max(-1.0, min(np.inner(zaxis, rijnorm[k]), 1.0))
                        thetak = acos(tmp)
                        xaxis = gramschmidt(rijnorm[k], zaxis)
                        if np.linalg.norm(xaxis) < very_small:
                            flag_xaxis = True
                        else:
                            xaxis = xaxis / np.linalg.norm(xaxis)
                            flag_xaxis = False
                        if self._comp_azi:
                            flag_yaxis = True
                            yaxis = np.cross(zaxis, xaxis)
                            if np.linalg.norm(yaxis) > very_small:
                                yaxis = yaxis / np.linalg.norm(yaxis)
                                flag_yaxis = False

                        # Contributions of j-i-k angles, where i represents the
                        # central atom and j and k two of the neighbors.
                        for i, t in enumerate(self._types):
                            if t in ["bent", "sq_pyr_legacy"]:
                                tmp = self._params[i]["IGW_TA"] * (thetak * ipi - self._params[i]["TA"])
                                qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                norms[i][j][kc] += 1
                            elif t in ["tri_plan", "tri_plan_max", "tet", "tet_max"]:
                                tmp = self._params[i]["IGW_TA"] * (thetak * ipi - self._params[i]["TA"])
                                gaussthetak[i] = exp(-0.5 * tmp * tmp)
                                if t in ["tri_plan_max", "tet_max"]:
                                    qsptheta[i][j][kc] += gaussthetak[i]
                                    norms[i][j][kc] += 1
                            elif t in ["T", "tri_pyr", "sq_pyr", "pent_pyr", "hex_pyr"]:
                                tmp = self._params[i]["IGW_EP"] * (thetak * ipi - 0.5)
                                qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                norms[i][j][kc] += 1
                            elif t in [
                                "sq_plan",
                                "oct",
                                "oct_legacy",
                                "cuboct",
                                "cuboct_max",
                            ]:
                                if thetak >= self._params[i]["min_SPP"]:
                                    tmp = self._params[i]["IGW_SPP"] * (thetak * ipi - 1.0)
                                    qsptheta[i][j][kc] += self._params[i]["w_SPP"] * exp(-0.5 * tmp * tmp)
                                    norms[i][j][kc] += self._params[i]["w_SPP"]
                            elif t in [
                                "see_saw_rect",
                                "tri_bipyr",
                                "sq_bipyr",
                                "pent_bipyr",
                                "hex_bipyr",
                                "oct_max",
                                "sq_plan_max",
                                "hex_plan_max",
                            ]:
                                if thetak < self._params[i]["min_SPP"]:
                                    tmp = (
                                        self._params[i]["IGW_EP"] * (thetak * ipi - 0.5)
                                        if t != "hex_plan_max"
                                        else self._params[i]["IGW_TA"]
                                        * (fabs(thetak * ipi - 0.5) - self._params[i]["TA"])
                                    )
                                    qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                    norms[i][j][kc] += 1
                            elif t in ["pent_plan", "pent_plan_max"]:
                                tmp = 0.4 if thetak <= self._params[i]["TA"] * pi else 0.8
                                tmp2 = self._params[i]["IGW_TA"] * (thetak * ipi - tmp)
                                gaussthetak[i] = exp(-0.5 * tmp2 * tmp2)
                                if t == "pent_plan_max":
                                    qsptheta[i][j][kc] += gaussthetak[i]
                                    norms[i][j][kc] += 1
                            elif t == "bcc" and j < k:
                                if thetak >= self._params[i]["min_SPP"]:
                                    tmp = self._params[i]["IGW_SPP"] * (thetak * ipi - 1.0)
                                    qsptheta[i][j][kc] += self._params[i]["w_SPP"] * exp(-0.5 * tmp * tmp)
                                    norms[i][j][kc] += self._params[i]["w_SPP"]
                            elif t == "sq_face_cap_trig_pris":
                                if thetak < self._params[i]["TA3"]:
                                    tmp = self._params[i]["IGW_TA1"] * (thetak * ipi - self._params[i]["TA1"])
                                    qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                    norms[i][j][kc] += 1

                        for m in range(nneigh):
                            if (m != j) and (m != k) and (not flag_xaxis):
                                tmp = max(-1.0, min(np.inner(zaxis, rijnorm[m]), 1.0))
                                thetam = acos(tmp)
                                xtwoaxistmp = gramschmidt(rijnorm[m], zaxis)
                                l = np.linalg.norm(xtwoaxistmp)
                                if l < very_small:
                                    flag_xtwoaxis = True
                                else:
                                    xtwoaxis = xtwoaxistmp / l
                                    phi = acos(max(-1.0, min(np.inner(xtwoaxis, xaxis), 1.0)))
                                    flag_xtwoaxis = False
                                    if self._comp_azi:
                                        phi2 = atan2(
                                            np.dot(xtwoaxis, yaxis),
                                            np.dot(xtwoaxis, xaxis),
                                        )
                                # South pole contributions of m.
                                if t in [
                                    "tri_bipyr",
                                    "sq_bipyr",
                                    "pent_bipyr",
                                    "hex_bipyr",
                                    "oct_max",
                                    "sq_plan_max",
                                    "hex_plan_max",
                                    "see_saw_rect",
                                ]:
                                    if thetam >= self._params[i]["min_SPP"]:
                                        tmp = self._params[i]["IGW_SPP"] * (thetam * ipi - 1.0)
                                        qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp)
                                        norms[i][j][kc] += 1

                                # Contributions of j-i-m angle and
                                # angles between plane j-i-k and i-m vector.
                                if not flag_xaxis and not flag_xtwoaxis:
                                    for i, t in enumerate(self._types):
                                        if t in [
                                            "tri_plan",
                                            "tri_plan_max",
                                            "tet",
                                            "tet_max",
                                        ]:
                                            tmp = self._params[i]["IGW_TA"] * (thetam * ipi - self._params[i]["TA"])
                                            tmp2 = cos(self._params[i]["fac_AA"] * phi) ** self._params[i]["exp_cos_AA"]
                                            tmp3 = 1 if t in ["tri_plan_max", "tet_max"] else gaussthetak[i]
                                            qsptheta[i][j][kc] += tmp3 * exp(-0.5 * tmp * tmp) * tmp2
                                            norms[i][j][kc] += 1
                                        elif t in ["pent_plan", "pent_plan_max"]:
                                            tmp = 0.4 if thetam <= self._params[i]["TA"] * pi else 0.8
                                            tmp2 = self._params[i]["IGW_TA"] * (thetam * ipi - tmp)
                                            tmp3 = cos(phi)
                                            tmp4 = 1 if t == "pent_plan_max" else gaussthetak[i]
                                            qsptheta[i][j][kc] += tmp4 * exp(-0.5 * tmp2 * tmp2) * tmp3 * tmp3
                                            norms[i][j][kc] += 1
                                        elif t in [
                                            "T",
                                            "tri_pyr",
                                            "sq_pyr",
                                            "pent_pyr",
                                            "hex_pyr",
                                        ]:
                                            tmp = cos(self._params[i]["fac_AA"] * phi) ** self._params[i]["exp_cos_AA"]
                                            tmp3 = self._params[i]["IGW_EP"] * (thetam * ipi - 0.5)
                                            qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp3 * tmp3)
                                            norms[i][j][kc] += 1
                                        elif t in ["sq_plan", "oct", "oct_legacy"]:
                                            if (
                                                thetak < self._params[i]["min_SPP"]
                                                and thetam < self._params[i]["min_SPP"]
                                            ):
                                                tmp = (
                                                    cos(self._params[i]["fac_AA"] * phi)
                                                    ** self._params[i]["exp_cos_AA"]
                                                )
                                                tmp2 = self._params[i]["IGW_EP"] * (thetam * ipi - 0.5)
                                                qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp2 * tmp2)
                                                if t == "oct_legacy":
                                                    qsptheta[i][j][kc] -= tmp * self._params[i][6] * self._params[i][7]
                                                norms[i][j][kc] += 1
                                        elif t in [
                                            "tri_bipyr",
                                            "sq_bipyr",
                                            "pent_bipyr",
                                            "hex_bipyr",
                                            "oct_max",
                                            "sq_plan_max",
                                            "hex_plan_max",
                                        ]:
                                            if thetam < self._params[i]["min_SPP"]:
                                                if thetak < self._params[i]["min_SPP"]:
                                                    tmp = (
                                                        cos(self._params[i]["fac_AA"] * phi)
                                                        ** self._params[i]["exp_cos_AA"]
                                                    )
                                                    tmp2 = (
                                                        self._params[i]["IGW_EP"] * (thetam * ipi - 0.5)
                                                        if t != "hex_plan_max"
                                                        else self._params[i]["IGW_TA"]
                                                        * (fabs(thetam * ipi - 0.5) - self._params[i]["TA"])
                                                    )
                                                    qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1
                                        elif t == "bcc" and j < k:
                                            if thetak < self._params[i]["min_SPP"]:
                                                if thetak > piover2:
                                                    fac = 1.0
                                                else:
                                                    fac = -1.0
                                                tmp = (thetam - piover2) / asin(1 / 3)
                                                qsptheta[i][j][kc] += (
                                                    fac * cos(3 * phi) * fac_bcc * tmp * exp(-0.5 * tmp * tmp)
                                                )
                                                norms[i][j][kc] += 1
                                        elif t == "see_saw_rect":
                                            if thetam < self._params[i]["min_SPP"]:
                                                if thetak < self._params[i]["min_SPP"] and phi < 0.75 * pi:
                                                    tmp = (
                                                        cos(self._params[i]["fac_AA"] * phi)
                                                        ** self._params[i]["exp_cos_AA"]
                                                    )
                                                    tmp2 = self._params[i]["IGW_EP"] * (thetam * ipi - 0.5)
                                                    qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1.0
                                        elif t in ["cuboct", "cuboct_max"]:
                                            if (
                                                thetam < self._params[i]["min_SPP"]
                                                and self._params[i][4] < thetak < self._params[i][2]
                                            ):
                                                if self._params[i][4] < thetam < self._params[i][2]:
                                                    tmp = cos(phi)
                                                    tmp2 = self._params[i][5] * (thetam * ipi - 0.5)
                                                    qsptheta[i][j][kc] += tmp * tmp * exp(-0.5 * tmp2 * tmp2)
                                                    norms[i][j][kc] += 1.0
                                                elif thetam < self._params[i][4]:
                                                    tmp = 0.0556 * (cos(phi - 0.5 * pi) - 0.81649658)
                                                    tmp2 = self._params[i][6] * (thetam * ipi - onethird)
                                                    qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp) * exp(
                                                        -0.5 * tmp2 * tmp2
                                                    )
                                                    norms[i][j][kc] += 1.0
                                                elif thetam > self._params[i][2]:
                                                    tmp = 0.0556 * (cos(phi - 0.5 * pi) - 0.81649658)
                                                    tmp2 = self._params[i][6] * (thetam * ipi - twothird)
                                                    qsptheta[i][j][kc] += exp(-0.5 * tmp * tmp) * exp(
                                                        -0.5 * tmp2 * tmp2
                                                    )
                                                    norms[i][j][kc] += 1.0
                                        elif t == "sq_face_cap_trig_pris" and not flag_yaxis:
                                            if thetak < self._params[i]["TA3"]:
                                                if thetam < self._params[i]["TA3"]:
                                                    tmp = (
                                                        cos(self._params[i]["fac_AA1"] * phi2)
                                                        ** self._params[i]["exp_cos_AA1"]
                                                    )
                                                    tmp2 = self._params[i]["IGW_TA1"] * (
                                                        thetam * ipi - self._params[i]["TA1"]
                                                    )
                                                else:
                                                    tmp = (
                                                        cos(
                                                            self._params[i]["fac_AA2"]
                                                            * (phi2 + self._params[i]["shift_AA2"])
                                                        )
                                                        ** self._params[i]["exp_cos_AA2"]
                                                    )
                                                    tmp2 = self._params[i]["IGW_TA2"] * (
                                                        thetam * ipi - self._params[i]["TA2"]
                                                    )

                                                qsptheta[i][j][kc] += tmp * exp(-0.5 * tmp2 * tmp2)
                                                norms[i][j][kc] += 1

                        kc += 1

            # Normalize Peters-style OPs.
            for i, t in enumerate(self._types):
                if t in [
                    "tri_plan",
                    "tet",
                    "bent",
                    "sq_plan",
                    "oct",
                    "oct_legacy",
                    "cuboct",
                    "pent_plan",
                ]:
                    ops[i] = tmp_norm = 0.0
                    for j in range(nneigh):
                        ops[i] += sum(qsptheta[i][j])
                        tmp_norm += float(sum(norms[i][j]))
                    ops[i] = ops[i] / tmp_norm if tmp_norm > 1.0e-12 else None  # type: ignore
                elif t in [
                    "T",
                    "tri_pyr",
                    "see_saw_rect",
                    "sq_pyr",
                    "tri_bipyr",
                    "sq_bipyr",
                    "pent_pyr",
                    "hex_pyr",
                    "pent_bipyr",
                    "hex_bipyr",
                    "oct_max",
                    "tri_plan_max",
                    "tet_max",
                    "sq_plan_max",
                    "pent_plan_max",
                    "cuboct_max",
                    "hex_plan_max",
                    "sq_face_cap_trig_pris",
                ]:
                    ops[i] = None  # type: ignore
                    if nneigh > 1:
                        for j in range(nneigh):
                            for k in range(len(qsptheta[i][j])):
                                qsptheta[i][j][k] = (
                                    qsptheta[i][j][k] / norms[i][j][k] if norms[i][j][k] > 1.0e-12 else 0.0
                                )
                            ops[i] = max(qsptheta[i][j]) if j == 0 else max(ops[i], max(qsptheta[i][j]))
                elif t == "bcc":
                    ops[i] = 0.0
                    for j in range(nneigh):
                        ops[i] += sum(qsptheta[i][j])
                    if nneigh > 3:
                        ops[i] = ops[i] / float(0.5 * float(nneigh * (6 + (nneigh - 2) * (nneigh - 3))))
                    else:
                        ops[i] = None  # type: ignore
                elif t == "sq_pyr_legacy":
                    if nneigh > 1:
                        dmean = np.mean(dist)
                        acc = 0.0
                        for d in dist:
                            tmp = self._params[i][2] * (d - dmean)
                            acc = acc + exp(-0.5 * tmp * tmp)
                        for j in range(nneigh):
                            ops[i] = max(qsptheta[i][j]) if j == 0 else max(ops[i], max(qsptheta[i][j]))
                        ops[i] = acc * ops[i] / float(nneigh)
                        # nneigh * (nneigh - 1))
                    else:
                        ops[i] = None  # type: ignore

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
            for neigh in neighsites:
                neighscent = neighscent + neigh.coords
            if nneigh > 0:
                neighscent = neighscent / float(nneigh)
            h = np.linalg.norm(neighscent - centvec)
            b = min(distjk_unique) if len(distjk_unique) > 0 else 0
            dhalf = max(distjk_unique) / 2.0 if len(distjk_unique) > 0 else 0

            for i, t in enumerate(self._types):
                if t in ("reg_tri", "sq"):
                    if nneigh < 3:
                        ops[i] = None  # type: ignore
                    else:
                        ops[i] = 1.0
                        if t == "reg_tri":
                            a = 2 * asin(b / (2 * sqrt(h * h + (b / (2 * cos(3 * pi / 18))) ** 2)))  # type: ignore
                            nmax = 3
                        elif t == "sq":
                            a = 2 * asin(b / (2 * sqrt(h * h + dhalf * dhalf)))  # type: ignore
                            nmax = 4
                        for j in range(min([nneigh, nmax])):
                            ops[i] = ops[i] * exp(-0.5 * ((aijs[j] - a) * self._params[i][0]) ** 2)

        return ops


class BrunnerNN_reciprocal(NearNeighbors):
    """
    Determine coordination number using Brunner's algorithm which counts the
    atoms that are within the largest gap in differences in real space
    interatomic distances. This algorithm uses Brunner's method of
    largest reciprocal gap in interatomic distances.
    """

    def __init__(self, tol: float = 1.0e-4, cutoff=8.0):
        """
        Args:
            tol (float): tolerance parameter for bond determination
                (default: 1E-4).
            cutoff (float): cutoff radius in Angstrom to look for near-neighbor
                atoms. Defaults to 8.0.
        """
        self.tol = tol
        self.cutoff = cutoff

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        ds = sorted(i.nn_distance for i in neighs_dists)

        ns = [1.0 / ds[i] - 1.0 / ds[i + 1] for i in range(len(ds) - 1)]

        d_max = ds[ns.index(max(ns))]
        siw = []
        for nn in neighs_dists:
            s, dist = nn, nn.nn_distance
            if dist < d_max + self.tol:
                w = ds[0] / dist
                siw.append(
                    {
                        "site": s,
                        "image": self._get_image(structure, s),
                        "weight": w,
                        "site_index": self._get_original_site(structure, s),
                    }
                )
        return siw


class BrunnerNN_relative(NearNeighbors):
    """
    Determine coordination number using Brunner's algorithm which counts the
    atoms that are within the largest gap in differences in real space
    interatomic distances. This algorithm uses Brunner's method of
    of largest relative gap in interatomic distances.
    """

    def __init__(self, tol: float = 1.0e-4, cutoff=8.0):
        """
        Args:
            tol (float): tolerance parameter for bond determination
                (default: 1E-4).
            cutoff (float): cutoff radius in Angstrom to look for near-neighbor
                atoms. Defaults to 8.0.
        """
        self.tol = tol
        self.cutoff = cutoff

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        ds = sorted(i.nn_distance for i in neighs_dists)

        ns = [ds[i + 1] / ds[i] for i in range(len(ds) - 1)]

        d_max = ds[ns.index(max(ns))]
        siw = []
        for nn in neighs_dists:
            s, dist = nn, nn.nn_distance
            if dist < d_max + self.tol:
                w = ds[0] / dist
                siw.append(
                    {
                        "site": s,
                        "image": self._get_image(structure, s),
                        "weight": w,
                        "site_index": self._get_original_site(structure, s),
                    }
                )
        return siw


class BrunnerNN_real(NearNeighbors):
    """
    Determine coordination number using Brunner's algorithm which counts the
    atoms that are within the largest gap in differences in real space
    interatomic distances. This algorithm uses Brunner's method of
    largest gap in interatomic distances.
    """

    def __init__(self, tol: float = 1.0e-4, cutoff=8.0):
        """
        Args:
            tol (float): tolerance parameter for bond determination
                (default: 1E-4).
            cutoff (float): cutoff radius in Angstrom to look for near-neighbor
                atoms. Defaults to 8.0.
        """
        self.tol = tol
        self.cutoff = cutoff

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        ds = sorted(i.nn_distance for i in neighs_dists)

        ns = [ds[i + 1] - ds[i] for i in range(len(ds) - 1)]

        d_max = ds[ns.index(max(ns))]
        siw = []
        for nn in neighs_dists:
            s, dist = nn, nn.nn_distance
            if dist < d_max + self.tol:
                w = ds[0] / dist
                siw.append(
                    {
                        "site": s,
                        "image": self._get_image(structure, s),
                        "weight": w,
                        "site_index": self._get_original_site(structure, s),
                    }
                )
        return siw


class EconNN(NearNeighbors):
    """
    Determines the average effective coordination number for each cation in a
    given structure using Hoppe's algorithm.

    This method follows the procedure outlined in:

    Hoppe, Rudolf. "Effective coordination numbers (ECoN) and mean fictive ionic
    radii (MEFIR)." Zeitschrift für Kristallographie-Crystalline Materials
    150.1-4 (1979): 23-52.
    """

    def __init__(
        self,
        tol: float = 0.2,
        cutoff: float = 10.0,
        cation_anion: bool = False,
        use_fictive_radius: bool = False,
    ):
        """
        Args:
            tol: Tolerance parameter for bond determination.
            cutoff: Cutoff radius in Angstrom to look for near-neighbor atoms.
            cation_anion: If set to True, will restrict bonding targets to
                sites with opposite or zero charge. Requires an oxidation states
                on all sites in the structure.
            use_fictive_radius: Whether to use the fictive radius in the
                EcoN calculation. If False, the bond distance will be used.
        """
        self.tol = tol
        self.cutoff = cutoff
        self.cation_anion = cation_anion
        self.use_fictive_radius = use_fictive_radius

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return True

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        site = structure[n]
        neighbors = structure.get_neighbors(site, self.cutoff)

        if self.cation_anion and hasattr(site.specie, "oxi_state"):
            # filter out neighbor of like charge (except for neutral sites)
            if site.specie.oxi_state >= 0:
                neighbors = [n for n in neighbors if n.oxi_state <= 0]
            elif site.specie.oxi_state <= 0:
                neighbors = [n for n in neighbors if n.oxi_state >= 0]

        if self.use_fictive_radius:
            # calculate fictive ionic radii
            firs = [_get_fictive_ionic_radius(site, neighbor) for neighbor in neighbors]
        else:
            # just use the bond distance
            firs = [neighbor.nn_distance for neighbor in neighbors]

        # calculate mean fictive ionic radius
        mefir = _get_mean_fictive_ionic_radius(firs)

        # # iteratively solve MEFIR; follows equation 4 in Hoppe's EconN paper
        prev_mefir = float("inf")
        while abs(prev_mefir - mefir) > 1e-4:
            # this is guaranteed to converge
            prev_mefir = mefir
            mefir = _get_mean_fictive_ionic_radius(firs, minimum_fir=mefir)

        siw = []
        for nn, fir in zip(neighbors, firs):
            if nn.nn_distance < self.cutoff:
                w = exp(1 - (fir / mefir) ** 6)
                if w > self.tol:
                    bonded_site = {
                        "site": nn,
                        "image": self._get_image(structure, nn),
                        "weight": w,
                        "site_index": self._get_original_site(structure, nn),
                    }
                    siw.append(bonded_site)
        return siw


def _get_fictive_ionic_radius(site: Site, neighbor: PeriodicNeighbor) -> float:
    """
    Get fictive ionic radius.

    Follows equation 1 of:

    Hoppe, Rudolf. "Effective coordination numbers (ECoN) and mean fictive ionic
    radii (MEFIR)." Zeitschrift für Kristallographie-Crystalline Materials
    150.1-4 (1979): 23-52.

    Args:
        site: The central site.
        neighbor neighboring site.

    Returns:
        Hoppe's fictive ionic radius.
    """
    r_h = _get_radius(site)
    if r_h == 0:
        r_h = _get_default_radius(site)

    r_i = _get_radius(neighbor)
    if r_i == 0:
        r_i = _get_default_radius(neighbor)

    return neighbor.nn_distance * (r_h / (r_h + r_i))


def _get_mean_fictive_ionic_radius(
    fictive_ionic_radii: list[float],
    minimum_fir: float | None = None,
) -> float:
    """
    Returns the mean fictive ionic radius.

    Follows equation 2:

    Hoppe, Rudolf. "Effective coordination numbers (ECoN) and mean fictive ionic
    radii (MEFIR)." Zeitschrift für Kristallographie-Crystalline Materials
    150.1-4 (1979): 23-52.

    Args:
        fictive_ionic_radii: List of fictive ionic radii for a center site
            and its neighbors.
        minimum_fir: Minimum fictive ionic radius to use.

    Returns:
        Hoppe's mean fictive ionic radius.
    """
    if not minimum_fir:
        minimum_fir = min(fictive_ionic_radii)

    weighted_sum = 0.0
    total_sum = 0.0
    for fir in fictive_ionic_radii:
        weighted_sum += fir * exp(1 - (fir / minimum_fir) ** 6)
        total_sum += exp(1 - (fir / minimum_fir) ** 6)

    return weighted_sum / total_sum


class CrystalNN(NearNeighbors):
    """
    This is a custom near-neighbor method intended for use in all kinds of periodic structures
    (metals, minerals, porous structures, etc). It is based on a Voronoi algorithm and uses the
    solid angle weights to determine the probability of various coordination environments. The
    algorithm can also modify probability using smooth distance cutoffs as well as Pauling
    electronegativity differences. The output can either be the most probable coordination
    environment or a weighted list of coordination environments.
    """

    NNData = namedtuple("NNData", ["all_nninfo", "cn_weights", "cn_nninfo"])

    def __init__(
        self,
        weighted_cn=False,
        cation_anion=False,
        distance_cutoffs=(0.5, 1),
        x_diff_weight=3.0,
        porous_adjustment=True,
        search_cutoff=7,
        fingerprint_length=None,
    ):
        """
        Initialize CrystalNN with desired parameters. Default parameters assume
        "chemical bond" type behavior is desired. For geometric neighbor
        finding (e.g., structural framework), set (i) distance_cutoffs=None,
        (ii) x_diff_weight=0.0 and (optionally) (iii) porous_adjustment=False
        which will disregard the atomic identities and perform best for a purely
        geometric match.

        Args:
            weighted_cn: (bool) if set to True, will return fractional weights
                for each potential near neighbor.
            cation_anion: (bool) if set True, will restrict bonding targets to
                sites with opposite or zero charge. Requires an oxidation states
                on all sites in the structure.
            distance_cutoffs: ([float, float]) - if not None, penalizes neighbor
                distances greater than sum of covalent radii plus
                distance_cutoffs[0]. Distances greater than covalent radii sum
                plus distance_cutoffs[1] are enforced to have zero weight.
            x_diff_weight: (float) - if multiple types of neighbor elements are
                possible, this sets preferences for targets with higher
                electronegativity difference.
            porous_adjustment: (bool) - if True, readjusts Voronoi weights to
                better describe layered / porous structures
            search_cutoff: (float) cutoff in Angstroms for initial neighbor
                search; this will be adjusted if needed internally
            fingerprint_length: (int) if a fixed_length CN "fingerprint" is
                desired from get_nn_data(), set this parameter
        """
        self.weighted_cn = weighted_cn
        self.cation_anion = cation_anion
        self.distance_cutoffs = distance_cutoffs
        self.x_diff_weight = x_diff_weight if x_diff_weight is not None else 0
        self.search_cutoff = search_cutoff
        self.porous_adjustment = porous_adjustment
        self.fingerprint_length = fingerprint_length

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_nn_info(self, structure: Structure, n: int) -> list[dict]:
        """
        Get all near-neighbor information.
        Args:
            structure: (Structure) pymatgen Structure
            n: (int) index of target site

        Returns:
            siw (list[dict]): each dictionary provides information
                about a single near neighbor, where key 'site' gives access to the
                corresponding Site object, 'image' gives the image location, and
                'weight' provides the weight that a given near-neighbor site contributes
                to the coordination number (1 or smaller), 'site_index' gives index of
                the corresponding site in the original structure.
        """
        nndata = self.get_nn_data(structure, n)

        if not self.weighted_cn:
            max_key = max(nndata.cn_weights, key=lambda k: nndata.cn_weights[k])
            nn = nndata.cn_nninfo[max_key]
            for entry in nn:
                entry["weight"] = 1
            return nn

        for entry in nndata.all_nninfo:
            weight = 0
            for cn in nndata.cn_nninfo:
                for cn_entry in nndata.cn_nninfo[cn]:
                    if entry["site"] == cn_entry["site"]:
                        weight += nndata.cn_weights[cn]

            entry["weight"] = weight

        return nndata.all_nninfo

    def get_nn_data(self, structure: Structure, n: int, length=None):
        """
        The main logic of the method to compute near neighbor.

        Args:
            structure: (Structure) enclosing structure object
            n: (int) index of target site to get NN info for
            length: (int) if set, will return a fixed range of CN numbers

        Returns:
            a namedtuple (NNData) object that contains:
                - all near neighbor sites with weights
                - a dict of CN -> weight
                - a dict of CN -> associated near neighbor sites
        """
        length = length or self.fingerprint_length

        # determine possible bond targets
        target = None
        if self.cation_anion:
            target = []
            m_oxi = structure[n].specie.oxi_state
            for site in structure:
                if site.specie.oxi_state * m_oxi <= 0:  # opposite charge
                    target.append(site.specie)
            if not target:
                raise ValueError("No valid targets for site within cation_anion constraint!")

        # get base VoronoiNN targets
        cutoff = self.search_cutoff
        vnn = VoronoiNN(weight="solid_angle", targets=target, cutoff=cutoff)
        nn = vnn.get_nn_info(structure, n)

        # solid angle weights can be misleading in open / porous structures
        # adjust weights to correct for this behavior
        if self.porous_adjustment:
            for x in nn:
                x["weight"] *= x["poly_info"]["solid_angle"] / x["poly_info"]["area"]

        # adjust solid angle weight based on electronegativity difference
        if self.x_diff_weight > 0:
            for entry in nn:
                X1 = structure[n].specie.X
                X2 = entry["site"].specie.X

                if math.isnan(X1) or math.isnan(X2):
                    chemical_weight = 1
                else:
                    # note: 3.3 is max deltaX between 2 elements
                    chemical_weight = 1 + self.x_diff_weight * math.sqrt(abs(X1 - X2) / 3.3)

                entry["weight"] = entry["weight"] * chemical_weight

        # sort nearest neighbors from highest to lowest weight
        nn = sorted(nn, key=lambda x: x["weight"], reverse=True)
        if nn[0]["weight"] == 0:
            return self.transform_to_length(self.NNData([], {0: 1.0}, {0: []}), length)

        # renormalize weights so the highest weight is 1.0
        highest_weight = nn[0]["weight"]
        for entry in nn:
            entry["weight"] = entry["weight"] / highest_weight

        # adjust solid angle weights based on distance
        if self.distance_cutoffs:
            r1 = _get_radius(structure[n])
            for entry in nn:
                r2 = _get_radius(entry["site"])
                if r1 > 0 and r2 > 0:
                    d = r1 + r2
                else:
                    warnings.warn(
                        "CrystalNN: cannot locate an appropriate radius, "
                        "covalent or atomic radii will be used, this can lead "
                        "to non-optimal results."
                    )
                    d = _get_default_radius(structure[n]) + _get_default_radius(entry["site"])

                dist = np.linalg.norm(structure[n].coords - entry["site"].coords)
                dist_weight: float = 0

                cutoff_low = d + self.distance_cutoffs[0]
                cutoff_high = d + self.distance_cutoffs[1]

                if dist <= cutoff_low:
                    dist_weight = 1
                elif dist < cutoff_high:
                    dist_weight = (math.cos((dist - cutoff_low) / (cutoff_high - cutoff_low) * math.pi) + 1) * 0.5
                entry["weight"] = entry["weight"] * dist_weight

        # sort nearest neighbors from highest to lowest weight
        nn = sorted(nn, key=lambda x: x["weight"], reverse=True)
        if nn[0]["weight"] == 0:
            return self.transform_to_length(self.NNData([], {0: 1.0}, {0: []}), length)

        for entry in nn:
            entry["weight"] = round(entry["weight"], 3)
            del entry["poly_info"]  # trim

        # remove entries with no weight
        nn = [x for x in nn if x["weight"] > 0]

        # get the transition distances, i.e. all distinct weights
        dist_bins: list[float] = []
        for entry in nn:
            if not dist_bins or dist_bins[-1] != entry["weight"]:
                dist_bins.append(entry["weight"])
        dist_bins.append(0)

        # main algorithm to determine fingerprint from bond weights
        cn_weights = {}  # CN -> score for that CN
        cn_nninfo = {}  # CN -> list of nearneighbor info for that CN
        for idx, val in enumerate(dist_bins):
            if val != 0:
                nn_info = []
                for entry in nn:
                    if entry["weight"] >= val:
                        nn_info.append(entry)
                cn = len(nn_info)
                cn_nninfo[cn] = nn_info
                cn_weights[cn] = self._semicircle_integral(dist_bins, idx)

        # add zero coord
        cn0_weight = 1.0 - sum(cn_weights.values())
        if cn0_weight > 0:
            cn_nninfo[0] = []
            cn_weights[0] = cn0_weight

        return self.transform_to_length(self.NNData(nn, cn_weights, cn_nninfo), length)

    def get_cn(self, structure: Structure, n: int, **kwargs) -> float:  # type: ignore
        """
        Get coordination number, CN, of site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine CN.
            use_weights (bool): flag indicating whether (True)
                to use weights for computing the coordination number
                or not (False, default: each coordinated site has equal
                weight).
            on_disorder ('take_majority_strict' | 'take_majority_drop' | 'take_max_species' | 'error'):
                What to do when encountering a disordered structure. 'error' will raise ValueError.
                'take_majority_strict' will use the majority specie on each site and raise
                ValueError if no majority exists. 'take_max_species' will use the first max specie
                on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, 'error' and 'take_majority_strict'
                will raise ValueError, while 'take_majority_drop' ignores this site altogether and
                'take_max_species' will use Fe as the site specie.

        Returns:
            cn (int or float): coordination number.
        """
        use_weights = kwargs.get("use_weights", False)
        if self.weighted_cn != use_weights:
            raise ValueError("The weighted_cn parameter and use_weights parameter should match!")

        return super().get_cn(structure, n, **kwargs)

    def get_cn_dict(self, structure: Structure, n: int, use_weights: bool = False, **kwargs):
        """
        Get coordination number, CN, of each element bonded to site with index n in structure

        Args:
            structure (Structure): input structure
            n (int): index of site for which to determine CN.
            use_weights (bool): flag indicating whether (True)
                to use weights for computing the coordination number
                or not (False, default: each coordinated site has equal
                weight).

        Returns:
            cn (dict): dictionary of CN of each element bonded to site
        """
        if self.weighted_cn != use_weights:
            raise ValueError("The weighted_cn parameter and use_weights parameter should match!")

        return super().get_cn_dict(structure, n, use_weights)

    @staticmethod
    def _semicircle_integral(dist_bins, idx):
        """
        An internal method to get an integral between two bounds of a unit
        semicircle. Used in algorithm to determine bond probabilities.
        Args:
            dist_bins: (float) list of all possible bond weights
            idx: (float) index of starting bond weight

        Returns:
            (float) integral of portion of unit semicircle
        """
        r = 1

        x1 = dist_bins[idx]
        x2 = dist_bins[idx + 1]

        if dist_bins[idx] == 1:
            area1 = 0.25 * math.pi * r**2
        else:
            area1 = 0.5 * ((x1 * math.sqrt(r**2 - x1**2)) + (r**2 * math.atan(x1 / math.sqrt(r**2 - x1**2))))

        area2 = 0.5 * ((x2 * math.sqrt(r**2 - x2**2)) + (r**2 * math.atan(x2 / math.sqrt(r**2 - x2**2))))

        return (area1 - area2) / (0.25 * math.pi * r**2)

    @staticmethod
    def transform_to_length(nndata, length):
        """
        Given NNData, transforms data to the specified fingerprint length
        Args:
            nndata: (NNData)
            length: (int) desired length of NNData
        """
        if length is None:
            return nndata

        if length:
            for cn in range(length):
                if cn not in nndata.cn_weights:
                    nndata.cn_weights[cn] = 0
                    nndata.cn_nninfo[cn] = []

        return nndata


def _get_default_radius(site):
    """
    An internal method to get a "default" covalent/element radius

    Args:
        site: (Site)

    Returns:
        Covalent radius of element on site, or Atomic radius if unavailable
    """
    try:
        return CovalentRadius.radius[site.specie.symbol]
    except Exception:
        return site.specie.atomic_radius


def _get_radius(site):
    """
    An internal method to get the expected radius for a site with
    oxidation state.
    Args:
        site: (Site)

    Returns:
        Oxidation-state dependent radius: ionic, covalent, or atomic.
        Returns 0 if no oxidation state or appropriate radius is found.
    """
    if hasattr(site.specie, "oxi_state"):
        el = site.specie.element
        oxi = site.specie.oxi_state

        if oxi == 0:
            return _get_default_radius(site)

        if oxi in el.ionic_radii:
            return el.ionic_radii[oxi]

        # e.g., oxi = 2.667, average together 2+ and 3+ radii
        if int(math.floor(oxi)) in el.ionic_radii and int(math.ceil(oxi)) in el.ionic_radii:
            oxi_low = el.ionic_radii[int(math.floor(oxi))]
            oxi_high = el.ionic_radii[int(math.ceil(oxi))]
            x = oxi - int(math.floor(oxi))
            return (1 - x) * oxi_low + x * oxi_high

        if oxi > 0 and el.average_cationic_radius > 0:
            return el.average_cationic_radius

        if el.average_anionic_radius > 0 > oxi:
            return el.average_anionic_radius

    else:
        warnings.warn(
            "No oxidation states specified on sites! For better results, set "
            "the site oxidation states in the structure."
        )
    return 0


class CutOffDictNN(NearNeighbors):
    """
    A basic NN class using a dictionary of fixed cut-off distances.
    Only pairs of elements listed in the cut-off dictionary are considered
    during construction of the neighbor lists.

    Omit passing a dictionary for a Null/Empty NN class.
    """

    def __init__(self, cut_off_dict=None):
        """
        Args:
            cut_off_dict (dict[str, float]): a dictionary
            of cut-off distances, e.g. {('Fe','O'): 2.0} for
            a maximum Fe-O bond length of 2.0 Angstroms.
            Bonds will only be created between pairs listed
            in the cut-off dictionary.
            If your structure is oxidation state decorated,
            the cut-off distances will have to explicitly include
            the oxidation state, e.g. {('Fe2+', 'O2-'): 2.0}
        """
        self.cut_off_dict = cut_off_dict or {}

        # for convenience
        self._max_dist = 0.0
        lookup_dict = defaultdict(dict)
        for (sp1, sp2), dist in self.cut_off_dict.items():
            lookup_dict[sp1][sp2] = dist
            lookup_dict[sp2][sp1] = dist
            if dist > self._max_dist:
                self._max_dist = dist
        self._lookup_dict = lookup_dict

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return True

    @staticmethod
    def from_preset(preset):
        """
        Initialise a CutOffDictNN according to a preset set of cut-offs.

        Args:
            preset (str): A preset name. The list of supported presets are:

                - "vesta_2019": The distance cut-offs used by the VESTA
                  visualisation program.

        Returns:
            A CutOffDictNN using the preset cut-off dictionary.
        """
        if preset == "vesta_2019":
            cut_offs = loadfn(os.path.join(_directory, "vesta_cutoffs.yaml"))
            return CutOffDictNN(cut_off_dict=cut_offs)

        raise ValueError(f"Unrecognised preset: {preset}")

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        site = structure[n]

        neighs_dists = structure.get_neighbors(site, self._max_dist)

        nn_info = []
        for nn in neighs_dists:
            n_site = nn
            dist = nn.nn_distance
            neigh_cut_off_dist = self._lookup_dict.get(site.species_string, {}).get(n_site.species_string, 0.0)

            if dist < neigh_cut_off_dist:
                nn_info.append(
                    {
                        "site": n_site,
                        "image": self._get_image(structure, n_site),
                        "weight": dist,
                        "site_index": self._get_original_site(structure, n_site),
                    }
                )

        return nn_info


class Critic2NN(NearNeighbors):
    """
    Performs a topological analysis using critic2 to obtain
    neighbor information, using a sum of atomic charge
    densities. If an actual charge density is available
    (e.g. from a VASP CHGCAR), see Critic2Caller directly
    instead.
    """

    def __init__(self):
        """
        Init for Critic2NN.
        """
        # we cache the last-used structure, in case user
        # calls get_nn_info() repeatedly for different
        # sites in the same structure to save redundant
        # computations
        self.__last_structure = None
        self.__last_bonded_structure = None

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return True

    @property
    def extend_structure_molecules(self):
        """
        Boolean property: Do Molecules need to be converted to Structures to use
        this NearNeighbors class? Note: this property is not defined for classes
        for which molecules_allowed == False.
        """
        return True

    def get_bonded_structure(self, structure: Structure, decorate: bool = False) -> StructureGraph:  # type: ignore
        """
        Args:
            structure (Structure): Input structure
            decorate (bool, optional): Whether to decorate the structure. Defaults to False.

        Returns:
            StructureGraph: Bonded structure
        """
        # not a top-level import because critic2 is an optional
        # dependency, only want to raise an import error if
        # Critic2NN() is used
        from pymatgen.command_line.critic2_caller import Critic2Caller

        if structure == self.__last_structure:
            sg = self.__last_bonded_structure
        else:
            c2_output = Critic2Caller(structure).output
            sg = c2_output.structure_graph()

            self.__last_structure = structure
            self.__last_bonded_structure = sg

        if decorate:
            order_parameters = [self.get_local_order_parameters(structure, n) for n in range(len(structure))]
            sg.structure.add_site_property("order_parameters", order_parameters)

        return sg

    def get_nn_info(self, structure: Structure, n: int):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (int): index of site for which to determine near-neighbor
                sites.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a coordinated site, its image location,
                and its weight.
        """
        sg = self.get_bonded_structure(structure)

        return [
            {
                "site": connected_site.site,
                "image": connected_site.jimage,
                "weight": connected_site.weight,
                "site_index": connected_site.index,
            }
            for connected_site in sg.get_connected_sites(n)
        ]


def metal_edge_extender(
    mol_graph,
    cutoff: float = 2.5,
    metals: list | tuple | None = ("Li", "Mg", "Ca", "Zn", "B", "Al"),
    coordinators: list | tuple = ("O", "N", "F", "S", "Cl"),
):
    """
    Function to identify and add missed coordinate bond edges for metals

    Args:
        mol_graph: pymatgen.analysis.graphs.MoleculeGraph object
        cutoff: cutoff in Angstrom. Metal-coordinator sites that are closer
            together than this value will be considered coordination bonds.
            If the MoleculeGraph contains a metal, but no coordination bonds are found
            with the chosen cutoff, the cutoff will be increased by 1 Angstrom
            and another attempt will be made to identify coordination bonds.
        metals: Species considered metals for the purpose of identifying
            missed coordinate bond edges. The set {"Li", "Mg", "Ca", "Zn", "B", "Al"}
            (default) corresponds to the settings used in the LIBE dataset.
            Alternatively, set to None to cause any Species classified as a metal
            by Specie.is_metal to be considered a metal.
        coordinators: Possible coordinating species to consider when identifying
            missed coordinate bonds. The default set {"O", "N", "F", "S", "Cl"} was
            used in the LIBE dataset.

    Returns:
        mol_graph: pymatgen.analysis.graphs.MoleculeGraph object with additional
            metal bonds (if any found) added

    """
    if metals is None:
        metals = []
        for idx in mol_graph.graph.nodes():
            if Species(mol_graph.graph.nodes()[idx]["specie"]).is_metal:
                metals.append(mol_graph.graph.nodes()[idx]["specie"])
    metal_sites: dict = {k: {} for k in metals}

    num_new_edges = 0
    for idx in mol_graph.graph.nodes():
        if mol_graph.graph.nodes()[idx]["specie"] in metals:
            metal_sites[mol_graph.graph.nodes()[idx]["specie"]][idx] = [
                site[2] for site in mol_graph.get_connected_sites(idx)
            ]
    for sites in metal_sites.values():
        for idx, indices in sites.items():
            for ii, site in enumerate(mol_graph.molecule):
                if ii != idx and ii not in indices:
                    if str(site.specie) in coordinators:
                        if site.distance(mol_graph.molecule[idx]) < cutoff:
                            mol_graph.add_edge(idx, ii)
                            num_new_edges += 1
                            indices.append(ii)
    # If no metal edges are found, increase cutoff by 1 Ang and repeat analysis
    total_metal_edges = 0
    for sites in metal_sites.values():
        for indices in sites.values():
            total_metal_edges += len(indices)
    if total_metal_edges == 0:
        for sites in metal_sites.values():
            for idx, indices in sites.items():
                for ii, site in enumerate(mol_graph.molecule):
                    if ii != idx and ii not in indices:
                        if str(site.specie) in coordinators:
                            if site.distance(mol_graph.molecule[idx]) < (cutoff + 1):
                                mol_graph.add_edge(idx, ii)
                                num_new_edges += 1
                                indices.append(ii)

    return mol_graph
