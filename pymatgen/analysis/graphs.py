# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Module for graph representations of crystals.
"""

import copy
import logging
import os.path
import subprocess
import warnings
from collections import defaultdict, namedtuple
from itertools import combinations
from operator import itemgetter

import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np
from monty.json import MSONable
from monty.os.path import which
from networkx.drawing.nx_agraph import write_dot
from networkx.readwrite import json_graph
from scipy.spatial import KDTree
from scipy.stats import describe

from pymatgen.core import Lattice, Molecule, PeriodicSite, Structure
from pymatgen.core.structure import FunctionalGroups
from pymatgen.util.coord import lattice_points_in_supercell
from pymatgen.vis.structure_vtk import EL_COLORS

try:
    import igraph

    IGRAPH_AVAILABLE = True
except ImportError:
    IGRAPH_AVAILABLE = False


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

__author__ = "Matthew Horton, Evan Spotte-Smith, Samuel Blau"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Production"
__date__ = "August 2017"

ConnectedSite = namedtuple("ConnectedSite", "site, jimage, index, weight, dist")


def _compare(g1, g2, i1, i2):
    """
    Helper function called by isomorphic to ensure comparison of node identities.
    """
    return g1.vs[i1]["species"] == g2.vs[i2]["species"]


def _igraph_from_nxgraph(graph):
    """
    Helper function that converts a networkx graph object into an igraph graph object.
    """
    nodes = graph.nodes(data=True)
    new_igraph = igraph.Graph()
    for node in nodes:
        new_igraph.add_vertex(name=str(node[0]), species=node[1]["specie"], coords=node[1]["coords"])
    new_igraph.add_edges([(str(edge[0]), str(edge[1])) for edge in graph.edges()])
    return new_igraph


def _isomorphic(frag1, frag2):
    """
    Internal function to check if two graph objects are isomorphic, using igraph if
    if is available and networkx if it is not.
    """
    f1_nodes = frag1.nodes(data=True)
    f2_nodes = frag2.nodes(data=True)
    if len(f1_nodes) != len(f2_nodes):
        return False
    f2_edges = frag2.edges()
    if len(f2_edges) != len(f2_edges):
        return False
    f1_comp_dict = {}
    f2_comp_dict = {}
    for node in f1_nodes:
        if node[1]["specie"] not in f1_comp_dict:
            f1_comp_dict[node[1]["specie"]] = 1
        else:
            f1_comp_dict[node[1]["specie"]] += 1
    for node in f2_nodes:
        if node[1]["specie"] not in f2_comp_dict:
            f2_comp_dict[node[1]["specie"]] = 1
        else:
            f2_comp_dict[node[1]["specie"]] += 1
    if f1_comp_dict != f2_comp_dict:
        return False
    if IGRAPH_AVAILABLE:
        ifrag1 = _igraph_from_nxgraph(frag1)
        ifrag2 = _igraph_from_nxgraph(frag2)
        return ifrag1.isomorphic_vf2(ifrag2, node_compat_fn=_compare)
    nm = iso.categorical_node_match("specie", "ERROR")
    return nx.is_isomorphic(frag1.to_undirected(), frag2.to_undirected(), node_match=nm)


class StructureGraph(MSONable):
    """
    This is a class for annotating a Structure with
    bond information, stored in the form of a graph. A "bond" does
    not necessarily have to be a chemical bond, but can store any
    kind of information that connects two Sites.
    """

    def __init__(self, structure, graph_data=None):
        """
        If constructing this class manually, use the `with_empty_graph`
        method or `with_local_env_strategy` method (using an algorithm
        provided by the `local_env` module, such as O'Keeffe).

        This class that contains connection information:
        relationships between sites represented by a Graph structure,
        and an associated structure object.

        This class uses the NetworkX package to store and operate
        on the graph itself, but contains a lot of helper methods
        to make associating a graph with a given crystallographic
        structure easier.

        Use cases for this include storing bonding information,
        NMR J-couplings, Heisenberg exchange parameters, etc.

        For periodic graphs, class stores information on the graph
        edges of what lattice image the edge belongs to.

        :param structure: a Structure object

        :param graph_data: dict containing graph information in
            dict format (not intended to be constructed manually,
        see as_dict method for format)
        """

        if isinstance(structure, StructureGraph):
            # just make a copy from input
            graph_data = structure.as_dict()["graphs"]

        self.structure = structure
        self.graph = nx.readwrite.json_graph.adjacency_graph(graph_data)

        # tidy up edge attr dicts, reading to/from json duplicates
        # information
        for u, v, k, d in self.graph.edges(keys=True, data=True):
            if "id" in d:
                del d["id"]
            if "key" in d:
                del d["key"]
            # ensure images are tuples (conversion to lists happens
            # when serializing back from json), it's important images
            # are hashable/immutable
            if "to_jimage" in d:
                d["to_jimage"] = tuple(d["to_jimage"])
            if "from_jimage" in d:
                d["from_jimage"] = tuple(d["from_jimage"])

    @classmethod
    def with_empty_graph(cls, structure, name="bonds", edge_weight_name=None, edge_weight_units=None):
        """
        Constructor for StructureGraph, returns a StructureGraph
        object with an empty graph (no edges, only nodes defined
        that correspond to Sites in Structure).

        :param structure (Structure):
        :param name (str): name of graph, e.g. "bonds"
        :param edge_weight_name (str): name of edge weights,
            e.g. "bond_length" or "exchange_constant"
        :param edge_weight_units (str): name of edge weight units
            e.g. "Å" or "eV"
        :return (StructureGraph):
        """

        if edge_weight_name and (edge_weight_units is None):
            raise ValueError(
                "Please specify units associated "
                "with your edge weights. Can be "
                "empty string if arbitrary or "
                "dimensionless."
            )

        # construct graph with one node per site
        # graph attributes don't change behavior of graph,
        # they're just for book-keeping
        graph = nx.MultiDiGraph(
            edge_weight_name=edge_weight_name,
            edge_weight_units=edge_weight_units,
            name=name,
        )
        graph.add_nodes_from(range(len(structure)))

        graph_data = json_graph.adjacency_data(graph)

        return cls(structure, graph_data=graph_data)

    @staticmethod
    def with_edges(structure, edges):
        """
        Constructor for MoleculeGraph, using pre-existing or pre-defined edges
        with optional edge parameters.

        :param molecule: Molecule object
        :param edges: dict representing the bonds of the functional
            group (format: {(from_index, to_index, from_image, to_image): props},
            where props is a dictionary of properties, including weight.
            Props should be None if no additional properties are to be
            specified.
        :return: sg, a StructureGraph
        """

        sg = StructureGraph.with_empty_graph(structure, name="bonds", edge_weight_name="weight", edge_weight_units="")

        for edge, props in edges.items():

            try:
                from_index = edge[0]
                to_index = edge[1]
                from_image = edge[2]
                to_image = edge[3]
            except TypeError:
                raise ValueError("Edges must be given as (from_index, to_index," " from_image, to_image) tuples")

            if props is not None:
                if "weight" in props.keys():
                    weight = props["weight"]
                    del props["weight"]
                else:
                    weight = None

                if len(props.items()) == 0:
                    props = None
            else:
                weight = None

            nodes = sg.graph.nodes
            if not (from_index in nodes and to_index in nodes):
                raise ValueError(
                    "Edges cannot be added if nodes are not" " present in the graph. Please check your" " indices."
                )

            sg.add_edge(
                from_index,
                to_index,
                from_jimage=from_image,
                to_jimage=to_image,
                weight=weight,
                edge_properties=props,
            )

        sg.set_node_attributes()
        return sg

    @staticmethod
    def with_local_env_strategy(structure, strategy, weights=False):
        """
        Constructor for StructureGraph, using a strategy
        from :Class: `pymatgen.analysis.local_env`.

        :param structure: Structure object
        :param strategy: an instance of a
            :Class: `pymatgen.analysis.local_env.NearNeighbors` object
        :param weights: if True, use weights from local_env class
            (consult relevant class for their meaning)
        :return:
        """

        if not strategy.structures_allowed:
            raise ValueError(
                "Chosen strategy is not designed for use with structures! " "Please choose another strategy."
            )

        sg = StructureGraph.with_empty_graph(structure, name="bonds")

        for n, neighbors in enumerate(strategy.get_all_nn_info(structure)):
            for neighbor in neighbors:
                # local_env will always try to add two edges
                # for any one bond, one from site u to site v
                # and another form site v to site u: this is
                # harmless, so warn_duplicates=False
                sg.add_edge(
                    from_index=n,
                    from_jimage=(0, 0, 0),
                    to_index=neighbor["site_index"],
                    to_jimage=neighbor["image"],
                    weight=neighbor["weight"] if weights else None,
                    warn_duplicates=False,
                )

        return sg

    @property
    def name(self):
        """
        :return: Name of graph
        """
        return self.graph.graph["name"]

    @property
    def edge_weight_name(self):
        """
        :return: Name of the edge weight property of graph
        """
        return self.graph.graph["edge_weight_name"]

    @property
    def edge_weight_unit(self):
        """
        :return: Units of the edge weight property of graph
        """
        return self.graph.graph["edge_weight_units"]

    def add_edge(
        self,
        from_index,
        to_index,
        from_jimage=(0, 0, 0),
        to_jimage=None,
        weight=None,
        warn_duplicates=True,
        edge_properties=None,
    ):
        """
        Add edge to graph.

        Since physically a 'bond' (or other connection
        between sites) doesn't have a direction, from_index,
        from_jimage can be swapped with to_index, to_jimage.

        However, images will always always be shifted so that
        from_index < to_index and from_jimage becomes (0, 0, 0).

        :param from_index: index of site connecting from
        :param to_index: index of site connecting to
        :param from_jimage (tuple of ints): lattice vector of periodic
            image, e.g. (1, 0, 0) for periodic image in +x direction
        :param to_jimage (tuple of ints): lattice vector of image
        :param weight (float): e.g. bond length
        :param warn_duplicates (bool): if True, will warn if
            trying to add duplicate edges (duplicate edges will not
            be added in either case)
        :param edge_properties (dict): any other information to
            store on graph edges, similar to Structure's site_properties
        :return:
        """

        # this is not necessary for the class to work, but
        # just makes it neater
        if to_index < from_index:
            to_index, from_index = from_index, to_index
            to_jimage, from_jimage = from_jimage, to_jimage

        # constrain all from_jimages to be (0, 0, 0),
        # initial version of this class worked even if
        # from_jimage != (0, 0, 0), but making this
        # assumption simplifies logic later
        if not np.array_equal(from_jimage, (0, 0, 0)):
            shift = from_jimage
            from_jimage = np.subtract(from_jimage, shift)
            to_jimage = np.subtract(to_jimage, shift)

        # automatic detection of to_jimage if user doesn't specify
        # will try and detect all equivalent images and add multiple
        # edges if appropriate
        if to_jimage is None:
            # assume we want the closest site
            warnings.warn("Please specify to_jimage to be unambiguous, trying to automatically detect.")
            dist, to_jimage = self.structure[from_index].distance_and_image(self.structure[to_index])
            if dist == 0:
                # this will happen when from_index == to_index,
                # typically in primitive single-atom lattices
                images = [1, 0, 0], [0, 1, 0], [0, 0, 1]
                dists = []
                for image in images:
                    dists.append(
                        self.structure[from_index].distance_and_image(self.structure[from_index], jimage=image)[0]
                    )
                dist = min(dists)
            equiv_sites = self.structure.get_neighbors_in_shell(
                self.structure[from_index].coords, dist, dist * 0.01, include_index=True
            )
            for nnsite in equiv_sites:
                to_jimage = np.subtract(nnsite.frac_coords, self.structure[from_index].frac_coords)
                to_jimage = np.round(to_jimage).astype(int)
                self.add_edge(
                    from_index=from_index,
                    from_jimage=(0, 0, 0),
                    to_jimage=to_jimage,
                    to_index=nnsite.index,
                )
            return

        # sanitize types
        from_jimage, to_jimage = (
            tuple(map(int, from_jimage)),
            tuple(map(int, to_jimage)),
        )
        from_index, to_index = int(from_index), int(to_index)

        # if edge is from site i to site i, constrain direction of edge
        # this is a convention to avoid duplicate hops
        if to_index == from_index:
            if to_jimage == (0, 0, 0):
                warnings.warn("Tried to create a bond to itself, " "this doesn't make sense so was ignored.")
                return

            # ensure that the first non-zero jimage index is positive
            # assumes that at least one non-zero index is present
            is_positive = [idx for idx in to_jimage if idx != 0][0] > 0

            if not is_positive:
                # let's flip the jimage,
                # e.g. (0, 1, 0) is equivalent to (0, -1, 0) in this case
                to_jimage = tuple(-idx for idx in to_jimage)

        # check we're not trying to add a duplicate edge
        # there should only ever be at most one edge
        # between a given (site, jimage) pair and another
        # (site, jimage) pair
        existing_edge_data = self.graph.get_edge_data(from_index, to_index)
        if existing_edge_data:
            for key, d in existing_edge_data.items():
                if d["to_jimage"] == to_jimage:
                    if warn_duplicates:
                        warnings.warn(
                            "Trying to add an edge that already exists from "
                            "site {} to site {} in {}.".format(from_index, to_index, to_jimage)
                        )
                    return

        # generic container for additional edge properties,
        # similar to site properties
        edge_properties = edge_properties or {}

        if weight:
            self.graph.add_edge(from_index, to_index, to_jimage=to_jimage, weight=weight, **edge_properties)
        else:
            self.graph.add_edge(from_index, to_index, to_jimage=to_jimage, **edge_properties)

    def insert_node(
        self,
        i,
        species,
        coords,
        coords_are_cartesian=False,
        validate_proximity=False,
        site_properties=None,
        edges=None,
    ):
        """
        A wrapper around Molecule.insert(), which also incorporates the new
        site into the MoleculeGraph.

        :param i: Index at which to insert the new site
        :param species: Species for the new site
        :param coords: 3x1 array representing coordinates of the new site
        :param coords_are_cartesian: Whether coordinates are cartesian.
            Defaults to False.
        :param validate_proximity: For Molecule.insert(); if True (default
            False), distance will be checked to ensure that site can be safely
            added.
        :param site_properties: Site properties for Molecule
        :param edges: List of dicts representing edges to be added to the
            MoleculeGraph. These edges must include the index of the new site i,
            and all indices used for these edges should reflect the
            MoleculeGraph AFTER the insertion, NOT before. Each dict should at
            least have a "to_index" and "from_index" key, and can also have a
            "weight" and a "properties" key.
        :return:
        """

        self.structure.insert(
            i,
            species,
            coords,
            coords_are_cartesian=coords_are_cartesian,
            validate_proximity=validate_proximity,
            properties=site_properties,
        )

        mapping = {}
        for j in range(len(self.structure) - 1):
            if j < i:
                mapping[j] = j
            else:
                mapping[j] = j + 1
        nx.relabel_nodes(self.graph, mapping, copy=False)

        self.graph.add_node(i)
        self.set_node_attributes()

        if edges is not None:
            for edge in edges:
                try:
                    self.add_edge(
                        edge["from_index"],
                        edge["to_index"],
                        from_jimage=(0, 0, 0),
                        to_jimage=edge["to_jimage"],
                        weight=edge.get("weight", None),
                        edge_properties=edge.get("properties", None),
                    )
                except KeyError:
                    raise RuntimeError("Some edges are invalid.")

    def set_node_attributes(self):
        """
        Gives each node a "specie" and a "coords" attribute, updated with the
        current species and coordinates.

        :return:
        """

        species = {}
        coords = {}
        properties = {}
        for node in self.graph.nodes():
            species[node] = self.structure[node].specie.symbol
            coords[node] = self.structure[node].coords
            properties[node] = self.structure[node].properties

        nx.set_node_attributes(self.graph, species, "specie")
        nx.set_node_attributes(self.graph, coords, "coords")
        nx.set_node_attributes(self.graph, properties, "properties")

    def alter_edge(
        self,
        from_index,
        to_index,
        to_jimage=None,
        new_weight=None,
        new_edge_properties=None,
    ):
        """
        Alters either the weight or the edge_properties of
        an edge in the StructureGraph.

        :param from_index: int
        :param to_index: int
        :param to_jimage: tuple
        :param new_weight: alter_edge does not require
            that weight be altered. As such, by default, this
            is None. If weight is to be changed, it should be a
            float.
        :param new_edge_properties: alter_edge does not require
            that edge_properties be altered. As such, by default,
            this is None. If any edge properties are to be changed,
            it should be a dictionary of edge properties to be changed.
        :return:
        """

        existing_edges = self.graph.get_edge_data(from_index, to_index)

        # ensure that edge exists before attempting to change it
        if not existing_edges:
            raise ValueError(
                "Edge between {} and {} cannot be altered;\
                                no edge exists between those sites.".format(
                    from_index, to_index
                )
            )

        if to_jimage is None:
            edge_index = 0
        else:
            for i, properties in existing_edges.items():
                if properties["to_jimage"] == to_jimage:
                    edge_index = i

        if new_weight is not None:
            self.graph[from_index][to_index][edge_index]["weight"] = new_weight

        if new_edge_properties is not None:
            for prop in list(new_edge_properties.keys()):
                self.graph[from_index][to_index][edge_index][prop] = new_edge_properties[prop]

    def break_edge(self, from_index, to_index, to_jimage=None, allow_reverse=False):
        """
        Remove an edge from the StructureGraph. If no image is given, this method will fail.

        :param from_index: int
        :param to_index: int
        :param to_jimage: tuple
        :param allow_reverse: If allow_reverse is True, then break_edge will
            attempt to break both (from_index, to_index) and, failing that,
            will attempt to break (to_index, from_index).
        :return:
        """

        # ensure that edge exists before attempting to remove it
        existing_edges = self.graph.get_edge_data(from_index, to_index)
        existing_reverse = None

        if to_jimage is None:
            raise ValueError("Image must be supplied, to avoid ambiguity.")

        if existing_edges:
            for i, properties in existing_edges.items():
                if properties["to_jimage"] == to_jimage:
                    edge_index = i

            self.graph.remove_edge(from_index, to_index, edge_index)

        else:
            if allow_reverse:
                existing_reverse = self.graph.get_edge_data(to_index, from_index)

            if existing_reverse:
                for i, properties in existing_reverse.items():
                    if properties["to_jimage"] == to_jimage:
                        edge_index = i

                self.graph.remove_edge(to_index, from_index, edge_index)
            else:
                raise ValueError(
                    "Edge cannot be broken between {} and {};\
                                no edge exists between those sites.".format(
                        from_index, to_index
                    )
                )

    def remove_nodes(self, indices):
        """
        A wrapper for Molecule.remove_sites().

        :param indices: list of indices in the current Molecule (and graph) to
            be removed.
        :return:
        """

        self.structure.remove_sites(indices)
        self.graph.remove_nodes_from(indices)

        mapping = {}
        for correct, current in enumerate(sorted(self.graph.nodes)):
            mapping[current] = correct

        nx.relabel_nodes(self.graph, mapping, copy=False)
        self.set_node_attributes()

    def substitute_group(
        self,
        index,
        func_grp,
        strategy,
        bond_order=1,
        graph_dict=None,
        strategy_params=None,
    ):
        """
        Builds off of Structure.substitute to replace an atom in self.structure
        with a functional group. This method also amends self.graph to
        incorporate the new functional group.

        NOTE: Care must be taken to ensure that the functional group that is
        substituted will not place atoms to close to each other, or violate the
        dimensions of the Lattice.

        :param index: Index of atom to substitute.
        :param func_grp: Substituent molecule. There are two options:

            1. Providing an actual Molecule as the input. The first atom
                must be a DummySpecies X, indicating the position of
                nearest neighbor. The second atom must be the next
                nearest atom. For example, for a methyl group
                substitution, func_grp should be X-CH3, where X is the
                first site and C is the second site. What the code will
                do is to remove the index site, and connect the nearest
                neighbor to the C atom in CH3. The X-C bond indicates the
                directionality to connect the atoms.
            2. A string name. The molecule will be obtained from the
                relevant template in func_groups.json.
        :param strategy: Class from pymatgen.analysis.local_env.
        :param bond_order: A specified bond order to calculate the bond
            length between the attached functional group and the nearest
            neighbor site. Defaults to 1.
        :param graph_dict: Dictionary representing the bonds of the functional
            group (format: {(u, v): props}, where props is a dictionary of
            properties, including weight. If None, then the algorithm
            will attempt to automatically determine bonds using one of
            a list of strategies defined in pymatgen.analysis.local_env.
        :param strategy_params: dictionary of keyword arguments for strategy.
            If None, default parameters will be used.
        :return:
        """

        def map_indices(grp):
            grp_map = {}

            # Get indices now occupied by functional group
            # Subtracting 1 because the dummy atom X should not count
            atoms = len(grp) - 1
            offset = len(self.structure) - atoms

            for i in range(atoms):
                grp_map[i] = i + offset

            return grp_map

        if isinstance(func_grp, Molecule):
            func_grp = copy.deepcopy(func_grp)
        else:
            try:
                func_grp = copy.deepcopy(FunctionalGroups[func_grp])
            except Exception:
                raise RuntimeError("Can't find functional group in list. " "Provide explicit coordinate instead")

        self.structure.substitute(index, func_grp, bond_order=bond_order)

        mapping = map_indices(func_grp)

        # Remove dummy atom "X"
        func_grp.remove_species("X")

        if graph_dict is not None:
            for (u, v) in graph_dict.keys():
                edge_props = graph_dict[(u, v)]
                if "to_jimage" in edge_props.keys():
                    to_jimage = edge_props["to_jimage"]
                    del edge_props["to_jimage"]
                else:
                    # By default, assume that all edges should stay remain
                    # inside the initial image
                    to_jimage = (0, 0, 0)
                if "weight" in edge_props.keys():
                    weight = edge_props["weight"]
                    del edge_props["weight"]
                self.add_edge(
                    mapping[u],
                    mapping[v],
                    to_jimage=to_jimage,
                    weight=weight,
                    edge_properties=edge_props,
                )

        else:
            if strategy_params is None:
                strategy_params = {}
            strat = strategy(**strategy_params)

            for site in mapping.values():
                neighbors = strat.get_nn_info(self.structure, site)

                for neighbor in neighbors:
                    self.add_edge(
                        from_index=site,
                        from_jimage=(0, 0, 0),
                        to_index=neighbor["site_index"],
                        to_jimage=neighbor["image"],
                        weight=neighbor["weight"],
                        warn_duplicates=False,
                    )

    def get_connected_sites(self, n, jimage=(0, 0, 0)):
        """
        Returns a named tuple of neighbors of site n:
        periodic_site, jimage, index, weight.
        Index is the index of the corresponding site
        in the original structure, weight can be
        None if not defined.
        :param n: index of Site in Structure
        :param jimage: lattice vector of site
        :return: list of ConnectedSite tuples,
            sorted by closest first
        """

        connected_sites = set()
        connected_site_images = set()

        out_edges = [(u, v, d, "out") for u, v, d in self.graph.out_edges(n, data=True)]
        in_edges = [(u, v, d, "in") for u, v, d in self.graph.in_edges(n, data=True)]

        for u, v, d, dir in out_edges + in_edges:

            to_jimage = d["to_jimage"]

            if dir == "in":
                u, v = v, u
                to_jimage = np.multiply(-1, to_jimage)

            to_jimage = tuple(map(int, np.add(to_jimage, jimage)))
            site_d = self.structure[v].as_dict()
            site_d["abc"] = np.add(site_d["abc"], to_jimage).tolist()
            site = PeriodicSite.from_dict(site_d)

            # from_site if jimage arg != (0, 0, 0)
            relative_jimage = np.subtract(to_jimage, jimage)
            dist = self.structure[u].distance(self.structure[v], jimage=relative_jimage)

            weight = d.get("weight", None)

            if (v, to_jimage) not in connected_site_images:
                connected_site = ConnectedSite(site=site, jimage=to_jimage, index=v, weight=weight, dist=dist)

                connected_sites.add(connected_site)
                connected_site_images.add((v, to_jimage))

        # return list sorted by closest sites first
        connected_sites = list(connected_sites)
        connected_sites.sort(key=lambda x: x.dist)

        return connected_sites

    def get_coordination_of_site(self, n):
        """
        Returns the number of neighbors of site n.
        In graph terms, simply returns degree
        of node corresponding to site n.
        :param n: index of site
        :return (int):
        """
        number_of_self_loops = sum([1 for n, v in self.graph.edges(n) if n == v])
        return self.graph.degree(n) - number_of_self_loops

    def draw_graph_to_file(
        self,
        filename="graph",
        diff=None,
        hide_unconnected_nodes=False,
        hide_image_edges=True,
        edge_colors=False,
        node_labels=False,
        weight_labels=False,
        image_labels=False,
        color_scheme="VESTA",
        keep_dot=False,
        algo="fdp",
    ):
        """
        Draws graph using GraphViz.

        The networkx graph object itself can also be drawn
        with networkx's in-built graph drawing methods, but
        note that this might give misleading results for
        multigraphs (edges are super-imposed on each other).

        If visualization is difficult to interpret,
        `hide_image_edges` can help, especially in larger
        graphs.

        :param filename: filename to output, will detect filetype
            from extension (any graphviz filetype supported, such as
            pdf or png)
        :param diff (StructureGraph): an additional graph to
            compare with, will color edges red that do not exist in diff
            and edges green that are in diff graph but not in the
            reference graph
        :param hide_unconnected_nodes: if True, hide unconnected
            nodes
        :param hide_image_edges: if True, do not draw edges that
            go through periodic boundaries
        :param edge_colors (bool): if True, use node colors to
            color edges
        :param node_labels (bool): if True, label nodes with
            species and site index
        :param weight_labels (bool): if True, label edges with
            weights
        :param image_labels (bool): if True, label edges with
            their periodic images (usually only used for debugging,
            edges to periodic images always appear as dashed lines)
        :param color_scheme (str): "VESTA" or "JMOL"
        :param keep_dot (bool): keep GraphViz .dot file for later
            visualization
        :param algo: any graphviz algo, "neato" (for simple graphs)
            or "fdp" (for more crowded graphs) usually give good outputs
        :return:
        """

        if not which(algo):
            raise RuntimeError("StructureGraph graph drawing requires " "GraphViz binaries to be in the path.")

        # Developer note: NetworkX also has methods for drawing
        # graphs using matplotlib, these also work here. However,
        # a dedicated tool like GraphViz allows for much easier
        # control over graph appearance and also correctly displays
        # mutli-graphs (matplotlib can superimpose multiple edges).

        g = self.graph.copy()

        g.graph = {"nodesep": 10.0, "dpi": 300, "overlap": "false"}

        # add display options for nodes
        for n in g.nodes():
            # get label by species name
            label = "{}({})".format(str(self.structure[n].specie), n) if node_labels else ""

            # use standard color scheme for nodes
            c = EL_COLORS[color_scheme].get(str(self.structure[n].specie.symbol), [0, 0, 0])

            # get contrasting font color
            # magic numbers account for perceived luminescence
            # https://stackoverflow.com/questions/1855884/determine-font-color-based-on-background-color
            fontcolor = "#000000" if 1 - (c[0] * 0.299 + c[1] * 0.587 + c[2] * 0.114) / 255 < 0.5 else "#ffffff"

            # convert color to hex string
            color = "#{:02x}{:02x}{:02x}".format(c[0], c[1], c[2])

            g.add_node(
                n,
                fillcolor=color,
                fontcolor=fontcolor,
                label=label,
                fontname="Helvetica-bold",
                style="filled",
                shape="circle",
            )

        edges_to_delete = []

        # add display options for edges
        for u, v, k, d in g.edges(keys=True, data=True):

            # retrieve from/to images, set as origin if not defined
            to_image = d["to_jimage"]

            # set edge style
            d["style"] = "solid"
            if to_image != (0, 0, 0):
                d["style"] = "dashed"
                if hide_image_edges:
                    edges_to_delete.append((u, v, k))

            # don't show edge directions
            d["arrowhead"] = "none"

            # only add labels for images that are not the origin
            if image_labels:
                d["headlabel"] = "" if to_image == (0, 0, 0) else "to {}".format((to_image))
                d["arrowhead"] = "normal" if d["headlabel"] else "none"

            # optionally color edges using node colors
            color_u = g.nodes[u]["fillcolor"]
            color_v = g.nodes[v]["fillcolor"]
            d["color_uv"] = "{};0.5:{};0.5".format(color_u, color_v) if edge_colors else "#000000"

            # optionally add weights to graph
            if weight_labels:
                units = g.graph.get("edge_weight_units", "")
                if d.get("weight"):
                    d["label"] = "{:.2f} {}".format(d["weight"], units)

            # update edge with our new style attributes
            g.edges[u, v, k].update(d)

        # optionally remove periodic image edges,
        # these can be confusing due to periodic boundaries
        if hide_image_edges:
            for edge_to_delete in edges_to_delete:
                g.remove_edge(*edge_to_delete)

        # optionally hide unconnected nodes,
        # these can appear when removing periodic edges
        if hide_unconnected_nodes:
            g = g.subgraph([n for n in g.degree() if g.degree()[n] != 0])

        # optionally highlight differences with another graph
        if diff:
            diff = self.diff(diff, strict=True)
            green_edges = []
            red_edges = []
            for u, v, k, d in g.edges(keys=True, data=True):
                if (u, v, d["to_jimage"]) in diff["self"]:
                    # edge has been deleted
                    red_edges.append((u, v, k))
                elif (u, v, d["to_jimage"]) in diff["other"]:
                    # edge has been added
                    green_edges.append((u, v, k))
            for u, v, k in green_edges:
                g.edges[u, v, k].update({"color_uv": "#00ff00"})
            for u, v, k in red_edges:
                g.edges[u, v, k].update({"color_uv": "#ff0000"})

        basename, extension = os.path.splitext(filename)
        extension = extension[1:]

        write_dot(g, basename + ".dot")

        with open(filename, "w") as f:

            args = [algo, "-T", extension, basename + ".dot"]
            rs = subprocess.Popen(args, stdout=f, stdin=subprocess.PIPE, close_fds=True)
            rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError("{} exited with return code {}.".format(algo, rs.returncode))

        if not keep_dot:
            os.remove(basename + ".dot")

    @property
    def types_and_weights_of_connections(self):
        """
        Extract a dictionary summarizing the types and weights
        of edges in the graph.

        :return: A dictionary with keys specifying the
            species involved in a connection in alphabetical order
            (e.g. string 'Fe-O') and values which are a list of
            weights for those connections (e.g. bond lengths).
        """

        def get_label(u, v):
            u_label = self.structure[u].species_string
            v_label = self.structure[v].species_string
            return "-".join(sorted((u_label, v_label)))

        types = defaultdict(list)
        for u, v, d in self.graph.edges(data=True):
            label = get_label(u, v)
            types[label].append(d["weight"])

        return dict(types)

    @property
    def weight_statistics(self):
        """
        Extract a statistical summary of edge weights present in
        the graph.

        :return: A dict with an 'all_weights' list, 'minimum',
            'maximum', 'median', 'mean', 'std_dev'
        """

        all_weights = [d.get("weight", None) for u, v, d in self.graph.edges(data=True)]
        stats = describe(all_weights, nan_policy="omit")

        return {
            "all_weights": all_weights,
            "min": stats.minmax[0],
            "max": stats.minmax[1],
            "mean": stats.mean,
            "variance": stats.variance,
        }

    def types_of_coordination_environments(self, anonymous=False):
        """
        Extract information on the different co-ordination environments
        present in the graph.

        :param anonymous: if anonymous, will replace specie names
            with A, B, C, etc.
        :return: a list of co-ordination environments,
            e.g. ['Mo-S(6)', 'S-Mo(3)']
        """

        motifs = set()
        for idx, site in enumerate(self.structure):

            centre_sp = site.species_string

            connected_sites = self.get_connected_sites(idx)
            connected_species = [connected_site.site.species_string for connected_site in connected_sites]

            labels = []
            for sp in set(connected_species):
                count = connected_species.count(sp)
                labels.append((count, sp))

            labels = sorted(labels, reverse=True)

            if anonymous:
                mapping = {centre_sp: "A"}
                available_letters = [chr(66 + i) for i in range(25)]
                for label in labels:
                    sp = label[1]
                    if sp not in mapping:
                        mapping[sp] = available_letters.pop(0)
                centre_sp = "A"
                labels = [(label[0], mapping[label[1]]) for label in labels]

            labels = ["{}({})".format(label[1], label[0]) for label in labels]
            motif = "{}-{}".format(centre_sp, ",".join(labels))
            motifs.add(motif)

        return sorted(list(motifs))

    def as_dict(self):
        """
        As in :Class: `pymatgen.core.Structure` except
        with using `to_dict_of_dicts` from NetworkX
        to store graph information.
        """

        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "structure": self.structure.as_dict(),
            "graphs": json_graph.adjacency_data(self.graph),
        }

        return d

    @classmethod
    def from_dict(cls, d):
        """
        As in :Class: `pymatgen.core.Structure` except
        restoring graphs using `from_dict_of_dicts`
        from NetworkX to restore graph information.
        """
        s = Structure.from_dict(d["structure"])
        return cls(s, d["graphs"])

    def __mul__(self, scaling_matrix):
        """
        Replicates the graph, creating a supercell,
        intelligently joining together
        edges that lie on periodic boundaries.
        In principle, any operations on the expanded
        graph could also be done on the original
        graph, but a larger graph can be easier to
        visualize and reason about.
        :param scaling_matrix: same as Structure.__mul__
        :return:
        """

        # Developer note: a different approach was also trialed, using
        # a simple Graph (instead of MultiDiGraph), with node indices
        # representing both site index and periodic image. Here, the
        # number of nodes != number of sites in the Structure. This
        # approach has many benefits, but made it more difficult to
        # keep the graph in sync with its corresponding Structure.

        # Broadly, it would be easier to multiply the Structure
        # *before* generating the StructureGraph, but this isn't
        # possible when generating the graph using critic2 from
        # charge density.

        # Multiplication works by looking for the expected position
        # of an image node, and seeing if that node exists in the
        # supercell. If it does, the edge is updated. This is more
        # computationally expensive than just keeping track of the
        # which new lattice images present, but should hopefully be
        # easier to extend to a general 3x3 scaling matrix.

        # code adapted from Structure.__mul__
        scale_matrix = np.array(scaling_matrix, np.int16)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
        else:
            # TODO: test __mul__ with full 3x3 scaling matrices
            raise NotImplementedError("Not tested with 3x3 scaling matrices yet.")
        new_lattice = Lattice(np.dot(scale_matrix, self.structure.lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        new_graphs = []

        for v in c_lat:

            # create a map of nodes from original graph to its image
            mapping = {n: n + len(new_sites) for n in range(len(self.structure))}

            for idx, site in enumerate(self.structure):
                s = PeriodicSite(
                    site.species,
                    site.coords + v,
                    new_lattice,
                    properties=site.properties,
                    coords_are_cartesian=True,
                    to_unit_cell=False,
                )

                new_sites.append(s)

            new_graphs.append(nx.relabel_nodes(self.graph, mapping, copy=True))

        new_structure = Structure.from_sites(new_sites)

        # merge all graphs into one big graph
        new_g = nx.MultiDiGraph()
        for new_graph in new_graphs:
            new_g = nx.union(new_g, new_graph)

        edges_to_remove = []  # tuple of (u, v, k)
        edges_to_add = []  # tuple of (u, v, attr_dict)

        # list of new edges inside supercell
        # for duplicate checking
        edges_inside_supercell = [{u, v} for u, v, d in new_g.edges(data=True) if d["to_jimage"] == (0, 0, 0)]
        new_periodic_images = []

        orig_lattice = self.structure.lattice

        # use k-d tree to match given position to an
        # existing Site in Structure
        kd_tree = KDTree(new_structure.cart_coords)

        # tolerance in Å for sites to be considered equal
        # this could probably be a lot smaller
        tol = 0.05

        for u, v, k, d in new_g.edges(keys=True, data=True):

            to_jimage = d["to_jimage"]  # for node v

            # reduce unnecessary checking
            if to_jimage != (0, 0, 0):

                # get index in original site
                n_u = u % len(self.structure)
                n_v = v % len(self.structure)

                # get fractional co-ordinates of where atoms defined
                # by edge are expected to be, relative to original
                # lattice (keeping original lattice has
                # significant benefits)
                v_image_frac = np.add(self.structure[n_v].frac_coords, to_jimage)
                u_frac = self.structure[n_u].frac_coords

                # using the position of node u as a reference,
                # get relative Cartesian co-ordinates of where
                # atoms defined by edge are expected to be
                v_image_cart = orig_lattice.get_cartesian_coords(v_image_frac)
                u_cart = orig_lattice.get_cartesian_coords(u_frac)
                v_rel = np.subtract(v_image_cart, u_cart)

                # now retrieve position of node v in
                # new supercell, and get asgolute Cartesian
                # co-ordinates of where atoms defined by edge
                # are expected to be
                v_expec = new_structure[u].coords + v_rel

                # now search in new structure for these atoms
                # query returns (distance, index)
                v_present = kd_tree.query(v_expec)
                v_present = v_present[1] if v_present[0] <= tol else None

                # check if image sites now present in supercell
                # and if so, delete old edge that went through
                # periodic boundary
                if v_present is not None:

                    new_u = u
                    new_v = v_present
                    new_d = d.copy()

                    # node now inside supercell
                    new_d["to_jimage"] = (0, 0, 0)

                    edges_to_remove.append((u, v, k))

                    # make sure we don't try to add duplicate edges
                    # will remove two edges for everyone one we add
                    if {new_u, new_v} not in edges_inside_supercell:

                        # normalize direction
                        if new_v < new_u:
                            new_u, new_v = new_v, new_u

                        edges_inside_supercell.append({new_u, new_v})
                        edges_to_add.append((new_u, new_v, new_d))

                else:

                    # want to find new_v such that we have
                    # full periodic boundary conditions
                    # so that nodes on one side of supercell
                    # are connected to nodes on opposite side

                    v_expec_frac = new_structure.lattice.get_fractional_coords(v_expec)

                    # find new to_jimage
                    # use np.around to fix issues with finite precision leading to incorrect image
                    v_expec_image = np.around(v_expec_frac, decimals=3)
                    v_expec_image = v_expec_image - v_expec_image % 1

                    v_expec_frac = np.subtract(v_expec_frac, v_expec_image)
                    v_expec = new_structure.lattice.get_cartesian_coords(v_expec_frac)
                    v_present = kd_tree.query(v_expec)
                    v_present = v_present[1] if v_present[0] <= tol else None

                    if v_present is not None:

                        new_u = u
                        new_v = v_present
                        new_d = d.copy()
                        new_to_jimage = tuple(map(int, v_expec_image))

                        # normalize direction
                        if new_v < new_u:
                            new_u, new_v = new_v, new_u
                            new_to_jimage = tuple(np.multiply(-1, d["to_jimage"]).astype(int))

                        new_d["to_jimage"] = new_to_jimage

                        edges_to_remove.append((u, v, k))

                        if (new_u, new_v, new_to_jimage) not in new_periodic_images:
                            edges_to_add.append((new_u, new_v, new_d))
                            new_periodic_images.append((new_u, new_v, new_to_jimage))

        logger.debug("Removing {} edges, adding {} new edges.".format(len(edges_to_remove), len(edges_to_add)))

        # add/delete marked edges
        for edges_to_remove in edges_to_remove:
            new_g.remove_edge(*edges_to_remove)
        for (u, v, d) in edges_to_add:
            new_g.add_edge(u, v, **d)

        # return new instance of StructureGraph with supercell
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "structure": new_structure.as_dict(),
            "graphs": json_graph.adjacency_data(new_g),
        }

        sg = StructureGraph.from_dict(d)

        return sg

    def __rmul__(self, other):
        return self.__mul__(other)

    @classmethod
    def _edges_to_string(cls, g):

        header = "from    to  to_image    "
        header_line = "----  ----  ------------"
        edge_weight_name = g.graph["edge_weight_name"]
        if edge_weight_name:
            print_weights = ["weight"]
            edge_label = g.graph["edge_weight_name"]
            edge_weight_units = g.graph["edge_weight_units"]
            if edge_weight_units:
                edge_label += " ({})".format(edge_weight_units)
            header += "  {}".format(edge_label)
            header_line += "  {}".format("-" * max([18, len(edge_label)]))
        else:
            print_weights = False

        s = header + "\n" + header_line + "\n"

        edges = list(g.edges(data=True))

        # sort edges for consistent ordering
        edges.sort(key=itemgetter(0, 1))

        if print_weights:
            for u, v, data in edges:
                s += "{:4}  {:4}  {:12}  {:.3e}\n".format(
                    u, v, str(data.get("to_jimage", (0, 0, 0))), data.get("weight", 0)
                )
        else:
            for u, v, data in edges:
                s += "{:4}  {:4}  {:12}\n".format(u, v, str(data.get("to_jimage", (0, 0, 0))))

        return s

    def __str__(self):
        s = "Structure Graph"
        s += "\nStructure: \n{}".format(self.structure.__str__())
        s += "\nGraph: {}\n".format(self.name)
        s += self._edges_to_string(self.graph)
        return s

    def __repr__(self):
        s = "Structure Graph"
        s += "\nStructure: \n{}".format(self.structure.__repr__())
        s += "\nGraph: {}\n".format(self.name)
        s += self._edges_to_string(self.graph)
        return s

    def __len__(self):
        """
        :return: length of Structure / number of nodes in graph
        """
        return len(self.structure)

    def sort(self, key=None, reverse=False):
        """
        Same as Structure.sort(), also remaps nodes in graph.
        :param key:
        :param reverse:
        :return:
        """

        old_structure = self.structure.copy()

        # sort Structure
        self.structure._sites = sorted(self.structure._sites, key=key, reverse=reverse)

        # apply Structure ordering to graph
        mapping = {idx: self.structure.index(site) for idx, site in enumerate(old_structure)}
        self.graph = nx.relabel_nodes(self.graph, mapping, copy=True)

        # normalize directions of edges
        edges_to_remove = []
        edges_to_add = []
        for u, v, k, d in self.graph.edges(keys=True, data=True):
            if v < u:
                new_v, new_u, new_d = u, v, d.copy()
                new_d["to_jimage"] = tuple(np.multiply(-1, d["to_jimage"]).astype(int))
                edges_to_remove.append((u, v, k))
                edges_to_add.append((new_u, new_v, new_d))

        # add/delete marked edges
        for edges_to_remove in edges_to_remove:
            self.graph.remove_edge(*edges_to_remove)
        for (u, v, d) in edges_to_add:
            self.graph.add_edge(u, v, **d)

    def __copy__(self):
        return StructureGraph.from_dict(self.as_dict())

    def __eq__(self, other):
        """
        Two StructureGraphs are equal if they have equal Structures,
        and have the same edges between Sites. Edge weights can be
        different and StructureGraphs can still be considered equal.

        :param other: StructureGraph
        :return (bool):
        """

        # sort for consistent node indices
        # PeriodicSite should have a proper __hash__() value,
        # using its frac_coords as a convenient key
        mapping = {tuple(site.frac_coords): self.structure.index(site) for site in other.structure}
        other_sorted = other.__copy__()
        other_sorted.sort(key=lambda site: mapping[tuple(site.frac_coords)])

        edges = {(u, v, d["to_jimage"]) for u, v, d in self.graph.edges(keys=False, data=True)}

        edges_other = {(u, v, d["to_jimage"]) for u, v, d in other_sorted.graph.edges(keys=False, data=True)}

        return (edges == edges_other) and (self.structure == other_sorted.structure)

    def diff(self, other, strict=True):
        """
        Compares two StructureGraphs. Returns dict with
        keys 'self', 'other', 'both' with edges that are
        present in only one StructureGraph ('self' and
        'other'), and edges that are present in both.

        The Jaccard distance is a simple measure of the
        dissimilarity between two StructureGraphs (ignoring
        edge weights), and is defined by 1 - (size of the
        intersection / size of the union) of the sets of
        edges. This is returned with key 'dist'.

        Important note: all node indices are in terms
        of the StructureGraph this method is called
        from, not the 'other' StructureGraph: there
        is no guarantee the node indices will be the
        same if the underlying Structures are ordered
        differently.

        :param other: StructureGraph
        :param strict: if False, will compare bonds
            from different Structures, with node indices
            replaced by Species strings, will not count
            number of occurrences of bonds
        :return:
        """

        if self.structure != other.structure and strict:
            return ValueError("Meaningless to compare StructureGraphs if " "corresponding Structures are different.")

        if strict:

            # sort for consistent node indices
            # PeriodicSite should have a proper __hash__() value,
            # using its frac_coords as a convenient key
            mapping = {tuple(site.frac_coords): self.structure.index(site) for site in other.structure}
            other_sorted = other.__copy__()
            other_sorted.sort(key=lambda site: mapping[tuple(site.frac_coords)])

            edges = {(u, v, d["to_jimage"]) for u, v, d in self.graph.edges(keys=False, data=True)}

            edges_other = {(u, v, d["to_jimage"]) for u, v, d in other_sorted.graph.edges(keys=False, data=True)}

        else:

            edges = {
                (str(self.structure[u].specie), str(self.structure[v].specie))
                for u, v, d in self.graph.edges(keys=False, data=True)
            }

            edges_other = {
                (str(other.structure[u].specie), str(other.structure[v].specie))
                for u, v, d in other.graph.edges(keys=False, data=True)
            }

        if len(edges) == 0 and len(edges_other) == 0:
            jaccard_dist = 0  # by definition
        else:
            jaccard_dist = 1 - len(edges.intersection(edges_other)) / len(edges.union(edges_other))

        return {
            "self": edges - edges_other,
            "other": edges_other - edges,
            "both": edges.intersection(edges_other),
            "dist": jaccard_dist,
        }

    def get_subgraphs_as_molecules(self, use_weights=False):
        """
        Retrieve subgraphs as molecules, useful for extracting
        molecules from periodic crystals.

        Will only return unique molecules, not any duplicates
        present in the crystal (a duplicate defined as an
        isomorphic subgraph).

        :param use_weights (bool): If True, only treat subgraphs
            as isomorphic if edges have the same weights. Typically,
            this means molecules will need to have the same bond
            lengths to be defined as duplicates, otherwise bond
            lengths can differ. This is a fairly robust approach,
            but will treat e.g. enantiomers as being duplicates.

        :return: list of unique Molecules in Structure
        """

        # creating a supercell is an easy way to extract
        # molecules (and not, e.g., layers of a 2D crystal)
        # without adding extra logic
        if getattr(self, "_supercell_sg", None) is None:
            self._supercell_sg = supercell_sg = self * (3, 3, 3)

        # make undirected to find connected subgraphs
        supercell_sg.graph = nx.Graph(supercell_sg.graph)

        # find subgraphs
        all_subgraphs = [supercell_sg.graph.subgraph(c) for c in nx.connected_components(supercell_sg.graph)]

        # discount subgraphs that lie across *supercell* boundaries
        # these will subgraphs representing crystals
        molecule_subgraphs = []
        for subgraph in all_subgraphs:
            intersects_boundary = any(d["to_jimage"] != (0, 0, 0) for u, v, d in subgraph.edges(data=True))
            if not intersects_boundary:
                molecule_subgraphs.append(nx.MultiDiGraph(subgraph))

        # add specie names to graph to be able to test for isomorphism
        for subgraph in molecule_subgraphs:
            for n in subgraph:
                subgraph.add_node(n, specie=str(supercell_sg.structure[n].specie))

        # now define how we test for isomorphism
        def node_match(n1, n2):
            return n1["specie"] == n2["specie"]

        def edge_match(e1, e2):
            if use_weights:
                return e1["weight"] == e2["weight"]
            return True

        # prune duplicate subgraphs
        unique_subgraphs = []
        for subgraph in molecule_subgraphs:

            already_present = [
                nx.is_isomorphic(subgraph, g, node_match=node_match, edge_match=edge_match) for g in unique_subgraphs
            ]

            if not any(already_present):
                unique_subgraphs.append(subgraph)

        # get Molecule objects for each subgraph
        molecules = []
        for subgraph in unique_subgraphs:
            coords = [supercell_sg.structure[n].coords for n in subgraph.nodes()]
            species = [supercell_sg.structure[n].specie for n in subgraph.nodes()]

            molecule = Molecule(species, coords)

            # shift so origin is at center of mass
            molecule = molecule.get_centered_molecule()

            molecules.append(molecule)

        return molecules


class MolGraphSplitError(Exception):
    """
    Raised when a molecule graph is failed to split into two disconnected
    subgraphs
    """

    pass


class MoleculeGraph(MSONable):
    """
    This is a class for annotating a Molecule with
    bond information, stored in the form of a graph. A "bond" does
    not necessarily have to be a chemical bond, but can store any
    kind of information that connects two Sites.
    """

    def __init__(self, molecule, graph_data=None):
        """
        If constructing this class manually, use the `with_empty_graph`
        method or `with_local_env_strategy` method (using an algorithm
        provided by the `local_env` module, such as O'Keeffe).

        This class that contains connection information:
        relationships between sites represented by a Graph structure,
        and an associated structure object.

        This class uses the NetworkX package to store and operate
        on the graph itself, but contains a lot of helper methods
        to make associating a graph with a given molecule easier.

        Use cases for this include storing bonding information,
        NMR J-couplings, Heisenberg exchange parameters, etc.

        :param molecule: Molecule object

        :param graph_data: dict containing graph information in
            dict format (not intended to be constructed manually,
            see as_dict method for format)
        """

        if isinstance(molecule, MoleculeGraph):
            # just make a copy from input
            graph_data = molecule.as_dict()["graphs"]

        self.molecule = molecule
        self.graph = nx.readwrite.json_graph.adjacency_graph(graph_data)

        # tidy up edge attr dicts, reading to/from json duplicates
        # information
        for u, v, k, d in self.graph.edges(keys=True, data=True):
            if "id" in d:
                del d["id"]
            if "key" in d:
                del d["key"]
            # ensure images are tuples (conversion to lists happens
            # when serializing back from json), it's important images
            # are hashable/immutable
            if "to_jimage" in d:
                d["to_jimage"] = tuple(d["to_jimage"])
            if "from_jimage" in d:
                d["from_jimage"] = tuple(d["from_jimage"])

        self.set_node_attributes()

    @classmethod
    def with_empty_graph(cls, molecule, name="bonds", edge_weight_name=None, edge_weight_units=None):
        """
        Constructor for MoleculeGraph, returns a MoleculeGraph
        object with an empty graph (no edges, only nodes defined
        that correspond to Sites in Molecule).

        :param molecule (Molecule):
        :param name (str): name of graph, e.g. "bonds"
        :param edge_weight_name (str): name of edge weights,
            e.g. "bond_length" or "exchange_constant"
        :param edge_weight_units (str): name of edge weight units
            e.g. "Å" or "eV"
        :return (MoleculeGraph):
        """

        if edge_weight_name and (edge_weight_units is None):
            raise ValueError(
                "Please specify units associated "
                "with your edge weights. Can be "
                "empty string if arbitrary or "
                "dimensionless."
            )

        # construct graph with one node per site
        # graph attributes don't change behavior of graph,
        # they're just for book-keeping
        graph = nx.MultiDiGraph(
            edge_weight_name=edge_weight_name,
            edge_weight_units=edge_weight_units,
            name=name,
        )
        graph.add_nodes_from(range(len(molecule)))

        graph_data = json_graph.adjacency_data(graph)

        return cls(molecule, graph_data=graph_data)

    @staticmethod
    def with_edges(molecule, edges):
        """
        Constructor for MoleculeGraph, using pre-existing or pre-defined edges
        with optional edge parameters.

        :param molecule: Molecule object
        :param edges: dict representing the bonds of the functional
            group (format: {(u, v): props}, where props is a dictionary of
            properties, including weight. Props should be None if no
            additional properties are to be specified.
        :return: mg, a MoleculeGraph
        """

        mg = MoleculeGraph.with_empty_graph(molecule, name="bonds", edge_weight_name="weight", edge_weight_units="")

        for edge, props in edges.items():

            try:
                from_index = edge[0]
                to_index = edge[1]
            except TypeError:
                raise ValueError("Edges must be given as (from_index, to_index)" "tuples")

            if props is not None:
                if "weight" in props.keys():
                    weight = props["weight"]
                    del props["weight"]
                else:
                    weight = None

                if len(props.items()) == 0:
                    props = None
            else:
                weight = None

            nodes = mg.graph.nodes
            if not (from_index in nodes and to_index in nodes):
                raise ValueError(
                    "Edges cannot be added if nodes are not" " present in the graph. Please check your" " indices."
                )

            mg.add_edge(from_index, to_index, weight=weight, edge_properties=props)

        mg.set_node_attributes()
        return mg

    @staticmethod
    def with_local_env_strategy(molecule, strategy):
        """
        Constructor for MoleculeGraph, using a strategy
        from :Class: `pymatgen.analysis.local_env`.

        :param molecule: Molecule object
        :param strategy: an instance of a
            :Class: `pymatgen.analysis.local_env.NearNeighbors` object
        :return: mg, a MoleculeGraph
        """

        if not strategy.molecules_allowed:
            raise ValueError(
                "Chosen strategy is not designed for use with molecules! " "Please choose another strategy."
            )
        extend_structure = strategy.extend_structure_molecules

        mg = MoleculeGraph.with_empty_graph(molecule, name="bonds", edge_weight_name="weight", edge_weight_units="")

        # NearNeighbor classes only (generally) work with structures
        # molecules have to be boxed first
        coords = molecule.cart_coords

        if extend_structure:
            a = max(coords[:, 0]) - min(coords[:, 0]) + 100
            b = max(coords[:, 1]) - min(coords[:, 1]) + 100
            c = max(coords[:, 2]) - min(coords[:, 2]) + 100

            structure = molecule.get_boxed_structure(a, b, c, no_cross=True, reorder=False)
        else:
            structure = None

        for n in range(len(molecule)):

            if structure is None:
                neighbors = strategy.get_nn_info(molecule, n)
            else:
                neighbors = strategy.get_nn_info(structure, n)
            for neighbor in neighbors:

                # all bonds in molecules should not cross
                # (artificial) periodic boundaries
                if not np.array_equal(neighbor["image"], [0, 0, 0]):
                    continue

                if n > neighbor["site_index"]:
                    from_index = neighbor["site_index"]
                    to_index = n
                else:
                    from_index = n
                    to_index = neighbor["site_index"]

                mg.add_edge(
                    from_index=from_index,
                    to_index=to_index,
                    weight=neighbor["weight"],
                    warn_duplicates=False,
                )

        duplicates = []
        for edge in mg.graph.edges:
            if edge[2] != 0:
                duplicates.append(edge)

        for duplicate in duplicates:
            mg.graph.remove_edge(duplicate[0], duplicate[1], key=duplicate[2])

        mg.set_node_attributes()
        return mg

    @property
    def name(self):
        """
        :return: Name of graph
        """
        return self.graph.graph["name"]

    @property
    def edge_weight_name(self):
        """
        :return: Name of the edge weight property of graph
        """
        return self.graph.graph["edge_weight_name"]

    @property
    def edge_weight_unit(self):
        """
        :return: Units of the edge weight property of graph
        """
        return self.graph.graph["edge_weight_units"]

    def add_edge(
        self,
        from_index,
        to_index,
        weight=None,
        warn_duplicates=True,
        edge_properties=None,
    ):
        """
        Add edge to graph.

        Since physically a 'bond' (or other connection
        between sites) doesn't have a direction, from_index,
        from_jimage can be swapped with to_index, to_jimage.

        However, images will always always be shifted so that
        from_index < to_index and from_jimage becomes (0, 0, 0).

        :param from_index: index of site connecting from
        :param to_index: index of site connecting to
        :param weight (float): e.g. bond length
        :param warn_duplicates (bool): if True, will warn if
            trying to add duplicate edges (duplicate edges will not
            be added in either case)
        :param edge_properties (dict): any other information to
            store on graph edges, similar to Structure's site_properties
        :return:
        """

        # this is not necessary for the class to work, but
        # just makes it neater
        if to_index < from_index:
            to_index, from_index = from_index, to_index

        # sanitize types
        from_index, to_index = int(from_index), int(to_index)

        # check we're not trying to add a duplicate edge
        # there should only ever be at most one edge
        # between two sites
        existing_edge_data = self.graph.get_edge_data(from_index, to_index)
        if existing_edge_data and warn_duplicates:
            warnings.warn(
                "Trying to add an edge that already exists from " "site {} to site {}.".format(from_index, to_index)
            )
            return

        # generic container for additional edge properties,
        # similar to site properties
        edge_properties = edge_properties or {}

        if weight:
            self.graph.add_edge(from_index, to_index, weight=weight, **edge_properties)
        else:
            self.graph.add_edge(from_index, to_index, **edge_properties)

    def insert_node(
        self,
        i,
        species,
        coords,
        validate_proximity=False,
        site_properties=None,
        edges=None,
    ):
        """
        A wrapper around Molecule.insert(), which also incorporates the new
        site into the MoleculeGraph.

        :param i: Index at which to insert the new site
        :param species: Species for the new site
        :param coords: 3x1 array representing coordinates of the new site
        :param validate_proximity: For Molecule.insert(); if True (default
            False), distance will be checked to ensure that site can be safely
            added.
        :param site_properties: Site properties for Molecule
        :param edges: List of dicts representing edges to be added to the
            MoleculeGraph. These edges must include the index of the new site i,
            and all indices used for these edges should reflect the
            MoleculeGraph AFTER the insertion, NOT before. Each dict should at
            least have a "to_index" and "from_index" key, and can also have a
            "weight" and a "properties" key.
        :return:
        """

        self.molecule.insert(
            i,
            species,
            coords,
            validate_proximity=validate_proximity,
            properties=site_properties,
        )

        mapping = {}
        for j in range(len(self.molecule) - 1):
            if j < i:
                mapping[j] = j
            else:
                mapping[j] = j + 1
        nx.relabel_nodes(self.graph, mapping, copy=False)

        self.graph.add_node(i)
        self.set_node_attributes()

        if edges is not None:
            for edge in edges:
                try:
                    self.add_edge(
                        edge["from_index"],
                        edge["to_index"],
                        weight=edge.get("weight", None),
                        edge_properties=edge.get("properties", None),
                    )
                except KeyError:
                    raise RuntimeError("Some edges are invalid.")

    def set_node_attributes(self):
        """
        Replicates molecule site properties (specie, coords, etc.) in the
        MoleculeGraph.

        :return:
        """

        species = {}
        coords = {}
        properties = {}
        for node in self.graph.nodes():
            species[node] = self.molecule[node].specie.symbol
            coords[node] = self.molecule[node].coords
            properties[node] = self.molecule[node].properties

        nx.set_node_attributes(self.graph, species, "specie")
        nx.set_node_attributes(self.graph, coords, "coords")
        nx.set_node_attributes(self.graph, properties, "properties")

    def alter_edge(self, from_index, to_index, new_weight=None, new_edge_properties=None):
        """
        Alters either the weight or the edge_properties of
        an edge in the MoleculeGraph.

        :param from_index: int
        :param to_index: int
        :param new_weight: alter_edge does not require
            that weight be altered. As such, by default, this
            is None. If weight is to be changed, it should be a
            float.
        :param new_edge_properties: alter_edge does not require
            that edge_properties be altered. As such, by default,
            this is None. If any edge properties are to be changed,
            it should be a dictionary of edge properties to be changed.
        :return:
        """

        existing_edge = self.graph.get_edge_data(from_index, to_index)

        # ensure that edge exists before attempting to change it
        if not existing_edge:
            raise ValueError(
                "Edge between {} and {} cannot be altered;\
                                no edge exists between those sites.".format(
                    from_index, to_index
                )
            )

        # Third index should always be 0 because there should only be one edge between any two nodes
        if new_weight is not None:
            self.graph[from_index][to_index][0]["weight"] = new_weight

        if new_edge_properties is not None:
            for prop in list(new_edge_properties.keys()):
                self.graph[from_index][to_index][0][prop] = new_edge_properties[prop]

    def break_edge(self, from_index, to_index, allow_reverse=False):
        """
        Remove an edge from the MoleculeGraph

        :param from_index: int
        :param to_index: int
        :param allow_reverse: If allow_reverse is True, then break_edge will
            attempt to break both (from_index, to_index) and, failing that,
            will attempt to break (to_index, from_index).
        :return:
        """

        # ensure that edge exists before attempting to remove it
        existing_edge = self.graph.get_edge_data(from_index, to_index)
        existing_reverse = None

        if existing_edge:
            self.graph.remove_edge(from_index, to_index)

        else:
            if allow_reverse:
                existing_reverse = self.graph.get_edge_data(to_index, from_index)

            if existing_reverse:
                self.graph.remove_edge(to_index, from_index)
            else:
                raise ValueError(
                    "Edge cannot be broken between {} and {};\
                                no edge exists between those sites.".format(
                        from_index, to_index
                    )
                )

    def remove_nodes(self, indices):
        """
        A wrapper for Molecule.remove_sites().

        :param indices: list of indices in the current Molecule (and graph) to
            be removed.
        :return:
        """

        self.molecule.remove_sites(indices)
        self.graph.remove_nodes_from(indices)

        mapping = {}
        for correct, current in enumerate(sorted(self.graph.nodes)):
            mapping[current] = correct

        nx.relabel_nodes(self.graph, mapping, copy=False)
        self.set_node_attributes()

    def get_disconnected_fragments(self):
        """
        Determine if the MoleculeGraph is connected. If it is not, separate the
        MoleculeGraph into different MoleculeGraphs, where each resulting
        MoleculeGraph is a disconnected subgraph of the original.
        Currently, this function naively assigns the charge
        of the total molecule to a single submolecule. A
        later effort will be to actually accurately assign
        charge.
        NOTE: This function does not modify the original
        MoleculeGraph. It creates a copy, modifies that, and
        returns two or more new MoleculeGraph objects.
        :return: list of MoleculeGraphs
        """

        if nx.is_weakly_connected(self.graph):
            return [copy.deepcopy(self)]

        original = copy.deepcopy(self)
        sub_mols = list()

        # Had to use nx.weakly_connected_components because of deprecation
        # of nx.weakly_connected_component_subgraphs
        subgraphs = [original.graph.subgraph(c) for c in nx.weakly_connected_components(original.graph)]

        for subg in subgraphs:

            nodes = sorted(list(subg.nodes))

            # Molecule indices are essentially list-based, so node indices
            # must be remapped, incrementing from 0
            mapping = {}
            for i, n in enumerate(nodes):
                mapping[n] = i

            # just give charge to whatever subgraph has node with index 0
            # TODO: actually figure out how to distribute charge
            if 0 in nodes:
                charge = self.molecule.charge
            else:
                charge = 0

            # relabel nodes in graph to match mapping
            new_graph = nx.relabel_nodes(subg, mapping)

            species = nx.get_node_attributes(new_graph, "specie")
            coords = nx.get_node_attributes(new_graph, "coords")
            raw_props = nx.get_node_attributes(new_graph, "properties")

            properties = {}
            for prop_set in raw_props.values():
                for prop in prop_set.keys():
                    if prop in properties:
                        properties[prop].append(prop_set[prop])
                    else:
                        properties[prop] = [prop_set[prop]]

            # Site properties must be present for all atoms in the molecule
            # in order to be used for Molecule instantiation
            for k, v in properties.items():
                if len(v) != len(species):
                    del properties[k]

            new_mol = Molecule(species, coords, charge=charge, site_properties=properties)
            graph_data = json_graph.adjacency_data(new_graph)

            # create new MoleculeGraph
            sub_mols.append(MoleculeGraph(new_mol, graph_data=graph_data))

        return sub_mols

    def split_molecule_subgraphs(self, bonds, allow_reverse=False, alterations=None):
        """
        Split MoleculeGraph into two or more MoleculeGraphs by
        breaking a set of bonds. This function uses
        MoleculeGraph.break_edge repeatedly to create
        disjoint graphs (two or more separate molecules).
        This function does not only alter the graph
        information, but also changes the underlying
        Molecules.
        If the bonds parameter does not include sufficient
        bonds to separate two molecule fragments, then this
        function will fail.
        Currently, this function naively assigns the charge
        of the total molecule to a single submolecule. A
        later effort will be to actually accurately assign
        charge.
        NOTE: This function does not modify the original
        MoleculeGraph. It creates a copy, modifies that, and
        returns two or more new MoleculeGraph objects.
        :param bonds: list of tuples (from_index, to_index)
            representing bonds to be broken to split the MoleculeGraph.
        :param alterations: a dict {(from_index, to_index): alt},
            where alt is a dictionary including weight and/or edge
            properties to be changed following the split.
        :param allow_reverse: If allow_reverse is True, then break_edge will
            attempt to break both (from_index, to_index) and, failing that,
            will attempt to break (to_index, from_index).
        :return: list of MoleculeGraphs
        """

        self.set_node_attributes()
        original = copy.deepcopy(self)

        for bond in bonds:
            original.break_edge(bond[0], bond[1], allow_reverse=allow_reverse)

        if nx.is_weakly_connected(original.graph):
            raise MolGraphSplitError(
                "Cannot split molecule; \
                                MoleculeGraph is still connected."
            )

        # alter any bonds before partition, to avoid remapping
        if alterations is not None:
            for (u, v) in alterations.keys():
                if "weight" in alterations[(u, v)]:
                    weight = alterations[(u, v)]["weight"]
                    del alterations[(u, v)]["weight"]
                    edge_properties = alterations[(u, v)] if len(alterations[(u, v)]) != 0 else None
                    original.alter_edge(u, v, new_weight=weight, new_edge_properties=edge_properties)
                else:
                    original.alter_edge(u, v, new_edge_properties=alterations[(u, v)])

        return original.get_disconnected_fragments()

    def build_unique_fragments(self):
        """
        Find all possible fragment combinations of the MoleculeGraphs (in other
        words, all connected induced subgraphs)

        :return:
        """
        self.set_node_attributes()

        graph = self.graph.to_undirected()

        # find all possible fragments, aka connected induced subgraphs
        frag_dict = {}
        for ii in range(1, len(self.molecule)):
            for combination in combinations(graph.nodes, ii):
                mycomp = []
                for idx in combination:
                    mycomp.append(str(self.molecule[idx].specie))
                mycomp = "".join(sorted(mycomp))
                subgraph = nx.subgraph(graph, combination)
                if nx.is_connected(subgraph):
                    mykey = mycomp + str(len(subgraph.edges()))
                    if mykey not in frag_dict:
                        frag_dict[mykey] = [copy.deepcopy(subgraph)]
                    else:
                        frag_dict[mykey].append(copy.deepcopy(subgraph))

        # narrow to all unique fragments using graph isomorphism
        unique_frag_dict = {}
        for key in frag_dict:
            unique_frags = []
            for frag in frag_dict[key]:
                found = False
                for f in unique_frags:
                    if _isomorphic(frag, f):
                        found = True
                        break
                if not found:
                    unique_frags.append(frag)
            unique_frag_dict[key] = copy.deepcopy(unique_frags)

        # convert back to molecule graphs
        unique_mol_graph_dict = {}
        for key in unique_frag_dict:
            unique_mol_graph_list = []
            for fragment in unique_frag_dict[key]:
                mapping = {e: i for i, e in enumerate(sorted(fragment.nodes))}
                remapped = nx.relabel_nodes(fragment, mapping)

                species = nx.get_node_attributes(remapped, "specie")
                coords = nx.get_node_attributes(remapped, "coords")

                edges = {}

                for from_index, to_index, key in remapped.edges:
                    edge_props = fragment.get_edge_data(from_index, to_index, key=key)

                    edges[(from_index, to_index)] = edge_props

                unique_mol_graph_list.append(
                    self.with_edges(
                        Molecule(species=species, coords=coords, charge=self.molecule.charge),
                        edges,
                    )
                )

            frag_key = (
                str(unique_mol_graph_list[0].molecule.composition.alphabetical_formula)
                + " E"
                + str(len(unique_mol_graph_list[0].graph.edges()))
            )
            unique_mol_graph_dict[frag_key] = copy.deepcopy(unique_mol_graph_list)
        return unique_mol_graph_dict

    def substitute_group(
        self,
        index,
        func_grp,
        strategy,
        bond_order=1,
        graph_dict=None,
        strategy_params=None,
    ):
        """
        Builds off of Molecule.substitute to replace an atom in self.molecule
        with a functional group. This method also amends self.graph to
        incorporate the new functional group.

        NOTE: using a MoleculeGraph will generally produce a different graph
        compared with using a Molecule or str (when not using graph_dict).

        :param index: Index of atom to substitute.
        :param func_grp: Substituent molecule. There are three options:

            1. Providing an actual molecule as the input. The first atom
                must be a DummySpecies X, indicating the position of
                nearest neighbor. The second atom must be the next
                nearest atom. For example, for a methyl group
                substitution, func_grp should be X-CH3, where X is the
                first site and C is the second site. What the code will
                do is to remove the index site, and connect the nearest
                neighbor to the C atom in CH3. The X-C bond indicates the
                directionality to connect the atoms.
            2. A string name. The molecule will be obtained from the
                relevant template in func_groups.json.
            3. A MoleculeGraph object.
        :param strategy: Class from pymatgen.analysis.local_env.
        :param bond_order: A specified bond order to calculate the bond
                length between the attached functional group and the nearest
                neighbor site. Defaults to 1.
        :param graph_dict: Dictionary representing the bonds of the functional
                group (format: {(u, v): props}, where props is a dictionary of
                properties, including weight. If None, then the algorithm
                will attempt to automatically determine bonds using one of
                a list of strategies defined in pymatgen.analysis.local_env.
        :param strategy_params: dictionary of keyword arguments for strategy.
                If None, default parameters will be used.
        :return:
        """

        def map_indices(grp):
            grp_map = {}

            # Get indices now occupied by functional group
            # Subtracting 1 because the dummy atom X should not count
            atoms = len(grp) - 1
            offset = len(self.molecule) - atoms

            for i in range(atoms):
                grp_map[i] = i + offset

            return grp_map

        # Work is simplified if a graph is already in place
        if isinstance(func_grp, MoleculeGraph):

            self.molecule.substitute(index, func_grp.molecule, bond_order=bond_order)

            mapping = map_indices(func_grp.molecule)

            for (u, v) in list(func_grp.graph.edges()):
                edge_props = func_grp.graph.get_edge_data(u, v)[0]
                weight = None
                if "weight" in edge_props.keys():
                    weight = edge_props["weight"]
                    del edge_props["weight"]
                self.add_edge(mapping[u], mapping[v], weight=weight, edge_properties=edge_props)

        else:
            if isinstance(func_grp, Molecule):
                func_grp = copy.deepcopy(func_grp)
            else:
                try:
                    func_grp = copy.deepcopy(FunctionalGroups[func_grp])
                except Exception:
                    raise RuntimeError("Can't find functional group in list. " "Provide explicit coordinate instead")

            self.molecule.substitute(index, func_grp, bond_order=bond_order)

            mapping = map_indices(func_grp)

            # Remove dummy atom "X"
            func_grp.remove_species("X")

            if graph_dict is not None:
                for (u, v) in graph_dict.keys():
                    edge_props = graph_dict[(u, v)]
                    if "weight" in edge_props.keys():
                        weight = edge_props["weight"]
                        del edge_props["weight"]
                    self.add_edge(
                        mapping[u],
                        mapping[v],
                        weight=weight,
                        edge_properties=edge_props,
                    )

            else:
                if strategy_params is None:
                    strategy_params = {}
                strat = strategy(**strategy_params)
                graph = self.with_local_env_strategy(func_grp, strat)

                for (u, v) in list(graph.graph.edges()):
                    edge_props = graph.graph.get_edge_data(u, v)[0]
                    weight = None
                    if "weight" in edge_props.keys():
                        weight = edge_props["weight"]
                        del edge_props["weight"]

                    if 0 not in list(graph.graph.nodes()):
                        # If graph indices have different indexing
                        u, v = (u - 1), (v - 1)

                    self.add_edge(
                        mapping[u],
                        mapping[v],
                        weight=weight,
                        edge_properties=edge_props,
                    )

    def replace_group(
        self,
        index,
        func_grp,
        strategy,
        bond_order=1,
        graph_dict=None,
        strategy_params=None,
    ):
        """
        Builds off of Molecule.substitute and MoleculeGraph.substitute_group
        to replace a functional group in self.molecule with a functional group.
        This method also amends self.graph to incorporate the new functional
        group.

        TODO: Figure out how to replace into a ring structure.

        :param index: Index of atom to substitute.
        :param func_grp: Substituent molecule. There are three options:

            1. Providing an actual molecule as the input. The first atom
               must be a DummySpecies X, indicating the position of
               nearest neighbor. The second atom must be the next
               nearest atom. For example, for a methyl group
               substitution, func_grp should be X-CH3, where X is the
               first site and C is the second site. What the code will
               do is to remove the index site, and connect the nearest
               neighbor to the C atom in CH3. The X-C bond indicates the
               directionality to connect the atoms.
            2. A string name. The molecule will be obtained from the
               relevant template in func_groups.json.
            3. A MoleculeGraph object.
        :param strategy: Class from pymatgen.analysis.local_env.
        :param bond_order: A specified bond order to calculate the bond
            length between the attached functional group and the nearest
            neighbor site. Defaults to 1.
        :param graph_dict: Dictionary representing the bonds of the functional
            group (format: {(u, v): props}, where props is a dictionary of
            properties, including weight. If None, then the algorithm
            will attempt to automatically determine bonds using one of
            a list of strategies defined in pymatgen.analysis.local_env.
        :param strategy_params: dictionary of keyword arguments for strategy.
            If None, default parameters will be used.
        :return:
        """

        self.set_node_attributes()
        neighbors = self.get_connected_sites(index)

        # If the atom at index is terminal
        if len(neighbors) == 1:
            self.substitute_group(
                index,
                func_grp,
                strategy,
                bond_order=bond_order,
                graph_dict=graph_dict,
                strategy_params=strategy_params,
            )

        else:
            rings = self.find_rings(including=[index])
            if len(rings) != 0:
                raise RuntimeError(
                    "Currently functional group replacement" "cannot occur at an atom within a ring" "structure."
                )

            to_remove = set()
            sizes = dict()
            disconnected = self.graph.to_undirected()
            disconnected.remove_node(index)
            for neighbor in neighbors:
                sizes[neighbor[2]] = len(nx.descendants(disconnected, neighbor[2]))

            keep = max(sizes, key=lambda x: sizes[x])
            for i in sizes.keys():
                if i != keep:
                    to_remove.add(i)

            self.remove_nodes(list(to_remove))
            self.substitute_group(
                index,
                func_grp,
                strategy,
                bond_order=bond_order,
                graph_dict=graph_dict,
                strategy_params=strategy_params,
            )

    def find_rings(self, including=None):
        """
        Find ring structures in the MoleculeGraph.

        :param including: list of site indices. If
            including is not None, then find_rings will
            only return those rings including the specified
            sites. By default, this parameter is None, and
            all rings will be returned.
        :return: dict {index:cycle}. Each
            entry will be a ring (cycle, in graph theory terms) including the index
            found in the Molecule. If there is no cycle including an index, the
            value will be an empty list.
        """

        # Copies self.graph such that all edges (u, v) matched by edges (v, u)
        undirected = self.graph.to_undirected()
        directed = undirected.to_directed()

        cycles_nodes = []
        cycles_edges = []

        # Remove all two-edge cycles
        all_cycles = [c for c in nx.simple_cycles(directed) if len(c) > 2]

        # Using to_directed() will mean that each cycle always appears twice
        # So, we must also remove duplicates
        unique_sorted = []
        unique_cycles = []
        for cycle in all_cycles:
            if sorted(cycle) not in unique_sorted:
                unique_sorted.append(sorted(cycle))
                unique_cycles.append(cycle)

        if including is None:
            cycles_nodes = unique_cycles
        else:
            for i in including:
                for cycle in unique_cycles:
                    if i in cycle and cycle not in cycles_nodes:
                        cycles_nodes.append(cycle)

        for cycle in cycles_nodes:
            edges = []
            for i, e in enumerate(cycle):
                edges.append((cycle[i - 1], e))
            cycles_edges.append(edges)

        return cycles_edges

    def get_connected_sites(self, n):
        """
        Returns a named tuple of neighbors of site n:
        periodic_site, jimage, index, weight.
        Index is the index of the corresponding site
        in the original structure, weight can be
        None if not defined.
        :param n: index of Site in Molecule
        :param jimage: lattice vector of site
        :return: list of ConnectedSite tuples,
            sorted by closest first
        """

        connected_sites = set()

        out_edges = list(self.graph.out_edges(n, data=True))
        in_edges = list(self.graph.in_edges(n, data=True))

        for u, v, d in out_edges + in_edges:

            weight = d.get("weight", None)

            if v == n:
                site = self.molecule[u]
                dist = self.molecule[v].distance(self.molecule[u])

                connected_site = ConnectedSite(site=site, jimage=(0, 0, 0), index=u, weight=weight, dist=dist)
            else:
                site = self.molecule[v]
                dist = self.molecule[u].distance(self.molecule[v])

                connected_site = ConnectedSite(site=site, jimage=(0, 0, 0), index=v, weight=weight, dist=dist)

            connected_sites.add(connected_site)

        # return list sorted by closest sites first
        connected_sites = list(connected_sites)
        connected_sites.sort(key=lambda x: x.dist)

        return connected_sites

    def get_coordination_of_site(self, n):
        """
        Returns the number of neighbors of site n.
        In graph terms, simply returns degree
        of node corresponding to site n.
        :param n: index of site
        :return (int):
        """
        number_of_self_loops = sum([1 for n, v in self.graph.edges(n) if n == v])
        return self.graph.degree(n) - number_of_self_loops

    def draw_graph_to_file(
        self,
        filename="graph",
        diff=None,
        hide_unconnected_nodes=False,
        hide_image_edges=True,
        edge_colors=False,
        node_labels=False,
        weight_labels=False,
        image_labels=False,
        color_scheme="VESTA",
        keep_dot=False,
        algo="fdp",
    ):
        """
        Draws graph using GraphViz.

        The networkx graph object itself can also be drawn
        with networkx's in-built graph drawing methods, but
        note that this might give misleading results for
        multigraphs (edges are super-imposed on each other).

        If visualization is difficult to interpret,
        `hide_image_edges` can help, especially in larger
        graphs.

        :param filename: filename to output, will detect filetype
            from extension (any graphviz filetype supported, such as
            pdf or png)
        :param diff (StructureGraph): an additional graph to
            compare with, will color edges red that do not exist in diff
            and edges green that are in diff graph but not in the
            reference graph
        :param hide_unconnected_nodes: if True, hide unconnected
            nodes
        :param hide_image_edges: if True, do not draw edges that
            go through periodic boundaries
        :param edge_colors (bool): if True, use node colors to
            color edges
        :param node_labels (bool): if True, label nodes with
            species and site index
        :param weight_labels (bool): if True, label edges with
            weights
        :param image_labels (bool): if True, label edges with
            their periodic images (usually only used for debugging,
            edges to periodic images always appear as dashed lines)
        :param color_scheme (str): "VESTA" or "JMOL"
        :param keep_dot (bool): keep GraphViz .dot file for later
            visualization
        :param algo: any graphviz algo, "neato" (for simple graphs)
            or "fdp" (for more crowded graphs) usually give good outputs
        :return:
        """

        if not which(algo):
            raise RuntimeError("StructureGraph graph drawing requires " "GraphViz binaries to be in the path.")

        # Developer note: NetworkX also has methods for drawing
        # graphs using matplotlib, these also work here. However,
        # a dedicated tool like GraphViz allows for much easier
        # control over graph appearance and also correctly displays
        # mutli-graphs (matplotlib can superimpose multiple edges).

        g = self.graph.copy()

        g.graph = {"nodesep": 10.0, "dpi": 300, "overlap": "false"}

        # add display options for nodes
        for n in g.nodes():
            # get label by species name
            label = "{}({})".format(str(self.molecule[n].specie), n) if node_labels else ""

            # use standard color scheme for nodes
            c = EL_COLORS[color_scheme].get(str(self.molecule[n].specie.symbol), [0, 0, 0])

            # get contrasting font color
            # magic numbers account for perceived luminescence
            # https://stackoverflow.com/questions/1855884/determine-font-color-based-on-background-color
            fontcolor = "#000000" if 1 - (c[0] * 0.299 + c[1] * 0.587 + c[2] * 0.114) / 255 < 0.5 else "#ffffff"

            # convert color to hex string
            color = "#{:02x}{:02x}{:02x}".format(c[0], c[1], c[2])

            g.add_node(
                n,
                fillcolor=color,
                fontcolor=fontcolor,
                label=label,
                fontname="Helvetica-bold",
                style="filled",
                shape="circle",
            )

        edges_to_delete = []

        # add display options for edges
        for u, v, k, d in g.edges(keys=True, data=True):

            # retrieve from/to images, set as origin if not defined
            if "to_image" in d:
                to_image = d["to_jimage"]
            else:
                to_image = (0, 0, 0)

            # set edge style
            d["style"] = "solid"
            if to_image != (0, 0, 0):
                d["style"] = "dashed"
                if hide_image_edges:
                    edges_to_delete.append((u, v, k))

            # don't show edge directions
            d["arrowhead"] = "none"

            # only add labels for images that are not the origin
            if image_labels:
                d["headlabel"] = "" if to_image == (0, 0, 0) else "to {}".format((to_image))
                d["arrowhead"] = "normal" if d["headlabel"] else "none"

            # optionally color edges using node colors
            color_u = g.node[u]["fillcolor"]
            color_v = g.node[v]["fillcolor"]
            d["color_uv"] = "{};0.5:{};0.5".format(color_u, color_v) if edge_colors else "#000000"

            # optionally add weights to graph
            if weight_labels:
                units = g.graph.get("edge_weight_units", "")
                if d.get("weight"):
                    d["label"] = "{:.2f} {}".format(d["weight"], units)

            # update edge with our new style attributes
            g.edges[u, v, k].update(d)

        # optionally remove periodic image edges,
        # these can be confusing due to periodic boundaries
        if hide_image_edges:
            for edge_to_delete in edges_to_delete:
                g.remove_edge(*edge_to_delete)

        # optionally hide unconnected nodes,
        # these can appear when removing periodic edges
        if hide_unconnected_nodes:
            g = g.subgraph([n for n in g.degree() if g.degree()[n] != 0])

        # optionally highlight differences with another graph
        if diff:
            diff = self.diff(diff, strict=True)
            green_edges = []
            red_edges = []
            for u, v, k, d in g.edges(keys=True, data=True):
                if (u, v, d["to_jimage"]) in diff["self"]:
                    # edge has been deleted
                    red_edges.append((u, v, k))
                elif (u, v, d["to_jimage"]) in diff["other"]:
                    # edge has been added
                    green_edges.append((u, v, k))
            for u, v, k in green_edges:
                g.edges[u, v, k].update({"color_uv": "#00ff00"})
            for u, v, k in red_edges:
                g.edges[u, v, k].update({"color_uv": "#ff0000"})

        basename, extension = os.path.splitext(filename)
        extension = extension[1:]

        write_dot(g, basename + ".dot")

        with open(filename, "w") as f:

            args = [algo, "-T", extension, basename + ".dot"]
            rs = subprocess.Popen(args, stdout=f, stdin=subprocess.PIPE, close_fds=True)
            rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError("{} exited with return code {}.".format(algo, rs.returncode))

        if not keep_dot:
            os.remove(basename + ".dot")

    def as_dict(self):
        """
        As in :Class: `pymatgen.core.Molecule` except
        with using `to_dict_of_dicts` from NetworkX
        to store graph information.
        """

        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "molecule": self.molecule.as_dict(),
            "graphs": json_graph.adjacency_data(self.graph),
        }

        return d

    @classmethod
    def from_dict(cls, d):
        """
        As in :Class: `pymatgen.core.Molecule` except
        restoring graphs using `from_dict_of_dicts`
        from NetworkX to restore graph information.
        """
        m = Molecule.from_dict(d["molecule"])
        return cls(m, d["graphs"])

    @classmethod
    def _edges_to_string(cls, g):

        header = "from    to  to_image    "
        header_line = "----  ----  ------------"
        edge_weight_name = g.graph["edge_weight_name"]
        if edge_weight_name:
            print_weights = ["weight"]
            edge_label = g.graph["edge_weight_name"]
            edge_weight_units = g.graph["edge_weight_units"]
            if edge_weight_units:
                edge_label += " ({})".format(edge_weight_units)
            header += "  {}".format(edge_label)
            header_line += "  {}".format("-" * max([18, len(edge_label)]))
        else:
            print_weights = False

        s = header + "\n" + header_line + "\n"

        edges = list(g.edges(data=True))

        # sort edges for consistent ordering
        edges.sort(key=itemgetter(0, 1))

        if print_weights:
            for u, v, data in edges:
                s += "{:4}  {:4}  {:12}  {:.3e}\n".format(
                    u, v, str(data.get("to_jimage", (0, 0, 0))), data.get("weight", 0)
                )
        else:
            for u, v, data in edges:
                s += "{:4}  {:4}  {:12}\n".format(u, v, str(data.get("to_jimage", (0, 0, 0))))

        return s

    def __str__(self):
        s = "Molecule Graph"
        s += "\nMolecule: \n{}".format(self.molecule.__str__())
        s += "\nGraph: {}\n".format(self.name)
        s += self._edges_to_string(self.graph)
        return s

    def __repr__(self):
        s = "Molecule Graph"
        s += "\nMolecule: \n{}".format(self.molecule.__repr__())
        s += "\nGraph: {}\n".format(self.name)
        s += self._edges_to_string(self.graph)
        return s

    def __len__(self):
        """
        :return: length of Molecule / number of nodes in graph
        """
        return len(self.molecule)

    def sort(self, key=None, reverse=False):
        """
        Same as Molecule.sort(), also remaps nodes in graph.
        :param key:
        :param reverse:
        :return:
        """

        old_molecule = self.molecule.copy()

        # sort Molecule
        self.molecule._sites = sorted(self.molecule._sites, key=key, reverse=reverse)

        # apply Molecule ordering to graph
        mapping = {idx: self.molecule.index(site) for idx, site in enumerate(old_molecule)}
        self.graph = nx.relabel_nodes(self.graph, mapping, copy=True)

        # normalize directions of edges
        edges_to_remove = []
        edges_to_add = []
        for u, v, k, d in self.graph.edges(keys=True, data=True):
            if v < u:
                new_v, new_u, new_d = u, v, d.copy()
                new_d["to_jimage"] = (0, 0, 0)
                edges_to_remove.append((u, v, k))
                edges_to_add.append((new_u, new_v, new_d))

        # add/delete marked edges
        for edges_to_remove in edges_to_remove:
            self.graph.remove_edge(*edges_to_remove)
        for (u, v, d) in edges_to_add:
            self.graph.add_edge(u, v, **d)

    def __copy__(self):
        return MoleculeGraph.from_dict(self.as_dict())

    def __eq__(self, other):
        """
        Two MoleculeGraphs are equal if they have equal Molecules,
        and have the same edges between Sites. Edge weights can be
        different and MoleculeGraphs can still be considered equal.

        :param other: MoleculeGraph
        :return (bool):
        """

        # sort for consistent node indices
        # PeriodicSite should have a proper __hash__() value,
        # using its frac_coords as a convenient key
        try:
            mapping = {tuple(site.coords): self.molecule.index(site) for site in other.molecule}
        except ValueError:
            return False
        other_sorted = other.__copy__()
        other_sorted.sort(key=lambda site: mapping[tuple(site.coords)])

        edges = {(u, v) for u, v, d in self.graph.edges(keys=False, data=True)}

        edges_other = {(u, v) for u, v, d in other_sorted.graph.edges(keys=False, data=True)}

        return (edges == edges_other) and (self.molecule == other_sorted.molecule)

    def isomorphic_to(self, other):
        """
        Checks if the graphs of two MoleculeGraphs are isomorphic to one
        another. In order to prevent problems with misdirected edges, both
        graphs are converted into undirected nx.Graph objects.

        :param other: MoleculeGraph object to be compared.
        :return: bool
        """
        if len(self.molecule) != len(other.molecule):
            return False
        if self.molecule.composition.alphabetical_formula != other.molecule.composition.alphabetical_formula:
            return False
        if len(self.graph.edges()) != len(other.graph.edges()):
            return False
        return _isomorphic(self.graph, other.graph)

    def diff(self, other, strict=True):
        """
        Compares two MoleculeGraphs. Returns dict with
        keys 'self', 'other', 'both' with edges that are
        present in only one MoleculeGraph ('self' and
        'other'), and edges that are present in both.

        The Jaccard distance is a simple measure of the
        dissimilarity between two MoleculeGraphs (ignoring
        edge weights), and is defined by 1 - (size of the
        intersection / size of the union) of the sets of
        edges. This is returned with key 'dist'.

        Important note: all node indices are in terms
        of the MoleculeGraph this method is called
        from, not the 'other' MoleculeGraph: there
        is no guarantee the node indices will be the
        same if the underlying Molecules are ordered
        differently.

        :param other: MoleculeGraph
        :param strict: if False, will compare bonds
            from different Molecules, with node indices
            replaced by Species strings, will not count
            number of occurrences of bonds
        :return:
        """

        if self.molecule != other.molecule and strict:
            return ValueError("Meaningless to compare MoleculeGraphs if " "corresponding Molecules are different.")

        if strict:
            # sort for consistent node indices
            # PeriodicSite should have a proper __hash__() value,
            # using its frac_coords as a convenient key
            mapping = {tuple(site.frac_coords): self.molecule.index(site) for site in other.molecule}
            other_sorted = other.__copy__()
            other_sorted.sort(key=lambda site: mapping[tuple(site.frac_coords)])

            edges = {(u, v, d.get("to_jimage", (0, 0, 0))) for u, v, d in self.graph.edges(keys=False, data=True)}

            edges_other = {
                (u, v, d.get("to_jimage", (0, 0, 0))) for u, v, d in other_sorted.graph.edges(keys=False, data=True)
            }

        else:

            edges = {
                (str(self.molecule[u].specie), str(self.molecule[v].specie))
                for u, v, d in self.graph.edges(keys=False, data=True)
            }

            edges_other = {
                (str(other.structure[u].specie), str(other.structure[v].specie))
                for u, v, d in other.graph.edges(keys=False, data=True)
            }

        if len(edges) == 0 and len(edges_other) == 0:
            jaccard_dist = 0  # by definition
        else:
            jaccard_dist = 1 - len(edges.intersection(edges_other)) / len(edges.union(edges_other))

        return {
            "self": edges - edges_other,
            "other": edges_other - edges,
            "both": edges.intersection(edges_other),
            "dist": jaccard_dist,
        }
