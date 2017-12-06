# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import warnings
import subprocess
import numpy as np
import os.path

from pymatgen.core import Structure, Lattice, PeriodicSite, Molecule
from pymatgen.util.coord import lattice_points_in_supercell
from pymatgen.vis.structure_vtk import EL_COLORS

from monty.json import MSONable
from monty.os.path import which
from operator import itemgetter
from collections import namedtuple
from scipy.spatial import KDTree

try:
    import networkx as nx
    from networkx.readwrite import json_graph
    from networkx.drawing.nx_agraph import write_dot
except ImportError:
    raise ImportError("This module requires the NetworkX "
                      "graph library to be installed.")

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

__author__ = "Matthew Horton"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Beta"
__date__ = "August 2017"

ConnectedSite = namedtuple('ConnectedSite', 'periodic_site, jimage, index, weight, dist')

class StructureGraph(MSONable):

    def __init__(self, structure, graph_data):
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

        :param *args: same as in :class: `pymatgen.core.Structure`
        :param graph_data: dict containing graph information in
        dict format, not intended to be constructed manually
        """

        self.structure = structure
        self.graph = nx.readwrite.json_graph.adjacency_graph(graph_data)

        # tidy up edge attr dicts, reading to/from json duplicates
        # information
        for u, v, k, d in self.graph.edges(keys=True, data=True):
            if 'id' in d:
                del d['id']
            if 'key' in d:
                del d['key']

    @classmethod
    def with_empty_graph(cls, structure, name="bonds",
                        edge_weight_name=None,
                         edge_weight_units=None):
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
            raise ValueError("Please specify units associated "
                             "with your edge weights. Can be "
                             "empty string if arbitrary or "
                             "dimensionless.")

        # construct graph with one node per site
        # graph attributes don't change behavior of graph,
        # they're just for book-keeping
        graph = nx.MultiDiGraph(edge_weight_name=edge_weight_name,
                                edge_weight_units=edge_weight_units,
                                name=name)
        graph.add_nodes_from(range(len(structure)))

        graph_data = json_graph.adjacency_data(graph)

        return cls(structure, graph_data=graph_data)

    @staticmethod
    def with_local_env_strategy(structure, strategy):
        """
        Constructor for StructureGraph, using a strategy
        from  :Class: `pymatgen.analysis.local_env`.

        :param structure: Structure object
        :param strategy: an instance of a
         :Class: `pymatgen.analysis.local_env.NearNeighbors`
         object
        :return:
        """

        sg = StructureGraph.with_empty_graph(structure, name="bonds",
                                             edge_weight_name="weight",
                                             edge_weight_units="")

        for n in range(len(structure)):
            neighbors = strategy.get_nn_info(structure, n)
            for neighbor in neighbors:

                # local_env will always try to add two edges
                # for any one bond, one from site u to site v
                # and another form site v to site u: this is
                # harmless, so warn_duplicates=False
                sg.add_edge(from_index=n,
                            from_jimage=(0, 0, 0),
                            to_index=neighbor['site_index'],
                            to_jimage=neighbor['image'],
                            weight=neighbor['weight'],
                            warn_duplicates=False)

        return sg

    @property
    def name(self):
        """
        :return: Name of graph
        """
        return self.graph.graph['name']

    @property
    def edge_weight_name(self):
        """
        :return: Name of the edge weight property of graph
        """
        return self.graph.graph['edge_weight_name']

    @property
    def edge_weight_unit(self):
        """
        :return: Units of the edge weight property of graph
        """
        return self.graph.graph['edge_weight_units']

    def add_edge(self, from_index, to_index,
                 from_jimage=(0, 0, 0), to_jimage=None,
                 weight=None, warn_duplicates=True):
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
            warnings.warn("Please specify to_jimage to be unambiguous, "
                          "trying to automatically detect.")
            dist, to_jimage = self.structure[from_index].distance_and_image(self.structure[to_index])
            if dist == 0:
                # this will happen when from_index == to_index,
                # typically in primitive single-atom lattices
                images = [1, 0, 0], [0, 1, 0], [0, 0, 1]
                dists = []
                for image in images:
                    dists.append(self.structure[from_index].distance_and_image(self.structure[from_index],
                                                                               jimage=image)[0])
                dist = min(dists)
            equiv_sites = self.structure.get_neighbors_in_shell(self.structure[from_index].coords,
                                                                dist,
                                                                dist*0.01,
                                                                include_index=True)
            for site, dist, to_index in equiv_sites:
                to_jimage = np.subtract(site.frac_coords, self.structure[from_index].frac_coords)
                to_jimage = to_jimage.astype(int)
                self.add_edge(from_index=from_index, from_jimage=(0, 0, 0),
                              to_jimage=to_jimage, to_index=to_index)
            return

        from_jimage, to_jimage = tuple(from_jimage), tuple(to_jimage)

        # check we're not trying to add a duplicate edge
        # there should only ever be at most one edge
        # between a given (site, jimage) pair and another
        # (site, jimage) pair
        existing_edge_data = self.graph.get_edge_data(from_index, to_index)
        if existing_edge_data:
            for key, d in existing_edge_data.items():
                if d["to_jimage"] == to_jimage:
                    if warn_duplicates:
                        warnings.warn("Trying to add an edge that already exists from "
                                      "site {} to site {} in {}.".format(from_index,
                                                                         to_index,
                                                                         to_jimage))
                    return

        if weight:
            self.graph.add_edge(from_index, to_index,
                                from_jimage=from_jimage, to_jimage=to_jimage,
                                weight=weight)
        else:
            self.graph.add_edge(from_index, to_index,
                                from_jimage=from_jimage, to_jimage=to_jimage)

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

        out_edges = [(u, v, d, 'out') for u, v, d in self.graph.out_edges(n, data=True)]
        in_edges = [(u, v, d, 'in') for u, v, d in self.graph.in_edges(n, data=True)]

        for u, v, d, dir in out_edges + in_edges:

            to_jimage = d['to_jimage']

            if dir == 'in':
                u, v = v, u
                to_jimage = np.multiply(-1, to_jimage)

            site_d = self.structure[v].as_dict()
            site_d['abc'] = np.add(site_d['abc'], to_jimage).tolist()
            to_jimage = tuple(map(int, np.add(to_jimage, jimage)))
            periodic_site = PeriodicSite.from_dict(site_d)

            weight = d.get('weight', None)

            # from_site if jimage arg != (0, 0, 0)
            relative_jimage = np.subtract(to_jimage, jimage)
            dist = self.structure[u].distance(self.structure[v], jimage=relative_jimage)

            connected_site = ConnectedSite(periodic_site=periodic_site,
                                           jimage=to_jimage,
                                           index=v,
                                           weight=weight,
                                           dist=dist)

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
        return self.graph.degree(n)

    def draw_graph_to_file(self, filename="graph",
                           diff=None,
                           hide_unconnected_nodes=False,
                           hide_image_edges=True,
                           edge_colors=False,
                           node_labels=False,
                           weight_labels=False,
                           image_labels=False,
                           color_scheme="VESTA",
                           keep_dot=False,
                           algo="fdp"):
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
            raise RuntimeError("StructureGraph graph drawing requires "
                               "GraphViz binaries to be in the path.")

        # Developer note: NetworkX also has methods for drawing
        # graphs using matplotlib, these also work here. However,
        # a dedicated tool like GraphViz allows for much easier
        # control over graph appearance and also correctly displays
        # mutli-graphs (matplotlib can superimpose multiple edges).

        g = self.graph.copy()

        g.graph = {'nodesep': 10.0, 'dpi': 300, 'overlap': "false"}

        # add display options for nodes
        for n in g.nodes():

            # get label by species name
            label = "{}({})".format(str(self.structure[n].specie), n) if node_labels else ""

            # use standard color scheme for nodes
            c = EL_COLORS[color_scheme].get(str(self.structure[n].specie.symbol), [0, 0, 0])

            # get contrasting font color
            # magic numbers account for perceived luminescence
            # https://stackoverflow.com/questions/1855884/determine-font-color-based-on-background-color
            fontcolor = '#000000' if 1 - (c[0] * 0.299 + c[1] * 0.587
                                          + c[2] * 0.114) / 255 < 0.5 else '#ffffff'

            # convert color to hex string
            color = "#{:02x}{:02x}{:02x}".format(c[0], c[1], c[2])

            g.add_node(n, fillcolor=color, fontcolor=fontcolor, label=label,
                       fontname="Helvetica-bold", style="filled", shape="circle")

        edges_to_delete = []

        # add display options for edges
        for u, v, k, d in g.edges(keys=True, data=True):

            # retrieve from/to images, set as origin if not defined
            to_image = d['to_jimage']

            # set edge style
            d['style'] = "solid"
            if to_image != (0, 0, 0):
                d['style'] = "dashed"
                if hide_image_edges:
                    edges_to_delete.append((u, v, k))

            # don't show edge directions
            d['arrowhead'] = "none"

            # only add labels for images that are not the origin
            if image_labels:
                d['headlabel'] = "" if to_image == (0, 0, 0) else "to {}".format((to_image))
                d['arrowhead'] = "normal" if d['headlabel'] else "none"

            # optionally color edges using node colors
            color_u = g.node[u]['fillcolor']
            color_v = g.node[v]['fillcolor']
            d['color_uv'] = "{};0.5:{};0.5".format(color_u, color_v) if edge_colors else "#000000"

            # optionally add weights to graph
            if weight_labels:
                units = g.graph.get('edge_weight_units', "")
                if d.get('weight'):
                    d['label'] = "{:.2f} {}".format(d['weight'], units)

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
                if (u, v, d['to_jimage']) in diff['self']:
                    # edge has been deleted
                    red_edges.append((u, v, k))
                elif (u, v, d['to_jimage']) in diff['other']:
                    # edge has been added
                    green_edges.append((u, v, k))
            for u, v, k in green_edges:
                g.edges[u, v, k].update({'color_uv': '#00ff00'})
            for u, v, k in red_edges:
                g.edges[u, v, k].update({'color_uv': '#ff0000'})

        basename, extension = os.path.splitext(filename)
        extension = extension[1:]

        write_dot(g, basename+".dot")

        with open(filename, "w") as f:

            args = [algo, "-T", extension, basename+".dot"]
            rs = subprocess.Popen(args,
                                  stdout=f,
                                  stdin=subprocess.PIPE, close_fds=True)
            rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError("{} exited with return code {}.".format(algo, rs.returncode))

        if not keep_dot:
            os.remove(basename+".dot")

    def as_dict(self):
        """
        As in :Class: `pymatgen.core.Structure` except
        with using `to_dict_of_dicts` from NetworkX
        to store graph information.
        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "structure": self.structure.as_dict(),
             "graphs": json_graph.adjacency_data(self.graph)}

        return d

    @classmethod
    def from_dict(cls, d):
        """
        As in :Class: `pymatgen.core.Structure` except
        restoring graphs using `from_dict_of_dicts`
        from NetworkX to restore graph information.
        """
        s = Structure.from_dict(d['structure'])
        return cls(s, d['graphs'])

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
            raise NotImplementedError('Not tested with 3x3 scaling matrices yet.')
        new_lattice = Lattice(np.dot(scale_matrix, self.structure.lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        new_graphs = []

        for v in c_lat:

            # create a map of nodes from original graph to its image
            mapping = {n: n + len(new_sites) for n in range(len(self.structure))}

            for idx, site in enumerate(self.structure):

                s = PeriodicSite(site.species_and_occu, site.coords + v,
                                 new_lattice, properties=site.properties,
                                 coords_are_cartesian=True, to_unit_cell=False)

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
        edges_inside_supercell = [{u, v} for u, v, d in new_g.edges(data=True)
                                  if d['to_jimage'] == (0, 0, 0)]
        new_periodic_images = []

        orig_lattice = self.structure.lattice

        # use k-d tree to match given position to an
        # existing Site in Structure
        kd_tree = KDTree(new_structure.cart_coords)

        # tolerance in Å for sites to be considered equal
        # this could probably be a lot smaller
        tol = 0.05

        for u, v, k, d in new_g.edges(keys=True, data=True):

            to_jimage = d['to_jimage']  # for node v

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
                # new supercell, and get absolute Cartesian
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
                    new_d['to_jimage'] = (0, 0, 0)

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
                    v_expec_image = v_expec_image - v_expec_image%1

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
                            new_to_jimage = tuple(np.multiply(-1, d['to_jimage']).astype(int))

                        new_d['to_jimage'] = new_to_jimage

                        edges_to_remove.append((u, v, k))

                        if (new_u, new_v, new_to_jimage) not in new_periodic_images:
                            edges_to_add.append((new_u, new_v, new_d))
                            new_periodic_images.append((new_u, new_v, new_to_jimage))

        logger.debug("Removing {} edges, adding {} new edges.".format(len(edges_to_remove),
                                                                      len(edges_to_add)))

        # add/delete marked edges
        for edges_to_remove in edges_to_remove:
            new_g.remove_edge(*edges_to_remove)
        for (u, v, d) in edges_to_add:
            new_g.add_edge(u, v, **d)

        # return new instance of StructureGraph with supercell
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "structure": new_structure.as_dict(),
             "graphs": json_graph.adjacency_data(new_g)}

        sg = StructureGraph.from_dict(d)

        return sg

    def __rmul__(self, other):
        return self.__mul__(other)

    def _edges_to_string(self, g):

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
            header_line += "  {}".format("-"*max([18, len(edge_label)]))
        else:
            print_weights = False

        s = header + "\n" + header_line + "\n"

        edges = list(g.edges(data=True))

        # sort edges for consistent ordering
        edges.sort(key=itemgetter(0,1))

        if print_weights:
            for u, v, data in edges:
                s += "{:4}  {:4}  {:12}  {:.12E}\n".format(u, v, str(data.get("to_jimage", (0, 0, 0))),
                                                           data.get("weight", 0))
        else:
            for u, v, data in edges:
                s += "{:4}  {:4}  {:12}\n".format(u, v,
                                                  str(data.get("to_jimage", (0, 0, 0))))

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
        mapping = {idx:self.structure.index(site) for idx, site in enumerate(old_structure)}
        self.graph = nx.relabel_nodes(self.graph, mapping, copy=True)

        # normalize directions of edges
        edges_to_remove = []
        edges_to_add = []
        for u, v, k, d in self.graph.edges(keys=True, data=True):
            if v < u:
                new_v, new_u, new_d = u, v, d.copy()
                new_d['to_jimage'] = tuple(np.multiply(-1, d['to_jimage']).astype(int))
                edges_to_remove.append((u,v,k))
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
        mapping = {tuple(site.frac_coords):self.structure.index(site) for site in other.structure}
        other_sorted = other.__copy__()
        other_sorted.sort(key=lambda site: mapping[tuple(site.frac_coords)])

        edges = {(u, v, d['to_jimage'])
                 for u, v, d in self.graph.edges(keys=False, data=True)}

        edges_other = {(u, v, d['to_jimage'])
                       for u, v, d in other_sorted.graph.edges(keys=False, data=True)}

        return (edges == edges_other) and \
               (self.structure == other_sorted.structure)

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
        replaced by Specie strings, will not count
        number of occurrences of bonds
        :return:
        """

        if self.structure != other.structure and strict:
            return ValueError("Meaningless to compare StructureGraphs if "
                              "corresponding Structures are different.")

        if strict:

            # sort for consistent node indices
            # PeriodicSite should have a proper __hash__() value,
            # using its frac_coords as a convenient key
            mapping = {tuple(site.frac_coords):self.structure.index(site) for site in other.structure}
            other_sorted = other.__copy__()
            other_sorted.sort(key=lambda site: mapping[tuple(site.frac_coords)])

            edges = {(u, v, d['to_jimage'])
                     for u, v, d in self.graph.edges(keys=False, data=True)}

            edges_other = {(u, v, d['to_jimage'])
                           for u, v, d in other_sorted.graph.edges(keys=False, data=True)}

        else:

            edges = {(str(self.structure[u].specie),
                      str(self.structure[v].specie))
                     for u, v, d in self.graph.edges(keys=False, data=True)}

            edges_other = {(str(other.structure[u].specie),
                            str(other.structure[v].specie))
                           for u, v, d in other.graph.edges(keys=False, data=True)}

        if len(edges) == 0 and len(edges_other) == 0:
            jaccard_dist = 0  # by definition
        else:
            jaccard_dist = 1 - len(edges.intersection(edges_other)) / len(edges.union(edges_other))

        return {
            'self': edges - edges_other,
            'other': edges_other - edges,
            'both': edges.intersection(edges_other),
            'dist': jaccard_dist
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
        if getattr(self, '_supercell_sg', None) is None:
            self._supercell_sg = supercell_sg = self*(3,3,3)

        # make undirected to find connected subgraphs
        supercell_sg.graph = nx.Graph(supercell_sg.graph)

        # find subgraphs
        all_subgraphs = list(nx.connected_component_subgraphs(supercell_sg.graph))

        # discount subgraphs that lie across *supercell* boundaries
        # these will subgraphs representing crystals
        molecule_subgraphs = []
        for subgraph in all_subgraphs:
            intersects_boundary = any([d['to_jimage'] != (0, 0, 0)
                                      for u, v, d in subgraph.edges(data=True)])
            if not intersects_boundary:
                molecule_subgraphs.append(subgraph)

        # add specie names to graph to be able to test for isomorphism
        for subgraph in molecule_subgraphs:
            for n in subgraph:
                subgraph.add_node(n, specie=str(supercell_sg.structure[n].specie))

        # now define how we test for isomorphism
        def node_match(n1, n2):
            return n1['specie'] == n2['specie']
        def edge_match(e1, e2):
            if use_weights:
                return e1['weight'] == e2['weight']
            else:
                return True

        # prune duplicate subgraphs
        unique_subgraphs = []
        for subgraph in molecule_subgraphs:

            already_present = [nx.is_isomorphic(subgraph, g,
                                                node_match=node_match,
                                                edge_match=edge_match)
                               for g in unique_subgraphs]

            if not any(already_present):
                unique_subgraphs.append(subgraph)

        # get Molecule objects for each subgraph
        molecules = []
        for subgraph in unique_subgraphs:

            coords = [supercell_sg.structure[n].coords for n
                      in subgraph.nodes()]
            species = [supercell_sg.structure[n].specie for n
                      in subgraph.nodes()]

            molecule = Molecule(species, coords)

            # shift so origin is at center of mass
            molecule = molecule.get_centered_molecule()

            molecules.append(molecule)

        return molecules
