# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import warnings
import subprocess
import numpy as np
import os.path

from pymatgen.core import Structure, Lattice, PeriodicSite
from pymatgen.util.coord_utils import lattice_points_in_supercell
from pymatgen.vis.structure_vtk import EL_COLORS
from monty.json import MSONable
from monty.os.path import which

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


class StructureGraph(MSONable):

    def __init__(self, structure, graph_data):
        """
        If constructing this class manually, use the `with_empty_graph`
        method or `with_local_env` method (using an algorithm provided
        by the `local_env` module, such as O'Keeffe).

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
        edges of what lattice image the edge belongs to: this is
        always relative to the original lattice, even when graph
        is multiplied to become a supercell. This is not the only
        possible approach, but has benefits, especially in
        terms of graph visualization (edges at boundaries
        become 'dangling bonds'). This approach may change in future.

        :param *args: same as in :class: `pymatgen.core.Structure`
        :param graph_data: dict containing graph information, not
        intended to be constructed manually, contains keys
        "graph_data", "edge_data" (NetworkX 'dict_of_dicts' format),
         and "node_data"
        :param **kwargs: same as in :Class: `pymatgen.core.Structure`
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

        if edge_weight_name and (not edge_weight_units):
            raise ValueError("Please specify units associated "
                             "with your edge weights. Can be "
                             "empty string if arbitrary or "
                             "dimensionless.")

        # construct graph with one node per site
        # graph attributes don't change behavior of graph,
        # they're just for book-keeping
        graph = nx.MultiDiGraph(edge_weight_name=edge_weight_name,
                                edge_weight_units=edge_weight_units,
                                lattice=tuple(structure.lattice.matrix),
                                name=name)
        graph.add_nodes_from(range(len(structure)))

        for n in graph:
            # add species name to graph, useful for visualization
            graph.node[n]['name'] = str(structure[n].specie)
            # add original node fractional coords, useful for creating supercells
            graph.node[n]['frac_coords'] = tuple(structure[n].frac_coords)

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

        sg = StructureGraph.with_empty_graph(structure, name="bonds")

        for n in range(len(structure)):
            neighbors = strategy.get_nn_info(structure, n)
            for neighbor in neighbors:
                sg.add_edge(from_index=n,
                            from_jimage=(0, 0, 0),
                            to_index=neighbor['site_index'],
                            to_jimage=neighbor['image'],
                            weight=neighbor['weight'])

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
                 weight=None):
        """
        Add edge to graph.

        Since physically a 'bond' (or other connection
        between sites) doesn't have a direction, from_index, from_jimage
        can be swapped with to_index, to_jimage. However, images will
        always always be shifted so from_jimage becomes (0, 0, 0).

        :param name: e.g. "bonds"
        :param from_index: index of site connecting from
        :param to_index: index of site connecting to
        :param from_jimage (tuple of ints): lattice vector of periodic
        image, e.g. (1, 0, 0) for periodic image in +x direction
        :param to_jimage (tuple of ints): lattice vector of image
        :param weight (float): e.g. bond length
        :param name: if you have multiple graphs defined, the name of
        the graph you want to use must be supplied
        :return:
        """

        # constrain all from_jimages to be (0, 0, 0),
        # simplifies logic later
        if from_jimage != (0, 0, 0):
            shift = from_jimage
            from_jimage = np.subtract(from_jimage, shift)
            to_jimage = np.subtract(to_jimage, shift)

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
                existing_from = d["from_jimage"]
                existing_to = d["to_jimage"]
                if existing_from == from_jimage and existing_to == to_jimage:
                    warnings.warn("Trying to add an edge that already exists from "
                                  "site {} {} to site {} {}.".format(from_index,
                                                                     from_jimage,
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

    def get_connected_sites(self, n):
        """
        Returns a list of neighbors of site n.
        :param n: index of Site in Structure
        :return: list of Sites that are connected
        to that site
        """
        sites = []
        edges = self.graph.out_edges(n, data=True) + self.graph.in_edges(n, data=True)
        for u, v, d in edges:
            site_d = self.structure[u].as_dict()
            site_d['abc'] = np.add(site_d['abc'], d['to_jimage']).tolist()
            sites.append(PeriodicSite.from_dict(site_d))
        return sites

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
                           hide_unconnected_nodes=False,
                           hide_image_edges=False,
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

        # Developer note: cannot use matplotlib-based
        # plotting since it does not handle MultiGraphs well
        # and can give misleading results. Also, drawing
        # graphs with matplotlib is a bit of a hack, using
        # a dedicated tool like GraphViz allows for much
        # easier control over graph appearance.

        g = self.graph.copy()

        if hide_unconnected_nodes:
            g = g.subgraph([n for n in g.degree() if g.degree()[n] != 0])

        g.graph = {'nodesep': 10.0, 'dpi': 300, 'overlap': "false"}

        # add display options for nodes
        for n in g.nodes():

            # get label by species name
            label = "{}({})".format(g.node[n]['name'], n) if node_labels else ""

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
            from_image = d['from_jimage']
            to_image = d['to_jimage']

            # set edge style
            d['style'] = "solid"
            if from_image != (0, 0, 0) or to_image != (0, 0, 0):
                d['style'] = "dashed"
                if hide_image_edges:
                    edges_to_delete.append([u, v, k])

            # don't show edge directions
            d['arrowhead'] = "none"

            # only add labels for images that are not the origin
            if image_labels:
                d['taillabel'] = "" if from_image == (0, 0, 0) else "from {}".format((from_image))
                d['headlabel'] = "" if to_image == (0, 0, 0) else "to {}".format((to_image))
                d['arrowhead'] = "normal" if (d['taillabel'] or d['headlabel']) else "none"

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
            g.edge[u][v][k].update(d)

        for edge_to_delete in edges_to_delete:
            g.remove_edge(*edge_to_delete)

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

        # code adapted from Structure.__mul__

        warnings.warn("StructureGraph.__mul__ in active development.")

        # TODO: faster implementation, initial implementation for correctness not speed

        scale_matrix = np.array(scaling_matrix, np.int16)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
        new_lattice = Lattice(np.dot(scale_matrix, self.structure.lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        new_graphs = []
        for v in c_lat:

            # create a map of nodes from original graph to its image
            mapping = {n: n + len(new_sites) for n in range(len(self.structure))}

            for site in self.structure:
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
        edges_inside_supercell = []  # sets of {u, v}
        for u, v, k, d in new_g.edges(keys=True, data=True):
            if d["to_jimage"] == (0, 0, 0):
                edges_inside_supercell.append({u, v})

        orig_lattice = Lattice(self.graph.graph['lattice'])
        new_coords = new_structure.cart_coords
        for u, v, k, d in new_g.edges(keys=True, data=True):

            from_jimage = d["from_jimage"]  # for node u
            to_jimage = d["to_jimage"]  # for node v

            # reduce unnecessary checking
            if to_jimage != (0, 0, 0):

                # get fractional co-ordinates of where atoms defined
                # by edge are expected to be, relative to original
                # lattice (keeping original lattice has
                # significant benefits)
                u_image_frac = np.add(new_g.node[u]['frac_coords'], from_jimage)
                v_image_frac = np.add(new_g.node[v]['frac_coords'], to_jimage)
                u_frac = new_g.node[u]['frac_coords']
                v_frac = new_g.node[v]['frac_coords']

                # using the position of node u as a reference,
                # get relative Cartesian co-ordinates of where
                # atoms defined by edge are expected to be
                u_image_cart = orig_lattice.get_cartesian_coords(u_image_frac)
                v_image_cart = orig_lattice.get_cartesian_coords(v_image_frac)
                u_cart = orig_lattice.get_cartesian_coords(u_frac)
                v_cart = orig_lattice.get_cartesian_coords(v_frac)
                u_rel = np.subtract(u_image_cart, u_cart)
                v_rel = np.subtract(v_image_cart, v_cart)

                # now retrieve position of node u (or v) in
                # new supercell, and get absolute Cartesian
                # co-ordinates of where atoms defined by edge
                # are expected to be
                u_expec = new_structure[u].coords + u_rel
                v_expec = new_structure[v].coords + v_rel

                # now search in new structure for these atoms
                # (these lines could/should be optimized)
                v_present = np.where([np.allclose(c, v_expec, atol=0.01) for c in new_coords])[0]

                # sanity check
                if len(v_present) > 1:
                    # could re-write to work in this instance,
                    # but this really shouldn't happen in practice
                    raise Exception("This shouldn't happen, do you have two atoms super-imposed?")

                v_is_present = len(v_present) == 1

                # check if image sites now present in supercell
                # and if so, delete old edge that went through
                # periodic boundary
                if v_is_present:

                    new_u = u
                    new_v = v_present[0]
                    new_d = d.copy()

                    # node now inside supercell
                    new_d['to_jimage'] = (0, 0, 0)

                    edges_to_remove.append((u, v, k))

                    # make sure we don't try to add duplicate edges
                    # will remove two edges for everyone one we add
                    if {new_u, new_v} not in edges_inside_supercell:
                        edges_inside_supercell.append({new_u, new_v})
                        edges_to_add.append((new_u, new_v, new_d))

        logger.debug("Removing {} edges, adding {} new edges.".format(len(edges_to_remove),
                                                                      len(edges_to_add)))

        # add/delete marked edges
        for edges_to_remove in edges_to_remove:
            new_g.remove_edge(*edges_to_remove)
        for (u, v, d) in edges_to_add:
            new_g.add_edge(u, v, **d)
            
        # TODO: re-connect periodic edges?

        # return new instance of StructureGraph with supercell
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "structure": new_structure.as_dict(),
             "graphs": json_graph.adjacency_data(new_g)}

        return StructureGraph.from_dict(d)

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
            header_line += "  {}".format("-"*max([10, len(edge_label)]))
        else:
            print_weights = False

        s = header + "\n" + header_line + "\n"

        if print_weights:
            for u, v, data in g.edges(data=True):
                s += "{:4}  {:4}  {:12}  {:.4E}\n".format(u, v, str(data.get("to_jimage", (0, 0, 0))),
                                                          data.get("weight", 0))
        else:
            for u, v, data in g.edges(data=True):
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
        return len(self.structure)
