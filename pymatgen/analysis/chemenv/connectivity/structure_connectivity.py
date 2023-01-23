"""
Structure connectivity class.
"""

from __future__ import annotations

import collections
import logging

import networkx as nx
import numpy as np
from monty.json import MSONable, jsanitize

from pymatgen.analysis.chemenv.connectivity.connected_components import (
    ConnectedComponent,
)
from pymatgen.analysis.chemenv.connectivity.environment_nodes import (
    get_environment_node,
)
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    LightStructureEnvironments,
)

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "1.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "June 25, 2019"


def get_delta_image(isite1, isite2, data1, data2):
    """
    Helper method to get the delta image between one environment and another
    from the ligand's delta images.
    """
    if data1["start"] == isite1:
        if data2["start"] == isite2:
            return np.array(data1["delta"]) - np.array(data2["delta"])
        return np.array(data1["delta"]) + np.array(data2["delta"])
    if data2["start"] == isite2:
        return -np.array(data1["delta"]) - np.array(data2["delta"])
    return -np.array(data1["delta"]) + np.array(data2["delta"])


class StructureConnectivity(MSONable):
    """
    Main class containing the connectivity of a structure.
    """

    def __init__(
        self,
        light_structure_environment,
        connectivity_graph=None,
        environment_subgraphs=None,
    ):
        """
        Constructor for the StructureConnectivity object.

        Args:
            light_structure_environment: a LightStructureEnvironments object
                                         containing the relevant local environments
                                         for the sites in the structure.
            connectivity_graph: the networkx MultiGraph if it has already been computed,
                                e.g. stored in a file or dict and StructureConnectivity
                                is reconstructed from that file or dict.
            environment_subgraphs: the different subgraphs of environments that have
                                   been computed if any (as for connectivity_graph, only
                                   if it is reconstructed from a file or dict).
        """
        self.light_structure_environments = light_structure_environment
        if connectivity_graph is None:
            self._graph = nx.MultiGraph()
        else:
            self._graph = connectivity_graph
        if environment_subgraphs is None:
            self.environment_subgraphs = {}
        else:
            self.environment_subgraphs = environment_subgraphs

    def environment_subgraph(self, environments_symbols=None, only_atoms=None):
        """
        Args:
            environments_symbols ():
            only_atoms ():

        Returns:
        """
        if environments_symbols is not None:
            self.setup_environment_subgraph(environments_symbols=environments_symbols, only_atoms=only_atoms)
        try:
            return self._environment_subgraph
        except AttributeError:
            all_envs = self.light_structure_environments.environments_identified()
            self.setup_environment_subgraph(environments_symbols=all_envs, only_atoms=only_atoms)
        return self._environment_subgraph

    def add_sites(self):
        """
        Add the sites in the structure connectivity graph.
        """
        self._graph.add_nodes_from(list(range(len(self.light_structure_environments.structure))))

    def add_bonds(self, isite, site_neighbors_set):
        """
        Add the bonds for a given site index to the structure connectivity graph.

        Args:
            isite: Index of the site for which the bonds have to be added.
            site_neighbors_set: site_neighbors_set: Neighbors set of the site
        """
        existing_edges = self._graph.edges(nbunch=[isite], data=True)
        for nb_index_and_image in site_neighbors_set.neighb_indices_and_images:
            nb_index_unitcell = nb_index_and_image["index"]
            nb_image_cell = nb_index_and_image["image_cell"]
            exists = False
            if np.allclose(nb_image_cell, np.zeros(3)):
                for _, ineighb1, data1 in existing_edges:
                    if np.allclose(data1["delta"], np.zeros(3)) and nb_index_unitcell == ineighb1:
                        exists = True
                        break
            else:
                if isite == nb_index_unitcell:
                    for (isite1, ineighb1, data1) in existing_edges:
                        if isite1 == ineighb1:
                            if np.allclose(data1["delta"], nb_image_cell) or np.allclose(
                                data1["delta"], -nb_image_cell
                            ):
                                exists = True
                                break
                else:
                    for _, ineighb1, data1 in existing_edges:
                        if nb_index_unitcell == ineighb1:
                            if data1["start"] == isite:
                                if np.allclose(data1["delta"], nb_image_cell):
                                    exists = True
                                    break
                            elif data1["end"] == isite:
                                if np.allclose(data1["delta"], -nb_image_cell):
                                    exists = True
                                    break
                            else:
                                raise ValueError("SHOULD NOT HAPPEN ???")
            if not exists:
                self._graph.add_edge(
                    isite,
                    nb_index_unitcell,
                    start=isite,
                    end=nb_index_unitcell,
                    delta=nb_image_cell,
                )

    def setup_environment_subgraph(self, environments_symbols, only_atoms=None):
        """
        Set up the graph for predefined environments and optionally atoms.

        Args:
            environments_symbols: Symbols of the environments for the environment subgraph.
            only_atoms: Atoms to be considered.
        """
        logging.info(f"Setup of environment subgraph for environments {', '.join(environments_symbols)}")
        if not isinstance(environments_symbols, collections.abc.Iterable):
            environments_symbols = [environments_symbols]
        environments_symbols = sorted(environments_symbols)
        envs_string = "-".join(environments_symbols)
        if only_atoms is not None:
            envs_string += "#" + "-".join(sorted(only_atoms))
        # Get it directly if it was already computed
        if envs_string in self.environment_subgraphs:
            self._environment_subgraph = self.environment_subgraphs[envs_string]
            return

        # Initialize graph for a subset of environments
        self._environment_subgraph = nx.MultiGraph()
        # Add the sites with the required environment(s)
        for isite, ce_this_site_all in enumerate(self.light_structure_environments.coordination_environments):
            if ce_this_site_all is None:
                continue
            if len(ce_this_site_all) == 0:
                continue
            ce_this_site = ce_this_site_all[0]["ce_symbol"]
            if ce_this_site in environments_symbols:
                if only_atoms is None:
                    env_node = get_environment_node(
                        self.light_structure_environments.structure[isite],
                        isite,
                        ce_this_site,
                    )
                    self._environment_subgraph.add_node(env_node)
                else:
                    if self.light_structure_environments.structure.is_ordered:
                        if self.light_structure_environments.structure[isite].specie.symbol in only_atoms:
                            env_node = get_environment_node(
                                self.light_structure_environments.structure[isite],
                                isite,
                                ce_this_site,
                            )
                            self._environment_subgraph.add_node(env_node)
                    else:
                        #  TODO: add the possibility of a "constraint" on the minimum percentage
                        #        of the atoms on the site
                        this_site_elements = [
                            sp.symbol for sp in self.light_structure_environments.structure[isite].species_and_occu
                        ]
                        for elem_symbol in this_site_elements:
                            if elem_symbol in only_atoms:
                                env_node = get_environment_node(
                                    self.light_structure_environments.structure[isite],
                                    isite,
                                    ce_this_site,
                                )
                                self._environment_subgraph.add_node(env_node)
                                break
        # Find the connections between the environments
        nodes = list(self._environment_subgraph.nodes())
        for inode1, node1 in enumerate(nodes):
            isite1 = node1.isite
            links_node1 = self._graph.edges(isite1, data=True)
            for node2 in nodes[inode1:]:
                isite2 = node2.isite
                links_node2 = self._graph.edges(isite2, data=True)
                # We look for ligands that are common to both site1 and site2
                connections_site1_site2 = {}
                for _, ilig_site1, d1 in links_node1:
                    for _, ilig_site2, d2 in links_node2:
                        if ilig_site1 == ilig_site2:
                            delta_image = get_delta_image(isite1, isite2, d1, d2)
                            if isite1 == isite2 and np.all(delta_image == 0):
                                continue
                            tuple_delta_image = tuple(delta_image)
                            if tuple_delta_image in connections_site1_site2:
                                connections_site1_site2[tuple_delta_image].append((ilig_site1, d1, d2))
                            else:
                                connections_site1_site2[tuple_delta_image] = [(ilig_site1, d1, d2)]
                # Remove the double self-loops ...
                if isite1 == isite2:
                    remove_deltas = []
                    alldeltas = list(connections_site1_site2)
                    alldeltas2 = list(connections_site1_site2)
                    if (0, 0, 0) in alldeltas:
                        alldeltas.remove((0, 0, 0))
                        alldeltas2.remove((0, 0, 0))
                    for current_delta in alldeltas:
                        opp_current_delta = tuple(-dd for dd in current_delta)
                        if opp_current_delta in alldeltas2:
                            remove_deltas.append(current_delta)
                            alldeltas2.remove(current_delta)
                            alldeltas2.remove(opp_current_delta)
                    for remove_delta in remove_deltas:
                        connections_site1_site2.pop(remove_delta)
                # Add all the edges
                for conn, ligands in list(connections_site1_site2.items()):
                    self._environment_subgraph.add_edge(
                        node1,
                        node2,
                        start=node1.isite,
                        end=node2.isite,
                        delta=conn,
                        ligands=ligands,
                    )
        self.environment_subgraphs[envs_string] = self._environment_subgraph

    def setup_connectivity_description(self):
        """
        Returns:
        """

    def get_connected_components(self, environments_symbols=None, only_atoms=None):
        """
        Args:
            environments_symbols ():
            only_atoms ():

        Returns:
        """
        connected_components = []
        env_subgraph = self.environment_subgraph(environments_symbols=environments_symbols, only_atoms=only_atoms)
        for component_nodes in nx.connected_components(env_subgraph):
            graph = env_subgraph.subgraph(component_nodes).copy()
            connected_components.append(ConnectedComponent.from_graph(graph))
        return connected_components

    def setup_atom_environment_subgraph(self, atom_environment):
        """
        Args:
            atom_environment ():

        Returns:
        """
        raise NotImplementedError()

    def setup_environments_subgraph(self, environments_symbols):
        """
        Args:
            environments_symbols ():

        Returns:
        """
        raise NotImplementedError()

    def setup_atom_environments_subgraph(self, atoms_environments):
        """
        Args:
            atoms_environments ():

        Returns:
        """
        raise NotImplementedError()

    def print_links(self):
        """
        Returns:
        """
        nodes = self.environment_subgraph().nodes()
        print("Links in graph :")
        for node in nodes:
            print(node.isite, " is connected with : ")
            for (n1, n2, data) in self.environment_subgraph().edges(node, data=True):
                if n1.isite == data["start"]:
                    print(
                        f"  - {n2.isite} by {len(data['ligands'])} ligands ({data['delta'][0]} "
                        f"{data['delta'][1]} {data['delta'][2]})"
                    )
                else:
                    print(
                        f"  - {n2.isite} by {len(data['ligands'])} ligands ({-data['delta'][0]} "
                        f"{-data['delta'][1]} {-data['delta'][2]})"
                    )

    def as_dict(self):
        """
        Returns:
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "light_structure_environments": self.light_structure_environments.as_dict(),
            "connectivity_graph": jsanitize(nx.to_dict_of_dicts(self._graph)),
            "environment_subgraphs": {
                env_key: jsanitize(nx.to_dict_of_dicts(subgraph))
                for env_key, subgraph in self.environment_subgraphs.items()
            },
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d ():

        Returns:
        """
        # Reconstructs the graph with integer as nodes (json's as_dict replaces integer keys with str keys)
        cgraph = nx.from_dict_of_dicts(d["connectivity_graph"], create_using=nx.MultiGraph, multigraph_input=True)
        cgraph = nx.relabel_nodes(cgraph, int)  # Just relabel the nodes using integer casting (maps str->int)
        # Relabel multiedges (removes multiedges with str keys and adds them back with int keys)
        edges = set(cgraph.edges())
        for n1, n2 in edges:
            new_edges = {int(iedge): edata for iedge, edata in cgraph[n1][n2].items()}
            cgraph.remove_edges_from([(n1, n2, iedge) for iedge, edata in cgraph[n1][n2].items()])
            cgraph.add_edges_from([(n1, n2, iedge, edata) for iedge, edata in new_edges.items()])
        return cls(
            LightStructureEnvironments.from_dict(d["light_structure_environments"]),
            connectivity_graph=cgraph,
            environment_subgraphs=None,
        )
        # TODO: also deserialize the environment_subgraphs
        #            environment_subgraphs={env_key: nx.from_dict_of_dicts(subgraph, multigraph_input=True)
        #                                   for env_key, subgraph in d['environment_subgraphs'].items()})
