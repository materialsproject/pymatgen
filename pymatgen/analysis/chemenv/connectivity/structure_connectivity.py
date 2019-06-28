import networkx as nx
import numpy as np
import collections
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.chemenv.connectivity.environment_nodes import get_environment_node
from pymatgen.analysis.chemenv.connectivity.connected_components import ConnectedComponent
from monty.json import MSONable
from monty.json import jsanitize
import logging

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "1.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "June 25, 2019"


def get_ligand_delta_image(isite1, isite2, data1, data2):
    if data1['start'] == isite1:
        if data2['start'] == isite2:
            return np.array(data1['delta']) - np.array(data2['delta'])
        else:
            return np.array(data1['delta']) + np.array(data2['delta'])
    else:
        if data2['start'] == isite2:
            return -np.array(data1['delta']) - np.array(data2['delta'])
        else:
            return -np.array(data1['delta']) + np.array(data2['delta'])


class StructureConnectivity(MSONable):
    """
    Main class containing the connectivity of a structure.
    """
    def __init__(self, light_structure_environment, connectivity_graph=None):
        """
        Constructore for the StructureConnectivity object.

        :param light_structure_environment: a LightStructureEnvironments object
            containing the relevant local environments for the sites in the
            structure
        :param connectivity_graph:
        """
        self.light_structure_environments = light_structure_environment
        if connectivity_graph is None:
            self._graph = nx.MultiGraph()
        else:
            self._graph = connectivity_graph
        self.environment_subgraphs = {}

    def add_sites(self):
        """
        Add the sites in the structure connectivity graph.
        :return: None
        """
        self._graph.add_nodes_from(list(range(len(self.light_structure_environments.structure))))

    def add_bonds(self, isite, site_neighbors_set):
        """
        Add the bonds for a given site index to the structure connectivity graph.

        :param isite: Index of the site for which the bonds have to be added.
        :param site_neighbors_set: Neighbors set of the site
        :return: None
        """
        existing_edges = self._graph.edges(nbunch=[isite], data=True)
        for nb_index_and_image in site_neighbors_set.neighb_indices_and_images:
            nb_index_unitcell = nb_index_and_image['index']
            nb_image_cell = nb_index_and_image['image_cell']
            exists = False
            if np.allclose(nb_image_cell, np.zeros(3)):
                for (isite1, ineighb1, data1) in existing_edges:
                    if np.allclose(data1['delta'], np.zeros(3)) and nb_index_unitcell == ineighb1:
                        exists = True
                        break
            else:
                if isite == nb_index_unitcell:
                    for (isite1, ineighb1, data1) in existing_edges:
                        if isite1 == ineighb1:
                            if np.allclose(data1['delta'], nb_image_cell) or np.allclose(data1['delta'], -nb_image_cell):
                                exists = True
                                break
                else:
                    for (isite1, ineighb1, data1) in existing_edges:
                        if nb_index_unitcell == ineighb1:
                            if data1['start'] == isite:
                                if np.allclose(data1['delta'], nb_image_cell):
                                    exists = True
                                    break
                            elif data1['end'] == isite:
                                if np.allclose(data1['delta'], -nb_image_cell):
                                    exists = True
                                    break
                            else:
                                raise ValueError('SHOULD NOT HAPPEN ???')
            if not exists:
                self._graph.add_edge(isite, nb_index_unitcell, start=isite, end=nb_index_unitcell, delta=nb_image_cell)

    def setup_environment_subgraph(self, environments_symbols, only_atoms=None):
        logging.info('Setup of environment subgraph for environments {}'.format(', '.join(environments_symbols)))
        if not isinstance(environments_symbols, collections.Iterable):
            environments_symbols = [environments_symbols]
        environments_symbols = sorted(environments_symbols)
        envs_string = '-'.join(environments_symbols)
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
            ce_this_site = ce_this_site_all[0]['ce_symbol']
            if ce_this_site in environments_symbols:
                if only_atoms is None:
                    env_node = get_environment_node(self.light_structure_environments.structure[isite], isite,
                                                    ce_this_site)
                    self._environment_subgraph.add_node(env_node)
                else:
                    if self.light_structure_environments.structure.is_ordered:
                        if self.light_structure_environments.structure[isite].specie.symbol in only_atoms:
                            env_node = get_environment_node(self.light_structure_environments.structure[isite], isite,
                                                            ce_this_site)
                            self._environment_subgraph.add_node(env_node)
                    else:
                        #TODO : add the possibility of a "constraint" on the minimum percentage of the atoms on the site
                        this_site_elements = [sp.symbol for sp in
                                              self.light_structure_environments.structure[isite].species_and_occu]
                        for elem_symbol in this_site_elements:
                            if elem_symbol in only_atoms:
                                env_node = get_environment_node(self.light_structure_environments.structure[isite],
                                                                isite, ce_this_site)
                                self._environment_subgraph.add_node(env_node)
                                break
        # Find the connections between the environments
        nodes = list(self._environment_subgraph.nodes())
        for inode1, node1 in enumerate(nodes):
            isite1 = node1.isite
            links_node1 = self._graph.edges(isite1, data=True)
            for inode2, node2 in enumerate(nodes[inode1:]):
                isite2 = node2.isite
                links_node2 = self._graph.edges(isite2, data=True)
                # We look for ligands that are common to both site1 and site2
                connections_site1_site2 = {}
                for (site1_1, ilig_site1, d1) in links_node1:
                    for (site2_1, ilig_site2, d2) in links_node2:
                        if ilig_site1 == ilig_site2:
                            delta_image = get_ligand_delta_image(isite1, isite2, d1, d2)
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
                    alldeltas = list(connections_site1_site2.keys())
                    alldeltas2 = list(connections_site1_site2.keys())
                    if (0, 0, 0) in alldeltas:
                        alldeltas.remove((0, 0, 0))
                        alldeltas2.remove((0, 0, 0))
                    for current_delta in alldeltas:
                        opp_current_delta = tuple([-dd for dd in current_delta])
                        if opp_current_delta in alldeltas2:
                            remove_deltas.append(current_delta)
                            alldeltas2.remove(current_delta)
                            alldeltas2.remove(opp_current_delta)
                    for remove_delta in remove_deltas:
                        connections_site1_site2.pop(remove_delta)
                # Add all the edges
                for conn, ligands in list(connections_site1_site2.items()):
                    self._environment_subgraph.add_edge(node1, node2, start=node1.isite, end=node2.isite,
                                                        delta=conn, ligands=ligands)
        self.environment_subgraphs[envs_string] = self._environment_subgraph

    def setup_connectivity_description(self):
        pass

    def get_connected_components(self):
        connected_components = []
        for graph in nx.connected_component_subgraphs(self._environment_subgraph):
            connected_components.append(ConnectedComponent.from_graph(graph))
        return connected_components

    def setup_atom_environment_subgraph(self, atom_environment):
        raise NotImplementedError()

    def setup_environments_subgraph(self, environments_symbols):
        raise NotImplementedError()

    def setup_atom_environments_subgraph(self, atoms_environments):
        raise NotImplementedError()

    def print_links(self):
        nodes = self._environment_subgraph.nodes()
        print('Links in graph :')
        for node in nodes:
            print(node.isite, ' is connected with : ')
            for (n1, n2, data) in self._environment_subgraph.edges(node, data=True):
                if n1.isite == data['start']:
                    print('  - {:d} by {:d} ligands ({:d} {:d} {:d})'.format(n2.isite, len(data['ligands']),
                                                                             data['delta'][0], data['delta'][1],
                                                                             data['delta'][2]))
                else:
                    print('  - {:d} by {:d} ligands ({:d} {:d} {:d})'.format(n2.isite, len(data['ligands']),
                                                                             -data['delta'][0], -data['delta'][1],
                                                                             -data['delta'][2]))

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "light_structure_environments": self.light_structure_environments.as_dict(),
                "connectivity_graph": jsanitize(nx.to_dict_of_dicts(self._graph))}

    @classmethod
    def from_dict(cls, d):
        return cls(LightStructureEnvironments.from_dict(d['light_structure_environments']),
                   connectivity_graph=nx.from_dict_of_dicts(d['connectivity_graph'], multigraph_input=True))

#     def add_uc_atom(self, site, isite):
#         """
#
#         :param site:
#         :param isite:
#         :return:
#         """
#         logger.debug('Adding Unit Cell Atom : #{} {}'.format(isite, str(site)))
#         self._graph.add_node(UnitCellNodeSite.from_psite(site, isite))
#
#     def add_uc_bond(self, site_1, site_2, isite_1, isite_2):
#         """
#
#         :param site_1:
#         :param site_2:
#         :param isite_1:
#         :param isite_2:
#         """
#         logger.debug('Adding Unit Cell Bond between sites #{} and #{}'.format(isite_1, isite_2))
#         self._graph.add_edge(UnitCellNodeSite.from_psite(site_1, isite_1), UnitCellNodeSite.from_psite(site_2, isite_2))
#
#     def add_fictive_node_and_bond(self, site, neighb, isite, neighb_index):
#         """
#
#         :param site:
#         :param neighb:
#         :param isite:
#         :param neighb_index:
#         """
#         logger.debug('Adding fictive boundary node between site #{} and '
#                      'its neighbor (#{} in unit cell)'.format(isite, neighb_index))
#         uc_neighb = UnitCellNodeSite.from_psite(neighb.to_unit_cell, neighb_index)
#         delta_uc = np.array(neighb.frac_coords - uc_neighb.frac_coords, np.int)
#         nodesite = UnitCellNodeSite.from_psite(site, isite)
#         fbn = FictiveBoundaryNode(nodesite, uc_neighb, delta_uc, isite, neighb_index)
#         #fbn = FictiveBoundaryNode(site, neighb, delta_uc, isite, neighb_index)
#         #translated_site = PeriodicSite(site._species, site._fcoords - delta_uc, site._lattice,
#         #                               properties=site._properties)
#         #fbn_dual = FictiveBoundaryNode(uc_neighb, translated_site, -delta_uc, neighb_index, isite)
#         self._graph.add_node(fbn)
#         self._graph.add_edge(nodesite, fbn)
#         self._graph.add_node(fbn.dual)
#         self._graph.add_edge(fbn.dual, uc_neighb)
#         #raw_input('number of boundary nodes : {}'.format(len(self.boundary_nodes)))
#
#     def join_dual_boundary_nodes(self):
#         for bn, attr in self.boundary_nodes:
#             self._graph.add_edge(bn, bn.dual)
#
#     # def check_connectivity(self, boundary_node_1, boundary_node_2, connectivity_condition):
#     #     from networkx.algorithms.simple_paths import all_simple_paths
#     #     print len(all_simple_paths(self._graph, boundary_node_1, boundary_node_2))
#
#     def check_connectivity(self):
#         """
#         Checks if all the atoms in the structure are connected to the rest of the structure. This can be used as a
#         check for oxygen atoms for example.
#         :return: True if all the atoms in the structure are connected, False otherwise.
#         """
#         return nx.is_connected(self._graph)
#
#     def setup_environment_subgraph(self, atom_symbol, environment_symbol):
#         """
#         Sets up a subgraph with nodes being the atoms with the given symbol in a given environment and the edges
#         being the ligands (with information on the number of ligands and their type)
#         :param atom_symbol:
#         :param environment_symbol:
#         :return:
#         """
#         logger.info('Setting up environment subgraph')
#         conditions = OrderedDict()
#         conditions['site'] = atom_symbol
#         conditions['environment'] = {'mp_symbol': environment_symbol}
#         self.environment_subgraph = nx.Graph()
#         # Get the sites satisfying the conditions of atom and environment
#         for node, attr in self.unit_cell_nodes:
#             print node, check_node_conditions(conditions, node, node.site_index,
#                                               light_structure_environments=self.light_structure_environment)
#         site_env_list = [node for node, attr in self.unit_cell_nodes
#                          if check_node_conditions(conditions, node, node.site_index,
#                                                   light_structure_environments=self.light_structure_environment)]
#         for node in site_env_list:
#             self.environment_subgraph.add_node(node)
#         #print 'len(environment_subgraph) before', len(self.environment_subgraph)
#         #raw_input('LEN OF ENVIRONMENT SUBGRAPH')
#         #self.environment_subgraph.add_nodes_from(site_env_list)
#         # First loop on atoms satisfying the conditions
#         #for inode, node in enumerate(site_env_list):
#         #    print inode, node, 'has ', len(self.graph.neighbors(node)), ' neighbors :'
#         #    print self.graph.neighbors(node)
#         #    print ''
#         #print 'NEIGHBORS OF EACH V O:6 ....'
#         for inode_1, node_1 in enumerate(site_env_list):
#             logger.info('INODE_1 {} : {}'.format(inode_1, str(node_1)))
#             nnn = 0
#             # Second loop on atoms satisfying the conditions
#             for inode_2, node_2 in enumerate(site_env_list):
#                 logger.debug('INODE_2 {} : {}'.format(inode_2, str(node_2)))
#                 tmpdict = {}
#                 # Loop on the neighbors of the first atom
#                 for inode_nb1, node_nb1 in enumerate(self.graph.neighbors(node_1)):
#                     if node_nb1.node_type == 'unit_cell_node':
#                         nb1 = node_nb1
#                         dnb1 = np.array([0, 0, 0])
#                     elif node_nb1.node_type == 'boundary_node':
#                         nb1 = UnitCellNodeSite.from_psite(node_nb1.neighb_uc, node_nb1.neighb_index)
#                         dnb1 = np.array(node_nb1.delta_uc)
#                         # Loop on the neighbors of the second atom
#                     for inode_nb2, node_nb2 in enumerate(self.graph.neighbors(node_2)):
#                         if node_nb2.node_type == 'unit_cell_node':
#                             nb2 = node_nb2
#                             dnb2 = np.array([0, 0, 0])
#                         elif node_nb2.node_type == 'boundary_node':
#                             nb2 = UnitCellNodeSite.from_psite(node_nb2.neighb_uc, node_nb2.neighb_index)
#                             dnb2 = np.array(node_nb2.delta_uc)
#                         #print '   inode_nb1', inode_nb1, nb1._fcoords, dnb1, '    and    inode_nb2', inode_nb2, nb2._fcoords, dnb2
#
#                         #if nb1 == nb2:
#                         #if np.allclose(nb1.frac_coords, nb2.frac_coords):
#                         #delta = tuple(dnb1 - dnb2)
#                         #print '        => delta : ', delta
#                         #if np.allclose(nb1._fcoords, nb2._fcoords):
#                         if nb1 == nb2:
#                             #print ' nb1 and nb2 : '
#                             #print nb1
#                             #print nb2
#                             #print nb1 == nb2
#                             #raw_input("we have ")
#                             delta = tuple(dnb1 - dnb2)
#                             if delta == (0, 0, 0) and inode_1 == inode_2:
#                                 #print inode_1, inode_2
#                                 #raw_input('These are the same nodes ...')
#                                 continue
#                             if delta in tmpdict:
#                                 #raw_input('delta is in the dict :-)')
#                                 tmpdict[delta].append(nb1)
#                             else:
#                                 #raw_input('adding a delta :-)')
#                                 tmpdict[delta] = [nb1]
#                 #if len(tmpdict) > 0:
#                 #    print 'these are the deltas between nodes {} and {}'.format(inode_1, inode_2)
#                 #    print tmpdict.keys()
#                     #raw_input('...')
#
#                 for delta in tmpdict:
#                     nnn += 1
#                     if delta == (0, 0, 0):
#                         #print 'Adding unit cell edge between nodes ', inode_1, ' and ', inode_2
#                         #print node_1, node_2
#                         self.environment_subgraph.add_edge(node_1, node_2)
#                         #print 'len(environment_subgraph) after', len(self.environment_subgraph)
#                         #raw_input('len...')
#                     else:
#                         npdelta = np.array(delta)
#                         #print 'Adding boundary edge between nodes ', inode_1, ' and ', inode_2, ' with delta = ', npdelta
#                         #print node_1, node_2
#                         fbn = FictiveBoundaryNode(node_1, node_2, npdelta, inode_1, inode_2)
#                         fbn_dual = fbn.dual
#                         self.environment_subgraph.add_node(fbn)
#                         self.environment_subgraph.add_node(fbn_dual)
#                         # Do we realy need the ligands here ?? I dont think so ...
#                         #translated_ligands = [PeriodicSite(lig._species, lig._fcoords - npdelta, lig._lattice,
#                         #                                   properties=lig._properties) for lig in tmpdict[delta]]
#                         self.environment_subgraph.add_edge(node_1, fbn)  #, ligands=tmpdict[delta])
#                         self.environment_subgraph.add_edge(node_2, fbn_dual) #, ligands=translated_ligands)
#                         #print 'len(environment_subgraph) after2', len(self.environment_subgraph)
#                         #raw_input('len...')
#             #print '{} links from previous nodes'.format(nnn)
#         #nx.draw(self.environment_subgraph)
#         #import matplotlib.pyplot as plt
#
#         #plt.axis('off')
#         #plt.show()
#
#     def setup_boundary_graph_and_cycles(self, graph_type='Graph'):
#         if graph_type == 'Graph':
#             self.setup_boundary_graph_and_cycles_graph()
#         elif graph_type == 'DiGraph':
#             self.setup_boundary_graph_and_cycles_digraph()
#         self.boundary_graph_cycles = [cycle for cycle in self.boundary_graph_cycles if self.check_cycle(cycle)]
#
#     def setup_boundary_graph_and_cycles_graph(self):
#         """
#         Sets up the boundary graph representing how the graph is connected from one side of the unit cell to another.
#         """
#         self.boundary_graph = nx.Graph()
#         self.boundary_graph.add_nodes_from([node for node in self.environment_subgraph.nodes()
#                                             if node.node_type == 'boundary_node'])
#         #for node in self.boundary_graph.nodes():
#         #    print node
#         #print 'Number of boundaries of envrionments before : ', len(self.boundary_graph)
#         #raw_input('length of boundary graph = number of nodes .?..')
#         self.boundary_shortest_paths = nx.shortest_path(self.environment_subgraph)
#         for inode_1, node_1 in enumerate(self.boundary_graph.nodes()):
#             for inode_2, node_2 in enumerate(self.boundary_graph.nodes()):
#             #if inode_2 <= inode_1:
#             #    continue
#             #if node_1 == node_2:
#                 #    continue
#                 if node_1 in self.boundary_shortest_paths and node_2 in self.boundary_shortest_paths[node_1]:
#                     # Adding the edges within the unit cell
#                     self.boundary_graph.add_edge(node_1, node_2)
#                     # Adding the edges between boundary nodes
#             fbcn = FictiveBoundaryConnectionNode(node_1, node_1.dual)
#             self.boundary_graph.add_node(fbcn)
#             self.boundary_graph.add_edge(node_1, fbcn)
#             self.boundary_graph.add_edge(fbcn, node_1.dual)
#             #self.boundary_graph.add_edge(node_1, attr_1['dual'], type='boundary_edge')
#         self.boundary_graph_cycles = list(nx.cycle_basis(self.boundary_graph))
#         #print 'Number of boundaries of envrionments after : ', len(self.boundary_graph)
#         #print len(self.boundary_graph_cycles)
#         #self.boundary_graph_cycles = list(nx.simple_cycles(self.boundary_graph))
#
#     def setup_boundary_graph_and_cycles_digraph(self):
#         """
#         Sets up the boundary graph representing how the graph is connected from one side of the unit cell to another.
#         """
#         self.boundary_graph = nx.DiGraph()
#         self.boundary_graph.add_nodes_from([node for node in self.environment_subgraph.nodes()
#                                             if node.node_type == 'boundary_node'])
#         for node in self.boundary_graph.nodes():
#             print node
#         print 'Number of boundaries of envrionments before : ', len(self.boundary_graph)
#         raw_input('length of boundary graph = number of nodes .?..')
#         self.boundary_shortest_paths = nx.shortest_path(self.environment_subgraph)
#         print 'after shortest paths'
#         for inode_1, node_1 in enumerate(self.boundary_graph.nodes()):
#             print 'INODE_1 {}'.format(inode_1)
#             for inode_2, node_2 in enumerate(self.boundary_graph.nodes()):
#                 print 'INODE_2 {}'.format(inode_2)
#             #if inode_2 <= inode_1:
#             #    continue
#             #if node_1 == node_2:
#                 #    continue
#                 if node_1 in self.boundary_shortest_paths and node_2 in self.boundary_shortest_paths[node_1]:
#                     # Adding the edges within the unit cell
#                     self.boundary_graph.add_edge(node_1, node_2)
#                     self.boundary_graph.add_edge(node_2, node_1)
#                     # Adding the edges between boundary nodes
#             fbcn = FictiveBoundaryConnectionNode(node_1, node_1.dual)
#             self.boundary_graph.add_node(fbcn)
#             self.boundary_graph.add_edge(node_1, fbcn)
#             self.boundary_graph.add_edge(fbcn, node_1.dual)
#             self.boundary_graph.add_edge(fbcn, node_1)
#             self.boundary_graph.add_edge(node_1.dual, fbcn)
#             #self.boundary_graph.add_edge(node_1, attr_1['dual'], type='boundary_edge')
#         self.boundary_graph_cycles = list(nx.simple_cycles(self.boundary_graph))
#         print 'Number of boundaries of envrionments after : ', len(self.boundary_graph)
#         print len(self.boundary_graph_cycles)
#         #self.boundary_graph_cycles = list(nx.simple_cycles(self.boundary_graph))
#
#     def get_fiber_unit(self, cycle, current_cell=None, start_index=0, to_periodic_site=True):
#         mycycle = [cycle[ii] for ii in range(start_index, len(cycle))]
#         mycycle.extend([cycle[ii] for ii in range(start_index)])
#         if current_cell is None:
#             current_cell = np.zeros(3, np.int)
#         site_env_list = []
#         for ib1 in range(len(mycycle)):
#             ib2 = np.mod(ib1 + 1, len(mycycle))
#             if mycycle[ib1].node_type == 'boundary_node' and mycycle[ib2].node_type == 'boundary_node':
#                 if (mycycle[ib1] in self.boundary_shortest_paths
#                     and mycycle[ib2] in self.boundary_shortest_paths[mycycle[ib1]]):
#                     if to_periodic_site:
#                         site_list = [PeriodicSite(ps._species, ps._fcoords + current_cell, ps._lattice,
#                                                   properties=ps._properties)
#                                      for ps in self.boundary_shortest_paths[mycycle[ib1]][mycycle[ib2]][1:-1]]
#                     else:
#                         site_list = [UnitCellNodeSite(ps._species, ps._fcoords + current_cell, ps._lattice,
#                                                       ps.site_index, properties=ps._properties)
#                                      for ps in self.boundary_shortest_paths[mycycle[ib1]][mycycle[ib2]][1:-1]]
#                     site_env_list.extend(site_list)
#                     current_cell = current_cell + mycycle[ib2].delta_uc
#         return site_env_list, current_cell
#
#     def get_fiber_unit_old(self, cycle, current_cell=None, inverse=False, start_index=0, to_periodic_site=True):
#         if inverse:
#             mycycle = [cycle[ii] for ii in range(start_index, 0, -1)]
#             mycycle.extend([cycle[ii] for ii in range(len(cycle) - 1, start_index, -1)])
#         else:
#             mycycle = [cycle[ii] for ii in range(start_index, len(cycle))]
#             mycycle.extend([cycle[ii] for ii in range(start_index)])
#         if current_cell is None:
#             current_cell = np.zeros(3, np.int)
#         site_env_list = []
#         for ib1 in range(len(mycycle)):
#             ib2 = np.mod(ib1 + 1, len(mycycle))
#             if mycycle[ib1].node_type == 'boundary_node' and mycycle[ib2].node_type == 'boundary_node':
#                 if (mycycle[ib1] in self.boundary_shortest_paths and
#                             mycycle[ib2] in self.boundary_shortest_paths[mycycle[ib1]]):
#                     if to_periodic_site:
#                         site_list = [PeriodicSite(ps._species, ps._fcoords + current_cell, ps._lattice,
#                                                   properties=ps._properties)
#                                      for ps in self.boundary_shortest_paths[mycycle[ib1]][mycycle[ib2]][1:-1]]
#                     else:
#                         site_list = [ps for ps in self.boundary_shortest_paths[mycycle[ib1]][mycycle[ib2]][1:-1]]
#                     site_env_list.extend(site_list)
#                     current_cell = current_cell + mycycle[ib2].delta_uc
#         return site_env_list, current_cell
#
#     def get_fiber_start_cell(self, cycle, start_index):
#         if start_index == 0:
#             return np.zeros(3, np.int)
#         mycycle = [cycle[ii] for ii in range(start_index, 0, -1)]
#         mycycle.extend([cycle[ii] for ii in range(len(cycle) - 1, start_index, -1)])
#         start_cell = np.zeros(3, np.int)
#         for ib1 in range(len(mycycle)):
#             ib2 = np.mod(ib1 + 1, len(mycycle))
#             if mycycle[ib1].node_type == 'boundary_node' and mycycle[ib2].node_type == 'boundary_node':
#                 if (mycycle[ib1] in self.boundary_shortest_paths
#                     and mycycle[ib2] in self.boundary_shortest_paths[mycycle[ib1]]):
#                     start_cell = start_cell + mycycle[ib2].delta_uc
#                 else:
#                     print 'DUH ?'
#         return start_cell
#
#     def check_cycle(self, cycle):
#         if np.mod(len(cycle), 3) != 0:
#             return False
#         first_node = -1
#         for inode in range(3):
#             if cycle[inode].node_type == 'boundary_connection_node':
#                 first_node = inode
#                 break
#         if first_node < 0:
#             return False
#         current_cell = np.zeros(3, np.int)
#         for ii, inode in enumerate(range(first_node, first_node + len(cycle))):
#             node = cycle[np.mod(inode, len(cycle))]
#             mymodulo = np.mod(ii, 3)
#             if mymodulo == 0:
#                 if node.node_type != 'boundary_connection_node':
#                     return False
#             else:
#                 if node.node_type != 'boundary_node':
#                     return False
#                 if mymodulo == 1:
#                     current_cell += node.delta_uc
#         if not np.allclose(current_cell, np.zeros(3, np.float)):
#             return [first_node]
#         return False
#
#         # found_first = []
#         # for inode in range(len(cycle) + 2):
#         #     node1 = cycle[np.mod(inode, len(cycle))]
#         #     node2 = cycle[np.mod(inode + 1, len(cycle))]
#         #     node3 = cycle[np.mod(inode + 2, len(cycle))]
#         #     if node1.node_type == 'boundary_node' and node2.node_type == 'boundary_node':
#         #         if not found_first:
#         #             found_first.append(inode)
#         #         if not node3.node_type == 'boundary_connection_node':
#         #             return False
#         #     if node1.node_type == 'boundary_connection_node':
#         #         found_fbcn = True
#         #     if inode < len(cycle) and node1.node_type == 'boundary_node' and node2.node_type == 'boundary_node':
#         #         current_cell += node2.delta_uc
#         # if found_fbcn and not np.allclose(current_cell, np.zeros(3, np.float)):
#         #     return found_first
#         # return False
#
#     def check_cycle_old(self, cycle):
#         found_fbcn = False
#         current_cell = np.zeros(3, np.int)
#         found_first = []
#         for inode in range(len(cycle) + 2):
#             node1 = cycle[np.mod(inode, len(cycle))]
#             node2 = cycle[np.mod(inode + 1, len(cycle))]
#             node3 = cycle[np.mod(inode + 2, len(cycle))]
#             if node1.node_type == 'boundary_node' and node2.node_type == 'boundary_node':
#                 if not found_first:
#                     found_first.append(inode)
#                 if not node3.node_type == 'boundary_connection_node':
#                     return False
#             if node1.node_type == 'boundary_connection_node':
#                 found_fbcn = True
#             if inode < len(cycle) and node1.node_type == 'boundary_node' and node2.node_type == 'boundary_node':
#                 current_cell += node2.delta_uc
#         if found_fbcn and not np.allclose(current_cell, np.zeros(3, np.float)):
#             return found_first
#         return False
#
#     def get_fiber_cycles(self, start_cycle=-1, end_cycle=1, to_periodic_site=True):
#         fibers = []
#         for cycle in self.boundary_graph_cycles:
#             checkcycle = self.check_cycle(cycle)
#             if not checkcycle:
#                 continue
#             start_cell = self.get_fiber_start_cell(cycle, start_cycle)
#             #print 'Starting cell is : ', start_cell
#             long_fiber = []
#             current_cell = start_cell
#             for ii in range(end_cycle - start_cycle + 1):
#                 #print ii, 'We are in this cell : ', current_cell
#                 fiber_unit, current_cell = self.get_fiber_unit(cycle, current_cell, start_index=checkcycle[0],
#                                                                to_periodic_site=to_periodic_site)
#                 long_fiber.extend(fiber_unit)
#             fibers.append(list(set(long_fiber)))
#         return fibers
#
#     def get_fiber_cycles_old(self, start_cycle=-1, end_cycle=1, to_periodic_site=True):
#         fibers = []
#         for cycle in self.boundary_graph_cycles:
#             checkcycle = self.check_cycle(cycle)
#             if not checkcycle:
#                 continue
#             start_cell = self.get_fiber_start_cell(cycle, start_cycle)
#             long_fiber = []
#             current_cell = start_cell
#             if start_cycle < 0:
#                 fiber_unit, current_cell = self.get_fiber_unit(cycle, current_cell, inverse=True,
#                                                                start_index=checkcycle[0],
#                                                                to_periodic_site=to_periodic_site)
#                 long_fiber.extend(fiber_unit)
#                 for ii in range(-start_cycle - 1):
#                     fiber_unit, current_cell = self.get_fiber_unit(cycle, current_cell, inverse=True,
#                                                                    start_index=checkcycle[0],
#                                                                    to_periodic_site=to_periodic_site)
#                     long_fiber.extend(fiber_unit)
#                 long_fiber = long_fiber[::-1]
#             if end_cycle > 0:
#                 current_cell = np.zeros(3, np.int)
#                 might_be_duplicate = True
#                 for ii in range(end_cycle):
#                     #print ii,
#                     fiber_unit, current_cell = self.get_fiber_unit(cycle, current_cell, start_index=checkcycle[0],
#                                                                    to_periodic_site=to_periodic_site)
#                     if might_be_duplicate:
#                         for site in fiber_unit:
#                             if not site in long_fiber:
#                                 might_be_duplicate = False
#                                 long_fiber.append(site)
#                     else:
#                         long_fiber.extend(fiber_unit)
#             fibers.append(list(set(long_fiber)))
#         return fibers
#
#     # def setup_boundary_graph_cycles(self):
#     #     """
#     #
#     #
#     #     """
#     #     self.boundary_graph_cycles = list(nx.cycle_basis(self.boundary_graph))
#
#     def draw_boundary_graph(self):
#         import matplotlib.pyplot as plt
#
#         pos = nx.spring_layout(self.boundary_graph)
#         ucnodes = [node for node in self.boundary_graph.nodes() if node.node_type == 'unit_cell_node']
#         boundarynodes = [node for node in self.boundary_graph.nodes() if node.node_type == 'boundary_node']
#         connectivitynodes = [node for node in self.boundary_graph.nodes()
#                              if node.node_type == 'boundary_connection_node']
#         #ucedges = [(e1, e2) for (e1, e2, attr) in self.boundary_graph.edges(data=True) if
#         #           attr['type'] == 'unit_cell_edge']
#         #boundaryedges = [(e1, e2) for (e1, e2, attr) in self.boundary_graph.edges(data=True) if
#         #                 attr['type'] == 'boundary_edge']
#         ucedges = [(e1, e2) for (e1, e2) in self.boundary_graph.edges() if
#                    (e1.node_type == 'boundary_node' and e2.node_type == 'boundary_node')]
#         boundaryedges = [(e1, e2) for (e1, e2) in self.boundary_graph.edges() if
#                          (e1.node_type == 'boundary_connection_node' or e2.node_type == 'boundary_connection_node')]
#         nx.draw_networkx_nodes(self.boundary_graph, pos, nodelist=ucnodes,
#                                node_color='b', node_size=500, alpha=0.8)
#         nx.draw_networkx_nodes(self.boundary_graph, pos, nodelist=boundarynodes,
#                                node_color='g', node_size=500, alpha=0.8)
#         nx.draw_networkx_nodes(self.boundary_graph, pos, nodelist=connectivitynodes,
#                                node_color='r', node_size=500, alpha=0.8)
#         nx.draw_networkx_edges(self.boundary_graph, pos, edgelist=ucedges, width=1.0, alpha=0.5, edge_color='b')
#         nx.draw_networkx_edges(self.boundary_graph, pos, edgelist=boundaryedges, width=1.0, alpha=0.5, edge_color='r')
#         plt.axis('off')
#         plt.show()
#         #raise NotImplementedError
#
#     #     #!/usr/bin/env python
#     # """
#     # Draw a graph with matplotlib, color by degree.
#     #
#     # You must have matplotlib for this to work.
#     # """
#     # __author__ = """Aric Hagberg (hagberg@lanl.gov)"""
#     # import matplotlib.pyplot as plt
#     #
#     # import networkx as nx
#     #
#     # G=nx.cubical_graph()
#     # pos=nx.spring_layout(G) # positions for all nodes
#     #
#     # # nodes
#     # nx.draw_networkx_nodes(G,pos,
#     #                        nodelist=[0,1,2,3],
#     #                        node_color='r',
#     #                        node_size=500,
#     #                    alpha=0.8)
#     # nx.draw_networkx_nodes(G,pos,
#     #                        nodelist=[4,5,6,7],
#     #                        node_color='b',
#     #                        node_size=500,
#     #                    alpha=0.8)
#     #
#     # # edges
#     # nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
#     # nx.draw_networkx_edges(G,pos,
#     #                        edgelist=[(0,1),(1,2),(2,3),(3,0)],
#     #                        width=8,alpha=0.5,edge_color='r')
#     # nx.draw_networkx_edges(G,pos,
#     #                        edgelist=[(4,5),(5,6),(6,7),(7,4)],
#     #                        width=8,alpha=0.5,edge_color='b')
#     #
#     #
#     # # some math labels
#     # labels={}
#     # labels[0]=r'$a$'
#     # labels[1]=r'$b$'
#     # labels[2]=r'$c$'
#     # labels[3]=r'$d$'
#     # labels[4]=r'$\alpha$'
#     # labels[5]=r'$\beta$'
#     # labels[6]=r'$\gamma$'
#     # labels[7]=r'$\delta$'
#     # nx.draw_networkx_labels(G,pos,labels,font_size=16)
#     #
#     # plt.axis('off')
#     # plt.savefig("labels_and_colors.png") # save as png
#     # plt.show() # display
#
#
#     """
#     def setup_directed_boundary_graph(self):
#
#         Sets up the boundary graph representing how the graph is connected from one side of the unit cell to another.
#
#         self.boundary_graph = nx.DiGraph()
#         self.boundary_graph.add_nodes_from([(node, attr) for node, attr in self.environment_subgraph.nodes(data=True)
#                                             if attr['type'] == 'boundary_node'])
#         shortest_paths = nx.shortest_path(self.environment_subgraph)
#         for node_1, attr_1 in self.boundary_graph.nodes(data=True):
#             for node_2 in self.boundary_graph.nodes():
#                 if node_1 == node_2:
#                     continue
#                 if node_1 in shortest_paths and node_2 in shortest_paths[node_1]:
#                     # Adding the edges within the unit cell
#                     self.boundary_graph.add_edge(node_1, node_2)
#                     self.boundary_graph.add_edge(node_2, node_1)
#             # Adding the edges between boundary nodes
#             self.boundary_graph.add_edge(node_1, attr_1['dual'])
#             self.boundary_graph.add_edge(attr_1['dual'], node_1)
#         cycles = list(nx.simple_cycles(self.boundary_graph))
#         for cycle in cycles:
#             print 'Here is a cycle of length ', len(cycle)
#             for nn in cycle:
#             #    if
#                 print nn
#             for nn in range(len(cycle)):
#                 for nn2 in range(len(cycle)):
#                     print nn,nn2, cycle[nn].is_dual_of(cycle[nn2])
#                 #print cycle[nn].is_dual_of(cycle[np.mod(len(cycle)-nn, len(cycle))])
#             #print cycle[0].is_dual_of(cycle[-1])
#             #print cycle[0] == cycle[-1]
#             print '----------------------------'
#     """
#
#     @property
#     def graph(self):
#         """
#         Returns the Graph object (networkx.Graph) representing the atoms and their bonds.
#         :return: networkx.Graph object representing the atoms and their bonds
#         """
#         return self._graph
#
#     @property
#     def boundary_nodes(self, data=True):
#         """
#         Returns a list of all the boundary nodes in the Graph object representing the atoms and their bonds.
#         :return: list of all the boundary nodes in the Graph object representing the atoms and their bonds.
#         """
#         if data:
#             return [(node, attr) for node, attr in self._graph.nodes(data=True) if node.node_type == "boundary_node"]
#         else:
#             return [node for node in self._graph.nodes() if node.node_type == "boundary_node"]
#
#     @property
#     def unit_cell_nodes(self, data=True):
#         """
#         Returns a list of all the unit cell nodes in the Graph object representing the atoms and their bonds.
#         :return: list of all the unit cell nodes in the Graph object representing the atoms and their bonds.
#         """
#         if data:
#             return [(node, attr) for node, attr in self.graph.nodes(data=True) if node.node_type == "unit_cell_node"]
#         else:
#             return [node for node in self.graph.nodes() if node.node_type == "unit_cell_node"]
#
#     @property
#     def edges(self, data=True):
#         if data:
#             return self.graph.edges(data=True)
#         else:
#             return self.graph.edges()
#
#
# class UnitCellNodeSite(PeriodicSite, MSONable):
#     def __init__(self, atoms_n_occu, coords, lattice, site_index, to_unit_cell=False,
#                  coords_are_cartesian=False, properties=None):
#         PeriodicSite.__init__(self, atoms_n_occu, coords, lattice, to_unit_cell,
#                               coords_are_cartesian, properties=properties)
#         self.node_type = 'unit_cell_node'
#         self.site_index = site_index
#
#     @classmethod
#     def from_psite(cls, psite, site_index=None):
#         #return cls(psite._species, psite._fcoords, psite._lattice, properties=psite._properties)
#         if site_index is not None:
#             return cls(psite._species, psite._fcoords, psite._lattice, site_index, properties=psite._properties)
#         else:
#             raise ValueError("site_index is None ...")
#
#     def to_periodic_site(self):
#         return PeriodicSite(self._species, self._fcoords, self._lattice, to_unit_cell=False,
#                             coords_are_cartesian=False, properties=self._properties)
#
#     def __eq__(self, other):
#         return self.node_type == other.node_type and self._site_eq(other) #super(UnitCellNodeSite, self).__eq__(other)
#
#     def __cmp__(self, other):
#         return self.__eq__(other) #super(UnitCellNodeSite, self).__eq__(other)
#
#     def __hash__(self):
#         return super(UnitCellNodeSite, self).__hash__()
#
#     def _site_eq(self, other):
#         return (self._species == other._species and self._lattice == other._lattice and
#                 np.allclose(self._fcoords, other._fcoords, atol=Site.position_atol))
#
#     def as_dict(self):
#         return {"psite": super(UnitCellNodeSite, self).as_dict(),
#                 "node_type": self.node_type}
#
#     @classmethod
#     def from_dict(cls, d):
#         return cls(PeriodicSite.from_dict(d["psite"]), d["node_type"])
#
#
# def check_node_conditions(conditions, node, isite, light_structure_environments=None):
#     """
#
#
#     :param isite:
#     :param conditions:
#     :param node:
#     :param structure_environments:
#     """
#     #if len(node.species_and_occu) != 1:
#     #    raise NotImplementedError('Does not allow unordered structures (yet)')
#     for cond in conditions:
#         #print 'checking CONDITION : ', cond
#         if cond == 'site':
#             species = node.species_and_occu.keys()
#             symbols = [sp.symbol for sp in species]
#             if conditions['site'] not in symbols:
#                 return False
#         elif cond == 'environment':
#             env = conditions['environment']['mp_symbol']
#             if not light_structure_environments.site_contains_environment(isite, env):
#                 return False
#     return True
#
#
# class FictiveBoundaryNode(MSONable):
#     """
#     Class representing the boundary nodes used in the connectivity identification algorithm
#     """
#
#     def __init__(self, site_uc, neighb_uc, delta_uc, isite_uc, neighb_index):
#         """
#
#         :param site_uc:
#         :param neighb:
#         :param delta_uc:
#         :param isite_uc:
#         :param neighb_index:
#         """
#         self.site_uc = site_uc
#         self.neighb_uc = neighb_uc
#         self.isite_uc = isite_uc
#         self.neighb_index = neighb_index
#         self.delta_uc = delta_uc
#         self.node_type = 'boundary_node'
#
#     def is_dual_of(self, other):
#         return (all(self.delta_uc == -other.delta_uc) and self.site_uc == other.neighb_uc and
#                 self.neighb_uc == other.site_uc)
#
#     @property
#     def dual(self):
#         return FictiveBoundaryNode(self.neighb_uc, self.site_uc, -self.delta_uc, self.neighb_index, self.isite_uc)
#
#     def __str__(self):
#         out = 'Fictive Boundary Node between the following sites :\n'
#         out += '   - {}'.format(self.site_uc.species_string)
#         out += ' ({}, {}, {})\n'.format(self.site_uc._fcoords[0], self.site_uc._fcoords[1], self.site_uc._fcoords[2])
#         out += '   - {}'.format(self.neighb_uc.species_string)
#         out += ' ({}, {}, {})\n'.format(self.neighb_uc._fcoords[0], self.neighb_uc._fcoords[1],
#                                         self.neighb_uc._fcoords[2])
#         out += ' with delta_uc = {} {} {}'.format(str(self.delta_uc[0]), str(self.delta_uc[1]), str(self.delta_uc[2]))
#         return out
#
#     def __hash__(self):
#         hashcode = int(sum((el.Z * occu for el, occu in self.site_uc._species.items())) + \
#                        sum((el.Z * occu for el, occu in self.neighb_uc._species.items())))
#         return hashcode
#
#     def __cmp__(self, other):
#         return self.__eq__(other)
#
#     def __eq__(self, other):
#         return (other.node_type == self.node_type and
#                 all(self.delta_uc == other.delta_uc) and
#                 self.site_uc == other.site_uc and
#                 self.neighb_uc == other.neighb_uc)
#
#     def as_dict(self):
#         return {"site_uc": self.site_uc.as_dict(),
#                 "neighb_uc": self.neighb_uc.as_dict(),
#                 "isite_uc": self.isite_uc,
#                 "neighb_index": self.neighb_index,
#                 "delta_uc": self.delta_uc}
#
#     @classmethod
#     def from_dict(cls, d):
#         dec = PMGJSONDecoder()
#         return cls(dec.process_decoded(d["site_uc"]), dec.process_decoded(d["neighb_uc"]),
#                    d["isite_uc"], d["neighb_index"], d["delta_uc"])
#
#
# class FictiveBoundaryConnectionNode(MSONable):
#     """
#     Class representing the fictive boundary connection nodes used in the connectivity identification algorithm. This
#     fictive node is needed to distinguish the unit cell edges in the boundary_graph from the boundary edges (periodic
#     edges)
#     """
#
#     def __init__(self, fbn_1, fbn_2):
#         """
#         Constructor for the FictiveBoundaryConnectionNode
#         :param fbn_1: FictiveBoundaryNode
#         :param fbn_2: FictiveBoundaryNode
#         """
#         self.fbn_1 = fbn_1
#         self.fbn_2 = fbn_2
#         self.node_type = 'boundary_connection_node'
#
#     def __hash__(self):
#         return hash(self.fbn_1) * 100 + hash(self.fbn_2)
#
#     def __cmp__(self, other):
#         return ((self.fbn_1 == other.fbn_1 and self.fbn_2 == other.fbn_2) or
#                 (self.fbn_1 == other.fbn_2 and self.fbn_2 == other.fbn_1))
#
#     def __eq__(self, other):
#         return (other.__class__.__name__ == self.__class__.__name__ and
#                 ((self.fbn_1 == other.fbn_1 and self.fbn_2 == other.fbn_2) or
#                  (self.fbn_1 == other.fbn_2 and self.fbn_2 == other.fbn_1)))
#
#     def as_dict(self):
#         return {"fbn_1": self.fbn_1.as_dict(),
#                 "fbn_2": self.fbn_2.as_dict()}
#
#     @classmethod
#     def from_dict(cls, d):
#         return cls(d["fbn_1"], d["fbn_2"])
#
#
# class AbstractConnectivityCondition():
#     """
#     Abstract class for connectivity conditions.
#     """
#     __metaclass__ = abc.ABCMeta
#
#     @abc.abstractmethod
#     def check(self, path):
#         """
#         Checks the condition on the path given. Should return True if the condition can be satisfied
#         :param path: path to be checked
#         """
#         return
#
#
# class LinkedEnvironmentsCondition(AbstractConnectivityCondition):
#     """
#     Class to check the conditions of linked environments (such as linked vanadium octahedrons)
#     """
#
#     def __init__(self, atom_symbol, environment_symbol, strict=True):
#         """
#
#         :param atom_symbol:
#         :param environment_symbol:
#         :param strict:
#         :raise:
#         """
#         if not strict:
#             raise NotImplementedError("One should allow to check for mixed environments, or mixed atoms")
#         self.atom_symbol = atom_symbol
#         self.environment_symbol = environment_symbol
#         self.conditions = deque([{'site': atom_symbol, 'environment': environment_symbol}, {}])
#
#     def check(self, path, structure, light_structure_environments):
#         """
#
#         :param path:
#         :param structure:
#         :param structure_environments:
#         :return:
#         """
#         if len(path) > len(self.conditions):
#             inodes_start = []
#             #Find the nodes that comply with the first condition
#             for inode, node in enumerate(path):
#                 if check_node_conditions(self.conditions[0], node,
#                                          light_structure_environments=light_structure_environments):
#                     inodes_start.append(inode)
#                     #Test the path
#             for inodestart in inodes_start:
#                 found = True
#                 #Check the head
#                 for inode in range(inodestart):
#                     if not check_node_conditions(self.conditions[np.mod(inode, len(self.conditions))]):
#                         found = False
#                         break
#                         #Check the tail
#                 if found:
#                     for inode in range(inodestart, len(path)):
#                         if not check_node_conditions(self.conditions[np.mod(inode - inodestart, len(self.conditions))]):
#                             found = False
#                             break
#                 if found: return True
#             pass
#         elif len(path) < len(self.conditions):
#             for nn in range(len(self.conditions)):
#                 found = True
#                 for inode, node in enumerate(path):
#                     if not check_node_conditions(self.conditions[np.mod(nn + inode, len(self.conditions))]):
#                         found = False
#                         break
#                 if found: return True
#             return False
#         else:
#             for nn in range(len(self.conditions)):
#                 found = True
#                 for inode, node in enumerate(path):
#                     if not check_node_conditions(self.conditions[inode], node,
#                                                  light_structure_environments=light_structure_environments):
#                         found = False
#                         break
#                 if found: return True
#                 self.conditions.rotate(1)
#             return False