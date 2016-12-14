from pymatgen.analysis.chemenv.utils.graph_utils import get_all_simple_paths_edges

__author__ = 'waroquiers'

from monty.json import MSONable
from pymatgen.analysis.chemenv.utils.chemenv_errors import ChemenvError
import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np
from pymatgen.analysis.chemenv.utils.math_utils import divisors, get_linearly_independent_vectors
import itertools
from matplotlib.patches import FancyArrowPatch, Circle


def draw_network(env_graph, pos, ax, sg=None, periodicity_vectors=None):
    for n in env_graph:
        c = Circle(pos[n], radius=0.02, alpha=0.5)
        ax.add_patch(c)
        env_graph.node[n]['patch'] = c
        x, y = pos[n]
        ax.annotate(n, pos[n], ha='center', va='center', xycoords='data')
    seen = {}
    for (u, v, d) in env_graph.edges(data=True):
        n1 = env_graph.node[u]['patch']
        n2 = env_graph.node[v]['patch']
        rad = 0.1
        if (u, v) in seen:
            rad = seen.get((u, v))
            rad = (rad + np.sign(rad) * 0.1) * -1
        alpha = 0.5
        color = 'k'
        periodicity_colors = ['r', 'g', 'b']

        delta = get_delta(u, v, d)

        #center = get_center_of_arc(n1.center, n2.center, rad)
        n1center = np.array(n1.center)
        n2center = np.array(n2.center)
        midpoint = (n1center + n2center) / 2
        dist = np.sqrt(np.power(n2.center[0] - n1.center[0], 2) + np.power(n2.center[1] - n1.center[1], 2))
        n1c_to_n2c = n2center - n1center
        vv = np.cross(np.array([n1c_to_n2c[0], n1c_to_n2c[1], 0], np.float), np.array([0, 0, 1], np.float))
        vv /= np.linalg.norm(vv)
        midarc = midpoint + rad * dist * np.array([vv[0], vv[1]], np.float)
        xytext_offset = 0.1 * dist * np.array([vv[0], vv[1]], np.float)

        if periodicity_vectors is not None and len(periodicity_vectors) == 1:
            if np.all(np.array(delta) ==
                    np.array(periodicity_vectors[0])) or np.all(np.array(delta) ==
                    -np.array(periodicity_vectors[0])):
                e = FancyArrowPatch(n1center, n2center, patchA=n1, patchB=n2,
                                    arrowstyle='-|>',
                                    connectionstyle='arc3,rad=%s' % rad,
                                    mutation_scale=15.0,
                                    lw=2,
                                    alpha=alpha,
                                    color='r',
                                    linestyle='dashed')
            else:
                e = FancyArrowPatch(n1center, n2center, patchA=n1, patchB=n2,
                                    arrowstyle='-|>',
                                    connectionstyle='arc3,rad=%s' % rad,
                                    mutation_scale=10.0,
                                    lw=2,
                                    alpha=alpha,
                                    color=color)
        else:
            e = FancyArrowPatch(n1center, n2center, patchA=n1, patchB=n2,
                                arrowstyle='-|>',
                                connectionstyle='arc3,rad=%s' % rad,
                                mutation_scale=10.0,
                                lw=2,
                                alpha=alpha,
                                color=color)

        ax.annotate(delta, midarc, ha='center', va='center', xycoords='data', xytext=xytext_offset,
                    textcoords='offset points')
        seen[(u, v)] = rad
        ax.add_patch(e)

    return e


def get_delta(node1, node2, edge_data):
    if node1.isite == edge_data['start'] and node2.isite == edge_data['end']:
        return np.array(edge_data['delta'])
    elif node2.isite == edge_data['start'] and node1.isite == edge_data['end']:
        return -np.array(edge_data['delta'])
    else:
        raise ValueError("Trying to find a delta between two nodes with an edge that seem not to link these nodes")

#
# def get_delta(node1, node2, edge_data):
#     if node1 == edge_data['start'] and node2 == edge_data['end']:
#         return np.array(edge_data['delta'])
#     elif node2 == edge_data['start'] and node1 == edge_data['end']:
#         return -np.array(edge_data['delta'])
#     else:
#         raise ValueError("Trying to find a delta between two nodes with an edge that seem not to link these nodes")


def get_ordered_path_isites(path):
    i_smallest = np.argmin(path)
    if path[np.mod(i_smallest + 1, len(path))] > path[np.mod(i_smallest - 1, len(path))]:
        return tuple([path[np.mod(ii, len(path))] for ii in range(i_smallest, i_smallest + len(path))])
    else:
        return tuple([path[np.mod(ii, len(path))] for ii in range(i_smallest, i_smallest - len(path), - 1)])


def get_ordered_node_group(node_group):
    min_groups = [np.min(gg) for gg in node_group]
    isorted = np.argsort(min_groups)
    return tuple([tuple(sorted(node_group[ii])) for ii in isorted])


def all_pairs_combinations(even_length_list, return_indices=False):
    indices_list = list(range(len(even_length_list)))
    pairs_combinations = []
    groups = []
    opposite_groups = []
    for group in itertools.combinations(indices_list, len(even_length_list) / 2):
        opposite_group = tuple(set(indices_list) - set(group))
        if group not in groups and opposite_group not in groups:
            groups.append(group)
            opposite_groups.append(opposite_group)
    for igroup, group in enumerate(groups):
        for group_perm in itertools.permutations(opposite_groups[igroup]):
            combination = tuple(tuple(sorted([group[ii], group_perm[ii]])) for ii in range(len(group)))
            if not combination in pairs_combinations:
                pairs_combinations.append(combination)
    if return_indices:
        return pairs_combinations
    return [[(even_length_list[pair[0]],
              even_length_list[pair[1]]) for pair in pair_combi] for pair_combi in pairs_combinations]


def cycle_contains_edge(cycle, edge):
    found = 0
    for cycle_edge in cycle:
        if cycle_edge[0] == edge[0] and cycle_edge[1] == edge[1] and cycle_edge[2] == edge[2]:
            found = 1
            break
        elif cycle_edge[0] == edge[1] and cycle_edge[1] == edge[0] and cycle_edge[2] == edge[2]:
            found = -1
            break
    if found == 0:
        return False
    delta = np.zeros(3, np.int)
    for n1, n2, key, data in cycle:
        delta += get_delta(n1, n2, data)
    return tuple(found*delta)


def make_supergraph(graph, multiplicity, periodicity_vectors):
    supergraph = nx.MultiGraph()
    print('peridoicity vectors :')
    print(periodicity_vectors)
    if isinstance(multiplicity, int) or len(multiplicity) == 1:
        mult = multiplicity if isinstance(multiplicity, int) else multiplicity[0]
        nodes = graph.nodes(data=True)
        inodes = [isite for isite, data in nodes]
        indices_nodes = {isite: inodes.index(isite) for isite in inodes}
        edges = graph.edges(data=True, keys=True)
        connecting_edges = []
        other_edges = []
        for (n1, n2, key, data) in edges:
            print(n1, n2, key, data)
            if np.all(np.array(data['delta']) == np.array(periodicity_vectors[0])):
                connecting_edges.append((n1, n2, key, data))
            elif np.all(np.array(data['delta']) == -np.array(periodicity_vectors[0])):
                new_data = dict(data)
                new_data['delta'] = tuple(-np.array(data['delta']))
                new_data['start'] = data['end']
                new_data['end'] = data['start']
                connecting_edges.append((n1, n2, key, new_data))
            else:
                if not np.all(np.array(data['delta']) == 0):
                    print('delta not equal to periodicity nor 0 ... : ', n1, n2, key, data['delta'], data)
                    input('Are we ok with this ?')
                other_edges.append((n1, n2, key, data))
        # for imult in range(mult):
        #     for n1, n2, key, data in other_edges:
        #         new_data = dict(data)
        #         new_data['start'] = (imult*len(nodes)) + indices_nodes[n1]
        #         new_data['end'] = (imult*len(nodes)) + indices_nodes[n2]
        #         supergraph.add_edge(new_data['start'], new_data['end'],
        #                             key=key, attr_dict=new_data)
        #     for n1, n2, key, data in connecting_edges:
        #         new_data = dict(data)
        #         new_data['start'] = (imult*len(nodes)) + indices_nodes[n1]
        #         new_data['end'] = (imult*len(nodes)) + indices_nodes[n2]
        #         #new_data['delta'] = (0, 0, 0)
        #         supergraph.add_edge(new_data['start'], new_data['end'],
        #                             key=key, attr_dict=new_data)
        print(periodicity_vectors)
        for imult in range(mult-1):
            for n1, n2, key, data in other_edges:
                new_data = dict(data)
                new_data['start'] = (imult*len(nodes)) + indices_nodes[n1]
                new_data['end'] = (imult*len(nodes)) + indices_nodes[n2]
                supergraph.add_edge(new_data['start'], new_data['end'],
                                    key=key, attr_dict=new_data)
            for n1, n2, key, data in connecting_edges:
                new_data = dict(data)
                new_data['start'] = (imult*len(nodes)) + indices_nodes[n1]
                new_data['end'] = np.mod(((imult+1)*len(nodes)) + indices_nodes[n2], len(nodes)*mult)
                new_data['delta'] = (0, 0, 0)
                supergraph.add_edge(new_data['start'], new_data['end'],
                                    key=key, attr_dict=new_data)
        imult = mult-1
        for n1, n2, key, data in other_edges:
            new_data = dict(data)
            new_data['start'] = (imult*len(nodes)) + indices_nodes[n1]
            new_data['end'] = (imult*len(nodes)) + indices_nodes[n2]
            supergraph.add_edge(new_data['start'], new_data['end'],
                                key=key, attr_dict=new_data)
        for n1, n2, key, data in connecting_edges:
            new_data = dict(data)
            new_data['start'] = (imult*len(nodes)) + indices_nodes[n1]
            new_data['end'] = indices_nodes[n2]
            supergraph.add_edge(new_data['start'], new_data['end'],
                                key=key, attr_dict=new_data)
        return supergraph
    else:
        raise NotImplementedError('make_supergraph not yet implemented for 2- and 3-periodic graphs')


class ConnectedComponent(MSONable):
    """
    Class used to describe the connected components in a structure in terms of coordination environments
    """

    def __init__(self, environments=None, links=None, environments_data=None, links_data=None, graph=None):
        """

        :param environments: list of environments
        :param links:
        :param environments_data:
        :param links_data:
        """
        self._periodicity_vectors = None
        self._primitive_reduced_connected_subgraph = None
        self._projected = False
        if graph is None:
            self._connected_subgraph = nx.MultiGraph()
            if environments_data is None:
                self._connected_subgraph.add_nodes_from(environments)
            else:
                self._connected_subgraph.add_nodes_from(environments, environments_data)
            for (env_node1, env_node2) in links:
                if ((not self._connected_subgraph.has_node(env_node1)) or
                        (not self._connected_subgraph.has_node(env_node2))):
                    raise ChemenvError(self.__class__, '__init__', 'Trying to add edge with some unexisting node ...')
                self._connected_subgraph.add_edge(env_node1, env_node2, attr_dict=links_data)
        else:
            self._connected_subgraph = graph

    def compute_periodicity(self):
        self_loop_nodes = self._connected_subgraph.nodes_with_selfloops()
        all_nodes_independent_cell_image_vectors = []
        for test_node in self._connected_subgraph.nodes():
            #TODO: do we need to go through all test nodes ?
            this_node_cell_img_vectors = []
            if test_node in self_loop_nodes:
                for key, edge_data in self._connected_subgraph[test_node][test_node].items():
                    this_node_cell_img_vectors.append(edge_data['delta'])
            # Here, we adopt a cutoff equal to the size of the graph, contrary to the default of networkX (size - 1),
            # because otherwise, the all_simple_paths algorithm fail when the source node is equal to the target node.
            paths = []
            #TODO: its probably possible to do just a dfs or bfs traversal instead of taking all simple paths!
            for path in nx.all_simple_paths(self._connected_subgraph, test_node, test_node,
                                            cutoff=len(self._connected_subgraph)):
                if path not in paths:
                    paths.append(path)
                else:
                    continue
                # TODO: there are some paths that appears twice for cycles, and there are some paths that should
                # probably not be considered
                this_path_deltas = [np.zeros(3, np.int)]
                for (node1, node2) in [(node1, path[inode1 + 1]) for inode1, node1 in enumerate(path[:-1])]:
                    this_path_deltas_new = []
                    for key, edge_data in self._connected_subgraph[node1][node2].items():
                        for current_delta in this_path_deltas:
                            delta = get_delta(node1, node2, edge_data)
                            this_path_deltas_new.append(current_delta + delta)
                    this_path_deltas = this_path_deltas_new
                this_node_cell_img_vectors.extend(this_path_deltas)
                this_node_cell_img_vectors = get_linearly_independent_vectors(this_node_cell_img_vectors)
                if len(this_node_cell_img_vectors) == 3:
                    break
            #independent_cell_img_vectors = get_linearly_independent_vectors(this_node_cell_img_vectors)
            independent_cell_img_vectors = this_node_cell_img_vectors
            all_nodes_independent_cell_image_vectors.append(independent_cell_img_vectors)
            #If we have found that the sub structure network is 3D-connected, we can stop ...
            if len(independent_cell_img_vectors) == 3:
                break
        self._periodicity_vectors = []
        if len(all_nodes_independent_cell_image_vectors) != 0:
            for independent_cell_img_vectors in all_nodes_independent_cell_image_vectors:
                if len(independent_cell_img_vectors) > len(self._periodicity_vectors):
                    self._periodicity_vectors = independent_cell_img_vectors
                if len(self._periodicity_vectors) == 3:
                    break

    def primitive_reduce(self):
        """
        Reduces the periodic subgraph to a primitive subgraph
        """
        if self.is_0d:
            return self._connected_subgraph
        elif self.is_1d:
            self.primitive_reduce_1d()
        else:
            print('WARNING: primitive reduce not yet available for 2d and 3d periodicities, the full (non reduced)' \
                  'graph will be provided')
            return self._connected_subgraph

    def primitive_reduce_1d(self, algo='node_connectivity'):
        if algo == 'subgraph_isomorphism':
            self.primitive_reduce_1d_subgraph_isomorphism()
        elif algo == 'node_connectivity':
            self.primitive_reduce_1d_node_connectivity()
        else:
            raise ValueError('"{}" is not a vaild algorithm for primitive_reduce_1d')

    def primitive_reduce_1d_subgraph_isomorphism(self):
        #TODO: look at the concept of "centrality" in networkx for the convention of the "primitive cell"
        self.project_1d()
        _divisors = divisors(len(self._connected_subgraph))
        _nodes_and_data = self._connected_subgraph.nodes(data=True)
        isites_nodes = sorted([nn for nn, dd in _nodes_and_data])
        all_node_groups = {}
        for divisor in _divisors[:-1]:
            all_node_groups[divisor] = []
            for node_bunch in itertools.combinations(_nodes_and_data, divisor):
                subgraph = self._connected_subgraph.subgraph([node for node, data in node_bunch])
                if not nx.is_connected(subgraph):  #This assumption should maybe be checked but I think its ok
                    continue
                gm = isomorphism.MultiGraphMatcher(self._connected_subgraph, subgraph)
                node_groups = []
                existing = []
                print('All node groups')
                for subgraph_isomorphism in gm.subgraph_isomorphisms_iter():
                    print([nn for nn in list(subgraph_isomorphism.keys())])
                    subgraph_isites = tuple(sorted([nn for nn in list(subgraph_isomorphism.keys())]))
                    if subgraph_isites in existing:
                        continue
                    existing.append(subgraph_isites)
                    node_groups.append(list(subgraph_isomorphism.keys()))
                print('Independent node groups : ')
                for node_group in node_groups:
                    print([nn for nn in node_group])
                for node, data in node_bunch:
                    print(node, node.central_site.frac_coords)
                print('')
                n_groups = len(self._connected_subgraph) / divisor
                if len(node_groups) < n_groups:
                    print('Not enough groups ... no possibility')
                else:
                    print('Number of groups equal or larger to len(graph) / divisor ... find all the ' \
                          'combinations of groups ' \
                          'for which the list of all isites gives the isites of the connected subgraph')
                    for node_groups_combination in itertools.combinations(node_groups, r=n_groups):
                        node_groups_isites = []
                        node_groups_isites_combinations = []
                        for ng in node_groups_combination:
                            node_groups_isites_combinations.append([])
                            for nn in ng:
                                node_groups_isites_combinations[-1].append(nn)
                                node_groups_isites.append(nn)
                        node_groups_isites = sorted(node_groups_isites)
                        if np.all(node_groups_isites == isites_nodes):
                            print(' => YES ! There might be something here :-) should check that this subgraph ' \
                                  'is not "collinear" to the periodicity')
                            print(node_groups_isites_combinations)
                            input()
                        else:
                            print(' => NO !')
        input()

    def _node_match(self):
        return None

    def _edge_match(self):
        return None

    def _get_node_groups(self):
        self.project_1d()
        _divisors = divisors(len(self._connected_subgraph))
        _nodes_and_data = self._connected_subgraph.nodes(data=True)
        isites_nodes = sorted([nn for nn, dd in _nodes_and_data])
        all_node_groups = {}
        matching = {}

        for divisor in _divisors[:-1]:
            all_node_groups[divisor] = {}
            for node_bunch in itertools.combinations(_nodes_and_data, divisor):
                subgraph = self._connected_subgraph.subgraph([node for node, data in node_bunch])
                if not nx.is_connected(subgraph):  #This assumption should maybe be checked but I think its ok
                    continue
                gm = isomorphism.MultiGraphMatcher(self._connected_subgraph, subgraph, node_match=self._node_match(),
                                                   edge_match=self._edge_match())
                matching_groups = {}
                for subgraph_isomorphism in gm.subgraph_isomorphisms_iter():
                    subgraph_isites = tuple(sorted([nn for nn in list(subgraph_isomorphism.keys())]))
                    if subgraph_isites in matching_groups:
                        matching_groups[subgraph_isites].append(subgraph_isomorphism)
                    else:
                        matching_groups[subgraph_isites] = [subgraph_isomorphism]
                node_groups = list(matching_groups.keys())
                n_groups = len(self._connected_subgraph) / divisor
                if len(node_groups) >= n_groups:
                    #Which combination of [n_groups] groups of [divisor] nodes fills the full graph
                    for node_groups_combination in itertools.combinations(node_groups, r=n_groups):
                        ord_node_group = get_ordered_node_group(node_groups_combination)
                        this_node_group_isites = sorted([isite for gg in ord_node_group for isite in gg])
                        if np.all(this_node_group_isites == isites_nodes):
                            if not ord_node_group in all_node_groups[divisor]:
                                all_node_groups[divisor][ord_node_group] = matching_groups
        return all_node_groups

    def _get_shortest_cycles(self):
        #TODO: should add a check on cycle duplicates (i.e. 1-(0)-6-(1)-1 is equivalent to 1-(1)-6-(0)-1 )
        shortest_cycles = []
        for (node, data) in self._connected_subgraph.nodes(data=True):
            for path in get_all_simple_paths_edges(self._connected_subgraph, node, node):
                #We only want the shortest paths
                if len(shortest_cycles) > 0 and len(path) > len(shortest_cycles[0]):
                    continue
                this_path_delta = np.zeros(3, np.int)
                for (node1, node2, key) in path:
                    this_path_delta += get_delta(node1, node2, self._connected_subgraph[node1][node2][key])
                if np.any(this_path_delta != 0):
                    if len(shortest_cycles) > 0 and len(path) < len(shortest_cycles[0]):
                        shortest_cycles = []
                    shortest_cycles.append(path)
        return shortest_cycles

    def _get_periodic_cycles(self):
        #TODO: should add a check on cycle duplicates (i.e. 1-(0)-6-(1)-1 is equivalent to 1-(1)-6-(0)-1 )
        periodic_cycles = {}
        for (node, data) in self._connected_subgraph.nodes(data=True):
            periodic_cycles[node] = []
            for path in get_all_simple_paths_edges(self._connected_subgraph, node, node):
                #We only want the shortest paths
                if len(periodic_cycles[node]) > 0 and len(path) > len(periodic_cycles[node][0]):
                    continue
                this_path_delta = np.zeros(3, np.int)
                for (node1, node2, key) in path:
                    this_path_delta += get_delta(node1, node2, self._connected_subgraph[node1][node2][key])
                if np.any(this_path_delta != 0):
                    if len(periodic_cycles[node]) > 0 and len(path) < len(periodic_cycles[node][0]):
                        periodic_cycles[node] = []
                    periodic_cycles[node].append(path)
        return periodic_cycles

    def primitive_reduce_1d_node_connectivity(self):
        self.project_1d()
        if self._connected_subgraph.number_of_selfloops() > 0:
            self._primitive_reduced_connected_subgraph = self._connected_subgraph
            return self._primitive_reduced_connected_subgraph
        shortest_cycles = self._get_shortest_cycles()
        periodic_cycles = self._get_periodic_cycles()
        for node, value in periodic_cycles.items():
            print('Node ', node)
            for cycle in value:
                print(cycle)
        input('this was the cycles for each node')

        if self._node_match() is None:
            def _check_node(nn1, nn2):
                return True
        else:
            def _check_node(nn1, nn2):
                nmatch = self._node_match()
                return nmatch(nn1, nn2)
        if self._edge_match() is None:
            def _check_edge(e1, e2):
                return True
        else:
            def _check_edge(e1, e2):
                ematch = self._edge_match()
                return ematch(e1, e2)

        all_node_groups = self._get_node_groups()
        found_reduced_candidate = False
        for icycle, shortest_cycle in enumerate(shortest_cycles):
            print('cycle #{:d}'.format(icycle))
            found_for_cycle = False
            cycle_nodes = [node1 for node1, node2, key in shortest_cycle]
            for divisor, node_groups_combinations in all_node_groups.items():
                found_for_divisor = False
                ngroups = len(self._connected_subgraph) / divisor
                n_cycle_nodes_per_group = len(shortest_cycle) / ngroups
                if ngroups > len(shortest_cycle):
                    continue
                for node_groups_combination, node_groups_matching in node_groups_combinations.items():
                    group_combination_coherent_with_cycle = True
                    indices_ordered_group = []
                    for group in node_groups_combination:
                        if len(set(cycle_nodes).intersection(group)) != n_cycle_nodes_per_group:
                            group_combination_coherent_with_cycle = False
                            break
                        for isite in group:
                            if isite in cycle_nodes:
                                indices_ordered_group.append(cycle_nodes.index(isite) / n_cycle_nodes_per_group)
                                break
                    if not group_combination_coherent_with_cycle:
                        continue
                    #Find the connections between the last group and the first group => that gives the reduced candidate
                    #Construct the n[G-] full graph from the reduced candidate
                    #Check for isomorphism
                    #Remove groups connecting edges that do not have the same delta as in the cycle

                    #We first make the subgraph from the first group, this gives the nodes and the "internal" edges
                    reduced_graph_0 = self._connected_subgraph.subgraph(node_groups_combination[0])
                    edges_rg0 = []
                    for node_rg0 in reduced_graph_0.nodes():
                        for edge in reduced_graph_0.edges(nbunch=node_rg0, keys=True):
                            edges_rg0.append(edge)
                    print('internal eges : ')
                    for edge in edges_rg0:
                        print(edge)
                    connecting_edges = []
                    connecting_edges_starts = []
                    connecting_edges_ends = []
                    neighbor_nodes = []
                    neighbor_nodes_starts = []
                    neighbor_nodes_ends = []
                    for node_rg0 in reduced_graph_0.nodes():
                        for n1, n2, key, data in self._connected_subgraph.edges(nbunch=node_rg0, keys=True, data=True):
                            if (n1, n2, key) not in edges_rg0:
                                connecting_edges.append((n1, n2, key, data))
                                neighbor_nodes.append(node_rg0)
                                if data['start'] == node_rg0:
                                    neighbor_nodes_starts.append(node_rg0)
                                    connecting_edges_starts.append((n1, n2, key, data))
                                elif data['end'] == node_rg0:
                                    neighbor_nodes_ends.append(node_rg0)
                                    connecting_edges_ends.append((n1, n2, key, data))
                    found_candidate = False
                    #Breadth-first search starting at some node,
                    print('neighbors starts : ', neighbor_nodes_starts)
                    for edge in connecting_edges_starts:
                        print(edge)
                    print('neighbors ends : ', neighbor_nodes_ends)
                    for edge in connecting_edges_ends:
                        print(edge)
                    input('is that normal ?')
                    #To add the correct connecting edges, test if the some edge is in a shortest path
                    if np.mod(len(connecting_edges), 2) != 0:
                        input('Number of edges is odd : {:d}'.format(len(connecting_edges)))
                        continue
                    for indices_pairs_combination in all_pairs_combinations(connecting_edges, return_indices=True):
                        connecting_edges_candidate = []
                        found_pair_combination = True
                        for pair in indices_pairs_combination:
                            node_edge_matching = True
                            if not _check_node(self._connected_subgraph.node[connecting_edges[pair[0]][0]],
                                               self._connected_subgraph.node[connecting_edges[pair[1]][0]]):
                                node_edge_matching = False
                                break
                            if not _check_edge(self._connected_subgraph
                                               [connecting_edges[pair[0]][0]]
                                               [connecting_edges[pair[0]][1]]
                                               [connecting_edges[pair[0]][2]],
                                               self._connected_subgraph
                                               [connecting_edges[pair[1]][0]]
                                               [connecting_edges[pair[1]][1]]
                                               [connecting_edges[pair[1]][2]]):
                                node_edge_matching = False
                                break
                            found_pair_cycle = False
                            for periodic_cycle in periodic_cycles[connecting_edges[pair[0]][0]]:
                                delta_candidate_1 = cycle_contains_edge(periodic_cycle, connecting_edges[pair[0]])
                                delta_candidate_2 = cycle_contains_edge(periodic_cycle, connecting_edges[pair[1]])
                                if delta_candidate_1 and delta_candidate_2:
                                    found_pair_cycle = True
                                    mydata = dict(connecting_edges[pair[0]][3])
                                    mydata['start'] = connecting_edges[pair[0]][0]
                                    mydata['end'] = connecting_edges[pair[1]][0]
                                    mydata['delta'] = delta_candidate_1
                                    connecting_edges_candidate.append((mydata['start'], mydata['end'], mydata))
                                    break
                            if not found_pair_cycle:
                                break
                        #check nodes and edges
                        #
                    #for edge in connecting_edges_starts:
                    #    node_start = edge[0]
                    #    for periodic_cycle in periodic_cycles[node_start]:
                    #        delta_candidate = cycle_contains_edge(periodic_cycle, edge)
                    #        if delta_candidate:
                    #            print 'duh'

                    # for end_node_indices_perm in itertools.permutations(range(len(neighbor_nodes_ends))):
                    #     node_edge_matching = True
                    #     connecting_edges_candidate = []
                    #     print 'TESTING PERMUTATION'
                    #     for i_node_couple in range(len(neighbor_nodes_ends)):
                    #         if not _check_node(neighbor_nodes_starts[i_node_couple],
                    #                            neighbor_nodes_ends[end_node_indices_perm[i_node_couple]]):
                    #             node_edge_matching = False
                    #             break
                    #         if not _check_edge(connecting_edges_starts[i_node_couple],
                    #                            connecting_edges_ends[end_node_indices_perm[i_node_couple]]):
                    #             node_edge_matching = False
                    #             break
                    #         if connecting_edges_starts[i_node_couple][3]['delta'] == (0, 0, 0):
                    #             if connecting_edges_ends[end_node_indices_perm[i_node_couple]][3]['delta'] == (0, 0, 0):
                    #                 raw_input('duuuhhh should not happen !')
                    #                 break
                    #             mydata = dict(connecting_edges_ends[end_node_indices_perm[i_node_couple]][3])
                    #             mydata['start'] = neighbor_nodes_starts[i_node_couple]
                    #             mydata['end'] = neighbor_nodes_ends[end_node_indices_perm[i_node_couple]]
                    #             print 'connecting : '
                    #             print connecting_edges_starts[i_node_couple]
                    #             print connecting_edges_ends[end_node_indices_perm[i_node_couple]]
                    #             print mydata['start'], mydata['end']
                    #             connecting_edges_candidate.append((mydata['start'], mydata['end'],
                    #                                                mydata))
                    #         elif connecting_edges_ends[end_node_indices_perm[i_node_couple]][3]['delta'] == (0, 0, 0):
                    #             mydata = dict(connecting_edges_starts[i_node_couple][3])
                    #             mydata['start'] = neighbor_nodes_starts[i_node_couple]
                    #             mydata['end'] = neighbor_nodes_ends[end_node_indices_perm[i_node_couple]]
                    #             connecting_edges_candidate.append((mydata['start'], mydata['end'],
                    #                                                mydata))
                    #             print 'connecting : '
                    #             print connecting_edges_starts[i_node_couple]
                    #             print connecting_edges_ends[end_node_indices_perm[i_node_couple]]
                    #             print mydata['start'], mydata['end']
                    #         else:
                    #             node_edge_matching = False
                    #             print connecting_edges_starts[i_node_couple]
                    #             print connecting_edges_ends[end_node_indices_perm[i_node_couple]]
                    #             raw_input('is that normal that None of the deltas are 0 ?????')
                    #             break
                    #     if not node_edge_matching:
                    #         continue
                    #
                    #     reduced_graph_candidate = self._connected_subgraph.subgraph(node_groups_combination[0])
                    #     reduced_graph_candidate.add_edges_from(connecting_edges_candidate)
                    #     supergraph = make_supergraph(reduced_graph_candidate, ngroups, self._periodicity)
                    #     if nx.is_isomorphic(supergraph, self._connected_subgraph, self._node_match(), self._edge_match()):
                    #         found_candidate = True
                    #         break
                    # if found_candidate:
                    #     found_for_divisor = True
                    #     break
                if found_for_divisor:
                    found_for_cycle = True
                    break
            if found_for_cycle:
                found_reduced_candidate = True
                break
        if found_reduced_candidate:
            print('WE FOUND A REDUCED !')
            self._reduced_connected_subgraph = reduced_graph_candidate

        else:
            print('no reduced found ...')
            self._reduced_connected_subgraph = self._connected_subgraph.copy()
        input('lets see it')

    def make_supergraph(self, multiplicity):
        supergraph = make_supergraph(self._connected_subgraph, multiplicity, self._periodicity_vectors)
        return supergraph

    def project_1d(self):
        """
        Projects all the delta's to the periodic vector of the 1d connected component
        :rtype : object
        """
        if not self._projected:
            if self._periodicity_vectors is None:
                self.compute_periodicity()
            for n1, n2, key, data in self._connected_subgraph.edges(data=True, keys=True):
                self._connected_subgraph[n1][n2][key]['unprojected_delta'] = data['delta']
                self._connected_subgraph[n1][n2][key]['delta'] = tuple(self._periodicity_vectors[0] * data['delta'])
                print(self._connected_subgraph[n1][n2][key]['ligands'])
                for iligand, (isite_ligand, e1_ligand, e2_ligand) in enumerate(self._connected_subgraph[n1][n2][key][
                    'ligands']):
                    new_e1_ligand = dict(e1_ligand)
                    new_e1_ligand['unprojected_delta'] = e1_ligand['delta']
                    new_e1_ligand['delta'] = tuple(self._periodicity_vectors[0] * e1_ligand['delta'])
                    new_e2_ligand = dict(e2_ligand)
                    new_e2_ligand['unprojected_delta'] = e2_ligand['delta']
                    new_e2_ligand['delta'] = tuple(self._periodicity_vectors[0] * e2_ligand['delta'])
                    self._connected_subgraph[n1][n2][key]['ligands'][iligand] = (isite_ligand, new_e1_ligand,
                                                                                 new_e2_ligand)
            self._projected = True
        return

    def unproject(self):
        """
        Recover all the initial (unprojected) delta's
        :return:
        """
        if self._projected:
            for n1, n2, key, data in self._connected_subgraph.edges(data=True, keys=True):
                self._connected_subgraph[n1][n2][key]['delta'] = self._connected_subgraph[n1][n2][key][
                    'unprojected_delta']
                self._connected_subgraph[n1][n2][key].pop('unprojected_delta')
                for iligand, (isite_ligand, e1_ligand, e2_ligand) in enumerate(
                        self._connected_subgraph[n1][n2][key]['ligands']):
                    old_e1_ligand = dict(e1_ligand)
                    old_e1_ligand['delta'] = e1_ligand['unprojected_delta']
                    old_e1_ligand.pop('unprojected_delta')
                    old_e2_ligand = dict(e2_ligand)
                    old_e2_ligand['delta'] = e2_ligand['unprojected_delta']
                    old_e2_ligand.pop('unprojected_delta')
                    self._connected_subgraph[n1][n2][key]['ligands'][iligand] = (isite_ligand, old_e1_ligand,
                                                                                 old_e2_ligand)
            self._projected = False
        return

    def show_graph(self, graph=None, save_file=None):
        import matplotlib.pyplot as plt

        if graph is None:
            shown_graph = self._connected_subgraph
        else:
            shown_graph = graph

        #pos = nx.spring_layout(shown_graph)
        pos = nx.shell_layout(shown_graph)
        ax = plt.gca()
        draw_network(shown_graph, pos, ax, periodicity_vectors=self._periodicity_vectors)
        ax.autoscale()
        plt.axis('equal')
        plt.axis('off')
        if save_file is not None:
            plt.savefig(save_file)
        #nx.draw(self._connected_subgraph)
        plt.show()

    @property
    def graph(self):
        return self._connected_subgraph

    @property
    def primitive_reduced_graph(self):
        if self._primitive_reduced_connected_subgraph is None:
            self.primitive_reduce()
        return self._primitive_reduced_connected_subgraph

    @property
    def is_periodic(self):
        return not self.is_0d

    @property
    def is_0d(self):
        if self._periodicity_vectors is None:
            self.compute_periodicity()
        return len(self._periodicity_vectors) == 0

    @property
    def is_1d(self):
        if self._periodicity_vectors is None:
            self.compute_periodicity()
        return len(self._periodicity_vectors) == 1

    @property
    def is_2d(self):
        if self._periodicity_vectors is None:
            self.compute_periodicity()
        return len(self._periodicity_vectors) == 2

    @property
    def is_3d(self):
        if self._periodicity_vectors is None:
            self.compute_periodicity()
        return len(self._periodicity_vectors) == 3

    @property
    def periodicity_vectors(self):
        if self._periodicity_vectors is None:
            self.compute_periodicity()
        return [np.array(pp) for pp in self._periodicity_vectors]

    @property
    def periodicity(self):
        if self._periodicity_vectors is None:
            self.compute_periodicity()
        return '{:d}D'.format(len(self._periodicity_vectors))

    def as_dict(self):
        """
        Bson-serializable dict representation of the ConnectedComponent object.
        :return: Bson-serializable dict representation of the ConnectedComponent object.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__, }

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs the ConnectedComponent object from a dict representation of the
        ConnectedComponent object created using the as_dict method.
        :param d: dict representation of the ConnectedComponent object
        :return: ConnectedComponent object
        """
        return cls()

    @classmethod
    def from_graph(cls, g):
        """
        Constructor for the ConnectedComponent object from a graph of the connected component
        :param g: Graph of the connected component
        :return: ConnectedComponent object
        """
        return cls(graph=g)
