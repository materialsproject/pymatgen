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
    e = None
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

    def make_supergraph(self, multiplicity):
        supergraph = make_supergraph(self._connected_subgraph, multiplicity, self._periodicity_vectors)
        return supergraph

    def show_graph(self, graph=None, save_file=None, drawing_type='internal'):
        import matplotlib.pyplot as plt

        if graph is None:
            shown_graph = self._connected_subgraph
        else:
            shown_graph = graph

        #pos = nx.spring_layout(shown_graph)
        if drawing_type == 'internal':
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
        elif drawing_type == 'draw_graphviz':
            import networkx
            networkx.nx_pydot.graphviz_layout(shown_graph)
            plt.show()
        elif drawing_type == 'draw_random':
            import networkx
            networkx.draw_random(shown_graph)
            plt.show()

    @property
    def graph(self):
        return self._connected_subgraph

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
