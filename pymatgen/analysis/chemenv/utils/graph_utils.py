import networkx as nx
import numpy as np
from collections import deque
from monty.json import MSONable

__author__ = 'waroquiers'


def get_delta(node1, node2, edge_data):
    if node1.isite == edge_data['start'] and node2.isite == edge_data['end']:
        return np.array(edge_data['delta'], dtype=np.int)
    elif node2.isite == edge_data['start'] and node1.isite == edge_data['end']:
        return -np.array(edge_data['delta'], dtype=np.int)
    else:
        raise ValueError('Trying to find a delta between two nodes with an edge '
                         'that seems not to link these nodes.')


def get_all_simple_paths_edges(graph, source, target, cutoff=None, data=True):
    edge_paths = []
    if not graph.is_multigraph():
        for path in nx.all_simple_paths(graph, source, target, cutoff=cutoff):
            edge_paths.append([(path[ii], path[ii+1]) for ii in range(len(path) - 1)])
        return edge_paths

    node_paths = []
    for path in nx.all_simple_paths(graph, source, target, cutoff=cutoff):
        exists = False
        for path2 in node_paths:
            if len(path2) == len(path) and np.all(path == path2):
                exists = True
                break
        if exists:
            continue
        node_paths.append(path)
        current_edge_paths = [[]]
        for (node1, node2) in [(node1, path[inode1 + 1]) for inode1, node1 in enumerate(path[:-1])]:
            new_edge_paths = []
            for key, edge_data in graph[node1][node2].items():
                for tmp_edge_path in current_edge_paths:
                    if data:
                        new_path = [(n1, n2, k, d) for (n1, n2, k, d) in tmp_edge_path]
                        new_path.append((node1, node2, key, edge_data))
                        new_edge_paths.append(new_path)
                    else:
                        new_path = [(n1, n2, k) for (n1, n2, k) in tmp_edge_path]
                        new_path.append((node1, node2, key))
                        new_edge_paths.append(new_path)
            current_edge_paths = new_edge_paths
        edge_paths.extend(current_edge_paths)
    return edge_paths


def _c2index_isreverse(c1, c2):
    """
    Private helper function to get the index c2_0_index of the first node of cycle c1
    in cycle c2 and whether the cycle c2 should be reversed or not.

    Returns None if the first node of cycle c1 is not found in cycle c2.
    The reverse value depends on the index c2_1_index of the second node of cycle c1 in
    cycle c2 : if it is *just after* the c2_0_index, reverse is False, if it is
    *just before* the c2_0_index, reverse is True, otherwise the function returns None).
    """

    c1_0 = c1.nodes[0]
    c1_1 = c1.nodes[1]
    if c1_0 not in c2.nodes:
        return None, None, 'First node of cycle c1 not found in cycle c2.'
    if c1_1 not in c2.nodes:
        return None, None, 'Second node of cycle c1 not found in cycle c2.'
    c2_0_index = c2.nodes.index(c1_0)
    c2_1_index = c2.nodes.index(c1_1)
    if c2_0_index == 0:
        if c2_1_index == 1:
            reverse = False
        elif c2_1_index == len(c2.nodes) - 1:
            reverse = True
        else:
            return None, None, 'Second node of cycle c1 is not second or last in cycle c2 ' \
                               '(first node of cycle c1 is first in cycle c2).'
    elif c2_0_index == len(c2.nodes) - 1:
        if c2_1_index == 0:
            reverse = False
        elif c2_1_index == c2_0_index - 1:
            reverse = True
        else:
            return None, None, 'Second node of cycle c1 is not first or before last in cycle c2 ' \
                               '(first node of cycle c1 is last in cycle c2).'
    else:
        if c2_1_index == c2_0_index + 1:
            reverse = False
        elif c2_1_index == c2_0_index - 1:
            reverse = True
        else:
            return None, None, 'Second node of cycle c1 in cycle c2 is not just after or ' \
                               'just before first node of cycle c1 in cycle c2.'
    return c2_0_index, reverse, ''


class SimpleGraphCycle(MSONable):
    """
    Class used to describe a cycle in a simple graph (graph without multiple edges).

    Note that the convention used here is the networkx convention for which simple graphs allow
    to have self-loops in a simple graph.
    """

    def __init__(self, nodes):
        self.nodes = deque(nodes)
        is_valid, msg = self._is_valid()
        if not is_valid:
            raise ValueError('SimpleGraphCycle is not valid : {}'.format(msg))

    def _is_valid(self):
        """Check if a cycle is valid.

        This method checks :
        - that there are no duplicate nodes,
        :return: True if the cycle is valid, False otherwise.
        """
        if len(self.nodes) != len(set(self.nodes)):  # Should not have duplicate nodes
            return False, 'Duplicate nodes.'
        return True, ''

    def __eq__(self, other):
        if len(self.nodes) != len(other.nodes):  # Trivial non equivalence
            return False
        if len(self.nodes) == 1:  # Self-loop case
            return self.nodes[0] == other.nodes[0]
        if len(self.nodes) == 2:  # Special case with two nodes
            return ((self.nodes[0] == other.nodes[0] and self.nodes[1] == other.nodes[1]) or
                    (self.nodes[0] == other.nodes[1] and self.nodes[1] == other.nodes[0]))
        other_0_index, reverse, errmsg = _c2index_isreverse(self, other)
        if other_0_index is None:
            return False
        other_nodes_copy = deque(other.nodes)
        if reverse:
            other_nodes_copy.reverse()
            n_rotation = len(other_nodes_copy) - other_0_index - 1
            other_nodes_copy.rotate(-n_rotation)
        else:
            other_nodes_copy.rotate(-other_0_index)
        return self.nodes == other_nodes_copy

    def __hash__(self):
        return len(self.nodes)

    def __len__(self):
        return len(self.nodes)

    @classmethod
    def from_edges(cls, edges):
        remaining_edges = list(edges)
        nodes = list(remaining_edges.pop())
        while remaining_edges:
            prev_node = nodes[-1]
            for ie, e in enumerate(remaining_edges):
                if prev_node == e[0]:
                    remaining_edges.pop(ie)
                    nodes.append(e[1])
                    break
                elif prev_node == e[1]:
                    remaining_edges.pop(ie)
                    nodes.append(e[0])
                    break
            else:  # did not find the next edge
                raise ValueError('Could not construct a cycle from edges.')
        if nodes[0] != nodes[-1]:
            raise ValueError('First node not equal to last node from edges.')
        nodes.pop()
        return cls(nodes)

    def as_dict(self):
        d = MSONable.as_dict(self)
        # Transforming deque object to a list to allow BSON and MongoDB
        d['nodes'] = list(d['nodes'])
        return d


class MultiGraphCycle(MSONable):
    """
    Class used to describe a cycle in a multigraph.

    nodes are the nodes of the cycle and edge_indices are the indices of the edges in the cycle.
    The nth index in edge_indices corresponds to the edge index between the nth node in nodes and
    the (n+1)th node in nodes with the exception of the last one being the edge index between
    the last node in nodes and the first node in nodes
    Examples :
      - A graph with one node (0) and two self-loops gives two cycles :
        - c1 = MultiGraphCycle([0], [0])
        - c2 = MultiGraphCycle([0], [1])
      - A graph with a simple cycle 0-1-2-3 with two additional edges (0-1) and (2-3)
        gives 6 cycles ("rotated" and "reversed" cycles are supposedly removed).
        - c1 = MultiGraphCycle([0, 1, 2, 3], [0, 0, 0, 0])
        - c2 = MultiGraphCycle([0, 1, 2, 3], [1, 0, 0, 0])
        - c3 = MultiGraphCycle([0, 1, 2, 3], [0, 0, 1, 0])
        - c4 = MultiGraphCycle([0, 1, 2, 3], [1, 0, 1, 0])
        - c5 = MultiGraphCycle([0, 1], [0, 1])
        - c6 = MultiGraphCycle([2, 3], [0, 1])

    Example: A cycle
        nodes:          1 - 3 - 4 - 0 - 2 - (1)
        edge_indices:     0 . 1 . 0 . 2 . 0 . (0)
    """
    def __init__(self, nodes, edge_indices):
        """

        :param nodes:
        :param edge_indices:
        """
        self.nodes = deque(nodes)
        self.edge_indices = deque(edge_indices)
        is_valid, msg = self._is_valid()
        if not is_valid:
            raise ValueError('MultiGraphCycle is not valid : {}'.format(msg))
        self.edge_deltas = None
        self.per = None

    def _is_valid(self):
        """Check if a cycle is valid.

        This method checks :
        - that the number of nodes and number of edges are equal,
        - that there are no duplicate nodes,
        - that edges in cycles with two nodes are different.
        :param cycle: Cycle object
        :return: True if the cycle is valid, False otherwise.
        """
        if len(self.nodes) != len(self.edge_indices):  # Should have the same number of nodes and edges
            return False, 'Number of nodes different from number of edge indices.'
        if len(self.nodes) != len(set(self.nodes)):  # Should not have duplicate nodes
            return False, 'Duplicate nodes.'
        if len(self.nodes) == 2:  # Cycles with two nodes cannot use the same edge for the cycle
            if self.edge_indices[0] == self.edge_indices[1]:
                return False, 'Cycles with two nodes cannot use the same edge for the cycle.'
        return True, ''

    def __eq__(self, other):
        if len(self.nodes) != len(other.nodes):  # Trivial non equivalence
            return False
        if len(self.nodes) == 1:  # Self-loop case
            return self.nodes[0] == other.nodes[0] and self.edge_indices[0] == other.edge_indices[0]
        if len(self.nodes) == 2:  # Special case with two nodes
            return (((self.nodes[0] == other.nodes[0] and self.nodes[1] == other.nodes[1]) or
                     (self.nodes[0] == other.nodes[1] and self.nodes[1] == other.nodes[0])) and
                    ((self.edge_indices[0] == other.edge_indices[0] and self.edge_indices[1] == other.edge_indices[1]) or
                     (self.edge_indices[0] == other.edge_indices[1] and self.edge_indices[1] == other.edge_indices[0])))
        other_0_index, reverse, errmsg = _c2index_isreverse(self, other)
        if other_0_index is None:
            return False
        other_nodes_copy = deque(other.nodes)
        other_edge_indices_copy = deque(other.edge_indices)
        if reverse:
            other_nodes_copy.reverse()
            other_edge_indices_copy.reverse()
            other_edge_indices_copy.rotate(-1)
            n_rotation = len(other_nodes_copy) - other_0_index - 1
            other_nodes_copy.rotate(-n_rotation)
            other_edge_indices_copy.rotate(-n_rotation)
        else:
            other_nodes_copy.rotate(-other_0_index)
            other_edge_indices_copy.rotate(-other_0_index)
        return self.nodes == other_nodes_copy and self.edge_indices == other_edge_indices_copy

    def __hash__(self):
        return len(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def periodicity(self, multigraph):
        if self.per is None:
            self.edge_deltas = []
            self.per = np.zeros(3, dtype=np.int)
            for inode, node in enumerate(self.nodes):
                n1 = node
                n2 = self.nodes[(inode+1) % len(self.nodes)]
                ie = self.edge_indices[inode]
                edge_data = multigraph[n1][n2][ie]
                delta = get_delta(node1=n1, node2=n2, edge_data=edge_data)
                self.edge_deltas.append(delta)
                self.per += delta
        return self.per

    def full_description(self, multigraph=None):
        if self.per is None and multigraph is None:
            raise ValueError('Should provide multigraph to get full description')
        out = ['Cycle [{}-({:d})] with periodicity ({}) :'.format('-'.join(['{:d}'.format(n) for n in self.nodes]),
                                                                  self.nodes[0],
                                                                  ' '.join(['{:d}'.format(ii) for ii in self.per]))]
        for inode, node in enumerate(self.nodes):
            out.append('  * Node {:d}'.format(node))
            out.append('     ===({})===>'.format(' '.join(['{:d}'.format(ii) for ii in self.edge_deltas[inode]])))
        out.append('  * Node {:d}'.format(self.nodes[0]))
        return '\n'.join(out)

    def as_dict(self):
        d = MSONable.as_dict(self)
        # Transforming deque objects to lists to allow BSON and MongoDB
        d['nodes'] = list(d['nodes'])
        d['edge_indices'] = list(d['edge_indices'])
        return d
