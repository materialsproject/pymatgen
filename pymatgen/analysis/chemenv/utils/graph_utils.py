"""
This module contains some graph utils that are used in the chemenv package.
"""

import itertools
import operator

import networkx as nx
import numpy as np
from monty.json import MSONable

__author__ = "waroquiers"


def get_delta(node1, node2, edge_data):
    """
    Get the delta.
    :param node1:
    :param node2:
    :param edge_data:
    :return:
    """
    if node1.isite == edge_data["start"] and node2.isite == edge_data["end"]:
        return np.array(edge_data["delta"], dtype=np.int_)
    if node2.isite == edge_data["start"] and node1.isite == edge_data["end"]:
        return -np.array(edge_data["delta"], dtype=np.int_)
    raise ValueError("Trying to find a delta between two nodes with an edge that seems not to link these nodes.")


def get_all_simple_paths_edges(graph, source, target, cutoff=None, data=True):
    """
    Get all the simple path and edges.
    :param graph:
    :param source:
    :param target:
    :param cutoff:
    :param data:
    :return:
    """
    edge_paths = []
    if not graph.is_multigraph():
        for path in nx.all_simple_paths(graph, source, target, cutoff=cutoff):
            edge_paths.append([(path[ii], path[ii + 1]) for ii in range(len(path) - 1)])
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
                        new_path = list(tmp_edge_path)
                        new_path.append((node1, node2, key, edge_data))
                        new_edge_paths.append(new_path)
                    else:
                        new_path = list(tmp_edge_path)
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
        return None, None, "First node of cycle c1 not found in cycle c2."
    if c1_1 not in c2.nodes:
        return None, None, "Second node of cycle c1 not found in cycle c2."
    c2_0_index = c2.nodes.index(c1_0)
    c2_1_index = c2.nodes.index(c1_1)
    if c2_0_index == 0:
        if c2_1_index == 1:
            reverse = False
        elif c2_1_index == len(c2.nodes) - 1:
            reverse = True
        else:
            return (
                None,
                None,
                "Second node of cycle c1 is not second or last in cycle c2 "
                "(first node of cycle c1 is first in cycle c2).",
            )
    elif c2_0_index == len(c2.nodes) - 1:
        if c2_1_index == 0:
            reverse = False
        elif c2_1_index == c2_0_index - 1:
            reverse = True
        else:
            return (
                None,
                None,
                "Second node of cycle c1 is not first or before last in cycle c2 "
                "(first node of cycle c1 is last in cycle c2).",
            )
    else:
        if c2_1_index == c2_0_index + 1:
            reverse = False
        elif c2_1_index == c2_0_index - 1:
            reverse = True
        else:
            return (
                None,
                None,
                "Second node of cycle c1 in cycle c2 is not just after or "
                "just before first node of cycle c1 in cycle c2.",
            )
    return c2_0_index, reverse, ""


class SimpleGraphCycle(MSONable):
    """
    Class used to describe a cycle in a simple graph (graph without multiple edges).

    Note that the convention used here is the networkx convention for which simple graphs allow
    to have self-loops in a simple graph.
    No simple graph cycle with two nodes is possible in a simple graph. The graph will not
    be validated if validate is set to False.
    By default, the "ordered" parameter is None, in which case the SimpleGraphCycle will be ordered.
    If the user explicitly sets ordered to False, the SimpleGraphCycle will not be ordered.
    """

    def __init__(self, nodes, validate=True, ordered=None):
        """
        :param nodes:
        :param validate:
        :param ordered:
        """
        self.nodes = tuple(nodes)
        if validate:
            self.validate()
        if ordered is not None:
            self.ordered = ordered
        else:
            self.order()

    def _is_valid(self, check_strict_ordering=False):
        """Check if a SimpleGraphCycle is valid.

        This method checks :
        - that there are no duplicate nodes,
        - that there are either 1 or more than 2 nodes
        :return: True if the SimpleGraphCycle is valid, False otherwise.
        """
        if len(self.nodes) == 1:
            return True, ""
        if len(self.nodes) == 2:
            return False, "Simple graph cycle with 2 nodes is not valid."
        if len(self.nodes) == 0:
            return False, "Empty cycle is not valid."
        if len(self.nodes) != len(set(self.nodes)):  # Should not have duplicate nodes
            return False, "Duplicate nodes."
        if check_strict_ordering:
            try:
                sorted_nodes = sorted(self.nodes)
            except TypeError as te:
                msg = te.args[0]
                if "'<' not supported between instances of" in msg:
                    return False, "The nodes are not sortable."
                raise
            res = all(i < j for i, j in zip(sorted_nodes, sorted_nodes[1:]))
            if not res:
                return (
                    False,
                    "The list of nodes in the cycle cannot be strictly ordered.",
                )
        return True, ""

    def validate(self, check_strict_ordering=False):
        """
        :param check_strict_ordering:
        :return:
        """
        is_valid, msg = self._is_valid(check_strict_ordering=check_strict_ordering)
        if not is_valid:
            raise ValueError("SimpleGraphCycle is not valid : {}".format(msg))

    def order(self, raise_on_fail=True):
        """Orders the SimpleGraphCycle.

        The ordering is performed such that the first node is the "lowest" one
        and the second node is the lowest one of the two neighbor nodes of the
        first node. If raise_on_fail is set to True a RuntimeError will be
        raised if the ordering fails.

        :param raise_on_fail: If set to True, will raise a RuntimeError if the
                              ordering fails.
        :return: None
        """
        # always validate the cycle if it needs to be ordered
        # also validates that the nodes can be strictly ordered
        try:
            self.validate(check_strict_ordering=True)
        except ValueError as ve:
            msg = ve.args[0]
            if "SimpleGraphCycle is not valid :" in msg and not raise_on_fail:
                self.ordered = False
                return
            raise

        if len(self.nodes) == 1:
            self.ordered = True
            return

        # Not sure whether the following should be checked here... if strict ordering was guaranteed by
        # the validate method, why would it be needed to have a unique class. One could have 2 subclasses of
        # the same parent class and things could be ok. To be checked what to do. (see also MultiGraphCycle)
        node_classes = set(n.__class__ for n in self.nodes)
        if len(node_classes) > 1:
            if raise_on_fail:
                raise ValueError("Could not order simple graph cycle as the nodes are of different classes.")
            self.ordered = False
            return

        min_index, min_node = min(enumerate(self.nodes), key=operator.itemgetter(1))
        reverse = self.nodes[(min_index - 1) % len(self.nodes)] < self.nodes[(min_index + 1) % len(self.nodes)]
        if reverse:
            self.nodes = self.nodes[min_index::-1] + self.nodes[:min_index:-1]
        else:
            self.nodes = self.nodes[min_index:] + self.nodes[:min_index]
        self.ordered = True

    def __hash__(self):
        return len(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def __str__(self):
        out = ["Simple cycle with nodes :"]
        out.extend([str(node) for node in self.nodes])
        return "\n".join(out)

    def __eq__(self, other):
        if not self.ordered or not other.ordered:
            raise RuntimeError("Simple cycles should be ordered in order to be compared.")
        return self.nodes == other.nodes

    @classmethod
    def from_edges(cls, edges, edges_are_ordered=True):
        """Constructs SimpleGraphCycle from a list edges.

        By default, the edges list is supposed to be ordered as it will be
        much faster to construct the cycle. If edges_are_ordered is set to
        False, the code will automatically try to find the corresponding edge
        order in the list.
        """
        if edges_are_ordered:
            nodes = [e[0] for e in edges]
            if not all(e1e2[0][1] == e1e2[1][0] for e1e2 in zip(edges, edges[1:])) or edges[-1][1] != edges[0][0]:
                raise ValueError("Could not construct a cycle from edges.")
        else:
            remaining_edges = list(edges)
            nodes = list(remaining_edges.pop())
            while remaining_edges:
                prev_node = nodes[-1]
                for ie, e in enumerate(remaining_edges):
                    if prev_node == e[0]:
                        remaining_edges.pop(ie)
                        nodes.append(e[1])
                        break
                    if prev_node == e[1]:
                        remaining_edges.pop(ie)
                        nodes.append(e[0])
                        break
                else:  # did not find the next edge
                    raise ValueError("Could not construct a cycle from edges.")
            if nodes[0] != nodes[-1]:
                raise ValueError("Could not construct a cycle from edges.")
            nodes.pop()
        return cls(nodes)

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        d = MSONable.as_dict(self)
        # Transforming tuple object to a list to allow BSON and MongoDB
        d["nodes"] = list(d["nodes"])
        return d

    @classmethod
    def from_dict(cls, d, validate=False):
        """
        Serialize from dict.
        :param d:
        :param validate:
        :return:
        """
        return cls(nodes=d["nodes"], validate=validate, ordered=d["ordered"])


class MultiGraphCycle(MSONable):
    """Class used to describe a cycle in a multigraph.

    nodes are the nodes of the cycle and edge_indices are the indices of the edges in the cycle.
    The nth index in edge_indices corresponds to the edge index between the nth node in nodes and
    the (n+1)th node in nodes with the exception of the last one being the edge index between
    the last node in nodes and the first node in nodes

    Example: A cycle
        nodes:          1 - 3 - 4 - 0 - 2 - (1)
        edge_indices:     0 . 1 . 0 . 2 . 0 . (0)
    """

    def __init__(self, nodes, edge_indices, validate=True, ordered=None):
        """
        :param nodes:
        :param edge_indices:
        :param validate:
        :param ordered:
        """
        self.nodes = tuple(nodes)
        self.edge_indices = tuple(edge_indices)
        if validate:
            self.validate()
        if ordered is not None:
            self.ordered = ordered
        else:
            self.order()
        self.edge_deltas = None
        self.per = None

    def _is_valid(self, check_strict_ordering=False):
        """Check if a MultiGraphCycle is valid.

        This method checks :
        - that there are no duplicate nodes,
        - that there are either 1 or more than 2 nodes
        :return: True if the SimpleGraphCycle is valid, False otherwise.
        """
        if len(self.nodes) != len(self.edge_indices):  # Should have the same number of nodes and edges
            return False, "Number of nodes different from number of edge indices."
        if len(self.nodes) == 0:
            return False, "Empty cycle is not valid."
        if len(self.nodes) != len(set(self.nodes)):  # Should not have duplicate nodes
            return False, "Duplicate nodes."
        if len(self.nodes) == 2:  # Cycles with two nodes cannot use the same edge for the cycle
            if self.edge_indices[0] == self.edge_indices[1]:
                return (
                    False,
                    "Cycles with two nodes cannot use the same edge for the cycle.",
                )
        if check_strict_ordering:
            try:
                sorted_nodes = sorted(self.nodes)
            except TypeError as te:
                msg = te.args[0]
                if "'<' not supported between instances of" in msg:
                    return False, "The nodes are not sortable."
                raise
            res = all(i < j for i, j in zip(sorted_nodes, sorted_nodes[1:]))
            if not res:
                return (
                    False,
                    "The list of nodes in the cycle cannot be strictly ordered.",
                )
        return True, ""

    def validate(self, check_strict_ordering=False):
        """
        :param check_strict_ordering:
        :return:
        """
        is_valid, msg = self._is_valid(check_strict_ordering=check_strict_ordering)
        if not is_valid:
            raise ValueError("MultiGraphCycle is not valid : {}".format(msg))

    def order(self, raise_on_fail=True):
        """Orders the SimpleGraphCycle.

        The ordering is performed such that the first node is the "lowest" one
        and the second node is the lowest one of the two neighbor nodes of the
        first node. If raise_on_fail is set to True a RuntimeError will be
        raised if the ordering fails.

        :param raise_on_fail: If set to True, will raise a RuntimeError if the
                              ordering fails.
        :return: None
        """
        # always validate the cycle if it needs to be ordered
        # also validates that the nodes can be strictly ordered
        try:
            self.validate(check_strict_ordering=True)
        except ValueError as ve:
            msg = ve.args[0]
            if "MultiGraphCycle is not valid :" in msg and not raise_on_fail:
                self.ordered = False
                return
            raise

        if len(self.nodes) == 1:
            self.ordered = True
            return

        # Not sure whether the following should be checked here... if strict ordering was guaranteed by
        # the validate method, why would it be needed to have a unique class. One could have 2 subclasses of
        # the same parent class and things could be ok. To be checked what to do. (see also SimpleGraphCycle)
        node_classes = set(n.__class__ for n in self.nodes)
        if len(node_classes) > 1:
            if raise_on_fail:
                raise ValueError("Could not order simple graph cycle as the nodes are of different classes.")
            self.ordered = False
            return

        min_index, min_node = min(enumerate(self.nodes), key=operator.itemgetter(1))

        # Special case when number of nodes is 2 because the two
        # edge_indices refer to the same pair of nodes
        if len(self.nodes) == 2:
            self.nodes = tuple(sorted(self.nodes))
            self.edge_indices = tuple(sorted(self.edge_indices))
            self.ordered = True
            return

        reverse = self.nodes[(min_index - 1) % len(self.nodes)] < self.nodes[(min_index + 1) % len(self.nodes)]
        if reverse:
            self.nodes = self.nodes[min_index::-1] + self.nodes[:min_index:-1]
            min_edge_index = (min_index - 1) % len(self.nodes)
            self.edge_indices = self.edge_indices[min_edge_index::-1] + self.edge_indices[:min_edge_index:-1]
        else:
            self.nodes = self.nodes[min_index:] + self.nodes[:min_index]
            self.edge_indices = self.edge_indices[min_index:] + self.edge_indices[:min_index]
        self.ordered = True

    def __hash__(self):
        return len(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def __str__(self):
        out = ["Multigraph cycle with nodes :"]
        cycle = []
        for inode, node1, node2 in zip(itertools.count(), self.nodes[:-1], self.nodes[1:]):
            cycle.append("{} -*{:d}*- {}".format(str(node1), self.edge_indices[inode], str(node2)))
        cycle.append("{} -*{:d}*- {}".format(str(self.nodes[-1]), self.edge_indices[-1], str(self.nodes[0])))
        # out.extend([str(node) for node in self.nodes])
        out.extend(cycle)
        return "\n".join(out)

    def __eq__(self, other):
        if not self.ordered or not other.ordered:
            raise RuntimeError("Multigraph cycles should be ordered in order to be compared.")
        return self.nodes == other.nodes and self.edge_indices == other.edge_indices


def get_all_elementary_cycles(graph):
    """
    :param graph:
    :return:
    """
    if not isinstance(graph, nx.Graph):
        raise TypeError("graph should be a networkx Graph object.")

    cycle_basis = nx.cycle_basis(graph)

    # print('CYCLE BASIS')
    # print(cycle_basis)

    if len(cycle_basis) < 2:
        return {SimpleGraphCycle(c) for c in cycle_basis}

    all_edges_dict = {}
    index2edge = []
    nedges = 0
    for n1, n2 in graph.edges:
        all_edges_dict[(n1, n2)] = nedges
        all_edges_dict[(n2, n1)] = nedges
        index2edge.append((n1, n2))
        nedges += 1
    cycles_matrix = np.zeros(shape=(len(cycle_basis), nedges), dtype=np.bool)
    for icycle, cycle in enumerate(cycle_basis):
        for in1, n1 in enumerate(cycle):
            n2 = cycle[(in1 + 1) % len(cycle)]
            iedge = all_edges_dict[(n1, n2)]
            cycles_matrix[icycle, iedge] = True

    # print(cycles_matrix)
    elementary_cycles_list = []

    for ncycles in range(1, len(cycle_basis) + 1):
        for cycles_combination in itertools.combinations(cycles_matrix, ncycles):
            edges_counts = np.array(np.mod(np.sum(cycles_combination, axis=0), 2), dtype=np.bool)
            myedges = [edge for iedge, edge in enumerate(index2edge) if edges_counts[iedge]]
            # print(myedges)
            try:
                sgc = SimpleGraphCycle.from_edges(myedges, edges_are_ordered=False)
                # print(sgc)
            except ValueError as ve:
                msg = ve.args[0]
                if msg == "SimpleGraphCycle is not valid : Duplicate nodes.":
                    continue
                if msg == "Could not construct a cycle from edges.":
                    continue
                raise
            elementary_cycles_list.append(sgc)

    return elementary_cycles_list
