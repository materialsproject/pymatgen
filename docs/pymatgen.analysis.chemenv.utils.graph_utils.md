---
layout: default
title: pymatgen.analysis.chemenv.utils.graph_utils.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.graph_utils module

This module contains some graph utils that are used in the chemenv package.


### _class_ pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle(nodes, edge_indices, validate=True, ordered=None)
Bases: `MSONable`

Class used to describe a cycle in a multigraph.

nodes are the nodes of the cycle and edge_indices are the indices of the edges in the cycle.
The nth index in edge_indices corresponds to the edge index between the nth node in nodes and
the (n+1)th node in nodes with the exception of the last one being the edge index between
the last node in nodes and the first node in nodes

Example: A cycle

    nodes:          1 - 3 - 4 - 0 - 2 - (1)
    edge_indices:     0 . 1 . 0 . 2 . 0 . (0)


* **Parameters**


    * **nodes** –


    * **edge_indices** –


    * **validate** –


    * **ordered** –



#### order(raise_on_fail=True)
Orders the SimpleGraphCycle.

The ordering is performed such that the first node is the “lowest” one
and the second node is the lowest one of the two neighbor nodes of the
first node. If raise_on_fail is set to True a RuntimeError will be
raised if the ordering fails.


* **Parameters**

    **raise_on_fail** – If set to True, will raise a RuntimeError if the
    ordering fails.



* **Returns**

    None



#### validate(check_strict_ordering=False)

* **Parameters**

    **check_strict_ordering** –



* **Returns**




### _class_ pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle(nodes, validate=True, ordered=None)
Bases: `MSONable`

Class used to describe a cycle in a simple graph (graph without multiple edges).

Note that the convention used here is the networkx convention for which simple graphs allow
to have self-loops in a simple graph.
No simple graph cycle with two nodes is possible in a simple graph. The graph will not
be validated if validate is set to False.
By default, the “ordered” parameter is None, in which case the SimpleGraphCycle will be ordered.
If the user explicitly sets ordered to False, the SimpleGraphCycle will not be ordered.


* **Parameters**


    * **nodes** –


    * **validate** –


    * **ordered** –



#### as_dict()

* **Returns**

    MSONable dict



#### _classmethod_ from_dict(d, validate=False)
Serialize from dict.
:param d:
:param validate:
:return:


#### _classmethod_ from_edges(edges, edges_are_ordered=True)
Constructs SimpleGraphCycle from a list edges.

By default, the edges list is supposed to be ordered as it will be
much faster to construct the cycle. If edges_are_ordered is set to
False, the code will automatically try to find the corresponding edge
order in the list.


#### order(raise_on_fail=True)
Orders the SimpleGraphCycle.

The ordering is performed such that the first node is the “lowest” one and the
second node is the lowest one of the two neighbor nodes of the first node. If
raise_on_fail is set to True a RuntimeError will be raised if the ordering fails.


* **Parameters**

    **raise_on_fail** (*bool*) – If set to True, will raise a RuntimeError if the ordering fails.



#### validate(check_strict_ordering=False)

* **Parameters**

    **check_strict_ordering** –



* **Returns**




### pymatgen.analysis.chemenv.utils.graph_utils.get_all_elementary_cycles(graph)

* **Parameters**

    **graph** –



* **Returns**




### pymatgen.analysis.chemenv.utils.graph_utils.get_all_simple_paths_edges(graph, source, target, cutoff=None, data=True)
Get all the simple path and edges.
:param graph:
:param source:
:param target:
:param cutoff:
:param data:
:return:


### pymatgen.analysis.chemenv.utils.graph_utils.get_delta(node1, node2, edge_data)
Get the delta.
:param node1:
:param node2:
:param edge_data:
:return: