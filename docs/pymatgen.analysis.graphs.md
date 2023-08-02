---
layout: default
title: pymatgen.analysis.graphs.md
nav_exclude: true
---

# pymatgen.analysis.graphs module

Module for graph representations of crystals and molecules.


### _class_ pymatgen.analysis.graphs.ConnectedSite(site, jimage, index, weight, dist)
Bases: `tuple`

Create new instance of ConnectedSite(site, jimage, index, weight, dist)


#### dist()
Alias for field number 4


#### index()
Alias for field number 2


#### jimage()
Alias for field number 1


#### site()
Alias for field number 0


#### weight()
Alias for field number 3


### _exception_ pymatgen.analysis.graphs.MolGraphSplitError()
Bases: `Exception`

Raised when a molecule graph is failed to split into two disconnected
subgraphs.


### _class_ pymatgen.analysis.graphs.MoleculeGraph(molecule, graph_data=None)
Bases: `MSONable`

This is a class for annotating a Molecule with
bond information, stored in the form of a graph. A “bond” does
not necessarily have to be a chemical bond, but can store any
kind of information that connects two Sites.

If constructing this class manually, use the with_empty_graph
method or with_local_env_strategy method (using an algorithm
provided by the local_env module, such as O’Keeffe).

This class that contains connection information:
relationships between sites represented by a Graph structure,
and an associated structure object.

This class uses the NetworkX package to store and operate
on the graph itself, but contains a lot of helper methods
to make associating a graph with a given molecule easier.

Use cases for this include storing bonding information,
NMR J-couplings, Heisenberg exchange parameters, etc.


* **Parameters**


    * **molecule** – Molecule object


    * **graph_data** – dict containing graph information in
    dict format (not intended to be constructed manually,
    see as_dict method for format)



#### add_edge(from_index, to_index, weight=None, warn_duplicates=True, edge_properties=None)
Add edge to graph.

Since physically a ‘bond’ (or other connection
between sites) doesn’t have a direction, from_index,
from_jimage can be swapped with to_index, to_jimage.

However, images will always be shifted so that
from_index < to_index and from_jimage becomes (0, 0, 0).


* **Parameters**


    * **from_index** – index of site connecting from


    * **to_index** – index of site connecting to


    * **(****float****)** (*weight*) – e.g. bond length


    * **(****bool****)** (*warn_duplicates*) – if True, will warn if
    trying to add duplicate edges (duplicate edges will not
    be added in either case)


    * **(****dict****)** (*edge_properties*) – any other information to
    store on graph edges, similar to Structure’s site_properties



* **Returns**




#### alter_edge(from_index, to_index, new_weight=None, new_edge_properties=None)
Alters either the weight or the edge_properties of
an edge in the MoleculeGraph.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **new_weight** – alter_edge does not require
    that weight be altered. As such, by default, this
    is None. If weight is to be changed, it should be a
    float.


    * **new_edge_properties** – alter_edge does not require
    that edge_properties be altered. As such, by default,
    this is None. If any edge properties are to be changed,
    it should be a dictionary of edge properties to be changed.



* **Returns**




#### as_dict()
As in `pymatgen.core.Molecule` except
with using to_dict_of_dicts from NetworkX
to store graph information.


#### break_edge(from_index, to_index, allow_reverse=False)
Remove an edge from the MoleculeGraph.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **allow_reverse** – If allow_reverse is True, then break_edge will
    attempt to break both (from_index, to_index) and, failing that,
    will attempt to break (to_index, from_index).



* **Returns**




#### build_unique_fragments()
Find all possible fragment combinations of the MoleculeGraphs (in other
words, all connected induced subgraphs).


* **Returns**




#### diff(other, strict=True)
Compares two MoleculeGraphs. Returns dict with
keys ‘self’, ‘other’, ‘both’ with edges that are
present in only one MoleculeGraph (‘self’ and
‘other’), and edges that are present in both.

The Jaccard distance is a simple measure of the
dissimilarity between two MoleculeGraphs (ignoring
edge weights), and is defined by 1 - (size of the
intersection / size of the union) of the sets of
edges. This is returned with key ‘dist’.

Important note: all node indices are in terms
of the MoleculeGraph this method is called
from, not the ‘other’ MoleculeGraph: there
is no guarantee the node indices will be the
same if the underlying Molecules are ordered
differently.


* **Parameters**


    * **other** – MoleculeGraph


    * **strict** – if False, will compare bonds
    from different Molecules, with node indices
    replaced by Species strings, will not count
    number of occurrences of bonds



* **Returns**




#### draw_graph_to_file(filename='graph', diff=None, hide_unconnected_nodes=False, hide_image_edges=True, edge_colors=False, node_labels=False, weight_labels=False, image_labels=False, color_scheme='VESTA', keep_dot=False, algo='fdp')
Draws graph using GraphViz.

The networkx graph object itself can also be drawn
with networkx’s in-built graph drawing methods, but
note that this might give misleading results for
multigraphs (edges are super-imposed on each other).

If visualization is difficult to interpret,
hide_image_edges can help, especially in larger
graphs.


* **Parameters**


    * **filename** – filename to output, will detect filetype
    from extension (any graphviz filetype supported, such as
    pdf or png)


    * **(****StructureGraph****)** (*diff*) – an additional graph to
    compare with, will color edges red that do not exist in diff
    and edges green that are in diff graph but not in the
    reference graph


    * **hide_unconnected_nodes** – if True, hide unconnected
    nodes


    * **hide_image_edges** – if True, do not draw edges that
    go through periodic boundaries


    * **(****bool****)** (*keep_dot*) – if True, use node colors to
    color edges


    * **(****bool****)** – if True, label nodes with
    species and site index


    * **(****bool****)** – if True, label edges with
    weights


    * **(****bool****)** – if True, label edges with
    their periodic images (usually only used for debugging,
    edges to periodic images always appear as dashed lines)


    * **(****str****)** (*color_scheme*) – “VESTA” or “JMOL”


    * **(****bool****)** – keep GraphViz .dot file for later
    visualization


    * **algo** – any graphviz algo, “neato” (for simple graphs)
    or “fdp” (for more crowded graphs) usually give good outputs



* **Returns**




#### _property_ edge_weight_name()
Name of the edge weight property of graph


* **Type**

    return



#### _property_ edge_weight_unit()
Units of the edge weight property of graph


* **Type**

    return



#### find_rings(including=None)
Find ring structures in the MoleculeGraph.


* **Parameters**


    * **including** (*list**[**int**]*) – list of site indices. If including is not None, then find_rings


    * **default** (*will only return those rings including the specified sites. By*) –


    * **parameter** (*this*) –


    * **None** (*is*) –


    * **returned.** (*and all rings will be*) –



* **Returns**

    Each entry will be a ring (cycle, in graph theory terms)

        including the index found in the Molecule. If there is no cycle including an index, the
        value will be an empty list.




* **Return type**

    list[list[tuple[int, int]]]



#### _classmethod_ from_dict(dct)
As in `pymatgen.core.Molecule` except
restoring graphs using from_dict_of_dicts
from NetworkX to restore graph information.


#### get_connected_sites(n)
Returns a named tuple of neighbors of site n:
periodic_site, jimage, index, weight.
Index is the index of the corresponding site
in the original structure, weight can be
None if not defined.
:param n: index of Site in Molecule
:param jimage: lattice vector of site
:return: list of ConnectedSite tuples,

> sorted by closest first.


#### get_coordination_of_site(n)
Returns the number of neighbors of site n.
In graph terms, simply returns degree
of node corresponding to site n.
:param n: index of site
:return (int):


#### get_disconnected_fragments(return_index_map: bool = False)
Determine if the MoleculeGraph is connected. If it is not, separate the
MoleculeGraph into different MoleculeGraphs, where each resulting MoleculeGraph is
a disconnected subgraph of the original. Currently, this function naively assigns
the charge of the total molecule to a single submolecule. A later effort will be
to actually accurately assign charge.


* **Parameters**

    **return_index_map** (*bool*) – If True, return a dictionary that maps the
    new indices to the original indices. Defaults to False.


NOTE: This function does not modify the original MoleculeGraph. It creates a copy,
modifies that, and returns two or more new MoleculeGraph objects.


* **Returns**

    Each MoleculeGraph is a disconnected subgraph of the original MoleculeGraph.



* **Return type**

    list[MoleculeGraph]



#### insert_node(i, species, coords, validate_proximity=False, site_properties=None, edges=None)
A wrapper around Molecule.insert(), which also incorporates the new
site into the MoleculeGraph.


* **Parameters**


    * **i** – Index at which to insert the new site


    * **species** – Species for the new site


    * **coords** – 3x1 array representing coordinates of the new site


    * **validate_proximity** – For Molecule.insert(); if True (default
    False), distance will be checked to ensure that site can be safely
    added.


    * **site_properties** – Site properties for Molecule


    * **edges** – List of dicts representing edges to be added to the
    MoleculeGraph. These edges must include the index of the new site i,
    and all indices used for these edges should reflect the
    MoleculeGraph AFTER the insertion, NOT before. Each dict should at
    least have a “to_index” and “from_index” key, and can also have a
    “weight” and a “properties” key.



* **Returns**




#### isomorphic_to(other)
Checks if the graphs of two MoleculeGraphs are isomorphic to one
another. In order to prevent problems with misdirected edges, both
graphs are converted into undirected nx.Graph objects.


* **Parameters**

    **other** – MoleculeGraph object to be compared.



* **Returns**

    bool



#### _property_ name()
Name of graph


* **Type**

    return



#### remove_nodes(indices)
A wrapper for Molecule.remove_sites().


* **Parameters**

    **indices** – list of indices in the current Molecule (and graph) to
    be removed.



* **Returns**




#### replace_group(index, func_grp, strategy, bond_order=1, graph_dict=None, strategy_params=None)
Builds off of Molecule.substitute and MoleculeGraph.substitute_group
to replace a functional group in self.molecule with a functional group.
This method also amends self.graph to incorporate the new functional
group.

TODO: Figure out how to replace into a ring structure.


* **Parameters**


    * **index** – Index of atom to substitute.


    * **func_grp** – Substituent molecule. There are three options:


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



    * **strategy** – Class from pymatgen.analysis.local_env.


    * **bond_order** – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.


    * **graph_dict** – Dictionary representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. If None, then the algorithm
    will attempt to automatically determine bonds using one of
    a list of strategies defined in pymatgen.analysis.local_env.


    * **strategy_params** – dictionary of keyword arguments for strategy.
    If None, default parameters will be used.



* **Returns**




#### set_node_attributes()
Replicates molecule site properties (specie, coords, etc.) in the
MoleculeGraph.


* **Returns**




#### sort(key: Callable[[[Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)], float] | None = None, reverse: bool = False)
Same as Molecule.sort(). Also remaps nodes in graph.


* **Parameters**


    * **key** (*callable**, **optional*) – Sort key. Defaults to None.


    * **reverse** (*bool**, **optional*) – Reverse sort order. Defaults to False.



#### split_molecule_subgraphs(bonds, allow_reverse=False, alterations=None)
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

> representing bonds to be broken to split the MoleculeGraph.


* **Parameters**


    * **alterations** – a dict {(from_index, to_index): alt},
    where alt is a dictionary including weight and/or edge
    properties to be changed following the split.


    * **allow_reverse** – If allow_reverse is True, then break_edge will
    attempt to break both (from_index, to_index) and, failing that,
    will attempt to break (to_index, from_index).



* **Returns**

    list of MoleculeGraphs.



#### substitute_group(index, func_grp, strategy, bond_order=1, graph_dict=None, strategy_params=None)
Builds off of Molecule.substitute to replace an atom in self.molecule
with a functional group. This method also amends self.graph to
incorporate the new functional group.

NOTE: using a MoleculeGraph will generally produce a different graph
compared with using a Molecule or str (when not using graph_dict).


* **Parameters**


    * **index** – Index of atom to substitute.


    * **func_grp** – Substituent molecule. There are three options:


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



    * **strategy** – Class from pymatgen.analysis.local_env.


    * **bond_order** – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.


    * **graph_dict** – Dictionary representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. If None, then the algorithm
    will attempt to automatically determine bonds using one of
    a list of strategies defined in pymatgen.analysis.local_env.


    * **strategy_params** – dictionary of keyword arguments for strategy.
    If None, default parameters will be used.



* **Returns**




#### _static_ with_edges(molecule, edges)
Constructor for MoleculeGraph, using pre-existing or pre-defined edges
with optional edge parameters.


* **Parameters**


    * **molecule** – Molecule object


    * **edges** – dict representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. Props should be None if no
    additional properties are to be specified.



* **Returns**

    mg, a MoleculeGraph



#### _classmethod_ with_empty_graph(molecule, name='bonds', edge_weight_name=None, edge_weight_units=None)
Constructor for MoleculeGraph, returns a MoleculeGraph
object with an empty graph (no edges, only nodes defined
that correspond to Sites in Molecule).


* **Parameters**


    * **(****Molecule****)** (*molecule*) –


    * **(****str****)** (*edge_weight_units*) – name of graph, e.g. “bonds”


    * **(****str****)** – name of edge weights,
    e.g. “bond_length” or “exchange_constant”


    * **(****str****)** – name of edge weight units
    e.g. “Å” or “eV”



* **Return (MoleculeGraph)**



#### _static_ with_local_env_strategy(molecule, strategy)
Constructor for MoleculeGraph, using a strategy
from [`pymatgen.analysis.local_env`](pymatgen.analysis.local_env.md#module-pymatgen.analysis.local_env).


* **Parameters**


    * **molecule** – Molecule object


    * **strategy** – an instance of a
    [`pymatgen.analysis.local_env.NearNeighbors`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors) object



* **Returns**

    mg, a MoleculeGraph



### _class_ pymatgen.analysis.graphs.StructureGraph(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), graph_data=None)
Bases: `MSONable`

This is a class for annotating a Structure with
bond information, stored in the form of a graph. A “bond” does
not necessarily have to be a chemical bond, but can store any
kind of information that connects two Sites.

If constructing this class manually, use the with_empty_graph
method or with_local_env_strategy method (using an algorithm
provided by the local_env module, such as O’Keeffe).

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


* **Parameters**


    * **structure** – a Structure object


    * **graph_data** – dict containing graph information in
    dict format (not intended to be constructed manually,


see as_dict method for format)


#### add_edge(from_index, to_index, from_jimage=(0, 0, 0), to_jimage=None, weight=None, warn_duplicates=True, edge_properties=None)
Add edge to graph.

Since physically a ‘bond’ (or other connection
between sites) doesn’t have a direction, from_index,
from_jimage can be swapped with to_index, to_jimage.

However, images will always be shifted so that
from_index < to_index and from_jimage becomes (0, 0, 0).


* **Parameters**


    * **from_index** – index of site connecting from


    * **to_index** – index of site connecting to


    * **ints****)** (*to_jimage** (**tuple of*) – lattice vector of periodic
    image, e.g. (1, 0, 0) for periodic image in +x direction


    * **ints****)** – lattice vector of image


    * **(****float****)** (*weight*) – e.g. bond length


    * **(****bool****)** (*warn_duplicates*) – if True, will warn if
    trying to add duplicate edges (duplicate edges will not
    be added in either case)


    * **(****dict****)** (*edge_properties*) – any other information to
    store on graph edges, similar to Structure’s site_properties



* **Returns**




#### alter_edge(from_index, to_index, to_jimage=None, new_weight=None, new_edge_properties=None)
Alters either the weight or the edge_properties of
an edge in the StructureGraph.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **to_jimage** – tuple


    * **new_weight** – alter_edge does not require
    that weight be altered. As such, by default, this
    is None. If weight is to be changed, it should be a
    float.


    * **new_edge_properties** – alter_edge does not require
    that edge_properties be altered. As such, by default,
    this is None. If any edge properties are to be changed,
    it should be a dictionary of edge properties to be changed.



* **Returns**




#### as_dict()
As in `pymatgen.core.Structure` except
with using to_dict_of_dicts from NetworkX
to store graph information.


#### break_edge(from_index, to_index, to_jimage=None, allow_reverse=False)
Remove an edge from the StructureGraph. If no image is given, this method will fail.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **to_jimage** – tuple


    * **allow_reverse** – If allow_reverse is True, then break_edge will
    attempt to break both (from_index, to_index) and, failing that,
    will attempt to break (to_index, from_index).



* **Returns**




#### diff(other, strict=True)
Compares two StructureGraphs. Returns dict with
keys ‘self’, ‘other’, ‘both’ with edges that are
present in only one StructureGraph (‘self’ and
‘other’), and edges that are present in both.

The Jaccard distance is a simple measure of the
dissimilarity between two StructureGraphs (ignoring
edge weights), and is defined by 1 - (size of the
intersection / size of the union) of the sets of
edges. This is returned with key ‘dist’.

Important note: all node indices are in terms
of the StructureGraph this method is called
from, not the ‘other’ StructureGraph: there
is no guarantee the node indices will be the
same if the underlying Structures are ordered
differently.


* **Parameters**


    * **other** – StructureGraph


    * **strict** – if False, will compare bonds
    from different Structures, with node indices
    replaced by Species strings, will not count
    number of occurrences of bonds



* **Returns**




#### draw_graph_to_file(filename='graph', diff=None, hide_unconnected_nodes=False, hide_image_edges=True, edge_colors=False, node_labels=False, weight_labels=False, image_labels=False, color_scheme='VESTA', keep_dot=False, algo='fdp')
Draws graph using GraphViz.

The networkx graph object itself can also be drawn
with networkx’s in-built graph drawing methods, but
note that this might give misleading results for
multigraphs (edges are super-imposed on each other).

If visualization is difficult to interpret,
hide_image_edges can help, especially in larger
graphs.


* **Parameters**


    * **filename** – filename to output, will detect filetype
    from extension (any graphviz filetype supported, such as
    pdf or png)


    * **(****StructureGraph****)** (*diff*) – an additional graph to
    compare with, will color edges red that do not exist in diff
    and edges green that are in diff graph but not in the
    reference graph


    * **hide_unconnected_nodes** – if True, hide unconnected
    nodes


    * **hide_image_edges** – if True, do not draw edges that
    go through periodic boundaries


    * **(****bool****)** (*keep_dot*) – if True, use node colors to
    color edges


    * **(****bool****)** – if True, label nodes with
    species and site index


    * **(****bool****)** – if True, label edges with
    weights


    * **(****bool****)** – if True, label edges with
    their periodic images (usually only used for debugging,
    edges to periodic images always appear as dashed lines)


    * **(****str****)** (*color_scheme*) – “VESTA” or “JMOL”


    * **(****bool****)** – keep GraphViz .dot file for later
    visualization


    * **algo** – any graphviz algo, “neato” (for simple graphs)
    or “fdp” (for more crowded graphs) usually give good outputs



* **Returns**




#### _property_ edge_weight_name()
Name of the edge weight property of graph


* **Type**

    return



#### _property_ edge_weight_unit()
Units of the edge weight property of graph


* **Type**

    return



#### _classmethod_ from_dict(d)
As in `pymatgen.core.Structure` except
restoring graphs using from_dict_of_dicts
from NetworkX to restore graph information.


#### get_connected_sites(n, jimage=(0, 0, 0))
Returns a named tuple of neighbors of site n:
periodic_site, jimage, index, weight.
Index is the index of the corresponding site
in the original structure, weight can be
None if not defined.
:param n: index of Site in Structure
:param jimage: lattice vector of site
:return: list of ConnectedSite tuples,

> sorted by closest first.


#### get_coordination_of_site(n)
Returns the number of neighbors of site n.
In graph terms, simply returns degree
of node corresponding to site n.
:param n: index of site
:return (int):


#### get_subgraphs_as_molecules(use_weights=False)
Retrieve subgraphs as molecules, useful for extracting
molecules from periodic crystals.

Will only return unique molecules, not any duplicates
present in the crystal (a duplicate defined as an
isomorphic subgraph).


* **Parameters**

    **(****bool****)** (*use_weights*) – If True, only treat subgraphs
    as isomorphic if edges have the same weights. Typically,
    this means molecules will need to have the same bond
    lengths to be defined as duplicates, otherwise bond
    lengths can differ. This is a fairly robust approach,
    but will treat e.g. enantiomers as being duplicates.



* **Returns**

    list of unique Molecules in Structure



#### insert_node(i, species, coords, coords_are_cartesian=False, validate_proximity=False, site_properties=None, edges=None)
A wrapper around Molecule.insert(), which also incorporates the new
site into the MoleculeGraph.


* **Parameters**


    * **i** – Index at which to insert the new site


    * **species** – Species for the new site


    * **coords** – 3x1 array representing coordinates of the new site


    * **coords_are_cartesian** – Whether coordinates are cartesian.
    Defaults to False.


    * **validate_proximity** – For Molecule.insert(); if True (default
    False), distance will be checked to ensure that site can be safely
    added.


    * **site_properties** – Site properties for Molecule


    * **edges** – List of dicts representing edges to be added to the
    MoleculeGraph. These edges must include the index of the new site i,
    and all indices used for these edges should reflect the
    MoleculeGraph AFTER the insertion, NOT before. Each dict should at
    least have a “to_index” and “from_index” key, and can also have a
    “weight” and a “properties” key.



* **Returns**




#### _property_ name()
Name of graph


* **Type**

    return



#### remove_nodes(indices)
A wrapper for Molecule.remove_sites().


* **Parameters**

    **indices** – list of indices in the current Molecule (and graph) to
    be removed.



* **Returns**




#### set_node_attributes()
Gives each node a “specie” and a “coords” attribute, updated with the
current species and coordinates.


* **Returns**




#### sort(key=None, reverse=False)
Same as Structure.sort(). Also remaps nodes in graph.


* **Parameters**


    * **key** – key to sort by


    * **reverse** – reverse sort order



#### substitute_group(index, func_grp, strategy, bond_order=1, graph_dict=None, strategy_params=None)
Builds off of Structure.substitute to replace an atom in self.structure
with a functional group. This method also amends self.graph to
incorporate the new functional group.

NOTE: Care must be taken to ensure that the functional group that is
substituted will not place atoms to close to each other, or violate the
dimensions of the Lattice.


* **Parameters**


    * **index** – Index of atom to substitute.


    * **func_grp** – Substituent molecule. There are two options:


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



    * **strategy** – Class from pymatgen.analysis.local_env.


    * **bond_order** – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.


    * **graph_dict** – Dictionary representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. If None, then the algorithm
    will attempt to automatically determine bonds using one of
    a list of strategies defined in pymatgen.analysis.local_env.


    * **strategy_params** – dictionary of keyword arguments for strategy.
    If None, default parameters will be used.



* **Returns**




#### _property_ types_and_weights_of_connections()
Extract a dictionary summarizing the types and weights
of edges in the graph.


* **Returns**

    A dictionary with keys specifying the
    species involved in a connection in alphabetical order
    (e.g. string ‘Fe-O’) and values which are a list of
    weights for those connections (e.g. bond lengths).



#### types_of_coordination_environments(anonymous=False)
Extract information on the different co-ordination environments
present in the graph.


* **Parameters**

    **anonymous** – if anonymous, will replace specie names
    with A, B, C, etc.



* **Returns**

    a list of co-ordination environments,
    e.g. [‘Mo-S(6)’, ‘S-Mo(3)’]



#### _property_ weight_statistics()
Extract a statistical summary of edge weights present in
the graph.


* **Returns**

    A dict with an ‘all_weights’ list, ‘minimum’,
    ‘maximum’, ‘median’, ‘mean’, ‘std_dev’



#### _static_ with_edges(structure, edges)
Constructor for MoleculeGraph, using pre-existing or pre-defined edges
with optional edge parameters.


* **Parameters**


    * **molecule** – Molecule object


    * **edges** – dict representing the bonds of the functional
    group (format: {(from_index, to_index, from_image, to_image): props},
    where props is a dictionary of properties, including weight.
    Props should be None if no additional properties are to be
    specified.



* **Returns**

    sg, a StructureGraph



#### _classmethod_ with_empty_graph(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), name='bonds', edge_weight_name=None, edge_weight_units=None)
Constructor for StructureGraph, returns a StructureGraph
object with an empty graph (no edges, only nodes defined
that correspond to Sites in Structure).


* **Parameters**


    * **(****Structure****)** (*structure*) –


    * **(****str****)** (*edge_weight_units*) – name of graph, e.g. “bonds”


    * **(****str****)** – name of edge weights,
    e.g. “bond_length” or “exchange_constant”


    * **(****str****)** – name of edge weight units
    e.g. “Å” or “eV”



* **Return (StructureGraph)**



#### _static_ with_local_env_strategy(structure, strategy, weights=False, edge_properties=False)
Constructor for StructureGraph, using a strategy
from [`pymatgen.analysis.local_env`](pymatgen.analysis.local_env.md#module-pymatgen.analysis.local_env).


* **Parameters**


    * **structure** – Structure object


    * **strategy** – an instance of a
    [`pymatgen.analysis.local_env.NearNeighbors`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors) object


    * **weights** – if True, use weights from local_env class
    (consult relevant class for their meaning)


    * **edge_properties** – if True, edge_properties from neighbors will be used



* **Returns**