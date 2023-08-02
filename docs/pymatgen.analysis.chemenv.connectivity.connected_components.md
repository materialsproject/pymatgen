---
layout: default
title: pymatgen.analysis.chemenv.connectivity.connected_components.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.connectivity.connected_components module

Connected components.


### _class_ pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent(environments=None, links=None, environments_data=None, links_data=None, graph=None)
Bases: `MSONable`

Class used to describe the connected components in a structure in terms of coordination environments.

Constructor for the ConnectedComponent object.


* **Parameters**


    * **environments** – Environments in the connected component.


    * **links** – Links between environments in the connected component.


    * **environments_data** – Data of environment nodes.


    * **links_data** – Data of links between environment nodes.


    * **graph** – Graph of the connected component.



* **Returns**

    Instance of this class



* **Return type**

    ConnectedComponent



#### as_dict()
Bson-serializable dict representation of the ConnectedComponent object.


* **Returns**

    Bson-serializable dict representation of the ConnectedComponent object.



* **Return type**

    dict



#### compute_periodicity(algorithm='all_simple_paths')

* **Parameters**

    **(****)** (*algorithm*) –


Returns:


#### compute_periodicity_all_simple_paths_algorithm()
Get the periodicity vectors of the connected component.


#### compute_periodicity_cycle_basis()
Compute periodicity vectors of the connected component.


#### coordination_sequence(source_node, path_size=5, coordination='number', include_source=False)
Get the coordination sequence for a given node.


* **Parameters**


    * **source_node** – Node for which the coordination sequence is computed.


    * **path_size** – Maximum length of the path for the coordination sequence.


    * **coordination** – Type of coordination sequence. The default (“number”) corresponds to the number
    of environment nodes that are reachable by following paths of sizes between 1 and path_size.
    For coordination “env:number”, this resulting coordination sequence is a sequence of dictionaries
    mapping the type of environment to the number of such environment reachable by following paths of
    sizes between 1 and path_size.


    * **include_source** – Whether to include the source_node in the coordination sequence.



* **Returns**

    Mapping between the nth “layer” of the connected component with the corresponding coordination.



* **Return type**

    dict


### Examples

The corner-sharing octahedral framework (as in perovskites) have the following coordination sequence (up to
a path of size 6) :
{1: 6, 2: 18, 3: 38, 4: 66, 5: 102, 6: 146}
Considering both the octahedrons and the cuboctahedrons of the typical BaTiO3 perovskite, the “env:number”
coordination sequence (up to a path of size 6) starting on the Ti octahedron and Ba cuboctahedron
are the following :
Starting on the Ti octahedron : {1: {‘O:6’: 6, ‘C:12’: 8}, 2: {‘O:6’: 26, ‘C:12’: 48},

> 3: {‘O:6’: 90, ‘C:12’: 128}, 4: {‘O:6’: 194, ‘C:12’: 248},
> 5: {‘O:6’: 338, ‘C:12’: 408}, 6: {‘O:6’: 522, ‘C:12’: 608}}

Starting on the Ba cuboctahedron

    3: {‘O:6’: 128, ‘C:12’: 170}, 4: {‘O:6’: 248, ‘C:12’: 306},
    5: {‘O:6’: 408, ‘C:12’: 482}, 6: {‘O:6’: 608, ‘C:12’: 698}}

If include_source is set to True, the source node is included in the sequence, e.g. for the corner-sharing
octahedral framework : {0: 1, 1: 6, 2: 18, 3: 38, 4: 66, 5: 102, 6: 146}. For the “env:number” coordination
starting on a Ba cuboctahedron (as shown above), the coordination sequence is then :
{0: {‘C:12’: 1}, 1: {‘O:6’: 8, ‘C:12’: 18}, 2: {‘O:6’: 48, ‘C:12’: 74}, 3: {‘O:6’: 128, ‘C:12’: 170},

> 4: {‘O:6’: 248, ‘C:12’: 306}, 5: {‘O:6’: 408, ‘C:12’: 482}, 6: {‘O:6’: 608, ‘C:12’: 698}}


#### description(full=False)

* **Parameters**

    **full** (*bool*) – Whether to return a short or full description.



* **Returns**

    A description of the connected component.



* **Return type**

    str



#### elastic_centered_graph(start_node=None)

* **Parameters**

    **(****)** (*start_node*) –



* **Returns**

    Elastic centered subgraph.



* **Return type**

    nx.MultiGraph



#### _classmethod_ from_dict(d)
Reconstructs the ConnectedComponent object from a dict representation of the
ConnectedComponent object created using the as_dict method.


* **Parameters**

    **d** (*dict*) – dict representation of the ConnectedComponent object



* **Returns**

    The connected component representing the links of a given set of environments.



* **Return type**

    ConnectedComponent



#### _classmethod_ from_graph(g)
Constructor for the ConnectedComponent object from a graph of the connected component.


* **Parameters**

    **g** (*MultiGraph*) – Graph of the connected component.



* **Returns**

    The connected component representing the links of a given set of environments.



* **Return type**

    ConnectedComponent



#### _property_ graph()
Return the graph of this connected component.


* **Returns**

    Networkx MultiGraph object with environment as nodes and links between these nodes as edges

        with information about the image cell difference if any.




* **Return type**

    MultiGraph



#### _property_ is_0d(_: boo_ )
Whether this connected component is 0-dimensional.


#### _property_ is_1d(_: boo_ )
Whether this connected component is 1-dimensional.


#### _property_ is_2d(_: boo_ )
Whether this connected component is 2-dimensional.


#### _property_ is_3d(_: boo_ )
Whether this connected component is 3-dimensional.


#### _property_ is_periodic(_: boo_ )
Whether this connected component is periodic.


#### make_supergraph(multiplicity)

* **Parameters**

    **(****)** (*multiplicity*) –


Returns:


#### _property_ periodicity()
Get periodicity of this connected component.


#### _property_ periodicity_vectors()
Get periodicity vectors of this connected component.


#### show_graph(graph: nx.MultiGraph | None = None, save_file: str | None = None, drawing_type: str = 'internal')
Displays the graph using the specified drawing type.


* **Parameters**


    * **graph** (*Graph**, **optional*) – The graph to display. If not provided, the current graph is used.


    * **save_file** (*str**, **optional*) – The file path to save the graph image to.
    If not provided, the graph is not saved.


    * **drawing_type** (*str*) – The type of drawing to use. Can be “internal” or “external”.



### pymatgen.analysis.chemenv.connectivity.connected_components.draw_network(env_graph, pos, ax, sg=None, periodicity_vectors=None)
Draw network of environments in a matplotlib figure axes.


* **Parameters**


    * **env_graph** – Graph of environments.


    * **pos** – Positions of the nodes of the environments in the 2D figure.


    * **ax** – Axes object in which the network should be drawn.


    * **sg** – Not used currently (drawing of supergraphs).


    * **periodicity_vectors** – List of periodicity vectors that should be drawn.


Returns: None


### pymatgen.analysis.chemenv.connectivity.connected_components.make_supergraph(graph, multiplicity, periodicity_vectors)
Make supergraph from a graph of environments.


* **Parameters**


    * **graph** – Graph of environments.


    * **multiplicity** – Multiplicity of the supergraph.


    * **periodicity_vectors** – Periodicity vectors needed to make the supergraph.


Returns: Super graph of the environments.