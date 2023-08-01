---
layout: default
title: pymatgen.analysis.chemenv.connectivity.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.connectivity package

Package for analyzing connectivity.

## Subpackages


* [pymatgen.analysis.chemenv.connectivity.tests package](pymatgen.analysis.chemenv.connectivity.tests.md)




    * [pymatgen.analysis.chemenv.connectivity.tests.test_connected_components module](pymatgen.analysis.chemenv.connectivity.tests.md#module-pymatgen.analysis.chemenv.connectivity.tests.test_connected_components)


        * [`ConnectedComponentTest`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest)


            * [`ConnectedComponentTest.test_coordination_sequences()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest.test_coordination_sequences)


            * [`ConnectedComponentTest.test_init()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest.test_init)


            * [`ConnectedComponentTest.test_periodicity()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest.test_periodicity)


            * [`ConnectedComponentTest.test_real_systems()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest.test_real_systems)


            * [`ConnectedComponentTest.test_serialization()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest.test_serialization)


            * [`ConnectedComponentTest.test_serialization_private_methods()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_connected_components.ConnectedComponentTest.test_serialization_private_methods)


    * [pymatgen.analysis.chemenv.connectivity.tests.test_environment_nodes module](pymatgen.analysis.chemenv.connectivity.tests.md#module-pymatgen.analysis.chemenv.connectivity.tests.test_environment_nodes)


        * [`EnvironmentNodesTest`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_environment_nodes.EnvironmentNodesTest)


            * [`EnvironmentNodesTest.test_as_dict()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_environment_nodes.EnvironmentNodesTest.test_as_dict)


            * [`EnvironmentNodesTest.test_equal()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_environment_nodes.EnvironmentNodesTest.test_equal)


            * [`EnvironmentNodesTest.test_str()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_environment_nodes.EnvironmentNodesTest.test_str)


    * [pymatgen.analysis.chemenv.connectivity.tests.test_structure_connectivity module](pymatgen.analysis.chemenv.connectivity.tests.md#module-pymatgen.analysis.chemenv.connectivity.tests.test_structure_connectivity)


        * [`StructureConnectivityTest`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_structure_connectivity.StructureConnectivityTest)


            * [`StructureConnectivityTest.test_serialization()`](pymatgen.analysis.chemenv.connectivity.tests.md#pymatgen.analysis.chemenv.connectivity.tests.test_structure_connectivity.StructureConnectivityTest.test_serialization)



## pymatgen.analysis.chemenv.connectivity.connected_components module

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

## pymatgen.analysis.chemenv.connectivity.connectivity_finder module

Module implementing connectivity finding.


### _class_ pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder(multiple_environments_choice=None)
Bases: `object`

Main class used to find the structure connectivity of a structure.

Constructor for the ConnectivityFinder.


* **Parameters**

    **multiple_environments_choice** – defines the procedure to apply when


the environment of a given site is described as a “mix” of more than one
coordination environments.


#### get_structure_connectivity(light_structure_environments)
Get the structure connectivity from the coordination environments provided
as an input.


* **Parameters**

    **light_structure_environments** – LightStructureEnvironments with the


relevant coordination environments in the structure
:return: a StructureConnectivity object describing the connectivity of
the environments in the structure


#### setup_parameters(multiple_environments_choice)
Setup of the parameters for the connectivity finder.

## pymatgen.analysis.chemenv.connectivity.environment_nodes module

Environment nodes module.


### _class_ pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode(central_site, i_central_site)
Bases: `MSONable`

Abstract class used to define an environment as a node in a graph.

Constructor for the AbstractEnvironmentNode object.


* **Parameters**


    * **central_site** ([*Site*](pymatgen.core.md#pymatgen.core.sites.Site)* or **subclass** of *[*Site*](pymatgen.core.md#pymatgen.core.sites.Site)) – central site as a pymatgen Site or
    subclass of Site (e.g. PeriodicSite, …).


    * **i_central_site** (*int*) – Index of the central site in the structure.



#### ATOM(_ = _ )

#### CE_NNBCES_NBCES_LIGANDS(_ = -_ )

#### COORDINATION_ENVIRONMENT(_ = _ )

#### DEFAULT_EXTENSIONS(_ = (6, 0_ )

#### LIGANDS_ARRANGEMENT(_ = _ )

#### NEIGHBORING_CES(_ = _ )

#### NEIGHBORING_COORDINATION_ENVIRONMENTS(_ = _ )

#### NEIGHBORS_LIGANDS_ARRANGEMENT(_ = _ )

#### NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_CE(_ = _ )

#### NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT(_ = _ )

#### NUMBER_OF_NEIGHBORING_CES(_ = _ )

#### NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS(_ = _ )

#### _property_ atom_symbol()
Symbol of the atom on the central site.


#### _property_ ce()
Coordination environment of this node.


#### _property_ ce_symbol()
Coordination environment of this node.


#### _abstract property_ coordination_environment()
Coordination environment of this node.


#### everything_equal(other)
Checks equality with respect to another AbstractEnvironmentNode using the index of the central site
as well as the central site itself.


#### _property_ isite()
Index of the central site.


#### _property_ mp_symbol()
Coordination environment of this node.


### _class_ pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode(central_site, i_central_site, ce_symbol)
Bases: `AbstractEnvironmentNode`

Class used to define an environment as a node in a graph.

Constructor for the EnvironmentNode object.


* **Parameters**


    * **central_site** ([*Site*](pymatgen.core.md#pymatgen.core.sites.Site)* or **subclass** of *[*Site*](pymatgen.core.md#pymatgen.core.sites.Site)) – central site as a pymatgen Site or
    subclass of Site (e.g. PeriodicSite, …).


    * **i_central_site** (*int*) – Index of the central site in the structure.


    * **ce_symbol** (*str*) – Symbol of the identified environment.



#### _property_ coordination_environment()
Coordination environment of this node.


#### everything_equal(other)
Compare with another environment node.


* **Returns**

    True if it is equal to the other node, False otherwise.



### pymatgen.analysis.chemenv.connectivity.environment_nodes.get_environment_node(central_site, i_central_site, ce_symbol)
Get the EnvironmentNode class or subclass for the given site and symbol.


* **Parameters**


    * **central_site** ([*Site*](pymatgen.core.md#pymatgen.core.sites.Site)* or **subclass** of *[*Site*](pymatgen.core.md#pymatgen.core.sites.Site)) – Central site of the environment.


    * **i_central_site** (*int*) – Index of the central site in the structure.


    * **ce_symbol** – Symbol of the environment.



* **Returns**

    An EnvironmentNode object.


## pymatgen.analysis.chemenv.connectivity.structure_connectivity module

Structure connectivity class.


### _class_ pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity(light_structure_environment, connectivity_graph=None, environment_subgraphs=None)
Bases: `MSONable`

Main class containing the connectivity of a structure.

Constructor for the StructureConnectivity object.


* **Parameters**


    * **light_structure_environment** – a LightStructureEnvironments object
    containing the relevant local environments
    for the sites in the structure.


    * **connectivity_graph** – the networkx MultiGraph if it has already been computed,
    e.g. stored in a file or dict and StructureConnectivity
    is reconstructed from that file or dict.


    * **environment_subgraphs** – the different subgraphs of environments that have
    been computed if any (as for connectivity_graph, only
    if it is reconstructed from a file or dict).



#### add_bonds(isite, site_neighbors_set)
Add the bonds for a given site index to the structure connectivity graph.


* **Parameters**


    * **isite** – Index of the site for which the bonds have to be added.


    * **site_neighbors_set** – site_neighbors_set: Neighbors set of the site



#### add_sites()
Add the sites in the structure connectivity graph.


#### as_dict()
Convert to MSONable dict.


#### environment_subgraph(environments_symbols=None, only_atoms=None)

* **Parameters**


    * **(****)** (*only_atoms*) –


    * **(****)** –


Returns:


#### _classmethod_ from_dict(d)

* **Parameters**

    **(****)** (*d*) –


Returns:


#### get_connected_components(environments_symbols=None, only_atoms=None)

#### print_links()
Print all links in the graph.


#### setup_atom_environment_subgraph(atom_environment)

#### setup_atom_environments_subgraph(atoms_environments)

#### setup_connectivity_description()

#### setup_environment_subgraph(environments_symbols, only_atoms=None)
Set up the graph for predefined environments and optionally atoms.


* **Parameters**


    * **environments_symbols** – Symbols of the environments for the environment subgraph.


    * **only_atoms** – Atoms to be considered.



#### setup_environments_subgraph(environments_symbols)

### pymatgen.analysis.chemenv.connectivity.structure_connectivity.get_delta_image(isite1, isite2, data1, data2)
Helper method to get the delta image between one environment and another
from the ligand’s delta images.