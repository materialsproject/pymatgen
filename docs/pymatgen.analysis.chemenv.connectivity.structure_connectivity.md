---
layout: default
title: pymatgen.analysis.chemenv.connectivity.structure_connectivity.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.connectivity.structure_connectivity module

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