---
layout: default
title: pymatgen.analysis.chemenv.connectivity.environment_nodes.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.connectivity.environment_nodes module

Environment nodes module.


### _class_ pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode(central_site, i_central_site)
Bases: `MSONable`

Abstract class used to define an environment as a node in a graph.

Constructor for the AbstractEnvironmentNode object.


* **Parameters**


    * **central_site** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)* or **subclass** of *[*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – central site as a pymatgen Site or
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


    * **central_site** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)* or **subclass** of *[*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – central site as a pymatgen Site or
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


    * **central_site** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)* or **subclass** of *[*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – Central site of the environment.


    * **i_central_site** (*int*) – Index of the central site in the structure.


    * **ce_symbol** – Symbol of the environment.



* **Returns**

    An EnvironmentNode object.