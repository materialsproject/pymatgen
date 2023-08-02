---
layout: default
title: pymatgen.analysis.chemenv.connectivity.connectivity_finder.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.connectivity.connectivity_finder module

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