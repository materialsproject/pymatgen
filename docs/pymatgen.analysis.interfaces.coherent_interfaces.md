---
layout: default
title: pymatgen.analysis.interfaces.coherent_interfaces.md
nav_exclude: true
---

# pymatgen.analysis.interfaces.coherent_interfaces module

This module provides classes to store, generate, and manipulate material interfaces.


### _class_ pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder(substrate_structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), film_structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), film_miller: tuple[int, int, int], substrate_miller: tuple[int, int, int], zslgen: [ZSLGenerator](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator) | None = None)
Bases: `object`

This class constructs the coherent interfaces between two crystalline slabs
Coherency is defined by matching lattices not sub-planes.


* **Parameters**


    * **substrate_structure** – structure of substrate


    * **film_structure** – structure of film


    * **film_miller** – miller index of the film layer


    * **substrate_miller** – miller index for the substrate layer


    * **zslgen** – BiDirectionalZSL if you want custom lattice matching tolerances for coherency.



#### get_interfaces(termination: tuple[str, str], gap: float = 2.0, vacuum_over_film: float = 20.0, film_thickness: float | int = 1, substrate_thickness: float | int = 1, in_layers: bool = True)
Generates interface structures given the film and substrate structure
as well as the desired terminations.


* **Parameters**


    * **termination** (*tuple**[**str**, **str**]*) – termination from self.termination list


    * **gap** (*float**, **optional*) – gap between film and substrate. Defaults to 2.0.


    * **vacuum_over_film** (*float**, **optional*) – vacuum over the top of the film. Defaults to 20.0.


    * **film_thickness** (*float** | **int**, **optional*) – the film thickness. Defaults to 1.


    * **substrate_thickness** (*float** | **int**, **optional*) – substrate thickness. Defaults to 1.


    * **in_layers** (*bool**, **optional*) – set the thickness in layer units. Defaults to True.



* **Yields**

    *Iterator[Interface]* – interfaces from slabs



### pymatgen.analysis.interfaces.coherent_interfaces.from_2d_to_3d(mat: ndarray)
Converts a 2D matrix to a 3D matrix.


### pymatgen.analysis.interfaces.coherent_interfaces.get_2d_transform(start: Sequence, end: Sequence)
Gets a 2d transformation matrix
that converts start to end.


### pymatgen.analysis.interfaces.coherent_interfaces.get_rot_3d_for_2d(film_matrix, sub_matrix)
Find transformation matrix that will rotate and strain the film to the substrate while preserving the c-axis.