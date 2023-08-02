---
layout: default
title: pymatgen.analysis.interfaces.substrate_analyzer.md
nav_exclude: true
---

# pymatgen.analysis.interfaces.substrate_analyzer module

This module provides classes to identify optimal substrates for film growth.


### _class_ pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer(film_max_miller=1, substrate_max_miller=1, \*\*kwargs)
Bases: [`ZSLGenerator`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator)

This class applies a set of search criteria to identify suitable
substrates for film growth. It first uses a topological search by Zur
and McGill to identify matching super-lattices on various faces of the
two materials. Additional criteria can then be used to identify the most
suitable substrate. Currently, the only additional criteria is the
elastic strain energy of the super-lattices.

Initializes the substrate analyzer
:param zslgen: Defaults to a ZSLGenerator with standard

> tolerances, but can be fed one with custom tolerances


* **Parameters**


    * **film_max_miller** (*int*) – maximum miller index to generate for film
    surfaces


    * **substrate_max_miller** (*int*) – maximum miller index to generate for
    substrate surfaces.



#### calculate(film, substrate, elasticity_tensor=None, film_millers=None, substrate_millers=None, ground_state_energy=0, lowest=False)
Finds all topological matches for the substrate and calculates elastic
strain energy and total energy for the film if elasticity tensor and
ground state energy are provided:


* **Parameters**


    * **film** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – conventional standard structure for the film


    * **substrate** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – conventional standard structure for the
    substrate


    * **elasticity_tensor** ([*ElasticTensor*](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor)) – elasticity tensor for the film
    in the IEEE orientation


    * **film_millers** (*array*) – film facets to consider in search as defined by
    miller indices


    * **substrate_millers** (*array*) – substrate facets to consider in search as
    defined by miller indices


    * **ground_state_energy** (*float*) – ground state energy for the film


    * **lowest** (*bool*) – only consider lowest matching area for each surface



#### generate_surface_vectors(film_millers, substrate_millers)
Generates the film/substrate slab combinations for a set of given
miller indices.


* **Parameters**


    * **film_millers** (*array*) – all miller indices to generate slabs for
    film


    * **substrate_millers** (*array*) – all miller indices to generate slabs
    for substrate



### _class_ pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch(film_sl_vectors: list, substrate_sl_vectors: list, film_vectors: list, substrate_vectors: list, film_transformation: list, substrate_transformation: list, film_miller: tuple[int, int, int], substrate_miller: tuple[int, int, int], strain: [Strain](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain), von_mises_strain: float, ground_state_energy: float, elastic_energy: float)
Bases: [`ZSLMatch`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch)

A substrate match building on the Zur and McGill algorithm. This match class includes the miller
planes of the film and substrate the full strain tensor, the Von Mises strain, the ground state
energy if provided, and the elastic energy.


#### elastic_energy(_: floa_ )

#### film_miller(_: tuple[int, int, int_ )

#### _classmethod_ from_zsl(match: [ZSLMatch](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch), film: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), film_miller, substrate_miller, elasticity_tensor=None, ground_state_energy=0)
Generate a substrate match from a ZSL match plus metadata.


#### ground_state_energy(_: floa_ )

#### strain(_: [Strain](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain_ )

#### substrate_miller(_: tuple[int, int, int_ )

#### _property_ total_energy()
Total energy of this match.


#### von_mises_strain(_: floa_ )