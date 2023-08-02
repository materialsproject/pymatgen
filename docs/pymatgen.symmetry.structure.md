---
layout: default
title: pymatgen.symmetry.structure.md
nav_exclude: true
---

# pymatgen.symmetry.structure module

This module implements symmetry-related structure forms.


### _class_ pymatgen.symmetry.structure.SymmetrizedStructure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), spacegroup: [SpacegroupOperations](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupOperations), equivalent_positions: Sequence[int], wyckoff_letters: Sequence[str])
Bases: [`Structure`](pymatgen.core.structure.md#pymatgen.core.structure.Structure)

This class represents a symmetrized structure, i.e. a structure
where the spacegroup and symmetry operations are defined. This class is
typically not called but instead is typically obtained by calling
pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.

<!-- attribute: equivalent_indices

indices of structure grouped by equivalency -->

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Original structure


    * **spacegroup** ([*SpacegroupOperations*](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupOperations)) – An input SpacegroupOperations from SpacegroupAnalyzer.


    * **equivalent_positions** (*list**[**int**]*) – Equivalent positions from SpacegroupAnalyzer.


    * **wyckoff_letters** (*list**[**str**]*) – Wyckoff letters.



#### as_dict()

* **Returns**

    MSONable dict



#### copy()

* **Returns**

    Copy of structure.



#### find_equivalent_sites(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite))
Finds all symmetrically equivalent sites for a particular site.


* **Parameters**

    **site** ([*PeriodicSite*](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)) – A site in the structure



* **Raises**

    **ValueError** – if site is not in the structure.



* **Returns**

    List of all symmetrically equivalent sites.



* **Return type**

    ([[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)])



#### _classmethod_ from_dict(dct)

* **Parameters**

    **d** – Dict representation



* **Returns**

    SymmetrizedStructure