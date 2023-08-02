---
layout: default
title: pymatgen.symmetry.site_symmetries.md
nav_exclude: true
---

# pymatgen.symmetry.site_symmetries module

Provides analysis of site symmetries.


### pymatgen.symmetry.site_symmetries.get_shared_symmetry_operations(struct: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), pointops: list[list[[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)]], tol: float = 0.1)
Get all the point group operations shared by a pair of atomic sites
in the form [[point operations of site index 1],[],…,[]].


* **Parameters**


    * **struct** – Pymatgen structure


    * **pointops** – list of point group operations from get_site_symmetries method


    * **tol** (*float*) – tolerance to find symmetry operations



* **Returns**

    list of lists of shared point operations for each pair of atomic sites



### pymatgen.symmetry.site_symmetries.get_site_symmetries(struct: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), precision: float = 0.1)
Get all the point group operations centered on each atomic site
in the form [[point operations of site index 1]…[[point operations of site index N]]].


* **Parameters**


    * **struct** – Pymatgen structure


    * **precision** (*float*) – tolerance to find symmetry operations



* **Returns**

    list of lists of point operations for each atomic site