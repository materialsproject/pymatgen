---
layout: default
title: pymatgen.analysis.disorder.md
nav_exclude: true
---

# pymatgen.analysis.disorder module

This module provides various methods to analyze order/disorder in materials.


### pymatgen.analysis.disorder.get_warren_cowley_parameters(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), r: float, dr: float)
Warren-Crowley parameters.


* **Parameters**


    * **structure** – Pymatgen Structure.


    * **r** – Radius


    * **dr** – Shell width



* **Returns**

    -1.0, …}



* **Return type**

    Warren-Crowley parameters in the form of a dict, e.g., {(Element Mo, Element W)