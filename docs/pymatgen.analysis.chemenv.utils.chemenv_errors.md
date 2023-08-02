---
layout: default
title: pymatgen.analysis.chemenv.utils.chemenv_errors.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.chemenv_errors module

This module contains the error classes for the chemenv package.


### _exception_ pymatgen.analysis.chemenv.utils.chemenv_errors.AbstractChemenvError(cls, method, msg)
Bases: `Exception`

Abstract class for Chemenv errors.


* **Parameters**


    * **cls** –


    * **method** –


    * **msg** –



### _exception_ pymatgen.analysis.chemenv.utils.chemenv_errors.ChemenvError(cls, method, msg)
Bases: `Exception`

Chemenv error.


* **Parameters**


    * **cls** –


    * **method** –


    * **msg** –



### _exception_ pymatgen.analysis.chemenv.utils.chemenv_errors.EquivalentSiteSearchError(site)
Bases: `AbstractChemenvError`

Equivalent site search error.


* **Parameters**

    **site** –



### _exception_ pymatgen.analysis.chemenv.utils.chemenv_errors.NeighborsNotComputedChemenvError(site)
Bases: `AbstractChemenvError`

Neighbors not computed error.


* **Parameters**

    **site** –



### _exception_ pymatgen.analysis.chemenv.utils.chemenv_errors.SolidAngleError(cosinus)
Bases: `AbstractChemenvError`

Solid angle error.


* **Parameters**

    **cosinus** –