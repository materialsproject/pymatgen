---
layout: default
title: pymatgen.entries.exp_entries.md
nav_exclude: true
---

# pymatgen.entries.exp_entries module

This module defines Entry classes for containing experimental data.


### _class_ pymatgen.entries.exp_entries.ExpEntry(composition, thermodata, temperature=298)
Bases: [`PDEntry`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry), `MSONable`

An lightweight ExpEntry object containing experimental data for a
composition for many purposes. Extends a PDEntry so that it can be used for
phase diagram generation and reaction calculation.

Current version works only with solid phases and at 298K. Further
extensions for temperature dependence are planned.


* **Parameters**


    * **composition** – Composition of the entry. For flexibility, this can take
    the form of all the typical input taken by a Composition, including
    a {symbol: amt} dict, a string formula, and others.


    * **thermodata** – A sequence of ThermoData associated with the entry.


    * **temperature** – A temperature for the entry in Kelvin. Defaults to 298K.



#### as_dict()

* **Returns**

    MSONable dict



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    ExpEntry