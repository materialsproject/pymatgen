---
layout: default
title: pymatgen.analysis.thermochemistry.md
nav_exclude: true
---

# pymatgen.analysis.thermochemistry module

A module to perform experimental thermochemical data analysis.


### _class_ pymatgen.analysis.thermochemistry.ThermoData(data_type, cpdname, phaseinfo, formula, value, ref='', method='', temp_range=(298, 298), uncertainty=None)
Bases: `object`

A object container for an experimental Thermochemical Data.


* **Parameters**


    * **data_type** – The thermochemical data type. Should be one of the
    following: fH - Formation enthalpy, S - Entropy,
    A, B, C, D, E, F, G, H - variables for use in the various
    equations for generating formation enthalpies or Cp at
    various temperatures.


    * **cpdname** (*str*) – A name for the compound. For example, hematite for
    Fe2O3.


    * **phaseinfo** (*str*) – Denoting the phase. For example, “solid”, “liquid”,
    “gas” or “tetragonal”.


    * **formula** (*str*) – A proper string formula, e.g., Fe2O3


    * **value** (*float*) – The value of the data.


    * **ref** (*str*) – A reference, if any, for the data.


    * **method** (*str*) – The method by which the data was determined,
    if available.


    * **temp_range** (*[**float**, **float**]*) – Temperature range of validity for the
    data in Kelvin. Defaults to 298 K only.


    * **uncertainty** (*float*) – An uncertainty for the data, if available.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    ThermoData