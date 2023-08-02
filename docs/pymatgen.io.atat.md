---
layout: default
title: pymatgen.io.atat.md
nav_exclude: true
---

# pymatgen.io.atat module

Classes for reading/writing mcsqs files following the rndstr.in format.


### _class_ pymatgen.io.atat.Mcsqs(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Bases: `object`

Handle input/output for the crystal definition format
used by mcsqs and other ATAT codes.


* **Parameters**

    **Structure** – input Structure.



#### _static_ structure_from_str(data)
Parses a rndstr.in, lat.in or bestsqs.out file into pymatgen’s
Structure format.


* **Parameters**

    **data** – contents of a rndstr.in, lat.in or bestsqs.out file



* **Returns**

    Structure object



#### structure_from_string(\*\*kwds)
structure_from_string is deprecated!
Use from_str instead


#### to_str()

* **Returns**

    a structure in mcsqs rndstr.in format.



* **Return type**

    str



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead