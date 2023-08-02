---
layout: default
title: pymatgen.io.cssr.md
nav_exclude: true
---

# pymatgen.io.cssr module

This module provides input and output from the CSSR file format.


### _class_ pymatgen.io.cssr.Cssr(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Bases: `object`

Basic object for working with Cssr file. Right now, only conversion from
a Structure to a Cssr file is supported.


* **Parameters**

    **structure** (*Structure/IStructure*) – A structure to create the Cssr object.



#### _static_ from_file(filename)
Reads a CSSR file to a Cssr object.


* **Parameters**

    **filename** (*str*) – Filename to read from.



* **Returns**

    Cssr object.



#### _static_ from_str(string)
Reads a string representation to a Cssr object.


* **Parameters**

    **string** (*str*) – A string representation of a CSSR.



* **Returns**

    Cssr object.



#### write_file(filename)
Write out a CSSR file.


* **Parameters**

    **filename** (*str*) – Filename to write to.