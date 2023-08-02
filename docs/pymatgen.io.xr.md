---
layout: default
title: pymatgen.io.xr.md
nav_exclude: true
---

# pymatgen.io.xr module

This module provides input and output mechanisms
for the xr file format, which is a modified CSSR
file format and, for example, used in GULP.
In particular, the module makes it easy
to remove shell positions from relaxations
that employed core-shell models.


### _class_ pymatgen.io.xr.Xr(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Bases: `object`

Basic object for working with xr files.


* **Parameters**

    **structure** (*Structure/IStructure*) – Structure object to create the
    Xr object.



#### _static_ from_file(filename, use_cores=True, thresh=0.0001)
Reads an xr-formatted file to create an Xr object.


* **Parameters**


    * **filename** (*str*) – name of file to read from.


    * **use_cores** (*bool*) – use core positions and discard shell
    positions if set to True (default). Otherwise,
    use shell positions and discard core positions.


    * **thresh** (*float*) – relative threshold for consistency check
    between cell parameters (lengths and angles) from
    header information and cell vectors, respectively.



* **Returns**

    Xr object corresponding to the input

        file.




* **Return type**

    xr (Xr)



#### _static_ from_str(string, use_cores=True, thresh=0.0001)
Creates an Xr object from a string representation.


* **Parameters**


    * **string** (*str*) – string representation of an Xr object.


    * **use_cores** (*bool*) – use core positions and discard shell
    positions if set to True (default). Otherwise,
    use shell positions and discard core positions.


    * **thresh** (*float*) – relative threshold for consistency check
    between cell parameters (lengths and angles) from
    header information and cell vectors, respectively.



* **Returns**

    Xr object corresponding to the input

        string representation.




* **Return type**

    xr (Xr)



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### write_file(filename)
Write out an xr file.


* **Parameters**

    **filename** (*str*) – name of the file to write to.