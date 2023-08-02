---
layout: default
title: pymatgen.io.xyz.md
nav_exclude: true
---

# pymatgen.io.xyz module

Module implementing an XYZ file object class.


### _class_ pymatgen.io.xyz.XYZ(mol: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), coord_precision: int = 6)
Bases: `object`

Basic class for importing and exporting Molecules or Structures in XYZ
format.

**NOTE**: Exporting periodic structures in the XYZ format will lose information
about the periodicity. Essentially, only Cartesian coordinates are
written in this format and no information is retained about the
lattice.


* **Parameters**


    * **mol** – Input molecule or list of molecules


    * **coord_precision** – Precision to be used for coordinates.



#### _property_ all_molecules()
Returns all the frames of molecule associated with this XYZ.


#### as_dataframe()
Generates a coordinates data frame with columns: atom, x, y, and z
In case of multiple frame XYZ, returns the last frame.


* **Returns**

    pandas.DataFrame



#### _static_ from_file(filename)
Creates XYZ object from a file.


* **Parameters**

    **filename** – XYZ filename



* **Returns**

    XYZ object



#### _static_ from_str(contents)
Creates XYZ object from a string.


* **Parameters**

    **contents** – String representing an XYZ file.



* **Returns**

    XYZ object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _property_ molecule(_: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule_ )
Returns molecule associated with this XYZ. In case multiple frame
XYZ, returns the last frame.


#### write_file(filename: str)
Writes XYZ to file.


* **Parameters**

    **filename** (*str*) – File name of output file.