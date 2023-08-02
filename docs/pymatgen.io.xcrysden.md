---
layout: default
title: pymatgen.io.xcrysden.md
nav_exclude: true
---

# pymatgen.io.xcrysden module

Support for reading XCrysDen files.


### _class_ pymatgen.io.xcrysden.XSF(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Bases: `object`

Class for parsing XCrysden files.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object.



#### _classmethod_ from_str(input_string, cls_=None)
Initialize a Structure object from a string with data in XSF format.


* **Parameters**


    * **input_string** – String with the structure in XSF format.
    See [http://www.xcrysden.org/doc/XSF.html](http://www.xcrysden.org/doc/XSF.html)


    * **cls** – Structure class to be created. default: pymatgen structure



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### to_str(atom_symbol=True)
Returns a string with the structure in XSF format
See [http://www.xcrysden.org/doc/XSF.html](http://www.xcrysden.org/doc/XSF.html).


* **Parameters**

    **atom_symbol** (*bool*) – Uses atom symbol instead of atomic number. Defaults to True.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead