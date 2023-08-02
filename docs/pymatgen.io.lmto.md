---
layout: default
title: pymatgen.io.lmto.md
nav_exclude: true
---

# pymatgen.io.lmto module

Module for implementing a CTRL file object class for the Stuttgart
LMTO-ASA code. It will primarily be used to generate a pymatgen
Structure object in the pymatgen.electronic_structure.cohp.py module.


### _class_ pymatgen.io.lmto.LMTOCopl(filename='COPL', to_eV=False)
Bases: `object`

Class for reading COPL files, which contain COHP data.

<!-- attribute: cohp_data

Dict that contains the COHP data of the form:
  {bond: {"COHP": {Spin.up: cohps, Spin.down:cohps},
          "ICOHP": {Spin.up: icohps, Spin.down: icohps},
          "length": bond length} -->
<!-- attribute: efermi

The Fermi energy in Ry or eV. -->
<!-- attribute: energies

Sequence of energies in Ry or eV. -->
<!-- attribute: is_spin_polarized

Boolean to indicate if the calculation is spin polarized. -->

* **Parameters**


    * **filename** – filename of the COPL file. Defaults to “COPL”.


    * **to_eV** – LMTO-ASA gives energies in Ry. To convert energies into
    eV, set to True. Defaults to False for energies in Ry.



### _class_ pymatgen.io.lmto.LMTOCtrl(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), header=None, version='LMASA-47')
Bases: `object`

Class for parsing CTRL files from the Stuttgart LMTO-ASA code.
Currently, only HEADER, VERS and the structure can be used.


* **Parameters**


    * **structure** – The structure as a pymatgen Structure object.


    * **header** – The header for the CTRL file .
    Defaults to None.


    * **version** – The LMTO version that is used for the VERS category.
    Defaults to the newest version (4.7).



#### as_dict()
Returns the CTRL as a dictionary. “SITE” and “CLASS” are of
the form {‘CATEGORY’: {‘TOKEN’: value}}, the rest is of the
form ‘TOKEN’/’CATEGORY’: value. It gets the conventional standard
structure because primitive cells use the conventional
a-lattice parameter as the scaling factor and not the a-lattice
parameter of the primitive cell.


#### _classmethod_ from_dict(dct)
Creates a CTRL file object from a dictionary. The dictionary
must contain the items “ALAT”, PLAT” and “SITE”.

Valid dictionary items are:

    ALAT: the a-lattice parameter
    PLAT: (3x3) array for the lattice vectors
    SITE: list of dictionaries: {‘ATOM’: class label, ‘POS’: (3x1) array of fractional coordinates}
    CLASS (optional): list of unique atom labels as str
    SPCGRP (optional): space group symbol (str) or number (int)
    HEADER (optional): HEADER text as a str
    VERS (optional): LMTO version as a str


* **Parameters**

    **dct** – The CTRL file as a dictionary.



* **Returns**

    An LMTOCtrl object.



#### _classmethod_ from_file(filename='CTRL', \*\*kwargs)
Creates a CTRL file object from an existing file.


* **Parameters**

    **filename** – The name of the CTRL file. Defaults to ‘CTRL’.



* **Returns**

    An LMTOCtrl object.



#### _classmethod_ from_str(data, sigfigs=8)
Creates a CTRL file object from a string. This will mostly be
used to read an LMTOCtrl object from a CTRL file. Empty spheres
are ignored.


* **Parameters**

    **data** – String representation of the CTRL file.



* **Returns**

    An LMTOCtrl object.



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_string(sigfigs=8)
Generates the string representation of the CTRL file. This is
the minimal CTRL file necessary to execute lmhart.run.


#### write_file(filename='CTRL', \*\*kwargs)
Writes a CTRL file with structure, HEADER, and VERS that can be
used as input for lmhart.run.