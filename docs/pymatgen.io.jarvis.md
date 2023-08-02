---
layout: default
title: pymatgen.io.jarvis.md
nav_exclude: true
---

# pymatgen.io.jarvis module

This module provides conversion between the JARVIS
Atoms object and pymatgen Structure objects.


### _class_ pymatgen.io.jarvis.JarvisAtomsAdaptor()
Bases: `object`

Adaptor serves as a bridge between JARVIS Atoms and pymatgen objects.


#### _static_ get_atoms(structure)
Returns JARVIS Atoms object from pymatgen structure.


* **Parameters**

    **structure** – pymatgen.core.structure.Structure



* **Returns**

    JARVIS Atoms object



#### _static_ get_structure(atoms)
Returns pymatgen structure from JARVIS Atoms.


* **Parameters**

    **atoms** – JARVIS Atoms object



* **Returns**

    Equivalent pymatgen.core.structure.Structure