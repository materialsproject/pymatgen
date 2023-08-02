---
layout: default
title: pymatgen.io.vasp.help.md
nav_exclude: true
---

# pymatgen.io.vasp.help module

Get help with VASP parameters from VASP wiki.


### _class_ pymatgen.io.vasp.help.VaspDoc()
Bases: `object`

A VASP documentation helper.

Init for VaspDoc.


#### _classmethod_ get_help(tag, fmt='text')
Get help on a VASP tag.


* **Parameters**

    **tag** (*str*) – VASP tag, e.g., ISYM.



* **Returns**

    Help text.



#### _classmethod_ get_incar_tags()
Returns: All incar tags.


#### print_help(tag)
Print the help for a TAG.


* **Parameters**

    **tag** (*str*) – Tag used in VASP.



#### print_jupyter_help(tag)
Display HTML help in ipython notebook.


* **Parameters**

    **tag** (*str*) – Tag used in VASP.