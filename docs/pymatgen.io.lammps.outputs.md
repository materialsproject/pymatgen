---
layout: default
title: pymatgen.io.lammps.outputs.md
nav_exclude: true
---

# pymatgen.io.lammps.outputs module

This module implements classes and methods for processing LAMMPS output
files (log and dump).


### _class_ pymatgen.io.lammps.outputs.LammpsDump(timestep, natoms, box, data)
Bases: `MSONable`

Object for representing dump data for a single snapshot.

Base constructor.


* **Parameters**


    * **timestep** (*int*) – Current time step.


    * **natoms** (*int*) – Total number of atoms in the box.


    * **box** ([*LammpsBox*](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsBox)) – Simulation box.


    * **data** (*pd.DataFrame*) – Dumped atomic data.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    LammpsDump



#### _classmethod_ from_str(string)
Constructor from string parsing.


* **Parameters**

    **string** (*str*) – Input string.



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### pymatgen.io.lammps.outputs.parse_lammps_dumps(file_pattern)
Generator that parses dump file(s).


* **Parameters**

    **file_pattern** (*str*) – Filename to parse. The timestep wildcard
    (e.g., dump.atom.’\*’) is supported and the files are parsed
    in the sequence of timestep.



* **Yields**

    LammpsDump for each available snapshot.



### pymatgen.io.lammps.outputs.parse_lammps_log(filename='log.lammps')
Parses log file with focus on thermo data. Both one and multi line
formats are supported. Any incomplete runs (no “Loop time” marker)
will not be parsed.

### Notes

SHAKE stats printed with thermo data are not supported yet.
They are ignored in multi line format, while they may cause
issues with dataframe parsing in one line format.


* **Parameters**

    **filename** (*str*) – Filename to parse.



* **Returns**

    [pd.DataFrame] containing thermo data for each completed run.