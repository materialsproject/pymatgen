---
layout: default
title: pymatgen.apps.borg.hive.md
nav_exclude: true
---

# pymatgen.apps.borg.hive module

This module define the various drones used to assimilate data.


### _class_ pymatgen.apps.borg.hive.AbstractDrone()
Bases: `MSONable`

Abstract drone class that defines the various methods that must be
implemented by drones. Because of the quirky nature of Python”s
multiprocessing, the intermediate data representations has to be in the
form of python primitives. So all objects that drones work with must be
MSONable. All drones must also implement the standard MSONable as_dict() and
from_dict API.


#### _abstract_ assimilate(path)
Assimilate data in a directory path into a pymatgen object. Because of
the quirky nature of Python’s multiprocessing, the object must support
pymatgen’s as_dict() for parallel processing.


* **Parameters**

    **path** – directory path



* **Returns**

    An assimilated object



#### _abstract_ get_valid_paths(path)
Checks if path contains valid data for assimilation, and then returns
the valid paths. The paths returned can be a list of directory or file
paths, depending on what kind of data you are assimilating. For
example, if you are assimilating VASP runs, you are only interested in
directories containing vasprun.xml files. On the other hand, if you are
interested converting all POSCARs in a directory tree to CIFs for
example, you will want the file paths.


* **Parameters**

    **path** – input path as a tuple generated from os.walk, i.e.,
    (parent, subdirs, files).



* **Returns**

    List of valid dir/file paths for assimilation



### _class_ pymatgen.apps.borg.hive.GaussianToComputedEntryDrone(inc_structure=False, parameters=None, data=None, file_extensions=('.log',))
Bases: `AbstractDrone`

GaussianToEntryDrone assimilates directories containing Gaussian output to
ComputedEntry/ComputedStructureEntry objects. By default, it is assumed
that Gaussian output files have a “.log” extension.

**NOTE**: Like the GaussianOutput class, this is still in early beta.


* **Parameters**


    * **inc_structure** (*bool*) – Set to True if you want
    ComputedStructureEntries to be returned instead of
    ComputedEntries.


    * **parameters** (*list*) – Input parameters to include. It has to be one of
    the properties supported by the GaussianOutput object. See
    `pymatgen.io.gaussianio GaussianOutput`. The parameters
    have to be one of python’s primitive types, i.e., list, dict of
    strings and integers. If parameters is None, a default set of
    parameters will be set.


    * **data** (*list*) – Output data to include. Has to be one of the properties
    supported by the GaussianOutput object. The parameters have to
    be one of python’s primitive types, i.e. list, dict of strings
    and integers. If data is None, a default set will be set.


    * **file_extensions** (*list*) – File extensions to be considered as Gaussian output files.
    Defaults to just the typical “log” extension.



#### as_dict()
Returns: MSONable dict.


#### assimilate(path)
Assimilate data in a directory path into a ComputedEntry object.


* **Parameters**

    **path** – directory path



* **Returns**

    ComputedEntry



#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict Representation.



* **Returns**

    GaussianToComputedEntryDrone



#### get_valid_paths(path)
Checks if path contains files with define extensions.


* **Parameters**

    **path** – input path as a tuple generated from os.walk, i.e.,
    (parent, subdirs, files).



* **Returns**

    List of valid dir/file paths for assimilation



### _class_ pymatgen.apps.borg.hive.SimpleVaspToComputedEntryDrone(inc_structure=False)
Bases: `VaspToComputedEntryDrone`

A simpler VaspToComputedEntryDrone. Instead of parsing vasprun.xml, it
parses only the INCAR, POTCAR, OSZICAR and KPOINTS files, which are much
smaller and faster to parse. However, much fewer properties are available
compared to the standard VaspToComputedEntryDrone.


* **Parameters**

    **inc_structure** (*bool*) – Set to True if you want
    ComputedStructureEntries to be returned instead of
    ComputedEntries. Structure will be parsed from the CONTCAR.



#### as_dict()
Returns: MSONable dict.


#### assimilate(path)
Assimilate data in a directory path into a ComputedEntry object.


* **Parameters**

    **path** – directory path



* **Returns**

    ComputedEntry



#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict Representation.



* **Returns**

    SimpleVaspToComputedEntryDrone



### _class_ pymatgen.apps.borg.hive.VaspToComputedEntryDrone(inc_structure=False, parameters=None, data=None)
Bases: `AbstractDrone`

VaspToEntryDrone assimilates directories containing VASP output to
ComputedEntry/ComputedStructureEntry objects.

There are some restrictions on the valid directory structures:


1. There can be only one vasp run in each directory.


2. Directories designated “relax1”, “relax2” are considered to be 2 parts
of an aflow style run, and only “relax2” is parsed.


3. The drone parses only the vasprun.xml file.


* **Parameters**


    * **inc_structure** (*bool*) – Set to True if you want
    ComputedStructureEntries to be returned instead of
    ComputedEntries.


    * **parameters** (*list*) – Input parameters to include. It has to be one of
    the properties supported by the Vasprun object. See
    `pymatgen.io.vasp.Vasprun`. If parameters is None,
    a default set of parameters that are necessary for typical
    post-processing will be set.


    * **data** (*list*) – Output data to include. Has to be one of the properties
    supported by the Vasprun object.



#### as_dict()
Returns: MSONABle dict.


#### assimilate(path)
Assimilate data in a directory path into a ComputedEntry object.


* **Parameters**

    **path** – directory path



* **Returns**

    ComputedEntry



#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict Representation.



* **Returns**

    VaspToComputedEntryDrone



#### get_valid_paths(path)
Checks if paths contains vasprun.xml or (POSCAR+OSZICAR).


* **Parameters**

    **path** – input path as a tuple generated from os.walk, i.e.,
    (parent, subdirs, files).



* **Returns**

    List of valid dir/file paths for assimilation