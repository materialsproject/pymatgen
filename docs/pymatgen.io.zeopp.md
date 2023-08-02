---
layout: default
title: pymatgen.io.zeopp.md
nav_exclude: true
---

# pymatgen.io.zeopp module

Module implementing classes and functions to use Zeo++
by Maciej Haranczyk.

If using this module, cite the following paper on Zeo++:
T.F. Willems, C.H. Rycroft, M. Kazi, J.C. Meza, and M. Haranczyk,
Algorithms and tools for high-throughput geometry-based analysis of crystalline porous materials,
Microporous and Mesoporous Materials, 149 (2012) 134-141.

## Zeo++ Installation Steps:

A stable version of Zeo++ can be obtained from [http://zeoplusplus.org](http://zeoplusplus.org).
Instructions can be found at [http://www.zeoplusplus.org/download.html](http://www.zeoplusplus.org/download.html)

## Zeo++ Post-Installation Checking:


1. Go to pymatgen/io/tests and run “python test_zeoio.py”
If Zeo++ python bindings are properly installed, the tests should
pass. One or two tests will be skipped.


1. Go to pymatgen/analysis/defects/tests and run
“python test_point_defects.py”. Lots of tests will be skipped if GULP
is not installed. But there should be no errors.


### _class_ pymatgen.io.zeopp.ZeoCssr(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Bases: [`Cssr`](pymatgen.io.cssr.md#pymatgen.io.cssr.Cssr)

ZeoCssr adds extra fields to CSSR sites to conform with Zeo++
input CSSR format. The coordinate system is rotated from xyz to zyx.
This change aligns the pivot axis of pymatgen (z-axis) to pivot axis
of Zeo++ (x-axis) for structural modifications.


* **Parameters**

    **structure** – A structure to create ZeoCssr object.



#### _static_ from_file(filename)
Reads a CSSR file to a ZeoCssr object.


* **Parameters**

    **filename** – Filename to read from.



* **Returns**

    ZeoCssr object.



#### _static_ from_str(string)
Reads a string representation to a ZeoCssr object.


* **Parameters**

    **string** – A string representation of a ZeoCSSR.



* **Returns**

    ZeoCssr object.



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### _class_ pymatgen.io.zeopp.ZeoVoronoiXYZ(mol)
Bases: [`XYZ`](pymatgen.io.xyz.md#pymatgen.io.xyz.XYZ)

Class to read Voronoi Nodes from XYZ file written by Zeo++.
The sites have an additional column representing the voronoi node radius.
The voronoi node radius is represented by the site property voronoi_radius.


* **Parameters**

    **mol** – Input molecule holding the voronoi node information.



#### _static_ from_file(filename)
Creates XYZ object from a file.


* **Parameters**

    **filename** – XYZ filename



* **Returns**

    XYZ object



#### _static_ from_str(contents)
Creates Zeo++ Voronoi XYZ object from a string.
from_string method of XYZ class is being redefined.


* **Parameters**

    **contents** – String representing Zeo++ Voronoi XYZ file.



* **Returns**

    ZeoVoronoiXYZ object



### pymatgen.io.zeopp.get_free_sphere_params(structure, rad_dict=None, probe_rad=0.1)
Analyze the void space in the input structure using voronoi decomposition
Calls Zeo++ for Voronoi decomposition.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **rad_dict** (*optional*) – Dictionary of radii of elements in structure.
    If not given, Zeo++ default values are used.
    Note: Zeo++ uses atomic radii of elements.
    For ionic structures, pass rad_dict with ionic radii


    * **probe_rad** (*optional*) – Sampling probe radius in Angstroms. Default is
    0.1 A



* **Returns**

    voronoi nodes as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure
    voronoi face centers as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure



### pymatgen.io.zeopp.get_high_accuracy_voronoi_nodes(structure, rad_dict, probe_rad=0.1)
Analyze the void space in the input structure using high accuracy
voronoi decomposition.
Calls Zeo++ for Voronoi decomposition.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **rad_dict** (*optional*) – Dictionary of radii of elements in structure.
    If not given, Zeo++ default values are used.
    Note: Zeo++ uses atomic radii of elements.
    For ionic structures, pass rad_dict with ionic radii


    * **probe_rad** (*optional*) – Sampling probe radius in Angstroms.
    Default is 0.1 A



* **Returns**

    voronoi nodes as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure
    voronoi face centers as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure



### pymatgen.io.zeopp.get_voronoi_nodes(structure, rad_dict=None, probe_rad=0.1)
Analyze the void space in the input structure using voronoi decomposition
Calls Zeo++ for Voronoi decomposition.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **rad_dict** (*optional*) – Dictionary of radii of elements in structure.
    If not given, Zeo++ default values are used.
    Note: Zeo++ uses atomic radii of elements.
    For ionic structures, pass rad_dict with ionic radii


    * **probe_rad** (*optional*) – Sampling probe radius in Angstroms. Default is
    0.1 A



* **Returns**

    voronoi nodes as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure
    voronoi face centers as pymatgen.core.structure.Structure within the
    unit cell defined by the lattice of input structure