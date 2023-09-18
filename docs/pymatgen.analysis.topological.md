---
layout: default
title: pymatgen.analysis.topological.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.topological package

Modules for prediciting topological properties.


## pymatgen.analysis.topological.spillage module

Code to calculate spin-orbit spillage.
Modified from JARVIS-Tools
[https://www.nature.com/articles/s41598-019-45028-y](https://www.nature.com/articles/s41598-019-45028-y)
[https://www.nature.com/articles/s41524-020-0319-4](https://www.nature.com/articles/s41524-020-0319-4).


### _class_ SOCSpillage(wf_noso='', wf_so='')
Bases: `object`

Spin-orbit spillage criteria to predict whether a material is topologically non-trival.
The spillage criteria physically signifies number of band-inverted electrons.
A non-zero, high value (generally >0.5) suggests non-trivial behavior.

Requires path to WAVECAR files with and without LSORBIT = .TRUE.


* **Parameters**


    * **wf_noso** – WAVECAR without spin-orbit coupling


    * **wf_so** – WAVECAR with spin-orbit coupling



#### _static_ isclose(n1, n2, rel_tol=1e-07)
Checking if the numbers are close enough.


#### _static_ orth(A)
Helper function to create orthonormal basis.


#### overlap_so_spinpol()
Main function to calculate SOC spillage.