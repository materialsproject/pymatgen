---
layout: default
title: pymatgen.io.xtb.outputs.md
nav_exclude: true
---

# pymatgen.io.xtb.outputs module

Parsers for XTB output files and directories.


### _class_ pymatgen.io.xtb.outputs.CRESTOutput(output_filename, path='.')
Bases: `MSONable`

Class to parse CREST output files.

Assumes runtype is iMTD-GC [default].


* **Parameters**


    * **output_filename** (*str*) – Filename to parse


    * **path** (*str*) – Path to directory including output_filename and all
    other xtb output files (crest_best.xyz, etc.)