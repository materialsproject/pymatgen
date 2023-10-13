---
layout: default
title: pymatgen.io.xtb.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.xtb package

This package implements modules for input and output to and from CREST.


## pymatgen.io.xtb.inputs module

Classes for writing XTB input files.


### _class_ CRESTInput(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), working_dir: str = '.', coords_filename: str | None = 'crest_in.xyz', constraints: dict[str, list[int] | float] | None = None)
Bases: `MSONable`

An object representing  CREST input files.
Because CREST is controlled through command line flags and external
files, the CRESTInput class mainly consists of methods for containing
and writing external files.


* **Parameters**


    * **object****)** (*molecule** (**pymatgen Molecule*) – Input molecule, the only required CREST input.


    * **(****str****)** (*coords_filename*) – Location to write input files, defaults to current directory


    * **(****str****)** – Name of input coordinates file


    * **(****dict****)** (*constraints*) – Dictionary of common editable parameters for .constrains file.
    {“atoms”: [List of 1-indexed atoms to fix], “force_constant”:
    float]



#### _static_ constrains_template(molecule, reference_fnm, constraints)

* **Parameters**


    * **Molecule****)** (*molecule** (**pymatgen*) – Molecule the constraints will be performed on


    * **reference_fnm** – Name of file containing reference structure in same directory


    * **constraints** – Dictionary of common editable parameters for .constrains file.
    {“atoms”: [List of 1-indexed atoms to fix], “force_constant”: float]



* **Returns**

    for .constrains file



* **Return type**

    str



#### write_input_files()
Write input files to working directory.

## pymatgen.io.xtb.outputs module

Parsers for XTB output files and directories.


### _class_ CRESTOutput(output_filename, path='.')
Bases: `MSONable`

Class to parse CREST output files.

Assumes runtype is iMTD-GC [default].


* **Parameters**


    * **output_filename** (*str*) – Filename to parse


    * **path** (*str*) – Path to directory including output_filename and all
    other xtb output files (crest_best.xyz, etc.)



#### _parse_crest_output()
Parse output file and directory to extract all command line inputs

    and output files.

Sets the attributes:

    cmd_options: Dict of type {flag: value}
    sorted_structrues_energies: n x m x 2 list, for n conformers,

    > m rotamers per conformer, and tuple of
    > [Molecule, energy]

    properly_terminated: True or False if run properly terminated.