---
layout: default
title: pymatgen.io.xtb.inputs.md
nav_exclude: true
---

# pymatgen.io.xtb.inputs module

Classes for writing XTB input files.


### _class_ pymatgen.io.xtb.inputs.CRESTInput(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), working_dir: str = '.', coords_filename: str | None = 'crest_in.xyz', constraints: dict[str, list[int] | float] | None = None)
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

        {“atoms”: [List of 1-indexed atoms to fix], “force_constant”:
        float]




* **Returns**

    String for .constrains file



#### write_input_files()
Write input files to working directory.