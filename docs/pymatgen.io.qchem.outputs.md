---
layout: default
title: pymatgen.io.qchem.outputs.md
nav_exclude: true
---

# pymatgen.io.qchem.outputs module

Parsers for Qchem output files.


### _class_ pymatgen.io.qchem.outputs.QCOutput(filename: str)
Bases: `MSONable`

Class to parse QChem output files.


* **Parameters**

    **filename** (*str*) – Filename to parse.



#### as_dict()

* **Returns**

    MSONable dict.



#### _static_ multiple_outputs_from_file(filename, keep_sub_files=True)
Parses a QChem output file with multiple calculations
# 1.) Separates the output into sub-files

> e.g. qcout -> qcout.0, qcout.1, qcout.2 … qcout.N
> a.) Find delimiter for multiple calculations
> b.) Make separate output sub-files

2.) Creates separate QCCalcs for each one from the sub-files.


### pymatgen.io.qchem.outputs.check_for_structure_changes(mol1: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), mol2: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Compares connectivity of two molecules (using MoleculeGraph w/ OpenBabelNN).
This function will work with two molecules with different atom orderings,

> but for proper treatment, atoms should be listed in the same order.

Possible outputs include:
- no_change: the bonding in the two molecules is identical
- unconnected_fragments: the MoleculeGraph of mol1 is connected, but the

> MoleculeGraph is mol2 is not connected


* fewer_bonds: the MoleculeGraph of mol1 has more bonds (edges) than the
MoleculeGraph of mol2


* more_bonds: the MoleculeGraph of mol2 has more bonds (edges) than the
MoleculeGraph of mol1


* bond_change: this case catches any other non-identical MoleculeGraphs


* **Parameters**


    * **mol1** – Pymatgen Molecule object to be compared.


    * **mol2** – Pymatgen Molecule object to be compared.



* **Returns**

    One of [“unconnected_fragments”, “fewer_bonds”, “more_bonds”,
    “bond_change”, “no_change”]



### pymatgen.io.qchem.outputs.get_percentage(line: str, orbital: str)
Retrieve the percent character of an orbital.


* **Parameters**


    * **line** – Line containing orbital and percentage.


    * **orbital** – Type of orbital (s, p, d, f).



* **Returns**

    Percentage of character.



* **Raises**

    **n/a** –



### pymatgen.io.qchem.outputs.jump_to_header(lines: list[str], header: str)
Given a list of lines, truncate the start of the list so that the first line
of the new list contains the header.


* **Parameters**


    * **lines** – List of lines.


    * **header** – Substring to match.



* **Returns**

    Truncated lines.



* **Raises**

    **RuntimeError** –



### pymatgen.io.qchem.outputs.nbo_parser(filename: str)
Parse all the important sections of NBO output.


* **Parameters**

    **filename** – Path to QChem NBO output.



* **Returns**

    Data frames of formatted output.



* **Raises**

    **RuntimeError** –



### pymatgen.io.qchem.outputs.parse_hybridization_character(lines: list[str])
Parse the hybridization character section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frames of formatted output.



* **Raises**

    **RuntimeError** –



### pymatgen.io.qchem.outputs.parse_hyperbonds(lines: list[str])
Parse the natural populations section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frame of formatted output.



* **Raises**

    **RuntimeError** –



### pymatgen.io.qchem.outputs.parse_natural_populations(lines: list[str])
Parse the natural populations section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frame of formatted output.



* **Raises**

    **RuntimeError** –



### pymatgen.io.qchem.outputs.parse_perturbation_energy(lines: list[str])
Parse the perturbation energy section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frame of formatted output.



* **Raises**

    **RuntimeError** –



### pymatgen.io.qchem.outputs.z_int(string: str)
Convert string to integer.
If string empty, return -1.


* **Parameters**

    **string** – Input to be cast to int.



* **Returns**

    Int representation.



* **Raises**

    **n/a** –