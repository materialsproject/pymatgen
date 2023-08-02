---
layout: default
title: pymatgen.io.lammps.sets.md
nav_exclude: true
---

# pymatgen.io.lammps.sets module

Input sets for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


### _class_ pymatgen.io.lammps.sets.LammpsInputSet(inputfile: [LammpsInputFile](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile) | str, data: [LammpsData](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData) | [CombinedData](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData), calc_type: str = '', template_file: str = '', keep_stages: bool = False)
Bases: [`InputSet`](pymatgen.io.core.md#pymatgen.io.core.InputSet)

Container class for all LAMMPS inputs. This class is intended to provide
general functionality that can be customized to many purposes.
InputGenerator classes elsewhere in this module are used to create
specific instances of LammpsInputSet that are tailored to specific purposes.

/!This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


* **Parameters**


    * **inputfile** – The input file containing settings.
    It can be a LammpsInputFile object or a string representation.


    * **data** – The data file containing structure and topology information.
    It can be a LammpsData or a CombinedData object.


    * **calc_type** – Human-readable string used to briefly describe the type of computations performed by LAMMPS.


    * **template_file** – Path (string) to the template file used to create the input file for LAMMPS.


    * **keep_stages** – Whether to keep the stage structure of the LammpsInputFile or not.



#### _classmethod_ from_directory(directory: str | Path, keep_stages: bool = False)
Construct a LammpsInputSet from a directory of two or more files.
TODO: accept directories with only the input file, that should include the structure as well.


* **Parameters**


    * **directory** – Directory to read input files from. It should contain at least two files:
    in.lammps for the LAMMPS input file, and system.data with the system information.


    * **keep_stages** – Whether to keep the stage structure of the LammpsInputFile or not.



#### validate()
A place to implement basic checks to verify the validity of an
input set. Can be as simple or as complex as desired.
Will raise a NotImplementedError unless overloaded by the inheriting class.