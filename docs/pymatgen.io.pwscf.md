---
layout: default
title: pymatgen.io.pwscf.md
nav_exclude: true
---

# pymatgen.io.pwscf module

This module implements input and output processing from PWSCF.


### _class_ pymatgen.io.pwscf.PWInput(structure, pseudo=None, control=None, system=None, electrons=None, ions=None, cell=None, kpoints_mode='automatic', kpoints_grid=(1, 1, 1), kpoints_shift=(0, 0, 0))
Bases: `object`

Base input file class. Right now, only supports no symmetry and is
very basic.

Initializes a PWSCF input file.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure. For spin-polarized calculation,
    properties (e.g. {“starting_magnetization”: -0.5,
    “pseudo”: “Mn.pbe-sp-van.UPF”}) on each site is needed instead of
    pseudo (dict).


    * **pseudo** (*dict*) – A dict of the pseudopotentials to use. Default to None.


    * **control** (*dict*) – Control parameters. Refer to official PWSCF doc
    on supported parameters. Default to {“calculation”: “scf”}


    * **system** (*dict*) – System parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **electrons** (*dict*) – Electron parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **ions** (*dict*) – Ions parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **cell** (*dict*) – Cell parameters. Refer to official PWSCF doc
    on supported parameters. Default to None, which means {}.


    * **kpoints_mode** (*str*) – Kpoints generation mode. Default to automatic.


    * **kpoints_grid** (*sequence*) – The kpoint grid. Default to (1, 1, 1).


    * **kpoints_shift** (*sequence*) – The shift for the kpoints. Defaults to
    (0, 0, 0).



#### as_dict()
Create a dictionary representation of a PWInput object.


* **Returns**

    dict



#### _classmethod_ from_dict(pwinput_dict)
Load a PWInput object from a dictionary.


* **Parameters**

    **pwinput_dict** (*dict*) – dictionary with PWInput data



* **Returns**

    PWInput object



#### _static_ from_file(filename)
Reads an PWInput object from a file.


* **Parameters**

    **filename** (*str*) – Filename for file



* **Returns**

    PWInput object



#### _static_ from_str(string)
Reads an PWInput object from a string.


* **Parameters**

    **string** (*str*) – PWInput string



* **Returns**

    PWInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _static_ proc_val(key, val)
Static helper method to convert PWINPUT parameters to proper type, e.g.,
integers, floats, etc.


* **Parameters**


    * **key** – PWINPUT parameter key


    * **val** – Actual value of PWINPUT parameter.



#### write_file(filename)
Write the PWSCF input file.


* **Parameters**

    **filename** (*str*) – The string filename to output to.



### _exception_ pymatgen.io.pwscf.PWInputError()
Bases: `BaseException`

Error for PWInput.


### _class_ pymatgen.io.pwscf.PWOutput(filename)
Bases: `object`

Parser for PWSCF output file.


* **Parameters**

    **filename** (*str*) – Filename.



#### _property_ final_energy()
Final energy.


* **Type**

    Returns



#### get_celldm(idx: int)

* **Parameters**

    **idx** (*int*) – index.



* **Returns**

    Cell dimension along index



#### _property_ lattice_type()
Lattice type.


* **Type**

    Returns



#### patterns(_ = {'celldm1': 'celldm\\\\(1\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm2': 'celldm\\\\(2\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm3': 'celldm\\\\(3\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm4': 'celldm\\\\(4\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm5': 'celldm\\\\(5\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'celldm6': 'celldm\\\\(6\\\\)=\\\\s+([\\\\d\\\\.]+)\\\\s', 'ecut': 'kinetic\\\\-energy cutoff\\\\s+=\\\\s+([\\\\d\\\\.\\\\-]+)\\\\s+Ry', 'energies': 'total energy\\\\s+=\\\\s+([\\\\d\\\\.\\\\-]+)\\\\sRy', 'lattice_type': 'bravais\\\\-lattice index\\\\s+=\\\\s+(\\\\d+)', 'nkpts': 'number of k points=\\\\s+([\\\\d]+)'_ )

#### read_pattern(patterns, reverse=False, terminate_on_match=False, postprocess=<class 'str'>)
General pattern reading. Uses monty’s regrep method. Takes the same
arguments.


* **Parameters**


    * **patterns** (*dict*) – A dict of patterns, e.g.,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”}.


    * **reverse** (*bool*) – Read files in reverse. Defaults to false. Useful for
    large files, esp OUTCARs, especially when used with
    terminate_on_match.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


Renders accessible:

    Any attribute in patterns. For example,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”} will set the
    value of self.data[“energy”] = [[-1234], [-3453], …], to the
    results from regex and postprocess. Note that the returned
    values are lists of lists, because you can grep multiple
    items on one line.