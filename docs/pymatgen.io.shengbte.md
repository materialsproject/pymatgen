---
layout: default
title: pymatgen.io.shengbte.md
nav_exclude: true
---

# pymatgen.io.shengbte module

This module implements reading and writing of ShengBTE CONTROL files.


### _class_ pymatgen.io.shengbte.Control(ngrid: list[int] | None = None, temperature: float | dict[str, float] = 300, \*\*kwargs)
Bases: `MSONable`, `dict`

Class for reading, updating, and writing ShengBTE CONTROL files.
See  [https://bitbucket.org/sousaw/shengbte/src/master/](https://bitbucket.org/sousaw/shengbte/src/master/) for more
detailed description and default values of CONTROL arguments.


* **Parameters**


    * **ngrid** – Reciprocal space grid density as a list of 3 ints.


    * **temperature** – The temperature to calculate the lattice thermal
    conductivity for. Can be given as a single float, or a dictionary
    with the keys “min”, “max”, “step”.


    * **\*\*kwargs** – Other ShengBTE parameters. Several parameters are required
    for ShengBTE to run - we have listed these parameters below:


        * nelements (int): number of different elements in the compound


        * natoms (int): number of atoms in the unit cell


        * lattvec (size 3x3 array): real-space lattice vectors, in units
    of lfactor


        * lfactor (float): unit of measurement for lattice vectors (nm).

        I.e., set to 0.1 if lattvec given in Angstrom.


        * types (size natom list): a vector of natom integers, ranging
    from 1 to nelements, assigning an element to each atom in the
    system


        * elements (size natom list): a vector of element names


        * positions (size natomx3 array): atomic positions in lattice
    coordinates


        * scell (size 3 list): supercell sizes along each crystal axis
    used for the 2nd-order force constant calculation




#### allocations_keys(_ = ('nelements', 'natoms', 'ngrid', 'norientations'_ )

#### as_dict()
Returns: MSONable dict.


#### crystal_keys(_ = ('lfactor', 'lattvec', 'types', 'elements', 'positions', 'masses', 'gfactors', 'epsilon', 'born', 'scell', 'orientations'_ )

#### data_keys(_ = ('nelements', 'natoms', 'ngrid', 'lattvec', 'types', 'elements', 'positions', 'scell'_ )

#### flags_keys(_ = ('nonanalytic', 'convergence', 'isotopes', 'autoisotopes', 'nanowires', 'onlyharmonic', 'espresso'_ )

#### _classmethod_ from_dict(control_dict: dict)
Write a CONTROL file from a Python dictionary. Description and default
parameters can be found at
[https://bitbucket.org/sousaw/shengbte/src/master/](https://bitbucket.org/sousaw/shengbte/src/master/).
Note some parameters are mandatory. Optional parameters default here to
None and will not be written to file.


* **Parameters**

    **control_dict** – A Python dictionary of ShengBTE input parameters.



#### _classmethod_ from_file(filepath: str)
Read a CONTROL namelist file and output a ‘Control’ object.


* **Parameters**

    **filepath** – Path of the CONTROL file.



* **Returns**

    ‘Control’ object with parameters instantiated.



#### _classmethod_ from_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), reciprocal_density: int | None = 50000, \*\*kwargs)
Get a ShengBTE control object from a structure.


* **Parameters**


    * **structure** – A structure object.


    * **reciprocal_density** – If not None, the q-point grid (“ngrid”) will be
    set using this density.


    * **kwargs** – Additional options to be passed to the Control constructor.
    See the docstring of the __init__ method for more details



* **Returns**

    A ShengBTE control object.



#### get_structure()
Get a pymatgen Structure from a ShengBTE control object.

The control object must have the “lattvec”, “types”, “elements”, and
“positions” settings otherwise an error will be thrown.


* **Returns**

    The structure.



#### params_keys(_ = ('t', 't_min', 't_max', 't_step', 'omega_max', 'scalebroad', 'rmin', 'rmax', 'dr', 'maxiter', 'nticks', 'eps'_ )

#### required_params(_ = ('nelements', 'natoms', 'ngrid', 'lattvec', 'types', 'elements', 'positions', 'scell'_ )

#### to_file(filename: str = 'CONTROL')
Writes ShengBTE CONTROL file from ‘Control’ object.


* **Parameters**

    **filename** – A file name.