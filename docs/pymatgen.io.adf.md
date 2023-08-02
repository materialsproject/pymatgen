---
layout: default
title: pymatgen.io.adf.md
nav_exclude: true
---

# pymatgen.io.adf module

IO for ADF files.


### _class_ pymatgen.io.adf.AdfInput(task)
Bases: `object`

A basic ADF input file writer.

Initialization method.


* **Parameters**

    **task** (*AdfTask*) – An ADF task.



#### write_file(molecule, inpfile)
Write an ADF input file.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – The molecule for this task.


    * **inpfile** (*str*) – The name where the input file will be saved.



### _exception_ pymatgen.io.adf.AdfInputError()
Bases: `Exception`

The default error class for ADF.


### _class_ pymatgen.io.adf.AdfKey(name, options=None, subkeys=None)
Bases: `MSONable`

The basic input unit for ADF. A key is a string of characters that does not
contain a delimiter (blank, comma or equal sign). A key may have multiple
subkeys and a set of options.

Initialization method.


* **Parameters**


    * **name** (*str*) – The name of this key.


    * **options** (*Sized*) – The options for this key. Each element can be a primitive object or
    a tuple/list with two elements: the first is the name and the second
    is a primitive object.


    * **subkeys** (*Sized*) – The subkeys for this key.


    * **Raises** –


    * **------** –


    * **ValueError** – If elements in `subkeys` are not `AdfKey` objects.



#### add_option(option)
Add a new option to this key.


* **Parameters**


    * **option** (*Sized** or **str** or **int** or **float*) – A new option to add. This must have the same format with existing
    options.


    * **Raises** –


    * **------** –


    * **TypeError** – If the format of the given `option` is different.



#### add_subkey(subkey)
Add a new subkey to this key.


* **Parameters**


    * **subkey** (*AdfKey*) – A new subkey.


    * **Notes** –


    * **-----** –


    * **block.** (*Duplicate check will not be performed if this is an 'Atoms'*) –



#### as_dict()
A JSON-serializable dict representation of self.


#### block_keys(_ = ('SCF', 'GEOMETRY', 'XC', 'UNITS', 'ATOMS', 'CHARGE', 'BASIS', 'SYMMETRY', 'RELATIVISTIC', 'OCCUPATIONS', 'SAVE', 'A1FIT', 'INTEGRATION', 'UNRESTRICTED', 'ZLMFIT', 'TITLE', 'EXACTDENSITY', 'TOTALENERGY', 'ANALYTICALFREQ'_ )

#### _classmethod_ from_dict(d)
Construct a MSONable AdfKey object from the JSON dict.


* **Parameters**


    * **d** (*dict*) – A dict of saved attributes.


    * **Returns** –


    * **-------** –


    * **adfkey** (*AdfKey*) – An AdfKey object recovered from the JSON dict `d`.



#### _static_ from_str(string)
Construct an AdfKey object from the string.


* **Parameters**


    * **string** (*str*) – A string.


    * **Returns** –


    * **-------** –


    * **adfkey** (*AdfKey*) – An AdfKey object recovered from the string.


    * **Raises** –


    * **------** –


    * **ValueError** – Currently nested subkeys are not supported. If `subend` was found
    a ValueError would be raised.


    * **Notes** –


    * **-----** –


    * **returned.** (*Only the first block key will be*) –



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### has_option(option)
Return True if the option is included in this key.


* **Parameters**


    * **option** (*str*) – The option.


    * **Returns** –


    * **-------** –


    * **has** (*bool*) – True if the option can be found. Otherwise False will be returned.



#### has_subkey(subkey)
Return True if this AdfKey contains the given subkey.


* **Parameters**


    * **subkey** (*str** or **AdfKey*) – A key name or an AdfKey object.


    * **Returns** –


    * **-------** –


    * **has** (*bool*) – True if this key contains the given key. Otherwise False.



#### is_block_key()
Return True if this key is a block key.


#### _property_ key()
Return the name of this key. If this is a block key, the name will be
converted to upper cases.


#### remove_option(option)
Remove an option.


* **Parameters**


    * **option** (*str** or **int*) – The name (str) or index (int) of the option to remove.


    * **Raises** –


    * **------** –


    * **TypeError** – If the option has a wrong type.



#### remove_subkey(subkey)
Remove the given subkey, if existed, from this AdfKey.


* **Parameters**

    **subkey** (*str** or **AdfKey*) – The subkey to remove.



#### sub_keys(_ = ('AtomDepQuality',_ )

### _class_ pymatgen.io.adf.AdfOutput(filename)
Bases: `object`

A basic ADF output file parser.

## Attributes:

is_failed

    True is the ADF job is terminated without success. Otherwise False.

is_internal_crash

    True if the job is terminated with internal crash. Please read ‘TAPE13’
    of the ADF manual for more detail.

error

    The error description.

run_type

    The RunType of this ADF job. Possible options are: ‘SinglePoint’,
    ‘GeometryOptimization’, ‘AnalyticalFreq’ and ‘NUmericalFreq’.

final_energy

    The final molecule energy (a.u).

final_structure

    The final structure of the molecule.

energies

    The energy of each cycle.

structures

    The structure of each cycle If geometry optimization is performed.

frequencies

    The frequencies of the molecule.

normal_modes

    The normal modes of the molecule.

freq_type

    Either ‘Analytical’ or ‘Numerical’.

Initialization method.


* **param filename**

    The ADF output file to parse.



* **type filename**

    str



### _exception_ pymatgen.io.adf.AdfOutputError()
Bases: `Exception`

The default error class for errors raised by `AdfOutput`.


### _class_ pymatgen.io.adf.AdfTask(operation='energy', basis_set=None, xc=None, title='ADF_RUN', units=None, geo_subkeys=None, scf=None, other_directives=None)
Bases: `MSONable`

Basic task for ADF. All settings in this class are independent of molecules.

## Notes:

Unlike other quantum chemistry packages (NWChem, Gaussian, …), ADF does
not support calculating force/gradient.

Initialization method.


* **param operation**

    The target operation.



* **type operation**

    str



* **param basis_set**

    The basis set definitions for this task. Defaults to ‘DZ/Large’.



* **type basis_set**

    AdfKey



* **param xc**

    The exchange-correlation functionals. Defaults to PBE.



* **type xc**

    AdfKey



* **param title**

    The title of this ADF task.



* **type title**

    str



* **param units**

    The units. Defaults to Angstroms/Degree.



* **type units**

    AdfKey



* **param geo_subkeys**

    The subkeys for the block key ‘GEOMETRY’.



* **type geo_subkeys**

    Sized



* **param scf**

    The scf options.



* **type scf**

    AdfKey



* **param other_directives**

    User-defined directives.



* **type other_directives**

    Sized



#### as_dict()
A JSON-serializable dict representation of self.


#### _classmethod_ from_dict(d)
Construct a MSONable AdfTask object from the JSON dict.


* **Parameters**


    * **d** (*dict*) – A dict of saved attributes.


    * **Returns** –


    * **-------** –


    * **task** (*AdfTask*) – An AdfTask object recovered from the JSON dict `d`.



#### _static_ get_default_basis_set()
Returns: Default basis set.


#### _static_ get_default_geo()
Returns: ADFKey using default geometry.


#### _static_ get_default_scf()
Returns: ADF using default SCF.


#### _static_ get_default_units()
Returns: Default units.


#### _static_ get_default_xc()
Returns: ADFKey using default XC.


#### operations(_ = {'energy': 'Evaluate the single point energy.', 'freq': 'Same as frequencies.', 'frequencies': 'Compute second derivatives and print out an analysis of molecular vibrations.', 'numerical_frequencies': 'Compute molecular frequencies using numerical method.', 'optimize': 'Minimize the energy by varying the molecular structure.'_ )

### pymatgen.io.adf.is_numeric(s)
Return True is the string `s` is a numeric string.


* **Parameters**


    * **s** (*str*) – A string.


    * **Returns** –


    * **-------** –


    * **res** (*bool*) – If True, `s` is a numeric string and can be converted to an int or a
    float. Otherwise False will be returned.



### pymatgen.io.adf.iterlines(s: str)
A generator form of s.split(’n’) for reducing memory overhead.


* **Parameters**

    **s** (*str*) – A multi-line string.



* **Yields**

    *str* – line