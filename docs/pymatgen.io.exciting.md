---
layout: default
title: pymatgen.io.exciting.md
nav_exclude: true
---

# pymatgen.io.exciting package

This package contains classes to parse input files from the exciting
code package.

## Subpackages


* [pymatgen.io.exciting.tests package](pymatgen.io.exciting.tests.md)




    * [pymatgen.io.exciting.tests.test_inputs module](pymatgen.io.exciting.tests.md#module-pymatgen.io.exciting.tests.test_inputs)


        * [`ExcitingInputTest`](pymatgen.io.exciting.tests.md#pymatgen.io.exciting.tests.test_inputs.ExcitingInputTest)


            * [`ExcitingInputTest.test_fromfile()`](pymatgen.io.exciting.tests.md#pymatgen.io.exciting.tests.test_inputs.ExcitingInputTest.test_fromfile)


            * [`ExcitingInputTest.test_paramdict()`](pymatgen.io.exciting.tests.md#pymatgen.io.exciting.tests.test_inputs.ExcitingInputTest.test_paramdict)


            * [`ExcitingInputTest.test_writebandstr()`](pymatgen.io.exciting.tests.md#pymatgen.io.exciting.tests.test_inputs.ExcitingInputTest.test_writebandstr)


            * [`ExcitingInputTest.test_writestring()`](pymatgen.io.exciting.tests.md#pymatgen.io.exciting.tests.test_inputs.ExcitingInputTest.test_writestring)



## pymatgen.io.exciting.inputs module

Classes for reading/manipulating/writing exciting input files.


### _class_ pymatgen.io.exciting.inputs.ExcitingInput(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), title=None, lockxyz=None)
Bases: `MSONable`

Object for representing the data stored in the structure part of the
exciting input.


#### structure()
Associated Structure.


#### title()
Optional title string.


#### lockxyz()
Lockxyz attribute for each site if available. A Nx3 array of
booleans.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object.


    * **title** (*str*) – Optional title for exciting input. Defaults to unit
    cell formula of structure. Defaults to None.


    * **lockxyz** (*Nx3 array*) – bool values for selective dynamics,
    where N is number of sites. Defaults to None.



#### bohr2ang(_ = 0.529169299821967_ )

#### _static_ from_file(filename)

* **Parameters**

    **filename** – Filename



* **Returns**

    ExcitingInput



#### _static_ from_str(data)
Reads the exciting input from a string.


#### _classmethod_ from_string(\*args, \*\*kwargs)
from_string is deprecated!
Use from_str instead


#### _property_ lockxyz()
Selective dynamics site properties.


* **Type**

    return



#### write_etree(celltype, cartesian=False, bandstr=False, symprec: float = 0.4, angle_tolerance=5, \*\*kwargs)
Writes the exciting input parameters to an xml object.


* **Parameters**


    * **celltype** (*str*) – Choice of unit cell. Can be either the unit cell


    * **self.structure** (*from*) –


    * **(****"conventional"****)** (*"primitive"*) –


    * **cell** (*or the primitive unit*) –


    * **cartesian** (*bool*) – Whether the atomic positions are provided in


    * **False.** (*celltype is set to "primitive". Default is*) –


    * **bandstr** (*bool*) – Whether the bandstructure path along the


    * **the** (*HighSymmKpath is included in the input file. Only supported if*) –


    * **False.** –


    * **symprec** (*float*) – Tolerance for the symmetry finding. Default is 0.4.


    * **angle_tolerance** (*float*) – Angle tolerance for the symmetry finding.


    * **5.** (*Default is*) –


    * **\*\*kwargs** – Additional parameters for the input file.



* **Returns**

    ET.Element containing the input XML structure



#### write_file(celltype, filename, cartesian=False, bandstr=False, symprec: float = 0.4, angle_tolerance=5, \*\*kwargs)
Writes exciting input file.


* **Parameters**


    * **celltype** (*str*) – Choice of unit cell. Can be either the unit cell


    * **self.structure** (*from*) –


    * **(****"conventional"****)** (*"primitive"*) –


    * **cell** (*or the primitive unit*) –


    * **filename** (*str*) – Filename for exciting input.


    * **cartesian** (*bool*) – Whether the atomic positions are provided in


    * **False.** (*celltype is set to "primitive". Default is*) –


    * **bandstr** (*bool*) – Whether the bandstructure path along the


    * **the** (*HighSymmKpath is included in the input file. Only supported if*) –


    * **False.** –


    * **symprec** (*float*) – Tolerance for the symmetry finding. Default is 0.4.


    * **angle_tolerance** (*float*) – Angle tolerance for the symmetry finding.


    * **5.** (*Default is*) –


    * **\*\*kwargs** – Additional parameters for the input file.



#### write_string(celltype, cartesian=False, bandstr=False, symprec: float = 0.4, angle_tolerance=5, \*\*kwargs)
Writes exciting input.xml as a string.


* **Parameters**


    * **celltype** (*str*) – Choice of unit cell. Can be either the unit cell


    * **self.structure** (*from*) –


    * **(****"conventional"****)** (*"primitive"*) –


    * **cell** (*or the primitive unit*) –


    * **cartesian** (*bool*) – Whether the atomic positions are provided in


    * **False.** (*celltype is set to "primitive". Default is*) –


    * **bandstr** (*bool*) – Whether the bandstructure path along the


    * **the** (*HighSymmKpath is included in the input file. Only supported if*) –


    * **False.** –


    * **symprec** (*float*) – Tolerance for the symmetry finding. Default is 0.4.


    * **angle_tolerance** (*float*) – Angle tolerance for the symmetry finding.


    * **5.** (*Default is*) –


    * **\*\*kwargs** – Additional parameters for the input file.



* **Returns**

    String