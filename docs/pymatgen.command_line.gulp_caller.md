---
layout: default
title: pymatgen.command_line.gulp_caller.md
nav_exclude: true
---

# pymatgen.command_line.gulp_caller module

Interface with command line GULP.
[http://projects.ivec.org](http://projects.ivec.org)
WARNING: you need to have GULP installed on your system.


### _class_ pymatgen.command_line.gulp_caller.BuckinghamPotential(bush_lewis_flag)
Bases: `object`

Generate the Buckingham Potential Table from the bush.lib and lewis.lib.

Ref:
T.S.Bush, J.D.Gale, C.R.A.Catlow and P.D. Battle,  J. Mater Chem.,
4, 831-837 (1994).
G.V. Lewis and C.R.A. Catlow, J. Phys. C: Solid State Phys., 18,
1149-1161 (1985)


* **Parameters**

    **bush_lewis_flag** (*str*) – Flag for using Bush or Lewis potential.



### _class_ pymatgen.command_line.gulp_caller.GulpCaller(cmd='gulp')
Bases: `object`

Class to run gulp from commandline.

Initialize with the executable if not in the standard path.


* **Parameters**

    **cmd** – Command. Defaults to gulp.



#### run(gin)
Run GULP using the gin as input.


* **Parameters**

    **gin** – GULP input string



* **Returns**

    GULP output string



* **Return type**

    gout



### _exception_ pymatgen.command_line.gulp_caller.GulpConvergenceError(msg='')
Bases: `Exception`

Exception class for GULP.
Raised when proper convergence is not reached in Mott-Littleton
defect energy optimization procedure in GULP.


* **Parameters**

    **msg** (*str*) – Message.



### _exception_ pymatgen.command_line.gulp_caller.GulpError(msg)
Bases: `Exception`

Exception class for GULP.
Raised when the GULP gives an error.


* **Parameters**

    **msg** (*str*) – Message.



### _class_ pymatgen.command_line.gulp_caller.GulpIO()
Bases: `object`

To generate GULP input and process output.


#### buckingham_input(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), keywords, library=None, uc=True, valence_dict=None)
Gets a GULP input for an oxide structure and buckingham potential
from library.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **keywords** – GULP first line keywords.


    * **library** (*Default=None*) – File containing the species and potential.


    * **uc** (*Default=True*) – Unit Cell Flag.


    * **valence_dict** – {El: valence}



#### _static_ buckingham_potential(structure, val_dict=None)
Generate species, buckingham, and spring options for an oxide structure
using the parameters in default libraries.

Ref:


    1. G.V. Lewis and C.R.A. Catlow, J. Phys. C: Solid State Phys.,
    18, 1149-1161 (1985)


    2. T.S.Bush, J.D.Gale, C.R.A.Catlow and P.D. Battle,
    J. Mater Chem., 4, 831-837 (1994)


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **val_dict** (*Needed if structure is not charge neutral*) – {El:valence}
    dict, where El is element.



#### _static_ get_energy(gout: str)

* **Parameters**

    **gout** (*str*) – GULP output string.



* **Returns**

    Energy



#### _static_ get_relaxed_structure(gout: str)

* **Parameters**

    **gout** (*str*) – GULP output string.



* **Returns**

    (Structure) relaxed structure.



#### _static_ keyword_line(\*args)
Checks if the input args are proper gulp keywords and
generates the 1st line of gulp input. Full keywords are expected.


* **Parameters**

    **args** – 1st line keywords



#### _static_ library_line(file_name)
Specifies GULP library file to read species and potential parameters.
If using library don’t specify species and potential
in the input file and vice versa. Make sure the elements of
structure are in the library file.


* **Parameters**

    **file_name** – Name of GULP library file



* **Returns**

    GULP input string specifying library option



#### _static_ specie_potential_lines(structure, potential, \*\*kwargs)
Generates GULP input specie and potential string for pymatgen
structure.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure object


    * **potential** – String specifying the type of potential used


    * **kwargs** – Additional parameters related to potential. For
    potential == “buckingham”,
    anion_shell_flg (default = False):
    If True, anions are considered polarizable.
    anion_core_chrg=float
    anion_shell_chrg=float
    cation_shell_flg (default = False):
    If True, cations are considered polarizable.
    cation_core_chrg=float
    cation_shell_chrg=float



* **Returns**

    string containing specie and potential specification for gulp
    input.



#### _static_ structure_lines(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), cell_flg: bool = True, frac_flg: bool = True, anion_shell_flg: bool = True, cation_shell_flg: bool = False, symm_flg: bool = True)
Generates GULP input string corresponding to pymatgen structure.


* **Parameters**


    * **structure** – pymatgen Structure object


    * **cell_flg** (*default = True*) – Option to use lattice parameters.


    * **frac_flg** (*default = True*) – If True, fractional coordinates
    are used. Else, Cartesian coordinates in Angstroms are used.
    **\*\***
    GULP convention is to use fractional coordinates for periodic
    structures and Cartesian coordinates for non-periodic
    structures.
    **\*\***


    * **anion_shell_flg** (*default = True*) – If True, anions are considered
    polarizable.


    * **cation_shell_flg** (*default = False*) – If True, cations are
    considered polarizable.


    * **symm_flg** (*default = True*) – If True, symmetry information is also
    written.



* **Returns**

    string containing structure for GULP input



#### tersoff_input(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), periodic=False, uc=True, \*keywords)
Gets a GULP input with Tersoff potential for an oxide structure.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **periodic** (*Default=False*) – Flag denoting whether periodic
    boundary conditions are used


    * **library** (*Default=None*) – File containing the species and potential.


    * **uc** (*Default=True*) – Unit Cell Flag.


    * **keywords** – GULP first line keywords.



#### _static_ tersoff_potential(structure)
Generate the species, Tersoff potential lines for an oxide structure.


* **Parameters**

    **structure** – pymatgen.core.structure.Structure



### _class_ pymatgen.command_line.gulp_caller.TersoffPotential()
Bases: `object`

Generate Tersoff Potential Table from “OxideTersoffPotentialentials” file.

Init TersoffPotential.


### pymatgen.command_line.gulp_caller.get_energy_buckingham(structure, gulp_cmd='gulp', keywords=('optimise', 'conp', 'qok'), valence_dict=None)
Compute the energy of a structure using Buckingham potential.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **gulp_cmd** – GULP command if not in standard place


    * **keywords** – GULP first line keywords


    * **valence_dict** – {El: valence}. Needed if the structure is not charge
    neutral.



### pymatgen.command_line.gulp_caller.get_energy_relax_structure_buckingham(structure, gulp_cmd='gulp', keywords=('optimise', 'conp'), valence_dict=None)
Relax a structure and compute the energy using Buckingham potential.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **gulp_cmd** – GULP command if not in standard place


    * **keywords** – GULP first line keywords


    * **valence_dict** – {El: valence}. Needed if the structure is not charge
    neutral.



### pymatgen.command_line.gulp_caller.get_energy_tersoff(structure, gulp_cmd='gulp')
Compute the energy of a structure using Tersoff potential.


* **Parameters**


    * **structure** – pymatgen.core.structure.Structure


    * **gulp_cmd** – GULP command if not in standard place