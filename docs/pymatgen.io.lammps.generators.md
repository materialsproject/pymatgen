---
layout: default
title: pymatgen.io.lammps.generators.md
nav_exclude: true
---

# pymatgen.io.lammps.generators module

Input set generators for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


### _class_ pymatgen.io.lammps.generators.BaseLammpsGenerator(template: str = <factory>, settings: dict = <factory>, calc_type: str = 'lammps', keep_stages: bool = True)
Bases: [`InputGenerator`](pymatgen.io.core.md#pymatgen.io.core.InputGenerator)

Base class to generate LAMMPS input sets.
Uses template files for the input. The variables that can be changed
in the input template file are those starting with a $ sign, e.g., $nsteps.
This generic class is specialized for each template in subclasses, e.g. LammpsMinimization.
You can create a template for your own task following those present in pymatgen/io/lammps/templates.
The parameters are then replaced based on the values found
in the settings dictionary that you provide, e.g., {“nsteps”: 1000}.


* **Parameters**


    * **template** – Path (string) to the template file used to create the InputFile for LAMMPS.


    * **calc_type** – Human-readable string used to briefly describe the type of computations performed by LAMMPS.


    * **settings** – Dictionary containing the values of the parameters to replace in the template.


    * **keep_stages** – If True, the string is formatted in a block structure with stage names
    and newlines that differentiate commands in the respective stages of the InputFile.
    If False, stage names are not printed and all commands appear in a single block.


/!This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


#### calc_type(_: st_ _ = 'lammps_ )

#### get_input_set(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [LammpsData](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData) | [CombinedData](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData) | None)
Generate a LammpsInputSet from the structure/data, tailored to the template file.


#### keep_stages(_: boo_ _ = Tru_ )

#### settings(_: dic_ )

#### template(_: st_ )

### _class_ pymatgen.io.lammps.generators.LammpsMinimization(template: str | None = None, units: str = 'metal', atom_style: str = 'full', dimension: int = 3, boundary: str = 'p p p', read_data: str = 'system.data', force_field: str = 'Unspecified force field!', keep_stages: bool = False)
Bases: `BaseLammpsGenerator`

Generator that yields a LammpsInputSet tailored for minimizing the energy of a system by iteratively
adjusting atom coordinates.
Example usage:
`\`
structure = Structure.from_file("mp-149.cif")
lmp_minimization = LammpsMinimization(units="atomic").get_input_set(structure)
\``.

Do not forget to specify the force field, otherwise LAMMPS will not be able to run!

/!This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


* **Parameters**


    * **template** – Path (string) to the template file used to create the InputFile for LAMMPS.


    * **units** – units to be used for the LAMMPS calculation (see official LAMMPS documentation).


    * **atom_style** – atom_style to be used for the LAMMPS calculation (see official LAMMPS documentation).


    * **dimension** – dimension to be used for the LAMMPS calculation (see official LAMMPS documentation).


    * **boundary** – boundary to be used for the LAMMPS calculation (see official LAMMPS documentation).


    * **read_data** – read_data to be used for the LAMMPS calculation (see official LAMMPS documentation).


    * **force_field** – force field to be used for the LAMMPS calculation (see official LAMMPS documentation).
    Note that you should provide all the required information as a single string.
    In case of multiple lines expected in the input file,
    separate them with ‘n’ in force_field.


    * **keep_stages** – If True, the string is formatted in a block structure with stage names
    and newlines that differentiate commands in the respective stages of the InputFile.
    If False, stage names are not printed and all commands appear in a single block.



#### _property_ atom_style()
Return the argument of the command ‘atom_style’ passed to the generator.


#### _property_ boundary()
Return the argument of the command ‘boundary’ passed to the generator.


#### _property_ dimension()
Return the argument of the command ‘dimension’ passed to the generator.


#### _property_ force_field()
Return the details of the force field commands passed to the generator.


#### _property_ read_data()
Return the argument of the command ‘read_data’ passed to the generator.


#### settings(_: dic_ )

#### template(_: st_ )

#### _property_ units()
Return the argument of the command ‘units’ passed to the generator.