---
layout: default
title: pymatgen.io.lammps.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.lammps package

IO for LAMMPS.


## pymatgen.io.lammps.data module

This module implements a core class LammpsData for generating/parsing
LAMMPS data file, and other bridging classes to build LammpsData from
molecules. This module also implements a subclass CombinedData for
merging LammpsData object.

Only point particle styles are supported for now (atom_style in angle,
atomic, bond, charge, full and molecular only). See the pages below for
more info.

> [https://docs.lammps.org/atom_style.html](https://docs.lammps.org/atom_style.html)
> [https://docs.lammps.org/read_data.html](https://docs.lammps.org/read_data.html)


### _class_ CombinedData(list_of_molecules, list_of_names, list_of_numbers, coordinates, atom_style='full')
Bases: `LammpsData`

Object for a collective set of data for a series of LAMMPS data file.
velocities not yet implemented.


* **Parameters**


    * **list_of_molecules** – A list of LammpsData objects of a chemical cluster.
    Each LammpsData object (cluster) may contain one or more molecule ID.


    * **list_of_names** – A list of name (string) for each cluster. The characters in each name are
    restricted to word characters ([

    ```
    a-zA-Z0-9_
    ```

    ]). If names with any non-word characters
    are passed in, the special characters will be substituted by ‘_’.


    * **list_of_numbers** – A list of Integer for counts of each molecule


    * **coordinates** (*pandas.DataFrame*) – DataFrame at least containing
    columns of [“x”, “y”, “z”] for coordinates of atoms.


    * **atom_style** (*str*) – Output atom_style. Default to “full”.



#### as_lammpsdata()
Convert a CombinedData object to a LammpsData object. attributes are deep-copied.

box (LammpsBox): Simulation box.
force_field (dict): Data for force field sections. Optional

> with default to None. Only keywords in force field and
> class 2 force field are valid keys, and each value is a
> DataFrame.

topology (dict): Data for topology sections. Optional with

    default to None. Only keywords in topology are valid
    keys, and each value is a DataFrame.


#### disassemble(atom_labels=None, guess_element=True, ff_label='ff_map')
Breaks down each LammpsData in CombinedData to building blocks
(LammpsBox, ForceField and a series of Topology).
RESTRICTIONS APPLIED:
1. No complex force field defined not just on atom

> types, where the same type or equivalent types of topology
> may have more than one set of coefficients.


1. No intermolecular topologies (with atoms from different

    molecule-ID) since a Topology object includes data for ONE
    molecule or structure only.


* **Parameters**


    * **atom_labels** (*[**str**]*) – List of strings (must be different
    from one another) for labelling each atom type found in
    Masses section. Default to None, where the labels are
    automatically added based on either element guess or
    dummy specie assignment.


    * **guess_element** (*bool*) – Whether to guess the element based on
    its atomic mass. Default to True, otherwise dummy
    species “Qa”, “Qb”, … will be assigned to various
    atom types. The guessed or assigned elements will be
    reflected on atom labels if atom_labels is None, as
    well as on the species of molecule in each Topology.


    * **ff_label** (*str*) – Site property key for labeling atoms of
    different types. Default to “ff_map”.



* **Returns**

    [(LammpsBox, ForceField, [Topology]), …]



#### _classmethod_ from_ff_and_topologies()
Unsupported constructor for CombinedData objects.


#### _classmethod_ from_files(coordinate_file, list_of_numbers, \*filenames)
Constructor that parse a series of data file.


* **Parameters**


    * **coordinate_file** (*str*) – The filename of xyz coordinates.


    * **list_of_numbers** (*list*) – A list of numbers specifying counts for each
    clusters parsed from files.


    * **filenames** (*str*) – A series of LAMMPS data filenames in string format.



#### _classmethod_ from_lammpsdata(mols, names, list_of_numbers, coordinates, atom_style=None)
Constructor that can infer atom_style.
The input LammpsData objects are used non-destructively.


* **Parameters**


    * **mols** – a list of LammpsData of a chemical cluster.Each LammpsData object (cluster)
    may contain one or more molecule ID.


    * **names** – a list of name for each cluster.


    * **list_of_numbers** – a list of Integer for counts of each molecule


    * **coordinates** (*pandas.DataFrame*) – DataFrame at least containing
    columns of [“x”, “y”, “z”] for coordinates of atoms.


    * **atom_style** (*str*) – Output atom_style. Default to “full”.



#### _classmethod_ from_structure()
Unsupported constructor for CombinedData objects.


#### get_str(distance=6, velocity=8, charge=4, hybrid=True)
Returns the string representation of CombinedData, essentially
the string to be written to a file. Combination info is included
as a comment. For single molecule ID data, the info format is:

> num name

For data with multiple molecule ID, the format is:

    num(mols_per_data) name.


* **Parameters**


    * **distance** (*int*) – No. of significant figures to output for
    box settings (bounds and tilt) and atomic coordinates.
    Default to 6.


    * **velocity** (*int*) – No. of significant figures to output for
    velocities. Default to 8.


    * **charge** (*int*) – No. of significant figures to output for
    charges. Default to 4.


    * **hybrid** (*bool*) – Whether to write hybrid coeffs types.
    Default to True. If the data object has no hybrid
    coeffs types and has large coeffs section, one may
    use False to speed up the process. Otherwise, the
    default is recommended.



* **Returns**

    String representation



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### _classmethod_ parse_xyz(filename)
Load xyz file generated from packmol (for those who find it hard to install openbabel).


* **Returns**

    pandas.DataFrame



#### _property_ structure()
Exports a periodic structure object representing the simulation
box.


* **Returns**

    Structure



### _class_ ForceField(mass_info, nonbond_coeffs=None, topo_coeffs=None)
Bases: `MSONable`

Class carrying most data in Masses and force field sections.


#### masses()
DataFrame for Masses section.


* **Type**

    pandas.DataFrame



#### force_field()
Force field section keywords (keys) and
data (values) as DataFrames.


* **Type**

    dict



#### maps()
Dict for labeling atoms and topologies.


* **Type**

    dict



* **Parameters**


    * **mass_info** (*list*) – List of atomic mass info. Elements,
    strings (symbols) and floats are all acceptable for the
    values, with the first two converted to the atomic mass
    of an element. It is recommended to use
    dict.items() to prevent key duplications.
    [(“C”, 12.01), (“H”, Element(“H”)), (“O”, “O”), …]


    * **nonbond_coeffs** (*list*) – List of Pair or PairIJ
    coefficients, of which the sequence must be sorted
    according to the species in mass_dict. Pair or PairIJ
    determined by the length of list. Optional with default
    to None.


    * **topo_coeffs** (*dict*) – Dict with force field coefficients for
    molecular topologies. Optional with default
    to None. All four valid keys listed below are optional.
    Each value is a list of dicts with non-optional keys
    “coeffs” and “types”, and related class2 force field
    keywords as optional keys.
    {

    > ”Bond Coeffs”:

    >     [{“coeffs”: [coeff],

    >         ”types”: [(“C”, “C”), …]}, …],

    > ”Angle Coeffs”:

    >     [{“coeffs”: [coeff],

    >         ”BondBond Coeffs”: [coeff],
    >         “types”: [(“H”, “C”, “H”), …]}, …],

    > ”Dihedral Coeffs”:

    >     [{“coeffs”: [coeff],

    >         ”BondBond13 Coeffs”: [coeff],
    >         “types”: [(“H”, “C”, “C”, “H”), …]}, …],

    > ”Improper Coeffs”:

    >     [{“coeffs”: [coeff],

    >         ”AngleAngle Coeffs”: [coeff],
    >         “types”: [(“H”, “C”, “C”, “H”), …]}, …],

    }
    Topology of same type or equivalent types (e.g.,
    (“C”, “H”) and (“H”, “C”) bonds) are NOT ALLOWED to
    be defined MORE THAN ONCE with DIFFERENT coefficients.




#### _static_ _is_valid(df)

#### _process_nonbond()

#### _process_topo(kw)

#### _classmethod_ from_dict(d)
Constructor that reads in a dictionary.


* **Parameters**

    **d** (*dict*) – Dictionary to read.



#### _classmethod_ from_file(filename)
Constructor that reads in a file in YAML format.


* **Parameters**

    **filename** (*str*) – Filename.



#### to_file(filename)
Saves object to a file in YAML format.


* **Parameters**

    **filename** (*str*) – Filename.



### _class_ LammpsBox(bounds, tilt=None)
Bases: `MSONable`

Object for representing a simulation box in LAMMPS settings.


* **Parameters**


    * **bounds** – A (3, 2) array/list of floats setting the
    boundaries of simulation box.


    * **tilt** – A (3,) array/list of floats setting the tilt of
    simulation box. Default to None, i.e., use an
    orthogonal box.



#### get_box_shift(i)
Calculates the coordinate shift due to PBC.


* **Parameters**


    * **i** – A (n, 3) integer array containing the labels for box


    * **entries.** (*images** of **n*) –



* **Returns**

    Coorindate shift array with the same shape of i



#### get_str(significant_figures: int = 6)
Returns the string representation of simulation box in LAMMPS
data file format.


* **Parameters**

    **significant_figures** (*int*) – No. of significant figures to
    output for box settings. Default to 6.



* **Returns**

    String representation



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### to_lattice()
Converts the simulation box to a more powerful Lattice backend.
Note that Lattice is always periodic in 3D space while a
simulation box is not necessarily periodic in all dimensions.


* **Returns**

    Lattice



#### _property_ volume()
Volume of simulation box.


### _class_ LammpsData(box, masses, atoms, velocities=None, force_field=None, topology=None, atom_style='full')
Bases: `MSONable`

Object for representing the data in a LAMMPS data file.

This is a low level constructor designed to work with parsed
data or other bridging objects (ForceField and Topology). Not
recommended to use directly.


* **Parameters**


    * **box** (*LammpsBox*) – Simulation box.


    * **masses** (*pandas.DataFrame*) – DataFrame with one column
    [“mass”] for Masses section.


    * **atoms** (*pandas.DataFrame*) – DataFrame with multiple columns
    for Atoms section. Column names vary with atom_style.


    * **velocities** (*pandas.DataFrame*) – DataFrame with three columns
    [“vx”, “vy”, “vz”] for Velocities section. Optional
    with default to None. If not None, its index should be
    consistent with atoms.


    * **force_field** (*dict*) – Data for force field sections. Optional
    with default to None. Only keywords in force field and
    class 2 force field are valid keys, and each value is a
    DataFrame.


    * **topology** (*dict*) – Data for topology sections. Optional with
    default to None. Only keywords in topology are valid
    keys, and each value is a DataFrame.


    * **atom_style** (*str*) – Output atom_style. Default to “full”.



#### disassemble(atom_labels=None, guess_element=True, ff_label='ff_map')
Breaks down LammpsData to building blocks
(LammpsBox, ForceField and a series of Topology).
RESTRICTIONS APPLIED:


1. No complex force field defined not just on atom

    types, where the same type or equivalent types of topology
    may have more than one set of coefficients.


2. No intermolecular topologies (with atoms from different

    molecule-ID) since a Topology object includes data for ONE
    molecule or structure only.


* **Parameters**


    * **atom_labels** (*[**str**]*) – List of strings (must be different
    from one another) for labelling each atom type found in
    Masses section. Default to None, where the labels are
    automatically added based on either element guess or
    dummy specie assignment.


    * **guess_element** (*bool*) – Whether to guess the element based on
    its atomic mass. Default to True, otherwise dummy
    species “Qa”, “Qb”, … will be assigned to various
    atom types. The guessed or assigned elements will be
    reflected on atom labels if atom_labels is None, as
    well as on the species of molecule in each Topology.


    * **ff_label** (*str*) – Site property key for labeling atoms of
    different types. Default to “ff_map”.



* **Returns**

    LammpsBox, ForceField, [Topology]



#### _classmethod_ from_ff_and_topologies(box, ff, topologies, atom_style='full')
Constructor building LammpsData from a ForceField object and a
list of Topology objects. Do not support intermolecular
topologies since a Topology object includes data for ONE
molecule or structure only.


* **Parameters**


    * **box** (*LammpsBox*) – Simulation box.


    * **ff** (*ForceField*) – ForceField object with data for Masses and
    force field sections.


    * **topologies** (*[**Topology**]*) – List of Topology objects with data
    for Atoms, Velocities and topology sections.


    * **atom_style** (*str*) – Output atom_style. Default to “full”.



#### _classmethod_ from_file(filename, atom_style='full', sort_id=False)
Constructor that parses a file.


* **Parameters**


    * **filename** (*str*) – Filename to read.


    * **atom_style** (*str*) – Associated atom_style. Default to “full”.


    * **sort_id** (*bool*) – Whether sort each section by id. Default to
    True.



#### _classmethod_ from_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), ff_elements=None, atom_style='charge', is_sort=False)
Simple constructor building LammpsData from a structure without
force field parameters and topologies.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **ff_elements** (*[**str**]*) – List of strings of elements that must
    be present due to force field settings but not
    necessarily in the structure. Default to None.


    * **atom_style** (*str*) – Choose between “atomic” (neutral) and
    “charge” (charged). Default to “charge”.


    * **is_sort** (*bool*) – whether to sort sites



#### get_str(distance: int = 6, velocity: int = 8, charge: int = 4, hybrid: bool = True)
Returns the string representation of LammpsData, essentially
the string to be written to a file. Support hybrid style
coeffs read and write.


* **Parameters**


    * **distance** (*int*) – No. of significant figures to output for
    box settings (bounds and tilt) and atomic coordinates.
    Default to 6.


    * **velocity** (*int*) – No. of significant figures to output for
    velocities. Default to 8.


    * **charge** (*int*) – No. of significant figures to output for
    charges. Default to 4.


    * **hybrid** (*bool*) – Whether to write hybrid coeffs types.
    Default to True. If the data object has no hybrid
    coeffs types and has large coeffs section, one may
    use False to speed up the process. Otherwise, the
    default is recommended.



* **Returns**

    String representation



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### set_charge_atom(charges: dict[int, float])
Set the charges of specific atoms of the data.


* **Parameters**

    **charges** – A dictionary with atom indexes as keys and
    charges as values, e.g., to set the charge
    of the atom with index 3 to -2, use {3: -2}.



#### set_charge_atom_type(charges: dict[str | int, float])
Add or modify charges of all atoms of a given type in the data.


* **Parameters**

    **charges** – Dict containing the charges for the atom types to set.
    The dict should contain atom types as integers or labels and charges.
    Example: change the charge of Li atoms to +3:

    > charges={“Li”: 3}
    > charges={1: 3} if Li atoms are of type 1




#### _property_ structure()
Exports a periodic structure object representing the simulation
box.


* **Returns**

    Structure



#### write_file(filename, distance=6, velocity=8, charge=4)
Writes LammpsData to file.


* **Parameters**


    * **filename** (*str*) – Filename.


    * **distance** (*int*) – No. of significant figures to output for
    box settings (bounds and tilt) and atomic coordinates.
    Default to 6.


    * **velocity** (*int*) – No. of significant figures to output for
    velocities. Default to 8.


    * **charge** (*int*) – No. of significant figures to output for
    charges. Default to 4.



### _class_ Topology(sites, ff_label=None, charges=None, velocities=None, topologies=None)
Bases: `MSONable`

Class carrying most data in Atoms, Velocities and molecular
topology sections for ONE SINGLE Molecule or Structure
object, or a plain list of Sites.


* **Parameters**


    * **sites** (*[*[*Site*](pymatgen.core.md#pymatgen.core.sites.Site)*] or *[*SiteCollection*](pymatgen.core.md#pymatgen.core.structure.SiteCollection)) – A group of sites in a
    list or as a Molecule/Structure.


    * **ff_label** (*str*) – Site property key for labeling atoms of
    different types. Default to None, i.e., use
    site.species_string.


    * **charges** (*[**q**, **...**]*) – Charge of each site in a (n,)
    array/list, where n is the No. of sites. Default to
    None, i.e., search site property for charges.


    * **velocities** (*[**[**vx**, **vy**, **vz**]**, **...**]*) – Velocity of each site in
    a (n, 3) array/list, where n is the No. of sites.
    Default to None, i.e., search site property for
    velocities.


    * **topologies** (*dict*) – Bonds, angles, dihedrals and improper
    dihedrals defined by site indices. Default to None,
    i.e., no additional topology. All four valid keys
    listed below are optional.
    {

    > ”Bonds”: [[i, j], …],
    > “Angles”: [[i, j, k], …],
    > “Dihedrals”: [[i, j, k, l], …],
    > “Impropers”: [[i, j, k, l], …]

    }.




#### _classmethod_ from_bonding(molecule, bond=True, angle=True, dihedral=True, tol: float = 0.1, \*\*kwargs)
Another constructor that creates an instance from a molecule.
Covalent bonds and other bond-based topologies (angles and
dihedrals) can be automatically determined. Cannot be used for
non bond-based topologies, e.g., improper dihedrals.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Input molecule.


    * **bond** (*bool*) – Whether find bonds. If set to False, angle and
    dihedral searching will be skipped. Default to True.


    * **angle** (*bool*) – Whether find angles. Default to True.


    * **dihedral** (*bool*) – Whether find dihedrals. Default to True.


    * **tol** (*float*) – Bond distance tolerance. Default to 0.1.
    Not recommended to alter.


    * **\*\*kwargs** – Other kwargs supported by Topology.



### lattice_2_lmpbox(lattice, origin=(0, 0, 0))
Converts a lattice object to LammpsBox, and calculates the symmetry
operation used.


* **Parameters**


    * **lattice** ([*Lattice*](pymatgen.core.md#pymatgen.core.lattice.Lattice)) – Input lattice.


    * **origin** – A (3,) array/list of floats setting lower bounds of
    simulation box. Default to (0, 0, 0).



* **Returns**

    LammpsBox, SymmOp


## pymatgen.io.lammps.generators module

Input set generators for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


### _class_ BaseLammpsGenerator(template: str = <factory>, settings: dict = <factory>, calc_type: str = 'lammps', keep_stages: bool = True)
Bases: [`InputGenerator`](pymatgen.io.md#pymatgen.io.core.InputGenerator)

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

#### get_input_set(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | LammpsData | CombinedData | None)
Generate a LammpsInputSet from the structure/data, tailored to the template file.


#### keep_stages(_: boo_ _ = Tru_ )

#### settings(_: dic_ )

#### template(_: st_ )

### _class_ LammpsMinimization(template: str | None = None, units: str = 'metal', atom_style: str = 'full', dimension: int = 3, boundary: str = 'p p p', read_data: str = 'system.data', force_field: str = 'Unspecified force field!', keep_stages: bool = False)
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

## pymatgen.io.lammps.inputs module

This module implements methods for reading/manupilating/writing LAMMPS input files.
It does not implement methods for automatically creating inputs based on a structure
and computation type. For this, see the InputSet and InputGenerator in sets.py, or
[https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps).


### _class_ LammpsInputFile(stages: list | None = None)
Bases: [`InputFile`](pymatgen.io.md#pymatgen.io.core.InputFile)

Class representing a LAMMPS input settings file, e.g. in.lammps.
Allows for LAMMPS input generation in line/stage wise manner. A stage
here is defined as a block of LAMMPS input commands usually performing a
specific task during the simulation such as energy minimization or
NPT/NVT runs. But more broadly, a stage can also be a block of LAMMPS
input where the simulation box is set up, a set of variables are declared or
quantities are computed.

The LammpsInputFile is defined by the attribute stages,
i.e. a list of dicts each with keys stage_name and commands,
defining the stage names and the corresponding LAMMPS input settings (list of tuples of two strings each).
The structure is the following:


```
``
```

\`
stages = [

> {“stage_name”: “Stage 1”, “commands”: [(cmd1, args1), (cmd2, args2)]},
> {“stage_name”: “Stage 2”, “commands”: [(cmd3, args3)]}

### ]

where cmd’s are the LAMMPS command names (e.g., “units”, or “pair_coeff”)
and the args are the corresponding arguments.
“Stage 1” and “Stage 2” are examples of stage names.


* **param stages**

    list of LAMMPS input settings.



#### _add_command(stage_name: str, command: str, args: str | float | None = None)
Helper method to add a single LAMMPS command and its arguments to
the LammpsInputFile. The stage name should be provided: a default behavior
is avoided here to avoid mistakes.

### Example

In order to add the command `pair_coeff 1 1 morse 0.0580 3.987 3.404`
to the stage “Definition of the potential”, simply use


```
``
```

\`
your_input_file._add_command(

> stage_name=”Definition of the potential”,
> command=”pair_coeff 1 1 morse 0.0580 3.987 3.404”

#### )

#### or

your_input_file._add_command(

    stage_name=”Definition of the potential”,
    command=”pair_coeff”,
    args=”1 1 morse 0.0580 3.987 3.404”

#### )


* **param stage_name**

    name of the stage to which the command should be added.



* **type stage_name**

    str



* **param command**

    LAMMPS command, with or without the arguments.



* **type command**

    str



* **param args**

    Arguments for the LAMMPS command.



* **type args**

    str



#### _add_comment(comment: str, inline: bool = False, stage_name: str | None = None, index_comment: bool = False)
Method to add a comment inside a stage (between actual commands)
or as a whole stage (which will do nothing when LAMMPS runs).


* **Parameters**


    * **comment** (*str*) – Comment string to be added. The comment will be
    preceded by ‘#’ in the generated input.


    * **inline** (*bool*) – True if the comment should be inline within a given block of commands.


    * **stage_name** (*str*) – set the stage_name to which the comment will be written. Required if inline is True.


    * **index_comment** (*bool*) – True if the comment should start with “Comment x” with x its number in the ordering.
    Used only for inline comments.



#### _static_ _check_stage_format(stage: dict)

#### _static_ _clean_lines(string_list: list, ignore_comments: bool = False)
Helper method to strips whitespaces, carriage returns and redundant empty
lines from a list of strings.
Transforms “& n” and “&n” into “” as the & symbol means the line continues.
Also removes lines with “# LAMMPS input generated from LammpsInputFile”
to avoid possible repetitions.


* **Parameters**


    * **string_list** (*list*) – List of strings.


    * **ignore_comments** (*bool*) – True if the strings starting with # should be ignored.



* **Returns**

    List of strings



#### _static_ _get_blocks(string_list: list[str], keep_stages: bool = False)
Helper method to return a list of blocks of LAMMPS commands,
separated from “” in a list of all commands.


* **Parameters**


    * **string_list** (*list*) – List of strings.


    * **keep_stages** (*bool*) – True if the block structure from the input file should be kept.
    If False, a single block is assumed.



* **Returns**

    List of list of strings containing the different blocks



#### _initialize_stage(stage_name: str | None = None, stage_index: int | None = None)
Initialize an empty stage with the given name in the LammpsInputFile.


* **Parameters**


    * **stage_name** (*str*) – If a stage name is mentioned, the command is added
    under that stage block, else the new stage is named from numbering.
    If given, stage_name cannot be one of the already present stage names.


    * **stage_index** (*int*) – Index of the stage where it should be added.
    If None, the stage is added at the end of the LammpsInputFile.



#### add_commands(stage_name: str, commands: str | list[str] | dict)
Method to add a LAMMPS commands and their arguments to a stage of
the LammpsInputFile. The stage name should be provided: a default behavior
is avoided here to avoid mistakes (e.g., the commands are added to the wrong stage).

### Example

In order to add the command `pair_coeff 1 1 morse 0.0580 3.987 3.404`
to the stage “Definition of the potential”, simply use


```
``
```

\`
your_input_file.add_commands(

> stage_name=”Definition of the potential”,
> commands=”pair_coeff 1 1 morse 0.0580 3.987 3.404”

#### )

To add multiple commands, use a dict or a list, e.g.,


```
``
```

\`
your_input_file.add_commands(

> stage_name=”Definition of the potential”,
> commands=[“pair_coeff 1 1 morse 0.0580 3.987 3.404”, “units atomic”]

)
your_input_file.add_commands(

> stage_name=”Definition of the potential”,
> commands={“pair_coeff”: “1 1 morse 0.0580 3.987 3.404”, “units”: “atomic”}

#### )


* **param stage_name**

    name of the stage to which the command should be added.



* **type stage_name**

    str



* **param commands**

    LAMMPS command, with or without the arguments.



* **type commands**

    str or list or dict



#### add_stage(stage: dict | None = None, commands: str | list[str] | dict[str, str | float] | None = None, stage_name: str | None = None, after_stage: str | int | None = None)
Adds a new stage to the LammpsInputFile, either from a whole stage (dict) or
from a stage_name and commands. Both ways are mutually exclusive.

### Examples

1) In order to add a stage defining the force field to be used, you can use:


```
``
```

\`
your_input_file.add_stage(

> commands=[“pair_coeff 1 1 morse 0.0580 3.987 3.404”, “pair_coeff 1 4 morse 0.0408 1.399 3.204”],
> stage_name=”Definition of the force field”

#### )

#### or

your_input_file.add_stage(

    {

        “stage_name”: “Definition of the force field”,
        “commands”: [

        > (“pair_coeff”, “1 1 morse 0.0580 3.987 3.404”),
        > (“pair_coeff”, “1 4 morse 0.0408 1.399 3.204”)

        ],

    }

#### )

2) Another stage could consist in an energy minimization. In that case, the commands could look like


```
``
```

\`
commands = [

> “thermo 1”,
> “thermo_style custom step lx ly lz press pxx pyy pzz pe”,
> “dump dmp all atom 5 run.dump”,
> “min_style cg”,
> “fix 1 all box/relax iso 0.0 vmax 0.001”,
> “minimize 1.0e-16 1.0e-16 5000 10000”,
> “write_data run.data”

#### ]

or a dictionary such as {“thermo”: 1, …}, or a string with a single command (e.g., “units atomic”).


* **param stage**

    if provided, this is the stage that will be added to the LammpsInputFile.stages



* **type stage**

    dict



* **param commands**

    if provided, these are the LAMMPS command(s)
    that will be included in the stage to add.
    Can pass a list of LAMMPS commands with their arguments.
    Also accepts a dictionary of LAMMPS commands and
    corresponding arguments as key, value pairs.
    A single string can also be passed (single command together with its arguments).
    Not used in case a whole stage is given.



* **type commands**

    str or list or dict



* **param stage_name**

    If a stage name is mentioned, the commands are added
    under that stage block, else the new stage is named from numbering.
    If given, stage_name cannot be one of the already present stage names.
    Not used in case a whole stage is given.



* **type stage_name**

    str



* **param after_stage**

    Name of the stage after which this stage should be added.
    If None, the stage is added at the end of the LammpsInputFile.



* **type after_stage**

    str



#### append(lmp_input_file: LammpsInputFile)
Appends a LammpsInputFile to another. The stages are merged,
and the numbering of stages/comments is either kept the same or updated.


* **Parameters**

    **lmp_input_file** (*LammpsInputFile*) – LammpsInputFile to append.



#### contains_command(command: str, stage_name: str | None = None)
Returns whether a given command is present in the LammpsInputFile.
A stage name can be given; in this case the search will happen only for this stage.


* **Parameters**


    * **command** (*str*) – String with the command to find in the input file (e.g., “units”).


    * **stage_name** (*str*) – String giving the stage name where the change should take place.



* **Returns**

    True if the command is present, False if not.



* **Return type**

    bool



#### _classmethod_ from_file(path: str | Path, ignore_comments: bool = False, keep_stages: bool = False)
Creates an InputFile object from a file.


* **Parameters**


    * **path** (*str** or **path*) – Filename to read, including path.


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the input file.


    * **keep_stages** (*bool*) – True if the block structure from the input file should be kept.
    If False, a single block is assumed.



* **Returns**

    LammpsInputFile



#### _classmethod_ from_str(contents: str, ignore_comments: bool = False, keep_stages: bool = False)
Helper method to parse string representation of LammpsInputFile.
If you created the input file by hand, there is no guarantee that the representation
will be perfect as it is difficult to account for all the cosmetic changes you
could have done on your input script. Always check that you have what you want !
By default, a single stage containing all the input settings is created.
If the block structure of your input file should be kept and stored as
different stages, set keep_stages to True.


* **Parameters**


    * **contents** (*str*) – String representation of LammpsInputFile.


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the input file.


    * **keep_stages** (*bool*) – True if the block structure from the input file should be kept.
    If False, a single block is assumed.



* **Returns**

    LammpsInputFile



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_args(command: str, stage_name: str | None = None)
Given a command, returns the corresponding arguments (or list of arguments) in the LammpsInputFile.
A stage name can be given; in this case the search will happen only for this stage.
If a stage name is not given, the command will be searched for through all of them.
If the command is not found, an empty list is returned.


* **Parameters**


    * **command** (*str*) – String with the command to find in the input file (e.g., “units”).


    * **stage_name** (*str*) – String giving the stage name where the change should take place.



* **Returns**

    Value of the argument corresponding to the command.
    List if the same command is used multiple times.



#### get_str(ignore_comments: bool = False, keep_stages: bool = True)
Generates and ² the string representation of the LammpsInputFile.
Stages are separated by empty lines.
The headers of the stages will be put in comments preceding each stage.
Other comments will be put inline within stages, where they have been added.


* **Parameters**


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the InputFile.


    * **keep_stages** (*bool*) – If True, the string is formatted in a block structure with stage names
    and newlines that differentiate commands in the respective stages of the InputFile.
    If False, stage names are not printed and all commands appear in a single block.



* **Returns**

    String representation of the LammpsInputFile.



* **Return type**

    str



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### merge_stages(stage_names: list[str])
Merges multiple stages of a LammpsInputFile together.
The merged stage will be at the same index as the first of the stages to be merged.
The others will appear in the same order as provided in the list. Other non-merged stages will follow.


* **Parameters**

    **stage_names** (*list*) – list of strings giving the names of the stages to be merged.



#### _property_ ncomments(_: in_ )
Returns the number of comments in the current LammpsInputFile. Includes the blocks of comments as well
as inline comments (comment lines within blocks of LAMMPS commands).


#### _property_ nstages(_: in_ )
Returns the number of stages in the current LammpsInputFile.


#### remove_command(command: str, stage_name: str | list[str] | None = None, remove_empty_stages: bool = True)
Removes a given command from a given stage. If no stage is given, removes all occurrences of the command.
In case removing a command completely empties a stage, the choice whether to keep this stage in the
LammpsInputFile is given by remove_empty_stages.


* **Parameters**


    * **command** (*str*) – command to be removed.


    * **stage_name** (*str** or **list*) – names of the stages where the command should be removed.


    * **remove_empty_stages** (*bool*) – whether to remove the stages emptied by removing the command or not.



#### remove_stage(stage_name: str)
Removes a whole stage from the LammpsInputFile.


* **Parameters**

    **stage_name** (*str*) – name of the stage to remove.



#### rename_stage(stage_name: str, new_name: str)
Renames a stage stage_name from LammpsInputFile into new_name.
First checks that the stage to rename is present, and that
the new name is not already a stage name.


* **Parameters**


    * **stage_name** (*str*) – name of the stage to rename.


    * **new_name** (*str*) – new name of the stage.



#### set_args(command: str, argument: str, stage_name: str | None = None, how: str | int | list[int] = 'all')
Sets the arguments for the given command to the given string.
If the command is not found, nothing is done. Use LammpsInputFile.add_commands instead.
If a stage name is specified, it will be replaced or set only for this stage.
If no stage name is given, it will apply the change in all of them that contain the given command.
If the command is set multiple times in the file/stage, it will be replaced based on “how”:
either the first occurrence, all of them, or the index of the occurrence.


* **Parameters**


    * **command** (*str*) – String representing the command to change, e.g., “units”.


    * **argument** (*str*) – String with the new value for the command, e.g., “atomic”.


    * **stage_name** (*str*) – String giving the stage name where the change should take place.


    * **how** (*str** or **int** or **list*) – “all” for changing all occurrences of the command within the stage_name
    or the whole InputFile, “first” for the first occurrence, int i for the i-th time the command
    is present in the stage_name or the whole InputFile, starting at 0. Can be a list of indexes as well.



#### _property_ stages_names(_: lis_ )
List of names for all the stages present in stages.


#### write_file(filename: str | Path, ignore_comments: bool = False, keep_stages: bool = True)
Writes the input file.


* **Parameters**


    * **filename** (*str** or **path*) – The filename to output to, including path.


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the InputFile.


    * **keep_stages** (*bool*) – True if the block structure from the InputFile should be kept.
    If False, a single block is assumed.



### _class_ LammpsRun(script_template, settings, data, script_filename)
Bases: `MSONable`

Examples for various simple LAMMPS runs with given simulation box,
force field and a few more settings. Experienced LAMMPS users should
consider using write_lammps_inputs method with more sophisticated
templates.

Base constructor.


* **Parameters**


    * **script_template** (*str*) – String template for input script
    with placeholders. The format for placeholders has to
    be ‘$variable_name’, e.g., ‘$temperature’


    * **settings** (*dict*) – Contains values to be written to the
    placeholders, e.g., {‘temperature’: 1}.


    * **data** (*LammpsData** or **str*) – Data file as a LammpsData
    instance or path to an existing data file. Default to
    None, i.e., no data file supplied. Useful only when
    read_data cmd is in the script.


    * **script_filename** (*str*) – Filename for the input script.



#### _classmethod_ md(data, force_field, temperature, nsteps, other_settings=None)
Example for a simple MD run based on template md.template.


* **Parameters**


    * **data** (*LammpsData** or **str*) – Data file as a LammpsData
    instance or path to an existing data file.


    * **force_field** (*str*) – Combined force field related cmds. For
    example, ‘pair_style eamnpair_coeff \* \* Cu_u3.eam’.


    * **temperature** (*float*) – Simulation temperature.


    * **nsteps** (*int*) – No. of steps to run.


    * **other_settings** (*dict*) – other settings to be filled into
    placeholders.



#### template_dir(_ = '/Users/shyue/repos/pymatgen/pymatgen/io/lammps/templates_ )

#### write_inputs(output_dir, \*\*kwargs)
Writes all input files (input script, and data if needed).
Other supporting files are not handled at this moment.


* **Parameters**


    * **output_dir** (*str*) – Directory to output the input files.


    * **\*\*kwargs** – kwargs supported by LammpsData.write_file.



### _class_ LammpsTemplateGen()
Bases: [`TemplateInputGen`](pymatgen.io.md#pymatgen.io.template.TemplateInputGen)

Creates an InputSet object for a LAMMPS run based on a template file.
The input script is constructed by substituting variables into placeholders
in the template file using python’s Template.safe_substitute() function.
The data file containing coordinates and topology information can be provided
as a LammpsData instance. Alternatively, you can include a read_data command
in the template file that points to an existing data file.
Other supporting files are not handled at the moment.

To write the input files to a directory, call LammpsTemplateSet.write_input()
See pymatgen.io.template.py for additional documentation of this method.


#### get_input_set(script_template: str | Path, settings: dict | None = None, script_filename: str = 'in.lammps', data: LammpsData | CombinedData | None = None, data_filename: str = 'system.data')

* **Parameters**


    * **script_template** – String template for input script with
    placeholders. The format for placeholders has to be
    ‘$variable_name’, e.g., ‘$temperature’


    * **settings** – Contains values to be written to the
    placeholders, e.g., {‘temperature’: 1}. Default to None.


    * **data** – Data file as a LammpsData instance. Default to None, i.e., no
    data file supplied. Note that a matching ‘read_data’ command
    must be provided in the script template in order for the data
    file to actually be read.


    * **script_filename** – Filename for the input file.


    * **data_filename** – Filename for the data file, if provided.


## pymatgen.io.lammps.outputs module

This module implements classes and methods for processing LAMMPS output
files (log and dump).


### _class_ LammpsDump(timestep, natoms, box, data)
Bases: `MSONable`

Object for representing dump data for a single snapshot.

Base constructor.


* **Parameters**


    * **timestep** (*int*) – Current time step.


    * **natoms** (*int*) – Total number of atoms in the box.


    * **box** (*LammpsBox*) – Simulation box.


    * **data** (*pd.DataFrame*) – Dumped atomic data.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    LammpsDump



#### _classmethod_ from_str(string)
Constructor from string parsing.


* **Parameters**

    **string** (*str*) – Input string.



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### parse_lammps_dumps(file_pattern)
Generator that parses dump file(s).


* **Parameters**

    **file_pattern** (*str*) – Filename to parse. The timestep wildcard
    (e.g., dump.atom.’\*’) is supported and the files are parsed
    in the sequence of timestep.



* **Yields**

    LammpsDump for each available snapshot.



### parse_lammps_log(filename='log.lammps')
Parses log file with focus on thermo data. Both one and multi line
formats are supported. Any incomplete runs (no “Loop time” marker)
will not be parsed.

### Notes

SHAKE stats printed with thermo data are not supported yet.
They are ignored in multi line format, while they may cause
issues with dataframe parsing in one line format.


* **Parameters**

    **filename** (*str*) – Filename to parse.



* **Returns**

    [pd.DataFrame] containing thermo data for each completed run.


## pymatgen.io.lammps.sets module

Input sets for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
([https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps)).


### _class_ LammpsInputSet(inputfile: LammpsInputFile | str, data: LammpsData | CombinedData, calc_type: str = '', template_file: str = '', keep_stages: bool = False)
Bases: [`InputSet`](pymatgen.io.md#pymatgen.io.core.InputSet)

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



#### _abc_impl(_ = <_abc._abc_data object_ )

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

## pymatgen.io.lammps.utils module

This module defines utility classes and functions.


### _class_ LammpsRunner(input_filename='lammps.in', bin='lammps')
Bases: `object`

LAMMPS wrapper.


* **Parameters**


    * **input_filename** (*str*) – input file name


    * **bin** (*str*) – command to run, excluding the input file name.



#### run()
Write the input/data files and run LAMMPS.


### _class_ Polymer(start_monomer, s_head, s_tail, monomer, head, tail, end_monomer, e_head, e_tail, n_units, link_distance=1.0, linear_chain=False)
Bases: `object`

Generate polymer chain via Random walk. At each position there are
a total of 5 possible moves(excluding the previous direction).


* **Parameters**


    * **start_monomer** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Starting molecule


    * **s_head** (*int*) – starting atom index of the start_monomer molecule


    * **s_tail** (*int*) – tail atom index of the start_monomer


    * **monomer** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – The monomer


    * **head** (*int*) – index of the atom in the monomer that forms the head


    * **tail** (*int*) – tail atom index. monomers will be connected from
    tail to head


    * **end_monomer** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Terminal molecule


    * **e_head** (*int*) – starting atom index of the end_monomer molecule


    * **e_tail** (*int*) – tail atom index of the end_monomer


    * **n_units** (*int*) – number of monomer units excluding the start and
    terminal molecules


    * **link_distance** (*float*) – distance between consecutive monomers


    * **linear_chain** (*bool*) – linear or random walk polymer chain.



#### _add_monomer(monomer, mon_vector, move_direction)
extend the polymer molecule by adding a monomer along mon_vector direction.


* **Parameters**


    * **monomer** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – monomer molecule


    * **mon_vector** (*numpy.array*) – monomer vector that points from head to tail.


    * **move_direction** (*numpy.array*) – direction along which the monomer
    will be positioned



#### _align_monomer(monomer, mon_vector, move_direction)
rotate the monomer so that it is aligned along the move direction.


* **Parameters**


    * **monomer** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) –


    * **mon_vector** (*numpy.array*) – molecule vector that starts from the
    start atom index to the end atom index


    * **move_direction** (*numpy.array*) – the direction of the polymer chain
    extension



#### _create(monomer, mon_vector)
create the polymer from the monomer.


* **Parameters**


    * **monomer** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) –


    * **mon_vector** (*numpy.array*) – molecule vector that starts from the
    start atom index to the end atom index



#### _next_move_direction()
Pick a move at random from the list of moves.