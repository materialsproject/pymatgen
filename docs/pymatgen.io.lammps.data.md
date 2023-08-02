---
layout: default
title: pymatgen.io.lammps.data.md
nav_exclude: true
---

# pymatgen.io.lammps.data module

This module implements a core class LammpsData for generating/parsing
LAMMPS data file, and other bridging classes to build LammpsData from
molecules. This module also implements a subclass CombinedData for
merging LammpsData object.

Only point particle styles are supported for now (atom_style in angle,
atomic, bond, charge, full and molecular only). See the pages below for
more info.

> [https://docs.lammps.org/atom_style.html](https://docs.lammps.org/atom_style.html)
> [https://docs.lammps.org/read_data.html](https://docs.lammps.org/read_data.html)


### _class_ pymatgen.io.lammps.data.CombinedData(list_of_molecules, list_of_names, list_of_numbers, coordinates, atom_style='full')
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
Convert a CombinedData object to a LammpsData object. attributes are deepcopied.

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


#### get_string(distance=6, velocity=8, charge=4, hybrid=True)
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



#### _classmethod_ parse_xyz(filename)
Load xyz file generated from packmol (for those who find it hard to install openbabel).


* **Returns**

    pandas.DataFrame



#### _property_ structure()
Exports a periodic structure object representing the simulation
box.


* **Returns**

    Structure



### _class_ pymatgen.io.lammps.data.ForceField(mass_info, nonbond_coeffs=None, topo_coeffs=None)
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



### _class_ pymatgen.io.lammps.data.LammpsBox(bounds, tilt=None)
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



#### get_string(significant_figures=6)
Returns the string representation of simulation box in LAMMPS
data file format.


* **Parameters**

    **significant_figures** (*int*) – No. of significant figures to
    output for box settings. Default to 6.



* **Returns**

    String representation



#### to_lattice()
Converts the simulation box to a more powerful Lattice backend.
Note that Lattice is always periodic in 3D space while a
simulation box is not necessarily periodic in all dimensions.


* **Returns**

    Lattice



#### _property_ volume()
Volume of simulation box.


### _class_ pymatgen.io.lammps.data.LammpsData(box, masses, atoms, velocities=None, force_field=None, topology=None, atom_style='full')
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



#### _classmethod_ from_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), ff_elements=None, atom_style='charge', is_sort=False)
Simple constructor building LammpsData from a structure without
force field parameters and topologies.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure.


    * **ff_elements** (*[**str**]*) – List of strings of elements that must
    be present due to force field settings but not
    necessarily in the structure. Default to None.


    * **atom_style** (*str*) – Choose between “atomic” (neutral) and
    “charge” (charged). Default to “charge”.


    * **is_sort** (*bool*) – whether to sort sites



#### get_string(distance=6, velocity=8, charge=4, hybrid=True)
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



### _class_ pymatgen.io.lammps.data.Topology(sites, ff_label=None, charges=None, velocities=None, topologies=None)
Bases: `MSONable`

Class carrying most data in Atoms, Velocities and molecular
topology sections for ONE SINGLE Molecule or Structure
object, or a plain list of Sites.


* **Parameters**


    * **sites** (*[*[*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)*] or *[*SiteCollection*](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection)) – A group of sites in a
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


    * **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Input molecule.


    * **bond** (*bool*) – Whether find bonds. If set to False, angle and
    dihedral searching will be skipped. Default to True.


    * **angle** (*bool*) – Whether find angles. Default to True.


    * **dihedral** (*bool*) – Whether find dihedrals. Default to True.


    * **tol** (*float*) – Bond distance tolerance. Default to 0.1.
    Not recommended to alter.


    * **\*\*kwargs** – Other kwargs supported by Topology.



### pymatgen.io.lammps.data.lattice_2_lmpbox(lattice, origin=(0, 0, 0))
Converts a lattice object to LammpsBox, and calculates the symmetry
operation used.


* **Parameters**


    * **lattice** ([*Lattice*](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice)) – Input lattice.


    * **origin** – A (3,) array/list of floats setting lower bounds of
    simulation box. Default to (0, 0, 0).



* **Returns**

    LammpsBox, SymmOp