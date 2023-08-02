---
layout: default
title: pymatgen.core.structure.md
nav_exclude: true
---

# pymatgen.core.structure module

This module provides classes used to define a non-periodic molecule and a
periodic structure.


### _class_ pymatgen.core.structure.IMolecule(species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float = 0.0, spin_multiplicity: int | None = None, validate_proximity: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, charge_spin_check: bool = True)
Bases: `SiteCollection`, `MSONable`

Basic immutable Molecule object without periodicity. Essentially a
sequence of sites. IMolecule is made to be immutable so that they can
function as keys in a dict. For a mutable molecule,
use the :class:Molecule.

Molecule extends Sequence and Hashable, which means that in many cases,
it can be used like any Python sequence. Iterating through a molecule is
equivalent to going through the sites in sequence.

Create a Molecule.


* **Parameters**


    * **species** – list of atomic species. Possible kinds of input include a
    list of dict of elements/species and occupancies, a List of
    elements/specie specified as actual Element/Species, Strings
    (“Fe”, “Fe2+”) or atomic numbers (1,56).


    * **coords** (*3x1 array*) – list of Cartesian coordinates of each species.


    * **charge** (*float*) – Charge for the molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multiplicity for molecule.
    Defaults to None, which means that the spin multiplicity is
    set to 1 if the molecule has no unpaired electrons and to 2
    if there are unpaired electrons.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 1 Ang apart. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as
    a dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The
    sequences have to be the same length as the atomic species
    and fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **charge_spin_check** (*bool*) – Whether to check that the charge and
    spin multiplicity are compatible with each other. Defaults
    to True.



#### as_dict()
JSON-serializable dict representation of Molecule.


#### break_bond(ind1: int, ind2: int, tol: float = 0.2)
Returns two molecules based on breaking the bond between atoms at index
ind1 and ind2.


* **Parameters**


    * **ind1** (*int*) – 1st site index


    * **ind2** (*int*) – 2nd site index


    * **tol** (*float*) – Relative tolerance to test. Basically, the code
    checks if the distance between the sites is less than (1 +
    tol) \* typical bond distances. Defaults to 0.2, i.e.,
    20% longer.



* **Returns**

    Two Molecule objects representing the two clusters formed from
    breaking the bond.



#### _property_ center_of_mass(_: ndarra_ )
Center of mass of molecule.


#### _property_ charge(_: floa_ )
Charge of molecule.


#### _classmethod_ from_dict(d)
Reconstitute a Molecule object from a dict representation created using
as_dict().


* **Parameters**

    **d** (*dict*) – dict representation of Molecule.



* **Returns**

    Molecule object



#### _classmethod_ from_file(filename)
Reads a molecule from a file. Supported formats include xyz,
gaussian input (gjf|g03|g09|com|inp), Gaussian output (.out|and
pymatgen’s JSON-serialized molecules. Using openbabel,
many more extensions are supported but requires openbabel to be
installed.


* **Parameters**

    **filename** (*str*) – The filename to read from.



* **Returns**

    Molecule



#### _classmethod_ from_sites(sites: Sequence[[Site](pymatgen.core.sites.md#pymatgen.core.sites.Site)], charge: float = 0, spin_multiplicity: int | None = None, validate_proximity: bool = False, charge_spin_check: bool = True)
Convenience constructor to make a Molecule from a list of sites.


* **Parameters**


    * **sites** (*[*[*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)*]*) – Sequence of Sites.


    * **charge** (*int*) – Charge of molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multicipity. Defaults to None,
    in which it is determined automatically.


    * **validate_proximity** (*bool*) – Whether to check that atoms are too
    close.


    * **charge_spin_check** (*bool*) – Whether to check that the charge and
    spin multiplicity are compatible with each other. Defaults
    to True.



#### _classmethod_ from_str(input_string: str, fmt: Literal['xyz', 'gjf', 'g03', 'g09', 'com', 'inp', 'json', 'yaml'])
Reads the molecule from a string.


* **Parameters**


    * **input_string** (*str*) – String to parse.


    * **fmt** (*str*) – Format to output to. Defaults to JSON unless filename
    is provided. If fmt is specifies, it overrides whatever the
    filename is. Options include “xyz”, “gjf”, “g03”, “json”. If
    you have OpenBabel installed, any of the formats supported by
    OpenBabel. Non-case sensitive.



* **Returns**

    IMolecule or Molecule.



#### get_boxed_structure(a: float, b: float, c: float, images: ArrayLike = (1, 1, 1), random_rotation: bool = False, min_dist: float = 1.0, cls=None, offset: ArrayLike | None = None, no_cross: bool = False, reorder: bool = True)
Creates a Structure from a Molecule by putting the Molecule in the
center of a orthorhombic box. Useful for creating Structure for
calculating molecules using periodic codes.


* **Parameters**


    * **a** (*float*) – a-lattice parameter.


    * **b** (*float*) – b-lattice parameter.


    * **c** (*float*) – c-lattice parameter.


    * **images** – No. of boxed images in each direction. Defaults to
    (1, 1, 1), meaning single molecule with 1 lattice parameter
    in each direction.


    * **random_rotation** (*bool*) – Whether to apply a random rotation to
    each molecule. This jumbles all the molecules so that they
    are not exact images of each other.


    * **min_dist** (*float*) – The minimum distance that atoms should be from
    each other. This is only used if random_rotation is True.
    The randomized rotations are searched such that no two atoms
    are less than min_dist from each other.


    * **cls** – The Structure class to instantiate (defaults to pymatgen
    structure)


    * **offset** – Translation to offset molecule from center of mass coords


    * **no_cross** – Whether to forbid molecule coords from extending beyond
    boundary of box.


    * **reorder** – Whether to reorder the sites to be in electronegativity
    order.



* **Returns**

    Structure containing molecule in a box.



#### get_centered_molecule()
Returns a Molecule centered at the center of mass.


* **Returns**

    Molecule centered with center of mass at origin.



#### get_covalent_bonds(tol: float = 0.2)
Determines the covalent bonds in a molecule.


* **Parameters**

    **tol** (*float*) – The tol to determine bonds in a structure. See
    CovalentBond.is_bonded.



* **Returns**

    List of bonds



#### get_distance(i: int, j: int)
Get distance between site i and j.


* **Parameters**


    * **i** (*int*) – 1st site index


    * **j** (*int*) – 2nd site index



* **Returns**

    Distance between the two sites.



#### get_neighbors(site: [Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), r: float)
Get all neighbors to a site within a sphere of radius r. Excludes the
site itself.


* **Parameters**


    * **site** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – Site at the center of the sphere.


    * **r** (*float*) – Radius of sphere.



* **Returns**

    Neighbor



#### get_neighbors_in_shell(origin: ArrayLike, r: float, dr: float)
Returns all sites in a shell centered on origin (coords) between radii
r-dr and r+dr.


* **Parameters**


    * **origin** (*3x1 array*) – Cartesian coordinates of center of sphere.


    * **r** (*float*) – Inner radius of shell.


    * **dr** (*float*) – Width of shell.



* **Returns**

    Neighbor



#### get_sites_in_sphere(pt: ArrayLike, r: float)
Find all sites within a sphere from a point.


* **Parameters**


    * **pt** (*3x1 array*) – Cartesian coordinates of center of sphere


    * **r** (*float*) – Radius of sphere.



* **Returns**

    Neighbor



#### get_zmatrix()
Returns a z-matrix representation of the molecule.


#### _property_ nelectrons(_: floa_ )
Number of electrons in the molecule.


#### _property_ spin_multiplicity(_: floa_ )
Spin multiplicity of molecule.


#### to(filename: str = '', fmt: str = '')
Outputs the molecule to a file or string.


* **Parameters**


    * **filename** (*str*) – If provided, output will be written to a file. If
    fmt is not specified, the format is determined from the
    filename. Defaults is None, i.e. string output.


    * **fmt** (*str*) – Format to output to. Defaults to JSON unless filename
    is provided. If fmt is specifies, it overrides whatever the
    filename is. Options include “xyz”, “gjf”, “g03”, “json”. If
    you have OpenBabel installed, any of the formats supported by
    OpenBabel. Non-case sensitive.



* **Returns**

    (str) if filename is None. None otherwise.



### _class_ pymatgen.core.structure.IStructure(lattice: ArrayLike | [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False, coords_are_cartesian: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None)
Bases: `SiteCollection`, `MSONable`

Basic immutable Structure object with periodicity. Essentially a sequence
of PeriodicSites having a common lattice. IStructure is made to be
(somewhat) immutable so that they can function as keys in a dict. To make
modifications, use the standard Structure object instead. Structure
extends Sequence and Hashable, which means that in many cases,
it can be used like any Python sequence. Iterating through a
structure is equivalent to going through the sites in sequence.

Create a periodic structure.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    [`pymatgen.core.lattice.Lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice) or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** (*[*[*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)*]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/Cartesian coordinates of
    each species.


    * **charge** (*int*) – overall charge of the structure. Defaults to behavior
    in SiteCollection where total charge is the sum of the oxidation
    states.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to map all sites into the unit cell,
    i.e. fractional coords between 0 and 1. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g. {“magmom”:[5, 5, 5, 5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.



#### as_dataframe()
Create a Pandas dataframe of the sites. Structure-level attributes are stored in DataFrame.attrs.

### Example

Species    a    b             c    x             y             z  magmom
0    (Si)  0.0  0.0  0.000000e+00  0.0  0.000000e+00  0.000000e+00       5
1    (Si)  0.0  0.0  1.000000e-07  0.0 -2.217138e-07  3.135509e-07      -5


#### as_dict(verbosity=1, fmt=None, \*\*kwargs)
Dict representation of Structure.


* **Parameters**


    * **verbosity** (*int*) – Verbosity level. Default of 1 includes both
    direct and Cartesian coordinates for all sites, lattice
    parameters, etc. Useful for reading and for insertion into a
    database. Set to 0 for an extremely lightweight version
    that only includes sufficient information to reconstruct the
    object.


    * **fmt** (*str*) – Specifies a format for the dict. Defaults to None,
    which is the default format used in pymatgen. Other options
    include “abivars”.


    * **\*\*kwargs** – Allow passing of other kwargs needed for certain


    * **formats** –


    * **e.g.** –


    * **"abivars".** –



* **Returns**

    JSON-serializable dict representation.



#### _property_ charge(_: floa_ )
Overall charge of the structure.


#### copy(site_properties=None, sanitize=False)
Convenience method to get a copy of the structure, with options to add
site properties.


* **Parameters**


    * **site_properties** (*dict*) – Properties to add or override. The
    properties are specified in the same way as the constructor,
    i.e., as a dict of the form {property: [values]}. The
    properties should be in the order of the *original* structure
    if you are performing sanitization.


    * **sanitize** (*bool*) – If True, this method will return a sanitized
    structure. Sanitization performs a few things: (i) The sites are
    sorted by electronegativity, (ii) a LLL lattice reduction is
    carried out to obtain a relatively orthogonalized cell,
    (iii) all fractional coords for sites are mapped into the
    unit cell.



* **Returns**

    A copy of the Structure, with optionally new site_properties and
    optionally sanitized.



#### _property_ density(_: floa_ )
Returns the density in units of g/cm^3.


#### _property_ distance_matrix(_: ndarra_ )
Returns the distance matrix between all sites in the structure. For
periodic structures, this should return the nearest image distance.


#### _property_ frac_coords()
Fractional coordinates as a Nx3 numpy array.


#### _classmethod_ from_dict(d: dict[str, Any], fmt: Literal['abivars'] | None = None)
Reconstitute a Structure object from a dict representation of Structure
created using as_dict().


* **Parameters**


    * **d** (*dict*) – Dict representation of structure.


    * **fmt** (*'abivars'** | **None*) – Use structure_from_abivars() to parse the dict. Defaults to None.



* **Returns**

    Structure object



#### _classmethod_ from_file(filename, primitive=False, sort=False, merge_tol=0.0, \*\*kwargs)
Reads a structure from a file. For example, anything ending in
a “cif” is assumed to be a Crystallographic Information Format file.
Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT,
vasprun.xml, CSSR, Netcdf and pymatgen’s JSON-serialized structures.


* **Parameters**


    * **filename** (*str*) – The filename to read from.


    * **primitive** (*bool*) – Whether to convert to a primitive cell. Only available for cifs. Defaults to False.


    * **sort** (*bool*) – Whether to sort sites. Default to False.


    * **merge_tol** (*float*) – If this is some positive number, sites that are within merge_tol from each other will be
    merged. Usually 0.01 should be enough to deal with common numerical issues.


    * **kwargs** – Passthrough to relevant reader. E.g. if the file has CIF format, the kwargs will be passed
    through to CifParser.



* **Returns**

    Structure.



#### _classmethod_ from_magnetic_spacegroup(msg: str | [MagneticSpaceGroup](pymatgen.symmetry.maggroups.md#pymatgen.symmetry.maggroups.MagneticSpaceGroup), lattice: list | np.ndarray | [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), species: Sequence[str | [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies) | [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition)], coords: Sequence[Sequence[float]], site_properties: dict[str, Sequence], coords_are_cartesian: bool = False, tol: float = 1e-05, labels: Sequence[str | None] | None = None)
Generate a structure using a magnetic spacegroup. Note that only
symmetrically distinct species, coords and magmoms should be provided.]
All equivalent sites are generated from the spacegroup operations.


* **Parameters**


    * **msg** (str/list/[`pymatgen.symmetry.maggroups.MagneticSpaceGroup`](pymatgen.symmetry.maggroups.md#pymatgen.symmetry.maggroups.MagneticSpaceGroup)) – The magnetic spacegroup.
    If a string, it will be interpreted as one of the notations
    supported by MagneticSymmetryGroup, e.g., “R-3’c” or “Fm’-3’m”.
    If a list of two ints, it will be interpreted as the number of
    the spacegroup in its Belov, Neronova and Smirnova (BNS) setting.


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    [`pymatgen.core.lattice.Lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice) or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
    Note that no attempt is made to check that the lattice is
    compatible with the spacegroup specified. This may be
    introduced in a future version.


    * **species** (*[*[*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)*]*) – Sequence of species on each site. Can take in
    flexible input, including:
    i.  A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        1. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Unlike Structure.from_spacegroup(),
    this argument is mandatory, since magnetic moment information
    has to be included. Note that the *direction* of the supplied
    magnetic moment relative to the crystal is important, even if
    the resulting structure is used for collinear calculations.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **tol** (*float*) – A fractional tolerance to deal with numerical
    precision issues in determining if orbits are the same.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.



* **Returns**

    Structure | IStructure



#### _classmethod_ from_sites(sites: list[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False)
Convenience constructor to make a Structure from a list of sites.


* **Parameters**


    * **sites** – Sequence of PeriodicSites. Sites must have the same
    lattice.


    * **charge** – Charge of structure.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to translate sites into the unit
    cell.



* **Returns**

    (Structure) Note that missing properties are set as None.



#### _classmethod_ from_spacegroup(sg: str | int, lattice: list | np.ndarray | [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), species: Sequence[str | [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies) | [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition)], coords: Sequence[Sequence[float]], site_properties: dict[str, Sequence] | None = None, coords_are_cartesian: bool = False, tol: float = 1e-05, labels: Sequence[str | None] | None = None)
Generate a structure using a spacegroup. Note that only symmetrically
distinct species and coords should be provided. All equivalent sites
are generated from the spacegroup operations.


* **Parameters**


    * **sg** (*str/int*) – The spacegroup. If a string, it will be interpreted
    as one of the notations supported by
    pymatgen.symmetry.groups.Spacegroup. E.g., “R-3c” or “Fm-3m”.
    If an int, it will be interpreted as an international number.


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    [`pymatgen.core.lattice.Lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice) or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
    Note that no attempt is made to check that the lattice is
    compatible with the spacegroup specified. This may be
    introduced in a future version.


    * **species** (*[*[*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)*]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **tol** (*float*) – A fractional tolerance to deal with numerical
    precision issues in determining if orbits are the same.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.



#### _classmethod_ from_str(input_string: str, fmt: Literal['cif', 'poscar', 'cssr', 'json', 'yaml', 'xsf', 'mcsqs', 'res'], primitive: bool = False, sort: bool = False, merge_tol: float = 0.0, \*\*kwargs)
Reads a structure from a string.


* **Parameters**


    * **input_string** (*str*) – String to parse.


    * **fmt** (*str*) – A file format specification. One of “cif”, “poscar”, “cssr”,
    “json”, “yaml”, “xsf”, “mcsqs”.


    * **primitive** (*bool*) – Whether to find a primitive cell. Defaults to
    False.


    * **sort** (*bool*) – Whether to sort the sites in accordance to the default
    ordering criteria, i.e., electronegativity.


    * **merge_tol** (*float*) – If this is some positive number, sites that
    are within merge_tol from each other will be merged. Usually
    0.01 should be enough to deal with common numerical issues.


    * **\*\*kwargs** – Passthrough to relevant parser.



* **Returns**

    IStructure | Structure



#### get_all_neighbors(r: float, include_index: bool = False, include_image: bool = False, sites: Sequence[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, numerical_tol: float = 1e-08)
Get neighbors for each atom in the unit cell, out to a distance r
Returns a list of list of neighbors for each site in structure.
Use this method if you are planning on looping over all sites in the
crystal. If you only want neighbors for a particular site, use the
method get_neighbors as it may not have to build such a large supercell
However if you are looping over all sites in the crystal, this method
is more efficient since it only performs one pass over a large enough
supercell to contain all possible atoms out to a distance r.
The return type is a [(site, dist) …] since most of the time,
subsequent processing requires the distance.

A note about periodic images: Before computing the neighbors, this
operation translates all atoms to within the unit cell (having
fractional coordinates within [0,1)). This means that the “image” of a
site does not correspond to how much it has been translates from its
current position, but which image of the unit cell it resides.


* **Parameters**


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.



* **Returns**

    [[`pymatgen.core.structure.PeriodicNeighbor`], ..]



#### get_all_neighbors_old(\*\*kwargs)

#### get_all_neighbors_py(r: float, include_index: bool = False, include_image: bool = False, sites: Sequence[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, numerical_tol: float = 1e-08)
Get neighbors for each atom in the unit cell, out to a distance r
Returns a list of list of neighbors for each site in structure.
Use this method if you are planning on looping over all sites in the
crystal. If you only want neighbors for a particular site, use the
method get_neighbors as it may not have to build such a large supercell
However if you are looping over all sites in the crystal, this method
is more efficient since it only performs one pass over a large enough
supercell to contain all possible atoms out to a distance r.
The return type is a [(site, dist) …] since most of the time,
subsequent processing requires the distance.

A note about periodic images: Before computing the neighbors, this
operation translates all atoms to within the unit cell (having
fractional coordinates within [0,1)). This means that the “image” of a
site does not correspond to how much it has been translates from its
current position, but which image of the unit cell it resides.


* **Parameters**


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.



* **Returns**

    list[list[PeriodicNeighbor]]



#### get_distance(i: int, j: int, jimage=None)
Get distance between site i and j assuming periodic boundary
conditions. If the index jimage of two sites atom j is not specified it
selects the jimage nearest to the i atom and returns the distance and
jimage indices in terms of lattice vector translations if the index
jimage of atom j is specified it returns the distance between the i
atom and the specified jimage atom.


* **Parameters**


    * **i** (*int*) – 1st site index


    * **j** (*int*) – 2nd site index


    * **jimage** – Number of lattice translations in each lattice direction.
    Default is None for nearest image.



* **Returns**

    distance



#### get_miller_index_from_site_indexes(site_ids, round_dp=4, verbose=True)
Get the Miller index of a plane from a set of sites indexes.

A minimum of 3 sites are required. If more than 3 sites are given
the best plane that minimises the distance to all points will be
calculated.


* **Parameters**


    * **site_ids** (*list** of **int*) – A list of site indexes to consider. A
    minimum of three site indexes are required. If more than three
    sites are provided, the best plane that minimises the distance
    to all sites will be calculated.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)



#### get_neighbor_list(r: float, sites: Sequence[[PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)] | None = None, numerical_tol: float = 1e-08, exclude_self: bool = True)
Get neighbor lists using numpy array representations without constructing
Neighbor objects. If the cython extension is installed, this method will
be orders of magnitude faster than get_all_neighbors_old and 2-3x faster
than get_all_neighbors.
The returned values are a tuple of numpy arrays
(center_indices, points_indices, offset_vectors, distances).
Atom center_indices[i] has neighbor atom points_indices[i] that is
translated by offset_vectors[i] lattice vectors, and the distance is
distances[i].


* **Parameters**


    * **r** (*float*) – Radius of sphere


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.


    * **exclude_self** (*bool*) – whether to exclude atom neighboring with itself within
    numerical tolerance distance, default to True


Returns: (center_indices, points_indices, offset_vectors, distances)


#### get_neighbors(site: [PeriodicSite](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite), r: float, include_index: bool = False, include_image: bool = False)
Get all neighbors to a site within a sphere of radius r. Excludes the
site itself.


* **Parameters**


    * **site** ([*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)) – Which is the center of the sphere.


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.



* **Returns**

    PeriodicNeighbor



#### get_neighbors_in_shell(origin: ArrayLike, r: float, dr: float, include_index: bool = False, include_image: bool = False)
Returns all sites in a shell centered on origin (coords) between radii
r-dr and r+dr.


* **Parameters**


    * **origin** (*3x1 array*) – Cartesian coordinates of center of sphere.


    * **r** (*float*) – Inner radius of shell.


    * **dr** (*float*) – Width of shell.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.



* **Returns**

    [NearestNeighbor] where Nearest Neighbor is a named tuple containing
    (site, distance, index, image).



#### get_neighbors_old(\*\*kwargs)

#### get_orderings(mode: Literal['enum', 'sqs'] = 'enum', \*\*kwargs)
Returns list of orderings for a disordered structure. If structure
does not contain disorder, the default structure is returned.


* **Parameters**


    * **mode** (*"enum"** | **"sqs"*) – Either “enum” or “sqs”. If enum,
    the enumlib will be used to return all distinct
    orderings. If sqs, mcsqs will be used to return
    an sqs structure.


    * **kwargs** – kwargs passed to either
    pymatgen.command_line..enumlib_caller.EnumlibAdaptor
    or pymatgen.command_line.mcsqs_caller.run_mcsqs.
    For run_mcsqs, a default cluster search of 2 cluster interactions
    with 1NN distance and 3 cluster interactions with 2NN distance
    is set.



* **Returns**

    List[Structure]



#### get_primitive_structure(tolerance: float = 0.25, use_site_props: bool = False, constrain_latt: list | dict | None = None)
This finds a smaller unit cell than the input. Sometimes it doesn”t
find the smallest possible one, so this method is recursively called
until it is unable to find a smaller cell.

NOTE: if the tolerance is greater than 1/2 the minimum inter-site
distance in the primitive cell, the algorithm will reject this lattice.


* **Parameters**


    * **tolerance** (*float*) – Tolerance for each coordinate of a
    particular site in Angstroms. For example, [0.1, 0, 0.1] in cartesian
    coordinates will be considered to be on the same coordinates
    as [0, 0, 0] for a tolerance of 0.25. Defaults to 0.25.


    * **use_site_props** (*bool*) – Whether to account for site properties in
    differentiating sites.


    * **constrain_latt** (*list/dict*) – List of lattice parameters we want to
    preserve, e.g. [“alpha”, “c”] or dict with the lattice
    parameter names as keys and values we want the parameters to
    be e.g. {“alpha”: 90, “c”: 2.5}.



* **Returns**

    The most primitive structure found.



#### get_reduced_structure(reduction_algo: Literal['niggli', 'LLL'] = 'niggli')
Get a reduced structure.


* **Parameters**

    **reduction_algo** (*"niggli"** | **"LLL"*) – The lattice reduction algorithm to use.
    Defaults to “niggli”.



#### get_sites_in_sphere(pt: ArrayLike, r: float, include_index: bool = False, include_image: bool = False)
Find all sites within a sphere from the point, including a site (if any)
sitting on the point itself. This includes sites in other periodic
images.

Algorithm:


1. place sphere of radius r in crystal and determine minimum supercell
(parallelepiped) which would contain a sphere of radius r. for this
we need the projection of a_1 on a unit vector perpendicular
to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
determine how many a_1”s it will take to contain the sphere.

Nxmax = r \* length_of_b_1 / (2 Pi)


2. keep points falling within r.


* **Parameters**


    * **pt** (*3x1 array*) – Cartesian coordinates of center of sphere.


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Whether the non-supercell site index
    is included in the returned data


    * **include_image** (*bool*) – Whether to include the supercell image
    is included in the returned data



* **Returns**

    PeriodicNeighbor



#### get_sorted_structure(key: Callable | None = None, reverse: bool = False)
Get a sorted copy of the structure. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### get_space_group_info(symprec=0.01, angle_tolerance=5.0)
Convenience method to quickly get the spacegroup of a structure.


* **Parameters**


    * **symprec** (*float*) – Same definition as in SpacegroupAnalyzer.
    Defaults to 1e-2.


    * **angle_tolerance** (*float*) – Same definition as in SpacegroupAnalyzer.
    Defaults to 5 degrees.



* **Returns**

    spacegroup_symbol, international_number



#### get_symmetric_neighbor_list(r: float, sg: str, unique: bool = False, numerical_tol: float = 1e-08, exclude_self: bool = True)
Similar to ‘get_neighbor_list’ with sites=None, but the neighbors are
grouped by symmetry. The returned values are a tuple of numpy arrays
(center_indices, points_indices, offset_vectors, distances,

> symmetry_indices). Atom center_indices[i] has neighbor atom

points_indices[i] that is translated by offset_vectors[i] lattice
vectors, and the distance is distances[i]. Symmetry_idx groups the bonds
that are related by a symmetry of the provided space group and symmetry_op
is the operation that relates the first bond of the same symmetry_idx to
the respective atom. The first bond maps onto itself via the Identity. The
output is sorted w.r.t. to symmetry_indices. If unique is True only one of the
two bonds connecting two points is given. Out of the two, the bond that does not
reverse the sites is chosen.


* **Parameters**


    * **r** (*float*) – Radius of sphere


    * **sg** (*str/int*) – The spacegroup the symmetry operations of which will be
    used to classify the neighbors. If a string, it will be interpreted
    as one of the notations supported by
    pymatgen.symmetry.groups.Spacegroup. E.g., “R-3c” or “Fm-3m”.
    If an int, it will be interpreted as an international number.
    If None, ‘get_space_group_info’ will be used to determine the
    space group, default to None.


    * **unique** (*bool*) – Whether a bond is given for both, or only a single
    direction is given. The default is False.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.


    * **exclude_self** (*bool*) – whether to exclude atom neighboring with itself within
    numerical tolerance distance, default to True


Returns: (center_indices, points_indices, offset_vectors, distances,

    symmetry_indices, symmetry_ops)


#### interpolate(end_structure: IStructure | Structure, nimages: int | Iterable = 10, interpolate_lattices: bool = False, pbc: bool = True, autosort_tol: float = 0)
Interpolate between this structure and end_structure. Useful for
construction of NEB inputs.


* **Parameters**


    * **end_structure** (*Structure*) – structure to interpolate between this
    structure and end.


    * **nimages** (*int**,**list*) – No. of interpolation images or a list of
    interpolation images. Defaults to 10 images.


    * **interpolate_lattices** (*bool*) – Whether to interpolate the lattices.
    Interpolates the lengths and angles (rather than the matrix)
    so orientation may be affected.


    * **pbc** (*bool*) – Whether to use periodic boundary conditions to find
    the shortest path between endpoints.


    * **autosort_tol** (*float*) – A distance tolerance in angstrom in
    which to automatically sort end_structure to match to the
    closest points in this particular structure. This is usually
    what you want in a NEB calculation. 0 implies no sorting.
    Otherwise, a 0.5 value usually works pretty well.



* **Returns**

    List of interpolated structures. The starting and ending
    structures included as the first and last structures respectively.
    A total of (nimages + 1) structures are returned.



#### _property_ is_3d_periodic(_: boo_ )
True if the Lattice is periodic in all directions.


#### _property_ lattice(_: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice_ )
Lattice of the structure.


#### matches(other: IStructure | Structure, anonymous: bool = False, \*\*kwargs)
Check whether this structure is similar to another structure.
Basically a convenience method to call structure matching.


* **Parameters**


    * **other** (*IStructure/Structure*) – Another structure.


    * **anonymous** (*bool*) – Whether to use anonymous structure matching which allows distinct
    species in one structure to map to another.


    * **\*\*kwargs** – Same

    ```
    **
    ```

    kwargs as in
    [`pymatgen.analysis.structure_matcher.StructureMatcher`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher).



* **Returns**

    True if the structures are similar under some affine transformation.



* **Return type**

    bool



#### _property_ pbc(_: tuple[bool, bool, bool_ )
Returns the periodicity of the structure.


#### to(filename: str = '', fmt: str = '', \*\*kwargs)
Outputs the structure to a file or string.


* **Parameters**


    * **filename** (*str*) – If provided, output will be written to a file. If
    fmt is not specified, the format is determined from the
    filename. Defaults is None, i.e. string output.


    * **fmt** (*str*) – Format to output to. Defaults to JSON unless filename
    is provided. If fmt is specifies, it overrides whatever the
    filename is. Options include “cif”, “poscar”, “cssr”, “json”,
    “xsf”, “mcsqs”, “prismatic”, “yaml”, “fleur-inpgen”.
    Non-case sensitive.


    * **\*\*kwargs** – Kwargs passthru to relevant methods. E.g., This allows
    the passing of parameters like symprec to the
    CifWriter.__init__ method for generation of symmetric cifs.



* **Returns**

    (str) if filename is None. None otherwise.



#### unset_charge()
Reset the charge to None, i.e., computed dynamically based on oxidation states.


#### _property_ volume(_: floa_ )
Returns the volume of the structure in Angstrom^3.


### _class_ pymatgen.core.structure.Molecule(species: Sequence[SpeciesLike], coords: Sequence[ArrayLike], charge: float = 0.0, spin_multiplicity: int | None = None, validate_proximity: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, charge_spin_check: bool = True)
Bases: `IMolecule`, `MutableSequence`

Mutable Molecule. It has all the methods in IMolecule, but in addition,
it allows a user to perform edits on the molecule.

Creates a MutableMolecule.


* **Parameters**


    * **species** – list of atomic species. Possible kinds of input include a
    list of dict of elements/species and occupancies, a List of
    elements/specie specified as actual Element/Species, Strings
    (“Fe”, “Fe2+”) or atomic numbers (1,56).


    * **coords** (*3x1 array*) – list of Cartesian coordinates of each species.


    * **charge** (*float*) – Charge for the molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multiplicity for molecule.
    Defaults to None, which means that the spin multiplicity is
    set to 1 if the molecule has no unpaired electrons and to 2
    if there are unpaired electrons.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 1 Ang apart. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as
    a dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The
    sequences have to be the same length as the atomic species
    and fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **charge_spin_check** (*bool*) – Whether to check that the charge and
    spin multiplicity are compatible with each other. Defaults
    to True.



#### append(species: CompositionLike, coords: ArrayLike, validate_proximity: bool = False, properties: dict | None = None)
Appends a site to the molecule.


* **Parameters**


    * **species** – Species of inserted site


    * **coords** – Coordinates of inserted site


    * **validate_proximity** (*bool*) – Whether to check if inserted site is
    too close to an existing site. Defaults to False.


    * **properties** (*dict*) – A dict of properties for the Site.



* **Returns**

    New molecule with inserted site.



#### apply_operation(symmop: [SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp))
Apply a symmetry operation to the molecule.


* **Parameters**

    **symmop** ([*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)) – Symmetry operation to apply.



#### calculate(calculator: str | Calculator = 'gfn2-xtb', verbose: bool = False)
Performs an ASE calculation.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “gfn2-xtb”.
    Defaults to ‘gfn2-xtb’.


    * **verbose** (*bool*) – whether to print stdout. Defaults to False.



* **Returns**

    ASE Calculator instance with a results attribute containing the output.



* **Return type**

    Calculator



#### copy()
Convenience method to get a copy of the molecule.


* **Returns**

    A copy of the Molecule.



#### insert(idx: int, species: CompositionLike, coords: ArrayLike, validate_proximity: bool = False, properties: dict | None = None, label: str | None = None)
Insert a site to the molecule.


* **Parameters**


    * **idx** (*int*) – Index to insert site


    * **species** – species of inserted site


    * **coords** (*3x1 array*) – coordinates of inserted site


    * **validate_proximity** (*bool*) – Whether to check if inserted site is
    too close to an existing site. Defaults to True.


    * **properties** (*dict*) – Dict of properties for the Site.


    * **label** (*str*) – Label of inserted site



* **Returns**

    New molecule with inserted site.



#### perturb(distance: float)
Performs a random perturbation of the sites in a structure to break
symmetries.


* **Parameters**

    **distance** (*float*) – Distance in angstroms by which to perturb each
    site.



#### relax(calculator: str | Calculator = 'gfn2-xtb', optimizer: str | Optimizer = 'FIRE', steps: int = 500, fmax: float = 0.1, opt_kwargs: dict | None = None, return_trajectory: bool = False, verbose: bool = False)
Performs a molecule relaxation using an ASE calculator.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “gfn2-xtb”.
    Defaults to ‘gfn2-xtb’.


    * **optimizer** (*str*) – name of the ASE optimizer class to use


    * **steps** (*int*) – max number of steps for relaxation. Defaults to 500.


    * **fmax** (*float*) – total force tolerance for relaxation convergence.
    Defaults to 0.1 eV/A.


    * **opt_kwargs** (*dict*) – kwargs for the ASE optimizer class.


    * **return_trajectory** (*bool*) – Whether to return the trajectory of relaxation.
    Defaults to False.


    * **verbose** (*bool*) – whether to print out relaxation steps. Defaults to False.



* **Returns**

    Relaxed Molecule or if return_trajectory=True,

        2-tuple of Molecule and ASE TrajectoryObserver.




* **Return type**

    Molecule | tuple[Molecule, [Trajectory](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory)]



#### remove_sites(indices: Sequence[int])
Delete sites with at indices.


* **Parameters**

    **indices** – Sequence of indices of sites to delete.



#### remove_species(species: Sequence[SpeciesLike])
Remove all occurrences of a species from a molecule.


* **Parameters**

    **species** – Species to remove.



#### rotate_sites(indices: Sequence[int] | None = None, theta: float = 0.0, axis: ArrayLike | None = None, anchor: ArrayLike | None = None)
Rotate specific sites by some angle around vector at anchor.


* **Parameters**


    * **indices** (*list*) – List of site indices on which to perform the
    translation.


    * **theta** (*float*) – Angle in radians


    * **axis** (*3x1 array*) – Rotation axis vector.


    * **anchor** (*3x1 array*) – Point of rotation.



#### set_charge_and_spin(charge: float, spin_multiplicity: int | None = None)
Set the charge and spin multiplicity.


* **Parameters**


    * **charge** (*int*) – Charge for the molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multiplicity for molecule.
    Defaults to None, which means that the spin multiplicity is
    set to 1 if the molecule has no unpaired electrons and to 2
    if there are unpaired electrons.



#### substitute(index: int, func_group: IMolecule | Molecule | str, bond_order: int = 1)
Substitute atom at index with a functional group.


* **Parameters**


    * **index** (*int*) – Index of atom to substitute.


    * **func_group** – Substituent molecule. There are two options:


        1. Providing an actual molecule as the input. The first atom
    must be a DummySpecies X, indicating the position of
    nearest neighbor. The second atom must be the next
    nearest atom. For example, for a methyl group
    substitution, func_group should be X-CH3, where X is the
    first site and C is the second site. What the code will
    do is to remove the index site, and connect the nearest
    neighbor to the C atom in CH3. The X-C bond indicates the
    directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the
    relevant template in func_groups.json.



    * **bond_order** (*int*) – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.



#### translate_sites(indices: Sequence[int] | None = None, vector: ArrayLike | None = None)
Translate specific sites by some vector, keeping the sites within the
unit cell.


* **Parameters**


    * **indices** (*list*) – List of site indices on which to perform the
    translation.


    * **vector** (*3x1 array*) – Translation vector for sites.



### _class_ pymatgen.core.structure.Neighbor(species: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), coords: np.ndarray, properties: dict | None = None, nn_distance: float = 0.0, index: int = 0, label: str | None = None)
Bases: [`Site`](pymatgen.core.sites.md#pymatgen.core.sites.Site)

Simple Site subclass to contain a neighboring atom that skips all the unnecessary checks for speed. Can be
used as a fixed-length tuple of size 3 to retain backwards compatibility with past use cases.

> (site, nn_distance, index).

In future, usage should be to call attributes, e.g., Neighbor.index, Neighbor.distance, etc.


* **Parameters**


    * **species** – Same as Site


    * **coords** – Same as Site, but must be fractional.


    * **properties** – Same as Site


    * **nn_distance** – Distance to some other Site.


    * **index** – Index within structure.


    * **label** – Label for the site. Defaults to None.



#### as_dict()
Note that method calls the super of Site, which is MSONable itself.

Returns: dict


#### coords(_: ndarra_ )

#### _classmethod_ from_dict(d: dict)
Returns a Neighbor from a dict.


* **Parameters**

    **d** – MSONable dict format.



* **Returns**

    Neighbor



#### properties(_: dic_ )

### _class_ pymatgen.core.structure.PeriodicNeighbor(species: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition), coords: np.ndarray, lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), properties: dict | None = None, nn_distance: float = 0.0, index: int = 0, image: tuple = (0, 0, 0), label: str | None = None)
Bases: [`PeriodicSite`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)

Simple PeriodicSite subclass to contain a neighboring atom that skips all
the unnecessary checks for speed. Can be used as a fixed-length tuple of
size 4 to retain backwards compatibility with past use cases.

> (site, distance, index, image).

In future, usage should be to call attributes, e.g., PeriodicNeighbor.index,
PeriodicNeighbor.distance, etc.


* **Parameters**


    * **species** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Same as PeriodicSite


    * **coords** (*np.ndarray*) – Same as PeriodicSite, but must be fractional.


    * **lattice** ([*Lattice*](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice)) – Same as PeriodicSite


    * **properties** (*dict**, **optional*) – Same as PeriodicSite. Defaults to None.


    * **nn_distance** (*float**, **optional*) – Distance to some other Site.. Defaults to 0.0.


    * **index** (*int**, **optional*) – Index within structure.. Defaults to 0.


    * **image** (*tuple**, **optional*) – PeriodicImage. Defaults to (0, 0, 0).


    * **label** (*str**, **optional*) – Label for the site. Defaults to None.



#### as_dict()
Note that method calls the super of Site, which is MSONable itself.

Returns: dict


#### _property_ coords(_: ndarra_ )
Cartesian coords.


* **Type**

    return



#### _classmethod_ from_dict(d: dict)
Returns a PeriodicNeighbor from a dict.


* **Parameters**

    **d** – MSONable dict format.



* **Returns**

    PeriodicNeighbor



#### properties(_: dic_ )

### _class_ pymatgen.core.structure.SiteCollection()
Bases: `Sequence`

Basic SiteCollection. Essentially a sequence of Sites or PeriodicSites.
This serves as a base class for Molecule (a collection of Site, i.e., no
periodicity) and Structure (a collection of PeriodicSites, i.e.,
periodicity). Not meant to be instantiated directly.


#### DISTANCE_TOLERANCE(_ = 0._ )

#### add_oxidation_state_by_element(oxidation_states: dict[str, float])
Add oxidation states.


* **Parameters**

    **oxidation_states** (*dict*) – Dict of oxidation states.
    E.g., {“Li”:1, “Fe”:2, “P”:5, “O”:-2}



* **Raises**

    **ValueError if oxidation states are not specified for all elements.** –



#### add_oxidation_state_by_guess(\*\*kwargs)
Decorates the structure with oxidation state, guessing
using Composition.oxi_state_guesses().


* **Parameters**

    **\*\*kwargs** – parameters to pass into oxi_state_guesses()



#### add_oxidation_state_by_site(oxidation_states: list[float])
Add oxidation states to a structure by site.


* **Parameters**

    **oxidation_states** (*list**[**float**]*) – List of oxidation states.
    E.g. [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2]



#### add_site_property(property_name: str, values: list)
Adds a property to a site. Note: This is the preferred method
for adding magnetic moments, selective dynamics, and related
site-specific properties to a structure/molecule object.

### Examples

structure.add_site_property(“magmom”, [1.0, 0.0])
structure.add_site_property(“selective_dynamics”, [[True, True, True], [False, False, False]])


* **Parameters**


    * **property_name** (*str*) – The name of the property to add.


    * **values** (*list*) – A sequence of values. Must be same length as
    number of sites.



#### add_spin_by_element(spins: dict[str, float])
Add spin states to structure.


* **Parameters**

    **spins** (*dict*) – Dict of spins associated with elements or species,
    e.g. {“Ni”:+5} or {“Ni2+”:5}



#### add_spin_by_site(spins: list[float])
Add spin states to structure by site.


* **Parameters**

    **spins** (*list*) – List of spins
    E.g., [+5, -5, 0, 0]



#### _property_ atomic_numbers(_: tuple[int, ..._ )
List of atomic numbers.


#### _property_ cart_coords(_: ndarra_ )
Returns an np.array of the Cartesian coordinates of sites in the
structure.


#### _property_ charge(_: floa_ )
Returns the net charge of the structure based on oxidation states. If
Elements are found, a charge of 0 is assumed.


#### _property_ composition(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
(Composition) Returns the composition.


#### _property_ distance_matrix(_: ndarra_ )
Returns the distance matrix between all sites in the structure. For
periodic structures, this is overwritten to return the nearest image
distance.


#### extract_cluster(target_sites: list[[pymatgen.core.sites.Site](pymatgen.core.sites.md#pymatgen.core.sites.Site)], \*\*kwargs)
Extracts a cluster of atoms based on bond lengths.


* **Parameters**


    * **target_sites** (*list**[*[*Site*](pymatgen.core.sites.md#pymatgen.core.sites.Site)*]*) – Initial sites from which to nucleate cluster.


    * **\*\*kwargs** – kwargs passed through to CovalentBond.is_bonded.



* **Returns**

    list[Site/PeriodicSite] Cluster of atoms.



#### _property_ formula(_: st_ )
(str) Returns the formula.


#### _abstract classmethod_ from_file(filename: str)
Reads in SiteCollection from a filename.


#### _abstract classmethod_ from_str(input_string: str, fmt: Any)
Reads in SiteCollection from a string.


#### get_angle(i: int, j: int, k: int)
Returns angle specified by three sites.


* **Parameters**


    * **i** – 1st site index


    * **j** – 2nd site index


    * **k** – 3rd site index



* **Returns**

    Angle in degrees.



#### get_dihedral(i: int, j: int, k: int, l: int)
Returns dihedral angle specified by four sites.


* **Parameters**


    * **i** – 1st site index


    * **j** – 2nd site index


    * **k** – 3rd site index


    * **l** – 4th site index



* **Returns**

    Dihedral angle in degrees.



#### _abstract_ get_distance(i: int, j: int)
Returns distance between sites at index i and j.


* **Parameters**


    * **i** – 1st site index


    * **j** – 2nd site index



* **Returns**

    Distance between sites at index i and index j.



#### group_by_types()
Iterate over species grouped by type.


#### indices_from_symbol(symbol: str)
Returns a tuple with the sequential indices of the sites
that contain an element with the given chemical symbol.


#### _property_ is_ordered(_: boo_ )
Checks if structure is ordered, meaning no partial occupancies in any
of the sites.


#### is_valid(tol: float = 0.5)
True if SiteCollection does not contain atoms that are too close
together. Note that the distance definition is based on type of
SiteCollection. Cartesian distances are used for non-periodic
Molecules, while PBC is taken into account for periodic structures.


* **Parameters**

    **tol** (*float*) – Distance tolerance. Default is 0.5 Angstrom, which is fairly large.



* **Returns**

    (bool) True if SiteCollection does not contain atoms that are too close together.



#### _property_ labels(_: list[str_ )
Return site labels as a list.


#### _property_ ntypesp(_: in_ )
Number of types of atoms.


#### _property_ num_sites(_: in_ )
Number of sites.


#### remove_oxidation_states()
Removes oxidation states from a structure.


#### remove_site_property(property_name: str)
Removes a property to a site.


* **Parameters**

    **property_name** (*str*) – The name of the property to remove.



#### remove_spin()
Remove spin states from structure.


#### replace_species(species_mapping: dict[SpeciesLike, SpeciesLike | dict[SpeciesLike, float]])
Swap species. Note that this method modifies the structure in place.


* **Parameters**

    **species_mapping** (*dict*) – dict of species to swap. Species can be elements too. E.g.,
    {Element(“Li”): Element(“Na”)} performs a Li for Na substitution. The second species can
    be a sp_and_occu dict. For example, a site with 0.5 Si that is passed the mapping
    {Element(‘Si): {Element(‘Ge’): 0.75, Element(‘C’): 0.25} } will have .375 Ge and .125 C.



#### _property_ site_properties(_: dict[str, Sequence_ )
(-4,4)}.


* **Type**

    Returns the site properties as a dict of sequences. E.g. {“magmom”



* **Type**

    (5,-5), “charge”



#### _property_ sites(_: list[[pymatgen.core.sites.Site](pymatgen.core.sites.md#pymatgen.core.sites.Site)_ )
Returns an iterator for the sites in the Structure.


#### _property_ species(_: list[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)_ )
Only works for ordered structures.


* **Raises**

    **AttributeError** – If structure is disordered.



* **Returns**

    ([Species]) List of species at each site of the structure.



#### _property_ species_and_occu(_: list[[pymatgen.core.composition.Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition)_ )
List of species and occupancies at each site of the structure.


#### _property_ symbol_set(_: tuple[str, ..._ )
Tuple with the set of chemical symbols.
Note that len(symbol_set) == len(types_of_specie).


#### _abstract_ to(filename: str = '', fmt: str = '')
Generates string representations (cif, json, poscar, ….) of SiteCollections (e.g.,
molecules / structures). Should return str or None if written to a file.


#### _property_ types_of_specie(_: tuple[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies)_ )
Specie->Species rename. Maintained for backwards compatibility.


#### _property_ types_of_species(_: tuple[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies)_ )
List of types of specie.


### _class_ pymatgen.core.structure.Structure(lattice: ArrayLike | [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False, coords_are_cartesian: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None)
Bases: `IStructure`, `MutableSequence`

Mutable version of structure.

Create a periodic structure.


* **Parameters**


    * **lattice** – The lattice, either as a pymatgen.core.lattice.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** – List of species on each site. Can take in flexible input,
    including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **charge** (*int*) – overall charge of the structure. Defaults to behavior
    in SiteCollection where total charge is the sum of the oxidation
    states.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to map all sites into the unit cell,
    i.e., fractional coords between 0 and 1. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.



#### append(species: CompositionLike, coords: ArrayLike, coords_are_cartesian: bool = False, validate_proximity: bool = False, properties: dict | None = None)
Append a site to the structure.


* **Parameters**


    * **species** – Species of inserted site


    * **coords** (*3x1 array*) – Coordinates of inserted site


    * **coords_are_cartesian** (*bool*) – Whether coordinates are cartesian.
    Defaults to False.


    * **validate_proximity** (*bool*) – Whether to check if inserted site is
    too close to an existing site. Defaults to False.


    * **properties** (*dict*) – Properties of the site.



* **Returns**

    New structure with inserted site.



#### apply_operation(symmop: [SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp), fractional: bool = False)
Apply a symmetry operation to the structure in place and return the modified
structure. The lattice is operated on by the rotation matrix only.
Coords are operated in full and then transformed to the new lattice.


* **Parameters**


    * **symmop** ([*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)) – Symmetry operation to apply.


    * **fractional** (*bool*) – Whether the symmetry operation is applied in
    fractional space. Defaults to False, i.e., symmetry operation
    is applied in Cartesian coordinates.



* **Returns**

    post-operation structure



* **Return type**

    Structure



#### apply_strain(strain: ArrayLike)
Apply a strain to the lattice.


* **Parameters**

    **strain** (*float** or **list*) – Amount of strain to apply. Can be a float,
    or a sequence of 3 numbers. E.g., 0.01 means all lattice
    vectors are increased by 1%. This is equivalent to calling
    modify_lattice with a lattice with lattice parameters that
    are 1% larger.



#### calculate(calculator: str | Calculator = 'm3gnet', verbose: bool = False)
Performs an ASE calculation.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “m3gnet”.
    Defaults to ‘m3gnet’, i.e. the M3GNet universal potential.


    * **verbose** (*bool*) – whether to print stdout. Defaults to False.



* **Returns**

    ASE Calculator instance with a results attribute containing the output.



* **Return type**

    Calculator



#### _classmethod_ from_prototype(prototype: str, species: Sequence, \*\*kwargs)
Method to rapidly construct common prototype structures.


* **Parameters**


    * **prototype** – Name of prototype. E.g., cubic, rocksalt, perovksite etc.


    * **species** – List of species corresponding to symmetrically distinct sites.


    * **\*\*kwargs** – Lattice parameters, e.g., a = 3.0, b = 4, c = 5. Only the required lattice parameters need to be
    specified. For example, if it is a cubic prototype, only a needs to be specified.



* **Returns**

    Structure



#### insert(idx: int, species: CompositionLike, coords: ArrayLike, coords_are_cartesian: bool = False, validate_proximity: bool = False, properties: dict | None = None, label: str | None = None)
Insert a site to the structure.


* **Parameters**


    * **idx** (*int*) – Index to insert site


    * **species** (*species-like*) – Species of inserted site


    * **coords** (*3x1 array*) – Coordinates of inserted site


    * **coords_are_cartesian** (*bool*) – Whether coordinates are cartesian.
    Defaults to False.


    * **validate_proximity** (*bool*) – Whether to check if inserted site is too close to
    an existing site. Controlled by self.DISTANCE_TOLERANCE. Defaults to False.


    * **properties** (*dict*) – Properties associated with the site.


    * **label** (*str*) – Label associated with the site.



* **Returns**

    New structure with inserted site.



#### _property_ lattice(_: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice_ )
Lattice associated with structure.


* **Type**

    return



#### make_supercell(scaling_matrix: ArrayLike, to_unit_cell: bool = True, in_place: bool = True)
Create a supercell.


* **Parameters**


    * **scaling_matrix** (*ArrayLike*) – A scaling matrix for transforming the lattice
    vectors. Has to be all integers. Several options are possible:


        1. A full 3x3 scaling matrix defining the linear combination
    the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
    1]] generates a new structure with lattice vectors a’ =
    2a + b, b’ = 3b, c’ = c where a, b, and c are the lattice
    vectors of the original structure.


        2. An sequence of three scaling factors. E.g., [2, 1, 1]
    specifies that the supercell should have dimensions 2a x b x
    c.


        3. A number, which simply scales all lattice vectors by the
    same factor.



    * **to_unit_cell** (*bool*) – Whether or not to fold sites back into the unit cell
    if they have fractional coords > 1. Defaults to True.


    * **in_place** (*bool*) – Whether to perform the operation in-place or to return
    a new Structure object. Defaults to True.



* **Returns**

    self if in_place is True else self.copy() after making supercell



* **Return type**

    Structure



#### merge_sites(tol: float = 0.01, mode: Literal['sum', 'delete', 'average'] = 'sum')
Merges sites (adding occupancies) within tol of each other.
Removes site properties.


* **Parameters**


    * **tol** (*float*) – Tolerance for distance to merge sites.


    * **mode** (*'sum'** | **'delete'** | **'average'*) – “delete” means duplicate sites are
    deleted. “sum” means the occupancies are summed for the sites.
    “average” means that the site is deleted but the properties are averaged
    Only first letter is considered.



#### perturb(distance: float, min_distance: float | None = None)
Performs a random perturbation of the sites in a structure to break
symmetries.


* **Parameters**


    * **distance** (*float*) – Distance in angstroms by which to perturb each
    site.


    * **min_distance** (*None**, **int**, or **float*) – if None, all displacements will
    be equal amplitude. If int or float, perturb each site a
    distance drawn from the uniform distribution between
    ‘min_distance’ and ‘distance’.



#### relax(calculator: str | Calculator = 'm3gnet', relax_cell: bool = True, optimizer: str | Optimizer = 'FIRE', steps: int = 500, fmax: float = 0.1, stress_weight: float = 0.01, opt_kwargs: dict | None = None, return_trajectory: bool = False, verbose: bool = False)
Performs a crystal structure relaxation using an ASE calculator.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “m3gnet”.
    Defaults to ‘m3gnet’, i.e. the M3GNet universal potential.


    * **relax_cell** (*bool*) – whether to relax the lattice cell. Defaults to True.


    * **optimizer** (*str*) – name of the ASE optimizer class to use


    * **steps** (*int*) – max number of steps for relaxation. Defaults to 500.


    * **fmax** (*float*) – total force tolerance for relaxation convergence.
    Here fmax is a sum of force and stress forces. Defaults to 0.1.


    * **stress_weight** (*float*) – the stress weight for relaxation with M3GNet.
    Defaults to 0.01.


    * **opt_kwargs** (*dict*) – kwargs for the ASE optimizer class.


    * **return_trajectory** (*bool*) – Whether to return the trajectory of relaxation.
    Defaults to False.


    * **verbose** (*bool*) – whether to print out relaxation steps. Defaults to False.



* **Returns**

    Relaxed structure or if return_trajectory=True,

        2-tuple of Structure and matgl TrajectoryObserver.




* **Return type**

    Structure | tuple[Structure, [Trajectory](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory)]



#### remove_sites(indices: Sequence[int])
Delete sites with at indices.


* **Parameters**

    **indices** – Sequence of indices of sites to delete.



#### remove_species(species: Sequence[SpeciesLike])
Remove all occurrences of several species from a structure.


* **Parameters**

    **species** – Sequence of species to remove, e.g., [“Li”, “Na”].



#### replace(idx: int, species: CompositionLike, coords: ArrayLike | None = None, coords_are_cartesian: bool = False, properties: dict | None = None, label: str | None = None)
Replace a single site. Takes either a species or a dict of species and
occupations.


* **Parameters**


    * **idx** (*int*) – Index of the site in the sites list.


    * **species** (*species-like*) – Species of replacement site


    * **coords** (*3x1 array*) – Coordinates of replacement site. If None,
    the current coordinates are assumed.


    * **coords_are_cartesian** (*bool*) – Whether coordinates are cartesian.
    Defaults to False.


    * **properties** (*dict*) – Properties associated with the site.


    * **label** (*str*) – Label associated with the site.



#### rotate_sites(indices: list[int] | None = None, theta: float = 0.0, axis: ArrayLike | None = None, anchor: ArrayLike | None = None, to_unit_cell: bool = True)
Rotate specific sites by some angle around vector at anchor.


* **Parameters**


    * **indices** (*list*) – List of site indices on which to perform the
    translation.


    * **theta** (*float*) – Angle in radians


    * **axis** (*3x1 array*) – Rotation axis vector.


    * **anchor** (*3x1 array*) – Point of rotation.


    * **to_unit_cell** (*bool*) – Whether new sites are transformed to unit
    cell



#### scale_lattice(volume: float)
Performs a scaling of the lattice vectors so that length proportions
and angles are preserved.


* **Parameters**

    **volume** (*float*) – New volume of the unit cell in A^3.



#### set_charge(new_charge: float = 0.0)
Sets the overall structure charge.


* **Parameters**

    **new_charge** (*float*) – new charge to set



#### sort(key: Callable | None = None, reverse: bool = False)
Sort a structure in place. The parameters have the same meaning as in
list.sort. By default, sites are sorted by the electronegativity of
the species. The difference between this method and
get_sorted_structure (which also works in IStructure) is that the
latter returns a new Structure, while this just sorts the Structure
in place.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### substitute(index: int, func_group: IMolecule | Molecule | str, bond_order: int = 1)
Substitute atom at index with a functional group.


* **Parameters**


    * **index** (*int*) – Index of atom to substitute.


    * **func_group** – Substituent molecule. There are two options:


        1. Providing an actual Molecule as the input. The first atom
    must be a DummySpecies X, indicating the position of
    nearest neighbor. The second atom must be the next
    nearest atom. For example, for a methyl group
    substitution, func_group should be X-CH3, where X is the
    first site and C is the second site. What the code will
    do is to remove the index site, and connect the nearest
    neighbor to the C atom in CH3. The X-C bond indicates the
    directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the
    relevant template in func_groups.json.



    * **bond_order** (*int*) – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.



#### translate_sites(indices: int | Sequence[int], vector: ArrayLike, frac_coords: bool = True, to_unit_cell: bool = True)
Translate specific sites by some vector, keeping the sites within the
unit cell.


* **Parameters**


    * **indices** – Integer or List of site indices on which to perform the
    translation.


    * **vector** – Translation vector for sites.


    * **frac_coords** (*bool*) – Whether the vector corresponds to fractional or
    Cartesian coordinates.


    * **to_unit_cell** (*bool*) – Whether new sites are transformed to unit
    cell



### _exception_ pymatgen.core.structure.StructureError()
Bases: `Exception`

Exception class for Structure.
Raised when the structure has problems, e.g., atoms that are too close.