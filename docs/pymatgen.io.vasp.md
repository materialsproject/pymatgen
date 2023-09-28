---
layout: default
title: pymatgen.io.vasp.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.vasp package

This package implements modules for input and output to and from VASP. It
imports the key classes form both vasp_input and vasp_output to allow most
classes to be simply called as pymatgen.io.vasp.Incar for example, to retain
backwards compatibility.


## pymatgen.io.vasp.help module

Get help with VASP parameters from VASP wiki.


### _class_ VaspDoc()
Bases: `object`

A VASP documentation helper.

Init for VaspDoc.


#### _classmethod_ get_help(tag, fmt='text')
Get help on a VASP tag.


* **Parameters**

    **tag** (*str*) – VASP tag, e.g., ISYM.



* **Returns**

    Help text.



#### _classmethod_ get_incar_tags()
Returns: All incar tags.


#### print_help(tag)
Print the help for a TAG.


* **Parameters**

    **tag** (*str*) – Tag used in VASP.



#### print_jupyter_help(tag)
Display HTML help in ipython notebook.


* **Parameters**

    **tag** (*str*) – Tag used in VASP.


## pymatgen.io.vasp.inputs module

Classes for reading/manipulating/writing VASP input files. All major VASP input
files.


### _exception_ BadIncarWarning()
Bases: `UserWarning`

Warning class for bad Incar parameters.


### _class_ Incar(params: dict[str, Any] | None = None)
Bases: `dict`, `MSONable`

INCAR object for reading and writing INCAR files. Essentially consists of
a dictionary with some helper functions.

Creates an Incar object.


* **Parameters**

    **params** (*dict*) – A set of input parameters as a dictionary.



#### as_dict()
MSONable dict.


#### check_params()
Raises a warning for nonsensical or non-existent INCAR tags and
parameters. If a keyword doesn’t exist (e.g. there’s a typo in a
keyword), your calculation will still run, however VASP will ignore the
parameter without letting you know, hence why we have this Incar method.


#### diff(other: Incar)
Diff function for Incar. Compares two Incars and indicates which
parameters are the same and which are not. Useful for checking whether
two runs were done using the same parameters.


* **Parameters**

    **other** (*Incar*) – The other Incar object to compare to.



* **Returns**

    {“Same” : parameters_that_are_the_same,
    “Different”: parameters_that_are_different}
    Note that the parameters are return as full dictionaries of values.
    E.g. {“ISIF”:3}



* **Return type**

    Dict of the following format



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Incar



#### _static_ from_file(filename: PathLike)
Reads an Incar object from a file.


* **Parameters**

    **filename** (*str*) – Filename for file



* **Returns**

    Incar object



#### _static_ from_str(string: str)
Reads an Incar object from a string.


* **Parameters**

    **string** (*str*) – Incar string



* **Returns**

    Incar object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_str(sort_keys: bool = False, pretty: bool = False)
Returns a string representation of the INCAR. The reason why this
method is different from the __str__ method is to provide options for
pretty printing.


* **Parameters**


    * **sort_keys** (*bool*) – Set to True to sort the INCAR parameters
    alphabetically. Defaults to False.


    * **pretty** (*bool*) – Set to True for pretty aligned output. Defaults
    to False.



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### _static_ proc_val(key: str, val: Any)
Static helper method to convert INCAR parameters to proper types, e.g.,
integers, floats, lists, etc.


* **Parameters**


    * **key** – INCAR parameter key


    * **val** – Actual value of INCAR parameter.



#### write_file(filename: PathLike)
Write Incar to a file.


* **Parameters**

    **filename** (*str*) – filename to write to.



### _class_ Kpoints(comment: str = 'Default gamma', num_kpts: int = 0, style: KpointsSupportedModes = KpointsSupportedModes.Gamma, kpts: Sequence[float | Sequence] = ((1, 1, 1),), kpts_shift: Vector3D = (0, 0, 0), kpts_weights=None, coord_type=None, labels=None, tet_number: int = 0, tet_weight: float = 0, tet_connections=None)
Bases: `MSONable`

KPOINT reader/writer.

Highly flexible constructor for Kpoints object. The flexibility comes
at the cost of usability and in general, it is recommended that you use
the default constructor only if you know exactly what you are doing and
requires the flexibility. For most usage cases, the three automatic
schemes can be constructed far more easily using the convenience static
constructors (automatic, gamma_automatic, monkhorst_automatic) and it
is recommended that you use those.


* **Parameters**


    * **comment** (*str*) – String comment for Kpoints. Defaults to “Default gamma”.


    * **num_kpts** – Following VASP method of defining the KPOINTS file, this
    parameter is the number of kpoints specified. If set to 0
    (or negative), VASP automatically generates the KPOINTS.


    * **style** – Style for generating KPOINTS. Use one of the
    Kpoints.supported_modes enum types.


    * **kpts** (*2D array*) – 2D array of kpoints. Even when only a single
    specification is required, e.g. in the automatic scheme,
    the kpts should still be specified as a 2D array. e.g.,
    [[20]] or [[2,2,2]].


    * **kpts_shift** (*3x1 array*) – Shift for Kpoints.


    * **kpts_weights** – Optional weights for kpoints. Weights should be
    integers. For explicit kpoints.


    * **coord_type** – In line-mode, this variable specifies whether the
    Kpoints were given in Cartesian or Reciprocal coordinates.


    * **labels** – In line-mode, this should provide a list of labels for
    each kpt. It is optional in explicit kpoint mode as comments for
    k-points.


    * **tet_number** – For explicit kpoints, specifies the number of
    tetrahedrons for the tetrahedron method.


    * **tet_weight** – For explicit kpoints, specifies the weight for each
    tetrahedron for the tetrahedron method.


    * **tet_connections** – For explicit kpoints, specifies the connections
    of the tetrahedrons for the tetrahedron method.
    Format is a list of tuples, [ (sym_weight, [tet_vertices]),
    …]


The default behavior of the constructor is for a Gamma centered,
1x1x1 KPOINTS with no shift.


#### as_dict()
MSONable dict.


#### _static_ automatic(subdivisions)
Convenient static constructor for a fully automatic Kpoint grid, with
gamma centered Monkhorst-Pack grids and the number of subdivisions
along each reciprocal lattice vector determined by the scheme in the
VASP manual.


* **Parameters**

    **subdivisions** – Parameter determining number of subdivisions along
    each reciprocal lattice vector.



* **Returns**

    Kpoints object



#### _static_ automatic_density(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), kppa: float, force_gamma: bool = False)
Returns an automatic Kpoint object based on a structure and a kpoint
density. Uses Gamma centered meshes for hexagonal cells and face-centered cells,
Monkhorst-Pack grids otherwise.

Algorithm:

    Uses a simple approach scaling the number of divisions along each
    reciprocal lattice vector proportional to its length.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **kppa** (*float*) – Grid density


    * **force_gamma** (*bool*) – Force a gamma centered mesh (default is to
    use gamma only for hexagonal cells or odd meshes)



* **Returns**

    Kpoints



#### _static_ automatic_density_by_lengths(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), length_densities: Sequence[float], force_gamma: bool = False)
Returns an automatic Kpoint object based on a structure and a k-point
density normalized by lattice constants.

Algorithm:

    For a given dimension, the # of k-points is chosen as
    length_density = # of kpoints \* lattice constant, e.g. [50.0, 50.0, 1.0] would
    have k-points of 50/a x 50/b x 1/c.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **length_densities** (*list**[**floats**]*) – Defines the density of k-points in each


    * **dimension** –


    * **[****50.0** (*e.g.*) –


    * **50.0** –


    * **1.0****]****.** –


    * **force_gamma** (*bool*) – Force a gamma centered mesh



* **Returns**

    Kpoints



#### _static_ automatic_density_by_vol(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), kppvol: int, force_gamma: bool = False)
Returns an automatic Kpoint object based on a structure and a kpoint
density per inverse Angstrom^3 of reciprocal cell.

Algorithm:

    Same as automatic_density()


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **kppvol** (*int*) – Grid density per Angstrom^(-3) of reciprocal cell


    * **force_gamma** (*bool*) – Force a gamma centered mesh



* **Returns**

    Kpoints



#### _static_ automatic_gamma_density(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), kppa: float)
Returns an automatic Kpoint object based on a structure and a kpoint
density. Uses Gamma centered meshes always. For GW.

Algorithm:

    Uses a simple approach scaling the number of divisions along each
    reciprocal lattice vector proportional to its length.


* **Parameters**


    * **structure** – Input structure


    * **kppa** – Grid density



#### _static_ automatic_linemode(divisions, ibz)
Convenient static constructor for a KPOINTS in mode line_mode.
gamma centered Monkhorst-Pack grids and the number of subdivisions
along each reciprocal lattice vector determined by the scheme in the
VASP manual.


* **Parameters**


    * **divisions** – Parameter determining the number of k-points along each high symmetry line.


    * **ibz** – HighSymmKpath object (pymatgen.symmetry.bandstructure)



* **Returns**

    Kpoints object



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Kpoints



#### _static_ from_file(filename)
Reads a Kpoints object from a KPOINTS file.


* **Parameters**

    **filename** (*str*) – filename to read from.



* **Returns**

    Kpoints object



#### _static_ from_str(string)
Reads a Kpoints object from a KPOINTS string.


* **Parameters**

    **string** (*str*) – KPOINTS string.



* **Returns**

    Kpoints object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### _static_ gamma_automatic(kpts: tuple[int, int, int] = (1, 1, 1), shift: Vector3D = (0, 0, 0))
Convenient static constructor for an automatic Gamma centered Kpoint
grid.


* **Parameters**


    * **kpts** – Subdivisions N_1, N_2 and N_3 along reciprocal lattice
    vectors. Defaults to (1,1,1)


    * **shift** – Shift to be applied to the kpoints. Defaults to (0,0,0).



* **Returns**

    Kpoints object



#### _static_ monkhorst_automatic(kpts: tuple[int, int, int] = (2, 2, 2), shift: Vector3D = (0, 0, 0))
Convenient static constructor for an automatic Monkhorst pack Kpoint
grid.


* **Parameters**


    * **kpts** – Subdivisions N_1, N_2, N_3 along reciprocal lattice
    vectors. Defaults to (2,2,2)


    * **shift** – Shift to be applied to the kpoints. Defaults to (0,0,0).



* **Returns**

    Kpoints object



#### _property_ style(_: KpointsSupportedMode_ )
Style for kpoint generation. One of Kpoints_supported_modes enum.


#### supported_modes()
alias of `KpointsSupportedModes`


#### write_file(filename)
Write Kpoints to a file.


* **Parameters**

    **filename** (*str*) – Filename to write to.



### _class_ KpointsSupportedModes(value)
Bases: `Enum`

Enum type of all supported modes for Kpoint generation.


#### Automatic(_ = _ )

#### Cartesian(_ = _ )

#### Gamma(_ = _ )

#### Line_mode(_ = _ )

#### Monkhorst(_ = _ )

#### Reciprocal(_ = _ )

#### _static_ from_str(s: str)

* **Parameters**

    **s** – String



* **Returns**

    Kpoints_supported_modes



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


### _class_ Orbital(n, l, j, E, occ)
Bases: `tuple`

Create new instance of Orbital(n, l, j, E, occ)


#### E()
Alias for field number 3


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('n', 'l', 'j', 'E', 'occ'_ )

#### _classmethod_ _make(iterable)
Make a new Orbital object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new Orbital object replacing specified fields with new values


#### j()
Alias for field number 2


#### l()
Alias for field number 1


#### n()
Alias for field number 0


#### occ()
Alias for field number 4


### _class_ OrbitalDescription(l, E, Type, Rcut, Type2, Rcut2)
Bases: `tuple`

Create new instance of OrbitalDescription(l, E, Type, Rcut, Type2, Rcut2)


#### E()
Alias for field number 1


#### Rcut()
Alias for field number 3


#### Rcut2()
Alias for field number 5


#### Type()
Alias for field number 2


#### Type2()
Alias for field number 4


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('l', 'E', 'Type', 'Rcut', 'Type2', 'Rcut2'_ )

#### _classmethod_ _make(iterable)
Make a new OrbitalDescription object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new OrbitalDescription object replacing specified fields with new values


#### l()
Alias for field number 0


### _class_ Poscar(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), comment: str | None = None, selective_dynamics: ArrayLike | None = None, true_names: bool = True, velocities: ArrayLike | None = None, predictor_corrector: ArrayLike | None = None, predictor_corrector_preamble: str | None = None, sort_structure: bool = False)
Bases: `MSONable`

Object for representing the data in a POSCAR or CONTCAR file.


#### structure()
Associated Structure.


#### comment()
Optional comment string.


#### true_names()
Boolean indication whether Poscar contains actual real names parsed
from either a POTCAR or the POSCAR itself.


#### selective_dynamics()
Selective dynamics attribute for each site if available.
A Nx3 array of booleans.


#### velocities()
Velocities for each site (typically read in from a CONTCAR).
A Nx3 array of floats.


#### predictor_corrector()
Predictor corrector coordinates and derivatives for each site;
i.e. a list of three 1x3 arrays for each site (typically read in from a MD CONTCAR).


#### predictor_corrector_preamble()
Predictor corrector preamble contains the predictor-corrector key,
POTIM, and thermostat parameters that precede the site-specific predictor corrector data in MD CONTCAR.


#### temperature()
Temperature of velocity Maxwell-Boltzmann initialization.
Initialized to -1 (MB hasn’t been performed).


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object.


    * **comment** (*str** | **None**, **optional*) – Optional comment line for POSCAR. Defaults to unit
    cell formula of structure. Defaults to None.


    * **selective_dynamics** (*ArrayLike** | **None**, **optional*) – Bool values for selective dynamics,
    where N is the number of sites. Defaults to None.


    * **true_names** (*bool**, **optional*) – Set to False if the names in the POSCAR are not
    well-defined and ambiguous. This situation arises commonly in
    VASP < 5 where the POSCAR sometimes does not contain element
    symbols. Defaults to True.


    * **velocities** (*ArrayLike** | **None**, **optional*) – Velocities for the POSCAR. Typically parsed
    in MD runs or can be used to initialize velocities. Defaults to None.


    * **predictor_corrector** (*ArrayLike** | **None**, **optional*) – Predictor corrector for the POSCAR.
    Typically parsed in MD runs. Defaults to None.


    * **predictor_corrector_preamble** (*str** | **None**, **optional*) – Preamble to the predictor
    corrector. Defaults to None.


    * **sort_structure** (*bool**, **optional*) – Whether to sort the structure. Useful if species
    are not grouped properly together. Defaults to False.



#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d: dict)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Poscar



#### _static_ from_file(filename, check_for_POTCAR=True, read_velocities=True)
Reads a Poscar from a file.

The code will try its best to determine the elements in the POSCAR in
the following order:

1. If check_for_POTCAR is True, the code will try to check if a POTCAR
is in the same directory as the POSCAR and use elements from that by
default. (This is the VASP default sequence of priority).
2. If the input file is VASP5-like and contains element symbols in the
6th line, the code will use that if check_for_POTCAR is False or there
is no POTCAR found.
3. Failing (2), the code will check if a symbol is provided at the end
of each coordinate.

If all else fails, the code will just assign the first n elements in
increasing atomic number, where n is the number of species, to the
Poscar. For example, H, He, Li, …. This will ensure at least a
unique element is assigned to each site and any analysis that does not
require specific elemental properties should work fine.


* **Parameters**


    * **filename** (*str*) – File name containing Poscar data.


    * **check_for_POTCAR** (*bool*) – Whether to check if a POTCAR is present
    in the same directory as the POSCAR. Defaults to True.


    * **read_velocities** (*bool*) – Whether to read or not velocities if they
    are present in the POSCAR. Default is True.



* **Returns**

    Poscar object.



#### _static_ from_str(data, default_names=None, read_velocities=True)
Reads a Poscar from a string.

The code will try its best to determine the elements in the POSCAR in
the following order:

1. If default_names are supplied and valid, it will use those. Usually,
default names comes from an external source, such as a POTCAR in the
same directory.

2. If there are no valid default names but the input file is VASP5-like
and contains element symbols in the 6th line, the code will use that.

3. Failing (2), the code will check if a symbol is provided at the end
of each coordinate.

If all else fails, the code will just assign the first n elements in
increasing atomic number, where n is the number of species, to the
Poscar. For example, H, He, Li, …. This will ensure at least a
unique element is assigned to each site and any analysis that does not
require specific elemental properties should work fine.


* **Parameters**


    * **data** (*str*) – String containing Poscar data.


    * **default_names** (*[**str**]*) – Default symbols for the POSCAR file,
    usually coming from a POTCAR in the same directory.


    * **read_velocities** (*bool*) – Whether to read or not velocities if they
    are present in the POSCAR. Default is True.



* **Returns**

    Poscar object.



#### _classmethod_ from_string(\*args, \*\*kwargs)

#### get_str(direct: bool = True, vasp4_compatible: bool = False, significant_figures: int = 16)
Returns a string to be written as a POSCAR file. By default, site
symbols are written, which means compatibility is for vasp >= 5.


* **Parameters**


    * **direct** (*bool*) – Whether coordinates are output in direct or
    Cartesian. Defaults to True.


    * **vasp4_compatible** (*bool*) – Set to True to omit site symbols on 6th
    line to maintain backward vasp 4.x compatibility. Defaults
    to False.


    * **significant_figures** (*int*) – No. of significant figures to
    output all quantities. Defaults to 16. Note that positions are
    output in fixed point, while velocities are output in
    scientific format.



* **Returns**

    String representation of POSCAR.



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### _property_ natoms()
Sequence of number of sites of each type associated with the Poscar.
Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.


#### _property_ predictor_corrector()
Predictor corrector in Poscar.


#### _property_ predictor_corrector_preamble()
Predictor corrector preamble in Poscar.


#### _property_ selective_dynamics()
Selective dynamics in Poscar.


#### set_temperature(temperature: float)
Initializes the velocities based on Maxwell-Boltzmann distribution.
Removes linear, but not angular drift (same as VASP).

Scales the energies to the exact temperature (microcanonical ensemble)
Velocities are given in A/fs. This is the vasp default when
direct/cartesian is not specified (even when positions are given in
direct coordinates)

Overwrites imported velocities, if any.


* **Parameters**

    **temperature** (*float*) – Temperature in Kelvin.



#### _property_ site_symbols()
Sequence of symbols associated with the Poscar. Similar to 6th line in
vasp 5+ POSCAR.


#### _property_ velocities()
Velocities in Poscar.


#### write_file(filename: PathLike, \*\*kwargs)
Writes POSCAR to a file. The supported kwargs are the same as those for
the Poscar.get_string method and are passed through directly.


### _class_ Potcar(symbols=None, functional=None, sym_potcar_map=None)
Bases: `list`, `MSONable`

Object for reading and writing POTCAR files for calculations. Consists of a
list of PotcarSingle.


* **Parameters**


    * **symbols** (*[**str**]*) – Element symbols for POTCAR. This should correspond
    to the symbols used by VASP. E.g., “Mg”, “Fe_pv”, etc.


    * **functional** (*str*) – Functional used. To know what functional options
    there are, use Potcar.FUNCTIONAL_CHOICES. Note that VASP has
    different versions of the same functional. By default, the old
    PBE functional is used. If you want the newer ones, use PBE_52 or
    PBE_54. Note that if you intend to compare your results with the
    Materials Project, you should use the default setting. You can also
    override the default by setting PMG_DEFAULT_FUNCTIONAL in your
    .pmgrc.yaml.


    * **sym_potcar_map** (*dict*) – Allows a user to specify a specific element
    symbol to raw POTCAR mapping.



#### FUNCTIONAL_CHOICES(_ = ('PBE', 'PBE_52', 'PBE_54', 'LDA', 'LDA_52', 'LDA_54', 'PW91', 'LDA_US', 'PW91_US', 'Perdew_Zunger81'_ )

#### as_dict()
MSONable dict representation


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    Potcar



#### _static_ from_file(filename: str)
Reads Potcar from file.


* **Parameters**

    **filename** – Filename



* **Returns**

    Potcar



#### set_symbols(symbols, functional=None, sym_potcar_map=None)
Initialize the POTCAR from a set of symbols. Currently, the POTCARs can
be fetched from a location specified in .pmgrc.yaml. Use pmg config
to add this setting.


* **Parameters**


    * **symbols** (*[**str**]*) – A list of element symbols


    * **functional** (*str*) – The functional to use. If None, the setting
    PMG_DEFAULT_FUNCTIONAL in .pmgrc.yaml is used, or if this is
    not set, it will default to PBE.


    * **sym_potcar_map** (*dict*) – A map of symbol:raw POTCAR string. If
    sym_potcar_map is specified, POTCARs will be generated from
    the given map data rather than the config file location.



#### _property_ spec()
Get the atomic symbols and hash of all the atoms in the POTCAR file.


#### _property_ symbols()
Get the atomic symbols of all the atoms in the POTCAR file.


#### write_file(filename: str)
Write Potcar to a file.


* **Parameters**

    **filename** (*str*) – filename to write to.



### _class_ PotcarSingle(data, symbol=None)
Bases: `object`

Object for a **single** POTCAR. The builder assumes the POTCAR contains
the complete untouched data in “data” as a string and a dict of keywords.


#### data()
POTCAR data as a string.


* **Type**

    str



#### keywords()
Keywords parsed from the POTCAR as a dict. All keywords are also
accessible as attributes in themselves. E.g., potcar.enmax, potcar.encut, etc.


* **Type**

    dict


md5 hashes of the entire POTCAR file and the actual data are validated
against a database of known good hashes. Appropriate warnings or errors
are raised if a POTCAR hash fails validation.


* **Parameters**


    * **data** – Complete and single potcar file as a string.


    * **symbol** – POTCAR symbol corresponding to the filename suffix
    e.g. “Tm_3” for POTCAR.TM_3”. If not given, pymatgen
    will attempt to extract the symbol from the file itself.
    However, this is not always reliable!



#### _property_ atomic_no(_: in_ )
Attempt to return the atomic number based on the VRHFIN keyword.


#### _property_ electron_configuration()
Electronic configuration of the PotcarSingle.


#### _property_ element(_: st_ )
Attempt to return the atomic symbol based on the VRHFIN keyword.


#### _static_ from_file(filename: str)
Reads PotcarSingle from file.


* **Parameters**

    **filename** – Filename.



* **Returns**

    PotcarSingle.



#### _static_ from_symbol_and_functional(symbol: str, functional: str | None = None)
Makes a PotcarSingle from a symbol and functional.


* **Parameters**


    * **symbol** – Symbol, e.g., Li_sv


    * **functional** – E.g., PBE



* **Returns**

    PotcarSingle



#### _property_ functional(_: str | Non_ )
Functional associated with PotcarSingle.


#### _property_ functional_class()
Functional class associated with PotcarSingle.


#### functional_dir(_ = {'LDA': 'POT_LDA_PAW', 'LDA_52': 'POT_LDA_PAW_52', 'LDA_54': 'POT_LDA_PAW_54', 'LDA_US': 'POT_LDA_US', 'PBE': 'POT_GGA_PAW_PBE', 'PBE_52': 'POT_GGA_PAW_PBE_52', 'PBE_54': 'POT_GGA_PAW_PBE_54', 'PW91': 'POT_GGA_PAW_PW91', 'PW91_US': 'POT_GGA_US_PW91', 'Perdew_Zunger81': 'POT_LDA_PAW'_ )

#### functional_tags(_ = {'91': {'class': 'GGA', 'name': 'PW91'}, 'am': {'class': 'GGA', 'name': 'AM05'}, 'ca': {'class': 'LDA', 'name': 'Perdew-Zunger81'}, 'hl': {'class': 'LDA', 'name': 'Hedin-Lundquist'}, 'lm': {'class': 'GGA', 'name': 'Langreth-Mehl-Hu'}, 'pb': {'class': 'GGA', 'name': 'Perdew-Becke'}, 'pe': {'class': 'GGA', 'name': 'PBE'}, 'ps': {'class': 'GGA', 'name': 'PBEsol'}, 'pw': {'class': 'GGA', 'name': 'PW86'}, 'rp': {'class': 'GGA', 'name': 'revPBE'}, 'wi': {'class': 'LDA', 'name': 'Wigner Interpolation'}_ )

#### get_potcar_file_hash()
Computes a md5 hash of the entire PotcarSingle.

This hash corresponds to the md5 hash of the POTCAR file itself.


* **Returns**

    Hash value.



#### get_potcar_hash()
Computes a md5 hash of the metadata defining the PotcarSingle.


* **Returns**

    Hash value.



#### get_sha256_file_hash()
Computes a SHA256 hash of the PotcarSingle EXCLUDING lines starting with ‘SHA256’ and ‘CPRY’.

This hash corresponds to the sha256 hash printed in the header of modern POTCAR files.


* **Returns**

    Hash value.



#### identify_potcar(mode: Literal['data', 'file'] = 'data')
Identify the symbol and compatible functionals associated with this PotcarSingle.

This method checks the md5 hash of either the POTCAR metadadata (PotcarSingle.hash)
or the entire POTCAR file (PotcarSingle.file_hash) against a database
of hashes for POTCARs distributed with VASP 5.4.4.


* **Parameters**

    **mode** (*'data'** | **'file'*) – ‘data’ mode checks the hash of the POTCAR metadata in self.PSCTR,
    while ‘file’ mode checks the hash of the entire POTCAR file.



* **Returns**

    List of symbols associated with the PotcarSingle
    potcar_functionals (list): List of potcar functionals associated with

    > the PotcarSingle




* **Return type**

    symbol (list)



#### _property_ nelectrons(_: floa_ )
Number of electrons


#### parse_functions(_ = {'COPYR': <function _parse_string>, 'DEXC': <function _parse_float>, 'EATOM': <function _parse_float>, 'EAUG': <function _parse_float>, 'EMMIN': <function _parse_float>, 'ENMAX': <function _parse_float>, 'ENMIN': <function _parse_float>, 'GGA': <function _parse_list>, 'ICORE': <function _parse_int>, 'IUNSCR': <function _parse_int>, 'LCOR': <function _parse_bool>, 'LEXCH': <function _parse_string>, 'LPAW': <function _parse_bool>, 'LULTRA': <function _parse_bool>, 'LUNSCR': <function _parse_bool>, 'NDATA': <function _parse_int>, 'POMASS': <function _parse_float>, 'QCUT': <function _parse_float>, 'QGAM': <function _parse_float>, 'RAUG': <function _parse_float>, 'RCLOC': <function _parse_float>, 'RCORE': <function _parse_float>, 'RDEP': <function _parse_float>, 'RDEPT': <function _parse_float>, 'RMAX': <function _parse_float>, 'RPACOR': <function _parse_float>, 'RRKJ': <function _parse_list>, 'RWIGS': <function _parse_float>, 'SHA256': <function _parse_string>, 'STEP': <function _parse_list>, 'TITEL': <function _parse_string>, 'VRHFIN': <function _parse_string>, 'ZVAL': <function _parse_float>_ )

#### _property_ potential_type(_: Literal['NC', 'PAW', 'US'_ )
Type of PSP. E.g., US, PAW, etc.


#### _property_ symbol(_: st_ )
The POTCAR symbol, e.g. W_pv


#### verify_potcar()
Attempts to verify the integrity of the POTCAR data.

This method checks the whole file (removing only the SHA256
metadata) against the SHA256 hash in the header if this is found.
If no SHA256 hash is found in the file, the file hash (md5 hash of the
whole file) is checked against all POTCAR file hashes known to pymatgen.

### Returns:

(bool, bool)

    has_sh256 and passed_hash_check are returned.


#### write_file(filename: str)
Write PotcarSingle to a file.


* **Parameters**

    **filename** (*str*) – Filename to write to.



### _exception_ UnknownPotcarWarning()
Bases: `UserWarning`

Warning raised when POTCAR hashes do not pass validation.


### _class_ VaspInput(incar, kpoints, poscar, potcar, optional_files=None, \*\*kwargs)
Bases: `dict`, `MSONable`

Class to contain a set of vasp input objects corresponding to a run.

Initializes a VaspInput object with the given input files.


* **Parameters**


    * **incar** (*Incar*) – The Incar object.


    * **kpoints** (*Kpoints*) – The Kpoints object.


    * **poscar** (*Poscar*) – The Poscar object.


    * **potcar** (*Potcar*) – The Potcar object.


    * **optional_files** (*dict*) – Other input files supplied as a dict of {filename: object}.
    The object should follow standard pymatgen conventions in implementing a
    as_dict() and from_dict method.


    * **\*\*kwargs** – Additional keyword arguments to be stored in the VaspInput object.



#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    VaspInput



#### _static_ from_directory(input_dir, optional_files=None)
Read in a set of VASP input from a directory. Note that only the
standard INCAR, POSCAR, POTCAR and KPOINTS files are read unless
optional_filenames is specified.


* **Parameters**


    * **input_dir** (*str*) – Directory to read VASP input from.


    * **optional_files** (*dict*) – Optional files to read in as well as a
    dict of {filename: Object type}. Object type must have a
    static method from_file.



#### run_vasp(run_dir: PathLike = '.', vasp_cmd: list | None = None, output_file: PathLike = 'vasp.out', err_file: PathLike = 'vasp.err')
Write input files and run VASP.


* **Parameters**


    * **run_dir** – Where to write input files and do the run.


    * **vasp_cmd** – Args to be supplied to run VASP. Otherwise, the
    PMG_VASP_EXE in .pmgrc.yaml is used.


    * **output_file** – File to write output.


    * **err_file** – File to write err.



#### write_input(output_dir='.', make_dir_if_not_present=True)
Write VASP input to a directory.


* **Parameters**


    * **output_dir** (*str*) – Directory to write to. Defaults to current
    directory (“.”).


    * **make_dir_if_not_present** (*bool*) – Create the directory if not
    present. Defaults to True.



### _parse_bool(s)

### _parse_float(s)

### _parse_int(s)

### _parse_list(s)

### _parse_string(s)
## pymatgen.io.vasp.optics module

Classes for parsing and manipulating VASP optical properties calculations.


### _class_ DielectricFunctionCalculator(cder_real: NDArray, cder_imag: NDArray, eigs: NDArray, kweights: NDArray, nedos: int, deltae: float, ismear: int, sigma: float, efermi: float, cshift: float, ispin: int, volume: float)
Bases: `MSONable`

Class for postprocessing VASP optical properties calculations.

This objects helps load the different parameters from the vasprun.xml file but allows users to override
them as needed.

The standard vasprun.xml from an `LOPTICS=.True.` calculation already contains
the complex frequency dependent dielectric functions.  However you have no way to decompose
the different contributions.  Since the `WAVEDER` file is also written during an optical calculation,
you can reconstruct the dielectric functions purely in Python and have full control over contribution
from different bands and k-points.

VASP’s linear optics follow these steps:


    * Calculate the imaginary part


    * Perform symmetry operations (this is not implemented here)


    * Calculate the real part

Currently, this Calculator only works for `ISYM=0` calculations since we cannot guarantee that our
externally defined symmetry operations are the same as VASP’s. This can be fixed by printing the
symmetry operators into the vasprun.xml file. If this happens in future versions of VASP,
we can dramatically speed up the calculations here by considering only the irreducible kpoints.


#### _property_ cder()
Complex CDER from WAVEDER.


#### cder_imag(_: NDArra_ )

#### cder_real(_: NDArra_ )

#### cshift(_: floa_ )

#### deltae(_: floa_ )

#### efermi(_: floa_ )

#### eigs(_: NDArra_ )

#### _classmethod_ from_directory(directory: Path | str)
Construct a DielectricFunction from a directory containing vasprun.xml and WAVEDER files.


#### _classmethod_ from_vasp_objects(vrun: Vasprun, waveder: Waveder)
Construct a DielectricFunction from Vasprun, Kpoint, and Waveder objects.


* **Parameters**


    * **vrun** – Vasprun object


    * **kpoint** – Kpoint object


    * **waveder** – Waveder object



#### get_epsilon(idir: int, jdir: int, efermi: float | None = None, nedos: int | None = None, deltae: float | None = None, ismear: int | None = None, sigma: float | None = None, cshift: float | None = None, mask: NDArray | None = None)
Compute the frequency dependent dielectric function.


* **Parameters**


    * **idir** – First direction of the dielectric tensor


    * **jdir** – Second direction of the dielectric tensor


    * **efermi** – Fermi energy


    * **nedos** – Number of points in the DOS


    * **deltae** – Energy step in the DOS


    * **ismear** – Smearing method (only has 0:gaussian, >0:Methfessel-Paxton)


    * **sigma** – Smearing width


    * **cshift** – Complex shift used for Kramer-Kronig transformation


    * **mask** – Mask for the bands/kpoint/spin index to include in the calculation



#### ismear(_: in_ )

#### ispin(_: in_ )

#### kweights(_: NDArra_ )

#### nedos(_: in_ )

#### plot_weighted_transition_data(idir: int, jdir: int, mask: NDArray | None = None, min_val: float = 0.0)
Data for plotting the weight matrix elements as a scatter plot.

Since the computation of the final spectrum (especially the smearing part)
is still fairly expensive.  This function can be used to check the values
of some portion of the spectrum (defined by the mask).
In a sense, we are lookin at the imaginary part of the dielectric function
before the smearing is applied.


* **Parameters**


    * **idir** – First direction of the dielectric tensor.


    * **jdir** – Second direction of the dielectric tensor.


    * **mask** – Mask to apply to the CDER for the bands/kpoint/spin
    index to include in the calculation


    * **min_val** – Minimum value below this value the matrix element will not be shown.



#### sigma(_: floa_ )

#### volume(_: floa_ )

### delta_func(x, ismear)
Replication of VASP’s delta function.


### delta_methfessel_paxton(x, n)
D_n (x) = exp -x^2 \* sum_i=0^n A_i H_2i(x)
where H is a Hermite polynomial and
A_i = (-1)^i / ( i! 4^i sqrt(pi) ).


### epsilon_imag(cder: NDArray, eigs: NDArray, kweights: ArrayLike, efermi: float, nedos: int, deltae: float, ismear: int, sigma: float, idir: int, jdir: int, mask: NDArray | None = None)
Replicate the EPSILON_IMAG function of VASP.


* **Parameters**


    * **cder** – The data written to the WAVEDER (nbands, nbands, nkpoints, nspin, diri, dirj)


    * **eigs** – The eigenvalues (nbands, nkpoints, nspin)


    * **kweights** – The kpoint weights (nkpoints)


    * **efermi** – The fermi energy


    * **nedos** – The sampling of the energy values


    * **deltae** – The energy grid spacing


    * **ismear** – The smearing parameter used by the `step_func`.


    * **sigma** – The width of the smearing


    * **idir** – The first direction of the dielectric tensor


    * **jdir** – The second direction of the dielectric tensor


    * **mask** – Mask for the bands/kpoint/spin index to include in the calculation



* **Returns**

    Array of size nedos with the imaginary part of the dielectric function.



* **Return type**

    np.array



### get_delta(x0: float, sigma: float, nx: int, dx: float, ismear: int = 3)
Get the smeared delta function to be added to form the spectrum.

This replaces the SLOT function from VASP. Uses finite differences instead of
evaluating the delta function since the step function is more likely to have analytic form.


* **Parameters**


    * **x0** – The center of the dielectric function.


    * **sigma** – The width of the smearing


    * **nx** – The number of grid points in the output grid.


    * **dx** – The gridspacing of the output grid.


    * **ismear** – The smearing parameter used by the `step_func`.



* **Returns**

    Array of size nx with delta function on the desired outputgrid.



* **Return type**

    np.array



### get_step(x0, sigma, nx, dx, ismear)
Get the smeared step function to be added to form the spectrum.

This replaces the SLOT function from VASP.


* **Parameters**


    * **x0** – The center of the dielectric function.


    * **sigma** – The width of the smearing


    * **nx** – The number of grid points in the output grid.


    * **dx** – The gridspacing of the output grid.


    * **ismear** – The smearing parameter used by the `step_func`.



* **Returns**

    Array of size nx with step function on the desired outputgrid.



* **Return type**

    np.array



### kramers_kronig(eps: np.ndarray, nedos: int, deltae: float, cshift: float = 0.1)
Perform the Kramers-Kronig transformation.

Perform the Kramers-Kronig transformation exactly as VASP does it.
The input eps should be complex and the imaginary part of the dielectric function
should be stored as the real part of the complex input array.
The output should be the complex dielectric function.


* **Parameters**


    * **eps** – The dielectric function with the imaginary part stored as the real part and nothing in the imaginary part.


    * **nedos** – The sampling of the energy values


    * **deltae** – The energy grid spacing


    * **cshift** – The shift of the imaginary part of the dielectric function.



* **Returns**

    Array of size nedos with the complex dielectric function.



* **Return type**

    np.array



### step_func(x, ismear)
Replication of VASP’s step function.


### step_methfessel_paxton(x, n)
S_n (x) = (1 + erf x)/2 - exp -x^2 \* sum_i=1^n A_i H_{2i-1}(x)
where H is a Hermite polynomial and
A_i = (-1)^i / ( i! 4^i sqrt(pi) ).

## pymatgen.io.vasp.outputs module

Classes for reading/manipulating/writing VASP output files.


### _class_ BSVasprun(filename: str, parse_projected_eigen: bool | str = False, parse_potcar_file: bool | str = False, occu_tol: float = 1e-08, separate_spins: bool = False)
Bases: `Vasprun`

A highly optimized version of Vasprun that parses only eigenvalues for
bandstructures. All other properties like structures, parameters,
etc. are ignored.


* **Parameters**


    * **filename** – Filename to parse


    * **parse_projected_eigen** – Whether to parse the projected
    eigenvalues. Defaults to False. Set to True to obtain projected
    eigenvalues. **Note that this can take an extreme amount of time
    and memory.** So use this wisely.


    * **parse_potcar_file** – Whether to parse the potcar file to read
    the potcar hashes for the potcar_spec attribute. Defaults to True,
    where no hashes will be determined and the potcar_spec dictionaries
    will read {“symbol”: ElSymbol, “hash”: None}. By Default, looks in
    the same directory as the vasprun.xml, with same extensions as

    > Vasprun.xml. If a string is provided, looks at that filepath.



    * **occu_tol** – Sets the minimum tol for the determination of the
    vbm and cbm. Usually the default of 1e-8 works well enough,
    but there may be pathological cases.


    * **separate_spins** (*bool*) – Whether the band gap, CBM, and VBM should be
    reported for each individual spin channel. Defaults to False,
    which computes the eigenvalue band properties independent of
    the spin orientation. If True, the calculation must be spin-polarized.



#### as_dict()
JSON-serializable dict representation.


### _class_ Chgcar(poscar, data, data_aug=None)
Bases: `VolumetricData`

Simple object for reading a CHGCAR file.


* **Parameters**


    * **poscar** (*Poscar** or *[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Object containing structure.


    * **data** – Actual data.


    * **data_aug** – Augmentation charge data.



#### _static_ from_file(filename: str)
Read a CHGCAR file.


* **Parameters**

    **filename** (*str*) – Path to CHGCAR file.



* **Returns**

    Chgcar



#### _property_ net_magnetization()
Net magnetization from Chgcar


### _class_ Dynmat(filename)
Bases: `object`

Object for reading a DYNMAT file.


#### data()
A nested dict containing the DYNMAT data of the form:
[atom <int>][disp <int>][‘dispvec’] =

> displacement vector (part of first line in dynmat block, e.g. “0.01 0 0”)

[atom <int>][disp <int>][‘dynmat’] =

    <list> list of dynmat lines for this atom and this displacement


* **Type**

    dict


Authors: Patrick Huck


* **Parameters**

    **filename** – Name of file containing DYNMAT.



#### get_phonon_frequencies()
Calculate phonon frequencies.


#### _property_ masses()
Returns the list of atomic masses.


#### _property_ natoms()
Returns the number of atoms.


#### _property_ ndisps()
Returns the number of displacements.


#### _property_ nspecs()
Returns the number of species.


### _class_ Eigenval(filename, occu_tol=1e-08, separate_spins=False)
Bases: `object`

Object for reading EIGENVAL file.


#### filename()
String containing input filename.


* **Type**

    str



#### occu_tol()
Tolerance for determining occupation in band properties.


* **Type**

    float



#### ispin()
Spin polarization tag.


* **Type**

    int



#### nelect()
Number of electrons.


* **Type**

    int



#### nkpt()
Number of kpoints.


* **Type**

    int



#### nbands()
Number of bands.


* **Type**

    int



#### kpoints()
List of kpoints.


* **Type**

    list



#### kpoints_weights()
Weights of each kpoint in the BZ, should sum to 1.


* **Type**

    list



#### eigenvalues()
Eigenvalues as a dict of {(spin): np.ndarray(shape=(nkpt, nbands, 2))}.
This representation is based on actual ordering in VASP and is meant as an intermediate representation
to be converted into proper objects. The kpoint index is 0-based (unlike the 1-based indexing in VASP).


* **Type**

    dict


Reads input from filename to construct Eigenval object.


* **Parameters**


    * **filename** (*str*) – filename of EIGENVAL to read in


    * **occu_tol** (*float*) – tolerance for determining band gap


    * **separate_spins** (*bool*) – whether the band gap, CBM, and VBM should be
    reported for each individual spin channel. Defaults to False,
    which computes the eigenvalue band properties independent of
    the spin orientation. If True, the calculation must be spin-polarized.



* **Returns**

    a pymatgen.io.vasp.outputs.Eigenval object



#### _property_ eigenvalue_band_properties()
Band properties from the eigenvalues as a tuple,
(band gap, cbm, vbm, is_band_gap_direct). In the case of separate_spins=True,
the band gap, cbm, vbm, and is_band_gap_direct are each lists of length 2,
with index 0 representing the spin-up channel and index 1 representing
the spin-down channel.


### _class_ Elfcar(poscar, data)
Bases: `VolumetricData`

Read an ELFCAR file which contains the Electron Localization Function (ELF)
as calculated by VASP.

For ELF, “total” key refers to Spin.up, and “diff” refers to Spin.down.

This also contains information on the kinetic energy density.


* **Parameters**


    * **poscar** (*Poscar** or *[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Object containing structure.


    * **data** – Actual data.



#### _classmethod_ from_file(filename)
Reads a ELFCAR file.


* **Parameters**

    **filename** – Filename



* **Returns**

    Elfcar



#### get_alpha()
Get the parameter alpha where ELF = 1/(1+alpha^2).


### _class_ Locpot(poscar, data)
Bases: `VolumetricData`

Simple object for reading a LOCPOT file.


* **Parameters**


    * **poscar** (*Poscar*) – Poscar object containing structure.


    * **data** – Actual data.



#### _classmethod_ from_file(filename, \*\*kwargs)
Read a LOCPOT file.


* **Parameters**

    **filename** (*str*) – Path to LOCPOT file.



* **Returns**

    Locpot



### _class_ Oszicar(filename)
Bases: `object`

A basic parser for an OSZICAR output from VASP. In general, while the
OSZICAR is useful for a quick look at the output from a VASP run, we
recommend that you use the Vasprun parser instead, which gives far richer
information about a run.


#### electronic_steps()
All electronic steps as a list of list of dict. e.g.,
[[{“rms”: 160.0, “E”: 4507.24605593, “dE”: 4507.2, “N”: 1, “deps”: -17777.0, “ncg”: 16576}, …], [….]
where electronic_steps[index] refers the list of electronic steps in one ionic_step,
electronic_steps[index][subindex] refers to a particular electronic step at subindex in ionic step at
index. The dict of properties depends on the type of VASP run, but in general, “E”, “dE” and “rms” should
be present in almost all runs.


* **Type**

    list



#### ionic_steps()
All ionic_steps as a list of dict, e.g.,
[{“dE”: -526.36, “E0”: -526.36024, “mag”: 0.0, “F”: -526.36024}, …]
This is the typical output from VASP at the end of each ionic step. The stored dict might be different
depending on the type of VASP run.


* **Type**

    list



* **Parameters**

    **filename** (*str*) – Filename of file to parse.



#### _property_ all_energies()
Compilation of all energies from all electronic steps and ionic steps
as a tuple of list of energies, e.g.,
((4507.24605593, 143.824705755, -512.073149912, …), …).


#### as_dict()
MSONable dict


#### _property_ final_energy()

### _class_ Outcar(filename)
Bases: `object`

Parser for data in OUTCAR that is not available in Vasprun.xml.

Note, this class works a bit differently than most of the other
VaspObjects, since the OUTCAR can be very different depending on which
“type of run” performed.

Creating the OUTCAR class with a filename reads “regular parameters” that
are always present.


#### magnetization()
Magnetization on each ion as a tuple of dict, e.g.,
({“d”: 0.0, “p”: 0.003, “s”: 0.002, “tot”: 0.005}, … )


* **Type**

    tuple



#### chemical_shielding()
Chemical shielding on each ion as a dictionary with core and valence contributions.


* **Type**

    dict



#### unsym_cs_tensor()
Unsymmetrized chemical shielding tensor matrixes on each ion as a list.
e.g., [[[sigma11, sigma12, sigma13], [sigma21, sigma22, sigma23], [sigma31, sigma32, sigma33]], …]


* **Type**

    list



#### cs_g0_contribution()
G=0 contribution to chemical shielding. 2D rank 3 matrix.


* **Type**

    numpy.ndarray



#### cs_core_contribution()
Core contribution to chemical shielding. dict. e.g.,
{‘Mg’: -412.8, ‘C’: -200.5, ‘O’: -271.1}


* **Type**

    dict



#### efg()
Electric Field Gradient (EFG) tensor on each ion as a tuple of dict, e.g.,
({“cq”: 0.1, “eta”, 0.2, “nuclear_quadrupole_moment”: 0.3}, {“cq”: 0.7, “eta”, 0.8,
“nuclear_quadrupole_moment”: 0.9}, …)


* **Type**

    tuple



#### charge()
Charge on each ion as a tuple of dict, e.g.,
({“p”: 0.154, “s”: 0.078, “d”: 0.0, “tot”: 0.232}, …)


* **Type**

    tuple



#### is_stopped()
True if OUTCAR is from a stopped run (using STOPCAR, see VASP Manual).


* **Type**

    bool



#### run_stats()
Various useful run stats as a dict including “System time (sec)”, “Total CPU time used (sec)”,
“Elapsed time (sec)”, “Maximum memory used (kb)”, “Average memory used (kb)”, “User time (sec)”, “cores”.


* **Type**

    dict



#### elastic_tensor()
Total elastic moduli (Kbar) is given in a 6x6 array matrix.


* **Type**

    numpy.ndarray



#### drift()
Total drift for each step in eV/Atom.


* **Type**

    numpy.ndarray



#### ngf()
Dimensions for the Augmentation grid.


* **Type**

    tuple



#### sampling_radii()
Size of the sampling radii in VASP for the test charges for the electrostatic
potential at each atom. Total array size is the number of elements present in the calculation.


* **Type**

    numpy.ndarray



#### electrostatic_potential()
Average electrostatic potential at each atomic position in order of
the atoms in POSCAR.


* **Type**

    numpy.ndarray



#### final_energy_contribs()
Individual contributions to the total final energy as a dictionary.
Include contributions from keys, e.g.:
{‘DENC’: -505778.5184347, ‘EATOM’: 15561.06492564, ‘EBANDS’: -804.53201231, ‘EENTRO’: -0.08932659,
‘EXHF’: 0.0, ‘Ediel_sol’: 0.0, ‘PAW double counting’: 664.6726974100002, ‘PSCENC’: 742.48691646,
‘TEWEN’: 489742.86847338, ‘XCENC’: -169.64189814}


* **Type**

    dict



#### efermi()
Fermi energy.


* **Type**

    float



#### filename()
Filename.


* **Type**

    str



#### final_energy()
Final energy after extrapolation of sigma back to 0, i.e. energy(sigma->0).


* **Type**

    float



#### final_energy_wo_entrp()
Final energy before extrapolation of sigma, i.e. energy without entropy.


* **Type**

    float



#### final_fr_energy()
Final “free energy”, i.e. free energy TOTEN.


* **Type**

    float



#### has_onsite_density_matrices()
Boolean for if onsite density matrices have been set.


* **Type**

    bool



#### lcalcpol()
If LCALCPOL has been set.


* **Type**

    bool



#### lepsilon()
If LEPSILON has been set.


* **Type**

    bool



#### nelect()
Returns the number of electrons in the calculation.


* **Type**

    float



#### spin()
If spin-polarization was enabled via ISPIN.


* **Type**

    bool



#### total_mag()
Total magnetization (in terms of the number of unpaired electrons).


* **Type**

    float


One can then call a specific reader depending on the type of run being
performed. These are currently: read_igpar(), read_lepsilon() and
read_lcalcpol(), read_core_state_eign(), read_avg_core_pot().

See the documentation of those methods for more documentation.

Authors: Rickard Armiento, Shyue Ping Ong


* **Parameters**

    **filename** (*str*) – OUTCAR filename to parse.



#### _static_ _parse_sci_notation(line)
Method to parse lines with values in scientific notation and potentially
without spaces in between the values. This assumes that the scientific
notation always lists two digits for the exponent, e.g. 3.535E-02


* **Parameters**

    **line** – line to parse.



* **Returns**

    numbers if found, empty ist if not



* **Return type**

    list[float]



#### as_dict()
MSONable dict.


#### read_avg_core_poten()
Read the core potential at each ionic step.


* **Returns**

    A list for each ionic step containing a list of the average core
    potentials for each atom: [[avg core pot]].


### Example

The average core potential of the 2nd atom of the structure at the
last ionic step is: [-1][1]


#### read_chemical_shielding()
Parse the NMR chemical shieldings data. Only the second part “absolute, valence and core”
will be parsed. And only the three right most field (ISO_SHIELDING, SPAN, SKEW) will be retrieved.


* **Returns**

    List of chemical shieldings in the order of atoms from the OUTCAR. Maryland notation is adopted.



#### read_core_state_eigen()
Read the core state eigenenergies at each ionic step.


* **Returns**

    [core state eig]}].
    The core state eigenenergie list for each AO is over all ionic
    step.



* **Return type**

    A list of dict over the atom such as [{“AO”


### Example

The core state eigenenergie of the 2s AO of the 6th atom of the
structure at the last ionic step is [5][“2s”][-1]


#### read_corrections(reverse=True, terminate_on_match=True)
Reads the dipol qudropol corrections into the
Outcar.data[“dipol_quadrupol_correction”].


* **Parameters**


    * **reverse** (*bool*) – Whether to start from end of OUTCAR. Defaults to True.


    * **terminate_on_match** (*bool*) – Whether to terminate once match is found. Defaults to True.



#### read_cs_core_contribution()
Parse the core contribution of NMR chemical shielding.


* **Returns**

    G0 contribution matrix.



* **Return type**

    list[list]



#### read_cs_g0_contribution()
Parse the  G0 contribution of NMR chemical shielding.


* **Returns**

    G0 contribution matrix as list of list.



#### read_cs_raw_symmetrized_tensors()
Parse the matrix form of NMR tensor before corrected to table.


* **Returns**

    nsymmetrized tensors list in the order of atoms.



#### read_elastic_tensor()
Parse the elastic tensor data.


* **Returns**

    6x6 array corresponding to the elastic tensor from the OUTCAR.



#### read_electrostatic_potential()
Parses the eletrostatic potential for the last ionic step.


#### read_fermi_contact_shift()
Output example:
Fermi contact (isotropic) hyperfine coupling parameter (MHz)
————————————————————-
ion      A_pw      A_1PS     A_1AE     A_1c      A_tot
————————————————————-

> 1      -0.002    -0.002    -0.051     0.000    -0.052
> 2      -0.002    -0.002    -0.051     0.000    -0.052
> 3       0.056     0.056     0.321    -0.048     0.321

> [-0.002, -0.002, -0.051, 0.0, -0.052],
> [0.056, 0.056, 0.321, -0.048, 0.321]] from ‘fch’ data.


#### read_freq_dielectric()
Parses the frequency dependent dielectric function (obtained with
LOPTICS). Frequencies (in eV) are in self.frequencies, and dielectric
tensor function is given as self.dielectric_tensor_function.


#### read_igpar()
Renders accessible:

    er_ev = e<r>_ev (dictionary with Spin.up/Spin.down as keys)
    er_bp = e<r>_bp (dictionary with Spin.up/Spin.down as keys)
    er_ev_tot = spin up + spin down summed
    er_bp_tot = spin up + spin down summed
    p_elc = spin up + spin down summed
    p_ion = spin up + spin down summed.

(See VASP section “LBERRY,  IGPAR,  NPPSTR,  DIPOL tags” for info on
what these are).


#### read_internal_strain_tensor()
Reads the internal strain tensor and populates self.internal_strain_tensor with an array of voigt notation

    tensors for each site.


#### read_lcalcpol()
Reads the lcalpol.

# TODO: Document the actual variables.


#### read_lepsilon()
Reads an LEPSILON run.

# TODO: Document the actual variables.


#### read_lepsilon_ionic()
Reads an LEPSILON run, the ionic component.

# TODO: Document the actual variables.


#### read_neb(reverse=True, terminate_on_match=True)
Reads NEB data. This only works with OUTCARs from both normal
VASP NEB calculations or from the CI NEB method implemented by
Henkelman et al.


* **Parameters**


    * **reverse** (*bool*) – Read files in reverse. Defaults to false. Useful for
    large files, esp OUTCARs, especially when used with
    terminate_on_match. Defaults to True here since we usually
    want only the final value.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern. Defaults to True here
    since we usually want only the final value.


Renders accessible:

    tangent_force - Final tangent force.
    energy - Final energy.
    These can be accessed under Outcar.data[key]


#### read_nmr_efg()
Parse the NMR Electric Field Gradient interpreted values.


* **Returns**

    Electric Field Gradient tensors as a list of dict in the order of atoms from OUTCAR.
    Each dict key/value pair corresponds to a component of the tensors.



#### read_nmr_efg_tensor()
Parses the NMR Electric Field Gradient Raw Tensors.


* **Returns**

    A list of Electric Field Gradient Tensors in the order of Atoms from OUTCAR



#### read_onsite_density_matrices()
Parse the onsite density matrices, returns list with index corresponding
to atom index in Structure.


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
    results from regex and postprocess. Note that the returned values
    are lists of lists, because you can grep multiple items on one line.


#### read_piezo_tensor()
Parse the piezo tensor data.


#### read_pseudo_zval()
Create pseudopotential ZVAL dictionary.


#### read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=<class 'str'>, attribute_name=None, last_one_only=True, first_one_only=False)
Parse table-like data. A table composes of three parts: header,
main body, footer. All the data matches “row pattern” in the main body
will be returned.


* **Parameters**


    * **header_pattern** (*str*) – The regular expression pattern matches the
    table header. This pattern should match all the text
    immediately before the main body of the table. For multiple
    sections table match the text until the section of
    interest. MULTILINE and DOTALL options are enforced, as a
    result, the “.” meta-character will also match “n” in this
    section.


    * **row_pattern** (*str*) – The regular expression matches a single line in
    the table. Capture interested field using regular expression
    groups.


    * **footer_pattern** (*str*) – The regular expression matches the end of the
    table. E.g. a long dash line.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


    * **attribute_name** (*str*) – Name of this table. If present the parsed data
    will be attached to “data. e.g. self.data[“efg”] = […]


    * **last_one_only** (*bool*) – All the tables will be parsed, if this option
    is set to True, only the last table will be returned. The
    enclosing list will be removed. i.e. Only a single table will
    be returned. Default to be True. Incompatible with first_one_only.


    * **first_one_only** (*bool*) – Only the first occurrence of the table will be
    parsed and the parsing procedure will stop. The enclosing list
    will be removed. i.e. Only a single table will be returned.
    Incompatible with last_one_only.



* **Returns**

    List of tables. 1) A table is a list of rows. 2) A row if either a list of
    attribute values in case the capturing group is defined without name in
    row_pattern, or a dict in case that named capturing groups are defined by
    row_pattern.



### _class_ Procar(filename)
Bases: `object`

Object for reading a PROCAR file.


#### data()
The PROCAR data of the form below. It should VASP uses 1-based indexing,
but all indices are converted to 0-based here.
{ spin: nd.array accessed with (k-point index, band index, ion index, orbital index) }


* **Type**

    dict



#### weights()
The weights associated with each k-point as an nd.array of length nkpoints.


* **Type**

    numpy.ndarray



#### phase_factors()
Phase factors, where present (e.g. LORBIT = 12). A dict of the form:
{ spin: complex nd.array accessed with (k-point index, band index, ion index, orbital index) }


* **Type**

    dict



#### nbands()
Number of bands.


* **Type**

    int



#### nkpoints()
Number of k-points.


* **Type**

    int



#### nions()
Number of ions.


* **Type**

    int



* **Parameters**

    **filename** – Name of file containing PROCAR.



#### get_occupation(atom_index, orbital)
Returns the occupation for a particular orbital of a particular atom.


* **Parameters**


    * **atom_num** (*int*) – Index of atom in the PROCAR. It should be noted
    that VASP uses 1-based indexing for atoms, but this is
    converted to 0-based indexing in this parser to be
    consistent with representation of structures in pymatgen.


    * **orbital** (*str*) – An orbital. If it is a single character, e.g., s,
    p, d or f, the sum of all s-type, p-type, d-type or f-type
    orbitals occupations are returned respectively. If it is a
    specific orbital, e.g., px, dxy, etc., only the occupation
    of that orbital is returned.



* **Returns**

    Sum occupation of orbital of atom.



#### get_projection_on_elements(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Method returning a dictionary of projections on elements.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.



* **Returns**

    [k index][b index][{Element:values}]]



* **Return type**

    a dictionary in the {Spin.up



### _exception_ UnconvergedVASPWarning()
Bases: `Warning`

Warning for unconverged vasp run.


### _exception_ VaspParseError()
Bases: [`ParseError`](pymatgen.io.md#pymatgen.io.core.ParseError)

Exception class for VASP parsing.


### _class_ Vasprun(filename, ionic_step_skip=None, ionic_step_offset=0, parse_dos=True, parse_eigen=True, parse_projected_eigen=False, parse_potcar_file=True, occu_tol=1e-08, separate_spins=False, exception_on_bad_xml=True)
Bases: `MSONable`

Vastly improved cElementTree-based parser for vasprun.xml files. Uses
iterparse to support incremental parsing of large files.
Speedup over Dom is at least 2x for smallish files (~1Mb) to orders of
magnitude for larger files (~10Mb).

**VASP results**


#### ionic_steps()
All ionic steps in the run as a list of {“structure”: structure at end of run,
“electronic_steps”: {All electronic step data in vasprun file}, “stresses”: stress matrix}.


* **Type**

    list



#### tdos()
Total dos calculated at the end of run.


* **Type**

    [Dos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.Dos)



#### idos()
Integrated dos calculated at the end of run.


* **Type**

    [Dos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.Dos)



#### pdos()
List of list of PDos objects. Access as pdos[atomindex][orbitalindex].


* **Type**

    list



#### efermi()
Fermi energy.


* **Type**

    float



#### eigenvalues()
Final eigenvalues as a dict of {(spin, kpoint index):[[eigenvalue, occu]]}.
The kpoint index is 0-based (unlike the 1-based indexing in VASP).


* **Type**

    dict



#### projected_eigenvalues()
Final projected eigenvalues as a dict of {spin: nd-array}.
To access a particular value, you need to do
Vasprun.projected_eigenvalues[spin][kpoint index][band index][atom index][orbital_index].
The kpoint, band and atom indices are 0-based (unlike the 1-based indexing in VASP).


* **Type**

    dict



#### projected_magnetisation()
Final projected magnetization as a numpy array with the
shape (nkpoints, nbands, natoms, norbitals, 3). Where the last axis is the contribution in the
3 Cartesian directions. This attribute is only set if spin-orbit coupling (LSORBIT = True) or
non-collinear magnetism (LNONCOLLINEAR = True) is turned on in the INCAR.


* **Type**

    numpy.ndarray



#### other_dielectric()
Dictionary, with the tag comment as key, containing other variants of
the real and imaginary part of the dielectric constant (e.g., computed by RPA) in function of
the energy (frequency). Optical properties (e.g. absorption coefficient) can be obtained through this.
The data is given as a tuple of 3 values containing each of them the energy, the real part tensor,
and the imaginary part tensor ([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,imag_partxy, imag_partyz, imag_partxz]]).


* **Type**

    dict



#### nionic_steps()
The total number of ionic steps. This number is always equal to the total number
of steps in the actual run even if ionic_step_skip is used.


* **Type**

    int



#### force_constants()
Force constants computed in phonon DFPT run(IBRION = 8).
The data is a 4D numpy array of shape (natoms, natoms, 3, 3).


* **Type**

    numpy.ndarray



#### normalmode_eigenvals()
Normal mode frequencies. 1D numpy array of size 3\*natoms.


* **Type**

    numpy.ndarray



#### normalmode_eigenvecs()
Normal mode eigen vectors. 3D numpy array of shape (3\*natoms, natoms, 3).


* **Type**

    numpy.ndarray



#### md_data()
Available only for ML MD runs, i.e., INCAR with ML_LMLFF = .TRUE. md_data is a list of
dict with the following format: [{‘energy’: {‘e_0_energy’: -525.07195568, ‘e_fr_energy’: -525.07195568,
‘e_wo_entrp’: -525.07195568, ‘kinetic’: 3.17809233, ‘lattice kinetic’: 0.0, ‘nosekinetic’: 1.323e-5,
‘nosepot’: 0.0, ‘total’: -521.89385012}, ‘forces’: [[0.17677989, 0.48309874, 1.85806696], …],
‘structure’: Structure object}].


* **Type**

    list



#### incar()
Incar object for parameters specified in INCAR file.


* **Type**

    Incar



#### parameters()
Incar object with parameters that vasp actually used, including all defaults.


* **Type**

    Incar



#### kpoints()
Kpoints object for KPOINTS specified in run.


* **Type**

    [Kpoints](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints)



#### actual_kpoints()
List of actual kpoints, e.g., [[0.25, 0.125, 0.08333333], [-0.25, 0.125, 0.08333333],
[0.25, 0.375, 0.08333333], ….].


* **Type**

    list



#### actual_kpoints_weights()
List of kpoint weights, E.g., [0.04166667, 0.04166667, 0.04166667, 0.04166667,
0.04166667, ….].


* **Type**

    list



#### atomic_symbols()
List of atomic symbols, e.g., [“Li”, “Fe”, “Fe”, “P”, “P”, “P”].


* **Type**

    list



#### potcar_symbols()
List of POTCAR symbols. e.g., [“PAW_PBE Li 17Jan2003”, “PAW_PBE Fe 06Sep2000”, ..].


* **Type**

    list


Author: Shyue Ping Ong


* **Parameters**


    * **filename** (*str*) – Filename to parse


    * **ionic_step_skip** (*int*) – If ionic_step_skip is a number > 1,
    only every ionic_step_skip ionic steps will be read for
    structure and energies. This is very useful if you are parsing
    very large vasprun.xml files and you are not interested in every
    single ionic step. Note that the final energies may not be the
    actual final energy in the vasprun.


    * **ionic_step_offset** (*int*) – Used together with ionic_step_skip. If set,
    the first ionic step read will be offset by the amount of
    ionic_step_offset. For example, if you want to start reading
    every 10th structure but only from the 3rd structure onwards,
    set ionic_step_skip to 10 and ionic_step_offset to 3. Main use
    case is when doing statistical structure analysis with
    extremely long time scale multiple VASP calculations of
    varying numbers of steps.


    * **parse_dos** (*bool*) – Whether to parse the dos. Defaults to True. Set
    to False to shave off significant time from the parsing if you
    are not interested in getting those data.


    * **parse_eigen** (*bool*) – Whether to parse the eigenvalues. Defaults to
    True. Set to False to shave off significant time from the
    parsing if you are not interested in getting those data.


    * **parse_projected_eigen** (*bool*) – Whether to parse the projected
    eigenvalues and magnetization. Defaults to False. Set to True to obtain
    projected eigenvalues and magnetization. **Note that this can take an
    extreme amount of time and memory.** So use this wisely.


    * **parse_potcar_file** (*bool/str*) – Whether to parse the potcar file to read
    the potcar hashes for the potcar_spec attribute. Defaults to True,
    where no hashes will be determined and the potcar_spec dictionaries
    will read {“symbol”: ElSymbol, “hash”: None}. By Default, looks in
    the same directory as the vasprun.xml, with same extensions as

    > Vasprun.xml. If a string is provided, looks at that filepath.



    * **occu_tol** (*float*) – Sets the minimum tol for the determination of the
    vbm and cbm. Usually the default of 1e-8 works well enough,
    but there may be pathological cases.


    * **separate_spins** (*bool*) – Whether the band gap, CBM, and VBM should be
    reported for each individual spin channel. Defaults to False,
    which computes the eigenvalue band properties independent of
    the spin orientation. If True, the calculation must be spin-polarized.


    * **exception_on_bad_xml** (*bool*) – Whether to throw a ParseException if a
    malformed XML is detected. Default to True, which ensures only
    proper vasprun.xml are parsed. You can set to False if you want
    partial results (e.g., if you are monitoring a calculation during a
    run), but use the results with care. A warning is issued.



#### _parse(stream, parse_dos, parse_eigen, parse_projected_eigen)

#### _static_ _parse_atominfo(elem)

#### _parse_calculation(elem)

#### _parse_chemical_shielding_calculation(elem)

#### _static_ _parse_diel(elem)

#### _static_ _parse_dos(elem)

#### _static_ _parse_dynmat(elem)

#### _static_ _parse_eigen(elem)

#### _static_ _parse_kpoints(elem)

#### _static_ _parse_optical_transition(elem)

#### _parse_params(elem)

#### _static_ _parse_projected_eigen(elem)

#### _parse_structure(elem)

#### as_dict()
JSON-serializable dict representation.


#### calculate_efermi(tol: float = 0.001)
Calculate the Fermi level using a robust algorithm.

Sometimes VASP can put the Fermi level just inside of a band due to issues in
the way band occupancies are handled. This algorithm tries to detect and correct
for this bug.

Slightly more details are provided here: [https://www.vasp.at/forum/viewtopic.php?f=4&t=17981](https://www.vasp.at/forum/viewtopic.php?f=4&t=17981)


#### _property_ complete_dos()
A complete dos object which incorporates the total dos and all
projected dos.


#### _property_ complete_dos_normalized(_: [CompleteDos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.CompleteDos_ )
A CompleteDos object which incorporates the total DOS and all
projected DOS. Normalized by the volume of the unit cell with
units of states/eV/unit cell volume.


#### _property_ converged()
Returns:
bool: True if a relaxation run is both ionically and electronically converged.


#### _property_ converged_electronic()
Returns:
bool: True if electronic step convergence has been reached in the final ionic step.


#### _property_ converged_ionic()
Returns:
bool: True if ionic step convergence has been reached, i.e. that vasp

> exited before reaching the max ionic steps for a relaxation run.


#### _property_ dielectric()
Returns:
The real and imaginary part of the dielectric constant (e.g., computed
by RPA) in function of the energy (frequency). Optical properties (e.g.
absorption coefficient) can be obtained through this.
The data is given as a tuple of 3 values containing each of them
the energy, the real part tensor, and the imaginary part tensor
([energies],[[real_partxx,real_partyy,real_partzz,real_partxy,
real_partyz,real_partxz]],[[imag_partxx,imag_partyy,imag_partzz,
imag_partxy, imag_partyz, imag_partxz]]).


#### _property_ eigenvalue_band_properties()
Band properties from the eigenvalues as a tuple,
(band gap, cbm, vbm, is_band_gap_direct). In the case of separate_spins=True,
the band gap, cbm, vbm, and is_band_gap_direct are each lists of length 2,
with index 0 representing the spin-up channel and index 1 representing
the spin-down channel.


#### _property_ epsilon_ionic()
Property only available for DFPT calculations and when IBRION=5, 6, 7 or 8.


* **Returns**

    The ionic part of the static dielectric constant. Present when it’s a
    DFPT run (LEPSILON=TRUE) and IBRION=5, 6, 7 or 8



#### _property_ epsilon_static()
Property only available for DFPT calculations.


* **Returns**

    The static part of the dielectric constant. Present when it’s a DFPT run
    (LEPSILON=TRUE)



#### _property_ epsilon_static_wolfe()
Property only available for DFPT calculations.


* **Returns**

    The static part of the dielectric constant without any local field
    effects. Present when it’s a DFPT run (LEPSILON=TRUE)



#### _property_ final_energy()

#### get_band_structure(kpoints_filename: str | None = None, efermi: float | Literal['smart'] | None = None, line_mode: bool = False, force_hybrid_mode: bool = False)
Get the band structure as a BandStructure object.


* **Parameters**


    * **kpoints_filename** – Full path of the KPOINTS file from which
    the band structure is generated.
    If none is provided, the code will try to intelligently
    determine the appropriate KPOINTS file by substituting the
    filename of the vasprun.xml with KPOINTS.
    The latter is the default behavior.


    * **efermi** – The Fermi energy associated with the bandstructure, in eV. By
    default (None), uses the value reported by VASP in vasprun.xml. To
    manually set the Fermi energy, pass a float. Pass ‘smart’ to use the
    calculate_efermi() method, which calculates the Fermi level by first
    checking whether it lies within a small tolerance (by default 0.001 eV)
    of a band edge) If it does, the Fermi level is placed in the center of
    the bandgap. Otherwise, the value is identical to the value reported by
    VASP.


    * **line_mode** – Force the band structure to be considered as
    a run along symmetry lines. (Default: False)


    * **force_hybrid_mode** – Makes it possible to read in self-consistent band
    structure calculations for every type of functional. (Default: False)



* **Returns**

    a BandStructure object (or more specifically a
    BandStructureSymmLine object if the run is detected to be a run
    along symmetry lines)

    Two types of runs along symmetry lines are accepted: non-sc with
    Line-Mode in the KPOINT file or hybrid, self-consistent with a
    uniform grid+a few kpoints along symmetry lines (explicit KPOINTS
    file) (it’s not possible to run a non-sc band structure with hybrid
    functionals). The explicit KPOINTS file needs to have data on the
    kpoint label as commentary.




#### get_computed_entry(inc_structure=True, parameters=None, data=None, entry_id: str | None = None)
Returns a ComputedEntry or ComputedStructureEntry from the Vasprun.


* **Parameters**


    * **inc_structure** (*bool*) – Set to True if you want
    ComputedStructureEntries to be returned instead of
    ComputedEntries.


    * **parameters** (*list*) – Input parameters to include. It has to be one of
    the properties supported by the Vasprun object. If
    parameters is None, a default set of parameters that are
    necessary for typical post-processing will be set.


    * **data** (*list*) – Output data to include. Has to be one of the properties
    supported by the Vasprun object.


    * **entry_id** (*str*) – Specify an entry id for the ComputedEntry. Defaults to
    “vasprun-{current datetime}”



* **Returns**

    ComputedStructureEntry/ComputedEntry



#### get_potcars(path: str | Path)
Returns the POTCAR from the specified path.


* **Parameters**

    **path** (*str*) – The path to search for POTCARs.



* **Returns**

    The POTCAR from the specified path.



* **Return type**

    Potcar | None



#### get_trajectory()
This method returns a Trajectory object, which is an alternative
representation of self.structures into a single object. Forces are
added to the Trajectory as site properties.


* **Returns**

    from pymatgen.core.trajectory



* **Return type**

    [Trajectory](pymatgen.core.md#pymatgen.core.trajectory.Trajectory)



#### _property_ hubbards()
Hubbard U values used if a vasprun is a GGA+U run. {} otherwise.


#### _property_ is_hubbard(_: boo_ )
True if run is a DFT+U run.


#### _property_ is_spin(_: boo_ )
True if run is spin-polarized.


#### _property_ optical_absorption_coeff()
Calculate the optical absorption coefficient
from the dielectric constants. Note that this method is only
implemented for optical properties calculated with GGA and BSE.


* **Returns**

    optical absorption coefficient in list



#### _property_ run_type()
Returns the run type. Currently detects GGA, metaGGA, HF, HSE, B3LYP,
and hybrid functionals based on relevant INCAR tags. LDA is assigned if
PAW POTCARs are used and no other functional is detected.

Hubbard U terms and vdW corrections are detected automatically as well.


#### _property_ structures()
Returns:
List of Structure objects for the structure at each ionic step.


#### update_charge_from_potcar(path)
Sets the charge of a structure based on the POTCARs found.


* **Parameters**

    **path** – Path to search for POTCARs



#### update_potcar_spec(path)

* **Parameters**

    **path** – Path to search for POTCARs



* **Returns**

    Potcar spec from path.



### _class_ VolumetricData(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), data, distance_matrix=None, data_aug=None)
Bases: [`VolumetricData`](pymatgen.io.md#pymatgen.io.common.VolumetricData)

Container for volumetric data that allows
for reading/writing with Poscar-type data.

Typically, this constructor is not used directly and the static
from_file constructor is used. This constructor is designed to allow
summation and other operations between VolumetricData objects.


* **Parameters**


    * **structure** – Structure associated with the volumetric data


    * **data** – Actual volumetric data. If the data is provided as in list format,
    it will be converted into an np.array automatically


    * **data_aug** – Any extra information associated with volumetric data
    (typically augmentation charges)


    * **distance_matrix** – A pre-computed distance matrix if available.
    Useful so pass distance_matrices between sums,
    short-circuiting an otherwise expensive operation.



#### _static_ parse_file(filename)
Convenience method to parse a generic volumetric data file in the vasp
like format. Used by subclasses for parsing file.


* **Parameters**

    **filename** (*str*) – Path of file to parse



* **Returns**

    (poscar, data)



#### write_file(file_name, vasp4_compatible=False)
Write the VolumetricData object to a vasp compatible file.


* **Parameters**


    * **file_name** (*str*) – Path to a file


    * **vasp4_compatible** (*bool*) – True if the format is vasp4 compatible



### _class_ WSWQ(nspin: int, nkpoints: int, nbands: int, me_real: ndarray, me_imag: ndarray)
Bases: `MSONable`

Class for reading a WSWQ file.
The WSWQ file is used to calculation the wave function overlaps between

>
> * W: Wavefunctions in the currenct directory’s WAVECAR file


> * WQ: Wavefunctions stored in a filed named the WAVECAR.qqq.

The overlap is computed using the overlap operator S
which make the PAW wavefunctions orthogonormal:

> <W_k,m| S | W_k,n> = delta_{mn}

The WSWQ file contains matrix elements of the overlap operator S evaluated
between the planewave wavefunctions W and WQ:

> COVL_k,mn = < W_s,k,m | S | WQ_s,k,n >

The indices of WSWQ.data are:

    [spin][kpoint][band_i][band_j]


#### nspin()
Number of spin channels


* **Type**

    int



#### nkpoints()
Number of k-points


* **Type**

    int



#### nbands()
Number of bands


* **Type**

    int



#### me_real()
Real part of the overlap matrix elements


* **Type**

    numpy.ndarray



#### me_imag()
Imaginary part of the overlap matrix elements


* **Type**

    numpy.ndarray



#### _property_ data()
Complex overlap matrix.


#### _classmethod_ from_file(filename: str)
Constructs a WSWQ object from a file.


* **Parameters**

    **filename** (*str*) – Name of WSWQ file.



* **Returns**

    WSWQ object.



#### me_imag(_: ndarra_ )

#### me_real(_: ndarra_ )

#### nbands(_: in_ )

#### nkpoints(_: in_ )

#### nspin(_: in_ )

### _class_ Wavecar(filename='WAVECAR', verbose=False, precision='normal', vasp_type=None)
Bases: `object`

This is a class that contains the (pseudo-) wavefunctions from VASP.

Coefficients are read from the given WAVECAR file and the corresponding
G-vectors are generated using the algorithm developed in WaveTrans (see
acknowledgments below). To understand how the wavefunctions are evaluated,
please see the evaluate_wavefunc docstring.

It should be noted that the pseudopotential augmentation is not included in
the WAVECAR file. As a result, some caution should be exercised when
deriving value from this information.

The usefulness of this class is to allow the user to do projections or band
unfolding style manipulations of the wavefunction. An example of this can
be seen in the work of Shen et al. 2017
([https://doi.org/10.1103/PhysRevMaterials.1.065001](https://doi.org/10.1103/PhysRevMaterials.1.065001)).


#### filename()
String of the input file (usually WAVECAR).


* **Type**

    str



#### vasp_type()
String that determines VASP type the WAVECAR was generated with.
One of ‘std’, ‘gam’, ‘ncl’.


* **Type**

    str



#### nk()
Number of k-points from the WAVECAR.


* **Type**

    int



#### nb()
Number of bands per k-point.


* **Type**

    int



#### encut()
Energy cutoff (used to define G_{cut}).


* **Type**

    float



#### efermi()
Fermi energy.


* **Type**

    float



#### a()
Primitive lattice vectors of the cell (e.g. a_1 = self.a[0, :]).


* **Type**

    numpy.ndarray



#### b()
Reciprocal lattice vectors of the cell (e.g. b_1 = self.b[0, :]).


* **Type**

    numpy.ndarray



#### vol()
The volume of the unit cell in real space.


* **Type**

    float



#### kpoints()
The list of k-points read from the WAVECAR file.


* **Type**

    numpy.ndarray



#### band_energy()
The list of band eigenenergies (and corresponding occupancies) for each kpoint,
where the first index corresponds to the index of the k-point (e.g. self.band_energy[kp]).


* **Type**

    list



#### Gpoints()
The list of generated G-points for each k-point (a double list), which
are used with the coefficients for each k-point and band to recreate
the wavefunction (e.g. self.Gpoints[kp] is the list of G-points for
k-point kp). The G-points depend on the k-point and reciprocal lattice
and therefore are identical for each band at the same k-point. Each
G-point is represented by integer multipliers (e.g. assuming
Gpoints[kp][n] == [n_1, n_2, n_3], then
G_n = n_1\*b_1 + n_2\*b_2 + n_3\*b_3)


* **Type**

    list



#### coeffs()
The list of coefficients for each k-point and band for reconstructing the wavefunction.
For non-spin-polarized, the first index corresponds to the kpoint and the second corresponds to the band
(e.g. self.coeffs[kp][b] corresponds to k-point kp and band b). For spin-polarized calculations,
the first index is for the spin. If the calculation was non-collinear, then self.coeffs[kp][b] will have
two columns (one for each component of the spinor).


* **Type**

    list


Acknowledgments:

    This code is based upon the Fortran program, WaveTrans, written by
    R. M. Feenstra and M. Widom from the Dept. of Physics at Carnegie
    Mellon University. To see the original work, please visit:
    [https://www.andrew.cmu.edu/user/feenstra/wavetrans/](https://www.andrew.cmu.edu/user/feenstra/wavetrans/)

Author: Mark Turiansky

Information is extracted from the given WAVECAR.


* **Parameters**


    * **filename** (*str*) – input file (default: WAVECAR)


    * **verbose** (*bool*) – determines whether processing information is shown


    * **precision** (*str*) – determines how fine the fft mesh is (normal or
    accurate), only the first letter matters


    * **vasp_type** (*str*) – determines the VASP type that is used, allowed
    values are [‘std’, ‘gam’, ‘ncl’] (only first letter is required)



#### _generate_G_points(kpoint: ndarray, gamma: bool = False)
Helper function to generate G-points based on nbmax.

This function iterates over possible G-point values and determines
if the energy is less than G_{cut}. Valid values are appended to
the output array. This function should not be called outside of
initialization.


* **Parameters**


    * **kpoint** (*np.array*) – the array containing the current k-point value


    * **gamma** (*bool*) – determines if G points for gamma-point only executable
    should be generated



* **Returns**

    a list containing valid G-points



#### _generate_nbmax()
Helper function that determines maximum number of b vectors for
each direction.

This algorithm is adapted from WaveTrans (see Class docstring). There
should be no reason for this function to be called outside of
initialization.


#### evaluate_wavefunc(kpoint: int, band: int, r: ndarray, spin: int = 0, spinor: int = 0)
Evaluates the wavefunction for a given position, r.

The wavefunction is given by the k-point and band. It is evaluated
at the given position by summing over the components. Formally,

psi_n^k (r) = sum_{i=1}^N c_i^{n,k} exp (i (k + G_i^{n,k}) cdot r)

where psi_n^k is the wavefunction for the nth band at k-point k, N is
the number of plane waves, c_i^{n,k} is the ith coefficient that
corresponds to the nth band and k-point k, and G_i^{n,k} is the ith
G-point corresponding to k-point k.

NOTE: This function is very slow; a discrete fourier transform is the
preferred method of evaluation (see Wavecar.fft_mesh).


* **Parameters**


    * **kpoint** (*int*) – the index of the kpoint where the wavefunction will be evaluated


    * **band** (*int*) – the index of the band where the wavefunction will be evaluated


    * **r** (*np.array*) – the position where the wavefunction will be evaluated


    * **spin** (*int*) – spin index for the desired wavefunction (only for
    ISPIN = 2, default = 0)


    * **spinor** (*int*) – component of the spinor that is evaluated (only used
    if vasp_type == ‘ncl’)



* **Returns**

    a complex value corresponding to the evaluation of the wavefunction



#### fft_mesh(kpoint: int, band: int, spin: int = 0, spinor: int = 0, shift: bool = True)
Places the coefficients of a wavefunction onto an fft mesh.

Once the mesh has been obtained, a discrete fourier transform can be
used to obtain real-space evaluation of the wavefunction. The output
of this function can be passed directly to numpy’s fft function. For
.. rubric:: Example

mesh = Wavecar(‘WAVECAR’).fft_mesh(kpoint, band)
evals = np.fft.ifftn(mesh)


* **Parameters**


    * **kpoint** (*int*) – the index of the kpoint where the wavefunction will be evaluated


    * **band** (*int*) – the index of the band where the wavefunction will be evaluated


    * **spin** (*int*) – the spin of the wavefunction for the desired
    wavefunction (only for ISPIN = 2, default = 0)


    * **spinor** (*int*) – component of the spinor that is evaluated (only used
    if vasp_type == ‘ncl’)


    * **shift** (*bool*) – determines if the zero frequency coefficient is
    placed at index (0, 0, 0) or centered



* **Returns**

    a numpy ndarray representing the 3D mesh of coefficients



#### get_parchg(poscar: Poscar, kpoint: int, band: int, spin: int | None = None, spinor: int | None = None, phase: bool = False, scale: int = 2)
Generates a Chgcar object, which is the charge density of the specified
wavefunction.

This function generates a Chgcar object with the charge density of the
wavefunction specified by band and kpoint (and spin, if the WAVECAR
corresponds to a spin-polarized calculation). The phase tag is a
feature that is not present in VASP. For a real wavefunction, the phase
tag being turned on means that the charge density is multiplied by the
sign of the wavefunction at that point in space. A warning is generated
if the phase tag is on and the chosen kpoint is not Gamma.

Note: Augmentation from the PAWs is NOT included in this function. The
maximal charge density will differ from the PARCHG from VASP, but the
qualitative shape of the charge density will match.


* **Parameters**


    * **poscar** (*pymatgen.io.vasp.inputs.Poscar*) – Poscar object that has the
    structure associated with the WAVECAR file


    * **kpoint** (*int*) – the index of the kpoint for the wavefunction


    * **band** (*int*) – the index of the band for the wavefunction


    * **spin** (*int*) – optional argument to specify the spin. If the Wavecar
    has ISPIN = 2, spin is None generates a Chgcar with total spin
    and magnetization, and spin == {0, 1} specifies just the spin
    up or down component.


    * **spinor** (*int*) – optional argument to specify the spinor component
    for noncollinear data wavefunctions (allowed values of None,
    0, or 1)


    * **phase** (*bool*) – flag to determine if the charge density is multiplied
    by the sign of the wavefunction. Only valid for real
    wavefunctions.


    * **scale** (*int*) – scaling for the FFT grid. The default value of 2 is at
    least as fine as the VASP default.



* **Returns**

    a pymatgen.io.vasp.outputs.Chgcar object



#### write_unks(directory: str)
Write the UNK files to the given directory.

Writes the cell-periodic part of the bloch wavefunctions from the
WAVECAR file to each of the UNK files. There will be one UNK file for
each of the kpoints in the WAVECAR file.

**NOTE**: wannier90 expects the full kpoint grid instead of the symmetry-
reduced one that VASP stores the wavefunctions on. You should run
a nscf calculation with ISYM=0 to obtain the correct grid.


* **Parameters**

    **directory** (*str*) – directory where the UNK files are written



### _class_ Waveder(cder_real: ndarray, cder_imag: ndarray)
Bases: `MSONable`

Representation of the WAVEDER file.

The LOPTICS tag produces a WAVEDER file which contains the derivative of the orbitals with respect to k.
Since the data is complex, we need to split it into the real and imaginary parts for JSON serialization.

**NOTE**: The way that VASP writes the WAVEDER and WAVEDERF has slightly different logic when indexing the bands.
This results in the formatted WAVDERF only indexing between filled bands. (i.e. all the matrix elements
are between the states i=1:8 and j=1:8 in a two atom Si calculation, which is likely a VASP bug).
As such, it is recommended to used the hidden `LVEL=.True.` flag in VASP which will force indexing over
all bands.

The order of the indices of the data are:
[

> band index1,
> band index2,
> kpoint index,
> spin index,
> cartesian direction,

]


#### cder_real()
Real part of the derivative of the orbitals with respect to k.


* **Type**

    numpy.ndarray



#### cder_imag()
Imaginary part of the derivative of the orbitals with respect to k.


* **Type**

    numpy.ndarray


Author: Miguel Dias Costa, Kamal Choudhary, Jimmy-Xuan Shen


#### _property_ cder()
Return the complex derivative of the orbitals with respect to k.


#### cder_imag(_: ndarra_ )

#### cder_real(_: ndarra_ )

#### _classmethod_ from_binary(filename, data_type='complex64')
Read the WAVEDER file and returns a Waveder object.


* **Parameters**


    * **filename** – Name of file containing WAVEDER.


    * **data_type** – Data type of the WAVEDER file. Default is complex64.
    If the file was generated with the “gamma” version of VASP,
    the data type can be either “float64” or “float32”.



* **Returns**

    Waveder object.



#### _classmethod_ from_formatted(filename)
Reads the WAVEDERF file and returns a Waveder object.

Note: This file is only produced when LOPTICS is true AND vasp has been
recompiled after uncommenting the line that calls
WRT_CDER_BETWEEN_STATES_FORMATTED in linear_optics.F
It is recommended to use from_binary instead since the binary file is
much smaller and contains the same information.


* **Parameters**

    **filename** (*str*) – The name of the WAVEDER file.



* **Returns**

    A Waveder object.



#### get_orbital_derivative_between_states(band_i, band_j, kpoint, spin, cart_dir)
Method returning a value
between bands band_i and band_j for k-point index, spin-channel and Cartesian direction.


* **Parameters**


    * **band_i** (*int*) – Index of band i


    * **band_j** (*int*) – Index of band j


    * **kpoint** (*int*) – Index of k-point


    * **spin** (*int*) – Index of spin-channel (0 or 1)


    * **cart_dir** (*int*) – Index of Cartesian direction (0,1,2)



* **Returns**

    a float value



#### _property_ nbands()
Returns the number of bands.


#### _property_ nkpoints()
Returns the number of k-points.


#### _property_ nspin()
Returns the number of spin channels.


### _class_ Xdatcar(filename, ionicstep_start=1, ionicstep_end=None, comment=None)
Bases: `object`

Class representing an XDATCAR file. Only tested with VASP 5.x files.


#### structures()
List of structures parsed from XDATCAR.


* **Type**

    list



#### comment()
Optional comment string.


* **Type**

    str


Authors: Ram Balachandran

Init a Xdatcar.


* **Parameters**


    * **filename** (*str*) – Filename of input XDATCAR file.


    * **ionicstep_start** (*int*) – Starting number of ionic step.


    * **ionicstep_end** (*int*) – Ending number of ionic step.


    * **comment** (*str*) – Optional comment attached to this set of structures.



#### concatenate(filename, ionicstep_start=1, ionicstep_end=None)
Concatenate structures in file to Xdatcar.


* **Parameters**


    * **filename** (*str*) – Filename of XDATCAR file to be concatenated.


    * **ionicstep_start** (*int*) – Starting number of ionic step.


    * **ionicstep_end** (*int*) – Ending number of ionic step.


TODO (rambalachandran): Requires a check to ensure if the new concatenating file

    has the same lattice structure and atoms as the Xdatcar class.


#### get_str(ionicstep_start: int = 1, ionicstep_end: int | None = None, significant_figures: int = 8)
Write  Xdatcar class to a string.


* **Parameters**


    * **ionicstep_start** (*int*) – Starting number of ionic step.


    * **ionicstep_end** (*int*) – Ending number of ionic step.


    * **significant_figures** (*int*) – Number of significant figures.



#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### _property_ natoms()
Sequence of number of sites of each type associated with the Poscar.
Similar to 7th line in vasp 5+ Xdatcar.


#### _property_ site_symbols()
Sequence of symbols associated with the Xdatcar. Similar to 6th line in
vasp 5+ Xdatcar.


#### write_file(filename, \*\*kwargs)
Write Xdatcar class into a file.


* **Parameters**


    * **filename** (*str*) – Filename of output XDATCAR file.


    * **\*\*kwargs** – Supported kwargs are the same as those for the
    Xdatcar.get_string method and are passed through directly.



### _parse_from_incar(filename, key)
Helper function to parse a parameter from the INCAR.


### _parse_parameters(val_type, val)
Helper function to convert a Vasprun parameter into the proper type.
Boolean, int and float types are converted.


* **Parameters**


    * **val_type** – Value type parsed from vasprun.xml.


    * **val** – Actual string value parsed for vasprun.xml.



### _parse_v_parameters(val_type, val, filename, param_name)
Helper function to convert a Vasprun array-type parameter into the proper
type. Boolean, int and float types are converted.


* **Parameters**


    * **val_type** – Value type parsed from vasprun.xml.


    * **val** – Actual string value parsed for vasprun.xml.


    * **filename** – Fullpath of vasprun.xml. Used for robust error handling.
    E.g., if vasprun.xml contains

    ```
    **
    ```

    \* for some Incar parameters,
    the code will try to read from an INCAR file present in the same
    directory.


    * **param_name** – Name of parameter.



* **Returns**

    Parsed value.



### _parse_varray(elem)

### _vasprun_float(f)
Large numbers are often represented as **\*\*\*\*\*** in the vasprun.
This function parses these values as np.nan.


### get_adjusted_fermi_level(efermi, cbm, band_structure)
When running a band structure computations the Fermi level needs to be
take from the static run that gave the charge density used for the non-self
consistent band structure run. Sometimes this Fermi level is however a
little too low because of the mismatch between the uniform grid used in
the static run and the band structure k-points (e.g., the VBM is on Gamma
and the Gamma point is not in the uniform mesh). Here we use a procedure
consisting in looking for energy levels higher than the static Fermi level
(but lower than the LUMO) if any of these levels make the band structure
appears insulating and not metallic anymore, we keep this adjusted fermi
level. This procedure has shown to detect correctly most insulators.


* **Parameters**


    * **efermi** (*float*) – The Fermi energy of the static run.


    * **cbm** (*float*) – The conduction band minimum of the static run.


    * **band_structure** ([*BandStructureSymmLine*](pymatgen.electronic_structure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine)) – A band structure object.



* **Returns**

    A new adjusted Fermi level.



* **Return type**

    float



### get_band_structure_from_vasp_multiple_branches(dir_name, efermi=None, projections=False)
This method is used to get band structure info from a VASP directory. It
takes into account that the run can be divided in several branches named
“branch_x”. If the run has not been divided in branches the method will
turn to parsing vasprun.xml directly.

The method returns None is there’s a parsing error


* **Parameters**


    * **dir_name** – Directory containing all bandstructure runs.


    * **efermi** – Efermi for bandstructure.


    * **projections** – True if you want to get the data on site projections if
    any. Note that this is sometimes very large



* **Returns**

    A BandStructure Object


## pymatgen.io.vasp.sets module

This module defines the VaspInputSet abstract base class and a concrete implementation for the parameters developed
and tested by the core team of pymatgen, including the Materials Virtual Lab, Materials Project and the MIT high
throughput project. The basic concept behind an input set is to specify a scheme to generate a consistent set of VASP
inputs from a structure without further user intervention. This ensures comparability across runs.

Read the following carefully before implementing new input sets:


1. 99% of what needs to be done can be done by specifying user_incar_settings to override some of the defaults of
various input sets. Unless there is an extremely good reason to add a new set, **do not** add one. E.g., if you want
to turn the Hubbard U off, just set “LDAU”: False as a user_incar_setting.


2. All derivative input sets should inherit appropriate configurations (e.g., from MPRelaxSet), and more often than
not, DictSet should be the superclass. Proper superclass delegation should be used where possible. In particular,
you are not supposed to implement your own as_dict or from_dict for derivative sets unless you know what you are
doing. Improper overriding the as_dict and from_dict protocols is the major cause of implementation headaches. If
you need an example, look at how the MPStaticSet is initialized.

The above are recommendations. The following are **UNBREAKABLE** rules:


1. All input sets must take in a structure, list of structures or None as the first argument. If None, the input set
should perform a stateless initialization and before any output can be written, a structure must be set.


2. user_incar_settings, user_kpoints_settings and user_<whatever>_settings are ABSOLUTE. Any new sets you implement
must obey this. If a user wants to override your settings, you assume he knows what he is doing. Do not
magically override user supplied settings. You can issue a warning if you think the user is wrong.


3. All input sets must save all supplied args and kwargs as instance variables. E.g., self.arg = arg and
self.kwargs = kwargs in the __init__. This ensures the as_dict and from_dict work correctly.


### _exception_ BadInputSetWarning()
Bases: `UserWarning`

Warning class for bad but legal VASP inputs.


### _class_ DictSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, config_dict: dict[str, Any] | None = None, files_to_transfer=None, user_incar_settings=None, user_kpoints_settings=None, user_potcar_settings=None, constrain_total_magmom: bool = False, sort_structure: bool = True, user_potcar_functional: UserPotcarFunctional | None = None, force_gamma: bool = False, reduce_structure=None, vdw=None, use_structure_charge: bool = False, standardize: bool = False, sym_prec=0.1, international_monoclinic: bool = True, validate_magmom: bool = True)
Bases: `VaspInputSet`

Concrete implementation of VaspInputSet that is initialized from a dict
settings. This allows arbitrary settings to be input. In general,
this is rarely used directly unless there is a source of settings in yaml
format (e.g., from a REST interface). It is typically used by other
VaspInputSets for initialization.

Special consideration should be paid to the way the MAGMOM initialization
for the INCAR is done. The initialization differs depending on the type of
structure and the configuration settings. The order in which the magmom is
determined is as follows:


1. If the site itself has a magmom setting (i.e. site.properties[“magmom”] = float),

    that is used. This can be set with structure.add_site_property().


2. If the species of the site has a spin setting, that is used. This can be set

    with structure.add_spin_by_element().


3. If the species itself has a particular setting in the config file, that
is used, e.g., Mn3+ may have a different magmom than Mn4+.


4. Lastly, the element symbol itself is checked in the config file. If
there are no settings, a default value of 0.6 is used.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **config_dict** (*dict*) – The config dictionary to use.


    * **files_to_transfer** (*dict*) – A dictionary of {filename: filepath}. This allows the transfer of files from a
    previous calculation.


    * **user_incar_settings** (*dict*) – User INCAR settings. This allows a user to override INCAR settings, e.g.,
    setting a different MAGMOM for various elements or species. Note that in the new scheme,
    ediff_per_atom and hubbard_u are no longer args. Instead, the config_dict supports EDIFF_PER_ATOM and
    EDIFF keys. The former scales with # of atoms, the latter does not. If both are present,
    EDIFF is preferred. To force such settings, just supply user_incar_settings={“EDIFF”: 1e-5,
    “LDAU”: False} for example. The keys ‘LDAUU’, ‘LDAUJ’, ‘LDAUL’ are special cases since
    pymatgen defines different values depending on what anions are present in the structure,
    so these keys can be defined in one of two ways, e.g. either {“LDAUU”:{“O”:{“Fe”:5}}} to set LDAUU
    for Fe to 5 in an oxide, or {“LDAUU”:{“Fe”:5}} to set LDAUU to 5 regardless of the input structure.
    If a None value is given, that key is unset. For example, {“ENCUT”: None} will remove ENCUT from the
    incar settings.


    * **user_kpoints_settings** (*dict** or *[*Kpoints*](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints)) – Allow user to override kpoints setting by supplying a dict.
    E.g., {“reciprocal_density”: 1000}. User can also supply Kpoints object. Default is None.


    * **(****dict** (*user_potcar_settings*) – Allow user to override POTCARs. E.g., {“Gd”: “Gd_3”}. This is generally not
    recommended. Default is None.


    * **constrain_total_magmom** (*bool*) – Whether to constrain the total magmom (NUPDOWN in INCAR) to be the sum of
    the expected MAGMOM for all species. Defaults to False.


    * **sort_structure** (*bool*) – Whether to sort the structure (using the default sort order of electronegativity)
    before generating input files. Defaults to True, the behavior you would want most of the time. This
    ensures that similar atomic species are grouped together.


    * **user_potcar_functional** (*str*) – Functional to use. Default (None) is to use the functional in the config
    dictionary. Valid values: “PBE”, “PBE_52”, “PBE_54”, “LDA”, “LDA_52”, “LDA_54”, “PW91”,
    “LDA_US”, “PW91_US”.


    * **force_gamma** (*bool*) – Force gamma centered kpoint generation. Default (False) is to use the Automatic
    Density kpoint scheme, which will use the Gamma centered generation scheme for hexagonal
    cells, and Monkhorst-Pack otherwise.


    * **reduce_structure** (*None/str*) – Before generating the input files, generate the reduced structure. Default (
    None), does not alter the structure. Valid values: None, “niggli”, “LLL”.


    * **vdw** – Adds default parameters for van-der-Waals functionals supported by VASP to INCAR. Supported
    functionals are: DFT-D2, undamped DFT-D3, DFT-D3 with Becke-Jonson damping, Tkatchenko-Scheffler,
    Tkatchenko-Scheffler with iterative Hirshfeld partitioning, [MBD@rSC](mailto:MBD@rSC), dDsC, Dion’s vdW-DF, DF2,
    optPBE, optB88, optB86b and rVV10.


    * **use_structure_charge** (*bool*) – If set to True, then the overall charge of the structure (structure.charge)
    is used to set the NELECT variable in the INCAR. Default is False.


    * **standardize** (*float*) – Whether to standardize to a primitive standard cell. Defaults to False.


    * **sym_prec** (*float*) – Tolerance for symmetry finding.


    * **international_monoclinic** (*bool*) – Whether to use international convention (vs Curtarolo) for monoclinic.
    Defaults True.


    * **validate_magmom** (*bool*) – Ensure that the missing magmom values are filled in with the VASP default value
    of 1.0.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### calculate_ng(max_prime_factor: int = 7, must_inc_2: bool = True, custom_encut: float | None = None, custom_prec: str | None = None)
Calculates the NGX, NGY, and NGZ values using the information available in the INCAR and POTCAR
This is meant to help with making initial guess for the FFT grid so we can interact with the Charge density API.


* **Parameters**


    * **max_prime_factor** (*int*) – the valid prime factors of the grid size in each direction
    VASP has many different setting for this to handle many compiling options.
    For typical MPI options all prime factors up to 7 are allowed


    * **must_inc_2** (*bool*) – Whether 2 must be a prime factor of the result. Defaults to True.


    * **custom_encut** (*float** | **None*) – Calculates the FFT grid parameters using a custom
    ENCUT that may be different from what is generated by the input set. Defaults to None.
    Do *not* use this unless you know what you are doing.


    * **custom_prec** (*str** | **None*) – Calculates the FFT grid parameters using a custom prec
    that may be different from what is generated by the input set. Defaults to None.
    Do *not* use this unless you know what you are doing.



#### estimate_nbands()
Estimate the number of bands that VASP will initialize a
calculation with by default. Note that in practice this
can depend on # of cores (if not set explicitly).
Note that this formula is slightly different than the formula on the VASP wiki
(as of July 2023). This is because the formula in the source code (main.F) is
slightly different than what is on the wiki.


#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints(_: [Kpoints](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints) | Non_ )
Returns a KPOINTS file using the fully automated grid method. Uses
Gamma centered meshes for hexagonal cells and Monk grids otherwise.

If KSPACING is set in user_incar_settings (or the INCAR file), no
file is created because VASP will automatically generate the kpoints.

Algorithm:

    Uses a simple approach scaling the number of divisions along each
    reciprocal lattice vector proportional to its length.


#### _property_ nelect(_: floa_ )
Gets the default number of electrons for a given structure.


#### _property_ poscar(_: Posca_ )
Poscar object.


#### _property_ potcar_functional(_: Literal['PBE', 'PBE_52', 'PBE_54', 'LDA', 'LDA_52', 'LDA_54', 'PW91', 'LDA_US', 'PW91_US'] | Non_ )
Returns the functional used for POTCAR generation.


#### _property_ structure(_: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure_ )
Structure


#### write_input(output_dir: str, make_dir_if_not_present: bool = True, include_cif: bool = False, potcar_spec: bool = False, zip_output: bool = False)
Writes out all input to a directory.


* **Parameters**


    * **output_dir** (*str*) – Directory to output the VASP input files


    * **make_dir_if_not_present** (*bool*) – Set to True if you want the
    directory (and the whole path) to be created if it is not
    present.


    * **include_cif** (*bool*) – Whether to write a CIF file in the output
    directory for easier opening by VESTA.


    * **potcar_spec** (*bool*) – Instead of writing the POTCAR, write a “POTCAR.spec”.
    This is intended to help sharing an input set with people who might
    not have a license to specific Potcar files. Given a “POTCAR.spec”,
    the specific POTCAR file can be re-generated using pymatgen with the
    “generate_potcar” function in the pymatgen CLI.


    * **zip_output** (*bool*) – Whether to zip each VASP input file written to the output directory.



### _class_ LobsterSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, isym: int = 0, ismear: int = -5, reciprocal_density: int | None = None, address_basis_file: str | None = None, user_supplied_basis: dict | None = None, \*\*kwargs)
Bases: `DictSet`

Input set to prepare VASP runs that can be digested by Lobster (See cohp.de).


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **isym** (*int*) – ISYM entry for INCAR, only isym=-1 and isym=0 are allowed


    * **ismear** (*int*) – ISMEAR entry for INCAR, only ismear=-5 and ismear=0 are allowed


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal volume


    * **user_supplied_basis** (*dict*) – dict including basis functions for all elements in structure,
    e.g. {“Fe”: “3d 3p 4s”, “O”: “2s 2p”}; if not supplied, a standard basis is used


    * **address_basis_file** (*str*) – address to a file similar to “BASIS_PBE_54_standaard.yaml”
    in pymatgen.io.lobster.lobster_basis


    * **user_potcar_settings** (*dict*) – dict including potcar settings for all elements in structure,
    e.g. {“Fe”: “Fe_pv”, “O”: “O”}; if not supplied, a standard basis is used.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _valid_potcars(_: Sequence[str] | Non_ _ = ('PBE_52', 'PBE_54'_ )

#### _property_ incar(_: Inca_ )
Incar


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MITMDSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, start_temp: float = 0.0, end_temp: float = 300.0, nsteps: int = 1000, time_step: float = 2, spin_polarized=False, \*\*kwargs)
Bases: `DictSet`

Class for writing a vasp md run. This DOES NOT do multiple stage
runs.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **start_temp** (*float*) – Starting temperature.


    * **end_temp** (*float*) – Final temperature.


    * **nsteps** (*int*) – Number of time steps for simulations. NSW parameter.


    * **time_step** (*float*) – The time step for the simulation. The POTIM
    parameter. Defaults to 2fs.


    * **spin_polarized** (*bool*) – Whether to do spin polarized calculations.
    The ISPIN parameter. Defaults to False.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ kpoints(_: Kpoint_ )
Kpoints


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MITNEBSet(structures, unset_encut=False, \*\*kwargs)
Bases: `DictSet`

Class for writing NEB inputs. Note that EDIFF is not on a per atom
basis for this input set.


* **Parameters**


    * **structures** – List of Structure objects.


    * **unset_encut** (*bool*) – Whether to unset ENCUT.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _process_structures(structures)
Remove any atom jumps across the cell.


#### _property_ poscar()
Poscar for structure of first end point.


#### _property_ poscars()
List of Poscars.


#### user_potcar_functional(_: UserPotcarFunctiona_ )

#### write_input(output_dir, make_dir_if_not_present=True, write_cif=False, write_path_cif=False, write_endpoint_inputs=False)
NEB inputs has a special directory structure where inputs are in 00,
01, 02, ….


* **Parameters**


    * **output_dir** (*str*) – Directory to output the VASP input files


    * **make_dir_if_not_present** (*bool*) – Set to True if you want the
    directory (and the whole path) to be created if it is not
    present.


    * **write_cif** (*bool*) – If true, writes a cif along with each POSCAR.


    * **write_path_cif** (*bool*) – If true, writes a cif for each image.


    * **write_endpoint_inputs** (*bool*) – If true, writes input files for
    running endpoint calculations.



### _class_ MITRelaxSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, \*\*kwargs)
Bases: `DictSet`

Standard implementation of VaspInputSet utilizing parameters in the MIT
High-throughput project.
The parameters are chosen specifically for a high-throughput project,
which means in general pseudopotentials with fewer electrons were chosen.

Please refer:

```default
A Jain, G. Hautier, C. Moore, S. P. Ong, C. Fischer, T. Mueller,
K. A. Persson, G. Ceder. A high-throughput infrastructure for density
functional theory calculations. Computational Materials Science,
2011, 50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023
```


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **\*\*kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF': 1e-05, 'ENCUT': 520, 'IBRION': 2, 'ICHARG': 1, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'ISYM': 0, 'LDAU': True, 'LDAUJ': {'F': {'Ag': 0, 'Co': 0, 'Cr': 0, 'Cu': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Nb': 0, 'Ni': 0, 'Re': 0, 'Ta': 0, 'V': 0, 'W': 0}, 'O': {'Ag': 0, 'Co': 0, 'Cr': 0, 'Cu': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Nb': 0, 'Ni': 0, 'Re': 0, 'Ta': 0, 'V': 0, 'W': 0}, 'S': {'Fe': 0, 'Mn': 0}}, 'LDAUL': {'F': {'Ag': 2, 'Co': 2, 'Cr': 2, 'Cu': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Nb': 2, 'Ni': 2, 'Re': 2, 'Ta': 2, 'V': 2, 'W': 2}, 'O': {'Ag': 2, 'Co': 2, 'Cr': 2, 'Cu': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Nb': 2, 'Ni': 2, 'Re': 2, 'Ta': 2, 'V': 2, 'W': 2}, 'S': {'Fe': 2, 'Mn': 2.5}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4, 'Fe': 4.0, 'Mn': 3.9, 'Mo': 4.38, 'Nb': 1.5, 'Ni': 6, 'Re': 2, 'Ta': 2, 'V': 3.1, 'W': 4.0}, 'O': {'Ag': 1.5, 'Co': 3.4, 'Cr': 3.5, 'Cu': 4, 'Fe': 4.0, 'Mn': 3.9, 'Mo': 4.38, 'Nb': 1.5, 'Ni': 6, 'Re': 2, 'Ta': 2, 'V': 3.1, 'W': 4.0}, 'S': {'Fe': 1.9, 'Mn': 2.5}}, 'LORBIT': '11', 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 200, 'NELMIN': 6, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'length': 25}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': {'hash': 'd6854224d20e3de6e6fd7399503791d1', 'symbol': 'Ac'}, 'Ag': {'hash': 'e8ffa02fe3f3a51338ac1ac91ae968b9', 'symbol': 'Ag'}, 'Al': {'hash': 'a6fd9a46aec185f4ad2acd0cbe4ae2fa', 'symbol': 'Al'}, 'Ar': {'hash': 'e782fc6292623b396091bf8b871c272f', 'symbol': 'Ar'}, 'As': {'hash': '8005364db225a254e52cba350bedd032', 'symbol': 'As'}, 'Au': {'hash': 'a9182d436a13194b744640ac940ab9b0', 'symbol': 'Au'}, 'B': {'hash': '18ed2875dfa6305324cec3d7d59273ae', 'symbol': 'B'}, 'Ba': {'hash': 'c0477913afb63dfae3439f3534fbf0ed', 'symbol': 'Ba_sv'}, 'Be': {'hash': 'fb974e44d56a8c62c6bbd1a1eb70c3a7', 'symbol': 'Be'}, 'Bi': {'hash': 'e29661c79d59abae3b3ba69eae24b1a5', 'symbol': 'Bi'}, 'Br': {'hash': '40f9594b4506684a69158c8975cfb9d6', 'symbol': 'Br'}, 'C': {'hash': 'c0a8167dbb174fe492a3db7f5006c0f8', 'symbol': 'C'}, 'Ca': {'hash': 'eb006721e214c04b3c13146e81b3a27d', 'symbol': 'Ca_sv'}, 'Cd': {'hash': '0506b2d0ac28d5fe2b5ced77a701aa86', 'symbol': 'Cd'}, 'Ce': {'hash': 'ff3a09f2ff91798e58eb4b9854e9be4a', 'symbol': 'Ce'}, 'Cl': {'hash': '779b9901046c78fe51c5d80224642aeb', 'symbol': 'Cl'}, 'Co': {'hash': 'b169bca4e137294d2ab3df8cbdd09083', 'symbol': 'Co'}, 'Cr': {'hash': '82c14307937c7509fda4e9bc023d243d', 'symbol': 'Cr'}, 'Cs': {'hash': '096b53a7d80cc0086976bcda50d536e5', 'symbol': 'Cs_sv'}, 'Cu': {'hash': '8ca4e43a30de0c397e51f16bbb20d678', 'symbol': 'Cu'}, 'Dy': {'hash': 'd4a05220ab0a2d4c03a76872ea724a1e', 'symbol': 'Dy_3'}, 'Er': {'hash': 'daa65a04877317f8c3c593ddeaa8a132', 'symbol': 'Er_3'}, 'Eu': {'hash': 'd466d046adf21f6146ee9644049ea268', 'symbol': 'Eu'}, 'F': {'hash': '180141c33d032bfbfff30b3bea9d23dd', 'symbol': 'F'}, 'Fe': {'hash': '9530da8244e4dac17580869b4adab115', 'symbol': 'Fe'}, 'Ga': {'hash': '6e0b9d58412b1bfcd7252aff13d476c2', 'symbol': 'Ga'}, 'Gd': {'hash': '1f0d42b1e5f6769d319d3f247992aeb9', 'symbol': 'Gd'}, 'Ge': {'hash': '79e788788c31e196a460553010512d3f', 'symbol': 'Ge'}, 'H': {'hash': 'bb43c666e3d36577264afe07669e9582', 'symbol': 'H'}, 'He': {'hash': '47f9434aa3db96c85d7c4b3e4c2df09b', 'symbol': 'He'}, 'Hf': {'hash': 'b113f150cbf9c736f8244a6c25b0482e', 'symbol': 'Hf'}, 'Hg': {'hash': 'c2f15dfb5fd53396c5427635e5019160', 'symbol': 'Hg'}, 'Ho': {'hash': '661891464a27e87cf7e1324dd1893b77', 'symbol': 'Ho_3'}, 'I': {'hash': 'f4ff16a495dd361ff5824ee61b418bb0', 'symbol': 'I'}, 'In': {'hash': '7df38c0cdb4e6d9a9b93f09d690bb3ae', 'symbol': 'In'}, 'Ir': {'hash': 'dbcf7dcc6f4fb40df7b3d26904f60a66', 'symbol': 'Ir'}, 'K': {'hash': '3e84f86d37f203a4fb01de36af57e430', 'symbol': 'K_sv'}, 'Kr': {'hash': '39b9b85ae3982e6c012fb549b2840ce5', 'symbol': 'Kr'}, 'La': {'hash': '9b3ce03d18f7c0b40471a817ff91b287', 'symbol': 'La'}, 'Li': {'hash': '65e83282d1707ec078c1012afbd05be8', 'symbol': 'Li'}, 'Lu': {'hash': 'd40a90babf1224b88ffb4c3273ac3848', 'symbol': 'Lu_3'}, 'Mg': {'hash': '1771eb72adbbfa6310d66e7517e49930', 'symbol': 'Mg'}, 'Mn': {'hash': 'd082dba29b57ab59b3165e605dbf71b8', 'symbol': 'Mn'}, 'Mo': {'hash': '84e18fd84a98e3d7fa8f055952410df0', 'symbol': 'Mo_pv'}, 'N': {'hash': 'b98fd027ddebc67da4063ff2cabbc04b', 'symbol': 'N'}, 'Na': {'hash': '1a89e79f7e21d99e8cf5788979f6a987', 'symbol': 'Na'}, 'Nb': {'hash': '7bcee99a4dc3094be0f9fd7961c02966', 'symbol': 'Nb_pv'}, 'Nd': {'hash': '0c64e63070cee837c967283fffa001df', 'symbol': 'Nd'}, 'Ne': {'hash': '52064eee378b9e37a295a674f1c278f0', 'symbol': 'Ne'}, 'Ni': {'hash': '653f5772e68b2c7fd87ffd1086c0d710', 'symbol': 'Ni'}, 'Np': {'hash': '20cb30b714200c4db870550b288ac4cd', 'symbol': 'Np'}, 'O': {'hash': '7a25bc5b9a5393f46600a4939d357982', 'symbol': 'O'}, 'Os': {'hash': '35c2cb48d48a9c38c40fb82bbe70626d', 'symbol': 'Os'}, 'P': {'hash': '7dc3393307131ae67785a0cdacb61d5f', 'symbol': 'P'}, 'Pa': {'hash': 'a1fdb1089d0727f415416ec8082246ba', 'symbol': 'Pa'}, 'Pb': {'hash': '704c2c967247d7f84090d2536c91877d', 'symbol': 'Pb'}, 'Pd': {'hash': 'a395eb3aaf2fcab12fac3030a1146f61', 'symbol': 'Pd'}, 'Pm': {'hash': 'a2c9485ea86b2a7cf175077e6e5c7b3e', 'symbol': 'Pm'}, 'Pr': {'hash': '92f191499bf5346ea652bb806350ad87', 'symbol': 'Pr'}, 'Pt': {'hash': 'a604ea3c6a9cc23c739b762f625cf449', 'symbol': 'Pt'}, 'Pu': {'hash': 'f1d01e845dccc52d448679911f301a73', 'symbol': 'Pu'}, 'Rb': {'hash': 'e447c648d870b066b3514e6b800727ab', 'symbol': 'Rb_pv'}, 'Re': {'hash': '72385e193c92a8acfe17ea49004c2be1', 'symbol': 'Re'}, 'Rh': {'hash': '2c3dba3fcc6058ca1b1cfa75e45084bc', 'symbol': 'Rh'}, 'Ru': {'hash': '7925f4d4b68076d70af7cd86eef9ba8d', 'symbol': 'Ru_pv'}, 'S': {'hash': 'd368db6899d8839859bbee4811a42a88', 'symbol': 'S'}, 'Sb': {'hash': 'd82c022b02fc5344e85bd1909f9ee3e7', 'symbol': 'Sb'}, 'Sc': {'hash': 'dc386f505ad0c43385a7715b4111cb75', 'symbol': 'Sc_sv'}, 'Se': {'hash': '67a8804ede9f1112726e3d136978ef19', 'symbol': 'Se'}, 'Si': {'hash': 'b2b0ea6feb62e7cde209616683b8f7f5', 'symbol': 'Si'}, 'Sm': {'hash': 'e5e274e7cd99602ca81d146155abdf88', 'symbol': 'Sm_3'}, 'Sn': {'hash': '849b0795e148f93113a06be8fd5f5001', 'symbol': 'Sn_d'}, 'Sr': {'hash': 'ca6a5429c120a0ab705824386a76fe5b', 'symbol': 'Sr_sv'}, 'Ta': {'hash': 'd4e2cfe9338ef80da592d5bb9dc782c7', 'symbol': 'Ta'}, 'Tb': {'hash': '0790955c547003956c0fd4f080f7f508', 'symbol': 'Tb_3'}, 'Tc': {'hash': '9592642886319309a39d55c5717c6f48', 'symbol': 'Tc'}, 'Te': {'hash': '72719856e22fb1d3032df6f96d98a0f2', 'symbol': 'Te'}, 'Th': {'hash': 'aea79f322180fa6f0bfa74cb2a156dcf', 'symbol': 'Th'}, 'Ti': {'hash': 'c617e8b539c3f44a0ab6e8da2a92d318', 'symbol': 'Ti'}, 'Tl': {'hash': '2aa0d5406aaab7ebfbc761da382f1352', 'symbol': 'Tl'}, 'Tm': {'hash': '94a07cb7949b01305cb161da0cbfb492', 'symbol': 'Tm_3'}, 'U': {'hash': '72702eabbb1bc02b4167590dc848ed5d', 'symbol': 'U'}, 'V': {'hash': '7f1297a2e1d963e2a4d81b61f85e4ded', 'symbol': 'V_pv'}, 'W': {'hash': '2a33e0d5c700640535f60ac0a12177ab', 'symbol': 'W_pv'}, 'Xe': {'hash': '338472e581f58b41d37c002a5e22353b', 'symbol': 'Xe'}, 'Y': {'hash': '4ed187e77cd54f198bb88020278b143d', 'symbol': 'Y_sv'}, 'Yb': {'hash': '9f472bd422f640710f7d93e2d9ce89f4', 'symbol': 'Yb'}, 'Zn': {'hash': 'e35ee27f8483a63bb68dbc236a343af3', 'symbol': 'Zn'}, 'Zr': {'hash': 'd221d2c0bac4f8e81af2f5c42a314274', 'symbol': 'Zr'}}, 'POTCAR_FUNCTIONAL': 'PBE'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPAbsorptionSet(structure, mode='IPA', copy_wavecar=True, nbands=None, nbands_factor=2, reciprocal_density=400, nkred=None, nedos=2001, prev_incar=None, \*\*kwargs)
Bases: `MPRelaxSet`

MP input set for generating frequency dependent dielectrics.
Two modes are supported: “IPA” or “RPA”.
A typical sequence is mode=”STATIC” -> mode=”IPA” -> mode=”RPA”(optional)
For all steps other than the first one (static), the
recommendation is to use from_prev_calculation on the preceding run in
the series. It is important to ensure Gamma centred kpoints for the RPA step.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **prev_incar** (*Incar/string*) – Incar file from previous run.


    * **mode** (*str*) – Supported modes are “IPA”, “RPA”


    * **copy_wavecar** (*bool*) – Whether to copy the WAVECAR from a previous run. Defaults to True.


    * **nbands** (*int*) – For subsequent calculations, it is generally
    recommended to perform NBANDS convergence starting from the
    NBANDS of the previous run for DIAG, and to use the exact same
    NBANDS for RPA. This parameter is used by
    from_previous_calculation to set nband.


    * **nbands_factor** (*int*) – Multiplicative factor for NBANDS when starting
    from a previous calculation. Only applies if mode==”IPA”.
    Need to be tested for convergence.


    * **reciprocal_density** – the k-points density


    * **nkred** – the reduced number of kpoints to calculate, equal to the k-mesh. Only applies in “RPA” mode
    because of the q->0 limit.


    * **nedos** – the density of DOS, default: 2001.


    * **\*\*kwargs** – All kwargs supported by DictSet. Typically, user_incar_settings is a commonly used option.



#### SUPPORTED_MODES(_ = ('IPA', 'RPA'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, mode, \*\*kwargs)
Generate a set of VASP input files for absorption calculation


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml of previous vasp run.


    * **mode** (*str*) – Supported modes are “IPA”, “RPA” (default)


    * **\*\*kwargs** – All kwargs supported by MPAbsorptionsSet, other than structure.



#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints(_: Kpoint_ )
Generate gamma center k-points mesh grid for optical calculation. It is not mandatory for ‘ALGO = Exact’,
but is requested by ‘ALGO = CHI’ calculation.


#### override_from_prev_calc(prev_calc_dir='.', \*\*kwargs)
Update the input set to include settings from a previous calculation.


* **Parameters**


    * **prev_calc_dir** (*str*) – The path to the previous calculation directory.


    * **\*\*kwargs** – unused



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPHSEBSSet(structure, user_incar_settings=None, added_kpoints=None, mode='Gap', reciprocal_density=None, copy_chgcar=True, kpoints_line_density=20, \*\*kwargs)
Bases: `MPHSERelaxSet`

Implementation of a VaspInputSet for HSE band structure computations.
Remember that HSE band structures must be self-consistent in VASP. A
band structure along symmetry lines for instance needs BOTH a uniform
grid with appropriate weights AND a path along the lines with weight 0.

Thus, the “Uniform” mode is just like regular static SCF but allows
adding custom kpoints (e.g., corresponding to known VBM/CBM) to the
uniform grid that have zero weight (e.g., for better gap estimate).

The “Gap” mode behaves just like the “Uniform” mode, however, if starting
from a previous calculation, the VBM and CBM k-points will automatically
be added to `added_kpoints`.

The “Line” mode is just like Uniform mode, but additionally adds
k-points along symmetry lines with zero weight.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to compute


    * **user_incar_settings** (*dict*) – A dict specifying additional incar
    settings


    * **added_kpoints** (*list*) – a list of kpoints (list of 3 number list)
    added to the run. The k-points are in fractional coordinates


    * **mode** (*str*) – “Line” - generate k-points along symmetry lines for
    bandstructure. “Uniform” - generate uniform k-points grid.


    * **reciprocal_density** (*int*) – k-point density to use for uniform mesh.


    * **copy_chgcar** (*bool*) – Whether to copy the CHGCAR of a previous run.


    * **kpoints_line_density** (*int*) – k-point density for high symmetry lines


    * **\*\*kwargs** (*dict*) – Any other parameters to pass into DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for HSE calculations from a
directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs
    (vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPHSEBSStaticSet, other than
    prev_structure which is determined from the previous calc dir.



#### _property_ kpoints(_: Kpoint_ )
Kpoints


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPHSERelaxSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, \*\*kwargs)
Bases: `DictSet`

Same as the MPRelaxSet, but with HSE parameters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **\*\*kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'All', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'HFSCREEN': 0.2, 'IBRION': 2, 'ICHARG': 1, 'ISIF': 3, 'ISMEAR': 0, 'ISPIN': 2, 'LHFCALC': True, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'PRECFOCK': 'Fast', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 50}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_52'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPMDSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, start_temp: float = 0.0, end_temp: float = 300.0, nsteps: int = 1000, time_step: float = 2, spin_polarized=False, \*\*kwargs)
Bases: `DictSet`

This a modified version of the old MITMDSet pre 2018/03/12.

This set serves as the basis for the amorphous skyline paper.


1. Aykol, M.; Dwaraknath, S. S.; Sun, W.; Persson, K. A. Thermodynamic
Limit for Synthesis of Metastable Inorganic Materials. Sci. Adv. 2018,
4 (4).

Class for writing a vasp md run. This DOES NOT do multiple stage runs.
Precision remains normal, to increase accuracy of stress tensor.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **start_temp** (*int*) – Starting temperature.


    * **end_temp** (*int*) – Final temperature.


    * **nsteps** (*int*) – Number of time steps for simulations. NSW parameter.


    * **time_step** (*int*) – The time step for the simulation. The POTIM
    parameter. Defaults to 2fs.


    * **spin_polarized** (*bool*) – Whether to do spin polarized calculations.
    The ISPIN parameter. Defaults to False.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints(_: Kpoint_ )
Kpoints


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPMetalRelaxSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, \*\*kwargs)
Bases: `DictSet`

Implementation of VaspInputSet utilizing parameters in the public
Materials Project, but with tuning for metals. Key things are a denser
k point density, and a.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **\*\*kwargs** – Same as those supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPNMRSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, mode: Literal['cs', 'efg'] = 'cs', isotopes: list | None = None, prev_incar: Incar = None, reciprocal_density: int = 100, \*\*kwargs)
Bases: `MPStaticSet`

Init a MPNMRSet.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to compute


    * **mode** (*str*) – The NMR calculation to run
    “cs”: for Chemical Shift
    “efg” for Electric Field Gradient


    * **isotopes** (*list*) – list of Isotopes for quadrupole moments


    * **prev_incar** (*Incar*) – Incar file from previous run.


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal volume. Defaults to 100.


    * **\*\*kwargs** – kwargs supported by MPStaticSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ incar(_: Inca_ )
Incar


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPNonSCFSet(structure, prev_incar=None, mode='line', nedos=2001, dedos=0.005, reciprocal_density=100, sym_prec=0.1, kpoints_line_density=20, optics=False, copy_chgcar=True, nbands_factor=1.2, small_gap_multiply=None, \*\*kwargs)
Bases: `DictSet`

Init a MPNonSCFSet. Typically, you would use the classmethod
from_prev_calc to initialize from a previous SCF run.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to compute


    * **prev_incar** (*Incar/string*) – Incar file from previous run.


    * **mode** (*str*) – Line, Uniform or Boltztrap mode supported.


    * **nedos** (*int*) – nedos parameter. Default to 2001.


    * **dedos** (*float*) – setting nedos=0 and uniform mode in from_prev_calc,
    an automatic nedos will be calculated using the total energy range
    divided by the energy step dedos


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal
    volume (defaults to 100)


    * **sym_prec** (*float*) – Symmetry precision (for Uniform mode).


    * **kpoints_line_density** (*int*) – Line density for Line mode.


    * **optics** (*bool*) – whether to add dielectric function


    * **copy_chgcar** – Whether to copy the old CHGCAR when starting from a
    previous calculation.


    * **nbands_factor** (*float*) – Multiplicative factor for NBANDS when starting
    from a previous calculation. Choose a higher number if you are
    doing an LOPTICS calculation.


    * **small_gap_multiply** (*[**float**, **float**]*) – When starting from a previous
    calculation, if the gap is less than 1st index, multiply the default
    reciprocal_density by the 2nd index.


    * **\*\*kwargs** – kwargs supported by MPRelaxSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for NonSCF calculations from a
directory of previous static VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPNonSCFSet, other than structure,
    prev_incar and prev_chgcar which are determined from the
    prev_calc_dir.



#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints(_: [Kpoints](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints) | Non_ )
Kpoints


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPRelaxSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, \*\*kwargs)
Bases: `DictSet`

Implementation of VaspInputSet utilizing parameters in the public
Materials Project. Typically, the pseudopotentials chosen contain more
electrons than the MIT parameters, and the k-point grid is ~50% more dense.
The LDAUU parameters are also different due to the different PSPs used,
which result in different fitted values.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **\*\*kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'IBRION': 2, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'LASPH': True, 'LDAU': True, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 64}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPSOCSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, saxis: tuple[int, int, int] = (0, 0, 1), copy_chgcar=True, nbands_factor=1.2, reciprocal_density=100, small_gap_multiply=None, magmom=None, \*\*kwargs)
Bases: `MPStaticSet`

An input set for running spin-orbit coupling (SOC) calculations.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – the structure must have the ‘magmom’ site
    property and each magnetic moment value must have 3
    components. eg: `magmom = [[0,0,2], ...]`


    * **saxis** (*tuple*) – magnetic moment orientation


    * **copy_chgcar** – Whether to copy the old CHGCAR. Defaults to True.


    * **nbands_factor** (*float*) – Multiplicative factor for NBANDS. Choose a
    higher number if you are doing an LOPTICS calculation.


    * **reciprocal_density** (*int*) – density of k-mesh by reciprocal volume.


    * **small_gap_multiply** (*[**float**, **float**]*) – If the gap is less than
    1st index, multiply the default reciprocal_density by the 2nd
    index.


    * **magmom** (*list**[**list**[**float**]**]*) – Override for the structure magmoms.


    * **\*\*kwargs** – kwargs supported by MPStaticSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for SOC calculations from a
directory of previous static VASP run. SOC calc requires all 3
components for MAGMOM for each atom in the structure.


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPSOCSet, other than structure,
    prev_incar and prev_chgcar which are determined from the
    prev_calc_dir.



#### _property_ incar(_: Inca_ )
Incar


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPScanRelaxSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, bandgap: float = 0, bandgap_tol: float = 0.0001, \*\*kwargs)
Bases: `DictSet`

Class for writing a relaxation input set using the accurate and numerically
efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
(SCAN) metaGGA density functional.

### Notes

1. This functional is officially supported in VASP 6.0.0 and above. On older version,
source code may be obtained by contacting the authors of the referenced manuscript.
The original SCAN functional, available from VASP 5.4.3 onwards, maybe used instead
by passing user_incar_settings={“METAGGA”: “SCAN”} when instantiating this InputSet.
r2SCAN and SCAN are expected to yield very similar results.

2. Meta-GGA calculations require POTCAR files that include
information on the kinetic energy density of the core-electrons,
i.e. “PBE_52” or “PBE_54”. Make sure the POTCARs include the
following lines (see VASP wiki for more details):

> $ grep kinetic POTCAR
> kinetic energy-density
> mkinetic energy-density pseudized
> kinetic energy density (partial)

### References

James W. Furness, Aaron D. Kaplan, Jinliang Ning, John P. Perdew, and Jianwei Sun.
Accurate and Numerically Efficient r2SCAN Meta-Generalized Gradient Approximation.
The Journal of Physical Chemistry Letters 0, 11 DOI: 10.1021/acs.jpclett.0c02405


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **bandgap** (*float*) – Bandgap of the structure in eV. The bandgap is used to
    compute the appropriate k-point density and determine the
    smearing settings.

    Metallic systems (default, bandgap = 0) use a KSPACING value of 0.22
    and Methfessel-Paxton order 2 smearing (ISMEAR=2, SIGMA=0.2).

    Non-metallic systems (bandgap > 0) use the tetrahedron smearing
    method (ISMEAR=-5, SIGMA=0.05). The KSPACING value is
    calculated from the bandgap via Eqs. 25 and 29 of Wisesa, McGill,
    and Mueller [1] (see References). Note that if ‘user_incar_settings’
    or ‘user_kpoints_settings’ override KSPACING, the calculation from
    bandgap is not performed.



    * **bandgap_tol** (*float*) – Tolerance for determining if a system is metallic.
    If the bandgap is less than this value, the system is considered
    metallic. Defaults to 1e-4 (eV).


    * **vdw** (*str*) – set “rVV10” to enable SCAN+rVV10, which is a versatile
    van der Waals density functional by combing the SCAN functional
    with the rVV10 non-local correlation functional. rvv10 is the only
    dispersion correction available for SCAN at this time.


    * **\*\*kwargs** – Same as those supported by DictSet.


### References

[1] P. Wisesa, K.A. McGill, T. Mueller, Efficient generation of
generalized Monkhorst-Pack grids through the use of informatics,
Phys. Rev. B. 93 (2016) 1-10. doi:10.1103/PhysRevB.93.155109.


#### CONFIG(_ = {'INCAR': {'ALGO': 'ALL', 'EDIFF': 1e-05, 'EDIFFG': -0.02, 'ENAUG': 1360, 'ENCUT': 680, 'IBRION': 2, 'ISIF': 3, 'ISPIN': 2, 'LAECHG': True, 'LASPH': True, 'LCHARG': True, 'LELF': False, 'LMIXTAU': True, 'LORBIT': 11, 'LREAL': 'Auto', 'LVTOT': True, 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'METAGGA': 'R2SCAN', 'NELM': 200, 'NSW': 99, 'PREC': 'Accurate'}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Am': 'Am', 'Ar': 'Ar', 'As': 'As', 'At': 'At', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cf': 'Cf', 'Cl': 'Cl', 'Cm': 'Cm', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Fr': 'Fr_sv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Po': 'Po_d', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Ra': 'Ra_sv', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Rn': 'Rn', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_sv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_3', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_54'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _valid_potcars(_: Sequence[str] | Non_ _ = ('PBE_52', 'PBE_54'_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPScanStaticSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, bandgap=0, prev_incar=None, lepsilon=False, lcalcpol=False, \*\*kwargs)
Bases: `MPScanRelaxSet`

Creates input files for a static calculation using the accurate and numerically
efficient r2SCAN variant of the Strongly Constrained and Appropriately Normed
(SCAN) metaGGA functional.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure from previous run.


    * **bandgap** (*float*) – Bandgap of the structure in eV. The bandgap is used to
    compute the appropriate k-point density and determine the smearing settings.


    * **prev_incar** (*Incar*) – Incar file from previous run.


    * **lepsilon** (*bool*) – Whether to add static dielectric calculation


    * **lcalcpol** (*bool*) – Whether to turn on evaluation of the Berry phase approximations
    for electronic polarization.


    * **\*\*kwargs** – kwargs supported by MPScanRelaxSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for static calculations from a
directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPScanStaticSet, other than prev_incar
    which is determined from the prev_calc_dir.



#### _property_ incar(_: Inca_ )
Incar


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MPStaticSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, prev_incar=None, prev_kpoints=None, lepsilon=False, lcalcpol=False, reciprocal_density=100, small_gap_multiply=None, \*\*kwargs)
Bases: `DictSet`

Creates input files for a static calculation.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure from previous run.


    * **prev_incar** (*Incar*) – Incar file from previous run.


    * **prev_kpoints** ([*Kpoints*](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints)) – Kpoints from previous run.


    * **lepsilon** (*bool*) – Whether to add static dielectric calculation


    * **lcalcpol** (*bool*) – Whether to turn on evaluation of the Berry phase approximations
    for electronic polarization


    * **reciprocal_density** (*int*) – For static calculations, we usually set the
    reciprocal density by volume. This is a convenience arg to change
    that, rather than using user_kpoints_settings. Defaults to 100,
    which is ~50% more than that of standard relaxation calculations.


    * **small_gap_multiply** (*[**float**, **float**]*) – If the gap is less than
    1st index, multiply the default reciprocal_density by the 2nd
    index.


    * **\*\*kwargs** – kwargs supported by MPRelaxSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for static calculations from a
directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MPStaticSet, other than prev_incar
    and prev_structure and prev_kpoints which are determined from
    the prev_calc_dir.



#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints(_: [Kpoints](pymatgen.io.cp2k.md#pymatgen.io.cp2k.inputs.Kpoints) | Non_ )
Kpoints


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLElasticSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, potim: float = 0.015, \*\*kwargs)
Bases: `DictSet`

MVL denotes VASP input sets that are implemented by the Materials Virtual
Lab ([http://materialsvirtuallab.org](http://materialsvirtuallab.org)) for various research.

This input set is used to calculate elastic constants in VASP. It is used
in the following work:

```default
Z. Deng, Z. Wang, I.-H. Chu, J. Luo, S. P. Ong.
“Elastic Properties of Alkali Superionic Conductor Electrolytes
from First Principles Calculations”, J. Electrochem. Soc.
2016, 163(2), A67-A74. doi: 10.1149/2.0061602jes
```

To read the elastic constants, you may use the Outcar class which parses the
elastic constants.


* **Parameters**


    * **structure** (*pymatgen.Structure*) – Input structure.


    * **potim** (*float*) – POTIM parameter. The default of 0.015 is usually fine,
    but some structures may require a smaller step.


    * **kwargs** – Parameters supported by MPRelaxSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLGBSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, k_product=40, slab_mode=False, is_metal=True, \*\*kwargs)
Bases: `DictSet`

Class for writing a vasp input files for grain boundary calculations, slab
or bulk.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – provide the structure


    * **k_product** – Kpoint number \* length for a & b directions, also for c
    direction in bulk calculations. Default to 40.


    * **slab_mode** (*bool*) – Defaults to False. Use default (False) for a
    bulk supercell. Use True if you are performing calculations on a
    slab-like (i.e., surface) of the GB, for example, when you are
    calculating the work of separation.


    * **is_metal** (*bool*) – Defaults to True. This determines whether an ISMEAR of
    1 is used (for metals) or not (for insulators and semiconductors)
    by default. Note that it does *not* override user_incar_settings,
    which can be set by the user to be anything desired.


    * **\*\*kwargs** – Other kwargs supported by MPRelaxSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints()
k_product, default to 40, is kpoint number \* length for a & b
directions, also for c direction in bulk calculations
Automatic mesh & Gamma is the default setting.


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLGWSet(structure, prev_incar=None, nbands=None, reciprocal_density=100, mode='STATIC', copy_wavecar=True, nbands_factor=5, ncores=16, \*\*kwargs)
Bases: `DictSet`

MVL denotes VASP input sets that are implemented by the Materials Virtual
Lab ([http://materialsvirtuallab.org](http://materialsvirtuallab.org)) for various research. This is a
flexible input set for GW calculations.

Note that unlike all other input sets in this module, the PBE_54 series of
functional is set as the default. These have much improved performance for
GW calculations.

A typical sequence is mode=”STATIC” -> mode=”DIAG” -> mode=”GW” ->
mode=”BSE”. For all steps other than the first one (static), the
recommendation is to use from_prev_calculation on the preceding run in
the series.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **prev_incar** (*Incar/string*) – Incar file from previous run.


    * **mode** (*str*) – Supported modes are “STATIC” (default), “DIAG”, “GW”,
    and “BSE”.


    * **nbands** (*int*) – For subsequent calculations, it is generally
    recommended to perform NBANDS convergence starting from the
    NBANDS of the previous run for DIAG, and to use the exact same
    NBANDS for GW and BSE. This parameter is used by
    from_previous_calculation to set nband.


    * **copy_wavecar** – Whether to copy the old WAVECAR, WAVEDER and associated
    files when starting from a previous calculation.


    * **nbands_factor** (*int*) – Multiplicative factor for NBANDS when starting
    from a previous calculation. Only applies if mode==”DIAG”.
    Need to be tested for convergence.


    * **reciprocal_density** (*int*) – Density of k-mesh by reciprocal atom. Only
    applies if mode==”STATIC”. Defaults to 100.


    * **ncores** (*int*) – Numbers of cores used for the calculation. VASP will alter
    NBANDS if it was not dividable by ncores. Only applies if
    mode==”DIAG”.


    * **\*\*kwargs** – All kwargs supported by DictSet. Typically,
    user_incar_settings is a commonly used option.



#### CONFIG(_ = {'INCAR': {'ALGO': 'Normal', 'EDIFF': 1e-08, 'IBRION': -1, 'ICHARG': 1, 'ISMEAR': 0, 'ISPIN': 2, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': True, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 100, 'PREC': 'Accurate', 'SIGMA': 0.01}, 'KPOINTS': {'reciprocal_density': 100}, 'PARENT': 'VASPIncarBase', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag_sv_GW', 'Al': 'Al_GW', 'Ar': 'Ar_GW', 'As': 'As_GW', 'At': 'At_d_GW', 'Au': 'Au_sv_GW', 'B': 'B_GW', 'Ba': 'Ba_sv_GW', 'Be': 'Be_sv_GW', 'Bi': 'Bi_d_GW', 'Br': 'Br_GW', 'C': 'C_GW', 'Ca': 'Ca_sv_GW', 'Cd': 'Cd_sv_GW', 'Ce': 'Ce_GW', 'Cl': 'Cl_GW', 'Co': 'Co_sv_GW', 'Cr': 'Cr_sv_GW', 'Cs': 'Cs_sv_GW', 'Cu': 'Cu_sv_GW', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F_GW', 'Fe': 'Fe_sv_GW', 'Ga': 'Ga_d_GW', 'Gd': 'Gd', 'Ge': 'Ge_d_GW', 'H': 'H_GW', 'He': 'He_GW', 'Hf': 'Hf_sv_GW', 'Hg': 'Hg_sv_GW', 'Ho': 'Ho_3', 'I': 'I_GW', 'In': 'In_d_GW', 'Ir': 'Ir_sv_GW', 'K': 'K_sv_GW', 'Kr': 'Kr_GW', 'La': 'La_GW', 'Li': 'Li_sv_GW', 'Lu': 'Lu_3', 'Mg': 'Mg_sv_GW', 'Mn': 'Mn_sv_GW', 'Mo': 'Mo_sv_GW', 'N': 'N_GW', 'Na': 'Na_sv_GW', 'Nb': 'Nb_sv_GW', 'Nd': 'Nd_3', 'Ne': 'Ne_GW', 'Ni': 'Ni_sv_GW', 'Np': 'Np', 'O': 'O_GW', 'Os': 'Os_sv_GW', 'P': 'P_GW', 'Pa': 'Pa', 'Pb': 'Pb_d_GW', 'Pd': 'Pd_sv_GW', 'Pm': 'Pm_3', 'Po': 'Po_d_GW', 'Pr': 'Pr_3', 'Pt': 'Pt_sv_GW', 'Pu': 'Pu', 'Rb': 'Rb_sv_GW', 'Re': 'Re_sv_GW', 'Rh': 'Rh_sv_GW', 'Rn': 'Rn_d_GW', 'Ru': 'Ru_sv_GW', 'S': 'S_GW', 'Sb': 'Sb_d_GW', 'Sc': 'Sc_sv_GW', 'Se': 'Se_GW', 'Si': 'Si_GW', 'Sm': 'Sm_3', 'Sn': 'Sn_d_GW', 'Sr': 'Sr_sv_GW', 'Ta': 'Ta_sv_GW', 'Tb': 'Tb_3', 'Tc': 'Tc_sv_GW', 'Te': 'Te_GW', 'Th': 'Th', 'Ti': 'Ti_sv_GW', 'Tl': 'Tl_d_GW', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_sv_GW', 'W': 'W_sv_GW', 'Xe': 'Xe_GW', 'Y': 'Y_sv_GW', 'Yb': 'Yb_3', 'Zn': 'Zn_sv_GW', 'Zr': 'Zr_sv_GW'}, 'POTCAR_FUNCTIONAL': 'PBE_54'_ )

#### SUPPORTED_MODES(_ = ('DIAG', 'GW', 'STATIC', 'BSE'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, mode='DIAG', \*\*kwargs)
Generate a set of VASP input files for GW or BSE calculations from a
directory of previous Exact Diag VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – The directory contains the outputs(
    vasprun.xml of previous vasp run.


    * **mode** (*str*) – Supported modes are “STATIC”, “DIAG” (default), “GW”,
    and “BSE”.


    * **\*\*kwargs** – All kwargs supported by MVLGWSet, other than structure,
    prev_incar and mode, which are determined from the
    prev_calc_dir.



#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints(_: Kpoint_ )
Generate gamma center k-points mesh grid for GW calc,
which is requested by GW calculation.


#### override_from_prev_calc(prev_calc_dir='.')
Update the input set to include settings from a previous calculation.


* **Parameters**

    **prev_calc_dir** (*str*) – The path to the previous calculation directory.



* **Returns**

    The input set with the settings (structure, k-points, incar, etc)
    updated using the previous VASP run.



#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLNPTMDSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, start_temp: float = 0.0, end_temp: float = 300.0, nsteps: int = 1000, time_step: float = 2, spin_polarized=False, \*\*kwargs)
Bases: `MITMDSet`

Class for writing a vasp md run in NPT ensemble.

### Notes

To eliminate Pulay stress, the default ENCUT is set to a rather large
value of ENCUT, which is 1.5 \* ENMAX.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **start_temp** (*int*) – Starting temperature.


    * **end_temp** (*int*) – Final temperature.


    * **nsteps** (*int*) – Number of time steps for simulations. NSW parameter.


    * **time_step** (*int*) – The time step for the simulation. The POTIM
    parameter. Defaults to 2fs.


    * **spin_polarized** (*bool*) – Whether to do spin polarized calculations.
    The ISPIN parameter. Defaults to False.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ incar(_: Inca_ )
Special processing of incar.


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLRelax52Set(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, \*\*kwargs)
Bases: `DictSet`

Implementation of VaspInputSet utilizing the public Materials Project
parameters for INCAR & KPOINTS and VASP’s recommended PAW potentials for
POTCAR.

Keynotes from VASP manual:


    1. Recommended potentials for calculations using vasp.5.2+


    2. If dimers with short bonds are present in the compound (O2, CO,

        N2, F2, P2, S2, Cl2), it is recommended to use the h potentials.
        Specifically, C_h, O_h, N_h, F_h, P_h, S_h, Cl_h


    3. Released on Oct 28, 2018 by VASP. Please refer to VASP

        Manual 1.2, 1.3 & 10.2.1 for more details.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **user_potcar_functional** (*str*) – choose from “PBE_52” and “PBE_54”.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'FAST', 'EDIFF_PER_ATOM': 5e-05, 'ENCUT': 520, 'IBRION': 2, 'ICHARG': 1, 'ISIF': 3, 'ISMEAR': -5, 'ISPIN': 2, 'LDAU': True, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUPRINT': 1, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LORBIT': 11, 'LREAL': 'AUTO', 'LWAVE': False, 'NELM': 100, 'NSW': 99, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'KPOINTS': {'reciprocal_density': 64}, 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Am': 'Am', 'Ar': 'Ar', 'As': 'As', 'At': 'At_d', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be', 'Bi': 'Bi_d', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Cm': 'Cm', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu_2', 'F': 'F', 'Fe': 'Fe', 'Fr': 'Fr_sv', 'Ga': 'Ga_d', 'Gd': 'Gd_3', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg', 'Mn': 'Mn_pv', 'Mo': 'Mo_sv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_sv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni', 'Np': 'Np', 'O': 'O', 'Os': 'Os', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Po': 'Po_d', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Ra': 'Ra_sv', 'Rb': 'Rb_sv', 'Re': 'Re', 'Rh': 'Rh_pv', 'Rn': 'Rn', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_sv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_sv', 'W': 'W_pv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_52'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _valid_potcars(_: Sequence[str] | Non_ _ = ('PBE_52', 'PBE_54'_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLScanRelaxSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, \*\*kwargs)
Bases: `DictSet`

Class for writing a relax input set using Strongly Constrained and
Appropriately Normed (SCAN) semilocal density functional.

### Notes


1. This functional is only available from VASP.5.4.3 upwards.

2. Meta-GGA calculations require POTCAR files that include
information on the kinetic energy density of the core-electrons,
i.e. “PBE_52” or “PBE_54”. Make sure the POTCAR including the
following lines (see VASP wiki for more details):

> $ grep kinetic POTCAR
> kinetic energy-density
> mkinetic energy-density pseudized
> kinetic energy density (partial)


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **vdw** (*str*) – set “rVV10” to enable SCAN+rVV10, which is a versatile
    van der Waals density functional by combing the SCAN functional
    with the rVV10 non-local correlation functional.


    * **\*\*kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _valid_potcars(_: Sequence[str] | Non_ _ = ('PBE_52', 'PBE_54'_ )

#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MVLSlabSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, k_product=50, bulk=False, auto_dipole=False, set_mix=True, sort_structure=True, \*\*kwargs)
Bases: `DictSet`

Class for writing a set of slab vasp runs,
including both slabs (along the c direction) and orient unit cells (bulk),
to ensure the same KPOINTS, POTCAR and INCAR criterion.


* **Parameters**


    * **structure** – Structure


    * **k_product** – default to 50, kpoint number \* length for a & b
    directions, also for c direction in bulk calculations


    * **bulk** –


    * **auto_dipole** –


    * **set_mix** –


    * **sort_structure** –


    * **kwargs** – Other kwargs supported by DictSet.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict(verbosity=2)

* **Parameters**

    **verbosity** – Verbosity of dict. E.g., whether to include Structure.



* **Returns**

    MSONable dict



#### _property_ incar(_: Inca_ )
Incar


#### _property_ kpoints()
k_product, default to 50, is kpoint number \* length for a & b

    directions, also for c direction in bulk calculations

Automatic mesh & Gamma is the default setting.


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ MatPESStaticSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None, xc_functional: Literal['R2SCAN', 'PBE', 'PBE+U'] = 'PBE', prev_incar: Incar | dict | None = None, \*\*kwargs: Any)
Bases: `DictSet`

Creates input files for a MatPES static calculation.

The goal of MatPES is to generate potential energy surface data. This is a distinctly different
from the objectives of the MP static calculations, which aims to obtain primarily accurate
energies and also electronic structure (DOS). For PES data, force accuracy (and to some extent,
stress accuracy) is of paramount importance.

The default POTCAR versions have been updated to PBE_54 from the old PBE set used in the
MPStaticSet. However, **U values** are still based on PBE. The implicit assumption here is that
the PBE_54 and PBE POTCARs are sufficiently similar that the U values fitted to the old PBE
functional still applies.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The Structure to create inputs for. If None, the input set is initialized without
    a Structure but one must be set separately before the inputs are generated.


    * **xc_functional** (*'R2SCAN'**|**'PBE'*) – Exchange-correlation functional to use. Defaults to ‘PBE’.


    * **prev_incar** (*Incar** | **dict*) – Incar file from previous run. Default settings of MatPESStaticSet
    are prioritized over inputs from previous runs. Defaults to None.


    * **\*\*kwargs** – Same as those supported by DictSet.



#### CONFIG(_ = {'INCAR': {'ALGO': 'Normal', 'EDIFF': 1e-05, 'ENAUG': 1360, 'ENCUT': 680, 'GGA': 'PE', 'ISMEAR': 0, 'ISPIN': 2, 'KSPACING': 0.22, 'LAECHG': True, 'LASPH': True, 'LCHARG': True, 'LDAU': False, 'LDAUJ': {'F': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}, 'O': {'Co': 0, 'Cr': 0, 'Fe': 0, 'Mn': 0, 'Mo': 0, 'Ni': 0, 'V': 0, 'W': 0}}, 'LDAUL': {'F': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}, 'O': {'Co': 2, 'Cr': 2, 'Fe': 2, 'Mn': 2, 'Mo': 2, 'Ni': 2, 'V': 2, 'W': 2}}, 'LDAUTYPE': 2, 'LDAUU': {'F': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}, 'O': {'Co': 3.32, 'Cr': 3.7, 'Fe': 5.3, 'Mn': 3.9, 'Mo': 4.38, 'Ni': 6.2, 'V': 3.25, 'W': 6.2}}, 'LMIXTAU': True, 'LORBIT': 11, 'LREAL': False, 'LWAVE': False, 'MAGMOM': {'Ce': 5, 'Ce3+': 1, 'Co': 0.6, 'Co3+': 0.6, 'Co4+': 1, 'Cr': 5, 'Dy3+': 5, 'Er3+': 3, 'Eu': 10, 'Eu2+': 7, 'Eu3+': 6, 'Fe': 5, 'Gd3+': 7, 'Ho3+': 4, 'La3+': 0.6, 'Lu3+': 0.6, 'Mn': 5, 'Mn3+': 4, 'Mn4+': 3, 'Mo': 5, 'Nd3+': 3, 'Ni': 5, 'Pm3+': 4, 'Pr3+': 2, 'Sm3+': 5, 'Tb3+': 6, 'Tm3+': 2, 'V': 5, 'W': 5, 'Yb3+': 1}, 'NELM': 200, 'NSW': 0, 'PREC': 'Accurate', 'SIGMA': 0.05}, 'PARENT': 'PBE54Base', 'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Am': 'Am', 'Ar': 'Ar', 'As': 'As', 'At': 'At', 'Au': 'Au', 'B': 'B', 'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv', 'Cd': 'Cd', 'Ce': 'Ce', 'Cf': 'Cf', 'Cl': 'Cl', 'Cm': 'Cm', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu_pv', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Fr': 'Fr_sv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He', 'Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Po': 'Po_d', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Ra': 'Ra_sv', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Rn': 'Rn', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_pv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W_sv', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_3', 'Zn': 'Zn', 'Zr': 'Zr_sv'}, 'POTCAR_FUNCTIONAL': 'PBE_54'_ )

#### INHERITED_INCAR_PARAMS(_ = ('LPEAD', 'NGX', 'NGY', 'NGZ', 'SYMPREC', 'IMIX', 'LMAXMIX', 'KGAMMA', 'ISYM', 'NCORE', 'NPAR', 'NELMIN', 'IOPT', 'NBANDS', 'KPAR', 'AMIN', 'NELMDL', 'BMIX', 'AMIX_MAG', 'BMIX_MAG'_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_prev_calc(prev_calc_dir, \*\*kwargs)
Generate a set of VASP input files for static calculations from a directory of previous VASP run.


* **Parameters**


    * **prev_calc_dir** (*str*) – Directory containing the outputs(
    vasprun.xml and OUTCAR) of previous vasp run.


    * **\*\*kwargs** – All kwargs supported by MatPESStaticSet, other than prev_incar
    and prev_structure and prev_kpoints which are determined from
    the prev_calc_dir.



#### _property_ incar(_: Inca_ )
Incar


#### user_potcar_functional(_: UserPotcarFunctiona_ )

### _class_ VaspInputSet()
Bases: `MSONable`

Base class representing a set of VASP input parameters with a structure
supplied as init parameters. Typically, you should not inherit from this
class. Start from DictSet or MPRelaxSet or MITRelaxSet.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _valid_potcars(_: Sequence[str] | Non_ _ = Non_ )

#### as_dict(verbosity=2)

* **Parameters**


    * **verbosity** – Verbosity for generated dict. If 1, structure is


    * **excluded.** –



* **Returns**

    MSONable dict



#### get_vasp_input(structure=None)

* **Returns**

    VaspInput.



#### _abstract property_ incar()
Incar object.


#### _abstract property_ kpoints()
Kpoints object.


#### _abstract property_ poscar()
Poscar object.


#### _property_ potcar(_: Potca_ )
Potcar object.


#### _property_ potcar_symbols()
List of POTCAR symbols.


#### write_input(output_dir: str, make_dir_if_not_present: bool = True, include_cif: bool = False, potcar_spec: bool = False, zip_output: bool = False)
Writes a set of VASP input to a directory.


* **Parameters**


    * **output_dir** (*str*) – Directory to output the VASP input files


    * **make_dir_if_not_present** (*bool*) – Set to True if you want the
    directory (and the whole path) to be created if it is not
    present.


    * **include_cif** (*bool*) – Whether to write a CIF file in the output
    directory for easier opening by VESTA.


    * **potcar_spec** (*bool*) – Instead of writing the POTCAR, write a “POTCAR.spec”.
    This is intended to help sharing an input set with people who might
    not have a license to specific Potcar files. Given a “POTCAR.spec”,
    the specific POTCAR file can be re-generated using pymatgen with the
    “generate_potcar” function in the pymatgen CLI.


    * **zip_output** (*bool*) – If True, output will be zipped into a file with the
    same name as the InputSet (e.g., MPStaticSet.zip)



### _load_yaml_config(fname)

### batch_write_input(structures, vasp_input_set=<class 'pymatgen.io.vasp.sets.MPRelaxSet'>, output_dir='.', make_dir_if_not_present=True, subfolder=None, sanitize=False, include_cif=False, potcar_spec=False, zip_output=False, \*\*kwargs)
Batch write vasp input for a sequence of structures to
output_dir, following the format output_dir/{group}/{formula}_{number}.


* **Parameters**


    * **structures** (*[*[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)*]*) – Sequence of Structures.


    * **vasp_input_set** (*VaspInputSet*) – VaspInputSet class that creates
    vasp input files from structures. Note that a class should be
    supplied. Defaults to MPRelaxSet.


    * **output_dir** (*str*) – Directory to output files. Defaults to current
    directory “.”.


    * **make_dir_if_not_present** (*bool*) – Create the directory if not present.
    Defaults to True.


    * **subfolder** (*callable*) – Function to create subdirectory name from
    structure. Defaults to simply “formula_count”.


    * **sanitize** (*bool*) – Boolean indicating whether to sanitize the
    structure before writing the VASP input files. Sanitized output
    are generally easier for viewing and certain forms of analysis.
    Defaults to False.


    * **include_cif** (*bool*) – Whether to output a CIF as well. CIF files are
    generally better supported in visualization programs.


    * **potcar_spec** (*bool*) – Instead of writing the POTCAR, write a “POTCAR.spec”.
    This is intended to help sharing an input set with people who might
    not have a license to specific Potcar files. Given a “POTCAR.spec”,
    the specific POTCAR file can be re-generated using pymatgen with the
    “generate_potcar” function in the pymatgen CLI.


    * **zip_output** (*bool*) – If True, output will be zipped into a file with the
    same name as the InputSet (e.g., MPStaticSet.zip)


    * **\*\*kwargs** – Additional kwargs are passed to the vasp_input_set class
    in addition to structure.



### get_structure_from_prev_run(vasprun, outcar=None)
Process structure from previous run.


* **Parameters**


    * **vasprun** (*Vasprun*) – Vasprun that contains the final structure
    from previous run.


    * **outcar** (*Outcar*) – Outcar that contains the magnetization info from
    previous run.



* **Returns**

    Returns the magmom-decorated structure that can be passed to get
    VASP input files, e.g. get_kpoints.



### get_valid_magmom_struct(structure, inplace=True, spin_mode='auto')
Make sure that the structure has valid magmoms based on the kind of calculation
Fill in missing Magmom values.


* **Parameters**


    * **structure** – The input structure


    * **inplace** – True - edit the magmom of the input structurel; False - return new structure


    * **spin_mode** – “scalar”/”vector”/”none”/”auto” only first letter (s/v/n) is needed.
    dictates how the spin configuration will be determined.


        * auto: read the existing magmom values and decide


        * scalar: use a single scalar value (for spin up/down)


        * vector: use a vector value for spin-orbit systems


        * none: Remove all the magmom information




* **Returns**

    New structure if inplace is False



### get_vasprun_outcar(path, parse_dos=True, parse_eigen=True)

* **Parameters**


    * **path** – Path to get the vasprun.xml and OUTCAR.


    * **parse_dos** – Whether to parse dos. Defaults to True.


    * **parse_eigen** – Whether to parse eigenvalue. Defaults to True.



### next_num_with_prime_factors(n: int, max_prime_factor: int, must_inc_2: bool = True)
Return the next number greater than or equal to n that only has the desired prime factors.


* **Parameters**


    * **n** (*int*) – Initial guess at the grid density


    * **max_prime_factor** (*int*) – the maximum prime factor


    * **must_inc_2** (*bool*) – 2 must be a prime factor of the result



* **Returns**

    first product of the prime_factors that is >= n



* **Return type**

    int



### primes_less_than(max_val: int)
Get the primes less than or equal to the max value.


### standardize_structure(structure, sym_prec=0.1, international_monoclinic=True)
Get the symmetrically standardized structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The structure.


    * **sym_prec** (*float*) – Tolerance for symmetry finding for standardization.


    * **international_monoclinic** (*bool*) – Whether to use international
    convention (vs Curtarolo) for monoclinic. Defaults True.



* **Returns**

    The symmetrized structure.