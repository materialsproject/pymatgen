---
layout: default
title: pymatgen.io.vasp.inputs.md
nav_exclude: true
---

# pymatgen.io.vasp.inputs module

Classes for reading/manipulating/writing VASP input files. All major VASP input
files.


### _exception_ pymatgen.io.vasp.inputs.BadIncarWarning()
Bases: `UserWarning`

Warning class for bad Incar parameters.


### _class_ pymatgen.io.vasp.inputs.Incar(params: dict[str, Any] | None = None)
Bases: `dict`, `MSONable`

INCAR object for reading and writing INCAR files. Essentially consists of
a dictionary with some helper functions.

Creates an Incar object.


* **Parameters**

    **params** (*dict*) – A set of input parameters as a dictionary.



#### as_dict()

* **Returns**

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


#### get_string(sort_keys: bool = False, pretty: bool = False)
Returns a string representation of the INCAR. The reason why this
method is different from the __str__ method is to provide options for
pretty printing.


* **Parameters**


    * **sort_keys** (*bool*) – Set to True to sort the INCAR parameters
    alphabetically. Defaults to False.


    * **pretty** (*bool*) – Set to True for pretty aligned output. Defaults
    to False.



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



### _class_ pymatgen.io.vasp.inputs.Kpoints(comment: str = 'Default gamma', num_kpts: int = 0, style: KpointsSupportedModes = KpointsSupportedModes.Gamma, kpts: Sequence[float | int | Sequence] = ((1, 1, 1),), kpts_shift: Vector3D = (0, 0, 0), kpts_weights=None, coord_type=None, labels=None, tet_number: int = 0, tet_weight: float = 0, tet_connections=None)
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

* **Returns**

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



#### _static_ automatic_density(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), kppa: float, force_gamma: bool = False)
Returns an automatic Kpoint object based on a structure and a kpoint
density. Uses Gamma centered meshes for hexagonal cells and
Monkhorst-Pack grids otherwise.

Algorithm:

    Uses a simple approach scaling the number of divisions along each
    reciprocal lattice vector proportional to its length.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **kppa** (*float*) – Grid density


    * **force_gamma** (*bool*) – Force a gamma centered mesh (default is to
    use gamma only for hexagonal cells or odd meshes)



* **Returns**

    Kpoints



#### _static_ automatic_density_by_lengths(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), length_densities: Sequence[float], force_gamma: bool = False)
Returns an automatic Kpoint object based on a structure and a k-point
density normalized by lattice constants.

Algorithm:

    For a given dimension, the # of k-points is chosen as
    length_density = # of kpoints \* lattice constant, e.g. [50.0, 50.0, 1.0] would
    have k-points of 50/a x 50/b x 1/c.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **length_densities** (*list**[**floats**]*) – Defines the density of k-points in each


    * **dimension** –


    * **[****50.0** (*e.g.*) –


    * **50.0** –


    * **1.0****]****.** –


    * **force_gamma** (*bool*) – Force a gamma centered mesh



* **Returns**

    Kpoints



#### _static_ automatic_density_by_vol(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), kppvol: int, force_gamma: bool = False)
Returns an automatic Kpoint object based on a structure and a kpoint
density per inverse Angstrom^3 of reciprocal cell.

Algorithm:

    Same as automatic_density()


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **kppvol** (*int*) – Grid density per Angstrom^(-3) of reciprocal cell


    * **force_gamma** (*bool*) – Force a gamma centered mesh



* **Returns**

    Kpoints



#### _static_ automatic_gamma_density(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), kppa: float)
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



#### _property_ style()
Style for kpoint generation. One of Kpoints_supported_modes
enum.


* **Type**

    return



#### supported_modes()
alias of `KpointsSupportedModes`


#### write_file(filename)
Write Kpoints to a file.


* **Parameters**

    **filename** (*str*) – Filename to write to.



### _class_ pymatgen.io.vasp.inputs.KpointsSupportedModes(value)
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


### _class_ pymatgen.io.vasp.inputs.Orbital(n, l, j, E, occ)
Bases: `tuple`

Create new instance of Orbital(n, l, j, E, occ)


#### E()
Alias for field number 3


#### j()
Alias for field number 2


#### l()
Alias for field number 1


#### n()
Alias for field number 0


#### occ()
Alias for field number 4


### _class_ pymatgen.io.vasp.inputs.OrbitalDescription(l, E, Type, Rcut, Type2, Rcut2)
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


#### l()
Alias for field number 0


### _class_ pymatgen.io.vasp.inputs.Poscar(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), comment: str | None = None, selective_dynamics: ArrayLike | None = None, true_names: bool = True, velocities: ArrayLike | None = None, predictor_corrector: ArrayLike | None = None, predictor_corrector_preamble: str | None = None, sort_structure: bool = False)
Bases: `MSONable`

Object for representing the data in a POSCAR or CONTCAR file.
Please note that this current implementation. Most attributes can be set
directly.


#### structure()
Associated Structure.


#### comment()
Optional comment string.


#### true_names()
Boolean indication whether Poscar contains actual real names parsed
from either a POTCAR or the POSCAR itself.


#### selective_dynamics()
Selective dynamics attribute for each site if available. A Nx3 array of
booleans.


#### velocities()
Velocities for each site (typically read in from a CONTCAR). A Nx3
array of floats.


#### predictor_corrector()
Predictor corrector coordinates and derivatives for each site; i.e.
a list of three 1x3 arrays for each site (typically read in from a MD
CONTCAR).


#### predictor_corrector_preamble()
Predictor corrector preamble contains the predictor-corrector key,
POTIM, and thermostat parameters that precede the site-specic predictor
corrector data in MD CONTCAR


#### temperature()
Temperature of velocity Maxwell-Boltzmann initialization. Initialized
to -1 (MB hasn”t been performed).


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object.


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

* **Returns**

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



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_string(direct: bool = True, vasp4_compatible: bool = False, significant_figures: int = 16)
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



#### _property_ natoms()
Sequence of number of sites of each type associated with the Poscar.
Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.


#### _property_ predictor_corrector()
Predictor corrector in Poscar.


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


### _class_ pymatgen.io.vasp.inputs.Potcar(symbols=None, functional=None, sym_potcar_map=None)
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

* **Returns**

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



### _class_ pymatgen.io.vasp.inputs.PotcarSingle(data, symbol=None)
Bases: `object`

Object for a **single** POTCAR. The builder assumes the POTCAR contains
the complete untouched data in “data” as a string and a dict of keywords.


#### data()
POTCAR data as a string.


#### keywords()
Keywords parsed from the POTCAR as a dict. All keywords are also
accessible as attributes in themselves. E.g., potcar.enmax,
potcar.encut, etc.

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


* **Type**

    return



#### _property_ element()
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



#### _property_ functional()
Functional associated with PotcarSingle.


* **Type**

    return



#### _property_ functional_class()
Functional class associated with PotcarSingle.


* **Type**

    return



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
    potcar_functionals (List): List of potcar functionals associated with

    > the PotcarSingle




* **Return type**

    symbol (List)



#### _property_ nelectrons()
Number of electrons


* **Type**

    return



#### parse_functions(_ = {'COPYR': <function _parse_string>, 'DEXC': <function _parse_float>, 'EATOM': <function _parse_float>, 'EAUG': <function _parse_float>, 'EMMIN': <function _parse_float>, 'ENMAX': <function _parse_float>, 'ENMIN': <function _parse_float>, 'GGA': <function _parse_list>, 'ICORE': <function _parse_int>, 'IUNSCR': <function _parse_int>, 'LCOR': <function _parse_bool>, 'LEXCH': <function _parse_string>, 'LPAW': <function _parse_bool>, 'LULTRA': <function _parse_bool>, 'LUNSCR': <function _parse_bool>, 'NDATA': <function _parse_int>, 'POMASS': <function _parse_float>, 'QCUT': <function _parse_float>, 'QGAM': <function _parse_float>, 'RAUG': <function _parse_float>, 'RCLOC': <function _parse_float>, 'RCORE': <function _parse_float>, 'RDEP': <function _parse_float>, 'RDEPT': <function _parse_float>, 'RMAX': <function _parse_float>, 'RPACOR': <function _parse_float>, 'RRKJ': <function _parse_list>, 'RWIGS': <function _parse_float>, 'SHA256': <function _parse_string>, 'STEP': <function _parse_list>, 'TITEL': <function _parse_string>, 'VRHFIN': <function _parse_string>, 'ZVAL': <function _parse_float>_ )

#### _property_ potential_type(_: st_ )
Type of PSP. E.g., US, PAW, etc.


* **Type**

    return



#### _property_ symbol()
The POTCAR symbol, e.g. W_pv


* **Type**

    return



#### verify_potcar()
Attempts to verify the integrity of the POTCAR data.

This method checks the whole file (removing only the SHA256
metadata) against the SHA256 hash in the header if this is found.
If no SHA256 hash is found in the file, the file hash (md5 hash of the
whole file) is checked against all POTCAR file hashes known to pymatgen.

## Returns:

(bool, bool)

    has_sh256 and passed_hash_check are returned.


#### write_file(filename: str)
Writes PotcarSingle to a file.
:param filename: Filename.


### _exception_ pymatgen.io.vasp.inputs.UnknownPotcarWarning()
Bases: `UserWarning`

Warning raised when POTCAR hashes do not pass validation.


### _class_ pymatgen.io.vasp.inputs.VaspInput(incar, kpoints, poscar, potcar, optional_files=None, \*\*kwargs)
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

* **Returns**

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