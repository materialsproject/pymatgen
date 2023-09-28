---
layout: default
title: pymatgen.io.cp2k.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.cp2k package

Module for CP2K input/output parsing as well as sets for standard calculations.


## pymatgen.io.cp2k.inputs module

This module defines the building blocks of a CP2K input file. The cp2k input structure is
essentially a collection of “sections” which are similar to dictionary objects that activate
modules of the cp2k executable, and then “keywords” which adjust variables inside of those
modules. For example, FORCE_EVAL section will activate CP2K’s ability to calculate forces,
and inside FORCE_EVAL, the Keyword “METHOD can be set to “QS” to set the method of force
evaluation to be the quickstep (DFT) module.

A quick overview of the module:

– Section class defines the basis of Cp2k input and contains methods for manipulating these

    objects similarly to Dicts.

– Keyword class defines the keywords used inside of Section objects that changes variables in

    Cp2k programs.

– SectionList and KeywordList classes are lists of Section and Keyword objects that have

    the same dictionary key. This deals with repeated sections and keywords.

– Cp2kInput class is special instantiation of Section that is used to represent the full cp2k

    calculation input.

– The rest of the classes are children of Section intended to make initialization of common

    sections easier.


### _class_ AtomicMetadata(info: BasisInfo | PotentialInfo | None = None, element: Element | None = None, potential: Literal['All Electron', 'Pseudopotential'] | None = None, name: str | None = None, alias_names: list = <factory>, filename: str | None = None, version: str | None = None)
Bases: `MSONable`

Metadata for basis sets and potentials in cp2k.


#### info()
Info about this object


* **Type**

    BasisInfo | PotentialInfo | None



#### element()
Element for this object


* **Type**

    [Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | None



#### potential()
The potential for this object


* **Type**

    Literal[‘All Electron’, ‘Pseudopotential’] | None



#### name()
Name of the object


* **Type**

    str | None



#### alias_names()
Optional aliases


* **Type**

    list



#### filename()
Name of the file containing this object


* **Type**

    str | None



#### version()
Version


* **Type**

    str | None



#### alias_names(_: lis_ )

#### element(_: [Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | Non_ _ = Non_ )

#### filename(_: str | Non_ _ = Non_ )

#### get_hash()
Get a hash of this object.


#### get_str()
Get string representation.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### info(_: BasisInfo | PotentialInfo | Non_ _ = Non_ )

#### name(_: str | Non_ _ = Non_ )

#### potential(_: Literal['All Electron', 'Pseudopotential'] | Non_ _ = Non_ )

#### softmatch(other)
Soft matching to see if a desired basis/potential matches requirements.

Does soft matching on the “info” attribute first. Then soft matches against the
element and name/aliases.


#### version(_: str | Non_ _ = Non_ )

### _class_ Band_Structure(kpoint_sets: Sequence[Kpoint_Set], filename: str = 'BAND.bs', added_mos: int = -1, keywords: dict | None = None, subsections: dict | None = None)
Bases: `Section`

Specifies high symmetry paths for outputting the band structure in CP2K.


* **Parameters**


    * **kpoint_sets** – Sequence of Kpoint_Set objects for the band structure calculation.


    * **filename** – Filename for the band structure output


    * **added_mos** – Added (unoccupied) molecular orbitals for the calculation.


    * **keywords** – additional keywords


    * **subsections** – additional subsections.



#### _static_ from_kpoints(kpoints: [Kpoints](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints), kpoints_line_density=20)
Initialize band structure section from a line-mode Kpoint object.


* **Parameters**


    * **kpoints** – a kpoint object from the vasp module, which was constructed in line mode


    * **kpoints_line_density** – Number of kpoints along each path



### _class_ BasisFile(objects: Sequence | None = None)
Bases: `DataFile`

Data file for basis sets only.


#### _classmethod_ from_str(string)
Initialize from a string representation.


### _class_ BasisInfo(electrons: int | None = None, core: int | None = None, valence: int | None = None, polarization: int | None = None, diffuse: int | None = None, cc: bool | None = False, pc: bool | None = False, sr: bool | None = False, molopt: bool | None = False, admm: bool | None = False, lri: bool | None = False, contracted: bool | None = None, xc: str | None = None)
Bases: `MSONable`

Summary info about a basis set.


#### electrons()
Number of electrons


* **Type**

    int | None



#### core()
Number of basis functions per core electron


* **Type**

    int | None



#### valence()
Number of basis functions per valence electron OR number of exp if it
is a FIT formatted admm basis


* **Type**

    int | None



#### polarization()
Number of polarization functions


* **Type**

    int | None



#### diffuse()
Number of added, diffuse/augmentation functions


* **Type**

    int | None



#### cc()
Correlation consistent


* **Type**

    bool | None



#### pc()
Polarization consistent


* **Type**

    bool | None



#### sr()
Short-range optimized


* **Type**

    bool | None



#### molopt()
Optimized for molecules/solids


* **Type**

    bool | None



#### admm()
Whether this is an auxiliary basis set for ADMM


* **Type**

    bool | None



#### lri()
Whether this is a local resolution of identity auxiliary basis


* **Type**

    bool | None



#### contracted()
Whether this basis set is contracted


* **Type**

    bool | None



#### xc()
Exchange correlation functional used for creating this potential


* **Type**

    str | None



#### admm(_: bool | Non_ _ = Fals_ )

#### cc(_: bool | Non_ _ = Fals_ )

#### contracted(_: bool | Non_ _ = Non_ )

#### core(_: int | Non_ _ = Non_ )

#### diffuse(_: int | Non_ _ = Non_ )

#### electrons(_: int | Non_ _ = Non_ )

#### _classmethod_ from_str(string: str)
Get summary info from a string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### lri(_: bool | Non_ _ = Fals_ )

#### molopt(_: bool | Non_ _ = Fals_ )

#### pc(_: bool | Non_ _ = Fals_ )

#### polarization(_: int | Non_ _ = Non_ )

#### softmatch(other)
Soft matching to see if two basis sets match.

Will only match those attributes which *are* defined for this basis info object (one way checking)


#### sr(_: bool | Non_ _ = Fals_ )

#### valence(_: int | Non_ _ = Non_ )

#### xc(_: str | Non_ _ = Non_ )

### _class_ BrokenSymmetry(l_alpha: Sequence = (-1,), n_alpha: Sequence = (0,), nel_alpha: Sequence = (-1,), l_beta: Sequence = (-1,), n_beta: Sequence = (0,), nel_beta: Sequence = (-1,))
Bases: `Section`

Define the required atomic orbital occupation assigned in initialization
of the density matrix, by adding or subtracting electrons from specific
angular momentum channels. It works only with GUESS ATOMIC.

Initialize the broken symmetry section.


* **Parameters**


    * **l_alpha** – Angular momentum quantum number of the orbitals whose occupation is changed


    * **n_alpha** – Principal quantum number of the orbitals whose occupation is changed.
    Default is the first not occupied


    * **nel_alpha** – Orbital occupation change per angular momentum quantum number. In
    unrestricted calculations applied to spin alpha


    * **l_beta** – Same as L_alpha for beta channel


    * **n_beta** – Same as N_alpha for beta channel


    * **nel_beta** – Same as NEL_alpha for beta channel



#### _classmethod_ from_el(el, oxi_state=0, spin=0)
Create section from element, oxidation state, and spin.


### _class_ Cell(lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice), keywords: dict | None = None, \*\*kwargs)
Bases: `Section`

Defines the simulation cell (lattice).

Initialize the cell section.


* **Parameters**


    * **lattice** – pymatgen lattice object


    * **keywords** – additional keywords



### _class_ Coord(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), aliases: dict | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Specifies the coordinates of the atoms using a pymatgen structure object.


* **Parameters**


    * **structure** – Pymatgen structure object


    * **alias** (*bool*) – whether or not to identify the sites by Element + number so you can do
    things like assign unique magnetization do different elements.


    * **keywords** – additional keywords


    * **subsections** – additional subsections.



### _class_ Cp2kInput(name: str = 'CP2K_INPUT', subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Special instance of ‘Section’ class that is meant to represent the overall cp2k input.
Distinguishes itself from Section by overriding get_str() to not print this section’s
title and by implementing the file i/o.

Initialize Cp2kInput by calling the super.


#### _classmethod_ _from_dict(d)
Initialize from a dictionary.


#### _from_lines(lines)
Helper method, reads lines of text to get a Cp2kInput.


#### _static_ from_file(file: str)
Initialize from a file.


#### _classmethod_ from_lines(lines: list | tuple)
Helper method to read lines of file.


#### _static_ from_str(s: str)
Initialize from a string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_str()
Get string representation of the Cp2kInput.


#### write_file(input_filename: str = 'cp2k.inp', output_dir: str = '.', make_dir_if_not_present: bool = True)
Write input to a file.


* **Parameters**


    * **input_filename** (*str**, **optional*) – Defaults to “cp2k.inp”.


    * **output_dir** (*str**, **optional*) – Defaults to “.”.


    * **make_dir_if_not_present** (*bool**, **optional*) – Defaults to True.



### _class_ DOS(ndigits: int = 6, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the density of states.

Initialize the DOS section.


* **Parameters**


    * **ndigits** – how many digits of precision to print. As of 2022.1,
    this is necessary to not lose information.


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ DataFile(objects: Sequence | None = None)
Bases: `MSONable`

A data file for a cp2k calc.


#### _classmethod_ from_file(fn)
Load from a file.


#### _classmethod_ from_str()
Initialize from a string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_str()
Get string representation.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### objects(_: Sequence | Non_ _ = Non_ )

#### write_file(fn)
Write to a file.


### _class_ Davidson(new_prec_each: int = 20, preconditioner: str = 'FULL_SINGLE_INVERSE', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Parameters for davidson diagonalization.


* **Parameters**


    * **new_prec_each** (*int*) – How often to recalculate the preconditioner.


    * **preconditioner** (*str*) – Preconditioner to use.
    “FULL_ALL”: Most effective state selective preconditioner based on diagonalization,

    > requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap.
    > This preconditioner is recommended for almost all systems, except very large
    > systems where make_preconditioner would dominate the total computational cost.

    ”FULL_KINETIC”: Cholesky inversion of S and T, fast construction, robust, use for

        very large systems.

    ”FULL_SINGLE”: Based on H-eS diagonalization, not as good as FULL_ALL, but

        somewhat cheaper to apply.

    ”FULL_SINGLE_INVERSE”: Based on H-eS cholesky inversion, similar to FULL_SINGLE

        in preconditioning efficiency but cheaper to construct, might be somewhat
        less robust. Recommended for large systems.

    ”FULL_S_INVERSE”: Cholesky inversion of S, not as good as FULL_KINETIC,

        yet equally expensive.

    ”NONE”: skip preconditioning



    * **keywords** – additional keywords


    * **subsections** – additional subsections.



### _class_ Dft(basis_set_filenames: Iterable = ('BASIS_MOLOPT',), potential_filename='GTH_POTENTIALS', uks: bool = True, wfn_restart_file_name: str | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the DFT parameters in Cp2k.

Initialize the DFT section.


* **Parameters**


    * **basis_set_filenames** – Name of the file that contains the basis set
    information. Defaults to “BASIS_MOLOPT”.


    * **potential_filename** – Name of the file that contains the pseudopotential
    information. Defaults to “GTH_POTENTIALS”.


    * **uks** – Whether to run unrestricted Kohn Sham (spin polarized).
    Defaults to True.


    * **wfn_restart_file_name** – Defaults to None.


    * **keywords** – additional keywords to add.


    * **subsections** – Any subsections to initialize with. Defaults to None.



### _class_ DftPlusU(eps_u_ramping=1e-05, init_u_ramping_each_scf=False, l=-1, u_minus_j=0, u_ramping=0, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls DFT+U for an atom kind.

Initialize the DftPlusU section.


* **Parameters**


    * **eps_u_ramping** – (float) SCF convergence threshold at which to start ramping the U value


    * **init_u_ramping_each_scf** – (bool) Whether or not to do u_ramping each scf cycle


    * **l** – (int) angular moment of the orbital to apply the +U correction


    * **u_minus_j** – (float) the effective U parameter, Ueff = U-J


    * **u_ramping** – (float) stepwise amount to increase during ramping until u_minus_j is reached


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ Diagonalization(eps_adapt: float = 0, eps_iter: float = 1e-08, eps_jacobi: float = 0, jacobi_threshold: float = 1e-07, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls diagonalization settings (if using traditional diagonalization).

Initialize the diagonalization section.


### _class_ E_Density_Cube(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the electron density cube file.

Basic object representing a CP2K Section. Sections activate different parts of the
calculation. For example, FORCE_EVAL section will activate CP2K’s ability to calculate
forces.


* **Parameters**


    * **name** – The name of the section (must match name in CP2K)


    * **subsections** – A dictionary of subsections that are nested in this section.
    Format is {‘NAME’: Section(*args, \*\*kwargs). The name you chose for ‘NAME’
    to index that subsection does not \*have* to be the same as the section’s true name,
    but we recommend matching them. You can specify a blank dictionary if there are
    no subsections, or if you want to insert the subsections later.


    * **repeats** – Whether or not this section can be repeated. Most sections cannot.
    Default=False.


    * **description** – Description of this section for easier readability


    * **keywords** – the keywords to be set for this section. Each element should be a
    Keyword object. This can be more cumbersome than simply using kwargs for building
    a class in a script, but is more convenient for the class instantiations of CP2K
    sections (see below).


    * **section_parameters** – the section parameters for this section. Section parameters
    are specialized keywords that modify the behavior of the section overall. Most
    sections do not have section parameters, but some do. Unlike normal Keywords,
    these are specified as strings and not as Keyword objects.


    * **location** – the path to the section in the form ‘SECTION/SUBSECTION1/SUBSECTION3’,
    example for QS module: ‘FORCE_EVAL/DFT/QS’. This location is used to automatically
    determine if a subsection requires a supersection to be activated.


    * **verbose** – Controls how much is printed to Cp2k input files (Also see Keyword).
    If True, then a description of the section will be printed with it as a comment
    (if description is set). Default=True.


    * **alias** – An alias for this class to use in place of the name.


    * **keyword** (*kwargs are interpreted as*) –


    * **as** (*value pairs and added to the keywords array*) –


    * **objects** (*Keyword*) –



### _class_ ForceEval(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the calculation of energy and forces in Cp2k.

Initialize the ForceEval section.


### _class_ GaussianTypeOrbitalBasisSet(info: BasisInfo | None = None, element: Element | None = None, potential: Literal['All Electron', 'Pseudopotential'] | None = None, name: str | None = None, alias_names: list = <factory>, filename: str | None = None, version: str | None = None, nset: int | None = None, n: list[int] | None = None, lmax: list[int] | None = None, lmin: list[int] | None = None, nshell: list[dict[int, int]] | None = None, exponents: list[list[float]] | None = None, coefficients: list[dict[int, dict[int, dict[int, float]]]] | None = None)
Bases: `AtomicMetadata`

Model definition of a GTO basis set.


#### info()
Cardinality of this basis


* **Type**

    BasisInfo | None



#### nset()
Number of exponent sets


* **Type**

    int | None



#### n()
Principle quantum number for each set


* **Type**

    list[int] | None



#### lmax()
Maximum angular momentum quantum number for each set


* **Type**

    list[int] | None



#### lmin()
Minimum angular momentum quantum number for each set


* **Type**

    list[int] | None



#### nshell()
Number of shells for angular momentum l for each set


* **Type**

    list[dict[int, int]] | None



#### exponents()
Exponents for each set


* **Type**

    list[list[float]] | None



#### coefficients()
Contraction coefficients for each set. Dict[exp->l->shell]


* **Type**

    list[dict[int, dict[int, dict[int, float]]]] | None



#### coefficients(_: list[dict[int, dict[int, dict[int, float]]]] | Non_ _ = Non_ )

#### exponents(_: list[list[float]] | Non_ _ = Non_ )

#### _classmethod_ from_str(string: str)
Read from standard cp2k GTO formatted string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_keyword()
Convert basis to keyword object.


#### get_str()
Get standard cp2k GTO formatted string.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### info(_: BasisInfo | Non_ _ = Non_ )

#### lmax(_: list[int] | Non_ _ = Non_ )

#### lmin(_: list[int] | Non_ _ = Non_ )

#### n(_: list[int] | Non_ _ = Non_ )

#### _property_ nexp()
Number of exponents.


#### nset(_: int | Non_ _ = Non_ )

#### nshell(_: list[dict[int, int]] | Non_ _ = Non_ )

### _class_ Global(project_name: str = 'CP2K', run_type: str = 'ENERGY_FORCE', keywords: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls ‘global’ settings for cp2k execution such as RUN_TYPE and PROJECT_NAME.

Initialize the global section.


* **Parameters**


    * **project_name** – Defaults to “CP2K”.


    * **run_type** – what type of calculation to run


    * **keywords** – Additional keywords to add



### _class_ GthPotential(info: PotentialInfo = None, element: Element | None = None, potential: Literal['All Electron', 'Pseudopotential'] | None = None, name: str | None = None, alias_names: list = <factory>, filename: str | None = None, version: str | None = None, n_elecs: dict[int, int] | None = None, r_loc: float | None = None, nexp_ppl: int | None = None, c_exp_ppl: Sequence | None = None, radii: dict[int, float] | None = None, nprj: int | None = None, nprj_ppnl: dict[int, int] | None = None, hprj_ppnl: dict[int, dict[int, dict[int, float]]] | None = None)
Bases: `AtomicMetadata`

Representation of GTH-type (pseudo)potential.


#### info()
Info about this potential


* **Type**

    PotentialInfo



#### n_elecs()
Number of electrons for each quantum number


* **Type**

    dict[int, int] | None



#### r_loc()
Radius of local projectors


* **Type**

    float | None



#### nexp_ppl()
Number of the local pseudopotential functions


* **Type**

    int | None



#### c_exp_ppl()
Sequence = field(None, description=”Coefficients of the local pseudopotential functions


* **Type**

    Sequence | None



#### radii()
Radius of the nonlocal part for angular momentum quantum number l defined by the Gaussian
function exponents alpha_prj_ppnl


* **Type**

    dict[int, float] | None



#### nprj()
Number of projectors


* **Type**

    int | None



#### nprj_ppnl()
Number of the non-local projectors for the angular momentum quantum number


* **Type**

    dict[int, int] | None



#### hprj_ppnl()
Coefficients of the non-local projector functions. Coeff ij for ang momentum l


* **Type**

    dict[int, dict[int, dict[int, float]]] | None



### )()

#### c_exp_ppl(_: Sequence | Non_ _ = Non_ )

#### _classmethod_ from_section(section: Section)
Extract GTH-formatted string from a section and convert it to model.


#### _classmethod_ from_str(string)
Initialize model from a GTH formatted string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_keyword()
Get keyword object for the potential.


#### get_section()
Convert model to a GTH-formatted section object for input files.


#### get_str()
Convert model to a GTH-formatted string.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### hprj_ppnl(_: dict[int, dict[int, dict[int, float]]] | Non_ _ = Non_ )

#### n_elecs(_: dict[int, int] | Non_ _ = Non_ )

#### nexp_ppl(_: int | Non_ _ = Non_ )

#### nprj(_: int | Non_ _ = Non_ )

#### nprj_ppnl(_: dict[int, int] | Non_ _ = Non_ )

#### r_loc(_: float | Non_ _ = Non_ )

#### radii(_: dict[int, float] | Non_ _ = Non_ )

### _class_ Keyword(name: str, \*values, description: str | None = None, units: str | None = None, verbose: bool | None = True, repeats: bool | None = False)
Bases: `MSONable`

Class representing a keyword argument in CP2K. Within CP2K Sections, which activate features
of the CP2K code, the keywords are arguments that control the functionality of that feature.
For example, the section “FORCE_EVAL” activates the evaluation of forces/energies, but within
“FORCE_EVAL” the keyword “METHOD” controls whether or not this will be done with, say,
“Quickstep” (DFT) or “EIP” (empirical interatomic potential).

Initializes a keyword. These Keywords and the value passed to them are sometimes as simple
as KEYWORD VALUE, but can also be more elaborate such as KEYWORD [UNITS] VALUE1 VALUE2,
which is why this class exists: to handle many values and control easy printing to an
input file.


* **Parameters**


    * **name** – The name of this keyword. Must match an acceptable keyword from CP2K


    * **values** – All non-keyword arguments after ‘name’ are interpreted as the values to set for
    this keyword. i.e: KEYWORD ARG1 ARG2 would provide two values to the keyword.


    * **description** – The description for this keyword. This can make readability of
    input files easier for some. Default=None.


    * **units** – The units for this keyword. If not specified, CP2K default units will be
    used. Consult manual for default units. Default=None.


    * **verbose** – Whether the description should be printed with the string of this keyword


    * **repeats** – Whether or not this keyword may be repeated. Default=False.



#### as_dict()
Get a dictionary representation of the Keyword.


#### _classmethod_ from_dict(d)
Initialize from dictionary.


#### _static_ from_str(s)
Initialize from a string.

Keywords must be labeled with strings. If the postprocessor finds
that the keywords is a number, then None is return (used by
the file reader).


* **Returns**

    Keyword or None



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_str()
String representation of Keyword.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### verbosity(v)
Change the printing of this keyword’s description.


### _class_ KeywordList(keywords: Sequence[Keyword])
Bases: `MSONable`

Some keywords can be repeated, which makes accessing them via the normal dictionary
methods a little unnatural. This class deals with this by defining a collection
of same-named keywords that are accessed by one name.

Initializes a keyword list given a sequence of keywords.


* **Parameters**

    **keywords** – A list of keywords. Must all have the same name (case-insensitive)



#### append(item)
Append the keyword list.


#### extend(lst: Sequence[Keyword])
Extend the keyword list.


#### get_str(indent: int = 0)
String representation of Keyword.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### verbosity(verbosity)
Silence all keywords in keyword list.


### _class_ Kind(specie: str, alias: str | None = None, magnetization: float = 0.0, basis_set: GaussianTypeOrbitalBasisSet | str | None = 'GTH_BASIS', potential: GthPotential | str | None = 'GTH_POTENTIALS', ghost: bool = False, aux_basis: GaussianTypeOrbitalBasisSet | str | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Specifies the information for the different atom types being simulated.

Initialize a KIND section.


* **Parameters**


    * **specie** – Object representing the atom.


    * **alias** – Alias for the atom, can be used for specifying modifications
    to certain atoms but not all, e.g. Mg_1 and Mg_2 to force difference
    oxidation states on the two atoms.


    * **magnetization** – From the CP2K Manual: The magnetization used
    in the atomic initial guess. Adds magnetization/2 spin-alpha
    electrons and removes magnetization/2 spin-beta electrons.


    * **basis_set** – Basis set for this atom, accessible from the
    basis set file specified


    * **potential** – Pseudopotential for this atom, accessible from the
    potential file


    * **ghost** – Turn this into ghost atom (disable the potential)


    * **aux_basis** – Auxiliary basis to use with ADMM


    * **keywords** – additional keywords


    * **subsections** – additional subsections


    * **kwargs** – Additional kwargs to pass to Section()



### _class_ Kpoint_Set(npoints: int, kpoints: Iterable, units: str = 'B_VECTOR')
Bases: `Section`

Specifies a kpoint line to be calculated between special points.


* **Parameters**


    * **npoints** (*int*) – Number of kpoints along the line.


    * **kpoints** – A dictionary of {label: kpoint} kpoints defining the path


    * **units** (*str*) – Units for the kpoint coordinates.
    Options: “B_VECTOR” (reciprocal coordinates)

    > ”CART_ANGSTROM” (units of 2\*Pi/Angstrom)
    > “CART_BOHR” (units of 2\*Pi/Bohr).




### _class_ Kpoints(kpts: Sequence | Sequence[Sequence[int]], weights: Sequence | None = None, eps_geo: float = 1e-06, full_grid: bool = False, parallel_group_size: int = -1, scheme: str = 'MONKHORST-PACK', symmetry: bool = False, units: str = 'B_VECTOR', verbose: bool = False, wavefunctions: str = 'COMPLEX')
Bases: `Section`

Description of the k-points to use for the calculation.


* **Parameters**


    * **kpts** (*list**, **tuple*) – a 2D array for the kpoints of the form
    [(1,1,1),]. If len(kpts) == 1. Then it is taken as subdivisions
    for automatic kpoint scheme. If it has more entries, it is
    taken as manual entries for kpoints.


    * **weights** (*list**, **tuple*) – a weight for each kpoint. Default is to
    weigh each by 1


    * **eps_geo** (*float*) – tolerance for symmetry. Default=1e-6


    * **full_grid** (*bool*) – use full (not reduced) kpoint grid. Default=False.


    * **parallel_group_size** (*int*) – from cp2k manual: Number of processors
    to be used for a single kpoint. This number must divide the
    total number of processes. The number of groups must divide
    the total number of kpoints. Value=-1 (smallest possible
    number of processes per group, satisfying the constraints).
    Value=0 (all processes). Value=n (exactly n processes).
    Default=-1.


    * **scheme** (*str*) – kpoint generation scheme. Default=’Monkhorst-Pack’


    * **symmetry** (*bool*) – Use symmetry to reduce the number of kpoints.
    Default=False.


    * **units** (*str*) – Units for the kpoint coordinates (reciprocal coordinates
    or cartesian). Default=’B_VECTOR’ (reciprocal)


    * **verbose** (*bool*) – verbose output for kpoints. Default=False


    * **wavefunctions** (*str*) – Whether to use complex or real valued wavefunctions
    (if available). Default=’complex’.



#### _classmethod_ from_kpoints(kpoints: [Kpoints](pymatgen.io.vasp.md#pymatgen.io.vasp.inputs.Kpoints), structure=None)
Initialize the section from a Kpoints object (pymatgen.io.vasp.inputs). CP2K
does not have an automatic gamma-point constructor, so this is generally used
to get the number of divisions from a kpoint static constructor and then
build a Monkhorst-Pack grid, which is sufficient for gamma-recommended systems
so long as the grid is fine enough.


* **Parameters**


    * **kpoints** – A pymatgen kpoints object.


    * **structure** – Pymatgen structure object. Required for automatically performing
    symmetry analysis and reducing the kpoint grid.


    * **reduce** – whether or not to reduce the grid using symmetry. CP2K itself cannot
    do this automatically without spglib present at execution time.



### _class_ LDOS(index: int = 1, alias: str | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the LDOS (List-Density of states). i.e. projects onto specific atoms.

Initialize the LDOS section.


* **Parameters**


    * **index** – Index of the atom to project onto


    * **alias** – section alias


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ MO_Cubes(write_cube: bool = False, nhomo: int = 1, nlumo: int = 1, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the molecular orbital eigenvalues.

Initialize the MO_CUBES section.


### _class_ Mgrid(cutoff: float = 1200, rel_cutoff: float = 80, ngrids: int = 5, progression_factor: int = 3, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the multigrid for numerical integration.

Initialize the MGRID section.


* **Parameters**


    * **cutoff** – Cutoff energy (in Rydbergs for historical reasons) defining how find of
    Gaussians will be used


    * **rel_cutoff** – The relative cutoff energy, which defines how to map the Gaussians onto
    the multigrid. If the value is too low then, even if you have a high cutoff
    with sharp Gaussians, they will be mapped to the course part of the multigrid


    * **ngrids** – number of grids to use


    * **progression_factor** – divisor that decides how to map Gaussians the multigrid after
    the highest mapping is decided by rel_cutoff


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ OrbitalTransformation(minimizer: str = 'CG', preconditioner: str = 'FULL_ALL', algorithm: str = 'STRICT', rotation: bool = False, occupation_preconditioner: bool = False, energy_gap: float = -1, linesearch: str = '2PNT', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Turns on the Orbital Transformation scheme for diagonalizing the Hamiltonian. Often faster
and with guaranteed convergence compared to normal diagonalization, but requires the system
to have a band gap.

NOTE: OT has poor convergence for metallic systems and cannot use SCF mixing or smearing.
Therefore, you should not use it for metals or systems with ‘small’ band gaps. In that
case, use normal diagonalization

Initialize the OT section.


* **Parameters**


    * **minimizer** – The minimizer to use with the OT method. Default is conjugate gradient
    method, which is more robust, but more well-behaved systems should use DIIS, which
    can be as much as 50% faster.


    * **preconditioner** – Preconditioner to use for OT, FULL_ALL tends to be most robust,
    but is not always most efficient. For difficult systems, FULL_SINGLE_INVERSE can be
    more robust, and is reasonably efficient with large systems. For huge, but well
    behaved, systems, where construction of the preconditioner can take a very long
    time, FULL_KINETIC can be a good choice.


    * **algorithm** – What algorithm to use for OT. ‘Strict’: Taylor or diagonalization
    based algorithm. IRAC: Orbital Transformation based Iterative Refinement of the
    Approximate Congruence transformation (OT/IR).


    * **rotation** – Introduce additional variables to allow subspace rotations (i.e fractional
    occupations)


    * **occupation_preconditioner** – include the fractional occupation in the preconditioning


    * **energy_gap** – Guess for the band gap. For FULL_ALL, should be smaller than the
    actual band gap, so simply using 0.01 is a robust value. Choosing a larger value
    will help if you start with a bad initial guess though. For FULL_SINGLE_INVERSE,
    energy_gap is treated as a lower bound. Values lower than 0.05 in this case can
    lead to stability issues.


    * **linesearch** (*str*) – From the manual: 1D line search algorithm to be used with the OT
    minimizer, in increasing order of robustness and cost. MINIMIZER CG combined with
    LINESEARCH GOLD should always find an electronic minimum. Whereas the 2PNT
    minimizer is almost always OK, 3PNT might be needed for systems in which successive
    OT CG steps do not decrease the total energy.


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ PBE(parameterization: str = 'ORIG', scale_c: float = 1, scale_x: float = 1, keywords: dict | None = None, subsections: dict | None = None)
Bases: `Section`

Info about the PBE functional.


* **Parameters**


    * **parameterization** (*str*) – ORIG: original PBE
    PBESOL: PBE for solids/surfaces
    REVPBE: revised PBE


    * **scale_c** (*float*) – scales the correlation part of the functional.


    * **scale_x** (*float*) – scales the exchange part of the functional.


    * **keywords** – additional keywords


    * **subsections** – additional subsections.



### _class_ PDOS(nlumo: int = -1, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of projected density of states onto the different atom KINDS
(elemental decomposed DOS).

Initialize the PDOS section.


* **Parameters**


    * **nlumo** – how many unoccupied orbitals to include (-1==ALL)


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ PotentialFile(objects: Sequence | None = None)
Bases: `DataFile`

Data file for potentials only.


#### _classmethod_ from_str(string)
Initialize from a string representation.


### _class_ PotentialInfo(electrons: int | None = None, potential_type: str | None = None, nlcc: bool | None = None, xc: str | None = None)
Bases: `MSONable`

Metadata for this potential.


#### electrons()
Total number of electrons


* **Type**

    int | None



#### potential_type()
Potential type (e.g. GTH)


* **Type**

    str | None



#### nlcc()
Nonlinear core corrected potential


* **Type**

    bool | None



#### xc()
Exchange correlation functional used for creating this potential


* **Type**

    str | None



#### electrons(_: int | Non_ _ = Non_ )

#### _classmethod_ from_str(string)
Get a cp2k formatted string representation.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### nlcc(_: bool | Non_ _ = Non_ )

#### potential_type(_: str | Non_ _ = Non_ )

#### softmatch(other)
Soft matching to see if two potentials match.

Will only match those attributes which *are* defined for this basis info object (one way checking)


#### xc(_: str | Non_ _ = Non_ )

### _class_ QS(method: str = 'GPW', eps_default: float = 1e-10, eps_pgf_orb: float | None = None, extrapolation: str = 'ASPC', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the quickstep settings (DFT driver).

Initialize the QS Section.


* **Parameters**


    * **method** (*"GPW"** | **"GAPW"*) – What DFT methodology to use. GPW (Gaussian Plane Waves) for
    DFT with pseudopotentials or GAPW (Gaussian Augmented Plane Waves) for all
    electron calculations.


    * **eps_default** (*float*) – The default level of convergence accuracy. NOTE: This is a
    global value for all the numerical value of all EPS_\* values in QS module.
    It is not the same as EPS_SCF, which sets convergence accuracy of the SCF cycle
    alone.


    * **eps_pgf_orb** – Precision for the overlap matrix. Default is to use sqrt(eps_default)


    * **extrapolation** (*"PS"** | **"ASPC"*) – Method use for extrapolation. If using
    gamma-point-only calculation, then one should either PS
    or ASPC (ASPC especially for MD runs). See the manual for other options.


    * **keywords** – Additional keywords to add


    * **subsections** – Subsections to initialize with.



### _class_ Scf(max_scf: int = 50, eps_scf: float = 1e-06, scf_guess: str = 'RESTART', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the self consistent field loop.

Initialize the Scf section.


* **Parameters**


    * **max_scf** (*int*) – Maximum number of SCF loops before terminating. Defaults to 50.


    * **eps_scf** (*float*) – Convergence criteria for SCF loop. Defaults to 1e-6.


    * **scf_guess** – Initial guess for SCF loop.
    “ATOMIC”: Generate an atomic density using the atomic code
    “CORE”: Diagonalize the core Hamiltonian for an initial guess.
    “HISTORY_RESTART”: Extrapolated from previous RESTART files.
    “MOPAC”: Use same guess as MOPAC for semi-empirical methods or a simple

    > diagonal density matrix for other methods.

    ”NONE”: Skip initial guess (only for NON-SCC DFTB).
    “RANDOM”: Use random wavefunction coefficients.
    “RESTART”: Use the RESTART file as an initial guess (and ATOMIC if not present).
    “SPARSE”: Generate a sparse wavefunction using the atomic code (for OT based

    > methods).



    * **keywords** – Additional keywords


    * **subsections** – Additional subsections



### _class_ Section(name: str, subsections: dict | None = None, repeats: bool = False, description: str | None = None, keywords: dict | None = None, section_parameters: list | tuple | None = None, location: str | None = None, verbose: bool | None = True, alias: str | None = None, \*\*kwargs)
Bases: `MSONable`

Basic input representation of input to Cp2k. Activates functionality inside of the
Cp2k executable.

Basic object representing a CP2K Section. Sections activate different parts of the
calculation. For example, FORCE_EVAL section will activate CP2K’s ability to calculate
forces.


* **Parameters**


    * **name** – The name of the section (must match name in CP2K)


    * **subsections** – A dictionary of subsections that are nested in this section.
    Format is {‘NAME’: Section(*args, \*\*kwargs). The name you chose for ‘NAME’
    to index that subsection does not \*have* to be the same as the section’s true name,
    but we recommend matching them. You can specify a blank dictionary if there are
    no subsections, or if you want to insert the subsections later.


    * **repeats** – Whether or not this section can be repeated. Most sections cannot.
    Default=False.


    * **description** – Description of this section for easier readability


    * **keywords** – the keywords to be set for this section. Each element should be a
    Keyword object. This can be more cumbersome than simply using kwargs for building
    a class in a script, but is more convenient for the class instantiations of CP2K
    sections (see below).


    * **section_parameters** – the section parameters for this section. Section parameters
    are specialized keywords that modify the behavior of the section overall. Most
    sections do not have section parameters, but some do. Unlike normal Keywords,
    these are specified as strings and not as Keyword objects.


    * **location** – the path to the section in the form ‘SECTION/SUBSECTION1/SUBSECTION3’,
    example for QS module: ‘FORCE_EVAL/DFT/QS’. This location is used to automatically
    determine if a subsection requires a supersection to be activated.


    * **verbose** – Controls how much is printed to Cp2k input files (Also see Keyword).
    If True, then a description of the section will be printed with it as a comment
    (if description is set). Default=True.


    * **alias** – An alias for this class to use in place of the name.


    * **keyword** (*kwargs are interpreted as*) –


    * **as** (*value pairs and added to the keywords array*) –


    * **objects** (*Keyword*) –



#### _static_ _get_str(d, indent=0)
Helper function to return a pretty string of the section. Includes indentation and
descriptions (if present).


#### _static_ _update(d1, d2, strict=False)
Helper method for self.update(d) method (see above).


#### add(other)
Add another keyword to the current section.


#### by_path(path: str)
Access a sub-section using a path. Used by the file parser.


* **Parameters**

    **path** (*str*) – Path to section of form ‘SUBSECTION1/SUBSECTION2/SUBSECTION_OF_INTEREST’



#### check(path: str)
Check if section exists within the current using a path. Can be useful for cross-checking
whether or not required dependencies have been satisfied, which CP2K does not enforce.


* **Parameters**

    **path** (*str*) – Path to section of form ‘SUBSECTION1/SUBSECTION2/SUBSECTION_OF_INTEREST’



#### get(d, default=None)
Similar to get for dictionaries. This will attempt to retrieve the
section or keyword matching d. Will not raise an error if d does not exist.


* **Parameters**


    * **d** – the key to retrieve, if present


    * **default** – what to return if d is not found



#### get_keyword(d, default=None)
Get function, only for subsections.


* **Parameters**


    * **d** – Name of keyword to get


    * **default** – return if d is not found in keyword list



#### get_section(d, default=None)
Get function, only for subsections.


* **Parameters**


    * **d** – Name of section to get


    * **default** – return if d is not found in subsections



#### get_str()
Get string representation of Section.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### inc(d: dict)
Mongo style dict modification. Include.


#### insert(d)
Insert a new section as a subsection of the current one.


#### safeset(d: dict)
Alias for update with strict (no insertions). Used by custodian.


#### set(d: dict)
Alias for update. Used by custodian.


#### setitem(key, value, strict=False)
Helper function for setting items. Kept separate from the double-underscore function so that
“strict” option can be made possible.

strict will only set values for items that already have a key entry (no insertion).


#### silence()
Recursively delete all print sections so that only defaults are printed out.


#### unset(d: dict)
Dict based deletion. Used by custodian.


#### update(d: dict, strict=False)
Update the Section according to a dictionary argument. This is most useful
for providing user-override settings to default parameters. As you pass a
dictionary the class variables like “description”, “location”, or “repeats”
are not included. Therefore, it is recommended that this be used to modify
existing Section objects to a user’s needs, but not used for the creation
of new Section child-classes.


* **Parameters**


    * **d** (*dict*) – A dictionary containing the update information. Should use nested dictionaries
    to specify the full path of the update. If a section or keyword does not exist, it
    will be created, but only with the values that are provided in “d”, not using
    default values from a Section object.
    Example: {

    > ’SUBSECTION1’: {

    >     ‘SUBSEC2’: {‘NEW_KEYWORD’: ‘NEW_VAL’},
    >     ‘NEW_SUBSEC’: {‘NEW_KWD’: ‘NEW_VAL’}
    >     }

    > }



    * **strict** (*bool*) – If true, only update existing sections and keywords. If false, allow
    new sections and keywords. Default: False



#### verbosity(verbosity: bool)
Change the section verbosity recursively by turning on/off the printing of descriptions.
Turning off descriptions may reduce the appealing documentation of input files, but also
helps de-clutter them.


### _class_ SectionList(sections: Sequence[Section])
Bases: `MSONable`

Section list.

Initializes a SectionList object using a sequence of sections.


* **Parameters**

    **sections** – A list of keywords. Must all have the same name (case-insensitive)



#### _static_ _get_str(d, indent=0)

#### append(item)
Append the section list.


#### extend(lst: list)
Extend the section list.


#### get(d, index=-1)
Get for section list. If index is specified, return the section at that index.
Otherwise, return a get on the last section.


#### get_str()
Return string representation of section list.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### verbosity(verbosity)
Silence all sections in section list.


### _class_ Smear(elec_temp: float = 300, method: str = 'FERMI_DIRAC', fixed_magnetic_moment: float = -100.0, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Control electron smearing.

Basic object representing a CP2K Section. Sections activate different parts of the
calculation. For example, FORCE_EVAL section will activate CP2K’s ability to calculate
forces.


* **Parameters**


    * **name** – The name of the section (must match name in CP2K)


    * **subsections** – A dictionary of subsections that are nested in this section.
    Format is {‘NAME’: Section(*args, \*\*kwargs). The name you chose for ‘NAME’
    to index that subsection does not \*have* to be the same as the section’s true name,
    but we recommend matching them. You can specify a blank dictionary if there are
    no subsections, or if you want to insert the subsections later.


    * **repeats** – Whether or not this section can be repeated. Most sections cannot.
    Default=False.


    * **description** – Description of this section for easier readability


    * **keywords** – the keywords to be set for this section. Each element should be a
    Keyword object. This can be more cumbersome than simply using kwargs for building
    a class in a script, but is more convenient for the class instantiations of CP2K
    sections (see below).


    * **section_parameters** – the section parameters for this section. Section parameters
    are specialized keywords that modify the behavior of the section overall. Most
    sections do not have section parameters, but some do. Unlike normal Keywords,
    these are specified as strings and not as Keyword objects.


    * **location** – the path to the section in the form ‘SECTION/SUBSECTION1/SUBSECTION3’,
    example for QS module: ‘FORCE_EVAL/DFT/QS’. This location is used to automatically
    determine if a subsection requires a supersection to be activated.


    * **verbose** – Controls how much is printed to Cp2k input files (Also see Keyword).
    If True, then a description of the section will be printed with it as a comment
    (if description is set). Default=True.


    * **alias** – An alias for this class to use in place of the name.


    * **keyword** (*kwargs are interpreted as*) –


    * **as** (*value pairs and added to the keywords array*) –


    * **objects** (*Keyword*) –



### _class_ Subsys(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the definition of the system to be simulated.

Initialize the subsys section.


### _class_ V_Hartree_Cube(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the hartree potential as a cube file.

Basic object representing a CP2K Section. Sections activate different parts of the
calculation. For example, FORCE_EVAL section will activate CP2K’s ability to calculate
forces.


* **Parameters**


    * **name** – The name of the section (must match name in CP2K)


    * **subsections** – A dictionary of subsections that are nested in this section.
    Format is {‘NAME’: Section(*args, \*\*kwargs). The name you chose for ‘NAME’
    to index that subsection does not \*have* to be the same as the section’s true name,
    but we recommend matching them. You can specify a blank dictionary if there are
    no subsections, or if you want to insert the subsections later.


    * **repeats** – Whether or not this section can be repeated. Most sections cannot.
    Default=False.


    * **description** – Description of this section for easier readability


    * **keywords** – the keywords to be set for this section. Each element should be a
    Keyword object. This can be more cumbersome than simply using kwargs for building
    a class in a script, but is more convenient for the class instantiations of CP2K
    sections (see below).


    * **section_parameters** – the section parameters for this section. Section parameters
    are specialized keywords that modify the behavior of the section overall. Most
    sections do not have section parameters, but some do. Unlike normal Keywords,
    these are specified as strings and not as Keyword objects.


    * **location** – the path to the section in the form ‘SECTION/SUBSECTION1/SUBSECTION3’,
    example for QS module: ‘FORCE_EVAL/DFT/QS’. This location is used to automatically
    determine if a subsection requires a supersection to be activated.


    * **verbose** – Controls how much is printed to Cp2k input files (Also see Keyword).
    If True, then a description of the section will be printed with it as a comment
    (if description is set). Default=True.


    * **alias** – An alias for this class to use in place of the name.


    * **keyword** (*kwargs are interpreted as*) –


    * **as** (*value pairs and added to the keywords array*) –


    * **objects** (*Keyword*) –



### _class_ Xc_Functional(functionals: Iterable | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Defines the XC functional(s) to use.

Basic object representing a CP2K Section. Sections activate different parts of the
calculation. For example, FORCE_EVAL section will activate CP2K’s ability to calculate
forces.


* **Parameters**


    * **name** – The name of the section (must match name in CP2K)


    * **subsections** – A dictionary of subsections that are nested in this section.
    Format is {‘NAME’: Section(*args, \*\*kwargs). The name you chose for ‘NAME’
    to index that subsection does not \*have* to be the same as the section’s true name,
    but we recommend matching them. You can specify a blank dictionary if there are
    no subsections, or if you want to insert the subsections later.


    * **repeats** – Whether or not this section can be repeated. Most sections cannot.
    Default=False.


    * **description** – Description of this section for easier readability


    * **keywords** – the keywords to be set for this section. Each element should be a
    Keyword object. This can be more cumbersome than simply using kwargs for building
    a class in a script, but is more convenient for the class instantiations of CP2K
    sections (see below).


    * **section_parameters** – the section parameters for this section. Section parameters
    are specialized keywords that modify the behavior of the section overall. Most
    sections do not have section parameters, but some do. Unlike normal Keywords,
    these are specified as strings and not as Keyword objects.


    * **location** – the path to the section in the form ‘SECTION/SUBSECTION1/SUBSECTION3’,
    example for QS module: ‘FORCE_EVAL/DFT/QS’. This location is used to automatically
    determine if a subsection requires a supersection to be activated.


    * **verbose** – Controls how much is printed to Cp2k input files (Also see Keyword).
    If True, then a description of the section will be printed with it as a comment
    (if description is set). Default=True.


    * **alias** – An alias for this class to use in place of the name.


    * **keyword** (*kwargs are interpreted as*) –


    * **as** (*value pairs and added to the keywords array*) –


    * **objects** (*Keyword*) –


## pymatgen.io.cp2k.outputs module

This module defines the Cp2k output parser along with a few other functions for parsing cp2k-related
outputs.


### _class_ Cp2kOutput(filename, verbose=False, auto_load=False)
Bases: `object`

Class for parsing output file from CP2K. The CP2K output file is very flexible in the way that
it is returned. This class will automatically parse parameters that should always be present,
but other parsing features may be called depending on the run type.

Initialize the Cp2kOutput object.


* **Parameters**


    * **filename** – (str) Name of the CP2K output file to parse


    * **verbose** – (bool) Whether or not to parse with verbosity (will parse lots of data that
    may not be useful)


    * **auto_load** (*bool*) – Whether or not to automatically load basic info like energies
    and structures.



#### _static_ _gauss_smear(densities, energies, npts, width)

#### as_dict()
Return dictionary representation of the output.


#### _property_ band_structure(_: [BandStructure](pymatgen.electronic_structure.md#pymatgen.electronic_structure.bandstructure.BandStructure_ )
Returns band structure object if it has been parsed.


#### _property_ calculation_type()
Returns the calculation type (what io.vasp.outputs calls run_type).


#### _property_ charge(_: floa_ )
Get charge from the input file.


#### _property_ complete_dos(_: [CompleteDos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.CompleteDos_ )
Returns complete dos object if it has been parsed.


#### _property_ completed()
Did the calculation complete.


#### convergence()
Check whether or not the SCF and geometry optimization cycles converged.


#### _property_ cp2k_version()
The cp2k version used in the calculation.


#### _property_ is_hubbard(_: boo_ )
Returns True if hubbard +U correction was used.


#### _property_ is_metal(_: boo_ )
Was a band gap found? i.e. is it a metal.


#### _property_ is_molecule(_: boo_ )
True if the cp2k output was generated for a molecule (i.e.
no periodicity in the cell).


#### _property_ multiplicity(_: in_ )
Get the spin multiplicity from input file.


#### _property_ num_warnings()
How many warnings showed up during the run.


#### parse_atomic_kind_info()
Parse info on what atomic kinds are present and what basis/pseudopotential is describing
each of them.


#### parse_bandstructure(bandstructure_filename=None)
Parse a CP2K bandstructure file.


* **Parameters**


    * **bandstructure_filename** – Filename containing bandstructure info. If


    * **provided** (*not*) –


    * **by** (*then the pmg name** of **"BAND.bs" will be assumed*) –


    * **parser.** (*the filename*) –



#### parse_cell_params()
Parse the lattice parameters (initial) from the output file.


#### parse_chi_tensor(chi_filename=None)
Parse the magnetic susceptibility tensor.


#### parse_cp2k_params()
Parse the CP2K general parameters from CP2K output file into a dictionary.


#### parse_dft_params()
Parse the DFT parameters (as well as functional, HF, vdW params).


#### parse_dos(dos_file=None, pdos_files=None, ldos_files=None)
Parse the dos files produced by cp2k calculation. CP2K produces different files based
on the input file rather than assimilating them all into one file.

One file type is the overall DOS file, which is used for k-point calculations. For
non-kpoint calculation, the overall DOS is generally not calculated, but the
element-projected pDOS is. Separate files are created for each spin channel and each
atom kind. If requested, cp2k can also do site/local projected dos (ldos). Each site
requested will have a separate file for each spin channel (if spin polarized calculation
is performed).

If possible, this function will assimilate the ldos files into a CompleteDos object.
Either provide a list of PDOS file paths, or use glob to find the .pdos_ALPHA extension
in the calculation directory.


* **Parameters**


    * **dos_file** (*str*) – Name of the dos file, otherwise will be inferred


    * **pdos_files** (*list*) – list of pdos file paths, otherwise they will be inferred


    * **ldos_files** (*list*) – list of ldos file paths, otherwise they will be inferred



#### parse_energies()
Get the total energy from a CP2K calculation. Presently, the energy reported in the
trajectory (pos.xyz) file takes presidence over the energy reported in the main output
file. This is because the trajectory file keeps track of energies in between restarts,
while the main output file may or may not depending on whether a particular machine
overwrites or appends it.


#### parse_files()
Identify files present in the directory with the cp2k output file. Looks for trajectories,
dos, and cubes.


#### parse_forces()
Get the forces from the forces file, or from the main output file.


#### parse_global_params()
Parse the GLOBAL section parameters from CP2K output file into a dictionary.


#### parse_gtensor(gtensor_filename=None)
Parse a file containing g tensor.


#### parse_hirshfeld()
Parse the hirshfeld population analysis for each step.


#### parse_homo_lumo()
Find the HOMO - LUMO gap in [eV]. Returns the last value. For gaps/eigenvalues decomposed
by spin up/spin down channel and over many ionic steps, see parse_mo_eigenvalues().


#### parse_hyperfine(hyperfine_filename=None)
Parse a file containing hyperfine coupling tensors for each atomic site.


#### parse_initial_structure()
Parse the initial structure from the main cp2k output file.


#### parse_input()
Load in the input set from the input file (if it can be found).


#### parse_ionic_steps()
Parse the ionic step info. If already parsed, this will just assimilate.


#### parse_mo_eigenvalues()
Parse the MO eigenvalues from the cp2k output file. Will get the eigenvalues (and band gap)
at each ionic step (if more than one exist).

Everything is decomposed by spin channel. If calculation was performed without spin
polarization, then only Spin.up will be present, which represents the average of up and
down.


#### parse_mulliken()
Parse the mulliken population analysis info for each step


#### parse_nmr_shift()
Parse NMR calculation.


#### parse_opt_steps()
Parse the geometry optimization information.


#### parse_overlap_condition()
Retrieve the overlap condition number


#### parse_plus_u_params()
Parse the DFT+U params.


#### parse_qs_params()
Parse the DFT parameters (as well as functional, HF, vdW params).


#### parse_raman()
Parse raman calculation.


#### parse_scf_opt()
Parse the SCF cycles (not usually important).


#### parse_scf_params()
Retrieve the most import SCF parameters: the max number of scf cycles (max_scf),
the convergence cutoff for scf (eps_scf),.


#### parse_stresses()
Get the stresses from stress file, or from the main output file.


#### parse_structures(trajectory_file=None, lattice_file=None)
Parses the structures from a cp2k calculation. Static calculations simply use the initial
structure. For calculations with ionic motion, the function will look for the appropriate
trajectory and lattice files based on naming convention. If no file is given, and no file
is found, it is assumed that the lattice/structure remained constant, and the initial
lattice/structure is used. Cp2k does not output the trajectory in the main output file by
default, so non static calculations have to reference the trajectory file.


#### parse_tddfpt()
Parse TDDFPT calculation.


#### parse_timing()
Parse the timing info (how long did the run take).


#### parse_total_numbers()
Parse total numbers (not usually important).


#### _property_ project_name(_: st_ )
What project name was used for this calculation.


#### ran_successfully()
Sanity checks that the program ran successfully. Looks at the bottom of the CP2K output
file for the “PROGRAM ENDED” line, which is printed when successfully ran. Also grabs
the number of warnings issued.


#### read_pattern(patterns, reverse=False, terminate_on_match=False, postprocess=<class 'str'>)
This function originally comes from pymatgen.io.vasp.outputs Outcar class.

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


#### read_table_pattern(header_pattern, row_pattern, footer_pattern, postprocess=<class 'str'>, attribute_name=None, last_one_only=True, strip=None)
This function originally comes from pymatgen.io.vasp.outputs Outcar class.

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
    be returned. Default to be True.


    * **strip** (*list*) – Whether or not to strip contents out of the file before
    reading for a table pattern. This is mainly used by parse_scf_opt(),
    to strip HFX info out of the SCF loop start or DFT+U warnings out
    of the SCF loop iterations.



* **Returns**

    List of tables. 1) A table is a list of rows. 2) A row if either a list of
    attribute values in case the the capturing group is defined without name in
    row_pattern, or a dict in case that named capturing groups are defined by
    row_pattern.



#### _property_ run_type()
What type of run (Energy, MD, etc.) was performed.


#### _property_ spin_polarized(_: boo_ )
Was the calculation spin polarized.


### parse_dos(dos_file=None)
Parse a dos file. This format is different from the pdos files.


### parse_energy_file(energy_file)
Parses energy file for calculations with multiple ionic steps.


### parse_pdos(dos_file=None, spin_channel=None, total=False)
Parse a single DOS file created by cp2k. Must contain one PDOS snapshot. i.e. you cannot
use this cannot deal with multiple concatenated dos files.


* **Parameters**


    * **dos_file** (*list*) – list of pdos_ALPHA file paths


    * **spin_channel** (*int*) – Which spin channel the file corresponds to. By default, CP2K will
    write the file with ALPHA or BETA in the filename (for spin up or down), but
    you can specify this here, in case you have a manual file name.
    spin_channel == 1 –> spin up, spin_channel == -1 –> spin down.


    * **total** (*bool*) – Whether to grab the total occupations, or the orbital decomposed ones.


    * **sigma** (*float*) – width for gaussian smearing, if desired



* **Returns**

    >
    > 1. orbital decomposed DOS dict:
    > i.e. pdoss = {specie: {orbital.s: {Spin.up: … }, orbital.px: {Spin.up: … } …}}


    > 2. energy levels of this dos file


    > 3. fermi energy (in eV).

    DOS object is not created here




* **Return type**

    Everything necessary to create a dos object, in dict format


## pymatgen.io.cp2k.sets module

This module defines input sets for CP2K and is a work in progress. The structure/philosophy
of this module is based on the Vasp input sets in Pymatgen. These sets are meant to contain
tested parameters that will result in successful, reproducible, consistent calculations without
need for intervention 99% of the time. 99% of the time, you only need to provide a pymatgen
structure object and let the defaults take over from there.

The sets are intended to be very general, e.g. a set for geometry relaxation, and so most of the
time, if you have specific needs, you can simply specify them via the keyword argument
override_default_params (see Section.update() method). If you have the need to create a new input
set (say for a standardized high throughput calculation) then you can create a new child of the
Cp2kInputSet class.

In order to implement a new Set within the current code structure, follow this 3 step flow:


    1. Inherit from Cp2kInputSet or one of its children and call the super() constructor


    2. Create the new sections and insert them into self and its subsections as needed


    3. Call self.update(override_default_params) in order to allow user settings.


### _class_ CellOptSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for cell optimization relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _exception_ Cp2kValidationError(message)
Bases: `Exception`

Cp2k Validation Exception. Not exhausted. May raise validation
errors for features which actually do work if using a newer version
of cp2k.


#### CP2K_VERSION(_ = 'v2022.1_ )

### _class_ DftSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), project_name: str = 'CP2K', basis_and_potential: dict | None = None, xc_functionals: list | str | None = None, multiplicity: int = 0, ot: bool = True, energy_gap: float = -1, qs_method: str = 'GPW', eps_default: float = 1e-12, eps_scf: float = 1e-06, max_scf: int | None = None, minimizer: str = 'DIIS', preconditioner: str = 'FULL_SINGLE_INVERSE', algorithm: str = 'IRAC', linesearch: str = '2PNT', rotation: bool = True, occupation_preconditioner: bool = False, cutoff: float | None = None, rel_cutoff: int = 50, ngrids: int = 5, progression_factor: int = 3, override_default_params: dict | None = None, wfn_restart_file_name: str | None = None, kpoints: VaspKpoints | None = None, smearing: bool = False, \*\*kwargs)
Bases: `Cp2kInput`

Base for an input set using the Quickstep module (i.e. a DFT calculation). The DFT section is
pretty vast in CP2K, so this set hopes to make the DFT setup fairly simple. The provided
parameters are pretty conservative, and so they should not need to be changed very often.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



#### activate_epr(\*\*kwargs)
Calculate g-tensor. Requires localize. Suggested with GAPW.


#### activate_fast_minimization(on)
Method to modify the set to use fast SCF minimization.


#### activate_hybrid(hybrid_functional: str = 'PBE0', hf_fraction: float = 0.25, gga_x_fraction: float = 0.75, gga_c_fraction: float = 1, max_memory: int = 2000, cutoff_radius: float = 8.0, potential_type: str | None = None, omega: float = 0.11, scale_coulomb: float = 1, scale_gaussian: float = 1, scale_longrange: float = 1, admm: bool = True, admm_method: str = 'BASIS_PROJECTION', admm_purification_method: str = 'NONE', admm_exch_correction_func: str = 'DEFAULT', eps_schwarz: float = 1e-07, eps_schwarz_forces: float = 1e-06, screen_on_initial_p: bool = True, screen_p_forces: bool = True)
Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.

Note 1: When running ADMM with cp2k, memory is very important. If the memory requirements
exceed what is available (see max_memory), then CP2K will have to calculate the 4-electron
integrals for HFX during each step of the SCF cycle. ADMM provides a huge speed up by
making the memory requirements *feasible* to fit into RAM, which means you only need to
calculate the integrals once each SCF cycle. But, this only works if it fits into memory.
When setting up ADMM calculations, we recommend doing whatever is possible to fit all the
4EI into memory.

Note 2: This set is designed for reliable high-throughput calculations, NOT for extreme
accuracy. Please review the in-line comments in this method if you want more control.


* **Parameters**


    * **hybrid_functional** (*str*) – Type of hybrid functional. This set supports HSE (screened)
    and PBE0 (truncated). Default is PBE0, which converges easier in the GPW basis
    used by cp2k.


    * **hf_fraction** (*float*) – fraction of exact HF exchange energy to mix. Default: 0.25


    * **gga_x_fraction** (*float*) – fraction of gga exchange energy to retain. Default: 0.75


    * **gga_c_fraction** (*float*) – fraction of gga correlation energy to retain. Default: 1.0


    * **max_memory** (*int*) – Maximum memory available to each MPI process (in Mb) in the
    calculation. Most modern computing nodes will have ~2Gb per core, or 2048 Mb,
    but check for your specific system. This value should be as large as possible
    while still leaving some memory for the other parts of cp2k. Important: If
    this value is set larger than the memory limits, CP2K will likely seg-fault.
    Default: 2000


    * **cutoff_radius** (*float*) – for truncated hybrid functional (i.e. PBE0), this is the cutoff
    radius. The default is selected as that which generally gives convergence, but
    maybe too low (if you want very high accuracy) or too high (if you want a quick
    screening). Default: 8 angstroms


    * **potential_type** (*str*) – what interaction potential to use for HFX. Available in CP2K are
    COULOMB, GAUSSIAN, IDENTITY, LOGRANGE, MIX_CL, MIX_CL_TRUNC, MIX_LG, SHORTRANGE,
    and TRUNCATED. Default is None, and it will be set automatically depending on the
    named hybrid_functional that you use, but setting it to one of the acceptable
    values will constitute a user-override.


    * **omega** (*float*) – For HSE, this specifies the screening parameter. HSE06 sets this as
    0.2, which is the default.


    * **scale_coulomb** – Scale for the coulomb operator if using a range separated functional


    * **scale_gaussian** – Scale for the gaussian operator (if applicable)


    * **scale_longrange** – Scale for the coulomb operator if using a range separated functional


    * **admm** – Whether or not to use the auxiliary density matrix method for the exact
    HF exchange contribution. Highly recommended. Speed ups between 10x and 1000x are
    possible when compared to non ADMM hybrid calculations.


    * **admm_method** – Method for constructing the auxiliary basis


    * **admm_purification_method** – Method for purifying the auxiliary density matrix so as to
    preserve properties, such as idempotency. May lead to shifts in the
    eigenvalues.


    * **admm_exch_correction_func** – Which functional to use to calculate the exchange correction
    E_x(primary) - E_x(aux)


    * **eps_schwarz** – Screening threshold for HFX, in Ha. Contributions smaller than
    this will be screened. The smaller the value, the more accurate, but also the more
    costly. Default value is 1e-7. 1e-6 works in a large number of cases, but is
    quite aggressive, which can lead to convergence issues.


    * **eps_schwarz_forces** – Same as for eps_schwarz, but for screening contributions to
    forces. Convergence is not as sensitive with respect to eps_schwarz forces as
    compared to eps_schwarz, and so 1e-6 should be good default.


    * **screen_on_initial_p** – If an initial density matrix is provided, in the form of a
    CP2K wfn restart file, then this initial density will be used for screening. This
    is generally very computationally efficient, but, as with eps_schwarz, can lead to
    instabilities if the initial density matrix is poor.


    * **screen_p_forces** – Same as screen_on_initial_p, but for screening of forces.



#### activate_hyperfine()
Print the hyperfine coupling constants.


#### activate_localize(states='OCCUPIED', preconditioner='FULL_ALL', restart=False)
Activate calculation of the maximally localized wannier functions.


* **Parameters**


    * **states** – Which states to calculate. occupied, unoccupied, mixed states, or all states. At
    present, unoccupied orbitals are only implemented for GPW.


    * **preconditioner** – Preconditioner to use for optimize


    * **restart** – Initialize from the localization restart file



#### activate_motion(max_drift: float = 0.003, rms_drift: float = 0.0015, max_force: float = 0.00045, rms_force: float = 0.0003, max_iter: int = 200, optimizer: str = 'BFGS', trust_radius: float = 0.25, line_search: str = '2PNT', ensemble: str = 'NVE', temperature: float = 300, timestep: float = 0.5, nsteps: int = 3, thermostat: str = 'NOSE', nproc_rep: int = 1)
Turns on the motion section for GEO_OPT, CELL_OPT, etc. calculations.
Will turn on the printing subsections and also bind any constraints
to their respective atoms.


#### activate_nmr(\*\*kwargs)
Calculate nmr shifts. Requires localize. Suggested with GAPW.


#### activate_nonperiodic(solver='ANALYTIC')
Activates a calculation with non-periodic calculations by turning of PBC and
changing the poisson solver. Still requires a CELL to put the atoms.


#### activate_polar(\*\*kwargs)
Calculate polarizations (including raman).


#### activate_robust_minimization()
Method to modify the set to use more robust SCF minimization technique.


#### activate_spinspin(\*\*kwargs)
Calculate spin-spin coupling tensor. Requires localize.


#### activate_tddfpt(\*\*kwargs)
Activate TDDFPT for calculating excited states. Only works with GPW. Supports hfx.


#### activate_vdw_potential(dispersion_functional: str, potential_type: str)
Activate van der Waals dispersion corrections.


* **Parameters**


    * **dispersion_functional** – Type of dispersion functional.
    Options: pair_potential or non_local


    * **potential_type** – What type of potential to use, given a dispersion functional type
    Options: DFTD2, DFTD3, DFTD3(BJ), DRSLL, LMKLL, RVV10



#### activate_very_strict_minimization()
Method to modify the set to use very strict SCF minimization scheme


#### create_subsys(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Create the structure for the input.


#### _static_ get_basis_and_potential(structure, basis_and_potential)
Get a dictionary of basis and potential info for constructing the input file.

data in basis_and_potential argument can be specified in several ways:

> Strategy 1: Element-specific info (takes precedence)

> >
> > 1. Provide a basis and potential object:

> > > el: {‘basis’: obj, ‘potential’: obj}


> > 2. Provide a hash of the object that matches the keys in the pmg configured cp2k data files.

> > > el: {‘basis’: hash, ‘potential’: hash}

> > 3. Provide the name of the basis and potential AND the basis_filenames and potential_filename
> > keywords specifying where to find these objects

> > > el: {

> > >     ‘basis’: name, ‘potential’: name, ‘basis_filenames’: [filenames],
> > >     ‘potential_filename’: filename

> > > }

> Strategy 2: global descriptors

> > In this case, any elements not present in the argument will be dealt with by searching the pmg
> > configured cp2k data files to find a objects matching your requirements.

> > >
> > > * functional: Find potential and basis that have been optimized for a specific functional like PBE.

> > >     Can be None if you do not require them to match.


> > > * basis_type: type of basis to search for (e.g. DZVP-MOLOPT).


> > > * aux_basis_type: type of basis to search for (e.g. pFIT). Some elements do not have all aux types

> > >     available. Use aux_basis_type that is universal to avoid issues, or avoid using this strategy.


> > > * potential_type: “Pseudopotential” or “All Electron”

> > **\*BE WARNED\*** CP2K data objects can have the same name, this will sort those and choose the first one
> > that matches.

Will raise an error if no basis/potential info can be found according to the input.


#### _static_ get_cutoff_from_basis(basis_sets, rel_cutoff)
Given a basis and a relative cutoff. Determine the ideal cutoff variable.


#### _static_ get_xc_functionals(xc_functionals: list | str | None = None)
Get XC functionals. If simplified names are provided in kwargs, they
will be expanded into their corresponding X and C names.


#### modify_dft_print_iters(iters, add_last='no')
Modify all DFT print iterations at once. Common use is to set iters to the max
number of iterations + 1 and then set add_last to numeric. This would have the
effect of printing only the first and last iteration, which might be useful for
speeding up/saving space on GEO_OPT or MD runs where you don’t need the intermediate
values.


* **Parameters**


    * **iters** (*int*) – print each “iters” iterations.


    * **add_last** (*str*) – Whether to explicitly include the last iteration, and how to mark it.
    numeric: mark last iteration with the iteration number
    symbolic: mark last iteration with the letter “l”
    no: do not explicitly include the last iteration



#### print_bandstructure(kpoints_line_density: int = 20)
Attaches a non-scf band structure calc the end of an SCF loop.

This requires a kpoint calculation, which is not always default in cp2k.


* **Parameters**

    **kpoints_line_density** – number of kpoints along each branch in line-mode calc.



#### print_dos(ndigits=6)
Activate printing of the overall DOS file.

Note: As of 2022.1, ndigits needs to be set to a sufficient value to ensure data is not lost.
Note: As of 2022.1, can only be used with a k-point calculation.


#### print_e_density(stride=(2, 2, 2))
Controls the printing of cube files with electronic density and, for UKS, the spin density.


#### print_forces()
Print out the forces and stress during calculation.


#### print_hirshfeld(on=True)
Activate or deactivate printing of Hirshfeld charges.


#### print_ldos(nlumo: int = -1)
Activate the printing of LDOS files, printing one for each atom kind by default.


* **Parameters**

    **nlumo** (*int*) – Number of virtual orbitals to be added to the MO set (-1=all).
    CAUTION: Setting this value to be higher than the number of states present may
    cause a Cholesky error.



#### print_mo()
Print molecular orbitals when running non-OT diagonalization.


#### print_mo_cubes(write_cube: bool = False, nlumo: int = -1, nhomo: int = -1)
Activate printing of molecular orbitals.


* **Parameters**


    * **write_cube** (*bool*) – whether to write cube file for the MOs instead of out file


    * **nlumo** (*int*) – Controls the number of lumos printed and dumped as a cube (-1=all)


    * **nhomo** (*int*) – Controls the number of homos printed and dumped as a cube (-1=all)



#### print_mulliken(on=False)
Activate or deactivate printing of Mulliken charges.


#### print_pdos(nlumo: int = -1)
Activate creation of the PDOS file.


* **Parameters**

    **nlumo** (*int*) – Number of virtual orbitals to be added to the MO set (-1=all).
    CAUTION: Setting this value to be higher than the number of states present may
    cause a Cholesky error.



#### print_v_hartree(stride=(2, 2, 2))
Controls the printing of a cube file with eletrostatic potential generated by the

    total density (electrons+ions). It is valid only for QS with GPW formalism.

Note that by convention the potential has opposite sign than the expected physical one.


#### set_charge(charge: int)
Set the overall charge of the simulation cell.


#### validate()
Implements a few checks for a valid input set.


#### write_basis_set_file(basis_sets, fn='BASIS')
Write the basis sets to a file.


#### write_potential_file(potentials, fn='POTENTIAL')
Write the potentials to a file.


### _class_ HybridCellOptSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for hybrid cell optimization relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ HybridRelaxSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for hybrid geometry relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ HybridStaticSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for static calculations.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ RelaxSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for geometry relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ StaticSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for static calculations.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** (*Kpoints*) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.


## pymatgen.io.cp2k.utils module

Utility functions for assisting with cp2k IO.


### chunk(string: str)
Chunk the string from a cp2k basis or potential file.


### get_truncated_coulomb_cutoff(inp_struct: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Get the truncated Coulomb cutoff for a given structure.


### get_unique_site_indices(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Get unique site indices for a structure according to site properties. Whatever site-property
has the most unique values is used for indexing.

For example, if you have magnetic CoO with half Co atoms having a positive moment, and the
other half having a negative moment. Then this function will create a dict of sites for
Co_1, Co_2, O. This function also deals with “Species” properties like oxi_state and spin by
pushing them to site properties.

This creates unique sites, based on site properties, but does not have anything to do with
turning those site properties into CP2K input parameters. This will only be done for properties
which can be turned into CP2K input parameters, which are stored in parsable_site_properties.


### natural_keys(text: str)
Sort text by numbers coming after an underscore with natural number
convention,
Ex: [file_1, file_12, file_2] becomes [file_1, file_2, file_12].


### postprocessor(data: str)
Helper function to post process the results of the pattern matching functions in Cp2kOutput
and turn them to Python types.


* **Parameters**

    **data** (*str*) – The data to be post processed.



* **Raises**

    **ValueError** – If the data cannot be parsed.



* **Returns**

    The post processed data.



* **Return type**

    str | float | bool | None



### preprocessor(data: str, dir: str = '.')
Cp2k contains internal preprocessor flags that are evaluated before execution. This helper
function recognizes those preprocessor flags and replaces them with an equivalent cp2k input
(this way everything is contained neatly in the cp2k input structure, even if the user preferred
to use the flags.

CP2K preprocessor flags (with arguments) are:

> @INCLUDE FILENAME: Insert the contents of FILENAME into the file at

>     this location.

> @SET VAR VALUE: set a variable, VAR, to have the value, VALUE.
> $VAR or ${VAR}: replace these with the value of the variable, as set

> > by the @SET flag.

> @IF/@ELIF: Not implemented yet.


* **Parameters**


    * **data** (*str*) – cp2k input to preprocess


    * **dir** (*str**, **optional*) – Path for include files. Default is ‘.’ (current directory).



* **Returns**

    Preprocessed string