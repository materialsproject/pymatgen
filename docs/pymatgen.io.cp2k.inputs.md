---
layout: default
title: pymatgen.io.cp2k.inputs.md
nav_exclude: true
---

# pymatgen.io.cp2k.inputs module

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


### _class_ pymatgen.io.cp2k.inputs.AtomicMetadata(info: BasisInfo | PotentialInfo | None = None, element: Element | None = None, potential: Literal['All Electron', 'Pseudopotential'] | None = None, name: str | None = None, alias_names: list = <factory>, filename: str | None = None, version: str | None = None)
Bases: `MSONable`

Metadata for basis sets and potentials in cp2k.


#### info()
Info about this object


* **Type**

    BasisInfo | PotentialInfo | None



#### element()
Element for this object


* **Type**

    [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | None



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

#### element(_: [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | Non_ _ = Non_ )

#### filename(_: str | Non_ _ = Non_ )

#### get_hash()
Get a hash of this object.


#### get_string()
Get string representation.


#### info(_: BasisInfo | PotentialInfo | Non_ _ = Non_ )

#### name(_: str | Non_ _ = Non_ )

#### potential(_: Literal['All Electron', 'Pseudopotential'] | Non_ _ = Non_ )

#### softmatch(other)
Soft matching to see if a desired basis/potential matches requirements.

Does soft matching on the “info” attribute first. Then soft matches against the
element and name/aliases.


#### version(_: str | Non_ _ = Non_ )

### _class_ pymatgen.io.cp2k.inputs.Band_Structure(kpoint_sets: Sequence[Kpoint_Set], filename: str = 'BAND.bs', added_mos: int = -1, keywords: dict | None = None, subsections: dict | None = None)
Bases: `Section`

Specifies high symmetry paths for outputting the band structure in CP2K.


* **Parameters**


    * **kpoint_sets** – Sequence of Kpoint_Set objects for the band structure calculation.


    * **filename** – Filename for the band structure output


    * **added_mos** – Added (unoccupied) molecular orbitals for the calculation.


    * **keywords** – additional keywords


    * **subsections** – additional subsections.



#### _static_ from_kpoints(kpoints: [Kpoints](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints), kpoints_line_density=20)
Initialize band structure section from a line-mode Kpoint object.


* **Parameters**


    * **kpoints** – a kpoint object from the vasp module, which was constructed in line mode


    * **kpoints_line_density** – Number of kpoints along each path



### _class_ pymatgen.io.cp2k.inputs.BasisFile(objects: Sequence | None = None)
Bases: `DataFile`

Data file for basis sets only.


#### _classmethod_ from_str(string)
Initialize from a string representation.


### _class_ pymatgen.io.cp2k.inputs.BasisInfo(electrons: int | None = None, core: int | None = None, valence: int | None = None, polarization: int | None = None, diffuse: int | None = None, cc: bool | None = False, pc: bool | None = False, sr: bool | None = False, molopt: bool | None = False, admm: bool | None = False, lri: bool | None = False, contracted: bool | None = None, xc: str | None = None)
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

### _class_ pymatgen.io.cp2k.inputs.BrokenSymmetry(l_alpha: Sequence = (-1,), n_alpha: Sequence = (0,), nel_alpha: Sequence = (-1,), l_beta: Sequence = (-1,), n_beta: Sequence = (0,), nel_beta: Sequence = (-1,))
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


### _class_ pymatgen.io.cp2k.inputs.Cell(lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), keywords: dict | None = None, \*\*kwargs)
Bases: `Section`

Defines the simulation cell (lattice).

Initialize the cell section.


* **Parameters**


    * **lattice** – pymatgen lattice object


    * **keywords** – additional keywords



### _class_ pymatgen.io.cp2k.inputs.Coord(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), aliases: dict | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Specifies the coordinates of the atoms using a pymatgen structure object.


* **Parameters**


    * **structure** – Pymatgen structure object


    * **alias** (*bool*) – whether or not to identify the sites by Element + number so you can do
    things like assign unique magnetization do different elements.


    * **keywords** – additional keywords


    * **subsections** – additional subsections.



### _class_ pymatgen.io.cp2k.inputs.Cp2kInput(name: str = 'CP2K_INPUT', subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Special instance of ‘Section’ class that is meant to represent the overall cp2k input.
Distinguishes itself from Section by overriding get_string() to not print this section’s
title and by implementing the file i/o.

Initialize Cp2kInput by calling the super.


#### _static_ from_file(file: str)
Initialize from a file.


#### _classmethod_ from_lines(lines: list | tuple)
Helper method to read lines of file.


#### _static_ from_str(s: str)
Initialize from a string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_string()
Get string representation of the Cp2kInput.


#### write_file(input_filename: str = 'cp2k.inp', output_dir: str = '.', make_dir_if_not_present: bool = True)
Write input to a file.


* **Parameters**


    * **input_filename** (*str**, **optional*) – Defaults to “cp2k.inp”.


    * **output_dir** (*str**, **optional*) – Defaults to “.”.


    * **make_dir_if_not_present** (*bool**, **optional*) – Defaults to True.



### _class_ pymatgen.io.cp2k.inputs.DOS(ndigits: int = 6, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the density of states.

Initialize the DOS section.


* **Parameters**


    * **ndigits** – how many digits of precision to print. As of 2022.1,
    this is necessary to not lose information.


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ pymatgen.io.cp2k.inputs.DataFile(objects: Sequence | None = None)
Bases: `MSONable`

A data file for a cp2k calc.


#### _classmethod_ from_file(fn)
Load from a file.


#### _classmethod_ from_str()
Initialize from a string.


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_string()
Get string representation.


#### objects(_: Sequence | Non_ _ = Non_ )

#### write_file(fn)
Write to a file.


### _class_ pymatgen.io.cp2k.inputs.Davidson(new_prec_each: int = 20, preconditioner: str = 'FULL_SINGLE_INVERSE', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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

    ”FULL_SINGLE”: Based on H-eS diagonalisation, not as good as FULL_ALL, but

        somewhat cheaper to apply.

    ”FULL_SINGLE_INVERSE”: Based on H-eS cholesky inversion, similar to FULL_SINGLE

        in preconditioning efficiency but cheaper to construct, might be somewhat
        less robust. Recommended for large systems.

    ”FULL_S_INVERSE”: Cholesky inversion of S, not as good as FULL_KINETIC,

        yet equally expensive.

    ”NONE”: skip preconditioning



    * **keywords** – additional keywords


    * **subsections** – additional subsections.



### _class_ pymatgen.io.cp2k.inputs.Dft(basis_set_filenames: Iterable = ('BASIS_MOLOPT',), potential_filename='GTH_POTENTIALS', uks: bool = True, wfn_restart_file_name: str | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.DftPlusU(eps_u_ramping=1e-05, init_u_ramping_each_scf=False, l=-1, u_minus_j=0, u_ramping=0, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.Diagonalization(eps_adapt: float = 0, eps_iter: float = 1e-08, eps_jacobi: float = 0, jacobi_threshold: float = 1e-07, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls diagonalization settings (if using traditional diagonalization).

Initialize the diagronalization section.


### _class_ pymatgen.io.cp2k.inputs.E_Density_Cube(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.ForceEval(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the calculation of energy and forces in Cp2k.

Initialize the ForceEval section.


### _class_ pymatgen.io.cp2k.inputs.GaussianTypeOrbitalBasisSet(info: BasisInfo | None = None, element: Element | None = None, potential: Literal['All Electron', 'Pseudopotential'] | None = None, name: str | None = None, alias_names: list = <factory>, filename: str | None = None, version: str | None = None, nset: int | None = None, n: list[int] | None = None, lmax: list[int] | None = None, lmin: list[int] | None = None, nshell: list[dict[int, int]] | None = None, exponents: list[list[float]] | None = None, coefficients: list[dict[int, dict[int, dict[int, float]]]] | None = None)
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


#### get_string()
Get standard cp2k GTO formatted string.


#### info(_: BasisInfo | Non_ _ = Non_ )

#### lmax(_: list[int] | Non_ _ = Non_ )

#### lmin(_: list[int] | Non_ _ = Non_ )

#### n(_: list[int] | Non_ _ = Non_ )

#### _property_ nexp()
Number of exponents.


#### nset(_: int | Non_ _ = Non_ )

#### nshell(_: list[dict[int, int]] | Non_ _ = Non_ )

### _class_ pymatgen.io.cp2k.inputs.Global(project_name: str = 'CP2K', run_type: str = 'ENERGY_FORCE', keywords: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls ‘global’ settings for cp2k execution such as RUN_TYPE and PROJECT_NAME.

Initialize the global section.


* **Parameters**


    * **project_name** – Defaults to “CP2K”.


    * **run_type** – what type of calculation to run


    * **keywords** – Additional keywords to add



### _class_ pymatgen.io.cp2k.inputs.GthPotential(info: PotentialInfo = None, element: Element | None = None, potential: Literal['All Electron', 'Pseudopotential'] | None = None, name: str | None = None, alias_names: list = <factory>, filename: str | None = None, version: str | None = None, n_elecs: dict[int, int] | None = None, r_loc: float | None = None, nexp_ppl: int | None = None, c_exp_ppl: Sequence | None = None, radii: dict[int, float] | None = None, nprj: int | None = None, nprj_ppnl: dict[int, int] | None = None, hprj_ppnl: dict[int, dict[int, dict[int, float]]] | None = None)
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


#### get_string()
Convert model to a GTH-formatted string.


#### hprj_ppnl(_: dict[int, dict[int, dict[int, float]]] | Non_ _ = Non_ )

#### n_elecs(_: dict[int, int] | Non_ _ = Non_ )

#### nexp_ppl(_: int | Non_ _ = Non_ )

#### nprj(_: int | Non_ _ = Non_ )

#### nprj_ppnl(_: dict[int, int] | Non_ _ = Non_ )

#### r_loc(_: float | Non_ _ = Non_ )

#### radii(_: dict[int, float] | Non_ _ = Non_ )

### _class_ pymatgen.io.cp2k.inputs.Keyword(name: str, \*values, description: str | None = None, units: str | None = None, verbose: bool | None = True, repeats: bool | None = False)
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


#### get_string()
String representation of Keyword.


#### verbosity(v)
Change the printing of this keyword’s description.


### _class_ pymatgen.io.cp2k.inputs.KeywordList(keywords: Sequence[Keyword])
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


#### get_string(indent=0)
String representation of Keyword.


#### verbosity(verbosity)
Silence all keywords in keyword list.


### _class_ pymatgen.io.cp2k.inputs.Kind(specie: str, alias: str | None = None, magnetization: float = 0.0, basis_set: GaussianTypeOrbitalBasisSet | str | None = 'GTH_BASIS', potential: GthPotential | str | None = 'GTH_POTENTIALS', ghost: bool = False, aux_basis: GaussianTypeOrbitalBasisSet | str | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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


    * **ghost** – Turn this into ghost atom (disaple the potential)


    * **aux_basis** – Auxiliary basis to use with ADMM


    * **keywords** – additional keywords


    * **subsections** – additional subsections


    * **kwargs** – Additional kwargs to pass to Section()



### _class_ pymatgen.io.cp2k.inputs.Kpoint_Set(npoints: int, kpoints: Iterable, units: str = 'B_VECTOR')
Bases: `Section`

Specifies a kpoint line to be calculated between special points.


* **Parameters**


    * **npoints** (*int*) – Number of kpoints along the line.


    * **kpoints** – A dictionary of {label: kpoint} kpoints defining the path


    * **units** (*str*) – Units for the kpoint coordinates.
    Options: “B_VECTOR” (reciprocal coordinates)

    > ”CART_ANGSTROM” (units of 2\*Pi/Angstrom)
    > “CART_BOHR” (units of 2\*Pi/Bohr).




### _class_ pymatgen.io.cp2k.inputs.Kpoints(kpts: Sequence | Sequence[Sequence[int]], weights: Sequence | None = None, eps_geo: float = 1e-06, full_grid: bool = False, parallel_group_size: int = -1, scheme: str = 'MONKHORST-PACK', symmetry: bool = False, units: str = 'B_VECTOR', verbose: bool = False, wavefunctions: str = 'COMPLEX')
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



#### _classmethod_ from_kpoints(kpoints: [Kpoints](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints), structure=None)
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



### _class_ pymatgen.io.cp2k.inputs.LDOS(index: int = 1, alias: str | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the LDOS (List-Density of states). i.e. projects onto specific atoms.

Initialize the LDOS section.


* **Parameters**


    * **index** – Index of the atom to project onto


    * **alias** – section alias


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ pymatgen.io.cp2k.inputs.MO_Cubes(write_cube: bool = False, nhomo: int = 1, nlumo: int = 1, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of the molecular orbital eigenvalues.

Initialize the MO_CUBES section.


### _class_ pymatgen.io.cp2k.inputs.Mgrid(cutoff: int | float = 1200, rel_cutoff: int | float = 80, ngrids: int = 5, progression_factor: int = 3, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the multigrid for numerical integration.

Initialize the MGRID section.


* **Parameters**


    * **cutoff** – Cutoff energy (in Rydbergs for historical reasons) defining how find of
    Gaussians will be used


    * **rel_cutoff** – The relative cutoff energy, which defines how to map the Gaussians onto
    the multigrid. If the the value is too low then, even if you have a high cutoff
    with sharp Gaussians, they will be mapped to the course part of the multigrid


    * **ngrids** – number of grids to use


    * **progression_factor** – divisor that decides how to map Gaussians the multigrid after
    the highest mapping is decided by rel_cutoff


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ pymatgen.io.cp2k.inputs.OrbitalTransformation(minimizer: str = 'CG', preconditioner: str = 'FULL_ALL', algorithm: str = 'STRICT', rotation: bool = False, occupation_preconditioner: bool = False, energy_gap: float = -1, linesearch: str = '2PNT', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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
    Approximative Congruence transformation (OT/IR).


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



### _class_ pymatgen.io.cp2k.inputs.PBE(parameterization: str = 'ORIG', scale_c: float | int = 1, scale_x: float | int = 1, keywords: dict | None = None, subsections: dict | None = None)
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



### _class_ pymatgen.io.cp2k.inputs.PDOS(nlumo: int = -1, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls printing of projected density of states onto the different atom KINDS
(elemental decomposed DOS).

Initialize the PDOS section.


* **Parameters**


    * **nlumo** – how many unoccupied orbitals to include (-1==ALL)


    * **keywords** – additional keywords


    * **subsections** – additional subsections



### _class_ pymatgen.io.cp2k.inputs.PotentialFile(objects: Sequence | None = None)
Bases: `DataFile`

Data file for potentials only.


#### _classmethod_ from_str(string)
Initialize from a string representation.


### _class_ pymatgen.io.cp2k.inputs.PotentialInfo(electrons: int | None = None, potential_type: str | None = None, nlcc: bool | None = None, xc: str | None = None)
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

### _class_ pymatgen.io.cp2k.inputs.QS(method: str = 'GPW', eps_default: float = 1e-10, eps_pgf_orb: float | None = None, extrapolation: str = 'ASPC', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.Scf(max_scf: int = 50, eps_scf: float = 1e-06, scf_guess: str = 'RESTART', keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.Section(name: str, subsections: dict | None = None, repeats: bool = False, description: str | None = None, keywords: dict | None = None, section_parameters: list | tuple | None = None, location: str | None = None, verbose: bool | None = True, alias: str | None = None, \*\*kwargs)
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
section or keyword matching d. Will not raise an error if d does not
exist.


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



#### get_string()
Get string representation of Section.


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



#### verbosity(verbosity)
Change the section verbosity recursively by turning on/off the printing of descriptions.
Turning off descriptions may reduce the appealing documentation of input files, but also
helps de-clutter them.


### _class_ pymatgen.io.cp2k.inputs.SectionList(sections: Sequence[Section])
Bases: `MSONable`

Section list.

Initializes a SectionList object using a sequence of sections.


* **Parameters**

    **sections** – A list of keywords. Must all have the same name (case-insensitive)



#### append(item)
Append the section list.


#### extend(lst: list)
Extend the section list.


#### get(d, index=-1)
Get for section list. If index is specified, return the section at that index.
Otherwise, return a get on the last section.


#### get_string()
Return string representation of section list.


#### verbosity(verbosity)
Silence all sections in section list.


### _class_ pymatgen.io.cp2k.inputs.Smear(elec_temp: int | float = 300, method: str = 'FERMI_DIRAC', fixed_magnetic_moment: float = -100.0, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.Subsys(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
Bases: `Section`

Controls the definition of the system to be simulated.

Initialize the subsys section.


### _class_ pymatgen.io.cp2k.inputs.V_Hartree_Cube(keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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



### _class_ pymatgen.io.cp2k.inputs.Xc_Functional(functionals: Iterable | None = None, keywords: dict | None = None, subsections: dict | None = None, \*\*kwargs)
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