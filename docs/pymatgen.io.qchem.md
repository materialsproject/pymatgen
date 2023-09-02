---
layout: default
title: pymatgen.io.qchem.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.qchem package

This package implements modules for input and output to and from Qchem.


## pymatgen.io.qchem.inputs module

Classes for reading/manipulating/writing QChem input files.


### _class_ QCInput(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule) | list[[Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule)] | Literal['read'], rem: dict, opt: dict[str, list] | None = None, pcm: dict | None = None, solvent: dict | None = None, smx: dict | None = None, scan: dict[str, list] | None = None, van_der_waals: dict[str, float] | None = None, vdw_mode: str = 'atomic', plots: dict | None = None, nbo: dict | None = None, geom_opt: dict | None = None, cdft: list[list[dict]] | None = None, almo_coupling: list[list[tuple[int, int]]] | None = None, svp: dict | None = None, pcm_nonels: dict | None = None)
Bases: [`InputFile`](pymatgen.io.md#pymatgen.io.core.InputFile)

An object representing a QChem input file. QCInput attributes represent different sections of a QChem input file.
To add a new section one needs to modify __init__, __str__, from_sting and add static methods
to read and write the new section i.e. section_template and read_section. By design, there is very little (or no)
checking that input parameters conform to the appropriate QChem format, this responsible lands on the user or a
separate error handling software.


* **Parameters**


    * **molecule** (*pymatgen Molecule object**, **list** of **Molecule objects**, or **"read"*) – Input molecule(s). molecule can be set as a pymatgen Molecule object, a list of such
    Molecule objects, or as the string “read”. “read” can be used in multi_job QChem input
    files where the molecule is read in from the previous calculation.


    * **rem** (*dict*) – A dictionary of all the input parameters for the rem section of QChem input file.
    Ex. rem = {‘method’: ‘rimp2’, ‘basis’: ‘6-31\*G++’ … }


    * **opt** (*dict** of **lists*) – A dictionary of opt sections, where each opt section is a key and the corresponding
    values are a list of strings. Strings must be formatted as instructed by the QChem manual.
    The different opt sections are: CONSTRAINT, FIXED, DUMMY, and CONNECT
    Ex. opt = {“CONSTRAINT”: [“tors 2 3 4 5 25.0”, “tors 2 5 7 9 80.0”], “FIXED”: [“2 XY”]}


    * **pcm** (*dict*) – A dictionary of the PCM section, defining behavior for use of the polarizable continuum model.
    Ex: pcm = {“theory”: “cpcm”, “hpoints”: 194}


    * **solvent** (*dict*) – A dictionary defining the solvent parameters used with PCM.
    Ex: solvent = {“dielectric”: 78.39, “temperature”: 298.15}


    * **smx** (*dict*) – A dictionary defining solvent parameters used with the SMD method, a solvent method that adds
    short-range terms to PCM.
    Ex: smx = {“solvent”: “water”}


    * **scan** (*dict** of **lists*) – A dictionary of scan variables. Because two constraints of the same type are allowed (for instance, two
    torsions or two bond stretches), each TYPE of variable (stre, bend, tors) should be its own key in the
    dict, rather than each variable. Note that the total number of variable (sum of lengths of all lists)
    CANNOT be
    more than two.
    Ex. scan = {“stre”: [“3 6 1.5 1.9 0.1”], “tors”: [“1 2 3 4 -180 180 15”]}


    * **van_der_waals** (*dict*) – A dictionary of custom van der Waals radii to be used when constructing cavities for the PCM
    model or when computing, e.g. Mulliken charges. They keys are strs whose meaning depends on
    the value of vdw_mode, and the values are the custom radii in angstroms.


    * **vdw_mode** (*str*) – Method of specifying custom van der Waals radii - ‘atomic’ or ‘sequential’.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., 12 = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **plots** (*dict*) – A dictionary of all the input parameters for the plots section of the QChem input file.


    * **nbo** (*dict*) – A dictionary of all the input parameters for the nbo section of the QChem input file.


    * **geom_opt** (*dict*) – A dictionary of input parameters for the geom_opt section of the QChem input file.
    This section is required when using the new libopt3 geometry optimizer.


    * **cdft** (*list** of **lists** of **dicts*) – A list of lists of dictionaries, where each dictionary represents a charge constraint in the
    cdft section of the QChem input file.

    Each entry in the main list represents one state (allowing for multi-configuration calculations
    using constrained density functional theory - configuration interaction (CDFT-CI).
    Each state is represented by a list, which itself contains some number of constraints
    (dictionaries).

    Ex:


        1. For a single-state calculation with two constraints:

    > cdft=[[

    >     {“value”: 1.0, “coefficients”: [1.0], “first_atoms”: [1], “last_atoms”: [2], “types”: [None]},
    >     {“value”: 2.0, “coefficients”: [1.0, -1.0], “first_atoms”: [1, 17], “last_atoms”: [3, 19],

    >     > ”types”: [“s”]}

    ]]

    Note that a type of None will default to a charge constraint (which can also be accessed by
    requesting a type of “c” or “charge”.

    2. For a multi-reference calculation:
    cdft=[

    > [

    >     {“value”: 1.0, “coefficients”: [1.0], “first_atoms”: [1], “last_atoms”: [27],

    >         ”types”: [“c”]},

    >     {“value”: 0.0, “coefficients”: [1.0], “first_atoms”: [1], “last_atoms”: [27],

    >         ”types”: [“s”]},

    > ],
    > [

    > > {“value”: 0.0, “coefficients”: [1.0], “first_atoms”: [1], “last_atoms”: [27],

    > >     ”types”: [“c”]},

    > > {“value”: -1.0, “coefficients”: [1.0], “first_atoms”: [1], “last_atoms”: [27],

    > >     ”types”: [“s”]},

    > ]

    ]



    * **almo_coupling** (*list** of **lists** of **int 2-tuples*) – A list of lists of int 2-tuples used for calculations of diabatization and state coupling calculations

        relying on the absolutely localized molecular orbitals (ALMO) methodology. Each entry in the main
        list represents a single state (two states are included in an ALMO calculation). Within a single
        state, each 2-tuple represents the charge and spin multiplicity of a single fragment.

    ex: almo=[

        > [

        >     (1, 2),
        >     (0, 1)

        > ],
        > [

        > > (0, 1),
        > > (1, 2)

        > ]

        ]




#### _static_ almo_template(almo_coupling: list[list[tuple[int, int]]])

* **Parameters**

    **almo** – list of lists of int 2-tuples.



* **Returns**

    (str)



#### _static_ cdft_template(cdft: list[list[dict]])

* **Parameters**

    **cdft** – list of lists of dicts.



* **Returns**

    (str)



#### _static_ find_sections(string: str)
Find sections in the string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    List of sections.



#### _static_ from_file(filename: str | Path)
Create QcInput from file.


* **Parameters**

    **filename** (*str*) – Filename



* **Returns**

    QcInput



#### _classmethod_ from_multi_jobs_file(filename: str)
Create list of QcInput from a file.


* **Parameters**

    **filename** (*str*) – Filename



* **Returns**

    List of QCInput objects



#### _classmethod_ from_str(string: str)
Read QcInput from string.


* **Parameters**

    **string** (*str*) – String input.



* **Returns**

    QcInput



#### _static_ geom_opt_template(geom_opt: dict)

* **Parameters**

    **(****)** (*geom_opt*) –



* **Returns**

    (str) geom_opt parameters.



#### get_str()
Return a string representation of an entire input file.


#### get_string(\*\*kwds)
get_string is deprecated!
Use get_str instead


#### _static_ molecule_template(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule) | list[[Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule)] | Literal['read'])

* **Parameters**

    **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)*, **list** of **Molecules**, or **"read"*) –



* **Returns**

    (str) Molecule template.



#### _static_ multi_job_string(job_list: list[QCInput])

* **Parameters**

    **(****)** (*job_list*) – List of jobs.



* **Returns**

    (str) String representation of multi job input file.



#### _static_ nbo_template(nbo: dict)

* **Parameters**

    **(****)** (*nbo*) –



* **Returns**

    (str)



#### _static_ opt_template(opt: dict[str, list])
Optimization template.


* **Parameters**

    **(****)** (*opt*) –



* **Returns**

    (str)



#### _static_ pcm_nonels_template(pcm_nonels: dict)
Template for the $pcm_nonels section.

Arg

    pcm_nonels: dict of CMIRS parameters, e.g.
    {

    > “a”: “-0.006736”,
    > “b”: “0.032698”,
    > “c”: “-1249.6”,
    > “d”: “-21.405”,
    > “gamma”: “3.7”,
    > “solvrho”: “0.05”,
    > “delta”: 7,
    > “gaulag_n”: 40,

    }


* **Returns**

    (str)



#### _static_ pcm_template(pcm: dict)
Pcm run template.


* **Parameters**

    **(****)** (*pcm*) –



* **Returns**

    (str)



#### _static_ plots_template(plots: dict)

* **Parameters**

    **(****)** (*plots*) –



* **Returns**

    (str)



#### _static_ read_almo(string: str)
Read ALMO coupling parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (list of lists of int 2-tuples) almo_coupling parameters



#### _static_ read_cdft(string: str)
Read cdft parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    cdft parameters



* **Return type**

    list[list[dict]]



#### _static_ read_geom_opt(string: str)
Read geom_opt parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) geom_opt parameters.



#### _static_ read_molecule(string: str)
Read molecule from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    Molecule



#### _static_ read_nbo(string: str)
Read nbo parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) nbo parameters.



#### _static_ read_opt(string: str)
Read opt section from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) Opt section



#### _static_ read_pcm(string: str)
Read pcm parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) PCM parameters



#### _static_ read_pcm_nonels(string: str)
Read pcm_nonels parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) PCM parameters



#### _static_ read_plots(string: str)
Read plots parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) plots parameters.



#### _static_ read_rem(string: str)
Parse rem from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) rem



#### _static_ read_scan(string: str)
Read scan section from a string.


* **Parameters**

    **string** – String to be parsed



* **Returns**

    Dict representing Q-Chem scan section



#### _static_ read_smx(string: str)
Read smx parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) SMX parameters.



#### _static_ read_solvent(string: str)
Read solvent parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (dict) Solvent parameters



#### _static_ read_svp(string: str)
Read svp parameters from string.


#### _static_ read_vdw(string: str)
Read van der Waals parameters from string.


* **Parameters**

    **string** (*str*) – String



* **Returns**

    (str, dict) vdW mode (‘atomic’ or ‘sequential’) and dict of van der Waals radii.



#### _static_ rem_template(rem: dict)

* **Parameters**

    **(****)** (*rem*) –



* **Returns**

    (str)



#### _static_ scan_template(scan: dict[str, list])

* **Parameters**

    **scan** (*dict*) – Dictionary with scan section information.
    Ex: {“stre”: [“3 6 1.5 1.9 0.1”], “tors”: [“1 2 3 4 -180 180 15”]}.



* **Returns**

    String representing Q-Chem input format for scan section



#### _static_ smx_template(smx: dict)

* **Parameters**

    **(****)** (*smx*) –



* **Returns**

    (str)



#### _static_ solvent_template(solvent: dict)
Solvent template.


* **Parameters**

    **(****)** (*solvent*) –



* **Returns**

    (str)



#### _static_ svp_template(svp: dict)
Template for the $svp section.


* **Parameters**


    * **svp** – dict of SVP parameters, e.g.


    * **{"rhoiso"** – “0.001”, “nptleb”: “1202”, “itrngr”: “2”, “irotgr”: “2”}



* **Returns**

    the $svp section. Note that all parameters will be concatenated onto

        a single line formatted as a FORTRAN namelist. This is necessary
        because the isodensity SS(V)PE model in Q-Chem calls a secondary code.




* **Return type**

    str



#### _static_ van_der_waals_template(radii: dict[str, float], mode: str = 'atomic')

* **Parameters**


    * **radii** (*dict*) – Dictionary with custom van der Waals radii, in
    Angstroms, keyed by either atomic number or sequential
    atom number (see ‘mode’ kwarg).
    Ex: {1: 1.20, 12: 1.70}


    * **mode** – ‘atomic’ or ‘sequential’. In ‘atomic’ mode (default), dict keys
    represent the atomic number associated with each radius (e.g., ‘12’ = carbon).
    In ‘sequential’ mode, dict keys represent the sequential position of
    a single specific atom in the input structure.
    **NOTE: keys must be given as strings even though they are numbers!**.



* **Returns**

    String representing Q-Chem input format for van_der_waals section



#### _static_ write_multi_job_file(job_list: list[QCInput], filename: str)
Write a multijob file.


* **Parameters**


    * **(****)** (*filename*) – List of jobs.


    * **(****)** – Filename


## pymatgen.io.qchem.outputs module

Parsers for Qchem output files.


### _class_ QCOutput(filename: str)
Bases: `MSONable`

Class to parse QChem output files.


* **Parameters**

    **filename** (*str*) – Filename to parse.



#### _check_completion_errors()
Parses potential errors that can cause jobs to crash.


#### _detect_general_warnings()

#### _get_grad_format_length(header)
Determines the maximum number of gradient entries printed on a line,
which changes for different versions of Q-Chem.


#### _read_SCF()
Parses both old and new SCFs.


#### _read_almo_msdft()
Parse output of ALMO(MSDFT) calculations for coupling between diabatic states.


#### _read_cdft()
Parses output from charge- or spin-constrained DFT (CDFT) calculations.


#### _read_charge_and_multiplicity()
Parses charge and multiplicity.


#### _read_charges_and_dipoles()
Parses Mulliken/ESP/RESP charges.
Parses associated dipoles.
Also parses spins given an unrestricted SCF.


#### _read_cmirs_information()
Parses information from CMIRS solvent calculations.

In addition to the 5 energies returned by ISOSVP (and read separately in
_read_isosvp_information), there are 4 additional energies reported, as shown
in the example below

### The Final SS(V)PE energies and Properties

#### Energies

The Final Solution-Phase Energy =     -40.4751881546
The Solute Internal Energy =          -40.4748568841
The Change in Solute Internal Energy =  0.0000089729  (   0.00563 KCAL/MOL)
The Reaction Field Free Energy =       -0.0003312705  (  -0.20788 KCAL/MOL)
The Dispersion Energy =                 0.6955550107  (  -2.27836 KCAL/MOL)
The Exchange Energy =                   0.2652679507  (   2.15397 KCAL/MOL)
Min. Negative Field Energy =            0.0005235850  (   0.00000 KCAL/MOL)
Max. Positive Field Energy =            0.0179866718  (   0.00000 KCAL/MOL)
The Total Solvation Free Energy =      -0.0005205275  (  -0.32664 KCAL/MOL)


#### _read_coefficient_matrix()
Parses the coefficient matrix from the output file. Done is much
the same was as the Fock matrix.


#### _read_eigenvalues()
Parse the orbital energies from the output file. An array of the
dimensions of the number of orbitals used in the calculation is stored.


#### _read_fock_matrix()
Parses the Fock matrix. The matrix is read in whole
from the output file and then transformed into the right dimensions.


#### _read_force_data()

#### _read_frequency_data()
Parses cpscf_nseg, frequencies, enthalpy, entropy, and mode vectors.


#### _read_geometries()
Parses all geometries from an optimization trajectory.


#### _read_gradients()
Parses all gradients obtained during an optimization trajectory.


#### _read_isosvp_information()
Parses information from ISOSVP solvent calculations.

There are 5 energies output, as in the example below

### The Final SS(V)PE energies and Properties

#### Energies

The Final Solution-Phase Energy =     -40.4850599390
The Solute Internal Energy =          -40.4846329759
The Change in Solute Internal Energy =  0.0000121970  (   0.00765 KCAL/MOL)
The Reaction Field Free Energy =       -0.0004269631  (  -0.26792 KCAL/MOL)
The Total Solvation Free Energy =      -0.0004147661  (  -0.26027 KCAL/MOL)

In addition, we need to parse the DIELST fortran variable to get the dielectric
constant used.


#### _read_nbo_data()
Parses NBO output.


#### _read_optimization_data()

#### _read_pcm_information()
Parses information from PCM solvent calculations.


#### _read_scan_data()

#### _read_smd_information()
Parses information from SMD solvent calculations.


#### _read_species_and_inital_geometry()
Parses species and initial geometry.


#### as_dict()

* **Returns**

    MSONable dict.



#### _static_ multiple_outputs_from_file(filename, keep_sub_files=True)
Parses a QChem output file with multiple calculations
# 1.) Separates the output into sub-files

> e.g. qcout -> qcout.0, qcout.1, qcout.2 … qcout.N
> a.) Find delimiter for multiple calculations
> b.) Make separate output sub-files

2.) Creates separate QCCalcs for each one from the sub-files.


### check_for_structure_changes(mol1: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), mol2: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Compares connectivity of two molecules (using MoleculeGraph w/ OpenBabelNN).
This function will work with two molecules with different atom orderings,

> but for proper treatment, atoms should be listed in the same order.

Possible outputs include:
- no_change: the bonding in the two molecules is identical
- unconnected_fragments: the MoleculeGraph of mol1 is connected, but the

> MoleculeGraph is mol2 is not connected


* fewer_bonds: the MoleculeGraph of mol1 has more bonds (edges) than the
MoleculeGraph of mol2


* more_bonds: the MoleculeGraph of mol2 has more bonds (edges) than the
MoleculeGraph of mol1


* bond_change: this case catches any other non-identical MoleculeGraphs


* **Parameters**


    * **mol1** – Pymatgen Molecule object to be compared.


    * **mol2** – Pymatgen Molecule object to be compared.



* **Returns**

    One of [“unconnected_fragments”, “fewer_bonds”, “more_bonds”,
    “bond_change”, “no_change”]



### get_percentage(line: str, orbital: str)
Retrieve the percent character of an orbital.


* **Parameters**


    * **line** – Line containing orbital and percentage.


    * **orbital** – Type of orbital (s, p, d, f).



* **Returns**

    Percentage of character.



* **Raises**

    **n/a** –



### jump_to_header(lines: list[str], header: str)
Given a list of lines, truncate the start of the list so that the first line
of the new list contains the header.


* **Parameters**


    * **lines** – List of lines.


    * **header** – Substring to match.



* **Returns**

    Truncated lines.



* **Raises**

    **RuntimeError** –



### nbo_parser(filename: str)
Parse all the important sections of NBO output.


* **Parameters**

    **filename** – Path to QChem NBO output.



* **Returns**

    Data frames of formatted output.



* **Raises**

    **RuntimeError** –



### parse_hybridization_character(lines: list[str])
Parse the hybridization character section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frames of formatted output.



* **Raises**

    **RuntimeError** –



### parse_hyperbonds(lines: list[str])
Parse the natural populations section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frame of formatted output.



* **Raises**

    **RuntimeError** –



### parse_natural_populations(lines: list[str])
Parse the natural populations section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frame of formatted output.



* **Raises**

    **RuntimeError** –



### parse_perturbation_energy(lines: list[str])
Parse the perturbation energy section of NBO output.


* **Parameters**

    **lines** – QChem output lines.



* **Returns**

    Data frame of formatted output.



* **Raises**

    **RuntimeError** –



### z_int(string: str)
Convert string to integer.
If string empty, return -1.


* **Parameters**

    **string** – Input to be cast to int.



* **Returns**

    Int representation.



* **Raises**

    **n/a** –


## pymatgen.io.qchem.sets module

Input sets for Qchem.


### _class_ ForceSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-tzvpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', cdft_constraints: list[list[dict]] | None = None, overwrite_inputs: dict | None = None)
Bases: `QChemDictSet`

QChemDictSet for a force (gradient) calculation.


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) –


    * **basis_set** (*str*) – Basis set to use. (Default: “def2-tzvpd”)


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.
    (Default: “diis”)


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2).

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **cdft_constraints** (*list** of **lists** of **dicts*) – A list of lists of dictionaries, where each dictionary represents a charge
    constraint in the cdft section of the QChem input file.

    Each entry in the main list represents one state (allowing for multi-configuration
    calculations using constrained density functional theory - configuration interaction
    (CDFT-CI). Each state is represented by a list, which itself contains some number of
    constraints (dictionaries).

    Ex:


        1. For a single-state calculation with two constraints:

    > cdft_constraints=[[

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [2],
    >         “types”: [None]

    >     },
    >     {

    >     > ”value”: 2.0,
    >     > “coefficients”: [1.0, -1.0],
    >     > “first_atoms”: [1, 17],
    >     > “last_atoms”: [3, 19],
    >     > “types”: [“s”]

    >     }

    ]]

    Note that a type of None will default to a charge constraint (which can also be
    accessed by requesting a type of “c” or “charge”).

    2. For a CDFT-CI multi-reference calculation:
    cdft_constraints=[

    > [

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [27],
    >         “types”: [“c”]

    >     },
    >     {

    >     > ”value”: 0.0,
    >     > “coefficients”: [1.0],
    >     > “first_atoms”: [1],
    >     > “last_atoms”: [27],
    >     > “types”: [“s”]

    >     },

    > ],
    > [

    > > {

    > >     “value”: 0.0,
    > >     “coefficients”: [1.0],
    > >     “first_atoms”: [1],
    > >     “last_atoms”: [27],
    > >     “types”: [“c”]

    > > },
    > > {

    > > > ”value”: -1.0,
    > > > “coefficients”: [1.0],
    > > > “first_atoms”: [1],
    > > > “last_atoms”: [27],
    > > > “types”: [“s”]

    > > },

    > ]

    ]



    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**




### _class_ FreqSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', cdft_constraints: list[list[dict]] | None = None, overwrite_inputs: dict | None = None)
Bases: `QChemDictSet`

QChemDictSet for a frequency calculation.


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) –


    * **basis_set** (*str*) – Basis set to use. (Default: “def2-svpd”)


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.
    (Default: “diis”)


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2).

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **cdft_constraints** (*list** of **lists** of **dicts*) – A list of lists of dictionaries, where each dictionary represents a charge
    constraint in the cdft section of the QChem input file.

    Each entry in the main list represents one state (allowing for multi-configuration
    calculations using constrained density functional theory - configuration interaction
    (CDFT-CI). Each state is represented by a list, which itself contains some number of
    constraints (dictionaries).

    Ex:


        1. For a single-state calculation with two constraints:

    > cdft_constraints=[[

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [2],
    >         “types”: [None]

    >     },
    >     {

    >     > ”value”: 2.0,
    >     > “coefficients”: [1.0, -1.0],
    >     > “first_atoms”: [1, 17],
    >     > “last_atoms”: [3, 19],
    >     > “types”: [“s”]

    >     }

    ]]

    Note that a type of None will default to a charge constraint (which can also be
    accessed by requesting a type of “c” or “charge”).

    2. For a CDFT-CI multi-reference calculation:
    cdft_constraints=[

    > [

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [27],
    >         “types”: [“c”]

    >     },
    >     {

    >     > ”value”: 0.0,
    >     > “coefficients”: [1.0],
    >     > “first_atoms”: [1],
    >     > “last_atoms”: [27],
    >     > “types”: [“s”]

    >     },

    > ],
    > [

    > > {

    > >     “value”: 0.0,
    > >     “coefficients”: [1.0],
    > >     “first_atoms”: [1],
    > >     “last_atoms”: [27],
    > >     “types”: [“c”]

    > > },
    > > {

    > > > ”value”: -1.0,
    > > > “coefficients”: [1.0],
    > > > “first_atoms”: [1],
    > > > “last_atoms”: [27],
    > > > “types”: [“s”]

    > > },

    > ]

    ]



    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**




### _class_ OptSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, opt_variables: dict[str, list] | None = None, geom_opt_max_cycles: int = 200, geom_opt: dict | None = None, cdft_constraints: list[list[dict]] | None = None, overwrite_inputs: dict | None = None)
Bases: `QChemDictSet`

QChemDictSet for a geometry optimization.


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) –


    * **job_type** (*str*) – QChem job type to run. Valid options are “opt” for optimization,
    “sp” for single point, “freq” for frequency calculation, or “force” for
    force evaluation.


    * **basis_set** (*str*) – Basis set to use. (Default: “def2-svpd”)


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.
    (Default: “diis”)


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2).

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **geom_opt_max_cycles** (*int*) – Maximum number of geometry optimization iterations. (Default: 200)


    * **geom_opt** (*dict*) – A dict containing parameters for the $geom_opt section of the Q-Chem input
    file, which control the new geometry optimizer available starting in version 5.4.2. The
    new optimizer remains under development but was officially released and became the default
    optimizer in Q-Chem version 6.0.0. Note that for version 5.4.2, the new optimizer must be
    explicitly requested by passing in a dictionary (empty or otherwise) for this input parameter.
    (Default: False)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **cdft_constraints** (*list** of **lists** of **dicts*) – A list of lists of dictionaries, where each dictionary represents a charge
    constraint in the cdft section of the QChem input file.

    Each entry in the main list represents one state (allowing for multi-configuration
    calculations using constrained density functional theory - configuration interaction
    (CDFT-CI). Each state is represented by a list, which itself contains some number of
    constraints (dictionaries).

    Ex:


        1. For a single-state calculation with two constraints:

    > cdft_constraints=[[

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [2],
    >         “types”: [None]

    >     },
    >     {

    >     > ”value”: 2.0,
    >     > “coefficients”: [1.0, -1.0],
    >     > “first_atoms”: [1, 17],
    >     > “last_atoms”: [3, 19],
    >     > “types”: [“s”]

    >     }

    ]]

    Note that a type of None will default to a charge constraint (which can also be
    accessed by requesting a type of “c” or “charge”).

    2. For a CDFT-CI multi-reference calculation:
    cdft_constraints=[

    > [

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [27],
    >         “types”: [“c”]

    >     },
    >     {

    >     > ”value”: 0.0,
    >     > “coefficients”: [1.0],
    >     > “first_atoms”: [1],
    >     > “last_atoms”: [27],
    >     > “types”: [“s”]

    >     },

    > ],
    > [

    > > {

    > >     “value”: 0.0,
    > >     “coefficients”: [1.0],
    > >     “first_atoms”: [1],
    > >     “last_atoms”: [27],
    > >     “types”: [“c”]

    > > },
    > > {

    > > > ”value”: -1.0,
    > > > “coefficients”: [1.0],
    > > > “first_atoms”: [1],
    > > > “last_atoms”: [27],
    > > > “types”: [“s”]

    > > },

    > ]

    ]



    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**



    * **vdw_mode** – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.



### _class_ PESScanSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, opt_variables: dict[str, list] | None = None, scan_variables: dict[str, list] | None = None, overwrite_inputs: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic')
Bases: `QChemDictSet`

QChemDictSet for a potential energy surface scan (PES_SCAN) calculation,
used primarily to identify possible transition states or to sample different
geometries.
Note: Because there are no defaults that can be used for a PES scan (the
variables are completely dependent on the molecular structure), by default
scan_variables = None. However, a PES Scan job should not be run with less
than one variable (or more than two variables).


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) –


    * **opt_variables** (*dict*) – A dictionary of opt sections, where each opt section is a key
    and the corresponding values are a list of strings. Strings must be formatted
    as instructed by the QChem manual. The different opt sections are: CONSTRAINT, FIXED,
    DUMMY, and CONNECT.

    Ex. opt = {“CONSTRAINT”: [“tors 2 3 4 5 25.0”, “tors 2 5 7 9 80.0”], “FIXED”: [“2 XY”]}



    * **scan_variables** (*dict*) – A dictionary of scan variables. Because two constraints of the
    same type are allowed (for instance, two torsions or two bond stretches), each TYPE of
    variable (stre, bend, tors) should be its own key in the dict, rather than each variable.
    Note that the total number of variable (sum of lengths of all lists) CANNOT be more than two.

    Ex. scan_variables = {“stre”: [“3 6 1.5 1.9 0.1”], “tors”: [“1 2 3 4 -180 180 15”]}



    * **basis_set** (*str*) – Basis set to use. (Default: “def2-svpd”)


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.
    (Default: “diis”)


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2)

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**



    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies only if
    you are using overwrite_inputs to add a $van_der_waals section to the input. In ‘atomic’ mode
    (default), dict keys represent the atomic number associated with each radius (e.g., ‘12’ = carbon).
    In ‘sequential’ mode, dict keys represent the sequential position of a single
    specific atom in the input structure.



### _class_ QChemDictSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), job_type: str, basis_set: str, scf_algorithm: str, qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, opt_variables: dict[str, list] | None = None, scan_variables: dict[str, list] | None = None, max_scf_cycles: int = 100, geom_opt_max_cycles: int = 200, plot_cubes: bool = False, nbo_params: dict | None = None, geom_opt: dict | None = None, cdft_constraints: list[list[dict]] | None = None, almo_coupling_states: list[list[tuple[int, int]]] | None = None, overwrite_inputs: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', extra_scf_print: bool = False)
Bases: `QCInput`

Build a QCInput given all the various input parameters. Can be extended by standard implementations below.


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) – Molecule to run QChem on.


    * **job_type** (*str*) – QChem job type to run. Valid options are “opt” for optimization,
    “sp” for single point, “freq” for frequency calculation, or “force” for
    force evaluation.


    * **basis_set** (*str*) – Basis set to use. For example, “def2-tzvpd”.


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2).

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **opt_variables** (*dict*) – A dictionary of opt sections, where each opt section is a key
    and the corresponding values are a list of strings. Strings must be formatted
    as instructed by the QChem manual. The different opt sections are: CONSTRAINT, FIXED,
    DUMMY, and CONNECT.

    Ex. opt = {“CONSTRAINT”: [“tors 2 3 4 5 25.0”, “tors 2 5 7 9 80.0”], “FIXED”: [“2 XY”]}



    * **scan_variables** (*dict*) – A dictionary of scan variables. Because two constraints of the
    same type are allowed (for instance, two torsions or two bond stretches), each TYPE of
    variable (stre, bend, tors) should be its own key in the dict, rather than each variable.
    Note that the total number of variable (sum of lengths of all lists) CANNOT be more than two.

    Ex. scan_variables = {“stre”: [“3 6 1.5 1.9 0.1”], “tors”: [“1 2 3 4 -180 180 15”]}



    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **geom_opt_max_cycles** (*int*) – Maximum number of geometry optimization iterations. (Default: 200)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **nbo_params** (*dict*) – A dict containing the desired NBO params. Note that a key:value pair of
    “version”:7 will trigger NBO7 analysis. Otherwise, NBO5 analysis will be performed,
    including if an empty dict is passed. Besides a key of “version”, all other key:value
    pairs will be written into the $nbo section of the QChem input file. (Default: False)


    * **geom_opt** (*dict*) – A dict containing parameters for the $geom_opt section of the Q-Chem input
    file, which control the new geometry optimizer available starting in version 5.4.2. The
    new optimizer remains under development but was officially released and became the default
    optimizer in Q-Chem version 6.0.0. Note that for version 5.4.2, the new optimizer must be
    explicitly requested by passing in a dictionary (empty or otherwise) for this input parameter.
    (Default: False)


    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.

    > cdft_constraints (list of lists of dicts):

    A list of lists of dictionaries, where each dictionary represents a charge
    constraint in the cdft section of the QChem input file.

    Each entry in the main list represents one state (allowing for multi-configuration
    calculations using constrained density functional theory - configuration interaction
    (CDFT-CI). Each state is represented by a list, which itself contains some number of
    constraints (dictionaries).

    Ex:


        1. For a single-state calculation with two constraints:

    > cdft_constraints=[[

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [2],
    >         “types”: [None]

    >     },
    >     {

    >     > ”value”: 2.0,
    >     > “coefficients”: [1.0, -1.0],
    >     > “first_atoms”: [1, 17],
    >     > “last_atoms”: [3, 19],
    >     > “types”: [“s”]

    >     }

    ]]

    Note that a type of None will default to a charge constraint (which can also be
    accessed by requesting a type of “c” or “charge”).

    2. For a CDFT-CI multi-reference calculation:
    cdft_constraints=[

    > [

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [27],
    >         “types”: [“c”]

    >     },
    >     {

    >     > ”value”: 0.0,
    >     > “coefficients”: [1.0],
    >     > “first_atoms”: [1],
    >     > “last_atoms”: [27],
    >     > “types”: [“s”]

    >     },

    > ],
    > [

    > > {

    > >     “value”: 0.0,
    > >     “coefficients”: [1.0],
    > >     “first_atoms”: [1],
    > >     “last_atoms”: [27],
    > >     “types”: [“c”]

    > > },
    > > {

    > > > ”value”: -1.0,
    > > > “coefficients”: [1.0],
    > > > “first_atoms”: [1],
    > > > “last_atoms”: [27],
    > > > “types”: [“s”]

    > > },

    > ]

    ]



    * **cdft_constraints** (*list**[**list**[**dict**]**]*) – A list of lists of dictionaries, where each


    * **almo_coupling_states** (*list** of **lists** of **int 2-tuples*) – A list of lists of int 2-tuples used for calculations of diabatization and state
    coupling calculations relying on the absolutely localized molecular orbitals (ALMO)
    methodology. Each entry in the main list represents a single state (two states are
    included in an ALMO calculation). Within a single state, each 2-tuple represents the
    charge and spin multiplicity of a single fragment.
    ex: almo_coupling_states=[

    > > [

    > >     (1, 2),
    > >     (0, 1)

    > > ],
    > > [

    > > > (0, 1),
    > > > (1, 2)

    > > ]

    > ]



    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**



    * **vdw_mode** – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **extra_scf_print** (*bool*) – Whether to store extra information generated from the SCF
    cycle. If switched on, the Fock Matrix, coefficients of MO and the density matrix
    will be stored.



#### write(input_file: str)

* **Parameters**

    **input_file** (*str*) – Filename.



### _class_ SinglePointSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-tzvpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', cdft_constraints: list[list[dict]] | None = None, almo_coupling_states: list[list[tuple[int, int]]] | None = None, extra_scf_print: bool = False, overwrite_inputs: dict | None = None)
Bases: `QChemDictSet`

QChemDictSet for a single point calculation.


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) –


    * **job_type** (*str*) – QChem job type to run. Valid options are “opt” for optimization,
    “sp” for single point, “freq” for frequency calculation, or “force” for
    force evaluation.


    * **basis_set** (*str*) – Basis set to use. (Default: “def2-tzvpd”)


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.
    (Default: “diis”)


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2).

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters. Note that due to limitations in Q-Chem, use of the ISOSVP
    or CMIRS solvent models will disable the GEN_SCFMAN algorithm, which may limit compatible choices
    for scf_algorithm.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics. Note also that due to limitations in Q-Chem, use of the ISOSVP
    or CMIRS solvent models will disable the GEN_SCFMAN algorithm, which may limit compatible choices
    for scf_algorithm.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **cdft_constraints** (*list** of **lists** of **dicts*) – A list of lists of dictionaries, where each dictionary represents a charge
    constraint in the cdft section of the QChem input file.

    Each entry in the main list represents one state (allowing for multi-configuration
    calculations using constrained density functional theory - configuration interaction
    (CDFT-CI). Each state is represented by a list, which itself contains some number of
    constraints (dictionaries).

    Ex:


        1. For a single-state calculation with two constraints:

    > cdft_constraints=[[

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [2],
    >         “types”: [None]

    >     },
    >     {

    >     > ”value”: 2.0,
    >     > “coefficients”: [1.0, -1.0],
    >     > “first_atoms”: [1, 17],
    >     > “last_atoms”: [3, 19],
    >     > “types”: [“s”]

    >     }

    ]]

    Note that a type of None will default to a charge constraint (which can also be
    accessed by requesting a type of “c” or “charge”).

    2. For a CDFT-CI multi-reference calculation:
    cdft_constraints=[

    > [

    >     {

    >         “value”: 1.0,
    >         “coefficients”: [1.0],
    >         “first_atoms”: [1],
    >         “last_atoms”: [27],
    >         “types”: [“c”]

    >     },
    >     {

    >     > ”value”: 0.0,
    >     > “coefficients”: [1.0],
    >     > “first_atoms”: [1],
    >     > “last_atoms”: [27],
    >     > “types”: [“s”]

    >     },

    > ],
    > [

    > > {

    > >     “value”: 0.0,
    > >     “coefficients”: [1.0],
    > >     “first_atoms”: [1],
    > >     “last_atoms”: [27],
    > >     “types”: [“c”]

    > > },
    > > {

    > > > ”value”: -1.0,
    > > > “coefficients”: [1.0],
    > > > “first_atoms”: [1],
    > > > “last_atoms”: [27],
    > > > “types”: [“s”]

    > > },

    > ]

    ]



    * **almo_coupling_states** (*list** of **lists** of **int 2-tuples*) – A list of lists of int 2-tuples used for calculations of diabatization and state
    coupling calculations relying on the absolutely localized molecular orbitals (ALMO)
    methodology. Each entry in the main list represents a single state (two states are
    included in an ALMO calculation). Within a single state, each 2-tuple represents the
    charge and spin multiplicity of a single fragment.
    ex: almo_coupling_states=[

    > > [

    > >     (1, 2),
    > >     (0, 1)

    > > ],
    > > [

    > > > (0, 1),
    > > > (1, 2)

    > > ]

    > ]



    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**



    * **vdw_mode** – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


    * **extra_scf_print** (*bool*) – Whether to store extra information generated from the SCF
    cycle. If switched on, the Fock Matrix, coefficients of MO and the density matrix
    will be stored.



### _class_ TransitionStateSet(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, opt_variables: dict[str, list] | None = None, geom_opt_max_cycles: int = 200, geom_opt: dict | None = None, overwrite_inputs: dict | None = None, vdw_mode='atomic')
Bases: `QChemDictSet`

QChemDictSet for a transition-state search.


* **Parameters**


    * **molecule** (*Pymatgen Molecule object*) –


    * **basis_set** (*str*) – Basis set to use. (Default: “def2-svpd”)


    * **scf_algorithm** (*str*) – Algorithm to use for converging the SCF. Recommended choices are
    “DIIS”, “GDM”, and “DIIS_GDM”. Other algorithms supported by Qchem’s GEN_SCFMAN
    module will also likely perform well. Refer to the QChem manual for further details.
    (Default: “diis”)


    * **qchem_version** (*int*) – Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)


    * **dft_rung** (*int*) – Select the rung on “Jacob’s Ladder of Density Functional Approximations” in
    order of increasing accuracy/cost. For each rung, we have prescribed one functional based
    on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
    1 (LSDA) = SPW92
    2 (GGA) = B97-D3(BJ)
    3 (metaGGA) = B97M-V
    4 (hybrid metaGGA) = ωB97M-V
    5 (double hybrid metaGGA) = ωB97M-(2).

    (Default: 4)

    To set a functional not given by one of the above, set the overwrite_inputs
    argument to {“method”:”<NAME OF FUNCTIONAL>”}



    * **pcm_dielectric** (*float*) – Dielectric constant to use for PCM implicit solvation model. (Default: None)
    If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
    Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
    custom keywords to overwrite_inputs, e.g.
    overwrite_inputs = {“pcm”: {“theory”: “ssvpe”}}
    Refer to the QChem manual for further details on the models available.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **isosvp_dielectric** (*float*) – Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
    (Default: None). If supplied, will set solvent_method to “isosvp” and populate the $svp section
    of the input file with appropriate parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **smd_solvent** (*str*) – Solvent to use for SMD implicit solvation model. (Default: None)
    Examples include “water”, “ethanol”, “methanol”, and “acetonitrile”. Refer to the QChem
    manual for a complete list of solvents available. To define a custom solvent, set this
    argument to “custom” and populate custom_smd with the necessary parameters.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **cmirs_solvent** (*str*) – Solvent to use for the CMIRS implicit solvation model. (Default: None).
    Only 5 solvents are presently available as of Q-Chem 6: “water”, “benzene”, “cyclohexane”,
    “dimethyl sulfoxide”, and “acetonitrile”. Note that selection of a solvent here will also
    populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
    to compute electrostatics.

    **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**



    * **custom_smd** (*str*) – List of parameters to define a custom solvent in SMD. (Default: None)
    Must be given as a string of seven comma separated values in the following order:
    “dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
    electronegative halogenicity”
    Refer to the QChem manual for further details.


    * **max_scf_cycles** (*int*) – Maximum number of SCF iterations. (Default: 100)


    * **geom_opt_max_cycles** (*int*) – Maximum number of geometry optimization iterations. (Default: 200)


    * **geom_opt** (*dict*) – A dict containing parameters for the $geom_opt section of the Q-Chem input
    file, which control the new geometry optimizer available starting in version 5.4.2. The
    new optimizer remains under development but was officially released and became the default
    optimizer in Q-Chem version 6.0.0. Note that for version 5.4.2, the new optimizer must be
    explicitly requested by passing in a dictionary (empty or otherwise) for this input parameter.
    (Default: False)


    * **plot_cubes** (*bool*) – Whether to write CUBE files of the electron density. (Default: False)


    * **overwrite_inputs** (*dict*) – Dictionary of QChem input sections to add or overwrite variables.
    The currently available sections (keys) are rem, pcm,
    solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
    dictionary of key value pairs relevant to that section. For example, to add
    a new variable to the rem section that sets symmetry to false, use

    overwrite_inputs = {“rem”: {“symmetry”: “false”}}

    **Note that if something like basis is added to the rem dict it will overwrite
    the default basis.**

    **Note that supplying a van_der_waals section here will automatically modify
    the PCM “radii” setting to “read”.**

    **Note that all keys must be given as strings, even when they are numbers!**



    * **vdw_mode** (*'atomic'** | **'sequential'*) – Method of specifying custom van der Waals radii. Applies
    only if you are using overwrite_inputs to add a $van_der_waals section to the input.
    In ‘atomic’ mode (default), dict keys represent the atomic number associated with each
    radius (e.g., ‘12’ = carbon). In ‘sequential’ mode, dict keys represent the sequential
    position of a single specific atom in the input structure.


## pymatgen.io.qchem.utils module

Utilities for Qchem io.


### lower_and_check_unique(dict_to_check)
Takes a dictionary and makes all the keys lower case. Also converts all numeric
values (floats, ints) to str and replaces “jobtype” with “job_type” just so that
key specifically can be called elsewhere without ambiguity. Finally, ensures that
multiple identical keys, that differed only due to different capitalizations, are not
present. If there are multiple equivalent keys, an Exception is raised.


* **Parameters**

    **dict_to_check** (*dict*) – The dictionary to check and standardize



* **Returns**

    An identical dictionary but with all keys made

        lower case and no identical keys present.




* **Return type**

    to_return (dict)



### process_parsed_HESS(hess_data)
Takes the information contained in a HESS file and converts it into
the format of the machine-readable 132.0 file which can be printed
out to be read into subsequent optimizations.


### process_parsed_coords(coords)
Takes a set of parsed coordinates, which come as an array of strings,
and returns a numpy array of floats.


### process_parsed_fock_matrix(fock_matrix)
The Fock matrix is parsed as a list, while it should actually be
a square matrix, this function takes the list of finds the right dimensions
in order to reshape the matrix.


### read_matrix_pattern(header_pattern, footer_pattern, elements_pattern, text, postprocess=<class 'str'>)
Parse a matrix to get the quantities in a numpy array.


### read_pattern(text_str, patterns, terminate_on_match=False, postprocess=<class 'str'>)
General pattern reading on an input string.


* **Parameters**


    * **text_str** (*str*) – the input string to search for patterns


    * **patterns** (*dict*) – A dict of patterns, e.g.,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”}.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


Renders accessible:

    Any attribute in patterns. For example,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”} will set the
    value of matches[“energy”] = [[-1234], [-3453], …], to the
    results from regex and postprocess. Note that the returned values
    are lists of lists, because you can grep multiple items on one line.


### read_table_pattern(text_str, header_pattern, row_pattern, footer_pattern, postprocess=<class 'str'>, attribute_name=None, last_one_only=False)
Parse table-like data. A table composes of three parts: header,
main body, footer. All the data matches “row pattern” in the main body
will be returned.


* **Parameters**


    * **text_str** (*str*) – the input string to search for patterns


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



* **Returns**

    List of tables. 1) A table is a list of rows. 2) A row if either a list of
    attribute values in case the capturing group is defined without name in
    row_pattern, or a dict in case that named capturing groups are defined by
    row_pattern.