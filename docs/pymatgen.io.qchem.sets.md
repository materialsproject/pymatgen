---
layout: default
title: pymatgen.io.qchem.sets.md
nav_exclude: true
---

# pymatgen.io.qchem.sets module

Input sets for Qchem.


### _class_ pymatgen.io.qchem.sets.ForceSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-tzvpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', cdft_constraints: list[list[dict]] | None = None, overwrite_inputs: dict | None = None)
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




### _class_ pymatgen.io.qchem.sets.FreqSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', cdft_constraints: list[list[dict]] | None = None, overwrite_inputs: dict | None = None)
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




### _class_ pymatgen.io.qchem.sets.OptSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, opt_variables: dict[str, list] | None = None, geom_opt_max_cycles: int = 200, geom_opt: dict | None = None, cdft_constraints: list[list[dict]] | None = None, overwrite_inputs: dict | None = None)
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



### _class_ pymatgen.io.qchem.sets.PESScanSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, opt_variables: dict[str, list] | None = None, scan_variables: dict[str, list] | None = None, overwrite_inputs: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic')
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



### _class_ pymatgen.io.qchem.sets.QChemDictSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), job_type: str, basis_set: str, scf_algorithm: str, qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, opt_variables: dict[str, list] | None = None, scan_variables: dict[str, list] | None = None, max_scf_cycles: int = 100, geom_opt_max_cycles: int = 200, plot_cubes: bool = False, nbo_params: dict | None = None, geom_opt: dict | None = None, cdft_constraints: list[list[dict]] | None = None, almo_coupling_states: list[list[tuple[int, int]]] | None = None, overwrite_inputs: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', extra_scf_print: bool = False)
Bases: [`QCInput`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput)

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



### _class_ pymatgen.io.qchem.sets.SinglePointSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-tzvpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, vdw_mode: Literal['atomic', 'sequential'] = 'atomic', cdft_constraints: list[list[dict]] | None = None, almo_coupling_states: list[list[tuple[int, int]]] | None = None, extra_scf_print: bool = False, overwrite_inputs: dict | None = None)
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



### _class_ pymatgen.io.qchem.sets.TransitionStateSet(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), basis_set: str = 'def2-svpd', scf_algorithm: str = 'diis', qchem_version: int = 5, dft_rung: int = 4, pcm_dielectric: float | None = None, isosvp_dielectric: float | None = None, smd_solvent: str | None = None, cmirs_solvent: Literal['water', 'acetonitrile', 'dimethyl sulfoxide', 'cyclohexane', 'benzene'] | None = None, custom_smd: str | None = None, max_scf_cycles: int = 100, plot_cubes: bool = False, nbo_params: dict | None = None, opt_variables: dict[str, list] | None = None, geom_opt_max_cycles: int = 200, geom_opt: dict | None = None, overwrite_inputs: dict | None = None, vdw_mode='atomic')
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