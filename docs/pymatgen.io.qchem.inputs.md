---
layout: default
title: pymatgen.io.qchem.inputs.md
nav_exclude: true
---

# pymatgen.io.qchem.inputs module

Classes for reading/manipulating/writing QChem input files.


### _class_ pymatgen.io.qchem.inputs.QCInput(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule) | list[[Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)] | Literal['read'], rem: dict, opt: dict[str, list] | None = None, pcm: dict | None = None, solvent: dict | None = None, smx: dict | None = None, scan: dict[str, list] | None = None, van_der_waals: dict[str, float] | None = None, vdw_mode: str = 'atomic', plots: dict | None = None, nbo: dict | None = None, geom_opt: dict | None = None, cdft: list[list[dict]] | None = None, almo_coupling: list[list[tuple[int, int]]] | None = None, svp: dict | None = None, pcm_nonels: dict | None = None)
Bases: [`InputFile`](pymatgen.io.core.md#pymatgen.io.core.InputFile)

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



#### get_string()
Return a string representation of an entire input file.


#### _static_ molecule_template(molecule: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule) | list[[Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)] | Literal['read'])

* **Parameters**

    **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)*, **list** of **Molecules**, or **"read"*) –



* **Returns**

    (str) Molecule template.



#### _static_ multi_job_string(job_list: list[pymatgen.io.qchem.inputs.QCInput])

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

    (list of lists of dicts) cdft parameters



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



#### _static_ write_multi_job_file(job_list: list[pymatgen.io.qchem.inputs.QCInput], filename: str)
Write a multijob file.


* **Parameters**


    * **(****)** (*filename*) – List of jobs.


    * **(****)** – Filename