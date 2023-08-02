---
layout: default
title: pymatgen.io.feff.sets.md
nav_exclude: true
---

# pymatgen.io.feff.sets module

This module defines the FeffInputSet abstract base class and a concrete
implementation for the Materials Project. The basic concept behind an input
set is to specify a scheme to generate a consistent set of Feff inputs from a
structure without further user intervention. This ensures comparability across
runs.


### _class_ pymatgen.io.feff.sets.AbstractFeffInputSet()
Bases: `MSONable`

Abstract base class representing a set of Feff input parameters.
The idea is that using a FeffInputSet, a complete set of input files
(feffPOT, feffXANES, feffEXAFS, ATOMS, feff.inp)set_
can be generated in an automated fashion for any structure.


#### all_input()
Returns all input files as a dict of {filename: feffio object}.


#### _abstract property_ atoms()
Returns Atoms string from a structure that goes in feff.inp file.


* **Returns**

    Atoms object.



#### _abstract_ header()
Returns header to be used in feff.inp file from a pymatgen structure.


#### _abstract property_ potential()
Returns POTENTIAL section used in feff.inp from a structure.


#### _abstract property_ tags()
Returns standard calculation parameters.


#### write_input(output_dir='.', make_dir_if_not_present=True)
Writes a set of FEFF input to a directory.


* **Parameters**


    * **output_dir** – Directory to output the FEFF input files


    * **make_dir_if_not_present** – Set to True if you want the directory (
    and the whole path) to be created if it is not present.



### _class_ pymatgen.io.feff.sets.FEFFDictSet(absorbing_atom: str | int, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), radius: float, config_dict: dict, edge: str = 'K', spectrum: str = 'EXAFS', nkpts=1000, user_tag_settings: dict | None = None, spacegroup_analyzer_settings: dict | None = None)
Bases: `AbstractFeffInputSet`

Standard implementation of FeffInputSet, which can be extended by specific
implementations.


* **Parameters**


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or site index


    * **structure** – Structure or Molecule object. If a Structure, SpaceGroupAnalyzer is used to
    determine symmetrically-equivalent sites. If a Molecule, there is no symmetry
    checking.


    * **radius** (*float*) – cluster radius


    * **config_dict** (*dict*) – control tag settings dict


    * **edge** (*str*) – absorption edge


    * **spectrum** (*str*) – type of spectrum to calculate, available options :
    EXAFS, XANES, DANES, XMCD, ELNES, EXELFS, FPRIME, NRIXS, XES.
    The default is EXAFS.


    * **nkpts** (*int*) – Total number of kpoints in the brillouin zone. Used
    only when feff is run in the reciprocal space mode.


    * **user_tag_settings** (*dict*) – override default tag settings. To delete
    tags, set the key ‘_del’ in the user_tag_settings.
    eg: user_tag_settings={“_del”: [“COREHOLE”, “EXCHANGE”]}
    To specify a net charge on the structure, pass an “IONS” tag containing a list

    > of tuples where the first element is the unique potential value (ipot value)
    > and the second element is the charge to be applied to atoms associated
    > with that potential, e.g. {“IONS”: [(0, 0.1), (1, 0.1), (2, 0.1)]}
    > will result in.

    > ION 0 0.1
    > ION 1 0.1
    > ION 2 0.1

    > being written to the input file.



    * **spacegroup_analyzer_settings** (*dict*) – parameters passed to SpacegroupAnalyzer.
    E.g., {“symprec”: 0.01, “angle_tolerance”: 4}



#### _property_ atoms(_: [Atoms](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms_ )
absorber + the rest.


* **Returns**

    Atoms



#### _static_ from_directory(input_dir)
Read in a set of FEFF input files from a directory, which is
useful when existing FEFF input needs some adjustment.


#### header(source: str = '', comment: str = '')
Creates header string from structure object.


* **Parameters**


    * **source** – Source identifier used to create structure, can be defined
    however user wants to organize structures, calculations, etc.
    example would be Materials Project material ID number.


    * **comment** – comment to include in header



* **Returns**

    Header



#### _property_ potential(_: [Potential](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Potential_ )
FEFF potential.


* **Returns**

    Potential



#### _property_ tags(_: [Tags](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags_ )
FEFF job parameters.


* **Returns**

    Tags



### _class_ pymatgen.io.feff.sets.MPEELSDictSet(absorbing_atom, structure, edge, spectrum, radius, beam_energy, beam_direction, collection_angle, convergence_angle, config_dict, user_eels_settings=None, nkpts: int = 1000, user_tag_settings: dict | None = None, \*\*kwargs)
Bases: `FEFFDictSet`

FeffDictSet for ELNES spectroscopy.


* **Parameters**


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or site index


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure


    * **edge** (*str*) – absorption edge


    * **spectrum** (*str*) – ELNES or EXELFS


    * **radius** (*float*) – cluster radius in Angstroms.


    * **beam_energy** (*float*) – Incident beam energy in keV


    * **beam_direction** (*list*) – Incident beam direction. If None, the
    cross section will be averaged.


    * **collection_angle** (*float*) – Detector collection angle in mrad.


    * **convergence_angle** (*float*) – Beam convergence angle in mrad.


    * **user_eels_settings** (*dict*) – override default EELS config.
    See MPELNESSet.yaml for supported keys.


    * **nkpts** (*int*) – Total number of kpoints in the brillouin zone. Used
    only when feff is run in the reciprocal space mode.


    * **user_tag_settings** (*dict*) – override default tag settings


    * **\*\*kwargs** – Passthrough to FEFFDictSet.



### _class_ pymatgen.io.feff.sets.MPELNESSet(absorbing_atom, structure, edge: str = 'K', radius: float = 10.0, beam_energy: float = 100, beam_direction=None, collection_angle: float = 1, convergence_angle: float = 1, user_eels_settings=None, nkpts: int = 1000, user_tag_settings: dict | None = None, \*\*kwargs)
Bases: `MPEELSDictSet`

FeffDictSet for ELNES spectroscopy.


* **Parameters**


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or site index


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure


    * **edge** (*str*) – absorption edge


    * **radius** (*float*) – cluster radius in Angstroms.


    * **beam_energy** (*float*) – Incident beam energy in keV


    * **beam_direction** (*list*) – Incident beam direction. If None, the
    cross section will be averaged.


    * **collection_angle** (*float*) – Detector collection angle in mrad.


    * **convergence_angle** (*float*) – Beam convergence angle in mrad.


    * **user_eels_settings** (*dict*) – override default EELS config.
    See MPELNESSet.yaml for supported keys.


    * **nkpts** (*int*) – Total number of kpoints in the brillouin zone. Used
    only when feff is run in the reciprocal space mode.


    * **user_tag_settings** (*dict*) – override default tag settings


    * **\*\*kwargs** – Passthrough to FEFFDictSet.



#### CONFIG(_ = {'CONTROL': '1 1 1 1 1 1', 'COREHOLE': 'FSR', 'EDGE': 'K', 'ELNES': {'ANGLES': '1 1', 'BEAM_DIRECTION': '0 1 0', 'BEAM_ENERGY': '100 0 1 1', 'ENERGY': '4 0.04 0.1', 'MESH': '50 1', 'POSITION': '0.0 0.0'}, 'EXCHANGE': '0 0.0 0.0 2', 'FMS': '7.5 0', 'LDOS': '-20.0 20.0 0.1', 'PRINT': '1 0 0 0 0 0', 'S02': 0.0, 'SCF': '6.0 0 30 0.2 1'_ )

### _class_ pymatgen.io.feff.sets.MPEXAFSSet(absorbing_atom, structure, edge: str = 'K', radius: float = 10.0, nkpts: int = 1000, user_tag_settings: dict | None = None, \*\*kwargs)
Bases: `FEFFDictSet`

FeffDictSet for EXAFS spectroscopy.


* **Parameters**


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or site index


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure


    * **edge** (*str*) – absorption edge


    * **radius** (*float*) – cluster radius in Angstroms.


    * **nkpts** (*int*) – Total number of kpoints in the brillouin zone. Used
    only when feff is run in the reciprocal space mode.


    * **user_tag_settings** (*dict*) – override default tag settings


    * **\*\*kwargs** – Passthrough to FEFFDictSet.



#### CONFIG(_ = {'CONTROL': '1 1 1 1 1 1', 'COREHOLE': 'FSR', 'EDGE': 'K', 'EXAFS': 20, 'PRINT': '1 0 0 0 0 0', 'RPATH': 10, 'S02': 0.0, 'SCF': '4.5 0 30 .2 1'_ )

### _class_ pymatgen.io.feff.sets.MPEXELFSSet(absorbing_atom, structure, edge='K', radius: float = 10.0, beam_energy: float = 100, beam_direction=None, collection_angle: float = 1, convergence_angle: float = 1, user_eels_settings=None, nkpts: int = 1000, user_tag_settings: dict | None = None, \*\*kwargs)
Bases: `MPEELSDictSet`

FeffDictSet for EXELFS spectroscopy.


* **Parameters**


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or site index


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure


    * **edge** (*str*) – absorption edge


    * **radius** (*float*) – cluster radius in Angstroms.


    * **beam_energy** (*float*) – Incident beam energy in keV


    * **beam_direction** (*list*) – Incident beam direction. If None, the
    cross section will be averaged.


    * **collection_angle** (*float*) – Detector collection angle in mrad.


    * **convergence_angle** (*float*) – Beam convergence angle in mrad.


    * **user_eels_settings** (*dict*) – override default EELS config.
    See MPEXELFSSet.yaml for supported keys.


    * **nkpts** (*int*) – Total number of kpoints in the brillouin zone. Used
    only when feff is run in the reciprocal space mode.


    * **user_tag_settings** (*dict*) – override default tag settings


    * **\*\*kwargs** – Passthrough to FEFFDictSet.



#### CONFIG(_ = {'CONTROL': '1 1 1 1 1 1', 'COREHOLE': 'FSR', 'EDGE': 'K', 'EXCHANGE': '0 0.0 0.0 2', 'EXELFS': {'ANGLES': '1 1', 'BEAM_DIRECTION': '0 1 0', 'BEAM_ENERGY': '100 0 1 1', 'ENERGY': 20, 'MESH': '50 1', 'POSITION': '0.0 0.0'}, 'PRINT': '1 0 0 0 0 0', 'RPATH': 10, 'S02': 0.0, 'SCF': '5.0 0 30 0.2 1'_ )

### _class_ pymatgen.io.feff.sets.MPXANESSet(absorbing_atom, structure, edge: str = 'K', radius: float = 10.0, nkpts: int = 1000, user_tag_settings: dict | None = None, \*\*kwargs)
Bases: `FEFFDictSet`

FeffDictSet for XANES spectroscopy.


* **Parameters**


    * **absorbing_atom** (*str/int*) – absorbing atom symbol or site index


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input


    * **edge** (*str*) – absorption edge


    * **radius** (*float*) – cluster radius in Angstroms.


    * **nkpts** (*int*) – Total number of kpoints in the brillouin zone. Used
    only when feff is run in the reciprocal space mode.


    * **user_tag_settings** (*dict*) – override default tag settings


    * **\*\*kwargs** – Passthrough to FEFFDictSet.



#### CONFIG(_ = {'CONTROL': '1 1 1 1 1 1', 'COREHOLE': 'FSR', 'EDGE': 'K', 'EXCHANGE': '0 0.0 0.0 2', 'FMS': '7.5 0', 'LDOS': '-30. 15. .1', 'PRINT': '1 0 0 0 0 0', 'S02': 0.0, 'SCF': '4.5 0 30 .2 1', 'XANES': '3.7 .04 .1'_ )