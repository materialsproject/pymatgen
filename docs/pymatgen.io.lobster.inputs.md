---
layout: default
title: pymatgen.io.lobster.inputs.md
nav_exclude: true
---

# pymatgen.io.lobster.inputs module

Module for reading Lobster input files. For more information
on LOBSTER see www.cohp.de.
If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
“Automated Bonding Analysis with Crystal Orbital Hamilton Populations”,
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.


### _class_ pymatgen.io.lobster.inputs.Lobsterin(settingsdict: dict)
Bases: `dict`, `MSONable`

This class can handle and generate lobsterin files
Furthermore, it can also modify INCAR files for lobster, generate KPOINT files for fatband calculations in Lobster,
and generate the standard primitive cells in a POSCAR file that are needed for the fatband calculations.
There are also several standard lobsterin files that can be easily generated.


* **Parameters**

    **settingsdict** – dict to initialize Lobsterin.



#### AVAILABLEKEYWORDS(_ = ('COHPstartEnergy', 'COHPendEnergy', 'gaussianSmearingWidth', 'useDecimalPlaces', 'COHPSteps', 'basisSet', 'cohpGenerator', 'realspaceHamiltonian', 'realspaceOverlap', 'printPAWRealSpaceWavefunction', 'printLCAORealSpaceWavefunction', 'kSpaceCOHP', 'EwaldSum', 'saveProjectionToFile', 'skipdos', 'skipcohp', 'skipcoop', 'skipcobi', 'skipMadelungEnergy', 'loadProjectionFromFile', 'forceEnergyRange', 'DensityOfEnergy', 'BWDF', 'BWDFCOHP', 'skipPopulationAnalysis', 'skipGrossPopulation', 'userecommendedbasisfunctions', 'skipProjection', 'writeBasisFunctions', 'writeMatricesToFile', 'noFFTforVisualization', 'RMSp', 'onlyReadVasprun.xml', 'noMemoryMappedFiles', 'skipPAWOrthonormalityTest', 'doNotIgnoreExcessiveBands', 'doNotUseAbsoluteSpilling', 'skipReOrthonormalization', 'forceV1HMatrix', 'useOriginalTetrahedronMethod', 'forceEnergyRange', 'bandwiseSpilling', 'kpointwiseSpilling', 'LSODOS', 'basisfunctions', 'cohpbetween', 'createFatband'_ )

#### BOOLEAN_KEYWORDS(_ = ('saveProjectionToFile', 'skipdos', 'skipcohp', 'skipcoop', 'skipcobi', 'skipMadelungEnergy', 'loadProjectionFromFile', 'forceEnergyRange', 'DensityOfEnergy', 'BWDF', 'BWDFCOHP', 'skipPopulationAnalysis', 'skipGrossPopulation', 'userecommendedbasisfunctions', 'skipProjection', 'writeBasisFunctions', 'writeMatricesToFile', 'noFFTforVisualization', 'RMSp', 'onlyReadVasprun.xml', 'noMemoryMappedFiles', 'skipPAWOrthonormalityTest', 'doNotIgnoreExcessiveBands', 'doNotUseAbsoluteSpilling', 'skipReOrthonormalization', 'forceV1HMatrix', 'useOriginalTetrahedronMethod', 'forceEnergyRange', 'bandwiseSpilling', 'kpointwiseSpilling', 'LSODOS'_ )

#### FLOAT_KEYWORDS(_ = ('COHPstartEnergy', 'COHPendEnergy', 'gaussianSmearingWidth', 'useDecimalPlaces', 'COHPSteps'_ )

#### LISTKEYWORDS(_ = ('basisfunctions', 'cohpbetween', 'createFatband'_ )

#### STRING_KEYWORDS(_ = ('basisSet', 'cohpGenerator', 'realspaceHamiltonian', 'realspaceOverlap', 'printPAWRealSpaceWavefunction', 'printLCAORealSpaceWavefunction', 'kSpaceCOHP', 'EwaldSum'_ )

#### as_dict()

* **Returns**

    MSONable dict



#### diff(other)
Diff function for lobsterin. Compares two lobsterin and indicates which parameters are the same.
Similar to the diff in INCAR.


* **Parameters**

    **other** (*Lobsterin*) – Lobsterin object to compare to



* **Returns**

    dict with differences and similarities



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    Lobsterin



#### _classmethod_ from_file(lobsterin: str)

* **Parameters**

    **lobsterin** (*str*) – path to lobsterin.



* **Returns**

    Lobsterin object



#### _static_ get_all_possible_basis_functions(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), potcar_symbols: list, address_basis_file_min: str | None = None, address_basis_file_max: str | None = None)

* **Parameters**


    * **structure** – Structure object


    * **potcar_symbols** – list of the potcar symbols


    * **address_basis_file_min** – path to file with the minimum required basis by the POTCAR


    * **address_basis_file_max** – path to file with the largest possible basis of the POTCAR.


Returns: List of dictionaries that can be used to create new Lobsterin objects in
standard_calculations_from_vasp_files as dict_for_basis


#### _static_ get_basis(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), potcar_symbols: list, address_basis_file: str | None = None)
Will get the basis from given potcar_symbols (e.g., [“Fe_pv”,”Si”]
#include this in lobsterin class.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object


    * **potcar_symbols** – list of potcar symbols



* **Returns**

    returns basis



#### _classmethod_ standard_calculations_from_vasp_files(POSCAR_input: str = 'POSCAR', INCAR_input: str = 'INCAR', POTCAR_input: str | None = None, Vasprun_output: str = 'vasprun.xml', dict_for_basis: dict | None = None, option: str = 'standard')
Will generate Lobsterin with standard settings.


* **Parameters**


    * **POSCAR_input** (*str*) – path to POSCAR


    * **INCAR_input** (*str*) – path to INCAR


    * **POTCAR_input** (*str*) – path to POTCAR


    * **dict_for_basis** (*dict*) – can be provided: it should look the following:
    dict_for_basis={“Fe”:’3p 3d 4s 4f’, “C”: ‘2s 2p’} and will overwrite all settings from POTCAR_input


    * **option** (*str*) – ‘standard’ will start a normal lobster run where COHPs, COOPs, DOS, CHARGE etc. will be
    calculated
    ‘standard_with_energy_range_from_vasprun’ will start a normal lobster run for entire energy range
    of VASP static run. vasprun.xml file needs to be in current directory.
    ‘standard_from_projection’ will start a normal lobster run from a projection
    ‘standard_with_fatband’ will do a fatband calculation, run over all orbitals
    ‘onlyprojection’ will only do a projection
    ‘onlydos’ will only calculate a projected dos
    ‘onlycohp’ will only calculate cohp
    ‘onlycoop’ will only calculate coop
    ‘onlycohpcoop’ will only calculate cohp and coop



* **Returns**

    Lobsterin Object with standard settings



#### write_INCAR(incar_input: str = 'INCAR', incar_output: str = 'INCAR.lobster', poscar_input: str = 'POSCAR', isym: int = -1, further_settings: dict | None = None)
Will only make the run static, insert nbands, make ISYM=-1, set LWAVE=True and write a new INCAR.
You have to check for the rest.


* **Parameters**


    * **incar_input** (*str*) – path to input INCAR


    * **incar_output** (*str*) – path to output INCAR


    * **poscar_input** (*str*) – path to input POSCAR


    * **isym** (*int*) – isym equal to -1 or 0 are possible. Current Lobster version only allow -1.


    * **further_settings** (*dict*) – A dict can be used to include further settings, e.g. {“ISMEAR”:-5}



#### _static_ write_KPOINTS(POSCAR_input: str = 'POSCAR', KPOINTS_output='KPOINTS.lobster', reciprocal_density: int = 100, isym: int = -1, from_grid: bool = False, input_grid: Sequence[int] = (5, 5, 5), line_mode: bool = True, kpoints_line_density: int = 20, symprec: float = 0.01)
Writes a KPOINT file for lobster (only ISYM=-1 and ISYM=0 are possible), grids are gamma centered.


* **Parameters**


    * **POSCAR_input** (*str*) – path to POSCAR


    * **KPOINTS_output** (*str*) – path to output KPOINTS


    * **reciprocal_density** (*int*) – Grid density


    * **isym** (*int*) – either -1 or 0. Current Lobster versions only allow -1.


    * **from_grid** (*bool*) – If True KPOINTS will be generated with the help of a grid given in input_grid. Otherwise,
    they will be generated from the reciprocal_density


    * **input_grid** (*list*) – grid to generate the KPOINTS file


    * **line_mode** (*bool*) – If True, band structure will be generated


    * **kpoints_line_density** (*int*) – density of the lines in the band structure


    * **symprec** (*float*) – precision to determine symmetry



#### _static_ write_POSCAR_with_standard_primitive(POSCAR_input='POSCAR', POSCAR_output='POSCAR.lobster', symprec: float = 0.01)
Writes a POSCAR with the standard primitive cell. This is needed to arrive at the correct kpath.


* **Parameters**


    * **POSCAR_input** (*str*) – filename of input POSCAR


    * **POSCAR_output** (*str*) – filename of output POSCAR


    * **symprec** (*float*) – precision to find symmetry



#### write_lobsterin(path='lobsterin', overwritedict=None)
Writes a lobsterin file.


* **Parameters**


    * **path** (*str*) – filename of the lobsterin file that will be written


    * **overwritedict** (*dict*) – dict that can be used to overwrite lobsterin, e.g. {“skipdos”: True}



### pymatgen.io.lobster.inputs.get_all_possible_basis_combinations(min_basis: list, max_basis: list)

* **Parameters**


    * **min_basis** – list of basis entries: e.g., [‘Si 3p 3s ‘]


    * **max_basis** – list of basis entries: e.g., [‘Si 3p 3s ‘].


Returns: all possible combinations of basis functions, e.g. [[‘Si 3p 3s’]]