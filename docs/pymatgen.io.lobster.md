---
layout: default
title: pymatgen.io.lobster.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.io.lobster package

This package implements modules for input and output to and from Lobster. It
imports the key classes form both lobster.inputs and lobster_outputs to allow most
classes to be simply called as pymatgen.io.lobster.Lobsterin for example, to retain
backwards compatibility.


## pymatgen.io.lobster.inputs module

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



#### _static_ get_all_possible_basis_functions(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), potcar_symbols: list, address_basis_file_min: str | None = None, address_basis_file_max: str | None = None)

* **Parameters**


    * **structure** – Structure object


    * **potcar_symbols** – list of the potcar symbols


    * **address_basis_file_min** – path to file with the minimum required basis by the POTCAR


    * **address_basis_file_max** – path to file with the largest possible basis of the POTCAR.


Returns: List of dictionaries that can be used to create new Lobsterin objects in
standard_calculations_from_vasp_files as dict_for_basis


#### _static_ get_basis(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), potcar_symbols: list, address_basis_file: str | None = None)
Will get the basis from given potcar_symbols (e.g., [“Fe_pv”,”Si”]
#include this in lobsterin class.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object


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

## pymatgen.io.lobster.lobsterenv module

This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures based on
bonding analysis with Lobster.
If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
“Automated Bonding Analysis with Crystal Orbital Hamilton Populations”,
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.


### _class_ pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo(total_icohp: float, list_icohps: list[float], n_bonds: int, labels: list[str], atoms: list[list[str]], central_isites: list[int] | None)
Bases: `NamedTuple`

Tuple to represent information on relevant bonds
:param total_icohp: sum of icohp values of neighbors to the selected sites [given by the id in structure]
:type total_icohp: float
:param list_icohps: list of summed icohp values for all identified interactions with neighbors
:type list_icohps: list
:param n_bonds: number of identified bonds to the selected sites
:type n_bonds: int
:param labels: labels (from ICOHPLIST) for all identified bonds
:type labels: list(str)
:param atoms: list of list describing the species present in the identified interactions

> (names from ICOHPLIST), e.g., [‘Ag3’, ‘O5’]


* **Parameters**

    **central_isites** (*list**(**int**)*) – list of the central isite for each identified interaction.


Create new instance of ICOHPNeighborsInfo(total_icohp, list_icohps, n_bonds, labels, atoms, central_isites)


#### atoms(_: list[list[str]_ )
Alias for field number 4


#### central_isites(_: list[int] | Non_ )
Alias for field number 5


#### labels(_: list[str_ )
Alias for field number 3


#### list_icohps(_: list[float_ )
Alias for field number 1


#### n_bonds(_: in_ )
Alias for field number 2


#### total_icohp(_: floa_ )
Alias for field number 0


### _class_ pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments(strategy, coordination_environments=None, all_nbs_sites=None, neighbors_sets=None, structure=None, valences=None, valences_origin=None)
Bases: [`LightStructureEnvironments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments)

Class to store LightStructureEnvironments based on Lobster outputs.

Constructor for the LightStructureEnvironments object.


* **Parameters**


    * **strategy** – ChemEnv strategy used to get the environments.


    * **coordination_environments** – The coordination environments identified.


    * **all_nbs_sites** – All the possible neighbors for each site in the structure.


    * **neighbors_sets** – The neighbors sets of each site in the structure.


    * **structure** – The structure.


    * **valences** – The valences used to get the environments (if needed).


    * **valences_origin** – How the valences were obtained (e.g. from the Bond-valence analysis or from the original
    structure).



#### as_dict()
Bson-serializable dict representation of the LightStructureEnvironments object.
:return: Bson-serializable dict representation of the LightStructureEnvironments object.


#### _classmethod_ from_Lobster(list_ce_symbol, list_csm, list_permutation, list_neighsite, list_neighisite, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), valences=None)
Will set up a LightStructureEnvironments from Lobster.


* **Parameters**


    * **structure** – Structure object


    * **list_ce_symbol** – list of symbols for coordination environments


    * **list_csm** – list of continuous symmetry measures


    * **list_permutation** – list of permutations


    * **list_neighsite** – list of neighboring sites


    * **list_neighisite** – list of neighboring isites (number of a site)


    * **valences** – list of valences


Returns: LobsterLightStructureEnvironments


#### _property_ uniquely_determines_coordination_environments()
True if the coordination environments are uniquely determined.


### _class_ pymatgen.io.lobster.lobsterenv.LobsterNeighbors(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), filename_ICOHP: str = 'ICOHPLIST.lobster', are_coops: bool = False, are_cobis: bool = False, valences: list[int | float] | None = None, limits: tuple[float, float] | None = None, additional_condition: int = 0, only_bonds_to: list[str] | None = None, perc_strength_ICOHP: float = 0.15, noise_cutoff: float = 0.1, valences_from_charges: bool = False, filename_CHARGE: str | None = None, which_charge: str = 'Mulliken', adapt_extremum_to_add_cond: bool = False, add_additional_data_sg: bool = False, filename_blist_sg1: str | None = None, filename_blist_sg2: str | None = None, id_blist_sg1: str = 'ICOOP', id_blist_sg2: str = 'ICOBI')
Bases: [`NearNeighbors`](pymatgen.analysis.md#pymatgen.analysis.local_env.NearNeighbors)

This class combines capabilities from LocalEnv and ChemEnv to determine coordination environments based on
bonding analysis.


* **Parameters**


    * **filename_ICOHP** – (str) Path to ICOHPLIST.lobster or ICOOPLIST.lobster or ICOBILIST.lobster


    * **structure** – (Structure) typically constructed by Structure.from_file(“POSCAR”)


    * **are_coops** – (bool) if True, the file is a ICOOPLIST.lobster and not a ICOHPLIST.lobster; only tested for
    ICOHPLIST.lobster so far


    * **are_cobis** – (bool) if True, the file is a ICOBILIST.lobster and not a ICOHPLIST.lobster


    * **valences** (*Mulliken**) **instead of*) – (list[int | float]): gives valence/charge for each element


    * **limits** (*tuple**[**float**, **float**] **| **None*) – limit to decide which ICOHPs (ICOOP or ICOBI) should be considered


    * **additional_condition** (*int*) – Additional condition that decides which kind of bonds will be considered
    NO_ADDITIONAL_CONDITION = 0
    ONLY_ANION_CATION_BONDS = 1
    NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
    ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
    ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
    DO_NOT_CONSIDER_ANION_CATION_BONDS=5
    ONLY_CATION_CATION_BONDS=6


    * **only_bonds_to** – (list[str]) will only consider bonds to certain elements (e.g. [“O”] for oxygen)


    * **perc_strength_ICOHP** – if no limits are given, this will decide which icohps will still be considered (


    * **to** (*relative*) –


    * **ICOHP** (*the strongest*) –


    * **noise_cutoff** – if provided hardcodes the lower limit of icohps considered


    * **valences_from_charges** – if True and path to CHARGE.lobster is provided, will use Lobster charges (


    * **valences** –


    * **filename_CHARGE** – (str) Path to Charge.lobster


    * **which_charge** – (str) “Mulliken” or “Loewdin”


    * **adapt_extremum_to_add_cond** – (bool) will adapt the limits to only focus on the bonds determined by the


    * **condition** (*additional*) –


    * **add_additional_data_sg** – (bool) will add the information from filename_add_bondinglist_sg1,


    * **filename_blist_sg1** – (str) Path to additional ICOOP, ICOBI data for structure graphs


    * **filename_blist_sg2** – (str) Path to dditional ICOOP, ICOBI data for structure graphs


    * **id_blist_sg1** – (str) Identity of data in filename_blist_sg1,
    e.g., “icoop” or “icobi”


    * **id_blist_sg2** – (str) Identity of data in filename_blist_sg2,
    e.g., “icoop” or “icobi”.



#### _property_ anion_types()
Return the types of anions present in crystal structure as a set
Returns: set of Element describing anions in the crystal structure.


#### get_anion_types(\*\*kwargs)

#### get_info_cohps_to_neighbors(path_to_COHPCAR='COHPCAR.lobster', isites=None, only_bonds_to=None, onlycation_isites=True, per_bond=True, summed_spin_channels=False)
Return info about the cohps (coops or cobis) as a summed cohp object and a label

    from all sites mentioned in isites with neighbors.


* **Parameters**


    * **path_to_COHPCAR** – str, path to COHPCAR or COOPCAR or COBICAR


    * **isites** – list of int that indicate the number of the site


    * **only_bonds_to** – list of str, e.g. [“O”] to only show cohps of anything to oxygen


    * **onlycation_isites** – if isites=None, only cation sites will be returned


    * **per_bond** – will normalize per bond


    * **summed_spin_channels** – will sum all spin channels


Returns: label for cohp (str), CompleteCohp object which describes all cohps (coops or cobis) of the sites
as given by isites and the other parameters


#### get_info_icohps_between_neighbors(isites=None, onlycation_isites=True)
Return infos about interactions between neighbors of a certain atom.


* **Parameters**


    * **isites** – list of site ids, if isite==None, all isites will be used


    * **onlycation_isites** – will only use cations, if isite==None


Returns: ICOHPNeighborsInfo


#### get_info_icohps_to_neighbors(isites=None, onlycation_isites=True)
This method returns information on the icohps of neighbors for certain sites as identified by their site id.
This is useful for plotting the relevant cohps of a site in the structure.
(could be ICOOPLIST.lobster or ICOHPLIST.lobster or ICOBILIST.lobster)


* **Parameters**


    * **isites** – list of site ids. If isite==None, all isites will be used to add the icohps of the neighbors


    * **onlycation_isites** – if True and if isite==None, it will only analyse the sites of the cations


Returns: ICOHPNeighborsInfo


#### get_light_structure_environment(only_cation_environments=False, only_indices=None)
Return a LobsterLightStructureEnvironments object
if the structure only contains coordination environments smaller 13.


* **Parameters**


    * **only_cation_environments** – only data for cations will be returned


    * **only_indices** – will only evaluate the list of isites in this list


Returns: LobsterLightStructureEnvironments Object


#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n, use_weights=False)
Get coordination number, CN, of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).
    True is not implemented for LobsterNeighbors



* **Returns**

    coordination number.



* **Return type**

    cn (integer or float)



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### plot_cohps_of_neighbors(path_to_COHPCAR='COHPCAR.lobster', isites=None, onlycation_isites=True, only_bonds_to=None, per_bond=False, summed_spin_channels=False, xlim=None, ylim=(-10, 6), integrated=False)
Will plot summed cohps or cobis or coops
(please be careful in the spin polarized case (plots might overlap (exactly!)).


* **Parameters**


    * **isites** – list of site ids, if isite==[], all isites will be used to add the icohps of the neighbors


    * **onlycation_isites** – bool, will only use cations, if isite==[]


    * **only_bonds_to** – list of str, only anions in this list will be considered


    * **per_bond** – bool, will lead to a normalization of the plotted COHP per number of bond if True,


    * **sum** (*otherwise the*) –


    * **plotted** (*will be*) –


    * **xlim** – list of float, limits of x values


    * **ylim** – list of float, limits of y values


    * **integrated** – bool, if true will show integrated cohp instead of cohp



* **Returns**

    plt of the cohps or coops or cobis



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property


## pymatgen.io.lobster.outputs module

Module for reading Lobster output files. For more information
on LOBSTER see www.cohp.de.
If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
“Automated Bonding Analysis with Crystal Orbital Hamilton Populations”,
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.


### _class_ pymatgen.io.lobster.outputs.Bandoverlaps(filename: str = 'bandOverlaps.lobster')
Bases: `object`

Class to read in bandOverlaps.lobster files. These files are not created during every Lobster run.
.. attribute: bandoverlapsdict is a dict of the following form:

> {spin:{“kpoint as string”: {“maxDeviation”: float that describes the max deviation, “matrix”: 2D
> array of the size number of bands times number of bands including the overlap matrices with } }}.

<!-- attribute: maxDeviation is a list of floats describing the maximal Deviation for each problematic kpoint -->

* **Parameters**

    **filename** – filename of the “bandOverlaps.lobster” file.



#### has_good_quality_check_occupied_bands(number_occ_bands_spin_up: int, number_occ_bands_spin_down: int | None = None, spin_polarized: bool = False, limit_deviation: float = 0.1)
Will check if the deviation from the ideal bandoverlap of all occupied bands is smaller or equal to
limit_deviation.

Args:
number_occ_bands_spin_up (int): number of occupied bands of spin up
number_occ_bands_spin_down (int): number of occupied bands of spin down
spin_polarized (bool):  If True, then it was a spin polarized calculation
limit_deviation (float): limit of the maxDeviation
:returns: Boolean that will give you information about the quality of the projection


#### has_good_quality_maxDeviation(limit_maxDeviation: float = 0.1)
Will check if the maxDeviation from the ideal bandoverlap is smaller or equal to limit_maxDeviation
:param limit_maxDeviation: limit of the maxDeviation


* **Returns**

    Boolean that will give you information about the quality of the projection.



### _class_ pymatgen.io.lobster.outputs.Charge(filename: str = 'CHARGE.lobster')
Bases: `object`

Class to read CHARGE files generated by LOBSTER.

<!-- attribute: atomlist
List of atoms in CHARGE.lobster -->
<!-- attribute: types
List of types of atoms in CHARGE.lobster -->
<!-- attribute: Mulliken
List of Mulliken charges of atoms in CHARGE.lobster -->
<!-- attribute: Loewdin
List of Loewdin charges of atoms in CHARGE.Loewdin -->
<!-- attribute: num_atoms
Number of atoms in CHARGE.lobster -->

* **Parameters**

    **filename** – filename for the CHARGE file, typically “CHARGE.lobster”.



#### get_structure_with_charges(structure_filename)
Get a Structure with Mulliken and Loewdin charges as site properties
:param structure_filename: filename of POSCAR


* **Returns**

    Structure Object with Mulliken and Loewdin charges as site properties.



### _class_ pymatgen.io.lobster.outputs.Cohpcar(are_coops: bool = False, are_cobis: bool = False, filename: str | None = None)
Bases: `object`

Class to read COHPCAR/COOPCAR files generated by LOBSTER.

<!-- attribute: cohp_data

Dict that contains the COHP data of the form:
  {bond: {"COHP": {Spin.up: cohps, Spin.down:cohps},
          "ICOHP": {Spin.up: icohps, Spin.down: icohps},
          "length": bond length,
          "sites": sites corresponding to the bond}
Also contains an entry for the average, which does not have
a "length" key. -->
<!-- attribute: efermi

The Fermi energy in eV. -->
<!-- attribute: energies

Sequence of energies in eV. Note that LOBSTER shifts the energies
so that the Fermi energy is at zero. -->
<!-- attribute: is_spin_polarized

Boolean to indicate if the calculation is spin polarized. -->
<!-- attribute: orb_cohp

orb_cohp[label] = {bond_data["orb_label"]: {"COHP": {Spin.up: cohps, Spin.down:cohps},
                                             "ICOHP": {Spin.up: icohps, Spin.down: icohps},
                                             "orbitals": orbitals,
                                             "length": bond lengths,
                                             "sites": sites corresponding to the bond}} -->

* **Parameters**


    * **are_coops** – Determines if the file is a list of COHPs or COOPs.
    Default is False for COHPs.


    * **are_cobis** – Determines if the file is a list of COHPs or COOPs.
    Default is False for COHPs.


    * **filename** – Name of the COHPCAR file. If it is None, the default
    file name will be chosen, depending on the value of are_coops.



### _class_ pymatgen.io.lobster.outputs.Doscar(doscar: str = 'DOSCAR.lobster', structure_file: str | None = 'POSCAR', structure: [IStructure](pymatgen.core.md#pymatgen.core.structure.IStructure) | [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None)
Bases: `object`

Class to deal with Lobster’s projected DOS and local projected DOS.
The beforehand quantum-chemical calculation was performed with VASP.


#### completedos()
LobsterCompleteDos Object


#### pdos()

### List of Dict including numpy arrays with pdos. Access as pdos[atomindex]['orbitalstring']['Spin.up/Spin.down']()

#### tdos()

### Dos Object of the total density of states()

#### energies()

### numpy array of the energies at which the DOS was calculated (in eV, relative to Efermi)()

#### tdensities()

### tdensities[Spin.up]: numpy array of the total density of states for the Spin.up contribution at each of the()

#### energies()

### tdensities[Spin.down]: numpy array of the total density of states for the Spin.down contribution at each of the()

#### energies()
if is_spin_polarized=False:
tdensities[Spin.up]: numpy array of the total density of states


### itdensities:()

### itdensities[Spin.up]: numpy array of the total density of states for the Spin.up contribution at each of the()

#### energies()

### itdensities[Spin.down]: numpy array of the total density of states for the Spin.down contribution at each of the()

#### energies()
if is_spin_polarized=False:
itdensities[Spin.up]: numpy array of the total density of states


#### is_spin_polarized()

### Boolean. Tells if the system is spin polarized()

* **Parameters**


    * **doscar** – DOSCAR filename, typically “DOSCAR.lobster”


    * **structure_file** – for vasp, this is typically “POSCAR”


    * **structure** – instead of a structure file, the structure can be given
    directly. structure_file will be preferred.



#### _property_ completedos(_: [LobsterCompleteDos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.LobsterCompleteDos_ )
CompleteDos


* **Type**

    return



#### _property_ energies(_: ndarra_ )
Energies


* **Type**

    return



#### _property_ is_spin_polarized(_: boo_ )
Whether run is spin polarized.


* **Type**

    return



#### _property_ itdensities(_: ndarra_ )
integrated total densities as a np.ndarray


* **Type**

    return



#### _property_ pdos(_: lis_ )
Projected DOS


* **Type**

    return



#### _property_ tdensities(_: ndarra_ )
total densities as a np.ndarray


* **Type**

    return



#### _property_ tdos(_: [Dos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.Dos_ )
Total DOS


* **Type**

    return



### _class_ pymatgen.io.lobster.outputs.Fatband(filenames='.', vasprun='vasprun.xml', Kpointsfile='KPOINTS')
Bases: `object`

Reads in FATBAND_x_y.lobster files.

<!-- attribute: efermi

efermi that was read in from vasprun.xml -->
<!-- attribute: eigenvals
{Spin.up:[][],Spin.down:[][]}, the first index of the array
    [][] refers to the band and the second to the index of the
    kpoint. The kpoints are ordered according to the order of the
    kpoints array. If the band structure is not spin polarized, we
    only store one data set under Spin.up. -->
<!-- attribute: is_spinpolarized

Boolean that tells you whether this was a spin-polarized calculation -->
<!-- attribute: kpoints_array

list of kpoint as numpy arrays, in frac_coords of the given lattice by default -->
<!-- attribute: label_dict

(dict) of {} this link a kpoint (in frac coords or Cartesian coordinates depending on the coords). -->
<!-- attribute: lattice

lattice object of reciprocal lattice as read in from vasprun.xml -->
<!-- attribute: nbands

number of bands used in the calculation -->
<!-- attribute: p_eigenvals

dict of orbital projections as {spin: array of dict}.
The indices of the array are [band_index, kpoint_index].
The dict is then built the following way:
{"string of element": "string of orbital as read in from FATBAND file"}
If the band structure is not spin polarized, we only store one data set under Spin.up. -->
<!-- attribute: structure

structure read in from vasprun.xml -->

* **Parameters**


    * **filenames** (*list** or **string*) – can be a list of file names or a path to a folder folder from which all
    “FATBAND_\*” files will be read


    * **vasprun** – corresponding vasprun file


    * **Kpointsfile** – KPOINTS file for bandstructure calculation, typically “KPOINTS”.



#### get_bandstructure()
Returns a LobsterBandStructureSymmLine object which can be plotted with a normal BSPlotter.


### _class_ pymatgen.io.lobster.outputs.Grosspop(filename: str = 'GROSSPOP.lobster')
Bases: `object`

Class to read in GROSSPOP.lobster files.

<!-- attribute: list_dict_grosspop
which is a list of dicts including all information about the grosspopulations, one sample dict looks like this:
 {'element': 'O', 'Mulliken GP': {'2s': '1.80', '2p_y': '1.83', '2p_z': '1.79', '2p_x': '1.75', 'total': '7.18'},
  'Loewdin GP': {'2s': '1.60', '2p_y': '1.82', '2p_z': '1.77', '2p_x': '1.73', 'total': '6.92'}}
 The 0. entry of the list refers to the first atom in GROSSPOP.lobster and so on. -->

* **Parameters**

    **filename** – filename of the “GROSSPOP.lobster” file.



#### get_structure_with_total_grosspop(structure_filename: str)
Get a Structure with Mulliken and Loewdin total grosspopulations as site properties
:param structure_filename: filename of POSCAR
:type structure_filename: str


* **Returns**

    Structure Object with Mulliken and Loewdin total grosspopulations as site properties.



### _class_ pymatgen.io.lobster.outputs.Icohplist(are_coops: bool = False, are_cobis: bool = False, filename: str | None = None)
Bases: `object`

Class to read ICOHPLIST/ICOOPLIST files generated by LOBSTER.

<!-- attribute: are_coops
Boolean to indicate if the populations are COOPs or COHPs. -->
<!-- attribute: is_spin_polarized
Boolean to indicate if the calculation is spin polarized. -->
<!-- attribute: Icohplist
Dict containing the listfile data of the form:
   {bond: "length": bond length,
          "number_of_bonds": number of bonds
          "icohp": {Spin.up: ICOHP(Ef) spin up, Spin.down: ...}} -->
<!-- attribute: IcohpCollection
IcohpCollection Object -->

* **Parameters**


    * **are_coops** – Determines if the file is a list of ICOOPs.
    Defaults to False for ICOHPs.


    * **are_cobis** – Determines if the file is a list of ICOBIs.
    Defaults to False for ICOHPs.


    * **filename** – Name of the ICOHPLIST file. If it is None, the default
    file name will be chosen, depending on the value of are_coops.



#### _property_ icohpcollection()
IcohpCollection object.


* **Type**

    Returns



#### _property_ icohplist(_: dict[Any, dict[str, Any]_ )
icohplist compatible with older version of this class.


* **Type**

    Returns



### _class_ pymatgen.io.lobster.outputs.Lobsterout(filename='lobsterout')
Bases: `object`

Class to read in the lobsterout and evaluate the spilling, save the basis, save warnings, save infos.

> <!-- attribute: basis_functions
> list of basis functions that were used in lobster run as strings -->
> <!-- attribute: basis_type
> list of basis type that were used in lobster run as strings -->
> <!-- attribute: charge_spilling
> list of charge spilling (first entry: result for spin 1, second entry: result for spin 2 or not present) -->
> <!-- attribute: dft_program
> string representing the dft program used for the calculation of the wave function -->
> <!-- attribute: elements
> list of strings of elements that were present in lobster calculation -->
> <!-- attribute: has_charge
> Boolean, indicates that CHARGE.lobster is present -->
> <!-- attribute: has_cohpcar
> Boolean, indicates that COHPCAR.lobster and ICOHPLIST.lobster are present -->
> <!-- attribute: has_madelung
> Boolean, indicates that SitePotentials.lobster and MadelungEnergies.lobster are present -->
> <!-- attribute: has_coopcar
> Boolean, indicates that COOPCAR.lobster and ICOOPLIST.lobster are present -->
> <!-- attribute: has_cobicar
> Boolean, indicates that COBICAR.lobster and ICOBILIST.lobster are present -->
> <!-- attribute: has_doscar
> Boolean, indicates that DOSCAR.lobster is present -->
> <!-- attribute: has_doscar_lso
> Boolean, indicates that DOSCAR.LSO.lobster is present -->
> <!-- attribute: has_projection
> Boolean, indicates that projectionData.lobster is present -->
> <!-- attribute: has_bandoverlaps
> Boolean, indicates that bandOverlaps.lobster is present -->
> <!-- attribute: has_density_of_energies
> Boolean, indicates that DensityOfEnergy.lobster is present -->
> <!-- attribute: has_fatbands
> Boolean, indicates that fatband calculation was performed -->
> <!-- attribute: has_grosspopulation
> Boolean, indicates that GROSSPOP.lobster is present -->
> <!-- attribute: info_lines
> string with additional infos on the run -->
> <!-- attribute: info_orthonormalization
> string with infos on orthonormalization -->
> <!-- attribute: is_restart_from_projection
> Boolean that indicates that calculation was restartet from existing projection file -->
> <!-- attribute: lobster_version
> string that indicates Lobster version -->
> <!-- attribute: number_of_spins
> Integer indicating the number of spins -->
> <!-- attribute: number_of_threads
> integer that indicates how many threads were used -->
> <!-- attribute: timing
> dict with infos on timing -->
> <!-- attribute: total_spilling
> list of values indicating the total spilling for spin channel 1 (and spin channel 2) -->
> <!-- attribute: warning_lines
> string with all warnings -->

* **Parameters**

    **filename** – filename of lobsterout.



#### get_doc()
Returns: LobsterDict with all the information stored in lobsterout.


### _class_ pymatgen.io.lobster.outputs.MadelungEnergies(filename: str = 'MadelungEnergies.lobster')
Bases: `object`

Class to read MadelungEnergies.lobster files generated by LOBSTER.

<!-- attribute: madelungenergies_Mulliken
float that gives the madelung energy based on the Mulliken approach -->
<!-- attribute: madelungenergies_Loewdin
float that gives the madelung energy based on the Loewdin approach -->
<!-- attribute: ewald_splitting
Ewald Splitting parameter to compute SitePotentials -->

* **Parameters**

    **filename** – filename of the “MadelungEnergies.lobster” file.



### _class_ pymatgen.io.lobster.outputs.SitePotential(filename: str = 'SitePotentials.lobster')
Bases: `object`

Class to read SitePotentials.lobster files generated by LOBSTER.

<!-- attribute: atomlist
List of atoms in SitePotentials.lobster -->
<!-- attribute: types
List of types of atoms in SitePotentials.lobster -->
<!-- attribute: num_atoms
Number of atoms in SitePotentials.lobster -->
<!-- attribute: sitepotentials_Mulliken
List of Mulliken potentials of sites in SitePotentials.lobster -->
<!-- attribute: sitepotentials_Loewdin
List of Loewdin potentials of sites in SitePotentials.lobster -->
<!-- attribute: madelung_Mulliken
float that gives the madelung energy based on the Mulliken approach -->
<!-- attribute: madelung_Loewdin
float that gives the madelung energy based on the Loewdin approach -->
<!-- attribute: ewald_splitting
Ewald Splitting parameter to compute SitePotentials -->

* **Parameters**

    **filename** – filename for the SitePotentials file, typically “SitePotentials.lobster”.



#### get_structure_with_site_potentials(structure_filename)
Get a Structure with Mulliken and Loewdin charges as site properties
:param structure_filename: filename of POSCAR


* **Returns**

    Structure Object with Mulliken and Loewdin charges as site properties.



### _class_ pymatgen.io.lobster.outputs.Wavefunction(filename, structure)
Bases: `object`

Class to read in wave function files from Lobster and transfer them into an object of the type VolumetricData.

<!-- attribute: grid

grid for the wave function [Nx+1,Ny+1,Nz+1] -->
<!-- attribute: points

list of points -->
<!-- attribute: real

list of real part of wave function -->
<!-- attribute: imaginary

list of imaginary part of wave function -->
<!-- attribute: distance

list of distance to first point in wave function file -->

* **Parameters**


    * **filename** – filename of wavecar file from Lobster


    * **structure** – Structure object (e.g., created by Structure.from_file(“”)).



#### get_volumetricdata_density()
Will return a VolumetricData object including the imaginary part of the wave function.

Returns: VolumetricData object


#### get_volumetricdata_imaginary()
Will return a VolumetricData object including the imaginary part of the wave function.

Returns: VolumetricData object


#### get_volumetricdata_real()
Will return a VolumetricData object including the real part of the wave function.

Returns: VolumetricData object


#### set_volumetric_data(grid, structure)
Will create the VolumetricData Objects.


* **Parameters**


    * **grid** – grid on which wavefunction was calculated, e.g. [1,2,2]


    * **structure** – Structure object



#### write_file(filename='WAVECAR.vasp', part='real')
Will save the wavefunction in a file format that can be read by VESTA
This will only work if the wavefunction from lobster was constructed with:
“printLCAORealSpaceWavefunction kpoint 1 coordinates 0.0 0.0 0.0 coordinates 1.0 1.0 1.0 box bandlist 1 2 3 4
5 6 ”
or similar (the whole unit cell has to be covered!).


* **Parameters**


    * **filename** – Filename for the output, e.g., WAVECAR.vasp


    * **part** – which part of the wavefunction will be saved (“real” or “imaginary”)



### pymatgen.io.lobster.outputs.get_orb_from_str(orbs)

* **Parameters**

    **orbs** – list of two str, e.g. [“2p_x”, “3s”].



* **Returns**

    list of tw Orbital objects