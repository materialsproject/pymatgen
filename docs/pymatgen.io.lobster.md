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


### _class_ Lobsterin(settingsdict: dict)
Bases: `dict`, `MSONable`

This class can handle and generate lobsterin files
Furthermore, it can also modify INCAR files for lobster, generate KPOINT files for fatband calculations in Lobster,
and generate the standard primitive cells in a POSCAR file that are needed for the fatband calculations.
There are also several standard lobsterin files that can be easily generated.


* **Parameters**

    **settingsdict** – dict to initialize Lobsterin.



#### AVAILABLE_KEYWORDS(_ = ('COHPstartEnergy', 'COHPendEnergy', 'gaussianSmearingWidth', 'useDecimalPlaces', 'COHPSteps', 'basisSet', 'cohpGenerator', 'realspaceHamiltonian', 'realspaceOverlap', 'printPAWRealSpaceWavefunction', 'printLCAORealSpaceWavefunction', 'kSpaceCOHP', 'EwaldSum', 'saveProjectionToFile', 'skipdos', 'skipcohp', 'skipcoop', 'skipcobi', 'skipMadelungEnergy', 'loadProjectionFromFile', 'forceEnergyRange', 'DensityOfEnergy', 'BWDF', 'BWDFCOHP', 'skipPopulationAnalysis', 'skipGrossPopulation', 'userecommendedbasisfunctions', 'skipProjection', 'writeBasisFunctions', 'writeMatricesToFile', 'noFFTforVisualization', 'RMSp', 'onlyReadVasprun.xml', 'noMemoryMappedFiles', 'skipPAWOrthonormalityTest', 'doNotIgnoreExcessiveBands', 'doNotUseAbsoluteSpilling', 'skipReOrthonormalization', 'forceV1HMatrix', 'useOriginalTetrahedronMethod', 'forceEnergyRange', 'bandwiseSpilling', 'kpointwiseSpilling', 'LSODOS', 'basisfunctions', 'cohpbetween', 'createFatband'_ )

#### BOOLEAN_KEYWORDS(_ = ('saveProjectionToFile', 'skipdos', 'skipcohp', 'skipcoop', 'skipcobi', 'skipMadelungEnergy', 'loadProjectionFromFile', 'forceEnergyRange', 'DensityOfEnergy', 'BWDF', 'BWDFCOHP', 'skipPopulationAnalysis', 'skipGrossPopulation', 'userecommendedbasisfunctions', 'skipProjection', 'writeBasisFunctions', 'writeMatricesToFile', 'noFFTforVisualization', 'RMSp', 'onlyReadVasprun.xml', 'noMemoryMappedFiles', 'skipPAWOrthonormalityTest', 'doNotIgnoreExcessiveBands', 'doNotUseAbsoluteSpilling', 'skipReOrthonormalization', 'forceV1HMatrix', 'useOriginalTetrahedronMethod', 'forceEnergyRange', 'bandwiseSpilling', 'kpointwiseSpilling', 'LSODOS'_ )

#### FLOAT_KEYWORDS(_ = ('COHPstartEnergy', 'COHPendEnergy', 'gaussianSmearingWidth', 'useDecimalPlaces', 'COHPSteps'_ )

#### LISTKEYWORDS(_ = ('basisfunctions', 'cohpbetween', 'createFatband'_ )

#### STRING_KEYWORDS(_ = ('basisSet', 'cohpGenerator', 'realspaceHamiltonian', 'realspaceOverlap', 'printPAWRealSpaceWavefunction', 'printLCAORealSpaceWavefunction', 'kSpaceCOHP', 'EwaldSum'_ )

#### _get_nbands(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Get number of bands.


#### _static_ _get_potcar_symbols(POTCAR_input: str)
Will return the name of the species in the POTCAR.


* **Parameters**

    **POTCAR_input** (*str*) – string to potcar file



* **Returns**

    list of the names of the species in string format



#### as_dict()
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



* **Returns**

    Can be used to create new Lobsterin objects in

        standard_calculations_from_vasp_files as dict_for_basis




* **Return type**

    list[dict]



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



### get_all_possible_basis_combinations(min_basis: list, max_basis: list)

* **Parameters**


    * **min_basis** – list of basis entries: e.g., [‘Si 3p 3s ‘]


    * **max_basis** – list of basis entries: e.g., [‘Si 3p 3s ‘].



* **Returns**

    all possible combinations of basis functions, e.g. [[‘Si 3p 3s’]]



* **Return type**

    list[list[str]]


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


### _class_ ICOHPNeighborsInfo(total_icohp: float, list_icohps: list[float], n_bonds: int, labels: list[str], atoms: list[list[str]], central_isites: list[int] | None)
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


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('total_icohp', 'list_icohps', 'n_bonds', 'labels', 'atoms', 'central_isites'_ )

#### _classmethod_ _make(iterable)
Make a new ICOHPNeighborsInfo object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new ICOHPNeighborsInfo object replacing specified fields with new values


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


### _class_ LobsterLightStructureEnvironments(strategy, coordination_environments=None, all_nbs_sites=None, neighbors_sets=None, structure=None, valences=None, valences_origin=None)
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


* **Returns**

    Bson-serializable dict representation of the LightStructureEnvironments object.



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



* **Returns**

    LobsterLightStructureEnvironments



#### _property_ uniquely_determines_coordination_environments()
True if the coordination environments are uniquely determined.


### _class_ LobsterNeighbors(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), filename_ICOHP: str = 'ICOHPLIST.lobster', are_coops: bool = False, are_cobis: bool = False, valences: list[float] | None = None, limits: tuple[float, float] | None = None, additional_condition: int = 0, only_bonds_to: list[str] | None = None, perc_strength_ICOHP: float = 0.15, noise_cutoff: float = 0.1, valences_from_charges: bool = False, filename_CHARGE: str | None = None, which_charge: str = 'Mulliken', adapt_extremum_to_add_cond: bool = False, add_additional_data_sg: bool = False, filename_blist_sg1: str | None = None, filename_blist_sg2: str | None = None, id_blist_sg1: str = 'ICOOP', id_blist_sg2: str = 'ICOBI')
Bases: [`NearNeighbors`](pymatgen.analysis.md#pymatgen.analysis.local_env.NearNeighbors)

This class combines capabilities from LocalEnv and ChemEnv to determine coordination environments based on
bonding analysis.


* **Parameters**


    * **filename_ICOHP** – (str) Path to ICOHPLIST.lobster or ICOOPLIST.lobster or ICOBILIST.lobster


    * **structure** – (Structure) typically constructed by Structure.from_file(“POSCAR”)


    * **are_coops** – (bool) if True, the file is a ICOOPLIST.lobster and not a ICOHPLIST.lobster; only tested for
    ICOHPLIST.lobster so far


    * **are_cobis** – (bool) if True, the file is a ICOBILIST.lobster and not a ICOHPLIST.lobster


    * **valences** (*Mulliken**) **instead of*) – (list[float]): gives valence/charge for each element


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



#### _adapt_extremum_to_add_cond(list_icohps, percentage)
Convinicence method for returning the extremum of the given icohps or icoops or icobis list


* **Parameters**

    **list_icohps** – can be a list of icohps or icobis or icobis



* **Returns**

    min value of input list of icohps / max value of input list of icobis or icobis



* **Return type**

    float



#### _static_ _determine_unit_cell(site)
Based on the site it will determine the unit cell, in which this site is based.


* **Parameters**

    **site** – site object



#### _evaluate_ce(lowerlimit, upperlimit, only_bonds_to=None, additional_condition=0, perc_strength_ICOHP=0.15, adapt_extremum_to_add_cond=False)

* **Parameters**


    * **lowerlimit** – lower limit which determines the ICOHPs that are considered for the determination of the


    * **neighbors** –


    * **upperlimit** – upper limit which determines the ICOHPs that are considered for the determination of the


    * **neighbors** –


    * **only_bonds_to** – restricts the types of bonds that will be considered


    * **additional_condition** – Additional condition for the evaluation


    * **perc_strength_ICOHP** – will be used to determine how strong the ICOHPs (percentage\*strongest ICOHP) will be


    * **evalulation** (*that are still considered for the*) –


    * **adapt_extremum_to_add_cond** – will recalculate the limit based on the bonding type and not on the overall


    * **extremum.** –


Returns:


#### _find_environments(additional_condition, lowerlimit, upperlimit, only_bonds_to)
Will find all relevant neighbors based on certain restrictions.


* **Parameters**


    * **additional_condition** (*int*) – additional condition (see above)


    * **lowerlimit** (*float*) – lower limit that tells you which ICOHPs are considered


    * **upperlimit** (*float*) – upper limit that tells you which ICOHPs are considered


    * **only_bonds_to** (*list*) – list of str, e.g. [“O”] that will ensure that only bonds to “O” will be considered


Returns:


#### _find_relevant_atoms_additional_condition(isite, icohps, additional_condition)
Will find all relevant atoms that fulfill the additional_conditions.


* **Parameters**


    * **isite** – number of site in structure (starts with 0)


    * **icohps** – icohps


    * **additional_condition** (*int*) – additional condition


Returns:


#### _static_ _get_atomnumber(atomstring)
Return the number of the atom within the initial POSCAR (e.g., Return 0 for “Na1”).


* **Parameters**

    **atomstring** – string such as “Na1”



* **Returns**

    indicating the position in the POSCAR



* **Return type**

    int



#### _static_ _get_icohps(icohpcollection, isite, lowerlimit, upperlimit, only_bonds_to)
Return icohp dict for certain site.


* **Parameters**


    * **icohpcollection** – Icohpcollection object


    * **isite** (*int*) – number of a site


    * **lowerlimit** (*float*) – lower limit that tells you which ICOHPs are considered


    * **upperlimit** (*float*) – upper limit that tells you which ICOHPs are considered


    * **only_bonds_to** (*list*) – list of str, e.g. [“O”] that will ensure that only bonds to “O” will be considered


Returns:


#### _get_limit_from_extremum(icohpcollection, percentage=0.15, adapt_extremum_to_add_cond=False, additional_condition=0)
Return limits for the evaluation of the icohp values from an icohpcollection
Return -float(‘inf’), min(max_icohp\*0.15,-0.1). Currently only works for ICOHPs.


* **Parameters**


    * **icohpcollection** – icohpcollection object


    * **percentage** – will determine which ICOHPs or ICOOP or ICOBI will be considered


    * **value****)** (*(**only 0.15 from the maximum*) –


    * **adapt_extremum_to_add_cond** – should the extrumum be adapted to the additional condition


    * **additional_condition** – additional condition to determine which bonds are relevant



* **Returns**

    [-inf, min(strongest_icohp\*0.15,-noise_cutoff)] / [max(strongest_icohp\*0.15,

        noise_cutoff), inf]




* **Return type**

    tuple[float, float]



#### _get_plot_label(atoms, per_bond)

#### _static_ _split_string(s)
Will split strings such as “Na1” in “Na” and “1” and return “1”.


* **Parameters**

    **s** (*str*) – string



#### _property_ anion_types()
Return the types of anions present in crystal structure as a set


* **Returns**

    describing anions in the crystal structure.



* **Return type**

    set[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)]



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



* **Returns**

    label for cohp (str), CompleteCohp object which describes all cohps (coops or cobis)

        of the sites as given by isites and the other parameters




* **Return type**

    str



#### get_info_icohps_between_neighbors(isites=None, onlycation_isites=True)
Return infos about interactions between neighbors of a certain atom.


* **Parameters**


    * **isites** – list of site ids, if isite==None, all isites will be used


    * **onlycation_isites** – will only use cations, if isite==None


Returns

    ICOHPNeighborsInfo


#### get_info_icohps_to_neighbors(isites=None, onlycation_isites=True)
This method returns information on the icohps of neighbors for certain sites as identified by their site id.
This is useful for plotting the relevant cohps of a site in the structure.
(could be ICOOPLIST.lobster or ICOHPLIST.lobster or ICOBILIST.lobster)


* **Parameters**


    * **isites** – list of site ids. If isite==None, all isites will be used to add the icohps of the neighbors


    * **onlycation_isites** – if True and if isite==None, it will only analyse the sites of the cations



* **Returns**

    ICOHPNeighborsInfo



#### get_light_structure_environment(only_cation_environments=False, only_indices=None)
Return a LobsterLightStructureEnvironments object
if the structure only contains coordination environments smaller 13.


* **Parameters**


    * **only_cation_environments** – only data for cations will be returned


    * **only_indices** – will only evaluate the list of isites in this list



* **Returns**

    LobsterLightStructureEnvironments



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
Whether this NearNeighbors class can be used with Molecule objects?


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
Whether this NearNeighbors class can be used with Structure objects?

## pymatgen.io.lobster.outputs module

Module for reading Lobster output files. For more information
on LOBSTER see www.cohp.de.
If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
“Automated Bonding Analysis with Crystal Orbital Hamilton Populations”,
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.


### _class_ Bandoverlaps(filename: str = 'bandOverlaps.lobster')
Bases: `object`

Class to read in bandOverlaps.lobster files. These files are not created during every Lobster run.
.. attribute:: bandoverlapsdict

> A dictionary
> containing the band overlap data of the form: {spin: {“kpoint as string”: {“maxDeviation”:
> float that describes the max deviation, “matrix”: 2D array of the size number of bands
> times number of bands including the overlap matrices with}}}.


> * **type**

>     dict[Spin, Dict[str, Dict[str, Union[float, np.ndarray]]]]



#### maxDeviation()
A list of floats describing the maximal deviation for each problematic kpoint.


* **Type**

    list[float]



* **Parameters**

    **filename** – filename of the “bandOverlaps.lobster” file.



#### _read(contents: list, spin_numbers: list)
Will read in all contents of the file


* **Parameters**


    * **contents** – list of strings


    * **spin_numbers** – list of spin numbers depending on Lobster version.



#### has_good_quality_check_occupied_bands(number_occ_bands_spin_up: int, number_occ_bands_spin_down: int | None = None, spin_polarized: bool = False, limit_deviation: float = 0.1)
Will check if the deviation from the ideal bandoverlap of all occupied bands
is smaller or equal to limit_deviation.


* **Parameters**


    * **number_occ_bands_spin_up** (*int*) – number of occupied bands of spin up


    * **number_occ_bands_spin_down** (*int*) – number of occupied bands of spin down


    * **spin_polarized** (*bool*) – If True, then it was a spin polarized calculation


    * **limit_deviation** (*float*) – limit of the maxDeviation



* **Returns**

    Boolean that will give you information about the quality of the projection



#### has_good_quality_maxDeviation(limit_maxDeviation: float = 0.1)
Will check if the maxDeviation from the ideal bandoverlap is smaller or equal to limit_maxDeviation


* **Parameters**

    **limit_maxDeviation** – limit of the maxDeviation



* **Returns**

    Boolean that will give you information about the quality of the projection.



### _class_ Charge(filename: str = 'CHARGE.lobster')
Bases: `object`

Class to read CHARGE files generated by LOBSTER.


#### atomlist()
List of atoms in CHARGE.lobster.


* **Type**

    list[str]



#### types()
List of types of atoms in CHARGE.lobster.


* **Type**

    list[str]



#### Mulliken()
List of Mulliken charges of atoms in CHARGE.lobster.


* **Type**

    list[float]



#### Loewdin()
List of Loewdin charges of atoms in CHARGE.Loewdin.


* **Type**

    list[float]



#### num_atoms()
Number of atoms in CHARGE.lobster.


* **Type**

    int



* **Parameters**

    **filename** – filename for the CHARGE file, typically “CHARGE.lobster”.



#### get_structure_with_charges(structure_filename)
Get a Structure with Mulliken and Loewdin charges as site properties


* **Parameters**

    **structure_filename** – filename of POSCAR



* **Returns**

    Structure Object with Mulliken and Loewdin charges as site properties.



### _class_ Cohpcar(are_coops: bool = False, are_cobis: bool = False, filename: str | None = None)
Bases: `object`

Class to read COHPCAR/COOPCAR files generated by LOBSTER.


#### cohp_data()
A dictionary containing the COHP data of the form:
{bond: {“COHP”: {Spin.up: cohps, Spin.down:cohps},

> “ICOHP”: {Spin.up: icohps, Spin.down: icohps},
> “length”: bond length,
> “sites”: sites corresponding to the bond}

Also contains an entry for the average, which does not have a “length” key.


* **Type**

    dict[str, Dict[str, Any]]



#### efermi()
The Fermi energy in eV.


* **Type**

    float



#### energies()
Sequence of energies in eV. Note that LOBSTER shifts the energies
so that the Fermi energy is at zero.


* **Type**

    Sequence[float]



#### is_spin_polarized()
Boolean to indicate if the calculation is spin polarized.


* **Type**

    bool



#### orb_cohp()
A dictionary containing the orbital-resolved COHPs of the form:
orb_cohp[label] = {bond_data[“orb_label”]: {

> “COHP”: {Spin.up: cohps, Spin.down:cohps},
> “ICOHP”: {Spin.up: icohps, Spin.down: icohps},
> “orbitals”: orbitals,
> “length”: bond lengths,
> “sites”: sites corresponding to the bond},

}


* **Type**

    dict[str, Dict[str, Dict[str, Any]]]



* **Parameters**


    * **are_coops** – Determines if the file is a list of COHPs or COOPs.
    Default is False for COHPs.


    * **are_cobis** – Determines if the file is a list of COHPs or COOPs.
    Default is False for COHPs.


    * **filename** – Name of the COHPCAR file. If it is None, the default
    file name will be chosen, depending on the value of are_coops.



#### _static_ _get_bond_data(line: str)
Subroutine to extract bond label, site indices, and length from
a LOBSTER header line. The site indices are zero-based, so they
can be easily used with a Structure object.

Example header line: No.4:Fe1->Fe9(2.4524893531900283)
Example header line for orbtial-resolved COHP:

> No.1:Fe1[3p_x]->Fe2[3d_x^2-y^2](2.456180552772262)


* **Parameters**

    **line** – line in the COHPCAR header describing the bond.



* **Returns**

    Dict with the bond label, the bond length, a tuple of the site
    indices, a tuple containing the orbitals (if orbital-resolved),
    and a label for the orbitals (if orbital-resolved).



### _class_ Doscar(doscar: str = 'DOSCAR.lobster', structure_file: str | None = 'POSCAR', structure: [IStructure](pymatgen.core.md#pymatgen.core.structure.IStructure) | [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None)
Bases: `object`

Class to deal with Lobster’s projected DOS and local projected DOS.
The beforehand quantum-chemical calculation was performed with VASP.


#### completedos()
LobsterCompleteDos Object.


* **Type**

    [LobsterCompleteDos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.LobsterCompleteDos)



#### pdos()
List of Dict including numpy arrays with pdos. Access as
pdos[atomindex][‘orbitalstring’][‘Spin.up/Spin.down’].


* **Type**

    list



#### tdos()
Dos Object of the total density of states.


* **Type**

    [Dos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.Dos)



#### energies()
Numpy array of the energies at which the DOS was calculated
(in eV, relative to Efermi).


* **Type**

    numpy.ndarray



#### tdensities()
tdensities[Spin.up]: numpy array of the total density of states for
the Spin.up contribution at each of the energies. tdensities[Spin.down]: numpy array
of the total density of states for the Spin.down contribution at each of the energies.
If is_spin_polarized=False, tdensities[Spin.up]: numpy array of the total density of states.


* **Type**

    dict



#### itdensities()
itdensities[Spin.up]: numpy array of the total density of states for
the Spin.up contribution at each of the energies. itdensities[Spin.down]: numpy array
of the total density of states for the Spin.down contribution at each of the energies.
If is_spin_polarized=False, itdensities[Spin.up]: numpy array of the total density of states.


* **Type**

    dict



#### is_spin_polarized()
Boolean. Tells if the system is spin polarized.


* **Type**

    bool



* **Parameters**


    * **doscar** – DOSCAR filename, typically “DOSCAR.lobster”


    * **structure_file** – for vasp, this is typically “POSCAR”


    * **structure** – instead of a structure file, the structure can be given
    directly. structure_file will be preferred.



#### _parse_doscar()

#### _property_ completedos(_: [LobsterCompleteDos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.LobsterCompleteDos_ )
LobsterCompleteDos


#### _property_ energies(_: ndarra_ )
Energies


#### _property_ is_spin_polarized(_: boo_ )
Whether run is spin polarized.


#### _property_ itdensities(_: ndarra_ )
integrated total densities as a np.ndarray


#### _property_ pdos(_: lis_ )
Projected DOS


#### _property_ tdensities(_: ndarra_ )
total densities as a np.ndarray


#### _property_ tdos(_: [Dos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.Dos_ )
Total DOS


### _class_ Fatband(filenames='.', vasprun='vasprun.xml', Kpointsfile='KPOINTS')
Bases: `object`

Reads in FATBAND_x_y.lobster files.


#### efermi()
Fermi energy read in from vasprun.xml.


* **Type**

    float



#### eigenvals()
Eigenvalues as a dictionary of numpy arrays of shape (nbands, nkpoints).
The first index of the array refers to the band and the second to the index of the kpoint.
The kpoints are ordered according to the order of the kpoints_array attribute.
If the band structure is not spin polarized, we only store one data set under Spin.up.


* **Type**

    dict[[Spin](pymatgen.electronic_structure.md#pymatgen.electronic_structure.core.Spin), np.ndarray]



#### is_spin_polarized()
Boolean that tells you whether this was a spin-polarized calculation.


* **Type**

    bool



#### kpoints_array()
List of kpoints as numpy arrays, in frac_coords of the given
lattice by default.


* **Type**

    list[np.ndarray]



#### label_dict()
Dictionary that links a kpoint (in frac coords or Cartesian
coordinates depending on the coords attribute) to a label.


* **Type**

    dict[str, Union[str, np.ndarray]]



#### lattice()
Lattice object of reciprocal lattice as read in from vasprun.xml.


* **Type**

    [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice)



#### nbands()
Number of bands used in the calculation.


* **Type**

    int



#### p_eigenvals()
Dictionary of orbital projections as {spin: array of dict}.
The indices of the array are [band_index, kpoint_index].
The dict is then built the following way: {“string of element”: “string of orbital as read in
from FATBAND file”}. If the band structure is not spin polarized, we only store one data set under Spin.up.


* **Type**

    dict[[Spin](pymatgen.electronic_structure.md#pymatgen.electronic_structure.core.Spin), np.ndarray]



#### structure()
Structure read in from vasprun.xml.


* **Type**

    [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)



* **Parameters**


    * **filenames** (*list** or **string*) – can be a list of file names or a path to a folder from which all
    “FATBAND_\*” files will be read


    * **vasprun** – corresponding vasprun file


    * **Kpointsfile** – KPOINTS file for bandstructure calculation, typically “KPOINTS”.



#### get_bandstructure()
Returns a LobsterBandStructureSymmLine object which can be plotted with a normal BSPlotter.


### _class_ Grosspop(filename: str = 'GROSSPOP.lobster')
Bases: `object`

Class to read in GROSSPOP.lobster files.


#### list_dict_grosspop()
List of dictionaries
including all information about the grosspopulations. Each dictionary contains the following keys:
- ‘element’: The element symbol of the atom.
- ‘Mulliken GP’: A dictionary of Mulliken gross populations, where the keys are the orbital labels and the

> values are the corresponding gross populations as strings.


* ‘Loewdin GP’: A dictionary of Loewdin gross populations, where the keys are the orbital labels and the

    values are the corresponding gross populations as strings.

The 0th entry of the list refers to the first atom in GROSSPOP.lobster and so on.


* **Type**

    list[dict[str, str| dict[str, str]]]



* **Parameters**

    **filename** – filename of the “GROSSPOP.lobster” file.



#### get_structure_with_total_grosspop(structure_filename: str)
Get a Structure with Mulliken and Loewdin total grosspopulations as site properties


* **Parameters**

    **structure_filename** (*str*) – filename of POSCAR



* **Returns**

    Structure Object with Mulliken and Loewdin total grosspopulations as site properties.



### _class_ Icohplist(are_coops: bool = False, are_cobis: bool = False, filename: str | None = None)
Bases: `object`

Class to read ICOHPLIST/ICOOPLIST files generated by LOBSTER.


#### are_coops()
Indicates whether the object is consisting of COOPs.


* **Type**

    bool



#### is_spin_polarized()
Boolean to indicate if the calculation is spin polarized.


* **Type**

    bool



#### Icohplist()
Dict containing the
listfile data of the form: {

> bond: “length”: bond length,
> “number_of_bonds”: number of bonds
> “icohp”: {Spin.up: ICOHP(Ef) spin up, Spin.down: …}

}


* **Type**

    dict[str, Dict[str, Union[float, int, Dict[[Spin](pymatgen.electronic_structure.md#pymatgen.electronic_structure.core.Spin), float]]]]



#### IcohpCollection()
IcohpCollection Object.


* **Type**

    [IcohpCollection](pymatgen.electronic_structure.md#pymatgen.electronic_structure.cohp.IcohpCollection)



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



### _class_ Lobsterout(filename='lobsterout')
Bases: `object`

Class to read in the lobsterout and evaluate the spilling, save the basis, save warnings, save infos.


#### basis_functions()
List of basis functions that were used in lobster run as strings.


* **Type**

    list[str]



#### basis_type()
List of basis type that were used in lobster run as strings.


* **Type**

    list[str]



#### charge_spilling()
List of charge spilling (first entry: result for spin 1,
second entry: result for spin 2 or not present).


* **Type**

    list[float]



#### dft_program()
String representing the DFT program used for the calculation of the wave function.


* **Type**

    str



#### elements()
List of strings of elements that were present in lobster calculation.


* **Type**

    list[str]



#### has_charge()
Whether CHARGE.lobster is present.


* **Type**

    bool



#### has_cohpcar()
Whether COHPCAR.lobster and ICOHPLIST.lobster are present.


* **Type**

    bool



#### has_madelung()
Whether SitePotentials.lobster and MadelungEnergies.lobster are present.


* **Type**

    bool



#### has_coopcar()
Whether COOPCAR.lobster and ICOOPLIST.lobster are present.


* **Type**

    bool



#### has_cobicar()
Whether COBICAR.lobster and ICOBILIST.lobster are present.


* **Type**

    bool



#### has_doscar()
Whether DOSCAR.lobster is present.


* **Type**

    bool



#### has_doscar_lso()
Whether DOSCAR.LSO.lobster is present.


* **Type**

    bool



#### has_projection()
Whether projectionData.lobster is present.


* **Type**

    bool



#### has_bandoverlaps()
Whether bandOverlaps.lobster is present.


* **Type**

    bool



#### has_density_of_energies()
Whether DensityOfEnergy.lobster is present.


* **Type**

    bool



#### has_fatbands()
Whether fatband calculation was performed.


* **Type**

    bool



#### has_grosspopulation()
Whether GROSSPOP.lobster is present.


* **Type**

    bool



#### info_lines()
String with additional infos on the run.


* **Type**

    str



#### info_orthonormalization()
String with infos on orthonormalization.


* **Type**

    str



#### is_restart_from_projection()
Boolean that indicates that calculation was restarted
from existing projection file.


* **Type**

    bool



#### lobster_version()
String that indicates Lobster version.


* **Type**

    str



#### number_of_spins()
Integer indicating the number of spins.


* **Type**

    int



#### number_of_threads()
Integer that indicates how many threads were used.


* **Type**

    int



#### timing()
Dictionary with infos on timing.


* **Type**

    dict[str, float]



#### total_spilling()
List of values indicating the total spilling for spin
channel 1 (and spin channel 2).


* **Type**

    list[float]



#### warning_lines()
String with all warnings.


* **Type**

    str



* **Parameters**

    **filename** – filename of lobsterout.



#### _static_ _get_all_info_lines(data)

#### _static_ _get_all_warning_lines(data)

#### _static_ _get_dft_program(data)

#### _static_ _get_elements_basistype_basisfunctions(data)

#### _static_ _get_lobster_version(data)

#### _static_ _get_number_of_spins(data)

#### _static_ _get_spillings(data, number_of_spins)

#### _static_ _get_threads(data)

#### _static_ _get_timing(data)

#### _static_ _get_warning_orthonormalization(data)

#### _static_ _has_fatband(data)

#### get_doc()
Returns: LobsterDict with all the information stored in lobsterout.


### _class_ MadelungEnergies(filename: str = 'MadelungEnergies.lobster')
Bases: `object`

Class to read MadelungEnergies.lobster files generated by LOBSTER.


#### madelungenergies_Mulliken()
Float that gives the Madelung energy based on the Mulliken approach.


* **Type**

    float



#### madelungenergies_Loewdin()
Float that gives the Madelung energy based on the Loewdin approach.


* **Type**

    float



#### ewald_splitting()
Ewald splitting parameter to compute SitePotentials.


* **Type**

    float



* **Parameters**

    **filename** – filename of the “MadelungEnergies.lobster” file.



### _class_ SitePotential(filename: str = 'SitePotentials.lobster')
Bases: `object`

Class to read SitePotentials.lobster files generated by LOBSTER.


#### atomlist()
List of atoms in SitePotentials.lobster.


* **Type**

    list[str]



#### types()
List of types of atoms in SitePotentials.lobster.


* **Type**

    list[str]



#### num_atoms()
Number of atoms in SitePotentials.lobster.


* **Type**

    int



#### sitepotentials_Mulliken()
List of Mulliken potentials of sites in SitePotentials.lobster.


* **Type**

    list[float]



#### sitepotentials_Loewdin()
List of Loewdin potentials of sites in SitePotentials.lobster.


* **Type**

    list[float]



#### madelung_Mulliken()
Float that gives the Madelung energy based on the Mulliken approach.


* **Type**

    float



#### madelung_Loewdin()
Float that gives the Madelung energy based on the Loewdin approach.


* **Type**

    float



#### ewald_splitting()
Ewald Splitting parameter to compute SitePotentials.


* **Type**

    float



* **Parameters**

    **filename** – filename for the SitePotentials file, typically “SitePotentials.lobster”.



#### get_structure_with_site_potentials(structure_filename)
Get a Structure with Mulliken and Loewdin charges as site properties


* **Parameters**

    **structure_filename** – filename of POSCAR



* **Returns**

    Structure Object with Mulliken and Loewdin charges as site properties.



### _class_ Wavefunction(filename, structure)
Bases: `object`

Class to read in wave function files from Lobster and transfer them into an object of the type VolumetricData.


#### grid()
Grid for the wave function [Nx+1,Ny+1,Nz+1].


* **Type**

    tuple[int, int, int]



#### points()
List of points.


* **Type**

    list[Tuple[float, float, float]]



#### real()
List of real part of wave function.


* **Type**

    list[float]



#### imaginary()
List of imaginary part of wave function.


* **Type**

    list[float]



#### distance()
List of distance to first point in wave function file.


* **Type**

    list[float]



* **Parameters**


    * **filename** – filename of wavecar file from Lobster


    * **structure** – Structure object (e.g., created by Structure.from_file(“”)).



#### _static_ _parse_file(filename)

#### get_volumetricdata_density()
Will return a VolumetricData object including the imaginary part of the wave function.


* **Returns**

    VolumetricData



#### get_volumetricdata_imaginary()
Will return a VolumetricData object including the imaginary part of the wave function.


* **Returns**

    VolumetricData



#### get_volumetricdata_real()
Will return a VolumetricData object including the real part of the wave function.


* **Returns**

    VolumetricData



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



### get_orb_from_str(orbs)

* **Parameters**

    **orbs** – list of two str, e.g. [“2p_x”, “3s”].



* **Returns**

    list of tw Orbital objects