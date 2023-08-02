---
layout: default
title: pymatgen.io.lobster.lobsterenv.md
nav_exclude: true
---

# pymatgen.io.lobster.lobsterenv module

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
Bases: [`LightStructureEnvironments`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments)

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


#### _classmethod_ from_Lobster(list_ce_symbol, list_csm, list_permutation, list_neighsite, list_neighisite, structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), valences=None)
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


### _class_ pymatgen.io.lobster.lobsterenv.LobsterNeighbors(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), filename_ICOHP: str = 'ICOHPLIST.lobster', are_coops: bool = False, are_cobis: bool = False, valences: list[int | float] | None = None, limits: tuple[float, float] | None = None, additional_condition: int = 0, only_bonds_to: list[str] | None = None, perc_strength_ICOHP: float = 0.15, noise_cutoff: float = 0.1, valences_from_charges: bool = False, filename_CHARGE: str | None = None, which_charge: str = 'Mulliken', adapt_extremum_to_add_cond: bool = False, add_additional_data_sg: bool = False, filename_blist_sg1: str | None = None, filename_blist_sg2: str | None = None, id_blist_sg1: str = 'ICOOP', id_blist_sg2: str = 'ICOBI')
Bases: [`NearNeighbors`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors)

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


#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n, use_weights=False)
Get coordination number, CN, of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


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