---
layout: default
title: pymatgen.analysis.local_env.md
nav_exclude: true
---

# pymatgen.analysis.local_env module

This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures.


### _class_ pymatgen.analysis.local_env.BrunnerNN_real(tol: float = 0.0001, cutoff=8.0)
Bases: `NearNeighbors`

Determine coordination number using Brunner’s algorithm which counts the
atoms that are within the largest gap in differences in real space
interatomic distances. This algorithm uses Brunner’s method of
largest gap in interatomic distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    (default: 1E-4).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 8.0.



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.BrunnerNN_reciprocal(tol: float = 0.0001, cutoff=8.0)
Bases: `NearNeighbors`

Determine coordination number using Brunner’s algorithm which counts the
atoms that are within the largest gap in differences in real space
interatomic distances. This algorithm uses Brunner’s method of
largest reciprocal gap in interatomic distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    (default: 1E-4).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 8.0.



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.BrunnerNN_relative(tol: float = 0.0001, cutoff=8.0)
Bases: `NearNeighbors`

Determine coordination number using Brunner’s algorithm which counts the
atoms that are within the largest gap in differences in real space
interatomic distances. This algorithm uses Brunner’s method of
of largest relative gap in interatomic distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    (default: 1E-4).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 8.0.



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.CovalentBondNN(tol: float = 0.2, order=True)
Bases: `NearNeighbors`

Determine near-neighbor sites and bond orders using built-in
pymatgen.Molecule CovalentBond functionality.

NOTE: This strategy is only appropriate for molecules, and not for
structures.


* **Parameters**


    * **tol** (*float*) – Tolerance for covalent bond checking.


    * **order** (*bool*) – If True (default), this class will compute bond
    orders. If False, bond lengths will be computed.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_bonded_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), decorate: bool = False)
Obtain a MoleculeGraph object using this NearNeighbor
class.


* **Parameters**


    * **structure** – Molecule object.


    * **decorate** (*bool*) – whether to annotate site properties


    * **by** (*with order parameters using neighbors determined*) –


    * **class** (*this NearNeighbor*) –


Returns: a pymatgen.analysis.graphs.MoleculeGraph object


#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites and weights (orders) of bonds for a given
atom.


* **Parameters**


    * **structure** – input Molecule.


    * **n** – index of site for which to determine near neighbors.



* **Returns**

    [dict] representing a neighboring site and the type of
    bond present between site n and the neighboring site.



#### get_nn_shell_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), site_idx, shell)
Get a certain nearest neighbor shell for a certain site.

Determines all non-backtracking paths through the neighbor network
computed by get_nn_info. The weight is determined by multiplying
the weight of the neighbor at each hop through the network. For
example, a 2nd-nearest-neighbor that has a weight of 1 from its
1st-nearest-neighbor and weight 0.5 from the original site will
be assigned a weight of 0.5.

As this calculation may involve computing the nearest neighbors of
atoms multiple times, the calculation starts by computing all of the
neighbor info and then calling _get_nn_shell_info. If you are likely
to call this method for more than one site, consider calling get_all_nn
first and then calling this protected method yourself.


* **Parameters**


    * **structure** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Input structure


    * **site_idx** (*int*) – index of site for which to determine neighbor
    information.


    * **shell** (*int*) – Which neighbor shell to retrieve (1 == 1st NN shell)



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info.




#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.Critic2NN()
Bases: `NearNeighbors`

Performs a topological analysis using critic2 to obtain
neighbor information, using a sum of atomic charge
densities. If an actual charge density is available
(e.g. from a VASP CHGCAR), see Critic2Caller directly
instead.

Init for Critic2NN.


#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_bonded_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), decorate: bool = False)

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **decorate** (*bool**, **optional*) – Whether to decorate the structure. Defaults to False.



* **Returns**

    Bonded structure



* **Return type**

    [StructureGraph](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.CrystalNN(weighted_cn=False, cation_anion=False, distance_cutoffs=(0.5, 1), x_diff_weight=3.0, porous_adjustment=True, search_cutoff=7, fingerprint_length=None)
Bases: `NearNeighbors`

This is a custom near-neighbor method intended for use in all kinds of periodic structures
(metals, minerals, porous structures, etc). It is based on a Voronoi algorithm and uses the
solid angle weights to determine the probability of various coordination environments. The
algorithm can also modify probability using smooth distance cutoffs as well as Pauling
electronegativity differences. The output can either be the most probable coordination
environment or a weighted list of coordination environments.

Initialize CrystalNN with desired parameters. Default parameters assume
“chemical bond” type behavior is desired. For geometric neighbor
finding (e.g., structural framework), set (i) distance_cutoffs=None,
(ii) x_diff_weight=0 and (optionally) (iii) porous_adjustment=False
which will disregard the atomic identities and perform best for a purely
geometric match.


* **Parameters**


    * **weighted_cn** – (bool) if set to True, will return fractional weights
    for each potential near neighbor.


    * **cation_anion** – (bool) if set True, will restrict bonding targets to
    sites with opposite or zero charge. Requires an oxidation states
    on all sites in the structure.


    * **distance_cutoffs** – ([float, float]) - if not None, penalizes neighbor
    distances greater than sum of covalent radii plus
    distance_cutoffs[0]. Distances greater than covalent radii sum
    plus distance_cutoffs[1] are enforced to have zero weight.


    * **x_diff_weight** – (float) - if multiple types of neighbor elements are
    possible, this sets preferences for targets with higher
    electronegativity difference.


    * **porous_adjustment** – (bool) - if True, readjusts Voronoi weights to
    better describe layered / porous structures


    * **search_cutoff** – (float) cutoff in Angstroms for initial neighbor
    search; this will be adjusted if needed internally


    * **fingerprint_length** – (int) if a fixed_length CN “fingerprint” is
    desired from get_nn_data(), set this parameter



#### _class_ NNData(all_nninfo, cn_weights, cn_nninfo)
Bases: `tuple`

Create new instance of NNData(all_nninfo, cn_weights, cn_nninfo)


#### all_nninfo()
Alias for field number 0


#### cn_nninfo()
Alias for field number 2


#### cn_weights()
Alias for field number 1


#### get_cn(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int, \*\*kwargs)
Get coordination number, CN, of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).


    * **on_disorder** (*'take_majority_strict'** | **'take_majority_drop'** | **'take_max_species'** | **'error'*) – What to do when encountering a disordered structure. ‘error’ will raise ValueError.
    ‘take_majority_strict’ will use the majority specie on each site and raise
    ValueError if no majority exists. ‘take_max_species’ will use the first max specie
    on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, ‘error’ and ‘take_majority_strict’
    will raise ValueError, while ‘take_majority_drop’ ignores this site altogether and
    ‘take_max_species’ will use Fe as the site specie.



* **Returns**

    coordination number.



* **Return type**

    cn (int or float)



#### get_cn_dict(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int, use_weights: bool = False, \*\*kwargs)
Get coordination number, CN, of each element bonded to site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).



* **Returns**

    dictionary of CN of each element bonded to site



* **Return type**

    cn (dict)



#### get_nn_data(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int, length=None)
The main logic of the method to compute near neighbor.


* **Parameters**


    * **structure** – (Structure) enclosing structure object


    * **n** – (int) index of target site to get NN info for


    * **length** – (int) if set, will return a fixed range of CN numbers



* **Returns**


    * all near neighbor sites with weights


    * a dict of CN -> weight


    * a dict of CN -> associated near neighbor sites




* **Return type**

    a namedtuple (NNData) object that contains



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor information.


* **Parameters**


    * **structure** – (Structure) pymatgen Structure


    * **n** – (int) index of target site



* **Returns**

    each dictionary provides information

        about a single near neighbor, where key ‘site’ gives access to the
        corresponding Site object, ‘image’ gives the image location, and
        ‘weight’ provides the weight that a given near-neighbor site contributes
        to the coordination number (1 or smaller), ‘site_index’ gives index of
        the corresponding site in the original structure.




* **Return type**

    siw (list[dict])



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



#### _static_ transform_to_length(nn_data, length)
Given NNData, transforms data to the specified fingerprint length
:param nn_data: (NNData)
:param length: (int) desired length of NNData.


### _class_ pymatgen.analysis.local_env.CutOffDictNN(cut_off_dict=None)
Bases: `NearNeighbors`

A basic NN class using a dictionary of fixed cut-off distances.
Only pairs of elements listed in the cut-off dictionary are considered
during construction of the neighbor lists.

Omit passing a dictionary for a Null/Empty NN class.


* **Parameters**


    * **cut_off_dict** (*dict**[**str**, **float**]*) – a dictionary


    * **distances** (*of cut-off*) – 2.0} for


    * **{** (*e.g.*) – 2.0} for


    * **Angstroms.** (*a maximum Fe-O bond length** of **2*) –


    * **listed** (*Bonds will only be created between pairs*) –


    * **dictionary.** (*in the cut-off*) –


    * **decorated** (*If your structure is oxidation state*) –


:param :
:param the cut-off distances will have to explicitly include:
:param the oxidation state: 2.0}.
:type the oxidation state: ‘Fe2+’, ‘O2-’
:param e.g. {: 2.0}.
:type e.g. {: ‘Fe2+’, ‘O2-’


#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### _static_ from_preset(preset)
Initialize a CutOffDictNN according to a preset set of cut-offs.


* **Parameters**

    **preset** (*str*) – A preset name. The list of supported presets are:


    * ”vesta_2019”: The distance cut-offs used by the VESTA
    visualisation program.




* **Returns**

    A CutOffDictNN using the preset cut-off dictionary.



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.EconNN(tol: float = 0.2, cutoff: float = 10.0, cation_anion: bool = False, use_fictive_radius: bool = False)
Bases: `NearNeighbors`

Determines the average effective coordination number for each cation in a
given structure using Hoppe’s algorithm.

This method follows the procedure outlined in:

Hoppe, Rudolf. “Effective coordination numbers (ECoN) and mean fictive ionic
radii (MEFIR).” Zeitschrift für Kristallographie-Crystalline Materials
150.1-4 (1979): 23-52.


* **Parameters**


    * **tol** – Tolerance parameter for bond determination.


    * **cutoff** – Cutoff radius in Angstrom to look for near-neighbor atoms.


    * **cation_anion** – If set to True, will restrict bonding targets to
    sites with opposite or zero charge. Requires an oxidation states
    on all sites in the structure.


    * **use_fictive_radius** – Whether to use the fictive radius in the
    EcoN calculation. If False, the bond distance will be used.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.IsayevNN(tol: float = 0.25, targets: [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | list[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)] | None = None, cutoff: float = 13.0, allow_pathological: bool = False, extra_nn_info: bool = True, compute_adj_neighbors: bool = True)
Bases: `VoronoiNN`

Uses the algorithm defined in 10.1038/ncomms15679.

Sites are considered neighbors if (i) they share a Voronoi facet and (ii) the
bond distance is less than the sum of the Cordero covalent radii + 0.25 Å.


* **Parameters**


    * **tol** – Tolerance in Å for bond distances that are considered coordinated.


    * **targets** – Target element(s).


    * **cutoff** – Cutoff radius in Angstrom to look for near-neighbor atoms.


    * **allow_pathological** – Whether to allow infinite vertices in Voronoi
    coordination.


    * **extra_nn_info** – Add all polyhedron info to get_nn_info.


    * **compute_adj_neighbors** – Whether to compute which neighbors are adjacent. Turn
    off for faster performance.



#### get_all_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.



* **Returns**

    List of near neighbor information for each site. See get_nn_info for the
    format of the data for each site.



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor site information.

Gets the associated image locations and weights of the site with index n
in structure using Voronoi decomposition and distance cutoff.


* **Parameters**


    * **structure** – Input structure.


    * **n** – Index of site for which to determine near-neighbor sites.



* **Returns**

    List of dicts containing the near-neighbor information. Each dict has the
    keys:


    * ”site”: The near-neighbor site.


    * ”image”: The periodic image of the near-neighbor site.


    * ”weight”: The face weight of the Voronoi decomposition.


    * ”site_index”: The index of the near-neighbor site in the original
    structure.




### _class_ pymatgen.analysis.local_env.JmolNN(tol: float = 0.45, min_bond_distance: float = 0.4, el_radius_updates: dict[SpeciesLike, float] | None = None)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using an emulation
of Jmol’s default autoBond() algorithm. This version of the algorithm
does not take into account any information regarding known charge
states.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    Defaults to 0.56.


    * **min_bond_distance** (*float*) – minimum bond distance for consideration
    Defaults to 0.4.


    * **el_radius_updates** – (dict) symbol->float to override default atomic
    radii table values.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_max_bond_distance(el1_sym, el2_sym)
Use Jmol algorithm to determine bond length from atomic parameters
:param el1_sym: (str) symbol of atom 1
:param el2_sym: (str) symbol of atom 2.

Returns: (float) max bond length


#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the bond identification
algorithm underlying Jmol.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    tuples, each one

        of which represents a neighbor site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.LocalStructOrderParams(types, parameters=None, cutoff=-10.0)
Bases: `object`

This class permits the calculation of various types of local
structure order parameters.


* **Parameters**


    * **types** (*[**string**]*) – list of strings representing the types of
    order parameters to be calculated. Note that multiple
    mentions of the same type may occur. Currently available
    types recognize following environments:

    > ”cn”: simple coordination number—normalized

    >     if desired;

    > ”sgl_bd”: single bonds;
    > “bent”: bent (angular) coordinations

    > > (Zimmermann & Jain, in progress, 2017);

    > ”T”: T-shape coordinations;
    > “see_saw_rect”: see saw-like coordinations;
    > “tet”: tetrahedra

    > > (Zimmermann et al., submitted, 2017);

    > ”oct”: octahedra

    >     (Zimmermann et al., submitted, 2017);

    > ”bcc”: body-centered cubic environments (Peters,


    >         1. Chem. Phys., 131, 244103, 2009);

    > ”tri_plan”: trigonal planar environments;
    > “sq_plan”: square planar environments;
    > “pent_plan”: pentagonal planar environments;
    > “tri_pyr”: trigonal pyramids (coordinated atom is in

    > > the center of the basal plane);

    > ”sq_pyr”: square pyramids;
    > “pent_pyr”: pentagonal pyramids;
    > “hex_pyr”: hexagonal pyramids;
    > “tri_bipyr”: trigonal bipyramids;
    > “sq_bipyr”: square bipyramids;
    > “pent_bipyr”: pentagonal bipyramids;
    > “hex_bipyr”: hexagonal bipyramids;
    > “cuboct”: cuboctahedra;
    > “q2”: motif-unspecific bond orientational order

    > > parameter (BOOP) of weight l=2 (Steinhardt
    > > et al., Phys. Rev. B, 28, 784-805, 1983);

    > ”q4”: BOOP of weight l=4;
    > “q6”: BOOP of weight l=6.
    > “reg_tri”: regular triangle with varying height

    > > to basal plane;

    > ”sq”: square coordination (cf., “reg_tri”);
    > “oct_legacy”: original Peters-style OP recognizing

    > > octahedral coordination environments
    > > (Zimmermann et al., J. Am. Chem. Soc.,
    > > 137, 13352-13361, 2015) that can, however,
    > > produce small negative values sometimes.

    > ”sq_pyr_legacy”: square pyramids (legacy);



    * **parameters** (*[**dict**]*) – list of dictionaries
    that store float-type parameters associated with the
    definitions of the different order parameters
    (length of list = number of OPs). If an entry
    is None, default values are used that are read from
    the op_params.yaml file. With few exceptions, 9 different
    parameters are used across all OPs:

    > ”norm”: normalizing constant (used in “cn”

    >     (default value: 1)).

    > ”TA”: target angle (TA) in fraction of 180 degrees

    >     (“bent” (1), “tet” (0.6081734479693927),
    >     “tri_plan” (0.66666666667), “pent_plan” (0.6),
    >     “sq_pyr_legacy” (0.5)).

    > ”IGW_TA”: inverse Gaussian width (IGW) for penalizing

    >     angles away from the target angle in inverse
    >     fractions of 180 degrees to (“bent” and “tet” (15),
    >     “tri_plan” (13.5), “pent_plan” (18),
    >     “sq_pyr_legacy” (30)).

    > ”IGW_EP”: IGW for penalizing angles away from the

    >     equatorial plane (EP) at 90 degrees (“T”, “see_saw_rect”,
    >     “oct”, “sq_plan”, “tri_pyr”, “sq_pyr”, “pent_pyr”,
    >     “hex_pyr”, “tri_bipyr”, “sq_bipyr”, “pent_bipyr”,
    >     “hex_bipyr”, and “oct_legacy” (18)).

    > ”fac_AA”: factor applied to azimuth angle (AA) in cosine

    >     term (“T”, “tri_plan”, and “sq_plan” (1), “tet”,
    >     “tri_pyr”, and “tri_bipyr” (1.5), “oct”, “sq_pyr”,
    >     “sq_bipyr”, and “oct_legacy” (2), “pent_pyr”
    >     and “pent_bipyr” (2.5), “hex_pyr” and
    >     “hex_bipyr” (3)).

    > ”exp_cos_AA”: exponent applied to cosine term of AA

    >     (“T”, “tet”, “oct”, “tri_plan”, “sq_plan”,
    >     “tri_pyr”, “sq_pyr”, “pent_pyr”, “hex_pyr”,
    >     “tri_bipyr”, “sq_bipyr”, “pent_bipyr”, “hex_bipyr”,
    >     and “oct_legacy” (2)).

    > ”min_SPP”: smallest angle (in radians) to consider

    >     a neighbor to be
    >     at South pole position (“see_saw_rect”, “oct”, “bcc”,
    >     “sq_plan”, “tri_bipyr”, “sq_bipyr”, “pent_bipyr”,
    >     “hex_bipyr”, “cuboct”, and “oct_legacy”
    >     (2.792526803190927)).

    > ”IGW_SPP”: IGW for penalizing angles away from South

    >     pole position (“see_saw_rect”, “oct”, “bcc”, “sq_plan”,
    >     “tri_bipyr”, “sq_bipyr”, “pent_bipyr”, “hex_bipyr”,
    >     “cuboct”, and “oct_legacy” (15)).

    > ”w_SPP”: weight for South pole position relative to

    >     equatorial positions (“see_saw_rect” and “sq_plan” (1),
    >     “cuboct” (1.8), “tri_bipyr” (2), “oct”,
    >     “sq_bipyr”, and “oct_legacy” (3), “pent_bipyr” (4),
    >     “hex_bipyr” (5), “bcc” (6)).



    * **cutoff** (*float*) – Cutoff radius to determine which nearest
    neighbors are supposed to contribute to the order
    parameters. If the value is negative the neighboring
    sites found by distance and cutoff radius are further
    pruned using the get_nn method from the
    VoronoiNN class.



#### compute_trigonometric_terms(thetas, phis)
Computes trigonometric terms that are required to
calculate bond orientational order parameters using
internal variables.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.
    The list of
    azimuth angles of all neighbors in radians. The list of
    azimuth angles is expected to have the same size as the
    list of polar angles; otherwise, a ValueError is raised.
    Also, the two lists of angles have to be coherent in
    order. That is, it is expected that the order in the list
    of azimuth angles corresponds to a distinct sequence of
    neighbors. And, this sequence has to equal the sequence
    of neighbors in the list of polar angles.



#### get_order_parameters(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int, indices_neighs: list[int] | None = None, tol: float = 0.0, target_spec=None)
Compute all order parameters of site n.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in input structure,
    for which OPs are to be
    calculated. Note that we do not use the sites iterator
    here, but directly access sites via struct[index].


    * **indices_neighs** (*list**[**int**]*) – list of indices of those neighbors
    in Structure object
    structure that are to be considered for OP computation.
    This optional argument overwrites the way neighbors are
    to be determined as defined in the constructor (i.e.,
    Voronoi coordination finder via negative cutoff radius
    vs constant cutoff radius if cutoff was positive).
    We do not use information about the underlying
    structure lattice if the neighbor indices are explicitly
    provided. This has two important consequences. First,
    the input Structure object can, in fact, be a
    simple list of Site objects. Second, no nearest images
    of neighbors are determined when providing an index list.
    Note furthermore that this neighbor
    determination type ignores the optional target_spec
    argument.


    * **tol** (*float*) – threshold of weight
    (= solid angle / maximal solid angle)
    to determine if a particular pair is
    considered neighbors; this is relevant only in the case
    when Voronoi polyhedra are used to determine coordination


    * **target_spec** ([*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – target species to be considered
    when calculating the order
    parameters of site n; None includes all species of input
    structure.



* **Returns**

    representing order parameters. Should it not be
    possible to compute a given OP for a conceptual reason, the
    corresponding entry is None instead of a float. For Steinhardt
    et al.’s bond orientational OPs and the other geometric OPs
    (“tet”, “oct”, “bcc”, etc.),
    this can happen if there is a single
    neighbor around site n in the structure because that
    does not permit calculation of angles between multiple
    neighbors.



* **Return type**

    [floats]



#### get_parameters(index)
Returns list of floats that represents
the parameters associated
with calculation of the order
parameter that was defined at the index provided.
Attention: the parameters do not need to equal those originally
inputted because of processing out of efficiency reasons.


* **Parameters**

    **index** (*int*) – index of order parameter for which associated parameters
    are to be returned.



* **Returns**

    parameters of a given OP.



* **Return type**

    [float]



#### get_q2(thetas=None, phis=None)
Calculates the value of the bond orientational order parameter of
weight l=2. If the function is called with non-empty lists of
polar and azimuthal angles the corresponding trigonometric terms
are computed afresh. Otherwise, it is expected that the
compute_trigonometric_terms function has been just called.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.



* **Returns**

    bond orientational order parameter of weight l=2

        corresponding to the input angles thetas and phis.




* **Return type**

    float



#### get_q4(thetas=None, phis=None)
Calculates the value of the bond orientational order parameter of
weight l=4. If the function is called with non-empty lists of
polar and azimuthal angles the corresponding trigonometric terms
are computed afresh. Otherwise, it is expected that the
compute_trigonometric_terms function has been just called.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.



* **Returns**

    bond orientational order parameter of weight l=4

        corresponding to the input angles thetas and phis.




* **Return type**

    float



#### get_q6(thetas=None, phis=None)
Calculates the value of the bond orientational order parameter of
weight l=6. If the function is called with non-empty lists of
polar and azimuthal angles the corresponding trigonometric terms
are computed afresh. Otherwise, it is expected that the
compute_trigonometric_terms function has been just called.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.



* **Returns**

    bond orientational order parameter of weight l=6

        corresponding to the input angles thetas and phis.




* **Return type**

    float



#### get_type(index)
Return type of order parameter at the index provided and
represented by a short string.


* **Parameters**

    **index** (*int*) – index of order parameter for which type is
    to be returned.



* **Returns**

    OP type.



* **Return type**

    str



#### _property_ last_nneigh()
Returns:
int: the number of neighbors encountered during the most

> recent order parameter calculation. A value of -1 indicates
> that no such calculation has yet been performed for this
> instance.


#### _property_ num_ops()
Returns:
int: the number of different order parameters that are targeted

> to be calculated.


### _class_ pymatgen.analysis.local_env.MinimumDistanceNN(tol: float = 0.1, cutoff=10, get_all_sites=False)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using the
nearest neighbor(s) at distance, d_min, plus all neighbors
within a distance (1 + tol) \* d_min, where tol is a
(relative) distance tolerance parameter.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for neighbor identification
    (default: 0.1).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10).


    * **get_all_sites** (*bool*) – If this is set to True then the neighbor
    sites are only determined by the cutoff radius, tol is ignored.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the closest neighbor
distance-based method.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    dicts with (Site, array, float) each one of which represents a

        neighbor site, its image location, and its weight.




* **Return type**

    siw (list[dict])



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.MinimumOKeeffeNN(tol: float = 0.1, cutoff=10)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using the
neighbor(s) at closest relative distance, d_min_OKeffee, plus some
relative tolerance, where bond valence parameters from O’Keeffe’s
bond valence method (J. Am. Chem. Soc. 1991, 3226-3229) are used
to calculate relative distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for neighbor identification
    (default: 0.1).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10).



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the closest relative
neighbor distance-based method with O’Keeffe parameters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    tuples, each one

        of which represents a neighbor site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.MinimumVIRENN(tol: float = 0.1, cutoff=10)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using the
neighbor(s) at closest relative distance, d_min_VIRE, plus some
relative tolerance, where atom radii from the
ValenceIonicRadiusEvaluator (VIRE) are used
to calculate relative distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for neighbor identification
    (default: 0.1).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10).



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the closest relative
neighbor distance-based method with VIRE atomic/ionic radii.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    tuples, each one

        of which represents a neighbor site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.NearNeighbors()
Bases: `object`

Base class to determine near neighbors that typically include nearest
neighbors and others that are within some tolerable distance.


#### _property_ extend_structure_molecules(_: boo_ )
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_all_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Get a listing of all neighbors for all sites in a structure.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure



* **Returns**

    List of NN site information for each site in the structure. Each

        entry has the same format as get_nn_info




#### get_bonded_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), decorate: bool = False, weights: bool = True, edge_properties: bool = False, on_disorder: on_disorder_options = 'take_majority_strict')
Obtain a StructureGraph object using this NearNeighbor
class. Requires the optional dependency networkx
(pip install networkx).


* **Parameters**


    * **structure** – Structure object.


    * **decorate** (*bool*) – whether to annotate site properties with order parameters using neighbors
    determined by this NearNeighbor class


    * **weights** (*bool*) – whether to include edge weights from NearNeighbor class in StructureGraph


    * **edge_properties** (*bool*) –


    * **on_disorder** (*'take_majority_strict'** | **'take_majority_drop'** | **'take_max_species'** | **'error'*) – What to do when encountering a disordered structure. ‘error’ will raise ValueError.
    ‘take_majority_strict’ will use the majority specie on each site and raise
    ValueError if no majority exists. ‘take_max_species’ will use the first max specie
    on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, ‘error’ and ‘take_majority_strict’
    will raise ValueError, while ‘take_majority_drop’ ignores this site altogether and
    ‘take_max_species’ will use Fe as the site specie.


Returns: a pymatgen.analysis.graphs.StructureGraph object


#### get_cn(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int, use_weights: bool = False, on_disorder: Literal['take_majority_strict', 'take_majority_drop', 'take_max_species', 'error'] = 'take_majority_strict')
Get coordination number, CN, of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True) to use weights for computing the coordination
    number or not (False, default: each coordinated site has equal weight).


    * **on_disorder** (*'take_majority_strict'** | **'take_majority_drop'** | **'take_max_species'** | **'error'*) – What to do when encountering a disordered structure. ‘error’ will raise ValueError.
    ‘take_majority_strict’ will use the majority specie on each site and raise
    ValueError if no majority exists. ‘take_max_species’ will use the first max specie
    on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, ‘error’ and ‘take_majority_strict’
    will raise ValueError, while ‘take_majority_drop’ ignores this site altogether and
    ‘take_max_species’ will use Fe as the site specie.



* **Returns**

    coordination number.



* **Return type**

    cn (int or float)



#### get_cn_dict(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int, use_weights: bool = False)
Get coordination number, CN, of each element bonded to site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).



* **Returns**

    dictionary of CN of each element bonded to site



* **Return type**

    cn (dict)



#### get_local_order_parameters(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Calculate those local structure order parameters for
the given site whose ideal CN corresponds to the
underlying motif (e.g., CN=4, then calculate the
square planar, tetrahedral, see-saw-like,
rectangular see-saw-like order parameters).


* **Parameters**


    * **structure** – Structure object


    * **n** (*int*) – site index.


Returns (dict[str, float]):

    A dict of order parameters (values) and the
    underlying motif type (keys; for example, tetrahedral).


#### get_nn(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get near neighbors of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in structure for which to determine
    neighbors.



* **Returns**

    near neighbors.



* **Return type**

    sites (list of Site objects)



#### get_nn_images(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get image location of all near neighbors of site with index n in
structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine the image
    location of near neighbors.



* **Returns**

    image locations of

        near neighbors.




* **Return type**

    images (list of 3D integer array)



#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    information.



* **Returns**

    each dictionary provides information

        about a single near neighbor, where key ‘site’ gives access to the
        corresponding Site object, ‘image’ gives the image location, and
        ‘weight’ provides the weight that a given near-neighbor site contributes
        to the coordination number (1 or smaller), ‘site_index’ gives index of
        the corresponding site in the original structure.




* **Return type**

    siw (list[dict])



#### get_nn_shell_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), site_idx, shell)
Get a certain nearest neighbor shell for a certain site.

Determines all non-backtracking paths through the neighbor network
computed by get_nn_info. The weight is determined by multiplying
the weight of the neighbor at each hop through the network. For
example, a 2nd-nearest-neighbor that has a weight of 1 from its
1st-nearest-neighbor and weight 0.5 from the original site will
be assigned a weight of 0.5.

As this calculation may involve computing the nearest neighbors of
atoms multiple times, the calculation starts by computing all of the
neighbor info and then calling _get_nn_shell_info. If you are likely
to call this method for more than one site, consider calling get_all_nn
first and then calling this protected method yourself.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **site_idx** (*int*) – index of site for which to determine neighbor
    information.


    * **shell** (*int*) – Which neighbor shell to retrieve (1 == 1st NN shell)



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info.




#### get_weights_of_nn_sites(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get weight associated with each near neighbor of site with
index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine the weights.



* **Returns**

    near-neighbor weights.



* **Return type**

    weights (list of floats)



#### _property_ molecules_allowed(_: boo_ )
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed(_: boo_ )
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.OpenBabelNN(order=True)
Bases: `NearNeighbors`

Determine near-neighbor sites and bond orders using OpenBabel API.

NOTE: This strategy is only appropriate for molecules, and not for
structures.


* **Parameters**


    * **order** (*bool*) – True if bond order should be returned as a weight, False


    * **weight.** (*if bond length should be used as a*) –



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_bonded_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), decorate: bool = False)
Obtain a MoleculeGraph object using this NearNeighbor
class. Requires the optional dependency networkx
(pip install networkx).


* **Parameters**


    * **structure** – Molecule object.


    * **decorate** (*bool*) – whether to annotate site properties


    * **by** (*with order parameters using neighbors determined*) –


    * **class** (*this NearNeighbor*) –


Returns: a pymatgen.analysis.graphs.MoleculeGraph object


#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites and weights (orders) of bonds for a given
atom.


* **Parameters**


    * **structure** – Molecule object.


    * **n** – index of site for which to determine near neighbors.



* **Returns**

    representing a neighboring site and the type of
    bond present between site n and the neighboring site.



* **Return type**

    (dict)



#### get_nn_shell_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), site_idx, shell)
Get a certain nearest neighbor shell for a certain site.

Determines all non-backtracking paths through the neighbor network
computed by get_nn_info. The weight is determined by multiplying
the weight of the neighbor at each hop through the network. For
example, a 2nd-nearest-neighbor that has a weight of 1 from its
1st-nearest-neighbor and weight 0.5 from the original site will
be assigned a weight of 0.5.

As this calculation may involve computing the nearest neighbors of
atoms multiple times, the calculation starts by computing all of the
neighbor info and then calling _get_nn_shell_info. If you are likely
to call this method for more than one site, consider calling get_all_nn
first and then calling this protected method yourself.


* **Parameters**


    * **structure** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Input structure


    * **site_idx** (*int*) – index of site for which to determine neighbor
    information.


    * **shell** (*int*) – Which neighbor shell to retrieve (1 == 1st NN shell)



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info.




#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ pymatgen.analysis.local_env.ValenceIonicRadiusEvaluator(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Bases: `object`

Computes site valences and ionic radii for a structure using bond valence
analyzer.


* **Parameters**

    **structure** – pymatgen.core.structure.Structure.



#### _property_ radii()
List of ionic radii of elements in the order of sites.


#### _property_ structure()
Returns oxidation state decorated structure.


#### _property_ valences()
List of oxidation states of elements in the order of sites.


### _class_ pymatgen.analysis.local_env.VoronoiNN(tol=0, targets=None, cutoff=13.0, allow_pathological=False, weight='solid_angle', extra_nn_info=True, compute_adj_neighbors=True)
Bases: `NearNeighbors`

Uses a Voronoi algorithm to determine near neighbors for each site in a
structure.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for near-neighbor finding. Faces that are
    smaller than tol fraction of the largest face are not included in the
    tessellation. (default: 0).


    * **targets** ([*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)* or **list** of **Elements*) – target element(s).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 13.0.


    * **allow_pathological** (*bool*) – whether to allow infinite vertices in
    determination of Voronoi coordination.


    * **weight** (*string*) – available in get_voronoi_polyhedra)


    * **extra_nn_info** (*bool*) –


    * **compute_adj_neighbors** (*bool*) – adjacent. Turn off for faster performance.



#### get_all_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.



* **Returns**

    All nn info for all sites.



#### get_all_voronoi_polyhedra(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Get the Voronoi polyhedra for all site in a simulation cell.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to be evaluated



* **Returns**

    A dict of sites sharing a common Voronoi facet with the site
    n mapped to a directory containing statistics about the facet:

    >
    > * solid_angle - Solid angle subtended by face


    > * angle_normalized - Solid angle normalized such that the

    >     faces with the largest


    > * area - Area of the facet


    > * face_dist - Distance between site n and the facet


    > * volume - Volume of Voronoi cell for this face


    > * n_verts - Number of vertices on the facet




#### get_nn_info(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure
using Voronoi decomposition.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.sites.md#pymatgen.core.sites.Site), array, float))



#### get_voronoi_polyhedra(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), n: int)
Gives a weighted polyhedra around a site.

See ref: A Proposed Rigorous Definition of Coordination Number,
M. O’Keeffe, Acta Cryst. (1979). A35, 772-775


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure for which to evaluate the
    coordination environment.


    * **n** (*int*) – site index.



* **Returns**

    A dict of sites sharing a common Voronoi facet with the site
    n mapped to a directory containing statistics about the facet:

    >
    > * solid_angle - Solid angle subtended by face


    > * angle_normalized - Solid angle normalized such that the

    >     faces with the largest


    > * area - Area of the facet


    > * face_dist - Distance between site n and the facet


    > * volume - Volume of Voronoi cell for this face


    > * n_verts - Number of vertices on the facet




#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### pymatgen.analysis.local_env.get_neighbors_of_site_with_index(struct, n, approach='min_dist', delta=0.1, cutoff=10)
Returns the neighbors of a given site using a specific neighbor-finding
method.


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in Structure object for which motif type
    is to be determined.


    * **approach** (*str*) – type of neighbor-finding approach, where
    “min_dist” will use the MinimumDistanceNN class,
    “voronoi” the VoronoiNN class, “min_OKeeffe” the
    MinimumOKeeffe class, and “min_VIRE” the MinimumVIRENN class.


    * **delta** (*float*) – tolerance involved in neighbor finding.


    * **cutoff** (*float*) – (large) radius to find tentative neighbors.


Returns: neighbor sites.


### pymatgen.analysis.local_env.get_okeeffe_distance_prediction(el1, el2)
Returns an estimate of the bond valence parameter (bond length) using
the derived parameters from ‘Atoms Sizes and Bond Lengths in Molecules
and Crystals’ (O’Keeffe & Brese, 1991). The estimate is based on two
experimental parameters: r and c. The value for r  is based off radius,
while c is (usually) the Allred-Rochow electronegativity. Values used
are *not* generated from pymatgen, and are found in
‘okeeffe_params.json’.


* **Parameters**


    * **el1** ([*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)) – two Element objects


    * **el2** ([*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)) – two Element objects



* **Returns**

    a float value of the predicted bond length



### pymatgen.analysis.local_env.get_okeeffe_params(el_symbol)
Returns the elemental parameters related to atom size and
electronegativity which are used for estimating bond-valence
parameters (bond length) of pairs of atoms on the basis of data
provided in ‘Atoms Sizes and Bond Lengths in Molecules and Crystals’
(O’Keeffe & Brese, 1991).


* **Parameters**

    **el_symbol** (*str*) – element symbol.



* **Returns**

    atom-size (‘r’) and electronegativity-related (‘c’)

        parameter.




* **Return type**

    (dict)



### pymatgen.analysis.local_env.gramschmidt(vin, uin)
Returns that part of the first input vector
that is orthogonal to the second input vector.
The output vector is not normalized.


* **Parameters**


    * **vin** (*numpy array*) – first input vector


    * **uin** (*numpy array*) – second input vector



### pymatgen.analysis.local_env.metal_edge_extender(mol_graph, cutoff: float = 2.5, metals: list | tuple | None = ('Li', 'Mg', 'Ca', 'Zn', 'B', 'Al'), coordinators: list | tuple = ('O', 'N', 'F', 'S', 'Cl'))
Function to identify and add missed coordinate bond edges for metals.


* **Parameters**


    * **mol_graph** – pymatgen.analysis.graphs.MoleculeGraph object


    * **cutoff** – cutoff in Angstrom. Metal-coordinator sites that are closer
    together than this value will be considered coordination bonds.
    If the MoleculeGraph contains a metal, but no coordination bonds are found
    with the chosen cutoff, the cutoff will be increased by 1 Angstrom
    and another attempt will be made to identify coordination bonds.


    * **metals** – Species considered metals for the purpose of identifying
    missed coordinate bond edges. The set {“Li”, “Mg”, “Ca”, “Zn”, “B”, “Al”}
    (default) corresponds to the settings used in the LIBE dataset.
    Alternatively, set to None to cause any Species classified as a metal
    by Specie.is_metal to be considered a metal.


    * **coordinators** – Possible coordinating species to consider when identifying
    missed coordinate bonds. The default set {“O”, “N”, “F”, “S”, “Cl”} was
    used in the LIBE dataset.



* **Returns**

    pymatgen.analysis.graphs.MoleculeGraph object with additional

        metal bonds (if any found) added




* **Return type**

    mol_graph



### pymatgen.analysis.local_env.oxygen_edge_extender(mol_graph: [MoleculeGraph](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph))
Identify and add missed O-C or O-H bonds. This is particularly
important when oxygen is forming three bonds, e.g. in H3O+ or XOH2+.
See [https://github.com/materialsproject/pymatgen/pull/2903](https://github.com/materialsproject/pymatgen/pull/2903) for details.


* **Parameters**

    **mol_graph** ([*MoleculeGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph)) – molecule graph to extend



* **Returns**

    object with additional O-C or O-H bonds added (if any found)



* **Return type**

    [MoleculeGraph](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph)



### pymatgen.analysis.local_env.site_is_of_motif_type(struct, n, approach='min_dist', delta=0.1, cutoff=10, thresh=None)
Returns the motif type of the site with index n in structure struct;
currently featuring “tetrahedral”, “octahedral”, “bcc”, and “cp”
(close-packed: fcc and hcp) as well as “square pyramidal” and
“trigonal bipyramidal”. If the site is not recognized,
“unrecognized” is returned. If a site should be assigned to two
different motifs, “multiple assignments” is returned.


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in Structure object for which motif type
    is to be determined.


    * **approach** (*str*) – type of neighbor-finding approach, where
    “min_dist” will use the MinimumDistanceNN class,
    “voronoi” the VoronoiNN class, “min_OKeeffe” the
    MinimumOKeeffe class, and “min_VIRE” the MinimumVIRENN class.


    * **delta** (*float*) – tolerance involved in neighbor finding.


    * **cutoff** (*float*) – (large) radius to find tentative neighbors.


    * **thresh** (*dict*) – thresholds for motif criteria (currently, required
    keys and their default values are “qtet”: 0.5,
    “qoct”: 0.5, “qbcc”: 0.5, “q6”: 0.4).


Returns: motif type (str).


### pymatgen.analysis.local_env.solid_angle(center, coords)
Helper method to calculate the solid angle of a set of coords from the
center.


* **Parameters**


    * **center** (*3x1 array*) – Center to measure solid angle from.


    * **coords** (*Nx3 array*) – List of coords to determine solid angle.



* **Returns**

    The solid angle.



### pymatgen.analysis.local_env.vol_tetra(vt1, vt2, vt3, vt4)
Calculate the volume of a tetrahedron, given the four vertices of vt1,
vt2, vt3 and vt4.


* **Parameters**


    * **vt1** (*array-like*) – coordinates of vertex 1.


    * **vt2** (*array-like*) – coordinates of vertex 2.


    * **vt3** (*array-like*) – coordinates of vertex 3.


    * **vt4** (*array-like*) – coordinates of vertex 4.



* **Returns**

    volume of the tetrahedron.



* **Return type**

    (float)