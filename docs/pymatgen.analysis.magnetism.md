---
layout: default
title: pymatgen.analysis.magnetism.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.magnetism package

Package for analysis of magnetic structures.


## pymatgen.analysis.magnetism.analyzer module

This module provides some useful functions for dealing with magnetic Structures
(e.g. Structures with associated magmom tags).


### _class_ CollinearMagneticStructureAnalyzer(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), overwrite_magmom_mode: OverwriteMagmomMode | str = 'none', round_magmoms: bool = False, detect_valences: bool = False, make_primitive: bool = True, default_magmoms: dict | None = None, set_net_positive: bool = True, threshold: float = 0, threshold_nonmag: float = 0.1)
Bases: `object`

A class which provides a few helpful methods to analyze
collinear magnetic structures.

If magnetic moments are not defined, moments will be
taken either from default_magmoms.yaml (similar to the
default magmoms in MPRelaxSet, with a few extra definitions)
or from a specie:magmom dict provided by the default_magmoms
kwarg.

Input magmoms can be replaced using the ‘overwrite_magmom_mode’
kwarg. This can be:
\* “none” to do nothing,
\* “respect_sign” which will overwrite existing magmoms with

> those from default_magmoms but will keep sites with positive magmoms
> positive, negative magmoms negative and zero magmoms zero,


* “respect_zeros”, which will give a ferromagnetic structure
(all positive magmoms from default_magmoms) but still keep sites with
zero magmoms as zero,


* “replace_all” which will try to guess initial magmoms for
all sites in the structure irrespective of input structure
(this is most suitable for an initial DFT calculation),


* “replace_all_if_undefined” is the same as “replace_all” but only if
no magmoms are defined in input structure, otherwise it will respect
existing magmoms.


* “normalize” will normalize magmoms to unity, but will respect sign
(used for comparing orderings), magmoms < threshold will be set to zero


* **Parameters**


    * **structure** – input Structure object


    * **overwrite_magmom_mode** – “respect_sign”, “respect_zeros”, “replace_all”,
    “replace_all_if_undefined”, “normalize” (default “none”)


    * **round_magmoms** – will round input magmoms to
    specified number of decimal places if integer is supplied, if set
    to a float will try and group magmoms together using a kernel density
    estimator of provided width, and extracting peaks of the estimator
    detect_valences: if True, will attempt to assign valences
    to input structure


    * **make_primitive** – if True, will transform to primitive
    magnetic cell


    * **default_magmoms** – (optional) dict specifying default magmoms


    * **set_net_positive** – if True, will change sign of magnetic
    moments such that the net magnetization is positive. Argument will be
    ignored if mode “respect_sign” is used.


    * **threshold** – number (in Bohr magnetons) below which magmoms
    will be rounded to zero


    * **threshold_nonmag** – number (in Bohr magneton)
    below which nonmagnetic ions (with no magmom specified
    in default_magmoms) will be rounded to zero



#### _static_ _round_magmoms(magmoms: ArrayLike, round_magmoms_mode: float)
If round_magmoms_mode is an integer, simply round to that number
of decimal places, else if set to a float will try and round
intelligently by grouping magmoms.


#### get_exchange_group_info(symprec: float = 0.01, angle_tolerance: float = 5)
Returns the information on the symmetry of the Hamiltonian
describing the exchange energy of the system, taking into
account relative direction of magnetic moments but not their
absolute direction.

This is not strictly accurate (e.g. some/many atoms will
have zero magnetic moments), but defining symmetry this
way is a useful way of keeping track of distinct magnetic
orderings within pymatgen.


* **Parameters**


    * **symprec** – same as SpacegroupAnalyzer (Default value = 1e-2)


    * **angle_tolerance** – same as SpacegroupAnalyzer (Default value = 5)



* **Returns**

    spacegroup_symbol, international_number



#### get_ferromagnetic_structure(make_primitive: bool = True)
Returns a Structure with all magnetic moments positive
or zero.


* **Parameters**

    **make_primitive** – Whether to make structure primitive after
    making all magnetic moments positive (Default value = True)



* **Returns**

    Structure



#### get_nonmagnetic_structure(make_primitive: bool = True)
Returns a Structure without magnetic moments defined.


* **Parameters**

    **make_primitive** – Whether to make structure primitive after
    removing magnetic information (Default value = True)



* **Returns**

    Structure



#### get_structure_with_only_magnetic_atoms(make_primitive: bool = True)
Returns a Structure with only magnetic atoms present.


* **Parameters**

    **make_primitive** – Whether to make structure primitive after
    removing non-magnetic atoms (Default value = True)



* **Returns**

    Structure



#### get_structure_with_spin()
Returns a Structure with species decorated with spin values instead
of using magmom site properties.


#### _property_ is_magnetic(_: boo_ )
Convenience property, returns True if any non-zero magmoms present.


#### _property_ magmoms(_: ndarra_ )
Convenience property, returns magmoms as a numpy array.


#### _property_ magnetic_species_and_magmoms(_: dict[str, Any_ )
Returns a dict of magnetic species and the magnitude of
their associated magmoms. Will return a list if there are
multiple magmoms per species.


* **Returns**

    dict of magnetic species and magmoms



#### matches_ordering(other: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Compares the magnetic orderings of one structure with another.


* **Parameters**

    **other** – Structure to compare



* **Returns**

    True if magnetic orderings match, False otherwise



* **Return type**

    bool



#### _property_ number_of_magnetic_sites(_: in_ )
Number of magnetic sites present in structure.


#### number_of_unique_magnetic_sites(symprec: float = 0.001, angle_tolerance: float = 5)

* **Parameters**


    * **symprec** – same as in SpacegroupAnalyzer (Default value = 1e-3)


    * **angle_tolerance** – same as in SpacegroupAnalyzer (Default value = 5).



* **Returns**

    Number of symmetrically-distinct magnetic sites present in structure.



* **Return type**

    int



#### _property_ ordering(_: Orderin_ )
Applies heuristics to return a magnetic ordering for a collinear
magnetic structure. Result is not guaranteed for correctness.


* **Returns**

    Enum (‘FiM’ is used as the abbreviation for ferrimagnetic)



* **Return type**

    Ordering



#### _property_ types_of_magnetic_specie(_: tuple[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.md#pymatgen.core.periodic_table.DummySpecies), ..._ )
Specie->Species rename. Used to maintain backwards compatibility.


#### _property_ types_of_magnetic_species(_: tuple[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.md#pymatgen.core.periodic_table.DummySpecies), ..._ )
Equivalent to Structure.types_of_specie but only returns
magnetic species.


* **Returns**

    types of Species



* **Return type**

    tuple



### _class_ MagneticDeformation(type, deformation)
Bases: `tuple`

Create new instance of MagneticDeformation(type, deformation)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('type', 'deformation'_ )

#### _classmethod_ _make(iterable)
Make a new MagneticDeformation object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new MagneticDeformation object replacing specified fields with new values


#### deformation()
Alias for field number 1


#### type()
Alias for field number 0


### _class_ MagneticStructureEnumerator(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), default_magmoms: dict[str, float] | None = None, strategies: list[str] | tuple[str, ...] = ('ferromagnetic', 'antiferromagnetic'), automatic: bool = True, truncate_by_symmetry: bool = True, transformation_kwargs: dict | None = None)
Bases: `object`

Combines MagneticStructureAnalyzer and MagOrderingTransformation to
automatically generate a set of transformations for a given structure
and produce a list of plausible magnetic orderings.

This class will try generated different collinear
magnetic orderings for a given input structure.

If the input structure has magnetic moments defined, it
is possible to use these as a hint as to which elements are
magnetic, otherwise magnetic elements will be guessed
(this can be changed using default_magmoms kwarg).


* **Parameters**


    * **structure** – input structure


    * **default_magmoms** – (optional, defaults provided) dict of
    magnetic elements to their initial magnetic moments in µB, generally
    these are chosen to be high-spin since they can relax to a low-spin
    configuration during a DFT electronic configuration


    * **strategies** – different ordering strategies to use, choose from:
    ferromagnetic, antiferromagnetic, antiferromagnetic_by_motif,
    ferrimagnetic_by_motif and ferrimagnetic_by_species (here, “motif”,
    means to use a different ordering parameter for symmetry inequivalent
    sites)


    * **automatic** – if True, will automatically choose sensible strategies


    * **truncate_by_symmetry** – if True, will remove very unsymmetrical
    orderings that are likely physically implausible


    * **transformation_kwargs** – keyword arguments to pass to
    MagOrderingTransformation, to change automatic cell size limits, etc.



#### _generate_ordered_structures(sanitized_input_structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), transformations: dict[str, [MagOrderingTransformation](pymatgen.transformations.md#pymatgen.transformations.advanced_transformations.MagOrderingTransformation)])
Apply our input structure to our list of transformations and output a list
of ordered structures that have been pruned for duplicates and for those
with low symmetry (optional). Sets self.ordered_structures
and self.ordered_structures_origins instance variables.


* **Parameters**


    * **sanitized_input_structure** – A sanitized input structure


    * **(****_sanitize_input_structure****)** –


    * **transformations** – A dict of transformations (values) and name of


    * **strategy** (*enumeration*) –


    * **keeping** (*for record*) –



* **Returns**

    list[Structures]



#### _generate_transformations(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
The central problem with trying to enumerate magnetic orderings is
that we have to enumerate orderings that might plausibly be magnetic
ground states, while not enumerating orderings that are physically
implausible. The problem is that it is not always obvious by e.g.
symmetry arguments alone which orderings to prefer. Here, we use a
variety of strategies (heuristics) to enumerate plausible orderings,
and later discard any duplicates that might be found by multiple
strategies. This approach is not ideal, but has been found to be
relatively robust over a wide range of magnetic structures.


* **Parameters**

    **structure** – A sanitized input structure (_sanitize_input_structure)



* **Returns**

    A dict of a transformation class instance (values) and name of enumeration strategy (keys)



* **Return type**

    dict



#### _static_ _sanitize_input_structure(struct: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Sanitize our input structure by removing magnetic information
and making primitive.


* **Parameters**

    **struct** – Structure



* **Returns**

    Structure



#### available_strategies(_ = ('ferromagnetic', 'antiferromagnetic', 'ferrimagnetic_by_motif', 'ferrimagnetic_by_species', 'antiferromagnetic_by_motif', 'nonmagnetic'_ )

### _class_ Ordering(value)
Bases: `Enum`

Enumeration defining possible magnetic orderings.


#### AFM(_ = 'AFM_ )

#### FM(_ = 'FM_ )

#### FiM(_ = 'FiM_ )

#### NM(_ = 'NM_ )

#### Unknown(_ = 'Unknown_ )

### _class_ OverwriteMagmomMode(value)
Bases: `Enum`

Enumeration defining different modes for analyzer.


#### none(_ = 'none_ )

#### normalize(_ = 'normalize_ )

#### replace_all(_ = 'replace_all_ )

#### respect_sign(_ = 'respect_sign_ )

#### respect_zero(_ = 'respect_zeros_ )

### magnetic_deformation(structure_A: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), structure_B: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Calculates ‘magnetic deformation proxy’,
a measure of deformation (norm of finite strain)
between ‘non-magnetic’ (non-spin-polarized) and
ferromagnetic structures.

Adapted from Bocarsly et al. 2017,
doi: 10.1021/acs.chemmater.6b04729


* **Parameters**


    * **structure_A** – Structure


    * **structure_B** – Structure



* **Returns**

    MagneticDeformation


## pymatgen.analysis.magnetism.heisenberg module

This module implements a simple algorithm for extracting nearest neighbor
exchange parameters by mapping low energy magnetic orderings to a Heisenberg
model.


### _class_ HeisenbergMapper(ordered_structures, energies, cutoff=0, tol: float = 0.02)
Bases: `object`

Class to compute exchange parameters from low energy magnetic orderings.

Exchange parameters are computed by mapping to a classical Heisenberg
model. Strategy is the scheme for generating neighbors. Currently only
MinimumDistanceNN is implemented.
n+1 unique orderings are required to compute n exchange
parameters.

First run a MagneticOrderingsWF to obtain low energy collinear magnetic
orderings and find the magnetic ground state. Then enumerate magnetic
states with the ground state as the input structure, find the subset
of supercells that map to the ground state, and do static calculations
for these orderings.


* **Parameters**


    * **ordered_structures** (*list*) – Structure objects with magmoms.


    * **energies** (*list*) – Total energies of each relaxed magnetic structure.


    * **cutoff** (*float*) – Cutoff in Angstrom for nearest neighbor search.
    Defaults to 0 (only NN, no NNN, etc.)


    * **tol** (*float*) – Tolerance (in Angstrom) on nearest neighbor distances
    being equal.


    * **strategy** (*object*) – Class from pymatgen.analysis.local_env for constructing graphs.


    * **sgraphs** (*list*) – StructureGraph objects.


    * **unique_site_ids** (*dict*) – Maps each site to its unique numerical identifier.


    * **wyckoff_ids** (*dict*) – Maps unique numerical identifier to wyckoff position.


    * **nn_interactions** (*dict*) – {i: j} pairs of NN interactions between unique sites.


    * **dists** (*dict*) – NN, NNN, and NNNN interaction distances


    * **ex_mat** (*DataFrame*) – Invertible Heisenberg Hamiltonian for each graph.


    * **ex_params** (*dict*) – Exchange parameter values (meV/atom)



#### _get_exchange_df()
Loop over all sites in a graph and count the number and types of
nearest neighbor interactions, computing +-

```
|S_i . S_j|
```

 to construct
a Heisenberg Hamiltonian for each graph. Sets self.ex_mat instance variable.

TODO Deal with large variance in

```
|S|
```

 across configs


#### _static_ _get_graphs(cutoff, ordered_structures)
Generate graph representations of magnetic structures with nearest
neighbor bonds. Right now this only works for MinimumDistanceNN.


* **Parameters**


    * **cutoff** (*float*) – Cutoff in Angstrom for nearest neighbor search.


    * **ordered_structures** (*list*) – Structure objects.



* **Returns**

    StructureGraph objects.



* **Return type**

    sgraphs (list)



#### _get_j_exc(i, j, dist)
Convenience method for looking up exchange parameter between two sites.


* **Parameters**


    * **i** (*int*) – index of ith site


    * **j** (*int*) – index of jth site


    * **dist** (*float*) – distance (Angstrom) between sites
    (10E-2 precision)



* **Returns**

    Exchange parameter in meV



* **Return type**

    j_exc (float)



#### _get_nn_dict()
Sets self.nn_interactions and self.dists instance variables describing unique
nearest neighbor interactions.


#### _static_ _get_unique_sites(structure)
Get dict that maps site indices to unique identifiers.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – ground state Structure object.



* **Returns**

    maps tuples of equivalent site indices to a

        unique int identifier

    wyckoff_ids (dict): maps tuples of equivalent site indices to their

        wyckoff symbols




* **Return type**

    unique_site_ids (dict)



#### estimate_exchange(fm_struct=None, afm_struct=None, fm_e=None, afm_e=None)
Estimate <J> for a structure based on low energy FM and AFM orderings.


* **Parameters**


    * **fm_struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – fm structure with ‘magmom’ site property


    * **afm_struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – afm structure with ‘magmom’ site property


    * **fm_e** (*float*) – fm energy/atom


    * **afm_e** (*float*) – afm energy/atom



* **Returns**

    Average exchange parameter (meV/atom)



* **Return type**

    j_avg (float)



#### get_exchange()
Take Heisenberg Hamiltonian and corresponding energy for each row and
solve for the exchange parameters.


* **Returns**

    Exchange parameter values (meV/atom).



* **Return type**

    ex_params (dict)



#### get_heisenberg_model()
Save results of mapping to a HeisenbergModel object.


* **Returns**

    MSONable object.



* **Return type**

    HeisenbergModel



#### get_interaction_graph(filename=None)
Get a StructureGraph with edges and weights that correspond to exchange
interactions and J_ij values, respectively.


* **Parameters**

    **filename** (*str*) – if not None, save interaction graph to filename.



* **Returns**

    Exchange interaction graph.



* **Return type**

    igraph ([StructureGraph](pymatgen.analysis.md#pymatgen.analysis.graphs.StructureGraph))



#### get_low_energy_orderings()
Find lowest energy FM and AFM orderings to compute E_AFM - E_FM.


* **Returns**

    fm structure with ‘magmom’ site property
    afm_struct (Structure): afm structure with ‘magmom’ site property
    fm_e (float): fm energy
    afm_e (float): afm energy



* **Return type**

    fm_struct ([Structure](pymatgen.core.md#pymatgen.core.structure.Structure))



#### get_mft_temperature(j_avg)
Crude mean field estimate of critical temperature based on <J> for
one sublattice, or solving the coupled equations for a multisublattice
material.


* **Parameters**

    **j_avg** (*float*) – j_avg (float): Average exchange parameter (meV/atom)



* **Returns**

    Critical temperature (K)



* **Return type**

    mft_t (float)



### _class_ HeisenbergModel(formula=None, structures=None, energies=None, cutoff=None, tol=None, sgraphs=None, unique_site_ids=None, wyckoff_ids=None, nn_interactions=None, dists=None, ex_mat=None, ex_params=None, javg=None, igraph=None)
Bases: `MSONable`

Store a Heisenberg model fit to low-energy magnetic orderings.
Intended to be generated by HeisenbergMapper.get_heisenberg_model().


* **Parameters**


    * **formula** (*str*) – Reduced formula of compound.


    * **structures** (*list*) – Structure objects with magmoms.


    * **energies** (*list*) – Energies of each relaxed magnetic structure.


    * **cutoff** (*float*) – Cutoff in Angstrom for nearest neighbor search.


    * **tol** (*float*) – Tolerance (in Angstrom) on nearest neighbor distances being equal.


    * **sgraphs** (*list*) – StructureGraph objects.


    * **unique_site_ids** (*dict*) – Maps each site to its unique numerical
    identifier.


    * **wyckoff_ids** (*dict*) – Maps unique numerical identifier to wyckoff
    position.


    * **nn_interactions** (*dict*) – {i: j} pairs of NN interactions
    between unique sites.


    * **dists** (*dict*) – NN, NNN, and NNNN interaction distances


    * **ex_mat** (*DataFrame*) – Invertible Heisenberg Hamiltonian for each
    graph.


    * **ex_params** (*dict*) – Exchange parameter values (meV/atom).


    * **javg** (*float*) – <J> exchange param (meV/atom).


    * **igraph** ([*StructureGraph*](pymatgen.analysis.md#pymatgen.analysis.graphs.StructureGraph)) – Exchange interaction graph.



#### _get_j_exc(i, j, dist)
Convenience method for looking up exchange parameter between two sites.


* **Parameters**


    * **i** (*int*) – index of ith site


    * **j** (*int*) – index of jth site


    * **dist** (*float*) – distance (Angstrom) between sites +- tol



* **Returns**

    Exchange parameter in meV



* **Return type**

    j_exc (float)



#### as_dict()
Because some dicts have tuple keys, some sanitization is required for json compatibility.


#### _classmethod_ from_dict(d)
Create a HeisenbergModel from a dict.


### _class_ HeisenbergScreener(structures, energies, screen=False)
Bases: `object`

Class to clean and screen magnetic orderings.

This class pre-processes magnetic orderings and energies for
HeisenbergMapper. It prioritizes low-energy orderings with large and
localized magnetic moments.


* **Parameters**


    * **structures** (*list*) – Structure objects with magnetic moments.


    * **energies** (*list*) – Energies/atom of magnetic orderings.


    * **screen** (*bool*) – Try to screen out high energy and low-spin configurations.



#### screened_structures()
Sorted structures.


* **Type**

    list



#### screened_energies()
Sorted energies.


* **Type**

    list



#### _static_ _do_cleanup(structures, energies)
Sanitize input structures and energies.

Takes magnetic structures and performs the following operations
- Erases nonmagnetic ions and gives all ions [‘magmom’] site prop
- Converts total energies -> energy / magnetic ion
- Checks for duplicate/degenerate orderings
- Sorts by energy


* **Parameters**


    * **structures** (*list*) – Structure objects with magmoms.


    * **energies** (*list*) – Corresponding energies.



* **Returns**

    Sanitized structures.
    ordered_energies (list): Sorted energies.



* **Return type**

    ordered_structures (list)



#### _static_ _do_screen(structures, energies)
Screen and sort magnetic orderings based on some criteria.

Prioritize low energy orderings and large, localized magmoms. do_clean should be run first to sanitize inputs.


* **Parameters**


    * **structures** (*list*) – At least three structure objects.


    * **energies** (*list*) – Energies.



* **Returns**

    Sorted structures.
    screened_energies (list): Sorted energies.



* **Return type**

    screened_structures (list)


## pymatgen.analysis.magnetism.jahnteller module

JahnTeller distortion analysis.


### _class_ JahnTellerAnalyzer()
Bases: `object`

Will attempt to classify if structure *may* be Jahn-Teller active.
Class currently uses datafile of hard-coded common Jahn-Teller
active ions.
If structure is annotated with magnetic moments, will estimate
if structure may be high-spin or low-spin.
Class aims for more false-positives than false-negatives.

Init for JahnTellerAnalyzer.


#### _static_ _estimate_spin_state(species: str | [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species), motif: Literal['oct', 'tet'], known_magmom: float)
Simple heuristic to estimate spin state. If magnetic moment
is sufficiently close to that predicted for a given spin state,
we assign it that state. If we only have data for one spin
state then that’s the one we use (e.g. we assume all tetrahedral
complexes are high-spin, since this is typically the case).


* **Parameters**


    * **species** – str or Species


    * **motif** (*"oct"** | **"tet"*) – Tetrahedron or octahedron crystal site coordination


    * **known_magmom** – magnetic moment in Bohr magnetons



* **Returns**

    “undefined” (if only one spin state possible), “low”, “high” or “unknown”



#### _static_ _get_number_of_d_electrons(species: [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species))
Get number of d electrons of a species.


* **Parameters**

    **species** – Species object



* **Returns**

    Number of d electrons.



* **Return type**

    int



#### get_analysis(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Convenience method, uses get_analysis_and_structure method.

Obtain an analysis of a given structure and if it may be Jahn-Teller
active or not. This is a heuristic, and may give false positives and
false negatives (false positives are preferred).


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    analysis of structure, with key ‘strength’ which may be ‘none’, ‘strong’,
    ‘weak’, or ‘unknown’ (Default value = 0.1)



#### get_analysis_and_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Obtain an analysis of a given structure and if it may be Jahn-Teller
active or not. This is a heuristic, and may give false positives and
false negatives (false positives are preferred).


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    analysis of structure, with key ‘strength’ which may be ‘none’, ‘strong’,
    ‘weak’, or ‘unknown’ (Default value = 0.1) and decorated structure



#### get_magnitude_of_effect_from_species(species: str | [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species), spin_state: str, motif: str)
Get magnitude of Jahn-Teller effect from provided species, spin state and motif.


* **Parameters**


    * **species** – e.g. Fe2+


    * **spin_state** – “high” or “low”


    * **motif** – “oct” or “tet”



* **Returns**

    “none”, “weak” or “strong”



* **Return type**

    str



#### _static_ get_magnitude_of_effect_from_spin_config(motif: str, spin_config: dict[str, float])
Roughly, the magnitude of Jahn-Teller distortion will be:
\* in octahedral environments, strong if e_g orbitals
unevenly occupied but weak if t_2g orbitals unevenly
occupied
\* in tetrahedral environments always weaker.


* **Parameters**


    * **motif** – “oct” or “tet”


    * **spin_config** – dict of ‘e’ (e_g) and ‘t’ (t2_g) with number of electrons in each state



* **Returns**

    “none”, “weak” or “strong”



* **Return type**

    str



#### is_jahn_teller_active(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Convenience method, uses get_analysis_and_structure method.
Check if a given structure and if it may be Jahn-Teller
active or not. This is a heuristic, and may give false positives and
false negatives (false positives are preferred).


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    boolean, True if might be Jahn-Teller active, False if not



#### _static_ mu_so(species: str | [Species](pymatgen.core.md#pymatgen.core.periodic_table.Species), motif: Literal['oct', 'tet'], spin_state: Literal['high', 'low'])
Calculates the spin-only magnetic moment for a
given species. Only supports transition metals.


* **Parameters**


    * **species** – Species


    * **motif** (*"oct"** | **"tet"*) – Tetrahedron or octahedron crystal site coordination


    * **spin_state** (*"low"** | **"high"*) – Whether the species is in a high or low spin state



* **Returns**

    Spin-only magnetic moment in Bohr magnetons or None if

        species crystal field not defined




* **Return type**

    float



#### tag_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Convenience method, uses get_analysis_and_structure method.
Add a “possible_jt_active” site property on Structure.


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    Decorated Structure, will be in primitive setting.