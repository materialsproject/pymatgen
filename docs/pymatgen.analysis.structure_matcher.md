---
layout: default
title: pymatgen.analysis.structure_matcher.md
nav_exclude: true
---

# pymatgen.analysis.structure_matcher module

This module provides classes to perform fitting of structures.


### _class_ pymatgen.analysis.structure_matcher.AbstractComparator()
Bases: `MSONable`

Abstract Comparator class. A Comparator defines how sites are compared in
a structure.


#### _abstract_ are_equal(sp1, sp2)
Defines how the species of two sites are considered equal. For
example, one can consider sites to have the same species only when
the species are exactly the same, i.e., Fe2+ matches Fe2+ but not
Fe3+. Or one can define that only the element matters,
and all oxidation state information are ignored.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are considered equal.



#### as_dict()

* **Returns**

    MSONable dict



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    Comparator.



#### _abstract_ get_hash(composition)
Defines a hash to group structures. This allows structures to be
grouped efficiently for comparison. The hash must be invariant under
supercell creation. (e.g. composition is not a good hash, but
fractional_composition might be). Reduced formula is not a good formula,
due to weird behavior with fractional occupancy.

Composition is used here instead of structure because for anonymous
matches it is much quicker to apply a substitution to a composition
object than a structure object.


* **Parameters**

    **composition** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – composition of the structure



* **Returns**

    A hashable object. Examples can be string formulas, integers etc.



### _class_ pymatgen.analysis.structure_matcher.ElementComparator()
Bases: `AbstractComparator`

A Comparator that matches elements. i.e. oxidation states are
ignored.


#### are_equal(sp1, sp2)
True if element:amounts are exactly the same, i.e.,
oxidation state is not considered.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are the same based on element
    and amounts.



#### get_hash(composition)
Returns: Fractional element composition.


### _class_ pymatgen.analysis.structure_matcher.FrameworkComparator()
Bases: `AbstractComparator`

A Comparator that matches sites, regardless of species.


#### are_equal(sp1, sp2)
True if there are atoms on both sites.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    True always



#### get_hash(composition)
No hash possible.


### _class_ pymatgen.analysis.structure_matcher.OccupancyComparator()
Bases: `AbstractComparator`

A Comparator that matches occupancies on sites,
irrespective of the species of those sites.


#### are_equal(sp1, sp2)

* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    True if sets of occupancies (amt) are equal on both sites.



#### get_hash(composition)

* **Parameters**

    **composition** – Composition.



* **Returns**


    1. Difficult to define sensible hash




### _class_ pymatgen.analysis.structure_matcher.OrderDisorderElementComparator()
Bases: `AbstractComparator`

A Comparator that matches sites, given some overlap in the element
composition.


#### are_equal(sp1, sp2)
True if there is some overlap in composition between the species.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    True always



#### get_hash(composition)
Returns: Fractional composition.


### _class_ pymatgen.analysis.structure_matcher.SpeciesComparator()
Bases: `AbstractComparator`

A Comparator that matches species exactly. The default used in StructureMatcher.


#### are_equal(sp1, sp2)
True if species are exactly the same, i.e., Fe2+ == Fe2+ but not Fe3+.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are equal.



#### get_hash(composition)
Returns: Fractional composition.


### _class_ pymatgen.analysis.structure_matcher.SpinComparator()
Bases: `AbstractComparator`

A Comparator that matches magnetic structures to their inverse spins.
This comparator is primarily used to filter magnetically ordered
structures with opposite spins, which are equivalent.


#### are_equal(sp1, sp2)
True if species are exactly the same, i.e., Fe2+ == Fe2+ but not
Fe3+. and the spins are reversed. i.e., spin up maps to spin down,
and vice versa.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are equal.



#### get_hash(composition)
Returns: Fractional composition.


### _class_ pymatgen.analysis.structure_matcher.StructureMatcher(ltol: float = 0.2, stol: float = 0.3, angle_tol: float = 5, primitive_cell: bool = True, scale: bool = True, attempt_supercell: bool = False, allow_subset: bool = False, comparator: AbstractComparator | None = None, supercell_size: Literal['num_sites', 'num_atoms', 'volume'] = 'num_sites', ignored_species: Sequence[SpeciesLike] = ())
Bases: `MSONable`

Class to match structures by similarity.

Algorithm:


1. Given two structures: s1 and s2


2. Optional: Reduce to primitive cells.


3. If the number of sites do not match, return False


4. Reduce to s1 and s2 to Niggli Cells


5. Optional: Scale s1 and s2 to same volume.


6. Optional: Remove oxidation states associated with sites


7. Find all possible lattice vectors for s2 within shell of ltol.


8. For s1, translate an atom in the smallest set to the origin


9. For s2: find all valid lattices from permutations of the list
of lattice vectors (invalid if: det(Lattice Matrix) < half
volume of original s2 lattice)


10. For each valid lattice:


    1. If the lattice angles of are within tolerance of s1,
basis change s2 into new lattice.


    2. For each atom in the smallest set of s2:

> i. Translate to origin and compare fractional sites in
> structure within a fractional tolerance.
> ii. If true:

> > ia. Convert both lattices to Cartesian and place
> > both structures on an average lattice
> > ib. Compute and return the average and max rms
> > displacement between the two structures normalized
> > by the average free length per atom

> > if fit function called:

> >     if normalized max rms displacement is less than
> >     stol. Return True

> > if get_rms_dist function called:

> >     if normalized average rms displacement is less
> >     than the stored rms displacement, store and
> >     continue. (This function will search all possible
> >     lattices for the smallest average rms displacement
> >     between the two structures)


* **Parameters**


    * **ltol** (*float*) – Fractional length tolerance. Default is 0.2.


    * **stol** (*float*) – Site tolerance. Defined as the fraction of the
    average free length per atom := ( V / Nsites ) \*\* (1/3)
    Default is 0.3.


    * **angle_tol** (*float*) – Angle tolerance in degrees. Default is 5 degrees.


    * **primitive_cell** (*bool*) – If true: input structures will be reduced to
    primitive cells prior to matching. Default to True.


    * **scale** (*bool*) – Input structures are scaled to equivalent volume if
    true; For exact matching, set to False.


    * **attempt_supercell** (*bool*) – If set to True and number of sites in
    cells differ after a primitive cell reduction (divisible by an
    integer) attempts to generate a supercell transformation of the
    smaller cell which is equivalent to the larger structure.


    * **allow_subset** (*bool*) – Allow one structure to match to the subset of
    another structure. Eg. Matching of an ordered structure onto a
    disordered one, or matching a delithiated to a lithiated
    structure. This option cannot be combined with
    attempt_supercell, or with structure grouping.


    * **comparator** (*Comparator*) – A comparator object implementing an equals
    method that declares equivalency of sites. Default is
    SpeciesComparator, which implies rigid species
    mapping, i.e., Fe2+ only matches Fe2+ and not Fe3+.

    Other comparators are provided, e.g., ElementComparator which
    matches only the elements and not the species.

    The reason why a comparator object is used instead of
    supplying a comparison function is that it is not possible to
    pickle a function, which makes it otherwise difficult to use
    StructureMatcher with Python’s multiprocessing.



    * **supercell_size** (*str** or **list*) – Method to use for determining the
    size of a supercell (if applicable). Possible values are
    ‘num_sites’, ‘num_atoms’, ‘volume’, or an element or list of elements
    present in both structures.


    * **ignored_species** (*list*) – A list of ions to be ignored in matching.
    Useful for matching structures that have similar frameworks
    except for certain ions, e.g., Li-ion intercalation frameworks.
    This is more useful than allow_subset because it allows better
    control over what species are ignored in the matching.



#### as_dict()

* **Returns**

    MSONable dict



#### fit(struct1: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symmetric: bool = False, skip_structure_reduction: bool = False)
Fit two structures.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 2nd structure


    * **symmetric** (*bool*) – Defaults to False
    If True, check the equality both ways.
    This only impacts a small percentage of structures


    * **skip_structure_reduction** (*bool*) – Defaults to False
    If True, skip to get a primitive structure and perform Niggli reduction for struct1 and struct2



* **Returns**

    True if the structures are equivalent



* **Return type**

    bool



#### fit_anonymous(struct1: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), niggli: bool = True, skip_structure_reduction: bool = False)
Performs an anonymous fitting, which allows distinct species in one structure to map
to another. E.g., to compare if the Li2O and Na2O structures are similar.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 2nd structure


    * **niggli** (*bool*) – If true, perform Niggli reduction for struct1 and struct2


    * **skip_structure_reduction** (*bool*) – Defaults to False
    If True, skip to get a primitive structure and perform Niggli reduction for struct1 and struct2



* **Returns**

    Whether a species mapping can map struct1 to stuct2



* **Return type**

    True/False



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    StructureMatcher



#### get_all_anonymous_mappings(struct1, struct2, niggli=True, include_dist=False)
Performs an anonymous fitting, which allows distinct species in one
structure to map to another. Returns a dictionary of species
substitutions that are within tolerance.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 2nd structure


    * **niggli** (*bool*) – Find niggli cell in preprocessing


    * **include_dist** (*bool*) – Return the maximin distance with each mapping



* **Returns**

    list of species mappings that map struct1 to struct2.



#### get_best_electronegativity_anonymous_mapping(struct1: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Performs an anonymous fitting, which allows distinct species in one
structure to map to another. E.g., to compare if the Li2O and Na2O
structures are similar. If multiple substitutions are within tolerance
this will return the one which minimizes the difference in
electronegativity between the matches species.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 2nd structure



* **Returns**

    Mapping of struct1 species to struct2 species



* **Return type**

    min_mapping (dict)



#### get_mapping(superset, subset)
Calculate the mapping from superset to subset.


* **Parameters**


    * **superset** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure containing at least the sites in
    subset (within the structure matching tolerance)


    * **subset** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure containing some of the sites in
    superset (within the structure matching tolerance)



* **Returns**

    numpy array such that superset.sites[mapping] is within matching
    tolerance of subset.sites or None if no such mapping is possible



#### get_rms_anonymous(struct1, struct2)
Performs an anonymous fitting, which allows distinct species in one
structure to map to another. E.g., to compare if the Li2O and Na2O
structures are similar.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 2nd structure



* **Returns**

    (min_rms, min_mapping)
    min_rms is the minimum rms distance, and min_mapping is the
    corresponding minimal species mapping that would map
    struct1 to struct2. (None, None) is returned if the minimax_rms
    exceeds the threshold.



#### get_rms_dist(struct1, struct2)
Calculate RMS displacement between two structures.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – 2nd structure



* **Returns**

    rms displacement normalized by (Vol / nsites) \*\* (1/3)
    and maximum distance between paired sites. If no matching
    lattice is found None is returned.



#### get_s2_like_s1(struct1, struct2, include_ignored_species=True)
Performs transformations on struct2 to put it in a basis similar to
struct1 (without changing any of the inter-site distances).


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Reference structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to transform.


    * **include_ignored_species** (*bool*) – Defaults to True,
    the ignored_species is also transformed to the struct1
    lattice orientation, though obviously there is no direct
    matching to existing sites.



* **Returns**

    A structure object similar to struct1, obtained by making a
    supercell, sorting, and translating struct2.



#### get_supercell_matrix(supercell, struct)
Returns the matrix for transforming struct to supercell. This
can be used for very distorted ‘supercells’ where the primitive cell
is impossible to find.


#### get_transformation(struct1, struct2)
Returns the supercell transformation, fractional translation vector,
and a mapping to transform struct2 to be similar to struct1.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Reference structure


    * **struct2** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure to transform.



* **Returns**

    supercell matrix
    vector (numpy.ndarray(3)): fractional translation vector
    mapping (list(int or None)):

    > The first len(struct1) items of the mapping vector are the
    > indices of struct1’s corresponding sites in struct2 (or None
    > if there is no corresponding site), and the other items are
    > the remaining site indices of struct2.




* **Return type**

    supercell (numpy.ndarray(3, 3))



#### group_structures(s_list, anonymous=False)
Given a list of structures, use fit to group
them by structural equality.


* **Parameters**


    * **s_list** (*[*[*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)*]*) – List of structures to be grouped


    * **anonymous** (*bool*) – Whether to use anonymous mode.



* **Returns**

    A list of lists of matched structures
    Assumption: if s1 == s2 but s1 != s3, than s2 and s3 will be put
    in different groups without comparison.