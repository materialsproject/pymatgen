---
layout: default
title: pymatgen.core.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.core package

This package contains core modules and classes for representing structures and operations on them.


### _load_pmg_settings()

## pymatgen.core.bonds module

This class implements definitions for various kinds of bonds. Typically used in
Molecule analysis.


### _class_ CovalentBond(site1: Site, site2: Site)
Bases: `object`

Defines a covalent bond between two sites.

Initializes a covalent bond between two sites.


* **Parameters**


    * **site1** (*Site*) – First site.


    * **site2** (*Site*) – Second site.



#### get_bond_order(tol: float = 0.2, default_bl: float | None = None)
The bond order according the distance between the two sites.


* **Parameters**


    * **tol** (*float*) – Relative tolerance to test.
    (1 + tol) \* the longest bond distance is considered
    to be the threshold length for a bond to exist.
    (1 - tol) \* the shortest bond distance is considered
    to be the shortest possible bond length
    Defaults to 0.2.


    * **default_bl** – If a particular type of bond does not exist,
    use this bond length as a default value
    (bond order = 1). If None, a ValueError will be thrown.



* **Returns**

    Float value of bond order. For example, for C-C bond in
    benzene, return 1.7.



#### _static_ is_bonded(site1, site2, tol: float = 0.2, bond_order: float | None = None, default_bl: float | None = None)
Test if two sites are bonded, up to a certain limit.


* **Parameters**


    * **site1** (*Site*) – First site


    * **site2** (*Site*) – Second site


    * **tol** (*float*) – Relative tolerance to test. Basically, the code
    checks if the distance between the sites is less than (1 +
    tol) \* typical bond distances. Defaults to 0.2, i.e.,
    20% longer.


    * **bond_order** – Bond order to test. If None, the code simply checks
    against all possible bond data. Defaults to None.


    * **default_bl** – If a particular type of bond does not exist, use this
    bond length. If None, a ValueError will be thrown.



* **Returns**

    Boolean indicating whether two sites are bonded.



#### _property_ length(_: floa_ )
Length of the bond.


### _load_bond_length_data()
Loads bond length data from json file.


### get_bond_length(sp1: SpeciesLike, sp2: SpeciesLike, bond_order: float = 1)
Get the bond length between two species.


* **Parameters**


    * **sp1** (*Species*) – First specie.


    * **sp2** (*Species*) – Second specie.


    * **bond_order** – For species with different possible bond orders,
    this allows one to obtain the bond length for a particular bond
    order. For example, to get the C=C bond length instead of the
    C-C bond length, this should be set to 2. Defaults to 1.



* **Returns**

    Bond length in Angstrom. If no data is available, the sum of the atomic
    radius is used.



### get_bond_order(sp1, sp2, dist: float, tol: float = 0.2, default_bl: float | None = None)
Calculate the bond order given the distance of 2 species.


* **Parameters**


    * **sp1** (*Species*) – First specie.


    * **sp2** (*Species*) – Second specie.


    * **dist** – Their distance in angstrom


    * **tol** (*float*) – Relative tolerance to test. Basically, the code
    checks if the distance between the sites is larger than
    (1 + tol) \* the longest bond distance or smaller than
    (1 - tol) \* the shortest bond distance to determine if
    they are bonded or the distance is too short.
    Defaults to 0.2.


    * **default_bl** – If a particular type of bond does not exist, use this
    bond length (bond order = 1). If None, a ValueError will be thrown.



* **Returns**

    Float value of bond order. For example, for C-C bond in benzene,
    return 1.7.



### obtain_all_bond_lengths(sp1, sp2, default_bl: float | None = None)
Obtain bond lengths for all bond orders from bond length database.


* **Parameters**


    * **sp1** (*Species*) – First specie.


    * **sp2** (*Species*) – Second specie.


    * **default_bl** – If a particular type of bond does not exist, use this
    bond length as a default value (bond order = 1).
    If None, a ValueError will be thrown.



* **Returns**

    A dict mapping bond order to bond length in angstrom


## pymatgen.core.composition module

This module implements a Composition class to represent compositions,
and a ChemicalPotential class to represent potentials.


### _class_ ChemicalPotential(\*args, \*\*kwargs)
Bases: `dict`, `MSONable`

Class to represent set of chemical potentials. Can be: multiplied/divided by a Number
multiplied by a Composition (returns an energy) added/subtracted with other ChemicalPotentials.


* **Parameters**


    * **\*args** – any valid dict init arguments


    * **\*\*kwargs** – any valid dict init arguments.



#### get_energy(composition: Composition, strict: bool = True)
Calculates the energy of a composition.


* **Parameters**


    * **composition** (*Composition*) – input composition


    * **strict** (*bool*) – Whether all potentials must be specified



### _class_ Composition(\*args, strict: bool = False, \*\*kwargs)
Bases: `Hashable`, `Mapping`, `MSONable`, [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)

Represents a Composition, which is essentially a {element:amount} mapping
type. Composition is written to be immutable and hashable,
unlike a standard Python dict.

Note that the key can be either an Element or a Species. Elements and Species
are treated differently. i.e., a Fe2+ is not the same as a Fe3+ Species and
would be put in separate keys. This differentiation is deliberate to
support using Composition to determine the fraction of a particular Species.

Works almost completely like a standard python dictionary, except that
__getitem__ is overridden to return 0 when an element is not found.
(somewhat like a defaultdict, except it is immutable).

Also adds more convenience methods relevant to compositions, e.g.,
get_fraction.

It should also be noted that many Composition related functionality takes
in a standard string as a convenient input. For example,
even though the internal representation of a Fe2O3 composition is
{Element(“Fe”): 2, Element(“O”): 3}, you can obtain the amount of Fe
simply by comp[“Fe”] instead of the more verbose comp[Element(“Fe”)].

```python
>>> comp = Composition("LiFePO4")
>>> comp.get_atomic_fraction(Element("Li"))
0.14285714285714285
>>> comp.num_atoms
7.0
>>> comp.reduced_formula
'LiFePO4'
>>> comp.formula
'Li1 Fe1 P1 O4'
>>> comp.get_wt_fraction(Element("Li"))
0.04399794666951898
>>> comp.num_atoms
7.0
```

Very flexible Composition construction, similar to the built-in Python
dict(). Also extended to allow simple string init.

Takes any inputs supported by the Python built-in dict function.


1. A dict of either {Element/Species: amount},

> {string symbol:amount}, or {atomic number:amount} or any mixture
> of these. E.g., {Element(“Li”): 2, Element(“O”): 1},
> {“Li”:2, “O”:1}, {3: 2, 8: 1} all result in a Li2O composition.


2. Keyword arg initialization, similar to a dict, e.g.,

> Composition(Li = 2, O = 1)

In addition, the Composition constructor also allows a single
string as an input formula. E.g., Composition(“Li2O”).


* **Parameters**


    * **\*args** – Any number of 2-tuples as key-value pairs.


    * **strict** (*bool*) – Only allow valid Elements and Species in the Composition. Defaults to False.


    * **allow_negative** (*bool*) – Whether to allow negative compositions. Defaults to False.


    * **\*\*kwargs** – Additional kwargs supported by the dict() constructor.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _comps_from_fuzzy_formula(fuzzy_formula: str, m_dict: dict[str, float] | None = None, m_points: int = 0, factor: float = 1)
A recursive helper method for formula parsing that helps in
interpreting and ranking indeterminate formulas.
Author: Anubhav Jain.


* **Parameters**


    * **fuzzy_formula** (*str*) – A formula string, such as “co2o3” or “MN”,
    that may or may not have multiple interpretations.


    * **m_dict** (*dict*) – A symbol:amt dictionary from the previously parsed
    formula.


    * **m_points** – Number of points gained from the previously parsed
    formula.


    * **factor** – Coefficient for this parse, e.g. (PO4)2 will feed in PO4
    as the fuzzy_formula with a coefficient of 2.



* **Returns**

    A list of tuples, with the first element being a Composition

        and the second element being the number of points awarded that Composition interpretation.




* **Return type**

    list[tuple[Composition, int]]



#### _get_oxi_state_guesses(all_oxi_states, max_sites, oxi_states_override, target_charge)
Utility operation for guessing oxidation states.

See oxi_state_guesses for full details. This operation does the
calculation of the most likely oxidation states


* **Parameters**


    * **oxi_states_override** (*dict*) – dict of str->list to override an
    element’s common oxidation states, e.g. {“V”: [2,3,4,5]}


    * **target_charge** (*int*) – the desired total charge on the structure.
    Default is 0 signifying charge balance.


    * **all_oxi_states** (*bool*) – if True, an element defaults to
    all oxidation states in pymatgen Element.icsd_oxidation_states.
    Otherwise, default is Element.common_oxidation_states. Note
    that the full oxidation state list is *very* inclusive and
    can produce nonsensical results.


    * **max_sites** (*int*) – if possible, will reduce Compositions to at most
    this many sites to speed up oxidation state guesses. If the
    composition cannot be reduced to this many sites a ValueError
    will be raised. Set to -1 to just reduce fully. If set to a
    number less than -1, the formula will be fully reduced but a
    ValueError will be thrown if the number of atoms in the reduced
    formula is greater than abs(max_sites).



* **Returns**

    Each dict maps the element symbol to a list of

        oxidation states for each site of that element. For example, Fe3O4 could
        return a list of [2,2,2,3,3,3] for the oxidation states of the 6 Fe sites.
        If the composition is not charge balanced, an empty list is returned.




* **Return type**

    list[dict]



#### _parse_formula(formula: str)

* **Parameters**

    **formula** (*str*) – A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.



* **Returns**

    Composition with that formula.


### Notes

In the case of Metallofullerene formula (e.g. [Y3N@C80](mailto:Y3N@C80)),
the @ mark will be dropped and passed to parser.


#### add_charges_from_oxi_state_guesses(oxi_states_override: dict | None = None, target_charge: float = 0, all_oxi_states: bool = False, max_sites: int | None = None)
Assign oxidation states based on guessed oxidation states.

See oxi_state_guesses for an explanation of how oxidation states are
guessed. This operation uses the set of oxidation states for each site
that were determined to be most likely from the oxidation state guessing
routine.


* **Parameters**


    * **oxi_states_override** (*dict*) – dict of str->list to override an
    element’s common oxidation states, e.g. {“V”: [2,3,4,5]}


    * **target_charge** (*int*) – the desired total charge on the structure.
    Default is 0 signifying charge balance.


    * **all_oxi_states** (*bool*) – if True, an element defaults to
    all oxidation states in pymatgen Element.icsd_oxidation_states.
    Otherwise, default is Element.common_oxidation_states. Note
    that the full oxidation state list is *very* inclusive and
    can produce nonsensical results.


    * **max_sites** (*int*) – if possible, will reduce Compositions to at most
    this many sites to speed up oxidation state guesses. If the
    composition cannot be reduced to this many sites a ValueError
    will be raised. Set to -1 to just reduce fully. If set to a
    number less than -1, the formula will be fully reduced but a
    ValueError will be thrown if the number of atoms in the reduced
    formula is greater than abs(max_sites).



* **Returns**

    Composition, where the elements are assigned oxidation states based
    on the results form guessing oxidation states. If no oxidation state
    is possible, returns a Composition where all oxidation states are 0.



#### almost_equals(other: Composition, rtol: float = 0.1, atol: float = 1e-08)
Returns true if compositions are equal within a tolerance.


* **Parameters**


    * **other** (*Composition*) – Other composition to check


    * **rtol** (*float*) – Relative tolerance


    * **atol** (*float*) – Absolute tolerance



#### _property_ alphabetical_formula(_: st_ )
Returns a formula string, with elements sorted by alphabetically
e.g., Fe4 Li4 O16 P4.


#### amount_tolerance(_ = 1e-0_ )

#### _property_ anonymized_formula(_: st_ )
An anonymized formula. Unique species are arranged in ordering of
increasing amounts and assigned ascending alphabets. Useful for
prototyping formulas. For example, all stoichiometric perovskites have
anonymized_formula ABC3.


#### as_dict()
Subtly different from get_el_amt_dict in that they keys here are str(Element)
instead of Element.symbol.


* **Returns**

    element symbol and (unreduced) amount. E.g.

        {“Fe”: 4.0, “O”:6.0} or {“Fe3+”: 4.0, “O2-“:6.0}




* **Return type**

    dict[str, float]



#### _property_ average_electroneg(_: floa_ )
Average electronegativity of the composition.


#### _property_ chemical_system(_: st_ )
Get the chemical system of a Composition, for example “O-Si” for
SiO2. Chemical system is a string of a list of elements
sorted alphabetically and joined by dashes, by convention for use
in database keys.


#### contains_element_type(category: str)
Check if Composition contains any elements matching a given category.


* **Parameters**

    **category** (*str*) – one of “noble_gas”, “transition_metal”,
    “post_transition_metal”, “rare_earth_metal”, “metal”, “metalloid”,
    “alkali”, “alkaline”, “halogen”, “chalcogen”, “lanthanoid”,
    “actinoid”, “quadrupolar”, “s-block”, “p-block”, “d-block”, “f-block”



* **Returns**

    True if any elements in Composition match category, otherwise False



* **Return type**

    bool



#### copy()
A copy of the composition.


#### _property_ element_composition(_: Compositio_ )
Returns the composition replacing any species by the corresponding element.


#### _property_ elements(_: list[Element | Species | DummySpecies_ )
Returns list of elements in Composition.


#### _property_ formula(_: st_ )
Returns a formula string, with elements sorted by electronegativity,
e.g., Li4 Fe4 P4 O16.


#### _property_ fractional_composition(_: Compositio_ )
Returns the normalized composition in which the amounts of each species sum to
1.
E.g. “Fe2 O3”.fractional_composition = “Fe0.4 O0.6”.


#### _classmethod_ from_dict(d)
Creates a composition from a dict generated by as_dict(). Strictly not
necessary given that the standard constructor already takes in such an
input, but this method preserves the standard pymatgen API of having
from_dict methods to reconstitute objects generated by as_dict(). Allows
for easier introspection.


* **Parameters**

    **d** (*dict*) – {symbol: amount} dict.



#### _classmethod_ from_weight_dict(weight_dict)
Creates a Composition based on a dict of atomic fractions calculated
from a dict of weight fractions. Allows for quick creation of the class
from weight-based notations commonly used in the industry, such as
Ti6V4Al and Ni60Ti40.


* **Parameters**

    **weight_dict** (*dict*) – {symbol: weight_fraction} dict.



* **Returns**

    Composition



#### get_atomic_fraction(el: str | Element | Species | DummySpecies)
Calculate atomic fraction of an Element or Species.


* **Parameters**

    **el** (*Element/Species*) – Element or Species to get fraction for.



* **Returns**

    Atomic fraction for element el in Composition



#### get_el_amt_dict()

* **Returns**

    element symbol and (unreduced) amount. E.g.
    {“Fe”: 4.0, “O”:6.0} or {“Fe3+”: 4.0, “O2-“:6.0}.



* **Return type**

    dict[str, float]



#### get_integer_formula_and_factor(max_denominator: int = 10000, iupac_ordering: bool = False)
Calculates an integer formula and factor.


* **Parameters**


    * **max_denominator** (*int*) – all amounts in the el:amt dict are
    first converted to a Fraction with this maximum denominator


    * **iupac_ordering** (*bool**, **optional*) – Whether to order the
    formula by the iupac “electronegativity” series, defined in
    Table VI of “Nomenclature of Inorganic Chemistry (IUPAC
    Recommendations 2005)”. This ordering effectively follows
    the groups and rows of the periodic table, except the
    Lanthanides, Actinides and hydrogen. Note that polyanions
    will still be determined based on the true electronegativity of
    the elements.



* **Returns**

    A pretty normalized formula and a multiplicative factor, i.e.,
    Li0.5O0.25 returns (Li2O, 0.25). O0.25 returns (O2, 0.125)



#### get_reduced_composition_and_factor()
Calculates a reduced composition and factor.


* **Returns**

    A normalized composition and a multiplicative factor, i.e.,
    Li4Fe4P4O16 returns (Composition(“LiFePO4”), 4).



#### get_reduced_formula_and_factor(iupac_ordering: bool = False)
Calculates a reduced formula and factor.


* **Parameters**

    **iupac_ordering** (*bool**, **optional*) – Whether to order the
    formula by the iupac “electronegativity” series, defined in
    Table VI of “Nomenclature of Inorganic Chemistry (IUPAC
    Recommendations 2005)”. This ordering effectively follows
    the groups and rows of the periodic table, except the
    Lanthanides, Actinides and hydrogen. Note that polyanions
    will still be determined based on the true electronegativity of
    the elements.



* **Returns**

    A pretty normalized formula and a multiplicative factor, i.e.,
    Li4Fe4P4O16 returns (LiFePO4, 4).



#### get_wt_fraction(el: str | Element | Species | DummySpecies)
Calculate weight fraction of an Element or Species.


* **Parameters**

    **el** (*Element** | **Species*) – Element or Species to get fraction for.



* **Returns**

    Weight fraction for element el in Composition.



* **Return type**

    float



#### _property_ hill_formula(_: st_ )
The Hill system (or Hill notation) is a system of writing empirical chemical
formulas, molecular chemical formulas and components of a condensed formula such
that the number of carbon atoms in a molecule is indicated first, the number of
hydrogen atoms next, and then the number of all other chemical elements
subsequently, in alphabetical order of the chemical symbols. When the formula
contains no carbon, all the elements, including hydrogen, are listed
alphabetically.


#### _property_ is_element(_: boo_ )
True if composition is an element.


#### _property_ iupac_formula(_: st_ )
Returns a formula string, with elements sorted by the iupac
electronegativity ordering defined in Table VI of “Nomenclature of
Inorganic Chemistry (IUPAC Recommendations 2005)”. This ordering
effectively follows the groups and rows of the periodic table, except
the Lanthanides, Actinides and hydrogen. Polyanions are still determined
based on the true electronegativity of the elements.
e.g. CH2(SO4)2.


#### _property_ num_atoms(_: floa_ )
Total number of atoms in Composition. For negative amounts, sum
of absolute values.


#### oxi_prob(_ = Non_ )

#### oxi_state_guesses(oxi_states_override: dict | None = None, target_charge: float = 0, all_oxi_states: bool = False, max_sites: int | None = None)
Checks if the composition is charge-balanced and returns back all
charge-balanced oxidation state combinations. Composition must have
integer values. Note that more num_atoms in the composition gives
more degrees of freedom. e.g., if possible oxidation states of
element X are [2,4] and Y are [-3], then XY is not charge balanced
but X2Y2 is. Results are returned from most to least probable based
on ICSD statistics. Use max_sites to improve performance if needed.


* **Parameters**


    * **oxi_states_override** (*dict*) – dict of str->list to override an
    element’s common oxidation states, e.g. {“V”: [2,3,4,5]}


    * **target_charge** (*int*) – the desired total charge on the structure.
    Default is 0 signifying charge balance.


    * **all_oxi_states** (*bool*) – if True, an element defaults to
    all oxidation states in pymatgen Element.icsd_oxidation_states.
    Otherwise, default is Element.common_oxidation_states. Note
    that the full oxidation state list is *very* inclusive and
    can produce nonsensical results.


    * **max_sites** (*int*) – if possible, will reduce Compositions to at most
    this many sites to speed up oxidation state guesses. If the
    composition cannot be reduced to this many sites a ValueError
    will be raised. Set to -1 to just reduce fully. If set to a
    number less than -1, the formula will be fully reduced but a
    ValueError will be thrown if the number of atoms in the reduced
    formula is greater than abs(max_sites).



* **Returns**

    each dict reports an element symbol and average

        oxidation state across all sites in that composition. If the
        composition is not charge balanced, an empty list is returned.




* **Return type**

    list[dict]



#### _static_ ranked_compositions_from_indeterminate_formula(fuzzy_formula: str, lock_if_strict: bool = True)
Takes in a formula where capitalization might not be correctly entered,
and suggests a ranked list of potential Composition matches.
Author: Anubhav Jain.


* **Parameters**


    * **fuzzy_formula** (*str*) – A formula string, such as “co2o3” or “MN”,
    that may or may not have multiple interpretations


    * **lock_if_strict** (*bool*) – If true, a properly entered formula will
    only return the one correct interpretation. For example,
    “Co1” will only return “Co1” if true, but will return both
    “Co1” and “C1 O1” if false.



* **Returns**

    A ranked list of potential Composition matches



#### _property_ reduced_composition(_: Compositio_ )
Returns the reduced composition, i.e. amounts normalized by greatest common denominator.
E.g. “Fe4 P4 O16”.reduced_composition = “Fe P O4”.


#### _property_ reduced_formula(_: st_ )
Returns a pretty normalized formula, i.e., LiFePO4 instead of
Li4Fe4P4O16.


#### remove_charges()
Returns a new Composition with charges from each Species removed.


* **Returns**

    Composition object without charge decoration, for example
    {“Fe3+”: 2.0, “O2-“:3.0} becomes {“Fe”: 2.0, “O”:3.0}



#### replace(elem_map: dict[str, str | dict[str, float]])
Replace elements in a composition. Returns a new Composition, leaving the old one unchanged.


* **Parameters**

    **elem_map** (*dict**[**str**, **str** | **dict**[**str**, **float**]**]*) – dict of elements or species to swap. E.g.
    {“Li”: “Na”} performs a Li for Na substitution. The target can be a {species: factor} dict. For
    example, in Fe2O3 you could map {“Fe”: {“Mg”: 0.5, “Cu”:0.5}} to obtain MgCuO3.



* **Returns**

    New object with elements remapped according to elem_map.



* **Return type**

    Composition



#### special_formulas(_ = {'Cl': 'Cl2', 'CsO': 'Cs2O2', 'F': 'F2', 'H': 'H2', 'HO': 'H2O2', 'KO': 'K2O2', 'LiO': 'Li2O2', 'N': 'N2', 'NaO': 'Na2O2', 'O': 'O2', 'RbO': 'Rb2O2'_ )

#### _property_ to_data_dict(_: dic_ )
Returns:
A dict with many keys and values relating to Composition/Formula,
including reduced_cell_composition, unit_cell_composition,
reduced_cell_formula, elements and nelements.


#### to_pretty_string()

* **Returns**

    Same output as __str__() but without spaces.



* **Return type**

    str



#### _property_ to_reduced_dict(_: dict[str, float_ )
Returns:
dict[str, float]: element symbols mapped to reduced amount e.g. {“Fe”: 2.0, “O”:3.0}.


#### _property_ to_weight_dict(_: dict[str, float_ )
Returns:
dict[str, float] with weight fraction of each component {“Ti”: 0.90, “V”: 0.06, “Al”: 0.04}.


#### _property_ total_electrons(_: floa_ )
Total number of electrons in composition.


#### _property_ valid(_: boo_ )
Returns True if Composition contains valid elements or species and
False if the Composition contains any dummy species.


#### _property_ weight(_: floa_ )
Total molecular weight of Composition.


### _exception_ CompositionError()
Bases: `Exception`

Exception class for composition errors.


### reduce_formula(sym_amt, iupac_ordering: bool = False)
Helper method to reduce a sym_amt dict to a reduced formula and factor.


* **Parameters**


    * **sym_amt** (*dict*) – {symbol: amount}.


    * **iupac_ordering** (*bool**, **optional*) – Whether to order the
    formula by the iupac “electronegativity” series, defined in
    Table VI of “Nomenclature of Inorganic Chemistry (IUPAC
    Recommendations 2005)”. This ordering effectively follows
    the groups and rows of the periodic table, except the
    Lanthanides, Actinides and hydrogen. Note that polyanions
    will still be determined based on the true electronegativity of
    the elements.



* **Returns**

    (reduced_formula, factor).


## pymatgen.core.interface module

This module provides classes to store, generate, and manipulate material interfaces.


### _class_ Interface(lattice, species, coords, site_properties, validate_proximity=False, to_unit_cell=False, coords_are_cartesian=False, in_plane_offset: tuple[float, float] = (0, 0), gap: float = 0, vacuum_over_film: float = 0, interface_properties: dict | None = None)
Bases: `Structure`

This class stores data for defining an interface between two structures.
It is a subclass of pymatgen.core.structure.Structure.

Makes an interface structure, a structure object with additional information
and methods pertaining to interfaces.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** (*[**Species**]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to translate sites into the unit cell. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **in_plane_offset** – fractional shift in plane for the film with respect
    to the substrate


    * **gap** – gap between substrate and film in Angstroms; zero corresponds to
    the original distance between substrate and film sites


    * **vacuum_over_film** – vacuum space above the film in Angstroms. Defaults to 0.


    * **interface_properties** – properties associated with the interface. Defaults to None.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### _update_c(new_c: float)
Modifies the c-direction of the lattice without changing the site Cartesian coordinates
Be careful you can mess up the interface by setting a c-length that can’t accommodate all the sites.


#### as_dict()
MSONable dict.


#### copy()

* **Returns**

    A copy of the Interface.



* **Return type**

    Interface



#### _property_ film(_: Structur_ )
A pymatgen Structure for just the film.


#### _property_ film_indices(_: list[int_ )
Site indices of the film sites.


#### _property_ film_layers(_: in_ )
Number of layers of the minimum element in the film composition.


#### _property_ film_sites(_: list[Site_ )
Return the film sites of the interface.


#### _property_ film_termination(_: st_ )
Label for the film termination chemistry.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    Creates slab from dict.



#### _classmethod_ from_slabs(substrate_slab: Slab, film_slab: Slab, in_plane_offset: tuple[float, float] = (0, 0), gap: float = 1.6, vacuum_over_film: float = 0, interface_properties: dict | None = None, center_slab: bool = True)
Makes an interface structure by merging a substrate and film slabs
The film a- and b-vectors will be forced to be the substrate slab’s
a- and b-vectors.

For now, it’s suggested to use a factory method that will ensure the
appropriate interface structure is already met.


* **Parameters**


    * **substrate_slab** (*Slab*) – slab for the substrate


    * **film_slab** (*Slab*) – slab for the film


    * **in_plane_offset** (*tuple*) – fractional shift in plane for the film with respect to the substrate.
    For example, (0.5, 0.5) will shift the film by half the substrate’s a- and b-vectors.
    Defaults to (0, 0).


    * **gap** (*float*) – gap between substrate and film in Angstroms. Defaults to 1.6.


    * **vacuum_over_film** (*float*) – vacuum space above the film in Angstroms. Defaults to 0.


    * **interface_properties** (*dict*) – misc properties to assign to the interface. Defaults to None.


    * **center_slab** (*bool*) – center the slab. Defaults to True.



#### _property_ gap(_: floa_ )
The gap in Cartesian units between the film and the substrate.


#### get_shifts_based_on_adsorbate_sites(tolerance: float = 0.1)
Computes possible in-plane shifts based on an adsorbate site  algorithm.


* **Parameters**

    **tolerance** – tolerance for “uniqueness” for shifts in Cartesian unit
    This is usually Angstroms.



#### get_sorted_structure(key=None, reverse=False)
Get a sorted structure for the interface. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### _property_ in_plane_offset(_: ndarra_ )
The shift between the film and substrate in fractional
coordinates.


#### _property_ substrate(_: Structur_ )
A pymatgen Structure for just the substrate.


#### _property_ substrate_indices(_: list[int_ )
Site indices for the substrate atoms.


#### _property_ substrate_layers(_: in_ )
Number of layers of the minimum element in the substrate composition.


#### _property_ substrate_sites(_: list[Site_ )
The site objects in the substrate.


#### _property_ substrate_termination(_: st_ )
Label for the substrate termination chemistry.


#### _property_ vacuum_over_film(_: floa_ )
The vacuum space over the film in Cartesian units.


### count_layers(struct: Structure, el=None)
Counts the number of ‘layers’ along the c-axis.


### label_termination(slab: Structure)
Labels the slab surface termination.

## pymatgen.core.ion module

Module containing class to create an ion.


### _class_ Ion(composition, charge=0.0, _properties=None)
Bases: `Composition`, `MSONable`, [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)

Ion object. Just a Composition object with an additional variable to store
charge.

The net charge can either be represented as Mn++, Mn+2, Mn[2+], Mn[++], or
Mn[+2]. Note the order of the sign and magnitude in each representation.

Flexible Ion construction, similar to Composition.
For more information, please see pymatgen.core.Composition.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ alphabetical_formula(_: st_ )
Returns a formula string, with elements sorted by alphabetically and
appended charge.


#### _property_ anonymized_formula(_: st_ )
An anonymized formula. Appends charge to the end
of anonymized composition.


#### as_dict()

* **Returns**

    dict with composition, as well as charge.



#### _property_ charge(_: floa_ )
Charge of the ion.


#### _property_ composition(_: Compositio_ )
Composition of ion.


#### _property_ formula(_: st_ )
Returns a formula string with appended charge. The
charge is written with the sign preceding the magnitude, e.g.,
‘Ca1 +2’. Uncharged species have “(aq)” appended, e.g. “O2 (aq)”.


#### _classmethod_ from_dict(dct)
Generates an ion object from a dict created by as_dict().


* **Parameters**

    **dct** – {symbol: amount} dict.



#### _classmethod_ from_formula(formula: str)
Creates Ion from formula. The net charge can either be represented as
Mn++, Mn+2, Mn[2+], Mn[++], or Mn[+2]. Note the order of the sign and
magnitude in each representation.

Also note that (aq) can be included in the formula, e.g. “NaOH (aq)”.


* **Parameters**

    **formula** –



* **Returns**

    Ion



#### get_reduced_formula_and_factor(iupac_ordering: bool = False, hydrates: bool = True)
Calculates a reduced formula and factor.

Similar to Composition.get_reduced_formula_and_factor except that O-H formulas
receive special handling to differentiate between hydrogen peroxide and OH-.
Formulas containing HO are written with oxygen first (e.g. ‘Fe(OH)2’ rather than
‘Fe(HO)2’), and special formulas that apply to solids (e.g. Li2O2 instead of LiO)
are not used.

Note that the formula returned by this method does not contain a charge.
To include charge, use formula or reduced_formula instead.


* **Parameters**


    * **iupac_ordering** (*bool**, **optional*) – Whether to order the
    formula by the iupac “electronegativity” series, defined in
    Table VI of “Nomenclature of Inorganic Chemistry (IUPAC
    Recommendations 2005)”. This ordering effectively follows
    the groups and rows of the periodic table, except the
    Lanthanides, Actinides and hydrogen. Note that polyanions
    will still be determined based on the true electronegativity of
    the elements.


    * **hydrates** – If True (default), attempt to recognize hydrated metal
    complexes and separate out the H2O in the reduced formula.
    For example, Zr(OH)4 becomes ZrO2.2H2O. Applies only to
    Ions containing metals.



* **Returns**

    A pretty normalized formula and a multiplicative factor, i.e.,
    H4O4 returns (‘H2O2’, 2.0).



#### oxi_state_guesses(oxi_states_override: dict | None = None, all_oxi_states: bool = False, max_sites: int | None = None)
Checks if the composition is charge-balanced and returns back all
charge-balanced oxidation state combinations. Composition must have
integer values. Note that more num_atoms in the composition gives
more degrees of freedom. e.g., if possible oxidation states of
element X are [2,4] and Y are [-3], then XY is not charge balanced
but X2Y2 is. Results are returned from most to least probable based
on ICSD statistics. Use max_sites to improve performance if needed.


* **Parameters**


    * **oxi_states_override** (*dict*) – dict of str->list to override an
    element’s common oxidation states, e.g. {“V”: [2,3,4,5]}


    * **all_oxi_states** (*bool*) – if True, an element defaults to
    all oxidation states in pymatgen Element.icsd_oxidation_states.
    Otherwise, default is Element.common_oxidation_states. Note
    that the full oxidation state list is *very* inclusive and
    can produce nonsensical results.


    * **max_sites** (*int*) – if possible, will reduce Compositions to at most
    this many sites to speed up oxidation state guesses. If the
    composition cannot be reduced to this many sites a ValueError
    will be raised. Set to -1 to just reduce fully. If set to a
    number less than -1, the formula will be fully reduced but a
    ValueError will be thrown if the number of atoms in the reduced
    formula is greater than abs(max_sites).



* **Returns**

    A list of dicts - each dict reports an element symbol and average

        oxidation state across all sites in that composition. If the
        composition is not charge balanced, an empty list is returned.




#### _property_ reduced_formula(_: st_ )
Returns a reduced formula string with appended charge. The
charge is placed in brackets with the sign preceding the magnitude, e.g.,
‘Ca[+2]’. Uncharged species have “(aq)” appended, e.g. “O2(aq)”.


#### to_pretty_string()
Pretty string with proper superscripts.


#### _property_ to_reduced_dict(_: dic_ )

* **Returns**

    dict with element symbol and reduced amount e.g.,


{“Fe”: 2.0, “O”:3.0}.

## pymatgen.core.lattice module

Defines the classes relating to 3D lattices.


### _class_ Lattice(matrix: ArrayLike, pbc: tuple[bool, bool, bool] = (True, True, True))
Bases: `MSONable`

A lattice object. Essentially a matrix with conversion matrices. In
general, it is assumed that length units are in Angstroms and angles are in
degrees unless otherwise stated.

Create a lattice from any sequence of 9 numbers. Note that the sequence
is assumed to be read one row at a time. Each row represents one
lattice vector.


* **Parameters**


    * **matrix** – Sequence of numbers in any form. Examples of acceptable
    input.
    i) An actual numpy array.
    ii) [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
    iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)
    Each row should correspond to a lattice vector.
    E.g., [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
    with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].


    * **pbc** – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



#### _calculate_lll(delta: float = 0.75)
Performs a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
c-reduced basis. This method returns a basis which is as “good” as
possible, with “good” defined by orthongonality of the lattice vectors.

This basis is used for all the periodic boundary condition calculations.


* **Parameters**

    **delta** (*float*) – Reduction parameter. Default of 0.75 is usually
    fine.



* **Returns**

    Reduced lattice matrix, mapping to get to that lattice.



#### _property_ a(_: floa_ )
*a* lattice parameter.


#### _property_ abc(_: Vector3_ )
Lengths of the lattice vectors, i.e. (a, b, c).


#### _property_ alpha(_: floa_ )
Angle alpha of lattice in degrees.


#### _property_ angles(_: Vector3_ )
Lattice angles.


* **Returns**

    The angles (alpha, beta, gamma) of the lattice.



#### as_dict(verbosity: int = 0)
MSONable dict representation of the Lattice.


* **Parameters**

    **verbosity** (*int*) – Default of 0 only includes the matrix representation.
    Set to 1 to include the lattice parameters.



#### _property_ b(_: floa_ )
*b* lattice parameter.


#### _property_ beta(_: floa_ )
Angle beta of lattice in degrees.


#### _property_ c(_: floa_ )
*c* lattice parameter.


#### copy()
Deep copy of self.


#### _static_ cubic(a: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a cubic lattice.


* **Parameters**


    * **a** (*float*) – The *a* lattice parameter of the cubic cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Cubic lattice of dimensions a x a x a.



#### d_hkl(miller_index: ArrayLike)
Returns the distance between the hkl plane and the origin.


* **Parameters**

    **miller_index** (*[**h**,**k**,**l**]*) – Miller index of plane



* **Returns**

    d_hkl (float)



#### dot(coords_a: ArrayLike, coords_b: ArrayLike, frac_coords: bool = False)
Compute the scalar product of vector(s).


* **Parameters**


    * **coords_a** – Array-like coordinates.


    * **coords_b** – Array-like coordinates.


    * **frac_coords** (*bool*) – Boolean stating whether the vector
    corresponds to fractional or Cartesian coordinates.



* **Returns**

    one-dimensional numpy array.



#### find_all_mappings(other_lattice: Lattice, ltol: float = 1e-05, atol: float = 1, skip_rotation_matrix: bool = False)
Finds all mappings between current lattice and another lattice.


* **Parameters**


    * **other_lattice** (*Lattice*) – Another lattice that is equivalent to this one.


    * **ltol** (*float*) – Tolerance for matching lengths. Defaults to 1e-5.


    * **atol** (*float*) – Tolerance for matching angles. Defaults to 1.


    * **skip_rotation_matrix** (*bool*) – Whether to skip calculation of the
    rotation matrix



* **Yields**

    (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
    found. aligned_lattice is a rotated version of other_lattice that
    has the same lattice parameters, but which is aligned in the
    coordinate system of this lattice so that translational points
    match up in 3D. rotation_matrix is the rotation that has to be
    applied to other_lattice to obtain aligned_lattice, i.e.,
    aligned_matrix = np.inner(other_lattice, rotation_matrix) and
    op = SymmOp.from_rotation_and_translation(rotation_matrix)
    aligned_matrix = op.operate_multi(latt.matrix)
    Finally, scale_matrix is the integer matrix that expresses
    aligned_matrix as a linear combination of this
    lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)

    None is returned if no matches are found.



#### find_mapping(other_lattice: Lattice, ltol: float = 1e-05, atol: float = 1, skip_rotation_matrix: bool = False)
Finds a mapping between current lattice and another lattice. There
are an infinite number of choices of basis vectors for two entirely
equivalent lattices. This method returns a mapping that maps
other_lattice to this lattice.


* **Parameters**


    * **other_lattice** (*Lattice*) – Another lattice that is equivalent to
    this one.


    * **ltol** (*float*) – Tolerance for matching lengths. Defaults to 1e-5.


    * **atol** (*float*) – Tolerance for matching angles. Defaults to 1.


    * **skip_rotation_matrix** (*bool*) – Whether to skip calculation of the rotation matrix.
    Defaults to False.



* **Returns**

    (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
    found. aligned_lattice is a rotated version of other_lattice that
    has the same lattice parameters, but which is aligned in the
    coordinate system of this lattice so that translational points
    match up in 3D. rotation_matrix is the rotation that has to be
    applied to other_lattice to obtain aligned_lattice, i.e.,
    aligned_matrix = np.inner(other_lattice, rotation_matrix) and
    op = SymmOp.from_rotation_and_translation(rotation_matrix)
    aligned_matrix = op.operate_multi(latt.matrix)
    Finally, scale_matrix is the integer matrix that expresses
    aligned_matrix as a linear combination of this
    lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)

    None is returned if no matches are found.




#### _classmethod_ from_dict(d: dict, fmt: str | None = None, \*\*kwargs)
Create a Lattice from a dictionary containing the a, b, c, alpha, beta,
and gamma parameters if fmt is None.

If fmt == “abivars”, the function build a Lattice object from a
dictionary with the Abinit variables acell and rprim in Bohr.
If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr

### Example

Lattice.from_dict(fmt=”abivars”, acell=3\*[10], rprim=np.eye(3))


#### _classmethod_ from_parameters(a: float, b: float, c: float, alpha: float, beta: float, gamma: float, vesta: bool = False, pbc: tuple[bool, bool, bool] = (True, True, True))
Create a Lattice using unit cell lengths (in Angstrom) and angles (in degrees).


* **Parameters**


    * **a** (*float*) – *a* lattice parameter.


    * **b** (*float*) – *b* lattice parameter.


    * **c** (*float*) – *c* lattice parameter.


    * **alpha** (*float*) – *alpha* angle in degrees.


    * **beta** (*float*) – *beta* angle in degrees.


    * **gamma** (*float*) – *gamma* angle in degrees.


    * **vesta** – True if you import Cartesian coordinates from VESTA.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Lattice with the specified lattice parameters.



#### _property_ gamma(_: floa_ )
Angle gamma of lattice in degrees.


#### get_all_distances(fcoords1: ArrayLike, fcoords2: ArrayLike)
Returns the distances between two lists of coordinates taking into
account periodic boundary conditions and the lattice. Note that this
computes an MxN array of distances (i.e. the distance between each
point in fcoords1 and every coordinate in fcoords2). This is
different functionality from pbc_diff.


* **Parameters**


    * **fcoords1** – First set of fractional coordinates. e.g., [0.5, 0.6,
    0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
    coord or any array of coords.


    * **fcoords2** – Second set of fractional coordinates.



* **Returns**

    2d array of Cartesian distances. E.g the distance between
    fcoords1[i] and fcoords2[j] is distances[i,j]



#### get_brillouin_zone()
Returns the Wigner-Seitz cell for the reciprocal lattice, aka the
Brillouin Zone.


* **Returns**

    A list of list of coordinates.
    Each element in the list is a “facet” of the boundary of the
    Brillouin Zone. For instance, a list of four coordinates will
    represent a square facet.



#### get_cartesian_coords(fractional_coords: ArrayLike)
Returns the Cartesian coordinates given fractional coordinates.


* **Parameters**

    **fractional_coords** (*3x1 array*) – Fractional coords.



* **Returns**

    Cartesian coordinates



#### get_distance_and_image(frac_coords1: ArrayLike, frac_coords2: ArrayLike, jimage: ArrayLike | None = None)
Gets distance between two frac_coords assuming periodic boundary
conditions. If the index jimage is not specified it selects the j
image nearest to the i atom and returns the distance and jimage
indices in terms of lattice vector translations. If the index jimage
is specified it returns the distance between the frac_coords1 and
the specified jimage of frac_coords2, and the given jimage is also
returned.


* **Parameters**


    * **frac_coords1** (*3x1 array*) – Reference fcoords to get distance from.


    * **frac_coords2** (*3x1 array*) – fcoords to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of
    lattice translations, e.g., [1,0,0] implies to take periodic
    image that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    distance and periodic lattice translations
    of the other site for which the distance applies. This means that
    the distance between frac_coords1 and (jimage + frac_coords2) is
    equal to distance.



* **Return type**

    (distance, jimage)



#### get_frac_coords_from_lll(lll_frac_coords: ArrayLike)
Given fractional coordinates in the lll basis, returns corresponding
fractional coordinates in the lattice basis.


#### get_fractional_coords(cart_coords: ArrayLike)
Returns the fractional coordinates given Cartesian coordinates.


* **Parameters**

    **cart_coords** (*3x1 array*) – Cartesian coords.



* **Returns**

    Fractional coordinates.



#### get_lll_frac_coords(frac_coords: ArrayLike)
Given fractional coordinates in the lattice basis, returns corresponding
fractional coordinates in the lll basis.


#### get_lll_reduced_lattice(delta: float = 0.75)

* **Parameters**

    **delta** – Delta parameter.



* **Returns**

    LLL reduced Lattice.



#### get_miller_index_from_coords(coords: ArrayLike, coords_are_cartesian: bool = True, round_dp: int = 4, verbose: bool = True)
Get the Miller index of a plane from a list of site coordinates.

A minimum of 3 sets of coordinates are required. If more than 3 sets of
coordinates are given, the best plane that minimises the distance to all
points will be calculated.


* **Parameters**


    * **coords** (*iterable*) – A list or numpy array of coordinates. Can be
    Cartesian or fractional coordinates. If more than three sets of
    coordinates are provided, the best plane that minimises the
    distance to all sites will be calculated.


    * **coords_are_cartesian** (*bool**, **optional*) – Whether the coordinates are
    in Cartesian space. If using fractional coordinates set to
    False.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)



#### get_niggli_reduced_lattice(tol: float = 1e-05)
Get the Niggli reduced lattice using the numerically stable algo
proposed by R. W. Grosse-Kunstleve, N. K. Sauter, & P. D. Adams,
Acta Crystallographica Section A Foundations of Crystallography, 2003,
60(1), 1-6. doi:10.1107/S010876730302186X.


* **Parameters**

    **tol** (*float*) – The numerical tolerance. The default of 1e-5 should
    result in stable behavior for most cases.



* **Returns**

    Niggli-reduced lattice.



#### get_points_in_sphere(frac_points: ArrayLike, center: ArrayLike, r: float, zip_results=True)
Find all points within a sphere from the point taking into account
periodic boundary conditions. This includes sites in other periodic
images.

Algorithm:


1. place sphere of radius r in crystal and determine minimum supercell
(parallelepiped) which would contain a sphere of radius r. for this
we need the projection of a_1 on a unit vector perpendicular
to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
determine how many a_1”s it will take to contain the sphere.

Nxmax = r \* length_of_b_1 / (2 Pi)


2. keep points falling within r.


* **Parameters**


    * **frac_points** – All points in the lattice in fractional coordinates.


    * **center** – Cartesian coordinates of center of sphere.


    * **r** – radius of sphere.


    * **zip_results** (*bool*) – Whether to zip the results together to group by
    point, or return the raw fcoord, dist, index arrays



* **Returns**

    [(fcoord, dist, index, supercell_image) …] since most of the time, subsequent

        processing requires the distance, index number of the atom, or index of the image

    else:

        fcoords, dists, inds, image




* **Return type**

    if zip_results



#### get_points_in_sphere_old(\*\*kwargs)

#### get_points_in_sphere_py(frac_points: ArrayLike, center: ArrayLike, r: float, zip_results=True)
Find all points within a sphere from the point taking into account
periodic boundary conditions. This includes sites in other periodic
images.

Algorithm:


1. place sphere of radius r in crystal and determine minimum supercell
(parallelpiped) which would contain a sphere of radius r. for this
we need the projection of a_1 on a unit vector perpendicular
to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
determine how many a_1”s it will take to contain the sphere.

Nxmax = r \* length_of_b_1 / (2 Pi)


2. keep points falling within r.


* **Parameters**


    * **frac_points** – All points in the lattice in fractional coordinates.


    * **center** – Cartesian coordinates of center of sphere.


    * **r** – radius of sphere.


    * **zip_results** (*bool*) – Whether to zip the results together to group by
    point, or return the raw fcoord, dist, index arrays



* **Returns**

    [(fcoord, dist, index, supercell_image) …] since most of the time, subsequent

        processing requires the distance, index number of the atom, or index of the image

    else:

        fcoords, dists, inds, image




* **Return type**

    if zip_results



#### get_recp_symmetry_operation(symprec: float = 0.01)
Find the symmetric operations of the reciprocal lattice,
to be used for hkl transformations.


* **Parameters**

    **symprec** – default is 0.001.



#### get_vector_along_lattice_directions(cart_coords: ArrayLike)
Returns the coordinates along lattice directions given Cartesian coordinates.

Note, this is different than a projection of the Cartesian vector along the
lattice parameters. It is simply the fractional coordinates multiplied by the
lattice vector magnitudes.

For example, this method is helpful when analyzing the dipole moment (in
units of electron Angstroms) of a ferroelectric crystal. See the Polarization
class in pymatgen.analysis.ferroelectricity.polarization.


* **Parameters**

    **cart_coords** (*3x1 array*) – Cartesian coords.



* **Returns**

    Lattice coordinates.



#### get_wigner_seitz_cell()
Returns the Wigner-Seitz cell for the given lattice.


* **Returns**

    A list of list of coordinates.
    Each element in the list is a “facet” of the boundary of the
    Wigner Seitz cell. For instance, a list of four coordinates will
    represent a square facet.



#### _static_ hexagonal(a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a hexagonal lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the hexagonal cell.


    * **c** (*float*) – *c* lattice parameter of the hexagonal cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Hexagonal lattice of dimensions a x a x c.



#### _property_ inv_matrix(_: ndarra_ )
Inverse of lattice matrix.


#### _property_ is_3d_periodic(_: boo_ )
True if the Lattice is periodic in all directions.


#### is_hexagonal(hex_angle_tol: float = 5, hex_length_tol: float = 0.01)

* **Parameters**


    * **hex_angle_tol** – Angle tolerance


    * **hex_length_tol** – Length tolerance



* **Returns**

    Whether lattice corresponds to hexagonal lattice.



#### _property_ is_orthogonal(_: boo_ )
Whether all angles are 90 degrees.


#### _property_ lengths(_: Vector3_ )
Lattice lengths.


* **Returns**

    The lengths (a, b, c) of the lattice.



#### _property_ lll_inverse(_: ndarra_ )
Inverse of self.lll_mapping.


#### _property_ lll_mapping(_: ndarra_ )
The mapping between the LLL reduced lattice and the original lattice.


#### _property_ lll_matrix(_: ndarra_ )
The matrix for LLL reduction.


#### _property_ matrix(_: ndarra_ )
Copy of matrix representing the Lattice.


#### _property_ metric_tensor(_: ndarra_ )
The metric tensor of the lattice.


#### _static_ monoclinic(a: float, b: float, c: float, beta: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a monoclinic lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the monoclinc cell.


    * **b** (*float*) – *b* lattice parameter of the monoclinc cell.


    * **c** (*float*) – *c* lattice parameter of the monoclinc cell.


    * **beta** (*float*) – *beta* angle between lattice vectors b and c in
    degrees.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Monoclinic lattice of dimensions a x b x c with non right-angle
    beta between lattice vectors a and c.



#### norm(coords: ArrayLike, frac_coords: bool = True)
Compute the norm of vector(s).


* **Parameters**


    * **coords** – Array-like object with the coordinates.


    * **frac_coords** – Boolean stating whether the vector corresponds to fractional or
    Cartesian coordinates.



* **Returns**

    one-dimensional numpy array.



#### _static_ orthorhombic(a: float, b: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for an orthorhombic lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the orthorhombic cell.


    * **b** (*float*) – *b* lattice parameter of the orthorhombic cell.


    * **c** (*float*) – *c* lattice parameter of the orthorhombic cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Orthorhombic lattice of dimensions a x b x c.



#### _property_ parameters(_: tuple[float, float, float, float, float, float_ )
Returns (a, b, c, alpha, beta, gamma).


#### _property_ params_dict(_: dict[str, float_ )
Dictionary of lattice parameters.


#### _property_ pbc(_: tuple[bool, bool, bool_ )
Tuple defining the periodicity of the Lattice.


#### _property_ reciprocal_lattice(_: Lattic_ )
Return the reciprocal lattice. Note that this is the standard
reciprocal lattice used for solid state physics with a factor of 2 \*
pi. If you are looking for the crystallographic reciprocal lattice,
use the reciprocal_lattice_crystallographic property.
The property is lazily generated for efficiency.


#### _property_ reciprocal_lattice_crystallographic(_: Lattic_ )
Returns the *crystallographic* reciprocal lattice, i.e. no factor of 2 \* pi.


#### _static_ rhombohedral(a: float, alpha: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a rhombohedral lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the rhombohedral cell.


    * **alpha** (*float*) – Angle for the rhombohedral lattice in degrees.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Rhombohedral lattice of dimensions a x a x a.



#### scale(new_volume: float)
Return a new Lattice with volume new_volume by performing a
scaling of the lattice vectors so that length proportions and angles
are preserved.


* **Parameters**

    **new_volume** – New volume to scale to.



* **Returns**

    New lattice with desired volume.



#### selling_dist(other)
Returns the minimum Selling distance between two lattices.


#### _property_ selling_vector(_: ndarra_ )
Returns the (1,6) array of Selling Scalars.


#### _static_ tetragonal(a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True))
Convenience constructor for a tetragonal lattice.


* **Parameters**


    * **a** (*float*) – *a* lattice parameter of the tetragonal cell.


    * **c** (*float*) – *c* lattice parameter of the tetragonal cell.


    * **pbc** (*tuple*) – a tuple defining the periodic boundary conditions along the three
    axis of the lattice. If None periodic in all directions.



* **Returns**

    Tetragonal lattice of dimensions a x a x c.



#### _property_ volume(_: floa_ )
Volume of the unit cell in Angstrom^3.


### _compute_cube_index(coords: ndarray, global_min: float, radius: float)
Compute the cube index from coordinates
:param coords: (nx3 array) atom coordinates
:param global_min: (float) lower boundary of coordinates
:param radius: (float) cutoff radius.


* **Returns**

    nx3 array int indices



* **Return type**

    np.ndarray



### _one_to_three(label1d: ndarray, ny: int, nz: int)
Convert a 1D index array to 3D index array.


* **Parameters**


    * **label1d** – (array) 1D index array


    * **ny** – (int) number of cells in y direction


    * **nz** – (int) number of cells in z direction



* **Returns**

    nx3 array int indices



* **Return type**

    np.ndarray



### _three_to_one(label3d: ndarray, ny: int, nz: int)
The reverse of _one_to_three.


### find_neighbors(label: ndarray, nx: int, ny: int, nz: int)
Given a cube index, find the neighbor cube indices.


* **Parameters**


    * **label** – (array) (n,) or (n x 3) indice array


    * **nx** – (int) number of cells in y direction


    * **ny** – (int) number of cells in y direction


    * **nz** – (int) number of cells in z direction



* **Returns**

    Neighbor cell indices.



### get_integer_index(miller_index: Sequence[float], round_dp: int = 4, verbose: bool = True)
Attempt to convert a vector of floats to whole numbers.


* **Parameters**


    * **miller_index** (*list** of **float*) – A list miller indexes.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)



### get_points_in_spheres(all_coords: np.ndarray, center_coords: np.ndarray, r: float, pbc: bool | list[bool] | tuple[bool, bool, bool] = True, numerical_tol: float = 1e-08, lattice: Lattice | None = None, return_fcoords: bool = False)
For each point in center_coords, get all the neighboring points in all_coords that are within the
cutoff radius r.


* **Parameters**


    * **all_coords** – (list of Cartesian coordinates) all available points


    * **center_coords** – (list of Cartesian coordinates) all centering points


    * **r** – (float) cutoff radius


    * **pbc** – (bool or a list of bool) whether to set periodic boundaries


    * **numerical_tol** – (float) numerical tolerance


    * **lattice** – (Lattice) lattice to consider when PBC is enabled


    * **return_fcoords** – (bool) whether to return fractional coords when pbc is set.



* **Returns**

    List[List[Tuple[coords, distance, index, image]]]


## pymatgen.core.libxcfunc module

Enumerator with the libxc identifiers.
This is a low level object, client code should not interact with LibxcFunc directly
but use the API provided by the Xcfunc object defined in core.xcfunc.py.
Part of this module is automatically generated so be careful when refactoring stuff.
Use the script ~pymatgen/dev_scripts/regen_libxcfunc.py to regenerate the enum values.


### _class_ LibxcFunc(value)
Bases: `Enum`

Enumerator with the identifiers. This object is used by Xcfunc
declared in xcfunc.py to create an internal representation of the XC functional.
This is a low level object, client code should not interact with LibxcFunc directly
but use the API provided by Xcfunc.


* **Parameters**

    **num** – Number for the xc.



#### GGA_C_AM05(_ = 13_ )

#### GGA_C_APBE(_ = 18_ )

#### GGA_C_BGCP(_ = 3_ )

#### GGA_C_FT97(_ = 8_ )

#### GGA_C_GAM(_ = 3_ )

#### GGA_C_HCTH_A(_ = 9_ )

#### GGA_C_LM(_ = 13_ )

#### GGA_C_LYP(_ = 13_ )

#### GGA_C_N12(_ = 8_ )

#### GGA_C_N12_SX(_ = 7_ )

#### GGA_C_OPTC(_ = 20_ )

#### GGA_C_OP_B88(_ = 8_ )

#### GGA_C_OP_G96(_ = 8_ )

#### GGA_C_OP_PBE(_ = 8_ )

#### GGA_C_OP_PW91(_ = 26_ )

#### GGA_C_OP_XALPHA(_ = 8_ )

#### GGA_C_P86(_ = 13_ )

#### GGA_C_PBE(_ = 13_ )

#### GGA_C_PBEFE(_ = 25_ )

#### GGA_C_PBEINT(_ = 6_ )

#### GGA_C_PBELOC(_ = 24_ )

#### GGA_C_PBE_JRGX(_ = 13_ )

#### GGA_C_PBE_SOL(_ = 13_ )

#### GGA_C_PW91(_ = 13_ )

#### GGA_C_Q2D(_ = 4_ )

#### GGA_C_REGTPSS(_ = 8_ )

#### GGA_C_REVTCA(_ = 9_ )

#### GGA_C_RGE2(_ = 14_ )

#### GGA_C_SOGGA11(_ = 15_ )

#### GGA_C_SOGGA11_X(_ = 15_ )

#### GGA_C_SPBE(_ = 8_ )

#### GGA_C_TCA(_ = 10_ )

#### GGA_C_WI(_ = 14_ )

#### GGA_C_WI0(_ = 15_ )

#### GGA_C_WL(_ = 14_ )

#### GGA_C_XPBE(_ = 13_ )

#### GGA_C_ZPBEINT(_ = 6_ )

#### GGA_C_ZPBESOL(_ = 6_ )

#### GGA_K_ABSP1(_ = 50_ )

#### GGA_K_ABSP2(_ = 50_ )

#### GGA_K_APBE(_ = 18_ )

#### GGA_K_APBEINT(_ = 5_ )

#### GGA_K_BALTIN(_ = 50_ )

#### GGA_K_DK(_ = 51_ )

#### GGA_K_ERNZERHOF(_ = 52_ )

#### GGA_K_FR_B88(_ = 51_ )

#### GGA_K_FR_PW86(_ = 51_ )

#### GGA_K_GE2(_ = 50_ )

#### GGA_K_GOLDEN(_ = 50_ )

#### GGA_K_GP85(_ = 51_ )

#### GGA_K_GR(_ = 50_ )

#### GGA_K_LC94(_ = 52_ )

#### GGA_K_LIEB(_ = 50_ )

#### GGA_K_LLP(_ = 52_ )

#### GGA_K_LUDENA(_ = 50_ )

#### GGA_K_MEYER(_ = 5_ )

#### GGA_K_OL1(_ = 51_ )

#### GGA_K_OL2(_ = 51_ )

#### GGA_K_PEARSON(_ = 51_ )

#### GGA_K_PERDEW(_ = 51_ )

#### GGA_K_REVAPBE(_ = 5_ )

#### GGA_K_REVAPBEINT(_ = 5_ )

#### GGA_K_TFVW(_ = 5_ )

#### GGA_K_THAKKAR(_ = 52_ )

#### GGA_K_TW1(_ = 18_ )

#### GGA_K_TW2(_ = 18_ )

#### GGA_K_TW3(_ = 18_ )

#### GGA_K_TW4(_ = 19_ )

#### GGA_K_VJKS(_ = 51_ )

#### GGA_K_VSK(_ = 51_ )

#### GGA_K_VW(_ = 50_ )

#### GGA_K_YT65(_ = 50_ )

#### GGA_XC_B97_D(_ = 17_ )

#### GGA_XC_B97_GGA1(_ = 9_ )

#### GGA_XC_EDF1(_ = 16_ )

#### GGA_XC_HCTH_120(_ = 16_ )

#### GGA_XC_HCTH_147(_ = 16_ )

#### GGA_XC_HCTH_407(_ = 16_ )

#### GGA_XC_HCTH_407P(_ = 9_ )

#### GGA_XC_HCTH_93(_ = 16_ )

#### GGA_XC_HCTH_P14(_ = 9_ )

#### GGA_XC_HCTH_P76(_ = 9_ )

#### GGA_XC_KT2(_ = 14_ )

#### GGA_XC_MOHLYP(_ = 19_ )

#### GGA_XC_MOHLYP2(_ = 19_ )

#### GGA_XC_MPWLYP1W(_ = 17_ )

#### GGA_XC_OBLYP_D(_ = 6_ )

#### GGA_XC_OPBE_D(_ = 6_ )

#### GGA_XC_OPWLYP_D(_ = 6_ )

#### GGA_XC_PBE1W(_ = 17_ )

#### GGA_XC_PBELYP1W(_ = 17_ )

#### GGA_XC_TH1(_ = 15_ )

#### GGA_XC_TH2(_ = 15_ )

#### GGA_XC_TH3(_ = 15_ )

#### GGA_XC_TH4(_ = 15_ )

#### GGA_XC_TH_FC(_ = 19_ )

#### GGA_XC_TH_FCFO(_ = 19_ )

#### GGA_XC_TH_FCO(_ = 19_ )

#### GGA_XC_TH_FL(_ = 19_ )

#### GGA_XC_VV10(_ = 25_ )

#### GGA_XC_XLYP(_ = 16_ )

#### GGA_X_2D_B86(_ = 12_ )

#### GGA_X_2D_B86_MGC(_ = 12_ )

#### GGA_X_2D_B88(_ = 12_ )

#### GGA_X_2D_PBE(_ = 12_ )

#### GGA_X_AIRY(_ = 19_ )

#### GGA_X_AK13(_ = 5_ )

#### GGA_X_AM05(_ = 12_ )

#### GGA_X_APBE(_ = 18_ )

#### GGA_X_B86(_ = 10_ )

#### GGA_X_B86_MGC(_ = 10_ )

#### GGA_X_B86_R(_ = 4_ )

#### GGA_X_B88(_ = 10_ )

#### GGA_X_BAYESIAN(_ = 12_ )

#### GGA_X_BGCP(_ = 3_ )

#### GGA_X_BPCCAC(_ = 9_ )

#### GGA_X_C09X(_ = 15_ )

#### GGA_X_CAP(_ = 27_ )

#### GGA_X_DK87_R1(_ = 11_ )

#### GGA_X_DK87_R2(_ = 11_ )

#### GGA_X_EV93(_ = 3_ )

#### GGA_X_FT97_A(_ = 11_ )

#### GGA_X_FT97_B(_ = 11_ )

#### GGA_X_G96(_ = 10_ )

#### GGA_X_GAM(_ = 3_ )

#### GGA_X_HCTH_A(_ = 3_ )

#### GGA_X_HERMAN(_ = 10_ )

#### GGA_X_HJS_B88(_ = 52_ )

#### GGA_X_HJS_B88_V2(_ = 4_ )

#### GGA_X_HJS_B97X(_ = 52_ )

#### GGA_X_HJS_PBE(_ = 52_ )

#### GGA_X_HJS_PBE_SOL(_ = 52_ )

#### GGA_X_HTBS(_ = 19_ )

#### GGA_X_ITYH(_ = 52_ )

#### GGA_X_KT1(_ = 14_ )

#### GGA_X_LAG(_ = 19_ )

#### GGA_X_LAMBDA_CH_N(_ = 4_ )

#### GGA_X_LAMBDA_LO_N(_ = 4_ )

#### GGA_X_LAMBDA_OC2_N(_ = 4_ )

#### GGA_X_LB(_ = 16_ )

#### GGA_X_LBM(_ = 18_ )

#### GGA_X_LG93(_ = 11_ )

#### GGA_X_LV_RPW86(_ = 5_ )

#### GGA_X_MB88(_ = 14_ )

#### GGA_X_MPBE(_ = 12_ )

#### GGA_X_MPW91(_ = 11_ )

#### GGA_X_N12(_ = 8_ )

#### GGA_X_OL2(_ = 18_ )

#### GGA_X_OPTB88_VDW(_ = 13_ )

#### GGA_X_OPTPBE_VDW(_ = 14_ )

#### GGA_X_OPTX(_ = 11_ )

#### GGA_X_PBE(_ = 10_ )

#### GGA_X_PBEA(_ = 12_ )

#### GGA_X_PBEFE(_ = 26_ )

#### GGA_X_PBEINT(_ = 6_ )

#### GGA_X_PBEK1_VDW(_ = 14_ )

#### GGA_X_PBE_JSJR(_ = 12_ )

#### GGA_X_PBE_MOL(_ = 4_ )

#### GGA_X_PBE_R(_ = 10_ )

#### GGA_X_PBE_SOL(_ = 11_ )

#### GGA_X_PBE_TCA(_ = 5_ )

#### GGA_X_PW86(_ = 10_ )

#### GGA_X_PW91(_ = 10_ )

#### GGA_X_Q2D(_ = 4_ )

#### GGA_X_RGE2(_ = 14_ )

#### GGA_X_RPBE(_ = 11_ )

#### GGA_X_RPW86(_ = 14_ )

#### GGA_X_SFAT(_ = 53_ )

#### GGA_X_SOGGA(_ = 15_ )

#### GGA_X_SOGGA11(_ = 15_ )

#### GGA_X_SSB(_ = 9_ )

#### GGA_X_SSB_D(_ = 9_ )

#### GGA_X_SSB_SW(_ = 9_ )

#### GGA_X_VMT84_GE(_ = 6_ )

#### GGA_X_VMT84_PBE(_ = 6_ )

#### GGA_X_VMT_GE(_ = 7_ )

#### GGA_X_VMT_PBE(_ = 7_ )

#### GGA_X_WC(_ = 11_ )

#### GGA_X_WPBEH(_ = 52_ )

#### GGA_X_XPBE(_ = 12_ )

#### HYB_GGA_XC_B1LYP(_ = 41_ )

#### HYB_GGA_XC_B1PW91(_ = 41_ )

#### HYB_GGA_XC_B1WC(_ = 41_ )

#### HYB_GGA_XC_B3LYP(_ = 40_ )

#### HYB_GGA_XC_B3LYP5(_ = 47_ )

#### HYB_GGA_XC_B3LYPs(_ = 45_ )

#### HYB_GGA_XC_B3P86(_ = 40_ )

#### HYB_GGA_XC_B3PW91(_ = 40_ )

#### HYB_GGA_XC_B97(_ = 40_ )

#### HYB_GGA_XC_B97_1(_ = 40_ )

#### HYB_GGA_XC_B97_1p(_ = 26_ )

#### HYB_GGA_XC_B97_2(_ = 41_ )

#### HYB_GGA_XC_B97_3(_ = 41_ )

#### HYB_GGA_XC_B97_K(_ = 41_ )

#### HYB_GGA_XC_BHANDH(_ = 43_ )

#### HYB_GGA_XC_BHANDHLYP(_ = 43_ )

#### HYB_GGA_XC_CAMY_B3LYP(_ = 47_ )

#### HYB_GGA_XC_CAMY_BLYP(_ = 45_ )

#### HYB_GGA_XC_CAM_B3LYP(_ = 43_ )

#### HYB_GGA_XC_CAP0(_ = 47_ )

#### HYB_GGA_XC_EDF2(_ = 47_ )

#### HYB_GGA_XC_HJS_B88(_ = 43_ )

#### HYB_GGA_XC_HJS_B97X(_ = 43_ )

#### HYB_GGA_XC_HJS_PBE(_ = 42_ )

#### HYB_GGA_XC_HJS_PBE_SOL(_ = 43_ )

#### HYB_GGA_XC_HPBEINT(_ = 47_ )

#### HYB_GGA_XC_HSE03(_ = 42_ )

#### HYB_GGA_XC_HSE06(_ = 42_ )

#### HYB_GGA_XC_LCY_BLYP(_ = 46_ )

#### HYB_GGA_XC_LCY_PBE(_ = 46_ )

#### HYB_GGA_XC_LC_VV10(_ = 46_ )

#### HYB_GGA_XC_LRC_WPBE(_ = 47_ )

#### HYB_GGA_XC_LRC_WPBEH(_ = 46_ )

#### HYB_GGA_XC_MB3LYP_RC04(_ = 43_ )

#### HYB_GGA_XC_MPW3LYP(_ = 41_ )

#### HYB_GGA_XC_MPW3PW(_ = 41_ )

#### HYB_GGA_XC_MPWLYP1M(_ = 45_ )

#### HYB_GGA_XC_O3LYP(_ = 40_ )

#### HYB_GGA_XC_PBE0_13(_ = 45_ )

#### HYB_GGA_XC_PBEH(_ = 40_ )

#### HYB_GGA_XC_REVB3LYP(_ = 45_ )

#### HYB_GGA_XC_SB98_1a(_ = 42_ )

#### HYB_GGA_XC_SB98_1b(_ = 42_ )

#### HYB_GGA_XC_SB98_1c(_ = 42_ )

#### HYB_GGA_XC_SB98_2a(_ = 42_ )

#### HYB_GGA_XC_SB98_2b(_ = 42_ )

#### HYB_GGA_XC_SB98_2c(_ = 42_ )

#### HYB_GGA_XC_TUNED_CAM_B3LYP(_ = 43_ )

#### HYB_GGA_XC_WB97(_ = 46_ )

#### HYB_GGA_XC_WB97X(_ = 46_ )

#### HYB_GGA_XC_WB97X_D(_ = 47_ )

#### HYB_GGA_XC_WB97X_V(_ = 46_ )

#### HYB_GGA_XC_X3LYP(_ = 41_ )

#### HYB_GGA_XC_mPW1K(_ = 40_ )

#### HYB_GGA_XC_mPW1PW(_ = 41_ )

#### HYB_GGA_X_N12_SX(_ = 8_ )

#### HYB_GGA_X_SOGGA11_X(_ = 42_ )

#### HYB_MGGA_XC_B86B95(_ = 44_ )

#### HYB_MGGA_XC_B88B95(_ = 44_ )

#### HYB_MGGA_XC_BB1K(_ = 44_ )

#### HYB_MGGA_XC_M05(_ = 43_ )

#### HYB_MGGA_XC_M05_2X(_ = 43_ )

#### HYB_MGGA_XC_M06(_ = 44_ )

#### HYB_MGGA_XC_M06_2X(_ = 45_ )

#### HYB_MGGA_XC_M06_HF(_ = 44_ )

#### HYB_MGGA_XC_M08_HX(_ = 46_ )

#### HYB_MGGA_XC_M08_SO(_ = 46_ )

#### HYB_MGGA_XC_M11(_ = 46_ )

#### HYB_MGGA_XC_MPW1B95(_ = 44_ )

#### HYB_MGGA_XC_MPWB1K(_ = 44_ )

#### HYB_MGGA_XC_PW6B95(_ = 45_ )

#### HYB_MGGA_XC_PW86B95(_ = 44_ )

#### HYB_MGGA_XC_PWB6K(_ = 45_ )

#### HYB_MGGA_XC_REVTPSSH(_ = 45_ )

#### HYB_MGGA_XC_TPSSH(_ = 45_ )

#### HYB_MGGA_XC_WB97M_V(_ = 53_ )

#### HYB_MGGA_XC_X1B95(_ = 44_ )

#### HYB_MGGA_XC_XB1K(_ = 44_ )

#### HYB_MGGA_X_DLDF(_ = 3_ )

#### HYB_MGGA_X_MN12_SX(_ = 24_ )

#### HYB_MGGA_X_MN15(_ = 26_ )

#### HYB_MGGA_X_MS2H(_ = 22_ )

#### HYB_MGGA_X_MVSH(_ = 47_ )

#### HYB_MGGA_X_SCAN0(_ = 26_ )

#### LDA_C_1D_CSC(_ = 1_ )

#### LDA_C_1D_LOOS(_ = 2_ )

#### LDA_C_2D_AMGB(_ = 1_ )

#### LDA_C_2D_PRM(_ = 1_ )

#### LDA_C_GL(_ = _ )

#### LDA_C_GOMBAS(_ = 2_ )

#### LDA_C_HL(_ = _ )

#### LDA_C_ML1(_ = 2_ )

#### LDA_C_ML2(_ = 2_ )

#### LDA_C_OB_PW(_ = 1_ )

#### LDA_C_OB_PZ(_ = 1_ )

#### LDA_C_PW(_ = 1_ )

#### LDA_C_PW_MOD(_ = 1_ )

#### LDA_C_PW_RPA(_ = 2_ )

#### LDA_C_PZ(_ = _ )

#### LDA_C_PZ_MOD(_ = 1_ )

#### LDA_C_RC04(_ = 2_ )

#### LDA_C_RPA(_ = _ )

#### LDA_C_VWN(_ = _ )

#### LDA_C_VWN_1(_ = 2_ )

#### LDA_C_VWN_2(_ = 2_ )

#### LDA_C_VWN_3(_ = 3_ )

#### LDA_C_VWN_4(_ = 3_ )

#### LDA_C_VWN_RPA(_ = _ )

#### LDA_C_WIGNER(_ = _ )

#### LDA_C_XALPHA(_ = _ )

#### LDA_C_vBH(_ = 1_ )

#### LDA_K_LP(_ = 5_ )

#### LDA_K_TF(_ = 5_ )

#### LDA_X(_ = _ )

#### LDA_XC_KSDT(_ = 25_ )

#### LDA_XC_TETER93(_ = 2_ )

#### LDA_XC_ZLP(_ = 4_ )

#### LDA_X_1D(_ = 2_ )

#### LDA_X_2D(_ = 1_ )

#### MGGA_C_BC95(_ = 24_ )

#### MGGA_C_CC06(_ = 22_ )

#### MGGA_C_CS(_ = 7_ )

#### MGGA_C_DLDF(_ = 3_ )

#### MGGA_C_M05(_ = 23_ )

#### MGGA_C_M05_2X(_ = 23_ )

#### MGGA_C_M06(_ = 23_ )

#### MGGA_C_M06_2X(_ = 23_ )

#### MGGA_C_M06_HF(_ = 23_ )

#### MGGA_C_M06_L(_ = 23_ )

#### MGGA_C_M08_HX(_ = 7_ )

#### MGGA_C_M08_SO(_ = 7_ )

#### MGGA_C_M11(_ = 7_ )

#### MGGA_C_M11_L(_ = 7_ )

#### MGGA_C_MN12_L(_ = 7_ )

#### MGGA_C_MN12_SX(_ = 7_ )

#### MGGA_C_MN15(_ = 26_ )

#### MGGA_C_MN15_L(_ = 26_ )

#### MGGA_C_PKZB(_ = 23_ )

#### MGGA_C_REVTPSS(_ = 24_ )

#### MGGA_C_SCAN(_ = 26_ )

#### MGGA_C_TPSS(_ = 23_ )

#### MGGA_C_TPSSLOC(_ = 24_ )

#### MGGA_C_VSXC(_ = 23_ )

#### MGGA_XC_B97M_V(_ = 25_ )

#### MGGA_XC_OTPSS_D(_ = 6_ )

#### MGGA_XC_TPSSLYP1W(_ = 24_ )

#### MGGA_XC_ZLP(_ = 4_ )

#### MGGA_X_2D_PRHG07(_ = 21_ )

#### MGGA_X_2D_PRHG07_PRP10(_ = 21_ )

#### MGGA_X_BJ06(_ = 20_ )

#### MGGA_X_BLOC(_ = 24_ )

#### MGGA_X_BR89(_ = 20_ )

#### MGGA_X_GVT4(_ = 20_ )

#### MGGA_X_LTA(_ = 20_ )

#### MGGA_X_M05(_ = 21_ )

#### MGGA_X_M05_2X(_ = 21_ )

#### MGGA_X_M06(_ = 21_ )

#### MGGA_X_M06_2X(_ = 21_ )

#### MGGA_X_M06_HF(_ = 21_ )

#### MGGA_X_M06_L(_ = 20_ )

#### MGGA_X_M08_HX(_ = 21_ )

#### MGGA_X_M08_SO(_ = 22_ )

#### MGGA_X_M11(_ = 22_ )

#### MGGA_X_M11_L(_ = 22_ )

#### MGGA_X_MBEEF(_ = 24_ )

#### MGGA_X_MBEEFVDW(_ = 25_ )

#### MGGA_X_MK00(_ = 23_ )

#### MGGA_X_MK00B(_ = 24_ )

#### MGGA_X_MN12_L(_ = 22_ )

#### MGGA_X_MN15_L(_ = 26_ )

#### MGGA_X_MODTPSS(_ = 24_ )

#### MGGA_X_MS0(_ = 22_ )

#### MGGA_X_MS1(_ = 22_ )

#### MGGA_X_MS2(_ = 22_ )

#### MGGA_X_MVS(_ = 25_ )

#### MGGA_X_PKZB(_ = 21_ )

#### MGGA_X_REVTPSS(_ = 21_ )

#### MGGA_X_RPP09(_ = 20_ )

#### MGGA_X_SCAN(_ = 26_ )

#### MGGA_X_TAU_HCTH(_ = 20_ )

#### MGGA_X_TB09(_ = 20_ )

#### MGGA_X_TPSS(_ = 20_ )

#### _static_ all_families()
List of strings with the libxc families.
Note that XC_FAMILY if removed from the string e.g. XC_FAMILY_LDA becomes LDA.


#### _static_ all_kinds()
List of strings with the libxc kinds.
Also in this case, the string is obtained by remove the

```
XC_
```

 prefix.
XC_CORRELATION –> CORRELATION.


#### as_dict()
Makes LibxcFunc obey the general json interface used in pymatgen for
easier serialization.


#### _classmethod_ from_dict(dct)
Makes LibxcFunc obey the general json interface used in pymatgen for
easier serialization.


#### _property_ info_dict()
Dictionary with metadata. see libxc_docs.json.


#### _property_ is_c_kind(_: boo_ )
True if this is a correlation-only functional.


#### _property_ is_gga_family(_: boo_ )
True if this functional belongs to the GGA family.


#### _property_ is_hyb_gga_family(_: boo_ )
True if this functional belongs to the hybrid + GGA family.


#### _property_ is_hyb_mgga_family(_: boo_ )
True if this functional belongs to the hybrid + meta-GGA family.


#### _property_ is_k_kind(_: boo_ )
True if this is a kinetic functional.


#### _property_ is_lda_family(_: boo_ )
True if this functional belongs to the LDA family.


#### _property_ is_mgga_family(_: boo_ )
True if this functional belongs to the meta-GGA family.


#### _property_ is_x_kind(_: boo_ )
True if this is an exchange-only functional.


#### _property_ is_xc_kind(_: boo_ )
True if this is a exchange+correlation functional.


#### to_json()
Returns a json string representation of the MSONable object.

## pymatgen.core.molecular_orbitals module

This module implements a MolecularOrbital class to represent band character in
solids. Useful for predicting PDOS character from structural information.


### _class_ MolecularOrbitals(formula)
Bases: `object`

Represents the character of bands in a solid. The input is a chemical
formula, since no structural characteristics are taken into account.

The band character of a crystal emerges from the atomic orbitals of the
constituent ions, hybridization/covalent bonds, and the spin-orbit
interaction (ex: Fe2O3). Right now the orbitals are only built from
the uncharged atomic species. Functionality can be improved by:
1) calculate charged ion orbital energies
2) incorporate the coordination environment to account for covalent bonds

The atomic orbital energies are stored in pymatgen.core.periodic_table.JSON

MOs = MolecularOrbitals(‘SrTiO3’)
MOs.band_edges
# gives {‘HOMO’:[‘O’,’2p’,-0.338381], ‘LUMO’:[‘Ti’,’3d’,-0.17001], ‘metal’:False}


* **Parameters**

    **formula** (*str*) – Chemical formula. Must have integer subscripts. Ex: ‘SrTiO3’.



#### composition()
the composition as a dictionary. Ex: {‘Sr’: 1, ‘Ti’: 1, ‘O’, 3}


#### elements()
the dictionary keys for the composition


#### elec_neg()
the maximum pairwise electronegativity difference


#### aos()
the constituent atomic orbitals for each element as a dictionary


#### band_edges()
dictionary containing the highest occupied molecular orbital (HOMO),
lowest unoccupied molecular orbital (LUMO), and whether the material is predicted
to be a metal


#### aos_as_list()
The orbitals energies in eV are represented as
[[‘O’, ‘1s’, -18.758245], [‘O’, ‘2s’, -0.871362], [‘O’, ‘2p’, -0.338381]]
Data is obtained from
[https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations](https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations).


* **Returns**

    A list of atomic orbitals, sorted from lowest to highest energy.



#### max_electronegativity()

* **Returns**

    The maximum pairwise electronegativity difference.



#### obtain_band_edges()
Fill up the atomic orbitals with available electrons.


* **Returns**

    HOMO, LUMO, and whether it’s a metal.


## pymatgen.core.operations module

This module provides classes that operate on points or vectors in 3D space.


### _class_ MagSymmOp(affine_transformation_matrix: ArrayLike, time_reversal: int, tol: float = 0.01)
Bases: `SymmOp`

Thin wrapper around SymmOp to extend it to support magnetic symmetry by including a time
reversal operator. Magnetic symmetry is similar to conventional crystal symmetry, except
symmetry is reduced by the addition of a time reversal operator which acts on an atom’s magnetic
moment.

Initializes the MagSymmOp from a 4x4 affine transformation matrix and time reversal
operator. In general, this constructor should not be used unless you are transferring
rotations. Use the static constructors instead to generate a SymmOp from proper rotations
and translation.


* **Parameters**


    * **affine_transformation_matrix** (*4x4 array*) – Representing an
    affine transformation.


    * **time_reversal** (*int*) – 1 or -1


    * **tol** (*float*) – Tolerance for determining if matrices are equal.



#### as_dict()
MSONable dict.


#### as_xyzt_str()
Returns a string of the form ‘x, y, z, +1’, ‘-x, -y, z, -1’,
‘-y+1/2, x+1/2, z+1/2, +1’, etc. Only works for integer rotation matrices.


#### as_xyzt_string(\*\*kwds)
as_xyzt_string is deprecated!
Use as_xyzt_str instead


#### _classmethod_ from_dict(d: dict)

* **Parameters**

    **d** – dict



* **Returns**

    MagneticSymmOp from dict representation.



#### _static_ from_rotation_and_translation_and_time_reversal(rotation_matrix: ArrayLike = ((1, 0, 0), (0, 1, 0), (0, 0, 1)), translation_vec: ArrayLike = (0, 0, 0), time_reversal: int = 1, tol: float = 0.1)
Creates a symmetry operation from a rotation matrix, translation
vector and time reversal operator.


* **Parameters**


    * **rotation_matrix** (*3x3 array*) – Rotation matrix.


    * **translation_vec** (*3x1 array*) – Translation vector.


    * **time_reversal** (*int*) – Time reversal operator, +1 or -1.


    * **tol** (*float*) – Tolerance to determine if rotation matrix is valid.



* **Returns**

    MagSymmOp object



#### _classmethod_ from_symmop(symmop: SymmOp, time_reversal)
Initialize a MagSymmOp from a SymmOp and time reversal operator.


* **Parameters**


    * **symmop** (*SymmOp*) – SymmOp


    * **time_reversal** (*int*) – Time reversal operator, +1 or -1.



* **Returns**

    MagSymmOp object



#### _classmethod_ from_xyzt_str(xyzt_string: str)

* **Parameters**

    **xyzt_string** (*str*) – of the form ‘x, y, z, +1’, ‘-x, -y, z, -1’,
    ‘-2y+1/2, 3x+1/2, z-y+1/2, +1’, etc.



* **Returns**

    MagSymmOp object



#### _classmethod_ from_xyzt_string(\*args, \*\*kwds)
from_xyzt_string is deprecated!
Use from_xyzt_str instead


#### operate_magmom(magmom)
Apply time reversal operator on the magnetic moment. Note that
magnetic moments transform as axial vectors, not polar vectors.

See ‘Symmetry and magnetic structures’, Rodríguez-Carvajal and
Bourée for a good discussion. DOI: 10.1051/epjconf/20122200010


* **Parameters**


    * **magmom** – Magnetic moment as electronic_structure.core.Magmom


    * **array-like** (*class** or **as list** or **np*) –



* **Returns**

    Magnetic moment after operator applied as Magmom class



### _class_ SymmOp(affine_transformation_matrix: ArrayLike, tol: float = 0.01)
Bases: `MSONable`

A symmetry operation in Cartesian space. Consists of a rotation plus a
translation. Implementation is as an affine transformation matrix of rank 4
for efficiency. Read: [http://en.wikipedia.org/wiki/Affine_transformation](http://en.wikipedia.org/wiki/Affine_transformation).


#### affine_matrix()
A 4x4 array representing the symmetry operation.


* **Type**

    np.ndarray


Initializes the SymmOp from a 4x4 affine transformation matrix.
In general, this constructor should not be used unless you are
transferring rotations. Use the static constructors instead to
generate a SymmOp from proper rotations and translation.


* **Parameters**


    * **affine_transformation_matrix** (*4x4 array*) – Representing an
    affine transformation.


    * **tol** (*float*) – Tolerance for determining if matrices are equal. Defaults to 0.01.



* **Raises**

    **ValueError** – if matrix is not 4x4.



#### apply_rotation_only(vector: ArrayLike)
Vectors should only be operated by the rotation matrix and not the
translation vector.


* **Parameters**

    **vector** (*3x1 array*) – A vector.



#### are_symmetrically_related(point_a: ArrayLike, point_b: ArrayLike, tol: float = 0.001)
Checks if two points are symmetrically related.


* **Parameters**


    * **point_a** (*3x1 array*) – First point.


    * **point_b** (*3x1 array*) – Second point.


    * **tol** (*float*) – Absolute tolerance for checking distance. Defaults to 0.001.



* **Returns**

    True if self.operate(point_a) == point_b or vice versa.



* **Return type**

    bool



#### are_symmetrically_related_vectors(from_a: ArrayLike, to_a: ArrayLike, r_a: ArrayLike, from_b: ArrayLike, to_b: ArrayLike, r_b: ArrayLike, tol: float = 0.001)
Checks if two vectors, or rather two vectors that connect two points
each are symmetrically related. r_a and r_b give the change of unit
cells. Two vectors are also considered symmetrically equivalent if starting
and end point are exchanged.


* **Parameters**


    * **from_a** (*3x1 array*) – Starting point of the first vector.


    * **to_a** (*3x1 array*) – Ending point of the first vector.


    * **from_b** (*3x1 array*) – Starting point of the second vector.


    * **to_b** (*3x1 array*) – Ending point of the second vector.


    * **r_a** (*3x1 array*) – Change of unit cell of the first vector.


    * **r_b** (*3x1 array*) – Change of unit cell of the second vector.


    * **tol** (*float*) – Absolute tolerance for checking distance.



* **Returns**

    (are_related, is_reversed)



#### as_dict()
MSONable dict.


#### as_xyz_str()
Returns a string of the form ‘x, y, z’, ‘-x, -y, z’, ‘-y+1/2, x+1/2, z+1/2’, etc.
Only works for integer rotation matrices.


#### as_xyz_string(\*\*kwds)
as_xyz_string is deprecated!
Use as_xyz_str instead


#### _static_ from_axis_angle_and_translation(axis: ArrayLike, angle: float, angle_in_radians: bool = False, translation_vec: ArrayLike = (0, 0, 0))
Generates a SymmOp for a rotation about a given axis plus translation.


* **Parameters**


    * **axis** – The axis of rotation in Cartesian space. For example,
    [1, 0, 0]indicates rotation about x-axis.


    * **angle** (*float*) – Angle of rotation.


    * **angle_in_radians** (*bool*) – Set to True if angles are given in
    radians. Or else, units of degrees are assumed.


    * **translation_vec** – A translation vector. Defaults to zero.



* **Returns**

    SymmOp for a rotation about given axis and translation.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    SymmOp from dict representation.



#### _static_ from_origin_axis_angle(origin: ArrayLike, axis: ArrayLike, angle: float, angle_in_radians: bool = False)
Generates a SymmOp for a rotation about a given axis through an
origin.


* **Parameters**


    * **origin** (*3x1 array*) – The origin which the axis passes through.


    * **axis** (*3x1 array*) – The axis of rotation in Cartesian space. For
    example, [1, 0, 0]indicates rotation about x-axis.


    * **angle** (*float*) – Angle of rotation.


    * **angle_in_radians** (*bool*) – Set to True if angles are given in
    radians. Or else, units of degrees are assumed.



* **Returns**

    SymmOp.



#### _static_ from_rotation_and_translation(rotation_matrix: ArrayLike = ((1, 0, 0), (0, 1, 0), (0, 0, 1)), translation_vec: ArrayLike = (0, 0, 0), tol: float = 0.1)
Creates a symmetry operation from a rotation matrix and a translation
vector.


* **Parameters**


    * **rotation_matrix** (*3x3 array*) – Rotation matrix.


    * **translation_vec** (*3x1 array*) – Translation vector.


    * **tol** (*float*) – Tolerance to determine if rotation matrix is valid.



* **Returns**

    SymmOp object



#### _classmethod_ from_xyz_str(xyz_str: str)

* **Parameters**

    **xyz_str** – string of the form ‘x, y, z’, ‘-x, -y, z’, ‘-2y+1/2, 3x+1/2, z-y+1/2’, etc.



* **Returns**

    SymmOp



#### _classmethod_ from_xyz_string(\*args, \*\*kwds)
from_xyz_string is deprecated!
Use from_xyz_str instead


#### _property_ inverse(_: SymmO_ )
Returns inverse of transformation.


#### _static_ inversion(origin: ArrayLike = (0, 0, 0))
Inversion symmetry operation about axis.


* **Parameters**

    **origin** (*3x1 array*) – Origin of the inversion operation. Defaults
    to [0, 0, 0].



* **Returns**

    SymmOp representing an inversion operation about the origin.



#### operate(point: ArrayLike)
Apply the operation on a point.


* **Parameters**

    **point** – Cartesian coordinate.



* **Returns**

    Coordinates of point after operation.



#### operate_multi(points: ArrayLike)
Apply the operation on a list of points.


* **Parameters**

    **points** – List of Cartesian coordinates



* **Returns**

    Numpy array of coordinates after operation



#### _static_ reflection(normal: ArrayLike, origin: ArrayLike = (0, 0, 0))
Returns reflection symmetry operation.


* **Parameters**


    * **normal** (*3x1 array*) – Vector of the normal to the plane of
    reflection.


    * **origin** (*3x1 array*) – A point in which the mirror plane passes
    through.



* **Returns**

    SymmOp for the reflection about the plane



#### _property_ rotation_matrix(_: ndarra_ )
A 3x3 numpy.array representing the rotation matrix.


#### _static_ rotoreflection(axis: ArrayLike, angle: float, origin: ArrayLike = (0, 0, 0))
Returns a roto-reflection symmetry operation.


* **Parameters**


    * **axis** (*3x1 array*) – Axis of rotation / mirror normal


    * **angle** (*float*) – Angle in degrees


    * **origin** (*3x1 array*) – Point left invariant by roto-reflection.
    Defaults to (0, 0, 0).



* **Returns**

    Roto-reflection operation



#### transform_tensor(tensor: ndarray)
Applies rotation portion to a tensor. Note that tensor has to be in
full form, not the Voigt form.


* **Parameters**

    **tensor** (*numpy array*) – A rank n tensor



* **Returns**

    Transformed tensor.



#### _property_ translation_vector(_: ndarra_ )
A rank 1 numpy.array of dim 3 representing the translation vector.

## pymatgen.core.periodic_table module

Classes representing Element, Species (Element + oxidation state) and PeriodicTable.


### _class_ DummySpecie(symbol: str = 'X', oxidation_state: float | None = 0, spin: float | None = None)
Bases: `DummySpecies`

This maps the historical grammatically inaccurate DummySpecie to DummySpecies
to maintain backwards compatibility.


* **Parameters**


    * **symbol** (*str*) – An assigned symbol for the dummy specie. Strict
    rules are applied to the choice of the symbol. The dummy
    symbol cannot have any part of first two letters that will
    constitute an Element symbol. Otherwise, a composition may
    be parsed wrongly. E.g., “X” is fine, but “Vac” is not
    because Vac contains V, a valid Element.


    * **oxidation_state** (*float*) – Oxidation state for dummy specie. Defaults to 0.
    deprecated and retained purely for backward compatibility.


    * **spin** – Spin associated with Species. Defaults to None.



### _class_ DummySpecies(symbol: str = 'X', oxidation_state: float | None = 0, spin: float | None = None)
Bases: `Species`

A special specie for representing non-traditional elements or species. For
example, representation of vacancies (charged or otherwise), or special
sites, etc.


#### oxi_state()
Oxidation state associated with Species.


* **Type**

    int



#### Z()
DummySpecies is always assigned an atomic number equal to the hash
number of the symbol. Obviously, it makes no sense whatsoever to use
the atomic number of a Dummy specie for anything scientific. The purpose
of this is to ensure that for most use cases, a DummySpecies behaves no
differently from an Element or Species.


* **Type**

    int



#### X()
DummySpecies is always assigned a Pauling electronegativity of 0.


* **Type**

    float



* **Parameters**


    * **symbol** (*str*) – An assigned symbol for the dummy specie. Strict
    rules are applied to the choice of the symbol. The dummy
    symbol cannot have any part of first two letters that will
    constitute an Element symbol. Otherwise, a composition may
    be parsed wrongly. E.g., “X” is fine, but “Vac” is not
    because Vac contains V, a valid Element.


    * **oxidation_state** (*float*) – Oxidation state for dummy specie. Defaults to 0.
    deprecated and retained purely for backward compatibility.


    * **spin** – Spin associated with Species. Defaults to None.



#### _property_ X(_: floa_ )
DummySpecies is always assigned a Pauling electronegativity of 0. The effect of
this is that DummySpecies are always sorted in front of actual Species.


#### _property_ Z(_: in_ )
DummySpecies is always assigned an atomic number equal to the hash of
the symbol. The expectation is that someone would be an actual dummy
to use atomic numbers for a Dummy specie.


#### as_dict()
MSONable dict representation.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    DummySpecies



#### _static_ from_str(species_string: str)
Returns a Dummy from a string representation.


* **Parameters**

    **species_string** (*str*) – A string representation of a dummy
    species, e.g., “X2+”, “X3+”.



* **Returns**

    A DummySpecies object.



* **Raises**

    **ValueError if species_string cannot be interpreted.** –



#### _property_ oxi_state(_: float | Non_ )
Oxidation state associated with DummySpecies.


#### _property_ symbol(_: st_ )
Symbol for DummySpecies.


### _class_ Element(value)
Bases: `ElementBase`

Enum representing an element in the periodic table.

Basic immutable element object with all relevant properties.

Only one instance of Element for each symbol is stored after creation,
ensuring that a particular element behaves like a singleton. For all
attributes, missing data (i.e., data for which is not available) is
represented by a None unless otherwise stated.


* **Parameters**

    **symbol** (*str*) – Element symbol, e.g., “H”, “Fe”



#### Z()
Atomic number.


* **Type**

    int



#### symbol()
Element symbol.


* **Type**

    str



#### long_name()
Long name for element. E.g., “Hydrogen”.


* **Type**

    str



#### atomic_radius_calculated()
Calculated atomic radius for the element. This is the empirical value.
Data is obtained from [http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page](http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)).


* **Type**

    float



#### van_der_waals_radius()
Van der Waals radius for the element. This is the empirical value determined
from critical reviews of X-ray diffraction, gas kinetic collision cross-section, and other experimental
data by Bondi and later workers. The uncertainty in these values is on the order of 0.1 Å.
Data are obtained from “Atomic Radii of the Elements” in CRC Handbook of Chemistry and Physics,
91st Ed.; Haynes, W.M., Ed.; CRC Press: Boca Raton, FL, 2010.


* **Type**

    float



#### mendeleev_no()
Mendeleev number from definition given by Pettifor, D. G. (1984). A chemical scale
for crystal-structure maps. Solid State Communications, 51 (1), 31-34.


* **Type**

    int



#### electrical_resistivity()
Electrical resistivity.


* **Type**

    float



#### velocity_of_sound()
Velocity of sound.


* **Type**

    float



#### reflectivity()
Reflectivity.


* **Type**

    float



#### refractive_index()
Refractive index.


* **Type**

    float



#### poissons_ratio()
Poisson’s ratio.


* **Type**

    float



#### molar_volume()
Molar volume.


* **Type**

    float



#### electronic_structure()
Electronic structure. E.g., The electronic structure for Fe is represented
as [Ar].3d6.4s2.


* **Type**

    str



#### atomic_orbitals()
Atomic Orbitals. Energy of the atomic orbitals as a dict. E.g., The orbitals
energies in eV are represented as {‘1s’: -1.0, ‘2s’: -0.1}. Data is obtained from
[https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations](https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations).
The LDA values for neutral atoms are used.


* **Type**

    dict



#### thermal_conductivity()
Thermal conductivity.


* **Type**

    float



#### boiling_point()
Boiling point.


* **Type**

    float



#### melting_point()
Melting point.


* **Type**

    float



#### critical_temperature()
Critical temperature.


* **Type**

    float



#### superconduction_temperature()
Superconduction temperature.


* **Type**

    float



#### liquid_range()
Liquid range.


* **Type**

    float



#### bulk_modulus()
Bulk modulus.


* **Type**

    float



#### youngs_modulus()
Young’s modulus.


* **Type**

    float



#### brinell_hardness()
Brinell hardness.


* **Type**

    float



#### rigidity_modulus()
Rigidity modulus.


* **Type**

    float



#### mineral_hardness()
Mineral hardness.


* **Type**

    float



#### vickers_hardness()
Vicker’s hardness.


* **Type**

    float



#### density_of_solid()
Density of solid phase.


* **Type**

    float



#### coefficient_of_linear_thermal_expansion()
Coefficient of linear thermal expansion.


* **Type**

    float



#### ground_level()
Ground level for element.


* **Type**

    float



#### ionization_energies()
List of ionization energies. First value is the first
ionization energy, second is the second ionization energy, etc. Note that this is zero-based indexing!
So Element.ionization_energies[0] refer to the 1st ionization energy. Values are from the NIST Atomic
Spectra Database. Missing values are None.


* **Type**

    list[Optional[float]]



#### Ac(_ = 'Ac_ )

#### Ag(_ = 'Ag_ )

#### Al(_ = 'Al_ )

#### Am(_ = 'Am_ )

#### Ar(_ = 'Ar_ )

#### As(_ = 'As_ )

#### At(_ = 'At_ )

#### Au(_ = 'Au_ )

#### B(_ = 'B_ )

#### Ba(_ = 'Ba_ )

#### Be(_ = 'Be_ )

#### Bh(_ = 'Bh_ )

#### Bi(_ = 'Bi_ )

#### Bk(_ = 'Bk_ )

#### Br(_ = 'Br_ )

#### C(_ = 'C_ )

#### Ca(_ = 'Ca_ )

#### Cd(_ = 'Cd_ )

#### Ce(_ = 'Ce_ )

#### Cf(_ = 'Cf_ )

#### Cl(_ = 'Cl_ )

#### Cm(_ = 'Cm_ )

#### Cn(_ = 'Cn_ )

#### Co(_ = 'Co_ )

#### Cr(_ = 'Cr_ )

#### Cs(_ = 'Cs_ )

#### Cu(_ = 'Cu_ )

#### Db(_ = 'Db_ )

#### Ds(_ = 'Ds_ )

#### Dy(_ = 'Dy_ )

#### Er(_ = 'Er_ )

#### Es(_ = 'Es_ )

#### Eu(_ = 'Eu_ )

#### F(_ = 'F_ )

#### Fe(_ = 'Fe_ )

#### Fl(_ = 'Fl_ )

#### Fm(_ = 'Fm_ )

#### Fr(_ = 'Fr_ )

#### Ga(_ = 'Ga_ )

#### Gd(_ = 'Gd_ )

#### Ge(_ = 'Ge_ )

#### H(_ = 'H_ )

#### He(_ = 'He_ )

#### Hf(_ = 'Hf_ )

#### Hg(_ = 'Hg_ )

#### Ho(_ = 'Ho_ )

#### Hs(_ = 'Hs_ )

#### I(_ = 'I_ )

#### In(_ = 'In_ )

#### Ir(_ = 'Ir_ )

#### K(_ = 'K_ )

#### Kr(_ = 'Kr_ )

#### La(_ = 'La_ )

#### Li(_ = 'Li_ )

#### Lr(_ = 'Lr_ )

#### Lu(_ = 'Lu_ )

#### Lv(_ = 'Lv_ )

#### Mc(_ = 'Mc_ )

#### Md(_ = 'Md_ )

#### Mg(_ = 'Mg_ )

#### Mn(_ = 'Mn_ )

#### Mo(_ = 'Mo_ )

#### Mt(_ = 'Mt_ )

#### N(_ = 'N_ )

#### Na(_ = 'Na_ )

#### Nb(_ = 'Nb_ )

#### Nd(_ = 'Nd_ )

#### Ne(_ = 'Ne_ )

#### Nh(_ = 'Nh_ )

#### Ni(_ = 'Ni_ )

#### No(_ = 'No_ )

#### Np(_ = 'Np_ )

#### O(_ = 'O_ )

#### Og(_ = 'Og_ )

#### Os(_ = 'Os_ )

#### P(_ = 'P_ )

#### Pa(_ = 'Pa_ )

#### Pb(_ = 'Pb_ )

#### Pd(_ = 'Pd_ )

#### Pm(_ = 'Pm_ )

#### Po(_ = 'Po_ )

#### Pr(_ = 'Pr_ )

#### Pt(_ = 'Pt_ )

#### Pu(_ = 'Pu_ )

#### Ra(_ = 'Ra_ )

#### Rb(_ = 'Rb_ )

#### Re(_ = 'Re_ )

#### Rf(_ = 'Rf_ )

#### Rg(_ = 'Rg_ )

#### Rh(_ = 'Rh_ )

#### Rn(_ = 'Rn_ )

#### Ru(_ = 'Ru_ )

#### S(_ = 'S_ )

#### Sb(_ = 'Sb_ )

#### Sc(_ = 'Sc_ )

#### Se(_ = 'Se_ )

#### Sg(_ = 'Sg_ )

#### Si(_ = 'Si_ )

#### Sm(_ = 'Sm_ )

#### Sn(_ = 'Sn_ )

#### Sr(_ = 'Sr_ )

#### Ta(_ = 'Ta_ )

#### Tb(_ = 'Tb_ )

#### Tc(_ = 'Tc_ )

#### Te(_ = 'Te_ )

#### Th(_ = 'Th_ )

#### Ti(_ = 'Ti_ )

#### Tl(_ = 'Tl_ )

#### Tm(_ = 'Tm_ )

#### Ts(_ = 'Ts_ )

#### U(_ = 'U_ )

#### V(_ = 'V_ )

#### W(_ = 'W_ )

#### Xe(_ = 'Xe_ )

#### Y(_ = 'Y_ )

#### Yb(_ = 'Yb_ )

#### Zn(_ = 'Zn_ )

#### Zr(_ = 'Zr_ )

### _class_ ElementBase(value)
Bases: `Enum`

Element class defined without any enum values so it can be subclassed.

This class is needed to get nested (as|from)_dict to work properly. All emmet classes that had
Element classes required custom construction whereas this definition behaves more like dataclasses
so serialization is less troublesome. There were many times where objects in as_dict serialized
only when they were top level. See [https://github.com/materialsproject/pymatgen/issues/2999](https://github.com/materialsproject/pymatgen/issues/2999).

Basic immutable element object with all relevant properties.

Only one instance of Element for each symbol is stored after creation,
ensuring that a particular element behaves like a singleton. For all
attributes, missing data (i.e., data for which is not available) is
represented by a None unless otherwise stated.


* **Parameters**

    **symbol** (*str*) – Element symbol, e.g., “H”, “Fe”



#### Z()
Atomic number.


* **Type**

    int



#### symbol()
Element symbol.


* **Type**

    str



#### long_name()
Long name for element. E.g., “Hydrogen”.


* **Type**

    str



#### atomic_radius_calculated()
Calculated atomic radius for the element. This is the empirical value.
Data is obtained from [http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page](http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)).


* **Type**

    float



#### van_der_waals_radius()
Van der Waals radius for the element. This is the empirical value determined
from critical reviews of X-ray diffraction, gas kinetic collision cross-section, and other experimental
data by Bondi and later workers. The uncertainty in these values is on the order of 0.1 Å.
Data are obtained from “Atomic Radii of the Elements” in CRC Handbook of Chemistry and Physics,
91st Ed.; Haynes, W.M., Ed.; CRC Press: Boca Raton, FL, 2010.


* **Type**

    float



#### mendeleev_no()
Mendeleev number from definition given by Pettifor, D. G. (1984). A chemical scale
for crystal-structure maps. Solid State Communications, 51 (1), 31-34.


* **Type**

    int



#### electrical_resistivity()
Electrical resistivity.


* **Type**

    float



#### velocity_of_sound()
Velocity of sound.


* **Type**

    float



#### reflectivity()
Reflectivity.


* **Type**

    float



#### refractive_index()
Refractive index.


* **Type**

    float



#### poissons_ratio()
Poisson’s ratio.


* **Type**

    float



#### molar_volume()
Molar volume.


* **Type**

    float



#### electronic_structure()
Electronic structure. E.g., The electronic structure for Fe is represented
as [Ar].3d6.4s2.


* **Type**

    str



#### atomic_orbitals()
Atomic Orbitals. Energy of the atomic orbitals as a dict. E.g., The orbitals
energies in eV are represented as {‘1s’: -1.0, ‘2s’: -0.1}. Data is obtained from
[https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations](https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations).
The LDA values for neutral atoms are used.


* **Type**

    dict



#### thermal_conductivity()
Thermal conductivity.


* **Type**

    float



#### boiling_point()
Boiling point.


* **Type**

    float



#### melting_point()
Melting point.


* **Type**

    float



#### critical_temperature()
Critical temperature.


* **Type**

    float



#### superconduction_temperature()
Superconduction temperature.


* **Type**

    float



#### liquid_range()
Liquid range.


* **Type**

    float



#### bulk_modulus()
Bulk modulus.


* **Type**

    float



#### youngs_modulus()
Young’s modulus.


* **Type**

    float



#### brinell_hardness()
Brinell hardness.


* **Type**

    float



#### rigidity_modulus()
Rigidity modulus.


* **Type**

    float



#### mineral_hardness()
Mineral hardness.


* **Type**

    float



#### vickers_hardness()
Vicker’s hardness.


* **Type**

    float



#### density_of_solid()
Density of solid phase.


* **Type**

    float



#### coefficient_of_linear_thermal_expansion()
Coefficient of linear thermal expansion.


* **Type**

    float



#### ground_level()
Ground level for element.


* **Type**

    float



#### ionization_energies()
List of ionization energies. First value is the first
ionization energy, second is the second ionization energy, etc. Note that this is zero-based indexing!
So Element.ionization_energies[0] refer to the 1st ionization energy. Values are from the NIST Atomic
Spectra Database. Missing values are None.


* **Type**

    list[Optional[float]]



#### _property_ X(_: floa_ )
Pauling electronegativity of element. Note that if an element does not
have an Pauling electronegativity, a NaN float is returned.


#### as_dict()
Makes Element obey the general json interface used in pymatgen for
easier serialization.


#### _property_ atomic_mass(_: FloatWithUni_ )
Returns:
float: The atomic mass of the element in amu.


#### _property_ atomic_radius(_: FloatWithUnit | Non_ )
Returns:
float | None: The atomic radius of the element in Ångstroms. Can be None for
some elements like noble gases.


#### _property_ average_anionic_radius(_: floa_ )
Average anionic radius for element (with units). The average is
taken over all negative oxidation states of the element for which
data is present.


#### _property_ average_cationic_radius(_: FloatWithUni_ )
Average cationic radius for element (with units). The average is
taken over all positive oxidation states of the element for which
data is present.


#### _property_ average_ionic_radius(_: FloatWithUni_ )
Average ionic radius for element (with units). The average is taken
over all oxidation states of the element for which data is present.


#### _property_ block(_: st_ )
Return the block character “s,p,d,f”.


#### _property_ common_oxidation_states(_: tuple[int, ..._ )
Tuple of common oxidation states.


#### _property_ data(_: dict[str, Any_ )
Returns dict of data for element.


#### _property_ electron_affinity(_: floa_ )
The amount of energy released when an electron is attached to a neutral atom.


#### _property_ electronic_structure(_: st_ )
Electronic structure as string, with only valence electrons.
E.g., The electronic structure for Fe is represented as ‘[Ar].3d6.4s2’.


#### _static_ from_Z(Z: int)
Get an element from an atomic number.


* **Parameters**

    **Z** (*int*) – Atomic number



* **Returns**

    Element with atomic number Z.



#### _static_ from_dict(d)
Makes Element obey the general json interface used in pymatgen for
easier serialization.


#### _static_ from_name(name: str)
Get an element from its long name.


* **Parameters**

    **name** – Long name of the element, e.g. ‘Hydrogen’ or ‘Iron’. Not case-sensitive.



* **Returns**

    Element with the name ‘name’



#### _static_ from_row_and_group(row: int, group: int)
Returns an element from a row and group number.
Important Note: For lanthanoids and actinoids, the row number must
be 8 and 9, respectively, and the group number must be
between 3 (La, Ac) and 17 (Lu, Lr). This is different than the
value for Element(symbol).row and Element(symbol).group for these
elements.


* **Parameters**


    * **row** (*int*) – (pseudo) row number. This is the
    standard row number except for the lanthanoids
    and actinoids for which it is 8 or 9, respectively.


    * **group** (*int*) – (pseudo) group number. This is the
    standard group number except for the lanthanoids
    and actinoids for which it is 3 (La, Ac) to 17 (Lu, Lr).


**NOTE**: The 18 group number system is used, i.e. noble gases are group 18.


#### _property_ full_electronic_structure(_: list[tuple[int, str, int]_ )
Full electronic structure as tuple.
E.g., The electronic structure for Fe is represented as:
[(1, “s”, 2), (2, “s”, 2), (2, “p”, 6), (3, “s”, 2), (3, “p”, 6),
(3, “d”, 6), (4, “s”, 2)].


#### _property_ ground_state_term_symbol()
Ground state term symbol
Selected based on Hund’s Rule.


#### _property_ group(_: in_ )
Returns the periodic table group of the element.
Note: For lanthanoids and actinoids, the group is always 3.


#### _property_ icsd_oxidation_states(_: tuple[int, ..._ )
Tuple of all oxidation states with at least 10 instances in
ICSD database AND at least 1% of entries for that element.


#### _property_ ionic_radii(_: dict[int, float_ )
All ionic radii of the element as a dict of
{oxidation state: ionic radii}. Radii are given in angstrom.


#### _property_ ionization_energy(_: float | Non_ )
First ionization energy of element.


#### _property_ is_actinoid(_: boo_ )
True if element is a actinoid.


#### _property_ is_alkali(_: boo_ )
True if element is an alkali metal.


#### _property_ is_alkaline(_: boo_ )
True if element is an alkaline earth metal (group II).


#### _property_ is_chalcogen(_: boo_ )
True if element is a chalcogen.


#### _property_ is_halogen(_: boo_ )
True if element is a halogen.


#### _property_ is_lanthanoid(_: boo_ )
True if element is a lanthanoid.


#### _property_ is_metal(_: boo_ )
True if is a metal.


#### _property_ is_metalloid(_: boo_ )
True if element is a metalloid.


#### _property_ is_noble_gas(_: boo_ )
True if element is noble gas.


#### _property_ is_post_transition_metal(_: boo_ )
True if element is a post-transition or poor metal.


#### _property_ is_quadrupolar(_: boo_ )
Checks if this element can be quadrupolar.


#### _property_ is_rare_earth_metal(_: boo_ )
True if element is a rare earth metal.


#### _property_ is_transition_metal(_: boo_ )
True if element is a transition metal.


#### _static_ is_valid_symbol(symbol: str)
Returns true if symbol is a valid element symbol.


* **Parameters**

    **symbol** (*str*) – Element symbol



* **Returns**

    True if symbol is a valid element (e.g., “H”).



* **Return type**

    bool



#### _property_ iupac_ordering()
Ordering according to Table VI of “Nomenclature of Inorganic Chemistry
(IUPAC Recommendations 2005)”. This ordering effectively follows the
groups and rows of the periodic table, except the Lanthanides, Actinides
and hydrogen.


#### _property_ max_oxidation_state(_: floa_ )
Maximum oxidation state for element.


#### _property_ min_oxidation_state(_: floa_ )
Minimum oxidation state for element.


#### _property_ nmr_quadrupole_moment(_: dict[str, FloatWithUnit_ )
Get a dictionary the nuclear electric quadrupole moment in units of
e\*millibarns for various isotopes.


#### _property_ number(_: in_ )
Alternative attribute for atomic number Z.


#### _property_ oxidation_states(_: tuple[int, ..._ )
Tuple of all known oxidation states.


#### _static_ print_periodic_table(filter_function: Callable | None = None)
A pretty ASCII printer for the periodic table, based on some
filter_function.


* **Parameters**

    **filter_function** – A filtering function taking an Element as input
    and returning a boolean. For example, setting
    filter_function = lambda el: el.X > 2 will print a periodic
    table containing only elements with Pauling electronegativity > 2.



#### _property_ row(_: in_ )
Returns the periodic table row of the element.
Note: For lanthanoids and actinoids, the row is always 6 or 7,
respectively.


#### _property_ term_symbols(_: list[list[str]_ )
All possible Russell-Saunders term symbol of the Element.
eg. L = 1, n_e = 2 (s2) returns [[‘1D2’], [‘3P0’, ‘3P1’, ‘3P2’], [‘1S0’]].


#### _property_ valence()
From full electron config obtain valence subshell angular moment (L) and number of valence e- (v_e).


### _class_ Specie(symbol: SpeciesLike, oxidation_state: float | None = None, spin: float | None = None)
Bases: `Species`

This maps the historical grammatically inaccurate Specie to Species
to maintain backwards compatibility.


* **Parameters**


    * **symbol** (*str*) – Element symbol optionally incl. oxidation state. E.g. Fe, Fe2+, O2-.


    * **oxidation_state** (*float*) – Explicit oxidation state of element, e.g. -2, -1, 0, 1, 2, …
    If oxidation state is present in symbol, this argument is ignored.


    * **spin** – Spin associated with Species. Defaults to None.



* **Raises**

    **ValueError** – If oxidation state passed both in symbol string and via
        oxidation_state kwarg.



### _class_ Species(symbol: SpeciesLike, oxidation_state: float | None = None, spin: float | None = None)
Bases: `MSONable`, [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)

An extension of Element with optional oxidation state and spin. Properties associated
with Species should be “idealized” values, not calculated values. For example,
high-spin Fe2+ may be assigned an idealized spin of +5, but an actual Fe2+ site may be
calculated to have a magmom of +4.5. Calculated properties should be assigned to Site
objects, and not Species.


* **Parameters**


    * **symbol** (*str*) – Element symbol optionally incl. oxidation state. E.g. Fe, Fe2+, O2-.


    * **oxidation_state** (*float*) – Explicit oxidation state of element, e.g. -2, -1, 0, 1, 2, …
    If oxidation state is present in symbol, this argument is ignored.


    * **spin** – Spin associated with Species. Defaults to None.



* **Raises**

    **ValueError** – If oxidation state passed both in symbol string and via
        oxidation_state kwarg.



#### STRING_MODE(_ = 'SUPERSCRIPT_ )

#### as_dict()
Json-able dictionary representation.


#### _property_ element(_: Elemen_ )
Underlying element object.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation.



* **Returns**

    Species.



#### _static_ from_str(species_string: str)
Returns a Species from a string representation.


* **Parameters**

    **species_string** (*str*) – A typical string representation of a
    species, e.g., “Mn2+”, “Fe3+”, “O2-“.



* **Returns**

    A Species object.



* **Raises**

    **ValueError if species_string cannot be interpreted.** –



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead

Use from_str instead.


#### get_crystal_field_spin(coordination: Literal['oct', 'tet'] = 'oct', spin_config: Literal['low', 'high'] = 'high')
Calculate the crystal field spin based on coordination and spin
configuration. Only works for transition metal species.


* **Parameters**


    * **coordination** (*"oct"** | **"tet"*) – Tetrahedron or octahedron crystal site coordination


    * **spin_config** (*"low"** | **"high"*) – Whether the species is in a high or low spin state



* **Returns**

    Crystal field spin in Bohr magneton.



* **Raises**


    * **AttributeError if species is not a valid transition metal**** or ****has** – an invalid oxidation state.


    * **ValueError if invalid coordination**** or ****spin_config.** –



#### get_nmr_quadrupole_moment(isotope: str | None = None)
Gets the nuclear electric quadrupole moment in units of e \* millibarns.


* **Parameters**

    **isotope** (*str*) – the isotope to get the quadrupole moment for
    default is None, which gets the lowest mass isotope



#### get_shannon_radius(cn: str, spin: Literal['', 'Low Spin', 'High Spin'] = '', radius_type: Literal['ionic', 'crystal'] = 'ionic')
Get the local environment specific ionic radius for species.


* **Parameters**


    * **cn** (*str*) – Coordination using roman letters. Supported values are
    I-IX, as well as IIIPY, IVPY and IVSQ.


    * **spin** (*str*) – Some species have different radii for different
    spins. You can get specific values using “High Spin” or
    “Low Spin”. Leave it as “” if not available. If only one spin
    data is available, it is returned and this spin parameter is
    ignored.


    * **radius_type** (*str*) – Either “crystal” or “ionic” (default).



* **Returns**

    Shannon radius for specie in the specified environment.



#### _property_ ionic_radius(_: float | Non_ )
Ionic radius of specie. Returns None if data is not present.


#### _property_ oxi_state(_: float | Non_ )
Oxidation state of Species.


#### _property_ spin(_: float | Non_ )
Spin of Species.


#### to_pretty_string()
String without properties.


### get_el_sp(obj: int | SpeciesLike)
Utility method to get an Element, Species or DummySpecies from any input.

If obj is in itself an element or a specie, it is returned automatically.
If obj is an int or a string representing an integer, the Element with the
atomic number obj is returned.
If obj is a string, Species parsing will be attempted (e.g. Mn2+). Failing that
Element parsing will be attempted (e.g. Mn). Failing that DummyElement parsing
will be attempted.


* **Parameters**

    **obj** (*Element/Species/str/int*) – An arbitrary object. Supported objects
    are actual Element/Species objects, integers (representing atomic
    numbers) or strings (element symbols or species strings).



* **Raises**

    **ValueError** – if obj cannot be converted into an Element or Species.



* **Returns**

    with a bias for the maximum number of properties

        that can be determined.




* **Return type**

    Species | Element


## pymatgen.core.sites module

This module defines classes representing non-periodic and periodic sites.


### _class_ PeriodicSite(species: SpeciesLike | CompositionLike, coords: ArrayLike, lattice: Lattice, to_unit_cell: bool = False, coords_are_cartesian: bool = False, properties: dict | None = None, label: str | None = None, skip_checks: bool = False)
Bases: `Site`, `MSONable`

Extension of generic Site object to periodic systems.
PeriodicSite includes a lattice system.

Create a periodic site.


* **Parameters**


    * **species** – Species on the site. Can be:
    i.  A Composition-type object (preferred)
    ii. An  element / species specified either as a string

    > symbols, e.g. “Li”, “Fe2+”, “P” or atomic numbers,
    > e.g., 3, 56, or actual Element or Species objects.

    iii.Dict of elements/species and occupancies, e.g.,

        {“Fe” : 0.5, “Mn”:0.5}. This allows the setup of
        disordered structures.



    * **coords** – Coordinates of site, fractional coordinates
    by default. See `coords_are_cartesian` for more details.


    * **lattice** – Lattice associated with the site.


    * **to_unit_cell** – Translates fractional coordinate to the
    basic unit cell, i.e. all fractional coordinates satisfy 0
    <= a < 1. Defaults to False.


    * **coords_are_cartesian** – Set to True if you are providing
    Cartesian coordinates. Defaults to False.


    * **properties** – Properties associated with the site as a dict, e.g.
    {“magmom”: 5}. Defaults to None.


    * **label** – Label for the site. Defaults to None.


    * **skip_checks** – Whether to ignore all the usual checks and just
    create the site. Use this if the PeriodicSite is created in a
    controlled manner and speed is desired.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ a(_: floa_ )
Fractional a coordinate.


#### as_dict(verbosity: int = 0)
JSON-serializable dict representation of PeriodicSite.


* **Parameters**

    **verbosity** (*int*) – Verbosity level. Default of 0 only includes the matrix
    representation. Set to 1 for more details such as Cartesian coordinates, etc.



#### _property_ b(_: floa_ )
Fractional b coordinate.


#### _property_ c(_: floa_ )
Fractional c coordinate.


#### _property_ coords(_: ndarra_ )
Cartesian coordinates.


#### distance(other: PeriodicSite, jimage: ArrayLike | None = None)
Get distance between two sites assuming periodic boundary conditions.


* **Parameters**


    * **other** (*PeriodicSite*) – Other site to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of lattice
    translations, e.g., [1,0,0] implies to take periodic image
    that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    Distance between the two sites



* **Return type**

    distance (float)



#### distance_and_image(other: PeriodicSite, jimage: ArrayLike | None = None)
Gets distance and instance between two sites assuming periodic boundary
conditions. If the index jimage of two sites atom j is not specified it
selects the j image nearest to the i atom and returns the distance and
jimage indices in terms of lattice vector translations. If the index
jimage of atom j is specified it returns the distance between the ith
atom and the specified jimage atom, the given jimage is also returned.


* **Parameters**


    * **other** (*PeriodicSite*) – Other site to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of lattice
    translations, e.g., [1,0,0] implies to take periodic image
    that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    distance and periodic lattice translations
    of the other site for which the distance applies.



* **Return type**

    (distance, jimage)



#### distance_and_image_from_frac_coords(fcoords: ArrayLike, jimage: ArrayLike | None = None)
Gets distance between site and a fractional coordinate assuming
periodic boundary conditions. If the index jimage of two sites atom j
is not specified it selects the j image nearest to the i atom and
returns the distance and jimage indices in terms of lattice vector
translations. If the index jimage of atom j is specified it returns the
distance between the i atom and the specified jimage atom, the given
jimage is also returned.


* **Parameters**


    * **fcoords** (*3x1 array*) – fcoords to get distance from.


    * **jimage** (*3x1 array*) – Specific periodic image in terms of
    lattice translations, e.g., [1,0,0] implies to take periodic
    image that is one a-lattice vector away. If jimage is None,
    the image that is nearest to the site is found.



* **Returns**

    distance and periodic lattice translations
    of the other site for which the distance applies.



* **Return type**

    (distance, jimage)



#### _property_ frac_coords(_: ndarra_ )
Fractional coordinates.


#### _classmethod_ from_dict(dct, lattice=None)
Create PeriodicSite from dict representation.


* **Parameters**


    * **dct** (*dict*) – dict representation of PeriodicSite


    * **lattice** – Optional lattice to override lattice specified in d.
    Useful for ensuring all sites in a structure share the same
    lattice.



* **Returns**

    PeriodicSite



#### is_periodic_image(other: PeriodicSite, tolerance: float = 1e-08, check_lattice: bool = True)
Returns True if sites are periodic images of each other.


* **Parameters**


    * **other** (*PeriodicSite*) – Other site


    * **tolerance** (*float*) – Tolerance to compare fractional coordinates


    * **check_lattice** (*bool*) – Whether to check if the two sites have the
    same lattice.



* **Returns**

    True if sites are periodic images of each other.



* **Return type**

    bool



#### _property_ lattice(_: Lattic_ )
Lattice associated with PeriodicSite.


#### to_unit_cell(in_place=False)
Move frac coords to within the unit cell.


#### _property_ x(_: floa_ )
Cartesian x coordinate.


#### _property_ y(_: floa_ )
Cartesian y coordinate.


#### _property_ z(_: floa_ )
Cartesian z coordinate.


### _class_ Site(species: SpeciesLike | CompositionLike, coords: ArrayLike, properties: dict | None = None, label: str | None = None, skip_checks: bool = False)
Bases: `Hashable`, `MSONable`

A generalized *non-periodic* site. This is essentially a composition
at a point in space, with some optional properties associated with it. A
Composition is used to represent the atoms and occupancy, which allows for
disordered site representation. Coords are given in standard Cartesian
coordinates.

Creates a non-periodic Site.


* **Parameters**


    * **species** – Species on the site. Can be:
    i.  A Composition-type object (preferred)
    ii. An  element / species specified either as a string

    > symbols, e.g. “Li”, “Fe2+”, “P” or atomic numbers,
    > e.g., 3, 56, or actual Element or Species objects.

    iii.Dict of elements/species and occupancies, e.g.,

        {“Fe” : 0.5, “Mn”:0.5}. This allows the setup of
        disordered structures.



    * **coords** – Cartesian coordinates of site.


    * **properties** – Properties associated with the site as a dict, e.g.
    {“magmom”: 5}. Defaults to None.


    * **label** – Label for the site. Defaults to None.


    * **skip_checks** – Whether to ignore all the usual checks and just
    create the site. Use this if the Site is created in a controlled
    manner and speed is desired.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
JSON-serializable dict representation for Site.


#### distance(other)
Get distance between two sites.


* **Parameters**

    **other** – Other site.



* **Returns**

    distance



* **Return type**

    float



#### distance_from_point(pt)
Returns distance between the site and a point in space.


* **Parameters**

    **pt** – Cartesian coordinates of point.



* **Returns**

    distance



* **Return type**

    float



#### _classmethod_ from_dict(dct: dict)
Create Site from dict representation.


#### _property_ is_ordered(_: boo_ )
True if site is an ordered site, i.e., with a single species with
occupancy 1.


#### _property_ label(_: st_ )
Site label.


#### position_atol(_ = 1e-0_ )

#### _property_ specie(_: Element | Species | DummySpecie_ )
The Species/Element at the site. Only works for ordered sites. Otherwise
an AttributeError is raised. Use this property sparingly.  Robust
design should make use of the property species instead. Note that the
singular of species is also species. So the choice of this variable
name is governed by programmatic concerns as opposed to grammar.


* **Raises**

    **AttributeError if Site is not ordered.** –



#### _property_ species(_: Compositio_ )
The species on the site as a composition, e.g., Fe0.5Mn0.5.


#### _property_ species_string(_: st_ )
String representation of species on the site.


#### _property_ x(_: floa_ )
Cartesian x coordinate.


#### _property_ y(_: floa_ )
Cartesian y coordinate.


#### _property_ z(_: floa_ )
Cartesian z coordinate.

## pymatgen.core.spectrum module

This module defines classes to represent any type of spectrum, essentially any
x y value pairs.


### _class_ Spectrum(x: ArrayLike, y: ArrayLike, \*args, \*\*kwargs)
Bases: `MSONable`

Base class for any type of xas, essentially just x, y values. Examples
include XRD patterns, XANES, EXAFS, NMR, DOS, etc.

Implements basic tools like application of smearing, normalization, addition
multiplication, etc.

Subclasses should extend this object and ensure that super is called with
ALL args and kwargs. That ensures subsequent things like add and mult work
properly.


* **Parameters**


    * **x** (*ndarray*) – A ndarray of N values.


    * **y** (*ndarray*) – A ndarray of N x k values. The first dimension must be
    the same as that of x. Each of the k values are interpreted as separate.


    * **\*args** – All subclasses should provide args other than x and y
    when calling super, e.g., super().__init__(
    x, y, arg1, arg2, kwarg1=val1, ..). This guarantees the +, -,

    ```
    *
    ```

    ,
    etc. operators work properly.


    * **\*\*kwargs** – Same as that for

    ```
    *
    ```

    args.



#### XLABEL(_ = 'x_ )

#### YLABEL(_ = 'y_ )

#### copy()

* **Returns**

    Copy of Spectrum object.



#### get_interpolated_value(x: float)
Returns an interpolated y value for a particular x value.


* **Parameters**

    **x** – x value to return the y value for



* **Returns**

    Value of y at x



#### normalize(mode: Literal['max', 'sum'] = 'max', value: float = 1.0)
Normalize the spectrum with respect to the sum of intensity.


* **Parameters**


    * **mode** (*"max"** | **"sum"*) – Normalization mode. “max” sets the max y value to value,
    e.g., in XRD patterns. “sum” sets the sum of y to a value, i.e., like a
    probability density.


    * **value** (*float*) – Value to normalize to. Defaults to 1.



#### smear(sigma: float = 0.0, func: str | Callable = 'gaussian')
Apply Gaussian/Lorentzian smearing to spectrum y value.


* **Parameters**


    * **sigma** – Std dev for Gaussian smear function


    * **func** – “gaussian” or “lorentzian” or a callable. If this is a callable, the sigma value is ignored. The
    callable should only take a single argument (a numpy array) and return a set of weights.



### lorentzian(x, x_0: float = 0, sigma: float = 1.0)

* **Parameters**


    * **x** – x values


    * **x_0** – Center


    * **sigma** – FWHM



* **Returns**

    Value of lorentzian at x.


## pymatgen.core.structure module

This module provides classes used to define a non-periodic molecule and a
periodic structure.


### _class_ IMolecule(species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float = 0.0, spin_multiplicity: int | None = None, validate_proximity: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, charge_spin_check: bool = True, properties: dict | None = None)
Bases: `SiteCollection`, `MSONable`

Basic immutable Molecule object without periodicity. Essentially a
sequence of sites. IMolecule is made to be immutable so that they can
function as keys in a dict. For a mutable object, use the Molecule class.

Molecule extends Sequence and Hashable, which means that in many cases,
it can be used like any Python sequence. Iterating through a molecule is
equivalent to going through the sites in sequence.

Create a Molecule.


* **Parameters**


    * **species** – list of atomic species. Possible kinds of input include a
    list of dict of elements/species and occupancies, a List of
    elements/specie specified as actual Element/Species, Strings
    (“Fe”, “Fe2+”) or atomic numbers (1,56).


    * **coords** (*3x1 array*) – list of Cartesian coordinates of each species.


    * **charge** (*float*) – Charge for the molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multiplicity for molecule.
    Defaults to None, which means that the spin multiplicity is
    set to 1 if the molecule has no unpaired electrons and to 2
    if there are unpaired electrons.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 1 Ang apart. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as
    a dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The
    sequences have to be the same length as the atomic species
    and fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **charge_spin_check** (*bool*) – Whether to check that the charge and
    spin multiplicity are compatible with each other. Defaults
    to True.


    * **properties** (*dict*) – dictionary containing properties associated
    with the whole molecule.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _find_nn_pos_before_site(site_idx)
Returns index of nearest neighbor atoms.


#### _properties(_: dic_ )

#### as_dict()
JSON-serializable dict representation of Molecule.


#### break_bond(ind1: int, ind2: int, tol: float = 0.2)
Returns two molecules based on breaking the bond between atoms at index
ind1 and ind2.


* **Parameters**


    * **ind1** (*int*) – 1st site index


    * **ind2** (*int*) – 2nd site index


    * **tol** (*float*) – Relative tolerance to test. Basically, the code
    checks if the distance between the sites is less than (1 +
    tol) \* typical bond distances. Defaults to 0.2, i.e.,
    20% longer.



* **Returns**

    Two Molecule objects representing the two clusters formed from
    breaking the bond.



#### _property_ center_of_mass(_: ndarra_ )
Center of mass of molecule.


#### _property_ charge(_: floa_ )
Charge of molecule.


#### copy()
Convenience method to get a copy of the molecule.


* **Returns**

    IMolecule | Molecule



#### _classmethod_ from_dict(dct)
Reconstitute a Molecule object from a dict representation created using as_dict().


* **Parameters**

    **dct** (*dict*) – dict representation of Molecule.



* **Returns**

    Molecule



#### _classmethod_ from_file(filename)
Reads a molecule from a file. Supported formats include xyz,
gaussian input (gjf|g03|g09|com|inp), Gaussian output (.out|and
pymatgen’s JSON-serialized molecules. Using openbabel,
many more extensions are supported but requires openbabel to be
installed.


* **Parameters**

    **filename** (*str*) – The filename to read from.



* **Returns**

    Molecule



#### _classmethod_ from_sites(sites: Sequence[Site], charge: float = 0, spin_multiplicity: int | None = None, validate_proximity: bool = False, charge_spin_check: bool = True, properties: dict | None = None)
Convenience constructor to make a Molecule from a list of sites.


* **Parameters**


    * **sites** (*[**Site**]*) – Sequence of Sites.


    * **charge** (*int*) – Charge of molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multicipity. Defaults to None,
    in which it is determined automatically.


    * **validate_proximity** (*bool*) – Whether to check that atoms are too
    close.


    * **charge_spin_check** (*bool*) – Whether to check that the charge and
    spin multiplicity are compatible with each other. Defaults
    to True.


    * **properties** (*dict*) – dictionary containing properties associated
    with the whole molecule.



* **Raises**

    **ValueError** – If sites is empty



* **Returns**

    Molecule



#### _classmethod_ from_str(input_string: str, fmt: Literal['xyz', 'gjf', 'g03', 'g09', 'com', 'inp', 'json', 'yaml'])
Reads the molecule from a string.


* **Parameters**


    * **input_string** (*str*) – String to parse.


    * **fmt** (*str*) – Format to output to. Defaults to JSON unless filename
    is provided. If fmt is specifies, it overrides whatever the
    filename is. Options include “xyz”, “gjf”, “g03”, “json”. If
    you have OpenBabel installed, any of the formats supported by
    OpenBabel. Non-case sensitive.



* **Returns**

    IMolecule or Molecule.



#### get_boxed_structure(a: float, b: float, c: float, images: ArrayLike = (1, 1, 1), random_rotation: bool = False, min_dist: float = 1.0, cls=None, offset: ArrayLike | None = None, no_cross: bool = False, reorder: bool = True)
Creates a Structure from a Molecule by putting the Molecule in the
center of a orthorhombic box. Useful for creating Structure for
calculating molecules using periodic codes.


* **Parameters**


    * **a** (*float*) – a-lattice parameter.


    * **b** (*float*) – b-lattice parameter.


    * **c** (*float*) – c-lattice parameter.


    * **images** – No. of boxed images in each direction. Defaults to
    (1, 1, 1), meaning single molecule with 1 lattice parameter
    in each direction.


    * **random_rotation** (*bool*) – Whether to apply a random rotation to
    each molecule. This jumbles all the molecules so that they
    are not exact images of each other.


    * **min_dist** (*float*) – The minimum distance that atoms should be from
    each other. This is only used if random_rotation is True.
    The randomized rotations are searched such that no two atoms
    are less than min_dist from each other.


    * **cls** – The Structure class to instantiate (defaults to pymatgen
    structure)


    * **offset** – Translation to offset molecule from center of mass coords


    * **no_cross** – Whether to forbid molecule coords from extending beyond
    boundary of box.


    * **reorder** – Whether to reorder the sites to be in electronegativity
    order.



* **Returns**

    Structure containing molecule in a box.



#### get_centered_molecule()
Returns a Molecule centered at the center of mass.


* **Returns**

    Molecule centered with center of mass at origin.



#### get_covalent_bonds(tol: float = 0.2)
Determines the covalent bonds in a molecule.


* **Parameters**

    **tol** (*float*) – The tol to determine bonds in a structure. See
    CovalentBond.is_bonded.



* **Returns**

    List of bonds



#### get_distance(i: int, j: int)
Get distance between site i and j.


* **Parameters**


    * **i** (*int*) – 1st site index


    * **j** (*int*) – 2nd site index



* **Returns**

    Distance between the two sites.



#### get_neighbors(site: Site, r: float)
Get all neighbors to a site within a sphere of radius r. Excludes the
site itself.


* **Parameters**


    * **site** (*Site*) – Site at the center of the sphere.


    * **r** (*float*) – Radius of sphere.



* **Returns**

    Neighbor



#### get_neighbors_in_shell(origin: ArrayLike, r: float, dr: float)
Returns all sites in a shell centered on origin (coords) between radii
r-dr and r+dr.


* **Parameters**


    * **origin** (*3x1 array*) – Cartesian coordinates of center of sphere.


    * **r** (*float*) – Inner radius of shell.


    * **dr** (*float*) – Width of shell.



* **Returns**

    Neighbor



#### get_sites_in_sphere(pt: ArrayLike, r: float)
Find all sites within a sphere from a point.


* **Parameters**


    * **pt** (*3x1 array*) – Cartesian coordinates of center of sphere


    * **r** (*float*) – Radius of sphere.



* **Returns**

    Neighbor



#### get_zmatrix()
Returns a z-matrix representation of the molecule.


#### _property_ nelectrons(_: floa_ )
Number of electrons in the molecule.


#### _property_ spin_multiplicity(_: floa_ )
Spin multiplicity of molecule.


#### to(filename: str = '', fmt: str = '')
Outputs the molecule to a file or string.


* **Parameters**


    * **filename** (*str*) – If provided, output will be written to a file. If
    fmt is not specified, the format is determined from the
    filename. Defaults is None, i.e. string output.


    * **fmt** (*str*) – Format to output to. Defaults to JSON unless filename
    is provided. If fmt is specifies, it overrides whatever the
    filename is. Options include “xyz”, “gjf”, “g03”, “json”. If
    you have OpenBabel installed, any of the formats supported by
    OpenBabel. Non-case sensitive.



* **Returns**

    String representation of molecule in given format. If a filename

        is provided, the same string is written to the file.




* **Return type**

    str



### _class_ IStructure(lattice: ArrayLike | Lattice, species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False, coords_are_cartesian: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, properties: dict | None = None)
Bases: `SiteCollection`, `MSONable`

Basic immutable Structure object with periodicity. Essentially a sequence
of PeriodicSites having a common lattice. IStructure is made to be
(somewhat) immutable so that they can function as keys in a dict. To make
modifications, use the standard Structure object instead. Structure
extends Sequence and Hashable, which means that in many cases,
it can be used like any Python sequence. Iterating through a
structure is equivalent to going through the sites in sequence.

Create a periodic structure.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** (*[**Species**]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/Cartesian coordinates of
    each species.


    * **charge** (*int*) – overall charge of the structure. Defaults to behavior
    in SiteCollection where total charge is the sum of the oxidation
    states.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to map all sites into the unit cell,
    i.e. fractional coords between 0 and 1. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g. {“magmom”:[5, 5, 5, 5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **properties** (*dict*) – Properties associated with the whole structure.
    Will be serialized when writing the structure to JSON or YAML but is
    lost when converting to other formats.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _get_neighbor_list_py(r: float, sites: list[PeriodicSite] | None = None, numerical_tol: float = 1e-08, exclude_self: bool = True)
A python version of getting neighbor_list. The returned values are a tuple of
numpy arrays (center_indices, points_indices, offset_vectors, distances).
Atom center_indices[i] has neighbor atom points_indices[i] that is
translated by offset_vectors[i] lattice vectors, and the distance is
distances[i].


* **Parameters**


    * **r** (*float*) – Radius of sphere


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.


    * **exclude_self** (*bool*) – whether to exclude atom neighboring with itself within
    numerical tolerance distance, default to True



* **Returns**

    (center_indices, points_indices, offset_vectors, distances)



* **Return type**

    tuple



#### _properties(_: dic_ )

#### as_dataframe()
Create a Pandas dataframe of the sites. Structure-level attributes are stored in DataFrame.attrs.

### Example

Species    a    b             c    x             y             z  magmom
0    (Si)  0.0  0.0  0.000000e+00  0.0  0.000000e+00  0.000000e+00       5
1    (Si)  0.0  0.0  1.000000e-7  0.0 -2.217138e-7  3.135509e-7      -5


#### as_dict(verbosity=1, fmt=None, \*\*kwargs)
Dict representation of Structure.


* **Parameters**


    * **verbosity** (*int*) – Verbosity level. Default of 1 includes both
    direct and Cartesian coordinates for all sites, lattice
    parameters, etc. Useful for reading and for insertion into a
    database. Set to 0 for an extremely lightweight version
    that only includes sufficient information to reconstruct the
    object.


    * **fmt** (*str*) – Specifies a format for the dict. Defaults to None,
    which is the default format used in pymatgen. Other options
    include “abivars”.


    * **\*\*kwargs** – Allow passing of other kwargs needed for certain


    * **formats** –


    * **e.g.** –


    * **"abivars".** –



* **Returns**

    JSON-serializable dict representation.



#### _property_ charge(_: floa_ )
Overall charge of the structure.


#### copy(site_properties=None, sanitize=False, properties=None)
Convenience method to get a copy of the structure, with options to add
site properties.


* **Parameters**


    * **site_properties** (*dict*) – Properties to add or override. The
    properties are specified in the same way as the constructor,
    i.e., as a dict of the form {property: [values]}. The
    properties should be in the order of the *original* structure
    if you are performing sanitization.


    * **sanitize** (*bool*) – If True, this method will return a sanitized
    structure. Sanitization performs a few things: (i) The sites are
    sorted by electronegativity, (ii) a LLL lattice reduction is
    carried out to obtain a relatively orthogonalized cell,
    (iii) all fractional coords for sites are mapped into the
    unit cell.


    * **properties** (*dict*) – General properties to add or override.



* **Returns**

    A copy of the Structure, with optionally new site_properties and
    optionally sanitized.



#### _property_ density(_: floa_ )
Returns the density in units of g/cm^3.


#### _property_ distance_matrix(_: ndarra_ )
Returns the distance matrix between all sites in the structure. For
periodic structures, this should return the nearest image distance.


#### _property_ frac_coords()
Fractional coordinates as a Nx3 numpy array.


#### _classmethod_ from_dict(d: dict[str, Any], fmt: Literal['abivars'] | None = None)
Reconstitute a Structure object from a dict representation of Structure
created using as_dict().


* **Parameters**


    * **d** (*dict*) – Dict representation of structure.


    * **fmt** (*'abivars'** | **None*) – Use structure_from_abivars() to parse the dict. Defaults to None.



* **Returns**

    Structure object



#### _classmethod_ from_file(filename, primitive=False, sort=False, merge_tol=0.0, \*\*kwargs)
Reads a structure from a file. For example, anything ending in
a “cif” is assumed to be a Crystallographic Information Format file.
Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT,
vasprun.xml, CSSR, Netcdf and pymatgen’s JSON-serialized structures.


* **Parameters**


    * **filename** (*str*) – The filename to read from.


    * **primitive** (*bool*) – Whether to convert to a primitive cell. Only available for CIFs. Defaults to False.


    * **sort** (*bool*) – Whether to sort sites. Default to False.


    * **merge_tol** (*float*) – If this is some positive number, sites that are within merge_tol from each other will be
    merged. Usually 0.01 should be enough to deal with common numerical issues.


    * **kwargs** – Passthrough to relevant reader. E.g. if the file has CIF format, the kwargs will be passed
    through to CifParser.



* **Returns**

    Structure.



#### _classmethod_ from_magnetic_spacegroup(msg: str | [MagneticSpaceGroup](pymatgen.symmetry.md#pymatgen.symmetry.maggroups.MagneticSpaceGroup), lattice: list | np.ndarray | Lattice, species: Sequence[str | Element | Species | DummySpecies | Composition], coords: Sequence[Sequence[float]], site_properties: dict[str, Sequence], coords_are_cartesian: bool = False, tol: float = 1e-05, labels: Sequence[str | None] | None = None)
Generate a structure using a magnetic spacegroup. Note that only
symmetrically distinct species, coords and magmoms should be provided.]
All equivalent sites are generated from the spacegroup operations.


* **Parameters**


    * **msg** (*str/list/pymatgen.symmetry.maggroups.MagneticSpaceGroup*) – The magnetic spacegroup.
    If a string, it will be interpreted as one of the notations
    supported by MagneticSymmetryGroup, e.g., “R-3’c” or “Fm’-3’m”.
    If a list of two ints, it will be interpreted as the number of
    the spacegroup in its Belov, Neronova and Smirnova (BNS) setting.


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
    Note that no attempt is made to check that the lattice is
    compatible with the spacegroup specified. This may be
    introduced in a future version.


    * **species** (*[**Species**]*) – Sequence of species on each site. Can take in
    flexible input, including:
    i.  A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        1. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Unlike Structure.from_spacegroup(),
    this argument is mandatory, since magnetic moment information
    has to be included. Note that the *direction* of the supplied
    magnetic moment relative to the crystal is important, even if
    the resulting structure is used for collinear calculations.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **tol** (*float*) – A fractional tolerance to deal with numerical
    precision issues in determining if orbits are the same.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.



* **Returns**

    Structure | IStructure



#### _classmethod_ from_sites(sites: list[PeriodicSite], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False, properties: dict | None = None)
Convenience constructor to make a Structure from a list of sites.


* **Parameters**


    * **sites** – Sequence of PeriodicSites. Sites must have the same
    lattice.


    * **charge** – Charge of structure.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to translate sites into the unit
    cell.


    * **properties** (*dict*) – Properties associated with the whole structure.
    Will be serialized when writing the structure to JSON or YAML but is
    lost when converting to other formats.



* **Raises**

    **ValueError** – If sites is empty or sites do not have the same lattice.



* **Returns**

    (Structure) Note that missing properties are set as None.



#### _classmethod_ from_spacegroup(sg: str | int, lattice: list | np.ndarray | Lattice, species: Sequence[str | Element | Species | DummySpecies | Composition], coords: Sequence[Sequence[float]], site_properties: dict[str, Sequence] | None = None, coords_are_cartesian: bool = False, tol: float = 1e-05, labels: Sequence[str | None] | None = None)
Generate a structure using a spacegroup. Note that only symmetrically
distinct species and coords should be provided. All equivalent sites
are generated from the spacegroup operations.


* **Parameters**


    * **sg** (*str/int*) – The spacegroup. If a string, it will be interpreted
    as one of the notations supported by
    pymatgen.symmetry.groups.Spacegroup. E.g., “R-3c” or “Fm-3m”.
    If an int, it will be interpreted as an international number.


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
    Note that no attempt is made to check that the lattice is
    compatible with the spacegroup specified. This may be
    introduced in a future version.


    * **species** (*[**Species**]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **tol** (*float*) – A fractional tolerance to deal with numerical
    precision issues in determining if orbits are the same.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.



#### _classmethod_ from_str(input_string: str, fmt: Literal['cif', 'poscar', 'cssr', 'json', 'yaml', 'xsf', 'mcsqs', 'res'], primitive: bool = False, sort: bool = False, merge_tol: float = 0.0, \*\*kwargs)
Reads a structure from a string.


* **Parameters**


    * **input_string** (*str*) – String to parse.


    * **fmt** (*str*) – A file format specification. One of “cif”, “poscar”, “cssr”,
    “json”, “yaml”, “xsf”, “mcsqs”.


    * **primitive** (*bool*) – Whether to find a primitive cell. Defaults to
    False.


    * **sort** (*bool*) – Whether to sort the sites in accordance to the default
    ordering criteria, i.e., electronegativity.


    * **merge_tol** (*float*) – If this is some positive number, sites that
    are within merge_tol from each other will be merged. Usually
    0.01 should be enough to deal with common numerical issues.


    * **\*\*kwargs** – Passthrough to relevant parser.



* **Returns**

    IStructure | Structure



#### get_all_neighbors(r: float, include_index: bool = False, include_image: bool = False, sites: Sequence[PeriodicSite] | None = None, numerical_tol: float = 1e-08)
Get neighbors for each atom in the unit cell, out to a distance r
Returns a list of list of neighbors for each site in structure.
Use this method if you are planning on looping over all sites in the
crystal. If you only want neighbors for a particular site, use the
method get_neighbors as it may not have to build such a large supercell
However if you are looping over all sites in the crystal, this method
is more efficient since it only performs one pass over a large enough
supercell to contain all possible atoms out to a distance r.
The return type is a [(site, dist) …] since most of the time,
subsequent processing requires the distance.

A note about periodic images: Before computing the neighbors, this
operation translates all atoms to within the unit cell (having
fractional coordinates within [0,1)). This means that the “image” of a
site does not correspond to how much it has been translates from its
current position, but which image of the unit cell it resides.


* **Parameters**


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.



* **Returns**

    [[pymatgen.core.structure.PeriodicNeighbor], ..]



#### get_all_neighbors_old(\*\*kwargs)

#### get_all_neighbors_py(r: float, include_index: bool = False, include_image: bool = False, sites: Sequence[PeriodicSite] | None = None, numerical_tol: float = 1e-08)
Get neighbors for each atom in the unit cell, out to a distance r
Returns a list of list of neighbors for each site in structure.
Use this method if you are planning on looping over all sites in the
crystal. If you only want neighbors for a particular site, use the
method get_neighbors as it may not have to build such a large supercell
However if you are looping over all sites in the crystal, this method
is more efficient since it only performs one pass over a large enough
supercell to contain all possible atoms out to a distance r.
The return type is a [(site, dist) …] since most of the time,
subsequent processing requires the distance.

A note about periodic images: Before computing the neighbors, this
operation translates all atoms to within the unit cell (having
fractional coordinates within [0,1)). This means that the “image” of a
site does not correspond to how much it has been translates from its
current position, but which image of the unit cell it resides.


* **Parameters**


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.



* **Returns**

    list[list[PeriodicNeighbor]]



#### get_distance(i: int, j: int, jimage=None)
Get distance between site i and j assuming periodic boundary
conditions. If the index jimage of two sites atom j is not specified it
selects the jimage nearest to the i atom and returns the distance and
jimage indices in terms of lattice vector translations if the index
jimage of atom j is specified it returns the distance between the i
atom and the specified jimage atom.


* **Parameters**


    * **i** (*int*) – 1st site index


    * **j** (*int*) – 2nd site index


    * **jimage** – Number of lattice translations in each lattice direction.
    Default is None for nearest image.



* **Returns**

    distance



#### get_miller_index_from_site_indexes(site_ids, round_dp=4, verbose=True)
Get the Miller index of a plane from a set of sites indexes.

A minimum of 3 sites are required. If more than 3 sites are given
the best plane that minimises the distance to all points will be
calculated.


* **Parameters**


    * **site_ids** (*list** of **int*) – A list of site indexes to consider. A
    minimum of three site indexes are required. If more than three
    sites are provided, the best plane that minimises the distance
    to all sites will be calculated.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)



#### get_neighbor_list(r: float, sites: Sequence[PeriodicSite] | None = None, numerical_tol: float = 1e-08, exclude_self: bool = True)
Get neighbor lists using numpy array representations without constructing
Neighbor objects. If the cython extension is installed, this method will
be orders of magnitude faster than get_all_neighbors_old and 2-3x faster
than get_all_neighbors.
The returned values are a tuple of numpy arrays
(center_indices, points_indices, offset_vectors, distances).
Atom center_indices[i] has neighbor atom points_indices[i] that is
translated by offset_vectors[i] lattice vectors, and the distance is
distances[i].


* **Parameters**


    * **r** (*float*) – Radius of sphere


    * **sites** (*list** of **Sites** or **None*) – sites for getting all neighbors,
    default is None, which means neighbors will be obtained for all
    sites. This is useful in the situation where you are interested
    only in one subspecies type, and makes it a lot faster.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.


    * **exclude_self** (*bool*) – whether to exclude atom neighboring with itself within
    numerical tolerance distance, default to True



* **Returns**

    (center_indices, points_indices, offset_vectors, distances)



* **Return type**

    tuple



#### get_neighbors(site: PeriodicSite, r: float, include_index: bool = False, include_image: bool = False)
Get all neighbors to a site within a sphere of radius r. Excludes the
site itself.


* **Parameters**


    * **site** (*Site*) – Which is the center of the sphere.


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.



* **Returns**

    PeriodicNeighbor



#### get_neighbors_in_shell(origin: ArrayLike, r: float, dr: float, include_index: bool = False, include_image: bool = False)
Returns all sites in a shell centered on origin (coords) between radii
r-dr and r+dr.


* **Parameters**


    * **origin** (*3x1 array*) – Cartesian coordinates of center of sphere.


    * **r** (*float*) – Inner radius of shell.


    * **dr** (*float*) – Width of shell.


    * **include_index** (*bool*) – Deprecated. Now, the non-supercell site index
    is always included in the returned data.


    * **include_image** (*bool*) – Deprecated. Now the supercell image
    is always included in the returned data.



* **Returns**

    [NearestNeighbor] where Nearest Neighbor is a named tuple containing
    (site, distance, index, image).



#### get_neighbors_old(\*\*kwargs)

#### get_orderings(mode: Literal['enum', 'sqs'] = 'enum', \*\*kwargs)
Returns list of orderings for a disordered structure. If structure
does not contain disorder, the default structure is returned.


* **Parameters**


    * **mode** (*"enum"** | **"sqs"*) – Either “enum” or “sqs”. If enum,
    the enumlib will be used to return all distinct
    orderings. If sqs, mcsqs will be used to return
    an sqs structure.


    * **kwargs** – kwargs passed to either
    pymatgen.command_line..enumlib_caller.EnumlibAdaptor
    or pymatgen.command_line.mcsqs_caller.run_mcsqs.
    For run_mcsqs, a default cluster search of 2 cluster interactions
    with 1NN distance and 3 cluster interactions with 2NN distance
    is set.



* **Returns**

    List[Structure]



#### get_primitive_structure(tolerance: float = 0.25, use_site_props: bool = False, constrain_latt: list | dict | None = None)
This finds a smaller unit cell than the input. Sometimes it doesn”t
find the smallest possible one, so this method is recursively called
until it is unable to find a smaller cell.

NOTE: if the tolerance is greater than 1/2 the minimum inter-site
distance in the primitive cell, the algorithm will reject this lattice.


* **Parameters**


    * **tolerance** (*float*) – Tolerance for each coordinate of a
    particular site in Angstroms. For example, [0.1, 0, 0.1] in cartesian
    coordinates will be considered to be on the same coordinates
    as [0, 0, 0] for a tolerance of 0.25. Defaults to 0.25.


    * **use_site_props** (*bool*) – Whether to account for site properties in
    differentiating sites.


    * **constrain_latt** (*list/dict*) – List of lattice parameters we want to
    preserve, e.g. [“alpha”, “c”] or dict with the lattice
    parameter names as keys and values we want the parameters to
    be e.g. {“alpha”: 90, “c”: 2.5}.



* **Returns**

    The most primitive structure found.



#### get_reduced_structure(reduction_algo: Literal['niggli', 'LLL'] = 'niggli')
Get a reduced structure.


* **Parameters**

    **reduction_algo** (*"niggli"** | **"LLL"*) – The lattice reduction algorithm to use.
    Defaults to “niggli”.



#### get_sites_in_sphere(pt: ArrayLike, r: float, include_index: bool = False, include_image: bool = False)
Find all sites within a sphere from the point, including a site (if any)
sitting on the point itself. This includes sites in other periodic
images.

Algorithm:


1. place sphere of radius r in crystal and determine minimum supercell
(parallelepiped) which would contain a sphere of radius r. for this
we need the projection of a_1 on a unit vector perpendicular
to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
determine how many a_1”s it will take to contain the sphere.

Nxmax = r \* length_of_b_1 / (2 Pi)


2. keep points falling within r.


* **Parameters**


    * **pt** (*3x1 array*) – Cartesian coordinates of center of sphere.


    * **r** (*float*) – Radius of sphere.


    * **include_index** (*bool*) – Whether the non-supercell site index
    is included in the returned data


    * **include_image** (*bool*) – Whether to include the supercell image
    is included in the returned data



* **Returns**

    PeriodicNeighbor



#### get_sorted_structure(key: Callable | None = None, reverse: bool = False)
Get a sorted copy of the structure. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### get_space_group_info(symprec=0.01, angle_tolerance=5.0)
Convenience method to quickly get the spacegroup of a structure.


* **Parameters**


    * **symprec** (*float*) – Same definition as in SpacegroupAnalyzer.
    Defaults to 1e-2.


    * **angle_tolerance** (*float*) – Same definition as in SpacegroupAnalyzer.
    Defaults to 5 degrees.



* **Returns**

    spacegroup_symbol, international_number



#### get_symmetric_neighbor_list(r: float, sg: str, unique: bool = False, numerical_tol: float = 1e-08, exclude_self: bool = True)
Similar to ‘get_neighbor_list’ with sites=None, but the neighbors are
grouped by symmetry. The returned values are a tuple of numpy arrays
(center_indices, points_indices, offset_vectors, distances, symmetry_indices).
Atom center_indices[i] has neighbor atom points_indices[i] that is translated
by offset_vectors[i] lattice vectors, and the distance is distances[i].
Symmetry_idx groups the bonds that are related by a symmetry of the provided space
group and symmetry_op is the operation that relates the first bond of the same
symmetry_idx to the respective atom. The first bond maps onto itself via the
Identity. The output is sorted w.r.t. to symmetry_indices. If unique is True only
one of the two bonds connecting two points is given. Out of the two, the bond that
does not reverse the sites is chosen.


* **Parameters**


    * **r** (*float*) – Radius of sphere


    * **sg** (*str/int*) – The spacegroup the symmetry operations of which will be
    used to classify the neighbors. If a string, it will be interpreted
    as one of the notations supported by
    pymatgen.symmetry.groups.Spacegroup. E.g., “R-3c” or “Fm-3m”.
    If an int, it will be interpreted as an international number.
    If None, ‘get_space_group_info’ will be used to determine the
    space group, default to None.


    * **unique** (*bool*) – Whether a bond is given for both, or only a single
    direction is given. The default is False.


    * **numerical_tol** (*float*) – This is a numerical tolerance for distances.
    Sites which are < numerical_tol are determined to be coincident
    with the site. Sites which are r + numerical_tol away is deemed
    to be within r from the site. The default of 1e-8 should be
    ok in most instances.


    * **exclude_self** (*bool*) – whether to exclude atom neighboring with itself within
    numerical tolerance distance, default to True



* **Returns**

    (center_indices, points_indices, offset_vectors, distances,

        symmetry_indices, symmetry_ops)




* **Return type**

    tuple



#### interpolate(end_structure: IStructure | Structure, nimages: int | Iterable = 10, interpolate_lattices: bool = False, pbc: bool = True, autosort_tol: float = 0)
Interpolate between this structure and end_structure. Useful for
construction of NEB inputs.


* **Parameters**


    * **end_structure** (*Structure*) – structure to interpolate between this
    structure and end.


    * **nimages** (*int**,**list*) – No. of interpolation images or a list of
    interpolation images. Defaults to 10 images.


    * **interpolate_lattices** (*bool*) – Whether to interpolate the lattices.
    Interpolates the lengths and angles (rather than the matrix)
    so orientation may be affected.


    * **pbc** (*bool*) – Whether to use periodic boundary conditions to find
    the shortest path between endpoints.


    * **autosort_tol** (*float*) – A distance tolerance in angstrom in
    which to automatically sort end_structure to match to the
    closest points in this particular structure. This is usually
    what you want in a NEB calculation. 0 implies no sorting.
    Otherwise, a 0.5 value usually works pretty well.



* **Returns**

    List of interpolated structures. The starting and ending
    structures included as the first and last structures respectively.
    A total of (nimages + 1) structures are returned.



#### _property_ is_3d_periodic(_: boo_ )
True if the Lattice is periodic in all directions.


#### _property_ lattice(_: Lattic_ )
Lattice of the structure.


#### matches(other: IStructure | Structure, anonymous: bool = False, \*\*kwargs)
Check whether this structure is similar to another structure.
Basically a convenience method to call structure matching.


* **Parameters**


    * **other** (*IStructure/Structure*) – Another structure.


    * **anonymous** (*bool*) – Whether to use anonymous structure matching which allows distinct
    species in one structure to map to another.


    * **\*\*kwargs** – Same

    ```
    **
    ```

    kwargs as in
    pymatgen.analysis.structure_matcher.StructureMatcher.



* **Returns**

    True if the structures are similar under some affine transformation.



* **Return type**

    bool



#### _property_ pbc(_: tuple[bool, bool, bool_ )
Returns the periodicity of the structure.


#### _property_ properties(_: dic_ )
Properties associated with the whole Structure. Note that this information is
only guaranteed to be saved if serializing to native pymatgen output formats (JSON/YAML).


#### to(filename: str = '', fmt: str = '', \*\*kwargs)
Outputs the structure to a file or string.


* **Parameters**


    * **filename** (*str*) – If provided, output will be written to a file. If
    fmt is not specified, the format is determined from the
    filename. Defaults is None, i.e. string output.


    * **fmt** (*str*) – Format to output to. Defaults to JSON unless filename
    is provided. If fmt is specifies, it overrides whatever the
    filename is. Options include “cif”, “poscar”, “cssr”, “json”,
    “xsf”, “mcsqs”, “prismatic”, “yaml”, “fleur-inpgen”.
    Non-case sensitive.


    * **\*\*kwargs** – Kwargs passthru to relevant methods. E.g., This allows
    the passing of parameters like symprec to the
    CifWriter.__init__ method for generation of symmetric CIFs.



* **Returns**

    String representation of molecule in given format. If a filename

        is provided, the same string is written to the file.




* **Return type**

    str



#### unset_charge()
Reset the charge to None, i.e., computed dynamically based on oxidation states.


#### _property_ volume(_: floa_ )
Returns the volume of the structure in Angstrom^3.


### _class_ Molecule(species: Sequence[SpeciesLike], coords: Sequence[ArrayLike], charge: float = 0.0, spin_multiplicity: int | None = None, validate_proximity: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, charge_spin_check: bool = True, properties: dict | None = None)
Bases: `IMolecule`, `MutableSequence`

Mutable Molecule. It has all the methods in IMolecule, but in addition,
it allows a user to perform edits on the molecule.

Creates a MutableMolecule.


* **Parameters**


    * **species** – list of atomic species. Possible kinds of input include a
    list of dict of elements/species and occupancies, a List of
    elements/specie specified as actual Element/Species, Strings
    (“Fe”, “Fe2+”) or atomic numbers (1,56).


    * **coords** (*3x1 array*) – list of Cartesian coordinates of each species.


    * **charge** (*float*) – Charge for the molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multiplicity for molecule.
    Defaults to None, which means that the spin multiplicity is
    set to 1 if the molecule has no unpaired electrons and to 2
    if there are unpaired electrons.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 1 Ang apart. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as
    a dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The
    sequences have to be the same length as the atomic species
    and fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **charge_spin_check** (*bool*) – Whether to check that the charge and
    spin multiplicity are compatible with each other. Defaults
    to True.


    * **properties** (*dict*) – dictionary containing properties associated
    with the whole molecule.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### append(species: CompositionLike, coords: ArrayLike, validate_proximity: bool = False, properties: dict | None = None)
Appends a site to the molecule.


* **Parameters**


    * **species** – Species of inserted site


    * **coords** – Coordinates of inserted site


    * **validate_proximity** (*bool*) – Whether to check if inserted site is
    too close to an existing site. Defaults to False.


    * **properties** (*dict*) – A dict of properties for the Site.



* **Returns**

    New molecule with inserted site.



#### apply_operation(symmop: SymmOp)
Apply a symmetry operation to the molecule.


* **Parameters**

    **symmop** (*SymmOp*) – Symmetry operation to apply.



#### calculate(calculator: str | Calculator = 'gfn2-xtb', verbose: bool = False)
Performs an ASE calculation.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “gfn2-xtb”.
    Defaults to ‘gfn2-xtb’.


    * **verbose** (*bool*) – whether to print stdout. Defaults to False.



* **Returns**

    ASE Calculator instance with a results attribute containing the output.



* **Return type**

    Calculator



#### insert(idx: int, species: CompositionLike, coords: ArrayLike, validate_proximity: bool = False, properties: dict | None = None, label: str | None = None)
Insert a site to the molecule.


* **Parameters**


    * **idx** (*int*) – Index to insert site


    * **species** – species of inserted site


    * **coords** (*3x1 array*) – coordinates of inserted site


    * **validate_proximity** (*bool*) – Whether to check if inserted site is
    too close to an existing site. Defaults to True.


    * **properties** (*dict*) – Dict of properties for the Site.


    * **label** (*str*) – Label of inserted site



* **Returns**

    New molecule with inserted site.



#### perturb(distance: float)
Performs a random perturbation of the sites in a structure to break
symmetries.


* **Parameters**

    **distance** (*float*) – Distance in angstroms by which to perturb each
    site.



#### relax(calculator: str | Calculator = 'gfn2-xtb', optimizer: str | Optimizer = 'FIRE', steps: int = 500, fmax: float = 0.1, opt_kwargs: dict | None = None, return_trajectory: bool = False, verbose: bool = False)
Performs a molecule relaxation using an ASE calculator.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “gfn2-xtb”.
    Defaults to ‘gfn2-xtb’.


    * **optimizer** (*str*) – name of the ASE optimizer class to use


    * **steps** (*int*) – max number of steps for relaxation. Defaults to 500.


    * **fmax** (*float*) – total force tolerance for relaxation convergence.
    Defaults to 0.1 eV/A.


    * **opt_kwargs** (*dict*) – kwargs for the ASE optimizer class.


    * **return_trajectory** (*bool*) – Whether to return the trajectory of relaxation.
    Defaults to False.


    * **verbose** (*bool*) – whether to print out relaxation steps. Defaults to False.



* **Returns**

    Relaxed Molecule or if return_trajectory=True,

        2-tuple of Molecule and ASE TrajectoryObserver.




* **Return type**

    Molecule | tuple[Molecule, Trajectory]



#### remove_sites(indices: Sequence[int])
Delete sites with at indices.


* **Parameters**

    **indices** – Sequence of indices of sites to delete.



#### remove_species(species: Sequence[SpeciesLike])
Remove all occurrences of a species from a molecule.


* **Parameters**

    **species** – Species to remove.



#### rotate_sites(indices: Sequence[int] | None = None, theta: float = 0.0, axis: ArrayLike | None = None, anchor: ArrayLike | None = None)
Rotate specific sites by some angle around vector at anchor.


* **Parameters**


    * **indices** (*list*) – List of site indices on which to perform the
    translation.


    * **theta** (*float*) – Angle in radians


    * **axis** (*3x1 array*) – Rotation axis vector.


    * **anchor** (*3x1 array*) – Point of rotation.



#### set_charge_and_spin(charge: float, spin_multiplicity: int | None = None)
Set the charge and spin multiplicity.


* **Parameters**


    * **charge** (*int*) – Charge for the molecule. Defaults to 0.


    * **spin_multiplicity** (*int*) – Spin multiplicity for molecule.
    Defaults to None, which means that the spin multiplicity is
    set to 1 if the molecule has no unpaired electrons and to 2
    if there are unpaired electrons.



#### substitute(index: int, func_group: IMolecule | Molecule | str, bond_order: int = 1)
Substitute atom at index with a functional group.


* **Parameters**


    * **index** (*int*) – Index of atom to substitute.


    * **func_group** – Substituent molecule. There are two options:


        1. Providing an actual molecule as the input. The first atom
    must be a DummySpecies X, indicating the position of
    nearest neighbor. The second atom must be the next
    nearest atom. For example, for a methyl group
    substitution, func_group should be X-CH3, where X is the
    first site and C is the second site. What the code will
    do is to remove the index site, and connect the nearest
    neighbor to the C atom in CH3. The X-C bond indicates the
    directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the
    relevant template in func_groups.json.



    * **bond_order** (*int*) – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.



#### translate_sites(indices: Sequence[int] | None = None, vector: ArrayLike | None = None)
Translate specific sites by some vector, keeping the sites within the
unit cell.


* **Parameters**


    * **indices** (*list*) – List of site indices on which to perform the
    translation.


    * **vector** (*3x1 array*) – Translation vector for sites.



### _class_ Neighbor(species: Composition, coords: np.ndarray, properties: dict | None = None, nn_distance: float = 0.0, index: int = 0, label: str | None = None)
Bases: `Site`

Simple Site subclass to contain a neighboring atom that skips all the unnecessary checks for speed. Can be
used as a fixed-length tuple of size 3 to retain backwards compatibility with past use cases.

> (site, nn_distance, index).

In future, usage should be to call attributes, e.g., Neighbor.index, Neighbor.distance, etc.


* **Parameters**


    * **species** – Same as Site


    * **coords** – Same as Site, but must be fractional.


    * **properties** – Same as Site


    * **nn_distance** – Distance to some other Site.


    * **index** – Index within structure.


    * **label** – Label for the site. Defaults to None.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _species(_: Compositio_ )

#### as_dict()
Note that method calls the super of Site, which is MSONable itself.


#### coords(_: ndarra_ )

#### _classmethod_ from_dict(dct: dict)
Returns a Neighbor from a dict.


* **Parameters**

    **dct** – MSONable dict format.



* **Returns**

    Neighbor



#### properties(_: dic_ )

### _class_ PeriodicNeighbor(species: Composition, coords: np.ndarray, lattice: Lattice, properties: dict | None = None, nn_distance: float = 0.0, index: int = 0, image: tuple = (0, 0, 0), label: str | None = None)
Bases: `PeriodicSite`

Simple PeriodicSite subclass to contain a neighboring atom that skips all
the unnecessary checks for speed. Can be used as a fixed-length tuple of
size 4 to retain backwards compatibility with past use cases.

> (site, distance, index, image).

In future, usage should be to call attributes, e.g., PeriodicNeighbor.index,
PeriodicNeighbor.distance, etc.


* **Parameters**


    * **species** (*Composition*) – Same as PeriodicSite


    * **coords** (*np.ndarray*) – Same as PeriodicSite, but must be fractional.


    * **lattice** (*Lattice*) – Same as PeriodicSite


    * **properties** (*dict**, **optional*) – Same as PeriodicSite. Defaults to None.


    * **nn_distance** (*float**, **optional*) – Distance to some other Site.. Defaults to 0.0.


    * **index** (*int**, **optional*) – Index within structure.. Defaults to 0.


    * **image** (*tuple**, **optional*) – PeriodicImage. Defaults to (0, 0, 0).


    * **label** (*str**, **optional*) – Label for the site. Defaults to None.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _coords(_: np.ndarray | Non_ )

#### _frac_coords(_: np.ndarra_ )

#### _lattice(_: Lattic_ )

#### _species(_: Compositio_ )

#### as_dict()
Note that method calls the super of Site, which is MSONable itself.


#### _property_ coords(_: ndarra_ )
Cartesian coords.


#### _classmethod_ from_dict(d: dict)
Returns a PeriodicNeighbor from a dict.


* **Parameters**

    **d** – MSONable dict format.



* **Returns**

    PeriodicNeighbor



#### properties(_: dic_ )

### _class_ SiteCollection()
Bases: `Sequence`

Basic SiteCollection. Essentially a sequence of Sites or PeriodicSites.
This serves as a base class for Molecule (a collection of Site, i.e., no
periodicity) and Structure (a collection of PeriodicSites, i.e.,
periodicity). Not meant to be instantiated directly.


#### DISTANCE_TOLERANCE(_ = 0._ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _calculate(calculator: str | Calculator, verbose: bool = False)
Performs an ASE calculation.


* **Parameters**


    * **calculator** (*str** | **Calculator*) – An ASE Calculator or a string from the following case-insensitive
    options: “m3gnet”, “gfn2-xtb”, “chgnet”.


    * **verbose** (*bool*) – whether to print stdout. Defaults to False.
    Has no effect when calculator==’chgnet’.



* **Returns**

    Structure or Molecule with new calc attribute containing

        result of ASE calculation.




* **Return type**

    Structure | Molecule



#### _prep_calculator(calculator: Literal['m3gnet', 'gfn2-xtb'] | Calculator, \*\*params)
Convert string name of special ASE calculators into ASE calculator objects.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “m3gnet”,
    “gfn2-xtb”.


    * **\*\*params** – Parameters for the calculator.



* **Returns**

    ASE calculator object.



* **Return type**

    Calculator



#### _properties(_: dic_ )

#### _relax(calculator: str | Calculator, relax_cell: bool = True, optimizer: str | Optimizer = 'FIRE', steps: int = 500, fmax: float = 0.1, stress_weight: float = 0.01, opt_kwargs: dict | None = None, return_trajectory: bool = False, verbose: bool = False)
Performs a structure relaxation using an ASE calculator.


* **Parameters**


    * **calculator** (*str** | **ase.Calculator*) – An ASE Calculator or a string from the following options: “M3GNet”,
    “gfn2-xtb”.


    * **relax_cell** (*bool*) – whether to relax the lattice cell. Defaults to True.


    * **optimizer** (*str*) – name of the ASE optimizer class to use


    * **steps** (*int*) – max number of steps for relaxation. Defaults to 500.


    * **fmax** (*float*) – total force tolerance for relaxation convergence.
    Here fmax is a sum of force and stress forces. Defaults to 0.1.


    * **stress_weight** (*float*) – the stress weight for relaxation with M3GNet.
    Defaults to 0.01.


    * **opt_kwargs** (*dict*) – kwargs for the ASE optimizer class.


    * **return_trajectory** (*bool*) – Whether to return the trajectory of relaxation.
    Defaults to False.


    * **verbose** (*bool*) – whether to print stdout. Defaults to False.



* **Returns**

    Relaxed structure or molecule



* **Return type**

    Structure | Molecule



#### add_oxidation_state_by_element(oxidation_states: dict[str, float])
Add oxidation states.


* **Parameters**

    **oxidation_states** (*dict*) – Dict of oxidation states.
    E.g., {“Li”:1, “Fe”:2, “P”:5, “O”:-2}



* **Raises**

    **ValueError if oxidation states are not specified for all elements.** –



#### add_oxidation_state_by_guess(\*\*kwargs)
Decorates the structure with oxidation state, guessing
using Composition.oxi_state_guesses().


* **Parameters**

    **\*\*kwargs** – parameters to pass into oxi_state_guesses()



#### add_oxidation_state_by_site(oxidation_states: list[float])
Add oxidation states to a structure by site.


* **Parameters**

    **oxidation_states** (*list**[**float**]*) – List of oxidation states.
    E.g. [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2]



#### add_site_property(property_name: str, values: Sequence | np.ndarray)
Adds a property to a site. Note: This is the preferred method
for adding magnetic moments, selective dynamics, and related
site-specific properties to a structure/molecule object.

### Examples

structure.add_site_property(“magmom”, [1.0, 0.0])
structure.add_site_property(“selective_dynamics”, [[True, True, True], [False, False, False]])


* **Parameters**


    * **property_name** (*str*) – The name of the property to add.


    * **values** (*list*) – A sequence of values. Must be same length as
    number of sites.



#### add_spin_by_element(spins: dict[str, float])
Add spin states to structure.


* **Parameters**

    **spins** (*dict*) – Dict of spins associated with elements or species,
    e.g. {“Ni”:+5} or {“Ni2+”:5}



#### add_spin_by_site(spins: Sequence[float])
Add spin states to structure by site.


* **Parameters**

    **spins** (*list*) – e.g. [+5, -5, 0, 0]



#### _property_ atomic_numbers(_: tuple[int, ..._ )
List of atomic numbers.


#### _property_ cart_coords(_: ndarra_ )
Returns an np.array of the Cartesian coordinates of sites in the structure.


#### _property_ charge(_: floa_ )
Returns the net charge of the structure based on oxidation states. If
Elements are found, a charge of 0 is assumed.


#### _property_ composition(_: Compositio_ )
Returns the structure’s corresponding Composition object.


#### _abstract_ copy()
Returns a copy of itself. Concrete subclasses should implement this
method.


#### _property_ distance_matrix(_: ndarra_ )
Returns the distance matrix between all sites in the structure. For
periodic structures, this is overwritten to return the nearest image
distance.


#### _property_ elements(_: list[Element | Species | DummySpecies_ )
Returns the elements in the structure as a list of Element objects.


#### extract_cluster(target_sites: list[Site], \*\*kwargs)
Extracts a cluster of atoms based on bond lengths.


* **Parameters**


    * **target_sites** (*list**[**Site**]*) – Initial sites from which to nucleate cluster.


    * **\*\*kwargs** – kwargs passed through to CovalentBond.is_bonded.



* **Returns**

    list[Site/PeriodicSite] Cluster of atoms.



#### _property_ formula(_: st_ )
Returns the formula as a string.


#### _abstract classmethod_ from_file(filename: str)
Reads in SiteCollection from a filename.


#### _abstract classmethod_ from_str(input_string: str, fmt: Any)
Reads in SiteCollection from a string.


#### get_angle(i: int, j: int, k: int)
Returns angle specified by three sites.


* **Parameters**


    * **i** – 1st site index


    * **j** – 2nd site index


    * **k** – 3rd site index



* **Returns**

    Angle in degrees.



#### get_dihedral(i: int, j: int, k: int, l: int)
Returns dihedral angle specified by four sites.


* **Parameters**


    * **i** – 1st site index


    * **j** – 2nd site index


    * **k** – 3rd site index


    * **l** – 4th site index



* **Returns**

    Dihedral angle in degrees.



#### _abstract_ get_distance(i: int, j: int)
Returns distance between sites at index i and j.


* **Parameters**


    * **i** – 1st site index


    * **j** – 2nd site index



* **Returns**

    Distance between sites at index i and index j.



#### group_by_types()
Iterate over species grouped by type.


#### indices_from_symbol(symbol: str)
Returns a tuple with the sequential indices of the sites
that contain an element with the given chemical symbol.


#### _property_ is_ordered(_: boo_ )
Checks if structure is ordered, meaning no partial occupancies in any
of the sites.


#### is_valid(tol: float = 0.5)
True if SiteCollection does not contain atoms that are too close
together. Note that the distance definition is based on type of
SiteCollection. Cartesian distances are used for non-periodic
Molecules, while PBC is taken into account for periodic structures.


* **Parameters**

    **tol** (*float*) – Distance tolerance. Default is 0.5 Angstrom, which is fairly large.



* **Returns**

    (bool) True if SiteCollection does not contain atoms that are too close together.



#### _property_ labels(_: list[str_ )
Return site labels as a list.


#### _property_ ntypesp(_: in_ )
Number of types of atoms.


#### _property_ num_sites(_: in_ )
Number of sites.


#### remove_oxidation_states()
Removes oxidation states from a structure.


#### remove_site_property(property_name: str)
Removes a property to a site.


* **Parameters**

    **property_name** (*str*) – The name of the property to remove.



#### remove_spin()
Remove spin states from structure.


#### replace_species(species_mapping: dict[SpeciesLike, SpeciesLike | dict[SpeciesLike, float]], in_place: bool = True)
Swap species.


* **Parameters**


    * **species_mapping** (*dict*) – Species to swap. Species can be elements too. E.g.,
    {Element(“Li”): Element(“Na”)} performs a Li for Na substitution. The second species can
    be a sp_and_occu dict. For example, a site with 0.5 Si that is passed the mapping
    {Element(‘Si): {Element(‘Ge’): 0.75, Element(‘C’): 0.25} } will have .375 Ge and .125 C.


    * **in_place** (*bool*) – Whether to perform the substitution in place or modify a copy.
    Defaults to True.



#### _property_ site_properties(_: dict[str, Sequence_ )
(-4, 4)}.


* **Type**

    Returns the site properties as a dict of sequences. E.g. {“magmom”



* **Type**

    (5, -5), “charge”



#### _property_ sites(_: list[Site_ )
Returns an iterator for the sites in the Structure.


#### _property_ species(_: list[Element | Species_ )
Only works for ordered structures.


* **Raises**

    **AttributeError** – If structure is disordered.



* **Returns**

    ([Species]) List of species at each site of the structure.



#### _property_ species_and_occu(_: list[Composition_ )
List of species and occupancies at each site of the structure.


#### _property_ symbol_set(_: tuple[str, ..._ )
Tuple with the set of chemical symbols.
Note that len(symbol_set) == len(types_of_specie).


#### _abstract_ to(filename: str = '', fmt: str = '')
Generates string representations (cif, json, poscar, ….) of SiteCollections (e.g.,
molecules / structures). Should return str or None if written to a file.


#### _property_ types_of_specie(_: tuple[Element | Species | DummySpecies_ )
Specie->Species rename. Maintained for backwards compatibility.


#### _property_ types_of_species(_: tuple[Element | Species | DummySpecies_ )
List of types of specie.


### _class_ Structure(lattice: ArrayLike | Lattice, species: Sequence[CompositionLike], coords: Sequence[ArrayLike], charge: float | None = None, validate_proximity: bool = False, to_unit_cell: bool = False, coords_are_cartesian: bool = False, site_properties: dict | None = None, labels: Sequence[str | None] | None = None, properties: dict | None = None)
Bases: `IStructure`, `MutableSequence`

Mutable version of structure.

Create a periodic structure.


* **Parameters**


    * **lattice** – The lattice, either as a pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** – List of species on each site. Can take in flexible input,
    including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of
    each species.


    * **charge** (*int*) – overall charge of the structure. Defaults to behavior
    in SiteCollection where total charge is the sum of the oxidation
    states.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **to_unit_cell** (*bool*) – Whether to map all sites into the unit cell,
    i.e., fractional coords between 0 and 1. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **labels** (*list**[**str**]*) – Labels associated with the sites as a
    list of strings, e.g. [‘Li1’, ‘Li2’]. Must have the same
    length as the species and fractional coords. Defaults to
    None for no labels.


    * **properties** (*dict*) – Properties associated with the whole structure.
    Will be serialized when writing the structure to JSON or YAML but is
    lost when converting to other formats.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### _sites(_: tuple[PeriodicSite, ..._ )

#### append(species: CompositionLike, coords: ArrayLike, coords_are_cartesian: bool = False, validate_proximity: bool = False, properties: dict | None = None)
Append a site to the structure.


* **Parameters**


    * **species** – Species of inserted site


    * **coords** (*3x1 array*) – Coordinates of inserted site


    * **coords_are_cartesian** (*bool*) – Whether coordinates are cartesian.
    Defaults to False.


    * **validate_proximity** (*bool*) – Whether to check if inserted site is
    too close to an existing site. Defaults to False.


    * **properties** (*dict*) – Properties of the site.



* **Returns**

    New structure with inserted site.



#### apply_operation(symmop: SymmOp, fractional: bool = False)
Apply a symmetry operation to the structure in place and return the modified
structure. The lattice is operated on by the rotation matrix only.
Coords are operated in full and then transformed to the new lattice.


* **Parameters**


    * **symmop** (*SymmOp*) – Symmetry operation to apply.


    * **fractional** (*bool*) – Whether the symmetry operation is applied in
    fractional space. Defaults to False, i.e., symmetry operation
    is applied in Cartesian coordinates.



* **Returns**

    post-operation structure



* **Return type**

    Structure



#### apply_strain(strain: ArrayLike)
Apply a strain to the lattice.


* **Parameters**

    **strain** (*float** or **list*) – Amount of strain to apply. Can be a float,
    or a sequence of 3 numbers. E.g., 0.01 means all lattice
    vectors are increased by 1%. This is equivalent to calling
    modify_lattice with a lattice with lattice parameters that
    are 1% larger.



#### calculate(calculator: str | Calculator = 'm3gnet', verbose: bool = False)
Performs an ASE calculation.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “m3gnet”.
    Defaults to ‘m3gnet’, i.e. the M3GNet universal potential.


    * **verbose** (*bool*) – whether to print stdout. Defaults to False.



* **Returns**

    ASE Calculator instance with a results attribute containing the output.



* **Return type**

    Calculator



#### _classmethod_ from_prototype(prototype: str, species: Sequence, \*\*kwargs)
Method to rapidly construct common prototype structures.


* **Parameters**


    * **prototype** – Name of prototype. E.g., cubic, rocksalt, perovksite etc.


    * **species** – List of species corresponding to symmetrically distinct sites.


    * **\*\*kwargs** – Lattice parameters, e.g., a = 3.0, b = 4, c = 5. Only the required lattice parameters need to be
    specified. For example, if it is a cubic prototype, only a needs to be specified.



* **Returns**

    Structure



#### insert(idx: int, species: CompositionLike, coords: ArrayLike, coords_are_cartesian: bool = False, validate_proximity: bool = False, properties: dict | None = None, label: str | None = None)
Insert a site to the structure.


* **Parameters**


    * **idx** (*int*) – Index to insert site


    * **species** (*species-like*) – Species of inserted site


    * **coords** (*3x1 array*) – Coordinates of inserted site


    * **coords_are_cartesian** (*bool*) – Whether coordinates are cartesian.
    Defaults to False.


    * **validate_proximity** (*bool*) – Whether to check if inserted site is too close to
    an existing site. Controlled by self.DISTANCE_TOLERANCE. Defaults to False.


    * **properties** (*dict*) – Properties associated with the site.


    * **label** (*str*) – Label associated with the site.



* **Returns**

    New structure with inserted site.



#### _property_ lattice(_: Lattic_ )
Lattice associated with structure.


#### make_supercell(scaling_matrix: ArrayLike, to_unit_cell: bool = True, in_place: bool = True)
Create a supercell.


* **Parameters**


    * **scaling_matrix** (*ArrayLike*) – A scaling matrix for transforming the lattice
    vectors. Has to be all integers. Several options are possible:


        1. A full 3x3 scaling matrix defining the linear combination
    the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
    1]] generates a new structure with lattice vectors a’ =
    2a + b, b’ = 3b, c’ = c where a, b, and c are the lattice
    vectors of the original structure.


        2. An sequence of three scaling factors. E.g., [2, 1, 1]
    specifies that the supercell should have dimensions 2a x b x
    c.


        3. A number, which simply scales all lattice vectors by the
    same factor.



    * **to_unit_cell** (*bool*) – Whether or not to fold sites back into the unit cell
    if they have fractional coords > 1. Defaults to True.


    * **in_place** (*bool*) – Whether to perform the operation in-place or to return
    a new Structure object. Defaults to True.



* **Returns**

    self if in_place is True else self.copy() after making supercell



* **Return type**

    Structure



#### merge_sites(tol: float = 0.01, mode: Literal['sum', 'delete', 'average'] = 'sum')
Merges sites (adding occupancies) within tol of each other.
Removes site properties.


* **Parameters**


    * **tol** (*float*) – Tolerance for distance to merge sites.


    * **mode** (*'sum'** | **'delete'** | **'average'*) – “delete” means duplicate sites are
    deleted. “sum” means the occupancies are summed for the sites.
    “average” means that the site is deleted but the properties are averaged
    Only first letter is considered.



#### perturb(distance: float, min_distance: float | None = None)
Performs a random perturbation of the sites in a structure to break
symmetries.


* **Parameters**


    * **distance** (*float*) – Distance in angstroms by which to perturb each
    site.


    * **min_distance** (*None**, **int**, or **float*) – if None, all displacements will
    be equal amplitude. If int or float, perturb each site a
    distance drawn from the uniform distribution between
    ‘min_distance’ and ‘distance’.



#### relax(calculator: str | Calculator = 'm3gnet', relax_cell: bool = True, optimizer: str | Optimizer = 'FIRE', steps: int = 500, fmax: float = 0.1, stress_weight: float = 0.01, opt_kwargs: dict | None = None, return_trajectory: bool = False, verbose: bool = False)
Performs a crystal structure relaxation using an ASE calculator.


* **Parameters**


    * **calculator** – An ASE Calculator or a string from the following options: “m3gnet”.
    Defaults to ‘m3gnet’, i.e. the M3GNet universal potential.


    * **relax_cell** (*bool*) – whether to relax the lattice cell. Defaults to True.


    * **optimizer** (*str*) – name of the ASE optimizer class to use


    * **steps** (*int*) – max number of steps for relaxation. Defaults to 500.


    * **fmax** (*float*) – total force tolerance for relaxation convergence.
    Here fmax is a sum of force and stress forces. Defaults to 0.1.


    * **stress_weight** (*float*) – the stress weight for relaxation with M3GNet.
    Defaults to 0.01.


    * **opt_kwargs** (*dict*) – kwargs for the ASE optimizer class.


    * **return_trajectory** (*bool*) – Whether to return the trajectory of relaxation.
    Defaults to False.


    * **verbose** (*bool*) – whether to print out relaxation steps. Defaults to False.



* **Returns**

    Relaxed structure or if return_trajectory=True,

        2-tuple of Structure and matgl TrajectoryObserver.




* **Return type**

    Structure | tuple[Structure, Trajectory]



#### remove_sites(indices: Sequence[int | None])
Delete sites with at indices.


* **Parameters**

    **indices** – Sequence of indices of sites to delete.



#### remove_species(species: Sequence[SpeciesLike])
Remove all occurrences of several species from a structure.


* **Parameters**

    **species** – Sequence of species to remove, e.g., [“Li”, “Na”].



#### replace(idx: int, species: CompositionLike, coords: ArrayLike | None = None, coords_are_cartesian: bool = False, properties: dict | None = None, label: str | None = None)
Replace a single site. Takes either a species or a dict of species and
occupations.


* **Parameters**


    * **idx** (*int*) – Index of the site in the sites list.


    * **species** (*species-like*) – Species of replacement site


    * **coords** (*3x1 array*) – Coordinates of replacement site. If None,
    the current coordinates are assumed.


    * **coords_are_cartesian** (*bool*) – Whether coordinates are cartesian.
    Defaults to False.


    * **properties** (*dict*) – Properties associated with the site.


    * **label** (*str*) – Label associated with the site.



#### rotate_sites(indices: list[int] | None = None, theta: float = 0.0, axis: ArrayLike | None = None, anchor: ArrayLike | None = None, to_unit_cell: bool = True)
Rotate specific sites by some angle around vector at anchor.


* **Parameters**


    * **indices** (*list*) – List of site indices on which to perform the
    translation.


    * **theta** (*float*) – Angle in radians


    * **axis** (*3x1 array*) – Rotation axis vector.


    * **anchor** (*3x1 array*) – Point of rotation.


    * **to_unit_cell** (*bool*) – Whether new sites are transformed to unit
    cell



#### scale_lattice(volume: float)
Performs a scaling of the lattice vectors so that length proportions
and angles are preserved.


* **Parameters**

    **volume** (*float*) – New volume of the unit cell in A^3.



#### set_charge(new_charge: float = 0.0)
Sets the overall structure charge.


* **Parameters**

    **new_charge** (*float*) – new charge to set



#### sort(key: Callable | None = None, reverse: bool = False)
Sort a structure in place. The parameters have the same meaning as in
list.sort. By default, sites are sorted by the electronegativity of
the species. The difference between this method and
get_sorted_structure (which also works in IStructure) is that the
latter returns a new Structure, while this just sorts the Structure
in place.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### substitute(index: int, func_group: IMolecule | Molecule | str, bond_order: int = 1)
Substitute atom at index with a functional group.


* **Parameters**


    * **index** (*int*) – Index of atom to substitute.


    * **func_group** – Substituent molecule. There are two options:


        1. Providing an actual Molecule as the input. The first atom
    must be a DummySpecies X, indicating the position of
    nearest neighbor. The second atom must be the next
    nearest atom. For example, for a methyl group
    substitution, func_group should be X-CH3, where X is the
    first site and C is the second site. What the code will
    do is to remove the index site, and connect the nearest
    neighbor to the C atom in CH3. The X-C bond indicates the
    directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the
    relevant template in func_groups.json.



    * **bond_order** (*int*) – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.



#### translate_sites(indices: int | Sequence[int], vector: ArrayLike, frac_coords: bool = True, to_unit_cell: bool = True)
Translate specific sites by some vector, keeping the sites within the
unit cell.


* **Parameters**


    * **indices** – Integer or List of site indices on which to perform the
    translation.


    * **vector** – Translation vector for sites.


    * **frac_coords** (*bool*) – Whether the vector corresponds to fractional or
    Cartesian coordinates.


    * **to_unit_cell** (*bool*) – Whether new sites are transformed to unit
    cell



### _exception_ StructureError()
Bases: `Exception`

Exception class for Structure.
Raised when the structure has problems, e.g., atoms that are too close.

## pymatgen.core.surface module

This module implements representations of slabs and surfaces + algorithms for generating them.

If you use this module, please consider citing the following work:

>
> 1. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,

> S. P. Ong, “Surface Energies of Elemental Crystals”, Scientific Data,
> 2016, 3:160080, doi: 10.1038/sdata.2016.80.

> Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
> Surface Science, 2013, 617, 53-59, doi:10.1016/j.susc.2013.05.016.


### _class_ ReconstructionGenerator(initial_structure, min_slab_size, min_vacuum_size, reconstruction_name)
Bases: `object`

This class takes in a pre-defined dictionary specifying the parameters
need to build a reconstructed slab such as the SlabGenerator parameters,
transformation matrix, sites to remove/add and slab/vacuum size. It will
then use the formatted instructions provided by the dictionary to build
the desired reconstructed slab from the initial structure.


#### slabgen_params()
Parameters for the SlabGenerator.


* **Type**

    dict



#### trans_matrix()
A 3x3 transformation matrix to generate the reconstructed
slab. Only the a and b lattice vectors are actually changed while the c vector remains
the same. This matrix is what the Wood’s notation is based on.


* **Type**

    np.ndarray



#### reconstruction_json()
The full json or dictionary containing the instructions for
building the reconstructed slab.


* **Type**

    dict



#### termination()
The index of the termination of the slab.


* **Type**

    int


Generates reconstructed slabs from a set of instructions

    specified by a dictionary or json file.


* **Parameters**


    * **initial_structure** (*Structure*) – Initial input structure. Note
    that to ensure that the miller indices correspond to usual
    crystallographic definitions, you should supply a conventional
    unit cell structure.


    * **min_slab_size** (*float*) – In Angstroms


    * **min_vacuum_size** (*float*) – In Angstroms


    * **reconstruction_name** (*str*) – Name of the dict containing the instructions
    for building a reconstructed slab. The dictionary can contain
    any item the creator deems relevant, however any instructions
    archived in pymatgen for public use needs to contain the
    following keys and items to ensure compatibility with the
    ReconstructionGenerator:

    > ”name” (str): A descriptive name for the type of

    >     reconstruction. Typically the name will have the type
    >     of structure the reconstruction is for, the Miller
    >     index, and Wood’s notation along with anything to
    >     describe the reconstruction: e.g.:
    >     “fcc_110_missing_row_1x2”

    > ”description” (str): A longer description of your

    >     reconstruction. This is to help future contributors who
    >     want to add other types of reconstructions to the
    >     archive on pymatgen to check if the reconstruction
    >     already exists. Please read the descriptions carefully
    >     before adding a new type of reconstruction to ensure it
    >     is not in the archive yet.

    > ”reference” (str): Optional reference to where the

    >     reconstruction was taken from or first observed.

    > ”spacegroup” (dict): e.g. {“symbol”: “Fm-3m”, “number”: 225}

    >     Indicates what kind of structure is this reconstruction.

    > ”miller_index” ([h,k,l]): Miller index of your reconstruction
    > “Woods_notation” (str): For a reconstruction, the a and b

    > > lattice may change to accommodate the symmetry of the
    > > reconstruction. This notation indicates the change in
    > > the vectors relative to the primitive (p) or
    > > conventional (c) slab cell. E.g. p(2x1):

    > > Wood, E. A. (1964). Vocabulary of surface
    > > crystallography. Journal of Applied Physics, 35(4),
    > > 1306-1312.

    > ”transformation_matrix” (numpy array): A 3x3 matrix to

    >     transform the slab. Only the a and b lattice vectors
    >     should change while the c vector remains the same.

    > ”SlabGenerator_parameters” (dict): A dictionary containing

    >     the parameters for the SlabGenerator class excluding the
    >     miller_index, min_slab_size and min_vac_size as the
    >     Miller index is already specified and the min_slab_size
    >     and min_vac_size can be changed regardless of what type
    >     of reconstruction is used. Having a consistent set of
    >     SlabGenerator parameters allows for the instructions to
    >     be reused to consistently build a reconstructed slab.

    > ”points_to_remove” (list of coords): A list of sites to

    >     remove where the first two indices are fraction (in a
    >     and b) and the third index is in units of 1/d (in c).

    > ”points_to_add” (list of frac_coords): A list of sites to add

    >     where the first two indices are fraction (in a an b) and
    >     the third index is in units of 1/d (in c).

    > ”base_reconstruction” (dict): Option to base a reconstruction on

    >     an existing reconstruction model also exists to easily build
    >     the instructions without repeating previous work. E.g. the
    >     alpha reconstruction of halites is based on the octopolar
    >     reconstruction but with the topmost atom removed. The dictionary
    >     for the alpha reconstruction would therefore contain the item
    >     “reconstruction_base”: “halite_111_octopolar_2x2”, and
    >     additional sites for “points_to_remove” and “points_to_add”
    >     can be added to modify this reconstruction.

    > For “points_to_remove” and “points_to_add”, the third index for

    >     the c vector is in units of 1/d where d is the spacing
    >     between atoms along hkl (the c vector) and is relative to
    >     the topmost site in the unreconstructed slab. e.g. a point
    >     of [0.5, 0.25, 1] corresponds to the 0.5 frac_coord of a,
    >     0.25 frac_coord of b and a distance of 1 atomic layer above
    >     the topmost site. [0.5, 0.25, -0.5] where the third index
    >     corresponds to a point half a atomic layer below the topmost
    >     site. [0.5, 0.25, 0] corresponds to a point in the same
    >     position along c as the topmost site. This is done because
    >     while the primitive units of a and b will remain constant,
    >     the user can vary the length of the c direction by changing
    >     the slab layer or the vacuum layer.



    * **NOTE** – THE DICTIONARY SHOULD ONLY CONTAIN “points_to_remove” AND


    * **ReconstructionGenerator** (*"points_to_add" FOR THE TOP SURFACE. THE*) –


    * **WITH** (*WILL MODIFY THE BOTTOM SURFACE ACCORDINGLY TO RETURN A SLAB*) –


    * **SURFACES.** (*EQUIVALENT*) –



#### build_slabs()
Builds the reconstructed slab by:


    1. Obtaining the unreconstructed slab using the specified
    parameters for the SlabGenerator.


    2. Applying the appropriate lattice transformation in the
    a and b lattice vectors.


    3. Remove any specified sites from both surfaces.


    4. Add any specified sites to both surfaces.


* **Returns**

    The reconstructed slab.



* **Return type**

    (Slab)



#### get_unreconstructed_slabs()
Generates the unreconstructed or pristine super slab.


### _class_ Slab(lattice, species, coords, miller_index, oriented_unit_cell, shift, scale_factor, reorient_lattice=True, validate_proximity=False, to_unit_cell=False, reconstruction=None, coords_are_cartesian=False, site_properties=None, energy=None, properties=None)
Bases: `Structure`

Subclass of Structure representing a Slab. Implements additional
attributes pertaining to slabs, but the init method does not
actually implement any algorithm that creates a slab. This is a
DUMMY class who’s init method only holds information about the
slab. Also has additional methods that returns other information
about a slab such as the surface area, normal, and atom adsorption.

Note that all Slabs have the surface normal oriented perpendicular to the a
and b lattice vectors. This means the lattice vectors a and b are in the
surface plane and the c vector is out of the surface plane (though not
necessarily perpendicular to the surface).


#### miller_index()
Miller index of plane parallel to surface.


* **Type**

    tuple



#### scale_factor()
Final computed scale factor that brings the parent cell to the surface cell.


* **Type**

    float



#### shift()
The shift value in Angstrom that indicates how much this slab has been shifted.


* **Type**

    float


Makes a Slab structure, a structure object with additional information
and methods pertaining to slabs.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as a
    pymatgen.core.Lattice or
    simply as any 2D array. Each row should correspond to a lattice
    vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
    lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].


    * **species** (*[**Species**]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates of each species.


    * **miller_index** (*[**h**, **k**, **l**]*) – Miller index of plane parallel to
    surface. Note that this is referenced to the input structure. If
    you need this to be based on the conventional cell,
    you should supply the conventional structure.


    * **oriented_unit_cell** (*Structure*) – The oriented_unit_cell from which
    this Slab is created (by scaling in the c-direction).


    * **shift** (*float*) – The shift in the c-direction applied to get the
    termination.


    * **scale_factor** (*np.ndarray*) – scale_factor Final computed scale factor
    that brings the parent cell to the surface cell.


    * **reorient_lattice** (*bool*) – reorients the lattice parameters such that
    the c direction is along the z axis.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **reconstruction** (*str*) – Type of reconstruction. Defaults to None if
    the slab is not reconstructed.


    * **to_unit_cell** (*bool*) – Translates fractional coordinates into the unit cell. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, e.g., {“magmom”:[5,5,5,5]}. The sequences
    have to be the same length as the atomic species and
    fractional_coords. Defaults to None for no properties.


    * **energy** (*float*) – A value for the energy.


    * **properties** (*dict*) – dictionary containing properties associated
    with the whole slab.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### _sites(_: tuple[PeriodicSite, ..._ )

#### add_adsorbate_atom(indices, specie, distance)
Gets the structure of single atom adsorption.
slab structure from the Slab class(in [0, 0, 1]).


* **Parameters**


    * **indices** (*[**int**]*) – Indices of sites on which to put the absorbate.
    Absorbed atom will be displaced relative to the center of
    these sites.


    * **specie** (*Species/Element/str*) – adsorbed atom species


    * **distance** (*float*) – between centers of the adsorbed atom and the
    given site in Angstroms.



#### as_dict()
MSONable dict.


#### _property_ center_of_mass()
Calculates the center of mass of the slab.


#### copy(site_properties=None, sanitize=False)
Convenience method to get a copy of the structure, with options to add
site properties.


* **Parameters**


    * **site_properties** (*dict*) – Properties to add or override. The
    properties are specified in the same way as the constructor,
    i.e., as a dict of the form {property: [values]}. The
    properties should be in the order of the *original* structure
    if you are performing sanitization.


    * **sanitize** (*bool*) – If True, this method will return a sanitized
    structure. Sanitization performs a few things: (i) The sites are
    sorted by electronegativity, (ii) a LLL lattice reduction is
    carried out to obtain a relatively orthogonalized cell,
    (iii) all fractional coords for sites are mapped into the
    unit cell.



* **Returns**

    A copy of the Structure, with optionally new site_properties and
    optionally sanitized.



#### _property_ dipole()
Calculates the dipole of the Slab in the direction of the surface
normal. Note that the Slab must be oxidation state-decorated for this
to work properly. Otherwise, the Slab will always have a dipole of 0.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    Creates slab from dict.



#### get_orthogonal_c_slab()
This method returns a Slab where the normal (c lattice vector) is
“forced” to be exactly orthogonal to the surface a and b lattice
vectors. **Note that this breaks inherent symmetries in the slab.**
It should be pointed out that orthogonality is not required to get good
surface energies, but it can be useful in cases where the slabs are
subsequently used for postprocessing of some kind, e.g. generating
GBs or interfaces.


#### get_sorted_structure(key=None, reverse=False)
Get a sorted copy of the structure. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species. Note that Slab has to override this
because of the different __init__ args.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### get_surface_sites(tag=False)
Returns the surface sites and their indices in a dictionary. The
oriented unit cell of the slab will determine the coordination number
of a typical site. We use VoronoiNN to determine the
coordination number of bulk sites and slab sites. Due to the
pathological error resulting from some surface sites in the
VoronoiNN, we assume any site that has this error is a surface
site as well. This will work for elemental systems only for now. Useful
for analysis involving broken bonds and for finding adsorption sites.


* **Parameters**

    **tag** (*bool*) – Option to adds site attribute “is_surfsite” (bool)
    to all sites of slab. Defaults to False



* **Returns**

    A dictionary grouping sites on top and bottom of the slab together.

        {“top”: [sites with indices], “bottom”: [sites with indices}




#### get_symmetric_site(point, cartesian=False)
This method uses symmetry operations to find equivalent sites on

    both sides of the slab. Works mainly for slabs with Laue
    symmetry. This is useful for retaining the non-polar and
    symmetric properties of a slab when creating adsorbed
    structures or symmetric reconstructions.

Arg:

    point: Fractional coordinate.


* **Returns**

    Fractional coordinate. A point equivalent to the

        parameter point, but on the other side of the slab




* **Return type**

    point



#### get_tasker2_slabs(tol: float = 0.01, same_species_only=True)
Get a list of slabs that have been Tasker 2 corrected.


* **Parameters**


    * **tol** (*float*) – Tolerance to determine if atoms are within same plane.
    This is a fractional tolerance, not an absolute one.


    * **same_species_only** (*bool*) – If True, only that are of the exact same
    species as the atom at the outermost surface are considered for
    moving. Otherwise, all atoms regardless of species that is
    within tol are considered for moving. Default is True (usually
    the desired behavior).



* **Returns**

    ([Slab]) List of tasker 2 corrected slabs.



#### is_polar(tol_dipole_per_unit_area=0.001)
Checks whether the surface is polar by computing the dipole per unit
area. Note that the Slab must be oxidation state-decorated for this
to work properly. Otherwise, the Slab will always be non-polar.


* **Parameters**

    **tol_dipole_per_unit_area** (*float*) – A tolerance. If the dipole
    magnitude per unit area is less than this value, the Slab is
    considered non-polar. Defaults to 1e-3, which is usually
    pretty good. Normalized dipole per unit area is used as it is
    more reliable than using the total, which tends to be larger for
    slabs with larger surface areas.



#### is_symmetric(symprec: float = 0.1)
Checks if surfaces are symmetric, i.e., contains inversion, mirror on (hkl) plane,

    or screw axis (rotation and translation) about [hkl].


* **Parameters**

    **symprec** (*float*) – Symmetry precision used for SpaceGroup analyzer.



* **Returns**

    Whether surfaces are symmetric.



* **Return type**

    bool



#### _property_ normal()
Calculates the surface normal vector of the slab.


#### _property_ surface_area()
Calculates the surface area of the slab.


#### symmetrically_add_atom(specie, point, coords_are_cartesian=False)
Class method for adding a site at a specified point in a slab.

    Will add the corresponding site on the other side of the
    slab to maintain equivalent surfaces.

Arg:

    specie (str): The specie to add
    point (coords): The coordinate of the site in the slab to add.
    coords_are_cartesian (bool): Is the point in Cartesian coordinates


* **Returns**

    The modified slab



* **Return type**

    (Slab)



#### symmetrically_remove_atoms(indices)
Class method for removing sites corresponding to a list of indices.

    Will remove the corresponding site on the other side of the
    slab to maintain equivalent surfaces.

Arg:

    indices ([indices]): The indices of the sites

        in the slab to remove.


### _class_ SlabGenerator(initial_structure, miller_index, min_slab_size, min_vacuum_size, lll_reduce=False, center_slab=False, in_unit_planes=False, primitive=True, max_normal_search=None, reorient_lattice=True)
Bases: `object`

This class generates different slabs using shift values determined by where
a unique termination can be found along with other criteria such as where a
termination doesn’t break a polyhedral bond. The shift value then indicates
where the slab layer will begin and terminate in the slab-vacuum system.


#### oriented_unit_cell()
A unit cell of the parent structure with the miller
index of plane parallel to surface.


* **Type**

    Structure



#### parent()
Parent structure from which Slab was derived.


* **Type**

    Structure



#### lll_reduce()
Whether or not the slabs will be orthogonalized.


* **Type**

    bool



#### center_slab()
Whether or not the slabs will be centered between the vacuum layer.


* **Type**

    bool



#### slab_scale_factor()
Final computed scale factor that brings the parent cell to the
surface cell.


* **Type**

    float



#### miller_index()
Miller index of plane parallel to surface.


* **Type**

    tuple



#### min_slab_size()
Minimum size in angstroms of layers containing atoms.


* **Type**

    float



#### min_vac_size()
Minimum size in angstroms of layers containing vacuum.


* **Type**

    float


Calculates the slab scale factor and uses it to generate a unit cell
of the initial structure that has been oriented by its miller index.
Also stores the initial information needed later on to generate a slab.


* **Parameters**


    * **initial_structure** (*Structure*) – Initial input structure. Note that to
    ensure that the miller indices correspond to usual
    crystallographic definitions, you should supply a conventional
    unit cell structure.


    * **miller_index** (*[**h**, **k**, **l**]*) – Miller index of plane parallel to
    surface. Note that this is referenced to the input structure. If
    you need this to be based on the conventional cell,
    you should supply the conventional structure.


    * **min_slab_size** (*float*) – In Angstroms or number of hkl planes


    * **min_vacuum_size** (*float*) – In Angstroms or number of hkl planes


    * **lll_reduce** (*bool*) – Whether to perform an LLL reduction on the
    eventual structure.


    * **center_slab** (*bool*) – Whether to center the slab in the cell with
    equal vacuum spacing from the top and bottom.


    * **in_unit_planes** (*bool*) – Whether to set min_slab_size and min_vac_size
    in units of hkl planes (True) or Angstrom (False/default).
    Setting in units of planes is useful for ensuring some slabs
    have a certain n_layer of atoms. e.g. for Cs (100), a 10 Ang
    slab will result in a slab with only 2 layer of atoms, whereas
    Fe (100) will have more layer of atoms. By using units of hkl
    planes instead, we ensure both slabs
    have the same number of atoms. The slab thickness will be in
    min_slab_size/math.ceil(self._proj_height/dhkl)
    multiples of oriented unit cells.


    * **primitive** (*bool*) – Whether to reduce any generated slabs to a
    primitive cell (this does **not** mean the slab is generated
    from a primitive cell, it simply means that after slab
    generation, we attempt to find shorter lattice vectors,
    which lead to less surface area and smaller cells).


    * **max_normal_search** (*int*) – If set to a positive integer, the code will
    conduct a search for a normal lattice vector that is as
    perpendicular to the surface as possible by considering
    multiples linear combinations of lattice vectors up to
    max_normal_search. This has no bearing on surface energies,
    but may be useful as a preliminary step to generating slabs
    for absorption and other sizes. It is typical that this will
    not be the smallest possible cell for simulation. Normality
    is not guaranteed, but the oriented cell will have the c
    vector as normal as possible (within the search range) to the
    surface. A value of up to the max absolute Miller index is
    usually sufficient.


    * **reorient_lattice** (*bool*) – reorients the lattice parameters such that
    the c direction is the third vector of the lattice matrix



#### _calculate_possible_shifts(tol: float = 0.1)

#### _get_c_ranges(bonds)

#### get_slab(shift=0, tol: float = 0.1, energy=None)
This method takes in shift value for the c lattice direction and
generates a slab based on the given shift. You should rarely use this
method. Instead, it is used by other generation algorithms to obtain
all slabs.

Arg:

    shift (float): A shift value in Angstrom that determines how much a

        slab should be shifted.

    tol (float): Tolerance to determine primitive cell.
    energy (float): An energy to assign to the slab.


* **Returns**

    (Slab) A Slab object with a particular shifted oriented unit cell.



#### get_slabs(bonds=None, ftol=0.1, tol=0.1, max_broken_bonds=0, symmetrize=False, repair=False)
This method returns a list of slabs that are generated using the list of
shift values from the method, _calculate_possible_shifts(). Before the
shifts are used to create the slabs however, if the user decides to take
into account whether or not a termination will break any polyhedral
structure (bonds is not None), this method will filter out any shift
values that do so.


* **Parameters**


    * **bonds** (*{**(**specie1**, **specie2*) – max_bond_dist}: bonds are
    specified as a dict of tuples: float of specie1, specie2
    and the max bonding distance. For example, PO4 groups may be
    defined as {(“P”, “O”): 3}.


    * **tol** (*float*) – General tolerance parameter for getting primitive
    cells and matching structures


    * **ftol** (*float*) – Threshold parameter in fcluster in order to check
    if two atoms are lying on the same plane. Default thresh set
    to 0.1 Angstrom in the direction of the surface normal.


    * **max_broken_bonds** (*int*) – Maximum number of allowable broken bonds
    for the slab. Use this to limit # of slabs (some structures
    may have a lot of slabs). Defaults to zero, which means no
    defined bonds must be broken.


    * **symmetrize** (*bool*) – Whether or not to ensure the surfaces of the
    slabs are equivalent.


    * **repair** (*bool*) – Whether to repair terminations with broken bonds
    or just omit them. Set to False as repairing terminations can
    lead to many possible slabs as oppose to just omitting them.



* **Returns**

    ([Slab]) List of all possible terminations of a particular surface.
    Slabs are sorted by the # of bonds broken.



#### move_to_other_side(init_slab, index_of_sites)
This method will Move a set of sites to the
other side of the slab (opposite surface).

Arg:

    init_slab (structure): A structure object representing a slab.
    index_of_sites (list of ints): The list of indices representing

    > the sites we want to move to the other side.


* **Returns**

    (Slab) A Slab object with a particular shifted oriented unit cell.



#### nonstoichiometric_symmetrized_slab(init_slab)
This method checks whether or not the two surfaces of the slab are
equivalent. If the point group of the slab has an inversion symmetry (
ie. belong to one of the Laue groups), then it is assumed that the
surfaces should be equivalent. Otherwise, sites at the bottom of the
slab will be removed until the slab is symmetric. Note the removal of sites
can destroy the stoichiometry of the slab. For non-elemental
structures, the chemical potential will be needed to calculate surface energy.

Arg:

    init_slab (Structure): A single slab structure


* **Returns**

    A symmetrized Slab object.



* **Return type**

    Slab (structure)



#### repair_broken_bonds(slab, bonds)
This method will find undercoordinated atoms due to slab
cleaving specified by the bonds parameter and move them
to the other surface to make sure the bond is kept intact.
In a future release of surface.py, the ghost_sites will be
used to tell us how the repair bonds should look like.

Arg:

    slab (structure): A structure object representing a slab.
    bonds ({(specie1, specie2): max_bond_dist}: bonds are

    > specified as a dict of tuples: float of specie1, specie2
    > and the max bonding distance. For example, PO4 groups may be
    > defined as {(“P”, “O”): 3}.


* **Returns**

    (Slab) A Slab object with a particular shifted oriented unit cell.



### _reduce_vector(vector)

### center_slab(slab)
The goal here is to ensure the center of the slab region

    is centered close to c=0.5. This makes it easier to
    find the surface sites and apply operations like doping.

There are three cases where the slab in not centered:

1. The slab region is completely between two vacuums in the
box but not necessarily centered. We simply shift the
slab by the difference in its center of mass and 0.5
along the c direction.

2. The slab completely spills outside the box from the bottom
and into the top. This makes it incredibly difficult to
locate surface sites. We iterate through all sites that
spill over (z>c) and shift all sites such that this specific
site is now on the other side. Repeat for all sites with z>c.

3. This is a simpler case of scenario 2. Either the top or bottom
slab sites are at c=0 or c=1. Treat as scenario 2.


* **Parameters**

    **slab** (*Slab*) – Slab structure to center



* **Returns**

    Returns a centered slab structure



### generate_all_slabs(structure, max_index, min_slab_size, min_vacuum_size, bonds=None, tol=0.1, ftol=0.1, max_broken_bonds=0, lll_reduce=False, center_slab=False, primitive=True, max_normal_search=None, symmetrize=False, repair=False, include_reconstructions=False, in_unit_planes=False)
A function that finds all different slabs up to a certain miller index.
Slabs oriented under certain Miller indices that are equivalent to other
slabs in other Miller indices are filtered out using symmetry operations
to get rid of any repetitive slabs. For example, under symmetry operations,
CsCl has equivalent slabs in the (0,0,1), (0,1,0), and (1,0,0) direction.


* **Parameters**


    * **structure** (*Structure*) – Initial input structure. Note that to
    ensure that the miller indices correspond to usual
    crystallographic definitions, you should supply a conventional
    unit cell structure.


    * **max_index** (*int*) – The maximum Miller index to go up to.


    * **min_slab_size** (*float*) – In Angstroms


    * **min_vacuum_size** (*float*) – In Angstroms


    * **bonds** (*{**(**specie1**, **specie2*) – max_bond_dist}: bonds are
    specified as a dict of tuples: float of specie1, specie2
    and the max bonding distance. For example, PO4 groups may be
    defined as {(“P”, “O”): 3}.


    * **tol** (*float*) – General tolerance parameter for getting primitive
    cells and matching structures


    * **ftol** (*float*) – Threshold parameter in fcluster in order to check
    if two atoms are lying on the same plane. Default thresh set
    to 0.1 Angstrom in the direction of the surface normal.


    * **max_broken_bonds** (*int*) – Maximum number of allowable broken bonds
    for the slab. Use this to limit # of slabs (some structures
    may have a lot of slabs). Defaults to zero, which means no
    defined bonds must be broken.


    * **lll_reduce** (*bool*) – Whether to perform an LLL reduction on the
    eventual structure.


    * **center_slab** (*bool*) – Whether to center the slab in the cell with
    equal vacuum spacing from the top and bottom.


    * **primitive** (*bool*) – Whether to reduce any generated slabs to a
    primitive cell (this does **not** mean the slab is generated
    from a primitive cell, it simply means that after slab
    generation, we attempt to find shorter lattice vectors,
    which lead to less surface area and smaller cells).


    * **max_normal_search** (*int*) – If set to a positive integer, the code will
    conduct a search for a normal lattice vector that is as
    perpendicular to the surface as possible by considering
    multiples linear combinations of lattice vectors up to
    max_normal_search. This has no bearing on surface energies,
    but may be useful as a preliminary step to generating slabs
    for absorption and other sizes. It is typical that this will
    not be the smallest possible cell for simulation. Normality
    is not guaranteed, but the oriented cell will have the c
    vector as normal as possible (within the search range) to the
    surface. A value of up to the max absolute Miller index is
    usually sufficient.


    * **symmetrize** (*bool*) – Whether or not to ensure the surfaces of the
    slabs are equivalent.


    * **repair** (*bool*) – Whether to repair terminations with broken bonds
    or just omit them


    * **include_reconstructions** (*bool*) – Whether to include reconstructed
    slabs available in the reconstructions_archive.json file. Defaults to False.


    * **in_unit_planes** (*bool*) – Whether to generate slabs in units of the primitive
    cell’s c lattice vector. This is useful for generating slabs with
    a specific number of layers, as the number of layers will be
    independent of the Miller index. Defaults to False.


    * **in_unit_planes** – Whether to set min_slab_size and min_vac_size
    in units of hkl planes (True) or Angstrom (False, the default). Setting in
    units of planes is useful for ensuring some slabs have a certain n_layer of
    atoms. e.g. for Cs (100), a 10 Ang slab will result in a slab with only 2
    layer of atoms, whereas Fe (100) will have more layer of atoms. By using units
    of hkl planes instead, we ensure both slabs have the same number of atoms. The
    slab thickness will be in min_slab_size/math.ceil(self._proj_height/dhkl)
    multiples of oriented unit cells.



### get_d(slab)
Determine the distance of space between
each layer of atoms along c.


### get_slab_regions(slab, blength=3.5)
Function to get the ranges of the slab regions. Useful for discerning where
the slab ends and vacuum begins if the slab is not fully within the cell
:param slab: Structure object modelling the surface
:type slab: Structure
:param blength: The bondlength between atoms. You generally

> want this value to be larger than the actual bondlengths in
> order to find atoms that are part of the slab.


### get_symmetrically_distinct_miller_indices(structure, max_index, return_hkil=False)
Returns all symmetrically distinct indices below a certain max-index for
a given structure. Analysis is based on the symmetry of the reciprocal
lattice of the structure.


* **Parameters**


    * **structure** (*Structure*) – input structure.


    * **max_index** (*int*) – The maximum index. For example, a max_index of 1
    means that (100), (110), and (111) are returned for the cubic
    structure. All other indices are equivalent to one of these.


    * **return_hkil** (*bool*) – If true, return hkil form of Miller
    index for hexagonal systems, otherwise return hkl



### get_symmetrically_equivalent_miller_indices(structure, miller_index, return_hkil=True, system: CrystalSystem | None = None)
Returns all symmetrically equivalent indices for a given structure. Analysis
is based on the symmetry of the reciprocal lattice of the structure.


* **Parameters**


    * **structure** (*Structure*) – Structure to analyze


    * **miller_index** (*tuple*) – Designates the family of Miller indices
    to find. Can be hkl or hkil for hexagonal systems


    * **return_hkil** (*bool*) – If true, return hkil form of Miller
    index for hexagonal systems, otherwise return hkl


    * **system** – If known, specify the crystal system of the structure
    so that it does not need to be re-calculated.



### hkl_transformation(transf, miller_index)
Returns the Miller index from setting
A to B using a transformation matrix
:param transf: The transformation matrix

> that transforms a lattice of A to B


* **Parameters**

    **miller_index** (*[**h**, **k**, **l**]*) – Miller index to transform to setting B.



### is_already_analyzed(miller_index: tuple, miller_list: list, symm_ops: list)
Helper function to check if a given Miller index is
part of the family of indices of any index in a list.


* **Parameters**


    * **miller_index** (*tuple*) – The Miller index to analyze


    * **miller_list** (*list*) – List of Miller indices. If the given
    Miller index belongs in the same family as any of the
    indices in this list, return True, else return False


    * **symm_ops** (*list*) – Symmetry operations of a
    lattice, used to define family of indices



### miller_index_from_sites(lattice, coords, coords_are_cartesian=True, round_dp=4, verbose=True)
Get the Miller index of a plane from a list of site coordinates.

A minimum of 3 sets of coordinates are required. If more than 3 sets of
coordinates are given, the best plane that minimises the distance to all
points will be calculated.


* **Parameters**


    * **lattice** (*list** or **Lattice*) – A 3x3 lattice matrix or Lattice object (for
    example obtained from Structure.lattice).


    * **coords** (*iterable*) – A list or numpy array of coordinates. Can be
    Cartesian or fractional coordinates. If more than three sets of
    coordinates are provided, the best plane that minimises the
    distance to all sites will be calculated.


    * **coords_are_cartesian** (*bool**, **optional*) – Whether the coordinates are
    in Cartesian space. If using fractional coordinates set to False.


    * **round_dp** (*int**, **optional*) – The number of decimal places to round the
    miller index to.


    * **verbose** (*bool**, **optional*) – Whether to print warnings.



* **Returns**

    The Miller index.



* **Return type**

    (tuple)


## pymatgen.core.tensors module

This module provides a base class for tensor-like objects and methods for
basic tensor manipulation. It also provides a class, SquareTensor,
that provides basic methods for creating and manipulating rank 2 tensors.


### _class_ SquareTensor(input_array, vscale=None)
Bases: `Tensor`

Base class for doing useful general operations on second rank tensors
(stress, strain etc.).

Create a SquareTensor object. Note that the constructor uses __new__ rather than
__init__ according to the standard method of subclassing numpy ndarrays. Error
is thrown when the class is initialized with non-square matrix.


* **Parameters**


    * **input_array** (*3x3 array-like*) – the 3x3 array-like
    representing the content of the tensor


    * **vscale** (*6x1 array-like*) – 6x1 array-like scaling the
    voigt-notation vector with the tensor entries



#### _property_ det()
Shorthand for the determinant of the SquareTensor.


#### get_scaled(scale_factor)
Scales the tensor by a certain multiplicative scale factor.


* **Parameters**

    **scale_factor** (*float*) – scalar multiplier to be applied to the
    SquareTensor object



#### _property_ inv()
Shorthand for matrix inverse on SquareTensor.


#### is_rotation(tol: float = 0.001, include_improper=True)
Test to see if tensor is a valid rotation matrix, performs a
test to check whether the inverse is equal to the transpose
and if the determinant is equal to one within the specified
tolerance.


* **Parameters**


    * **tol** (*float*) – tolerance to both tests of whether the
    the determinant is one and the inverse is equal
    to the transpose


    * **include_improper** (*bool*) – whether to include improper
    rotations in the determination of validity



#### polar_decomposition(side='right')
Calculates matrices for polar decomposition.


#### _property_ principal_invariants()
Returns a list of principal invariants for the tensor,
which are the values of the coefficients of the characteristic
polynomial for the matrix.


#### refine_rotation()
Helper method for refining rotation matrix by ensuring
that second and third rows are perpendicular to the first.
Gets new y vector from an orthogonal projection of x onto y
and the new z vector from a cross product of the new x and y.


* **Parameters**

    **rotation** (*tol to test for*) –



* **Returns**

    new rotation matrix



#### _property_ trans()
Shorthand for transpose on SquareTensor.


### _class_ Tensor(input_array, vscale=None, check_rank=None)
Bases: `ndarray`, `MSONable`

Base class for doing useful general operations on Nth order tensors,
without restrictions on the type (stress, elastic, strain, piezo, etc.).

Create a Tensor object. Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **input_array** – (array-like with shape 3^N): array-like representing
    a tensor quantity in standard (i. e. non-voigt) notation


    * **vscale** – (N x M array-like): a matrix corresponding
    to the coefficients of the voigt-notation tensor


    * **check_rank** – (int): If not None, checks that input_array’s rank == check_rank.
    Defaults to None.



#### as_dict(voigt: bool = False)
Serializes the tensor object.


* **Parameters**

    **voigt** (*bool*) – flag for whether to store entries in
    voigt-notation. Defaults to false, as information
    may be lost in conversion.


Returns (dict):

    serialized format tensor object


#### average_over_unit_sphere(quad=None)
Method for averaging the tensor projection over the unit
with option for custom quadrature.


* **Parameters**

    **quad** (*dict*) – quadrature for integration, should be
    dictionary with “points” and “weights” keys defaults
    to quadpy.sphere.Lebedev(19) as read from file



* **Returns**

    Average of tensor projected into vectors on the unit sphere



#### convert_to_ieee(structure: Structure, initial_fit=True, refine_rotation=True)
Given a structure associated with a tensor, attempts a
calculation of the tensor in IEEE format according to
the 1987 IEEE standards.


* **Parameters**


    * **structure** (*Structure*) – a structure associated with the
    tensor to be converted to the IEEE standard


    * **initial_fit** (*bool*) – flag to indicate whether initial
    tensor is fit to the symmetry of the structure.
    Defaults to true. Note that if false, inconsistent
    results may be obtained due to symmetrically
    equivalent, but distinct transformations
    being used in different versions of spglib.


    * **refine_rotation** (*bool*) – whether to refine the rotation
    produced by the ieee transform generator, default True



#### einsum_sequence(other_arrays, einsum_string=None)
Calculates the result of an einstein summation expression.


#### fit_to_structure(structure: Structure, symprec: float = 0.1)
Returns a tensor that is invariant with respect to symmetry
operations corresponding to a structure.


* **Parameters**


    * **structure** (*Structure*) – structure from which to generate
    symmetry operations


    * **symprec** (*float*) – symmetry tolerance for the Spacegroup Analyzer
    used to generate the symmetry operations



#### _classmethod_ from_dict(d)
Instantiate Tensors from dicts (using MSONable API).


* **Returns**

    hydrated tensor object



* **Return type**

    Tensor



#### _classmethod_ from_values_indices(values, indices, populate=False, structure=None, voigt_rank=None, vsym=True, verbose=False)
Creates a tensor from values and indices, with options
for populating the remainder of the tensor.


* **Parameters**


    * **values** (*floats*) – numbers to place at indices


    * **indices** (*array-likes*) – indices to place values at


    * **populate** (*bool*) – whether to populate the tensor


    * **structure** (*Structure*) – structure to base population
    or fit_to_structure on


    * **voigt_rank** (*int*) – full tensor rank to indicate the
    shape of the resulting tensor. This is necessary
    if one provides a set of indices more minimal than
    the shape of the tensor they want, e.g.
    Tensor.from_values_indices((0, 0), 100)


    * **vsym** (*bool*) – whether to voigt symmetrize during the
    optimization procedure


    * **verbose** (*bool*) – whether to populate verbosely



#### _classmethod_ from_voigt(voigt_input)
Constructor based on the voigt notation vector or matrix.


* **Parameters**

    **voigt_input** (*array-like*) – voigt input for a given tensor



#### get_grouped_indices(voigt=False, \*\*kwargs)
Gets index sets for equivalent tensor values.


* **Parameters**


    * **voigt** (*bool*) – whether to get grouped indices
    of voigt or full notation tensor, defaults
    to false


    * **\*\*kwargs** – keyword args for np.isclose. Can take atol
    and rtol for absolute and relative tolerance, e. g.

    ```python
    >>> tensor.group_array_indices(atol=1e-8)
    ```

    or

    ```python
    >>> tensor.group_array_indices(rtol=1e-5)
    ```




* **Returns**

    list of index groups where tensor values are equivalent to
    within tolerances



#### _static_ get_ieee_rotation(structure, refine_rotation=True)
Given a structure associated with a tensor, determines
the rotation matrix for IEEE conversion according to
the 1987 IEEE standards.


* **Parameters**


    * **structure** (*Structure*) – a structure associated with the
    tensor to be converted to the IEEE standard


    * **refine_rotation** (*bool*) – whether to refine the rotation
    using SquareTensor.refine_rotation



#### get_symbol_dict(voigt=True, zero_index=False, \*\*kwargs)
Creates a summary dict for tensor with associated symbol.


* **Parameters**


    * **voigt** (*bool*) – whether to get symbol dict for voigt
    notation tensor, as opposed to full notation,
    defaults to true


    * **zero_index** (*bool*) – whether to set initial index to zero,
    defaults to false, since tensor notations tend to use
    one-indexing, rather than zero indexing like python


    * **\*\*kwargs** – keyword args for np.isclose. Can take atol
    and rtol for absolute and relative tolerance, e. g.

    ```python
    >>> tensor.get_symbol_dict(atol=1e-8)
    ```

    or

    ```python
    >>> tensor.get_symbol_dict(rtol=1e-5)
    ```




* **Returns**

    list of index groups where tensor values are equivalent to
    within tolerances



#### _static_ get_voigt_dict(rank)
Returns a dictionary that maps indices in the tensor to those
in a voigt representation based on input rank.


* **Parameters**

    **rank** (*int*) – Tensor rank to generate the voigt map



#### is_fit_to_structure(structure: Structure, tol: float = 0.01)
Tests whether a tensor is invariant with respect to the
symmetry operations of a particular structure by testing
whether the residual of the symmetric portion is below a
tolerance.


* **Parameters**


    * **structure** (*Structure*) – structure to be fit to


    * **tol** (*float*) – tolerance for symmetry testing



#### is_symmetric(tol: float = 1e-05)
Tests whether a tensor is symmetric or not based on the residual
with its symmetric part, from self.symmetrized.


* **Parameters**

    **tol** (*float*) – tolerance to test for symmetry



#### is_voigt_symmetric(tol: float = 1e-06)
Tests symmetry of tensor to that necessary for voigt-conversion
by grouping indices into pairs and constructing a sequence of
possible permutations to be used in a tensor transpose.


#### populate(structure: Structure, prec: float = 1e-05, maxiter: int = 200, verbose: bool = False, precond: bool = True, vsym: bool = True)
Takes a partially populated tensor, and populates the non-zero
entries according to the following procedure, iterated until
the desired convergence (specified via prec) is achieved.


1. Find non-zero entries


2. Symmetrize the tensor with respect to crystal symmetry and
(optionally) voigt symmetry


3. Reset the non-zero entries of the original tensor


* **Parameters**


    * **structure** (*Structure*) – structure to base population on


    * **prec** (*float*) – precision for determining a non-zero value


    * **maxiter** (*int*) – maximum iterations for populating the tensor


    * **verbose** (*bool*) – whether to populate verbosely


    * **precond** (*bool*) – whether to precondition by cycling through
    all symmops and storing new nonzero values, default True


    * **vsym** (*bool*) – whether to enforce voigt symmetry, defaults
    to True



* **Returns**

    Populated tensor



* **Return type**

    Tensor



#### project(n)
Convenience method for projection of a tensor into a
vector. Returns the tensor dotted into a unit vector
along the input n.


* **Parameters**

    **n** (*3x1 array-like*) – direction to project onto


Returns (float):

    scalar value corresponding to the projection of
    the tensor into the vector


#### rotate(matrix, tol: float = 0.001)
Applies a rotation directly, and tests input matrix to ensure a valid
rotation.


* **Parameters**


    * **matrix** (*3x3 array-like*) – rotation matrix to be applied to tensor


    * **tol** (*float*) – tolerance for testing rotation matrix validity



#### round(decimals=0)
Wrapper around numpy.round to ensure object
of same type is returned.


* **Parameters**

    **decimals** – Number of decimal places to round to (default: 0).
    If decimals is negative, it specifies the number of
    positions to the left of the decimal point.


Returns (Tensor):

    rounded tensor of same type


#### structure_transform(original_structure, new_structure, refine_rotation=True)
Transforms a tensor from one basis for an original structure
into a new basis defined by a new structure.


* **Parameters**


    * **original_structure** (*Structure*) – structure corresponding
    to the basis of the current tensor


    * **new_structure** (*Structure*) – structure corresponding to the
    desired basis


    * **refine_rotation** (*bool*) – whether to refine the rotations
    generated in get_ieee_rotation



* **Returns**

    Tensor that has been transformed such that its basis
    corresponds to the new_structure’s basis



#### symbol(_ = 'T_ )

#### _property_ symmetrized()
Returns a generally symmetrized tensor, calculated by taking
the sum of the tensor and its transpose with respect to all
possible permutations of indices.


#### transform(symm_op)
Applies a transformation (via a symmetry operation) to a tensor.


* **Parameters**

    **symm_op** (*SymmOp*) – a symmetry operation to apply to the tensor



#### _property_ voigt()
Returns the tensor in Voigt notation.


#### _property_ voigt_symmetrized()
Returns a “voigt”-symmetrized tensor, i. e. a voigt-notation
tensor such that it is invariant wrt permutation of indices.


#### zeroed(tol: float = 0.001)
Returns the matrix with all entries below a certain threshold (i.e. tol) set to zero.


### _class_ TensorCollection(tensor_list, base_class=<class 'pymatgen.core.tensors.Tensor'>)
Bases: `Sequence`, `MSONable`

A sequence of tensors that can be used for fitting data
or for having a tensor expansion.


* **Parameters**


    * **tensor_list** – List of tensors.


    * **base_class** – Class to be used.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict(voigt=False)

* **Parameters**

    **voigt** – Whether to use voight form.



* **Returns**

    Dict representation of TensorCollection.



#### convert_to_ieee(structure: Structure, initial_fit=True, refine_rotation=True)
Convert all tensors to IEEE.


* **Parameters**


    * **structure** – Structure


    * **initial_fit** – Whether to perform an initial fit.


    * **refine_rotation** – Whether to refine the rotation.



* **Returns**

    TensorCollection.



#### fit_to_structure(structure: Structure, symprec: float = 0.1)
Fits all tensors to a Structure.


* **Parameters**


    * **structure** – Structure


    * **symprec** – symmetry precision.



* **Returns**

    TensorCollection.



#### _classmethod_ from_dict(d)
Creates TensorCollection from dict.


* **Parameters**

    **d** – dict



* **Returns**

    TensorCollection



#### _classmethod_ from_voigt(voigt_input_list, base_class=<class 'pymatgen.core.tensors.Tensor'>)
Creates TensorCollection from voigt form.


* **Parameters**


    * **voigt_input_list** – List of voigt tensors


    * **base_class** – Class for tensor.



* **Returns**

    TensorCollection.



#### is_fit_to_structure(structure: Structure, tol: float = 0.01)

* **Parameters**


    * **structure** – Structure


    * **tol** – tolerance



* **Returns**

    Whether all tensors are fitted to Structure.



#### is_symmetric(tol: float = 1e-05)

* **Parameters**

    **tol** – tolerance



* **Returns**

    Whether all tensors are symmetric.



#### is_voigt_symmetric(tol: float = 1e-06)

* **Parameters**

    **tol** – tolerance



* **Returns**

    Whether all tensors are voigt symmetric.



#### _property_ ranks()
Ranks for all tensors.


#### rotate(matrix, tol: float = 0.001)
Rotates TensorCollection.


* **Parameters**


    * **matrix** – Rotation matrix.


    * **tol** – tolerance.



* **Returns**

    TensorCollection.



#### round(\*args, \*\*kwargs)
Round all tensors.


* **Parameters**


    * **args** – Passthrough to Tensor.round


    * **kwargs** – Passthrough to Tensor.round



* **Returns**

    TensorCollection.



#### _property_ symmetrized()
TensorCollection where all tensors are symmetrized.


#### transform(symm_op)
Transforms TensorCollection with a symmetry operation.


* **Parameters**

    **symm_op** – SymmetryOperation.



* **Returns**

    TensorCollection.



#### _property_ voigt()
TensorCollection where all tensors are in voight form.


#### _property_ voigt_symmetrized()
TensorCollection where all tensors are voigt symmetrized.


#### zeroed(tol: float = 0.001)

* **Parameters**

    **tol** – Tolerance



* **Returns**

    TensorCollection where small values are set to 0.



### _class_ TensorMapping(tensors: Sequence[Tensor] = (), values: Sequence = (), tol: float = 1e-05)
Bases: `MutableMapping`

Base class for tensor mappings, which function much like
a dictionary, but use numpy routines to determine approximate
equality to keys for getting and setting items.

This is intended primarily for convenience with things like
stress-strain pairs and fitting data manipulation. In general,
it is significantly less robust than a typical hashing
and should be used with care.

Initialize a TensorMapping.


* **Parameters**


    * **tensors** (*Sequence**[**Tensor**]**, **optional*) – Defaults to (,).


    * **values** (*Sequence**, **optional*) – Values to be associated with tensors. Defaults to (,).


    * **tol** (*float**, **optional*) – an absolute tolerance for getting and setting items in the mapping.
    Defaults to 1e-5.



* **Raises**

    **ValueError** – if tensors and values are not the same length



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _get_item_index(item)

#### items()
Items in mapping.


#### values()
Values in mapping.


### get_uvec(vec)
Gets a unit vector parallel to input vector.


### symmetry_reduce(tensors, structure: Structure, tol: float = 1e-08, \*\*kwargs)
Function that converts a list of tensors corresponding to a structure
and returns a dictionary consisting of unique tensor keys with symmop
values corresponding to transformations that will result in derivative
tensors from the original list.


* **Parameters**


    * **tensors** (*list** of **tensors*) – list of Tensor objects to test for
    symmetrically-equivalent duplicates


    * **structure** (*Structure*) – structure from which to get symmetry


    * **tol** (*float*) – tolerance for tensor equivalence


    * **kwargs** – keyword arguments for the SpacegroupAnalyzer



* **Returns**

    dictionary consisting of unique tensors with symmetry operations
    corresponding to those which will reconstruct the remaining
    tensors as values


## pymatgen.core.trajectory module

This module provides classes to define a simulation trajectory, which could come from
either relaxation or molecular dynamics.


### _class_ Trajectory(species: list[str | Element | Species | DummySpecies | Composition], coords: list[list[Vector3D]] | np.ndarray | list[np.ndarray], charge: float | None = None, spin_multiplicity: float | None = None, lattice: Lattice | Matrix3D | list[Lattice] | list[Matrix3D] | np.ndarray | None = None, \*, site_properties: SitePropsType | None = None, frame_properties: list[dict] | None = None, constant_lattice: bool | None = True, time_step: float | None = None, coords_are_displacement: bool = False, base_positions: list[list[Vector3D]] | np.ndarray | None = None)
Bases: `MSONable`

Trajectory of a geometry optimization or molecular dynamics simulation.

Provides basic functions such as slicing trajectory, combining trajectories, and
obtaining displacements.

In below, N denotes the number of sites in the structure, and M denotes the
number of frames in the trajectory.


* **Parameters**


    * **species** – shape (N,). List of species on each site. Can take in flexible
    input, including:
    i.  A sequence of element / species specified either as string

    > symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    > e.g., (3, 56, …) or actual Element or Species objects.


        1. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** – shape (M, N, 3). fractional coordinates of the sites.


    * **charge** – int or float. Charge of the system. This is only used for Molecule-based
    trajectories.


    * **spin_multiplicity** – int or float. Spin multiplicity of the system. This is only
    used for Molecule-based trajectories.


    * **lattice** – shape (3, 3) or (M, 3, 3). Lattice of the structures in the
    trajectory; should be used together with constant_lattice.
    If constant_lattice=True, this should be a single lattice that is
    common for all structures in the trajectory (e.g. in an NVT run).
    If constant_lattice=False, this should be a list of lattices,
    each for one structure in the trajectory (e.g. in an NPT run or a
    relaxation that allows changing the cell size). This is only used for
    Structure-based trajectories.


    * **site_properties** – Properties associated with the sites. This should be a
    list of M dicts for a single dict. If a list of dicts, each provides
    the site properties for a frame. Each value in a dict should be a
    sequence of length N, giving the properties of the N sites.
    For example, for a trajectory with M=2 and N=4, the
    site_properties can be: [{“magmom”:[5,5,5,5]}, {“magmom”:[5,5,5,5]}].
    If a single dict, the site properties in the dict apply to all frames
    in the trajectory. For example, for a trajectory with M=2 and N=4,
    {“magmom”:[2,2,2,2]} means that, through the entire trajectory,
    the magmom are kept constant at 2 for all four atoms.


    * **frame_properties** – Properties associated with the structure (e.g. total
    energy). This should be a sequence of M dicts, with each dict
    providing the properties for a frame. For example, for a trajectory with
    M=2, the frame_properties can be [{‘energy’:1.0}, {‘energy’:2.0}].


    * **constant_lattice** – Whether the lattice changes during the simulation.
    Should be used together with lattice. See usage there. This is only
    used for Structure-based trajectories.


    * **time_step** – Time step of MD simulation in femto-seconds. Should be None
    for a trajectory representing a geometry optimization.


    * **coords_are_displacement** – Whether coords are given in displacements
    (True) or positions (False). Note, if this is True, coords
    of a frame (say i) should be relative to the previous frame (i.e.
    i-1), but not relative to the base_position.


    * **base_positions** – shape (N, 3). The starting positions of all atoms in the
    trajectory. Used to reconstruct positions when converting from
    displacements to positions. Only needs to be specified if
    coords_are_displacement=True. Defaults to the first index of
    coords when coords_are_displacement=False.



#### _check_frame_props(frame_props: list[dict] | None)
Check data shape of site properties.


#### _check_site_props(site_props: SitePropsType | None)
Check data shape of site properties.


* **Parameters**

    **site_props** (*dict** | **list**[**dict**] **| **None*) – Returns immediately if None.



* **Raises**

    **AssertionError** – If the size of the site properties does not match
        the number of sites in the structure.



#### _static_ _combine_frame_props(prop1: list[dict] | None, prop2: list[dict] | None, len1: int, len2: int)
Combine frame properties.


#### _static_ _combine_lattice(lat1: ndarray, lat2: ndarray, len1: int, len2: int)
Helper function to combine trajectory lattice.


#### _static_ _combine_site_props(prop1: SitePropsType | None, prop2: SitePropsType | None, len1: int, len2: int)
Combine site properties.

Either one of prop1 or prop2 can be None, dict, or a list of dict. All
possibilities of combining them are considered.


#### _get_site_props(frames: int | list[int])
Slice site properties.


#### as_dict()
Return the trajectory as a MSONable dict.


#### extend(trajectory: Trajectory)
Append a trajectory to the current one.

The lattice, coords, and all other properties are combined.


* **Parameters**

    **trajectory** – Trajectory to append.



#### _classmethod_ from_file(filename: str | Path, constant_lattice: bool = True, \*\*kwargs)
Create trajectory from XDATCAR or vasprun.xml file.


* **Parameters**


    * **filename** – Path to the file to read from.


    * **constant_lattice** – Whether the lattice changes during the simulation,
    such as in an NPT MD simulation.


    * **\*\*kwargs** – Additional kwargs passed to Trajectory constructor.



* **Returns**

    A trajectory from the file.



#### _classmethod_ from_molecules(molecules: list[Molecule], \*\*kwargs)
Create trajectory from a list of molecules.

Note: Assumes no atoms removed during simulation.


* **Parameters**


    * **molecules** – pymatgen Molecule objects.


    * **\*\*kwargs** – Additional kwargs passed to Trajectory constructor.



* **Returns**

    A trajectory from the structures.



#### _classmethod_ from_structures(structures: list[Structure], constant_lattice: bool = True, \*\*kwargs)
Create trajectory from a list of structures.

Note: Assumes no atoms removed during simulation.


* **Parameters**


    * **structures** – pymatgen Structure objects.


    * **constant_lattice** – Whether the lattice changes during the simulation,
    such as in an NPT MD simulation.


    * **\*\*kwargs** – Additional kwargs passed to Trajectory constructor.



* **Returns**

    A trajectory from the structures.



#### get_molecule(idx: int)
Get molecule at specified index.


* **Parameters**

    **idx** – Index of molecule.



* **Returns**

    A pymatgen Molecule object.



#### get_structure(idx: int)
Get structure at specified index.


* **Parameters**

    **idx** – Index of structure.



* **Returns**

    A pymatgen Structure object.



#### to_displacements()
Converts positions of trajectory into displacements between consecutive frames.

base_positions and coords should both be in fractional coords. Does
not work for absolute coords because the atoms are to be wrapped into the
simulation box.

This is the opposite operation of to_positions().


#### to_positions()
Convert displacements between consecutive frames into positions.

base_positions and coords should both be in fractional coords or
absolute coords.

This is the opposite operation of to_displacements().


#### write_Xdatcar(filename: str | Path = 'XDATCAR', system: str | None = None, significant_figures: int = 6)
Writes to Xdatcar file.

The supported kwargs are the same as those for the
Xdatcar_from_structs.get_string method and are passed through directly.


* **Parameters**


    * **filename** – Name of file to write.  It’s prudent to end the filename with
    ‘XDATCAR’, as most visualization and analysis software require this
    for autodetection.


    * **system** – Description of system (e.g. 2D MoS2).


    * **significant_figures** – Significant figures in the output file.


## pymatgen.core.units module

This module implements a FloatWithUnit, which is a subclass of float. It
also defines supported units for some commonly used units for energy, length,
temperature, time and charge. FloatWithUnit also support conversion to one
another, and additions and subtractions perform automatic conversion if
units are detected. An ArrayWithUnit is also implemented, which is a subclass
of numpy’s ndarray with similar unit features.


### _class_ ArrayWithUnit(input_array, unit, unit_type=None)
Bases: `ndarray`

Subclasses numpy.ndarray to attach a unit type. Typically, you should
use the pre-defined unit type subclasses such as EnergyArray,
LengthArray, etc. instead of using ArrayWithFloatWithUnit directly.

Supports conversion, addition and subtraction of the same unit type. E.g.,
1 m + 20 cm will be automatically converted to 1.2 m (units follow the
leftmost quantity).

```python
>>> a = EnergyArray([1, 2], "Ha")
>>> b = EnergyArray([1, 2], "eV")
>>> c = a + b
>>> print(c)
[ 1.03674933  2.07349865] Ha
>>> c.to("eV")
array([ 28.21138386,  56.42276772]) eV
```

Override __new__.


#### Error()
alias of `UnitError`


#### _property_ as_base_units()
Returns this ArrayWithUnit in base SI units, including derived units.


* **Returns**

    An ArrayWithUnit object in base SI units



#### conversions()
Returns a string showing the available conversions.
Useful tool in interactive mode.


#### _property_ supported_units()
Supported units for specific unit type.


#### to(new_unit)
Conversion to a new_unit.


* **Parameters**

    **new_unit** – New unit type.



* **Returns**

    A ArrayWithFloatWithUnit object in the new units.


Example usage:
>>> e = EnergyArray([1, 1.1], “Ha”)
>>> e.to(“eV”)
array([ 27.21138386,  29.93252225]) eV


#### _property_ unit(_: st_ )
The unit, e.g., “eV”.


#### _property_ unit_type(_: st_ )
The type of unit. Energy, Charge, etc.


### Charge(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='charge'_ )
A float with a charge unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., C, e (electron charge). Must be valid unit or UnitError
    is raised.



### Energy(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='energy'_ )
A float with an energy unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., eV, kJ, etc. Must be valid unit or UnitError is raised.



### _class_ FloatWithUnit(val, unit, unit_type=None)
Bases: `float`

Subclasses float to attach a unit type. Typically, you should use the
pre-defined unit type subclasses such as Energy, Length, etc. instead of
using FloatWithUnit directly.

Supports conversion, addition and subtraction of the same unit type. E.g.,
1 m + 20 cm will be automatically converted to 1.2 m (units follow the
leftmost quantity). Note that FloatWithUnit does not override the eq
method for float, i.e., units are not checked when testing for equality.
The reason is to allow this class to be used transparently wherever floats
are expected.

```python
>>> e = Energy(1.1, "Ha")
>>> a = Energy(1.1, "Ha")
>>> b = Energy(3, "eV")
>>> c = a + b
>>> print(c)
1.2102479761938871 Ha
>>> c.to("eV")
32.932522246000005 eV
```

Initializes a float with unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – A unit. E.g., “C”.


    * **unit_type** (*str*) – A type of unit. E.g., “charge”



#### Error()
alias of `UnitError`


#### _property_ as_base_units()
Returns this FloatWithUnit in base SI units, including derived units.


* **Returns**

    A FloatWithUnit object in base SI units



#### _classmethod_ from_str(s)
Parse string to FloatWithUnit.

Example: Memory.from_str(“1. Mb”)


#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead

Use from_str instead.


#### _property_ supported_units()
Supported units for specific unit type.


#### to(new_unit)
Conversion to a new_unit. Right now, only supports 1 to 1 mapping of
units of each type.


* **Parameters**

    **new_unit** – New unit type.



* **Returns**

    A FloatWithUnit object in the new units.


Example usage:
>>> e = Energy(1.1, “eV”)
>>> e = Energy(1.1, “Ha”)
>>> e.to(“eV”)
29.932522246 eV


#### _property_ unit(_: st_ )
The unit, e.g., “eV”.


#### _property_ unit_type(_: st_ )
The type of unit. Energy, Charge, etc.


### Length(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='length'_ )
A float with a length unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., m, ang, bohr, etc. Must be valid unit or UnitError is
    raised.



### Mass(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='mass'_ )
A float with a mass unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., amu, kg, etc. Must be valid unit or UnitError is
    raised.



### Memory(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='memory'_ )
A float with a memory unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., Kb, Mb, Gb, Tb. Must be valid unit or UnitError
    is raised.



### Temp(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='temperature'_ )
A float with a temperature unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., K. Only K (kelvin) is supported.



### Time(_ = functools.partial(<class 'pymatgen.core.units.FloatWithUnit'>, unit_type='time'_ )
A float with a time unit.


* **Parameters**


    * **val** (*float*) – Value


    * **unit** (*Unit*) – E.g., s, min, h. Must be valid unit or UnitError is
    raised.



### _class_ Unit(unit_def)
Bases: `Mapping`

Represents a unit, e.g., “m” for meters, etc. Supports compound units.
Only integer powers are supported for units.

Constructs a unit.


* **Parameters**

    **unit_def** – A definition for the unit. Either a mapping of unit to
    powers, e.g., {“m”: 2, “s”: -1} represents “m^2 s^-1”,
    or simply as a string “kg m^2 s^-1”. Note that the supported
    format uses “^” as the power operator and all units must be
    space-separated.



#### Error()
alias of `UnitError`


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ as_base_units()
Converts all units to base SI units, including derived units.


* **Returns**

    (base_units_dict, scaling factor). base_units_dict will not
    contain any constants, which are gathered in the scaling factor.



#### get_conversion_factor(new_unit)
Returns a conversion factor between this unit and a new unit.
Compound units are supported, but must have the same powers in each
unit type.


* **Parameters**

    **new_unit** – The new unit.



### _exception_ UnitError()
Bases: `BaseException`

Exception class for unit errors.


### _check_mappings(u)

### _get_si_unit(unit)

### _my_partial(func, \*args, \*\*kwargs)
Partial returns a partial object and therefore we cannot inherit class
methods defined in FloatWithUnit. This function calls partial and patches
the new class before returning.


### kb(_ = 8.617333262e-0_ )
Definitions of supported units. Values below are essentially scaling and
conversion factors. What matters is the relative values, not the absolute.
The SI units must have factor 1.


### obj_with_unit(obj: Any, unit: str)
Returns a FloatWithUnit instance if obj is scalar, a dictionary of
objects with units if obj is a dict, else an instance of
ArrayWithFloatWithUnit.


* **Parameters**


    * **obj** (*Any*) – Object to be given a unit.


    * **unit** (*str*) – Specific units (eV, Ha, m, ang, etc.).



### unitized(unit)
Useful decorator to assign units to the output of a function. You can also
use it to standardize the output units of a function that already returns
a FloatWithUnit or ArrayWithUnit. For sequences, all values in the sequences
are assigned the same unit. It works with Python sequences only. The creation
of numpy arrays loses all unit information. For mapping types, the values
are assigned units.


* **Parameters**

    **unit** – Specific unit (eV, Ha, m, ang, etc.).


Example usage:

```default
@unitized(unit="kg")
def get_mass():
    return 123.45
```

## pymatgen.core.xcfunc module

This module provides.


### _class_ XcFunc(xc=None, x=None, c=None)
Bases: `MSONable`

This object stores information about the XC correlation functional.

Client code usually creates the object by calling the class methods:

>
> * from_name


> * from_type_name

or code-specific methods such as:

>
> * from_abinit_ixc

Ax XcFunc instance is hashable and can therefore be used as key in dictionaries.

The implementation is based on the libxc conventions
and is inspired to the XML specification for atomic PAW datasets documented at:

> [https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html](https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html)

For convenience, part of the pawxml documentation is reported here.

The xc_functional element defines the exchange-correlation functional used for
generating the dataset. It has the two attributes type and name.

The type attribute can be LDA, GGA, MGGA or HYB.
The name attribute designates the exchange-correlation functional
and can be specified in the following ways:

[1] Taking the names from the LibXC library. The correlation and exchange names

    are stripped from their

    ```
    XC_
    ```

     part and combined with a + sign.
    Here is an example for an LDA functional:

    <xc_functional type=”LDA”, name=”LDA_X+LDA_C_PW”/>

    and this is what PBE will look like:

    <xc_functional type=”GGA”, name=”GGA_X_PBE+GGA_C_PBE”/>

[2] Using one of the following pre-defined aliases:

type    name    LibXC equivalent             Reference
LDA     PW      LDA_X+LDA_C_PW               LDA exchange; Perdew, Wang, PRB 45, 13244 (1992)
GGA     PW91    GGA_X_PW91+GGA_C_PW91        Perdew et al PRB 46, 6671 (1992)
GGA     PBE     GGA_X_PBE+GGA_C_PBE          Perdew, Burke, Ernzerhof, PRL 77, 3865 (1996)
GGA     RPBE    GGA_X_RPBE+GGA_C_PBE         Hammer, Hansen, Nørskov, PRB 59, 7413 (1999)
GGA     revPBE  GGA_X_PBE_R+GGA_C_PBE        Zhang, Yang, PRL 80, 890 (1998)
GGA     PBEsol  GGA_X_PBE_SOL+GGA_C_PBE_SOL  Perdew et al, PRL 100, 136406 (2008)
GGA     AM05    GGA_X_AM05+GGA_C_AM05        Armiento, Mattsson, PRB 72, 085108 (2005)
GGA     BLYP    GGA_X_B88+GGA_C_LYP          Becke, PRA 38, 3098 (1988); Lee, Yang, Parr, PRB 37, 785


* **Parameters**


    * **xc** – LibxcFunc for XC functional.


    * **x** – LibxcFunc for exchange part. Mutually exclusive with xc.


    * **c** – LibxcFunc for correlation part. Mutually exclusive with xc.



#### abinitixc_to_libxc(_ = {1: {'xc': LibxcFunc.LDA_XC_TETER93}, 2: {'c': LibxcFunc.LDA_C_PZ, 'x': LibxcFunc.LDA_X}, 4: {'c': LibxcFunc.LDA_C_WIGNER, 'x': LibxcFunc.LDA_X}, 5: {'c': LibxcFunc.LDA_C_HL, 'x': LibxcFunc.LDA_X}, 7: {'c': LibxcFunc.LDA_C_PW, 'x': LibxcFunc.LDA_X}, 11: {'c': LibxcFunc.GGA_C_PBE, 'x': LibxcFunc.GGA_X_PBE}, 14: {'c': LibxcFunc.GGA_C_PBE, 'x': LibxcFunc.GGA_X_PBE_R}, 15: {'c': LibxcFunc.GGA_C_PBE, 'x': LibxcFunc.GGA_X_RPBE}_ )

#### _classmethod_ aliases()
List of registered names.


#### as_dict()
Makes XcFunc obey the general json interface used in pymatgen for easier serialization.


#### _classmethod_ asxc(obj)
Convert object into Xcfunc.


#### defined_aliases(_ = {(LibxcFunc.GGA_X_AM05, LibxcFunc.GGA_C_AM05): ('GGA', 'AM05'), (LibxcFunc.GGA_X_B88, LibxcFunc.GGA_C_LYP): ('GGA', 'BLYP'), (LibxcFunc.GGA_X_PBE, LibxcFunc.GGA_C_PBE): ('GGA', 'PBE'), (LibxcFunc.GGA_X_PBE_R, LibxcFunc.GGA_C_PBE): ('GGA', 'revPBE'), (LibxcFunc.GGA_X_PBE_SOL, LibxcFunc.GGA_C_PBE_SOL): ('GGA', 'PBEsol'), (LibxcFunc.GGA_X_PW91, LibxcFunc.GGA_C_PW91): ('GGA', 'PW91'), (LibxcFunc.GGA_X_RPBE, LibxcFunc.GGA_C_PBE): ('GGA', 'RPBE'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_GL): ('LDA', 'GL'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_HL): ('LDA', 'HL'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_PW): ('LDA', 'PW'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_PW_MOD): ('LDA', 'PW_MOD'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_PZ): ('LDA', 'PZ'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_VWN): ('LDA', 'VWN'), (LibxcFunc.LDA_X, LibxcFunc.LDA_C_WIGNER): ('LDA', 'W')_ )

#### _classmethod_ from_abinit_ixc(ixc)
Build the object from Abinit ixc (integer).


#### _classmethod_ from_dict(d)
Makes XcFunc obey the general json interface used in pymatgen for easier serialization.


#### _classmethod_ from_name(name)
Build the object from one of the registered names.


#### _classmethod_ from_type_name(typ, name)
Build the object from (type, name).


#### name()
The name of the functional. If the functional is not found in the aliases,
the string has the form X_NAME+C_NAME.


#### type()
The type of the functional.