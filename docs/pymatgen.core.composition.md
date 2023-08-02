---
layout: default
title: pymatgen.core.composition.md
nav_exclude: true
---

# pymatgen.core.composition module

This module implements a Composition class to represent compositions,
and a ChemicalPotential class to represent potentials.


### _class_ pymatgen.core.composition.ChemicalPotential(\*args, \*\*kwargs)
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



### _class_ pymatgen.core.composition.Composition(\*args, strict: bool = False, \*\*kwargs)
Bases: `Hashable`, `Mapping`, `MSONable`, [`Stringify`](pymatgen.util.string.md#pymatgen.util.string.Stringify)

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
Note: Subtly different from get_el_amt_dict in that they keys here are str(Element) instead of Element.symbol.


* **Returns**

    element symbol and (unreduced) amount. E.g.

        {“Fe”: 4.0, “O”:6.0} or {“Fe3+”: 4.0, “O2-“:6.0}




* **Return type**

    dict[str, float]



#### _property_ average_electroneg(_: floa_ )
Average electronegativity of the composition.


* **Type**

    return



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



#### copy()

* **Returns**

    A copy of the composition.



#### _property_ element_composition(_: Compositio_ )
Returns the composition replacing any species by the corresponding
element.


#### _property_ elements(_: list[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies)_ )
Returns view of elements in Composition.


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



#### get_atomic_fraction(el: str | [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies))
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



#### get_wt_fraction(el: str | [Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species) | [DummySpecies](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies))
Calculate weight fraction of an Element or Species.


* **Parameters**

    **el** ([*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)* | *[*Species*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)) – Element or Species to get fraction for.



* **Returns**

    Weight fraction for element el in Composition.



* **Return type**

    float



#### _property_ hill_formula(_: st_ )
Hill formula. The Hill system (or Hill notation) is a system
of writing empirical chemical formulas, molecular chemical formulas and
components of a condensed formula such that the number of carbon atoms
in a molecule is indicated first, the number of hydrogen atoms next,
and then the number of all other chemical elements subsequently, in
alphabetical order of the chemical symbols. When the formula contains
no carbon, all the elements, including hydrogen, are listed
alphabetically.


* **Type**

    return



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

    A list of dicts - each dict reports an element symbol and average

        oxidation state across all sites in that composition. If the
        composition is not charge balanced, an empty list is returned.




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



#### replace(elem_map: dict[str, str | dict[str, int | float]])
Replace elements in a composition. Returns a new Composition, leaving the old one unchanged.


* **Parameters**

    **elem_map** (*dict**[**str**, **str** | **dict**[**str**, **int** | **float**]**]*) – dict of elements or species to swap. E.g.
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

    Same as output __str__() but without spaces.



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


* **Type**

    return



#### _property_ valid(_: boo_ )
Returns True if Composition contains valid elements or species and
False if the Composition contains any dummy species.


#### _property_ weight(_: floa_ )
Total molecular weight of Composition.


### _exception_ pymatgen.core.composition.CompositionError()
Bases: `Exception`

Exception class for composition errors.


### pymatgen.core.composition.reduce_formula(sym_amt, iupac_ordering: bool = False)
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