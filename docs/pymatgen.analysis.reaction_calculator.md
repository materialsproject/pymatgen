---
layout: default
title: pymatgen.analysis.reaction_calculator.md
nav_exclude: true
---

# pymatgen.analysis.reaction_calculator module

This module provides classes that define a chemical reaction.


### _class_ pymatgen.analysis.reaction_calculator.BalancedReaction(reactants_coeffs, products_coeffs)
Bases: `MSONable`

An object representing a complete chemical reaction.

Reactants and products to be specified as dict of {Composition: coeff}.


* **Parameters**


    * **(****{Composition** (*products_coeffs*) – float}): Reactants as dict of
    {Composition: amt}.


    * **(****{Composition** – float}): Products as dict of
    {Composition: amt}.



#### TOLERANCE(_ = 1e-0_ )

#### _property_ all_comp()
List of all compositions in the reaction.


#### as_dict()

* **Returns**

    A dictionary representation of BalancedReaction.



#### as_entry(energies)
Returns a ComputedEntry representation of the reaction.
:return:


#### calculate_energy(energies)
Calculates the energy of the reaction.


* **Parameters**

    **(****{Composition** (*energies*) – float}): Energy for each composition.
    E.g ., {comp1: energy1, comp2: energy2}.



* **Returns**

    reaction energy as a float.



#### _property_ coeffs()
Final coefficients of the calculated reaction.


#### _property_ elements()
List of elements in the reaction.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – from as_dict().



* **Returns**

    A BalancedReaction object.



#### _static_ from_str(rxn_string)
Generates a balanced reaction from a string. The reaction must
already be balanced.


* **Parameters**

    **rxn_string** – The reaction string. For example, “4 Li + O2-> 2Li2O”



* **Returns**

    BalancedReaction



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_coeff(comp)
Returns coefficient for a particular composition.


#### get_el_amount(element)
Returns the amount of the element in the reaction.


* **Parameters**

    **element** (*Element/Species*) – Element in the reaction



* **Returns**

    Amount of that element in the reaction.



#### normalize_to(comp, factor=1)
Normalizes the reaction to one of the compositions.
By default, normalizes such that the composition given has a
coefficient of 1. Another factor can be specified.


* **Parameters**


    * **comp** ([*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)) – Composition to normalize to


    * **factor** (*float*) – Factor to normalize to. Defaults to 1.



#### normalize_to_element(element, factor=1)
Normalizes the reaction to one of the elements.
By default, normalizes such that the amount of the element is 1.
Another factor can be specified.


* **Parameters**


    * **element** (*Element/Species*) – Element to normalize to.


    * **factor** (*float*) – Factor to normalize to. Defaults to 1.



#### _property_ normalized_repr()
A normalized representation of the reaction. All factors are converted
to lowest common factors.


#### normalized_repr_and_factor()
Normalized representation for a reaction
For example, `4 Li + 2 O -> 2Li2O` becomes `2 Li + O -> Li2O`.


#### _property_ products()
List of products.


#### _property_ reactants()
List of reactants.


### _class_ pymatgen.analysis.reaction_calculator.ComputedReaction(reactant_entries, product_entries)
Bases: `Reaction`

Convenience class to generate a reaction from ComputedEntry objects, with
some additional attributes, such as a reaction energy based on computed
energies.


* **Parameters**


    * **reactant_entries** (*[*[*ComputedEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)*]*) – List of reactant_entries.


    * **product_entries** (*[*[*ComputedEntry*](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)*]*) – List of product_entries.



#### _property_ all_entries()
Equivalent of all_comp but returns entries, in the same order as the
coefficients.


#### as_dict()

* **Returns**

    A dictionary representation of ComputedReaction.



#### _property_ calculated_reaction_energy()
Returns (float):
The calculated reaction energy.


#### _property_ calculated_reaction_energy_uncertainty()
Calculates the uncertainty in the reaction energy based on the uncertainty in the
energies of the products and reactants.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – from as_dict().



* **Returns**

    A ComputedReaction object.



### _class_ pymatgen.analysis.reaction_calculator.Reaction(reactants, products)
Bases: `BalancedReaction`

A more flexible class representing a Reaction. The reaction amounts will
be automatically balanced. Reactants and products can swap sides so that
all coefficients are positive, however this class will find the solution
with the minimum number of swaps and coefficients of 0. Normalizes so that
the *FIRST* product (or products, if underdetermined) has a coefficient of one.

Reactants and products to be specified as list of
pymatgen.core.structure.Composition. e.g., [comp1, comp2].


* **Parameters**


    * **reactants** (*[*[*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)*]*) – List of reactants.


    * **products** (*[*[*Composition*](pymatgen.core.composition.md#pymatgen.core.composition.Composition)*]*) – List of products.



#### as_dict()

* **Returns**

    A dictionary representation of Reaction.



#### copy()
Returns a copy of the Reaction object.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – from as_dict().



* **Returns**

    A Reaction object.



### _exception_ pymatgen.analysis.reaction_calculator.ReactionError(msg)
Bases: `Exception`

Exception class for Reactions. Allows more information in exception
messages to cover situations not covered by standard exception classes.

Create a ReactionError.


* **Parameters**

    **msg** (*str*) – More information about the ReactionError.