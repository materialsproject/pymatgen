---
layout: default
title: pymatgen.core.ion.md
nav_exclude: true
---

# pymatgen.core.ion module

Module containing class to create an ion.


### _class_ pymatgen.core.ion.Ion(composition, charge=0.0, _properties=None)
Bases: [`Composition`](pymatgen.core.composition.md#pymatgen.core.composition.Composition), `MSONable`, [`Stringify`](pymatgen.util.string.md#pymatgen.util.string.Stringify)

Ion object. Just a Composition object with an additional variable to store
charge.

The net charge can either be represented as Mn++, Mn+2, Mn[2+], Mn[++], or
Mn[+2]. Note the order of the sign and magnitude in each representation.

Flexible Ion construction, similar to Composition.
For more information, please see pymatgen.core.Composition.


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


#### _property_ composition(_: [Composition](pymatgen.core.composition.md#pymatgen.core.composition.Composition_ )
Composition of ion.


#### _property_ formula(_: st_ )
Returns a formula string, with elements sorted by electronegativity,
e.g., Li4 Fe4 P4 O16.


#### _classmethod_ from_dict(d)
Generates an ion object from a dict created by as_dict().


* **Parameters**

    **d** – {symbol: amount} dict.



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
‘Ca[+2]’.


#### to_pretty_string()

* **Returns**

    Pretty string with proper superscripts.



#### _property_ to_reduced_dict(_: dic_ )
Returns:
dict with element symbol and reduced amount e.g.,
{“Fe”: 2.0, “O”:3.0}.