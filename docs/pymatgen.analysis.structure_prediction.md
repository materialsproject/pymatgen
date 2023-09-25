---
layout: default
title: pymatgen.analysis.structure_prediction.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.structure_prediction package

Utilities to predict new structures.


## pymatgen.analysis.structure_prediction.dopant_predictor module

Predicting potential dopants.


### _get_dopants(substitutions, num_dopants, match_oxi_sign)
Utility method to get n- and p-type dopants from a list of substitutions.


### _int_to_roman(number)
Utility method to convert an int (less than 20) to a roman numeral.


### _shannon_radii_from_cn(species_list, cn_roman, radius_to_compare=0)
Utility func to get Shannon radii for a particular coordination number.

As the Shannon radii depends on charge state and coordination number,
species without an entry for a particular coordination number will
be skipped.


* **Parameters**


    * **species_list** (*list*) – A list of Species to get the Shannon radii for.


    * **cn_roman** (*str*) – The coordination number as a roman numeral. See
    Species.get_shannon_radius for more details.


    * **radius_to_compare** (*float**, **optional*) – If set, the data will be returned
    with a “radii_diff” key, containing the difference between the
    shannon radii and this radius.



* **Returns**

    The Shannon radii for all Species in species. Formatted
    as a list of dictionaries, with the keys:


    * ”species”: The species with charge state.


    * ”radius”: The Shannon radius for the species.


    * ”radius_diff”: The difference between the Shannon radius and the

        radius_to_compare optional argument.




* **Return type**

    (list of dict)



### get_dopants_from_shannon_radii(bonded_structure, num_dopants=5, match_oxi_sign=False)
Get dopant suggestions based on Shannon radii differences.


* **Parameters**


    * **bonded_structure** ([*StructureGraph*](pymatgen.analysis.md#pymatgen.analysis.graphs.StructureGraph)) – A pymatgen structure graph
    decorated with oxidation states. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **num_dopants** (*int*) – The number of suggestions to return for
    n- and p-type dopants.


    * **match_oxi_sign** (*bool*) – Whether to force the dopant and original species
    to have the same sign of oxidation state. E.g. If the original site
    is in a negative charge state, then only negative dopants will be
    returned.



* **Returns**

    Dopant suggestions, given as a dictionary with keys “n_type” and
    “p_type”. The suggestions for each doping type are given as a list of
    dictionaries, each with they keys:


    * ”radii_diff”: The difference between the Shannon radii of the species.


    * ”dopant_spcies”: The dopant species.


    * ”original_species”: The substituted species.




* **Return type**

    (dict)



### get_dopants_from_substitution_probabilities(structure, num_dopants=5, threshold=0.001, match_oxi_sign=False)
Get dopant suggestions based on substitution probabilities.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – A pymatgen structure decorated with
    oxidation states.


    * **num_dopants** (*int*) – The number of suggestions to return for
    n- and p-type dopants.


    * **threshold** (*float*) – Probability threshold for substitutions.


    * **match_oxi_sign** (*bool*) – Whether to force the dopant and original species
    to have the same sign of oxidation state. E.g. If the original site
    is in a negative charge state, then only negative dopants will be
    returned.



* **Returns**

    Dopant suggestions, given as a dictionary with keys “n_type” and
    “p_type”. The suggestions for each doping type are given as a list of
    dictionaries, each with they keys:


    * ”probability”: The probability of substitution.


    * ”dopant_species”: The dopant species.


    * ”original_species”: The substituted species.




* **Return type**

    (dict)


## pymatgen.analysis.structure_prediction.substitution_probability module

This module provides classes for representing species substitution
probabilities.


### _class_ SubstitutionPredictor(lambda_table=None, alpha=-5, threshold=0.001)
Bases: `object`

Predicts likely substitutions either to or from a given composition
or species list using the SubstitutionProbability.


* **Parameters**


    * **(****)** (*lambda_table*) – Input lambda table.


    * **alpha** (*float*) – weight function for never observed substitutions


    * **threshold** (*float*) – Threshold to use to identify high probability structures.



#### composition_prediction(composition, to_this_composition=True)
Returns charged balanced substitutions from a starting or ending
composition.


* **Parameters**


    * **composition** – starting or ending composition


    * **to_this_composition** – If true, substitutions with this as a final composition
    will be found. If false, substitutions with this as a
    starting composition will be found (these are slightly
    different)



* **Returns**

    List of predictions in the form of dictionaries.
    If to_this_composition is true, the values of the dictionary
    will be from the list species. If false, the keys will be
    from that list.



#### list_prediction(species, to_this_composition=True)

* **Parameters**


    * **species** – list of species


    * **to_this_composition** – If true, substitutions with this as a final composition
    will be found. If false, substitutions with this as a
    starting composition will be found (these are slightly
    different).



* **Returns**

    List of predictions in the form of dictionaries.
    If to_this_composition is true, the values of the dictionary
    will be from the list species. If false, the keys will be
    from that list.



### _class_ SubstitutionProbability(\*args, \*\*kwargs)
Bases: `SubstitutionProbability`

This class finds substitution probabilities given lists of atoms
to substitute. The inputs make more sense if you look through the
from_defaults static method.

The substitution prediction algorithm is presented in:
Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011)
Data Mined Ionic Substitutions for the Discovery of New Compounds.
Inorganic Chemistry, 50(2), 656-663. doi:10.1021/ic102031h


* **Parameters**


    * **lambda_table** – json table of the weight functions lambda if None,
    will use the default lambda.json table


    * **alpha** – weight function for never observed substitutions.


## pymatgen.analysis.structure_prediction.substitutor module

This module provides classes for predicting new structures from existing ones.


### _class_ Substitutor(threshold=0.001, symprec: float = 0.1, \*\*kwargs)
Bases: `MSONable`

This object uses a data mined ionic substitution approach to propose
compounds likely to be stable. It relies on an algorithm presented in
Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011).
Data Mined Ionic Substitutions for the Discovery of New Compounds.
Inorganic Chemistry, 50(2), 656-663. doi:10.1021/ic102031h.

This substitutor uses the substitution probability class to
find good substitutions for a given chemistry or structure.


* **Parameters**


    * **threshold** – probability threshold for predictions


    * **symprec** – symmetry precision to determine if two structures
    are duplicates


    * **kwargs** – kwargs for the SubstitutionProbability object
    lambda_table, alpha



#### _static_ _is_charge_balanced(struct)
Checks if the structure object is charge balanced.


#### _static_ _is_from_chemical_system(chemical_system, struct)
Checks if the structure object is from the given chemical system.


#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    Class



#### get_allowed_species()
Returns the species in the domain of the probability function
any other specie will not work.


#### pred_from_comp(composition)
Similar to pred_from_list except this method returns a list after
checking that compositions are charge balanced.


#### pred_from_list(species_list)
There are an exceptionally large number of substitutions to
look at (260^n), where n is the number of species in the
list. We need a more efficient than brute force way of going
through these possibilities. The brute force method would be:

```default
output = []
for p in itertools.product(self._sp.species_list
                           , repeat = len(species_list)):
    if self._sp.conditional_probability_list(p, species_list)
                           > self._threshold:
        output.append(dict(zip(species_list,p)))
return output
```

Instead of that we do a branch and bound.


* **Parameters**

    **species_list** – list of species in the starting structure



* **Returns**

    list of dictionaries, each including a substitutions
    dictionary, and a probability value



#### pred_from_structures(target_species, structures_list, remove_duplicates=True, remove_existing=False)
Performs a structure prediction targeting compounds containing all of
the target_species, based on a list of structure (those structures
can for instance come from a database like the ICSD). It will return
all the structures formed by ionic substitutions with a probability
higher than the threshold.

### Notes

If the default probability model is used, input structures must
be oxidation state decorated. See AutoOxiStateDecorationTransformation

This method does not change the number of species in a structure. i.e
if the number of target species is 3, only input structures containing
3 species will be considered.


* **Parameters**


    * **target_species** – a list of species with oxidation states
    e.g., [Species(‘Li+’), Species(‘Ni2+’), Species(‘O-2’)]


    * **structures_list** – list of dictionary of the form {‘structure’: Structure object, ‘id’: some id where it comes from}
    The id can for instance refer to an ICSD id.


    * **remove_duplicates** – if True, the duplicates in the predicted structures will
    be removed


    * **remove_existing** – if True, the predicted structures that already exist in the
    structures_list will be removed



* **Returns**

    a list of TransformedStructure objects.


## pymatgen.analysis.structure_prediction.volume_predictor module

Predict volumes of crystal structures.


### _class_ DLSVolumePredictor(cutoff=4.0, min_scaling=0.5, max_scaling=1.5)
Bases: `object`

Data-mined lattice scaling (DLS) scheme that relies on data-mined bond
lengths to predict the crystal volume of a given structure.

As of 2/12/19, we suggest this method be used in conjunction with
min_scaling and max_scaling to prevent instances of very large, unphysical
predicted volumes found in a small subset of structures.


* **Parameters**


    * **cutoff** (*float*) – cutoff radius added to site radius for finding
    site pairs. Necessary to increase only if your initial
    structure guess is extremely bad (atoms way too far apart). In
    all other instances, increasing cutoff gives same answer
    but takes more time.


    * **min_scaling** (*float*) – if not None, this will ensure that the new
    volume is at least this fraction of the original (preventing
    too-small volumes)


    * **max_scaling** (*float*) – if not None, this will ensure that the new
    volume is at most this fraction of the original (preventing
    too-large volumes).



#### get_predicted_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), icsd_vol=False)
Given a structure, returns back the structure scaled to predicted
volume.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure w/unknown volume



* **Returns**

    a Structure object with predicted volume



#### predict(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), icsd_vol=False)
Given a structure, returns the predicted volume.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – a crystal structure with an unknown volume.


    * **icsd_vol** (*bool*) – True if the input structure’s volume comes from ICSD.



* **Returns**

    a float value of the predicted volume.



### _class_ RLSVolumePredictor(check_isostructural=True, radii_type='ionic-atomic', use_bv=True)
Bases: `object`

Reference lattice scaling (RLS) scheme that predicts the volume of a
structure based on a known crystal structure.


* **Parameters**


    * **check_isostructural** – Whether to test that the two structures are
    isostructural. This algo works best for isostructural compounds.
    Defaults to True.


    * **radii_type** (*str*) – Types of radii to use. You can specify “ionic”
    (only uses ionic radii), “atomic” (only uses atomic radii) or
    “ionic-atomic” (uses either ionic or atomic radii, with a
    preference for ionic where possible).


    * **use_bv** (*bool*) – Whether to use BVAnalyzer to determine oxidation
    states if not present.



#### get_predicted_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), ref_structure)
Given a structure, returns back the structure scaled to predicted
volume.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure w/unknown volume


    * **ref_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – A reference structure with a similar
    structure but different species.



* **Returns**

    a Structure object with predicted volume



#### predict(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), ref_structure)
Given a structure, returns the predicted volume.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure w/unknown volume


    * **ref_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – A reference structure with a similar
    structure but different species.



* **Returns**

    a float value of the predicted volume



### _is_ox(structure)