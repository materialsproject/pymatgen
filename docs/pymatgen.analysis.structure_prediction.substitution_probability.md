---
layout: default
title: pymatgen.analysis.structure_prediction.substitution_probability.md
nav_exclude: true
---

# pymatgen.analysis.structure_prediction.substitution_probability module

This module provides classes for representing species substitution
probabilities.


### _class_ pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor(lambda_table=None, alpha=-5, threshold=0.001)
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



### _class_ pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionProbability(\*args, \*\*kwargs)
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