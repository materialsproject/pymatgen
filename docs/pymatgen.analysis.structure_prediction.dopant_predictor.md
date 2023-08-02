---
layout: default
title: pymatgen.analysis.structure_prediction.dopant_predictor.md
nav_exclude: true
---

# pymatgen.analysis.structure_prediction.dopant_predictor module

Predicting potential dopants.


### pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_shannon_radii(bonded_structure, num_dopants=5, match_oxi_sign=False)
Get dopant suggestions based on Shannon radii differences.


* **Parameters**


    * **bonded_structure** ([*StructureGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)) – A pymatgen structure graph
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



### pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_substitution_probabilities(structure, num_dopants=5, threshold=0.001, match_oxi_sign=False)
Get dopant suggestions based on substitution probabilities.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – A pymatgen structure decorated with
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