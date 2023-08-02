---
layout: default
title: pymatgen.analysis.structure_prediction.substitutor.md
nav_exclude: true
---

# pymatgen.analysis.structure_prediction.substitutor module

This module provides classes for predicting new structures from existing ones.


### _class_ pymatgen.analysis.structure_prediction.substitutor.Substitutor(threshold=0.001, symprec: float = 0.1, \*\*kwargs)
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

Notes:
If the default probability model is used, input structures must
be oxidation state decorated. See AutoOxiStateDecorationTransformation

This method does not change the number of species in a structure. i.e
if the number of target species is 3, only input structures containing
3 species will be considered.


* **Parameters**


    * **target_species** – a list of species with oxidation states
    e.g., [Species(‘Li’,1),Species(‘Ni’,2), Species(‘O’,-2)]


    * **structures_list** – a list of dictionary of the form {‘structure’:Structure object
    ,’id’:some id where it comes from}
    the id can for instance refer to an ICSD id.


    * **remove_duplicates** – if True, the duplicates in the predicted structures will
    be removed


    * **remove_existing** – if True, the predicted structures that already exist in the
    structures_list will be removed



* **Returns**

    a list of TransformedStructure objects.