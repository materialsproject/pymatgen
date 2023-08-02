---
layout: default
title: pymatgen.analysis.molecule_matcher.md
nav_exclude: true
---

# pymatgen.analysis.molecule_matcher module

This module provides classes to perform fitting of molecule with arbitrary
atom orders.
This module is supposed to perform exact comparisons without the atom order
correspondence prerequisite, while molecule_structure_comparator is supposed
to do rough comparisons with the atom order correspondence prerequisite.

The implementation is based on an excellent python package called rmsd that
you can find at [https://github.com/charnley/rmsd](https://github.com/charnley/rmsd).


### _class_ pymatgen.analysis.molecule_matcher.AbstractMolAtomMapper()
Bases: `MSONable`

Abstract molecular atom order mapping class. A mapping will be able to
find the uniform atom order of two molecules that can pair the
geometrically equivalent atoms.


#### _classmethod_ from_dict(d)

* **Parameters**

    **(****)** (*d*) – Dict.



* **Returns**

    AbstractMolAtomMapper



#### _abstract_ get_molecule_hash(mol)
Defines a hash for molecules. This allows molecules to be grouped
efficiently for comparison.


* **Parameters**

    **mol** – The molecule. OpenBabel OBMol or pymatgen Molecule object



* **Returns**

    A hashable object. Examples can be string formulas, etc.



#### _abstract_ uniform_labels(mol1, mol2)
Pair the geometrically equivalent atoms of the molecules.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol or pymatgen Molecule object.


    * **mol2** – Second molecule. OpenBabel OBMol or pymatgen Molecule object.



* **Returns**

    (list1, list2) if uniform atom order is found. list1 and list2
    are for mol1 and mol2, respectively. Their length equal
    to the number of atoms. They represents the uniform atom order
    of the two molecules. The value of each element is the original
    atom index in mol1 or mol2 of the current atom in uniform atom
    order.
    (None, None) if unform atom is not available.



### _class_ pymatgen.analysis.molecule_matcher.BruteForceOrderMatcher(target: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Bases: `KabschMatcher`

Finding the best match between molecules by selecting molecule order
with the smallest RMSD from all the possible order combinations.

### Notes

When aligning molecules, the atoms of the two molecules **must** have same number
of atoms from the same species.

Constructor of the matcher object.


* **Parameters**

    **target** – a Molecule object used as a target during the alignment



#### fit(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), ignore_warning=False)
Order, rotate and transform p molecule according to the best match.

A ValueError will be raised when the total number of possible combinations
become unfeasible (more than a million combinations).


* **Parameters**


    * **p** – a Molecule object what will be matched with the target one.


    * **ignore_warning** – ignoring error when the number of combination is too large



* **Returns**

    Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    p_prime



#### match(mol: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), ignore_warning: bool = False)
Similar as KabschMatcher.match but this method also finds the order of
atoms which belongs to the best match.

A ValueError will be raised when the total number of possible combinations
become unfeasible (more than a million combination).


* **Parameters**


    * **mol** – a Molecule object what will be matched with the target one.


    * **ignore_warning** – ignoring error when the number of combination is too large



* **Returns**

    The indices of atoms
    U: 3x3 rotation matrix
    V: Translation vector
    rmsd: Root mean squared deviation between P and Q



* **Return type**

    inds



#### _static_ permutations(atoms)
Generates all the possible permutations of atom order. To achieve better
performance all the cases where the atoms are different has been ignored.


### _class_ pymatgen.analysis.molecule_matcher.GeneticOrderMatcher(target: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), threshold: float)
Bases: `KabschMatcher`

This method was inspired by genetic algorithms and tries to match molecules
based on their already matched fragments.

It uses the fact that when two molecule is matching their sub-structures have to match as well.
The main idea here is that in each iteration (generation) we can check the match of all possible
fragments and ignore those which are not feasible.

Although in the worst case this method has N! complexity (same as the brute force one),
in practice it performs much faster because many of the combination can be eliminated
during the fragment matching.

### Notes

This method very robust and returns with all the possible orders.

There is a well known weakness/corner case: The case when there is
a outlier with large deviation with a small index might be ignored.
This happens due to the nature of the average function
used to calculate the RMSD for the fragments.

When aligning molecules, the atoms of the two molecules **must** have the
same number of atoms from the same species.

Constructor of the matcher object.


* **Parameters**


    * **target** – a Molecule object used as a target during the alignment


    * **threshold** – value used to match fragments and prune configuration



#### fit(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Order, rotate and transform all of the matched p molecule
according to the given threshold.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    p_prime: Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    Array of the possible matches where the elements are



#### match(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Similar as KabschMatcher.match but this method also finds all of the
possible atomic orders according to the threshold.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    inds: The indices of atoms
    U: 3x3 rotation matrix
    V: Translation vector
    rmsd: Root mean squared deviation between P and Q



* **Return type**

    Array of the possible matches where the elements are



#### permutations(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Generates all of possible permutations of atom order according the threshold.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Array of index arrays



### _class_ pymatgen.analysis.molecule_matcher.HungarianOrderMatcher(target: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Bases: `KabschMatcher`

This method pre-aligns the molecules based on their principal inertia
axis and then re-orders the input atom list using the Hungarian method.

### Notes

This method cannot guarantee the best match but is very fast.

When aligning molecules, the atoms of the two molecules **must** have same number
of atoms from the same species.

Constructor of the matcher object.


* **Parameters**

    **target** – a Molecule object used as a target during the alignment



#### fit(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Order, rotate and transform p molecule according to the best match.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    p_prime



#### _static_ get_principal_axis(coords, weights)
Get the molecule’s principal axis.


* **Parameters**


    * **coords** – coordinates of atoms


    * **weights** – the weight use for calculating the inertia tensor



* **Returns**

    Array of dim 3 containing the principal axis



#### match(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Similar as KabschMatcher.match but this method also finds the order of
atoms which belongs to the best match.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    The indices of atoms
    U: 3x3 rotation matrix
    V: Translation vector
    rmsd: Root mean squared deviation between P and Q



* **Return type**

    inds



#### _static_ permutations(p_atoms, p_centroid, p_weights, q_atoms, q_centroid, q_weights)
Generates two possible permutations of atom order. This method uses the principle component
of the inertia tensor to prealign the molecules and hungarian method to determine the order.
There are always two possible permutation depending on the way to pre-aligning the molecules.


* **Parameters**


    * **p_atoms** – atom numbers


    * **p_centroid** – array of atom positions


    * **p_weights** – array of atom weights


    * **q_atoms** – atom numbers


    * **q_centroid** – array of atom positions


    * **q_weights** – array of atom weights



* **Yields**

    *perm_inds* – array of atoms’ order



#### _static_ rotation_matrix_vectors(v1, v2)
Returns the rotation matrix that rotates v1 onto v2 using
Rodrigues’ rotation formula.

See more: [https://math.stackexchange.com/a/476311](https://math.stackexchange.com/a/476311)


* **Parameters**


    * **v1** – initial vector


    * **v2** – target vector



* **Returns**

    3x3 rotation matrix



### _class_ pymatgen.analysis.molecule_matcher.InchiMolAtomMapper(angle_tolerance=10.0)
Bases: `AbstractMolAtomMapper`

Pair atoms by inchi labels.


* **Parameters**

    **angle_tolerance** (*float*) – Angle threshold to assume linear molecule. In degrees.



#### as_dict()

* **Returns**

    MSONable dict.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict Representation.



* **Returns**

    InchiMolAtomMapper



#### get_molecule_hash(mol)
Return inchi as molecular hash.


#### uniform_labels(mol1, mol2)

* **Parameters**


    * **mol1** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Molecule 1


    * **mol2** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Molecule 2.



* **Returns**

    Labels



### _class_ pymatgen.analysis.molecule_matcher.IsomorphismMolAtomMapper()
Bases: `AbstractMolAtomMapper`

Pair atoms by isomorphism permutations in the OpenBabel::OBAlign class.


#### as_dict()

* **Returns**

    Jsonable dict.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    IsomorphismMolAtomMapper



#### get_molecule_hash(mol)
Return inchi as molecular hash.


#### uniform_labels(mol1, mol2)
Pair the geometrically equivalent atoms of the molecules.
Calculate RMSD on all possible isomorphism mappings and return mapping
with the least RMSD.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol or pymatgen Molecule object.


    * **mol2** – Second molecule. OpenBabel OBMol or pymatgen Molecule object.



* **Returns**

    (list1, list2) if uniform atom order is found. list1 and list2
    are for mol1 and mol2, respectively. Their length equal
    to the number of atoms. They represents the uniform atom order
    of the two molecules. The value of each element is the original
    atom index in mol1 or mol2 of the current atom in uniform atom
    order.
    (None, None) if unform atom is not available.



### _class_ pymatgen.analysis.molecule_matcher.KabschMatcher(target: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Bases: `MSONable`

Molecule matcher using Kabsch algorithm.

The Kabsch algorithm capable aligning two molecules by finding the parameters
(translation, rotation) which minimize the root-mean-square-deviation (RMSD) of
two molecules which are topologically (atom types, geometry) similar two each other.

### Notes

When aligning molecules, the atoms of the two molecules **must** be in the same
order for the results to be sensible.

Constructor of the matcher object.


* **Parameters**

    **target** – a Molecule object used as a target during the alignment



#### fit(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Rotate and transform p molecule according to the best match.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    p_prime



#### _static_ kabsch(P: ndarray, Q: ndarray)
The Kabsch algorithm is a method for calculating the optimal rotation matrix
that minimizes the root mean squared deviation (RMSD) between two paired sets of points
P and Q, centered around the their centroid.

For more info see:
- [http://en.wikipedia.org/wiki/Kabsch_algorithm](http://en.wikipedia.org/wiki/Kabsch_algorithm) and
- [https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures](https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures)


* **Parameters**


    * **P** – Nx3 matrix, where N is the number of points.


    * **Q** – Nx3 matrix, where N is the number of points.



* **Returns**

    3x3 rotation matrix



* **Return type**

    U



#### match(p: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Using the Kabsch algorithm the alignment of two molecules (P, Q)
happens in three steps:
- translate the P and Q into their centroid
- compute of the optimal rotation matrix (U) using Kabsch algorithm
- compute the translation (V) and rmsd.

The function returns the rotation matrix (U), translation vector (V),
and RMSD between Q and P’, where P’ is:

> P’ = P \* U + V


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Rotation matrix (D,D)
    V: Translation vector (D)
    RMSD : Root mean squared deviation between P and Q



* **Return type**

    U



### _class_ pymatgen.analysis.molecule_matcher.MoleculeMatcher(tolerance: float = 0.01, mapper=None)
Bases: `MSONable`

Class to match molecules and identify whether molecules are the same.


* **Parameters**


    * **tolerance** (*float*) – RMSD difference threshold whether two molecules are
    different


    * **mapper** (*AbstractMolAtomMapper*) – MolAtomMapper object that is able to map the atoms of two
    molecule to uniform order.



#### as_dict()

* **Returns**

    MSONable dict.



#### fit(mol1, mol2)
Fit two molecules.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol or pymatgen Molecule object


    * **mol2** – Second molecule. OpenBabel OBMol or pymatgen Molecule object



* **Returns**

    A boolean value indicates whether two molecules are the same.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    MoleculeMatcher



#### get_rmsd(mol1, mol2)
Get RMSD between two molecule with arbitrary atom order.


* **Returns**

    RMSD if topology of the two molecules are the same
    Infinite if  the topology is different



#### group_molecules(mol_list)
Group molecules by structural equality.


* **Parameters**

    **mol_list** – List of OpenBabel OBMol or pymatgen objects



* **Returns**

    A list of lists of matched molecules
    Assumption: if s1=s2 and s2=s3, then s1=s3
    This may not be true for small tolerances.