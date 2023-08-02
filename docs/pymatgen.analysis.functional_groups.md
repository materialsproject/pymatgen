---
layout: default
title: pymatgen.analysis.functional_groups.md
nav_exclude: true
---

# pymatgen.analysis.functional_groups module

Determine functional groups present in a Molecule.


### _class_ pymatgen.analysis.functional_groups.FunctionalGroupExtractor(molecule, optimize=False)
Bases: `object`

This class is used to algorithmically parse a molecule (represented by an
instance of pymatgen.analysis.graphs.MoleculeGraph) and determine arbitrary
functional groups.

Instantiation method for FunctionalGroupExtractor.


* **Parameters**


    * **molecule** – Either a filename, a pymatgen.core.structure.Molecule
    object, or a pymatgen.analysis.graphs.MoleculeGraph object.


    * **optimize** – Default False. If True, then the input molecule will be
    modified, adding Hydrogens, performing a simple conformer search,
    etc.



#### categorize_functional_groups(groups)
Determine classes of functional groups present in a set.


* **Parameters**

    **groups** – Set of functional groups.



* **Returns**

    dict containing representations of the groups, the indices of
    where the group occurs in the MoleculeGraph, and how many of each
    type of group there is.



#### get_all_functional_groups(elements=None, func_groups=None, catch_basic=True)
Identify all functional groups (or all within a certain subset) in the
molecule, combining the methods described above.


* **Parameters**


    * **elements** – List of elements that will qualify a carbon as special
    (if only certain functional groups are of interest).
    Default None.


    * **func_groups** – List of strs representing the functional groups of
    interest. Default to None, meaning that all of the functional groups
    defined in this function will be sought.


    * **catch_basic** – bool. If True, use get_basic_functional_groups and
    other methods



* **Returns**

    list of sets of ints, representing groups of connected atoms



#### get_basic_functional_groups(func_groups=None)
Identify functional groups that cannot be identified by the Ertl method
of get_special_carbon and get_heteroatoms, such as benzene rings, methyl
groups, and ethyl groups.

TODO: Think of other functional groups that are important enough to be
added (ex: do we need ethyl, butyl, propyl?)


* **Parameters**

    **func_groups** – List of strs representing the functional groups of
    interest. Default to None, meaning that all of the functional groups
    defined in this function will be sought.



* **Returns**

    list of sets of ints, representing groups of connected atoms



#### get_heteroatoms(elements=None)
Identify non-H, non-C atoms in the MoleculeGraph, returning a list of
their node indices.


* **Parameters**

    **elements** – List of elements to identify (if only certain
    functional groups are of interest).



* **Returns**

    set of ints representing node indices



#### get_special_carbon(elements=None)
Identify Carbon atoms in the MoleculeGraph that fit the characteristics
defined Ertl (2017), returning a list of their node indices.

The conditions for marking carbon atoms are (quoted from Ertl):

    “- atoms connected by non-aromatic double or triple bond to any
    heteroatom
    - atoms in nonaromatic carbon-carbon double or triple bonds
    - acetal carbons, i.e. sp3 carbons connected to two or more oxygens,
    nitrogens or sulfurs; these O, N or S atoms must have only single bonds
    - all atoms in oxirane, aziridine and thiirane rings”


* **Parameters**

    **elements** – List of elements that will qualify a carbon as special
    (if only certain functional groups are of interest).
    Default None.



* **Returns**

    set of ints representing node indices



#### link_marked_atoms(atoms)
Take a list of marked “interesting” atoms (heteroatoms, special carbons)
and attempt to connect them, returning a list of disjoint groups of
special atoms (and their connected hydrogens).


* **Parameters**

    **atoms** – set of marked “interesting” atoms, presumably identified
    using other functions in this class.



* **Returns**

    list of sets of ints, representing groups of connected atoms