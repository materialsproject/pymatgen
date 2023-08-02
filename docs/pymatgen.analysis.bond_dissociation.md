---
layout: default
title: pymatgen.analysis.bond_dissociation.md
nav_exclude: true
---

# pymatgen.analysis.bond_dissociation module

Module for BondDissociationEnergies.


### _class_ pymatgen.analysis.bond_dissociation.BondDissociationEnergies(molecule_entry: dict[str, str | dict[str, str | int]], fragment_entries: list[dict[str, str | dict[str, str | int]]], allow_additional_charge_separation: bool = False, multibreak: bool = False)
Bases: `MSONable`

Standard constructor for bond dissociation energies. All bonds in the principle molecule are
looped through and their dissociation energies are calculated given the energies of the resulting
fragments, or, in the case of a ring bond, from the energy of the molecule obtained from breaking
the bond and opening the ring. This class should only be called after the energies of the optimized
principle molecule and all relevant optimized fragments have been determined, either from quantum
chemistry or elsewhere. It was written to provide the analysis after running an Atomate fragmentation
workflow.

Note that the entries passed by the user must have the following keys: formula_pretty, initial_molecule,
final_molecule. If a PCM is present, all entries should also have a pcm_dielectric key.


* **Parameters**


    * **molecule_entry** (*dict*) – Entry for the principle molecule. Should have the keys mentioned above.


    * **fragment_entries** (*list** of **dicts*) – List of fragment entries. Each should have the keys mentioned above.


    * **allow_additional_charge_separation** (*bool*) – If True, consider larger than normal charge separation
    among fragments. Defaults to False. See the definition of self.expected_charges below for more
    specific information.


    * **multibreak** (*bool*) – If True, additionally attempt to break pairs of bonds. Defaults to False.



#### build_new_entry(frags, bonds)
Simple function to format a bond dissociation entry that will eventually be returned to the user.


* **Parameters**


    * **frags** –


    * **bonds** –



* **Returns**




#### filter_fragment_entries(fragment_entries)
Filter the fragment entries.


* **Parameters**

    **fragment_entries** –



* **Returns**




#### fragment_and_process(bonds)
Fragment and process bonds.


* **Parameters**

    **bonds** – Bonds to process.



* **Returns**




#### search_fragment_entries(frag)
Search all fragment entries for those isomorphic to the given fragment.
We distinguish between entries where both initial and final molgraphs are isomorphic to the
given fragment (entries) vs those where only the initial molgraph is isomorphic to the given
fragment (initial_entries) vs those where only the final molgraph is isomorphic (final_entries).


* **Parameters**

    **frag** – Fragment