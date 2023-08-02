---
layout: default
title: pymatgen.entries.entry_tools.md
nav_exclude: true
---

# pymatgen.entries.entry_tools module

This module implements functions to perform various useful operations on
entries, such as grouping entries by structure.


### _class_ pymatgen.entries.entry_tools.EntrySet(entries: Iterable[[PDEntry](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry) | [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry) | [ComputedStructureEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)])
Bases: `MutableSet`, `MSONable`

A convenient container for manipulating entries. Allows for generating
subsets, dumping into files, etc.


* **Parameters**

    **entries** – All the entries.



#### add(element)
Add an entry.


* **Parameters**

    **element** – Entry



#### as_dict()
Returns MSONable dict.


#### _property_ chemsys(_: se_ )
Returns:
set representing the chemical system, e.g., {“Li”, “Fe”, “P”, “O”}.


#### discard(element)
Discard an entry.


* **Parameters**

    **element** – Entry



#### _classmethod_ from_csv(filename: str)
Imports PDEntries from a csv.


* **Parameters**

    **filename** – Filename to import from.



* **Returns**

    List of Elements, List of PDEntries



#### get_subset_in_chemsys(chemsys: list[str])
Returns an EntrySet containing only the set of entries belonging to
a particular chemical system (in this definition, it includes all sub
systems). For example, if the entries are from the
Li-Fe-P-O system, and chemsys=[“Li”, “O”], only the Li, O,
and Li-O entries are returned.


* **Parameters**

    **chemsys** – Chemical system specified as list of elements. E.g.,
    [“Li”, “O”]



* **Returns**

    EntrySet



#### _property_ ground_states(_: se_ )
A set containing only the entries that are ground states, i.e., the lowest energy
per atom entry at each composition.


#### is_ground_state(entry)
Boolean indicating whether a given Entry is a ground state.


#### remove_non_ground_states()
Removes all non-ground state entries, i.e., only keep the lowest energy
per atom entry at each composition.


#### to_csv(filename: str, latexify_names: bool = False)
Exports PDEntries to a csv.


* **Parameters**


    * **filename** – Filename to write to.


    * **entries** – PDEntries to export.


    * **latexify_names** – Format entry names to be LaTex compatible,
    e.g., Li_{2}O



### pymatgen.entries.entry_tools.group_entries_by_composition(entries, sort_by_e_per_atom=True)
Given a sequence of Entry-like objects, group them by composition and

    optionally sort by energy above hull.


* **Parameters**


    * **entries** (*List*) – Sequence of Entry-like objects.


    * **sort_by_e_per_atom** (*bool*) – Whether to sort the grouped entries by
    energy per atom (lowest energy first). Default True.



* **Returns**

    Sequence of sequence of entries by composition. e.g,
    [[ entry1, entry2], [entry3, entry4, entry5]]



### pymatgen.entries.entry_tools.group_entries_by_structure(entries, species_to_remove=None, ltol=0.2, stol=0.4, angle_tol=5, primitive_cell=True, scale=True, comparator=None, ncpus=None)
Given a sequence of ComputedStructureEntries, use structure fitter to group
them by structural similarity.


* **Parameters**


    * **entries** – Sequence of ComputedStructureEntries.


    * **species_to_remove** – Sometimes you want to compare a host framework
    (e.g., in Li-ion battery analysis). This allows you to specify
    species to remove before structural comparison.


    * **ltol** (*float*) – Fractional length tolerance. Default is 0.2.


    * **stol** (*float*) – Site tolerance in Angstrom. Default is 0.4 Angstrom.


    * **angle_tol** (*float*) – Angle tolerance in degrees. Default is 5 degrees.


    * **primitive_cell** (*bool*) – If true: input structures will be reduced to
    primitive cells prior to matching. Defaults to True.


    * **scale** – Input structures are scaled to equivalent volume if true;
    For exact matching, set to False.


    * **comparator** – A comparator object implementing an equals method that
    declares equivalency of sites. Default is SpeciesComparator,
    which implies rigid species mapping.


    * **ncpus** – Number of cpus to use. Use of multiple cpus can greatly improve
    fitting speed. Default of None means serial processing.



* **Returns**

    Sequence of sequence of entries by structural similarity. e.g,
    [[ entry1, entry2], [entry3, entry4, entry5]]