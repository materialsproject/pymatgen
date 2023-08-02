---
layout: default
title: pymatgen.analysis.cost.md
nav_exclude: true
---

# pymatgen.analysis.cost module

This module is used to estimate the cost of various compounds. Costs are taken
from the a CostDB instance, for example a CSV file via CostDBCSV.
For compounds with no cost listed, a Phase Diagram style convex hull
optimization is performed to determine a set of compositions that can be mixed
to give the desired compound with lowest total cost.


### _class_ pymatgen.analysis.cost.CostAnalyzer(costdb)
Bases: `object`

Given a CostDB, figures out the minimum cost solutions via convex hull.


* **Parameters**

    **(****)** (*costdb*) – Cost database.



#### get_cost_per_kg(comp)
Get best estimate of minimum cost/kg based on known data.


* **Parameters**

    **comp** – Composition as a pymatgen.core.structure.Composition



* **Returns**

    float of cost/kg



#### get_cost_per_mol(comp)
Get best estimate of minimum cost/mol based on known data.


* **Parameters**

    **comp** – Composition as a pymatgen.core.structure.Composition



* **Returns**

    float of cost/mol



#### get_lowest_decomposition(composition)
Get the decomposition leading to lowest cost.


* **Parameters**

    **composition** – Composition as a pymatgen.core.structure.Composition



* **Returns**

    amount}



* **Return type**

    Decomposition as a dict of {Entry



### _class_ pymatgen.analysis.cost.CostDB()
Bases: `object`

Abstract class for representing a Cost database.
Can be extended, e.g. for file-based or REST-based databases.


#### _abstract_ get_entries(chemsys)
For a given chemical system, return an array of CostEntries.


* **Parameters**

    **chemsys** – array of Elements defining the chemical system.



* **Returns**

    array of CostEntries



### _class_ pymatgen.analysis.cost.CostDBCSV(filename)
Bases: `CostDB`

Read a CSV file to get costs
Format is formula,cost_per_kg,name,BibTeX.


* **Parameters**

    **filename** (*str*) – Filename of cost database.



#### get_entries(chemsys)
For a given chemical system, return an array of CostEntries.


* **Parameters**

    **chemsys** – array of Elements defining the chemical system.



* **Returns**

    array of CostEntries



### _class_ pymatgen.analysis.cost.CostEntry(composition, cost, name, reference)
Bases: [`PDEntry`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)

Extends PDEntry to include a BibTeX reference and include language about cost.


* **Parameters**


    * **composition** – Composition as a pymatgen.core.structure.Composition


    * **cost** – Cost (per mol, NOT per kg) of the full Composition


    * **name** – Optional parameter to name the entry. Defaults to the reduced
    chemical formula as in PDEntry.


    * **reference** – Reference data as BiBTeX string.