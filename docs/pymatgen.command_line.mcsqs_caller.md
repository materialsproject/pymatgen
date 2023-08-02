---
layout: default
title: pymatgen.command_line.mcsqs_caller.md
nav_exclude: true
---

# pymatgen.command_line.mcsqs_caller module

Module to call mcsqs, distributed with AT-AT
[https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/).


### _class_ pymatgen.command_line.mcsqs_caller.Sqs(bestsqs, objective_function, allsqs, clusters, directory)
Bases: `tuple`

Return type for run_mcsqs.
bestsqs: Structure
objective_function: float | str
allsqs: List
clusters: List
directory: str


#### allsqs()
Alias for field number 2


#### bestsqs()
Alias for field number 0


#### clusters()
Alias for field number 3


#### directory()
Alias for field number 4


#### objective_function()
Alias for field number 1


### pymatgen.command_line.mcsqs_caller.run_mcsqs(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), clusters: dict[int, float], scaling: int | list[int] = 1, search_time: float = 60, directory: str | None = None, instances: int | None = None, temperature: int | float = 1, wr: float = 1, wn: float = 1, wd: float = 0.5, tol: float = 0.001)
Helper function for calling mcsqs with different arguments
:param structure: Disordered pymatgen Structure object
:type structure: Structure
:param clusters: Dictionary of cluster interactions with entries in the form

> number of atoms: cutoff in angstroms


* **Parameters**


    * **scaling** (*int** or **list*) – Scaling factor to determine supercell. Two options are possible:


            1. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,
        scaling=4 would lead to a 32 atom supercell


            2. A sequence of three scaling factors, e.g., [2, 1, 1], which
        specifies that the supercell should have dimensions 2a x b x c

    Defaults to 1.



    * **search_time** (*float*) – Time spent looking for the ideal SQS in minutes (default: 60)


    * **directory** (*str*) – Directory to run mcsqs calculation and store files (default: None
    runs calculations in a temp directory)


    * **instances** (*int*) – Specifies the number of parallel instances of mcsqs to run
    (default: number of cpu cores detected by Python)


    * **temperature** (*int** or **float*) – Monte Carlo temperature (default: 1), “T” in atat code


    * **wr** (*int** or **float*) – Weight assigned to range of perfect correlation match in objective
    function (default = 1)


    * **wn** (*int** or **float*) – Multiplicative decrease in weight per additional point in cluster (default: 1)


    * **wd** (*int** or **float*) – Exponent of decay in weight as function of cluster diameter (default: 0.5)


    * **tol** (*int** or **float*) – Tolerance for matching correlations (default: 1e-3).



* **Returns**

    Tuple of Pymatgen structure SQS of the input structure, the mcsqs objective function,

        list of all SQS structures, and the directory where calculations are run