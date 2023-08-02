---
layout: default
title: pymatgen.cli.pmg_analyze.md
nav_exclude: true
---

# pymatgen.cli.pmg_analyze module

Implementation for pmg analyze CLI.


### pymatgen.cli.pmg_analyze.analyze(args)
Master function controlling which analysis to call.


* **Parameters**

    **args** (*dict*) – args from argparse.



### pymatgen.cli.pmg_analyze.get_energies(rootdir, reanalyze, verbose, quick, sort, fmt)
Get energies of all vaspruns in directory (nested).


* **Parameters**


    * **rootdir** (*str*) – Root directory.


    * **reanalyze** (*bool*) – Whether to ignore saved results and reanalyze


    * **verbose** (*bool*) – Verbose mode or not.


    * **quick** (*bool*) – Whether to perform a quick analysis (using OSZICAR instead
    of vasprun.xml


    * **sort** (*bool*) – Whether to sort the results in ascending order.


    * **fmt** (*str*) – tablefmt passed to tabulate.



### pymatgen.cli.pmg_analyze.get_magnetizations(dir: str, ion_list: list[int])
Get magnetization info from OUTCARs.


* **Parameters**


    * **dir** (*str*) – Directory name


    * **ion_list** (*list**[**int**]*) – List of ions to obtain magnetization information for.



* **Returns**

    0 if successful.



* **Return type**

    int