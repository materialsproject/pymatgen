---
layout: default
title: pymatgen.cli.pmg_potcar.md
nav_exclude: true
---

# pymatgen.cli.pmg_potcar module

Implementation for pmg potcar CLI.


### pymatgen.cli.pmg_potcar.gen_potcar(dirname, filename)
Generate POTCAR from POTCAR.spec in directories.


* **Parameters**


    * **dirname** (*str*) – Directory name.


    * **filename** (*str*) – Filename in directory.



### pymatgen.cli.pmg_potcar.generate_potcar(args)
Generate POTCAR.


* **Parameters**

    **args** (*dict*) – Args from argparse.



### pymatgen.cli.pmg_potcar.proc_dir(dirname, procfilefunction)
Process a directory.


* **Parameters**


    * **dirname** (*str*) – Directory name.


    * **procfilefunction** (*callable*) – Callable to execute on directory.