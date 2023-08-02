---
layout: default
title: pymatgen.io.lammps.utils.md
nav_exclude: true
---

# pymatgen.io.lammps.utils module

This module defines utility classes and functions.


### _class_ pymatgen.io.lammps.utils.LammpsRunner(input_filename='lammps.in', bin='lammps')
Bases: `object`

LAMMPS wrapper.


* **Parameters**


    * **input_filename** (*str*) – input file name


    * **bin** (*str*) – command to run, excluding the input file name.



#### run()
Write the input/data files and run LAMMPS.


### _class_ pymatgen.io.lammps.utils.Polymer(start_monomer, s_head, s_tail, monomer, head, tail, end_monomer, e_head, e_tail, n_units, link_distance=1.0, linear_chain=False)
Bases: `object`

Generate polymer chain via Random walk. At each position there are
a total of 5 possible moves(excluding the previous direction).


* **Parameters**


    * **start_monomer** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Starting molecule


    * **s_head** (*int*) – starting atom index of the start_monomer molecule


    * **s_tail** (*int*) – tail atom index of the start_monomer


    * **monomer** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – The monomer


    * **head** (*int*) – index of the atom in the monomer that forms the head


    * **tail** (*int*) – tail atom index. monomers will be connected from
    tail to head


    * **end_monomer** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Terminal molecule


    * **e_head** (*int*) – starting atom index of the end_monomer molecule


    * **e_tail** (*int*) – tail atom index of the end_monomer


    * **n_units** (*int*) – number of monomer units excluding the start and
    terminal molecules


    * **link_distance** (*float*) – distance between consecutive monomers


    * **linear_chain** (*bool*) – linear or random walk polymer chain.