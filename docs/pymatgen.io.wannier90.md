---
layout: default
title: pymatgen.io.wannier90.md
nav_exclude: true
---

# pymatgen.io.wannier90 module

Modules for working with wannier90 input and output.


### _class_ pymatgen.io.wannier90.Unk(ik: int, data: ndarray)
Bases: `object`

Object representing the data in a UNK file.


#### ik()
int index of kpoint for this file


#### data()
numpy.ndarray that contains the wavefunction data for in the UNK file.
The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and
(nbnd, 2, ngx, ngy, ngz) for noncollinear calculations.


#### is_noncollinear()
bool that specifies if data is from a noncollinear calculation


#### nbnd()
int number of bands in data


#### ng()
sequence of three integers that correspond to the grid size of the
given data. The definition is ng = (ngx, ngy, ngz).

Initialize Unk class.


* **Parameters**


    * **ik** (*int*) – index of the kpoint UNK file is for


    * **data** (*np.ndarray*) – data from the UNK file that has shape (nbnd,
    ngx, ngy, ngz) or (nbnd, 2, ngx, ngy, ngz) if noncollinear



#### _property_ data(_: ndarra_ )
contains the wavefunction data for in the UNK file.
The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and
(nbnd, 2, ngx, ngy, ngz) for noncollinear calculations.


* **Type**

    np.ndarray



#### _static_ from_file(filename: str)
Reads the UNK data from file.


* **Parameters**

    **filename** (*str*) – path to UNK file to read



* **Returns**

    Unk object



#### ik(_: in_ )

#### is_noncollinear(_: boo_ )

#### nbnd(_: in_ )

#### ng(_: Sequence[int_ )

#### write_file(filename: str)
Write the UNK file.


* **Parameters**

    **filename** (*str*) – path to UNK file to write, the name should have the
    form ‘UNKXXXXX.YY’ where XXXXX is the kpoint index (Unk.ik) and
    YY is 1 or 2 for the spin index or NC if noncollinear