---
layout: default
title: pymatgen.core.molecular_orbitals.md
nav_exclude: true
---

# pymatgen.core.molecular_orbitals module

This module implements a MolecularOrbital class to represent band character in
solids. Useful for predicting PDOS character from structural information.


### _class_ pymatgen.core.molecular_orbitals.MolecularOrbitals(formula)
Bases: `object`

Represents the character of bands in a solid. The input is a chemical
formula, since no structural characteristics are taken into account.

The band character of a crystal emerges from the atomic orbitals of the
constituent ions, hybridization/covalent bonds, and the spin-orbit
interaction (ex: Fe2O3). Right now the orbitals are only built from
the uncharged atomic species. Functionality can be improved by:
1) calculate charged ion orbital energies
2) incorporate the coordination environment to account for covalent bonds

The atomic orbital energies are stored in pymatgen.core.periodic_table.JSON

```python
>>> MOs = MolecularOrbitals('SrTiO3')
>>> MOs.band_edges
{'HOMO':['O','2p',-0.338381], 'LUMO':['Ti','3d',-0.17001], 'metal':False}
```


* **Parameters**

    **formula** (*str*) – Chemical formula. Must have integer subscripts. Ex: ‘SrTiO3’.



#### composition()
the composition as a dictionary. Ex: {‘Sr’: 1, ‘Ti’: 1, ‘O’, 3}


#### elements()
the dictionary keys for the composition


#### elec_neg()
the maximum pairwise electronegativity difference


#### aos()
the constituent atomic orbitals for each element as a dictionary


#### band_edges()
dictionary containing the highest occupied molecular orbital (HOMO),
lowest unoccupied molecular orbital (LUMO), and whether the material is predicted
to be a metal


#### aos_as_list()

* **Returns**

    A list of atomic orbitals, sorted from lowest to highest energy.

    The orbitals energies in eV are represented as

        [[‘O’, ‘1s’, -18.758245], [‘O’, ‘2s’, -0.871362], [‘O’, ‘2p’, -0.338381]]

    Data is obtained from
    [https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations](https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations)




#### max_electronegativity()

* **Returns**

    The maximum pairwise electronegativity difference.



#### obtain_band_edges()
Fill up the atomic orbitals with available electrons.


* **Returns**

    HOMO, LUMO, and whether it’s a metal.