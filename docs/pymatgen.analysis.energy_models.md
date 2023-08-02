---
layout: default
title: pymatgen.analysis.energy_models.md
nav_exclude: true
---

# pymatgen.analysis.energy_models module

This module implements a EnergyModel abstract class and some basic
implementations. Basically, an EnergyModel is any model that returns an
“energy” for any given structure.


### _class_ pymatgen.analysis.energy_models.EnergyModel()
Bases: `MSONable`

Abstract structure filter class.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    EnergyModel



#### _abstract_ get_energy(structure)

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ pymatgen.analysis.energy_models.EwaldElectrostaticModel(real_space_cut=None, recip_space_cut=None, eta=None, acc_factor=8.0)
Bases: `EnergyModel`

Wrapper around EwaldSum to calculate the electrostatic energy.

Initializes the model. Args have the same definitions as in
[`pymatgen.analysis.ewald.EwaldSummation`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation).


* **Parameters**


    * **real_space_cut** (*float*) – Real space cutoff radius dictating how
    many terms are used in the real space sum. Defaults to None,
    which means determine automagically using the formula given
    in gulp 3.1 documentation.


    * **recip_space_cut** (*float*) – Reciprocal space cutoff radius.
    Defaults to None, which means determine automagically using
    the formula given in gulp 3.1 documentation.


    * **eta** (*float*) – Screening parameter. Defaults to None, which means
    determine automatically.


    * **acc_factor** (*float*) – No. of significant figures each sum is
    converged to.



#### as_dict()

* **Returns**

    MSONable dict



#### get_energy(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ pymatgen.analysis.energy_models.IsingModel(j, max_radius)
Bases: `EnergyModel`

A very simple Ising model, with r^2 decay.


* **Parameters**


    * **j** (*float*) – The interaction parameter. E = J \* spin1 \* spin2.


    * **radius** (*float*) – max_radius for the interaction.



#### as_dict()

* **Returns**

    MSONable dict



#### get_energy(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ pymatgen.analysis.energy_models.NsitesModel()
Bases: `EnergyModel`

Sets the energy to the number of sites. More sites => higher “energy”.
Used to rank structures from smallest number of sites to largest number
of sites after enumeration.


#### as_dict()

* **Returns**

    MSONable dict



#### get_energy(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ pymatgen.analysis.energy_models.SymmetryModel(symprec: float = 0.1, angle_tolerance=5)
Bases: `EnergyModel`

Sets the energy to the -ve of the spacegroup number. Higher symmetry =>
lower “energy”.

Args have same meaning as in
`pymatgen.symmetry.finder.SpacegroupAnalyzer`.


* **Parameters**


    * **symprec** (*float*) – Symmetry tolerance. Defaults to 0.1.


    * **angle_tolerance** (*float*) – Tolerance for angles. Defaults to 5 degrees.



#### as_dict()

* **Returns**

    MSONable dict



#### get_energy(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value