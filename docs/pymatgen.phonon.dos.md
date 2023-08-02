---
layout: default
title: pymatgen.phonon.dos.md
nav_exclude: true
---

# pymatgen.phonon.dos module

This module defines classes to represent the phonon density of states, etc.


### _class_ pymatgen.phonon.dos.CompletePhononDos(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), total_dos, pdoss)
Bases: `PhononDos`

This wrapper class defines a total dos, and also provides a list of PDos.


#### pdos()
Dict of partial densities of the form {Site:Densities}


* **Parameters**


    * **structure** – Structure associated with this particular DOS.


    * **total_dos** – total Dos for structure


    * **pdoss** – The pdoss are supplied as an {Site: Densities}.



#### as_dict()
JSON-serializable dict representation of CompletePhononDos.


#### _classmethod_ from_dict(d)
Returns CompleteDos object from dict representation.


#### get_element_dos()
Get element projected Dos.


* **Returns**

    Dos}



* **Return type**

    dict of {Element



#### get_site_dos(site)
Get the Dos for a site.


* **Parameters**

    **site** – Site in Structure associated with CompletePhononDos.



* **Returns**

    PhononDos containing summed orbital densities for site.



### _class_ pymatgen.phonon.dos.PhononDos(frequencies, densities)
Bases: `MSONable`

Basic DOS object. All other DOS objects are extended versions of this
object.


* **Parameters**


    * **frequencies** – A sequences of frequencies in THz


    * **densities** – A list representing the density of states.



#### as_dict()
JSON-serializable dict representation of PhononDos.


#### cv(t, structure=None)
Constant volume specific heat C_v at temperature T obtained from the integration of the DOS.
Only positive frequencies will be used.
Result in J/(K\*mol-c). A mol-c is the abbreviation of a mole-cell, that is, the number
of Avogadro times the atoms in a unit cell. To compare with experimental data the result
should be divided by the number of unit formulas in the cell. If the structure is provided
the division is performed internally and the result is in J/(K\*mol).


* **Parameters**


    * **t** – a temperature in K


    * **structure** – the structure of the system. If not None it will be used to determine the number of
    formula units



* **Returns**

    Constant volume specific heat C_v



#### entropy(t, structure=None)
Vibrational entropy at temperature T obtained from the integration of the DOS.
Only positive frequencies will be used.
Result in J/(K\*mol-c). A mol-c is the abbreviation of a mole-cell, that is, the number
of Avogadro times the atoms in a unit cell. To compare with experimental data the result
should be divided by the number of unit formulas in the cell. If the structure is provided
the division is performed internally and the result is in J/(K\*mol).


* **Parameters**


    * **t** – a temperature in K


    * **structure** – the structure of the system. If not None it will be used to determine the number of
    formula units



* **Returns**

    Vibrational entropy



#### _classmethod_ from_dict(d)
Returns PhononDos object from dict representation of PhononDos.


#### get_interpolated_value(frequency)
Returns interpolated density for a particular frequency.


* **Parameters**

    **frequency** – frequency to return the density for.



#### get_smeared_densities(sigma)
Returns the densities, but with a Gaussian smearing of
std dev sigma applied.


* **Parameters**

    **sigma** – Std dev of Gaussian smearing function.



* **Returns**

    Gaussian-smeared densities.



#### helmholtz_free_energy(t, structure=None)
Phonon contribution to the Helmholtz free energy at temperature T obtained from the integration of the DOS.
Only positive frequencies will be used.
Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
of Avogadro times the atoms in a unit cell. To compare with experimental data the result
should be divided by the number of unit formulas in the cell. If the structure is provided
the division is performed internally and the result is in J/mol.


* **Parameters**


    * **t** – a temperature in K


    * **structure** – the structure of the system. If not None it will be used to determine the number of
    formula units



* **Returns**

    Phonon contribution to the Helmholtz free energy



#### ind_zero_freq()
Index of the first point for which the frequencies are equal or greater than zero.


#### internal_energy(t, structure=None)
Phonon contribution to the internal energy at temperature T obtained from the integration of the DOS.
Only positive frequencies will be used.
Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
of Avogadro times the atoms in a unit cell. To compare with experimental data the result
should be divided by the number of unit formulas in the cell. If the structure is provided
the division is performed internally and the result is in J/mol.


* **Parameters**


    * **t** – a temperature in K


    * **structure** – the structure of the system. If not None it will be used to determine the number of
    formula units



* **Returns**

    Phonon contribution to the internal energy



#### zero_point_energy(structure=None)
Zero point energy of the system. Only positive frequencies will be used.
Result in J/mol-c. A mol-c is the abbreviation of a mole-cell, that is, the number
of Avogadro times the atoms in a unit cell. To compare with experimental data the result
should be divided by the number of unit formulas in the cell. If the structure is provided
the division is performed internally and the result is in J/mol.


* **Parameters**

    **structure** – the structure of the system. If not None it will be used to determine the number of
    formula units



* **Returns**

    Phonon contribution to the internal energy



### pymatgen.phonon.dos.coth(x)
Coth function.


* **Parameters**

    **(****)** (*x*) – value



* **Returns**

    coth(x)