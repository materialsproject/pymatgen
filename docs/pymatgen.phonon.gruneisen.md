---
layout: default
title: pymatgen.phonon.gruneisen.md
nav_exclude: true
---

# pymatgen.phonon.gruneisen module

This module provides classes to define a Grueneisen band structure.


### _class_ pymatgen.phonon.gruneisen.GruneisenParameter(qpoints, gruneisen, frequencies, multiplicities=None, structure=None, lattice=None)
Bases: `MSONable`

Class for Grueneisen parameters on a regular grid.


* **Parameters**


    * **qpoints** – list of qpoints as numpy arrays, in frac_coords of the given lattice by default


    * **gruneisen** – list of gruneisen parameters as numpy arrays, shape: (3\*len(structure), len(qpoints))


    * **frequencies** – list of phonon frequencies in THz as a numpy array with shape (3\*len(structure), len(qpoints))


    * **multiplicities** – list of multiplicities


    * **structure** – The crystal structure (as a pymatgen Structure object) associated with the gruneisen parameters.


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object. Pymatgen uses the physics convention of
    reciprocal lattice vectors WITH a 2\*pi coefficient.



#### _property_ acoustic_debye_temp()
Acoustic Debye temperature in K, i.e. the Debye temperature divided by nsites\*\*(1/3).
Adapted from abipy.


#### average_gruneisen(t=None, squared=True, limit_frequencies=None)
Calculates the average of the Gruneisen based on the values on the regular grid.
If squared is True the average will use the squared value of the Gruneisen and a squared root
is performed on the final result.
Values associated to negative frequencies will be ignored.
See Scripta Materialia 129, 88 for definitions.
Adapted from classes in abipy that have been written by Guido Petretto (UCLouvain).


* **Parameters**


    * **t** – the temperature at which the average Gruneisen will be evaluated. If None the acoustic Debye
    temperature is used (see acoustic_debye_temp).


    * **squared** – if True the average is performed on the squared values of the Grueneisen.


    * **limit_frequencies** – if None (default) no limit on the frequencies will be applied.
    Possible values are “debye” (only modes with frequencies lower than the acoustic Debye
    temperature) and “acoustic” (only the acoustic modes, i.e. the first three modes).



* **Returns**

    The average Gruneisen parameter



#### _property_ debye_temp_limit()
Debye temperature in K. Adapted from apipy.


#### debye_temp_phonopy(freq_max_fit=None)
Get Debye temperature in K as implemented in phonopy.


* **Parameters**

    **freq_max_fit** – Maximum frequency to include for fitting.
    Defaults to include first quartile of frequencies.



* **Returns**

    Debye temperature in K.



#### _property_ phdos()
PhononDos object.


* **Type**

    Returns



#### _property_ tdos()
The total DOS (re)constructed from the gruneisen.yaml file.


#### thermal_conductivity_slack(squared=True, limit_frequencies=None, theta_d=None, t=None)
Calculates the thermal conductivity at the acoustic Debye temperature with the Slack formula,
using the average Gruneisen.
Adapted from abipy.


* **Parameters**


    * **squared** (*bool*) – if True the average is performed on the squared values of the Gruenisen


    * **limit_frequencies** – if None (default) no limit on the frequencies will be applied.
    Possible values are “debye” (only modes with frequencies lower than the acoustic Debye
    temperature) and “acoustic” (only the acoustic modes, i.e. the first three modes).


    * **theta_d** – the temperature used to estimate the average of the Gruneisen used in the
    Slack formula. If None the acoustic Debye temperature is used (see
    acoustic_debye_temp). Will also be considered as the Debye temperature in the
    Slack formula.


    * **t** – temperature at which the thermal conductivity is estimated. If None the value at
    the calculated acoustic Debye temperature is given. The value is obtained as a
    simple rescaling of the value at the Debye temperature.



* **Returns**

    The value of the thermal conductivity in W/(m\*K)



### _class_ pymatgen.phonon.gruneisen.GruneisenPhononBandStructure(qpoints, frequencies, gruneisenparameters, lattice, eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure=None)
Bases: [`PhononBandStructure`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure)

This is the most generic phonon band structure data possible
it’s defined by a list of qpoints + frequencies for each of them.
Additional information may be given for frequencies at Gamma, where
non-analytical contribution may be taken into account.


* **Parameters**


    * **qpoints** – list of qpoint as numpy arrays, in frac_coords of the
    given lattice by default


    * **frequencies** – list of phonon frequencies in THz as a numpy array with shape
    (3\*len(structure), len(qpoints)). The First index of the array
    refers to the band and the second to the index of the qpoint.


    * **gruneisenparameters** – list of Grueneisen parameters with the same structure
    frequencies.


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient.


    * **eigendisplacements** – the phonon eigendisplacements associated to the
    frequencies in Cartesian coordinates. A numpy array of complex
    numbers with shape (3\*len(structure), len(qpoints), len(structure), 3).
    The first index of the array refers to the band, the second to the index
    of the qpoint, the third to the atom in the structure and the fourth
    to the Cartesian coordinates.


    * **labels_dict** – (dict) of {} this links a qpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether the qpoint coordinates are Cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.



#### as_dict()

* **Returns**

    MSONable (dict).



#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    Phonon band structure with Grueneisen parameters.



* **Return type**

    GruneisenPhononBandStructure



### _class_ pymatgen.phonon.gruneisen.GruneisenPhononBandStructureSymmLine(qpoints, frequencies, gruneisenparameters, lattice, eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure=None)
Bases: `GruneisenPhononBandStructure`, [`PhononBandStructureSymmLine`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine)

This object stores a GruneisenPhononBandStructureSymmLine together with Grueneisen parameters
for every frequency.


* **Parameters**


    * **qpoints** – list of qpoints as numpy arrays, in frac_coords of the
    given lattice by default


    * **frequencies** – list of phonon frequencies in eV as a numpy array with shape
    (3\*len(structure), len(qpoints))


    * **gruneisenparameters** – list of Grueneisen parameters as a numpy array with the
    shape (3\*len(structure), len(qpoints))


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient


    * **eigendisplacements** – the phonon eigendisplacements associated to the
    frequencies in Cartesian coordinates. A numpy array of complex
    numbers with shape (3\*len(structure), len(qpoints), len(structure), 3).
    The first index of the array refers to the band, the second to the index
    of the qpoint, the third to the atom in the structure and the fourth
    to the Cartesian coordinates.


    * **labels_dict** – (dict) of {} this links a qpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether the qpoint coordinates are cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.



#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** – Dict representation.


Returns: GruneisenPhononBandStructureSymmLine