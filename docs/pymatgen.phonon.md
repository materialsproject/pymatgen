---
layout: default
title: pymatgen.phonon.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.phonon package

Phonon DOS and bandstructure analysis package.


## pymatgen.phonon.bandstructure module

This module provides classes to define a phonon band structure.


### _class_ PhononBandStructure(qpoints: list[[Kpoint](pymatgen.electronic_structure.md#pymatgen.electronic_structure.bandstructure.Kpoint)], frequencies: np.ndarray, lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice), nac_frequencies=None, eigendisplacements=None, nac_eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | None = None)
Bases: `MSONable`

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


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient.


    * **nac_frequencies** – Frequencies with non-analytical contributions at Gamma in THz.
    A list of tuples. The first element of each tuple should be a list
    defining the direction (not necessarily a versor, will be normalized
    internally). The second element containing the 3\*len(structure)
    phonon frequencies with non-analytical correction for that direction.


    * **eigendisplacements** – the phonon eigendisplacements associated to the
    frequencies in Cartesian coordinates. A numpy array of complex
    numbers with shape (3\*len(structure), len(qpoints), len(structure), 3).
    he First index of the array refers to the band, the second to the index
    of the qpoint, the third to the atom in the structure and the fourth
    to the Cartesian coordinates.


    * **nac_eigendisplacements** – the phonon eigendisplacements associated to the
    non-analytical frequencies in nac_frequencies in Cartesian coordinates.
    A list of tuples. The first element of each tuple should be a list
    defining the direction. The second element containing a numpy array of
    complex numbers with shape (3\*len(structure), len(structure), 3).


    * **labels_dict** – (dict) of {} this links a qpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether the qpoint coordinates are Cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.



#### as_dict()
MSONable dict.


#### asr_breaking(tol_eigendisplacements=1e-05)
Returns the breaking of the acoustic sum rule for the three acoustic modes,
if Gamma is present. None otherwise.
If eigendisplacements are available they are used to determine the acoustic
modes: selects the bands corresponding  to the eigendisplacements that
represent to a translation within tol_eigendisplacements. If these are not
identified or eigendisplacements are missing the first 3 modes will be used
(indices [0:3]).


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** – Dict representation



* **Returns**

    PhononBandStructure



#### get_nac_eigendisplacements_along_dir(direction)
Returns the nac_eigendisplacements for the given direction (not necessarily a versor).
None if the direction is not present or nac_eigendisplacements has not been calculated.


* **Parameters**

    **direction** – the direction as a list of 3 elements



* **Returns**

    the eigendisplacements as a numpy array of complex numbers with shape
    (3\*len(structure), len(structure), 3). None if not found.



#### get_nac_frequencies_along_dir(direction)
Returns the nac_frequencies for the given direction (not necessarily a versor).
None if the direction is not present or nac_frequencies has not been calculated.


* **Parameters**

    **direction** – the direction as a list of 3 elements



* **Returns**

    the frequencies as a numpy array o(3\*len(structure), len(qpoints)).
    None if not found.



#### _property_ has_eigendisplacements(_: boo_ )
True if eigendisplacements are present.


#### has_imaginary_freq(tol: float = 1e-05)
True if imaginary frequencies are present in the BS.


#### _property_ has_nac(_: boo_ )
True if nac_frequencies are present.


#### min_freq()
Returns the point where the minimum frequency is reached and its value.


### _class_ PhononBandStructureSymmLine(qpoints, frequencies, lattice, has_nac=False, eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure=None)
Bases: `PhononBandStructure`

This object stores phonon band structures along selected (symmetry) lines in the
Brillouin zone. We call the different symmetry lines (ex: \\Gamma to Z)
“branches”.


* **Parameters**


    * **qpoints** – list of qpoints as numpy arrays, in frac_coords of the
    given lattice by default


    * **frequencies** – list of phonon frequencies in eV as a numpy array with shape
    (3\*len(structure), len(qpoints))


    * **lattice** – The reciprocal lattice as a pymatgen Lattice object.
    Pymatgen uses the physics convention of reciprocal lattice vectors
    WITH a 2\*pi coefficient


    * **has_nac** – specify if the band structure has been produced taking into account
    non-analytical corrections at Gamma. If True frequencies at Gamma from
    different directions will be stored in naf. Default False.


    * **eigendisplacements** – the phonon eigendisplacements associated to the
    frequencies in Cartesian coordinates. A numpy array of complex
    numbers with shape (3\*len(structure), len(qpoints), len(structure), 3).
    he First index of the array refers to the band, the second to the index
    of the qpoint, the third to the atom in the structure and the fourth
    to the Cartesian coordinates.


    * **labels_dict** – (dict) of {} this links a qpoint (in frac coords or
    Cartesian coordinates depending on the coords) to a label.


    * **coords_are_cartesian** – Whether the qpoint coordinates are cartesian.


    * **structure** – The crystal structure (as a pymatgen Structure object)
    associated with the band structure. This is needed if we
    provide projections to the band structure.



#### _reuse_init(eigendisplacements, frequencies, has_nac, qpoints)

#### as_dict()
Returns: MSONable dict.


#### as_phononwebsite()
Return a dictionary with the phononwebsite format:
[http://henriquemiranda.github.io/phononwebsite](http://henriquemiranda.github.io/phononwebsite).


#### band_reorder()
Re-order the eigenvalues according to the similarity of the eigenvectors.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** – Dict representation.



* **Returns**

    PhononBandStructureSymmLine



#### get_branch(index)
Returns in what branch(es) is the qpoint. There can be several
branches.


* **Parameters**

    **index** – the qpoint index



* **Returns**

    A list of dictionaries [{“name”,”start_index”,”end_index”,”index”}]
    indicating all branches in which the qpoint is. It takes into
    account the fact that one qpoint (e.g., \\Gamma) can be in several
    branches



#### get_equivalent_qpoints(index)
Returns the list of qpoint indices equivalent (meaning they are the
same frac coords) to the given one.


* **Parameters**

    **index** – the qpoint index



* **Returns**

    a list of equivalent indices


TODO: now it uses the label we might want to use coordinates instead
(in case there was a mislabel)


#### write_phononwebsite(filename)
Write a json file for the phononwebsite:
[http://henriquemiranda.github.io/phononwebsite](http://henriquemiranda.github.io/phononwebsite).


### eigenvectors_from_displacements(disp, masses)
Calculate the eigenvectors from the atomic displacements.


### estimate_band_connection(prev_eigvecs, eigvecs, prev_band_order)
A function to order the phonon eigenvectors taken from phonopy.


### get_reasonable_repetitions(n_atoms: int)
Choose the number of repetitions in a supercell
according to the number of atoms in the system.

## pymatgen.phonon.dos module

This module defines classes to represent the phonon density of states, etc.


### _class_ CompletePhononDos(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), total_dos, pdoses)
Bases: `PhononDos`

This wrapper class defines a total dos, and also provides a list of PDos.


#### pdos()
Dict of partial densities of the form {Site:Densities}.
Densities are a dict of {Orbital:Values} where Values are a list of floats.
Site is a pymatgen.core.sites.Site object.


* **Type**

    dict



* **Parameters**


    * **structure** – Structure associated with this particular DOS.


    * **total_dos** – total Dos for structure


    * **pdoses** – The pdoses are supplied as a dict of {Site: Densities}.



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



### _class_ PhononDos(frequencies, densities)
Bases: `MSONable`

Basic DOS object. All other DOS objects are extended versions of this
object.


* **Parameters**


    * **frequencies** – A sequences of frequencies in THz


    * **densities** – A list representing the density of states.



#### _positive_densities()
Numpy array containing the list of densities corresponding to positive frequencies.


#### _positive_frequencies()
Numpy array containing the list of positive frequencies.


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



### coth(x)
Coth function.


* **Parameters**

    **(****)** (*x*) – value



* **Returns**

    coth(x)


## pymatgen.phonon.gruneisen module

This module provides classes to define a Grueneisen band structure.


### _class_ GruneisenParameter(qpoints, gruneisen, frequencies, multiplicities=None, structure=None, lattice=None)
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
Acoustic Debye temperature in K, i.e. the Debye temperature divided by n_sites\*\*(1/3).
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



### _class_ GruneisenPhononBandStructure(qpoints, frequencies, gruneisenparameters, lattice, eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure=None)
Bases: `PhononBandStructure`

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



### _class_ GruneisenPhononBandStructureSymmLine(qpoints, frequencies, gruneisenparameters, lattice, eigendisplacements=None, labels_dict=None, coords_are_cartesian=False, structure=None)
Bases: `GruneisenPhononBandStructure`, `PhononBandStructureSymmLine`

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



* **Returns**

    GruneisenPhononBandStructureSymmLine


## pymatgen.phonon.ir_spectra module

This module provides classes to handle the calculation of the IR spectra
This implementation is adapted from Abipy
[https://github.com/abinit/abipy](https://github.com/abinit/abipy)
where it was originally done by Guido Petretto and Matteo Giantomassi.


### _class_ IRDielectricTensor(oscillator_strength, ph_freqs_gamma, epsilon_infinity, structure)
Bases: `MSONable`

Class to handle the Ionic Dielectric Tensor
The implementation is adapted from Abipy
See the definitions Eq.(53-54) in

```
:cite:`Gonze1997`
```

 PRB55, 10355 (1997).


* **Parameters**


    * **oscillator_strength** – IR oscillator strengths as defined in Eq. 54 in


    ```
    :cite:`Gonze1997`
    ```

     PRB55, 10355 (1997).


    * **ph_freqs_gamma** – Phonon frequencies at the Gamma point


    * **epsilon_infinity** – electronic susceptibility as defined in Eq. 29.


    * **structure** – A Structure object corresponding to the structure used for the calculation.



#### as_dict()
JSON-serializable dict representation of IRDielectricTensor.


#### _classmethod_ from_dict(d)
Returns IRDielectricTensor from dict representation.


#### get_ir_spectra(broad=5e-05, emin=0, emax=None, divs=500)
The IR spectra is obtained for the different directions.


* **Parameters**


    * **broad** – a list of broadenings or a single broadening for the phonon peaks


    * **emin** (*float*) – minimum energy in which to obtain the spectra. Defaults to 0.


    * **emax** (*float*) – maximum energy in which to obtain the spectra. Defaults to None.


    * **divs** – number of frequency samples between emin and emax



* **Returns**

    divs array with the frequencies at which the

        dielectric tensor is calculated

    dielectric_tensor: divsx3x3 numpy array with the dielectric tensor

        for the range of frequencies




* **Return type**

    frequencies



#### get_plotter(components=('xx',), reim='reim', broad=5e-05, emin=0, emax=None, divs=500, \*\*kwargs)
Return an instance of the Spectrum plotter containing the different requested components.


* **Parameters**


    * **components** – A list with the components of the dielectric tensor to plot.
    Can be either two indexes or a string like ‘xx’ to plot the (0,0) component


    * **reim** – If ‘re’ (im) is present in the string plots the real (imaginary) part of the dielectric tensor


    * **broad** (*float*) – a list of broadenings or a single broadening for the phonon peaks. Defaults to 0.00005.


    * **emin** (*float*) – minimum energy in which to obtain the spectra. Defaults to 0.


    * **emax** (*float*) – maximum energy in which to obtain the spectra. Defaults to None.


    * **divs** – number of frequency samples between emin and emax


    * **\*\*kwargs** – Passed to IRDielectricTensor.get_spectrum()



#### get_spectrum(component, reim, broad=5e-05, emin=0, emax=None, divs=500, label=None)
component: either two indexes or a string like ‘xx’ to plot the (0,0) component
reim: only “re” or “im”
broad: a list of broadenings or a single broadening for the phonon peaks.


#### _property_ max_phfreq()
Maximum phonon frequency.


#### _property_ nph_freqs()
Number of phonon frequencies.


#### plot(components=('xx',), reim='reim', show_phonon_frequencies=True, xlim=None, ylim=None, \*\*kwargs)
Helper function to generate the Spectrum plotter and directly plot the results.


* **Parameters**


    * **components** – A list with the components of the dielectric tensor to plot.
    Can be either two indexes or a string like ‘xx’ to plot the (0,0) component


    * **reim** – If ‘re’ (im) is present in the string plots the real (imaginary) part of the dielectric tensor


    * **show_phonon_frequencies** – plot a dot where the phonon frequencies are to help identify IR inactive modes


    * **xlim** – x-limits of the plot. Defaults to None for automatic determination.


    * **ylim** – y-limits of the plot. Defaults to None for automatic determination.


    * **kwargs** – keyword arguments passed to the plotter


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |  |  |  |  |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### write_json(filename)
Save a json file with this data.

## pymatgen.phonon.plotter module

This module implements plotter for DOS and band structure.


### _class_ FreqUnits(factor, label)
Bases: `tuple`

Create new instance of FreqUnits(factor, label)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('factor', 'label'_ )

#### _classmethod_ _make(iterable)
Make a new FreqUnits object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new FreqUnits object replacing specified fields with new values


#### factor()
Alias for field number 0


#### label()
Alias for field number 1


### _class_ GruneisenPhononBSPlotter(bs)
Bases: `PhononBSPlotter`

Class to plot or get data to facilitate the plot of band structure objects.


* **Parameters**

    **bs** – A GruneisenPhononBandStructureSymmLine object.



#### bs_plot_data()
Get the data nicely formatted for a plot.


* **Returns**

    ticks: A dict with the ‘distances’ at which there is a qpoint (the
    x axis) and the labels (None if no label)
    frequencies: A list (one element for each branch) of frequencies for
    each qpoint: [branch][qpoint][mode]. The data is
    stored by branch to facilitate the plotting
    gruneisen: GruneisenPhononBandStructureSymmLine
    lattice: The reciprocal lattice.



* **Return type**

    A dict of the following format



#### get_plot_gs(ylim=None)
Get a matplotlib object for the gruneisen bandstructure plot.


* **Parameters**

    **ylim** – Specify the y-axis (gruneisen) limits; by default None let
    the code choose.



#### plot_compare_gs(other_plotter: GruneisenPhononBSPlotter)
Plot two band structure for comparison. One is in red the other in blue.
The two band structures need to be defined on the same symmetry lines!
and the distance between symmetry lines is
the one of the band structure used to build the PhononBSPlotter.


* **Parameters**

    **other_plotter** (*GruneisenPhononBSPlotter*) – another phonon DOS plotter defined along
    the same symmetry lines.



* **Raises**

    **ValueError** – if the two plotters are incompatible (due to different data lengths)



* **Returns**

    a matplotlib object with both band structures



#### save_plot_gs(filename, img_format='eps', ylim=None)
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **ylim** – Specifies the y-axis limits.



#### show_gs(ylim=None)
Show the plot using matplotlib.


* **Parameters**

    **ylim** – Specifies the y-axis limits.



### _class_ GruneisenPlotter(gruneisen)
Bases: `object`

Class to plot Gruneisenparameter Object.

Class to plot information from Gruneisenparameter Object.


* **Parameters**

    **gruneisen** – GruneisenParameter Object.



#### get_plot(marker='o', markersize=None, units='thz')
Will produce a plot.


* **Parameters**


    * **marker** – marker for the depiction


    * **markersize** – size of the marker


    * **units** – unit for the plots, accepted units: thz, ev, mev, ha, cm-1, cm^-1.



* **Returns**

    matplotlib axes object



* **Return type**

    plt.Axes



#### save_plot(filename, img_format='pdf', units='thz')
Will save the plot to a file.


* **Parameters**


    * **filename** – name of the filename


    * **img_format** – format of the saved plot


    * **units** – accepted units: thz, ev, mev, ha, cm-1, cm^-1.



#### show(units='thz')
Will show the plot.


* **Parameters**

    **units** – units for the plot, accepted units: thz, ev, mev, ha, cm-1, cm^-1.



### _class_ PhononBSPlotter(bs)
Bases: `object`

Class to plot or get data to facilitate the plot of band structure objects.


* **Parameters**

    **bs** – A PhononBandStructureSymmLine object.



#### _get_weight(vec: ndarray, indices: list[list[int]])
Compute the weight for each combination of sites according to the
eigenvector.


#### _static_ _make_color(colors: list[int])
Convert the eigendisplacements to rgb colors.


#### _make_ticks(ax: Axes)
Utility private method to add ticks to a band structure.


#### bs_plot_data()
Get the data nicely formatted for a plot.


* **Returns**

    ticks: A dict with the ‘distances’ at which there is a qpoint (the
    x axis) and the labels (None if no label)
    frequencies: A list (one element for each branch) of frequencies for
    each qpoint: [branch][qpoint][mode]. The data is
    stored by branch to facilitate the plotting
    lattice: The reciprocal lattice.



* **Return type**

    A dict of the following format



#### get_plot(ylim=None, units='thz')
Get a matplotlib object for the bandstructure plot.


* **Parameters**


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### get_proj_plot(site_comb: str | list[list[int]] = 'element', ylim: tuple[None | float, None | float] | None = None, units: str = 'thz', rgb_labels: tuple[None | str] | None = None)
Get a matplotlib object for the bandstructure plot projected along atomic
sites.


* **Parameters**


    * **site_comb** – a list of list, for example, [[0],[1],[2,3,4]];
    the numbers in each sublist represents the indices of atoms;
    the atoms in a same sublist will be plotted in a same color;
    if not specified, unique elements are automatically grouped.


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
    Defaults to “thz”.


    * **rgb_labels** – a list of rgb colors for the labels; if not specified,
    the colors will be automatically generated.



#### get_ticks()
Get all ticks and labels for a band structure plot.


* **Returns**

    a list of distance at which ticks should
    be set and ‘label’: a list of label for each of those ticks.



* **Return type**

    A dict with ‘distance’



#### plot_brillouin()
Plot the Brillouin zone.


#### plot_compare(other_plotter, units='thz')
Plot two band structure for comparison. One is in red the other in blue.
The two band structures need to be defined on the same symmetry lines!
and the distance between symmetry lines is the one of the band structure
used to build the PhononBSPlotter.


* **Parameters**


    * **other_plotter** – another PhononBSPlotter object defined along the same symmetry lines


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
    Defaults to ‘thz’.



* **Returns**

    a matplotlib object with both band structures



#### save_plot(filename, img_format='eps', ylim=None, units='thz')
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### show(ylim=None, units='thz')
Show the plot using matplotlib.


* **Parameters**


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### show_proj(site_comb: str | list[list[int]] = 'element', ylim: tuple[None | float, None | float] | None = None, units: str = 'thz', rgb_labels: tuple[str] | None = None)
Show the projected plot using matplotlib.


* **Parameters**


    * **site_comb** – A list of list of indices of sites to combine. For example,
    [[0, 1], [2, 3]] will combine the projections of sites 0 and 1,
    and sites 2 and 3. Defaults to “element”, which will combine
    sites by element.


    * **ylim** – Specify the y-axis (frequency) limits; by default None let
    the code choose.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.
    Defaults to “thz”.


    * **rgb_labels** – A list of labels for the rgb triangle. Defaults to None,
    which will use the element symbols.



### _class_ PhononDosPlotter(stack=False, sigma=None)
Bases: `object`

Class for plotting phonon DOSs. The interface is extremely flexible given there are many
different ways in which people want to view DOS.
Typical usage is:

> # Initializes plotter with some optional args. Defaults are usually fine
> plotter = PhononDosPlotter().

> # Add DOS with a label
> plotter.add_dos(“Total DOS”, dos)

> # Alternatively, you can add a dict of DOSes. This is the typical form
> # returned by CompletePhononDos.get_element_dos().


* **Parameters**


    * **stack** – Whether to plot the DOS as a stacked area graph


    * **sigma** – A float specifying a standard deviation for Gaussian smearing
    the DOS for nicer looking plots. Defaults to None for no smearing.



#### add_dos(label, dos)
Adds a dos for plotting.


* **Parameters**


    * **label** – label for the DOS. Must be unique.


    * **dos** – PhononDos object



#### add_dos_dict(dos_dict, key_sort_func=None)
Add a dictionary of doses, with an optional sorting function for the
keys.


* **Parameters**


    * **dos_dict** – dict of {label: Dos}


    * **key_sort_func** – function used to sort the dos_dict keys.



#### get_dos_dict()
Returns the added doses as a json-serializable dict. Note that if you
have specified smearing for the DOS plot, the densities returned will
be the smeared densities, not the original densities.


* **Returns**

    {‘frequencies’:..,
    ‘densities’: …}}



* **Return type**

    Dict of dos data. Generally of the form, {label



#### get_plot(xlim=None, ylim=None, units='thz')
Get a matplotlib plot showing the DOS.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



#### save_plot(filename, img_format='eps', xlim=None, ylim=None, units='thz')
Save matplotlib plot to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1



#### show(xlim=None, ylim=None, units='thz')
Show the plot using matplotlib.


* **Parameters**


    * **xlim** – Specifies the x-axis limits. Set to None for automatic
    determination.


    * **ylim** – Specifies the y-axis limits.


    * **units** – units for the frequencies. Accepted values thz, ev, mev, ha, cm-1, cm^-1.



### _class_ ThermoPlotter(dos, structure=None)
Bases: `object`

Plotter for thermodynamic properties obtained from phonon DOS.
If the structure corresponding to the DOS, it will be used to extract the formula unit and provide
the plots in units of mol instead of mole-cell.


* **Parameters**


    * **dos** – A PhononDos object.


    * **structure** – A Structure object corresponding to the structure used for the calculation.



#### _plot_thermo(func, temperatures, factor=1, ax: Axes | None = None, ylabel=None, label=None, ylim=None, \*\*kwargs)
Plots a thermodynamic property for a generic function from a PhononDos instance.


* **Parameters**


    * **func** – the thermodynamic function to be used to calculate the property


    * **temperatures** – a list of temperatures


    * **factor** – a multiplicative factor applied to the thermodynamic property calculated. Used to change
    the units.


    * **ax** – matplotlib Axes or None if a new figure should be created.


    * **ylabel** – label for the y axis


    * **label** – label of the plot


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure



* **Return type**

    plt.figure



#### plot_cv(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the constant volume specific heat C_v in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure



* **Return type**

    plt.figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_entropy(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the vibrational entrpy in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure



* **Return type**

    plt.figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_helmholtz_free_energy(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the vibrational contribution to the Helmoltz free energy in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure



* **Return type**

    plt.figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_internal_energy(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots the vibrational internal energy in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure



* **Return type**

    plt.figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_thermodynamic_properties(tmin, tmax, ntemp, ylim=None, \*\*kwargs)
Plots all the thermodynamic properties in a temperature range.


* **Parameters**


    * **tmin** – minimum temperature


    * **tmax** – maximum temperature


    * **ntemp** – number of steps


    * **ylim** – tuple specifying the y-axis limits.


    * **kwargs** – kwargs passed to the matplotlib function ‘plot’.



* **Returns**

    matplotlib figure



* **Return type**

    plt.figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

### freq_units(units)

* **Parameters**

    **units** – str, accepted values: thz, ev, mev, ha, cm-1, cm^-1.



* **Returns**

    Returns conversion factor from THz to the required units and the label in the form of a namedtuple


## pymatgen.phonon.thermal_displacements module

This module provides classes to handle thermal displacement matrices (anisotropic displacement parameters).


### _class_ ThermalDisplacementMatrices(thermal_displacement_matrix_cart, structure, temperature, thermal_displacement_matrix_cif=None)
Bases: `MSONable`

Class to handle thermal displacement matrices
This class stores thermal displacement matrices in Ucart format.

An earlier implementation based on Matlab can be found here:
[https://github.com/JaGeo/MolecularToolbox](https://github.com/JaGeo/MolecularToolbox)
( J. George, A. Wang, V. L. Deringer, R. Wang, R. Dronskowski, U. Englert, CrystEngComm, 2015, 17, 7414-7422.)


* **Parameters**


    * **thermal_displacement_matrix_cart** – 2D numpy array including the thermal_displacement matrix Ucart
    1st dimension atom types, then compressed thermal displacement matrix will follow
    U11, U22, U33, U23, U13, U12 (xx, yy, zz, yz, xz, xy)
    convention similar to “thermal_displacement_matrices.yaml” in phonopy


    * **structure** – A pymatgen Structure object


    * **temperature** – temperature at which thermal displacement matrix was determined


    * **thermal_displacement_matrix_cif** – 2D numpy array including the thermal_displacement matrix Ucif format
    1st dimension atom types, then compressed thermal displacement matrix will follow
    U11, U22, U33, U23, U13, U12 (xx, yy, zz, yz, xz, xy)
    convention similar to “thermal_displacement_matrices.yaml” in phonopy.



#### _property_ B()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.


* **Returns**

    First dimension are the atoms in the structure.



* **Return type**

    np.array



#### _property_ U1U2U3()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.


* **Returns**

    eigenvalues of Ucart. First dimension are the atoms in the structure.



* **Return type**

    np.array



#### _property_ Ucif()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.


* **Returns**

    Ucif as array. First dimension are the atoms in the structure.



* **Return type**

    np.array



#### _property_ Ustar()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.


* **Returns**

    Ustar as array. First dimension are the atoms in the structure.



* **Return type**

    np.array



#### _static_ _angle_dot(a, b)

#### _property_ beta()
Computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477-480.


* **Returns**

    First dimension are the atoms in the structure.



* **Return type**

    np.array



#### compute_directionality_quality_criterion(other)
Will compute directionality of prolate displacement ellipsoids as described in
[https://doi.org/10.1039/C9CE00794F](https://doi.org/10.1039/C9CE00794F) with the earlier implementation: [https://github.com/damMroz/Angle/](https://github.com/damMroz/Angle/).


* **Parameters**


    * **other** – ThermalDisplacementMatrix


    * **compared** (*please make sure that the order** of **the atoms in both objects that are*) –


    * **Otherwise** (*is the same.*) –


    * **results** (*this analysis will deliver wrong*) –



* **Returns**

    will return a list including dicts for each atom that include “vector0”
    (largest principal axes of self object),

    > ”vector1” (largest principal axes of the other object), “angle” between both axes,

    >     These vectors can then, for example, be drawn into the structure with VESTA.
    >     Vectors are given in Cartesian coordinates




#### _static_ from_Ucif(thermal_displacement_matrix_cif, structure, temperature)
Starting from a numpy array, it will convert Ucif values into Ucart values and initialize the class.


* **Parameters**


    * **thermal_displacement_matrix_cif** – np.array,
    first dimension are the atoms,
    then reduced form of thermal displacement matrix will follow
    Order as above: U11, U22, U33, U23, U13, U12


    * **structure** – Structure object


    * **temperature** – float
    Corresponding temperature



* **Returns**

    ThermalDisplacementMatrices



#### _static_ from_cif_P1(filename: str)
Reads a cif with P1 symmetry including positions and ADPs.
Currently, no check of symmetry is performed as CifParser methods cannot be easily reused.


* **Parameters**

    **filename** – Filename of the CIF.



* **Returns**

    ThermalDisplacementMatrices



#### _static_ from_structure_with_site_properties_Ucif(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), temperature: float | None = None)
Will create this object with the help of a structure with site properties.


* **Parameters**


    * **structure** – Structure object including U11_cif, U22_cif, U33_cif, U23_cif, U13_cif, U12_cif as site


    * **properties** –


    * **temperature** – temperature for Ucif data



* **Returns**

    ThermalDisplacementMatrices



#### _static_ get_full_matrix(thermal_displacement)
Transfers the reduced matrix to the full matrix (order of reduced matrix U11, U22, U33, U23, U13, U12).


* **Parameters**

    **thermal_displacement** – 2d numpy array, first dimension are the atoms



* **Returns**

    3d numpy array including thermal displacements, first dimensions are the atoms



#### _static_ get_reduced_matrix(thermal_displacement)
Transfers the full matrix to reduced matrix (order of reduced matrix U11, U22, U33, U23, U13, U12).


* **Parameters**

    **thermal_displacement** – 2d numpy array, first dimension are the atoms



* **Returns**

    3d numpy array including thermal displacements, first dimensions are the atoms



#### _property_ ratio_prolate()
This will compute ratio between largest and smallest eigenvalue of Ucart.


#### to_structure_with_site_properties_Ucif()
Transfers this object into a structure with site properties (Ucif).
This is useful for sorting the atoms in the structure including site properties.
E.g., with code like this:
def sort_order(site):

> return [site.specie.X, site.frac_coords[0], site.frac_coords[1], site.frac_coords[2]]

new_structure0 = Structure.from_sites(sorted(structure0, key=sort_order)).


* **Returns**

    Structure



#### visualize_directionality_quality_criterion(other, filename: str = 'visualization.vesta', which_structure: int = 0)
Will create a VESTA file for visualization of the directionality criterion.


* **Parameters**


    * **other** – ThermalDisplacementMatrices


    * **filename** – Filename of the VESTA file


    * **which_structure** – 0 means structure of the self object will be used, 1 means structure of the other
    object will be used



#### write_cif(filename)
Writes a cif including thermal displacements.


* **Parameters**

    **filename** – name of the cif file