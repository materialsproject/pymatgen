---
layout: default
title: pymatgen.analysis.excitation.md
nav_exclude: true
---

# pymatgen.analysis.excitation module

This module defines an excitation spectrum class.


### _class_ pymatgen.analysis.excitation.ExcitationSpectrum(x, y)
Bases: [`Spectrum`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum)

Basic excitation spectrum object.

<!-- attribute: x
The sequence of energies -->
<!-- attribute: y
The sequence of mu(E) -->

* **Parameters**


    * **x** – A sequence of x-ray energies in eV


    * **y** – A sequence of intensity values.



#### XLABEL(_ = 'Energy (eV)_ )

#### YLABEL(_ = 'Intensity_ )