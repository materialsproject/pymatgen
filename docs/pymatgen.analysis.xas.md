---
layout: default
title: pymatgen.analysis.xas.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.xas package

Package for analysis of X-ray Absorption Spectroscopy.


## pymatgen.analysis.xas.spectrum module

This module defines classes to represent all xas and stitching methods.


### _class_ XAS(x, y, structure, absorbing_element, edge='K', spectrum_type='XANES', absorbing_index=None)
Bases: [`Spectrum`](pymatgen.core.md#pymatgen.core.spectrum.Spectrum)

Basic XAS object.


* **Parameters**


    * **x** – A sequence of x-ray energies in eV


    * **y** – A sequence of mu(E)


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure associated with the spectrum


    * **absorbing_element** ([*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)) – Element associated with the spectrum


    * **edge** (*str*) – Absorption edge associated with the spectrum


    * **spectrum_type** (*str*) – ‘XANES’ or ‘EXAFS’


    * **absorbing_index** (*None** or **int*) – If None, the spectrum is assumed to be a
    site-weighted spectrum, which is comparable to experimental one.
    Otherwise, it indicates that the absorbing_index for a site-wise spectrum.



#### x()
The sequence of energies.


* **Type**

    Sequence[float]



#### y()
The sequence of mu(E).


* **Type**

    Sequence[float]



#### absorbing_element()
The absorbing element of the spectrum.


* **Type**

    str



#### edge()
The edge of the spectrum.


* **Type**

    str



#### spectrum_type()
The type of the spectrum (XANES or EXAFS).


* **Type**

    str



#### absorbing_index()
The absorbing index of the spectrum.


* **Type**

    int


Initializes a spectrum object.


#### XLABEL(_ = 'Energy_ )

#### YLABEL(_ = 'Intensity_ )

#### stitch(other: XAS, num_samples: int = 500, mode: Literal['XAFS', 'L23'] = 'XAFS')
Stitch XAS objects to get the full XAFS spectrum or L23 edge XANES
spectrum depending on the mode.


1. Use XAFS mode for stitching XANES and EXAFS with same absorption edge.

    The stitching will be performed based on wavenumber, k.
    for k <= 3, XAS(k) = XAS[XANES(k)]
    for 3 < k < max(xanes_k), will interpolate according to

    > XAS(k)=f(k)\*mu[XANES(k)]+(1-f(k))\*mu[EXAFS(k)]
    > where f(k)=cos^2((pi/2) (k-3)/(max(xanes_k)-3)

    for k > max(xanes_k), XAS(k) = XAS[EXAFS(k)]


2. Use L23 mode for stitching L2 and L3 edge XANES for elements with

    atomic number <=30.


* **Parameters**


    * **other** – Another XAS object.


    * **num_samples** (*int*) – Number of samples for interpolation.


    * **mode** (*"XAFS"** | **"L23"*) – Either XAFS mode for stitching XANES and EXAFS
    or L23 mode for stitching L2 and L3.



* **Returns**

    The stitched spectrum.



* **Return type**

    XAS object



### site_weighted_spectrum(xas_list: list[XAS], num_samples: int = 500)
Obtain site-weighted XAS object based on site multiplicity for each
absorbing index and its corresponding site-wise spectrum.


* **Parameters**


    * **xas_list** (*[**XAS**]*) – List of XAS object to be weighted


    * **num_samples** (*int*) – Number of samples for interpolation



* **Returns**

    The site-weighted spectrum



* **Return type**

    XAS object