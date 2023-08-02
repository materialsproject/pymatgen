---
layout: default
title: pymatgen.analysis.diffraction.neutron.md
nav_exclude: true
---

# pymatgen.analysis.diffraction.neutron module

This module implements a neutron diffraction (ND) pattern calculator.


### _class_ pymatgen.analysis.diffraction.neutron.NDCalculator(wavelength=1.54184, symprec: float = 0, debye_waller_factors=None)
Bases: [`AbstractDiffractionPatternCalculator`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator)

Computes the powder neutron diffraction pattern of a crystal structure.
This code is a slight modification of XRDCalculator in
pymatgen.analysis.diffraction.xrd. See it for details of the algorithm.
Main changes by using neutron instead of X-ray are as follows:


1. Atomic scattering length is a constant.


2. Polarization correction term of Lorentz factor is unnecessary.

Reference:
Marc De Graef and Michael E. McHenry, Structure of Materials 2nd ed,
Chapter13, Cambridge University Press 2003.

Initializes the ND calculator with a given radiation.


* **Parameters**


    * **wavelength** (*float*) – The wavelength of neutron in angstroms.
    Defaults to 1.54, corresponds to Cu K_alpha x-ray radiation.


    * **symprec** (*float*) – Symmetry precision for structure refinement. If
    set to 0, no refinement is done. Otherwise, refinement is
    performed using spglib with provided precision.


    * **symbol** (*debye_waller_factors** (**{element*) – float}): Allows the
    specification of Debye-Waller factors. Note that these
    factors are temperature dependent.



#### get_pattern(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), scaled=True, two_theta_range=(0, 90))
Calculates the powder neutron diffraction pattern for a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **scaled** (*bool*) – Whether to return scaled intensities. The maximum
    peak is set to a value of 100. Defaults to True. Use False if
    you need the absolute values to combine ND plots.


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.



* **Returns**

    (NDPattern)