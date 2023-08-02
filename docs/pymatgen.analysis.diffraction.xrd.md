---
layout: default
title: pymatgen.analysis.diffraction.xrd.md
nav_exclude: true
---

# pymatgen.analysis.diffraction.xrd module

This module implements an XRD pattern calculator.


### _class_ pymatgen.analysis.diffraction.xrd.XRDCalculator(wavelength='CuKa', symprec: float = 0, debye_waller_factors=None)
Bases: [`AbstractDiffractionPatternCalculator`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator)

Computes the XRD pattern of a crystal structure.

This code is implemented by Shyue Ping Ong as part of UCSD’s NANO106 -
Crystallography of Materials. The formalism for this code is based on
that given in Chapters 11 and 12 of Structure of Materials by Marc De
Graef and Michael E. McHenry. This takes into account the atomic
scattering factors and the Lorentz polarization factor, but not
the Debye-Waller (temperature) factor (for which data is typically not
available). Note that the multiplicity correction is not needed since
this code simply goes through all reciprocal points within the limiting
sphere, which includes all symmetrically equivalent facets. The algorithm
is as follows


1. Calculate reciprocal lattice of structure. Find all reciprocal points
within the limiting sphere given by :math:\` frac{2}{ lambda}\`.


2. For each reciprocal point :math:\` mathbf{g_{hkl}}\` corresponding to
lattice plane $(hkl)$, compute the Bragg condition
:math:\` sin( theta) =  frac{ lambda}{2d_{hkl}}\`


3. Compute the structure factor as the sum of the atomic scattering
factors. The atomic scattering factors are given by

f(s) = Z - 41.78214 \\times s^2 \\times \\sum \\limits_{i=1}^n a_i \\
 \\exp(-b_is^2)where $s = \ frac{\ sin(\ theta)}{\ lambda}$ and $a_i$
and $b_i$ are the fitted parameters for each element. The
structure factor is then given by

F_{hkl} =  \\sum \\limits_{j=1}^N f_j  \\exp(2 \\pi i  \\mathbf{g_{hkl}}
 \\cdot  \\mathbf{r})
4. The intensity is then given by the modulus square of the structure
factor.

I_{hkl} = F_{hkl}F_{hkl}^\*
5. Finally, the Lorentz polarization correction factor is applied. This
factor is given by:

P( \\theta) = \\frac{1 +  \\cos^2(2 \\theta)}
{ \\sin^2( \\theta) \\cos( \\theta)}Initializes the XRD calculator with a given radiation.


* **Parameters**


    * **wavelength** (*str/float*) – The wavelength can be specified as either a
    float or a string. If it is a string, it must be one of the
    supported definitions in the AVAILABLE_RADIATION class
    variable, which provides useful commonly used wavelengths.
    If it is a float, it is interpreted as a wavelength in
    angstroms. Defaults to “CuKa”, i.e, Cu K_alpha radiation.


    * **symprec** (*float*) – Symmetry precision for structure refinement. If
    set to 0, no refinement is done. Otherwise, refinement is
    performed using spglib with provided precision.


    * **symbol** (*debye_waller_factors** (**{element*) – float}): Allows the
    specification of Debye-Waller factors. Note that these
    factors are temperature dependent.



#### AVAILABLE_RADIATION(_ = ('CuKa', 'CuKa2', 'CuKa1', 'CuKb1', 'MoKa', 'MoKa2', 'MoKa1', 'MoKb1', 'CrKa', 'CrKa2', 'CrKa1', 'CrKb1', 'FeKa', 'FeKa2', 'FeKa1', 'FeKb1', 'CoKa', 'CoKa2', 'CoKa1', 'CoKb1', 'AgKa', 'AgKa2', 'AgKa1', 'AgKb1'_ )

#### get_pattern(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), scaled=True, two_theta_range=(0, 90))
Calculates the diffraction pattern for a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **scaled** (*bool*) – Whether to return scaled intensities. The maximum
    peak is set to a value of 100. Defaults to True. Use False if
    you need the absolute values to combine XRD plots.


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.



* **Returns**

    (XRDPattern)