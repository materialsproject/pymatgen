---
layout: default
title: pymatgen.analysis.diffraction.core.md
nav_exclude: true
---

# pymatgen.analysis.diffraction.core module

This module implements core classes for calculation of diffraction patterns.


### _class_ pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator()
Bases: `ABC`

Abstract base class for computing the diffraction pattern of a crystal.


#### SCALED_INTENSITY_TOL(_ = 0.00_ )

#### TWO_THETA_TOL(_ = 1e-0_ )

#### _abstract_ get_pattern(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), scaled=True, two_theta_range=(0, 90))
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

    (DiffractionPattern)



#### get_plot(structure, two_theta_range=(0, 90), annotate_peaks='compact', ax=None, with_labels=True, fontsize=16)
Returns the diffraction plot as a matplotlib.pyplot.


* **Parameters**


    * **structure** – Input structure


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.


    * **annotate_peaks** (*str** or **None*) – Whether and how to annotate the peaks
    with hkl indices. Default is ‘compact’, i.e. show short
    version (oriented vertically), e.g. 100. If ‘full’, show
    long version, e.g. (1, 0, 0). If None, do not show anything.


    * **ax** – matplotlib `Axes` or None if a new figure should be
    created.


    * **with_labels** – True to add xlabels and ylabels to the plot.


    * **fontsize** – (int) fontsize for peak labels.



* **Returns**

    (matplotlib.pyplot)



#### plot_structures(structures, fontsize=6, \*\*kwargs)
Plot diffraction patterns for multiple structures on the same figure.


* **Parameters**


    * **structures** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – List of structures


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.


    * **annotate_peaks** (*str** or **None*) – Whether and how to annotate the peaks
    with hkl indices. Default is ‘compact’, i.e. show short
    version (oriented vertically), e.g. 100. If ‘full’, show
    long version, e.g. (1, 0, 0). If None, do not show anything.


    * **fontsize** – (int) fontsize for peak labels.


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------ | ------- |
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

#### show_plot(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), \*\*kwargs)
Shows the diffraction plot.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Input structure


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.


    * **annotate_peaks** (*str** or **None*) – Whether and how to annotate the peaks
    with hkl indices. Default is ‘compact’, i.e. show short
    version (oriented vertically), e.g. 100. If ‘full’, show
    long version, e.g. (1, 0, 0). If None, do not show anything.



### _class_ pymatgen.analysis.diffraction.core.DiffractionPattern(x, y, hkls, d_hkls)
Bases: [`Spectrum`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum)

A representation of a diffraction pattern.


* **Parameters**


    * **x** – Two theta angles.


    * **y** – Intensities


    * **hkls** – [{“hkl”: (h, k, l), “multiplicity”: mult}],
    where {“hkl”: (h, k, l), “multiplicity”: mult}
    is a dict of Miller
    indices for all diffracted lattice facets contributing to each
    intensity.


    * **d_hkls** – List of interplanar spacings.



#### XLABEL(_ = '$2\\\\Theta$_ )

#### YLABEL(_ = 'Intensity_ )

### pymatgen.analysis.diffraction.core.get_unique_families(hkls)
Returns unique families of Miller indices. Families must be permutations
of each other.


* **Parameters**

    **hkls** (*[**h**, **k**, **l**]*) – List of Miller indices.



* **Returns**

    multiplicity}: A dict with unique hkl and multiplicity.



* **Return type**

    {hkl