---
layout: default
title: pymatgen.analysis.diffraction.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.diffraction package

This package implements various diffraction analyses.


## pymatgen.analysis.diffraction.core module

This module implements core classes for calculation of diffraction patterns.


### _class_ AbstractDiffractionPatternCalculator()
Bases: `ABC`

Abstract base class for computing the diffraction pattern of a crystal.


#### SCALED_INTENSITY_TOL(_ = 0.00_ )

#### TWO_THETA_TOL(_ = 1e-0_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ get_pattern(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), scaled=True, two_theta_range=(0, 90))
Calculates the diffraction pattern for a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **scaled** (*bool*) – Whether to return scaled intensities. The maximum
    peak is set to a value of 100. Defaults to True. Use False if
    you need the absolute values to combine XRD plots.


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.



* **Returns**

    (DiffractionPattern)



#### get_plot(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), two_theta_range: tuple[float, float] = (0, 90), annotate_peaks='compact', ax: plt.Axes = None, with_labels=True, fontsize=16)
Returns the diffraction plot as a matplotlib Axes.


* **Parameters**


    * **structure** – Input structure


    * **two_theta_range** (*tuple**[**float**, **float**]*) – Range of two_thetas to calculate in degrees.
    Defaults to (0, 90). Set to None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.


    * **annotate_peaks** (*str** or **None*) – Whether and how to annotate the peaks
    with hkl indices. Default is ‘compact’, i.e. show short
    version (oriented vertically), e.g. 100. If ‘full’, show
    long version, e.g. (1, 0, 0). If None, do not show anything.


    * **ax** – matplotlib Axes or None if a new figure should be
    created.


    * **with_labels** – True to add xlabels and ylabels to the plot.


    * **fontsize** – (int) fontsize for peak labels.



* **Returns**

    matplotlib Axes object



* **Return type**

    plt.Axes



#### plot_structures(structures, fontsize=6, \*\*kwargs)
Plot diffraction patterns for multiple structures on the same figure.


* **Parameters**


    * **structures** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – List of structures


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
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |
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

#### show_plot(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), \*\*kwargs)
Shows the diffraction plot.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.


    * **annotate_peaks** (*str** or **None*) – Whether and how to annotate the peaks
    with hkl indices. Default is ‘compact’, i.e. show short
    version (oriented vertically), e.g. 100. If ‘full’, show
    long version, e.g. (1, 0, 0). If None, do not show anything.



### _class_ DiffractionPattern(x, y, hkls, d_hkls)
Bases: [`Spectrum`](pymatgen.core.md#pymatgen.core.spectrum.Spectrum)

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

### get_unique_families(hkls)
Returns unique families of Miller indices. Families must be permutations
of each other.


* **Parameters**

    **hkls** (*[**h**, **k**, **l**]*) – List of Miller indices.



* **Returns**

    multiplicity}: A dict with unique hkl and multiplicity.



* **Return type**

    {hkl


## pymatgen.analysis.diffraction.neutron module

This module implements a neutron diffraction (ND) pattern calculator.


### _class_ NDCalculator(wavelength=1.54184, symprec: float = 0, debye_waller_factors=None)
Bases: `AbstractDiffractionPatternCalculator`

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



#### _abc_impl(_ = <_abc._abc_data object_ )

#### get_pattern(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), scaled=True, two_theta_range=(0, 90))
Calculates the powder neutron diffraction pattern for a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **scaled** (*bool*) – Whether to return scaled intensities. The maximum
    peak is set to a value of 100. Defaults to True. Use False if
    you need the absolute values to combine ND plots.


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.



* **Returns**

    (NDPattern)


## pymatgen.analysis.diffraction.tem module

This module implements a TEM pattern calculator.


### _class_ TEMCalculator(symprec: float | None = None, voltage: float = 200, beam_direction: tuple[int, int, int] = (0, 0, 1), camera_length: int = 160, debye_waller_factors: dict[str, float] | None = None, cs: float = 1)
Bases: `AbstractDiffractionPatternCalculator`

Computes the TEM pattern of a crystal structure for multiple Laue zones.
Code partially inspired from XRD calculation implementation. X-ray factor to electron factor

> conversion based on the International Table of Crystallography.

#TODO: Could add “number of iterations”, “magnification”, “critical value of beam”,

    “twin direction” for certain materials, “sample thickness”, and “excitation error s”.


* **Parameters**


    * **symprec** (*float*) – Symmetry precision for structure refinement. If
    set to 0, no refinement is done. Otherwise, refinement is
    performed using spglib with provided precision.


    * **voltage** (*float*) – The wavelength is a function of the TEM microscope’s
    voltage (in kV). Defaults to 200.


    * **beam_direction** (*tuple*) – The direction of the electron beam fired onto the sample.
    By default, set to [0,0,1], which corresponds to the normal direction
    of the sample plane.


    * **camera_length** (*int*) – The distance from the sample to the projected diffraction pattern.
    By default, set to 160 cm. Units in cm.


    * **symbol** (*debye_waller_factors** (**{element*) – float}): Allows the
    specification of Debye-Waller factors. Note that these
    factors are temperature dependent.


    * **cs** (*float*) – The chromatic aberration coefficient (in mm). Defaults to 1.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### bragg_angles(interplanar_spacings: dict[tuple[int, int, int], float])
Gets the Bragg angles for every hkl point passed in (where n = 1).


* **Parameters**

    **interplanar_spacings** (*dict*) – dictionary of hkl to interplanar spacing



* **Returns**

    dict of hkl plane (3-tuple) to Bragg angle in radians (float)



#### cell_intensity(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates cell intensity for each hkl plane. For simplicity’s sake, take I =

```
|
```

F|\*\*2.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict of hkl plane to cell intensity



#### cell_scattering_factors(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates the scattering factor for the whole cell.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict of hkl plane (3-tuple) to scattering factor (in angstroms).



#### electron_scattering_factors(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates atomic scattering factors for electrons using the Mott-Bethe formula (1st order Born approximation).


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict from atomic symbol to another dict of hkl plane to factor (in angstroms)



#### _static_ generate_points(coord_left: int = -10, coord_right: int = 10)
Generates a bunch of 3D points that span a cube.


* **Parameters**


    * **coord_left** (*int*) – The minimum coordinate value.


    * **coord_right** (*int*) – The maximum coordinate value.



* **Returns**

    2d array



* **Return type**

    np.array



#### get_first_point(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), points: list)
Gets the first point to be plotted in the 2D DP, corresponding to maximum d/minimum R.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **points** (*list*) – All points to be checked.



* **Returns**

    dict of a hkl plane to max interplanar distance.



#### _static_ get_interplanar_angle(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), p1: tuple[int, int, int], p2: tuple[int, int, int])
Returns the interplanar angle (in degrees) between the normal of two crystal planes.
Formulas from International Tables for Crystallography Volume C pp. 2-9.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **p1** (*3-tuple*) – plane 1


    * **p2** (*3-tuple*) – plane 2



* **Returns**

    float



#### get_interplanar_spacings(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), points: list[tuple[int, int, int]] | np.ndarray)

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – the input structure.


    * **points** (*tuple*) – the desired hkl indices.



* **Returns**

    Dict of hkl to its interplanar spacing, in angstroms (float).



#### get_pattern(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), scaled: bool | None = None, two_theta_range: tuple[float, float] | None = None)
Returns all relevant TEM DP info in a pandas dataframe.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **scaled** (*bool*) – Required value for inheritance, does nothing in TEM pattern


    * **two_theta_range** (*tuple**[**float**, **float**]*) – Required value for inheritance, does nothing in TEM pattern



* **Returns**

    pd.DataFrame



#### get_plot_2d(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Generates the 2D diffraction pattern of the input structure.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.



* **Returns**

    Figure



#### get_plot_2d_concise(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Generates the concise 2D diffraction pattern of the input structure of a smaller size and without layout.
Does not display.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.



* **Returns**

    Figure



#### _static_ get_plot_coeffs(p1: tuple[int, int, int], p2: tuple[int, int, int], p3: tuple[int, int, int])
Calculates coefficients of the vector addition required to generate positions for each DP point
by the Moore-Penrose inverse method.


* **Parameters**


    * **p1** (*3-tuple*) – The first point. Fixed.


    * **p2** (*3-tuple*) – The second point. Fixed.


    * **p3** (*3-tuple*) – The point whose coefficients are to be calculted.



* **Returns**

    Numpy array



#### get_positions(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), points: list)
Calculates all the positions of each hkl point in the 2D diffraction pattern by vector addition.
Distance in centimeters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **points** (*list*) – All points to be checked.



* **Returns**

    dict of hkl plane to xy-coordinates.



#### get_s2(bragg_angles: dict[tuple[int, int, int], float])
Calculates the s squared parameter (= square of sin theta over lambda) for each hkl plane.


* **Parameters**

    **bragg_angles** (*dict*) – The bragg angles for each hkl plane.



* **Returns**

    Dict of hkl plane to s2 parameter, calculates the s squared parameter

        (= square of sin theta over lambda).




#### is_parallel(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), plane: tuple[int, int, int], other_plane: tuple[int, int, int])
Checks if two hkl planes are parallel in reciprocal space.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **plane** (*3-tuple*) – The first plane to be compared.


    * **other_plane** (*3-tuple*) – The other plane to be compared.



* **Returns**

    boolean



#### normalized_cell_intensity(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Normalizes the cell_intensity dict to 1, for use in plotting.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict of hkl plane to normalized cell intensity



#### tem_dots(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), points)
Generates all TEM_dot as named tuples that will appear on the 2D diffraction pattern.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **points** (*list*) – All points to be checked.



* **Returns**

    list of TEM_dots



#### wavelength_rel()
Calculates the wavelength of the electron beam with relativistic kinematic effects taken

    into account.


* **Returns**

    Relativistic Wavelength (in angstroms)



* **Return type**

    float



#### x_ray_factors(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates x-ray factors, which are required to calculate atomic scattering factors. Method partially inspired
by the equivalent process in the xrd module.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict*) – Dictionary of hkl plane to Bragg angle.



* **Returns**

    dict of atomic symbol to another dict of hkl plane to x-ray factor (in angstroms).



#### zone_axis_filter(points: list[tuple[int, int, int]] | np.ndarray, laue_zone: int = 0)
Filters out all points that exist within the specified Laue zone according to the zone axis rule.


* **Parameters**


    * **points** (*np.ndarray*) – The list of points to be filtered.


    * **laue_zone** (*int*) – The desired Laue zone.



* **Returns**

    list of 3-tuples


## pymatgen.analysis.diffraction.xrd module

This module implements an XRD pattern calculator.


### _class_ XRDCalculator(wavelength='CuKa', symprec: float = 0, debye_waller_factors=None)
Bases: `AbstractDiffractionPatternCalculator`

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
within the limiting sphere given by frac{2}{lambda}.


2. For each reciprocal point mathbf{g_{hkl}} corresponding to
lattice plane (hkl), compute the Bragg condition
sin(theta) = frac{ lambda}{2d_{hkl}}


3. Compute the structure factor as the sum of the atomic scattering
factors. The atomic scattering factors are given by

> f(s) = Z - 41.78214 times s^2 times sum limits_{i=1}^n a_i exp(-b_is^2)

where s = frac{sin(theta)}{lambda} and a_i
and b_i are the fitted parameters for each element. The
structure factor is then given by

> F_{hkl} = sum limits_{j=1}^N f_j  exp(2 pi i  mathbf{g_{hkl}} cdot  mathbf{r})


4. The intensity is then given by the modulus square of the structure factor.

> I_{hkl} = F_{hkl}F_{hkl}^\*


5. Finally, the Lorentz polarization correction factor is applied. This
factor is given by:

> P(theta) = frac{1 + cos^2(2 theta)}{sin^2(theta) cos(theta)}

Initializes the XRD calculator with a given radiation.


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

#### _abc_impl(_ = <_abc._abc_data object_ )

#### get_pattern(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), scaled=True, two_theta_range=(0, 90))
Calculates the diffraction pattern for a structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **scaled** (*bool*) – Whether to return scaled intensities. The maximum
    peak is set to a value of 100. Defaults to True. Use False if
    you need the absolute values to combine XRD plots.


    * **two_theta_range** (*[**float** of **length 2**]*) – Tuple for range of
    two_thetas to calculate in degrees. Defaults to (0, 90). Set to
    None if you want all diffracted beams within the limiting
    sphere of radius 2 / wavelength.



* **Returns**

    (XRDPattern)