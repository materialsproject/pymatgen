---
layout: default
title: pymatgen.analysis.diffraction.tem.md
nav_exclude: true
---

# pymatgen.analysis.diffraction.tem module

This module implements a TEM pattern calculator.


### _class_ pymatgen.analysis.diffraction.tem.TEMCalculator(symprec: float | None = None, voltage: float = 200, beam_direction: tuple[int, int, int] = (0, 0, 1), camera_length: int = 160, debye_waller_factors: dict[str, float] | None = None, cs: float = 1)
Bases: [`AbstractDiffractionPatternCalculator`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator)

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
    voltage. By default, set to 200 kV. Units in kV.


    * **beam_direction** (*tuple*) – The direction of the electron beam fired onto the sample.
    By default, set to [0,0,1], which corresponds to the normal direction
    of the sample plane.


    * **camera_length** (*int*) – The distance from the sample to the projected diffraction pattern.
    By default, set to 160 cm. Units in cm.


    * **symbol** (*debye_waller_factors** (**{element*) – float}): Allows the
    specification of Debye-Waller factors. Note that these
    factors are temperature dependent.


    * **cs** (*float*) – the chromatic aberration coefficient. set by default to 1 mm.



#### bragg_angles(interplanar_spacings: dict[tuple[int, int, int], float])
Gets the Bragg angles for every hkl point passed in (where n = 1).


* **Parameters**

    **interplanar_spacings** (*dict*) – dictionary of hkl to interplanar spacing



* **Returns**

    dict of hkl plane (3-tuple) to Bragg angle in radians (float)



#### cell_intensity(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates cell intensity for each hkl plane. For simplicity’s sake, take I =

```
|
```

F|\*\*2.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict of hkl plane to cell intensity



#### cell_scattering_factors(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates the scattering factor for the whole cell.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict of hkl plane (3-tuple) to scattering factor (in angstroms).



#### electron_scattering_factors(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates atomic scattering factors for electrons using the Mott-Bethe formula (1st order Born approximation).


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict from atomic symbol to another dict of hkl plane to factor (in angstroms)



#### _static_ generate_points(coord_left: int = -10, coord_right: int = 10)
Generates a bunch of 3D points that span a cube.


* **Parameters**


    * **coord_left** (*int*) – The minimum coordinate value.


    * **coord_right** (*int*) – The maximum coordinate value.



* **Returns**

    Numpy 2d array



#### get_first_point(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), points: list)
Gets the first point to be plotted in the 2D DP, corresponding to maximum d/minimum R.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **points** (*list*) – All points to be checked.



* **Returns**

    dict of a hkl plane to max interplanar distance.



#### _static_ get_interplanar_angle(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), p1: tuple[int, int, int], p2: tuple[int, int, int])
Returns the interplanar angle (in degrees) between the normal of two crystal planes.
Formulas from International Tables for Crystallography Volume C pp. 2-9.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **p1** (*3-tuple*) – plane 1


    * **p2** (*3-tuple*) – plane 2



* **Returns**

    float



#### get_interplanar_spacings(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), points: list[tuple[int, int, int]] | np.ndarray)

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – the input structure.


    * **points** (*tuple*) – the desired hkl indices.



* **Returns**

    Dict of hkl to its interplanar spacing, in angstroms (float).



#### get_pattern(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), scaled: bool | None = None, two_theta_range: tuple[float, float] | None = None)
Returns all relevant TEM DP info in a pandas dataframe.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **scaled** (*bool*) – Required value for inheritance, does nothing in TEM pattern


    * **two_theta_range** (*Tuple*) – Required value for inheritance, does nothing in TEM pattern



* **Returns**

    PandasDataFrame



#### get_plot_2d(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Generates the 2D diffraction pattern of the input structure.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.



* **Returns**

    Figure



#### get_plot_2d_concise(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure))
Generates the concise 2D diffraction pattern of the input structure of a smaller size and without layout.
Does not display.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.



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



#### get_positions(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), points: list)
Calculates all the positions of each hkl point in the 2D diffraction pattern by vector addition.
Distance in centimeters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


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




#### is_parallel(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), plane: tuple[int, int, int], other_plane: tuple[int, int, int])
Checks if two hkl planes are parallel in reciprocal space.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **plane** (*3-tuple*) – The first plane to be compared.


    * **other_plane** (*3-tuple*) – The other plane to be compared.



* **Returns**

    boolean



#### normalized_cell_intensity(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Normalizes the cell_intensity dict to 1, for use in plotting.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


    * **bragg_angles** (*dict** of **3-tuple to float*) – The Bragg angles for each hkl plane.



* **Returns**

    dict of hkl plane to normalized cell intensity



#### tem_dots(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), points)
Generates all TEM_dot as named tuples that will appear on the 2D diffraction pattern.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


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



#### x_ray_factors(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bragg_angles: dict[tuple[int, int, int], float])
Calculates x-ray factors, which are required to calculate atomic scattering factors. Method partially inspired
by the equivalent process in the xrd module.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – The input structure.


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