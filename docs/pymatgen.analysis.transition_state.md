---
layout: default
title: pymatgen.analysis.transition_state.md
nav_exclude: true
---

# pymatgen.analysis.transition_state module

Some reimplementation of Henkelman’s Transition State Analysis utilities,
which are originally in Perl. Additional features beyond those offered by
Henkelman’s utilities will be added.

This allows the usage and customization in Python.


### _class_ pymatgen.analysis.transition_state.NEBAnalysis(r, energies, forces, structures, spline_options=None)
Bases: `MSONable`

An NEBAnalysis class.

Initializes an NEBAnalysis from the cumulative root mean squared distances
between structures, the energies, the forces, the structures and the
interpolation_order for the analysis.


* **Parameters**


    * **r** – Root mean square distances between structures


    * **energies** – Energies of each structure along reaction coordinate


    * **forces** – Tangent forces along the reaction coordinate.


    * **structures** (*[*[*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)*]*) – List of Structures along reaction
    coordinate.


    * **spline_options** (*dict*) – Options for cubic spline. For example,
    {“saddle_point”: “zero_slope”} forces the slope at the saddle to
    be zero.



#### as_dict()
Dict representation of NEBAnalysis.


* **Returns**

    JSON-serializable dict representation.



#### _classmethod_ from_dir(root_dir, relaxation_dirs=None, \*\*kwargs)
Initializes a NEBAnalysis object from a directory of a NEB run.
Note that OUTCARs must be present in all image directories. For the
terminal OUTCARs from relaxation calculations, you can specify the
locations using relaxation_dir. If these are not specified, the code
will attempt to look for the OUTCARs in 00 and 0n directories,
followed by subdirs “start”, “end” or “initial”, “final” in the
root_dir. These are just some typical conventions used
preferentially in Shyue Ping’s MAVRL research group. For the
non-terminal points, the CONTCAR is read to obtain structures. For
terminal points, the POSCAR is used. The image directories are
assumed to be the only directories that can be resolved to integers.
E.g., “00”, “01”, “02”, “03”, “04”, “05”, “06”. The minimum
sub-directory structure that can be parsed is of the following form (
a 5-image example is shown):

00:
- POSCAR
- OUTCAR
01, 02, 03, 04, 05:
- CONTCAR
- OUTCAR
06:
- POSCAR
- OUTCAR


* **Parameters**


    * **root_dir** (*str*) – Path to the root directory of the NEB calculation.


    * **relaxation_dirs** (*tuple*) – This specifies the starting and ending
    relaxation directories from which the OUTCARs are read for the
    terminal points for the energies.



* **Returns**

    NEBAnalysis object.



#### _classmethod_ from_outcars(outcars, structures, \*\*kwargs)
Initializes an NEBAnalysis from Outcar and Structure objects. Use
the static constructors, e.g., `from_dir` instead if you
prefer to have these automatically generated from a directory of NEB
calculations.


* **Parameters**


    * **outcars** (*[*[*Outcar*](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar)*]*) – List of Outcar objects. Note that these have
    to be ordered from start to end along reaction coordinates.


    * **structures** (*[*[*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)*]*) – List of Structures along reaction
    coordinate. Must be same length as outcar.


    * **interpolation_order** (*int*) – Order of polynomial to use to
    interpolate between images. Same format as order parameter in
    scipy.interplotate.PiecewisePolynomial.



#### get_extrema(normalize_rxn_coordinate=True)
Returns the positions of the extrema along the MEP. Both local
minimums and maximums are returned.


* **Parameters**

    **normalize_rxn_coordinate** (*bool*) – Whether to normalize the
    reaction coordinate to between 0 and 1. Defaults to True.



* **Returns**

    (min_extrema, max_extrema), where the extrema are given as
    [(x1, y1), (x2, y2), …].



#### get_plot(normalize_rxn_coordinate=True, label_barrier=True)
Returns the NEB plot. Uses Henkelman’s approach of spline fitting
each section of the reaction path based on tangent force and energies.


* **Parameters**


    * **normalize_rxn_coordinate** (*bool*) – Whether to normalize the
    reaction coordinate to between 0 and 1. Defaults to True.


    * **label_barrier** (*bool*) – Whether to label the maximum barrier.



* **Returns**

    matplotlib.pyplot object.



#### setup_spline(spline_options=None)
Setup of the options for the spline interpolation.


* **Parameters**

    **spline_options** (*dict*) – Options for cubic spline. For example,
    {“saddle_point”: “zero_slope”} forces the slope at the saddle to
    be zero.



### pymatgen.analysis.transition_state.combine_neb_plots(neb_analyses, arranged_neb_analyses=False, reverse_plot=False)
neb_analyses: a list of NEBAnalysis objects.

arranged_neb_analyses: The code connects two end points with the
smallest-energy difference. If all end points have very close energies, it’s
likely to result in an inaccurate connection. Manually arrange neb_analyses
if the combined plot is not as expected compared with all individual plots.
E.g., if there are two NEBAnalysis objects to combine, arrange in such a
way that the end-point energy of the first NEBAnalysis object is the
start-point energy of the second NEBAnalysis object.
Note that the barrier labeled in y-axis in the combined plot might be
different from that in the individual plot due to the reference energy used.
reverse_plot: reverse the plot or percolation direction.
return: a NEBAnalysis object