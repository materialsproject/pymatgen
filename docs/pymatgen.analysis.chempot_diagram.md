---
layout: default
title: pymatgen.analysis.chempot_diagram.md
nav_exclude: true
---

# pymatgen.analysis.chempot_diagram module

This module implements the construction and plotting of chemical potential diagrams
from a list of entries within a chemical system containing 2 or more elements. The
chemical potential diagram is the mathematical dual to the traditional compositional
phase diagram.

For more information, please cite/reference the paper below:

> Todd, P. K., McDermott, M. J., Rom, C. L., Corrao, A. A., Denney, J. J., Dwaraknath,
> S. S.,  Khalifah, P. G., Persson, K. A., & Neilson, J. R. (2021). Selectivity in
> Yttrium Manganese Oxide Synthesis via Local Chemical Potentials in Hyperdimensional
> Phase Space. Journal of the American Chemical Society, 143(37), 15185-15194.
> [https://doi.org/10.1021/jacs.1c06229](https://doi.org/10.1021/jacs.1c06229)

Please also consider referencing the original 1999 paper by H. Yokokawa,
who outlined many of its possible uses:

> Yokokawa, H. “Generalized chemical potential diagram and its applications to
> chemical reactions at interfaces between dissimilar materials.” JPE 20,
> 258 (1999). [https://doi.org/10.1361/105497199770335794](https://doi.org/10.1361/105497199770335794)


### _class_ pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram(entries: list[[PDEntry](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)], limits: dict[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element), tuple[float, float]] | None = None, default_min_limit: float = -50.0, formal_chempots: bool = True)
Bases: `MSONable`

The chemical potential diagram is the mathematical dual to the compositional
phase diagram. To create the diagram, convex minimization is
performed in energy (E) vs. chemical potential (μ) space by taking the lower convex
envelope of hyperplanes. Accordingly, “points” on the compositional phase diagram
become N-dimensional convex polytopes (domains) in chemical potential space.

For more information on this specific implementation of the algorithm,
please cite/reference the paper below:

> Todd, P. K., McDermott, M. J., Rom, C. L., Corrao, A. A., Denney, J. J., Dwaraknath,
> S. S.,  Khalifah, P. G., Persson, K. A., & Neilson, J. R. (2021). Selectivity in
> Yttrium Manganese Oxide Synthesis via Local Chemical Potentials in Hyperdimensional
> Phase Space. Journal of the American Chemical Society, 143(37), 15185-15194.
> [https://doi.org/10.1021/jacs.1c06229](https://doi.org/10.1021/jacs.1c06229)


* **Parameters**


    * **entries** (*list**[*[*PDEntry*](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)*]*) – PDEntry-like objects containing a composition and
    energy. Must contain elemental references and be suitable for typical
    phase diagram construction. Entries must be within a chemical system
    of with 2+ elements.


    * **limits** (*dict**[*[*Element*](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)*, **float**] **| **None*) – Bounds of elemental chemical potentials (min, max),
    which are used to construct the border hyperplanes used in the HalfSpaceIntersection
    algorithm; these constrain the space over which the domains are calculated and also
    determine the size of the plotted diagram. Any elemental limits not specified are
    covered in the default_min_limit argument. e.g., {Element(“Li”): [-12.0, 0.0], …}


    * **default_min_limit** (*float*) – Default minimum chemical potential limit (i.e.,
    lower bound) for unspecified elements within the “limits” argument.


    * **formal_chempots** (*bool*) – Whether to plot the formal (‘reference’) chemical potentials
    (i.e. μ_X - μ_X^0) or the absolute DFT reference energies (i.e. μ_X(DFT)).
    Default is True (i.e. plot formal chemical potentials).



#### _property_ border_hyperplanes(_: ndarra_ )
Returns bordering hyperplanes.


#### _property_ chemical_system(_: st_ )
Returns the chemical system (A-B-C-…) of diagram object.


#### _property_ domains(_: dict[str, numpy.ndarray_ )
Mapping of formulas to array of domain boundary points.


#### _property_ el_refs(_: dict[[pymatgen.core.periodic_table.Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element), [pymatgen.analysis.phase_diagram.PDEntry](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)_ )
Returns a dictionary of elements and reference entries.


#### _property_ entry_dict(_: dict[str, [ComputedEntry](pymatgen.entries.computed_entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )
Mapping between reduced formula and ComputedEntry.


#### get_plot(elements: list[[Element](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element) | str] | None = None, label_stable: bool | None = True, formulas_to_draw: list[str] | None = None, draw_formula_meshes: bool | None = True, draw_formula_lines: bool | None = True, formula_colors: list[str] = ['rgb(27,158,119)', 'rgb(217,95,2)', 'rgb(117,112,179)', 'rgb(231,41,138)', 'rgb(102,166,30)', 'rgb(230,171,2)', 'rgb(166,118,29)', 'rgb(102,102,102)'], element_padding: float | None = 1.0)
Plot the 2-dimensional or 3-dimensional chemical potential diagram using an
interactive Plotly interface.

Elemental axes can be specified; if none provided, will automatically default
to first 2-3 elements within the “elements” attribute.

In 3D, this method also allows for plotting of lower-dimensional “slices” of
hyperdimensional polytopes (e.g., the LiMnO2 domain within a Y-Mn-O diagram).
This allows for visualization of some of the phase boundaries that can only
be seen fully in high dimensional space; see the “formulas_to_draw” argument.


* **Parameters**


    * **elements** – list of elements to use as axes in the diagram. If None,
    automatically defaults to the first 2 or elements within the
    object’s “elements” attribute.


    * **label_stable** – whether to label stable phases by their reduced
    formulas. Defaults to True.


    * **formulas_to_draw** – for 3-dimensional diagrams, an optional list of
    formulas to plot on the diagram; if these are from a different
    chemical system a 3-d polyhedron “slice” will be plotted. Defaults to None.


    * **draw_formula_meshes** – whether to draw a colored mesh for the
    optionally specified formulas_to_draw. Defaults to True.


    * **draw_formula_lines** – whether to draw bounding lines for the
    optionally specified formulas_to_draw. Defaults to True.


    * **formula_colors** – a list of colors to use in the plotting of the optionally
    specified formulas_to-draw. Defaults to the Plotly Dark2 color scheme.


    * **element_padding** – if provided, automatically adjusts chemical potential axis
    limits of the plot such that elemental domains have the specified padding
    (in eV/atom), helping provide visual clarity. Defaults to 1.0.



* **Returns**

    A Plotly Figure object



#### _property_ hyperplane_entries(_: list[[pymatgen.analysis.phase_diagram.PDEntry](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)_ )
Returns list of entries corresponding to hyperplanes.


#### _property_ hyperplanes(_: ndarra_ )
Returns array of hyperplane data.


#### _property_ lims(_: ndarra_ )
Returns array of limits used in constructing hyperplanes.


### pymatgen.analysis.chempot_diagram.get_2d_orthonormal_vector(line_pts: ndarray)
Calculates a vector that is orthonormal to a line given by a set of points. Used
for determining the location of an annotation on a 2-d chemical potential diagram.


* **Parameters**

    **line_pts** – a 2x2 array in the form of [[x0, y0], [x1, y1]] giving the
    coordinates of a line



* **Returns**

    A length-2 vector that is orthonormal to the line.



* **Return type**

    np.ndarray



### pymatgen.analysis.chempot_diagram.get_centroid_2d(vertices: ndarray)
A bare-bones implementation of the formula for calculating the centroid of a 2D
polygon. Useful for calculating the location of an annotation on a chemical
potential domain within a 3D chemical potential diagram.

**NOTE**: vertices must be ordered circumferentially!


* **Parameters**

    **vertices** – array of 2-d coordinates corresponding to a polygon, ordered
    circumferentially



* **Returns**

    Array giving 2-d centroid coordinates



### pymatgen.analysis.chempot_diagram.simple_pca(data: ndarray, k: int = 2)
A bare-bones implementation of principal component analysis (PCA) used in the
ChemicalPotentialDiagram class for plotting.


* **Parameters**


    * **data** – array of observations


    * **k** – Number of principal components returned



* **Returns**

    Tuple of projected data, eigenvalues, eigenvectors