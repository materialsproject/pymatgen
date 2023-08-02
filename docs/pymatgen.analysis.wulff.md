---
layout: default
title: pymatgen.analysis.wulff.md
nav_exclude: true
---

# pymatgen.analysis.wulff module

This module define a WulffShape class to generate the Wulff shape from
a lattice, a list of indices and their corresponding surface energies,
and the total area and volume of the Wulff shape, the weighted surface energy,
the anisotropy and shape_factor can also be calculated.
In support of plotting from a given view in terms of miller index.

The lattice is from the conventional unit cell, and (hkil) for hexagonal
lattices.

If you use this code extensively, consider citing the following:

Tran, R.; Xu, Z.; Radhakrishnan, B.; Winston, D.; Persson, K. A.; Ong, S. P.
(2016). Surface energies of elemental crystals. Scientific Data.


### _class_ pymatgen.analysis.wulff.WulffFacet(normal, e_surf, normal_pt, dual_pt, index, m_ind_orig, miller)
Bases: `object`

Helper container for each Wulff plane.


* **Parameters**


    * **normal** –


    * **e_surf** –


    * **normal_pt** –


    * **dual_pt** –


    * **index** –


    * **m_ind_orig** –


    * **miller** –



### _class_ pymatgen.analysis.wulff.WulffShape(lattice: [Lattice](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice), miller_list, e_surf_list, symprec=1e-05)
Bases: `object`

Generate Wulff Shape from list of miller index and surface energies,
with given conventional unit cell.
surface energy (Jm^2) is the length of normal.

Wulff shape is the convex hull.
Based on:
[http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html](http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html)

Process:


    1. get wulff simplices


    2. label with color


    3. get wulff_area and other properties


#### debug(bool)

#### alpha()

#### transparency()

#### color_set()

#### grid_off(bool)

#### axis_off(bool)

#### show_area()

#### off_color()

### color of facets off wulff()

#### structure()

### Structure object, input conventional unit cell (with H ) from lattice()

#### miller_list()

### list of input miller index, for hcp in the form of hkil()

#### hkl_list()

### modify hkill to hkl, in the same order with input_miller()

#### e_surf_list()

### list of input surface energies, in the same order with input_miller()

#### lattice()

### Lattice object, the input lattice for the conventional unit cell()

#### facets()

### [WulffFacet] for all facets considering symm()

#### dual_cv_simp()

### simplices from the dual convex hull (dual_pt)()

#### wulff_pt_list()

#### wulff_cv_simp()

### simplices from the convex hull of wulff_pt_list()

#### on_wulff()

### list for all input_miller, True is on wulff.()

#### color_area()

### list for all input_miller, total area on wulff, off_wulff = 0.()

#### miller_area()

### ($hkl$): area for all input_miller()

* **Parameters**


    * **lattice** – Lattice object of the conventional unit cell


    * **miller_list** (*[**(**hkl*) – list of hkl or hkil for hcp


    * **e_surf_list** (*[**float**]*) – list of corresponding surface energies


    * **symprec** (*float*) – for recp_operation, default is 1e-5.



#### _property_ anisotropy()
Returns:
(float) Coefficient of Variation from weighted surface energy
The ideal sphere is 0.


#### _property_ area_fraction_dict()
Returns:
(dict): {hkl: area_hkl/total area on wulff}.


#### _property_ effective_radius()
Radius of the Wulffshape when the
Wulffshape is approximated as a sphere.


* **Returns**

    (float) radius.



#### get_line_in_facet(facet)
Returns the sorted pts in a facet used to draw a line.


#### get_plot(color_set='PuBu', grid_off=True, axis_off=True, show_area=False, alpha=1, off_color='red', direction=None, bar_pos=(0.75, 0.15, 0.05, 0.65), bar_on=False, units_in_JPERM2=True, legend_on=True, aspect_ratio=(8, 8), custom_colors=None)
Get the Wulff shape plot.


* **Parameters**


    * **color_set** – default is ‘PuBu’


    * **grid_off** (*bool*) – default is True


    * **axis_off** (*bool*) – default is True


    * **show_area** (*bool*) – default is False


    * **alpha** (*float*) – chosen from 0 to 1 (float), default is 1


    * **off_color** – Default color for facets not present on the Wulff shape.


    * **direction** – default is (1, 1, 1)


    * **bar_pos** – default is [0.75, 0.15, 0.05, 0.65]


    * **bar_on** (*bool*) – default is False


    * **legend_on** (*bool*) – default is True


    * **aspect_ratio** – default is (8, 8)


    * **(****{****(****h** (*custom_colors*) – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **k** – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **l}** – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **units_in_JPERM2** (*bool*) – Units of surface energy, defaults to
    Joules per square meter (True)



* **Returns**

    (matplotlib.pyplot)



#### get_plotly(color_set='PuBu', off_color='red', alpha=1, custom_colors=None, units_in_JPERM2=True)
Get the Wulff shape as a plotly Figure object.


* **Parameters**


    * **color_set** – default is ‘PuBu’


    * **alpha** (*float*) – chosen from 0 to 1 (float), default is 1


    * **off_color** – Default color for facets not present on the Wulff shape.


    * **(****{****(****h** (*custom_colors*) – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **k** – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **l}** – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **units_in_JPERM2** (*bool*) – Units of surface energy, defaults to
    Joules per square meter (True)



* **Returns**

    (plotly.graph_objs.Figure)



#### _property_ miller_area_dict()
area_hkl on wulff}.


* **Type**

    Returns {hkl



#### _property_ miller_energy_dict()
surface energy_hkl}.


* **Type**

    Returns {hkl



#### _property_ shape_factor()
This is useful for determining the critical nucleus size.
A large shape factor indicates great anisotropy.
See Ballufi, R. W., Allen, S. M. & Carter, W. C. Kinetics

> of Materials. (John Wiley & Sons, 2005), p.461.


* **Returns**

    (float) Shape factor.



#### show(\*args, \*\*kwargs)
Show the Wulff plot.


* **Parameters**


    * **\*args** – Passed to get_plot.


    * **\*\*kwargs** – Passed to get_plot.



#### _property_ surface_area()
Total surface area of Wulff shape.


#### _property_ tot_corner_sites()
Returns the number of vertices in the convex hull.
Useful for identifying catalytically active sites.


#### _property_ tot_edges()
Returns the number of edges in the convex hull.
Useful for identifying catalytically active sites.


#### _property_ total_surface_energy()
Total surface energy of the Wulff shape.


* **Returns**

    (float) sum(surface_energy_hkl \* area_hkl)



#### _property_ volume()
Volume of the Wulff shape.


#### _property_ weighted_surface_energy()
Returns:
sum(surface_energy_hkl \* area_hkl)/ sum(area_hkl).


### pymatgen.analysis.wulff.get_tri_area(pts)
Given a list of coords for 3 points,
Compute the area of this triangle.


* **Parameters**

    **pts** – [a, b, c] three points



### pymatgen.analysis.wulff.hkl_tuple_to_str(hkl)
Prepare for display on plots “(hkl)” for surfaces


* **Parameters**

    **hkl** – in the form of [h, k, l] or (h, k, l).