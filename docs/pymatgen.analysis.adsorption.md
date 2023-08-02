---
layout: default
title: pymatgen.analysis.adsorption.md
nav_exclude: true
---

# pymatgen.analysis.adsorption module

This module provides classes used to enumerate surface sites and to find
adsorption sites on slabs.


### _class_ pymatgen.analysis.adsorption.AdsorbateSiteFinder(slab, selective_dynamics: bool = False, height: float = 0.9, mi_vec: ArrayLike | None = None)
Bases: `object`

This class finds adsorbate sites on slabs and generates adsorbate
structures according to user-defined criteria.

The algorithm for finding sites is essentially as follows:


    1. Determine “surface sites” by finding those within

        a height threshold along the miller index of the
        highest site


    2. Create a network of surface sites using the Delaunay

        triangulation of the surface sites


    3. Assign on-top, bridge, and hollow adsorption sites

        at the nodes, edges, and face centers of the Del.
        Triangulation


    4. Generate structures from a molecule positioned at

        these sites

Create an AdsorbateSiteFinder object.


* **Parameters**


    * **slab** ([*Slab*](pymatgen.core.surface.md#pymatgen.core.surface.Slab)) – slab object for which to find adsorbate sites


    * **selective_dynamics** (*bool*) – flag for whether to assign
    non-surface sites as fixed for selective dynamics


    * **height** (*float*) – height criteria for selection of surface sites


    * **mi_vec** (*3-D array-like*) – vector corresponding to the vector
    concurrent with the miller index, this enables use with
    slabs that have been reoriented, but the miller vector
    must be supplied manually



#### add_adsorbate(molecule, ads_coord, repeat=None, translate=True, reorient=True)
Adds an adsorbate at a particular coordinate. Adsorbate represented
by a Molecule object and is translated to (0, 0, 0) if translate is
True, or positioned relative to the input adsorbate coordinate if
translate is False.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – molecule object representing the adsorbate


    * **ads_coord** (*array*) – coordinate of adsorbate position


    * **repeat** (*3-tuple** or **list*) – input for making a supercell of slab
    prior to placing the adsorbate


    * **translate** (*bool*) – flag on whether to translate the molecule so
    that its CoM is at the origin prior to adding it to the surface


    * **reorient** (*bool*) – flag on whether to reorient the molecule to
    have its z-axis concurrent with miller index



#### adsorb_both_surfaces(molecule, repeat=None, min_lw=5.0, translate=True, reorient=True, find_args=None)
Function that generates all adsorption structures for a given
molecular adsorbate on both surfaces of a slab. This is useful for
calculating surface energy where both surfaces need to be equivalent or
if we want to calculate nonpolar systems.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – molecule corresponding to adsorbate


    * **repeat** (*3-tuple** or **list*) – repeat argument for supercell generation


    * **min_lw** (*float*) – minimum length and width of the slab, only used
    if repeat is None


    * **reorient** (*bool*) – flag on whether or not to reorient adsorbate
    along the miller index


    * **find_args** (*dict*) – dictionary of arguments to be passed to the
    call to self.find_adsorption_sites, e.g. {“distance”:2.0}



#### _classmethod_ assign_selective_dynamics(slab)
Helper function to assign selective dynamics site_properties based
on surface, subsurface site properties.


* **Parameters**

    **slab** ([*Slab*](pymatgen.core.surface.md#pymatgen.core.surface.Slab)) – slab for which to assign selective dynamics



#### assign_site_properties(slab, height=0.9)
Assigns site properties.


#### _classmethod_ ensemble_center(site_list, indices, cartesian=True)
Finds the center of an ensemble of sites selected from a list of
sites. Helper method for the find_adsorption_sites algorithm.


* **Parameters**


    * **site_list** (*list** of **sites*) – list of sites


    * **indices** (*list** of **ints*) – list of ints from which to select
    sites from site list


    * **cartesian** (*bool*) – whether to get average fractional or
    Cartesian coordinate



#### find_adsorption_sites(distance=2.0, put_inside=True, symm_reduce=0.01, near_reduce=0.01, positions=('ontop', 'bridge', 'hollow'), no_obtuse_hollow=True)
Finds surface sites according to the above algorithm. Returns a list
of corresponding Cartesian coordinates.


* **Parameters**


    * **distance** (*float*) – distance from the coordinating ensemble
    of atoms along the miller index for the site (i. e.
    the distance from the slab itself)


    * **put_inside** (*bool*) – whether to put the site inside the cell


    * **symm_reduce** (*float*) – symm reduction threshold


    * **near_reduce** (*float*) – near reduction threshold


    * **positions** (*list*) – which positions to include in the site finding
    “ontop”: sites on top of surface sites
    “bridge”: sites at edges between surface sites in Delaunay

    > triangulation of surface sites in the miller plane

    ”hollow”: sites at centers of Delaunay triangulation faces
    “subsurface”: subsurface positions projected into miller plane



    * **no_obtuse_hollow** (*bool*) – flag to indicate whether to include
    obtuse triangular ensembles in hollow sites



#### find_surface_sites_by_height(slab, height=0.9, xy_tol=0.05)
This method finds surface sites by determining which sites are
within a threshold value in height from the topmost site in a list of
sites.


* **Parameters**


    * **site_list** (*list*) – list of sites from which to select surface sites


    * **height** (*float*) – threshold in angstroms of distance from topmost
    site in slab along the slab c-vector to include in surface
    site determination


    * **xy_tol** (*float*) – if supplied, will remove any sites which are
    within a certain distance in the miller plane.



* **Returns**

    list of sites selected to be within a threshold of the highest



#### _classmethod_ from_bulk_and_miller(structure, miller_index, min_slab_size=8.0, min_vacuum_size=10.0, max_normal_search=None, center_slab=True, selective_dynamics=False, undercoord_threshold=0.09)
This method constructs the adsorbate site finder from a bulk
structure and a miller index, which allows the surface sites to be
determined from the difference in bulk and slab coordination, as
opposed to the height threshold.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure from which slab
    input to the ASF is constructed


    * **miller_index** (*3-tuple** or **list*) – miller index to be used


    * **min_slab_size** (*float*) – min slab size for slab generation


    * **min_vacuum_size** (*float*) – min vacuum size for slab generation


    * **max_normal_search** (*int*) – max normal search for slab generation


    * **center_slab** (*bool*) – whether to center slab in slab generation


    * **dynamics** (*selective*) – whether to assign surface sites
    to selective dynamics


    * **undercoord_threshold** (*float*) – threshold of “undercoordation”
    to use for the assignment of surface sites. Default is
    0.1, for which surface sites will be designated if they
    are 10% less coordinated than their bulk counterpart



#### generate_adsorption_structures(molecule, repeat=None, min_lw=5.0, translate=True, reorient=True, find_args=None)
Function that generates all adsorption structures for a given
molecular adsorbate. Can take repeat argument or minimum length/width
of precursor slab as an input.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – molecule corresponding to adsorbate


    * **repeat** (*3-tuple** or **list*) – repeat argument for supercell generation


    * **min_lw** (*float*) – minimum length and width of the slab, only used
    if repeat is None


    * **translate** (*bool*) – flag on whether to translate the molecule so
    that its CoM is at the origin prior to adding it to the surface


    * **reorient** (*bool*) – flag on whether or not to reorient adsorbate
    along the miller index


    * **find_args** (*dict*) – dictionary of arguments to be passed to the
    call to self.find_adsorption_sites, e.g. {“distance”:2.0}



#### generate_substitution_structures(atom, target_species=None, sub_both_sides=False, range_tol=0.01, dist_from_surf=0)
Function that performs substitution-type doping on the surface and
returns all possible configurations where one dopant is substituted per
surface. Can substitute one surface or both.


* **Parameters**


    * **atom** (*str*) – atom corresponding to substitutional dopant


    * **sub_both_sides** (*bool*) – If true, substitute an equivalent
    site on the other surface


    * **target_species** (*list*) – List of specific species to substitute


    * **range_tol** (*float*) – Find viable substitution sites at a specific
    distance from the surface +- this tolerance


    * **dist_from_surf** (*float*) – Distance from the surface to find viable
    substitution sites, defaults to 0 to substitute at the surface



#### get_extended_surface_mesh(repeat=(5, 5, 1))
Gets an extended surface mesh for to use for adsorption site finding
by constructing supercell of surface sites.


* **Parameters**

    **repeat** (*3-tuple*) – repeat for getting extended surface mesh



#### near_reduce(coords_set, threshold=0.0001)
Prunes coordinate set for coordinates that are within threshold.


* **Parameters**


    * **coords_set** (*Nx3 array-like*) – list or array of coordinates


    * **threshold** (*float*) – threshold value for distance



#### subsurface_sites()
Convenience method to return list of subsurface sites.


#### _property_ surface_sites()
Convenience method to return a list of surface sites.


#### symm_reduce(coords_set, threshold=1e-06)
Reduces the set of adsorbate sites by finding removing symmetrically
equivalent duplicates.


* **Parameters**


    * **coords_set** – coordinate set in Cartesian coordinates


    * **threshold** – tolerance for distance equivalence, used
    as input to in_coord_list_pbc for dupl. checking



### pymatgen.analysis.adsorption.get_mi_vec(slab)
Convenience function which returns the unit vector aligned with the
miller index.


### pymatgen.analysis.adsorption.get_rot(slab)
Gets the transformation to rotate the z axis into the miller index.


### pymatgen.analysis.adsorption.plot_slab(slab, ax, scale=0.8, repeat=5, window=1.5, draw_unit_cell=True, decay=0.2, adsorption_sites=True, inverse=False)
Function that helps visualize the slab in a 2-D plot, for convenient
viewing of output of AdsorbateSiteFinder.


* **Parameters**


    * **slab** (*slab*) – Slab object to be visualized


    * **ax** (*axes*) – matplotlib axes with which to visualize


    * **scale** (*float*) – radius scaling for sites


    * **repeat** (*int*) – number of repeating unit cells to visualize


    * **window** (*float*) – window for setting the axes limits, is essentially
    a fraction of the unit cell limits


    * **draw_unit_cell** (*bool*) – flag indicating whether or not to draw cell


    * **decay** (*float*) – how the alpha-value decays along the z-axis


    * **inverse** (*bool*) – invert z axis to plot opposite surface



### pymatgen.analysis.adsorption.put_coord_inside(lattice, cart_coordinate)
Converts a Cartesian coordinate such that it is inside the unit cell.


### pymatgen.analysis.adsorption.reorient_z(structure)
Reorients a structure such that the z axis is concurrent with the normal
to the A-B plane.