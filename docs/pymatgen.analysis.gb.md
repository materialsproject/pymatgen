---
layout: default
title: pymatgen.analysis.gb.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.gb package

This package implements various grain boundary analyses.


## pymatgen.analysis.gb.grain module

Module containing classes to generate grain boundaries.


### _class_ GrainBoundary(lattice: np.ndarray | [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice), species: Sequence[CompositionLike], coords: Sequence[ArrayLike], rotation_axis: Vector3D, rotation_angle: float, gb_plane: Vector3D, join_plane: Vector3D, init_cell: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), vacuum_thickness: float, ab_shift: tuple[float, float], site_properties: dict[str, Any], oriented_unit_cell: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), validate_proximity: bool = False, coords_are_cartesian: bool = False, properties: dict | None = None)
Bases: [`Structure`](pymatgen.core.md#pymatgen.core.structure.Structure)

Subclass of Structure representing a GrainBoundary (GB) object. Implements additional
attributes pertaining to gbs, but the init method does not actually implement any
algorithm that creates a GB. This is a DUMMY class who’s init method only holds
information about the GB. Also has additional methods that returns other information
about a GB such as sigma value.

Note that all gbs have the GB surface normal oriented in the c-direction. This means
the lattice vectors a and b are in the GB surface plane (at least for one grain) and
the c vector is out of the surface plane (though not necessarily perpendicular to the
surface).

Makes a GB structure, a structure object with additional information
and methods pertaining to gbs.


* **Parameters**


    * **lattice** (*Lattice/3x3 array*) – The lattice, either as an instance or
    any 2D array. Each row should correspond to a lattice vector.


    * **species** (*[*[*Species*](pymatgen.core.md#pymatgen.core.periodic_table.Species)*]*) – Sequence of species on each site. Can take in
    flexible input, including:


        1. A sequence of element / species specified either as string
    symbols, e.g. [“Li”, “Fe2+”, “P”, …] or atomic numbers,
    e.g., (3, 56, …) or actual Element or Species objects.


        2. List of dict of elements/species and occupancies, e.g.,
    [{“Fe” : 0.5, “Mn”:0.5}, …]. This allows the setup of
    disordered structures.



    * **coords** (*Nx3 array*) – list of fractional/cartesian coordinates for each species.


    * **rotation_axis** (*list**[**int**]*) – Rotation axis of GB in the form of a list of integers, e.g. [1, 1, 0].


    * **rotation_angle** (*float**, **in unit** of **degree*) – rotation angle of GB.


    * **gb_plane** (*list*) – Grain boundary plane in the form of a list of integers
    e.g.: [1, 2, 3].


    * **join_plane** (*list*) – Joining plane of the second grain in the form of a list of
    integers. e.g.: [1, 2, 3].


    * **init_cell** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – initial bulk structure to form the GB.


    * **site_properties** (*dict*) – Properties associated with the sites as a
    dict of sequences, The sequences have to be the same length as
    the atomic species and fractional_coords. For GB, you should
    have the ‘grain_label’ properties to classify the sites as ‘top’,
    ‘bottom’, ‘top_incident’, or ‘bottom_incident’.


    * **vacuum_thickness** (*float in angstrom*) – The thickness of vacuum inserted
    between two grains of the GB.


    * **ab_shift** (*list** of **float**, **in unit** of **crystal vector a**, **b*) – The relative
    shift along a, b vectors.


    * **oriented_unit_cell** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – oriented unit cell of the bulk init_cell.
    Helps to accurately calculate the bulk properties that are consistent
    with GB calculations.


    * **validate_proximity** (*bool*) – Whether to check if there are sites
    that are less than 0.01 Ang apart. Defaults to False.


    * **coords_are_cartesian** (*bool*) – Set to True if you are providing
    coordinates in Cartesian coordinates. Defaults to False.


    * **properties** (*dict*) – dictionary containing properties associated
    with the whole GrainBoundary.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### as_dict()

* **Returns**

    Dictionary representation of GrainBoundary object.



#### _property_ bottom_grain(_: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure_ )
Return the bottom grain (Structure) of the GB.


#### _property_ coincidents(_: list[[Site](pymatgen.core.md#pymatgen.core.sites.Site)_ )
Return the a list of coincident sites.


#### copy()
Convenience method to get a copy of the structure, with options to add
site properties.


* **Returns**

    A copy of the Structure, with optionally new site_properties and
    optionally sanitized.



#### _classmethod_ from_dict(d)
Generates a GrainBoundary object from a dictionary created by as_dict().


* **Parameters**

    **d** – dict



* **Returns**

    GrainBoundary object



#### get_sorted_structure(key=None, reverse=False)
Get a sorted copy of the structure. The parameters have the same
meaning as in list.sort. By default, sites are sorted by the
electronegativity of the species. Note that Slab has to override this
because of the different __init__ args.


* **Parameters**


    * **key** – Specifies a function of one argument that is used to extract
    a comparison key from each list element: key=str.lower. The
    default value is None (compare the elements directly).


    * **reverse** (*bool*) – If set to True, then the list elements are sorted
    as if each comparison were reversed.



#### _property_ sigma(_: in_ )
This method returns the sigma value of the GB.
If using ‘quick_gen’ to generate GB, this value is not valid.


#### _property_ sigma_from_site_prop(_: in_ )
This method returns the sigma value of the GB from site properties.
If the GB structure merge some atoms due to the atoms too closer with
each other, this property will not work.


#### _property_ top_grain(_: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure_ )
Return the top grain (Structure) of the GB.


### _class_ GrainBoundaryGenerator(initial_structure, symprec: float = 0.1, angle_tolerance=1)
Bases: `object`

This class is to generate grain boundaries (GBs) from bulk
conventional cell (fcc, bcc can from the primitive cell), and works for Cubic,
Tetragonal, Orthorhombic, Rhombohedral, and Hexagonal systems.
It generate GBs from given parameters, which includes
GB plane, rotation axis, rotation angle.

This class works for any general GB, including twist, tilt and mixed GBs.
The three parameters, rotation axis, GB plane and rotation angle, are
sufficient to identify one unique GB. While sometimes, users may not be able
to tell what exactly rotation angle is but prefer to use sigma as an parameter,
this class also provides the function that is able to return all possible
rotation angles for a specific sigma value.
The same sigma value (with rotation axis fixed) can correspond to
multiple rotation angles.
Users can use structure matcher in pymatgen to get rid of the redundant structures.

initial_structure (Structure): Initial input structure. It can

    be conventional or primitive cell (primitive cell works for bcc and fcc).
    For fcc and bcc, using conventional cell can lead to a non-primitive
    grain boundary structure.
    This code supplies Cubic, Tetragonal, Orthorhombic, Rhombohedral, and
    Hexagonal systems.

symprec (float): Tolerance for symmetry finding. Defaults to 0.1 (the value used

    in Materials Project), which is for structures with slight deviations
    from their proper atomic positions (e.g., structures relaxed with
    electronic structure codes).
    A smaller value of 0.01 is often used for properly refined
    structures with atoms in the proper symmetry coordinates.
    User should make sure the symmetry is what you want.

angle_tolerance (float): Angle tolerance for symmetry finding.


#### _static_ enum_possible_plane_cubic(plane_cutoff, r_axis, r_angle)
Find all possible plane combinations for GBs given a rotation axis and angle for
cubic system, and classify them to different categories, including ‘Twist’,
‘Symmetric tilt’, ‘Normal tilt’, ‘Mixed’ GBs.


* **Parameters**


    * **plane_cutoff** (*int*) – the cutoff of plane miller index.


    * **r_axis** (*list** of **three integers**, **e.g. u**, **v**, **w*) – the rotation axis of the grain boundary, with the format of [u,v,w].


    * **r_angle** (*float*) – rotation angle of the GBs.



* **Returns**

    dictionary with keys as GB type, e.g. ‘Twist’,’Symmetric tilt’,etc.
    and values as the combination of the two plane miller index

    > (GB plane and joining plane).




* **Return type**

    all_combinations (dict)



#### _static_ enum_sigma_cubic(cutoff, r_axis)
Find all possible sigma values and corresponding rotation angles
within a sigma value cutoff with known rotation axis in cubic system.
The algorithm for this code is from reference, Acta Cryst, A40,108(1984).


* **Parameters**


    * **cutoff** (*int*) – the cutoff of sigma values.


    * **r_axis** (*list** of **three integers**, **e.g. u**, **v**, **w*) – the rotation axis of the grain boundary, with the format of [u,v,w].



* **Returns**

    dictionary with keys as the possible integer sigma values
    and values as list of the possible rotation angles to the
    corresponding sigma values.
    e.g. the format as
    {sigma1: [angle11,angle12,…], sigma2: [angle21, angle22,…],…}
    Note: the angles are the rotation angles of one grain respect to
    the other grain.
    When generate the microstructures of the grain boundary using these angles,
    you need to analyze the symmetry of the structure. Different angles may
    result in equivalent microstructures.



* **Return type**

    sigmas (dict)



#### _static_ enum_sigma_hex(cutoff, r_axis, c2_a2_ratio)
Find all possible sigma values and corresponding rotation angles
within a sigma value cutoff with known rotation axis in hexagonal system.
The algorithm for this code is from reference, Acta Cryst, A38,550(1982).


* **Parameters**


    * **cutoff** (*int*) – the cutoff of sigma values.


    * **integers** (*r_axis** (**list** of **three*) – or four integers, e.g. u, v, t, w):
    the rotation axis of the grain boundary.


    * **u** (*e.g.*) – or four integers, e.g. u, v, t, w):
    the rotation axis of the grain boundary.


    * **v** – or four integers, e.g. u, v, t, w):
    the rotation axis of the grain boundary.


    * **w** – or four integers, e.g. u, v, t, w):
    the rotation axis of the grain boundary.


    * **c2_a2_ratio** (*list** of **two integers**, **e.g. mu**, **mv*) – mu/mv is the square of the hexagonal axial ratio, which is rational
    number. If irrational, set c2_a2_ratio = None



* **Returns**

    dictionary with keys as the possible integer sigma values
    and values as list of the possible rotation angles to the
    corresponding sigma values.
    e.g. the format as
    {sigma1: [angle11,angle12,…], sigma2: [angle21, angle22,…],…}
    Note: the angles are the rotation angle of one grain respect to the
    other grain.
    When generate the microstructure of the grain boundary using these
    angles, you need to analyze the symmetry of the structure. Different
    angles may result in equivalent microstructures.



* **Return type**

    sigmas (dict)



#### _static_ enum_sigma_ort(cutoff, r_axis, c2_b2_a2_ratio)
Find all possible sigma values and corresponding rotation angles
within a sigma value cutoff with known rotation axis in orthorhombic system.
The algorithm for this code is from reference, Scipta Metallurgica 27, 291(1992).


* **Parameters**


    * **cutoff** (*int*) – the cutoff of sigma values.


    * **r_axis** (*list** of **three integers**, **e.g. u**, **v**, **w*) – the rotation axis of the grain boundary, with the format of [u,v,w].


    * **c2_b2_a2_ratio** (*list** of **three integers**, **e.g. mu**,**lambda**, **mv*) – mu:lam:mv is the square of the orthorhombic axial ratio with rational
    numbers. If irrational for one axis, set it to None.
    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.



* **Returns**

    dictionary with keys as the possible integer sigma values
    and values as list of the possible rotation angles to the
    corresponding sigma values.
    e.g. the format as
    {sigma1: [angle11,angle12,…], sigma2: [angle21, angle22,…],…}
    Note: the angles are the rotation angle of one grain respect to the
    other grain.
    When generate the microstructure of the grain boundary using these
    angles, you need to analyze the symmetry of the structure. Different
    angles may result in equivalent microstructures.



* **Return type**

    sigmas (dict)



#### _static_ enum_sigma_rho(cutoff, r_axis, ratio_alpha)
Find all possible sigma values and corresponding rotation angles
within a sigma value cutoff with known rotation axis in rhombohedral system.
The algorithm for this code is from reference, Acta Cryst, A45,505(1989).


* **Parameters**


    * **cutoff** (*int*) – the cutoff of sigma values.


    * **r_axis** (*list**[**int**]*) – of three integers, e.g. u, v, w
    or four integers, e.g. u, v, t, w):
    the rotation axis of the grain boundary, with the format of [u,v,w]
    or Weber indices [u, v, t, w].


    * **ratio_alpha** (*list** of **two integers**, **e.g. mu**, **mv*) – mu/mv is the ratio of (1+2\*cos(alpha))/cos(alpha) with rational number.
    If irrational, set ratio_alpha = None.



* **Returns**

    dictionary with keys as the possible integer sigma values
    and values as list of the possible rotation angles to the
    corresponding sigma values.
    e.g. the format as
    {sigma1: [angle11,angle12,…], sigma2: [angle21, angle22,…],…}
    Note: the angles are the rotation angle of one grain respect to the
    other grain.
    When generate the microstructure of the grain boundary using these
    angles, you need to analyze the symmetry of the structure. Different
    angles may result in equivalent microstructures.



* **Return type**

    sigmas (dict)



#### _static_ enum_sigma_tet(cutoff, r_axis, c2_a2_ratio)
Find all possible sigma values and corresponding rotation angles
within a sigma value cutoff with known rotation axis in tetragonal system.
The algorithm for this code is from reference, Acta Cryst, B46,117(1990).


* **Parameters**


    * **cutoff** (*int*) – the cutoff of sigma values.


    * **r_axis** (*list** of **three integers**, **e.g. u**, **v**, **w*) – the rotation axis of the grain boundary, with the format of [u,v,w].


    * **c2_a2_ratio** (*list** of **two integers**, **e.g. mu**, **mv*) – mu/mv is the square of the tetragonal axial ratio with rational number.
    if irrational, set c2_a2_ratio = None



* **Returns**

    dictionary with keys as the possible integer sigma values
    and values as list of the possible rotation angles to the
    corresponding sigma values.
    e.g. the format as
    {sigma1: [angle11,angle12,…], sigma2: [angle21, angle22,…],…}
    Note: the angles are the rotation angle of one grain respect to the
    other grain.
    When generate the microstructure of the grain boundary using these
    angles, you need to analyze the symmetry of the structure. Different
    angles may result in equivalent microstructures.



* **Return type**

    sigmas (dict)



#### gb_from_parameters(rotation_axis, rotation_angle, expand_times=4, vacuum_thickness=0.0, ab_shift: tuple[float, float] = (0, 0), normal=False, ratio=None, plane=None, max_search=20, tol_coi=1e-08, rm_ratio=0.7, quick_gen=False)

* **Parameters**


    * **rotation_axis** (*list*) – Rotation axis of GB in the form of a list of integer
    e.g.: [1, 1, 0]


    * **rotation_angle** (*float**, **in unit** of **degree*) – rotation angle used to generate GB.
    Make sure the angle is accurate enough. You can use the enum\* functions
    in this class to extract the accurate angle.
    e.g.: The rotation angle of sigma 3 twist GB with the rotation axis
    [1, 1, 1] and GB plane (1, 1, 1) can be 60 degree.
    If you do not know the rotation angle, but know the sigma value, we have
    provide the function get_rotation_angle_from_sigma which is able to return
    all the rotation angles of sigma value you provided.


    * **expand_times** (*int*) – The multiple times used to expand one unit grain to larger grain.
    This is used to tune the grain length of GB to warrant that the two GBs in one
    cell do not interact with each other. Default set to 4.


    * **vacuum_thickness** (*float**, **in angstrom*) – The thickness of vacuum that you want to insert
    between two grains of the GB. Default to 0.


    * **ab_shift** (*list** of **float**, **in unit** of **a**, **b vectors** of **Gb*) – in plane shift of two grains


    * **normal** (*logic*) – determine if need to require the c axis of top grain (first transformation matrix)
    perpendicular to the surface or not.
    default to false.


    * **ratio** (*list** of **integers*) – lattice axial ratio.
    For cubic system, ratio is not needed.
    For tetragonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.
    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
    For rhombohedral system, ratio = [mu, mv], list of two integers,
    that is, mu/mv is the ratio of (1+2\*cos(alpha))/cos(alpha).
    If irrational, set it to None.
    For hexagonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.
    This code also supplies a class method to generate the ratio from the
    structure (get_ratio). User can also make their own approximation and
    input the ratio directly.


    * **plane** (*list*) – Grain boundary plane in the form of a list of integers
    e.g.: [1, 2, 3]. If none, we set it as twist GB. The plane will be perpendicular
    to the rotation axis.


    * **max_search** (*int*) – max search for the GB lattice vectors that give the smallest GB
    lattice. If normal is true, also max search the GB c vector that perpendicular
    to the plane. For complex GB, if you want to speed up, you can reduce this value.
    But too small of this value may lead to error.


    * **tol_coi** (*float*) – tolerance to find the coincidence sites. When making approximations to
    the ratio needed to generate the GB, you probably need to increase this tolerance to
    obtain the correct number of coincidence sites. To check the number of coincidence
    sites are correct or not, you can compare the generated Gb object’s sigma_from_site_prop
    with enum\* sigma values (what user expected by input).


    * **rm_ratio** (*float*) – the criteria to remove the atoms which are too close with each other.
    rm_ratio\*bond_length of bulk system is the criteria of bond length, below which the atom
    will be removed. Default to 0.7.


    * **quick_gen** (*bool*) – whether to quickly generate a supercell, if set to true, no need to
    find the smallest cell.



* **Returns**

    Grain boundary structure (GB object).



#### get_ratio(max_denominator=5, index_none=None)
find the axial ratio needed for GB generator input.


* **Parameters**


    * **max_denominator** (*int*) – the maximum denominator for
    the computed ratio, default to be 5.


    * **index_none** (*int*) – specify the irrational axis.
    0-a, 1-b, 2-c. Only may be needed for orthorhombic system.



* **Returns**

    axial ratio needed for GB generator (list of integers).



#### _static_ get_rotation_angle_from_sigma(sigma, r_axis, lat_type='C', ratio=None)
Find all possible rotation angle for the given sigma value.


* **Parameters**


    * **sigma** (*int*) – sigma value provided


    * **integers** (*r_axis** (**list** of **three*) – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **u** (*e.g.*) – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **v** – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **w** – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **lat_type** (*one character*) – ‘c’ or ‘C’: cubic system

        ’t’ or ‘T’: tetragonal system
        ‘o’ or ‘O’: orthorhombic system
        ‘h’ or ‘H’: hexagonal system
        ‘r’ or ‘R’: rhombohedral system
        default to cubic system



    * **ratio** (*list** of **integers*) – lattice axial ratio.
    For cubic system, ratio is not needed.
    For tetragonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.
    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
    For rhombohedral system, ratio = [mu, mv], list of two integers,
    that is, mu/mv is the ratio of (1+2\*cos(alpha)/cos(alpha).
    If irrational, set it to None.
    For hexagonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.



* **Returns**

    rotation_angles corresponding to the provided sigma value.
    If the sigma value is not correct, return the rotation angle corresponding
    to the correct possible sigma value right smaller than the wrong sigma value provided.



#### _static_ get_trans_mat(r_axis, angle, normal=False, trans_cry=None, lat_type='c', ratio=None, surface=None, max_search=20, quick_gen=False)
Find the two transformation matrix for each grain from given rotation axis,
GB plane, rotation angle and corresponding ratio (see explanation for ratio
below).
The structure of each grain can be obtained by applying the corresponding
transformation matrix to the conventional cell.
The algorithm for this code is from reference, Acta Cryst, A32,783(1976).


* **Parameters**


    * **integers** (*surface** (**list** of **three*) – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **u** (*e.g.*) – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **v** – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **w** – or four integers, e.g. u, v, t, w for hex/rho system only):
    the rotation axis of the grain boundary.


    * **angle** (*float**, **in unit** of **degree*) – the rotation angle of the grain boundary


    * **normal** (*logic*) – determine if need to require the c axis of one grain associated with
    the first transformation matrix perpendicular to the surface or not.
    default to false.


    * **trans_cry** (*3 by 3 array*) – if the structure given are primitive cell in cubic system, e.g.
    bcc or fcc system, trans_cry is the transformation matrix from its
    conventional cell to the primitive cell.


    * **lat_type** (*one character*) – ‘c’ or ‘C’: cubic system

        ’t’ or ‘T’: tetragonal system
        ‘o’ or ‘O’: orthorhombic system
        ‘h’ or ‘H’: hexagonal system
        ‘r’ or ‘R’: rhombohedral system
        default to cubic system



    * **ratio** (*list** of **integers*) – lattice axial ratio.
    For cubic system, ratio is not needed.
    For tetragonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.
    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
    For rhombohedral system, ratio = [mu, mv], list of two integers,
    that is, mu/mv is the ratio of (1+2\*cos(alpha)/cos(alpha).
    If irrational, set it to None.
    For hexagonal system, ratio = [mu, mv], list of two integers,
    that is, mu/mv = c2/a2. If it is irrational, set it to none.


    * **integers** – or four integers, e.g. h, k, i, l for hex/rho system only):
    the miller index of grain boundary plane, with the format of [h,k,l]
    if surface is not given, the default is perpendicular to r_axis, which is
    a twist grain boundary.


    * **h** (*e.g.*) – or four integers, e.g. h, k, i, l for hex/rho system only):
    the miller index of grain boundary plane, with the format of [h,k,l]
    if surface is not given, the default is perpendicular to r_axis, which is
    a twist grain boundary.


    * **k** – or four integers, e.g. h, k, i, l for hex/rho system only):
    the miller index of grain boundary plane, with the format of [h,k,l]
    if surface is not given, the default is perpendicular to r_axis, which is
    a twist grain boundary.


    * **l** – or four integers, e.g. h, k, i, l for hex/rho system only):
    the miller index of grain boundary plane, with the format of [h,k,l]
    if surface is not given, the default is perpendicular to r_axis, which is
    a twist grain boundary.


    * **max_search** (*int*) – max search for the GB lattice vectors that give the smallest GB
    lattice. If normal is true, also max search the GB c vector that perpendicular
    to the plane.


    * **quick_gen** (*bool*) – whether to quickly generate a supercell, if set to true, no need to
    find the smallest cell.



* **Returns**

    The transformation array for one grain.
    t2 (3 by 3 integer array):

    > The transformation array for the other grain




* **Return type**

    t1 (3 by 3 integer array)



#### _static_ reduce_mat(mat, mag, r_matrix)
Reduce integer array mat’s determinant mag times by linear combination
of its row vectors, so that the new array after rotation (r_matrix) is
still an integer array.


* **Parameters**


    * **mat** (*3 by 3 array*) – input matrix


    * **mag** (*int*) – reduce times for the determinant


    * **r_matrix** (*3 by 3 array*) – rotation matrix



* **Returns**

    the reduced integer array



#### _static_ slab_from_csl(csl, surface, normal, trans_cry, max_search=20, quick_gen=False)
By linear operation of csl lattice vectors to get the best corresponding
slab lattice. That is the area of a,b vectors (within the surface plane)
is the smallest, the c vector first, has shortest length perpendicular
to surface [h,k,l], second, has shortest length itself.


* **Parameters**


    * **csl** (*3 by 3 integer array*) – input csl lattice.


    * **surface** (*list** of **three integers**, **e.g. h**, **k**, **l*) – the miller index of the surface, with the format of [h,k,l]


    * **normal** (*logic*) – determine if the c vector needs to perpendicular to surface


    * **trans_cry** (*3 by 3 array*) – transform matrix from crystal system to orthogonal system


    * **max_search** (*int*) – max search for the GB lattice vectors that give the smallest GB
    lattice. If normal is true, also max search the GB c vector that perpendicular
    to the plane.


    * **quick_gen** (*bool*) – whether to quickly generate a supercell, no need to find the smallest
    cell if set to true.



* **Returns**

    a slab lattice ( 3 by 3 integer array):



* **Return type**

    t_matrix



#### _static_ vec_to_surface(vec)
Transform a float vector to a surface miller index with integers.


* **Parameters**

    **vec** (*1 by 3 array float vector*) – input float vector



* **Returns**

    the surface miller index of the input vector.



### fix_pbc(structure, matrix=None)
Set all frac_coords of the input structure within [0,1].


* **Parameters**


    * **structure** (*pymatgen structure object*) – input structure


    * **matrix** (*lattice matrix**, **3 by 3 array/matrix*) – new structure’s lattice matrix,
    If None, use input structure’s matrix.



* **Returns**

    new structure with fixed frac_coords and lattice matrix



### symm_group_cubic(mat)
> obtain cubic symmetric equivalents of the list of vectors.


* **Parameters**

    **matrix** (*lattice matrix**, **n by 3 array/matrix*) –



* **Returns**

    cubic symmetric equivalents of the list of vectors.