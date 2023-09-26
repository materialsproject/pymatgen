---
layout: default
title: pymatgen.symmetry.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.symmetry package

The symmetry package implements symmetry tools like spacegroup determination, etc.


## pymatgen.symmetry.analyzer module

An interface to the excellent spglib library by Atsushi Togo
([http://spglib.sourceforge.net/](http://spglib.sourceforge.net/)) for pymatgen.

v1.0 - Now works with both ordered and disordered structure.
v2.0 - Updated for spglib 1.6.
v3.0 - pymatgen no longer ships with spglib. Instead, spglib (the python

> version) is now a dependency and the SpacegroupAnalyzer merely serves
> as an interface to spglib for pymatgen Structures.


### _class_ PointGroupAnalyzer(mol, tolerance=0.3, eigen_tolerance=0.01, matrix_tolerance=0.1)
Bases: `object`

A class to analyze the point group of a molecule.

The general outline of the algorithm is as follows:


1. Center the molecule around its center of mass.


2. Compute the inertia tensor and the eigenvalues and eigenvectors.


3. Handle the symmetry detection based on eigenvalues.

>
>     1. Linear molecules have one zero eigenvalue. Possible symmetry
> operations are C\*v or D\*v


>     2. Asymmetric top molecules have all different eigenvalues. The
> maximum rotational symmetry in such molecules is 2


>     3. Symmetric top molecules have 1 unique eigenvalue, which gives a
> unique rotation axis. All axial point groups are possible
> except the cubic groups (T & O) and I.


>     4. Spherical top molecules have all three eigenvalues equal. They
> have the rare T, O or I point groups.

Attribute:

    sch_symbol (str): Schoenflies symbol of the detected point group.

The default settings are usually sufficient.


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Molecule to determine point group for.


    * **tolerance** (*float*) – Distance tolerance to consider sites as
    symmetrically equivalent. Defaults to 0.3 Angstrom.


    * **eigen_tolerance** (*float*) – Tolerance to compare eigen values of
    the inertia tensor. Defaults to 0.01.


    * **matrix_tolerance** (*float*) – Tolerance used to generate the full set of
    symmetry operations of the point group.



#### _analyze()

#### _check_R2_axes_asym()
Test for 2-fold rotation along the principal axes.

Used to handle asymmetric top molecules.


#### _check_perpendicular_r2_axis(axis)
Checks for R2 axes perpendicular to unique axis.

For handling symmetric top molecules.


#### _check_rot_sym(axis)
Determines the rotational symmetry about supplied axis.

Used only for symmetric top molecules which has possible rotational symmetry
operations > 2.


#### _static_ _combine_eq_sets(equiv_sets, sym_ops)
Combines the dicts of _get_equivalent_atom_dicts into one.


* **Parameters**


    * **equiv_sets** (*dict*) – Map of equivalent atoms onto each other (i.e. indices to indices).


    * **sym_ops** (*dict*) – Map of symmetry operations that map atoms onto each other.



* **Returns**

    The returned dictionary has two possible keys:

    `eq_sets`:
    A dictionary of indices mapping to sets of indices,
    each key maps to indices of all equivalent atoms.
    The keys are guaranteed to be not equivalent.

    `sym_ops`:
    Twofold nested dictionary.
    `operations[i][j]` gives the symmetry operation
    that maps atom `i` unto `j`.




* **Return type**

    dict



#### _find_mirror(axis)
Looks for mirror symmetry of specified type about axis.

Possible types are “h” or “vd”. Horizontal (h) mirrors are perpendicular to the
axis while vertical (v) or diagonal (d) mirrors are parallel. v mirrors has atoms
lying on the mirror plane while d mirrors do not.


#### _find_spherical_axes()
Looks for R5, R4, R3 and R2 axes in spherical top molecules.

Point group T molecules have only one unique 3-fold and one unique 2-fold axis. O
molecules have one unique 4, 3 and 2-fold axes. I molecules have a unique 5-fold
axis.


#### _get_eq_sets()
Calculates the dictionary for mapping equivalent atoms onto each other.


* **Parameters**

    **None** –



* **Returns**

    The returned dictionary has two possible keys:

    `eq_sets`:
    A dictionary of indices mapping to sets of indices,
    each key maps to indices of all equivalent atoms.
    The keys are guaranteed to be not equivalent.

    `sym_ops`:
    Twofold nested dictionary.
    `operations[i][j]` gives the symmetry operation
    that maps atom `i` unto `j`.




* **Return type**

    dict



#### _get_smallest_set_not_on_axis(axis)
Returns the smallest list of atoms with the same species and distance from
origin AND does not lie on the specified axis.

This maximal set limits the possible rotational symmetry operations, since atoms
lying on a test axis is irrelevant in testing rotational symmetryOperations.


#### _proc_asym_top()
Handles asymmetric top molecules, which cannot contain rotational symmetry
larger than 2.


#### _proc_cyclic()
Handles cyclic group molecules.


#### _proc_dihedral()
Handles dihedral group molecules, i.e those with intersecting R2 axes and a
main axis.


#### _proc_linear()

#### _proc_no_rot_sym()
Handles molecules with no rotational symmetry.

Only possible point groups are C1, Cs and Ci.


#### _proc_sph_top()
Handles Spherical Top Molecules, which belongs to the T, O or I point
groups.


#### _proc_sym_top()
Handles symmetric top molecules which has one unique eigenvalue whose
corresponding principal axis is a unique rotational axis.

More complex handling required to look for R2 axes perpendicular to this unique
axis.


#### get_equivalent_atoms()
Returns sets of equivalent atoms with symmetry operations.


* **Parameters**

    **None** –



* **Returns**

    The returned dictionary has two possible keys:

    `eq_sets`:
    A dictionary of indices mapping to sets of indices,
    each key maps to indices of all equivalent atoms.
    The keys are guaranteed to be not equivalent.

    `sym_ops`:
    Twofold nested dictionary.
    `operations[i][j]` gives the symmetry operation
    that maps atom `i` unto `j`.




* **Return type**

    dict



#### get_pointgroup()
Returns a PointGroup object for the molecule.


#### get_rotational_symmetry_number()
Return the rotational symmetry number.


#### get_symmetry_operations()
Return symmetry operations as a list of SymmOp objects. Returns Cartesian coord
symmops.


* **Returns**

    List of symmetry operations.



* **Return type**

    ([[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)])



#### inversion_op(_ = SymmOp(affine_matrix=array([[-1., -0., -0.,  0.],        [-0., -1., -0.,  0.],        [-0., -0., -1.,  0.],        [-0., -0., -0.,  1.]])_ )

#### is_valid_op(symmop)
Check if a particular symmetry operation is a valid symmetry operation for a
molecule, i.e., the operation maps all atoms to another equivalent atom.


* **Parameters**

    **symmop** ([*SymmOp*](pymatgen.core.md#pymatgen.core.operations.SymmOp)) – Symmetry operation to test.



* **Returns**

    Whether SymmOp is valid for Molecule.



* **Return type**

    (bool)



#### symmetrize_molecule()
Returns a symmetrized molecule.

The equivalent atoms obtained via
`get_equivalent_atoms()`
are rotated, mirrored… unto one position.
Then the average position is calculated.
The average position is rotated, mirrored… back with the inverse
of the previous symmetry operations, which gives the
symmetrized molecule


* **Parameters**

    **None** –



* **Returns**

    The returned dictionary has three possible keys:

    `sym_mol`:
    A symmetrized molecule instance.

    `eq_sets`:
    A dictionary of indices mapping to sets of indices,
    each key maps to indices of all equivalent atoms.
    The keys are guaranteed to be not equivalent.

    `sym_ops`:
    Twofold nested dictionary.
    `operations[i][j]` gives the symmetry operation
    that maps atom `i` unto `j`.




* **Return type**

    dict



### _class_ PointGroupOperations(sch_symbol, operations, tol: float = 0.1)
Bases: `list`

Defines a point group, which is essentially a sequence of symmetry operations.


#### sch_symbol()
Schoenflies symbol of the point group.


* **Type**

    str



* **Parameters**


    * **sch_symbol** (*str*) – Schoenflies symbol of the point group.


    * **operations** (*[*[*SymmOp*](pymatgen.core.md#pymatgen.core.operations.SymmOp)*]*) – Initial set of symmetry operations. It is
    sufficient to provide only just enough operations to generate
    the full set of symmetries.


    * **tol** (*float*) – Tolerance to generate the full set of symmetry
    operations.



### _class_ SpacegroupAnalyzer(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), symprec: float | None = 0.01, angle_tolerance=5.0)
Bases: `object`

Takes a pymatgen.core.structure.Structure object and a symprec.

Uses spglib to perform various symmetry finding operations.


* **Parameters**


    * **structure** (*Structure/IStructure*) – Structure to find symmetry


    * **symprec** (*float*) – Tolerance for symmetry finding. Defaults to 0.01,
    which is fairly strict and works well for properly refined
    structures with atoms in the proper symmetry coordinates. For
    structures with slight deviations from their proper atomic
    positions (e.g., structures relaxed with electronic structure
    codes), a looser tolerance of 0.1 (the value used in Materials
    Project) is often needed.


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.



#### _get_symmetry()
Get the symmetry operations associated with the structure.


* **Returns**

    Symmetry operations as a tuple of two equal length sequences.
    (rotations, translations). “rotations” is the numpy integer array
    of the rotation matrices for scaled positions
    “translations” gives the numpy float64 array of the translation
    vectors in scaled positions.



#### find_primitive(keep_site_properties=False)
Find a primitive version of the unit cell.


* **Parameters**

    **keep_site_properties** (*bool*) – Whether to keep the input site properties (including
    magnetic moments) on the sites that are still present after the refinement. Note:
    This is disabled by default because the magnetic moments are not always directly
    transferable between unit cell definitions. For instance, long-range magnetic
    ordering or antiferromagnetic character may no longer be present (or exist in
    the same way) in the returned structure. If keep_site_properties is True,
    each site retains the same site property as in the original structure without
    further adjustment.



* **Returns**

    A primitive cell in the input cell is searched and returned
    as an Structure object. If no primitive cell is found, None is
    returned.



#### get_conventional_standard_structure(international_monoclinic=True, keep_site_properties=False)
Gives a structure with a conventional cell according to certain standards. The
standards are defined in Setyawan, W., & Curtarolo, S. (2010). High-throughput
electronic band structure calculations: Challenges and tools. Computational
Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010 They
basically enforce as much as possible norm(a1)<norm(a2)<norm(a3). NB This is not
necessarily the same as the standard settings within the International Tables of
Crystallography, for which get_refined_structure should be used instead.


* **Parameters**


    * **international_monoclinic** (*bool*) – Whether to convert to proper international convention
    such that beta is the non-right angle.


    * **keep_site_properties** (*bool*) – Whether to keep the input site properties (including
    magnetic moments) on the sites that are still present after the refinement. Note:
    This is disabled by default because the magnetic moments are not always directly
    transferable between unit cell definitions. For instance, long-range magnetic
    ordering or antiferromagnetic character may no longer be present (or exist in
    the same way) in the returned structure. If keep_site_properties is True,
    each site retains the same site property as in the original structure without
    further adjustment.



* **Returns**

    The structure in a conventional standardized cell



#### get_conventional_to_primitive_transformation_matrix(international_monoclinic=True)
Gives the transformation matrix to transform a conventional unit cell to a
primitive cell according to certain standards the standards are defined in
Setyawan, W., & Curtarolo, S. (2010). High-throughput electronic band structure
calculations: Challenges and tools. Computational Materials Science, 49(2),
299-312. doi:10.1016/j.commatsci.2010.05.010.


* **Parameters**

    **international_monoclinic** (*bool*) – Whether to convert to proper international convention
    such that beta is the non-right angle.



* **Returns**

    Transformation matrix to go from conventional to primitive cell



#### get_crystal_system()
Get the crystal system for the structure, e.g., (triclinic, orthorhombic,
cubic, etc.).


* **Raises**

    **ValueError** – on invalid space group numbers < 1 or > 230.



* **Returns**

    Crystal system for structure



* **Return type**

    (str)



#### get_hall()
Returns Hall symbol for structure.


* **Returns**

    Hall symbol



* **Return type**

    (str)



#### get_ir_reciprocal_mesh(mesh=(10, 10, 10), is_shift=(0, 0, 0))
k-point mesh of the Brillouin zone generated taken into account symmetry.The
method returns the irreducible kpoints of the mesh and their weights.


* **Parameters**


    * **mesh** (*3x1 array*) – The number of kpoint for the mesh needed in
    each direction


    * **is_shift** (*3x1 array*) – Whether to shift the kpoint grid. (1, 1,


    * **0.5** (*1**) **means all points are shifted by*) –


    * **0.5** –


    * **0.5.** –



* **Returns**

    A list of irreducible kpoints and their weights as a list of
    tuples [(ir_kpoint, weight)], with ir_kpoint given
    in fractional coordinates



#### get_ir_reciprocal_mesh_map(mesh=(10, 10, 10), is_shift=(0, 0, 0))
Same as ‘get_ir_reciprocal_mesh’ but the full grid together with the mapping
that maps a reducible to an irreducible kpoint is returned.


* **Parameters**


    * **mesh** (*3x1 array*) – The number of kpoint for the mesh needed in
    each direction


    * **is_shift** (*3x1 array*) – Whether to shift the kpoint grid. (1, 1,


    * **0.5** (*1**) **means all points are shifted by*) –


    * **0.5** –


    * **0.5.** –



* **Returns**

    A tuple containing two numpy.ndarray. The first is the mesh in
    fractional coordinates and the second is an array of integers
    that maps all the reducible kpoints from to irreducible ones.



#### get_kpoint_weights(kpoints, atol=1e-05)
Calculate the weights for a list of kpoints.


* **Parameters**


    * **kpoints** (*Sequence*) – Sequence of kpoints. np.arrays is fine. Note
    that the code does not check that the list of kpoints
    provided does not contain duplicates.


    * **atol** (*float*) – Tolerance for fractional coordinates comparisons.



* **Returns**

    List of weights, in the SAME order as kpoints.



#### get_lattice_type()
Get the lattice for the structure, e.g., (triclinic, orthorhombic, cubic,
etc.).This is the same as the crystal system with the exception of the
hexagonal/rhombohedral lattice.


* **Raises**

    **ValueError** – on invalid space group numbers < 1 or > 230.



* **Returns**

    Lattice type for structure



* **Return type**

    (str)



#### get_point_group_operations(cartesian=False)
Return symmetry operations as a list of SymmOp objects. By default returns
fractional coord symmops. But Cartesian can be returned too.


* **Parameters**

    **cartesian** (*bool*) – Whether to return SymmOps as Cartesian or
    direct coordinate operations.



* **Returns**

    Point group symmetry operations.



* **Return type**

    list[[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)]



#### get_point_group_symbol()
Get the point group associated with the structure.


* **Returns**

    Point group for structure.



* **Return type**

    (Pointgroup)



#### get_primitive_standard_structure(international_monoclinic=True, keep_site_properties=False)
Gives a structure with a primitive cell according to certain standards the
standards are defined in Setyawan, W., & Curtarolo, S. (2010). High-throughput
electronic band structure calculations: Challenges and tools. Computational
Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010.


* **Parameters**


    * **international_monoclinic** (*bool*) – Whether to convert to proper international convention
    such that beta is the non-right angle.


    * **keep_site_properties** (*bool*) – Whether to keep the input site properties (including
    magnetic moments) on the sites that are still present after the refinement. Note:
    This is disabled by default because the magnetic moments are not always directly
    transferable between unit cell definitions. For instance, long-range magnetic
    ordering or antiferromagnetic character may no longer be present (or exist in
    the same way) in the returned structure. If keep_site_properties is True,
    each site retains the same site property as in the original structure without
    further adjustment.



* **Returns**

    The structure in a primitive standardized cell



#### get_refined_structure(keep_site_properties=False)
Get the refined structure based on detected symmetry. The refined structure is
a *conventional* cell setting with atoms moved to the expected symmetry positions.


* **Parameters**

    **keep_site_properties** (*bool*) – Whether to keep the input site properties (including
    magnetic moments) on the sites that are still present after the refinement. Note:
    This is disabled by default because the magnetic moments are not always directly
    transferable between unit cell definitions. For instance, long-range magnetic
    ordering or antiferromagnetic character may no longer be present (or exist in
    the same way) in the returned structure. If keep_site_properties is True,
    each site retains the same site property as in the original structure without
    further adjustment.



* **Returns**

    Refined structure.



#### get_space_group_number()
Get the international spacegroup number (e.g., 62) for structure.


* **Returns**

    International spacegroup number for structure.



* **Return type**

    (int)



#### get_space_group_operations()
Get the SpacegroupOperations for the Structure.


* **Returns**

    SpacegroupOperations object.



#### get_space_group_symbol()
Get the spacegroup symbol (e.g., Pnma) for structure.


* **Returns**

    Spacegroup symbol for structure.



* **Return type**

    (str)



#### get_symmetrized_structure()
Get a symmetrized structure. A symmetrized structure is one where the sites
have been grouped into symmetrically equivalent groups.


* **Returns**

    pymatgen.symmetry.structure.SymmetrizedStructure object.



#### get_symmetry_dataset()
Returns the symmetry dataset as a dict.


* **Returns**

    With the following properties:

        number: International space group number
        international: International symbol
        hall: Hall symbol
        transformation_matrix: Transformation matrix from lattice of
        input cell to Bravais lattice L^bravais = L^original \* Tmat
        origin shift: Origin shift in the setting of “Bravais lattice”
        rotations, translations: Rotation matrices and translation
        vectors. Space group operations are obtained by
        [(r,t) for r, t in zip(rotations, translations)]
        wyckoffs: Wyckoff letters




* **Return type**

    (dict)



#### get_symmetry_operations(cartesian=False)
Return symmetry operations as a list of SymmOp objects. By default returns
fractional coord symmops. But Cartesian can be returned too.


* **Returns**

    List of symmetry operations.



* **Return type**

    ([[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)])



#### is_laue()
Check if the point group of the structure has Laue symmetry (centrosymmetry).


### _class_ SpacegroupOperations(int_symbol, int_number, symmops)
Bases: `list`

Represents a space group, which is a collection of symmetry operations.


* **Parameters**


    * **int_symbol** (*str*) – International symbol of the spacegroup.


    * **int_number** (*int*) – International number of the spacegroup.


    * **symmops** (*[*[*SymmOp*](pymatgen.core.md#pymatgen.core.operations.SymmOp)*]*) – Symmetry operations associated with the
    spacegroup.



#### are_symmetrically_equivalent(sites1, sites2, symm_prec=0.001)
Given two sets of PeriodicSites, test if they are actually symmetrically
equivalent under this space group. Useful, for example, if you want to test if
selecting atoms 1 and 2 out of a set of 4 atoms are symmetrically the same as
selecting atoms 3 and 4, etc.

One use is in PartialRemoveSpecie transformation to return only
symmetrically distinct arrangements of atoms.


* **Parameters**


    * **sites1** (*[*[*PeriodicSite*](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)*]*) – 1st set of sites


    * **sites2** (*[*[*PeriodicSite*](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)*]*) – 2nd set of sites


    * **symm_prec** (*float*) – Tolerance in atomic distance to test if atoms
    are symmetrically similar.



* **Returns**

    Whether the two sets of sites are symmetrically
    equivalent.



* **Return type**

    (bool)



### _get_symmetry_dataset(cell, symprec, angle_tolerance)
Simple wrapper to cache results of spglib.get_symmetry_dataset since this call is
expensive.


### cluster_sites(mol: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), tol: float, give_only_index: bool = False)
Cluster sites based on distance and species type.


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Molecule **with origin at center of mass**.


    * **tol** (*float*) – Tolerance to use.


    * **give_only_index** (*bool*) – Whether to return only the index of the
    origin site, instead of the site itself. Defaults to False.



* **Returns**

    origin_site is a site at the center
    of mass (None if there are no origin atoms). clustered_sites is a
    dict of {(avg_dist, species_and_occu): [list of sites]}



* **Return type**

    (origin_site, clustered_sites)



### generate_full_symmops(symmops: Sequence[[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)], tol: float)
Recursive algorithm to permute through all possible combinations of the initially
supplied symmetry operations to arrive at a complete set of operations mapping a
single atom to all other equivalent atoms in the point group. This assumes that the
initial number already uniquely identifies all operations.


* **Parameters**


    * **symmops** (*list**[*[*SymmOp*](pymatgen.core.md#pymatgen.core.operations.SymmOp)*]*) – Initial set of symmetry operations.


    * **tol** (*float*) – Tolerance for detecting symmetry.



* **Returns**

    Full set of symmetry operations.



* **Return type**

    list[[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)]



### iterative_symmetrize(mol, max_n=10, tolerance=0.3, epsilon=0.01)
Returns a symmetrized molecule.

The equivalent atoms obtained via
`get_equivalent_atoms()`
are rotated, mirrored… unto one position.
Then the average position is calculated.
The average position is rotated, mirrored… back with the inverse
of the previous symmetry operations, which gives the
symmetrized molecule


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – A pymatgen Molecule instance.


    * **max_n** (*int*) – Maximum number of iterations.


    * **tolerance** (*float*) – Tolerance for detecting symmetry.
    Gets passed as Argument into
    ~pymatgen.analyzer.symmetry.PointGroupAnalyzer.


    * **epsilon** (*float*) – If the elementwise absolute difference of two
    subsequently symmetrized structures is smaller epsilon,
    the iteration stops before `max_n` is reached.



* **Returns**

    The returned dictionary has three possible keys:

    `sym_mol`:
    A symmetrized molecule instance.

    `eq_sets`:
    A dictionary of indices mapping to sets of indices,
    each key maps to indices of all equivalent atoms.
    The keys are guaranteed to be not equivalent.

    `sym_ops`:
    Twofold nested dictionary.
    `operations[i][j]` gives the symmetry operation
    that maps atom `i` unto `j`.




* **Return type**

    dict


## pymatgen.symmetry.bandstructure module

Provides a class for interacting with KPath classes to
generate high-symmetry k-paths using different conventions.


### _class_ HighSymmKpath(structure, has_magmoms=False, magmom_axis=None, path_type='setyawan_curtarolo', symprec=0.01, angle_tolerance=5, atol=1e-05)
Bases: `KPathBase`

This class generates path along high symmetry lines in the
Brillouin zone according to different conventions.
The class is designed to be used with a specific primitive
cell setting. The definitions for the primitive cell
used can be found in: Computational Materials Science,
49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010.
The space group analyzer can be used to produce the correct
primitive structure
(method get_primitive_standard_structure(international_monoclinic=False)).
Ensure input structure is correct before ‘get_kpoints()’ method is used.
See individual KPath classes for details on specific conventions.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object


    * **has_magmoms** (*bool*) – Whether the input structure contains
    magnetic moments as site properties with the key ‘magmom.’
    Values may be in the form of 3-component vectors given in
    the basis of the input lattice vectors, in
    which case the spin axis will default to a_3, the third
    real-space lattice vector (this triggers a warning).


    * **magmom_axis** (*list** or **numpy array*) – 3-component vector specifying
    direction along which magnetic moments given as scalars
    should point. If all magnetic moments are provided as
    vectors then this argument is not used.


    * **path_type** (*str*) – Chooses which convention to use to generate
    the high symmetry path. Options are: ‘setyawan_curtarolo’, ‘hinuma’,
    ‘latimer_munro’ for the Setyawan & Curtarolo, Hinuma et al., and
    Latimer & Munro conventions. Choosing ‘all’ will generate one path
    with points from all three conventions. Equivalent labels between
    each will also be generated. Order will always be Latimer & Munro,
    Setyawan & Curtarolo, and Hinuma et al. Lengths for each of the paths
    will also be generated and output as a list. Note for ‘all’ the user
    will have to alter the labels on their own for plotting.


    * **symprec** (*float*) – Tolerance for symmetry finding


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to determine symmetric
    equivalence of points and lines on the BZ.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _get_hin_kpath(symprec, angle_tolerance, atol, tri)

* **Returns**

    Hinuma et al. k-path with labels.



#### _get_klabels(lm_bs, sc_bs, hin_bs, rpg)

* **Returns**

    Dictionary of equivalent labels for paths if ‘all’ is chosen.
    If an exact kpoint match cannot be found, symmetric equivalency will be
    searched for and indicated with an asterisk in the equivalent label.
    If an equivalent label can still not be found, or the point is not in
    the explicit kpath, its equivalent label will be set to itself in the output.



* **Return type**

    labels (dict)



#### _get_lm_kpath(has_magmoms, magmom_axis, symprec, angle_tolerance, atol)

* **Returns**

    Latimer and Munro k-path with labels.



#### _get_sc_kpath(symprec, angle_tolerance, atol)

* **Returns**

    Setyawan and Curtarolo k-path with labels.



#### _property_ equiv_labels()
Returns:
The correspondence between the kpoint symbols in the Latimer and
Munro convention, Setyawan and Curtarolo, and Hinuma
conventions respectively. Only generated when path_type = ‘all’.


#### _static_ get_continuous_path(bandstructure)
Obtain a continuous version of an inputted path using graph theory.
This routine will attempt to add connections between nodes of
odd-degree to ensure a Eulerian path can be formed. Initial
k-path must be able to be converted to a connected graph. See
npj Comput Mater 6, 112 (2020). 10.1038/s41524-020-00383-7
for more details.


* **Parameters**

    **bandstructure** (*BandstructureSymmLine*) – BandstructureSymmLine object.



* **Returns**

    New BandstructureSymmLine object with continuous path.



* **Return type**

    bandstructure (BandstructureSymmLine)



#### _property_ label_index()

* **Returns**

    The correspondence between numbers and kpoint symbols for the


combined kpath generated when path_type = ‘all’. None otherwise.


#### _property_ path_lengths()
Returns:
List of lengths of the Latimer and Munro, Setyawan and Curtarolo, and Hinuma
conventions in the combined HighSymmKpath object when path_type = ‘all’ respectively.
None otherwise.


#### _property_ path_type()
Returns:
The type of kpath chosen.

## pymatgen.symmetry.groups module

Defines SymmetryGroup parent class and PointGroup and SpaceGroup classes.
Shyue Ping Ong thanks Marc De Graef for his generous sharing of his
SpaceGroup data as published in his textbook “Structure of Materials”.


### _class_ PointGroup(\*args, \*\*kwargs)
Bases: `PointGroup`

Class representing a Point Group, with generators and symmetry operations.


#### symbol()
Full International or Hermann-Mauguin Symbol.


* **Type**

    str



#### generators()
List of generator matrices. Note that 3x3 matrices are used for Point Groups.


* **Type**

    list



#### symmetry_ops()
Full set of symmetry operations as matrices.


* **Type**

    list


Initializes a Point Group from its international symbol.


* **Parameters**

    **int_symbol** (*str*) – International or Hermann-Mauguin Symbol.



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ SpaceGroup(\*args, \*\*kwargs)
Bases: `SpaceGroup`

Class representing a SpaceGroup.


#### symbol()
Full International or Hermann-Mauguin Symbol.


* **Type**

    str



#### int_number()
International number.


* **Type**

    int



#### generators()
List of generator matrices. Note that 4x4 matrices are used for Space Groups.


* **Type**

    list



#### order()
Order of Space Group.


* **Type**

    int


Initializes a Space Group from its full or abbreviated international
symbol. Only standard settings are supported.


* **Parameters**

    **int_symbol** (*str*) – Full International (e.g., “P2/m2/m2/m”) or
    Hermann-Mauguin Symbol (“Pmmm”) or abbreviated symbol. The
    notation is a LaTeX-like string, with screw axes being
    represented by an underscore. For example, “P6_3/mmc”.
    Alternative settings can be accessed by adding a “:identifier”.
    For example, the hexagonal setting  for rhombohedral cells can be
    accessed by adding a “:H”, e.g., “R-3m:H”. To find out all
    possible settings for a spacegroup, use the get_settings()
    classmethod. Alternative origin choices can be indicated by a
    translation vector, e.g., ‘Fm-3m(a-1/4,b-1/4,c-1/4)’.



#### _abc_impl(_ = <_abc._abc_data object_ )

### _class_ SymmetryGroup()
Bases: `Sequence`, [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)

Abstract class representing a symmetry group.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### is_subgroup(supergroup: SymmetryGroup)
True if this group is a subgroup of the supplied group.


* **Parameters**

    **supergroup** (*SymmetryGroup*) – Supergroup to test.



* **Returns**

    True if this group is a subgroup of the supplied group.



* **Return type**

    bool



#### is_supergroup(subgroup: SymmetryGroup)
True if this group is a supergroup of the supplied group.


* **Parameters**

    **subgroup** (*SymmetryGroup*) – Subgroup to test.



* **Returns**

    True if this group is a supergroup of the supplied group.



* **Return type**

    bool



#### _abstract property_ symmetry_ops(_: set[[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)_ )
Returns:
List of symmetry operations associated with the group.


#### to_latex_string()

* **Returns**

    A latex formatted group symbol with proper subscripts and overlines.



### in_array_list(array_list: list[np.ndarray] | np.ndarray, arr: np.ndarray, tol: float = 1e-05)
Extremely efficient nd-array comparison using numpy’s broadcasting. This
function checks if a particular array a, is present in a list of arrays.
It works for arrays of any size, e.g., even matrix searches.


* **Parameters**


    * **array_list** (*[**array**]*) – A list of arrays to compare to.


    * **arr** (*array*) – The test array for comparison.


    * **tol** (*float*) – The tolerance. Defaults to 1e-5. If 0, an exact match is done.



* **Returns**

    (bool)



### sg_symbol_from_int_number(int_number: int, hexagonal: bool = True)
Obtains a SpaceGroup name from its international number.


* **Parameters**


    * **int_number** (*int*) – International number.


    * **hexagonal** (*bool*) – For rhombohedral groups, whether to return the
    hexagonal setting (default) or rhombohedral setting.



* **Returns**

    Spacegroup symbol



* **Return type**

    str


## pymatgen.symmetry.kpath module

Provides classes for generating high-symmetry k-paths using different conventions.


### _class_ KPathBase(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), symprec: float = 0.01, angle_tolerance=5, atol=1e-05, \*args, \*\*kwargs)
Bases: `object`

This is the base class for classes used to generate high-symmetry
paths in reciprocal space (k-paths) for band structure calculations.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object.


    * **symprec** (*float*) – Tolerance for symmetry finding.


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to compare structures
    and determine symmetric equivalence of points and lines in the BZ.


    * **\*args** – Other arguments supported by subclasses.


    * **\*\*kwargs** – Other keyword arguments supported by subclasses.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### get_kpoints(line_density=20, coords_are_cartesian=True)

* **Returns**

    kpoints along the path in Cartesian coordinates


together with the critical-point labels.


#### _property_ kpath()
Returns:
The symmetry line path in reciprocal space.


#### _property_ lattice()
Returns:
The real space lattice.


#### _property_ rec_lattice()
Returns:
The reciprocal space lattice.


#### _property_ structure()
Returns:
The input structure.


### _class_ KPathLatimerMunro(structure, has_magmoms=False, magmom_axis=None, symprec=0.01, angle_tolerance=5, atol=1e-05)
Bases: `KPathBase`

This class looks for a path along high-symmetry lines in the
Brillouin zone. It is based on the method outlined in:
npj Comput Mater 6, 112 (2020). 10.1038/s41524-020-00383-7
The user should ensure that the unit cell of the input structure
is as reduced as possible, i.e. that there is no linear
combination of lattice vectors which can produce a vector of
lesser magnitude than the given set (this is required to
obtain the correct Brillouin zone within the current
implementation). This is checked during initialization and a
warning is issued if the condition is not fulfilled.
In the case of magnetic structures, care must also be taken to
provide the magnetic primitive cell (i.e. that which reproduces
the entire crystal, including the correct magnetic ordering,
upon application of lattice translations). There is no algorithm to

> check for this, so if the input structure is

incorrect, the class will output the incorrect k-path without
any warning being issued.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object


    * **has_magmoms** (*bool*) – Whether the input structure contains
    magnetic moments as site properties with the key ‘magmom.’
    Values may be in the form of 3-component vectors given in
    the basis of the input lattice vectors, or as scalars, in
    which case the spin axis will default to a_3, the third
    real-space lattice vector (this triggers a warning).


    * **magmom_axis** (*list** or **numpy array*) – 3-component vector specifying
    direction along which magnetic moments given as scalars
    should point. If all magnetic moments are provided as
    vectors then this argument is not used.


    * **symprec** (*float*) – Tolerance for symmetry finding


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to determine symmetric
    equivalence of points and lines in the BZ.



#### _static_ LabelPoints(index)
Axes used in generating labels for Latimer-Munro convention.


#### _static_ LabelSymbol(index)
Letters used in generating labels for the Latimer-Munro convention.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _all_ints(arr, atol)

#### _static_ _apply_op_to_magmom(r, magmom)

#### _choose_path(key_points, key_points_inds_orbits, key_lines_inds_orbits, little_groups_points, little_groups_lines)

#### _static_ _closewrapped(pos1, pos2, tolerance)

#### _convert_all_magmoms_to_vectors(magmom_axis, axis_specified)

#### _get_IRBZ(recip_point_group, W, key_points, face_center_inds, atol)

#### _get_coset_factor(G, H)

#### _get_key_line_orbits(key_points, key_lines, key_points_inds_orbits)

#### _static_ _get_key_lines(key_points, bz_as_key_point_inds)

#### _get_key_point_orbits(key_points)

#### _get_key_points()

#### _get_ksymm_kpath(has_magmoms, magmom_axis, axis_specified, symprec, angle_tolerance, atol)

#### _get_little_groups(key_points, key_points_inds_orbits, key_lines_inds_orbits)

#### _get_magnetic_symmetry_operations(struct, grey_ops, atol)

#### _get_max_cosine_labels(max_cosine_orbits_orig, key_points_inds_orbits, atol)

#### _get_orbit_labels(orbit_cosines_orig, key_points_inds_orbits, atol)

#### _static_ _get_reciprocal_point_group(ops, R, A)

#### _static_ _get_reciprocal_point_group_dict(recip_point_group, atol)

#### _kpath(_: dict[str, Any] | Non_ )

#### _static_ _op_maps_IRBZ_to_self(op, IRBZ_points, atol)

#### _static_ _reduce_IRBZ(IRBZ_points, boundaries, g, atol)

#### _static_ _reduce_cosines_array(orbit_cosines, pop_orbits, pop_labels)

#### _property_ mag_type()
Returns:
The type of magnetic space group as a string. Current implementation does not
distinguish between types 3 and 4, so return value is ‘3/4’. If has_magmoms is
False, returns ‘0’.


### _class_ KPathSeek(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), symprec: float = 0.01, angle_tolerance=5, atol=1e-05, system_is_tri=True)
Bases: `KPathBase`

This class looks for a path along high-symmetry lines in the Brillouin zone. It is based on
Hinuma, Y., Pizzi, G., Kumagai, Y., Oba, F., & Tanaka, I. (2017). Band structure diagram paths
based on crystallography. Computational Materials Science, 128, 140-184.
[https://doi.org/10.1016/j.commatsci.2016.10.015](https://doi.org/10.1016/j.commatsci.2016.10.015). It should be used with primitive structures that
comply with the definition given in the paper. The symmetry is determined by spglib using the
SpacegroupAnalyzer class. k-points are generated using the get_kpoints() method for the
reciprocal cell basis defined in the paper.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object


    * **symprec** (*float*) – Tolerance for symmetry finding


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to determine edge cases
    for settings of structures.


    * **system_is_tri** (*bool*) – Indicates if the system is time-reversal
    invariant.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _kpath(_: dict[str, Any] | Non_ )

#### _static_ _trans_sc_to_Hin(sub_class)

### _class_ KPathSetyawanCurtarolo(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), symprec: float = 0.01, angle_tolerance=5, atol=1e-05)
Bases: `KPathBase`

This class looks for a path along high-symmetry lines in
the Brillouin zone.
It is based on Setyawan, W., & Curtarolo, S. (2010).
High-throughput electronic band structure calculations:
Challenges and tools. Computational Materials Science,
49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
It should be used with primitive structures that
comply with the definition given in the paper.
The symmetry is determined by spglib using the
SpacegroupAnalyzer class. The analyzer can be used to
produce the correct primitive structure with the method
get_primitive_standard_structure(international_monoclinic=False).
A warning will signal possible compatibility problems
with the given structure. k-points generated using the get_kpoints() method
are returned for the reciprocal cell basis defined in the paper.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object.


    * **symprec** (*float*) – Tolerance for symmetry finding.


    * **angle_tolerance** (*float*) – Angle tolerance for symmetry finding.


    * **atol** (*float*) – Absolute tolerance used to compare the input
    structure with the one expected as primitive standard.
    A warning will be issued if the cells don’t match.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _kpath(_: dict[str, Any] | Non_ )

#### bcc()
BCC Path.


#### bctet1(c, a)
BCT1 Path.


#### bctet2(c, a)
BCT2 Path.


#### _property_ conventional()
Returns:
The conventional cell structure.


#### cubic()
CUB Path.


#### fcc()
FCC Path.


#### hex()
HEX Path.


#### mcl(b, c, beta)
MCL Path.


#### mclc1(a, b, c, alpha)
MCLC1 Path.


#### mclc2(a, b, c, alpha)
MCLC2 Path.


#### mclc3(a, b, c, alpha)
MCLC3 Path.


#### mclc4(a, b, c, alpha)
MCLC4 Path.


#### mclc5(a, b, c, alpha)
MCLC5 Path.


#### orc()
ORC Path.


#### orcc(a, b, c)
ORCC Path.


#### orcf1(a, b, c)
ORFC1 Path.


#### orcf2(a, b, c)
ORFC2 Path.


#### orcf3(a, b, c)
ORFC3 Path.


#### orci(a, b, c)
ORCI Path.


#### _property_ prim()
Returns:
The primitive cell structure.


#### _property_ prim_rec()
Returns:
The primitive reciprocal cell structure.


#### rhl1(alpha)
RHL1 Path.


#### rhl2(alpha)
RHL2 Path.


#### tet()
TET Path.


#### tria()
TRI1a Path.


#### trib()
TRI1b Path.

## pymatgen.symmetry.maggroups module

Magnetic space groups.


### _class_ MagneticSpaceGroup(\*args, \*\*kwargs)
Bases: `MagneticSpaceGroup`

Representation of a magnetic space group.

Initializes a MagneticSpaceGroup from its Belov, Neronova and
Smirnova (BNS) number supplied as a list or its label supplied
as a string. To create a magnetic structure in pymatgen, the
Structure.from_magnetic_spacegroup() method can be used, which
relies on this class.

The main difference between magnetic space groups and normal
crystallographic space groups is the inclusion of a time reversal
operator that acts on an atom’s magnetic moment. This is
indicated by a prime symbol (’) next to the respective symmetry
operation in its label, e.g. the standard crystallographic
space group Pnma has magnetic subgroups Pn’ma, Pnm’a, Pnma’,
Pn’m’a, Pnm’a’, Pn’ma’, Pn’m’a’.

The magnetic space groups are classified as one of 4 types
where G = magnetic space group, and F = parent crystallographic
space group:


1. G=F no time reversal, i.e. the same as corresponding

    crystallographic group


2. G=F+F1’, “grey” groups, where avg. magnetic moment is zero,

    e.g. a paramagnet in zero ext. mag. field


3. G=D+(F-D)1’, where D is an equi-translation subgroup of F of

    index 2, lattice translations do not include time reversal


4. G=D+(F-D)1’, where D is an equi-class subgroup of F of index 2

There are two common settings for magnetic space groups, BNS
and OG. In case 4, the BNS setting != OG setting, and so a
transformation to go between the two settings is required:
specifically, the BNS setting is derived from D, and the OG
setting is derived from F.

This means that the OG setting refers to the unit cell if magnetic
order is neglected, and requires multiple unit cells to reproduce
the full crystal periodicity when magnetic moments are present.
This does not make the OG setting, in general, useful for
electronic structure calculations and the BNS setting is preferred.
However, this class does contain information on the OG setting and
can be initialized from OG labels or numbers if required.

Conventions: ITC monoclinic unique axis b, monoclinic cell choice 1,
hexagonal axis for trigonal groups, origin choice 2 for groups with
more than one origin choice (ISO-MAG).

Raw data comes from ISO-MAG, ISOTROPY Software Suite, iso.byu.edu
[http://stokes.byu.edu/iso/magnetic_data.txt](http://stokes.byu.edu/iso/magnetic_data.txt)
with kind permission from Professor Branton Campbell, BYU

Data originally compiled from:
(1) Daniel B. Litvin, Magnetic Group Tables (International Union

> of Crystallography, 2013) www.iucr.org/publ/978-0-9553602-2-0.


1. C. J. Bradley and A. P. Cracknell, The Mathematical Theory of
Symmetry in Solids (Clarendon Press, Oxford, 1972).

See [http://stokes.byu.edu/iso/magneticspacegroupshelp.php](http://stokes.byu.edu/iso/magneticspacegroupshelp.php) for more
information on magnetic symmetry.


* **Parameters**

    **id** – BNS number supplied as list of 2 ints or BNS label as
    str or index as int (1-1651) to iterate over all space groups



#### _abc_impl(_ = <_abc._abc_data object_ )

### _write_all_magnetic_space_groups_to_file(filename)
Write all magnetic space groups to a human-readable text file.
Should contain same information as text files provided by ISO-MAG.

## pymatgen.symmetry.settings module

This module provides classes for non-standard space-group settings.


### _class_ JonesFaithfulTransformation(P, p)
Bases: `object`

Transformation for space-groups defined in a non-standard setting.

Transform between settings using matrix P and origin shift vector p,
using same notation as reference.

Should initialize using from_transformation_string in Jones
faithful notation, given by a string specifying both a
transformation matrix and an origin shift, with parts delimited
by a semi-colon. Best shown by example:


* a,b,c;0,0,0 is the identity (no change)


* -b+c,a+c,-a+b+c;0,0,0 is R3:r to R3:h (rhombohedral to
hexagonal setting)


* a,b,c;-1/4,-1/4,-1/4 is Pnnn:1 to Pnnn:2 (change in origin
choice)


* b,c,a;-1/2,-1/2,-1/2 is Bbab:1 to Ccca:2 (change settin
and origin)

Can transform points (coords), lattices and symmetry operations.

Used for transforming magnetic space groups since these are
commonly used in multiple settings, due to needing to transform
between magnetic and non-magnetic settings.

See: International Tables for Crystallography (2016). Vol. A,
Chapter 1.5, pp. 75-106.


#### _property_ P(_: list[list[float]_ )
Transformation matrix.


#### _static_ _get_transformation_string_from_Pp(P: list[list[float]] | np.ndarray, p: list[float])

#### _classmethod_ from_origin_shift(origin_shift='0,0,0')
Construct SpaceGroupTransformation from its origin shift string.


* **Parameters**

    **origin_shift** (*str**, **optional*) – Defaults to “0,0,0”.



* **Returns**

    JonesFaithfulTransformation



#### _classmethod_ from_transformation_str(transformation_string='a,b,c;0,0,0')
Construct SpaceGroupTransformation from its transformation string.


* **Parameters**

    **transformation_string** (*str**, **optional*) – Defaults to “a,b,c;0,0,0”.



* **Returns**

    JonesFaithfulTransformation



#### _classmethod_ from_transformation_string(\*args, \*\*kwds)
from_transformation_string is deprecated!
Use from_transformation_str instead


#### _property_ inverse(_: JonesFaithfulTransformatio_ )
JonesFaithfulTransformation.


#### _property_ p(_: list[float_ )
Translation vector.


#### _static_ parse_transformation_string(transformation_string: str = 'a,b,c;0,0,0')

* **Parameters**

    **transformation_string** (*str**, **optional*) – Defaults to “a,b,c;0,0,0”.



* **Raises**

    **ValueError** – When transformation string fails to parse.



* **Returns**

    transformation matrix & vector



* **Return type**

    tuple[list[list[float]] | np.ndarray, list[float]]



#### transform_coords(coords: list[list[float]] | np.ndarray)
Takes a list of coordinates and transforms them.


#### transform_lattice(lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice))
Transforms a lattice.


#### transform_symmop(symmop: [SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp) | [MagSymmOp](pymatgen.core.md#pymatgen.core.operations.MagSymmOp))
Takes a symmetry operation and transforms it.


#### _property_ transformation_string(_: st_ )
Transformation string.

## pymatgen.symmetry.site_symmetries module

Provides analysis of site symmetries.


### get_shared_symmetry_operations(struct: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), pointops: list[list[[SymmOp](pymatgen.core.md#pymatgen.core.operations.SymmOp)]], tol: float = 0.1)
Get all the point group operations shared by a pair of atomic sites
in the form [[point operations of site index 1],[],…,[]].


* **Parameters**


    * **struct** – Pymatgen structure


    * **pointops** – list of point group operations from get_site_symmetries method


    * **tol** (*float*) – tolerance to find symmetry operations



* **Returns**

    list of lists of shared point operations for each pair of atomic sites



### get_site_symmetries(struct: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), precision: float = 0.1)
Get all the point group operations centered on each atomic site
in the form [[point operations of site index 1]…[[point operations of site index N]]].


* **Parameters**


    * **struct** – Pymatgen structure


    * **precision** (*float*) – tolerance to find symmetry operations



* **Returns**

    list of lists of point operations for each atomic site


## pymatgen.symmetry.structure module

This module implements symmetry-related structure forms.


### _class_ SymmetrizedStructure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), spacegroup: SpacegroupOperations, equivalent_positions: Sequence[int], wyckoff_letters: Sequence[str])
Bases: [`Structure`](pymatgen.core.md#pymatgen.core.structure.Structure)

This class represents a symmetrized structure, i.e. a structure
where the spacegroup and symmetry operations are defined. This class is
typically not called but instead is typically obtained by calling
pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.


#### equivalent_indices()
A list of lists of indices of the sites in the structure that are
considered equivalent based on the symmetry operations of the space group.


* **Type**

    list[List[int]]



* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Original structure


    * **spacegroup** (*SpacegroupOperations*) – An input SpacegroupOperations from SpacegroupAnalyzer.


    * **equivalent_positions** (*list**[**int**]*) – Equivalent positions from SpacegroupAnalyzer.


    * **wyckoff_letters** (*list**[**str**]*) – Wyckoff letters.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _properties(_: dic_ )

#### _sites(_: tuple[[PeriodicSite](pymatgen.core.md#pymatgen.core.sites.PeriodicSite), ..._ )

#### as_dict()
MSONable dict.


#### copy()
Copy of structure.


#### find_equivalent_sites(site: [PeriodicSite](pymatgen.core.md#pymatgen.core.sites.PeriodicSite))
Finds all symmetrically equivalent sites for a particular site.


* **Parameters**

    **site** ([*PeriodicSite*](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)) – A site in the structure



* **Raises**

    **ValueError** – if site is not in the structure.



* **Returns**

    List of all symmetrically equivalent sites.



* **Return type**

    ([[PeriodicSite](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)])



#### _classmethod_ from_dict(dct)

* **Parameters**

    **d** – Dict representation



* **Returns**

    SymmetrizedStructure