---
layout: default
title: pymatgen.symmetry.analyzer.md
nav_exclude: true
---

# pymatgen.symmetry.analyzer module

An interface to the excellent spglib library by Atsushi Togo
([http://spglib.sourceforge.net/](http://spglib.sourceforge.net/)) for pymatgen.

v1.0 - Now works with both ordered and disordered structure.
v2.0 - Updated for spglib 1.6.
v3.0 - pymatgen no longer ships with spglib. Instead, spglib (the python

> version) is now a dependency and the SpacegroupAnalyzer merely serves
> as an interface to spglib for pymatgen Structures.


### _class_ pymatgen.symmetry.analyzer.PointGroupAnalyzer(mol, tolerance=0.3, eigen_tolerance=0.01, matrix_tolerance=0.1)
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


    * **mol** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Molecule to determine point group for.


    * **tolerance** (*float*) – Distance tolerance to consider sites as
    symmetrically equivalent. Defaults to 0.3 Angstrom.


    * **eigen_tolerance** (*float*) – Tolerance to compare eigen values of
    the inertia tensor. Defaults to 0.01.


    * **matrix_tolerance** (*float*) – Tolerance used to generate the full set of
    symmetry operations of the point group.



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

    ([[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)])



#### inversion_op(_ = Rot: [[-1. -0. -0.]  [-0. -1. -0.]  [-0. -0. -1.]] tau [0. 0. 0._ )

#### is_valid_op(symmop)
Check if a particular symmetry operation is a valid symmetry operation for a
molecule, i.e., the operation maps all atoms to another equivalent atom.


* **Parameters**

    **symmop** ([*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)) – Symmetry operation to test.



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



### _class_ pymatgen.symmetry.analyzer.PointGroupOperations(sch_symbol, operations, tol: float = 0.1)
Bases: `list`

Defines a point group, which is essentially a sequence of symmetry operations.


#### sch_symbol()
Schoenflies symbol of the point group.


* **Parameters**


    * **sch_symbol** (*str*) – Schoenflies symbol of the point group.


    * **operations** (*[*[*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)*]*) – Initial set of symmetry operations. It is
    sufficient to provide only just enough operations to generate
    the full set of symmetries.


    * **tol** (*float*) – Tolerance to generate the full set of symmetry
    operations.



### _class_ pymatgen.symmetry.analyzer.SpacegroupAnalyzer(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), symprec: float | None = 0.01, angle_tolerance=5.0)
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

    list[[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)]



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

    [`pymatgen.symmetry.structure.SymmetrizedStructure`](pymatgen.symmetry.structure.md#pymatgen.symmetry.structure.SymmetrizedStructure) object.



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

    ([[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)])



#### is_laue()
Check if the point group of the structure has Laue symmetry (centrosymmetry).


### _class_ pymatgen.symmetry.analyzer.SpacegroupOperations(int_symbol, int_number, symmops)
Bases: `list`

Represents a space group, which is a collection of symmetry operations.


* **Parameters**


    * **int_symbol** (*str*) – International symbol of the spacegroup.


    * **int_number** (*int*) – International number of the spacegroup.


    * **symmops** (*[*[*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)*]*) – Symmetry operations associated with the
    spacegroup.



#### are_symmetrically_equivalent(sites1, sites2, symm_prec=0.001)
Given two sets of PeriodicSites, test if they are actually symmetrically
equivalent under this space group. Useful, for example, if you want to test if
selecting atoms 1 and 2 out of a set of 4 atoms are symmetrically the same as
selecting atoms 3 and 4, etc.

One use is in PartialRemoveSpecie transformation to return only
symmetrically distinct arrangements of atoms.


* **Parameters**


    * **sites1** (*[*[*PeriodicSite*](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)*]*) – 1st set of sites


    * **sites2** (*[*[*PeriodicSite*](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)*]*) – 2nd set of sites


    * **symm_prec** (*float*) – Tolerance in atomic distance to test if atoms
    are symmetrically similar.



* **Returns**

    Whether the two sets of sites are symmetrically
    equivalent.



* **Return type**

    (bool)



### pymatgen.symmetry.analyzer.cluster_sites(mol: [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), tol: float, give_only_index: bool = False)
Cluster sites based on distance and species type.


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – Molecule **with origin at center of mass**.


    * **tol** (*float*) – Tolerance to use.


    * **give_only_index** (*bool*) – Whether to return only the index of the
    origin site, instead of the site itself. Defaults to False.



* **Returns**

    origin_site is a site at the center
    of mass (None if there are no origin atoms). clustered_sites is a
    dict of {(avg_dist, species_and_occu): [list of sites]}



* **Return type**

    (origin_site, clustered_sites)



### pymatgen.symmetry.analyzer.generate_full_symmops(symmops: Sequence[[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)], tol: float)
Recursive algorithm to permute through all possible combinations of the initially
supplied symmetry operations to arrive at a complete set of operations mapping a
single atom to all other equivalent atoms in the point group. This assumes that the
initial number already uniquely identifies all operations.


* **Parameters**


    * **symmops** (*list**[*[*SymmOp*](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)*]*) – Initial set of symmetry operations.


    * **tol** (*float*) – Tolerance for detecting symmetry.



* **Returns**

    Full set of symmetry operations.



* **Return type**

    list[[SymmOp](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)]



### pymatgen.symmetry.analyzer.iterative_symmetrize(mol, max_n=10, tolerance=0.3, epsilon=0.01)
Returns a symmetrized molecule.

The equivalent atoms obtained via
`get_equivalent_atoms()`
are rotated, mirrored… unto one position.
Then the average position is calculated.
The average position is rotated, mirrored… back with the inverse
of the previous symmetry operations, which gives the
symmetrized molecule


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – A pymatgen Molecule instance.


    * **max_n** (*int*) – Maximum number of iterations.


    * **tolerance** (*float*) – Tolerance for detecting symmetry.
    Gets passed as Argument into
    `PointGroupAnalyzer`.


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