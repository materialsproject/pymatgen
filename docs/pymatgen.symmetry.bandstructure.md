---
layout: default
title: pymatgen.symmetry.bandstructure.md
nav_exclude: true
---

# pymatgen.symmetry.bandstructure module

Provides a class for interacting with KPath classes to
generate high-symmetry k-paths using different conventions.


### _class_ pymatgen.symmetry.bandstructure.HighSymmKpath(structure, has_magmoms=False, magmom_axis=None, path_type='setyawan_curtarolo', symprec=0.01, angle_tolerance=5, atol=1e-05)
Bases: [`KPathBase`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase)

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


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – Structure object


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

Args:
bandstructure (BandstructureSymmLine): BandstructureSymmLine object.

Returns:
bandstructure (BandstructureSymmLine): New BandstructureSymmLine object with continuous path.


#### _property_ label_index()
Returns:
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