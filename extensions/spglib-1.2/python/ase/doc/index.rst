pyspglib for ASE
==========================================

This is document for Pyspglib for ASE (Atomic Simulation Environment).
Pyspglib is the python module to use spglib library.

How to build spglib python module
=================================

The C sources of spglib and interface for the python C/API are
compiled. The development environment for python and gcc are required
before starting to build.

1. Go to the :file:`python/ase` directory
2. Type the command::

    % python setup.py install --home=<my-directory>

The :file:`{<my-directory>}` is possibly current directory, :file:`.`.

3. Put ``lib/python`` path into :envvar:`$PYTHONPATH`, e.g., in your .bashrc.

How to use it
=============
1. Import spglib::

    from pyspglib import spglib

2. Call the methods with ASE Atoms object.

Methods
=======

The tolerance is given in Cartesian coordinates.

``get_spacegroup``
------------------
::

    spacegroup = get_spacegroup(atoms, symprec=1e-5)

``atoms`` is the object of ASE Atoms class. ``symprec`` is the float
variable, which is used as tolerance in symmetry search.

International space group symbol and the number are obtained as a string.

``get_symmetry``
----------------
::

    symmetry = get_symmetry(atoms, symprec=1e-5)

``atoms`` is the object of ASE Atoms class. ``symprec`` is the float
variable, which is used as tolerance in symmetry search.

Symmetry operations are obtained as a dictionary. The key ``rotation``
contains a numpy array of integer, which is "number of symmetry
operations" x "3x3 matrices". The key ``translation`` contains a numpy
array of float, which is "number of symmetry operations" x
"vectors". The orders of the rotation matrices and the translation
vectors correspond with each other, e.g. , the second symmetry
operation is organized by the second rotation matrix and second
translation vector in the respective arrays. The operations are
applied for the fractional coordinates (not for Cartesian
coordinates).

The rotation matrix and translation vector are used as follows::

    new_vector[3x1] = rotation[3x3] * vector[3x1] + translation[3x1]

The three values in the vector are given for the a, b, and c axes,
respectively.

``refine_cell``
-------------------------------
::

    lattice, scaled_positions, numbers = refine_cell(atoms, symprec=1e-5)

``atoms`` is the object of ASE Atoms class. ``symprec`` is the float
variable, which is used as tolerance in symmetry search. 

Bravais lattice (3x3 numpy array), atomic scaled positions (a numpy
array of [number_of_atoms,3]), and atomic numbers (a 1D numpy array)
that are symmetrized following space group type are returned.

``get_symmetry_dataset``
----------------------------
::

    dataset = get_symmetry_dataset(atoms, symprec=1e-5)

``dataset`` is a dictionary. The keys are:

* ``number``: International space group number
* ``international``: International symbol
* ``hall``: Hall symbol
* ``transformation_matrix``: Transformation matrix from lattice of input cell to Bravais lattice :math:`L^{bravais} = L^{original} * T`
* ``origin shift``: Origin shift in the setting of Bravais lattice
* ``wyckoffs``: Wyckoff letters
* ``equivalent_atoms``: Mapping table to equivalent atoms
* ``rotations`` and ``translations``: Rotation matrices and translation vectors. Space group operations are obtained by::

    [ ( r,t ) for r, t in zip( dataset['rotations'], dataset['translations'] ) ]


``find_primitive``
------------------
::

   lattice, scaled_positions, numbers = find_primitive(atoms, symprec=1e-5)

``atoms`` is the object of ASE Atoms class. ``symprec`` is the float
variable, which is used as tolerance in symmetry search.

When a primitive cell is found, lattice parameters (3x3 numpy array),
scaled positions (a numpy array of [number_of_atoms,3]), and atomic
numbers (a 1D numpy array) is returned. When no primitive cell is
found, (``None``, ``None``, ``None``) is returned.

``get_ir_reciprocal_mesh``
--------------------------

::

   mapping, grid = get_ir_reciprocal_mesh( mesh, atoms, is_shift=[0,0,0] )

Irreducible k-points are obtained from a sampling mesh of k-points.
``mesh`` is given by three integers by array and specifies mesh
numbers along reciprocal primitive axis. ``atoms`` is an Atoms object
of ASE. ``is_shift`` is given by the three integers by array. When
``is_shift`` is set for each reciprocal primitive axis, the mesh is
shifted along the axis in half of adjacent mesh points irrespective of
the mesh numbers. When the value is not 0, ``is_shift`` is set.

``mapping`` and ``grid`` are returned. ``grid`` gives the mesh points in
fractional coordinates in reciprocal space. ``mapping`` gives mapping to
the irreducible k-point indices that are obtained by ::

   np.unique( mapping )

Here ``np`` is the imported numpy module. The grid point is accessed
by ``grid[ index ]``.

For example, the irreducible k-points in fractional coordinates are
obtained by ::

   ir_grid = []
   mapping, grid = get_ir_reciprocal_mesh( [ 8, 8, 8 ], atoms, [ 1, 1, 1 ] )
   for i in np.unique( mapping ):
     ir_grid.append( grid[ i ] )

Example
=============

Examples are found in ``examples`` directory. An example code is :ref:`this <examples>`.

Instead of ASE's ``Atoms`` class, ``Atoms`` class in ``atoms.py`` in the ``examples`` directory may be used. To use this ``atoms.py``, ::

   from atoms import Atoms


.. |sflogo| image:: http://sflogo.sourceforge.net/sflogo.php?group_id=161614&type=1
            :target: http://sourceforge.net



|sflogo|
