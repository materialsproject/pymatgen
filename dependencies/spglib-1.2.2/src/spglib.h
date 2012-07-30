/* spglib.h version 1.2.1 */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __spglib_H__
#define __spglib_H__

/* SPGCONST is used instead of 'const' so to avoid gcc warning. */
/* However there should be better way than this way.... */
#define SPGCONST

/*
  ------------------------------------------------------------------

  lattice: Lattice vectors (in Cartesian)

  [ [ a_x, b_x, c_x ],
  [ a_y, b_y, c_y ],
  [ a_z, b_z, c_z ] ]

  position: Atomic positions (in fractional coordinates)
  
  [ [ x1_a, x1_b, x1_c ], 
  [ x2_a, x2_b, x2_c ], 
  [ x3_a, x3_b, x3_c ], 
  ...                   ]

  types: Atom types, i.e., species identified by number

  [ type_1, type_2, type_3, ... ]

  rotation: Rotation matricies of symmetry operations

  each rotation is:
  [ [ r_aa, r_ab, r_ac ],
  [ r_ba, r_bb, r_bc ],
  [ r_ca, r_cb, r_cc ] ]

  translation: Translation vectors of symmetry operations

  each translation is:
  [ t_a, t_b, t_c ]

  symprec: Tolerance of atomic positions (in fractional coordinate)
  in finding symmetry operations

  ------------------------------------------------------------------

  Definitio of the operation:
  r : rotation     3x3 matrix
  t : translation  vector

  x_new = r * x + t:
  [ x_new_a ]   [ r_aa, r_ab, r_ac ]   [ x_a ]   [ t_a ]
  [ x_new_b ] = [ r_ba, r_bb, r_bc ] * [ x_b ] + [ t_b ]
  [ x_new_c ]   [ r_ca, r_cb, r_cc ]   [ x_c ]   [ t_c ]

  ------------------------------------------------------------------
*/

typedef struct {
  int spacegroup_number;
  int hall_number;
  char international_symbol[11];
  char hall_symbol[17];
  double transformation_matrix[3][3]; /* bravais_lattice = T * original_lattice */
  double origin_shift[3]; /* Origin shift in Bravais lattice */
  int n_operations; /* Symmetry operations from database */
  int (*rotations)[3][3];
  double (*translations)[3];
  int n_atoms;
  int *wyckoffs; /* Wyckoff letters */
  int *equivalent_atoms;
} SpglibDataset;

typedef struct {
  int size;
  int (*triplets)[3];
  int *weights;
  int mesh[3];
  int (*mesh_points)[3];
} SpglibTriplets;

SpglibDataset * spg_get_dataset(SPGCONST double lattice[3][3],
				SPGCONST double position[][3],
				const int types[],
				const int num_atom,
				const double symprec);

SpglibDataset * spgat_get_dataset(SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec,
				  const double angle_tolerance);

void spg_free_dataset(SpglibDataset *dataset);

/* Find symmetry operations. The operations are stored in */
/* ``rotatiion`` and ``translation``. The number of operations is */
/* return as the return value. Rotations and translations are */
/* given in fractional coordinates, and ``rotation[i]`` and */
/* ``translation[i]`` with same index give a symmetry oprations, */
/* i.e., these have to be used togather. */
int spg_get_symmetry(int rotation[][3][3],
		     double translation[][3],
		     const int max_size,
		     SPGCONST double lattice[3][3],
		     SPGCONST double position[][3],
		     const int types[],
		     const int num_atom,
		     const double symprec);

int spgat_get_symmetry(int rotation[][3][3],
		       double translation[][3],
		       const int max_size,
		       SPGCONST double lattice[3][3],
		       SPGCONST double position[][3],
		       const int types[],
		       const int num_atom,
		       const double symprec,
		       const double angle_tolerance);

/* Find symmetry operations with collinear spins on atoms. */
int spg_get_symmetry_with_collinear_spin(int rotation[][3][3],
					 double translation[][3],
					 const int max_size,
					 SPGCONST double lattice[3][3],
					 SPGCONST double position[][3],
					 const int types[],
					 const double spins[],
					 const int num_atom,
					 const double symprec);

int spgat_get_symmetry_with_collinear_spin(int rotation[][3][3],
					   double translation[][3],
					   const int max_size,
					   SPGCONST double lattice[3][3],
					   SPGCONST double position[][3],
					   const int types[],
					   const double spins[],
					   const int num_atom,
					   const double symprec,
					   const double angle_tolerance);

/* Return exact number of symmetry operations. This function may */
/* be used in advance to allocate memoery space for symmetry */
/* operations. Only upper bound is required, */
/* ``spg_get_max_multiplicity`` can be used instead of this */
/* function and ``spg_get_max_multiplicity`` is faster than this */
/* function. */
int spg_get_multiplicity(SPGCONST double lattice[3][3],
			 SPGCONST double position[][3],
			 const int types[],
			 const int num_atom,
			 const double symprec);

int spgat_get_multiplicity(SPGCONST double lattice[3][3],
			   SPGCONST double position[][3],
			   const int types[],
			   const int num_atom,
			   const double symprec,
			   const double angle_tolerance);

/* Considering periodicity of crystal, one of the possible smallest */
/* lattice is searched. The lattice is stored in ``smallest_lattice``. */
int spg_get_smallest_lattice(double smallest_lattice[3][3],
			     SPGCONST double lattice[3][3],
			     const double symprec);

/* A primitive cell is found from an input cell. Be careful that  */
/* ``lattice``, ``position``, and ``types`` are overwritten. */
/* ``num_atom`` is returned as return value. */
/* When any primitive cell is not found, 0 is returned. */
int spg_find_primitive(double lattice[3][3],
		       double position[][3],
		       int types[],
		       const int num_atom,
		       const double symprec);

int spgat_find_primitive(double lattice[3][3],
			 double position[][3],
			 int types[],
			 const int num_atom,
			 const double symprec,
			 const double angle_tolerance);

/* Space group is found in international table symbol (``symbol``) and */
/* number (return value). 0 is returned when it fails. */
int spg_get_international(char symbol[11],
			  SPGCONST double lattice[3][3],
			  SPGCONST double position[][3],
			  const int types[],
			  const int num_atom,
			  const double symprec);

int spgat_get_international(char symbol[11],
			    SPGCONST double lattice[3][3],
			    SPGCONST double position[][3],
			    const int types[],
			    const int num_atom,
			    const double symprec,
			    const double angle_tolerance);

/* Space group is found in schoenflies (``symbol``) and as number (return */
/* value).  0 is returned when it fails. */
int spg_get_schoenflies(char symbol[10],
			SPGCONST double lattice[3][3],
			SPGCONST double position[][3],
			const int types[],
			const int num_atom,
			const double symprec);

int spgat_get_schoenflies(char symbol[10],
			  SPGCONST double lattice[3][3],
			  SPGCONST double position[][3],
			  const int types[],
			  const int num_atom,
			  const double symprec,
			  const double angle_tolerance);

/* Point group symbol is obtained from the rotation part of */
/* symmetry operations */
int spg_get_pointgroup(char symbol[6],
		       int trans_mat[3][3],
		       SPGCONST int rotations[][3][3],
		       const int num_rotations);

/* Bravais lattice with internal atomic points are returned. */
/* The arrays are require to have 4 times larger memory space */
/* those of input cell. */
/* When bravais lattice could not be found, or could not be */
/* symmetrized, 0 is returned. */
int spg_refine_cell(double lattice[3][3],
		    double position[][3],
		    int types[],
		    const int num_atom,
		    const double symprec);

int spgat_refine_cell(double lattice[3][3],
		      double position[][3],
		      int types[],
		      const int num_atom,
		      const double symprec,
		      const double angle_tolerance);

/* Irreducible k-points are searched from the input k-points */
/* (``kpoints``).  The result is returned as a map of */
/* numbers (``map``), where ``kpoints`` and ``map`` have to have */
/* the same number of elements.  The array index of ``map`` */
/* corresponds to the reducible k-point numbering.  After finding */
/* irreducible k-points, the indices of the irreducible k-points */
/* are mapped to the elements of ``map``, i.e., number of unique */
/* values in ``map`` is the number of the irreducible k-points. */
/* The number of the irreducible k-points is also returned as the */
/* return value. */
int spg_get_ir_kpoints(int map[],
		       SPGCONST double kpoints[][3],
		       const int num_kpoints,
		       SPGCONST double lattice[3][3],
		       SPGCONST double position[][3],
		       const int types[],
		       const int num_atom,
		       const int is_time_reversal,
		       const double symprec);

/* Irreducible reciprocal grid points are searched from uniform */
/* mesh grid points specified by ``mesh`` and ``is_shift``. */
/* ``mesh`` stores three integers. Reciprocal primitive vectors */
/* are divided by the number stored in ``mesh`` with (0,0,0) point */
/* centering. The centering can be shifted only half of one mesh */
/* by setting 1 for each ``is_shift`` element. If 0 is set for */
/* ``is_shift``, it means there is no shift. This limitation of */
/* shifting enables the irreducible k-point search significantly */
/* faster when the mesh is very dense. */

/* The reducible uniform grid points are returned in reduced */
/* coordinates as ``grid_point``. A map between reducible and */
/* irreducible points are returned as ``map`` as in the indices of */
/* ``grid_point``. The number of the irreducible k-points are */
/* returned as the return value.  The time reversal symmetry is */
/* imposed by setting ``is_time_reversal`` 1. */
int spg_get_ir_reciprocal_mesh(int grid_point[][3],
			       int map[],
			       const int mesh[3],
			       const int is_shift[3],
			       const int is_time_reversal,
			       SPGCONST double lattice[3][3],
			       SPGCONST double position[][3],
			       const int types[],
			       const int num_atom,
			       const double symprec);

/* The irreducible k-points are searched from unique k-point mesh */
/* grids from real space lattice vectors and rotation matrices of */
/* symmetry operations in real space with stabilizers. The */
/* stabilizers are written in reduced coordinates. Number of the */
/* stabilizers are given by ``num_q``. Reduced k-points are stored */
/* in ``map`` as indices of ``grid_point``. The number of the */
/* reduced k-points with stabilizers are returned as the return */
/* value. */
int spg_get_stabilized_reciprocal_mesh(int grid_point[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       SPGCONST double lattice[3][3],
				       const int num_rot,
				       SPGCONST int rotations[][3][3],
				       const int num_q,
				       SPGCONST double qpoints[][3],
				       const double symprec);

/* Irreducible triplets of k-points are searched under conservation of */
/* :math:``\mathbf{k}_1 + \mathbf{k}_2 + \mathbf{k}_3 = \mathbf{G}``. */
/* Don't forget to free memory space of triplets using spg_free_triplets */
SpglibTriplets * spg_get_triplets_reciprocal_mesh(const int mesh[3],
						  const int is_time_reversal,
						  SPGCONST double lattice[3][3],
						  const int num_rot,
						  SPGCONST int rotations[][3][3],
						  const double symprec);

void spg_free_triplets(SpglibTriplets * triplets);

int spg_get_triplets_reciprocal_mesh_at_q(int weights[],
					  int grid_points[][3],
					  int third_q[],
					  const int grid_point,
					  const int mesh[3],
					  const int is_time_reversal,
					  SPGCONST double lattice[3][3],
					  const int num_rot,
					  SPGCONST int rotations[][3][3],
					  const double symprec);

int spg_extract_triplets_reciprocal_mesh_at_q(int triplets_at_q[][3],
					      int weight_triplets_at_q[],
					      const int fixed_grid_number,
					      const int num_triplets,
					      SPGCONST int triplets[][3],
					      const int mesh[3],
					      const int is_time_reversal,
					      SPGCONST double lattice[3][3],
					      const int num_rot,
					      SPGCONST int rotations[][3][3],
					      const double symprec);

#endif
