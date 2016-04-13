/* Copyright (C) 2008 Atsushi Togo */
/* All rights reserved. */

/* This file is part of spglib. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the phonopy project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#ifndef __spglib_H__
#define __spglib_H__

/* SPGCONST is used instead of 'const' so to avoid gcc warning. */
/* However there should be better way than this way.... */
#ifndef SPGCONST
#define SPGCONST
#endif

int spg_get_major_version(void);
int spg_get_minor_version(void);
int spg_get_micro_version(void);

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
  char setting[6];
  double transformation_matrix[3][3];
  double origin_shift[3];
  int n_operations;
  int (*rotations)[3][3];
  double (*translations)[3];
  int n_atoms;
  int *wyckoffs;
  int *equivalent_atoms;
  int n_std_atoms;
  double std_lattice[3][3];
  int *std_types;
  double (*std_positions)[3];
  int pointgroup_number;
  char pointgroup_symbol[6];
} SpglibDataset;

typedef struct {
  int number;
  char schoenflies[7];
  char hall_symbol[17];
  char international[32];
  char international_full[20];
  char international_short[11];
  int pointgroup_number;
} SpglibSpacegroupType;

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

/* hall_number = 0 gives the same as spg_get_dataset. */
SpglibDataset * spg_get_dataset_with_hall_number(SPGCONST double lattice[3][3],
						 SPGCONST double position[][3],
						 const int types[],
						 const int num_atom,
						 const int hall_number,
						 const double symprec);

/* hall_number = 0 gives the same as spgat_get_dataset. */
SpglibDataset *
spgat_get_dataset_with_hall_number(SPGCONST double lattice[3][3],
				   SPGCONST double position[][3],
				   const int types[],
				   const int num_atom,
				   const int hall_number,
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

/* This is only used to check consistensy with spg_get_symmetry. */
int spg_get_symmetry_numerical(int rotation[][3][3],
			       double translation[][3],
			       const int max_size,
			       SPGCONST double lattice[3][3],
			       SPGCONST double position[][3],
			       const int types[],
			       const int num_atom,
			       const double symprec);

int spgat_get_symmetry_numerical(int rotation[][3][3],
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
					 int equivalent_atoms[],
					 const int max_size,
					 SPGCONST double lattice[3][3],
					 SPGCONST double position[][3],
					 const int types[],
					 const double spins[],
					 const int num_atom,
					 const double symprec);

int spgat_get_symmetry_with_collinear_spin(int rotation[][3][3],
					   double translation[][3],
					   int equivalent_atoms[],
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
/* operations. */
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

/* Space-group operations in built-in database are accessed by index */
/* of hall symbol. The index is defined as number from 1 to 530. */
/* The muximum number of symmetry operations is 192. */
int spg_get_symmetry_from_database(int rotations[192][3][3],
				   double translations[192][3],
				   const int hall_number);

/* Space-group type information is accessed by index of hall symbol. */
/* The index is defined as number from 1 to 530. */
SpglibSpacegroupType spg_get_spacegroup_type(const int hall_number);


int spg_standardize_cell(double lattice[3][3],
			 double position[][3],
			 int types[],
			 const int num_atom,
			 const int to_primitive,
			 const int no_idealize,
			 const double symprec);

int spgat_standardize_cell(double lattice[3][3],
			   double position[][3],
			   int types[],
			   const int num_atom,
			   const int to_primitive,
			   const int no_idealize,
			   const double symprec,
			   const double angle_tolerance);

/************/
/* Obsolete */
/************/
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

/************/
/* Obsolete */
/************/
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


/*---------*/
/* kpoints */
/*---------*/

/* Translate grid address to grid point index in the kspclib definition */
/* (see the comment in kgrid.h.) */
/* A q-point in fractional coordinates is given as */
/* ((grid_address * 2 + (shift != 0)) / (mesh * 2)). */
/* Each element of shift[] is 0 or non-zero. */
int spg_get_grid_point_from_address(const int grid_address[3],
				    const int mesh[3]);

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
/* coordinates as ``grid_address``. A map between reducible and */
/* irreducible points are returned as ``map`` as in the indices of */
/* ``grid_address``. The number of the irreducible k-points are */
/* returned as the return value.  The time reversal symmetry is */
/* imposed by setting ``is_time_reversal`` 1. */
int spg_get_ir_reciprocal_mesh(int grid_address[][3],
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
/* in ``map`` as indices of ``grid_address``. The number of the */
/* reduced k-points with stabilizers are returned as the return */
/* value. */
int spg_get_stabilized_reciprocal_mesh(int grid_address[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const int num_rot,
				       SPGCONST int rotations[][3][3],
				       const int num_q,
				       SPGCONST double qpoints[][3]);

/* Rotation operations in reciprocal space ``rot_reciprocal`` are applied */
/* to a grid address ``address_orig`` and resulting grid points are stored in */
/* ``rot_grid_points``. Return 0 if failed. */
int spg_get_grid_points_by_rotations(int rot_grid_points[],
				     const int address_orig[3],
				     const int num_rot,
				     SPGCONST int rot_reciprocal[][3][3],
				     const int mesh[3],
				     const int is_shift[3]);

int spg_get_BZ_grid_points_by_rotations(int rot_grid_points[],
					const int address_orig[3],
					const int num_rot,
					SPGCONST int rot_reciprocal[][3][3],
					const int mesh[3],
					const int is_shift[3],
					const int bz_map[]);

/* Grid addresses are relocated inside Brillouin zone. */
/* Number of ir-grid-points inside Brillouin zone is returned. */
/* It is assumed that the following arrays have the shapes of */
/*   bz_grid_address[prod(mesh + 1)][3] */
/*   bz_map[prod(mesh * 2)] */
/* where grid_address[prod(mesh)][3]. */
/* Each element of grid_address is mapped to each element of */
/* bz_grid_address with keeping element order. bz_grid_address has */
/* larger memory space to represent BZ surface even if some points */
/* on a surface are translationally equivalent to the other points */
/* on the other surface. Those equivalent points are added successively */
/* as grid point numbers to bz_grid_address. Those added grid points */
/* are stored after the address of end point of grid_address, i.e. */
/*                                                                       */
/* |-----------------array size of bz_grid_address---------------------| */
/* |--grid addresses similar to grid_address--|--newly added ones--|xxx| */
/*                                                                       */
/* where xxx means the memory space that may not be used. Number of grid */
/* points stored in bz_grid_address is returned. */
/* bz_map is used to recover grid point index expanded to include BZ */
/* surface from grid address. The grid point indices are mapped to */
/* (mesh[0] * 2) x (mesh[1] * 2) x (mesh[2] * 2) space (bz_map). */
int spg_relocate_BZ_grid_address(int bz_grid_address[][3],
				 int bz_map[],
				 SPGCONST int grid_address[][3],
				 const int mesh[3],
				 SPGCONST double rec_lattice[3][3],
				 const int is_shift[3]);

void spg_get_neighboring_grid_points(int relative_grid_points[],
				     const int grid_point,
				     SPGCONST int relative_grid_address[][3],
				     const int num_relative_grid_address,
				     const int mesh[3],
				     SPGCONST int bz_grid_address[][3],
				     const int bz_map[]);

/*--------*/
/* Niggli */
/*--------*/
/* Return 0 if failed */
int spg_niggli_reduce(double lattice[3][3], const double symprec);

#endif

