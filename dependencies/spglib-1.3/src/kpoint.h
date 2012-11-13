/* kpoints.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __kpoints_H__
#define __kpoints_H__

#include "symmetry.h"
#include "mathfunc.h"

typedef struct {
  int size;
  int (*triplets)[3];
  int *weights;
  int mesh[3];
  int (*mesh_points)[3];
} Triplets;

int kpt_get_irreducible_kpoints( int map[],
				 SPGCONST double kpoints[][3], 
				 const int num_kpoint,
				 SPGCONST double lattice[3][3],
				 const Symmetry * symmetry,
				 const int is_time_reversal,
				 const double symprec );
int kpt_get_irreducible_reciprocal_mesh( int grid_points[][3],
					 int map[],
					 const int mesh[3],
					 const int is_shift[3],
					 const int is_time_reversal,
					 SPGCONST double lattice[3][3],
					 const Symmetry * symmetry,
					 const double symprec );
int kpt_get_stabilized_reciprocal_mesh( int grid_points[][3],
					int map[],
					const int mesh[3],
					const int is_shift[3],
					const int is_time_reversal,
					SPGCONST double lattice[3][3],
					const MatINT * pointgroup_real,
					const int num_q,
					SPGCONST double qpoints[][3],
					const double symprec );
Triplets * kpt_get_triplets_reciprocal_mesh( const int mesh[3],
					     const int is_time_reversal,
					     SPGCONST double lattice[3][3],
					     const MatINT * pointgroup_real,
					     const double symprec );
void kpt_free_triplets( Triplets * t );
int kpt_get_ir_triplets_at_q( int weights[],
			      int grid_points[][3],
			      int third_q[],
			      const int grid_point,
			      const int mesh[3],
			      const int is_time_reversal,
			      SPGCONST double lattice[3][3],
			      const MatINT * rotations,
			      const double symprec );
int kpt_extract_triplets_reciprocal_mesh_at_q( int triplets_at_q[][3],
					       int weight_at_q[],
					       const int fixed_grid_number,
					       const int num_triplets,
					       SPGCONST int triplets[][3],
					       const int mesh[3],
					       const int is_time_reversal,
					       SPGCONST double lattice[3][3],
					       const MatINT * pointgroup_real,
					       const double symprec );
#endif
