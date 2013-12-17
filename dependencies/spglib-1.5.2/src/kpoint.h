/* kpoints.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __kpoints_H__
#define __kpoints_H__

#include "symmetry.h"
#include "mathfunc.h"

int kpt_get_irreducible_kpoints(int map[],
				SPGCONST double kpoints[][3], 
				const int num_kpoint,
				const Symmetry * symmetry,
				const int is_time_reversal,
				const double symprec);
int kpt_get_irreducible_reciprocal_mesh(int grid_address[][3],
					int map[],
					const int mesh[3],
					const int is_shift[3],
					const int is_time_reversal,
					const Symmetry * symmetry);
int kpt_get_stabilized_reciprocal_mesh(int grid_address[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const MatINT * pointgroup_real,
				       const int num_q,
				       SPGCONST double qpoints[][3]);
int kpt_relocate_BZ_grid_address(int bz_grid_address[][3],
				 int bz_map[],
				 SPGCONST int grid_address[][3],
				 const int mesh[3],
				 SPGCONST double rec_lattice[3][3],
				 const int is_shift[3]);
int kpt_get_ir_triplets_at_q(int weights[],
			     int grid_address[][3],
			     int third_q[],
			     const int grid_point,
			     const int mesh[3],
			     const int is_time_reversal,
			     const MatINT * rotations);
int kpt_get_BZ_triplets_at_q(int triplets[][3],
			     const int grid_point,
			     SPGCONST int bz_grid_address[][3],
			     const int bz_map[],
			     const int weights[],
			     const int mesh[3]);
#endif
