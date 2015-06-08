/* kpoint.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __kpoint_H__
#define __kpoint_H__

#include "mathfunc.h"

/* #define GRID_ORDER_XYZ */
/* The addressing order of mesh grid is defined as running left */
/* element first. But when GRID_ORDER_XYZ is defined, it is changed to right */ 
/* element first. */

#define KPT_NUM_BZ_SEARCH_SPACE 125
extern const int kpt_bz_search_space[KPT_NUM_BZ_SEARCH_SPACE][3];

int kpt_get_grid_point_double_mesh(const int address_double[3],
				   const int mesh[3]);
int kpt_get_irreducible_reciprocal_mesh(int grid_address[][3],
					int map[],
					const int mesh[3],
					const int is_shift[3],
					const MatINT *rot_reciprocal);
int kpt_get_stabilized_reciprocal_mesh(int grid_address[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const MatINT * rotations,
				       const int num_q,
				       SPGCONST double qpoints[][3]);
void kpt_get_grid_points_by_rotations(int rot_grid_points[],
				      const int address_orig[3],
				      const MatINT * rot_reciprocal,
				      const int mesh[3],
				      const int is_shift[3]);
void kpt_get_BZ_grid_points_by_rotations(int rot_grid_points[],
					 const int address_orig[3],
					 const MatINT * rot_reciprocal,
					 const int mesh[3],
					 const int is_shift[3],
					 const int bz_map[]);
int kpt_relocate_BZ_grid_address(int bz_grid_address[][3],
				 int bz_map[],
				 SPGCONST int grid_address[][3],
				 const int mesh[3],
				 SPGCONST double rec_lattice[3][3],
				 const int is_shift[3]);
MatINT *kpt_get_point_group_reciprocal(const MatINT * rotations,
				       const int is_time_reversal);
MatINT *kpt_get_point_group_reciprocal_with_q(const MatINT * rot_reciprocal,
					      const double symprec,
					      const int num_q,
					      SPGCONST double qpoints[][3]);

#endif
