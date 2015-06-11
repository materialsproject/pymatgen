/* triplet_kpoint.h */
/* Copyright (C) 2015 Atsushi Togo */

#ifndef __triplet_kpoint_H__
#define __triplet_kpoint_H__

#include "mathfunc.h"

int tpk_get_ir_triplets_at_q(int map_triplets[],
			     int map_q[],
			     int grid_address[][3],
			     const int grid_point,
			     const int mesh[3],
			     const int is_time_reversal,
			     const MatINT * rotations);
int tpk_get_BZ_triplets_at_q(int triplets[][3],
			     const int grid_point,
			     SPGCONST int bz_grid_address[][3],
			     const int bz_map[],
			     const int map_triplets[],
			     const int num_map_triplets,
			     const int mesh[3]);

#endif
