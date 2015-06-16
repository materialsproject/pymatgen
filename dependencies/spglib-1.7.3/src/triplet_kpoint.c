/* triplet_kpoint.c */
/* Copyright (C) 2015 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include "mathfunc.h"
#include "kpoint.h"
#include "triplet_kpoint.h"


static void grid_point_to_address_double(int address_double[3],
					 const int grid_point,
					 const int mesh[3],
					 const int is_shift[3]);
static int get_ir_triplets_at_q(int map_triplets[],
				int map_q[],
				int grid_address[][3],
				const int grid_point,
				const int mesh[3],
				const MatINT * rot_reciprocal);
static int get_BZ_triplets_at_q(int triplets[][3],
				const int grid_point,
				SPGCONST int bz_grid_address[][3],
				const int bz_map[],
				const int map_triplets[],
				const int num_map_triplets,
				const int mesh[3]);
static int get_third_q_of_triplets_at_q(int bz_address[3][3],
					const int q_index,
					const int bz_map[],
					const int mesh[3],
					const int bzmesh[3]);

int tpk_get_ir_triplets_at_q(int map_triplets[],
			     int map_q[],
			     int grid_address[][3],
			     const int grid_point,
			     const int mesh[3],
			     const int is_time_reversal,
			     const MatINT * rotations)
{
  int num_ir;
  MatINT *rot_reciprocal;

  rot_reciprocal = kpt_get_point_group_reciprocal(rotations, is_time_reversal);
  num_ir = get_ir_triplets_at_q(map_triplets,
				map_q,
				grid_address,
				grid_point,
				mesh,
				rot_reciprocal);
  mat_free_MatINT(rot_reciprocal);
  return num_ir;
}

int tpk_get_BZ_triplets_at_q(int triplets[][3],
			     const int grid_point,
			     SPGCONST int bz_grid_address[][3],
			     const int bz_map[],
			     const int map_triplets[],
			     const int num_map_triplets,
			     const int mesh[3])
{
  return get_BZ_triplets_at_q(triplets,
			      grid_point,
			      bz_grid_address,
			      bz_map,
			      map_triplets,
			      num_map_triplets,
			      mesh);
}

static int get_ir_triplets_at_q(int map_triplets[],
				int map_q[],
				int grid_address[][3],
				const int grid_point,
				const int mesh[3],
				const MatINT * rot_reciprocal)
{
  int i, j, num_grid, q_2, num_ir_q, num_ir_triplets, ir_grid_point;
  int mesh_double[3], is_shift[3];
  int address_double0[3], address_double1[3], address_double2[3];
  int *ir_grid_points, *third_q;
  double tolerance;
  double stabilizer_q[1][3];
  MatINT *rot_reciprocal_q;

  tolerance = 0.01 / (mesh[0] + mesh[1] + mesh[2]);
  num_grid = mesh[0] * mesh[1] * mesh[2];

  for (i = 0; i < 3; i++) {
    /* Only consider the gamma-point */
    is_shift[i] = 0;
    mesh_double[i] = mesh[i] * 2;
  }

  /* Search irreducible q-points (map_q) with a stabilizer */
  /* q */  
  grid_point_to_address_double(address_double0, grid_point, mesh, is_shift);
  for (i = 0; i < 3; i++) {
    stabilizer_q[0][i] =
      (double)address_double0[i] / mesh_double[i] - (address_double0[i] > mesh[i]);
  }

  rot_reciprocal_q = kpt_get_point_group_reciprocal_with_q(rot_reciprocal,
							   tolerance,
							   1,
							   stabilizer_q);
  num_ir_q = kpt_get_irreducible_reciprocal_mesh(grid_address,
						 map_q,
						 mesh,
						 is_shift,
						 rot_reciprocal_q);
  mat_free_MatINT(rot_reciprocal_q);

  third_q = (int*) malloc(sizeof(int) * num_ir_q);
  ir_grid_points = (int*) malloc(sizeof(int) * num_ir_q);
  num_ir_q = 0;
  for (i = 0; i < num_grid; i++) {
    if (map_q[i] == i) {
      ir_grid_points[num_ir_q] = i;
      num_ir_q++;
    }
    map_triplets[i] = -1;
  }

#pragma omp parallel for private(j, address_double1, address_double2)
  for (i = 0; i < num_ir_q; i++) {
    grid_point_to_address_double(address_double1,
				 ir_grid_points[i],
				 mesh,
				 is_shift); /* q' */
    for (j = 0; j < 3; j++) { /* q'' */
      address_double2[j] = - address_double0[j] - address_double1[j];
    }
    third_q[i] = kpt_get_grid_point_double_mesh(address_double2, mesh);
  }

  num_ir_triplets = 0;
  for (i = 0; i < num_ir_q; i++) {
    ir_grid_point = ir_grid_points[i];
    q_2 = third_q[i];
    if (map_triplets[map_q[q_2]] > -1) {
      map_triplets[ir_grid_point] = map_q[q_2];
    } else {
      map_triplets[ir_grid_point] = ir_grid_point;
      num_ir_triplets++;
    }
  }

#pragma omp parallel for
  for (i = 0; i < num_grid; i++) {
    map_triplets[i] = map_triplets[map_q[i]];
  }
  
  free(third_q);
  third_q = NULL;
  free(ir_grid_points);
  ir_grid_points = NULL;

  return num_ir_triplets;
}

static int get_BZ_triplets_at_q(int triplets[][3],
				const int grid_point,
				SPGCONST int bz_grid_address[][3],
				const int bz_map[],
				const int map_triplets[],
				const int num_map_triplets,
				const int mesh[3])
{
  int i, j, k, num_ir;
  int bz_address[3][3], bz_address_double[3], bzmesh[3];
  int *ir_grid_points;

  for (i = 0; i < 3; i++) {
    bzmesh[i] = mesh[i] * 2;
  }

  num_ir = 0;
  ir_grid_points = (int*) malloc(sizeof(int) * num_map_triplets);
  for (i = 0; i < num_map_triplets; i++) {
    if (map_triplets[i] == i) {
      ir_grid_points[num_ir] = i;
      num_ir++;
    }
  }
 
#pragma omp parallel for private(j, k, bz_address, bz_address_double)
  for (i = 0; i < num_ir; i++) {
    for (j = 0; j < 3; j++) {
      bz_address[0][j] = bz_grid_address[grid_point][j];
      bz_address[1][j] = bz_grid_address[ir_grid_points[i]][j];
      bz_address[2][j] = - bz_address[0][j] - bz_address[1][j];
    }
    for (j = 2; j > -1; j--) {
      if (get_third_q_of_triplets_at_q(bz_address,
    				       j,
    				       bz_map,
    				       mesh,
    				       bzmesh) == 0) {
    	break;
      }
    }
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	bz_address_double[k] = bz_address[j][k] * 2;
      }
      triplets[i][j] =
	bz_map[kpt_get_grid_point_double_mesh(bz_address_double, bzmesh)];
    }
  }

  free(ir_grid_points);
  
  return num_ir;
}

static int get_third_q_of_triplets_at_q(int bz_address[3][3],
					const int q_index,
					const int bz_map[],
					const int mesh[3],
					const int bzmesh[3])
{
  int i, j, smallest_g, smallest_index, sum_g, delta_g[3];
  int bzgp[KPT_NUM_BZ_SEARCH_SPACE], bz_address_double[3];

  mat_modulo_i3(bz_address[q_index], mesh);
  for (i = 0; i < 3; i++) {
    delta_g[i] = 0;
    for (j = 0; j < 3; j++) {
      delta_g[i] += bz_address[j][i];
    }
    delta_g[i] /= mesh[i];
  }
  
  for (i = 0; i < KPT_NUM_BZ_SEARCH_SPACE; i++) {
    for (j = 0; j < 3; j++) {
      bz_address_double[j] = (bz_address[q_index][j] +
			      kpt_bz_search_space[i][j] * mesh[j]) * 2;
    }
    bzgp[i] = bz_map[kpt_get_grid_point_double_mesh(bz_address_double, bzmesh)];
  }

  for (i = 0; i < KPT_NUM_BZ_SEARCH_SPACE; i++) {
    if (bzgp[i] != -1) {
      goto escape;
    }
  }

 escape:

  smallest_g = 4;
  smallest_index = 0;

  for (i = 0; i < KPT_NUM_BZ_SEARCH_SPACE; i++) {
    if (bzgp[i] > -1) { /* q'' is in BZ */
      sum_g = (abs(delta_g[0] + kpt_bz_search_space[i][0]) +
	       abs(delta_g[1] + kpt_bz_search_space[i][1]) +
	       abs(delta_g[2] + kpt_bz_search_space[i][2]));
      if (sum_g < smallest_g) {
	smallest_index = i;
	smallest_g = sum_g;
      }
    }
  }

  for (i = 0; i < 3; i++) {
    bz_address[q_index][i] += kpt_bz_search_space[smallest_index][i] * mesh[i];
  }

  return smallest_g;
}

static void grid_point_to_address_double(int address_double[3],
					 const int grid_point,
					 const int mesh[3],
					 const int is_shift[3])
{
  int i;
  int address[3];

#ifndef GRID_ORDER_XYZ
  address[2] = grid_point / (mesh[0] * mesh[1]);
  address[1] = (grid_point - address[2] * mesh[0] * mesh[1]) / mesh[0];
  address[0] = grid_point % mesh[0];
#else
  address[0] = grid_point / (mesh[1] * mesh[2]);
  address[1] = (grid_point - address[0] * mesh[1] * mesh[2]) / mesh[2];
  address[2] = grid_point % mesh[2];
#endif

  for (i = 0; i < 3; i++) {
    address_double[i] = address[i] * 2 + is_shift[i];
  }
}

