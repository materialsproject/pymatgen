/* kpoint.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include "mathfunc.h"
#include "kpoint.h"

const int kpt_bz_search_space[KPT_NUM_BZ_SEARCH_SPACE][3] = {
  { 0,  0,  0},
  { 0,  0,  1},
  { 0,  0,  2},
  { 0,  0, -2},
  { 0,  0, -1},
  { 0,  1,  0},
  { 0,  1,  1},
  { 0,  1,  2},
  { 0,  1, -2},
  { 0,  1, -1},
  { 0,  2,  0},
  { 0,  2,  1},
  { 0,  2,  2},
  { 0,  2, -2},
  { 0,  2, -1},
  { 0, -2,  0},
  { 0, -2,  1},
  { 0, -2,  2},
  { 0, -2, -2},
  { 0, -2, -1},
  { 0, -1,  0},
  { 0, -1,  1},
  { 0, -1,  2},
  { 0, -1, -2},
  { 0, -1, -1},
  { 1,  0,  0},
  { 1,  0,  1},
  { 1,  0,  2},
  { 1,  0, -2},
  { 1,  0, -1},
  { 1,  1,  0},
  { 1,  1,  1},
  { 1,  1,  2},
  { 1,  1, -2},
  { 1,  1, -1},
  { 1,  2,  0},
  { 1,  2,  1},
  { 1,  2,  2},
  { 1,  2, -2},
  { 1,  2, -1},
  { 1, -2,  0},
  { 1, -2,  1},
  { 1, -2,  2},
  { 1, -2, -2},
  { 1, -2, -1},
  { 1, -1,  0},
  { 1, -1,  1},
  { 1, -1,  2},
  { 1, -1, -2},
  { 1, -1, -1},
  { 2,  0,  0},
  { 2,  0,  1},
  { 2,  0,  2},
  { 2,  0, -2},
  { 2,  0, -1},
  { 2,  1,  0},
  { 2,  1,  1},
  { 2,  1,  2},
  { 2,  1, -2},
  { 2,  1, -1},
  { 2,  2,  0},
  { 2,  2,  1},
  { 2,  2,  2},
  { 2,  2, -2},
  { 2,  2, -1},
  { 2, -2,  0},
  { 2, -2,  1},
  { 2, -2,  2},
  { 2, -2, -2},
  { 2, -2, -1},
  { 2, -1,  0},
  { 2, -1,  1},
  { 2, -1,  2},
  { 2, -1, -2},
  { 2, -1, -1},
  {-2,  0,  0},
  {-2,  0,  1},
  {-2,  0,  2},
  {-2,  0, -2},
  {-2,  0, -1},
  {-2,  1,  0},
  {-2,  1,  1},
  {-2,  1,  2},
  {-2,  1, -2},
  {-2,  1, -1},
  {-2,  2,  0},
  {-2,  2,  1},
  {-2,  2,  2},
  {-2,  2, -2},
  {-2,  2, -1},
  {-2, -2,  0},
  {-2, -2,  1},
  {-2, -2,  2},
  {-2, -2, -2},
  {-2, -2, -1},
  {-2, -1,  0},
  {-2, -1,  1},
  {-2, -1,  2},
  {-2, -1, -2},
  {-2, -1, -1},
  {-1,  0,  0},
  {-1,  0,  1},
  {-1,  0,  2},
  {-1,  0, -2},
  {-1,  0, -1},
  {-1,  1,  0},
  {-1,  1,  1},
  {-1,  1,  2},
  {-1,  1, -2},
  {-1,  1, -1},
  {-1,  2,  0},
  {-1,  2,  1},
  {-1,  2,  2},
  {-1,  2, -2},
  {-1,  2, -1},
  {-1, -2,  0},
  {-1, -2,  1},
  {-1, -2,  2},
  {-1, -2, -2},
  {-1, -2, -1},
  {-1, -1,  0},
  {-1, -1,  1},
  {-1, -1,  2},
  {-1, -1, -2},
  {-1, -1, -1}
};

static MatINT *get_point_group_reciprocal(const MatINT * rotations,
					  const int is_time_reversal);
static MatINT *get_point_group_reciprocal_with_q(const MatINT * rot_reciprocal,
						 const double symprec,
						 const int num_q,
						 SPGCONST double qpoints[][3]);
static int get_ir_reciprocal_mesh(int grid_address[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  const MatINT * rot_reciprocal);
static int
get_ir_reciprocal_mesh_openmp(int grid_address[][3],
			      int map[],
			      const int mesh[3],
			      const int is_shift[3],
			      const MatINT* rot_reciprocal);
static int relocate_BZ_grid_address(int bz_grid_address[][3],
				    int bz_map[],
				    SPGCONST int grid_address[][3],
				    const int mesh[3],
				    SPGCONST double rec_lattice[3][3],
				    const int is_shift[3]);
static double get_tolerance_for_BZ_reduction(SPGCONST double rec_lattice[3][3],
					     const int mesh[3]);
static int get_grid_point_double_mesh(const int address_double[3],
				      const int mesh[3]);
static int get_grid_point_single_mesh(const int address[3],
				      const int mesh[3]);
static void reduce_grid_address(int address[3],
				const int address_double[3],
				const int mesh[3]);

int kpt_get_grid_point_double_mesh(const int address_double[3],
				   const int mesh[3])
{
  return get_grid_point_double_mesh(address_double, mesh);
}

/* grid_address (e.g. 4x4x4 mesh, unless GRID_ORDER_XYZ is defined) */
/*    [[ 0  0  0]                                                   */
/*     [ 1  0  0]                                                   */
/*     [ 2  0  0]                                                   */
/*     [-1  0  0]                                                   */
/*     [ 0  1  0]                                                   */
/*     [ 1  1  0]                                                   */
/*     [ 2  1  0]                                                   */
/*     [-1  1  0]                                                   */
/*     ....      ]                                                  */
/*                                                                  */
/* Each value of 'map' correspnds to the index of grid_point.       */
int kpt_get_irreducible_reciprocal_mesh(int grid_address[][3],
					int map[],
					const int mesh[3],
					const int is_shift[3],
					const MatINT *rot_reciprocal)
{
  int num_ir;

#ifdef _OPENMP
  num_ir = get_ir_reciprocal_mesh_openmp(grid_address,
					 map,
					 mesh,
					 is_shift,
					 rot_reciprocal);
  
#else
  num_ir = get_ir_reciprocal_mesh(grid_address,
				  map,
				  mesh,
				  is_shift,
				  rot_reciprocal);
#endif
  
  return num_ir;
}

int kpt_get_stabilized_reciprocal_mesh(int grid_address[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const MatINT * rotations,
				       const int num_q,
				       SPGCONST double qpoints[][3])
{
  int num_ir;
  MatINT *rot_reciprocal, *rot_reciprocal_q;
  double tolerance;
  
  rot_reciprocal = get_point_group_reciprocal(rotations, is_time_reversal);
  tolerance = 0.01 / (mesh[0] + mesh[1] + mesh[2]);
  rot_reciprocal_q = get_point_group_reciprocal_with_q(rot_reciprocal,
						       tolerance,
						       num_q,
						       qpoints);

#ifdef _OPENMP
  num_ir = get_ir_reciprocal_mesh_openmp(grid_address,
					 map,
					 mesh,
					 is_shift,
					 rot_reciprocal_q);
#else
  num_ir = get_ir_reciprocal_mesh(grid_address,
				  map,
				  mesh,
				  is_shift,
				  rot_reciprocal_q);
#endif
  mat_free_MatINT(rot_reciprocal_q);
  mat_free_MatINT(rot_reciprocal);
  return num_ir;
}

void kpt_get_grid_points_by_rotations(int rot_grid_points[],
				      const int address_orig[3],
				      const MatINT * rot_reciprocal,
				      const int mesh[3],
				      const int is_shift[3])
{
  int i;
  int address_double_orig[3], address_double[3];

  for (i = 0; i < 3; i++) {
    address_double_orig[i] = address_orig[i] * 2 + is_shift[i];
  }
  for (i = 0; i < rot_reciprocal->size; i++) {
    mat_multiply_matrix_vector_i3(address_double,
				  rot_reciprocal->mat[i],
				  address_double_orig);
    rot_grid_points[i] = get_grid_point_double_mesh(address_double, mesh);
  }
}

void kpt_get_BZ_grid_points_by_rotations(int rot_grid_points[],
					 const int address_orig[3],
					 const MatINT * rot_reciprocal,
					 const int mesh[3],
					 const int is_shift[3],
					 const int bz_map[])
{
  int i;
  int address_double_orig[3], address_double[3], bzmesh[3];

  for (i = 0; i < 3; i++) {
    bzmesh[i] = mesh[i] * 2;
    address_double_orig[i] = address_orig[i] * 2 + is_shift[i];
  }
  for (i = 0; i < rot_reciprocal->size; i++) {
    mat_multiply_matrix_vector_i3(address_double,
				  rot_reciprocal->mat[i],
				  address_double_orig);
    rot_grid_points[i] =
      bz_map[get_grid_point_double_mesh(address_double, bzmesh)];
  }
}

int kpt_relocate_BZ_grid_address(int bz_grid_address[][3],
				 int bz_map[],
				 SPGCONST int grid_address[][3],
				 const int mesh[3],
				 SPGCONST double rec_lattice[3][3],
				 const int is_shift[3])
{
  return relocate_BZ_grid_address(bz_grid_address,
				  bz_map,
				  grid_address,
				  mesh,
				  rec_lattice,
				  is_shift);
}


MatINT *kpt_get_point_group_reciprocal(const MatINT * rotations,
				       const int is_time_reversal)
{
  return get_point_group_reciprocal(rotations, is_time_reversal);
}

MatINT *kpt_get_point_group_reciprocal_with_q(const MatINT * rot_reciprocal,
					      const double symprec,
					      const int num_q,
					      SPGCONST double qpoints[][3])
{
  return get_point_group_reciprocal_with_q(rot_reciprocal,
					   symprec,
					   num_q,
					   qpoints);
}

static MatINT *get_point_group_reciprocal(const MatINT * rotations,
					  const int is_time_reversal)
{
  int i, j, num_rot;
  MatINT *rot_reciprocal, *rot_return;
  int *unique_rot;
  SPGCONST int inversion[3][3] = {
    {-1, 0, 0 },
    { 0,-1, 0 },
    { 0, 0,-1 }
  };
  
  if (is_time_reversal) {
    rot_reciprocal = mat_alloc_MatINT(rotations->size * 2);
  } else {
    rot_reciprocal = mat_alloc_MatINT(rotations->size);
  }
  unique_rot = (int*)malloc(sizeof(int) * rot_reciprocal->size);
  for (i = 0; i < rot_reciprocal->size; i++) {
    unique_rot[i] = -1;
  }

  for (i = 0; i < rotations->size; i++) {
    mat_transpose_matrix_i3(rot_reciprocal->mat[i], rotations->mat[i]);
    
    if (is_time_reversal) {
      mat_multiply_matrix_i3(rot_reciprocal->mat[rotations->size+i],
			     inversion,
			     rot_reciprocal->mat[i]);
    }
  }

  num_rot = 0;
  for (i = 0; i < rot_reciprocal->size; i++) {
    for (j = 0; j < num_rot; j++) {
      if (mat_check_identity_matrix_i3(rot_reciprocal->mat[unique_rot[j]],
				       rot_reciprocal->mat[i])) {
	goto escape;
      }
    }
    unique_rot[num_rot] = i;
    num_rot++;
  escape:
    ;
  }

  rot_return = mat_alloc_MatINT(num_rot);
  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot_return->mat[i], rot_reciprocal->mat[unique_rot[i]]);    }
  free(unique_rot);
  mat_free_MatINT(rot_reciprocal);

  return rot_return;
}

static MatINT *get_point_group_reciprocal_with_q(const MatINT * rot_reciprocal,
						 const double symprec,
						 const int num_q,
						 SPGCONST double qpoints[][3])
{
  int i, j, k, l, is_all_ok, num_rot;
  int *ir_rot;
  double q_rot[3], diff[3];
  MatINT * rot_reciprocal_q;

  is_all_ok = 0;
  num_rot = 0;
  ir_rot = (int*)malloc(sizeof(int) * rot_reciprocal->size);
  for (i = 0; i < rot_reciprocal->size; i++) {
    ir_rot[i] = -1;
  }
  for (i = 0; i < rot_reciprocal->size; i++) {
    for (j = 0; j < num_q; j++) {
      is_all_ok = 0;
      mat_multiply_matrix_vector_id3(q_rot,
				     rot_reciprocal->mat[i],
				     qpoints[j]);

      for (k = 0; k < num_q; k++) {
	for (l = 0; l < 3; l++) {
	  diff[l] = q_rot[l] - qpoints[k][l];
	  diff[l] -= mat_Nint(diff[l]);
	}
	
	if (mat_Dabs(diff[0]) < symprec && 
	    mat_Dabs(diff[1]) < symprec &&
	    mat_Dabs(diff[2]) < symprec) {
	  is_all_ok = 1;
	  break;
	}
      }
      
      if (! is_all_ok) {
	break;
      }
    }

    if (is_all_ok) {
      ir_rot[num_rot] = i;
      num_rot++;
    }
  }

  rot_reciprocal_q = mat_alloc_MatINT(num_rot);
  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot_reciprocal_q->mat[i],
		       rot_reciprocal->mat[ir_rot[i]]);  
  }

  free(ir_rot);

  return rot_reciprocal_q;
}

static int get_ir_reciprocal_mesh(int grid_address[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  const MatINT *rot_reciprocal)
{
  /* In the following loop, mesh is doubled. */
  /* Even and odd mesh numbers correspond to */
  /* is_shift[i] are 0 or 1, respectively. */
  /* is_shift = [0,0,0] gives Gamma center mesh. */
  /* grid: reducible grid points */
  /* map: the mapping from each point to ir-point. */

  int i, j, k, l, grid_point, grid_point_rot, num_ir = 0;
  int address[3], address_double[3], address_double_rot[3];

  /* "-1" means the element is not touched yet. */
  for (i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++) {
    map[i] = -1;
  }

#ifndef GRID_ORDER_XYZ
  for (i = 0; i < mesh[2]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[0]; k++) {
	address[0] = k;
	address[1] = j;
	address[2] = i;
#else
  for (i = 0; i < mesh[0]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[2]; k++) {
	address[0] = i;
	address[1] = j;
	address[2] = k;
#endif	
	for (l = 0; l < 3; l++) {
	  address_double[l] = address[l] * 2 + is_shift[l];
	}
	grid_point = get_grid_point_double_mesh(address_double, mesh);
	reduce_grid_address(grid_address[grid_point], address, mesh);

	for (l = 0; l < rot_reciprocal->size; l++) {
	  mat_multiply_matrix_vector_i3(address_double_rot,
					rot_reciprocal->mat[l],
					address_double);
	  grid_point_rot = get_grid_point_double_mesh(address_double_rot, mesh);

	  if (grid_point_rot > -1) { /* Invalid if even --> odd or odd --> even */
	    if (map[grid_point_rot] > -1) {
	      map[grid_point] = map[grid_point_rot];
	      break;
	    }
	  }
	}
	
	if (map[grid_point] == -1) {
	  map[grid_point] = grid_point;
	  num_ir++;
	}
      }
    }
  }

  return num_ir;
}

static int
get_ir_reciprocal_mesh_openmp(int grid_address[][3],
			      int map[],
			      const int mesh[3],
			      const int is_shift[3],
			      const MatINT * rot_reciprocal)
{
  int i, j, k, l, grid_point, grid_point_rot, num_ir;
  int address[3], address_double[3], address_double_rot[3];

#ifndef GRID_ORDER_XYZ
#pragma omp parallel for private(j, k, l, grid_point, grid_point_rot, address_double, address_double_rot)
  for (i = 0; i < mesh[2]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[0]; k++) {
	address[0] = k;
	address[1] = j;
	address[2] = i;
#else
#pragma omp parallel for private(j, k, l, grid_point, grid_point_rot, address_double, address_double_rot)
  for (i = 0; i < mesh[0]; i++) {
    for (j = 0; j < mesh[1]; j++) {
      for (k = 0; k < mesh[2]; k++) {
	address[0] = i;
	address[1] = j;
	address[2] = k;
#endif	
	for (l = 0; l < 3; l++) {
	  address_double[l] = address[l] * 2 + is_shift[l];
	}

	grid_point = get_grid_point_double_mesh(address_double, mesh);
	map[grid_point] = grid_point;
	reduce_grid_address(grid_address[grid_point], address, mesh);

	for (l = 0; l < rot_reciprocal->size; l++) {
	  mat_multiply_matrix_vector_i3(address_double_rot,
					rot_reciprocal->mat[l],
					address_double);
	  grid_point_rot = get_grid_point_double_mesh(address_double_rot, mesh);

	  if (grid_point_rot > -1) { /* Invalid if even --> odd or odd --> even */
	    if (grid_point_rot < map[grid_point]) {
	      map[grid_point] = grid_point_rot;
	    }
	  }
	}
      }
    }
  }

  num_ir = 0;

#pragma omp parallel for reduction(+:num_ir)
  for (i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++) {
    if (map[i] == i) {
      num_ir++;
    }
  }
  
  return num_ir;
}

/* Relocate grid addresses to first Brillouin zone */
/* bz_grid_address[prod(mesh + 1)][3] */
/* bz_map[prod(mesh * 2)] */
static int relocate_BZ_grid_address(int bz_grid_address[][3],
				    int bz_map[],
				    SPGCONST int grid_address[][3],
				    const int mesh[3],
				    SPGCONST double rec_lattice[3][3],
				    const int is_shift[3])
{
  double tolerance, min_distance;
  double q_vector[3], distance[KPT_NUM_BZ_SEARCH_SPACE];
  int bzmesh[3], bz_address_double[3];
  int i, j, k, min_index, boundary_num_gp, total_num_gp, bzgp, gp;

  tolerance = get_tolerance_for_BZ_reduction(rec_lattice, mesh);
  for (i = 0; i < 3; i++) {
    bzmesh[i] = mesh[i] * 2;
  }
  for (i = 0; i < bzmesh[0] * bzmesh[1] * bzmesh[2]; i++) {
    bz_map[i] = -1;
  }
  
  boundary_num_gp = 0;
  total_num_gp = mesh[0] * mesh[1] * mesh[2];
  for (i = 0; i < total_num_gp; i++) {
    for (j = 0; j < KPT_NUM_BZ_SEARCH_SPACE; j++) {
      for (k = 0; k < 3; k++) {
	q_vector[k] = 
	  ((grid_address[i][k] + kpt_bz_search_space[j][k] * mesh[k]) * 2 +
	   is_shift[k]) / ((double)mesh[k]) / 2;
      }
      mat_multiply_matrix_vector_d3(q_vector, rec_lattice, q_vector);
      distance[j] = mat_norm_squared_d3(q_vector);
    }
    min_distance = distance[0];
    min_index = 0;
    for (j = 1; j < KPT_NUM_BZ_SEARCH_SPACE; j++) {
      if (distance[j] < min_distance) {
	min_distance = distance[j];
	min_index = j;
      }
    }

    for (j = 0; j < KPT_NUM_BZ_SEARCH_SPACE; j++) {
      if (distance[j] < min_distance + tolerance) {
	if (j == min_index) {
	  gp = i;
	} else {
	  gp = boundary_num_gp + total_num_gp;
	}
	
	for (k = 0; k < 3; k++) {
	  bz_grid_address[gp][k] = 
	    grid_address[i][k] + kpt_bz_search_space[j][k] * mesh[k];
	  bz_address_double[k] = bz_grid_address[gp][k] * 2 + is_shift[k];
	}
	bzgp = get_grid_point_double_mesh(bz_address_double, bzmesh);
	bz_map[bzgp] = gp;
	if (j != min_index) {
	  boundary_num_gp++;
	}
      }
    }
  }

  return boundary_num_gp + total_num_gp;
}

static double get_tolerance_for_BZ_reduction(SPGCONST double rec_lattice[3][3],
					     const int mesh[3])
{
  int i, j;
  double tolerance;
  double length[3];
  
  for (i = 0; i < 3; i++) {
    length[i] = 0;
    for (j = 0; j < 3; j++) {
      length[i] += rec_lattice[j][i] * rec_lattice[j][i];
    }
    length[i] /= mesh[i] * mesh[i];
  }
  tolerance = length[0];
  for (i = 1; i < 3; i++) {
    if (tolerance < length[i]) {
      tolerance = length[i];
    }
  }
  tolerance *= 0.01;
  
  return tolerance;
}
 
static int get_grid_point_double_mesh(const int address_double[3],
				      const int mesh[3])
{
  int i, address[3];

  for (i = 0; i < 3; i++) {
    if (address_double[i] % 2 == 0) {
      address[i] = address_double[i] / 2;
    } else {
      address[i] = (address_double[i] - 1) / 2;
    }
  }
  mat_modulo_i3(address, mesh);
  return get_grid_point_single_mesh(address, mesh);
}

static int get_grid_point_single_mesh(const int address[3],
				      const int mesh[3])
{  
#ifndef GRID_ORDER_XYZ
  return address[2] * mesh[0] * mesh[1] + address[1] * mesh[0] + address[0];
#else
  return address[0] * mesh[1] * mesh[2] + address[1] * mesh[2] + address[2];
#endif  
}

static void reduce_grid_address(int reduced_address[3],
				const int address[3],
				const int mesh[3])
{
  int i;

  for (i = 0; i < 3; i++) {
#ifndef GRID_BOUNDARY_AS_NEGATIVE
    reduced_address[i] = address[i] - mesh[i] * (address[i] > mesh[i] / 2);
#else
    reduced_address[i] = address[i] - mesh[i] * (address[i] >= mesh[i] / 2);
#endif
  }  
}
