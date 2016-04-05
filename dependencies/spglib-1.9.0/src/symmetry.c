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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "lattice.h"
#include "mathfunc.h"
#include "symmetry.h"

#include "debug.h"

#define NUM_ATOMS_CRITERION_FOR_OPENMP 1000
#define REDUCE_RATE 0.95
#define PI 3.14159265358979323846
/* Tolerance of angle between lattice vectors in degrees */
/* Negative value invokes converter from symprec. */
static double angle_tolerance = -1.0; 

static int relative_axes[][3] = {
  { 1, 0, 0},
  { 0, 1, 0},
  { 0, 0, 1},
  {-1, 0, 0},
  { 0,-1, 0}, /* 5 */
  { 0, 0,-1},
  { 0, 1, 1},
  { 1, 0, 1},
  { 1, 1, 0},
  { 0,-1,-1}, /* 10 */
  {-1, 0,-1},
  {-1,-1, 0},
  { 0, 1,-1},
  {-1, 0, 1},
  { 1,-1, 0}, /* 15 */
  { 0,-1, 1},
  { 1, 0,-1},
  {-1, 1, 0},
  { 1, 1, 1},
  {-1,-1,-1}, /* 20 */
  {-1, 1, 1},
  { 1,-1, 1},
  { 1, 1,-1},
  { 1,-1,-1},
  {-1, 1,-1}, /* 25 */
  {-1,-1, 1},
};

static int identity[3][3] = {{1, 0, 0},
			     {0, 1, 0},
			     {0, 0, 1}};

static int get_index_with_least_atoms(const Cell *cell);
static VecDBL * get_translation(SPGCONST int rot[3][3],
				SPGCONST Cell *cell,
				const double symprec,
				const int is_identity);
static Symmetry * get_operations(SPGCONST Cell *primitive,
				 const double symprec);
static Symmetry * reduce_operation(SPGCONST Cell * primitive,
				   SPGCONST Symmetry * symmetry,
				   const double symprec);
static int search_translation_part(int lat_point_atoms[],
				   SPGCONST Cell * cell,
				   SPGCONST int rot[3][3],
				   const int min_atom_index,
				   const double origin[3],
				   const double symprec,
				   const int is_identity);
static int is_overlap_all_atoms(const double test_trans[3],
				SPGCONST int rot[3][3],
				SPGCONST Cell * cell,
				const double symprec,
				const int is_identity);
static PointSymmetry
transform_pointsymmetry(SPGCONST PointSymmetry * point_sym_prim,
			SPGCONST double new_lattice[3][3],
			SPGCONST double original_lattice[3][3]);
static Symmetry *
get_space_group_operations(SPGCONST PointSymmetry *lattice_sym,
			   SPGCONST Cell *primitive,
			   const double symprec);
static void set_axes(int axes[3][3],
		     const int a1, const int a2, const int a3);
static PointSymmetry get_lattice_symmetry(SPGCONST double cell_lattice[3][3],
					  const double symprec);
static int is_identity_metric(SPGCONST double metric_rotated[3][3],
			      SPGCONST double metric_orig[3][3],
			      const double symprec);
static double get_angle(SPGCONST double metric[3][3],
			const int i,
			const int j);

/* Return NULL if failed */
Symmetry * sym_alloc_symmetry(const int size)
{
  Symmetry *symmetry;

  symmetry = NULL;

  if (size < 1) {
    return NULL;
  }

  if ((symmetry = (Symmetry*) malloc(sizeof(Symmetry))) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  symmetry->size = size;
  symmetry->rot = NULL;
  symmetry->trans = NULL;

  if ((symmetry->rot =
       (int (*)[3][3]) malloc(sizeof(int[3][3]) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    free(symmetry);
    symmetry = NULL;
    return NULL;
  }
  if ((symmetry->trans =
       (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    free(symmetry->rot);
    symmetry->rot = NULL;
    free(symmetry);
    symmetry = NULL;
    return NULL;
  } 

  return symmetry;
}

void sym_free_symmetry(Symmetry *symmetry)
{
  if (symmetry->size > 0) {
    free(symmetry->rot);
    symmetry->rot = NULL;
    free(symmetry->trans);
    symmetry->trans = NULL;
  }
  free(symmetry);
  symmetry = NULL;
}

/* Return NULL if failed */
Symmetry * sym_get_operation(SPGCONST Cell * primitive,
			     const double symprec)
{

  debug_print("sym_get_operations:\n");

  return get_operations(primitive, symprec);
}

/* Return NULL if failed */
Symmetry * sym_reduce_operation(SPGCONST Cell * primitive,
				SPGCONST Symmetry * symmetry,
				const double symprec)
{
  return reduce_operation(primitive, symmetry, symprec);
}

/* Return NULL if failed */
VecDBL * sym_get_pure_translation(SPGCONST Cell *cell,
				  const double symprec)
{
  int multi;
  VecDBL * pure_trans;

  debug_print("sym_get_pure_translation (tolerance = %f):\n", symprec);

  multi = 0;
  pure_trans = NULL;

  if ((pure_trans = get_translation(identity, cell, symprec, 1)) == NULL) {
    return NULL;
  }

  multi = pure_trans->size;
  if ((cell->size / multi) * multi == cell->size) {
    debug_print("  sym_get_pure_translation: pure_trans->size = %d\n", multi);
  } else {
    ;
    warning_print("spglib: Finding pure translation failed (line %d, %s).\n", __LINE__, __FILE__);
    warning_print("        cell->size %d, multi %d\n", cell->size, multi);
  }

  return pure_trans;
}

/* Return NULL if failed */
VecDBL * sym_reduce_pure_translation(SPGCONST Cell * cell,
				     const VecDBL * pure_trans,
				     const double symprec)
{
  int i, multi;
  Symmetry *symmetry, *symmetry_reduced;
  VecDBL * pure_trans_reduced;

  symmetry = NULL;
  symmetry_reduced = NULL;
  pure_trans_reduced = NULL;
  
  multi = pure_trans->size;

  if ((symmetry = sym_alloc_symmetry(multi)) == NULL) {
    return NULL;
  }

  for (i = 0; i < multi; i++) {
    mat_copy_matrix_i3(symmetry->rot[i], identity);
    mat_copy_vector_d3(symmetry->trans[i], pure_trans->vec[i]);
  }

  if ((symmetry_reduced = reduce_operation(cell, symmetry, symprec)) == NULL) {
    sym_free_symmetry(symmetry);
    return NULL;
  }

  sym_free_symmetry(symmetry);
  multi = symmetry_reduced->size;

  if ((pure_trans_reduced = mat_alloc_VecDBL(multi)) == NULL) {
    sym_free_symmetry(symmetry_reduced);
    return NULL;
  }
  
  for (i = 0; i < multi; i++) {
    mat_copy_vector_d3(pure_trans_reduced->vec[i], symmetry_reduced->trans[i]);
  }
  sym_free_symmetry(symmetry_reduced);
  
  return pure_trans_reduced;
}

void sym_set_angle_tolerance(double tolerance)
{
  angle_tolerance = tolerance;
}

double sym_get_angle_tolerance(void)
{
  return angle_tolerance;
}


/* 1) Pointgroup operations of the primitive cell are obtained. */
/*    These are constrained by the input cell lattice pointgroup, */
/*    i.e., even if the lattice of the primitive cell has higher */
/*    symmetry than that of the input cell, it is not considered. */
/* 2) Spacegroup operations are searched for the primitive cell */
/*    using the constrained point group operations. */
/* 3) The spacegroup operations for the primitive cell are */
/*    transformed to those of original input cells, if the input cell */
/*    was not a primitive cell. */
static Symmetry * get_operations(SPGCONST Cell *primitive,
				 const double symprec)
{
  PointSymmetry lattice_sym;
  Symmetry *symmetry;

  debug_print("get_operations:\n");

  symmetry = NULL;

  lattice_sym = get_lattice_symmetry(primitive->lattice, symprec);
  if (lattice_sym.size == 0) {
    debug_print("get_lattice_symmetry failed.\n");
    goto end;
  }

  if ((symmetry = get_space_group_operations(&lattice_sym,
					     primitive,
					     symprec)) == NULL) {
    goto end;
  }

 end:
  return symmetry;
}

/* Return NULL if failed */
static Symmetry * reduce_operation(SPGCONST Cell * primitive,
				   SPGCONST Symmetry * symmetry,
				   const double symprec)
{
  int i, j, num_sym;
  Symmetry * sym_reduced;
  PointSymmetry point_symmetry;
  MatINT *rot;
  VecDBL *trans;

  debug_print("reduce_operation:\n");

  sym_reduced = NULL;
  rot = NULL;
  trans = NULL;

  point_symmetry = get_lattice_symmetry(primitive->lattice, symprec);

  if ((rot = mat_alloc_MatINT(symmetry->size)) == NULL) {
    return NULL;
  }

  if ((trans = mat_alloc_VecDBL(symmetry->size)) == NULL) {
    mat_free_MatINT(rot);
    return NULL;
  }
  
  num_sym = 0;
  for (i = 0; i < point_symmetry.size; i++) {
    for (j = 0; j < symmetry->size; j++) {
      if (mat_check_identity_matrix_i3(point_symmetry.rot[i],
				       symmetry->rot[j])) {
	if (is_overlap_all_atoms(symmetry->trans[j],
				 symmetry->rot[j],
				 primitive,
				 symprec,
				 0)) {
	  mat_copy_matrix_i3(rot->mat[num_sym], symmetry->rot[j]);
	  mat_copy_vector_d3(trans->vec[num_sym], symmetry->trans[j]);
	  num_sym++;
	}
      }
    }
  }

  if ((sym_reduced = sym_alloc_symmetry(num_sym)) != NULL) {
    for (i = 0; i < num_sym; i++) {
      mat_copy_matrix_i3(sym_reduced->rot[i], rot->mat[i]);
      mat_copy_vector_d3(sym_reduced->trans[i], trans->vec[i]);
    }
  }

  mat_free_MatINT(rot);
  mat_free_VecDBL(trans);

  return sym_reduced;
}

/* Look for the translations which satisfy the input symmetry operation. */
/* This function is heaviest in this code. */
/* Return NULL if failed */
static VecDBL * get_translation(SPGCONST int rot[3][3],
				SPGCONST Cell *cell,
				const double symprec,
				const int is_identity)
{
  int i, j, k, min_atom_index, num_trans;
  int *is_found;
  double origin[3];
  VecDBL *trans;

  debug_print("get_translation (tolerance = %f):\n", symprec);

  num_trans = 0;
  is_found = NULL;
  trans = NULL;

#ifdef _OPENMP
  int num_min_type_atoms;
  int *min_type_atoms;
  double vec[3];

  min_type_atoms = NULL;
#endif

  if ((is_found = (int*) malloc(sizeof(int)*cell->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  for (i = 0; i < cell->size; i++) {
    is_found[i] = 0;
  }

  /* Look for the atom index with least number of atoms within same type */
  min_atom_index = get_index_with_least_atoms(cell);
  if (min_atom_index == -1) {
    debug_print("spglib: get_index_with_least_atoms failed.\n");
    goto ret;
  }

  /* Set min_atom_index as the origin to measure the distance between atoms. */
  mat_multiply_matrix_vector_id3(origin, rot, cell->position[min_atom_index]);

#ifdef _OPENMP
  if (cell->size < NUM_ATOMS_CRITERION_FOR_OPENMP) {
    num_trans = search_translation_part(is_found,
					cell,
					rot,
					min_atom_index,
					origin,
					symprec,
					is_identity);
    if (num_trans == 0) {
      goto ret;
    }
  } else {
    /* Collect indices of atoms with the type where the minimum number */
    /* of atoms belong. */
    if ((min_type_atoms = (int*) malloc(sizeof(int)*cell->size)) == NULL) {
      warning_print("spglib: Memory could not be allocated ");
      goto ret;
    }

    num_min_type_atoms = 0;
    for (i = 0; i < cell->size; i++) {
      if (cell->types[i] == cell->types[min_atom_index]) {
	min_type_atoms[num_min_type_atoms] = i;
	num_min_type_atoms++;
      }
    }
#pragma omp parallel for private(j, vec)
    for (i = 0; i < num_min_type_atoms; i++) {
      for (j = 0; j < 3; j++) {
	vec[j] = cell->position[min_type_atoms[i]][j] - origin[j];
      }
      if (is_overlap_all_atoms(vec,
			       rot,
			       cell,
			       symprec,
			       is_identity)) {
	is_found[min_type_atoms[i]] = 1;
      }
    }

    free(min_type_atoms);
    min_type_atoms = NULL;

    for (i = 0; i < cell->size; i++) {
      num_trans += is_found[i];
    }
  }
#else
  num_trans = search_translation_part(is_found,
				      cell,
				      rot,
				      min_atom_index,
				      origin,
				      symprec,
				      is_identity);
  if (num_trans == 0) {
    goto ret;
  }
#endif

  if ((trans = mat_alloc_VecDBL(num_trans)) == NULL) {
    goto ret;
  }

  k = 0;
  for (i = 0; i < cell->size; i++) {
    if (is_found[i]) {
      for (j = 0; j < 3; j++) {
	trans->vec[k][j] = cell->position[i][j] - origin[j];
      }
      k++;
    }
  }

 ret:
  free(is_found);
  is_found = NULL;
  
  return trans;
}

static int search_translation_part(int lat_point_atoms[],
				   SPGCONST Cell * cell,
				   SPGCONST int rot[3][3],
				   const int min_atom_index,
				   const double origin[3],
				   const double symprec,
				   const int is_identity)
{
  int i, j, num_trans;
  double vec[3];

  num_trans = 0;

  for (i = 0; i < cell->size; i++) {
    if (cell->types[i] != cell->types[min_atom_index]) {
      continue;
    }

    for (j = 0; j < 3; j++) {
      vec[j] = cell->position[i][j] - origin[j];
    }
    if (is_overlap_all_atoms(vec,
			     rot,
			     cell,
			     symprec,
			     is_identity)) {
      lat_point_atoms[i] = 1;
      num_trans++;
    }
  }

  return num_trans;
}

static int is_overlap_all_atoms(const double trans[3],
				SPGCONST int rot[3][3],
				SPGCONST Cell * cell,
				const double symprec,
				const int is_identity)
{
  int i, j, k, is_found;
  double symprec2;
  double pos_rot[3], d_frac[3], d[3];

  symprec2 = symprec * symprec;
  
  for (i = 0; i < cell->size; i++) {
    if (is_identity) { /* Identity matrix is treated as special for speed-up. */
      for (j = 0; j < 3; j++) {
	pos_rot[j] = cell->position[i][j] + trans[j];
      }
    } else {
      mat_multiply_matrix_vector_id3(pos_rot,
				     rot,
				     cell->position[i]);
      for (j = 0; j < 3; j++) {
	pos_rot[j] += trans[j];
      }
    }

    is_found = 0;
    for (j = 0; j < cell->size; j++) {
      if (cell->types[i] == cell->types[j]) {
	/* here cel_is_overlap can be used, but for the tuning */
	/* purpose, write it again */
	for (k = 0; k < 3; k++) {
	  d_frac[k] = pos_rot[k] - cell->position[j][k];
	  d_frac[k] -= mat_Nint(d_frac[k]);
	}
	mat_multiply_matrix_vector_d3(d, cell->lattice, d_frac);
	if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < symprec2) {
	  is_found = 1;
	  break;
	}
      }
    }

    if (! is_found) {
      goto not_found;
    }
  }

  return 1;  /* found */

 not_found:
  return 0;
}

static int get_index_with_least_atoms(const Cell *cell)
{
  int i, j, min, min_index;
  int *mapping;

  mapping = NULL;

  if ((mapping = (int *) malloc(sizeof(int) * cell->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return -1;
  }
  
  for (i = 0; i < cell->size; i++) {
    mapping[i] = 0;
  }
  
  for (i = 0; i < cell->size; i++) {
    for (j = 0; j < cell->size; j++) {
      if (cell->types[i] == cell->types[j]) {
	mapping[j]++;
	break;
      }
    }
  }
  
  min = mapping[0];
  min_index = 0;
  for (i = 0; i < cell->size; i++) {
    if (min > mapping[i] && mapping[i] >0) {
      min = mapping[i];
      min_index = i;
    }
  }

  free(mapping);
  mapping = NULL;

  return min_index;
}

/* Return NULL if failed */
static Symmetry *
get_space_group_operations(SPGCONST PointSymmetry *lattice_sym,
			   SPGCONST Cell *primitive,
			   const double symprec)
{
  int i, j, num_sym, total_num_sym;
  VecDBL **trans;
  Symmetry *symmetry;

  debug_print("get_space_group_operations (tolerance = %f):\n", symprec);

  trans = NULL;
  symmetry = NULL;
  
  if ((trans = (VecDBL**) malloc(sizeof(VecDBL*) * lattice_sym->size))
      == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  for (i = 0; i < lattice_sym->size; i++) {
    trans[i] = NULL;
  }

  total_num_sym = 0;
  for (i = 0; i < lattice_sym->size; i++) {

    if ((trans[i] = get_translation(lattice_sym->rot[i], primitive, symprec, 0))
	!= NULL) {

      debug_print("  match translation %d/%d; tolerance = %f\n",
		  i + 1, lattice_sym->size, symprec);

      total_num_sym += trans[i]->size;
    }
  }

  if ((symmetry = sym_alloc_symmetry(total_num_sym)) == NULL) {
    goto ret;
  }

  num_sym = 0;
  for (i = 0; i < lattice_sym->size; i++) {
    if (trans[i] == NULL) {
      continue;
    }
    for (j = 0; j < trans[i]->size; j++) {
      mat_copy_vector_d3(symmetry->trans[num_sym + j], trans[i]->vec[j]);
      mat_copy_matrix_i3(symmetry->rot[num_sym + j], lattice_sym->rot[i]);
    }
    num_sym += trans[i]->size;
  }

 ret:
  for (i = 0; i < lattice_sym->size; i++) {
    if (trans[i] != NULL) {
      mat_free_VecDBL(trans[i]);
    }
  }
  free(trans);
  trans = NULL;

  return symmetry;
}

static PointSymmetry get_lattice_symmetry(SPGCONST double cell_lattice[3][3],
					  const double symprec)
{
  int i, j, k, num_sym;
  int axes[3][3];
  double lattice[3][3], min_lattice[3][3];
  double metric[3][3], metric_orig[3][3];
  PointSymmetry lattice_sym;

  debug_print("get_lattice_symmetry:\n");

  if (! lat_smallest_lattice_vector(min_lattice,
				    cell_lattice,
				    symprec)) {
    goto err;
  }

  mat_get_metric(metric_orig, min_lattice);

  num_sym = 0;
  for (i = 0; i < 26; i++) {
    for (j = 0; j < 26; j++) {
      for (k = 0; k < 26; k++) {
	set_axes(axes, i, j, k);
	if (! ((mat_get_determinant_i3(axes) == 1) ||
	       (mat_get_determinant_i3(axes) == -1))) {
	  continue;
	}
	mat_multiply_matrix_di3(lattice, min_lattice, axes);
	mat_get_metric(metric, lattice);
	
	if (is_identity_metric(metric, metric_orig, symprec)) {
	  mat_copy_matrix_i3(lattice_sym.rot[num_sym], axes);
	  num_sym++;
	}
	  
	if (num_sym > 48) {
	  warning_print("spglib: Too many lattice symmetries was found.\n");
	  warning_print("        Tolerance may be too large ");
	  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
	  goto err;
	}
      }
    }
  }

  lattice_sym.size = num_sym;
  return transform_pointsymmetry(&lattice_sym,
				 cell_lattice,
				 min_lattice);
  
 err:
  lattice_sym.size = 0;
  return lattice_sym;
}

static int is_identity_metric(SPGCONST double metric_rotated[3][3],
			      SPGCONST double metric_orig[3][3],
			      const double symprec)
{
  int i, j, k;
  int elem_sets[3][2] = {{0, 1},
			 {0, 2},
			 {1, 2}};
  double cos1, cos2, x, length_ave2, sin_dtheta2;
  double length_orig[3], length_rot[3];

  for (i = 0; i < 3; i++) {
    length_orig[i] = sqrt(metric_orig[i][i]);
    length_rot[i] = sqrt(metric_rotated[i][i]);
    if (mat_Dabs(length_orig[i] - length_rot[i]) > symprec) {
      goto fail;
    }
  }
  
  for (i = 0; i < 3; i++) {
    j = elem_sets[i][0];
    k = elem_sets[i][1];
    if (angle_tolerance > 0) {
      if (mat_Dabs(get_angle(metric_orig, j, k) -
		   get_angle(metric_rotated, j, k)) > angle_tolerance) {
	goto fail;
      }
    } else {
      /* dtheta = arccos(cos(theta1) - arccos(cos(theta2))) */
      /*        = arccos(c1) - arccos(c2) */
      /*        = arccos(c1c2 + sqrt((1-c1^2)(1-c2^2))) */
      /* sin(dtheta) = sin(arccos(x)) = sqrt(1 - x^2) */
      cos1 = metric_orig[j][k] / length_orig[j] / length_orig[k];
      cos2 = metric_rotated[j][k] / length_rot[j] / length_rot[k];
      x = cos1 * cos2 + sqrt(1 - cos1 * cos1) * sqrt(1 - cos2 * cos2);
      sin_dtheta2 = 1 - x * x;
      length_ave2 = ((length_orig[j] + length_rot[j]) *
		     (length_orig[k] + length_rot[k])) / 4;
      if (sin_dtheta2 > 1e-12) {
	if (sin_dtheta2 * length_ave2 > symprec * symprec) {
	  goto fail;
	}
      }
    }
  }

  return 1;

 fail:
  return 0;
}

static double get_angle(SPGCONST double metric[3][3],
			const int i,
			const int j)
{
  double length_i, length_j;

  length_i = sqrt(metric[i][i]);
  length_j = sqrt(metric[j][j]);

  return acos(metric[i][j] / length_i / length_j) / PI * 180;
}

static PointSymmetry
transform_pointsymmetry(SPGCONST PointSymmetry * lat_sym_orig,
			SPGCONST double new_lattice[3][3],
			SPGCONST double original_lattice[3][3])
{
  int i, size;
  double trans_mat[3][3], inv_mat[3][3], drot[3][3];
  PointSymmetry lat_sym_new;

  mat_inverse_matrix_d3(inv_mat, original_lattice, 0);
  mat_multiply_matrix_d3(trans_mat, inv_mat, new_lattice);

  size = 0;
  for (i = 0; i < lat_sym_orig->size; i++) {
    mat_cast_matrix_3i_to_3d(drot, lat_sym_orig->rot[i]);
    mat_get_similar_matrix_d3(drot, drot, trans_mat, 0);

    /* new_lattice may have lower point symmetry than original_lattice.*/
    /* The operations that have non-integer elements are not counted. */
    if (mat_is_int_matrix(drot, mat_Dabs(mat_get_determinant_d3(trans_mat)) / 10)) {
      mat_cast_matrix_3d_to_3i(lat_sym_new.rot[size], drot);
      if (abs(mat_get_determinant_i3(lat_sym_new.rot[size])) != 1) {
	warning_print("spglib: A point symmetry operation is not unimodular.");
	warning_print("(line %d, %s).\n", __LINE__, __FILE__);
	goto err;
      }
      size++;
    }
  }

#ifdef SPGWARNING
  if (! (lat_sym_orig->size == size)) {
    warning_print("spglib: Some of point symmetry operations were dropped.");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }
#endif

  lat_sym_new.size = size;
  return lat_sym_new;

 err:
  lat_sym_new.size = 0;
  return lat_sym_new;
}

static void set_axes(int axes[3][3],
		     const int a1, const int a2, const int a3)
{
  int i;
  for (i = 0; i < 3; i++) {axes[i][0] = relative_axes[a1][i]; }
  for (i = 0; i < 3; i++) {axes[i][1] = relative_axes[a2][i]; }
  for (i = 0; i < 3; i++) {axes[i][2] = relative_axes[a3][i]; }
}
