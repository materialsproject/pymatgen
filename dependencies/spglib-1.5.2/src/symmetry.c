/* symmetry.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "debug.h"
#include "lattice.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "primitive.h"
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
static Symmetry * get_operations(SPGCONST Cell * cell,
				 const double symprec);
static Symmetry * reduce_operation(SPGCONST Cell * cell,
				   SPGCONST Symmetry * symmetry,
				   const double symprec);
static void search_translation_part(int lat_point_atoms[],
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
static Symmetry * recover_operations_original(SPGCONST Symmetry *symmetry,
					      const VecDBL * pure_trans,
					      SPGCONST Cell *cell,
					      SPGCONST Cell *primitive);
static void set_axes(int axes[3][3],
		     const int a1, const int a2, const int a3);
static PointSymmetry get_lattice_symmetry(SPGCONST Cell *cell,
					  const double symprec);
static int is_identity_metric(SPGCONST double metric_rotated[3][3],
			      SPGCONST double metric_orig[3][3],
			      const double symprec);
static double get_angle(SPGCONST double metric[3][3],
			const int i,
			const int j);

Symmetry * sym_alloc_symmetry(const int size)
{
  Symmetry *symmetry;

  symmetry = (Symmetry*) malloc(sizeof(Symmetry));
  symmetry->size = size;
  if (size > 0) {
    if ((symmetry->rot =
	 (int (*)[3][3]) malloc(sizeof(int[3][3]) * size)) == NULL) {
      warning_print("spglib: Memory could not be allocated ");
      warning_print("(line %d, %s).\n", __LINE__, __FILE__);
      exit(1);
    }
    if ((symmetry->trans =
	 (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
      warning_print("spglib: Memory could not be allocated ");
      warning_print("(line %d, %s).\n", __LINE__, __FILE__);
      exit(1);
    } 
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

Symmetry * sym_get_operation(SPGCONST Cell *cell,
			     const double symprec) {
  Symmetry *symmetry;
  
  symmetry = get_operations(cell, symprec);

  return symmetry;
}

/* Number of operations may be reduced with smaller symprec. */
Symmetry * sym_reduce_operation(SPGCONST Cell * cell,
				SPGCONST Symmetry * symmetry,
				const double symprec)
{
  return reduce_operation(cell, symmetry, symprec);
}

int sym_get_multiplicity(SPGCONST Cell *cell,
			 const double symprec)
{
  int multi;
  VecDBL * trans;

  trans = get_translation(identity, cell, symprec, 1);
  multi = trans->size;
  mat_free_VecDBL(trans);
  return multi;
}

VecDBL * sym_get_pure_translation(SPGCONST Cell *cell,
				  const double symprec)
{
  int multi;
  VecDBL * pure_trans;

  pure_trans = get_translation(identity, cell, symprec, 1);
  multi = pure_trans->size;
  if ((cell->size / multi) * multi == cell->size) {
    debug_print("sym_get_pure_translation: pure_trans->size = %d\n", multi);
  } else {
    ;
    warning_print("spglib: Finding pure translation failed (line %d, %s).\n", __LINE__, __FILE__);
    warning_print("        cell->size %d, multi %d\n", cell->size, multi);
  }

  return pure_trans;
}

VecDBL * sym_reduce_pure_translation(SPGCONST Cell * cell,
				     const VecDBL * pure_trans,
				     const double symprec)
{
  int i, multi;
  Symmetry *symmetry, *symmetry_reduced;
  VecDBL * pure_trans_reduced;
  
  multi = pure_trans->size;
  symmetry = sym_alloc_symmetry(multi);
  for (i = 0; i < multi; i++) {
    mat_copy_matrix_i3(symmetry->rot[i], identity);
    mat_copy_vector_d3(symmetry->trans[i], pure_trans->vec[i]);
  }

  symmetry_reduced = reduce_operation(cell, symmetry, symprec);
  sym_free_symmetry(symmetry);

  multi = symmetry_reduced->size;
  pure_trans_reduced = mat_alloc_VecDBL(multi);
  
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


/* 1) A primitive cell of the input cell is searched. */
/* 2) Pointgroup operations of the primitive cell are obtained. */
/*    These are constrained by the input cell lattice pointgroup, */
/*    i.e., even if the lattice of the primitive cell has higher */
/*    symmetry than that of the input cell, it is not considered. */
/* 3) Spacegroup operations are searched for the primitive cell */
/*    using the constrained point group operations. */
/* 4) The spacegroup operations for the primitive cell are */
/*    transformed to those of original input cells, if the input cell */
/*    was not a primitive cell. */
static Symmetry * get_operations(SPGCONST Cell *cell,
				 const double symprec)
{
  int i, j, attempt;
  double tolerance;
  PointSymmetry lattice_sym;
  Symmetry *symmetry, *symmetry_orig, *symmetry_reduced;
  Primitive primitive;

  debug_print("get_operations:\n");

  symmetry_orig = NULL;

  lattice_sym = get_lattice_symmetry(cell, symprec);
  if (lattice_sym.size == 0) {
    debug_print("get_lattice_symmetry failed.\n");
    goto end;
  }

  primitive = prm_get_primitive_and_pure_translations(cell, symprec);
  if (primitive.cell->size == 0) {
    goto deallocate_and_end;
  }

  lattice_sym = transform_pointsymmetry(&lattice_sym,
					primitive.cell->lattice,
					cell->lattice);
  if (lattice_sym.size == 0) {
    goto deallocate_and_end;
  }

  
  symmetry = get_space_group_operations(&lattice_sym,
					primitive.cell,
					symprec);
  if (symmetry->size > 48) {
    tolerance = symprec;
    for (attempt = 0; attempt < 100; attempt++) {
      tolerance *= REDUCE_RATE;
      warning_print("spglib: number of symmetry operations for primitive cell > 48 was found. (line %d, %s).\n", __LINE__, __FILE__);
      warning_print("tolerance is reduced to %f\n", tolerance);
      symmetry_reduced = reduce_operation(primitive.cell,
					  symmetry,
					  tolerance);
      sym_free_symmetry(symmetry);
      symmetry = symmetry_reduced;
      if (symmetry_reduced->size > 48) {
	;
      } else {
	break;
      }
    }
  }

  symmetry_orig = recover_operations_original(symmetry,
					      primitive.pure_trans,
					      cell,
					      primitive.cell);
  sym_free_symmetry(symmetry);

  for (i = 0; i < symmetry_orig->size; i++) {
    for (j = 0; j < 3; j++) {
      symmetry_orig->trans[i][j] -= mat_Nint(symmetry_orig->trans[i][j]);
    }
  }

 deallocate_and_end:
  cel_free_cell(primitive.cell);
  mat_free_VecDBL(primitive.pure_trans);

 end:
  if (! symmetry_orig) {
    symmetry_orig = sym_alloc_symmetry(0);
  }
  return symmetry_orig;
}

static Symmetry * reduce_operation(SPGCONST Cell * cell,
				   SPGCONST Symmetry * symmetry,
				   const double symprec)
{
  int i, j, num_sym;
  Symmetry * sym_reduced;
  PointSymmetry point_symmetry;
  MatINT *rot;
  VecDBL *trans;

  debug_print("reduce_operation:\n");

  point_symmetry = get_lattice_symmetry(cell, symprec);
  rot = mat_alloc_MatINT(symmetry->size);
  trans = mat_alloc_VecDBL(symmetry->size);

  num_sym = 0;
  for (i = 0; i < point_symmetry.size; i++) {
    for (j = 0; j < symmetry->size; j++) {
      if (mat_check_identity_matrix_i3(point_symmetry.rot[i],
				       symmetry->rot[j])) {
	if (is_overlap_all_atoms(symmetry->trans[j],
				 symmetry->rot[j],
				 cell,
				 symprec,
				 0)) {
	  mat_copy_matrix_i3(rot->mat[num_sym], symmetry->rot[j]);
	  mat_copy_vector_d3(trans->vec[num_sym], symmetry->trans[j]);
	  num_sym++;
	}
      }
    }
  }

  sym_reduced = sym_alloc_symmetry(num_sym);
  for (i = 0; i < num_sym; i++) {
    mat_copy_matrix_i3(sym_reduced->rot[i], rot->mat[i]);
    mat_copy_vector_d3(sym_reduced->trans[i], trans->vec[i]);
  }

  mat_free_MatINT(rot);
  mat_free_VecDBL(trans);

  debug_print("  num_sym %d -> %d\n", symmetry->size, num_sym);

  return sym_reduced;
}

/* Look for the translations which satisfy the input symmetry operation. */
/* This function is heaviest in this code. */
static VecDBL * get_translation(SPGCONST int rot[3][3],
				SPGCONST Cell *cell,
				const double symprec,
				const int is_identity)
{
  int i, j, min_atom_index, num_trans = 0;
  int *is_found;
  double origin[3];
  VecDBL *trans;

#ifdef _OPENMP
  int num_min_type_atoms;
  int *min_type_atoms;
  double vec[3];
#endif

  is_found = (int*) malloc(sizeof(int)*cell->size);
  for (i = 0; i < cell->size; i++) {
    is_found[i] = 0;
  }

  /* Look for the atom index with least number of atoms within same type */
  min_atom_index = get_index_with_least_atoms(cell);

  /* Set min_atom_index as the origin to measure the distance between atoms. */
  mat_multiply_matrix_vector_id3(origin, rot, cell->position[min_atom_index]);

#ifdef _OPENMP
  if (cell->size < NUM_ATOMS_CRITERION_FOR_OPENMP) {
    search_translation_part(is_found,
			    cell,
			    rot,
			    min_atom_index,
			    origin,
			    symprec,
			    is_identity);
  } else {
    /* Collect indices of atoms with the type where the minimum number */
    /* of atoms belong. */
    min_type_atoms = (int*) malloc(sizeof(int)*cell->size);
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
  }
#else
  search_translation_part(is_found,
			  cell,
			  rot,
			  min_atom_index,
			  origin,
			  symprec,
			  is_identity);
#endif

  
  for (i = 0; i < cell->size; i++) {
    num_trans += is_found[i];
  }
  trans = mat_alloc_VecDBL(num_trans);
  num_trans = 0;
  for (i = 0; i < cell->size; i++) {
    if (is_found[i]) {
      for (j = 0; j < 3; j++) {
	trans->vec[num_trans][j] = cell->position[i][j] - origin[j];
      }
      num_trans++;
    }
  }

  free(is_found);
  is_found = NULL;
  
  return trans;
}

static void search_translation_part(int lat_point_atoms[],
				    SPGCONST Cell * cell,
				    SPGCONST int rot[3][3],
				    const int min_atom_index,
				    const double origin[3],
				    const double symprec,
				    const int is_identity)
{
  int i, j;
  double vec[3];

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
    }
  }
}

static int is_overlap_all_atoms(const double trans[3],
				SPGCONST int rot[3][3],
				SPGCONST Cell * cell,
				const double symprec,
				const int is_identity)
{
  int i, j, k, is_found;
  double symprec2;
  double pos_rot[3], d[3];

  symprec2 = symprec*symprec;
  
  for (i = 0; i < cell->size; i++) {
    if (is_identity) { /* Identity matrix is treated as special for speed. */
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
	  d[k] = pos_rot[k] - cell->position[j][k];
	  d[k] -= mat_Nint(d[k]);
	}
	mat_multiply_matrix_vector_d3(d, cell->lattice, d);
	if (d[0]*d[0]+d[1]*d[1]+d[2]*d[2] < symprec2) {
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
  mapping = (int *) malloc(sizeof(int) * cell->size);
  
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

static Symmetry *
get_space_group_operations(SPGCONST PointSymmetry *lattice_sym,
			   SPGCONST Cell *cell,
			   const double symprec)
{
  int i, j, num_sym, total_num_sym;
  VecDBL **trans;
  Symmetry *symmetry;

  debug_print("get_space_group_operations:\n");
  
  trans = (VecDBL**) malloc(sizeof(VecDBL*) * lattice_sym->size);
  total_num_sym = 0;
  for (i = 0; i < lattice_sym->size; i++) {
    trans[i] = get_translation(lattice_sym->rot[i], cell, symprec, 0);
    total_num_sym += trans[i]->size;
  }

  symmetry = sym_alloc_symmetry(total_num_sym);
  num_sym = 0;
  for (i = 0; i < lattice_sym->size; i++) {
    for (j = 0; j < trans[i]->size; j++) {
      mat_copy_vector_d3(symmetry->trans[num_sym + j], trans[i]->vec[j]);
      mat_copy_matrix_i3(symmetry->rot[num_sym + j], lattice_sym->rot[i]);
    }
    num_sym += trans[i]->size;
  }

  for (i = 0; i < lattice_sym->size; i++) {
    mat_free_VecDBL(trans[i]);
  }
  free(trans);
  trans = NULL;

  return symmetry;
}

static Symmetry * recover_operations_original(SPGCONST Symmetry *symmetry,
					      const VecDBL * pure_trans,
					      SPGCONST Cell *cell,
					      SPGCONST Cell *primitive)
{
  int i, j, k, multi;
  double inv_prim_lat[3][3], drot[3][3], trans_mat[3][3], trans_mat_inv[3][3];
  Symmetry *symmetry_orig, *sym_tmp;

  debug_print("recover_operations_original:\n");

  multi = pure_trans->size;
  sym_tmp = sym_alloc_symmetry(symmetry->size);
  symmetry_orig = sym_alloc_symmetry(symmetry->size * multi);

  mat_inverse_matrix_d3(inv_prim_lat, primitive->lattice, 0);
  mat_multiply_matrix_d3(trans_mat, inv_prim_lat, cell->lattice);
  mat_inverse_matrix_d3(trans_mat_inv, trans_mat, 0);

  for(i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(sym_tmp->rot[i], symmetry->rot[i]);
    mat_copy_vector_d3(sym_tmp->trans[i], symmetry->trans[i]);
  }

  for(i = 0; i < symmetry->size; i++) {
    mat_cast_matrix_3i_to_3d(drot, sym_tmp->rot[i]);
    mat_get_similar_matrix_d3(drot, drot, trans_mat, 0);
    mat_cast_matrix_3d_to_3i(sym_tmp->rot[i], drot);

    mat_multiply_matrix_vector_d3(sym_tmp->trans[i],
				  trans_mat_inv,
				  sym_tmp->trans[i]);
  }


  for(i = 0; i < symmetry->size; i++) {
    for(j = 0; j < multi; j++) {
      mat_copy_matrix_i3(symmetry_orig->rot[i * multi + j], sym_tmp->rot[i]);
      for (k = 0; k < 3; k++) {
	symmetry_orig->trans[i * multi + j][k] =
	  sym_tmp->trans[i][k] + pure_trans->vec[j][k];
      }
    }
  }

  sym_free_symmetry(sym_tmp);

  return symmetry_orig;
}

static PointSymmetry get_lattice_symmetry(SPGCONST Cell *cell,
					  const double symprec)
{
  int i, j, k, num_sym;
  int axes[3][3];
  double lattice[3][3], min_lattice[3][3];
  double metric[3][3], metric_orig[3][3];
  PointSymmetry lattice_sym;

  debug_print("get_lattice_symmetry:\n");

  if (! lat_smallest_lattice_vector(min_lattice,
				    cell->lattice,
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
				 cell->lattice,
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
      if (! abs(mat_get_determinant_i3(lat_sym_new.rot[size])) == 1) {
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
