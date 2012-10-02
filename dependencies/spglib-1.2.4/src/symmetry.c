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
static int get_operation(int rot[][3][3],
			 double trans[][3],
			 SPGCONST Cell * cell,
			 const double symprec);
static Symmetry * reduce_operation(SPGCONST Cell * cell,
				   SPGCONST Symmetry * symmetry,
				   const double symprec);
static int is_overlap_all_atoms(const double test_trans[3],
				SPGCONST int rot[3][3],
				SPGCONST Cell * cell,
				const double symprec,
				const int is_identity);
static PointSymmetry transform_pointsymmetry(SPGCONST PointSymmetry * point_sym_prim,
					     SPGCONST double new_lattice[3][3],
					     SPGCONST double original_lattice[3][3]);
static int get_space_group_operation(int rot[][3][3],
				     double trans[][3],
				     SPGCONST PointSymmetry *lattice_sym,
				     SPGCONST Cell *primitive,
				     const double symprec);
static int get_operation_supercell(int rot[][3][3],
				   double trans[][3],
				   const int num_sym, 
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
  int i, j, num_sym;
  MatINT *rot;
  VecDBL *trans;
  Symmetry *symmetry;
  
  rot = mat_alloc_MatINT(cell->size * 48);
  trans = mat_alloc_VecDBL(cell->size * 48);

  num_sym = get_operation(rot->mat, trans->vec, cell, symprec);


  symmetry = sym_alloc_symmetry(num_sym);
  for (i = 0; i < num_sym; i++) {
    mat_copy_matrix_i3(symmetry->rot[i], rot->mat[i]);
    for (j = 0; j < 3; j++) {
      symmetry->trans[i][j] = trans->vec[i][j] - mat_Nint(trans->vec[i][j]);
    }
  }

  mat_free_MatINT(rot);
  mat_free_VecDBL(trans);

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
  int attempt, multi;
  double tolerance;
  VecDBL * pure_trans;

  tolerance = symprec;
  
  for (attempt = 0; attempt < 100; attempt++) {
    pure_trans = get_translation(identity, cell, tolerance, 1);
    multi = pure_trans->size;
    if ((cell->size / multi) * multi == cell->size) {
      goto found;
    } else {
      mat_free_VecDBL(pure_trans);
      tolerance *= REDUCE_RATE;
    }
  }

  pure_trans = mat_alloc_VecDBL(0);
  warning_print("spglib: Finding pure translation failed (line %d, %s).\n", __LINE__, __FILE__);
  warning_print("        cell->size %d, multi %d\n", cell->size, multi);

 found:

#ifdef SPGWARNING
  if (attempt > 0) {
    printf("spglib: Tolerance to search pure_translation is changed to %f\n", tolerance);
  }
#endif
  
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
static int get_operation(int rot[][3][3],
			 double trans[][3],
			 SPGCONST Cell *cell,
			 const double symprec)
{
  int num_sym, attempt, is_found;
  double tolerance;
  PointSymmetry lattice_sym;
  Cell *primitive;
  VecDBL *pure_trans;
#ifdef DEBUG
  int i;
#endif
  debug_print("*** get_symmetry (found symmetry operations) *** \n");

  num_sym = 0;

  /* Lattice symmetry for input cell*/
  lattice_sym = get_lattice_symmetry(cell, symprec);
  if (lattice_sym.size == 0) {
    debug_print("get_lattice_symmetry failed.\n");
    goto end;
  }

  tolerance = symprec;
  for (attempt = 0; attempt < 100; attempt++) {
    is_found = 0;
    pure_trans = sym_get_pure_translation(cell, tolerance);
    if (pure_trans->size == 0) {
      mat_free_VecDBL(pure_trans);
      goto end;
    }

    /* Obtain primitive cell */
    primitive = prm_get_primitive_with_pure_translations(cell,
							 pure_trans,
							 tolerance);
    if (primitive->size > 0) {
      is_found = 1;
      break;
    }
    cel_free_cell(primitive);
    mat_free_VecDBL(pure_trans);
    tolerance *= REDUCE_RATE;
    warning_print("spglib: Reduce tolerance to %f\n", tolerance);
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }

  if (! is_found) {goto end;}

  lattice_sym = transform_pointsymmetry(&lattice_sym,
					primitive->lattice,
					cell->lattice);
  if (lattice_sym.size == 0) {goto deallocate_and_end;}

  /* Symmetry operation search for primitive cell */
  num_sym = get_space_group_operation(rot, trans, &lattice_sym,
				      primitive, symprec);

  /* Recover symmetry operation for the input structure (overwritten) */
  num_sym = get_operation_supercell(rot,
				    trans,
				    num_sym,
				    pure_trans,
				    cell,
				    primitive);

 deallocate_and_end:
  cel_free_cell(primitive);
  mat_free_VecDBL(pure_trans);

 end:
#ifdef DEBUG
  debug_print("num_sym = %d\n", num_sym);
  debug_print("num_lattice_sym = %d\n", lattice_sym.size);
  debug_print("Lattice \n");
  debug_print_matrix_d3(cell->lattice);
  for (i = 0; i < num_sym; i++) {
    debug_print("--- %d ---\n", i + 1);
    debug_print_matrix_i3(rot[i]);
    debug_print("%f %f %f\n", trans[i][0], trans[i][1], trans[i][2]);
  }
#endif
  
  return num_sym;
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

  debug_print("reduce_operation: num_sym = %d\n", num_sym);

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
  VecDBL *tmp_trans, *trans;

  tmp_trans = mat_alloc_VecDBL(cell->size);
  is_found = (int*) malloc(sizeof(int)*cell->size);
  for (i = 0; i < cell->size; i++) {
    is_found[i] = 0;
  }

  /* Look for the atom index with least number of atoms within same type */
  min_atom_index = get_index_with_least_atoms(cell);

  /* Set min_atom_index as the origin to measure the distance between atoms. */
  mat_multiply_matrix_vector_id3(origin, rot, cell->position[min_atom_index]);

  
#pragma omp parallel for private(j)
  for (i = 0; i < cell->size; i++) {	/* test translation */
    if (cell->types[i] != cell->types[min_atom_index]) {
      continue;
    }

    for (j = 0; j < 3; j++) {
      tmp_trans->vec[i][j] = cell->position[i][j] - origin[j];
    }
    if (is_overlap_all_atoms(tmp_trans->vec[i],
			     rot,
			     cell,
			     symprec,
			     is_identity)) {
      is_found[i] = 1;
    }
  }

  for (i = 0; i < cell->size; i++) {
    num_trans += is_found[i];
  }
  trans = mat_alloc_VecDBL(num_trans);
  num_trans = 0;
  for (i = 0; i < cell->size; i++) {
    if (is_found[i]) {
      mat_copy_vector_d3(trans->vec[num_trans], tmp_trans->vec[i]);
      num_trans++;
    }
  }

  mat_free_VecDBL(tmp_trans);
  free(is_found);
  is_found = NULL;
  
  return trans;
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
    if (is_identity) {
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
	
	/* if (cel_is_overlap(cell->position[j], */
	/* 		     pos_rot, */
	/* 		     cell->lattice, */
	/* 		     symprec)) {*/
	/* is_found = 1; */
	/* break; */
	/* } */
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

static int get_space_group_operation(int rot[][3][3],
				     double trans[][3],
				     SPGCONST PointSymmetry *lattice_sym,
				     SPGCONST Cell *cell,
				     const double symprec)
{
  int i, j, k, num_sym;
  VecDBL **tmp_trans;

  num_sym = 0;
  tmp_trans = (VecDBL**) malloc(sizeof(VecDBL*) * lattice_sym->size);

#pragma omp parallel for
  for (i = 0; i < lattice_sym->size; i++) {
    /* get translation corresponding to a rotation */
    tmp_trans[i] = get_translation(lattice_sym->rot[i], cell, symprec, 0);
  }


  for (i = 0; i < lattice_sym->size; i++) {
    for (j = 0; j < tmp_trans[i]->size; j++) {
      for (k = 0; k < 3; k++) {
	trans[num_sym + j][k] = tmp_trans[i]->vec[j][k];
      }
      mat_copy_matrix_i3(rot[num_sym + j], lattice_sym->rot[i]);
    }
    num_sym += tmp_trans[i]->size;
    mat_free_VecDBL(tmp_trans[i]);
  }

  free(tmp_trans);
  tmp_trans = NULL;

  return num_sym;
}

static int get_operation_supercell(int rot[][3][3],
				   double trans[][3],
				   const int num_sym, 
				   const VecDBL * pure_trans,
				   SPGCONST Cell *cell,
				   SPGCONST Cell *primitive)
{
  int i, j, k, multi;
  double inv_prim_lat[3][3], drot[3][3], trans_mat[3][3], trans_mat_inv[3][3];
  MatINT *rot_prim;
  VecDBL *trans_prim;

  rot_prim = mat_alloc_MatINT(num_sym);
  trans_prim = mat_alloc_VecDBL(num_sym);
  multi = pure_trans->size;

  debug_print("get_operation_supercell\n");

  mat_inverse_matrix_d3(inv_prim_lat, primitive->lattice, 0);
  mat_multiply_matrix_d3(trans_mat, inv_prim_lat, cell->lattice);
  mat_inverse_matrix_d3(trans_mat_inv, trans_mat, 0);

  for(i = 0; i < num_sym; i++) {

    /* Translations  */
    mat_multiply_matrix_vector_d3(trans[i], trans_mat_inv, trans[i]);

    /* Rotations */
    mat_cast_matrix_3i_to_3d(drot, rot[i]);
    mat_get_similar_matrix_d3(drot, drot, trans_mat, 0);
    mat_cast_matrix_3d_to_3i(rot[i], drot);
  }

  for(i = 0; i < num_sym; i++) {
    mat_copy_matrix_i3(rot_prim->mat[i], rot[i]);
    for(j = 0; j < 3; j++)
      trans_prim->vec[i][j] = trans[i][j];
  }

  /* Rotations and translations are copied with the set of */
  /* pure translations. */
  for(i = 0; i < num_sym; i++) {
    for(j = 0; j < multi; j++) {
      mat_copy_matrix_i3(rot[ i * multi + j ], rot_prim->mat[i]);
      for (k = 0; k < 3; k++) {
	trans[i * multi + j][k] =
	  mat_Dmod1(trans_prim->vec[i][k] + pure_trans->vec[j][k]);
      }
    }
  }

  mat_free_MatINT(rot_prim);
  mat_free_VecDBL(trans_prim);

  /* return number of symmetry operation of supercell */
  return num_sym * multi;
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

	debug_print_matrix_i3(axes);
	debug_print_matrix_d3(metric);
	debug_print_matrix_d3(metric_orig);
	
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
      debug_print("eliminate by length\n");
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
	  debug_print("eliminate by angle %f %f %20.17f %e\n", cos1, cos2, x, sin_dtheta2);
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

static PointSymmetry transform_pointsymmetry(SPGCONST PointSymmetry * lat_sym_orig,
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

#ifdef DEBUG
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
