/* primitive.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "lattice.h"
#include "mathfunc.h"
#include "primitive.h"
#include "symmetry.h"

#include "debug.h"

#define INCREASE_RATE 2.0
#define REDUCE_RATE 0.95

static Primitive * get_primitive(SPGCONST Cell * cell, const double symprec);
static int set_primitive_positions(Cell * primitive_cell,
				   const VecDBL * position,
				   const Cell * cell,
				   int * const * table);
static VecDBL * get_positions_primitive(SPGCONST Cell * cell,
					SPGCONST double prim_lat[3][3]);
static int get_overlap_table(int ** table,
			     SPGCONST Cell *primitive_cell,
			     const VecDBL * position,
			     const int *types,
			     const double symprec);
static int check_overlap_table(SPGCONST int **overlap_table,
			       const int cell_size,
			       const int ratio);
static void free_overlap_table(int ** table, const int size);
static int ** allocate_overlap_table(const int size);
static Cell * get_cell_with_smallest_lattice(SPGCONST Cell * cell,
					     const double symprec);
static Cell * get_primitive_cell(int * mapping_table,
				 SPGCONST Cell * cell,
				 const VecDBL * pure_trans,
				 const double symprec);
static int trim_cell(Cell * primitive_cell,
		     int * mapping_table,
		     SPGCONST Cell * cell,
		     const double symprec);
static int get_primitive_lattice_vectors_iterative(double prim_lattice[3][3],
						   SPGCONST Cell * cell,
						   const VecDBL * pure_trans,
						   const double symprec);
static int get_primitive_lattice_vectors(double prim_lattice[3][3],
					 const VecDBL * vectors,
					 SPGCONST Cell * cell,
					 const double symprec);
static VecDBL * get_translation_candidates(const VecDBL * pure_trans);

Primitive * prm_alloc_primitive(const int size)
{
  Primitive *primitive;

  primitive = (Primitive*) malloc(sizeof(Primitive));
  primitive->size = size;
  if (size > 0) {
    if ((primitive->mapping_table = (int*) malloc(sizeof(int) * size)) == NULL)
      {
	warning_print("spglib: Memory could not be allocated ");
	warning_print("(Primitive, line %d, %s).\n", __LINE__, __FILE__);
	exit(1);
      }
  }
  primitive->tolerance = 0;
  primitive->pure_trans = NULL;
  primitive->cell = NULL;

  return primitive;
}

void prm_free_primitive(Primitive * primitive)
{
  if (primitive->size > 0) {
    free(primitive->mapping_table);
    primitive->mapping_table = NULL;
  }
  primitive->size = 0;
  mat_free_VecDBL(primitive->pure_trans);
  cel_free_cell(primitive->cell);
  free(primitive);
  primitive = NULL;
}

Cell * prm_get_primitive_cell(SPGCONST Cell * cell, const double symprec)
{
  Cell *primitive_cell;
  Primitive *primitive;

  primitive = get_primitive(cell, symprec);
  primitive_cell = cel_copy_cell(primitive->cell);
  prm_free_primitive(primitive);

  return primitive_cell;
}

Primitive * prm_get_primitive(SPGCONST Cell * cell, const double symprec)
{
  return get_primitive(cell, symprec);
}

/* If primitive could not be found, primitive->size = 0 is returned. */
static Primitive * get_primitive(SPGCONST Cell * cell, const double symprec)
{
  int i, attempt, is_found = 0;
  double tolerance;
  Primitive *primitive;

  primitive = prm_alloc_primitive(cell->size);

  tolerance = symprec;
  for (attempt = 0; attempt < 100; attempt++) {
    primitive->pure_trans = sym_get_pure_translation(cell, tolerance);
    if (primitive->pure_trans->size == 0) {
      mat_free_VecDBL(primitive->pure_trans);
      continue;
    }

    if (primitive->pure_trans->size == 1) {
      primitive->cell = get_cell_with_smallest_lattice(cell, tolerance);
      for (i = 0; i < cell->size; i++) {
	primitive->mapping_table[i] = i;
      }
    } else {
      primitive->cell = get_primitive_cell(primitive->mapping_table,
					   cell,
					   primitive->pure_trans,
					   tolerance);
    }

    if (primitive->cell->size > 0) {
      is_found = 1;
      break;
    }

    cel_free_cell(primitive->cell);
    mat_free_VecDBL(primitive->pure_trans);
    
    tolerance *= REDUCE_RATE;
    warning_print("spglib: Reduce tolerance to %f ", tolerance);
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }

  if (is_found) {
    primitive->tolerance = tolerance;
  } else {
    primitive->cell = cel_alloc_cell(0);
    primitive->pure_trans = mat_alloc_VecDBL(0);
  }

  return primitive;
}

static Cell * get_cell_with_smallest_lattice(SPGCONST Cell * cell,
					     const double symprec)
{
  int i, j;
  double min_lat[3][3], trans_mat[3][3], inv_lat[3][3];
  Cell * smallest_cell;

  debug_print("get_cell_with_smallest_lattice:\n");
  
  if (lat_smallest_lattice_vector(min_lat,
				  cell->lattice,
				  symprec)) {
    mat_inverse_matrix_d3(inv_lat, min_lat, 0);
    mat_multiply_matrix_d3(trans_mat, inv_lat, cell->lattice);
    smallest_cell = cel_alloc_cell(cell->size);
    mat_copy_matrix_d3(smallest_cell->lattice, min_lat);
    for (i = 0; i < cell->size; i++) {
      smallest_cell->types[i] = cell->types[i];
      mat_multiply_matrix_vector_d3(smallest_cell->position[i],
				    trans_mat, cell->position[i]);
      for (j = 0; j < 3; j++) {
	cell->position[i][j] -= mat_Nint(cell->position[i][j]);
      }
    }
  } else {
    smallest_cell = cel_alloc_cell(0);
  }

  return smallest_cell;
}

/* If primitive could not be found, primitive->size = 0 is returned. */
static Cell * get_primitive_cell(int * mapping_table,
				 SPGCONST Cell * cell,
				 const VecDBL * pure_trans,
				 const double symprec)
{
  int multi;
  double prim_lattice[3][3];
  Cell * primitive_cell;

  debug_print("get_primitive:\n");

  /* Primitive lattice vectors are searched. */
  /* To be consistent, sometimes tolerance is decreased iteratively. */
  /* The descreased tolerance is stored in 'static double tolerance'. */
  multi = get_primitive_lattice_vectors_iterative(prim_lattice,
						  cell,
						  pure_trans,
						  symprec);
  if (! multi) {
    goto not_found;
  }

  primitive_cell = cel_alloc_cell(cell->size / multi);

  if (! lat_smallest_lattice_vector(primitive_cell->lattice,
				    prim_lattice,
				    symprec)) {
    cel_free_cell(primitive_cell);
    goto not_found;
  }

  /* Fit atoms into new primitive cell */
  if (! trim_cell(primitive_cell, mapping_table, cell, symprec)) {
    cel_free_cell(primitive_cell);
    goto not_found;
  }

  /* found */
  return primitive_cell;

 not_found:
  primitive_cell = cel_alloc_cell(0);
  warning_print("spglib: Primitive cell could not be found ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return primitive_cell;
}


static int trim_cell(Cell * primitive_cell,
		     int * mapping_table,
		     SPGCONST Cell * cell,
		     const double symprec)
{
  int i, index_prim_atom;
  VecDBL * position;
  int **overlap_table;

  overlap_table = allocate_overlap_table(cell->size);
  position = get_positions_primitive(cell, primitive_cell->lattice);

  if (! get_overlap_table(overlap_table,
			  primitive_cell,
			  position,
			  cell->types,
			  symprec)) {goto err;}

  index_prim_atom = 0;
  for (i = 0; i < cell->size; i++) {
    if (overlap_table[i][0] == i) {
      mapping_table[i] = index_prim_atom;
      index_prim_atom++;
    } else {
      mapping_table[i] = mapping_table[overlap_table[i][0]];
    }
  }

  if (! set_primitive_positions(primitive_cell,
				position,
				cell,
				overlap_table)) {goto err;}

  mat_free_VecDBL(position);
  free_overlap_table(overlap_table, cell->size);
  return 1;

 err:
  mat_free_VecDBL(position);
  free_overlap_table(overlap_table, cell->size);
  return 0;
}

static int set_primitive_positions(Cell * primitive_cell,
				   const VecDBL * position,
				   const Cell * cell,
				   int * const * table)
{
  int i, j, k, ratio, index_prim_atom;
  int *is_equivalent;

  is_equivalent = (int*)malloc(cell->size * sizeof(int));
  for (i = 0; i < cell->size; i++) {
    is_equivalent[i] = 0;
  }
  ratio = cell->size / primitive_cell->size;

  /* Copy positions. Positions of overlapped atoms are averaged. */
  index_prim_atom = 0;
  for (i = 0; i < cell->size; i++) {

    if (! is_equivalent[i]) {
      primitive_cell->types[index_prim_atom] = cell->types[i];

      for (j = 0; j < 3; j++) {
	primitive_cell->position[index_prim_atom][j] = 0;
      }

      for (j = 0; j < ratio; j++) { /* Loop for averaging positions */
	is_equivalent[table[i][j]] = 1;

	for (k = 0; k < 3; k++) {
	  /* boundary treatment */
	  /* One is at right and one is at left or vice versa. */
	  if (mat_Dabs(position->vec[table[i][0]][k] -
		       position->vec[table[i][j]][k]) > 0.5) {
	    if (position->vec[table[i][j]][k] < 0) {
	      primitive_cell->position[index_prim_atom][k] =
		primitive_cell->position[index_prim_atom][k] +
		position->vec[table[i][j]][k] + 1;
	    } else {
	      primitive_cell->position[index_prim_atom][k] =
		primitive_cell->position[index_prim_atom][k] +
		position->vec[table[i][j]][k] - 1;
	    }

	  } else {
	    primitive_cell->position[index_prim_atom][k] =
	      primitive_cell->position[index_prim_atom][k] +
	      position->vec[table[i][j]][k];
	  }
	}
	
      }

      for (j = 0; j < 3; j++) {	/* take average and reduce */
	primitive_cell->position[index_prim_atom][j] =
	  primitive_cell->position[index_prim_atom][j] / ratio;
	primitive_cell->position[index_prim_atom][j] =
	  primitive_cell->position[index_prim_atom][j] -
	  mat_Nint(primitive_cell->position[index_prim_atom][j]);
      }
      index_prim_atom++;
    }
  }

  free(is_equivalent);
  is_equivalent = NULL;

  if (! (index_prim_atom == primitive_cell->size)) {
    warning_print("spglib: Atomic positions of primitive cell could not be determined ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    goto err;
  }

  return 1;

 err:
  return 0;
}

static VecDBL * get_positions_primitive(SPGCONST Cell * cell,
					SPGCONST double prim_lat[3][3])
{
  int i, j;
  double tmp_matrix[3][3], axis_inv[3][3];
  VecDBL * position;

  position = mat_alloc_VecDBL(cell->size);

  mat_inverse_matrix_d3(tmp_matrix, prim_lat, 0);
  mat_multiply_matrix_d3(axis_inv, tmp_matrix, cell->lattice);

  /* Send atoms into the primitive cell */
  for (i = 0; i < cell->size; i++) {
    mat_multiply_matrix_vector_d3(position->vec[i],
				  axis_inv, cell->position[i]);
    for (j = 0; j < 3; j++) {
      position->vec[i][j] -= mat_Nint(position->vec[i][j]);
    }
  }

  return position;
}


/* If overlap_table is correctly obtained, */
/* shape of overlap_table will be (cell->size, cell->size / primitive->size). */
static int get_overlap_table(int **overlap_table,
			     SPGCONST Cell *primitive_cell,
			     const VecDBL * position,
			     const int *types,
			     const double symprec)
{
  int i, j, attempt, num_overlap, ratio, cell_size;
  double trim_tolerance;

  cell_size = position->size;
  ratio = cell_size / primitive_cell->size;
  trim_tolerance = symprec;
  
  for (attempt = 0; attempt < 100; attempt++) {
    /* Each value of -1 has to be overwritten by 0 or positive numbers. */
    for (i = 0; i < cell_size; i++) {
      for (j = 0; j < cell_size; j++) {
        overlap_table[i][j] = -1;
      }
    }

    for (i = 0; i < cell_size; i++) {
      num_overlap = 0;
      for (j = 0; j < cell_size; j++) {
	if (types[i] == types[j]) {
	  if (cel_is_overlap(position->vec[i],
			     position->vec[j],
			     primitive_cell->lattice,
			     trim_tolerance)) {
	    overlap_table[i][num_overlap] = j;
	    num_overlap++;
	  }
	}
      }
    }

    if (check_overlap_table(overlap_table, cell_size, ratio)) {
      goto found;
    }

    if (num_overlap < ratio) {
      trim_tolerance *= INCREASE_RATE;
      warning_print("spglib: Increase tolerance to %f ", trim_tolerance);
    } else {
      trim_tolerance *= REDUCE_RATE;
      warning_print("spglib: Reduce tolerance to %f ", trim_tolerance);
    }
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }

  warning_print("spglib: Could not trim cell into primitive ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return 0;

 found:
  return 1;

}


static int check_overlap_table(SPGCONST int **overlap_table,
			       const int cell_size,
			       const int ratio) {
  int i, j, index_compared, all_ok;

  all_ok = 1;
  for (i = 0; i < cell_size; i++) {
    index_compared = overlap_table[i][0];
    for (j = 0; j < cell_size; j++) {
      if (! (overlap_table[i][j] == overlap_table[index_compared][j])) {
	all_ok = 0;
	break;
      }
      if (j < ratio) {
	if (overlap_table[i][j] == -1) {
	  all_ok = 0;
	  break;
	}
      } else {
	if (overlap_table[i][j] > -1) {
	  all_ok = 0;
	  break;
	}
      }
    }
    if (! all_ok) {
      break;
    }
  }
  return all_ok;
}

static void free_overlap_table(int **table, const int size)
{
  int i;
  for (i = 0; i < size; i++) {
    free(table[i]);
    table[i] = NULL;
  }
  free(table);
  table = NULL;
}

static int ** allocate_overlap_table(const int size)
{
  int i;
  int **table = (int**)malloc(size * sizeof(int*));
  for (i = 0; i < size; i++) {
    table[i] = (int*)malloc(size * sizeof(int));
  }
  return table;
}


static int get_primitive_lattice_vectors_iterative(double prim_lattice[3][3],
						   SPGCONST Cell * cell,
						   const VecDBL * pure_trans,
						   const double symprec)
{
  int i, multi, attempt;
  double tolerance;
  VecDBL * vectors, * pure_trans_reduced, *tmp_vec;

  tolerance = symprec;
  pure_trans_reduced = mat_alloc_VecDBL(pure_trans->size);
  for (i = 0; i < pure_trans->size; i++) {
    mat_copy_vector_d3(pure_trans_reduced->vec[i], pure_trans->vec[i]);
  }
  
  for (attempt = 0; attempt < 100; attempt++) {
    multi = pure_trans_reduced->size;
    vectors = get_translation_candidates(pure_trans_reduced);

    /* Lattice of primitive cell is found among pure translation vectors */
    if (get_primitive_lattice_vectors(prim_lattice,
				      vectors,
				      cell,
				      tolerance)) {

      mat_free_VecDBL(vectors);
      mat_free_VecDBL(pure_trans_reduced);

      goto found;
    } else {

      tmp_vec = mat_alloc_VecDBL(multi);
      for (i = 0; i < multi; i++) {
	mat_copy_vector_d3(tmp_vec->vec[i], pure_trans_reduced->vec[i]);
      }
      mat_free_VecDBL(pure_trans_reduced);
      pure_trans_reduced = sym_reduce_pure_translation(cell,
						       tmp_vec,
						       tolerance);
      warning_print("Tolerance is reduced to %f (%d), size = %d\n",
		    tolerance, attempt, pure_trans_reduced->size);

      mat_free_VecDBL(tmp_vec);
      mat_free_VecDBL(vectors);

      tolerance *= REDUCE_RATE;
    }
  }

  /* Not found */
  return 0;

 found:
#ifdef SPGWARNING
  if (attempt > 0) {
    printf("spglib: Tolerance to find primitive lattice vectors was changed to %f\n", tolerance);
  }
#endif
  return multi;
}

static int get_primitive_lattice_vectors(double prim_lattice[3][3],
					 const VecDBL * vectors,
					 SPGCONST Cell * cell,
					 const double symprec)
{
  int i, j, k, size;
  double initial_volume, volume;
  double relative_lattice[3][3], min_vectors[3][3], tmp_lattice[3][3];
  double inv_mat_dbl[3][3];
  int inv_mat_int[3][3];

  debug_print("get_primitive_lattice_vectors:\n");

  size = vectors->size;
  initial_volume = mat_Dabs(mat_get_determinant_d3(cell->lattice));

  /* check volumes of all possible lattices, find smallest volume */
  for (i = 0; i < size; i++) {
    for (j = i + 1; j < size; j++) {
      for (k = j + 1; k < size; k++) {
	mat_multiply_matrix_vector_d3(tmp_lattice[0],
				      cell->lattice,
				      vectors->vec[i]);
	mat_multiply_matrix_vector_d3(tmp_lattice[1],
				      cell->lattice,
				      vectors->vec[j]);
	mat_multiply_matrix_vector_d3(tmp_lattice[2],
				      cell->lattice,
				      vectors->vec[k]);
	volume = mat_Dabs(mat_get_determinant_d3(tmp_lattice));
	if (volume > symprec) {
	  if (mat_Nint(initial_volume / volume) == size-2) {
	    mat_copy_vector_d3(min_vectors[0], vectors->vec[i]);
	    mat_copy_vector_d3(min_vectors[1], vectors->vec[j]);
	    mat_copy_vector_d3(min_vectors[2], vectors->vec[k]);
	    goto ret;
	  }
	}
      }
    }
  }

  /* Not found */
  warning_print("spglib: Primitive lattice vectors cound not be found ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return 0;

  /* Found */
 ret:
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      relative_lattice[j][i] = min_vectors[i][j];
    }
  }

  mat_inverse_matrix_d3(inv_mat_dbl, relative_lattice, 0);
  mat_cast_matrix_3d_to_3i(inv_mat_int, inv_mat_dbl);
  if (abs(mat_get_determinant_i3(inv_mat_int)) == size-2) {
    mat_cast_matrix_3i_to_3d(inv_mat_dbl, inv_mat_int);
    mat_inverse_matrix_d3(relative_lattice, inv_mat_dbl, 0);
  } else {
    warning_print("spglib: Primitive lattice cleaning is incomplete ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }
  mat_multiply_matrix_d3(prim_lattice, cell->lattice, relative_lattice);

  return 1;  
}

static VecDBL * get_translation_candidates(const VecDBL * pure_trans)
{
  int i, j, multi;
  VecDBL * vectors;

  multi = pure_trans->size;
  vectors = mat_alloc_VecDBL(multi+2);

  /* store pure translations in original cell */ 
  /* as trial primitive lattice vectors */
  for (i = 0; i < multi - 1; i++) {
    mat_copy_vector_d3(vectors->vec[i], pure_trans->vec[i + 1]);
  }

  /* store lattice translations of original cell */
  /* as trial primitive lattice vectors */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      if (i == j) {
	vectors->vec[i+multi-1][j] = 1;
      } else {
	vectors->vec[i+multi-1][j] = 0;
      }
    }
  }

  return vectors;
}

