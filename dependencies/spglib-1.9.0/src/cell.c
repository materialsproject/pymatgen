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

#include <stdlib.h>
#include <stdio.h>
#include "cell.h"
#include "mathfunc.h"

#include "debug.h"

#include <assert.h>


#define INCREASE_RATE 2.0
#define REDUCE_RATE 0.95


static Cell * trim_cell(int * mapping_table,
			SPGCONST double trimmed_lattice[3][3],
			SPGCONST Cell * cell,
			const double symprec);
static void set_positions(Cell * trim_cell,
			  const VecDBL * position,
			  const int * mapping_table,
			  const int * overlap_table);
static VecDBL *
translate_atoms_in_trimmed_lattice(SPGCONST Cell * cell,
				   SPGCONST double prim_lat[3][3]);
static int * get_overlap_table(const VecDBL * position,
			       SPGCONST Cell * cell,
			       SPGCONST Cell * trimmed_cell,
			       const double symprec);

/* NULL is returned if faied */
Cell * cel_alloc_cell(const int size)
{
  Cell *cell;
  int i, j;
  
  cell = NULL;

  if ((cell = (Cell*) malloc(sizeof(Cell))) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      cell->lattice[i][j] = 0;
    }
  }
  cell->size = size;
  
  if (size > 0) {
    if ((cell->types = (int *) malloc(sizeof(int) * size)) == NULL) {
      warning_print("spglib: Memory could not be allocated.");
      free(cell);
      cell = NULL;
      return NULL;
    }
    if ((cell->position =
	 (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
      warning_print("spglib: Memory could not be allocated.");
      free(cell->types);
      cell->types = NULL;
      free(cell);
      cell = NULL;
      return NULL;
    }
  }

  return cell;
}

void cel_free_cell(Cell * cell)
{
  if (cell->size > 0) {
    free(cell->position);
    cell->position = NULL;
    free(cell->types);
    cell->types = NULL;
  }
  free (cell);
  cell = NULL;
}

void cel_set_cell(Cell * cell,
		  SPGCONST double lattice[3][3],
		  SPGCONST double position[][3],
		  const int types[])
{
  int i, j;
  mat_copy_matrix_d3(cell->lattice, lattice);
  for (i = 0; i < cell->size; i++) {
    for (j = 0; j < 3; j++) {
      cell->position[i][j] = position[i][j];
    }
    cell->types[i] = types[i];
  }
}

Cell * cel_copy_cell(SPGCONST Cell * cell)
{
  Cell * cell_new;

  cell_new = NULL;
  
  if ((cell_new = cel_alloc_cell(cell->size)) == NULL) {
    return NULL;
  }

  cel_set_cell(cell_new,
	       cell->lattice,
	       cell->position,
	       cell->types);

  return cell_new;
}

int cel_is_overlap(const double a[3],
		   const double b[3],
		   SPGCONST double lattice[3][3],
		   const double symprec)
{
  int i;
  double v_diff[3];

  for (i = 0; i < 3; i++) {
    v_diff[i] = a[i] - b[i];
    v_diff[i] -= mat_Nint(v_diff[i]);
  }

  mat_multiply_matrix_vector_d3(v_diff, lattice, v_diff);
  if ( mat_norm_squared_d3(v_diff) < symprec * symprec) {
    return 1;
  } else {
    return 0;
  }
}

Cell * cel_trim_cell(int * mapping_table,
		     SPGCONST double trimmed_lattice[3][3],
		     SPGCONST Cell * cell,
		     const double symprec)
{
  return trim_cell(mapping_table,
		   trimmed_lattice,
		   cell,
		   symprec);
}


/* Return NULL if failed */
static Cell * trim_cell(int * mapping_table,
			SPGCONST double trimmed_lattice[3][3],
			SPGCONST Cell * cell,
			const double symprec)
{
  int i, index_atom, ratio;
  Cell *trimmed_cell;
  VecDBL * position;
  int *overlap_table;

  position = NULL;
  overlap_table = NULL;
  trimmed_cell = NULL;

  ratio = mat_Nint(mat_get_determinant_d3(cell->lattice) /
		   mat_get_determinant_d3(trimmed_lattice));

  /* Check if cell->size is dividable by ratio */
  if ((cell->size / ratio) * ratio != cell->size) {
    return NULL;
  }

  if ((trimmed_cell = cel_alloc_cell(cell->size / ratio)) == NULL) {
    return NULL;
  }

  if ((position = translate_atoms_in_trimmed_lattice(cell,
						     trimmed_lattice))
      == NULL) {
    cel_free_cell(trimmed_cell);
    goto err;
  }

  mat_copy_matrix_d3(trimmed_cell->lattice, trimmed_lattice);

  if ((overlap_table = get_overlap_table(position,
					 cell,
					 trimmed_cell,
					 symprec)) == NULL) {
    mat_free_VecDBL(position);
    cel_free_cell(trimmed_cell);
    goto err;
  }

  index_atom = 0;
  for (i = 0; i < cell->size; i++) {
    if (overlap_table[i] == i) {
      mapping_table[i] = index_atom;
      trimmed_cell->types[index_atom] = cell->types[i];
      index_atom++;
    } else {
      mapping_table[i] = mapping_table[overlap_table[i]];
    }
  }

  set_positions(trimmed_cell,
		position,
		mapping_table,
		overlap_table);

  mat_free_VecDBL(position);
  free(overlap_table);

  return trimmed_cell;

 err:
  return NULL;
}

static void set_positions(Cell * trimmed_cell,
			  const VecDBL * position,
			  const int * mapping_table,
			  const int * overlap_table)
{
  int i, j, k, l, multi;

  for (i = 0; i < trimmed_cell->size; i++) {
    for (j = 0; j < 3; j++) {
      trimmed_cell->position[i][j] = 0;
    }
  }

  /* Positions of overlapped atoms are averaged. */
  for (i = 0; i < position->size; i++) {
    j = mapping_table[i];
    k = overlap_table[i];
    for (l = 0; l < 3; l++) {
      /* boundary treatment */
      /* One is at right and one is at left or vice versa. */
      if (mat_Dabs(position->vec[k][l] - position->vec[i][l]) > 0.5) {
	if (position->vec[i][l] < position->vec[k][l]) {
	  trimmed_cell->position[j][l] += position->vec[i][l] + 1;
	} else {
	  trimmed_cell->position[j][l] += position->vec[i][l] - 1;
	}
      } else {
	trimmed_cell->position[j][l] += position->vec[i][l];
      }
    }
	
  }

  multi = position->size / trimmed_cell->size;
  for (i = 0; i < trimmed_cell->size; i++) {
    for (j = 0; j < 3; j++) {
      trimmed_cell->position[i][j] /= multi;
      trimmed_cell->position[i][j] = mat_Dmod1(trimmed_cell->position[i][j]);
    }
  }
}

/* Return NULL if failed */
static VecDBL *
translate_atoms_in_trimmed_lattice(SPGCONST Cell * cell,
				   SPGCONST double trimmed_lattice[3][3])
{
  int i, j;
  double tmp_matrix[3][3], axis_inv[3][3];
  VecDBL * position;

  position = NULL;

  if ((position = mat_alloc_VecDBL(cell->size)) == NULL) {
    return NULL;
  }

  mat_inverse_matrix_d3(tmp_matrix, trimmed_lattice, 0);
  mat_multiply_matrix_d3(axis_inv, tmp_matrix, cell->lattice);

  /* Send atoms into the trimmed cell */
  for (i = 0; i < cell->size; i++) {
    mat_multiply_matrix_vector_d3(position->vec[i],
				  axis_inv,
				  cell->position[i]);
    for (j = 0; j < 3; j++) {
      position->vec[i][j] = mat_Dmod1(position->vec[i][j]);
    }
  }

  return position;
}


/* Return NULL if failed */
static int * get_overlap_table(const VecDBL * position,
			       SPGCONST Cell * cell,
			       SPGCONST Cell * trimmed_cell,
			       const double symprec)
{
  int i, j, attempt, num_overlap, ratio, count;
  double trim_tolerance;
  int *overlap_table;

  trim_tolerance = symprec;

  ratio = cell->size / trimmed_cell->size;

  if ((overlap_table = (int*)malloc(sizeof(int) * cell->size)) == NULL) {
    return NULL;
  }
  
  for (attempt = 0; attempt < 100; attempt++) {
    for (i = 0; i < cell->size; i++) {
      overlap_table[i] = -1;
      num_overlap = 0;
      for (j = 0; j < cell->size; j++) {
	if (cell->types[i] == cell->types[j]) {
	  if (cel_is_overlap(position->vec[i],
			     position->vec[j],
			     trimmed_cell->lattice,
			     trim_tolerance)) {
	    num_overlap++;
	    if (overlap_table[i] == -1) {
	      overlap_table[i] = j;
	      assert(j <= i);
	    }
	  }
	}
      }

      if (num_overlap == ratio)	{
	continue;
      }
      if (num_overlap < ratio) {
	trim_tolerance *= INCREASE_RATE;
	warning_print("spglib: Increase tolerance to %f ", trim_tolerance);
	warning_print("(line %d, %s).\n", __LINE__, __FILE__);
	goto cont;
      }
      if (num_overlap > ratio) {
	trim_tolerance *= REDUCE_RATE;
	warning_print("spglib: Reduce tolerance to %f ", trim_tolerance);
	warning_print("(line %d, %s).\n", __LINE__, __FILE__);
	goto cont;
      }
    }

    for (i = 0; i < cell->size; i++) {
      if (overlap_table[i] != i) {
	continue;
      }
      count = 0;
      for (j = 0; j < cell->size; j++) {
	if (i == overlap_table[j]) {
	  count++;
	}
      }
      assert(count == ratio);
    }

    goto found;

  cont:
    ;
  }

  warning_print("spglib: Could not trim cell well ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return NULL;

found:
  return overlap_table;
}
