/* Copyright (C) 2012 Atsushi Togo */
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

#include <stdio.h>
#include <stdlib.h>
#include "mathfunc.h"
#include "symmetry.h"
#include "cell.h"
#include "debug.h"

static Symmetry * get_collinear_operations(SPGCONST Symmetry *sym_nonspin,
					   SPGCONST Cell *cell,
					   const double spins[],
					   const double symprec);
static int set_equivalent_atoms(int * equiv_atoms,
				SPGCONST Symmetry *symmetry,
				SPGCONST Cell * cell,
				const double symprec);
static int * get_mapping_table(SPGCONST Symmetry *symmetry,
			       SPGCONST Cell * cell,
			       const double symprec);

/* Return NULL if failed */
Symmetry * spn_get_collinear_operations(int equiv_atoms[],
					SPGCONST Symmetry *sym_nonspin,
					SPGCONST Cell *cell,
					const double spins[],
					const double symprec)
{
  Symmetry *symmetry;

  symmetry = NULL;

  if ((symmetry = get_collinear_operations(sym_nonspin,
					   cell,
					   spins,
					   symprec)) == NULL) {
    return NULL;
  }

  if ((set_equivalent_atoms(equiv_atoms,
			    symmetry,
			    cell,
			    symprec)) == 0) {
    sym_free_symmetry(symmetry);
    symmetry = NULL;
  }

  return symmetry;
}

static Symmetry * get_collinear_operations(SPGCONST Symmetry *sym_nonspin,
					   SPGCONST Cell *cell,
					   const double spins[],
					   const double symprec)
{
  Symmetry *symmetry;
  int i, j, k, sign, is_found, num_sym;
  double pos[3];
  MatINT * rot;
  VecDBL * trans;

  rot = mat_alloc_MatINT(sym_nonspin->size);
  trans = mat_alloc_VecDBL(sym_nonspin->size);
  num_sym = 0;
  
  for (i = 0; i < sym_nonspin->size; i++) {
    sign = 0; /* Set sign as undetermined */
    is_found = 1;
    for (j = 0; j < cell->size; j++) {
      mat_multiply_matrix_vector_id3(pos, sym_nonspin->rot[i], cell->position[j]);
      for (k = 0; k < 3; k++) {
	pos[k] += sym_nonspin->trans[i][k];
      }
      for (k = 0; k < cell->size; k++) {
	if (cel_is_overlap(cell->position[k],
			   pos,
			   cell->lattice,
			   symprec)) {
	  if (sign == 0) {
	    if (mat_Dabs(spins[j] - spins[k]) < symprec) {
	      sign = 1;
	      break;
	    }
	    if (mat_Dabs(spins[j] + spins[k]) < symprec) {
	      sign = -1;
	      break;
	    }
	    is_found = 0;
	    break;
	  } else {
	    if (mat_Dabs(spins[j] - spins[k] * sign) < symprec) {
	      break;
	    } else {
	      is_found = 0;
	      break;
	    }
	  }
	}
      }
      if (! is_found) {
	break;
      }
    }
    if (is_found) {
      mat_copy_matrix_i3(rot->mat[num_sym], sym_nonspin->rot[i]);
      mat_copy_vector_d3(trans->vec[num_sym], sym_nonspin->trans[i]);
      num_sym++;
    }
  }

  symmetry = sym_alloc_symmetry(num_sym);
  for (i = 0; i < num_sym; i++) {
    mat_copy_matrix_i3(symmetry->rot[i], rot->mat[ i ]);
    mat_copy_vector_d3(symmetry->trans[i], trans->vec[ i ]);
  }

  mat_free_MatINT(rot);
  mat_free_VecDBL(trans);

  return symmetry;
}

/* Return 0 if failed */
static int set_equivalent_atoms(int * equiv_atoms,
				SPGCONST Symmetry *symmetry,
				SPGCONST Cell * cell,
				const double symprec)
{
  int i, j, k, is_found;
  double pos[3];
  int *mapping_table;

  mapping_table = NULL;

  if ((mapping_table = get_mapping_table(symmetry, cell, symprec)) == NULL) {
    return 0;
  }
  
  for (i = 0; i < cell->size; i++) {
    if (mapping_table[i] != i) {
      continue;
    }
    is_found = 0;
    for (j = 0; j < symmetry->size; j++) {
      mat_multiply_matrix_vector_id3(pos,
				     symmetry->rot[j],
				     cell->position[i]);
      for (k = 0; k < 3; k++) {
	pos[k] += symmetry->trans[j][k];
      }
      for (k = 0; k < cell->size; k++) {
	if (cel_is_overlap(pos, cell->position[k], cell->lattice, symprec)) {
	  if (mapping_table[k] < i) {
	    equiv_atoms[i] = equiv_atoms[mapping_table[k]];
	    is_found = 1;
	    break;
	  }
	}
      }
      if (is_found) {
	break;
      }
    }
    if (!is_found) {
      equiv_atoms[i] = i;
    }
  }

  for (i = 0; i < cell->size; i++) {
    if (mapping_table[i] == i) {
      continue;
    }
    equiv_atoms[i] = equiv_atoms[mapping_table[i]];
  }

  free(mapping_table);
  mapping_table = NULL;

  return 1;
}

/* Return NULL if failed */
static int * get_mapping_table(SPGCONST Symmetry *symmetry,
			       SPGCONST Cell * cell,
			       const double symprec)
{
  int i, j, k, is_found;
  double pos[3];
  int *mapping_table;
  SPGCONST int I[3][3] = {{ 1, 0, 0},
			  { 0, 1, 0},
			  { 0, 0, 1}};

  mapping_table = NULL;

  if ((mapping_table = (int*) malloc(sizeof(int) * cell->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  for (i = 0; i < cell->size; i++) {
    is_found = 0;
    for (j = 0; j < symmetry->size; j++) {
      if (mat_check_identity_matrix_i3(symmetry->rot[j], I)) {
	for (k = 0; k < 3; k++) {
	  pos[k] = cell->position[i][k] + symmetry->trans[j][k];
	}
	for (k = 0; k < cell->size; k++) {
	  if (cel_is_overlap(pos, cell->position[k], cell->lattice, symprec)) {
	    if (k < i) {
	      mapping_table[i] = mapping_table[k];
	      is_found = 1;
	      break;
	    }
	  }
	}
      }
      if (is_found) {
	break;
      }
    }
    if (!is_found) {
      mapping_table[i] = i;
    }
  }

  return mapping_table;
}
