/* Copyright (C) 2011 Atsushi Togo */
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
#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"
#include "sitesym_database.h"

#include "debug.h"

static VecDBL * get_exact_positions(int * wyckoffs,
				    int * equiv_atoms,
				    SPGCONST Cell * bravais,
				    SPGCONST Symmetry * conv_sym,
				    const int hall_number,
				    const double symprec);
static void get_exact_location(double position[3],
			       SPGCONST Symmetry * conv_sym,
			       SPGCONST double bravais_lattice[3][3],
			       const double symprec);
static int get_Wyckoff_notation(double position[3],
				SPGCONST Symmetry * conv_sym,
				SPGCONST double bravais_lattice[3][3],
				const int hall_number,
				const double symprec);

/* Return if failed */
VecDBL * ssm_get_exact_positions(int *wyckoffs,
				 int *equiv_atoms,
				 SPGCONST Cell * bravais,
				 SPGCONST Symmetry * conv_sym,
				 const int hall_number,
				 const double symprec)
{
  return get_exact_positions(wyckoffs,
			     equiv_atoms,
			     bravais,
			     conv_sym,
			     hall_number,
			     symprec);
}

/* Return NULL if failed */
static VecDBL * get_exact_positions(int * wyckoffs,
				    int * equiv_atoms,
				    SPGCONST Cell * bravais,
				    SPGCONST Symmetry * conv_sym,
				    const int hall_number,
				    const double symprec)
{
  int i, j, k, l, num_indep_atoms;
  double pos[3];
  int *indep_atoms;
  VecDBL *positions;

  debug_print("get_exact_positions\n");

  indep_atoms = NULL;
  positions = NULL;

  if ((indep_atoms = (int*) malloc(sizeof(int) * bravais->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  if ((positions = mat_alloc_VecDBL(bravais->size)) == NULL) {
    free(indep_atoms);
    indep_atoms = NULL;
    return NULL;
  }

  num_indep_atoms = 0;

  for (i = 0; i < bravais->size; i++) {
    /* Check if atom_i overlap to an atom already set at the exact position. */
    for (j = 0; j < num_indep_atoms; j++) {
      for (k = 0; k < conv_sym->size; k++) {
	mat_multiply_matrix_vector_id3(pos,
				       conv_sym->rot[k],
				       positions->vec[indep_atoms[j]]);
	for (l = 0; l < 3; l++) {
	  pos[l] += conv_sym->trans[k][l];
	}
	if (cel_is_overlap(pos,
			   bravais->position[i],
			   bravais->lattice,
			   symprec)) {
	  /* Equivalent atom was found. */
	  for (l = 0; l < 3; l++) {
	    positions->vec[i][l] = mat_Dmod1(pos[l]);
	  }
	  wyckoffs[i] = wyckoffs[indep_atoms[j]];
	  equiv_atoms[i] = indep_atoms[j];
	  goto escape;
	}
      }
    }
    
    /* No equivalent atom was found. */
    indep_atoms[num_indep_atoms] = i;
    num_indep_atoms++;
    mat_copy_vector_d3(positions->vec[i], bravais->position[i]);
    get_exact_location(positions->vec[i],
		       conv_sym,
		       bravais->lattice,
		       symprec);
    wyckoffs[i] = get_Wyckoff_notation(positions->vec[i],
				       conv_sym,
				       bravais->lattice,
				       hall_number,
				       symprec);
    equiv_atoms[i] = i;
  escape:
    ;
  }

  free(indep_atoms);
  indep_atoms = NULL;

  for (i = 0; i < bravais->size; i++) {
    if (wyckoffs[i] == -1) {
      mat_free_VecDBL(positions);
      positions = NULL;
      break;
    }
  }

  return positions;
}


/* Site-symmetry is used to determine exact location of an atom */
/* R. W. Grosse-Kunstleve and P. D. Adams */
/* Acta Cryst. (2002). A58, 60-65 */
static void get_exact_location(double position[3],
			       SPGCONST Symmetry * conv_sym,
			       SPGCONST double bravais_lattice[3][3],
			       const double symprec)
{
  int i, j, k, num_sum;
  double sum_rot[3][3];
  double pos[3], sum_trans[3];

  debug_print("get_exact_location\n");

  num_sum = 0;
  for (i = 0; i < 3; i++) {
    sum_trans[i] = 0.0;
    for (j = 0; j < 3; j++) {
      sum_rot[i][j] = 0;
    }
  }
  
  for (i = 0; i < conv_sym->size; i++) {
    mat_multiply_matrix_vector_id3(pos,
				   conv_sym->rot[i],
				   position);
    for (j = 0; j < 3; j++) {
      pos[j] += conv_sym->trans[i][j];
    }

    if (cel_is_overlap(pos,
		       position,
		       bravais_lattice,
		       symprec)) {
      for (j = 0; j < 3; j++) {
	sum_trans[j] += conv_sym->trans[i][j] - 
	  mat_Nint(pos[j] - position[j]);
	for (k = 0; k < 3; k++) {
	  sum_rot[j][k] += conv_sym->rot[i][j][k];
	}
      }
      num_sum++;
    }
  }

  for (i = 0; i < 3; i++) {
    sum_trans[i] /= num_sum;
    for (j = 0; j < 3; j++) {
      sum_rot[i][j] /= num_sum;
    }
  }

  /* (sum_rot|sum_trans) is the special-position operator. */
  /* Elements of sum_rot can be fractional values. */
  mat_multiply_matrix_vector_d3(position,
				sum_rot,
				position);

  for (i = 0; i < 3; i++) {
    position[i] += sum_trans[i];
  }

}

/* Return -1 if failed */
static int get_Wyckoff_notation(double position[3],
				SPGCONST Symmetry * conv_sym,
				SPGCONST double bravais_lattice[3][3],
				const int hall_number,
				const double symprec)
{
  int i, j, k, l, at_orbit, num_sitesym, wyckoff_letter;
  int indices_wyc[2];
  int rot[3][3];
  double trans[3], orbit[3];
  VecDBL *pos_rot;

  debug_print("get_Wyckoff_notation\n");

  wyckoff_letter = -1;
  pos_rot = NULL;

  if ((pos_rot = mat_alloc_VecDBL(conv_sym->size)) == NULL) {
    return -1;
  }

  for (i = 0; i < conv_sym->size; i++) {
    mat_multiply_matrix_vector_id3(pos_rot->vec[i], conv_sym->rot[i], position);
    for (j = 0; j < 3; j++) {
      pos_rot->vec[i][j] += conv_sym->trans[i][j];
    }
  }

  ssmdb_get_wyckoff_indices(indices_wyc, hall_number);
  for (i = 0; i < indices_wyc[1]; i++) {
    num_sitesym = ssmdb_get_coordinate(rot, trans, i + indices_wyc[0]);
    for (j = 0; j < pos_rot->size; j++) {
      at_orbit = 0;
      for (k = 0; k < pos_rot->size; k++) {
	if (cel_is_overlap(pos_rot->vec[j],
			   pos_rot->vec[k],
			   bravais_lattice,
			   symprec)) {
	  mat_multiply_matrix_vector_id3(orbit, rot, pos_rot->vec[k]);
	  for (l = 0; l < 3; l++) {
	    orbit[l] += trans[l];
	  }
	  if (cel_is_overlap(pos_rot->vec[k],
			     orbit,
			     bravais_lattice,
			     symprec)) {
	    at_orbit++;
	  }
	}
      }
      if (at_orbit == conv_sym->size / num_sitesym) {
	/* Database is made reversed order, e.g., gfedcba. */
	/* wyckoff is set 0 1 2 3 4... for a b c d e..., respectively. */
	wyckoff_letter = indices_wyc[1] - i - 1;
	goto end;
      }
    }
  }

 end:
  mat_free_VecDBL(pos_rot);
  return wyckoff_letter;
}
