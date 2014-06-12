/* spin.c */
/* Copyright (C) 2012 Atsushi Togo */

#include <stdio.h>
#include "mathfunc.h"
#include "symmetry.h"
#include "cell.h"

Symmetry  * spn_get_collinear_operation_with_symmetry(SPGCONST Symmetry *sym_nonspin,
						      SPGCONST Cell *cell,
						      const double spins[],
						      const double symprec) {
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
      mat_copy_matrix_i3(rot->mat[ num_sym ], sym_nonspin->rot[i]);
      mat_copy_vector_d3(trans->vec[ num_sym ], sym_nonspin->trans[i]);
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

Symmetry * spn_get_collinear_operation(SPGCONST Cell *cell,
				       const double spins[],
				       const double symprec) {
  Symmetry *sym_nonspin, *symmetry;

  sym_nonspin = sym_get_operation(cell, symprec);
  symmetry = spn_get_collinear_operation_with_symmetry(sym_nonspin,
						       cell,
						       spins,
						       symprec);
  sym_free_symmetry(sym_nonspin);
  
  return symmetry;
}
