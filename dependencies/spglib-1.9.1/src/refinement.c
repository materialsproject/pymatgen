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
/* refinement.c */
/* Copyright (C) 2011 Atsushi Togo */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "refinement.h"
#include "cell.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "spg_database.h"
#include "site_symmetry.h"
#include "symmetry.h"

#include "debug.h"

#define REDUCE_RATE 0.95

static Cell * get_Wyckoff_positions(int * wyckoffs,
				    int * equiv_atoms,
				    SPGCONST Cell * primitive,
				    SPGCONST Cell * cell,
				    SPGCONST Spacegroup * spacegroup,
				    SPGCONST Symmetry * symmetry,
				    const int * mapping_table,
				    const double symprec);
static Cell *
get_bravais_exact_positions_and_lattice(int * wyckoffs,
					int * equiv_atoms,
					SPGCONST Spacegroup * spacegroup,
					SPGCONST Cell * primitive,
					const double symprec);
static Cell * expand_positions(int * wyckoffs,
			       int * equiv_atoms,
			       SPGCONST Cell * conv_prim,
			       SPGCONST Symmetry * conv_sym,
			       const int * wyckoffs_prim,
			       const int * equiv_atoms_prim);
static Cell * get_conventional_primitive(SPGCONST Spacegroup * spacegroup,
					 SPGCONST Cell * primitive);
static int get_number_of_pure_translation(SPGCONST Symmetry * conv_sym);
static void get_conventional_lattice(double lattice[3][3],
				     SPGCONST Spacegroup *spacegroup);
static void set_tricli(double lattice[3][3],
		       SPGCONST double metric[3][3]);
static void set_monocli(double lattice[3][3],
			SPGCONST double metric[3][3]);
static void set_ortho(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_tetra(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_trigo(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_rhomb(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_cubic(double lattice[3][3],
		      SPGCONST double metric[3][3]);

static Symmetry *
get_refined_symmetry_operations(SPGCONST Cell * cell,
				SPGCONST Cell * primitive,
				SPGCONST Spacegroup * spacegroup,
				const double symprec);
static void set_translation_with_origin_shift(Symmetry *conv_sym,
					      const double origin_shift[3]);
static Symmetry * get_primitive_db_symmetry(SPGCONST double t_mat[3][3],
					    const Symmetry *conv_sym);
static void get_corners(int corners[3][8],
			SPGCONST int t_mat[3][3]);
static void get_surrounding_frame(int frame[3],
				  SPGCONST int t_mat[3][3]);
static int set_equivalent_atoms(int * equiv_atoms_cell,
				SPGCONST Cell * primitive,
				SPGCONST Cell * cell,
				const int * equiv_atoms_prim,
				const int * mapping_table);
static void set_equivalent_atoms_broken_symmetry(int * equiv_atoms_cell,
						 SPGCONST Cell * cell,
						 const Symmetry *symmetry,
						 const int * mapping_table,
						 const double symprec);
static int search_equivalent_atom(const int atom_index,
				  SPGCONST Cell * cell,
				  const Symmetry *symmetry,
				  const double symprec);
static Symmetry * reduce_symmetry_in_frame(const int frame[3],
					   SPGCONST Symmetry *prim_sym,
					   SPGCONST int t_mat[3][3],
					   SPGCONST double lattice[3][3],
					   const int multiplicity,
					   const double symprec);
static VecDBL * get_lattice_translations(const int frame[3],
					 SPGCONST double inv_tmat[3][3]);
static VecDBL *
remove_overlapping_lattice_points(SPGCONST double lattice[3][3],
				  const VecDBL *lattice_trans,
				  const double symprec);
static Symmetry *
get_symmetry_in_original_cell(SPGCONST int t_mat[3][3],
			      SPGCONST double inv_tmat[3][3],
			      SPGCONST double lattice[3][3],
			      SPGCONST Symmetry *prim_sym,
			      const double symprec);
static Symmetry *
copy_symmetry_upon_lattice_points(const VecDBL *pure_trans,
				  SPGCONST Symmetry *t_sym);


static SPGCONST int identity[3][3] = {
  { 1, 0, 0},
  { 0, 1, 0},
  { 0, 0, 1},
};


/* Return NULL if failed */
Symmetry *
ref_get_refined_symmetry_operations(SPGCONST Cell * cell,
				    SPGCONST Cell * primitive,
				    SPGCONST Spacegroup * spacegroup,
				    const double symprec)
{
  return get_refined_symmetry_operations(cell,
					 primitive,
					 spacegroup,
					 symprec);
}

/* Return NULL if failed */
Cell * ref_get_Wyckoff_positions(int * wyckoffs,
				 int * equiv_atoms,
				 SPGCONST Cell * primitive,
				 SPGCONST Cell * cell,
				 SPGCONST Spacegroup * spacegroup,
				 SPGCONST Symmetry * symmetry,
				 const int * mapping_table,
				 const double symprec)
{
  return get_Wyckoff_positions(wyckoffs,
			       equiv_atoms,
			       primitive,
			       cell,
			       spacegroup,
			       symmetry,
			       mapping_table,
			       symprec);
}

Cell * get_Wyckoff_positions(int * wyckoffs,
			     int * equiv_atoms,
			     SPGCONST Cell * primitive,
			     SPGCONST Cell * cell,
			     SPGCONST Spacegroup * spacegroup,
			     SPGCONST Symmetry * symmetry,
			     const int * mapping_table,
			     const double symprec)
{
  Cell *bravais;
  int i, num_prim_sym;
  int *wyckoffs_bravais, *equiv_atoms_bravais;
  int operation_index[2];

  debug_print("get_Wyckoff_positions\n");

  bravais = NULL;
  wyckoffs_bravais = NULL;
  equiv_atoms_bravais = NULL;

  if ((wyckoffs_bravais = (int*)malloc(sizeof(int) * primitive->size * 4))
      == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  if ((equiv_atoms_bravais = (int*)malloc(sizeof(int) * primitive->size * 4))
      == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    free(wyckoffs_bravais);
    wyckoffs_bravais = NULL;
    return NULL;
  }
  
  if ((bravais = get_bravais_exact_positions_and_lattice
       (wyckoffs_bravais,
	equiv_atoms_bravais,
	spacegroup,
	primitive,
	symprec)) == NULL) {
    goto ret;
  }

  for (i = 0; i < cell->size; i++) {
    wyckoffs[i] = wyckoffs_bravais[mapping_table[i]];
  }
  
  spgdb_get_operation_index(operation_index, spacegroup->hall_number);
  num_prim_sym = operation_index[0] / (bravais->size / primitive->size);

  /* Check symmetry breaking by unusual multiplicity of primitive cell. */
  if (cell->size * num_prim_sym != symmetry->size * primitive->size) {
    set_equivalent_atoms_broken_symmetry(equiv_atoms,
					 cell,
					 symmetry,
					 mapping_table,
					 symprec);
  } else {
    if (set_equivalent_atoms(equiv_atoms,
			     primitive,
			     cell,
			     equiv_atoms_bravais,
			     mapping_table) == 0) {
      cel_free_cell(bravais);
      bravais = NULL;
    }
  }

 ret:  
  free(equiv_atoms_bravais);
  equiv_atoms_bravais = NULL;
  free(wyckoffs_bravais);
  wyckoffs_bravais = NULL;
  
  return bravais;
}

/* Only the atoms corresponding to those in primitive are returned. */
/* Return NULL if failed */
static Cell *
get_bravais_exact_positions_and_lattice(int * wyckoffs,
					int * equiv_atoms,
					SPGCONST Spacegroup *spacegroup,
					SPGCONST Cell * primitive,
					const double symprec)
{
  int i;
  int *wyckoffs_prim, *equiv_atoms_prim;
  Symmetry *conv_sym;
  Cell *bravais, *conv_prim;
  VecDBL *exact_positions;

  debug_print("get_bravais_exact_positions_and_lattice\n");

  wyckoffs_prim = NULL;
  equiv_atoms_prim = NULL;
  conv_prim = NULL;
  bravais = NULL;
  conv_sym = NULL;
  exact_positions = NULL;

  /* Symmetrize atomic positions of conventional unit cell */
  if ((wyckoffs_prim = (int*)malloc(sizeof(int) * primitive->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return NULL;
  }

  if ((equiv_atoms_prim = (int*)malloc(sizeof(int) * primitive->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    free(wyckoffs_prim);
    wyckoffs_prim = NULL;
    return NULL;
  }

  for (i = 0; i < primitive->size; i++) {
    wyckoffs_prim[i] = -1;
    equiv_atoms_prim[i] = -1;
  }

  /* Positions of primitive atoms are represented wrt Bravais lattice */
  if ((conv_prim = get_conventional_primitive(spacegroup, primitive)) == NULL) {
    free(wyckoffs_prim);
    wyckoffs_prim = NULL;
    free(equiv_atoms_prim);
    equiv_atoms_prim = NULL;
    return NULL;
  }

  /* Symmetries in database (wrt Bravais lattice) */
  if ((conv_sym = spgdb_get_spacegroup_operations(spacegroup->hall_number))
      == NULL) {
    goto err;
  }

  /* Lattice vectors are set. */
  get_conventional_lattice(conv_prim->lattice, spacegroup);

  if ((exact_positions = ssm_get_exact_positions(wyckoffs_prim,
						 equiv_atoms_prim,
						 conv_prim,
						 conv_sym,
						 spacegroup->hall_number,
						 symprec)) == NULL) {
    sym_free_symmetry(conv_sym);
    goto err;
  }

  for (i = 0; i < conv_prim->size; i++) {
    mat_copy_vector_d3(conv_prim->position[i], exact_positions->vec[i]);
  }

  bravais = expand_positions(wyckoffs,
			     equiv_atoms,
			     conv_prim,
			     conv_sym,
			     wyckoffs_prim,
			     equiv_atoms_prim);

  mat_free_VecDBL(exact_positions);
  sym_free_symmetry(conv_sym);
 err:
  free(wyckoffs_prim);
  wyckoffs_prim = NULL;
  free(equiv_atoms_prim);
  equiv_atoms_prim = NULL;
  cel_free_cell(conv_prim);

  return bravais;
}

/* Return NULL if failed */
static Cell * expand_positions(int * wyckoffs,
			       int * equiv_atoms,
			       SPGCONST Cell * conv_prim,
			       SPGCONST Symmetry * conv_sym,
			       const int * wyckoffs_prim,
			       const int * equiv_atoms_prim)
{
  int i, j, k, num_pure_trans;
  int num_atom;
  Cell * bravais;

  bravais = NULL;

  num_pure_trans = get_number_of_pure_translation(conv_sym);

  if ((bravais = cel_alloc_cell(conv_prim->size * num_pure_trans)) == NULL) {
    return NULL;
  }

  num_atom = 0;
  for (i = 0; i < conv_sym->size; i++) {
    /* Referred atoms in Bravais lattice */
    if (mat_check_identity_matrix_i3(identity, conv_sym->rot[i])) {
      for (j = 0; j < conv_prim->size; j++) {
	bravais->types[num_atom] = conv_prim->types[j];
	mat_copy_vector_d3(bravais->position[ num_atom ],
			   conv_prim->position[j]);
	for (k = 0; k < 3; k++) {
	  bravais->position[num_atom][k] += conv_sym->trans[i][k];
	  bravais->position[num_atom][k] = 
	    mat_Dmod1(bravais->position[num_atom][k]);
	}
	wyckoffs[num_atom] = wyckoffs_prim[j];
	equiv_atoms[num_atom] = equiv_atoms_prim[j];
	num_atom++;
      }
    }
  }

  mat_copy_matrix_d3(bravais->lattice, conv_prim->lattice);

  return bravais;
}


static int get_number_of_pure_translation(SPGCONST Symmetry * conv_sym)
{
  int i, num_pure_trans = 0;
  
  for (i = 0; i < conv_sym->size; i++) {
    if (mat_check_identity_matrix_i3(identity, conv_sym->rot[i])) {
      num_pure_trans++;
    }
  }

  return num_pure_trans;
}

static Cell * get_conventional_primitive(SPGCONST Spacegroup * spacegroup,
					 SPGCONST Cell * primitive)
{
  int i, j;
  double inv_brv[3][3], trans_mat[3][3];
  Cell * conv_prim;

  conv_prim = NULL;

  if ((conv_prim = cel_alloc_cell(primitive->size)) == NULL) {
    return NULL;
  }

  mat_inverse_matrix_d3(inv_brv, spacegroup->bravais_lattice, 0);
  mat_multiply_matrix_d3(trans_mat, inv_brv, primitive->lattice);
  
  for (i = 0; i < primitive->size; i++) {
    conv_prim->types[i] = primitive->types[i];
    mat_multiply_matrix_vector_d3(conv_prim->position[i],
				  trans_mat,
				  primitive->position[i]);
    for (j = 0; j < 3; j++) {
      conv_prim->position[i][j] += spacegroup->origin_shift[j];
      conv_prim->position[i][j] = mat_Dmod1(conv_prim->position[i][j]);
    }
  }

  return conv_prim;
}

static void get_conventional_lattice(double lattice[3][3],
				     SPGCONST Spacegroup *spacegroup)
{
  int i, j;
  double metric[3][3];
  Pointgroup pointgroup;

  pointgroup = ptg_get_pointgroup(spacegroup->pointgroup_number);
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lattice[i][j] = 0;
    }
  }

  mat_get_metric(metric, spacegroup->bravais_lattice);

  debug_print("bravais lattice\n");
  debug_print_matrix_d3(spacegroup->bravais_lattice);
  debug_print("%s\n", spacegroup->setting);

  switch (pointgroup.holohedry) {
  case TRICLI:
    set_tricli(lattice, metric);
    break;
  case MONOCLI: /* b-axis is the unique axis. */
    set_monocli(lattice, metric);
    break;
  case ORTHO:
    set_ortho(lattice, metric);
    break;
  case TETRA:
    set_tetra(lattice, metric);
    break;
  case TRIGO:
    if (spacegroup->setting[0] == 'R') {
      set_rhomb(lattice, metric);
    } else {
      set_trigo(lattice, metric);
    }
    break;
  case HEXA:
    set_trigo(lattice, metric);
    break;
  case CUBIC:
    set_cubic(lattice, metric);
    break;
  case HOLOHEDRY_NONE:
    break;
  }
}

/* The conversion refers the wikipedia, */
/* http://en.wikipedia.org/wiki/Fractional_coordinates */
static void set_tricli(double lattice[3][3],
		       SPGCONST double metric[3][3])
{
  double a, b, c, alpha, beta, gamma, cg, cb, ca, sg;

  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  alpha = acos(metric[1][2] / b / c);
  beta = acos(metric[0][2] / a / c);
  gamma = acos(metric[0][1] / a / b);

  cg = cos(gamma);
  cb = cos(beta);
  ca = cos(alpha);
  sg = sin(gamma);

  lattice[0][0] = a;
  lattice[0][1] = b * cg;
  lattice[0][2] = c * cb;
  lattice[1][1] = b * sg;
  lattice[1][2] = c * (ca - cb * cg) / sg;
  lattice[2][2] = c * sqrt(1 - ca * ca - cb * cb - cg * cg +
			   2 * ca * cb * cg) / sg;
}

static void set_monocli(double lattice[3][3],
			SPGCONST double metric[3][3])
{
  /* Lattice is expected to be C centring */
  double a, b, c, beta;

  debug_print("set_monocli:\n");
  debug_print_matrix_d3(metric);

  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  lattice[0][0] = a;
  lattice[1][1] = b;
  beta = acos(metric[0][2] / a / c);
  lattice[0][2] = c * cos(beta);
  lattice[2][2] = c * sin(beta);

  debug_print("beta %f\n", beta);
  debug_print_matrix_d3(lattice);
}
			 
static void set_ortho(double lattice[3][3],
		      SPGCONST double metric[3][3])
{
  double a, b, c;
  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  lattice[0][0] = a;
  lattice[1][1] = b;
  lattice[2][2] = c;
}

static void set_tetra(double lattice[3][3],
		      SPGCONST double metric[3][3])
{
  double a, b, c;
  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  lattice[0][0] = (a + b) / 2;
  lattice[1][1] = (a + b) / 2;
  lattice[2][2] = c;
}

static void set_rhomb(double lattice[3][3],
		      SPGCONST double metric[3][3])
{
  double a, b, c, angle, ahex, chex;

  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  angle = acos((metric[0][1] / a / b +
		metric[0][2] / a / c +
		metric[1][2] / b / c) / 3);

  /* Reference, http://cst-www.nrl.navy.mil/lattice/struk/rgr.html */
  ahex = 2 * (a+b+c)/3 * sin(angle / 2);
  chex = (a+b+c)/3 * sqrt(3 * (1 + 2 * cos(angle))) ;
  lattice[0][0] = ahex / 2;
  lattice[1][0] = -ahex / (2 * sqrt(3));
  lattice[2][0] = chex / 3;
  lattice[1][1] = ahex / sqrt(3);
  lattice[2][1] = chex / 3;
  lattice[0][2] = -ahex / 2;
  lattice[1][2] = -ahex / (2 * sqrt(3));
  lattice[2][2] = chex / 3;


#ifdef SPGDEBUG
  debug_print("Rhombo lattice: %f %f %f %f %f %f\n", a, b, c,
	      acos(metric[0][1] / a / b) / 3.14 * 180,
	      acos(metric[0][2] / a / c) / 3.14 * 180,
	      acos(metric[1][2] / b / c) / 3.14 * 180);
  double dmetric[3][3];
  mat_get_metric(dmetric, lattice);
  a = sqrt(dmetric[0][0]);
  b = sqrt(dmetric[1][1]);
  c = sqrt(dmetric[2][2]);
  debug_print("Rhombo lattice symmetrized: %f %f %f %f %f %f\n",
	      a, b, c,
	      acos(dmetric[0][1] / a / b) / 3.14 * 180,
	      acos(dmetric[0][2] / a / c) / 3.14 * 180,
	      acos(dmetric[1][2] / b / c) / 3.14 * 180);
#endif
}

static void set_trigo(double lattice[3][3],
		      SPGCONST double metric[3][3])
{
  double a, b, c;

  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  lattice[0][0] = (a + b) / 2;
  lattice[0][1] = - (a + b) / 4;
  lattice[1][1] = (a + b) / 4 * sqrt(3);
  lattice[2][2] = c;
}

static void set_cubic(double lattice[3][3],
		      SPGCONST double metric[3][3])
{
  double a, b, c;

  a = sqrt(metric[0][0]);
  b = sqrt(metric[1][1]);
  c = sqrt(metric[2][2]);
  lattice[0][0] = (a + b + c) / 3;
  lattice[1][1] = (a + b + c) / 3;
  lattice[2][2] = (a + b + c) / 3;
}

/* Return NULL if failed */
static Symmetry *
get_refined_symmetry_operations(SPGCONST Cell * cell,
				SPGCONST Cell * primitive,
				SPGCONST Spacegroup * spacegroup,
				const double symprec)
{
  int t_mat_int[3][3];
  int frame[3];
  double inv_prim_lat[3][3], t_mat[3][3];
  Symmetry *conv_sym, *prim_sym, *symmetry;

  conv_sym = NULL;
  prim_sym = NULL;
  symmetry = NULL;

  /* Primitive symmetry from database */
  if ((conv_sym = spgdb_get_spacegroup_operations(spacegroup->hall_number))
      == NULL) {
    return NULL;
  }

  mat_inverse_matrix_d3(inv_prim_lat, primitive->lattice, 0);
  mat_multiply_matrix_d3(t_mat, inv_prim_lat, spacegroup->bravais_lattice);

  set_translation_with_origin_shift(conv_sym, spacegroup->origin_shift);

  if ((prim_sym = get_primitive_db_symmetry(t_mat, conv_sym)) == NULL) {
    sym_free_symmetry(conv_sym);
    return NULL;
  }

  sym_free_symmetry(conv_sym);

  /* Input cell symmetry from primitive symmetry */
  mat_multiply_matrix_d3(t_mat, inv_prim_lat, cell->lattice);
  mat_cast_matrix_3d_to_3i(t_mat_int, t_mat);
  get_surrounding_frame(frame, t_mat_int);

  symmetry = reduce_symmetry_in_frame(frame,
				      prim_sym,
				      t_mat_int,
				      cell->lattice,
				      cell->size / primitive->size,
				      symprec);

  sym_free_symmetry(prim_sym);

  return symmetry;
}

static int set_equivalent_atoms(int * equiv_atoms_cell,
				SPGCONST Cell * primitive,
				SPGCONST Cell * cell,
				const int * equiv_atoms_prim,
				const int * mapping_table)
{
  int i, j;
  int *equiv_atoms;

  equiv_atoms = NULL;

  if ((equiv_atoms = (int*) malloc(sizeof(int) * primitive->size)) == NULL) {
    warning_print("spglib: Memory could not be allocated ");
    return 0;
  }

  for (i = 0; i < primitive->size; i++) {
    for (j = 0; j < cell->size; j++) {
      if (mapping_table[j] == equiv_atoms_prim[i]) {
	equiv_atoms[i] = j;
	break;
      }
    }
  }
  for (i = 0; i < cell->size; i++) {
    equiv_atoms_cell[i] = equiv_atoms[mapping_table[i]];
  }
  free(equiv_atoms);
  equiv_atoms = NULL;

  return 1;
}

static void set_equivalent_atoms_broken_symmetry(int * equiv_atoms_cell,
						 SPGCONST Cell * cell,
						 const Symmetry *symmetry,
						 const int * mapping_table,
						 const double symprec)
{
  int i, j;

  for (i = 0; i < cell->size; i++) {
    equiv_atoms_cell[i] = i;
    for (j = 0; j < cell->size; j++) {
      if (mapping_table[i] == mapping_table[j]) {
	if (i == j) {
	  equiv_atoms_cell[i] =
	    equiv_atoms_cell[search_equivalent_atom(i,
						    cell,
						    symmetry,
						    symprec)];
	} else {
	  equiv_atoms_cell[i] = equiv_atoms_cell[j];
	}
	break;
      }
    }
  }
}

static int search_equivalent_atom(const int atom_index,
				  SPGCONST Cell * cell,
				  const Symmetry *symmetry,
				  const double symprec)
{
  int i, j;
  double pos_rot[3];

  for (i = 0; i < symmetry->size; i++) {
    mat_multiply_matrix_vector_id3(pos_rot,
				   symmetry->rot[i],
				   cell->position[atom_index]);
    for (j = 0; j < 3; j++) {
      pos_rot[j] += symmetry->trans[i][j];
    }
    for (j = 0; j < atom_index; j++) {
      if (cel_is_overlap(cell->position[j], pos_rot, cell->lattice, symprec)) {
	return j;
      }
    }
  }
  return atom_index;
}


static void set_translation_with_origin_shift(Symmetry *conv_sym,
					      const double origin_shift[3])
{
  int i, j;
  double tmp_vec[3];
  int tmp_mat[3][3];

  /* t' = t - (R-E)w, w is the origin shift */
  for (i = 0; i < conv_sym->size; i++) {
    mat_copy_matrix_i3(tmp_mat, conv_sym->rot[i]);
    tmp_mat[0][0]--;
    tmp_mat[1][1]--;
    tmp_mat[2][2]--;
    mat_multiply_matrix_vector_id3(tmp_vec, tmp_mat, origin_shift);
    for (j = 0; j < 3; j++) {
      conv_sym->trans[i][j] += tmp_vec[j];
    }
  }
}

static Symmetry * get_primitive_db_symmetry(SPGCONST double t_mat[3][3],
					    const Symmetry *conv_sym)
{
  int i, j, num_op;
  double inv_mat[3][3], tmp_mat[3][3];
  MatINT *r_prim;
  VecDBL *t_prim;
  Symmetry *prim_sym;

  r_prim = NULL;
  t_prim = NULL;
  prim_sym = NULL;
  
  if ((r_prim = mat_alloc_MatINT(conv_sym->size)) == NULL) {
    return NULL;
  }

  if ((t_prim = mat_alloc_VecDBL(conv_sym->size)) == NULL) {
    mat_free_MatINT(r_prim);
    return NULL;
  }

  mat_inverse_matrix_d3(inv_mat, t_mat, 0);

  num_op = 0;
  for (i = 0; i < conv_sym->size; i++) {
    for (j = 0; j < i; j++) {
      if (mat_check_identity_matrix_i3(conv_sym->rot[i],
				       conv_sym->rot[j])) {
	goto pass;
      }
    }

    /* R' = T*R*T^-1 */
    mat_multiply_matrix_di3(tmp_mat, t_mat, conv_sym->rot[i]);
    mat_multiply_matrix_d3(tmp_mat, tmp_mat, inv_mat);
    mat_cast_matrix_3d_to_3i(r_prim->mat[ num_op ], tmp_mat);
    /* t' = T*t */
    mat_multiply_matrix_vector_d3(t_prim->vec[num_op],
				  t_mat,
				  conv_sym->trans[i]);
    num_op++;

  pass:
    ;
  }

  if ((prim_sym = sym_alloc_symmetry(num_op)) == NULL) {
    goto ret;
  }

  for (i = 0; i < num_op; i++) {
    mat_copy_matrix_i3(prim_sym->rot[i], r_prim->mat[i]);
    for (j = 0; j < 3; j++) {
      prim_sym->trans[i][j] = mat_Dmod1(t_prim->vec[i][j]);
    }
  }

 ret:
  mat_free_MatINT(r_prim);
  mat_free_VecDBL(t_prim);

  return prim_sym;
}

static void get_surrounding_frame(int frame[3],
				  SPGCONST int t_mat[3][3])
{
  int i, j, max, min;
  int corners[3][8];

  get_corners(corners, t_mat);

  for (i = 0; i < 3; i++) {
    max = corners[i][0];
    min = corners[i][0];
    for (j = 1; j < 8; j++) {
      if (max < corners[i][j]) {
	max = corners[i][j];
      }
      if (min > corners[i][j]) {
	min = corners[i][j];
      }
    }
    frame[i] = max - min;
  }
}

static void get_corners(int corners[3][8],
			SPGCONST int t_mat[3][3])
{
  int i, j;

  /* O */
  for (i = 0; i < 3; i++) {
    corners[i][0] = 0;
  }

  /* a,b,c */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      corners[j][i+1] = t_mat[j][i];
    }
  }

  /* b+c, c+a, a+b */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      corners[j][i+4] = t_mat[j][(i+1) % 3] + t_mat[j][(i+2) % 3];
    }
  }

  /* a+b+c */
  for (i = 0; i < 3; i++) {
    corners[i][7] = t_mat[i][0] + t_mat[i][1] + t_mat[i][2];
  }
}

static Symmetry * reduce_symmetry_in_frame(const int frame[3],
					   SPGCONST Symmetry *prim_sym,
					   SPGCONST int t_mat[3][3],
					   SPGCONST double lattice[3][3],
					   const int multiplicity,
					   const double symprec)
{
  Symmetry *symmetry, *t_sym;
  double inv_tmat[3][3], tmp_mat[3][3];
  VecDBL *pure_trans, *lattice_trans;

  symmetry = NULL;
  t_sym = NULL;
  pure_trans = NULL;
  lattice_trans = NULL;

  mat_cast_matrix_3i_to_3d(tmp_mat, t_mat);
  mat_inverse_matrix_d3(inv_tmat, tmp_mat, 0);

  if ((lattice_trans = get_lattice_translations(frame, inv_tmat)) == NULL) {
    return NULL;
  }

  if ((pure_trans = remove_overlapping_lattice_points(lattice,
						      lattice_trans,
						      symprec)) == NULL) {
    mat_free_VecDBL(lattice_trans);
    return NULL;
  }

  if ((t_sym = get_symmetry_in_original_cell(t_mat,
					     inv_tmat,
					     lattice,
					     prim_sym,
					     symprec)) == NULL) {
    mat_free_VecDBL(pure_trans);
    mat_free_VecDBL(lattice_trans);
    return NULL;
  }

  if (pure_trans->size == multiplicity) {
    symmetry = copy_symmetry_upon_lattice_points(pure_trans, t_sym);
  }

  mat_free_VecDBL(lattice_trans);
  mat_free_VecDBL(pure_trans);
  sym_free_symmetry(t_sym);

  return symmetry;
}

/* Return NULL if failed */
static VecDBL * get_lattice_translations(const int frame[3],
					 SPGCONST double inv_tmat[3][3])
{
  int i, j, k, l, num_trans;
  VecDBL * lattice_trans;

  lattice_trans = NULL;

  if ((lattice_trans = mat_alloc_VecDBL(frame[0] * frame[1] * frame[2]))
      == NULL) {
    return NULL;
  }

  num_trans = 0;
  for (i = 0; i < frame[0]; i++) {
    for (j = 0; j < frame[1]; j++) {
      for (k = 0; k < frame[2]; k++) {
	lattice_trans->vec[num_trans][0] = i;
	lattice_trans->vec[num_trans][1] = j;
	lattice_trans->vec[num_trans][2] = k;

	/* t' = T^-1*t */
	mat_multiply_matrix_vector_d3(lattice_trans->vec[num_trans],
				      inv_tmat,
				      lattice_trans->vec[num_trans]);
	for (l = 0; l < 3; l++) {
	  lattice_trans->vec[num_trans][l] =
	    mat_Dmod1(lattice_trans->vec[num_trans][l]); 
	}
	num_trans++;
      }
    }
  }

  return lattice_trans;
}

static VecDBL *
remove_overlapping_lattice_points(SPGCONST double lattice[3][3],
				  const VecDBL *lattice_trans,
				  const double symprec)
{
  int i, j, is_found, num_pure_trans;
  VecDBL *pure_trans, *t;
  
  pure_trans = NULL;
  t = NULL;

  num_pure_trans = 0;

  if ((t = mat_alloc_VecDBL(lattice_trans->size)) == NULL) {
    return NULL;
  }

  for (i = 0; i < lattice_trans->size; i++) {
    is_found = 0;
    for (j = 0; j < num_pure_trans; j++) {
      if (cel_is_overlap(lattice_trans->vec[i], t->vec[j], lattice, symprec)) {
	is_found = 1;
	break;
      }
    }
    if (! is_found) {
      mat_copy_vector_d3(t->vec[num_pure_trans], lattice_trans->vec[i]);
      num_pure_trans++;
    }
  }

  if ((pure_trans = mat_alloc_VecDBL(num_pure_trans)) == NULL) {
    mat_free_VecDBL(t);
    return NULL;
  }

  for (i = 0; i < num_pure_trans; i++) {
    mat_copy_vector_d3(pure_trans->vec[i], t->vec[i]);
  }
  mat_free_VecDBL(t);

  return pure_trans;
}

/* Return NULL if failed */
static Symmetry *
get_symmetry_in_original_cell(SPGCONST int t_mat[3][3],
			      SPGCONST double inv_tmat[3][3],
			      SPGCONST double lattice[3][3],
			      SPGCONST Symmetry *prim_sym,
			      const double symprec)
{				    
  int i, size_sym_orig;
  double tmp_rot_d[3][3], tmp_lat_d[3][3], tmp_lat_i[3][3], tmp_mat[3][3];
  int tmp_rot_i[3][3];
  Symmetry *t_sym, *t_red_sym;

  t_sym = NULL;
  t_red_sym = NULL;

  if ((t_sym = sym_alloc_symmetry(prim_sym->size)) == NULL) {
    return NULL;
  }

  /* transform symmetry operations of primitive cell to those of original */
  size_sym_orig = 0;
  for (i = 0; i < prim_sym->size; i++) {
    /* R' = T^-1*R*T */
    mat_multiply_matrix_di3(tmp_mat, inv_tmat, prim_sym->rot[i]);
    mat_multiply_matrix_di3(tmp_rot_d, tmp_mat, t_mat);

    /* In spglib, symmetry of supercell is defined by the set of symmetry */
    /* operations that are searched among supercell lattice point group */
    /* operations. The supercell lattice may be made by breaking the */
    /* unit cell lattice symmetry. In this case, a part of symmetry */
    /* operations is discarded. */
    mat_cast_matrix_3d_to_3i(tmp_rot_i, tmp_rot_d);
    mat_multiply_matrix_di3(tmp_lat_i, lattice, tmp_rot_i);
    mat_multiply_matrix_d3(tmp_lat_d, lattice, tmp_rot_d);
    if (mat_check_identity_matrix_d3(tmp_lat_i, tmp_lat_d, symprec)) {
      mat_copy_matrix_i3(t_sym->rot[size_sym_orig], tmp_rot_i);
      /* t' = T^-1*t */
      mat_multiply_matrix_vector_d3(t_sym->trans[size_sym_orig],
				    inv_tmat,
				    prim_sym->trans[i]);
      size_sym_orig++;
    }
  }

  /* Broken symmetry due to supercell multiplicity */
  if (size_sym_orig != prim_sym->size) {

    if ((t_red_sym = sym_alloc_symmetry(size_sym_orig)) == NULL) {
      sym_free_symmetry(t_sym);
      return NULL;
    }

    for (i = 0; i < size_sym_orig; i++) {
      mat_copy_matrix_i3(t_red_sym->rot[i], t_sym->rot[i]);
      mat_copy_vector_d3(t_red_sym->trans[i], t_sym->trans[i]);
    }

    sym_free_symmetry(t_sym);

    t_sym = t_red_sym;
    t_red_sym = NULL;
  }

  return t_sym;
}

/* Return NULL if failed */
static Symmetry *
copy_symmetry_upon_lattice_points(const VecDBL *pure_trans,
				  SPGCONST Symmetry *t_sym)
{
  int i, j, k, size_sym_orig;
  Symmetry *symmetry;

  symmetry = NULL;

  size_sym_orig = t_sym->size;

  if ((symmetry = sym_alloc_symmetry(pure_trans->size * size_sym_orig))
      == NULL) {
    return NULL;
  }

  for (i = 0; i < pure_trans->size; i++) {
    for (j = 0; j < size_sym_orig; j++) {
      mat_copy_matrix_i3(symmetry->rot[size_sym_orig * i + j], t_sym->rot[j]);
      mat_copy_vector_d3(symmetry->trans[size_sym_orig * i + j],
			 t_sym->trans[j]);
      for (k = 0; k < 3; k++) {
	symmetry->trans[size_sym_orig * i + j][k] += pure_trans->vec[i][k];
	symmetry->trans[size_sym_orig * i + j][k] =
	  mat_Dmod1(symmetry->trans[size_sym_orig * i + j][k]);
      }
    }
  }

  return symmetry;
}
