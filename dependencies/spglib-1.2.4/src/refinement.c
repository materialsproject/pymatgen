/* bravais.c */
/* Copyright (C) 2011 Atsushi Togo */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "refinement.h"
#include "cell.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "primitive.h"
#include "spacegroup.h"
#include "spg_database.h"
#include "site_symmetry.h"
#include "symmetry.h"

#include "debug.h"

#define REDUCE_RATE 0.95

static Cell * refine_cell(SPGCONST Cell * cell,
			  const double symprec);
static Cell * get_bravais_exact_positions_and_lattice(int * wyckoffs,
						      int * equiv_atoms,
						      SPGCONST Spacegroup * spacegroup,
						      SPGCONST Cell * primitive,
						      const double symprec);
static Cell * expand_positions(SPGCONST Cell * conv_prim,
			       SPGCONST Symmetry * conv_sym);
static Cell * get_conventional_primitive(SPGCONST Spacegroup * spacegroup,
					 SPGCONST Cell * primitive);
static Symmetry * get_db_symmetry(const int hall_number);
static int get_number_of_pure_translation(SPGCONST Symmetry * conv_sym);
static int get_conventional_lattice(double lattice[3][3],
				    const Holohedry holohedry,
				    SPGCONST double bravais_lattice[3][3]);
static void set_monocli(double lattice[3][3],
			SPGCONST double metric[3][3]);
static void set_ortho(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_tetra(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_rhomb(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_trigo(double lattice[3][3],
		      SPGCONST double metric[3][3]);
static void set_cubic(double lattice[3][3],
		      SPGCONST double metric[3][3]);

static Symmetry * get_refined_symmetry_operations(SPGCONST Cell * cell,
						  SPGCONST Cell * primitive,
						  SPGCONST Spacegroup * spacegroup,
						  const double symprec);
static void set_translation_with_origin_shift(Symmetry *conv_sym,
					      const double origin_shift[3]);
static Symmetry * get_primitive_db_symmetry(SPGCONST double t_mat[3][3],
					    const Symmetry *conv_sym,
					    const double symprec);
static void get_corners(int corners[3][8],
			SPGCONST int t_mat[3][3]);
static void get_surrounding_frame(int frame[3],
				  SPGCONST int t_mat[3][3]);
static Symmetry * reduce_symmetry_in_frame(const int frame[3],
					   SPGCONST Symmetry *prim_sym,
					   SPGCONST int t_mat[3][3],
					   SPGCONST double lattice[3][3],
					   const int multiplicity,
					   const double symprec);
static VecDBL * reduce_lattice_points(SPGCONST double lattice[3][3],
				      const VecDBL *lattice_trans,
				      const double symprec);


static SPGCONST int identity[3][3] = {
  { 1, 0, 0},
  { 0, 1, 0},
  { 0, 0, 1},
};


Cell * ref_refine_cell(SPGCONST Cell * cell,
		       const double symprec)
{
  return refine_cell(cell, symprec);
}

/* symmetry->size = 0 is returned when it failed. */
Symmetry * ref_get_refined_symmetry_operations(SPGCONST Cell * cell,
					       SPGCONST Cell * primitive,
					       SPGCONST Spacegroup * spacegroup,
					       const double symprec)
{
  return get_refined_symmetry_operations(cell, primitive, spacegroup, symprec);
}

void ref_get_Wyckoff_positions(int * wyckoffs,
			       int * equiv_atoms,
			       SPGCONST Cell * primitive,
			       SPGCONST Spacegroup * spacegroup,
			       const double symprec)
{
  Cell *conv_prim;

  conv_prim = get_bravais_exact_positions_and_lattice(wyckoffs,
						      equiv_atoms,
						      spacegroup,
						      primitive,
						      symprec);
  cel_free_cell(conv_prim);
}

static Cell * refine_cell(SPGCONST Cell * cell,
			  const double symprec)
{
  int *wyckoffs, *equiv_atoms;
  double tolerance;
  Cell *primitive, *bravais, *conv_prim;
  Symmetry *conv_sym;
  Spacegroup spacegroup;

  debug_print("refine_cell:\n");
  
  primitive = prm_get_primitive(cell, symprec);
  
  if (primitive->size == 0) {
    cel_free_cell(primitive);
    bravais = cel_alloc_cell(0);
    goto end;
  }

  tolerance = prm_get_current_tolerance();
  spacegroup = spa_get_spacegroup_with_primitive(primitive, tolerance);

  wyckoffs = (int*)malloc(sizeof(int) * primitive->size);
  equiv_atoms = (int*)malloc(sizeof(int) * primitive->size);
  conv_prim = get_bravais_exact_positions_and_lattice(wyckoffs,
						      equiv_atoms,
						      &spacegroup,
						      primitive,
						      tolerance);
  free(equiv_atoms);
  equiv_atoms = NULL;
  free(wyckoffs);
  wyckoffs = NULL;

  conv_sym = get_db_symmetry(spacegroup.hall_number);
  bravais = expand_positions(conv_prim, conv_sym);

  debug_print("primitive cell in refine_cell:\n");
  debug_print_matrix_d3(primitive->lattice);
  debug_print("conventional lattice in refine_cell:\n");
  debug_print_matrix_d3(conv_prim->lattice);
  debug_print("bravais lattice in refine_cell:\n");
  debug_print_matrix_d3(bravais->lattice);
  
  cel_free_cell(conv_prim);
  sym_free_symmetry(conv_sym);
  cel_free_cell(primitive);

 end:  /* Return bravais->size = 0, if the bravais could not be found. */
  return bravais;
}

/* Only the atoms corresponding to those in primitive are returned. */
static Cell * get_bravais_exact_positions_and_lattice(int * wyckoffs,
						      int * equiv_atoms,
						      SPGCONST Spacegroup *spacegroup,
						      SPGCONST Cell * primitive,
						      const double symprec)
{
  int i;
  Symmetry *conv_sym;
  Cell *bravais;
  VecDBL *exact_positions;

  /* Positions of primitive atoms are represented wrt Bravais lattice */
  bravais = get_conventional_primitive(spacegroup, primitive);
  /* Symmetries in database (wrt Bravais lattice) */
  conv_sym = get_db_symmetry(spacegroup->hall_number);
  /* Lattice vectors are set. */
  get_conventional_lattice(bravais->lattice,
			   spacegroup->holohedry,
			   spacegroup->bravais_lattice);

  /* Symmetrize atomic positions of conventional unit cell */
  exact_positions = ssm_get_exact_positions(wyckoffs,
					    equiv_atoms,
					    bravais,
					    conv_sym,
					    spacegroup->hall_number,
					    symprec);
  sym_free_symmetry(conv_sym);

  if (exact_positions->size > 0) {
    for (i = 0; i < bravais->size; i++) {
      mat_copy_vector_d3(bravais->position[i], exact_positions->vec[i]);
    }
  } else {
    cel_free_cell(bravais);
    bravais = cel_alloc_cell(0);
  }

  mat_free_VecDBL(exact_positions);

  return bravais;
}

static Cell * expand_positions(SPGCONST Cell * conv_prim,
			       SPGCONST Symmetry * conv_sym)
{
  int i, j, k, num_pure_trans;
  int num_atom;
  Cell * bravais;

  num_pure_trans = get_number_of_pure_translation(conv_sym);
  bravais = cel_alloc_cell(conv_prim->size * num_pure_trans);

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

  conv_prim = cel_alloc_cell(primitive->size);

  mat_inverse_matrix_d3(inv_brv, spacegroup->bravais_lattice, 0);
  mat_multiply_matrix_d3(trans_mat, inv_brv, primitive->lattice);
  
  for (i = 0; i < primitive->size; i++) {
    conv_prim->types[i] = primitive->types[i];
    mat_multiply_matrix_vector_d3(conv_prim->position[i],
				  trans_mat,
				  primitive->position[i]);
    for (j = 0; j < 3; j++) {
      conv_prim->position[i][j] -= spacegroup->origin_shift[j];
      conv_prim->position[i][j] -= mat_Nint(conv_prim->position[i][j]);
    }
  }

  return conv_prim;
}

static Symmetry * get_db_symmetry(const int hall_number)
{
  int i;
  int operation_index[2];
  int rot[3][3];
  double trans[3];
  Symmetry *symmetry;

  spgdb_get_operation_index(operation_index, hall_number);
  symmetry = sym_alloc_symmetry(operation_index[0]);

  for (i = 0; i < operation_index[0]; i++) {
    /* rotation matrix matching and set difference of translations */
    spgdb_get_operation(rot, trans, operation_index[1] + i);
    mat_copy_matrix_i3(symmetry->rot[i], rot);
    mat_copy_vector_d3(symmetry->trans[i], trans);
  }

  return symmetry;
}

static int get_conventional_lattice(double lattice[3][3],
				    const Holohedry holohedry,
				    SPGCONST double bravais_lattice[3][3])
{
  int i, j;
  double metric[3][3];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      lattice[i][j] = 0;
    }
  }

  mat_get_metric(metric, bravais_lattice);

  switch (holohedry) {
  case TRICLI:
    mat_copy_matrix_d3(lattice, bravais_lattice);
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
  case RHOMB:
    set_rhomb(lattice, metric);
    break;
  case TRIGO:
    set_trigo(lattice, metric);
    break;
  case HEXA:
    set_trigo(lattice, metric);
    break;
  case CUBIC:
    set_cubic(lattice, metric);
    break;
  case NONE:
    break;
  }

  return 1;
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
  lattice[1][1] = b;
  lattice[2][2] = c;
  beta = acos(metric[0][2] / a / c);
  lattice[2][0] = a * cos(beta);
  lattice[0][0] = a * sin(beta);

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


#ifdef DEBUG
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





static Symmetry * get_refined_symmetry_operations(SPGCONST Cell * cell,
						  SPGCONST Cell * primitive,
						  SPGCONST Spacegroup * spacegroup,
						  const double symprec)
{
  int t_mat_int[3][3];
  int frame[3];
  double inv_mat[3][3], t_mat[3][3];
  Symmetry *conv_sym, *prim_sym, *symmetry;

  /* Primitive symmetry from database */
  conv_sym = get_db_symmetry(spacegroup->hall_number);
  set_translation_with_origin_shift(conv_sym, spacegroup->origin_shift);
  mat_inverse_matrix_d3(inv_mat, primitive->lattice, symprec);
  mat_multiply_matrix_d3(t_mat, inv_mat, spacegroup->bravais_lattice);
  prim_sym = get_primitive_db_symmetry(t_mat, conv_sym, symprec);

  /* Input cell symmetry from primitive symmetry */
  mat_inverse_matrix_d3(inv_mat, primitive->lattice, symprec);
  mat_multiply_matrix_d3(t_mat, inv_mat, cell->lattice);
  mat_cast_matrix_3d_to_3i(t_mat_int, t_mat);
  get_surrounding_frame(frame, t_mat_int);
  symmetry = reduce_symmetry_in_frame(frame,
				      prim_sym,
				      t_mat_int,
				      cell->lattice,
				      cell->size / primitive->size,
				      symprec);

  /* sym_free_symmetry(f_sym); */
  sym_free_symmetry(prim_sym);
  sym_free_symmetry(conv_sym);

  return symmetry;
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
      conv_sym->trans[i][j] -= tmp_vec[j];
    }
  }
}

static Symmetry * get_primitive_db_symmetry(SPGCONST double t_mat[3][3],
					    const Symmetry *conv_sym,
					    const double symprec)
{
  int i, j, num_op;
  double inv_mat[3][3], tmp_mat[3][3];
  MatINT *r_prim;
  VecDBL *t_prim;
  Symmetry *prim_sym;
  
  r_prim = mat_alloc_MatINT(conv_sym->size);
  t_prim = mat_alloc_VecDBL(conv_sym->size);

  mat_inverse_matrix_d3(inv_mat, t_mat, symprec);

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
    mat_multiply_matrix_vector_d3(t_prim->vec[ num_op ],
				  t_mat,
				  conv_sym->trans[ i ]);
    num_op++;

  pass:
    ;
  }

  prim_sym = sym_alloc_symmetry(num_op);
  for (i = 0; i < num_op; i++) {
    mat_copy_matrix_i3(prim_sym->rot[i], r_prim->mat[i]);
    for (j = 0; j < 3; j++) {
      prim_sym->trans[i][j] = t_prim->vec[i][j] - mat_Nint(t_prim->vec[i][j]);
    }
  }

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
  int i, j, k, l, num_trans, size_sym_orig;
  Symmetry *symmetry, *t_sym;
  double inv_tmat[3][3], tmp_mat[3][3], tmp_rot_d[3][3], tmp_lat_d[3][3], tmp_lat_i[3][3];
  int tmp_rot_i[3][3];
  VecDBL *pure_trans, *lattice_trans;

  mat_cast_matrix_3i_to_3d(tmp_mat, t_mat);
  mat_inverse_matrix_d3(inv_tmat, tmp_mat, symprec);

  /* transformed lattice points */
  lattice_trans = mat_alloc_VecDBL(frame[0]*frame[1]*frame[2]);
  num_trans = 0;
  for (i = 0; i < frame[0]; i++) {
    for (j = 0; j < frame[1]; j++) {
      for (k = 0; k < frame[2]; k++) {
	lattice_trans->vec[num_trans][0] = i;
	lattice_trans->vec[num_trans][1] = j;
	lattice_trans->vec[num_trans][2] = k;

	mat_multiply_matrix_vector_d3(lattice_trans->vec[num_trans],
				      inv_tmat,
				      lattice_trans->vec[num_trans]);
	for (l = 0; l < 3; l++) {
	  /* t' = T^-1*t */
	  lattice_trans->vec[num_trans][l] = \
	    mat_Dmod1(lattice_trans->vec[num_trans][l]); 
	}
	num_trans++;
      }
    }
  }

  /* transformed symmetry operations of primitive cell */
  t_sym = sym_alloc_symmetry(prim_sym->size);
  size_sym_orig = 0;
  for (i = 0; i < prim_sym->size; i++) {
    /* R' = T^-1*R*T */
    mat_multiply_matrix_di3(tmp_mat, inv_tmat, prim_sym->rot[i]);
    mat_multiply_matrix_di3(tmp_rot_d, tmp_mat, t_mat);
    mat_cast_matrix_3d_to_3i(tmp_rot_i, tmp_rot_d);
    mat_multiply_matrix_di3(tmp_lat_i, lattice, tmp_rot_i);
    mat_multiply_matrix_d3(tmp_lat_d, lattice, tmp_rot_d);
    if (mat_check_identity_matrix_d3(tmp_lat_i, tmp_lat_d, symprec)) {
      mat_copy_matrix_i3(t_sym->rot[size_sym_orig], tmp_rot_i);
      /* t' = T^-1*t */
      mat_multiply_matrix_vector_d3(t_sym->trans[size_sym_orig],
				    inv_tmat, prim_sym->trans[i]);
      size_sym_orig++;
    }
  }

  /* reduce lattice points */
  pure_trans = reduce_lattice_points(lattice,
				     lattice_trans,
				     symprec);

  if (! (pure_trans->size == multiplicity)) {
    symmetry = sym_alloc_symmetry(0);
    goto ret;
  }

  /* copy symmetry operations upon lattice points */
  symmetry = sym_alloc_symmetry(pure_trans->size * size_sym_orig);
  for (i = 0; i < pure_trans->size; i++) {
    for (j = 0; j < size_sym_orig; j++) {
      mat_copy_matrix_i3(symmetry->rot[size_sym_orig * i + j],
			 t_sym->rot[j]);
      mat_copy_vector_d3(symmetry->trans[size_sym_orig * i + j],
			 t_sym->trans[j]);
      for (k = 0; k < 3; k++) {
	symmetry->trans[size_sym_orig * i + j][k] += pure_trans->vec[i][k];
	symmetry->trans[size_sym_orig * i + j][k] = \
	  mat_Dmod1(symmetry->trans[size_sym_orig * i + j][k]);
      }
    }
  }
  

 ret:
  mat_free_VecDBL(lattice_trans);
  mat_free_VecDBL(pure_trans);
  sym_free_symmetry(t_sym);

  return symmetry;
}

static VecDBL * reduce_lattice_points(SPGCONST double lattice[3][3],
				      const VecDBL *lattice_trans,
				      const double symprec)
{
  int i, j, is_found, num_pure_trans;
  VecDBL *pure_trans, *t;
  
  num_pure_trans = 0;
  t = mat_alloc_VecDBL(lattice_trans->size);
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

  pure_trans = mat_alloc_VecDBL(num_pure_trans);
  for (i = 0; i < num_pure_trans; i++) {
    mat_copy_vector_d3(pure_trans->vec[i], t->vec[i]);
  }
  mat_free_VecDBL(t);

  return pure_trans;
}
