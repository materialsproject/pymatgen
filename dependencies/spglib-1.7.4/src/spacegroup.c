/* spacegroup.c */
/* Copyright (C) 2010 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cell.h"
#include "hall_symbol.h"
#include "lattice.h"
#include "mathfunc.h"
#include "niggli.h"
#include "pointgroup.h"
#include "primitive.h"
#include "spacegroup.h"
#include "spg_database.h"
#include "symmetry.h"

#include "debug.h"

#define REDUCE_RATE 0.95

static double change_of_basis_monocli[18][3][3] = {{{ 1, 0, 0 },
						    { 0, 1, 0 },
						    { 0, 0, 1 }},
						   {{ 0, 0, 1 },
						    { 0,-1, 0 },
						    { 1, 0, 0 }},
						   {{ 0, 0, 1 },
						    { 1, 0, 0 },
						    { 0, 1, 0 }},
						   {{ 1, 0, 0 },
						    { 0, 0, 1 },
						    { 0,-1, 0 }},
						   {{ 0, 1, 0 },
						    { 0, 0, 1 },
						    { 1, 0, 0 }},
						   {{ 0,-1, 0 },
						    { 1, 0, 0 },
						    { 0, 0, 1 }},
						   {{-1, 0, 1 },
						    { 0, 1, 0 },
						    {-1, 0, 0 }},
						   {{ 1, 0,-1 },
						    { 0,-1, 0 },
						    { 0, 0,-1 }},
						   {{ 0, 1,-1 },
						    { 1, 0, 0 },
						    { 0, 0,-1 }},
						   {{-1,-1, 0 },
						    { 0, 0, 1 },
						    {-1, 0, 0 }},
						   {{ 1,-1, 0 },
						    { 0, 0, 1 },
						    { 0,-1, 0 }},
						   {{ 0, 1, 1 },
						    { 1, 0, 0 },
						    { 0, 1, 0 }},
						   {{ 0, 0,-1 },
						    { 0, 1, 0 },
						    { 1, 0,-1 }},
						   {{-1, 0, 0 },
						    { 0,-1, 0 },
						    {-1, 0, 1 }},
						   {{ 0,-1, 0 },
						    { 1, 0, 0 },
						    { 0,-1, 1 }},
						   {{ 0, 1, 0 },
						    { 0, 0, 1 },
						    { 1, 1, 0 }},
						   {{-1, 0, 0 },
						    { 0, 0, 1 },
						    {-1, 1, 0 }},
						   {{ 0, 0,-1 },
						    { 1, 0, 0 },
						    { 0,-1,-1 }}};

static Centering change_of_centering_monocli[18] = {C_FACE,
						    A_FACE,
						    B_FACE,
						    B_FACE,
						    A_FACE,
						    C_FACE,
						    BASE,
						    BASE,
						    BASE,
						    BASE,
						    BASE,
						    BASE,
						    A_FACE,
						    C_FACE,
						    C_FACE,
						    A_FACE,
						    B_FACE,
						    B_FACE};

static int change_of_unique_axis_monocli[18] =
  {1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0};

static double change_of_basis_ortho[6][3][3] = {{{ 1, 0, 0 },
						 { 0, 1, 0 },
						 { 0, 0, 1 }},
						{{ 0, 0, 1 },
						 { 1, 0, 0 },
						 { 0, 1, 0 }},
						{{ 0, 1, 0 },
						 { 0, 0, 1 },
						 { 1, 0, 0 }},
						{{ 0, 1, 0 },
						 { 1, 0, 0 },
						 { 0, 0,-1 }},
						{{ 1, 0, 0 },
						 { 0, 0, 1 },
						 { 0,-1, 0 }},
						{{ 0, 0, 1 },
						 { 0, 1, 0 },
						 {-1, 0, 0 }}};

static Centering change_of_centering_ortho[6] = {C_FACE,
						 B_FACE,
						 A_FACE,
						 C_FACE,
						 B_FACE,
						 A_FACE};
static int change_of_unique_axis_ortho[6] = {2, 1, 0, 2, 1, 0};

static double hR_to_hP[3][3] = {{ 1, 0, 1 },
				{-1, 1, 1 },
				{ 0,-1, 1 }};
static double change_of_basis_501[3][3] = {{ 0, 0, 1},
					   { 0,-1, 0},
					   { 1, 0, 0}};

static int spacegroup_to_hall_number[230] = {
    1,   2,   3,   6,   9,  18,  21,  30,  39,  57,
   60,  63,  72,  81,  90, 108, 109, 112, 115, 116,
  119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
  155, 161, 164, 170, 173, 176, 182, 185, 191, 197,
  203, 209, 212, 215, 218, 221, 227, 228, 230, 233,
  239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
  290, 292, 298, 304, 310, 313, 316, 322, 334, 335,
  337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
  355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
  368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
  378, 379, 380, 381, 382, 383, 384, 385, 386, 387,
  388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
  398, 399, 400, 401, 402, 404, 406, 407, 408, 410,
  412, 413, 414, 416, 418, 419, 420, 422, 424, 425,
  426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
  440, 441, 442, 443, 444, 446, 447, 448, 449, 450,
  452, 454, 455, 456, 457, 458, 460, 462, 463, 464,
  465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
  475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
  485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
  495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
  507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
  517, 518, 520, 521, 523, 524, 525, 527, 529, 530,
};

static Spacegroup search_spacegroup(SPGCONST Cell * primitive,
				    const int candidates[],
				    const int num_candidates,
				    const double symprec);
static Spacegroup get_spacegroup(const int hall_number,
				 const double origin_shift[3],
				 SPGCONST double conv_lattice[3][3]);
static int iterative_search_hall_number(double origin_shift[3],
					double conv_lattice[3][3],
					const int candidates[],
					const int num_candidates,
					SPGCONST Cell * primitive,
					SPGCONST Symmetry * symmetry,
					const double symprec);
static Symmetry * get_symmetry_settings(double conv_lattice[3][3],
					Pointgroup *pointgroup,
					Centering *centering,
					SPGCONST double primitive_lattice[3][3],
					SPGCONST Symmetry * symmetry,
					const double symprec);
static int search_hall_number(double origin_shift[3],
			      double conv_lattice[3][3],
			      const int candidates[],
			      const int num_candidates,
			      SPGCONST double primitive_lattice[3][3],
			      SPGCONST Symmetry * symmetry,
			      const double symprec);
static int match_hall_symbol_db(double origin_shift[3],
				double lattice[3][3],
				const int hall_number,
				const Pointgroup *pointgroup,
				const Centering centering,
				SPGCONST Symmetry *symmetry,
				const double symprec);
static int match_hall_symbol_db_monocli(double origin_shift[3],
					double lattice[3][3],
					const int hall_number,
					const int num_hall_types,
					const Centering centering,
					SPGCONST Symmetry *symmetry,
					const double symprec);
static int match_hall_symbol_db_ortho(double origin_shift[3],
				      double lattice[3][3],
				      const int hall_number,
				      const Centering centering,
				      SPGCONST Symmetry *symmetry,
				      const int num_free_axes,
				      const double symprec);
static Symmetry * get_conventional_symmetry(SPGCONST double transform_mat[3][3],
					    const Centering centering,
					    const Symmetry *primitive_sym);

Primitive * spa_get_spacegroup(Spacegroup * spacegroup,
			       SPGCONST Cell * cell,
			       const double symprec)
{
  int attempt;
  double tolerance;
  Primitive *primitive;

  tolerance = symprec;

  for (attempt = 0; attempt < 100; attempt++) {
    primitive = prm_get_primitive(cell, tolerance);
    if (primitive->size > 0) {
      *spacegroup = search_spacegroup(primitive->cell,
				      spacegroup_to_hall_number,
				      230,
				      primitive->tolerance);
      if (spacegroup->number > 0) {
	break;
      }
    }
    
    warning_print("spglib: Attempt %d tolerance = %f failed.", attempt, tolerance);
    warning_print(" (line %d, %s).\n", __LINE__, __FILE__);

    tolerance *= REDUCE_RATE;
    prm_free_primitive(primitive);
  }

  if (primitive->size == 0) {
    warning_print("spglib: Space group could not be found ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    primitive = prm_alloc_primitive(0);
    primitive->cell = cel_alloc_cell(0);
    primitive->pure_trans = mat_alloc_VecDBL(0);
  }

  return primitive;
}

Spacegroup spa_get_spacegroup_with_hall_number(SPGCONST Cell * primitive,
					       const int hall_number,
					       const double symprec)
{
  int num_candidates;
  int candidate[1];
  Spacegroup spacegroup;
  
  if (hall_number < 1 || hall_number > 530 || primitive->size < 1) {
    goto err;
  }
    
  num_candidates = 1;
  candidate[0] = hall_number;
  spacegroup = search_spacegroup(primitive,
				 candidate,
				 num_candidates,
				 symprec);
  if (spacegroup.number > 0) {
    goto ret;
  }

 err:
  spacegroup.number = 0;
  warning_print("spglib: Space group with the input setting could not be found ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);

 ret:
  return spacegroup;
}

static Spacegroup search_spacegroup(SPGCONST Cell * primitive,
				    const int candidates[],
				    const int num_candidates,
				    const double symprec)
{
  int hall_number;
  double conv_lattice[3][3];
  double origin_shift[3];
  Symmetry *symmetry;

  hall_number = 0;
  symmetry = sym_get_operation(primitive, symprec);
  if (symmetry->size > 0) {
    hall_number = iterative_search_hall_number(origin_shift,
					       conv_lattice,
					       candidates,
					       num_candidates,
					       primitive,
					       symmetry,
					       symprec);
  }
  sym_free_symmetry(symmetry);

  return get_spacegroup(hall_number, origin_shift, conv_lattice);
}

static Spacegroup get_spacegroup(const int hall_number,
				 const double origin_shift[3],
				 SPGCONST double conv_lattice[3][3])
{
  Spacegroup spacegroup;
  SpacegroupType spacegroup_type;
  
  if (hall_number == 0) {
    spacegroup.number = 0;
    warning_print("spglib: Find Hall symbol failed.\n");
    warning_print("spglib: Space group could not be found ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
    goto ret;
  }

  spacegroup_type = spgdb_get_spacegroup_type(hall_number);

  if (spacegroup_type.number > 0) {
    mat_copy_matrix_d3(spacegroup.bravais_lattice, conv_lattice);
    mat_copy_vector_d3(spacegroup.origin_shift, origin_shift);
    spacegroup.number = spacegroup_type.number;
    spacegroup.hall_number = hall_number;
    spacegroup.pointgroup_number = spacegroup_type.pointgroup_number;
    strcpy(spacegroup.schoenflies,
	   spacegroup_type.schoenflies);
    strcpy(spacegroup.hall_symbol,
	   spacegroup_type.hall_symbol);
    strcpy(spacegroup.international,
	   spacegroup_type.international);
    strcpy(spacegroup.international_long,
	   spacegroup_type.international_full);
    strcpy(spacegroup.international_short,
	   spacegroup_type.international_short);
    strcpy(spacegroup.setting,
	   spacegroup_type.setting);
  } else {
    spacegroup.number = 0;
    warning_print("spglib: Space group could not be found ");
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  }

 ret:
  /* spacegroup.number = 0 when space group was not found. */
  return spacegroup;
}

static int iterative_search_hall_number(double origin_shift[3],
					double conv_lattice[3][3],
					const int candidates[],
					const int num_candidates,
					SPGCONST Cell * primitive,
					SPGCONST Symmetry * symmetry,
					const double symprec)
{
  int i, attempt, hall_number=0;
  double tolerance;
  Symmetry * sym_reduced;

  debug_print("iterative_search_hall_number:\n");

  sym_reduced = sym_alloc_symmetry(symmetry->size);
  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(sym_reduced->rot[i], symmetry->rot[i]);
    mat_copy_vector_d3(sym_reduced->trans[i], symmetry->trans[i]);
  }
  
  tolerance = symprec;
  for (attempt = 0; attempt < 100; attempt++) {
    hall_number = search_hall_number(origin_shift,
				     conv_lattice,
				     candidates,
				     num_candidates,
				     primitive->lattice,
				     sym_reduced,
				     symprec);
    if (hall_number > 0) {
      break;
    }

    warning_print("spglib: Attempt %d tolerance = %f failed", attempt, tolerance);
    warning_print("(line %d, %s).\n", __LINE__, __FILE__);

    sym_free_symmetry(sym_reduced);
    tolerance *= REDUCE_RATE;
    sym_reduced = sym_reduce_operation(primitive, symmetry, tolerance);
  }

#ifdef SPGWARNING
  if (hall_number == 0) {
    warning_print("spglib: Iterative attempt with sym_reduce_operation to find Hall symbol failed.\n");
  }
#endif

  sym_free_symmetry(sym_reduced);
  
  return hall_number;
}

static int search_hall_number(double origin_shift[3],
			      double conv_lattice[3][3],
			      const int candidates[],
			      const int num_candidates,
			      SPGCONST double primitive_lattice[3][3],
			      SPGCONST Symmetry * symmetry,
			      const double symprec)
{
  int i, hall_number;
  Centering centering;
  Pointgroup pointgroup;
  Symmetry * conv_symmetry;

  debug_print("get_hall_number:\n");

  conv_symmetry = get_symmetry_settings(conv_lattice,
					&pointgroup,
					&centering,
					primitive_lattice,
					symmetry,
					symprec);
  if (conv_symmetry->size == 0) {
    hall_number = 0;
    goto ret;
  }

  for (i = 0; i < num_candidates; i++) {
    hall_number = candidates[i];
    if (match_hall_symbol_db(origin_shift,
			     conv_lattice,
			     hall_number,
			     &pointgroup,
			     centering,
			     conv_symmetry,
			     symprec)) {goto ret;}
  }
  
  hall_number = 0;

 ret:
  sym_free_symmetry(conv_symmetry);
  return hall_number;
}

static Symmetry * get_symmetry_settings(double conv_lattice[3][3],
					Pointgroup *pointgroup,
					Centering *centering,
					SPGCONST double primitive_lattice[3][3],
					SPGCONST Symmetry * symmetry,
					const double symprec)
{
  int i, j;
  int int_transform_mat[3][3];
  double correction_mat[3][3], transform_mat[3][3], inv_lattice[3][3], smallest_lattice[3][3];
  double niggli_cell[9];
  Symmetry * conv_symmetry;
  
  *pointgroup = ptg_get_transformation_matrix(int_transform_mat,
					      symmetry->rot,
					      symmetry->size);

  if (pointgroup->number < 1) {
    conv_symmetry = sym_alloc_symmetry(0);
    *centering = NO_CENTER;
    goto ret;
  }

  mat_multiply_matrix_di3(conv_lattice,
			  primitive_lattice,
			  int_transform_mat);

  /* Triclinic: Niggli cell reduction */
  if (pointgroup->laue == LAUE1) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
  	niggli_cell[i * 3 + j] = conv_lattice[i][j];
      }
    }
    niggli_reduce(niggli_cell, symprec * symprec);
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
  	smallest_lattice[i][j] = niggli_cell[i * 3 + j];
      }
    }
    if (mat_get_determinant_d3(smallest_lattice) < 0) {
      for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	  smallest_lattice[i][j] = -smallest_lattice[i][j];
	}
      }
    }
    mat_inverse_matrix_d3(inv_lattice, primitive_lattice, 0);
    mat_multiply_matrix_d3(transform_mat, inv_lattice, smallest_lattice);
    mat_cast_matrix_3d_to_3i(int_transform_mat, transform_mat);
  }

  /* Monoclinic: choose shortest a, c lattice vectors (|a| < |c|) */
  if (pointgroup->laue == LAUE2M) {
    if (lat_smallest_lattice_vector_2D(smallest_lattice,
				       conv_lattice,
				       1, /* unique axis of b */
				       symprec)) {
      mat_inverse_matrix_d3(inv_lattice, primitive_lattice, 0);
      mat_multiply_matrix_d3(transform_mat, inv_lattice, smallest_lattice);
      mat_cast_matrix_3d_to_3i(int_transform_mat, transform_mat);
    }
  }

  *centering = lat_get_centering(correction_mat,
				 int_transform_mat,
				 pointgroup->laue);
  
  mat_multiply_matrix_id3(transform_mat, int_transform_mat, correction_mat);
  mat_multiply_matrix_d3(conv_lattice, primitive_lattice, transform_mat);

  if (*centering == R_CENTER) {
    /* hP for rhombohedral */
    conv_symmetry = get_conventional_symmetry(transform_mat,
					      NO_CENTER,
					      symmetry);
  } else {
    conv_symmetry = get_conventional_symmetry(transform_mat,
					      *centering,
					      symmetry);
  }

 ret:
  return conv_symmetry;
}

static int match_hall_symbol_db(double origin_shift[3],
				double lattice[3][3],
				const int hall_number,
				const Pointgroup *pointgroup,
				const Centering centering,
				SPGCONST Symmetry *symmetry,
				const double symprec)
{
  int is_found, num_hall_types;
  SpacegroupType spacegroup_type;
  Symmetry * changed_symmetry;
  double changed_lattice[3][3], inv_lattice[3][3], transform_mat[3][3];
  
  spacegroup_type = spgdb_get_spacegroup_type(hall_number);
  num_hall_types = (spacegroup_to_hall_number[spacegroup_type.number] -
		    spacegroup_to_hall_number[spacegroup_type.number - 1]);

  if (pointgroup->number == spacegroup_type.pointgroup_number) {
    switch (pointgroup->holohedry) {
    case MONOCLI:
      if (match_hall_symbol_db_monocli(origin_shift,
				       lattice,
				       hall_number,
				       num_hall_types,
				       centering,
				       symmetry,
				       symprec)) {return 1;}
      break;
      
    case ORTHO:
      if (spacegroup_type.number == 48 ||
	  spacegroup_type.number == 50 ||
	  spacegroup_type.number == 59 ||
	  spacegroup_type.number == 68 ||
	  spacegroup_type.number == 70) { /* uncount origin shift */
	num_hall_types /= 2;
      }
      
      if (num_hall_types == 1) {
	if (match_hall_symbol_db_ortho(origin_shift,
				       lattice,
				       hall_number,
				       centering,
				       symmetry,
				       6,
				       symprec)) {return 1;}
	break;
      }

      if (num_hall_types == 2) {
	if (match_hall_symbol_db_ortho(origin_shift,
				       lattice,
				       hall_number,
				       centering,
				       symmetry,
				       3,
				       symprec)) {return 1;}
	break;
      }

      if (num_hall_types == 3) {
	mat_copy_matrix_d3(changed_lattice, lattice);
	if (! match_hall_symbol_db_ortho
	    (origin_shift,
	     changed_lattice,
	     spacegroup_to_hall_number[spacegroup_type.number - 1],
	     centering,
	     symmetry,
	     0,
	     symprec)) {break;}
	mat_inverse_matrix_d3(inv_lattice, lattice, 0);
	mat_multiply_matrix_d3(transform_mat, inv_lattice, changed_lattice);
	changed_symmetry = get_conventional_symmetry(transform_mat,
						     NO_CENTER,
						     symmetry);
	is_found = match_hall_symbol_db_ortho(origin_shift,
					      changed_lattice,
					      hall_number,
					      centering,
					      changed_symmetry,
					      2,
					      symprec);
	sym_free_symmetry(changed_symmetry);	
	if (is_found) {
	  mat_copy_matrix_d3(lattice, changed_lattice);
	  return 1;
	}
	break;
      }

      if (num_hall_types == 6) {
	if (match_hall_symbol_db_ortho(origin_shift,
				       lattice,
				       hall_number,
				       centering,
				       symmetry,
				       1,
				       symprec)) {return 1;}
	break;
      }

      break;

    case CUBIC:
      if (hal_match_hall_symbol_db(origin_shift,
				   lattice,
				   hall_number,
				   centering,
				   symmetry,
				   symprec)) {return 1;}
      
      if (hall_number == 501) { /* Try another basis for No.205 */
	mat_multiply_matrix_d3(changed_lattice,
			       lattice,
			       change_of_basis_501);
	changed_symmetry =
	  get_conventional_symmetry(change_of_basis_501,
				    NO_CENTER,
				    symmetry);
	is_found = hal_match_hall_symbol_db(origin_shift,
					    changed_lattice,
					    hall_number,
					    NO_CENTER,
					    changed_symmetry,
					    symprec);
	sym_free_symmetry(changed_symmetry);
	if (is_found) {
	  mat_copy_matrix_d3(lattice, changed_lattice);
	  return 1;
	}
      }
      break;
      
    case TRIGO:
      if (centering == R_CENTER) {
	if (hall_number == 433 ||
	    hall_number == 436 ||
	    hall_number == 444 ||
	    hall_number == 450 ||
	    hall_number == 452 ||
	    hall_number == 458 ||
	    hall_number == 460) {
	  mat_multiply_matrix_d3(changed_lattice,
				 lattice,
				 hR_to_hP);
	  changed_symmetry =
	    get_conventional_symmetry(hR_to_hP,
				      R_CENTER,
				      symmetry);
	  is_found = hal_match_hall_symbol_db(origin_shift,
					      changed_lattice,
					      hall_number,
					      centering,
					      changed_symmetry,
					      symprec);
	  sym_free_symmetry(changed_symmetry);
	  if (is_found) {
	    mat_copy_matrix_d3(lattice, changed_lattice);
	    return 1;
	  }
	}
      }
      /* Do not break for other trigonal cases */
    default: /* HEXA, TETRA, TRICLI and rest of TRIGO */
      if (hal_match_hall_symbol_db(origin_shift,
				   lattice,
				   hall_number,
				   centering,
				   symmetry,
				   symprec)) {
	return 1;
      }
      break;
    }
  }
  return 0;
}

static int match_hall_symbol_db_monocli(double origin_shift[3],
					double lattice[3][3],
					const int hall_number,
					const int num_hall_types,
					const Centering centering,
					SPGCONST Symmetry *symmetry,
					const double symprec)
{
  int i, j, k, l, is_found;
  double vec[3], norms[3];
  Centering changed_centering;
  Symmetry * changed_symmetry;
  double changed_lattice[3][3];
  
  for (i = 0; i < 18; i++) {
    if (centering == C_FACE) {
      changed_centering = change_of_centering_monocli[i];
    } else { /* suppose NO_CENTER */
      changed_centering = centering;
    }

    mat_multiply_matrix_d3(changed_lattice,
			   lattice,
			   change_of_basis_monocli[i]);

    /* Choose |a| < |b| < |c| if there are freedom. */
    if (num_hall_types == 3) {
      l = 0;
      for (j = 0; j < 3; j++) {
	if (j == change_of_unique_axis_monocli[i]) {continue;}
	for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
	norms[l] = mat_norm_squared_d3(vec);
	l++;
      }
      if (norms[0] > norms[1]) {continue;}
    }

    changed_symmetry =
      get_conventional_symmetry(change_of_basis_monocli[i],
				NO_CENTER,
				symmetry);
    is_found = hal_match_hall_symbol_db(origin_shift,
					changed_lattice,
					hall_number,
					changed_centering,
					changed_symmetry,
					symprec);
    sym_free_symmetry(changed_symmetry);
    if (is_found) {
      mat_copy_matrix_d3(lattice, changed_lattice);
      return 1;
    }
  }
  return 0;
}

static int match_hall_symbol_db_ortho(double origin_shift[3],
				      double lattice[3][3],
				      const int hall_number,
				      const Centering centering,
				      SPGCONST Symmetry *symmetry,
				      const int num_free_axes,
				      const double symprec)
{
  int i, j, k, l, is_found;
  double vec[3], norms[3];
  Centering changed_centering;
  Symmetry * changed_symmetry;
  double changed_lattice[3][3];
  
  for (i = 0; i < 6; i++) {
    if (centering == C_FACE) {
      changed_centering = change_of_centering_ortho[i];
    } else {
      changed_centering = centering;
    }
    
    mat_multiply_matrix_d3(changed_lattice,
			   lattice,
			   change_of_basis_ortho[i]);

    if (num_free_axes == 2) {
      l = 0;
      for (j = 0; j < 3; j++) {
    	if (j == change_of_unique_axis_ortho[i]) {continue;}
    	for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
    	norms[l] = mat_norm_squared_d3(vec);
    	l++;
      }
      if (norms[0] > norms[1]) {continue;}
    }

    if (num_free_axes == 3) {
      for (j = 0; j < 3; j++) {
    	for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
    	norms[j] = mat_norm_squared_d3(vec);
      }
      if (norms[0] > norms[1] || norms[0] > norms[2]) {continue;}
    }

    if (num_free_axes == 6) {
      for (j = 0; j < 3; j++) {
    	for (k = 0; k < 3; k++) {vec[k] = changed_lattice[k][j];}
    	norms[j] = mat_norm_squared_d3(vec);
      }
      if (norms[0] > norms[1] || norms[1] > norms[2]) {continue;}
    }
    
    changed_symmetry = get_conventional_symmetry(change_of_basis_ortho[i],
						 NO_CENTER,
						 symmetry);

    is_found = hal_match_hall_symbol_db(origin_shift,
					changed_lattice,
					hall_number,
					changed_centering,
					changed_symmetry,
					symprec);
    sym_free_symmetry(changed_symmetry);
    if (is_found) {
      mat_copy_matrix_d3(lattice, changed_lattice);
      return 1;
    }
  }
  return 0;
}

static Symmetry * get_conventional_symmetry(SPGCONST double transform_mat[3][3],
					    const Centering centering,
					    const Symmetry *primitive_sym)
{
  int i, j, k, multi, size;
  double tmp_trans;
  double tmp_matrix_d3[3][3], shift[4][3];
  double symmetry_rot_d3[3][3], primitive_sym_rot_d3[3][3];
  Symmetry *symmetry;

  size = primitive_sym->size;

  switch (centering) {
  case FACE:
    symmetry = sym_alloc_symmetry(size * 4);
    break;
  case R_CENTER:
    symmetry = sym_alloc_symmetry(size * 3);
    break;
  case BODY:
  case A_FACE:
  case B_FACE:
  case C_FACE:
    symmetry = sym_alloc_symmetry(size * 2);
    break;
  default:
    symmetry = sym_alloc_symmetry(size);
    break;
  }

  for (i = 0; i < size; i++) {
    mat_cast_matrix_3i_to_3d(primitive_sym_rot_d3, primitive_sym->rot[i]);

    /* C*S*C^-1: recover conventional cell symmetry operation */
    mat_get_similar_matrix_d3(symmetry_rot_d3,
			      primitive_sym_rot_d3,
			      transform_mat,
			      0);
    mat_cast_matrix_3d_to_3i(symmetry->rot[i], symmetry_rot_d3);

    /* translation in conventional cell: C = B^-1*P */
    mat_inverse_matrix_d3(tmp_matrix_d3,
			  transform_mat,
			  0);
    mat_multiply_matrix_vector_d3(symmetry->trans[i],
				  tmp_matrix_d3,
				  primitive_sym->trans[i]);
  }

  multi = 1;

  if (centering != NO_CENTER) {
    if (centering != FACE && centering != R_CENTER) {
      for (i = 0; i < 3; i++) {	shift[0][i] = 0.5; } /* BASE */
      if (centering == A_FACE) { shift[0][0] = 0; }
      if (centering == B_FACE) { shift[0][1] = 0; }
      if (centering == C_FACE) { shift[0][2] = 0; }
      multi = 2;
    }

    if (centering == R_CENTER) {
      shift[0][0] = 2. / 3;
      shift[0][1] = 1. / 3;
      shift[0][2] = 1. / 3;
      shift[1][0] = 1. / 3;
      shift[1][1] = 2. / 3;
      shift[1][2] = 2. / 3;
      multi = 3;
    }
    
    if (centering == FACE) {
      shift[0][0] = 0;
      shift[0][1] = 0.5;
      shift[0][2] = 0.5;
      shift[1][0] = 0.5;
      shift[1][1] = 0;
      shift[1][2] = 0.5;
      shift[2][0] = 0.5;
      shift[2][1] = 0.5;
      shift[2][2] = 0;
      multi = 4;
    }

    for (i = 0; i < multi - 1; i++) {
      for (j = 0; j < size; j++) {
	mat_copy_matrix_i3(symmetry->rot[(i+1) * size + j],
			   symmetry->rot[j]);
	for (k = 0; k < 3; k++) {
	  tmp_trans = symmetry->trans[j][k] + shift[i][k];
	  symmetry->trans[(i+1) * size + j][k] = tmp_trans;
	}
      }
    }
  }


  /* Reduce translations into -0 < trans < 1.0 */
  for (i = 0; i < multi; i++) {
    for (j = 0; j < size; j++) {
      for (k = 0; k < 3; k++) {
  	tmp_trans = symmetry->trans[i * size + j][k];
  	tmp_trans -= mat_Nint(tmp_trans);
  	if (tmp_trans < 0) {
  	  tmp_trans += 1.0;
  	}
  	symmetry->trans[i * size + j][k] = tmp_trans;
      }
    }
  }

  return symmetry;
}

