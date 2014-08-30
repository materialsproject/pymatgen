/* lattice.c */
/* Copyright (C) 2010 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include "lattice.h"
#include "mathfunc.h"
#include "debug.h"

#define INT_PREC 0.1

static double identity[3][3] = {{ 1, 0, 0 },
				{ 0, 1, 0 },
				{ 0, 0, 1 }};
static double monocli_i2c[3][3] = {{ 1, 0, 0 },
				   { 0, 1, 0 },
				   { 1, 0,-1 }};
static double monocli_a2c[3][3] = {{ 0, 0, 1 },
				   { 0,-1, 0 },
				   { 1, 0, 0 }};
#ifdef DEBUG
static double tetra_f2i_c2p[3][3] = {{ 0.5,-0.5, 0.0 },
				     { 0.5, 0.5, 0.0 },
				     { 0.0, 0.0, 1.0 }};
static double hexa_h2p[3][3] = {{ 2./3,-1./3, 0.0 },
				{ 1./3, 1./3, 0.0 },
				{  0.0,  0.0, 1.0 }};
#endif
static double rhombo_obverse[3][3] = {{ 2./3,-1./3,-1./3 },
				      { 1./3, 1./3,-2./3 },
				      { 1./3, 1./3, 1./3 }};
static double rhomb_reverse[3][3] = {{ 1./3,-2./3, 1./3 },
				     { 2./3,-1./3,-1./3 },
				     { 1./3, 1./3, 1./3 }};
static double a2c[3][3] = {{ 0, 0, 1 },
			   { 1, 0, 0 },
			   { 0, 1, 0 }};
static double b2c[3][3] = {{ 0, 1, 0 },
			   { 0, 0, 1 },
			   { 1, 0, 0 }};

static Centering get_centering(double correction_mat[3][3],
			       SPGCONST int transform_mat[3][3],
			       const Laue laue);
static int get_Delaunay_reduction(double red_lattice[3][3], 
				  SPGCONST double lattice[3][3],
				  SPGCONST double symprec);
static Centering get_base_center(SPGCONST int transform_mat[3][3]);
static int get_Delaunay_reduction_basis(double basis[4][3],
					const double symprec);
static void get_Delaunay_shortest_vectors(double basis[4][3],
					  const double symprec);
static void get_exteneded_basis(double basis[4][3],
				SPGCONST double lattice[3][3]);

int lat_smallest_lattice_vector(double min_lattice[3][3],
				SPGCONST double lattice[3][3],
				const double symprec)
{
  debug_print("lat_smallest_lattice_vector:\n");
  return get_Delaunay_reduction(min_lattice, lattice, symprec);
}

Centering lat_get_centering(double correction_mat[3][3],
			    SPGCONST int transform_mat[3][3],
			    const Laue laue)
{
  return get_centering(correction_mat,
		       transform_mat,
		       laue);
}

static Centering get_centering(double correction_mat[3][3],
			       SPGCONST int transform_mat[3][3],
			       const Laue laue)
{
  int det;
  double trans_corr_mat[3][3];
  Centering centering;

  mat_copy_matrix_d3(correction_mat, identity);
  det = abs(mat_get_determinant_i3(transform_mat));
  debug_print("laue class: %d\n", laue);
  debug_print("multiplicity: %d\n", det);

  if (det == 1) { centering = NO_CENTER; }
  if (det == 2) { centering = get_base_center(transform_mat);
    if (centering == A_FACE) {
      if (laue == LAUE2M) {
	debug_print("Monocli A to C\n");
	mat_copy_matrix_d3(correction_mat, monocli_a2c);
      } else {
	mat_copy_matrix_d3(correction_mat, a2c);
      }
      centering = C_FACE;
    }
    if (centering == B_FACE) {
      mat_copy_matrix_d3(correction_mat, b2c);
      centering = C_FACE;
    }
    if (laue == LAUE2M && centering == BODY) {
      debug_print("Monocli I to C\n");
      mat_copy_matrix_d3(correction_mat, monocli_i2c);
      centering = C_FACE;
    }
  }
  if (det == 3) {
    centering = NO_CENTER;
    mat_multiply_matrix_id3(trans_corr_mat,
			    transform_mat, rhombo_obverse);
    if (mat_is_int_matrix(trans_corr_mat, INT_PREC)) {
      mat_copy_matrix_d3(correction_mat, rhombo_obverse);
      debug_print("R-center observe setting\n");
      debug_print_matrix_d3(trans_corr_mat);
    }
    mat_multiply_matrix_id3(trans_corr_mat,
			    transform_mat, rhomb_reverse);
    if (mat_is_int_matrix(trans_corr_mat, INT_PREC)) {
      mat_copy_matrix_d3(correction_mat, rhomb_reverse);
      debug_print("R-center reverse setting\n");
      debug_print_matrix_d3(trans_corr_mat);
    }
  }
  if (det == 4) { centering = FACE; }

  return centering;
}

static Centering get_base_center(SPGCONST int transform_mat[3][3])
{
  int i;
  Centering centering = NO_CENTER;

  debug_print("lat_get_base_center\n");

  /* C center */
  for (i = 0; i < 3; i++) {
    if (transform_mat[i][0] == 0 &&
	transform_mat[i][1] == 0 &&
	abs(transform_mat[i][2]) == 1) {
      centering = C_FACE;
      goto end;
    }
  }

  /* A center */
  for (i = 0; i < 3; i++) {
    if (abs(transform_mat[i][0]) == 1 && 
	transform_mat[i][1] == 0 &&
	transform_mat[i][2] == 0) {
      centering = A_FACE;
      goto end;
    }
  }

  /* B center */
  for (i = 0; i < 3; i++) {
    if (transform_mat[i][0] == 0 &&
	abs(transform_mat[i][1]) == 1 && 
	transform_mat[i][2] == 0) {
      centering = B_FACE;
      goto end;
    }
  }

  /* body center */
  if (abs(transform_mat[0][0]) +
      abs(transform_mat[0][1]) + 
      abs(transform_mat[0][2]) == 2) {
    centering = BODY;
    goto end;
  }

  /* This should not happen. */
  warning_print("spglib: No centring was found (line %d, %s).\n", __LINE__, __FILE__);
  return NO_CENTER;

 end:
  debug_print("centering: %d\n", centering);
  return centering;
}

/* Delaunay reduction */
/* Reference can be found in International table A. */
static int get_Delaunay_reduction(double red_lattice[3][3], 
				  SPGCONST double lattice[3][3],
				  const double symprec)
{
  int i, j;
  double volume, sum, red_sum;
  double basis[4][3];

  get_exteneded_basis(basis, lattice);
  
  sum = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      sum += basis[i][j] * basis[i][j];
    }
  }

  while (1) {
    if (get_Delaunay_reduction_basis(basis, symprec)) {
      break;
    }
  }

  get_Delaunay_shortest_vectors(basis, symprec);

  red_sum = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      red_sum += basis[i][j] * basis[i][j];
    }
  }

  if (sum - red_sum > symprec * symprec) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
	red_lattice[i][j] = basis[j][i];
      }
    }
  } else {
    mat_copy_matrix_d3(red_lattice, lattice);
  }

  volume = mat_get_determinant_d3(red_lattice);
  if (mat_Dabs(volume) < symprec) {
    warning_print("spglib: Minimum lattice has no volume (line %d, %s).\n", __LINE__, __FILE__);
    goto err;
  }

  if (volume  < 0) {
    /* Flip axes */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
	red_lattice[i][j] = -red_lattice[i][j];
      }
    }
  }

  return 1;

 err:
  return 0;
}

static void get_Delaunay_shortest_vectors(double basis[4][3],
					  const double symprec)
{
  int i, j;
  double tmpmat[3][3], b[7][3], tmpvec[3];
  
  /* Search in the set {b1, b2, b3, b4, b1+b2, b2+b3, b3+b1} */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      b[i][j] = basis[i][j];
    }
  }
  
  for (i = 0; i < 3; i++) {
    b[4][i] = basis[0][i] + basis[1][i];
  }
  for (i = 0; i < 3; i++) {
    b[5][i] = basis[1][i] + basis[2][i];
  }
  for (i = 0; i < 3; i++) {
    b[6][i] = basis[2][i] + basis[0][i];
  }
  
  /* Bubble sort */
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      if (mat_norm_squared_d3(b[j]) > mat_norm_squared_d3(b[j+1])) {
	mat_copy_vector_d3(tmpvec, b[j]);
	mat_copy_vector_d3(b[j], b[j+1]);
	mat_copy_vector_d3(b[j+1], tmpvec);
      }
    }
  }

  for (i = 2; i < 7; i++) {
    for (j = 0; j < 3; j++) {
      tmpmat[j][0] = b[0][j];
      tmpmat[j][1] = b[1][j];
      tmpmat[j][2] = b[i][j];
    }
    if (mat_Dabs(mat_get_determinant_d3(tmpmat)) > symprec) {
      for (j = 0; j < 3; j++) {
	basis[0][j] = b[0][j];
	basis[1][j] = b[1][j];
	basis[2][j] = b[i][j];
      }
      break;
    }
  }
}

static int get_Delaunay_reduction_basis(double basis[4][3],
					const double symprec)
{
  int i, j, k, l;
  double dot_product;

  for (i = 0; i < 4; i++) {
    for (j = i+1; j < 4; j++) {
      dot_product = 0.0;
      for (k = 0; k < 3; k++) {
	dot_product += basis[i][k] * basis[j][k];
      }
      if (dot_product > symprec) {
	for (k = 0; k < 4; k++) {
	  if (! (k == i || k == j)) {
	    for (l = 0; l < 3; l++) {
	      basis[k][l] += basis[i][l];
	    }
	  }
	}
	for (k = 0; k < 3; k++) {
	  basis[i][k] = -basis[i][k];
	}
	return 0;
      }
    }
  }

  return 1;
}

static void get_exteneded_basis(double basis[4][3],
				SPGCONST double lattice[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      basis[i][j] = lattice[j][i];
    }
  }

  for (i = 0; i < 3; i++) {
    basis[3][i] = -lattice[i][0] -lattice[i][1] -lattice[i][2];
  }
}

