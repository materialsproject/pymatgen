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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "lattice.h"
#include "pointgroup.h"
#include "symmetry.h"
#include "mathfunc.h"

#include "debug.h"

#define NUM_ROT_AXES 73

typedef struct {
  int table[10];
  char symbol[6];
  char schoenflies[4];
  Holohedry holohedry;
  Laue laue;
} PointgroupType;

static PointgroupType pointgroup_data[33] = {
  { /* 0 */
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    "     ",
    "   ",
    HOLOHEDRY_NONE,
    LAUE_NONE,
  },
  { /* 1 */
    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
    "1    ",
    "C1 ",
    TRICLI,
    LAUE1,
  },
  { /* 2 */
    {0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
    "-1   ",
    "Ci ",
    TRICLI,
    LAUE1,
  },
  { /* 3 */
    {0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
    "2    ",
    "C2 ",
    MONOCLI,
    LAUE2M,
  },
  { /* 4 */
    {0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
    "m    ",
    "Cs ",
    MONOCLI,
    LAUE2M,
  },
  { /* 5 */
    {0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
    "2/m  ",
    "C2h",
    MONOCLI,
    LAUE2M,
  },
  { /* 6 */
    {0, 0, 0, 0, 0, 1, 3, 0, 0, 0},
    "222  ",
    "D2 ",
    ORTHO,
    LAUEMMM,
  },
  { /* 7 */
    {0, 0, 0, 2, 0, 1, 1, 0, 0, 0},
    "mm2  ",
    "C2v",
    ORTHO,
    LAUEMMM,
  },
  { /* 8 */
    {0, 0, 0, 3, 1, 1, 3, 0, 0, 0},
    "mmm  ",
    "D2h",
    ORTHO,
    LAUEMMM,
  },
  { /* 9 */
    {0, 0, 0, 0, 0, 1, 1, 0, 2, 0},
    "4    ",
    "C4 ",
    TETRA,
    LAUE4M,
  },
  { /* 10 */
    {0, 2, 0, 0, 0, 1, 1, 0, 0, 0},
    "-4   ",
    "S4 ",
    TETRA,
    LAUE4M,
  },
  { /* 11 */
    {0, 2, 0, 1, 1, 1, 1, 0, 2, 0},
    "4/m  ",
    "C4h",
    TETRA,
    LAUE4M,
  },
  { /* 12 */
    {0, 0, 0, 0, 0, 1, 5, 0, 2, 0},
    "422  ",
    "D4 ",
    TETRA,
    LAUE4MMM,
  },
  { /* 13 */
    {0, 0, 0, 4, 0, 1, 1, 0, 2, 0},
    "4mm  ",
    "C4v",
    TETRA,
    LAUE4MMM,
  },
  { /* 14 */
    {0, 2, 0, 2, 0, 1, 3, 0, 0, 0},
    "-42m ",
    "D2d",
    TETRA,
    LAUE4MMM,
  },
  { /* 15 */
    {0, 2, 0, 5, 1, 1, 5, 0, 2, 0},
    "4/mmm",
    "D4h",
    TETRA,
    LAUE4MMM,
  },
  { /* 16 */
    {0, 0, 0, 0, 0, 1, 0, 2, 0, 0},
    "3    ",
    "C3 ",
    TRIGO,
    LAUE3,
  },
  { /* 17 */
    {0, 0, 2, 0, 1, 1, 0, 2, 0, 0},
    "-3   ",
    "C3i",
    TRIGO,
    LAUE3,
  },
  { /* 18 */
    {0, 0, 0, 0, 0, 1, 3, 2, 0, 0},
    "32   ",
    "D3 ",
    TRIGO,
    LAUE3M,
  },
  { /* 19 */
    {0, 0, 0, 3, 0, 1, 0, 2, 0, 0},
    "3m   ",
    "C3v",
    TRIGO,
    LAUE3M,
  },
  { /* 20 */
    {0, 0, 2, 3, 1, 1, 3, 2, 0, 0},
    "-3m  ",
    "D3d",
    TRIGO,
    LAUE3M,
  },
  { /* 21 */
    {0, 0, 0, 0, 0, 1, 1, 2, 0, 2},
    "6    ",
    "C6 ",
    HEXA,
    LAUE6M,
  },
  { /* 22 */
    {2, 0, 0, 1, 0, 1, 0, 2, 0, 0},
    "-6   ",
    "C3h",
    HEXA,
    LAUE6M,
  },
  { /* 23 */
    {2, 0, 2, 1, 1, 1, 1, 2, 0, 2},
    "6/m  ",
    "C6h",
    HEXA,
    LAUE6M,
  },
  { /* 24 */
    {0, 0, 0, 0, 0, 1, 7, 2, 0, 2},
    "622  ",
    "D6 ",
    HEXA,
    LAUE6MMM,
  },
  { /* 25 */
    {0, 0, 0, 6, 0, 1, 1, 2, 0, 2},
    "6mm  ",
    "C6v",
    HEXA,
    LAUE6MMM,
  },
  { /* 26 */
    {2, 0, 0, 4, 0, 1, 3, 2, 0, 0},
    "-6m2 ",
    "D3h",
    HEXA,
    LAUE6MMM,
  },
  { /* 27 */
    {2, 0, 2, 7, 1, 1, 7, 2, 0, 2},
    "6/mmm",
    "D6h",
    HEXA,
    LAUE6MMM,
  },
  { /* 28 */
    {0, 0, 0, 0, 0, 1, 3, 8, 0, 0},
    "23   ",
    "T  ",
    CUBIC,
    LAUEM3,
  },
  { /* 29 */
    {0, 0, 8, 3, 1, 1, 3, 8, 0, 0},
    "m-3  ",
    "Th ",
    CUBIC,
    LAUEM3,
  },
  { /* 30 */
    {0, 0, 0, 0, 0, 1, 9, 8, 6, 0},
    "432  ",
    "O  ",
    CUBIC,
    LAUEM3M,
  },
  { /* 31 */
    {0, 6, 0, 6, 0, 1, 3, 8, 0, 0},
    "-43m ",
    "Td ",
    CUBIC,
    LAUEM3M,
  },
  { /* 32 */
    {0, 6, 8, 9, 1, 1, 9, 8, 6, 0},
    "m-3m ",
    "Oh ",
    CUBIC,
    LAUEM3M,
  }
};

static int identity[3][3] = {
  { 1, 0, 0},
  { 0, 1, 0},
  { 0, 0, 1},
};

static int inversion[3][3] = {
  {-1, 0, 0},
  { 0,-1, 0},
  { 0, 0,-1},
};

static int rot_axes[][3] = {
  { 1, 0, 0},
  { 0, 1, 0},
  { 0, 0, 1},
  { 0, 1, 1},
  { 1, 0, 1},
  { 1, 1, 0},
  { 0, 1,-1},
  {-1, 0, 1},
  { 1,-1, 0}, 
  { 1, 1, 1}, /* 10 */
  {-1, 1, 1},
  { 1,-1, 1},
  { 1, 1,-1},
  { 0, 1, 2},
  { 2, 0, 1},
  { 1, 2, 0},
  { 0, 2, 1},
  { 1, 0, 2},
  { 2, 1, 0},
  { 0,-1, 2}, /* 20 */
  { 2, 0,-1},
  {-1, 2, 0},
  { 0,-2, 1},
  { 1, 0,-2},
  {-2, 1, 0},
  { 2, 1, 1},
  { 1, 2, 1},
  { 1, 1, 2},
  { 2,-1,-1},
  {-1, 2,-1}, /* 30 */
  {-1,-1, 2},
  { 2, 1,-1},
  {-1, 2, 1},
  { 1,-1, 2},
  { 2,-1, 1}, /* 35 */
  { 1, 2,-1},
  {-1, 1, 2},
  { 3, 1, 2},
  { 2, 3, 1},
  { 1, 2, 3}, /* 40 */
  { 3, 2, 1},
  { 1, 3, 2},
  { 2, 1, 3},
  { 3,-1, 2},
  { 2, 3,-1}, /* 45 */
  {-1, 2, 3},
  { 3,-2, 1},
  { 1, 3,-2},
  {-2, 1, 3},
  { 3,-1,-2}, /* 50 */
  {-2, 3,-1},
  {-1,-2, 3},
  { 3,-2,-1},
  {-1, 3,-2},
  {-2,-1, 3}, /* 55 */
  { 3, 1,-2},
  {-2, 3, 1},
  { 1,-2, 3},
  { 3, 2,-1},
  {-1, 3, 2}, /* 60 */
  { 2,-1, 3},
  { 1, 1, 3},
  {-1, 1, 3},
  { 1,-1, 3},
  {-1,-1, 3}, /* 65 */
  { 1, 3, 1},
  {-1, 3, 1},
  { 1, 3,-1},
  {-1, 3,-1},
  { 3, 1, 1}, /* 70 */
  { 3, 1,-1},
  { 3,-1, 1},
  { 3,-1,-1},
};

static int get_pointgroup_number_by_rotations(SPGCONST int rotations[][3][3],
					      const int num_rotations);
static PointSymmetry get_pointsymmetry(SPGCONST int rotations[][3][3],
				       const int num_rotations);
static int get_pointgroup_number(SPGCONST PointSymmetry * pointsym);
static int get_pointgroup_class_table(int table[10],
				      SPGCONST PointSymmetry * pointsym);
static int get_rotation_type(SPGCONST int rot[3][3]);
static int get_rotation_axis(SPGCONST int rot[3][3]);
static int get_orthogonal_axis(int ortho_axes[],
			       SPGCONST int proper_rot[3][3],
			       const int rot_order);
static int laue2m(int axes[3],
		  SPGCONST PointSymmetry * pointsym);

#ifdef SPGDEBUG
static int lauemmm(int axes[3],
		   SPGCONST PointSymmetry * pointsym);
static int laue4m(int axes[3],
		  SPGCONST PointSymmetry * pointsym);
static int laue4mmm(int axes[3],
		    SPGCONST PointSymmetry * pointsym);
static int laue3(int axes[3],
		 SPGCONST PointSymmetry * pointsym);
static int laue3m(int axes[3],
		  SPGCONST PointSymmetry * pointsym);
static int lauem3m(int axes[3],
		   SPGCONST PointSymmetry * pointsym);
#endif

static int laue_one_axis(int axes[3],
			 SPGCONST PointSymmetry * pointsym,
			 const int rot_order);
static int lauennn(int axes[3],
		   SPGCONST PointSymmetry * pointsym,
		   const int rot_order);
static int get_axes(int axes[3],
		    const Laue laue,
		    SPGCONST PointSymmetry * pointsym);
static void get_proper_rotation(int prop_rot[3][3],
				SPGCONST int rot[3][3]);
static void set_transformation_matrix(int tmat[3][3],
				      const int axes[3]);
static int is_exist_axis(const int axis_vec[3], const int axis_index);
static void sort_axes(int axes[3]);

/* Retrun pointgroup.number = 0 if failed */
Pointgroup ptg_get_transformation_matrix(int transform_mat[3][3],
					 SPGCONST int rotations[][3][3],
					 const int num_rotations)
{
  int i, j, pg_num;
  int axes[3];
  PointSymmetry pointsym;
  Pointgroup pointgroup;

  debug_print("ptg_get_transformation_matrix:\n");

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      transform_mat[i][j] = 0;
    }
  }
  
  pg_num = get_pointgroup_number_by_rotations(rotations, num_rotations);
  
  if (pg_num > 0) {
    pointgroup = ptg_get_pointgroup(pg_num);
    pointsym = get_pointsymmetry(rotations, num_rotations);
    get_axes(axes, pointgroup.laue, &pointsym);
    set_transformation_matrix(transform_mat, axes);
  } else {
    pointgroup = ptg_get_pointgroup(0);
  }    

  return pointgroup;
}

Pointgroup ptg_get_pointgroup(const int pointgroup_number)
{
  int i;
  Pointgroup pointgroup;
  PointgroupType pointgroup_type;

  pointgroup.number = pointgroup_number;
  pointgroup_type = pointgroup_data[pointgroup_number];
  strcpy(pointgroup.symbol, pointgroup_type.symbol);
  for (i = 0; i < 5; i++) {
    if (pointgroup.symbol[i] == ' ') {pointgroup.symbol[i] = '\0';}
  }
  pointgroup.holohedry = pointgroup_type.holohedry;
  pointgroup.laue = pointgroup_type.laue;

  debug_print("ptg_get_pointgroup: %s\n", pointgroup_type.symbol);

  return pointgroup;
}

static int get_pointgroup_number_by_rotations(SPGCONST int rotations[][3][3],
					      const int num_rotations)
{
  PointSymmetry pointsym;

  pointsym = get_pointsymmetry(rotations, num_rotations);
  return get_pointgroup_number(&pointsym);
}

static PointSymmetry get_pointsymmetry(SPGCONST int rotations[][3][3],
				       const int num_rotations)
{
  int i, j;
  PointSymmetry pointsym;

  pointsym.size = 0;
  for (i = 0; i < num_rotations; i++) {
    for (j = 0; j < pointsym.size; j++) {
      if (mat_check_identity_matrix_i3(rotations[i], pointsym.rot[j])) {
	goto escape;
      }
    }
    mat_copy_matrix_i3(pointsym.rot[pointsym.size], rotations[i]);
    pointsym.size++;
  escape:
    ;
  }

  return pointsym;
}

static int get_pointgroup_number(SPGCONST PointSymmetry * pointsym)
{
  int i, j, pg_num, counter;
  int table[10];
  PointgroupType pointgroup_type;

  debug_print("get_pointgroup_number:");


  pg_num = 0;
  
  /* Get list of point symmetry operations */
  if (! get_pointgroup_class_table(table, pointsym)) {
    goto end;
  }

  for (i = 1; i < 33; i++) {
    counter = 0;
    pointgroup_type = pointgroup_data[i];
    for (j = 0; j < 10; j++) {
      if (pointgroup_type.table[j] == table[j]) {counter++;}
    }
    if (counter == 10) {
      pg_num = i;
      break;
    }
  }

 end:
  debug_print(" %d\n", pg_num);
  return pg_num;
}

static int get_pointgroup_class_table(int table[10],
				      SPGCONST PointSymmetry * pointsym)
{
  /* Look-up table */
  /* Operation   -6 -4 -3 -2 -1  1  2  3  4  6 */
  /* Trace     -  2 -1  0  1 -3  3 -1  0  1  2 */
  /* Determinant -1 -1 -1 -1 -1  1  1  1  1  1 */

  /* table[0] = -6 axis */
  /* table[1] = -4 axis */
  /* table[2] = -3 axis */
  /* table[3] = -2 axis */
  /* table[4] = -1 axis */
  /* table[5] =  1 axis */
  /* table[6] =  2 axis */
  /* table[7] =  3 axis */
  /* table[8] =  4 axis */
  /* table[9] =  6 axis */

  int i, rot_type;

  for (i = 0; i < 10; i++) { table[i] = 0; }
  for (i = 0; i < pointsym->size; i++) {
    rot_type = get_rotation_type(pointsym->rot[i]);
    if (rot_type == -1) {
      goto err;
    } else {
      table[rot_type]++;
    }
  }
  
  return 1;

 err:
  warning_print("spglib: No point group symbol found ");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return 0;
}

static int get_rotation_type(SPGCONST int rot[3][3])
{
  int rot_type;

  if (mat_get_determinant_i3(rot) == -1) {
    switch (mat_get_trace_i3(rot)) {
    case -2: /* -6 */
      rot_type = 0;
      break;
    case -1: /* -4 */
      rot_type = 1;
      break;
    case 0:  /* -3 */
      rot_type = 2;
      break;
    case 1:  /* -2 */
      rot_type = 3;
      break;
    case -3: /* -1 */
      rot_type = 4;
      break;
    default:
      rot_type = -1;
      break;
    }
  } else {
    switch (mat_get_trace_i3(rot)) {
    case 3:  /* 1 */
      rot_type = 5;
      break;
    case -1: /* 2 */
      rot_type = 6;
      break;
    case 0:  /* 3 */
      rot_type = 7;
      break;
    case 1:  /* 4 */
      rot_type = 8;
      break;
    case 2:  /* 6 */
      rot_type = 9;
      break;
    default:
      rot_type = -1;
      break;
    }	
  }

  return rot_type;
}

static int get_axes(int axes[3],
		    const Laue laue,
		    SPGCONST PointSymmetry * pointsym)
{
  switch (laue) {
  case LAUE1:
    axes[0] = 0;
    axes[1] = 1;
    axes[2] = 2;
    break;
  case LAUE2M:
    laue2m(axes, pointsym);
    break;
  case LAUEMMM:
    lauennn(axes, pointsym, 2);
    break;
  case LAUE4M:
    laue_one_axis(axes, pointsym, 4);
    break;
  case LAUE4MMM:
    laue_one_axis(axes, pointsym, 4);
    break;
  case LAUE3:
    laue_one_axis(axes, pointsym, 3);
    break;
  case LAUE3M:
    laue_one_axis(axes, pointsym, 3);
    break;
  case LAUE6M:
    laue_one_axis(axes, pointsym, 3);
    break;
  case LAUE6MMM:
    laue_one_axis(axes, pointsym, 3);
    break;
  case LAUEM3:
    lauennn(axes, pointsym, 2);
    break;
  case LAUEM3M:
    lauennn(axes, pointsym, 4);
    break;
  default:
    break;
  }

  return 1;
}

static int laue2m(int axes[3],
		  SPGCONST PointSymmetry * pointsym)
{
  int i, num_ortho_axis, norm, min_norm, is_found, tmpval;
  int prop_rot[3][3], t_mat[3][3];
  int ortho_axes[NUM_ROT_AXES];

  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search two-fold rotation */
    if (! (mat_get_trace_i3(prop_rot) == -1)) {
      continue;
    }

    /* The first axis */
    axes[1] = get_rotation_axis(prop_rot);
    break;
  }

  /* The second axis */
  num_ortho_axis = get_orthogonal_axis(ortho_axes, prop_rot, 2);
  if (! num_ortho_axis) { goto err; }
  
  min_norm = 8;
  is_found = 0;
  for (i = 0; i < num_ortho_axis; i++) {
    norm = mat_norm_squared_i3(rot_axes[ortho_axes[i]]);
    if (norm < min_norm) {
      min_norm = norm;
      axes[0] = ortho_axes[i];
      is_found = 1;
    }
  }
  if (! is_found) { goto err; }
  
  /* The third axis */
  min_norm = 8;
  is_found = 0;
  for (i = 0; i < num_ortho_axis; i++) {
    norm = mat_norm_squared_i3(rot_axes[ortho_axes[i]]);
    if (norm < min_norm && (! (ortho_axes[i] == axes[0]))) {
      min_norm = norm;
      axes[2] = ortho_axes[i];
      is_found = 1;
    }
  }
  if (! is_found) { goto err; }

  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[2];
    axes[2] = tmpval;
  }

  return 1;

 err:
  return 0;
}

#ifdef SPGDEBUG
static int lauemmm(int axes[3],
		   SPGCONST PointSymmetry * pointsym)
{
  int i, count, axis;
  int prop_rot[3][3];


  for (i = 0; i < 3; i++) { axes[i] = -1; }
  count = 0;
  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search two-fold rotation */
    if (! (mat_get_trace_i3(prop_rot) == -1)) {
      continue;
    }

    axis = get_rotation_axis(prop_rot);
    if (! ((axis == axes[0]) ||
	   (axis == axes[1]) ||
	   (axis == axes[2]))) {
      axes[count] = axis;
      count++;
    }
  }

  sort_axes(axes);

  return 1;
}

static int laue4m(int axes[3],
		  SPGCONST PointSymmetry * pointsym)
{
  int i, num_ortho_axis, norm, min_norm, is_found, tmpval;
  int axis_vec[3];
  int prop_rot[3][3], t_mat[3][3];
  int ortho_axes[NUM_ROT_AXES];

  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search foud-fold rotation */
    if ( mat_get_trace_i3(prop_rot) == 1) {
      /* The first axis */
      axes[2] = get_rotation_axis(prop_rot);
      break;
    }
  }

  /* The second axis */
  num_ortho_axis = get_orthogonal_axis(ortho_axes, prop_rot, 4);
  if (! num_ortho_axis) { goto err; }
  
  min_norm = 8;
  is_found = 0;
  for (i = 0; i < num_ortho_axis; i++) {
    norm = mat_norm_squared_i3(rot_axes[ortho_axes[i]]);
    if (norm < min_norm) {
      min_norm = norm;
      axes[0] = ortho_axes[i];
      is_found = 1;
    }
  }
  if (! is_found) { goto err; }
  
  /* The third axis */
  mat_multiply_matrix_vector_i3(axis_vec, prop_rot, rot_axes[axes[0]]);
  is_found = 0;
  for (i = 0; i < NUM_ROT_AXES; i++) {
    if (is_exist_axis(axis_vec, i)) {
      is_found = 1;
      axes[1] = i;
      break;
    }
  }
  if (! is_found) { goto err; }

  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[1];
    axes[1] = tmpval;
  }

  return 1;

 err:
  return 0;
}

static int laue4mmm(int axes[3],
		    SPGCONST PointSymmetry * pointsym)
{
  int i, is_found, tmpval, axis;
  int prop_rot[3][3], prop_rot2[3][3], t_mat[3][3];
  int axis_vec[3];

  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search foud-fold rotation */
    if (mat_get_trace_i3(prop_rot) == 1) {
      /* The first axis */
      axes[2] = get_rotation_axis(prop_rot);
      break;
    }
  }

  is_found = 0;
  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot2, pointsym->rot[i]);

    /* Search two-fold rotation */
    if (! (mat_get_trace_i3(prop_rot2) == -1)) {
      continue;
    }

    /* The second axis */
    axis = get_rotation_axis(prop_rot2);
    if (! (axis == axes[2])) {
      axes[0] = axis;
      is_found = 1;
      break;
    }
  }
  if (! is_found) { goto err; }


  /* The third axis */
  mat_multiply_matrix_vector_i3(axis_vec, prop_rot, rot_axes[axes[0]]);
  is_found = 0;
  for (i = 0; i < NUM_ROT_AXES; i++) {
    if (is_exist_axis(axis_vec, i)) {
      is_found = 1;
      axes[1] = i;
      break;
    }
  }
  if (! is_found) { goto err; }

  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[1];
    axes[1] = tmpval;
  }

  return 1;

 err:
  return 0;
}

static int laue3(int axes[3],
		 SPGCONST PointSymmetry * pointsym)
{
  int i, num_ortho_axis, norm, min_norm, is_found, tmpval;
  int prop_rot[3][3], t_mat[3][3];
  int axis_vec[3];
  int ortho_axes[NUM_ROT_AXES];

  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search thee-fold rotation */
    if (mat_get_trace_i3(prop_rot) == 0) {
      /* The first axis */
      axes[2] = get_rotation_axis(prop_rot);
      break;
    }
  }

  /* The second axis */
  num_ortho_axis = get_orthogonal_axis(ortho_axes, prop_rot, 3);
  if (! num_ortho_axis) { goto err; }
  min_norm = 8;
  is_found = 0;
  for (i = 0; i < num_ortho_axis; i++) {
    norm = mat_norm_squared_i3(rot_axes[ortho_axes[i]]);
    if (norm < min_norm) {
      min_norm = norm;
      axes[0] = ortho_axes[i];
      is_found = 1;
    }
  }
  if (! is_found) { goto err; }
  
  /* The third axis */
  mat_multiply_matrix_vector_i3(axis_vec, prop_rot, rot_axes[axes[0]]);
  is_found = 0;
  for (i = 0; i < NUM_ROT_AXES; i++) {
    is_found = is_exist_axis(axis_vec, i);
    if (is_found == 1) {
      axes[1] = i;
      break;
    }
    if (is_found == -1) {
      axes[1] = i + NUM_ROT_AXES;
      break;
    }
  }
  if (! is_found) { goto err; }

  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[1];
    axes[1] = tmpval;
  }

  return 1;

 err:
  return 0;
}

static int laue3m(int axes[3],
		  SPGCONST PointSymmetry * pointsym)
{
  int i, is_found, tmpval, axis;
  int prop_rot[3][3], prop_rot2[3][3], t_mat[3][3];
  int axis_vec[3];

  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search three-fold rotation */
    if (mat_get_trace_i3(prop_rot) == 0) {
      /* The first axis */
      axes[2] = get_rotation_axis(prop_rot);
      debug_print("laue3m prop_rot\n");
      debug_print_matrix_i3(prop_rot);
      break;
    }
  }

  is_found = 0;
  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot2, pointsym->rot[i]);

    /* Search two-fold rotation */
    if (! (mat_get_trace_i3(prop_rot2) == -1)) {
      continue;
    }

    /* The second axis */
    axis = get_rotation_axis(prop_rot2);
    if (! (axis == axes[2])) {
      axes[0] = axis;
      is_found = 1;
      break;
    }
  }
  if (! is_found) { goto err; }

  /* The third axis */
  mat_multiply_matrix_vector_i3(axis_vec, prop_rot, rot_axes[axes[0]]);
  is_found = 0;
  for (i = 0; i < NUM_ROT_AXES; i++) {
    is_found = is_exist_axis(axis_vec, i);
    if (is_found == 1) {
      axes[1] = i;
      break;
    }
    if (is_found == -1) {
      axes[1] = i + NUM_ROT_AXES;
      break;
    }
  }
  if (! is_found) { goto err; }

  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[1];
    axes[1] = tmpval;
  }

  return 1;

 err:
  return 0;
}

static int lauem3m(int axes[3],
		   SPGCONST PointSymmetry * pointsym)
{
  int i, count, axis;
  int prop_rot[3][3];

  for (i = 0; i < 3; i++) { axes[i] = -1; }
  count = 0;
  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search four-fold rotation */
    if (! (mat_get_trace_i3(prop_rot) == 1)) {
      continue;
    }

    axis = get_rotation_axis(prop_rot);
    if (! ((axis == axes[0]) ||
	   (axis == axes[1]) ||
	   (axis == axes[2]))) {
      axes[count] = axis;
      count++;
    }
  }

  sort_axes(axes);

  return 1;
}
#endif

static int laue_one_axis(int axes[3],
			 SPGCONST PointSymmetry * pointsym,
			 const int rot_order)
{
  int i, j, num_ortho_axis, det, is_found, tmpval;
  int axis_vec[3], tmp_axes[3];
  int prop_rot[3][3], t_mat[3][3];
  int ortho_axes[NUM_ROT_AXES];

  debug_print("laue_one_axis with rot_order %d\n", rot_order);
  
  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search foud-fold rotation */
    if (rot_order == 4) {
      if (mat_get_trace_i3(prop_rot) == 1) {
	/* The first axis */
	axes[2] = get_rotation_axis(prop_rot);
	break;
      }
    }

    /* Search three-fold rotation */
    if (rot_order == 3) {
      if (mat_get_trace_i3(prop_rot) == 0) {
	/* The first axis */
	axes[2] = get_rotation_axis(prop_rot);
	break;
      }
    }
  }

  /* Candidates of the second axis */
  num_ortho_axis = get_orthogonal_axis(ortho_axes, prop_rot, rot_order);
  if (! num_ortho_axis) { goto err; }

  tmp_axes[1] = -1;
  tmp_axes[2] = axes[2];
  for (i = 0; i < num_ortho_axis; i++) {
    is_found = 0;
    tmp_axes[0] = ortho_axes[i];
    mat_multiply_matrix_vector_i3(axis_vec,
				  prop_rot,
				  rot_axes[tmp_axes[0]]);
    for (j = 0; j < num_ortho_axis; j++) {
      is_found = is_exist_axis(axis_vec, ortho_axes[j]);
      if (is_found == 1) {
	tmp_axes[1] = ortho_axes[j];
	break;
      }
      if (is_found == -1) {
	tmp_axes[1] = ortho_axes[j] + NUM_ROT_AXES;
	break;
      }
    }

    if (!is_found) { continue; }
    
    set_transformation_matrix(t_mat, tmp_axes);
    det = abs(mat_get_determinant_i3(t_mat));
    if (det < 4) { /* to avoid F-center choice det=4 */
      axes[0] = tmp_axes[0];
      axes[1] = tmp_axes[1];
      goto end;
    }
  }

 err: /* axes are not correctly found. */
  warning_print("spglib: Secondary axis is not found.");
  warning_print("(line %d, %s).\n", __LINE__, __FILE__);
  return 0;

 end:
  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[1];
    axes[1] = tmpval;
  }

  debug_print("axes[0] = %d\n", axes[0]);
  debug_print("axes[1] = %d\n", axes[1]);
  debug_print("axes[2] = %d\n", axes[2]);

  return 1;

}

static int lauennn(int axes[3],
		   SPGCONST PointSymmetry * pointsym,
		   const int rot_order)
{
  int i, count, axis;
  int prop_rot[3][3];

  for (i = 0; i < 3; i++) {
    axes[i] = -1;
  }

  count = 0;
  for (i = 0; i < pointsym->size; i++) {
    get_proper_rotation(prop_rot, pointsym->rot[i]);

    /* Search two- or four-fold rotation */
    if ((mat_get_trace_i3(prop_rot) == -1 && rot_order == 2) ||
	(mat_get_trace_i3(prop_rot) == 1 && rot_order == 4)) {
      axis = get_rotation_axis(prop_rot);
      if (! ((axis == axes[0]) ||
	     (axis == axes[1]) ||
	     (axis == axes[2]))) {
	axes[count] = axis;
	count++;
      }
    }
  }

  sort_axes(axes);

  return 1;
}

static int get_rotation_axis(SPGCONST int proper_rot[3][3])
{
  int i, axis = -1;
  int vec[3];

  /* No specific axis for I and -I */
  if (mat_check_identity_matrix_i3(proper_rot, identity)) {
    goto end;
  }

  /* Look for eigenvector = rotation axis */
  for (i = 0; i < NUM_ROT_AXES; i++) {
    mat_multiply_matrix_vector_i3(vec, proper_rot, rot_axes[i]);
    if (vec[0] == rot_axes[i][0] &&
	vec[1] == rot_axes[i][1] &&
	vec[2] == rot_axes[i][2]) {
      axis = i;
      break;
    }
  }
  
 end:
#ifdef SPGDEBUG
  if (axis == -1) {
    printf("rotation axis cound not found.\n");
  }
#endif
  
  return axis;
}

static int get_orthogonal_axis(int ortho_axes[],
			       SPGCONST int proper_rot[3][3],
			       const int rot_order)
{
  int i, num_ortho_axis;
  int vec[3];
  int sum_rot[3][3], rot[3][3];

  num_ortho_axis = 0;

  mat_copy_matrix_i3(sum_rot, identity);
  mat_copy_matrix_i3(rot, identity);
  for (i = 0; i < rot_order - 1; i++) {
    mat_multiply_matrix_i3(rot, proper_rot, rot);
    mat_add_matrix_i3(sum_rot, rot, sum_rot);
  }
  
  for (i = 0; i < NUM_ROT_AXES; i++) {
    mat_multiply_matrix_vector_i3(vec, sum_rot, rot_axes[i]);
    if (vec[0] == 0 &&
	vec[1] == 0 &&
	vec[2] == 0) {
      ortho_axes[num_ortho_axis] = i;
      num_ortho_axis++;
    }
  }

  return num_ortho_axis;
}

static void get_proper_rotation(int prop_rot[3][3],
				SPGCONST int rot[3][3])
{
  if (mat_get_determinant_i3(rot) == -1) {
    mat_multiply_matrix_i3(prop_rot, inversion, rot);
  } else {
    mat_copy_matrix_i3(prop_rot, rot);
  }
}

static void set_transformation_matrix(int tmat[3][3],
				      const int axes[3])
{
  int i, j, s[3];
  
  for (i = 0; i < 3; i++) {
    if (axes[i] < NUM_ROT_AXES) {
      s[i] = 1;
    } else {
      s[i] = -1; /* axes[i] + NUM_ROT_AXES means improper rotation. */
    }
  }
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tmat[i][j] = s[j] * rot_axes[axes[j]%NUM_ROT_AXES][i];
    }
  }
}

static int is_exist_axis(const int axis_vec[3], const int axis_index)
{
  if ((axis_vec[0] == rot_axes[axis_index][0]) &&
      (axis_vec[1] == rot_axes[axis_index][1]) &&
      (axis_vec[2] == rot_axes[axis_index][2])) { return 1; }
  if ((axis_vec[0] == -rot_axes[axis_index][0]) &&
      (axis_vec[1] == -rot_axes[axis_index][1]) &&
      (axis_vec[2] == -rot_axes[axis_index][2])) { return -1; }
  return 0;
}

static void sort_axes(int axes[3])
{
  int axis;
  int t_mat[3][3];

  if (axes[1] > axes[2]) {
    axis = axes[1];
    axes[1] = axes[2];
    axes[2] = axis;
  }

  if (axes[0] > axes[1]) {
    axis = axes[0];
    axes[0] = axes[1];
    axes[1] = axis;
  }

  if (axes[1] > axes[2]) {
    axis = axes[1];
    axes[1] = axes[2];
    axes[2] = axis;
  }

  set_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    axis = axes[1];
    axes[1] = axes[2];
    axes[2] = axis;
  }
}

