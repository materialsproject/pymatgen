/* pointgroup.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <string.h>
#include <stdio.h>
#include "lattice.h"
#include "pointgroup.h"
#include "symmetry.h"
#include "mathfunc.h"

#include "debug.h"

#define NUM_ROT_AXES 73

typedef struct {
  int table[10];
  char symbol[6];
  Holohedry holohedry;
  Laue laue;
} PointgroupType;

static PointgroupType pointgroup_data[32] = {
  {
    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
    "1    ",
    TRICLI,
    LAUE1,
  },
  {
    {0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
    "-1   ",
    TRICLI,
    LAUE1,
  },
  {
    {0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
    "2    ",
    MONOCLI,
    LAUE2M,
  },
  {
    {0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
    "m    ",
    MONOCLI,
    LAUE2M,
  },
  {
    {0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
    "2/m  ",
    MONOCLI,
    LAUE2M,
  },
  {
    {0, 0, 0, 0, 0, 1, 3, 0, 0, 0},
    "222  ",
    ORTHO,
    LAUEMMM,
  },
  {
    {0, 0, 0, 2, 0, 1, 1, 0, 0, 0},
    "mm2  ",
    ORTHO,
    LAUEMMM,
  },
  {
    {0, 0, 0, 3, 1, 1, 3, 0, 0, 0},
    "mmm  ",
    ORTHO,
    LAUEMMM,
  },
  {
    {0, 0, 0, 0, 0, 1, 1, 0, 2, 0},
    "4    ",
    TETRA,
    LAUE4M,
  },
  {
    {0, 2, 0, 0, 0, 1, 1, 0, 0, 0},
    "-4   ",
    TETRA,
    LAUE4M,
  },
  {
    {0, 2, 0, 1, 1, 1, 1, 0, 2, 0},
    "4/m  ",
    TETRA,
    LAUE4M,
  },
  {
    {0, 0, 0, 0, 0, 1, 5, 0, 2, 0},
    "422  ",
    TETRA,
    LAUE4MMM,
  },
  {
    {0, 0, 0, 4, 0, 1, 1, 0, 2, 0},
    "4mm  ",
    TETRA,
    LAUE4MMM,
  },
  {
    {0, 2, 0, 2, 0, 1, 3, 0, 0, 0},
    "-42m ",
    TETRA,
    LAUE4MMM,
  },
  {
    {0, 2, 0, 5, 1, 1, 5, 0, 2, 0},
    "4/mmm",
    TETRA,
    LAUE4MMM,
  },
  {
    {0, 0, 0, 0, 0, 1, 0, 2, 0, 0},
    "3    ",
    TRIGO,
    LAUE3,
  },
  {
    {0, 0, 2, 0, 1, 1, 0, 2, 0, 0},
    "-3   ",
    TRIGO,
    LAUE3,
  },
  {
    {0, 0, 0, 0, 0, 1, 3, 2, 0, 0},
    "32   ",
    TRIGO,
    LAUE3M,
  },
  {
    {0, 0, 0, 3, 0, 1, 0, 2, 0, 0},
    "3m   ",
    TRIGO,
    LAUE3M,
  },
  {
    {0, 0, 2, 3, 1, 1, 3, 2, 0, 0},
    "-3m  ",
    TRIGO,
    LAUE3M,
  },
  {
    {0, 0, 0, 0, 0, 1, 1, 2, 0, 2},
    "6    ",
    HEXA,
    LAUE6M,
  },
  {
    {2, 0, 0, 1, 0, 1, 0, 2, 0, 0},
    "-6   ",
    HEXA,
    LAUE6M,
  },
  {
    {2, 0, 2, 1, 1, 1, 1, 2, 0, 2},
    "6/m  ",
    HEXA,
    LAUE6M,
  },
  {
    {0, 0, 0, 0, 0, 1, 7, 2, 0, 2},
    "622  ",
    HEXA,
    LAUE6MMM,
  },
  {
    {0, 0, 0, 6, 0, 1, 1, 2, 0, 2},
    "6mm  ",
    HEXA,
    LAUE6MMM,
  },
  {
    {2, 0, 0, 4, 0, 1, 3, 2, 0, 0},
    "-6m2 ",
    HEXA,
    LAUE6MMM,
  },
  {
    {2, 0, 2, 7, 1, 1, 7, 2, 0, 2},
    "6/mmm",
    HEXA,
    LAUE6MMM,
  },
  {
    {0, 0, 0, 0, 0, 1, 3, 8, 0, 0},
    "23   ",
    CUBIC,
    LAUEM3,
  },
  {
    {0, 0, 8, 3, 1, 1, 3, 8, 0, 0},
    "m-3  ",
    CUBIC,
    LAUEM3,
  },
  {
    {0, 0, 0, 0, 0, 1, 9, 8, 6, 0},
    "432  ",
    CUBIC,
    LAUEM3M,
  },
  {
    {0, 6, 0, 6, 0, 1, 3, 8, 0, 0},
    "-43m ",
    CUBIC,
    LAUEM3M,
  },
  {
    {0, 6, 8, 9, 1, 1, 9, 8, 6, 0},
    "m-3m ",
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

static void set_transformation_matrix(Pointgroup * pointgroup,
				      SPGCONST int rotations[][3][3],
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

#ifdef DEBUG
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
static void get_transformation_matrix(int tmat[3][3],
				      const int axes[3]);
static int is_exist_axis(const int axis_vec[3], const int axis_index);
static void sort_axes(int axes[3]);


int ptg_get_pointgroup_number(const Symmetry * symmetry)
{
  PointSymmetry pointsym;

  pointsym = get_pointsymmetry(symmetry->rot,
			       symmetry->size);
  return get_pointgroup_number(&pointsym);
}

int ptg_get_pointgroup_number_by_rotations(SPGCONST int rotations[][3][3],
					   const int num_rotations)
{
  PointSymmetry pointsym;

  pointsym = get_pointsymmetry(rotations, num_rotations);
  return get_pointgroup_number(&pointsym);
}

Pointgroup ptg_get_pointgroup(const int pointgroup_number)
{
  Pointgroup pointgroup;
  PointgroupType pointgroup_type;
  
  pointgroup_type = pointgroup_data[ pointgroup_number ];
  strcpy(pointgroup.symbol, pointgroup_type.symbol);
  pointgroup.holohedry = pointgroup_type.holohedry;
  pointgroup.laue = pointgroup_type.laue;

  debug_print("ptg_get_pointgroup: %s\n", pointgroup_type.symbol);

  return pointgroup;
}

Centering ptg_get_transformation_matrix(double trans_mat[3][3],
					SPGCONST int rotations[][3][3],
					const int num_rotations)
{
  int pg_num;
  double correction_mat[3][3];
  Centering centering;
  Pointgroup pointgroup;

  debug_print("ptg_get_transformation_matrix:\n");

  pg_num = ptg_get_pointgroup_number_by_rotations(rotations,
						  num_rotations);
  pointgroup = ptg_get_pointgroup(pg_num);
  set_transformation_matrix(&pointgroup, rotations, num_rotations);

  debug_print("transformation matrix:\n");
  debug_print_matrix_i3(pointgroup.transform_mat);

  /* Centering is not determined only from symmetry operations */
  /* sometimes. Therefore centering and transformation matrix are */
  /* related. */
  centering = lat_get_centering(correction_mat,
				pointgroup.transform_mat,
				pointgroup.laue);

  mat_multiply_matrix_id3(trans_mat,
			  pointgroup.transform_mat,
			  correction_mat);

  debug_print("correction matrix:\n");
  debug_print_matrix_d3(correction_mat);

  return centering;
}

/* pointgroup is modified. */
static void set_transformation_matrix(Pointgroup * pointgroup,
				      SPGCONST int rotations[][3][3],
				      const int num_rotations)
{
  int axes[3];
  int transform_mat[3][3];
  PointSymmetry pointsym;

  pointsym = get_pointsymmetry(rotations, num_rotations);
  get_axes(axes, pointgroup->laue, &pointsym);
  get_transformation_matrix(transform_mat, axes);
  mat_copy_matrix_i3(pointgroup->transform_mat, transform_mat);
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
  
  /* Get list of point symmetry operations */
  if (! get_pointgroup_class_table(table, pointsym)) {
    pg_num = -1;
    goto end;
  }

  pg_num = -1;
  for (i = 0; i < 32; i++) {
    counter = 0;
    pointgroup_type = pointgroup_data[ i ];
    for (j = 0; j < 10; j++) {
      if (pointgroup_type.table[j] == table[j]) { counter++; }
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

  get_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    tmpval = axes[0];
    axes[0] = axes[2];
    axes[2] = tmpval;
  }

  return 1;

 err:
  return 0;
}

#ifdef DEBUG
static int lauemmm(int axes[3],
		   SPGCONST PointSymmetry * pointsym)
{
  int i, count, axis, tmpval;
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

  get_transformation_matrix(t_mat, axes);
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

  get_transformation_matrix(t_mat, axes);
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

  get_transformation_matrix(t_mat, axes);
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

  get_transformation_matrix(t_mat, axes);
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
  int i, count, axis, tmpval;
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
  int i, j, num_ortho_axis, det, min_det, is_found, tmpval;
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

  tmp_axes[2] = axes[2];
  min_det = 4;
  is_found = 0;
  for (i = 0; i < num_ortho_axis; i++) {
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

    get_transformation_matrix(t_mat, tmp_axes);
    det = mat_get_determinant_i3(t_mat);
    if (det < 0) { det = -det; }
    if (det < min_det) {
      min_det = det;
      axes[0] = tmp_axes[0];
      axes[1] = tmp_axes[1];
    }

    if (is_found) { goto end; }
  }

 err: /* axes are not correctly found. */
  return 0;

 end:
  get_transformation_matrix(t_mat, axes);
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

  for (i = 0; i < 3; i++) { axes[i] = -1; }
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
#ifdef DEBUG
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

static void get_transformation_matrix(int tmat[3][3],
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

  get_transformation_matrix(t_mat, axes);
  if (mat_get_determinant_i3(t_mat) < 0) {
    axis = axes[1];
    axes[1] = axes[2];
    axes[2] = axis;
  }
}

