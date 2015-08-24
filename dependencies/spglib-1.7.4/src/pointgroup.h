/* pointgroup.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __pointgroup_H__
#define __pointgroup_H__

#include "symmetry.h"
#include "lattice.h"

typedef enum {
  HOLOHEDRY_NONE,
  TRICLI,
  MONOCLI,
  ORTHO,
  TETRA,
  TRIGO,
  HEXA,
  CUBIC,
} Holohedry;

typedef struct {
  int number;
  char symbol[6];
  Holohedry holohedry;
  Laue laue;
} Pointgroup;

Pointgroup ptg_get_transformation_matrix(int transform_mat[3][3],
					 SPGCONST int rotations[][3][3],
					 const int num_rotations);
Pointgroup ptg_get_pointgroup(const int pointgroup_number);
#endif
