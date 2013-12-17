/* pointgroup.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __pointgroup_H__
#define __pointgroup_H__

#include "symmetry.h"
#include "lattice.h"

typedef enum {
  NONE,
  TRICLI,
  MONOCLI,
  ORTHO,
  TETRA,
  RHOMB,
  TRIGO,
  HEXA,
  CUBIC,
} Holohedry;

typedef struct {
  char symbol[6];
  Holohedry holohedry;
  Laue laue;
  int transform_mat[3][3];
} Pointgroup;

int ptg_get_pointgroup_number( const Symmetry * symmetry );
int ptg_get_pointgroup_number_by_rotations( SPGCONST int rotations[][3][3],
					    const int num_rotations );
Pointgroup ptg_get_pointgroup( const int pointgroup_number );
Centering ptg_get_transformation_matrix( double trans_mat[3][3],
					 SPGCONST int rotations[][3][3],
					 const int num_rotations );
#endif
