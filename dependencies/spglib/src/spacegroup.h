/* spacegroup.h */
/* Copyright (C) 2010 Atsushi Togo */

#ifndef __spacegroup_H__
#define __spacegroup_H__

#include "cell.h"
#include "lattice.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "symmetry.h"

typedef struct {
  int number;
  int hall_number;
  char schoenflies[7];
  char hall_symbol[17];
  char international[32];
  char international_long[20];
  char international_short[11];
  Holohedry holohedry;
  double bravais_lattice[3][3];
  double origin_shift[3];
} Spacegroup;

Spacegroup spa_get_spacegroup(SPGCONST Cell * cell,
			      const double symprec);
Spacegroup spa_get_spacegroup_with_primitive(SPGCONST Cell * primitive,
					     const double symprec);
Symmetry * spa_get_conventional_symmetry(SPGCONST double transform_mat[3][3],
					 const Centering centering,
					 const Symmetry *primitive_sym);

#endif
