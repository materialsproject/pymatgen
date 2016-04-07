/* spacegroup.h */
/* Copyright (C) 2010 Atsushi Togo */

#ifndef __spacegroup_H__
#define __spacegroup_H__

#include "cell.h"
#include "lattice.h"
#include "mathfunc.h"
#include "primitive.h"
#include "symmetry.h"

typedef struct {
  int number;
  int hall_number;
  int pointgroup_number;
  char schoenflies[7];
  char hall_symbol[17];
  char international[32];
  char international_long[20];
  char international_short[11];
  char setting[6];
  double bravais_lattice[3][3];
  double origin_shift[3];
} Spacegroup;

Primitive * spa_get_spacegroup(Spacegroup * spacegroup,
			       SPGCONST Cell * cell,
			       const double symprec);
Spacegroup spa_get_spacegroup_with_hall_number(SPGCONST Cell * primitive,
					       const int hall_number,
					       const double symprec);
#endif
