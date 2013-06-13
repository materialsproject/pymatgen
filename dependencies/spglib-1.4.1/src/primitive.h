/* primitive.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __primitive_H__
#define __primitive_H__

#include "symmetry.h"
#include "cell.h"
#include "mathfunc.h"

typedef struct {
  Cell *cell;
  VecDBL *pure_trans;
} Primitive;

Cell * prm_get_primitive(SPGCONST Cell * cell,
			 const double symprec);
Cell * prm_get_primitive_with_mapping_table(int * mapping_table,
					    SPGCONST Cell * cell,
					    const double symprec);
Primitive prm_get_primitive_and_pure_translations(SPGCONST Cell * cell,
						  const double symprec);
double prm_get_current_tolerance(void);
#endif
