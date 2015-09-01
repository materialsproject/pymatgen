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
  int * mapping_table;
  int size;
  double tolerance;
} Primitive;

Primitive * prm_alloc_primitive(const int size);
void prm_free_primitive(Primitive * primitive);
Cell * prm_get_primitive_cell(SPGCONST Cell * cell, const double symprec);
Primitive * prm_get_primitive(SPGCONST Cell * cell, const double symprec);
#endif
