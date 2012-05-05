/* primitive.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __primitive_H__
#define __primitive_H__

#include "symmetry.h"
#include "cell.h"
#include "mathfunc.h"

Cell * prm_get_primitive(SPGCONST Cell * cell,
			 const double symprec);
Cell * prm_get_primitive_with_mapping_table(int * mapping_table,
					    SPGCONST Cell * cell,
					    const double symprec);
Cell * prm_get_primitive_with_pure_translations(SPGCONST Cell * cell,
						const VecDBL *pure_trans,
						const double symprec);
Cell * prm_get_primitive_with_all(int * mapping_table,
				  SPGCONST Cell * cell,
				  const VecDBL *pure_trans,
				  const double symprec);
double prm_get_current_tolerance(void);
#endif
