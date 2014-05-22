/* spin.h */
/* Copyright (C) 2012 Atsushi Togo */

#ifndef __spin_H__
#define __spin_H__

#include "mathfunc.h"
#include "symmetry.h"
#include "cell.h"

Symmetry  * spn_get_collinear_operation_with_symmetry(SPGCONST Symmetry *sym_nonspin,
						      SPGCONST Cell *cell,
						      const double spins[],
						      const double symprec);
Symmetry * spn_get_collinear_operation(SPGCONST Cell *cell,
				       const double spins[],
				       const double symprec);

#endif
