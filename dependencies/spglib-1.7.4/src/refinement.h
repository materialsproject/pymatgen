/* bravais.h */
/* Copyright (C) 2011 Atsushi Togo */

#ifndef __refinement_H__
#define __refinement_H__

#include "cell.h"
#include "mathfunc.h"
#include "spacegroup.h"
#include "symmetry.h"

Symmetry *
ref_get_refined_symmetry_operations(SPGCONST Cell * cell,
				    SPGCONST Cell * primitive,
				    SPGCONST Spacegroup * spacegroup,
				    const double symprec);
Cell * ref_get_Wyckoff_positions(int * wyckoffs,
				 int * equiv_atoms,
				 SPGCONST Cell * primitive,
				 SPGCONST Cell * cell,
				 SPGCONST Spacegroup * spacegroup,
				 SPGCONST Symmetry * symmetry,
				 const int * mapping_table,
				 const double symprec);

#endif
