/* site_symmetry.h */
/* Copyright (C) 2011 Atsushi Togo */

#ifndef __site_symmetry_H__
#define __site_symmetry_H__

#include "cell.h"
#include "mathfunc.h"
#include "symmetry.h"

VecDBL * ssm_get_exact_positions(int * wyckoffs,
				 int * equiv_atoms,
				 SPGCONST Cell * bravais,
				 SPGCONST Symmetry * conv_sym,
				 const int hall_number,
				 const double symprec);

#endif
