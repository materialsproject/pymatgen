/* hall_symbol.h */
/* Copyright (C) 2010 Atsushi Togo */

#ifndef __hall_symbol_H__
#define __hall_symbol_H__

#include "lattice.h"
#include "symmetry.h"
#include "mathfunc.h"

int hal_match_hall_symbol_db(double origin_shift[3],
			     SPGCONST double bravais_lattice[3][3],
			     const int hall_number,
			     const Centering centering,
			     SPGCONST Symmetry *symmetry,
			     const double symprec);

#endif
