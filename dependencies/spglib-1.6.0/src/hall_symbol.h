/* hall_symbol.h */
/* Copyright (C) 2010 Atsushi Togo */

#ifndef __hall_symbol_H__
#define __hall_symbol_H__

#include "lattice.h"
#include "symmetry.h"
#include "mathfunc.h"

int hal_get_hall_symbol( double origin_shift[3],
			 const Centering centering,
			 double bravais_lattice[3][3],
			 SPGCONST Symmetry *symmetry,
			 const double symprec );
#endif
