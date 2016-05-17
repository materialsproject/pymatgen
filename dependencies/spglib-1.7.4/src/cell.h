/* cell.h */
/* Copyright (C) 2008 Atsushi Togo */

#ifndef __cell_H__
#define __cell_H__

#include "mathfunc.h"

typedef struct {
    int size;
    double lattice[3][3];
    int *types;
    double (*position)[3];
} Cell;

Cell *cel_alloc_cell( const int size );
void cel_free_cell( Cell * cell );
void cel_set_cell( Cell * cell,
		   SPGCONST double lattice[3][3],
		   SPGCONST double position[][3],
		   const int types[] );
Cell * cel_copy_cell( SPGCONST Cell * cell );
int cel_is_overlap( const double a[3],
		    const double b[3],
		    SPGCONST double lattice[3][3],
		    const double symprec );

#endif
