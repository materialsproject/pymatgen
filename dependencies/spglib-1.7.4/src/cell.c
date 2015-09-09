/* cell.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdlib.h>
#include <stdio.h>
#include "cell.h"
#include "mathfunc.h"

#include "debug.h"

/* cell->size = 0 is a sign of False */
Cell * cel_alloc_cell( const int size )
{
    Cell *cell;
    int i, j;
    
    cell = (Cell*) malloc( sizeof( Cell ) );

    for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ ) {
	cell->lattice[i][j] = 0;
      }
    }
    cell->size = size;
    
    if ( size > 0 ) {
      if ((cell->types = (int *) malloc(sizeof(int) * size)) == NULL) {
        warning_print("spglib: Memory of cell could not be allocated.");
        exit(1);
      }
      if ((cell->position =
	   (double (*)[3]) malloc(sizeof(double[3]) * size)) == NULL) {
        warning_print("spglib: Memory of cell could not be allocated.");
        exit(1);
      }
    }

    return cell;
}

void cel_free_cell( Cell * cell )
{
  if ( cell->size > 0 ) {
    free( cell->position );
    cell->position = NULL;
    free( cell->types );
    cell->types = NULL;
  }
  free ( cell );
  cell = NULL;
}

void cel_set_cell( Cell * cell,
		   SPGCONST double lattice[3][3],
		   SPGCONST double position[][3],
		   const int types[] )
{
  int i, j;
  mat_copy_matrix_d3(cell->lattice, lattice);
  for (i = 0; i < cell->size; i++) {
    for (j = 0; j < 3; j++) {
      cell->position[i][j] = position[i][j];
    }
    cell->types[i] = types[i];
  }
}

Cell * cel_copy_cell( SPGCONST Cell * cell )
{
  Cell * cell_new;
  
  cell_new = cel_alloc_cell( cell->size );
  cel_set_cell( cell_new,
		cell->lattice,
		cell->position,
		cell->types );
  return cell_new;
}

int cel_is_overlap( const double a[3],
		    const double b[3],
		    SPGCONST double lattice[3][3],
		    const double symprec )
{
  int i;
  double v_diff[3];

  for ( i = 0; i < 3; i++ ) {
    v_diff[i] = a[i] - b[i];
    v_diff[i] -= mat_Nint( v_diff[i] );
  }

  mat_multiply_matrix_vector_d3( v_diff, lattice, v_diff );
  if ( mat_norm_squared_d3( v_diff ) < symprec*symprec ) {
    return 1;
  } else {
    return 0;
  }
}

