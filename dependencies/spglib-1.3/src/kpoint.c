/* kpoint.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include "mathfunc.h"
#include "symmetry.h"
#include "kpoint.h"

#include "debug.h"

/* #define GRID_ORDER_XYZ */
/* The addressing order of mesh grid is defined as running left */
/* element first. But when GRID_ORDER_XYZ is defined, it is changed to right */ 
/* element first. */

static PointSymmetry get_point_group_reciprocal( SPGCONST double lattice[3][3], 
						 const MatINT * rotations,
						 const int is_time_reversal,
						 const double symprec );
static PointSymmetry get_point_group_reciprocal_with_q( SPGCONST PointSymmetry * pointgroup,
							const double symprec,
							const int num_q,
							SPGCONST double qpoints[][3] );
static int get_ir_kpoints( int map[],
			   SPGCONST double kpoints[][3],
			   const int num_kpoint,
			   SPGCONST PointSymmetry * point_symmetry,
			   const double symprec );
static int get_ir_reciprocal_mesh( int grid_point[][3],
				   int map[],
				   const int mesh[3],
				   const int is_shift[3],
				   SPGCONST PointSymmetry * point_symmetry );
static Triplets * get_ir_triplets( const int mesh[3],
				   const int is_time_reversal,
				   SPGCONST double lattice[3][3],
				   const MatINT * rotations,
				   const double symprec );
static int get_ir_triplets_with_q( int weights[],
				   int grid_points[][3],
				   int third_q[],
				   const int grid_point,
				   const int mesh[3],
				   PointSymmetry * pointgroup,
				   const double symprec );
static int extract_ir_triplets_with_q( int triplets_with_q[][3], 
				       int weight_with_q[],
				       const int fixed_grid_number,
				       SPGCONST int triplets[][3],
				       const int num_triplets,
				       const int mesh[3],
				       SPGCONST PointSymmetry * point_symmetry );
static void get_grid_mapping_table( int **map_sym,
				    SPGCONST PointSymmetry * point_symmetry,
				    const int mesh[3],
				    const int is_shift[3] );
static void address_to_grid( int grid_double[3],
			     const int address,
			     const int mesh[3],
			     const int is_shift[3] );
static void get_grid_points( int grid_point[3],
			     const int grid[3],
			     const int mesh[3] );
static void get_vector_modulo( int v[3],
			       const int m[3] );
static int grid_to_address( const int grid[3],
			    const int mesh[3],
			    const int is_shift[3] );
static void free_array2D_int( int **array,
			      const int num_row );
static int ** allocate_array2d_int( const int num_row,
				    const int num_column );
static Triplets * allocate_triplets( const int num_triplets, const int mesh[3] );




int kpt_get_irreducible_kpoints( int map[],
				 SPGCONST double kpoints[][3],
				 const int num_kpoint,
				 SPGCONST double lattice[3][3],
				 const Symmetry * symmetry,
				 const int is_time_reversal,
				 const double symprec )
{
  int i;
  PointSymmetry point_symmetry;
  MatINT *rotations;
  
  rotations = mat_alloc_MatINT( symmetry->size );
  for ( i = 0; i < symmetry->size; i++ ) {
    mat_copy_matrix_i3( rotations->mat[i], symmetry->rot[i] );
  }

  point_symmetry = get_point_group_reciprocal( lattice,
					       rotations,
					       is_time_reversal,
					       symprec );
  mat_free_MatINT( rotations );

  return get_ir_kpoints(map, kpoints, num_kpoint, &point_symmetry, symprec);
}

/* grid_point (e.g. 4x4x4 mesh)                               */
/*    [[ 0  0  0]                                             */
/*     [ 1  0  0]                                             */
/*     [ 2  0  0]                                             */
/*     [-1  0  0]                                             */
/*     [ 0  1  0]                                             */
/*     [ 1  1  0]                                             */
/*     [ 2  1  0]                                             */
/*     [-1  1  0]                                             */
/*     ....      ]                                            */
/*                                                            */
/* Each value of 'map' correspnds to the index of grid_point. */
int kpt_get_irreducible_reciprocal_mesh( int grid_points[][3],
					 int map[],
					 const int mesh[3],
					 const int is_shift[3],
					 const int is_time_reversal,
					 SPGCONST double lattice[3][3],
					 const Symmetry * symmetry,
					 const double symprec )
{
  int i;
  PointSymmetry point_symmetry;
  MatINT *rotations;
  
  rotations = mat_alloc_MatINT( symmetry->size );
  for ( i = 0; i < symmetry->size; i++ ) {
    mat_copy_matrix_i3( rotations->mat[i], symmetry->rot[i] );
  }

  point_symmetry = get_point_group_reciprocal( lattice,
					       rotations,
					       is_time_reversal,
					       symprec );
  mat_free_MatINT( rotations );

  return get_ir_reciprocal_mesh( grid_points,
				 map,
				 mesh,
				 is_shift,
				 &point_symmetry );
}

void kpt_free_triplets( Triplets * t )
{
  free( t->triplets );
  t->triplets = NULL;
  free( t->weights );
  t->weights = NULL;
  free( t->mesh_points );
  t->mesh_points = NULL;
  free( t );
  t = NULL;
}

int kpt_get_stabilized_reciprocal_mesh( int grid_points[][3],
					int map[],
					const int mesh[3],
					const int is_shift[3],
					const int is_time_reversal,
					SPGCONST double lattice[3][3],
					const MatINT * rotations,
					const int num_q,
					SPGCONST double qpoints[][3],
					const double symprec )
{
  PointSymmetry pointgroup, pointgroup_q;
  
  pointgroup = get_point_group_reciprocal( lattice,
					   rotations,
					   is_time_reversal,
					   symprec );

  pointgroup_q = get_point_group_reciprocal_with_q( &pointgroup,
						    symprec,
						    num_q,
						    qpoints );

  return get_ir_reciprocal_mesh( grid_points,
				 map,
				 mesh,
				 is_shift,
				 &pointgroup_q );
}

Triplets * kpt_get_triplets_reciprocal_mesh( const int mesh[3],
					     const int is_time_reversal,
					     SPGCONST double lattice[3][3],
					     const MatINT * rotations,
					     const double symprec )
{
  return get_ir_triplets( mesh,
			  is_time_reversal,
			  lattice,
			  rotations,
			  symprec );
}

int kpt_get_ir_triplets_at_q( int weights[],
			      int grid_points[][3],
			      int third_q[],
			      const int grid_point,
			      const int mesh[3],
			      const int is_time_reversal,
			      SPGCONST double lattice[3][3],
			      const MatINT * rotations,
			      const double symprec )
{
  PointSymmetry pointgroup;

  pointgroup = get_point_group_reciprocal( lattice,
					   rotations,
					   is_time_reversal,
					   symprec );
  return get_ir_triplets_with_q( weights,
				 grid_points,
				 third_q,
				 grid_point,
				 mesh,
				 &pointgroup,
				 symprec );
}

int kpt_extract_triplets_reciprocal_mesh_at_q( int triplets_with_q[][3],
					       int weight_with_q[],
					       const int fixed_grid_number,
					       const int num_triplets,
					       SPGCONST int triplets[][3],
					       const int mesh[3],
					       const int is_time_reversal,
					       SPGCONST double lattice[3][3],
					       const MatINT * rotations,
					       const double symprec )
{
  PointSymmetry point_group;

  point_group = get_point_group_reciprocal( lattice,
					    rotations,
					    is_time_reversal,
					    symprec );

  return extract_ir_triplets_with_q( triplets_with_q,
				     weight_with_q,
				     fixed_grid_number,
				     triplets,
				     num_triplets,
				     mesh,
				     &point_group );
}

/* qpoints are used to find stabilizers (operations). */
/* num_q is the number of the qpoints. */
static PointSymmetry get_point_group_reciprocal( SPGCONST double lattice[3][3], 
						 const MatINT * rotations,
						 const int is_time_reversal,
						 const double symprec )
{
  int i, j, num_pt = 0;
  double volume;
  double rot_d[3][3], lat_inv[3][3], glat[3][3], tmp_mat[3][3], grot_d[3][3];
  MatINT *rot_reciprocal;
  PointSymmetry point_symmetry;
  SPGCONST int inversion[3][3] = {
    {-1, 0, 0 },
    { 0,-1, 0 },
    { 0, 0,-1 }
  };
  
  if ( is_time_reversal ) {
    rot_reciprocal = mat_alloc_MatINT( rotations->size * 2 );
  } else {
    rot_reciprocal = mat_alloc_MatINT( rotations->size );
  }

  volume = mat_get_determinant_d3(lattice);

  mat_inverse_matrix_d3(lat_inv, lattice, symprec);
  mat_transpose_matrix_d3(glat, lat_inv);
  mat_multiply_matrix_d3(tmp_mat, lat_inv, glat);
  
  for ( i = 0; i < rotations->size; i++ ) {
    mat_cast_matrix_3i_to_3d( rot_d, rotations->mat[ i ] );
    mat_get_similar_matrix_d3( grot_d, rot_d, tmp_mat,
			       symprec / volume / volume );
    mat_cast_matrix_3d_to_3i( rot_reciprocal->mat[ i ], grot_d );
    
    if ( is_time_reversal ) {
      mat_multiply_matrix_i3( rot_reciprocal->mat[ rotations->size+i ],
			      inversion,
			      rot_reciprocal->mat[ i ] );
    }
  }


  for ( i = 0; i < rot_reciprocal->size; i++ ) {
    for ( j = 0; j < num_pt; j++ ) {
      if ( mat_check_identity_matrix_i3( point_symmetry.rot[ j ],
					 rot_reciprocal->mat[ i ] ) ) {
	goto escape;
      }
    }
    
    mat_copy_matrix_i3( point_symmetry.rot[ num_pt ],
			rot_reciprocal->mat[ i ] );
    num_pt++;
  escape:
    ;
  }

  point_symmetry.size = num_pt;

  mat_free_MatINT( rot_reciprocal );

  return point_symmetry;
}

static PointSymmetry get_point_group_reciprocal_with_q( SPGCONST PointSymmetry * pointgroup,
							const double symprec,
							const int num_q,
							SPGCONST double qpoints[][3] )
{
  int i, j, k, l, is_all_ok=0, num_ptq = 0;
  double q_rot[3], diff[3];
  PointSymmetry pointgroup_q;

  for ( i = 0; i < pointgroup->size; i++ ) {
    for ( j = 0; j < num_q; j++ ) {
      is_all_ok = 0;
      mat_multiply_matrix_vector_id3( q_rot,
				      pointgroup->rot[ i ],
				      qpoints[ j ] );

      for ( k = 0; k < num_q; k++ ) {
	for ( l = 0; l < 3; l++ ) {
	  diff[l] = q_rot[l] - qpoints[k][l];
	  diff[l] -= mat_Nint( diff[l] );
	}
	
	if ( mat_Dabs( diff[0] ) < symprec && 
	     mat_Dabs( diff[1] ) < symprec &&
	     mat_Dabs( diff[2] ) < symprec ) {
	  is_all_ok = 1;
	  break;
	}
      }
      
      if ( ! is_all_ok ) {
	break;
      }
    }

    if ( is_all_ok ) {
      mat_copy_matrix_i3( pointgroup_q.rot[ num_ptq ], pointgroup->rot[ i ] );
      num_ptq++;
    }
  }
  pointgroup_q.size = num_ptq;

  return pointgroup_q;
}


static int get_ir_kpoints( int map[],
			   SPGCONST double kpoints[][3],
			   const int num_kpoint,
			   SPGCONST PointSymmetry * point_symmetry,
			   const double symprec )
{
  int i, j, k, l, num_ir_kpoint = 0, is_found;
  int *ir_map;
  double kpt_rot[3], diff[3];

  ir_map = (int*)malloc(num_kpoint*sizeof(int));

  for ( i = 0; i < num_kpoint; i++ ) {

    map[i] = i;

    is_found = 1;

    for ( j = 0; j < point_symmetry->size; j++ ) {
      mat_multiply_matrix_vector_id3(kpt_rot, point_symmetry->rot[j], kpoints[i]);

      for ( k = 0; k < 3; k++ ) {
	diff[k] = kpt_rot[k] - kpoints[i][k];
	diff[k] = diff[k] - mat_Nint(diff[k]);
      }

      if ( mat_Dabs(diff[0]) < symprec && 
	   mat_Dabs(diff[1]) < symprec && 
	   mat_Dabs(diff[2]) < symprec ) {
	continue;
      }
      
      for ( k = 0; k < num_ir_kpoint; k++ ) {
	mat_multiply_matrix_vector_id3(kpt_rot, point_symmetry->rot[j], kpoints[i]);

	for ( l = 0; l < 3; l++ ) {
	  diff[l] = kpt_rot[l] - kpoints[ir_map[k]][l];
	  diff[l] = diff[l] - mat_Nint(diff[l]);
	}

	if ( mat_Dabs(diff[0]) < symprec && 
	     mat_Dabs(diff[1]) < symprec && 
	     mat_Dabs(diff[2]) < symprec ) {
	  is_found = 0;
	  map[i] = ir_map[k];
	  break;
	}
      }

      if ( ! is_found )
	break;
    }

    if ( is_found ) {
      ir_map[num_ir_kpoint] = i;
      num_ir_kpoint++;
    }
  }

  free( ir_map );
  ir_map = NULL;

  return num_ir_kpoint;
}

static int get_ir_reciprocal_mesh( int grid[][3],
				   int map[],
				   const int mesh[3],
				   const int is_shift[3],
				   SPGCONST PointSymmetry * point_symmetry )
{
  /* In the following loop, mesh is doubled. */
  /* Even and odd mesh numbers correspond to */
  /* is_shift[i] = 0 and 1, respectively. */
  /* is_shift = [0,0,0] gives Gamma center mesh. */
  /* grid: reducible grid points */
  /* map: the mapping from each point to ir-point. */
  int i, j, k, l, address, address_rot, num_ir = 0;
  int grid_double[3], grid_rot[3], mesh_double[3];

  for ( i = 0; i < 3; i++ ) {
    mesh_double[i] = mesh[i] * 2;
  }

  /* "-1" means the element is not touched yet. */
  for ( i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++ ) {
    map[i] = -1;
  }

#ifndef GRID_ORDER_XYZ
  for ( i = 0; i < mesh_double[2]; i++ ) {
    if ( ( is_shift[2] && i % 2 == 0 ) ||
	 ( is_shift[2] == 0 && i % 2 != 0 ) ) 
      continue;

    for ( j = 0; j < mesh_double[1]; j++ ) {
      if ( ( is_shift[1] && j % 2 == 0 ) ||
	   ( is_shift[1] == 0 && j % 2 != 0 ) ) 
	continue;
      
      for ( k = 0; k < mesh_double[0]; k++ ) {
	if ( ( is_shift[0] && k % 2 == 0 ) ||
	     ( is_shift[0] == 0 && k % 2 != 0 ) ) 
	  continue;

	grid_double[0] = k;
	grid_double[1] = j;
	grid_double[2] = i;
#else
  for ( i = 0; i < mesh_double[0]; i++ ) {
    if ( ( is_shift[0] && i % 2 == 0 ) ||
  	 ( is_shift[0] == 0 && i % 2 != 0 ) )
      continue;

    for ( j = 0; j < mesh_double[1]; j++ ) {
      if ( ( is_shift[1] && j % 2 == 0 ) ||
  	   ( is_shift[1] == 0 && j % 2 != 0 ) )
  	continue;
      
      for ( k = 0; k < mesh_double[2]; k++ ) {
  	if ( ( is_shift[2] && k % 2 == 0 ) ||
  	     ( is_shift[2] == 0 && k % 2 != 0 ) )
  	  continue;

  	grid_double[0] = i;
  	grid_double[1] = j;
  	grid_double[2] = k;
#endif	

	address = grid_to_address( grid_double, mesh, is_shift );
	get_grid_points(grid[ address ], grid_double, mesh);

	for ( l = 0; l < point_symmetry->size; l++ ) {

	  mat_multiply_matrix_vector_i3( grid_rot, point_symmetry->rot[l], grid_double );
	  get_vector_modulo(grid_rot, mesh_double);
	  address_rot = grid_to_address( grid_rot, mesh, is_shift );

	  if ( address_rot > -1 ) { /* Invalid if even --> odd or odd --> even */
	    if ( map[ address_rot ] > -1 ) {
	      map[ address ] = map[ address_rot ];
	      break;
	    }
	  }
	}
	
	/* Set itself to the map when equivalent point */
	/* with smaller numbering could not be found. */
	if ( map[ address ] == -1 ) {
	  map[ address ] = address;
	  num_ir++;
	}
      }
    }
  }

  return num_ir;
}


/* Unique q-point triplets that conserve the momentum,  */
/* q+q'+q''=G, are obtained.                            */
/*                                                      */
/* The first q-point is selected among the ir-q-points. */
/* The second q-point is selected among the ir-q-points */
/* constrained by the first q-point (stabilizer)        */
/* The third q-point is searched through the all grid   */
/* points and is checked if it satisfies q+q'+q''=G,    */
/* here q, q', and q'' can be exchanged one another.    */
static Triplets * get_ir_triplets( const int mesh[3],
				   const int is_time_reversal,
				   SPGCONST double lattice[3][3],
				   const MatINT * rotations,
				   const double symprec )
{
  int i, j, k, l, num_ir, num_grid, weight, weight_q, count, q_2;
  int num_triplets, num_unique_q;
  int mesh_double[3], address[3], is_shift[3];
  int grid_double[3][3];
  int (*grid)[3], (*grid_local)[3];
  int *map, *map_q, *unique_q;
  int **map_sym = NULL;
  int **weight_counts;
  double stabilizer_q[1][3];
  PointSymmetry point_symmetry, point_symmetry_q;
  Triplets * tps;

  const int index_exchange[6][3] = { { 0, 1, 2 },
				     { 2, 0, 1 },
				     { 1, 2, 0 },
				     { 2, 1, 0 },
				     { 0, 2, 1 },
				     { 1, 0, 2 } };
  
  num_grid = mesh[0] * mesh[1] * mesh[2];
  map = (int*) malloc( num_grid * sizeof(int) );
  unique_q = (int*) malloc( num_grid * sizeof(int) );
  grid = (int (*)[3]) malloc( sizeof(int[3]) * num_grid );

  point_symmetry = get_point_group_reciprocal( lattice,
					       rotations,
					       is_time_reversal,
					       symprec );

  /* Only consider the gamma-point */
  for ( i = 0; i < 3; i++ ) {
    is_shift[i] = 0;
  }

  num_ir = get_ir_reciprocal_mesh( grid,
				   map,
				   mesh,
				   is_shift,
				   &point_symmetry );

  weight_counts = allocate_array2d_int( num_ir, num_grid );
  for ( i = 0; i < num_ir; i++ ) {
    for ( j = 0; j < num_grid; j++ ) {
      weight_counts[i][j] = 0;
    }
  }

  for ( i = 0; i < 3; i++ ) {
    mesh_double[i] = mesh[i] * 2;
  }

  /* Prepare triplet mapping table to enhance speed of query */
  /* 'unique_q' numbering is prepared for saving memory space */
  num_unique_q = 0;
  for ( i = 0; i < num_grid; i++ ) {
    if ( i == map[i] ) {
      unique_q[i] = num_unique_q;
      num_unique_q++;
    } 
    else {
      unique_q[i] = unique_q[map[i]];
    }
  }

  /* Prepare grid point mapping table */
  map_sym = allocate_array2d_int( point_symmetry.size, num_grid );
  get_grid_mapping_table( map_sym,
			  &point_symmetry,
			  mesh,
			  is_shift );

  /* Search triplets without considersing combination */
/* #pragma omp parallel for private( j, k, l, grid_double, point_symmetry_q, stabilizer_q, weight_q, grid_local, address, map_q, weight  ) */
  for ( i = 0; i < num_grid; i++ ) {
    if ( ! ( i == map[ i ] ) ) {
      continue;
    }

    weight = 0;
    for ( j = 0; j < num_grid; j++ ) {
      if ( i == map[j] ) {
	weight++;
      }
    }

    /* Search irreducible q-points (map_q) with a stabilizer */
    address_to_grid( grid_double[0], i, mesh, is_shift ); /* q */
    for ( j = 0; j < 3; j++ ) {
      stabilizer_q[0][j] = (double)grid_double[0][j] / mesh_double[j];
    }

    point_symmetry_q = get_point_group_reciprocal_with_q( &point_symmetry,
							  symprec,
							  1,
							  stabilizer_q );

    grid_local = (int (*)[3]) malloc( sizeof(int[3]) * num_grid );
    map_q = (int*) malloc( num_grid * sizeof(int) );
    get_ir_reciprocal_mesh( grid_local,
			    map_q,
			    mesh,
			    is_shift,
			    &point_symmetry_q);
    free( grid_local );
    grid_local = NULL;

    for ( j = 0; j < num_grid; j++ ) {
      if ( ! ( j == map_q[ j ] ) ) {
	continue;
      }

      weight_q = 0;
      for ( k = 0; k < num_grid; k++ ) {
	if ( j == map_q[k] ) {
	  weight_q++;
	}
      }

      address_to_grid( grid_double[1], j, mesh, is_shift ); /* q' */

      for ( k = 0; k < 3; k++ ) { /* q'' */
	grid_double[2][k] = - grid_double[0][k] - grid_double[1][k];
      }
      get_vector_modulo( grid_double[2], mesh_double );
      q_2 = grid_to_address( grid_double[2], mesh, is_shift );

      /* Look for irreducible triplets exchanging three q-points */
      /* and equivalent by symmetry rotations */
      for ( k = 0; k < point_symmetry.size; k++ ) {
	/* Index exchange */
	for ( l = 0; l < 6; l++ ) {
	  /* Rotated grid point addresses with index exchange */
	  address[index_exchange[l][0]] = map_sym[k][ i ];
	  address[index_exchange[l][1]] = map_sym[k][ j ];
	  address[index_exchange[l][2]] = map_sym[k][ q_2 ];

	  /* address[0] has to be one of ir-q-points. */
	  if ( address[0] == map[ address[0] ] ) {
	    /* Is the set of ddress[0] and address[1] already found? */
	    if ( weight_counts[ unique_q[ address[0] ] ][ address[1] ] ) {
	      weight_counts[ unique_q[ address[0] ] ][ address[1] ] += \
		weight * weight_q;
	      goto escape;
	    }
	  }
	}
      }

      /* Not found, then this is an irreducible triplet. */
      weight_counts[ unique_q[i] ][j] = weight * weight_q;

    escape:
      ;
    }

    free( map_q );
    map_q = NULL;

  }

  num_triplets = 0;
  for ( i = 0; i < num_grid; i++ ) {
    if ( ! ( i == map[i] ) ) {
      continue;
    }
    for ( j = 0; j < num_grid; j++ ) {
      if ( weight_counts[ unique_q[ i ] ][ j ]  ) {
	num_triplets++;
      }
    }
  }

  tps = allocate_triplets( num_triplets, mesh );
  for ( i = 0; i < num_grid; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      tps->mesh_points[i][j] = grid[i][j];
    }
  }

  count = 0;
  for ( i = 0; i < num_grid; i++ ) {
    if ( ! ( i == map[i] ) ) {
      continue;
    }
    for ( j = 0; j < num_grid; j++ ) {
      if ( weight_counts[ unique_q[ i ] ][ j ]  ) {
	tps->triplets[count][0] = i;
	tps->triplets[count][1] = j;
	address_to_grid( grid_double[0], i, mesh, is_shift ); /* q */
	address_to_grid( grid_double[1], j, mesh, is_shift ); /* q' */
	for ( l = 0; l < 3; l++ ) { /* q'' */
	  grid_double[2][l] = - grid_double[0][l] - grid_double[1][l];
	}
	get_vector_modulo( grid_double[2], mesh_double );
	tps->triplets[count][2] = grid_to_address( grid_double[2], mesh, is_shift );
	tps->weights[count] = weight_counts[ unique_q[ i ] ][ j ];
	count++;
      }
    }
  }

  free_array2D_int( map_sym, point_symmetry.size );
  free_array2D_int( weight_counts, num_ir );
  free( map );
  map = NULL;
  free( unique_q );
  unique_q = NULL;
  free( grid );
  grid = NULL;

  return tps;
}

static int get_ir_triplets_with_q( int weights[],
				   int grid_points[][3],
				   int third_q[],
				   const int grid_point,
				   const int mesh[3],
				   PointSymmetry * pointgroup,
				   const double symprec )
{
  int i, j, k, num_grid, weight_q, q_2, num_ir;
  int mesh_double[3], address[3], is_shift[3];
  int grid_double[3][3];
  int *map_q;
  int **map_sym = NULL;
  double stabilizer_q[1][3];
  PointSymmetry pointgroup_q;

  const int index_exchange[6][3] = { { 0, 1, 2 },
				     { 2, 0, 1 },
				     { 1, 2, 0 },
				     { 2, 1, 0 },
				     { 0, 2, 1 },
				     { 1, 0, 2 } };

  num_grid = mesh[0] * mesh[1] * mesh[2];

  for ( i = 0; i < 3; i++ ) {
    /* Only consider the gamma-point */
    is_shift[i] = 0;
    mesh_double[i] = mesh[i] * 2;
  }

  /* Search irreducible q-points (map_q) with a stabilizer */
  address_to_grid( grid_double[0], grid_point, mesh, is_shift ); /* q */
  for ( i = 0; i < 3; i++ ) {
    stabilizer_q[0][i] = (double)grid_double[0][i] / mesh_double[i];
  }

  pointgroup_q = get_point_group_reciprocal_with_q( pointgroup,
						    symprec,
						    1,
						    stabilizer_q );

  map_sym = allocate_array2d_int( pointgroup->size, num_grid );
  get_grid_mapping_table( map_sym,
			  pointgroup,
			  mesh,
			  is_shift );

  map_q = (int*) malloc( sizeof(int) * num_grid );
  get_ir_reciprocal_mesh( grid_points,
			  map_q,
			  mesh,
			  is_shift,
			  &pointgroup_q );

  for ( i = 0; i < num_grid; i++ ) {
    weights[i] = 0;
    third_q[i] = -1;
  }
  num_ir = 0;

  for ( i = 0; i < num_grid; i++ ) {
    if ( ! ( i == map_q[ i ] ) ) {
      continue;
    }

    weight_q = 0;
    for ( j = 0; j < num_grid; j++ ) {
      if ( i == map_q[j] ) {
	weight_q++;
      }
    }

    address_to_grid( grid_double[1], i, mesh, is_shift ); /* q' */
    for ( j = 0; j < 3; j++ ) { /* q'' */
      grid_double[2][j] = - grid_double[0][j] - grid_double[1][j];
    }
    get_vector_modulo( grid_double[2], mesh_double );
    q_2 = grid_to_address( grid_double[2], mesh, is_shift );
    third_q[i] = q_2;

    /* Look for irreducible triplets exchanging three q-points */
    /* and equivalent by symmetry rotations */
    for ( j = 0; j < pointgroup->size; j++ ) {
      /* Index exchange */
      for ( k = 0; k < 6; k++ ) {
	/* Rotated grid point addresses with index exchange */
	address[index_exchange[k][0]] = map_sym[j][ grid_point ];
	address[index_exchange[k][1]] = map_sym[j][ i ];
	address[index_exchange[k][2]] = map_sym[j][ q_2 ];

	if ( address[0] == grid_point ) {
	  /* Is the set of ddress[0] and address[1] already found? */
	  if ( weights[ address[1] ] ) {
	    weights[ address[1] ] += weight_q;
	    goto escape;
	  }
	}
      }
    }

    /* Not found, then this is an irreducible triplet. */
    weights[ i ] = weight_q;
    num_ir++;

  escape:
    ;
  }

  free( map_q );
  map_q = NULL;
  free_array2D_int( map_sym, pointgroup->size );

  return num_ir;
}

static int extract_ir_triplets_with_q( int triplets_with_q[][3], 
				       int weight_with_q[],
				       const int fixed_grid_number,
				       SPGCONST int triplets[][3],
				       const int num_triplets,
				       const int mesh[3],
				       SPGCONST PointSymmetry *point_symmetry )
{
  int i, j, k, sym_num, rest_index, num_triplets_with_q;
  int address0, address1, address1_orig, found;
  int is_shift[3];
  int num_grid;
  int **map_sym;

  num_grid = mesh[0] * mesh[1] * mesh[2];
  map_sym = allocate_array2d_int( point_symmetry->size, num_grid );

  /* Only consider the gamma-point */
  for ( i = 0; i < 3; i++ ) {
    is_shift[i] = 0;
  }

  /* Prepare mapping tables */
  get_grid_mapping_table( map_sym,
			  point_symmetry,
			  mesh,
			  is_shift );

  num_triplets_with_q = 0;

  for ( i = 0; i < num_triplets; i++ ) {
    sym_num = -1;
    for ( j = 0; j < point_symmetry->size; j++ ) {
      address0 = map_sym[j][fixed_grid_number];
      if ( triplets[i][0] == address0 ||
	   triplets[i][1] == address0 ||
	   triplets[i][2] == address0 ) {
	for ( k = 0; k < num_grid; k++ ) {
	  address1 = map_sym[j][k];
	  /* Matching indices 0 and 1 */
	  if ( ( triplets[i][0] == address0 && triplets[i][1] == address1 ) ||
	       ( triplets[i][1] == address0 && triplets[i][0] == address1 ) ) {
	    sym_num = j;
	    rest_index = 2;
	    address1_orig = k;
	    break;
	  }
	  /* Matching indices 1 and 2 */
	  if ( ( triplets[i][1] == address0 && triplets[i][2] == address1 ) ||
	       ( triplets[i][2] == address0 && triplets[i][1] == address1 ) ) {
	    sym_num = j;
	    rest_index = 0;
	    address1_orig = k;
	    break;
	  }
	  /* Matching indices 2 and 0 */
	  if ( ( triplets[i][2] == address0 && triplets[i][0] == address1 ) ||
	       ( triplets[i][0] == address0 && triplets[i][2] == address1 ) ) {
	    sym_num = j;
	    rest_index = 1;
	    address1_orig = k;
	    break;
	  }
	}
	if ( sym_num > -1 ) {
	  break;
	}
      }
    }

    /* Found? */
    if ( sym_num > -1 ) {
      for ( j = 0; j < num_grid; j++ ) {
	if ( map_sym[sym_num][j] == triplets[i][rest_index] ) {
	  triplets_with_q[num_triplets_with_q][0] = fixed_grid_number;
	  if ( j > address1_orig ) {
	    triplets_with_q[num_triplets_with_q][1] = address1_orig;
	    triplets_with_q[num_triplets_with_q][2] = j;
	  } else {
	    triplets_with_q[num_triplets_with_q][2] = address1_orig;
	    triplets_with_q[num_triplets_with_q][1] = j;
	  }
	  num_triplets_with_q++;
	  break;
	}
      }
    }
  }

  for ( i = 0; i < num_triplets_with_q; i++ ) {
    weight_with_q[i] = 0;
  }

  for ( i = 0; i < num_grid; i++ ) {
    found = 0;
    for ( j = 0; j < num_triplets_with_q; j++ ) {
      for ( k = 0; k < point_symmetry->size; k++ ) {

	if ( map_sym[k][fixed_grid_number] == triplets_with_q[j][0] ) {
	  if ( map_sym[k][i] == triplets_with_q[j][1] ||
	       map_sym[k][i] == triplets_with_q[j][2] ) {
	    weight_with_q[j]++;
	    found = 1;
	    break;
	  }	  
	}
	if ( map_sym[k][fixed_grid_number] == triplets_with_q[j][1] ) {
	  if ( map_sym[k][i] == triplets_with_q[j][2] ||
	       map_sym[k][i] == triplets_with_q[j][0] ) {
	    weight_with_q[j]++;
	    found = 1;
	    break;
	  }	  
	}
	if ( map_sym[k][fixed_grid_number] == triplets_with_q[j][2] ) {
	  if ( map_sym[k][i] == triplets_with_q[j][0] ||
	       map_sym[k][i] == triplets_with_q[j][1] ) {
	    weight_with_q[j]++;
	    found = 1;
	    break;
	  }
	}
      }
      if ( found ) {
	break;
      }
    }
    if ( ! found ) {
      warning_print("spglib: Unexpected behavior in extract_ir_triplets_with_q ");
      warning_print("(line %d, %s).\n", __LINE__, __FILE__);
      num_triplets_with_q = 0;
      break;
    }
  }

  free_array2D_int( map_sym, point_symmetry->size );
  return num_triplets_with_q;
}

static void get_grid_mapping_table( int **map_sym,
				    SPGCONST PointSymmetry *point_symmetry,
				    const int mesh[3],
				    const int is_shift[3] )
{
  int i, j;
  int grid_rot[3], grid_double[3], mesh_double[3];

  for ( i = 0; i < 3; i++ ) {
    mesh_double[i] = mesh[i] * 2;
  }

  for ( i = 0; i < point_symmetry->size; i++ ) {
    for ( j = 0; j < mesh[0]*mesh[1]*mesh[2]; j++ ) {
      address_to_grid( grid_double, j, mesh, is_shift );
      mat_multiply_matrix_vector_i3( grid_rot,
				     point_symmetry->rot[i],
				     grid_double );
      get_vector_modulo( grid_rot, mesh_double );
      map_sym[i][j] = grid_to_address( grid_rot, mesh, is_shift );
    }
  }
}  


static int grid_to_address( const int grid_double[3],
			    const int mesh[3],
			    const int is_shift[3] )
{
  int i, grid[3];

  for ( i = 0; i < 3; i++ ) {
    if ( grid_double[i] % 2 == 0 && (! is_shift[i])  ) {
      grid[i] = grid_double[i] / 2;
    } else {
      if ( grid_double[i] % 2 != 0 && is_shift[i] ) {
	grid[i] = ( grid_double[i] - 1 ) / 2;
      } else {
	return -1;
      }
    }
  }

#ifndef GRID_ORDER_XYZ
  return grid[2] * mesh[0] * mesh[1] + grid[1] * mesh[0] + grid[0];
#else
  return grid[0] * mesh[1] * mesh[2] + grid[1] * mesh[2] + grid[2];
#endif  
}

static void address_to_grid( int grid_double[3],
			     const int address,
			     const int mesh[3],
			     const int is_shift[3] )
{
  int i;
  int grid[3];

#ifndef GRID_ORDER_XYZ
  grid[2] = address / ( mesh[0] * mesh[1] );
  grid[1] = ( address - grid[2] * mesh[0] * mesh[1] ) / mesh[0];
  grid[0] = address % mesh[0];
#else
  grid[0] = address / ( mesh[1] * mesh[2] );
  grid[1] = ( address - grid[0] * mesh[1] * mesh[2] ) / mesh[2];
  grid[2] = address % mesh[2];
#endif

  for ( i = 0; i < 3; i++ ) {
    grid_double[i] = grid[i] * 2 + is_shift[i];
  }
}

static void get_grid_points( int grid[3],
			     const int grid_double[3],
			     const int mesh[3] )
{
  int i;

  for ( i = 0; i < 3; i++ ) {
    if ( grid_double[i] % 2 == 0 ) {
      grid[i] = grid_double[i] / 2;
    } else {
      grid[i] = ( grid_double[i] - 1 ) / 2;
    }

#ifndef GRID_BOUNDARY_AS_NEGATIVE
    grid[i] = grid[i] - mesh[i] * ( grid[i] > mesh[i] / 2 );
#else
    grid[i] = grid[i] - mesh[i] * ( grid[i] >= mesh[i] / 2 );
#endif
  }  
}

static void get_vector_modulo( int v[3],
			       const int m[3] )
{
  int i;

  for ( i = 0; i < 3; i++ ) {
    v[i] = v[i] % m[i];

    if ( v[i] < 0 )
      v[i] += m[i];
  }
}

static void free_array2D_int( int **array,
			      const int num_row )
{
  int i;
  for ( i = 0; i < num_row; i++ ) {
    free( array[i] );
    array[i] = NULL;
  }
  free( array );
  array = NULL;
}

static int ** allocate_array2d_int( const int num_row,
				    const int num_column )
{
  int i;
  int **array;
  
  array = (int**) malloc( num_row * sizeof(int*) );
  for (i = 0; i < num_row; i++) {
    array[i] = (int*) malloc( num_column * sizeof(int) );
  }
  return array;
}

static Triplets * allocate_triplets( const int num_triplets, const int mesh[3] )
{
  int i, num_grid;
  Triplets * tps;

  num_grid = mesh[0] * mesh[1] * mesh[2];
  tps = (Triplets*) malloc( sizeof( Triplets ) );
  tps->size = num_triplets;
  tps->triplets = (int (*)[3]) malloc( sizeof(int[3]) * num_triplets );
  tps->weights = (int*) malloc( sizeof(int) * num_triplets );
  tps->mesh_points = (int (*)[3]) malloc( sizeof(int[3]) * num_grid );
  for ( i = 0; i < 3; i++ ) {
    tps->mesh[i] = mesh[i];
  }
 
  return tps;
}


