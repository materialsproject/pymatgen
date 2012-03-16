#include "spglib.h"
#include <stdio.h>
#include <math.h>

void test_spg_get_symmetry();
void test_spg_get_symmetry_with_collinear_spin();
void test_spg_get_multiplicity();
void test_spg_find_primitive();
void test_spg_get_international();
void test_spg_get_schoenflies();
void test_spg_refine_cell();
void test_spg_get_dataset();
void test_spg_get_ir_reciprocal_mesh();

int main()
{
  test_spg_find_primitive();
  test_spg_get_multiplicity();
  test_spg_get_symmetry();
  test_spg_get_symmetry_with_collinear_spin();
  test_spg_get_international();
  test_spg_get_schoenflies();
  test_spg_refine_cell();
  test_spg_get_dataset();
  test_spg_get_ir_reciprocal_mesh();

  return 0;
}

void test_spg_find_primitive()
{
  double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} };
  double position[][3] = {
    {0, 0, 0},
    {0.5, 0.5, 0.5}
  };
  int types[] = { 1, 1 };
  int i, num_atom = 2, num_primitive_atom;
  double symprec = 1e-5;
  
  /* lattice, position, and types are overwirtten. */
  printf("*** Example of spg_find_primitive (BCC unitcell --> primitive) ***:\n");
  num_primitive_atom = spg_find_primitive(lattice, position, types, num_atom, symprec);
  if ( num_primitive_atom == 0 ) {
    printf("Primitive cell was not found.\n");
  } else { 
    printf("Lattice parameter:\n");
    for (i = 0; i < 3; i++) {
      printf("%f %f %f\n", lattice[0][i], lattice[1][i], lattice[2][i]);
    }
    printf("Atomic positions:\n");
    for (i=0; i<num_primitive_atom; i++) {
      printf("%d: %f %f %f\n", types[i], position[i][0], position[i][1],
	     position[i][2]);
    }
  }
}

void test_spg_refine_cell()
{
  double lattice[3][3] = { {0, 2, 2}, {2, 0, 2}, {2, 2, 0} };

  /* 4 times larger memory space must be prepared. */
  double position[4][3];
  int types[4];
  
  int i, num_atom_bravais, num_atom = 1;
  double symprec = 1e-5;

  position[0][0] = 0;
  position[0][1] = 0;
  position[0][2] = 0;
  types[0] = 1;
  
  /* lattice, position, and types are overwirtten. */
  printf("*** Example of spg_refine_cell ***:\n");
  num_atom_bravais = spg_refine_cell( lattice,
				      position,
				      types,
				      num_atom,
				      symprec );
  printf("Lattice parameter:\n");
  for ( i = 0; i < 3; i++ ) {
    printf("%f %f %f\n", lattice[0][i], lattice[1][i], lattice[2][i]);
  }
  printf("Atomic positions:\n");
  for ( i = 0; i<num_atom_bravais; i++ ) {
    printf("%d: %f %f %f\n", types[i], position[i][0], position[i][1],
	   position[i][2]);
  }
}

void test_spg_get_international()
{
  double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
  double position[][3] =
    {
      {0,0,0},
      {0.5,0.5,0.5},
      {0.3,0.3,0},
      {0.7,0.7,0},
      {0.2,0.8,0.5},
      {0.8,0.2,0.5},
    };
  int types[] = {1,1,2,2,2,2};
  int num_spg, num_atom = 6;
  char symbol[21];
  
  num_spg = spg_get_international(symbol, lattice, position, types, num_atom, 1e-5);
  printf("*** Example of spg_get_international ***:\n");
  if ( num_spg > 0 ) {
    printf("%s (%d)\n", symbol, num_spg);
  } else {
    printf("Space group could not be found.\n");
  }
}

void test_spg_get_schoenflies()
{
  double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
  double position[][3] =
    {
      {0,0,0},
      {0.5,0.5,0.5},
      {0.3,0.3,0},
      {0.7,0.7,0},
      {0.2,0.8,0.5},
      {0.8,0.2,0.5},
    };
  int types[] = {1,1,2,2,2,2};
  int num_atom = 6;
  char symbol[10];
  
  spg_get_schoenflies(symbol, lattice, position, types, num_atom, 1e-5);
  printf("*** Example of spg_get_schoenflies ***:\n");
  printf("Schoenflies: %s\n", symbol);
}

void test_spg_get_multiplicity()
{
  double lattice[3][3] = { {4, 0, 0}, {0, 4, 0}, {0, 0, 4} };
  double position[][3] = {
    {0, 0, 0},
    {0.5, 0.5, 0.5}
  };
  int types[] = { 1, 2 };
  int num_atom = 2;
  int size;

  size = spg_get_multiplicity(lattice, position, types, num_atom, 1e-5);

  printf("*** Example of spg_get_multiplicity ***:\n");
  printf("Number of symmetry operations: %d\n", size);
}


void test_spg_get_symmetry()
{
  double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
  double position[][3] =
    {
      {0,0,0},
      {0.5,0.5,0.25},
      {0.3,0.3,0},
      {0.7,0.7,0},
      {0.2,0.8,0.25},
      {0.8,0.2,0.25},
      {0,0,0.5},
      {0.5,0.5,0.75},
      {0.3,0.3,0.5},
      {0.7,0.7,0.5},
      {0.2,0.8,0.75},
      {0.8,0.2,0.75}
    };
  int types[] = {1,1,2,2,2,2,1,1,2,2,2,2};
  int num_atom = 12;
  int max_size = 50;
  int i, j, size;
  int rotation[max_size][3][3];
  double translation[max_size][3];

  double origin_shift[3] = { 0.1, 0.1, 0 };
  for ( i = 0; i < num_atom; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      position[i][j] += origin_shift[j];
    }
  }

  printf("*** Example of spg_get_symmetry (Rutile two unit cells) ***:\n");
  size = spg_get_symmetry( rotation,
			   translation,
			   max_size,
			   lattice,
			   position,
			   types,
			   num_atom,
			   1e-5 );
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++)
      printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    printf("%f %f %f\n", translation[i][0], translation[i][1],
	   translation[i][2]);
  }

}

void test_spg_get_symmetry_with_collinear_spin() {
  double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,4}};
  double position[][3] =
    {
      {0,0,0},
      {0.5,0.5,0.5}
    };	
  int types[] = { 1, 1 };
  double spins[2];
  int num_atom = 2;
  int max_size = 300;
  int i, j, size;
  int rotation[max_size][3][3];
  double translation[max_size][3];

  printf("*** Example of spg_get_symmetry_with_spin (BCC ferro) ***:\n");
  spins[0] = 1;
  spins[1] = 1;
  size = spg_get_symmetry_with_collinear_spin( rotation,
					       translation,
					       max_size,
					       lattice,
					       position,
					       types,
					       spins,
					       num_atom,
					       1e-5 );
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++)
      printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    printf("%f %f %f\n", translation[i][0], translation[i][1],
           translation[i][2]);
  }

  printf("*** Example of spg_get_symmetry_with_spin (BCC antiferro) ***:\n");
  spins[0] = 1;
  spins[1] = -1;
  size = spg_get_symmetry_with_collinear_spin( rotation,
					       translation,
					       max_size,
					       lattice,
					       position,
					       types,
					       spins,
					       num_atom,
					       1e-5 );
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++)
      printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    printf("%f %f %f\n", translation[i][0], translation[i][1],
           translation[i][2]);
  }

  printf("*** Example of spg_get_symmetry_with_spin (BCC broken spin) ***:\n");
  spins[0] = 1;
  spins[1] = 2;
  size = spg_get_symmetry_with_collinear_spin( rotation,
					       translation,
					       max_size,
					       lattice,
					       position,
					       types,
					       spins,
					       num_atom,
					       1e-5 );
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++)
      printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    printf("%f %f %f\n", translation[i][0], translation[i][1],
           translation[i][2]);
  }
}


void test_spg_get_dataset()
{
  SpglibDataset *dataset;
  double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
  double origin_shift[3] = { 0.1, 0.1, 0 };
  double position[][3] =
    {
      {0,0,0},
      {0.5,0.5,0.25},
      {0.3,0.3,0},
      {0.7,0.7,0},
      {0.2,0.8,0.25},
      {0.8,0.2,0.25},
      {0,0,0.5},
      {0.5,0.5,0.75},
      {0.3,0.3,0.5},
      {0.7,0.7,0.5},
      {0.2,0.8,0.75},
      {0.8,0.2,0.75}
    };
  int types[] = {1,1,2,2,2,2,1,1,2,2,2,2};
  int num_atom = 12;
  int i, j, size;
  const char *wl = "abcdefghijklmnopqrstuvwxyz";

  for ( i = 0; i < num_atom; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      position[i][j] += origin_shift[j];
    }
  }

  printf("*** Example of spg_get_dataset (Rutile two unit cells) ***:\n");
  dataset = spg_get_dataset( lattice,
			     position,
			     types,
			     num_atom,
			     1e-5 );
  
  printf("International: %s (%d)\n", dataset->international_symbol, dataset->spacegroup_number );
  printf("Hall symbol:   %s\n", dataset->hall_symbol );
  printf("Transformation matrix:\n");
  for ( i = 0; i < 3; i++ ) {
    printf("%f %f %f\n",
	   dataset->transformation_matrix[i][0],
	   dataset->transformation_matrix[i][1],
	   dataset->transformation_matrix[i][2]);
  }
  printf("Wyckoff letters:\n");
  for ( i = 0; i < dataset->n_atoms; i++ ) {
    printf("%c ", wl[dataset->wyckoffs[i]]);
  }
  printf("\n");
  printf("Equivalent atoms:\n");
  for ( i = 0; i < dataset->n_atoms; i++ ) {
    printf("%d ", dataset->equivalent_atoms[i]);
  }
  printf("\n");
  
  for ( i = 0; i < dataset->n_operations; i++ ) {
    printf("--- %d ---\n", i + 1);
    for ( j = 0; j < 3; j++ ) {
      printf("%2d %2d %2d\n",
	     dataset->rotations[i][j][0],
	     dataset->rotations[i][j][1],
	     dataset->rotations[i][j][2]);
    }
    printf("%f %f %f\n",
	   dataset->translations[i][0],
	   dataset->translations[i][1],
	   dataset->translations[i][2]);
  }

  spg_free_dataset( dataset );

}

void test_spg_get_ir_reciprocal_mesh()
{
  double lattice[3][3] = {{4,0,0},{0,4,0},{0,0,3}};
  double position[][3] =
    {
      {0,0,0},
      {0.5,0.5,0.5},
      {0.3,0.3,0},
      {0.7,0.7,0},
      {0.2,0.8,0.5},
      {0.8,0.2,0.5},
    };
  int types[] = {1,1,2,2,2,2};
  int num_atom = 6;
  int m = 100;
  int mesh[] = { m, m, m };
  int is_shift[] = { 1, 1, 1 };
  int grid_point[m*m*m][3];
  int map[m*m*m];

  printf("*** Example of spg_get_ir_reciprocal_mesh of Rutile structure ***:\n");

  int num_ir = 
    spg_get_ir_reciprocal_mesh( grid_point,
				map,
				mesh,
				is_shift,
				1,
				lattice,
				position,
				types,
				num_atom,
				1e-5 );

  printf("Number of irreducible k-points of Rutile with\n");
  printf("100x100x100 Monkhorst-Pack mesh is %d.\n", num_ir);
}

