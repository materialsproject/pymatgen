#include "spglib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int test_spg_get_symmetry(void);
static int test_spg_get_symmetry_with_collinear_spin(void);
static int test_spg_get_multiplicity(void);
static int test_spg_find_primitive_BCC(void);
static int test_spg_find_primitive_corundum(void);
static int test_spg_standardize_cell_BCC(void);
static int test_spg_standardize_cell_corundum(void);
static int sub_spg_standardize_cell(double lattice[3][3],
				    double position[][3],
				    int types[],
				    const int num_atom,
				    const double symprec,
				    const int to_primitive,
				    const int no_idealize);
static int test_spg_get_international(void);
static int test_spg_get_schoenflies(void);
static int test_spg_get_spacegroup_type(void);
static int test_spg_get_symmetry_from_database(void);
static int test_spg_refine_cell_BCC(void);
static int test_spg_get_dataset(void);
static int test_spg_get_ir_reciprocal_mesh(void);
static int test_spg_get_stabilized_reciprocal_mesh(void);
static int test_spg_relocate_BZ_grid_address(void);
static void show_spg_dataset(double lattice[3][3],
			     const double origin_shift[3],
			     double position[][3],
			     const int num_atom,
			     const int types[]);
static void show_cell(double lattice[3][3],
		      double position[][3],
		      const int types[],
		      const int num_atom);

int main(void)
{

  int (*funcs[])(void) = {
    test_spg_find_primitive_BCC,
    test_spg_find_primitive_corundum,
    test_spg_standardize_cell_BCC,
    test_spg_standardize_cell_corundum,
    test_spg_get_multiplicity,
    test_spg_get_symmetry,
    test_spg_get_symmetry_with_collinear_spin,
    test_spg_get_international,
    test_spg_get_schoenflies,
    test_spg_get_spacegroup_type,
    test_spg_get_symmetry_from_database,
    test_spg_refine_cell_BCC,
    test_spg_get_dataset,
    test_spg_get_ir_reciprocal_mesh,
    test_spg_get_stabilized_reciprocal_mesh,
    test_spg_relocate_BZ_grid_address,
    NULL};

  int i, result;

  for (i = 0; i < 30; i++) {
    if (*funcs[i] == NULL) {
      break;
    } else {
      result = (funcs[i])();
    }
    if (result) {
      return 1;
    }
  }

  return 0;
}

static int test_spg_find_primitive_BCC(void)
{
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
  double position[][3] = {
    {0, 0, 0},
    {0.5, 0.5, 0.5}
  };
  int types[] = {1, 1};
  int i, num_atom = 2, num_primitive_atom;
  double symprec = 1e-5;
  
  /* lattice, position, and types are overwirtten. */
  printf("*** spg_find_primitive (BCC unitcell --> primitive) ***:\n");
  num_primitive_atom = spg_find_primitive(lattice, position, types, num_atom, symprec);
  if (num_primitive_atom == 0) {
    printf("Primitive cell was not found.\n");
    return 1;
  } else { 
    show_cell(lattice, position, types, num_primitive_atom);
    return 0;
  }
}

static int test_spg_find_primitive_corundum(void)
{
  double lattice[3][3] = {{4.8076344022756095, -2.4038172011378047, 0},
			  {0, 4.1635335244786962, 0},
			  {0, 0, 13.1172699198127543}};
  double position[][3] = {
    {0.0000000000000000, 0.0000000000000000, 0.3521850942289043},
    {0.6666666666666643, 0.3333333333333357, 0.6855184275622400},
    {0.3333333333333357, 0.6666666666666643, 0.0188517608955686},
    {0.0000000000000000, 0.0000000000000000, 0.6478149057711028},
    {0.6666666666666643, 0.3333333333333357, 0.9811482391044314},
    {0.3333333333333357, 0.6666666666666643, 0.3144815724377600},
    {0.0000000000000000, 0.0000000000000000, 0.1478149057710957},
    {0.6666666666666643, 0.3333333333333357, 0.4811482391044314},
    {0.3333333333333357, 0.6666666666666643, 0.8144815724377600},
    {0.0000000000000000, 0.0000000000000000, 0.8521850942288972},
    {0.6666666666666643, 0.3333333333333357, 0.1855184275622400},
    {0.3333333333333357, 0.6666666666666643, 0.5188517608955686},
    {0.3061673906454899, 0.0000000000000000, 0.2500000000000000},
    {0.9728340573121541, 0.3333333333333357, 0.5833333333333357},
    {0.6395007239788255, 0.6666666666666643, 0.9166666666666643},
    {0.6938326093545102, 0.0000000000000000, 0.7500000000000000},
    {0.3604992760211744, 0.3333333333333357, 0.0833333333333357},
    {0.0271659426878458, 0.6666666666666643, 0.4166666666666643},
    {0.0000000000000000, 0.3061673906454899, 0.2500000000000000},
    {0.6666666666666643, 0.6395007239788255, 0.5833333333333357},
    {0.3333333333333357, 0.9728340573121541, 0.9166666666666643},
    {0.0000000000000000, 0.6938326093545102, 0.7500000000000000},
    {0.6666666666666643, 0.0271659426878458, 0.0833333333333357},
    {0.3333333333333357, 0.3604992760211744, 0.4166666666666643},
    {0.6938326093545102, 0.6938326093545102, 0.2500000000000000},
    {0.3604992760211744, 0.0271659426878458, 0.5833333333333357},
    {0.0271659426878458, 0.3604992760211744, 0.9166666666666643},
    {0.3061673906454899, 0.3061673906454899, 0.7500000000000000},
    {0.9728340573121541, 0.6395007239788255, 0.0833333333333357},
    {0.6395007239788255, 0.9728340573121541, 0.4166666666666643},
  };
  int types[30];
  int i, num_atom = 30, num_primitive_atom;
  double symprec = 1e-5;
  
  for (i = 0; i < 12; i++) {
    types[i] = 1;
  }
  for (i = 12; i < 30; i++) {
    types[i] = 2;
  }

  /* lattice, position, and types are overwirtten. */
  printf("*** spg_find_primitive (Corundum) ***:\n");
  num_primitive_atom = spg_find_primitive(lattice, position, types, num_atom, symprec);
  if (num_primitive_atom == 0) {
    printf("Primitive cell was not found.\n");
    return 1;
  } else { 
    show_cell(lattice, position, types, num_primitive_atom);
    return 0;
  }
}

static int test_spg_refine_cell_BCC(void)
{
  double lattice[3][3] = {{0, 2, 2}, {2, 0, 2}, {2, 2, 0}};

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
  printf("*** spg_refine_cell ***:\n");
  num_atom_bravais = spg_refine_cell( lattice,
				      position,
				      types,
				      num_atom,
				      symprec );
  show_cell(lattice, position, types, num_atom_bravais);
  return 0;
}

static int test_spg_standardize_cell_BCC(void)
{
  double lattice[3][3] = {{3.97, 0, 0}, {0, 4.03, 0}, {0, 0, 4.0}};
  double position[][3] = {
    {0.002, 0, 0},
    {0.5, 0.5001, 0.5}
  };
  int types[] = {1, 1};
  int i, j, k, num_atom = 2, num_primitive_atom;
  double symprec = 1e-1;
  
  /* lattice, position, and types are overwirtten. */
  printf("*** spg_standardize_cell (BCC unitcell --> primitive) ***:\n");
  printf("------------------------------------------------------\n");
  for (j = 0; j < 2; j++) {
    for (k = 0; k < 2; k++) {
      sub_spg_standardize_cell(lattice,
			       position,
			       types,
			       num_atom,
			       symprec,
			       j,
			       k);
      printf("------------------------------------------------------\n");
    }
  }

  return 0;
}

static int test_spg_standardize_cell_corundum(void)
{
  double lattice[3][3] = {{4.8076344022756095, -2.4038172011378047, 0},
			  {0, 4.1635335244786962, 0},
			  {0, 0, 13.1172699198127543}};
  double position[][3] = {
    {0.0000000000000000, 0.0000000000000000, 0.3521850942289043},
    {0.6666666666666643, 0.3333333333333357, 0.6855184275622400},
    {0.3333333333333357, 0.6666666666666643, 0.0188517608955686},
    {0.0000000000000000, 0.0000000000000000, 0.6478149057711028},
    {0.6666666666666643, 0.3333333333333357, 0.9811482391044314},
    {0.3333333333333357, 0.6666666666666643, 0.3144815724377600},
    {0.0000000000000000, 0.0000000000000000, 0.1478149057710957},
    {0.6666666666666643, 0.3333333333333357, 0.4811482391044314},
    {0.3333333333333357, 0.6666666666666643, 0.8144815724377600},
    {0.0000000000000000, 0.0000000000000000, 0.8521850942288972},
    {0.6666666666666643, 0.3333333333333357, 0.1855184275622400},
    {0.3333333333333357, 0.6666666666666643, 0.5188517608955686},
    {0.3061673906454899, 0.0000000000000000, 0.2500000000000000},
    {0.9728340573121541, 0.3333333333333357, 0.5833333333333357},
    {0.6395007239788255, 0.6666666666666643, 0.9166666666666643},
    {0.6938326093545102, 0.0000000000000000, 0.7500000000000000},
    {0.3604992760211744, 0.3333333333333357, 0.0833333333333357},
    {0.0271659426878458, 0.6666666666666643, 0.4166666666666643},
    {0.0000000000000000, 0.3061673906454899, 0.2500000000000000},
    {0.6666666666666643, 0.6395007239788255, 0.5833333333333357},
    {0.3333333333333357, 0.9728340573121541, 0.9166666666666643},
    {0.0000000000000000, 0.6938326093545102, 0.7500000000000000},
    {0.6666666666666643, 0.0271659426878458, 0.0833333333333357},
    {0.3333333333333357, 0.3604992760211744, 0.4166666666666643},
    {0.6938326093545102, 0.6938326093545102, 0.2500000000000000},
    {0.3604992760211744, 0.0271659426878458, 0.5833333333333357},
    {0.0271659426878458, 0.3604992760211744, 0.9166666666666643},
    {0.3061673906454899, 0.3061673906454899, 0.7500000000000000},
    {0.9728340573121541, 0.6395007239788255, 0.0833333333333357},
    {0.6395007239788255, 0.9728340573121541, 0.4166666666666643},
  };
  int types[30];
  int i, j, k, num_atom = 30;
  double symprec = 1e-5;
  
  for (i = 0; i < 12; i++) {
    types[i] = 1;
  }
  for (i = 12; i < 30; i++) {
    types[i] = 2;
  }
  
  /* lattice, position, and types are overwirtten. */
  printf("*** spg_standardize_cell (Corundum) ***:\n");
  printf("------------------------------------------------------\n");
  for (j = 0; j < 2; j++) {
    for (k = 0; k < 2; k++) {
      sub_spg_standardize_cell(lattice,
			       position,
			       types,
			       num_atom,
			       symprec,
			       j,
			       k);
      printf("------------------------------------------------------\n");
    }
  }

  return 0;
}

static int sub_spg_standardize_cell(double lattice[3][3],
				    double position[][3],
				    int types[],
				    const int num_atom,
				    const double symprec,
				    const int to_primitive,
				    const int no_idealize)
{
  int i, num_primitive_atom;
  double lat[3][3], pos[num_atom][3];
  int typ[num_atom];

  for (i = 0; i < 3; i++) {
    lat[i][0] = lattice[i][0];
    lat[i][1] = lattice[i][1];
    lat[i][2] = lattice[i][2];
  }

  for (i = 0; i < num_atom; i++) {
    pos[i][0] = position[i][0];
    pos[i][1] = position[i][1];
    pos[i][2] = position[i][2];
    typ[i] = types[i];
  }
  
  /* lattice, position, and types are overwirtten. */
  num_primitive_atom = spg_standardize_cell(lat,
					    pos,
					    typ,
					    num_atom,
					    to_primitive,
					    no_idealize,
					    symprec);
  printf("VASP POSCAR format: ");
  if (to_primitive == 0) {
    printf("to_primitive=0 and ");
  } else {
    printf("to_primitive=1 and ");
  }

  if (no_idealize == 0) {
    printf("no_idealize=0\n");
  } else {
    printf("no_idealize=1\n");
  }
  printf("1.0\n");
  for (i = 0; i < 3; i++) {
    printf("%f %f %f\n", lat[0][i], lat[1][i], lat[2][i]);
  }
  printf("%d\n", num_primitive_atom);
  printf("Direct\n");
  for (i = 0; i < num_primitive_atom; i++) {
    printf("%f %f %f\n", pos[i][0], pos[i][1], pos[i][2]);
  }

  return 0;
}

static int test_spg_get_international(void)
{
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.5},
      {0.3, 0.3, 0},
      {0.7, 0.7, 0},
      {0.2, 0.8, 0.5},
      {0.8, 0.2, 0.5},
    };
  int types[] = {1, 1, 2, 2, 2, 2};
  int num_spg, num_atom = 6;
  char symbol[21];
  
  num_spg = spg_get_international(symbol, lattice, position, types, num_atom, 1e-5);
  printf("*** spg_get_international ***:\n");
  if ( num_spg > 0 ) {
    printf("%s (%d)\n", symbol, num_spg);
    return 0;
  } else {
    printf("Space group could not be found.\n");
    return 1;
  }
}

static int test_spg_get_schoenflies(void)
{
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.5},
      {0.3, 0.3, 0},
      {0.7, 0.7, 0},
      {0.2, 0.8, 0.5},
      {0.8, 0.2, 0.5},
    };
  int types[] = {1, 1, 2, 2, 2, 2};
  int num_atom = 6;
  char symbol[7];
  
  spg_get_schoenflies(symbol, lattice, position, types, num_atom, 1e-5);
  printf("*** spg_get_schoenflies ***:\n");
  printf("Schoenflies: %s\n", symbol);
  
  return 0;
}

static int test_spg_get_spacegroup_type(void)
{
  SpglibSpacegroupType spgtype;
  spgtype = spg_get_spacegroup_type(446);

  printf("*** spg_get_spacegroup_type ***:\n");
  printf("Number:        %d\n", spgtype.number);
  printf("Schoenflies:   %s\n", spgtype.schoenflies);
  printf("International: %s\n", spgtype.international);
  printf("International: %s\n", spgtype.international_full);
  printf("International: %s\n", spgtype.international_short);
  printf("Hall symbol:   %s\n", spgtype.hall_symbol);
 
  return 0;
}

static int test_spg_get_symmetry_from_database(void)
{
  int rotations[192][3][3];
  double translations[192][3];
  int max_size, i, j, size;
  
  size = spg_get_symmetry_from_database(rotations,
					translations,
					460);
  
  printf("*** spg_get_symmetry_from_database ***:\n");
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
      printf("%2d %2d %2d\n",
	     rotations[i][j][0], rotations[i][j][1], rotations[i][j][2]);
    }
    printf("%f %f %f\n",
	   translations[i][0], translations[i][1], translations[i][2]);
  }

  return 0;
}

static int test_spg_get_multiplicity(void)
{
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
  double position[][3] = {
    {0, 0, 0},
    {0.5, 0.5, 0.5}
  };
  int types[] = {1, 2};
  int num_atom = 2;
  int size;

  size = spg_get_multiplicity(lattice, position, types, num_atom, 1e-5);

  printf("*** spg_get_multiplicity ***:\n");
  printf("Number of symmetry operations: %d\n", size);

  return 0;
}


static int test_spg_get_symmetry(void)
{
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.25},
      {0.3, 0.3, 0},
      {0.7, 0.7, 0},
      {0.2, 0.8, 0.25},
      {0.8, 0.2, 0.25},
      {0, 0, 0.5},
      {0.5, 0.5, 0.75},
      {0.3, 0.3, 0.5},
      {0.7, 0.7, 0.5},
      {0.2, 0.8, 0.75},
      {0.8, 0.2, 0.75}
    };
  int types[] = {1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2};
  int num_atom = 12;
  int max_size = 50;
  int i, j, size;
  int rotation[max_size][3][3];
  double translation[max_size][3];

  double origin_shift[3] = {0.1, 0.1, 0};
  for (i = 0; i < num_atom; i++) {
    for (j = 0; j < 3; j++) {
      position[i][j] += origin_shift[j];
    }
  }

  printf("*** spg_get_symmetry (Rutile two unit cells) ***:\n");
  size = spg_get_symmetry(rotation,
			  translation,
			  max_size,
			  lattice,
			  position,
			  types,
			  num_atom,
			  1e-5);
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++)
      printf("%2d %2d %2d\n", rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    printf("%f %f %f\n",
	   translation[i][0], translation[i][1], translation[i][2]);
  }

  return 0;
}

static int test_spg_get_symmetry_with_collinear_spin(void) {
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 4}};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.5}
    };	
  int types[] = {1, 1};
  int equivalent_atoms[2];
  double spins[2];
  int num_atom = 2;
  int max_size = 300;
  int i, j, size;
  int rotation[max_size][3][3];
  double translation[max_size][3];

  printf("*** spg_get_symmetry_with_spin (BCC ferro) ***:\n");
  spins[0] = 1;
  spins[1] = 1;
  size = spg_get_symmetry_with_collinear_spin(rotation,
					      translation,
					      equivalent_atoms,
					      max_size,
					      lattice,
					      position,
					      types,
					      spins,
					      num_atom,
					      1e-5);
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
      printf("%2d %2d %2d\n",
	     rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    }
    printf("%f %f %f\n", translation[i][0], translation[i][1],
           translation[i][2]);
  }

  printf("*** Example of spg_get_symmetry_with_spin (BCC antiferro) ***:\n");
  spins[0] = 1;
  spins[1] = -1;
  size = spg_get_symmetry_with_collinear_spin(rotation,
					      translation,
					      equivalent_atoms,
					      max_size,
					      lattice,
					      position,
					      types,
					      spins,
					      num_atom,
					      1e-5);
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
      printf("%2d %2d %2d\n",
	     rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    }
    printf("%f %f %f\n", translation[i][0], translation[i][1],
           translation[i][2]);
  }

  printf("*** spg_get_symmetry_with_spin (BCC broken spin) ***:\n");
  spins[0] = 1;
  spins[1] = 2;
  size = spg_get_symmetry_with_collinear_spin(rotation,
					      translation,
					      equivalent_atoms,
					      max_size,
					      lattice,
					      position,
					      types,
					      spins,
					      num_atom,
					      1e-5);
  for (i = 0; i < size; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
      printf("%2d %2d %2d\n",
	     rotation[i][j][0], rotation[i][j][1], rotation[i][j][2]);
    }
    printf("%f %f %f\n", translation[i][0], translation[i][1],
           translation[i][2]);
  }

  return 0;
}


static int test_spg_get_dataset(void)
{
  printf("*** spg_get_dataset (Rutile two unit cells) ***:\n");
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double origin_shift[3] = {0.1, 0.1, 0};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.25},
      {0.3, 0.3, 0},
      {0.7, 0.7, 0},
      {0.2, 0.8, 0.25},
      {0.8, 0.2, 0.25},
      {0, 0, 0.5},
      {0.5, 0.5, 0.75},
      {0.3, 0.3, 0.5},
      {0.7, 0.7, 0.5},
      {0.2, 0.8, 0.75},
      {0.8, 0.2, 0.75}
    };
  int types[] = {1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2};
  int num_atom = 12;

  show_spg_dataset(lattice, origin_shift, position, num_atom, types);

  double lattice_2[3][3] = {{3.7332982433264039, -1.8666491216632011, 0},
			    {0, 3.2331311186244847, 0},
			    {0, 0, 6.0979971306362799}};
  double origin_shift_2[3] = {0.1, 0.1, 0};
  double position_2[][3] =
    {
      {0, 0, 0},
      {1.0 / 3, 2.0 / 3, 0.4126},
      {1.0 / 3, 2.0 / 3, 0.776},
      {2.0 / 3, 1.0 / 3, 0.2542},
    };
  int types_2[] = {1, 2, 3, 3};
  int num_atom_2 = 4;

  show_spg_dataset(lattice_2, origin_shift_2, position_2, num_atom_2, types_2);

  return 0;
}

static int test_spg_get_ir_reciprocal_mesh(void)
{
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.5},
      {0.3, 0.3, 0},
      {0.7, 0.7, 0},
      {0.2, 0.8, 0.5},
      {0.8, 0.2, 0.5},
    };
  int types[] = {1, 1, 2, 2, 2, 2};
  int num_atom = 6;
  int m = 40;
  int mesh[] = {m, m, m};
  int is_shift[] = {1, 1, 1};
  int grid_address[m * m * m][3];
  int grid_mapping_table[m * m * m];

  printf("*** spg_get_ir_reciprocal_mesh of Rutile structure ***:\n");

  int num_ir = spg_get_ir_reciprocal_mesh(grid_address,
					  grid_mapping_table,
					  mesh,
					  is_shift,
					  1,
					  lattice,
					  position,
					  types,
					  num_atom,
					  1e-5);

  printf("Number of irreducible k-points of Rutile with\n");
  printf("40x40x40 Monkhorst-Pack mesh is %d (4200).\n", num_ir);

  return 0;
}

static int test_spg_get_stabilized_reciprocal_mesh(void)
{
  SpglibDataset *dataset;
  double lattice[3][3] = {{4, 0, 0}, {0, 4, 0}, {0, 0, 3}};
  double position[][3] =
    {
      {0, 0, 0},
      {0.5, 0.5, 0.5},
      {0.3, 0.3, 0},
      {0.7, 0.7, 0},
      {0.2, 0.8, 0.5},
      {0.8, 0.2, 0.5},
    };
  int types[] = {1, 1, 2, 2, 2, 2};
  int num_atom = 6;

  dataset = spg_get_dataset(lattice,
			    position,
			    types,
			    num_atom,
			    1e-5);

  
  int m = 40;
  int mesh[] = {m, m, m};
  int is_shift[] = {1, 1, 1};
  int grid_address[m * m * m][3];
  int grid_mapping_table[m * m * m];
  double q[] = {0, 0.5, 0.5};

  printf("*** spg_get_stabilized_reciprocal_mesh of Rutile structure ***:\n");

  int num_ir = spg_get_stabilized_reciprocal_mesh(grid_address,
						  grid_mapping_table,
						  mesh,
						  is_shift,
						  1,
						  dataset->n_operations,
						  dataset->rotations,
						  1,
						  (double(*)[3])q);
  spg_free_dataset(dataset);

  printf("Number of irreducible k-points stabilized by q=(0, 1/2, 1/2) of Rutile with\n");
  printf("40x40x40 Monkhorst-Pack mesh is %d (8000).\n", num_ir);

  return 0;
}

static int test_spg_relocate_BZ_grid_address(void)
{
  double rec_lattice[3][3] = {{-0.17573761,  0.17573761,  0.17573761},
			      { 0.17573761, -0.17573761,  0.17573761},
			      { 0.17573761,  0.17573761, -0.17573761}};
  int rotations[][3][3] = {{{1, 0, 0},
			    {0, 1, 0},
			    {0, 0, 1}}};

  
  int m = 40;
  int mesh[] = {m, m, m};
  int is_shift[] = {0, 0, 0};
  int bz_grid_address[(m + 1) * (m + 1) * (m + 1)][3];
  int bz_map[m * m * m * 8];
  int grid_address[m * m * m][3];
  int grid_mapping_table[m * m * m];
  double q[] = {0, 0, 0};

  int num_ir = spg_get_stabilized_reciprocal_mesh(grid_address,
						  grid_mapping_table,
						  mesh,
						  is_shift,
						  1,
						  1,
						  rotations,
						  1,
						  (double(*)[3])q);
  

  printf("*** spg_relocate_BZ_grid_address of NaCl structure ***:\n");

  int num_q = spg_relocate_BZ_grid_address(bz_grid_address,
					   bz_map,
					   grid_address,
					   mesh,
					   rec_lattice,
					   is_shift);

  printf("Number of k-points of NaCl Brillouin zone\n");
  printf("with Gamma-centered 40x40x40 Monkhorst-Pack mesh is %d (65861).\n", num_q);

  return 0;
}

static void show_spg_dataset(double lattice[3][3],
			     const double origin_shift[3],
			     double position[][3],
			     const int num_atom,
			     const int types[])
{
  SpglibDataset *dataset;
  char ptsymbol[6];
  int pt_trans_mat[3][3];

  int i, j, size;
  const char *wl = "abcdefghijklmnopqrstuvwxyz";

  for ( i = 0; i < num_atom; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      position[i][j] += origin_shift[j];
    }
  }

  dataset = spg_get_dataset(lattice,
			    position,
			    types,
			    num_atom,
			    1e-5);
  
  printf("International: %s (%d)\n", dataset->international_symbol, dataset->spacegroup_number );
  printf("Hall symbol:   %s\n", dataset->hall_symbol );
  spg_get_pointgroup(ptsymbol,
		     pt_trans_mat,
		     dataset->rotations,
		     dataset->n_operations);
  printf("Point group:   %s\n", ptsymbol);
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
  for (i = 0; i < dataset->n_atoms; i++) {
    printf("%d ", dataset->equivalent_atoms[i]);
  }
  printf("\n");
  
  for (i = 0; i < dataset->n_operations; i++) {
    printf("--- %d ---\n", i + 1);
    for (j = 0; j < 3; j++) {
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

  spg_free_dataset(dataset);
}

static void show_cell(double lattice[3][3],
		      double position[][3],
		      const int types[],
		      const int num_atom)
{
  int i;

  printf("Lattice parameter:\n");
  for (i = 0; i < 3; i++) {
    printf("%f %f %f\n", lattice[0][i], lattice[1][i], lattice[2][i]);
  }
  printf("Atomic positions:\n");
  for (i = 0; i < num_atom; i++) {
    printf("%d: %f %f %f\n",
	   types[i], position[i][0], position[i][1], position[i][2]);
  }
}
