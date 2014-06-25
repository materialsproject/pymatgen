#include "spglib.h"
#include <string.h>

void spg_get_multiplicity_(int *size,
			   double lattice[3][3],
			   double position[][3],
			   int types[],
			   int *num_atom,
			   double *symprec);
void spg_get_symmetry_(int *nsym,
		       int rot[][3][3],
		       double trans[][3],
		       int *size,
		       double lattice[3][3],
		       double position[][3],
		       int types[],
		       int *num_atom,
		       double *symprec);
void spg_get_smallest_lattice_(double smallest_lattice[3][3],
			       double lattice[3][3],
			       double *symprec);
void spg_get_international_(int *spacegroup,
			    char symbol[11],
			    double lattice[3][3],
			    double position[][3],
			    int types[],
			    int *num_atom,
			    double *symprec);
void spg_refine_cell_(double lattice[3][3],
		      double position[][3],
		      int types[],
		      int *num_atom,
		      double *symprec);
void spg_get_schoenflies_(int *spacegroup,
			  char symbol[7],
			  double lattice[3][3],
			  double position[][3],
			  int types[],
			  int *num_atom,
			  double *symprec);
void spg_find_primitive_(double lattice[3][3],
			 double position[][3],
			 int types[],
			 int *num_atom,
			 double *symprec);
void spg_get_ir_reciprocal_mesh_(int *num_ir_grid,
				 int grid_point[][3],
				 int map[],
				 int mesh[3],
				 int is_shift[3],
				 int *is_time_reversal,
				 double lattice[3][3],
				 double position[][3],
				 int types[],
				 int *num_atom,
				 double *symprec);





void spg_get_multiplicity_(int *size,
			   double lattice[3][3],
			   double position[][3],
			   int types[],
			   int *num_atom,
			   double *symprec)
{
  *size = spg_get_multiplicity(lattice, position, types, *num_atom, *symprec);
}

void spg_get_symmetry_(int *nsym,
		       int rot[][3][3],
		       double trans[][3],
		       int *size,
		       double lattice[3][3],
		       double position[][3],
		       int types[],
		       int *num_atom,
		       double *symprec)
{
  *nsym = spg_get_symmetry(rot, trans, *size, lattice, position,
			   types, *num_atom, *symprec);
}

void spg_get_smallest_lattice_(double smallest_lattice[3][3],
			       double lattice[3][3],
			       double *symprec)
{
  spg_get_smallest_lattice(smallest_lattice, lattice, *symprec);
}

void spg_get_international_(int *spacegroup,
			    char symbol[11],
			    double lattice[3][3],
			    double position[][3],
			    int types[],
			    int *num_atom,
			    double *symprec)
{
  char symbol_c[11];
  int i, length;

  *spacegroup = spg_get_international(symbol_c, lattice, position, types,
				      *num_atom, *symprec);
  if (*spacegroup > 0) {
    length = strlen(symbol_c);
    strncpy(symbol, symbol_c, length);
  } else {
    length = 0;
  }

  for (i = length; i < 11; i++) {
    symbol[i] = ' ';
  }
}

void spg_refine_cell_(double lattice[3][3],
		      double position[][3],
		      int types[],
		      int *num_atom,
		      double *symprec)
{
  int num_atom_bravais;

  num_atom_bravais = spg_refine_cell(lattice,
				     position,
				     types,
				     *num_atom,
				     *symprec);
  *num_atom = num_atom_bravais;
}


void spg_get_schoenflies_(int *spacegroup,
			  char symbol[7],
			  double lattice[3][3],
			  double position[][3],
			  int types[],
			  int *num_atom,
			  double *symprec)
{
  char symbol_c[7];
  int i, length;
    
  *spacegroup = spg_get_schoenflies(symbol_c, lattice, position, types,
				    *num_atom, *symprec);
  if (*spacegroup > 0) {
    length = strlen(symbol_c);
    strncpy(symbol, symbol_c, length);
  } else {
    length = 0;
  }

  for (i = length; i < 7; i++) {
    symbol[i] = ' ';
  }
}

void spg_find_primitive_(double lattice[3][3],
			 double position[][3],
			 int types[],
			 int *num_atom,
			 double *symprec)
{

  *num_atom = spg_find_primitive(lattice, position, types, *num_atom,
				 *symprec);
}

void spg_get_ir_reciprocal_mesh_(int *num_ir_grid,
				 int grid_point[][3],
				 int map[],
				 int mesh[3],
				 int is_shift[3],
				 int *is_time_reversal,
				 double lattice[3][3],
				 double position[][3],
				 int types[],
				 int *num_atom,
				 double *symprec)
{
  int i;
  *num_ir_grid = spg_get_ir_reciprocal_mesh(grid_point,
					    map,
					    mesh,
					    is_shift,
					    *is_time_reversal,
					    lattice,
					    position,
					    types,
					    *num_atom,
					    *symprec);

  for (i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++) {
    map[i]++;
  }
}
