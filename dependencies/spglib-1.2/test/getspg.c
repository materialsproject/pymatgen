#include <stdio.h>
#include "ruby.h"
#include "spglib.h"

VALUE Getspg = Qnil;
void Init_getspg(void);
VALUE method_getspg(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec);
VALUE method_getptg(VALUE self, VALUE r_rotations);
VALUE method_refine_cell(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec);
VALUE method_get_operations(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec);
VALUE method_get_dataset(VALUE self, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec);

void Init_getspg(void)
{
  Getspg = rb_define_module("Getspg");
  rb_define_method(Getspg, "getspg", method_getspg, 6);
  rb_define_method(Getspg, "getptg", method_getptg, 1);
  rb_define_method(Getspg, "refine_cell", method_refine_cell, 6);
  rb_define_method(Getspg, "get_operations", method_get_operations, 6);
  rb_define_method(Getspg, "get_dataset", method_get_dataset, 5);
}

VALUE method_getspg(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec)
{
  int i, j, size, spgroup;
  double symprec, lattice[3][3];
  VALUE array;
     
  size = NUM2INT(r_size);
     
  double position[size*4][3];
  int types[size*4];
  char symbol[21];
  char output[21];

  spg_set_angle_tolerance();
  symprec = NUM2DBL(r_symprec);

  for (i=0; i<size; i++)
    for (j=0; j<3; j++) {
      position[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_position, i), j));
      types[i] = NUM2DBL(rb_ary_entry(r_types, i));
    }

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      lattice[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_lattice, i), j));

  /* Space group */
  spgroup = spgat_get_international(symbol,
				    lattice,
				    position,
				    types,
				    size,
				    symprec,
				    NUM2DBL(r_angle_symprec));
  sprintf(output, "%s", symbol);

  array = rb_ary_new();
  rb_ary_push(array, rb_str_new2(output));
  rb_ary_push(array, INT2NUM(spgroup));

  return array;
}

VALUE method_getptg(VALUE self, VALUE r_rotations)
{
  int i, j, k, size, ptg_num;
  char symbol[6];
  int trans_mat[3][3];
  VALUE array, matrix, vector;

  size = RARRAY( r_rotations )->len;
  int rotations[size][3][3];

  for ( i = 0; i < size; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      for ( k = 0; k < 3; k++ ) {
	rotations[i][j][k] = NUM2INT( rb_ary_entry( rb_ary_entry( rb_ary_entry( r_rotations, i ), j ), k ) );
      }
    }				      
  }
  
  ptg_num = spg_get_pointgroup( symbol, trans_mat, rotations, size );
  array = rb_ary_new();
  rb_ary_push( array, rb_str_new2( symbol ) );
  rb_ary_push( array, INT2NUM( ptg_num ) );
  matrix = rb_ary_new();
  for ( i = 0; i < 3; i++) {
    vector = rb_ary_new();
    for ( j = 0; j < 3; j++ ) {
      rb_ary_push( vector, INT2NUM( trans_mat[i][j] ) );
    }
    rb_ary_push( matrix, vector);
  }
  rb_ary_push( array, matrix );
  
  return array;
}

VALUE method_refine_cell(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec)
{
  int i, j, size, spgroup, num_atom_bravais ;
  double symprec, angle_tolerance, lattice[3][3];
  VALUE vector, array, r_brv_lattice, r_brv_positions, r_brv_types;
     
  size = NUM2INT(r_size);
     
  double position[size*4][3];
  int types[size*4];
  char symbol[21];
  char output[21];

  angle_tolerance = NUM2DBL(r_angle_symprec);
  symprec = NUM2DBL(r_symprec);

  for (i=0; i<size; i++)
    for (j=0; j<3; j++) {
      position[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_position, i), j));
      types[i] = NUM2DBL(rb_ary_entry(r_types, i));
    }

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      lattice[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_lattice, i), j));

  /* Space group */
  spgroup = spgat_get_international(symbol,
				    lattice,
				    position,
				    types,
				    size,
				    symprec,
				    angle_tolerance);
  sprintf(output, "%s", symbol);

  num_atom_bravais = spgat_refine_cell(lattice,
				       position,
				       types,
				       size,
				       symprec,
				       angle_tolerance);

  r_brv_lattice = rb_ary_new();
  for (i=0; i<3; i++) {
    vector = rb_ary_new();
    for (j=0; j<3; j++) {
      rb_ary_push(vector, rb_float_new(lattice[i][j]));
    }
    rb_ary_push(r_brv_lattice, vector);
  }

  r_brv_positions = rb_ary_new();
  r_brv_types = rb_ary_new();
  for (i=0; i<num_atom_bravais; i++) {
    vector = rb_ary_new();
    for (j=0; j<3; j++) {
      rb_ary_push(vector, rb_float_new(position[i][j]));
    }
    rb_ary_push(r_brv_positions, vector);
    rb_ary_push(r_brv_types, INT2NUM(types[i]));
  }

  array = rb_ary_new();
  rb_ary_push(array, rb_str_new2(output));
  rb_ary_push(array, r_brv_lattice);
  rb_ary_push(array, r_brv_positions);
  rb_ary_push(array, r_brv_types);
  rb_ary_push(array, INT2NUM(spgroup));

  return array;
}

VALUE method_get_operations(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec)
{
  int i, j, k, num_atom, num_sym;
  double symprec, lattice[3][3];
  VALUE matrix, matrix_row, vector, array, r_rotation, r_translation;
     
  num_atom = NUM2INT(r_size);
     
  double position[num_atom][3];
  int types[num_atom];
  int rotation[num_atom*48][3][3];
  double translation[num_atom*48][3];

  symprec = NUM2DBL(r_symprec);

  for (i=0; i<num_atom; i++) {
    for (j=0; j<3; j++) {
      position[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_position, i), j));
      types[i] = NUM2DBL(rb_ary_entry(r_types, i));
    }
  }

  for (i=0; i<3; i++) 
    for (j=0; j<3; j++)
      lattice[i][j] =
	NUM2DBL(rb_ary_entry(rb_ary_entry(r_lattice, i), j));

  num_sym = spgat_get_symmetry(rotation,
			       translation,
			       num_atom*192,
			       lattice,
			       position,
			       types,
			       num_atom,
			       symprec,
			       NUM2DBL(r_angle_symprec));

  
  r_rotation = rb_ary_new();
  r_translation = rb_ary_new();

  for ( i = 0; i < num_sym; i++ ) {
    vector = rb_ary_new();
    matrix = rb_ary_new();
    for ( j = 0; j < 3 ; j++ ) {
      rb_ary_push( vector, rb_float_new( translation[i][j] ) );
      matrix_row = rb_ary_new();
      for ( k = 0; k < 3; k++ ) {
	rb_ary_push( matrix_row, INT2NUM( rotation[i][j][k] ) );
      }
      rb_ary_push( matrix, matrix_row );
    }
    rb_ary_push( r_rotation, matrix );
    rb_ary_push( r_translation, vector );
  }

  array = rb_ary_new();
  rb_ary_push(array, r_rotation);
  rb_ary_push(array, r_translation);

  return array;
}

VALUE method_get_dataset(VALUE self, VALUE r_lattice, VALUE r_position, VALUE r_types, VALUE r_symprec, VALUE r_angle_symprec)
{
  int i, j, k, num_atom;
  double symprec, lattice[3][3];
  SpglibDataset *dataset;
  VALUE mat, vec, row, array, r_tmat, r_oshift, r_rot, r_trans, r_wyckoffs;

  num_atom = RARRAY( r_types )->len;

  double position[num_atom][3];
  int types[num_atom];
  int rotation[num_atom*48][3][3];
  double translation[num_atom*48][3];

  symprec = NUM2DBL(r_symprec);

  for ( i = 0; i < num_atom; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      position[i][j] =
	NUM2DBL( rb_ary_entry( rb_ary_entry( r_position, i ), j ) );
      types[i] = NUM2DBL( rb_ary_entry( r_types, i ) );
    }
  }

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      lattice[i][j] =
	NUM2DBL( rb_ary_entry( rb_ary_entry( r_lattice, i), j ) );
    }
  }

  dataset = spgat_get_dataset(lattice,
			      position,
			      types,
			      num_atom,
			      symprec,
			      NUM2DBL(r_angle_symprec));

  array = rb_ary_new();

  rb_ary_push( array, INT2NUM( dataset->spacegroup_number ) );
  rb_ary_push( array, rb_str_new2( dataset->international_symbol ) );
  rb_ary_push( array, rb_str_new2( dataset->hall_symbol ) );
  
  /* Transformation_matrix */
  r_tmat = rb_ary_new();
  for ( i = 0; i < 3; i++ ) {
    row = rb_ary_new();
    for ( j = 0; j < 3; j++ ) {
      rb_ary_push( row,
		   rb_float_new( dataset->transformation_matrix[i][j] ) );
    }
    rb_ary_push( r_tmat, row );
  }
  rb_ary_push( array, r_tmat );

  /* Origin shift */
  r_oshift = rb_ary_new();
  for ( i = 0; i < 3; i++ ) {
    rb_ary_push( r_oshift, rb_float_new( dataset->origin_shift[i] ) );
  }
  rb_ary_push( array, r_oshift );
 
  /* Rotations, translations */
  r_rot = rb_ary_new();
  r_trans = rb_ary_new();
  for ( i = 0; i < dataset->n_operations; i++ ) {
    mat = rb_ary_new();
    vec = rb_ary_new();
    for ( j = 0; j < 3; j++ ) {
      rb_ary_push( vec, rb_float_new( dataset->translations[i][j] ) );
      row = rb_ary_new();
      for ( k = 0; k < 3; k++ ) {
	rb_ary_push( row, rb_float_new( dataset->rotations[i][j][k] ) );
      }
      rb_ary_push( mat, row );
    }
    rb_ary_push( r_trans, vec );
    rb_ary_push( r_rot, mat );
  }
  rb_ary_push( array, r_rot );
  rb_ary_push( array, r_trans );

  /* Wyckoff letters */
  r_wyckoffs = rb_ary_new();
  for ( i = 0; i < dataset->n_atoms; i++ ) {
    r_wyckoffs = rb_ary_push( r_wyckoffs, INT2NUM( dataset->wyckoffs[i] ) );
  }
  rb_ary_push( array, r_wyckoffs );

  spg_free_dataset( dataset );
  
  return array;
}
