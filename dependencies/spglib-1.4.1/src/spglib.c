/* spglib.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spglib.h"
#include "refinement.h"
#include "cell.h"
#include "lattice.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "primitive.h"
#include "spacegroup.h"
#include "symmetry.h"
#include "kpoint.h"
#include "spin.h"

static SpglibDataset * get_dataset(SPGCONST double lattice[3][3],
				   SPGCONST double position[][3],
				   const int types[],
				   const int num_atom,
				   const double symprec);

static int get_symmetry(int rotation[][3][3],
			double translation[][3],
			const int max_size,
			SPGCONST double lattice[3][3],
			SPGCONST double position[][3],
			const int types[],
			const int num_atom,
			const double symprec);

static int get_symmetry_with_collinear_spin(int rotation[][3][3],
					    double translation[][3],
					    const int max_size,
					    SPGCONST double lattice[3][3],
					    SPGCONST double position[][3],
					    const int types[],
					    const double spins[],
					    const int num_atom,
					    const double symprec);

static int get_multiplicity(SPGCONST double lattice[3][3],
			    SPGCONST double position[][3],
			    const int types[],
			    const int num_atom,
			    const double symprec);

static int find_primitive(double lattice[3][3],
			  double position[][3],
			  int types[],
			  const int num_atom,
			  const double symprec);

static int get_international(char symbol[11],
			     SPGCONST double lattice[3][3],
			     SPGCONST double position[][3],
			     const int types[],
			     const int num_atom,
			     const double symprec);

static int get_schoenflies(char symbol[10],
			   SPGCONST double lattice[3][3],
			   SPGCONST double position[][3],
			   const int types[], const int num_atom,
			   const double symprec);

static int refine_cell(double lattice[3][3],
		       double position[][3],
		       int types[],
		       const int num_atom,
		       const double symprec);

static int get_ir_kpoints(int map[],
			  SPGCONST double kpoints[][3],
			  const int num_kpoint,
			  SPGCONST double lattice[3][3],
			  SPGCONST double position[][3],
			  const int types[],
			  const int num_atom,
			  const int is_time_reversal,
			  const double symprec);

static int get_ir_reciprocal_mesh(int grid_point[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  const int is_time_reversal,
				  SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec);

static int get_stabilized_reciprocal_mesh(int grid_point[][3],
					  int map[],
					  const int mesh[3],
					  const int is_shift[3],
					  const int is_time_reversal,
					  const int num_rot,
					  SPGCONST int rotations[][3][3],
					  const int num_q,
					  SPGCONST double qpoints[][3]);

static SpglibTriplets * get_triplets_reciprocal_mesh(const int mesh[3],
						     const int is_time_reversal,
						     const int num_rot,
						     SPGCONST int rotations[][3][3]);

static int get_triplets_reciprocal_mesh_at_q(int weights[],
					     int grid_points[][3],
					     int third_q[],
					     const int grid_point,
					     const int mesh[3],
					     const int is_time_reversal,
					     const int num_rot,
					     SPGCONST int rotations[][3][3]);

static int extract_triplets_reciprocal_mesh_at_q(int triplets_at_q[][3],
						 int weight_triplets_at_q[],
						 const int fixed_grid_number,
						 const int num_triplets,
						 SPGCONST int triplets[][3],
						 const int mesh[3],
						 const int is_time_reversal,
						 const int num_rot,
						 SPGCONST int rotations[][3][3]);


SpglibDataset * spg_get_dataset(SPGCONST double lattice[3][3],
				SPGCONST double position[][3],
				const int types[],
				const int num_atom,
				const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_dataset(lattice,
		     position,
		     types,
		     num_atom,
		     symprec);
}

SpglibDataset * spgat_get_dataset(SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec,
				  const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_dataset(lattice,
		     position,
		     types,
		     num_atom,
		     symprec);
}

void spg_free_dataset(SpglibDataset *dataset)
{
  if (dataset->n_operations > 0) {
    free(dataset->rotations);
    dataset->rotations = NULL;
    free(dataset->translations);
    dataset->translations = NULL;
  }
  
  if (! (dataset->wyckoffs == NULL)) {
    free(dataset->wyckoffs);
    dataset->wyckoffs = NULL;
  }
  
  if (! (dataset->equivalent_atoms == NULL)) {
    free(dataset->equivalent_atoms);
    dataset->equivalent_atoms = NULL;
  }

  free(dataset);
  dataset = NULL;
}

int spg_get_symmetry(int rotation[][3][3],
		     double translation[][3],
		     const int max_size,
		     SPGCONST double lattice[3][3],
		     SPGCONST double position[][3],
		     const int types[],
		     const int num_atom,
		     const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_symmetry(rotation,
		      translation,
		      max_size,
		      lattice,
		      position,
		      types,
		      num_atom,
		      symprec);
}

int spgat_get_symmetry(int rotation[][3][3],
		       double translation[][3],
		       const int max_size,
		       SPGCONST double lattice[3][3],
		       SPGCONST double position[][3],
		       const int types[],
		       const int num_atom,
		       const double symprec,
		       const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_symmetry(rotation,
		      translation,
		      max_size,
		      lattice,
		      position,
		      types,
		      num_atom,
		      symprec);
}

int spg_get_symmetry_with_collinear_spin(int rotation[][3][3],
					 double translation[][3],
					 const int max_size,
					 SPGCONST double lattice[3][3],
					 SPGCONST double position[][3],
					 const int types[],
					 const double spins[],
					 const int num_atom,
					 const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_symmetry_with_collinear_spin(rotation,
					  translation,
					  max_size,
					  lattice,
					  position,
					  types,
					  spins,
					  num_atom,
					  symprec);
}

int spgat_get_symmetry_with_collinear_spin(int rotation[][3][3],
					   double translation[][3],
					   const int max_size,
					   SPGCONST double lattice[3][3],
					   SPGCONST double position[][3],
					   const int types[],
					   const double spins[],
					   const int num_atom,
					   const double symprec,
					   const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_symmetry_with_collinear_spin(rotation,
					  translation,
					  max_size,
					  lattice,
					  position,
					  types,
					  spins,
					  num_atom,
					  symprec);
}

int spg_get_multiplicity(SPGCONST double lattice[3][3],
			 SPGCONST double position[][3],
			 const int types[],
			 const int num_atom,
			 const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_multiplicity(lattice,
			  position,
			  types,
			  num_atom,
			  symprec);
}

int spgat_get_multiplicity(SPGCONST double lattice[3][3],
			   SPGCONST double position[][3],
			   const int types[],
			   const int num_atom,
			   const double symprec,
			   const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_multiplicity(lattice,
			  position,
			  types,
			  num_atom,
			  symprec);
}

int spg_get_smallest_lattice(double smallest_lattice[3][3],
			     SPGCONST double lattice[3][3],
			     const double symprec)
{
  return lat_smallest_lattice_vector(smallest_lattice, lattice, symprec);
}

int spg_find_primitive(double lattice[3][3],
		       double position[][3],
		       int types[],
		       const int num_atom,
		       const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return find_primitive(lattice,
			position,
			types,
			num_atom,
			symprec);
}

int spgat_find_primitive(double lattice[3][3],
			 double position[][3],
			 int types[],
			 const int num_atom,
			 const double symprec,
			 const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return find_primitive(lattice,
			position,
			types,
			num_atom,
			symprec);
}

int spg_get_international(char symbol[11],
			  SPGCONST double lattice[3][3],
			  SPGCONST double position[][3],
			  const int types[],
			  const int num_atom,
			  const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_international(symbol,
			   lattice,
			   position,
			   types,
			   num_atom,
			   symprec);
}

int spgat_get_international(char symbol[11],
			    SPGCONST double lattice[3][3],
			    SPGCONST double position[][3],
			    const int types[],
			    const int num_atom,
			    const double symprec,
			    const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_international(symbol,
			   lattice,
			   position,
			   types,
			   num_atom,
			   symprec);
}

int spg_get_schoenflies(char symbol[10],
			SPGCONST double lattice[3][3],
			SPGCONST double position[][3],
			const int types[],
			const int num_atom,
			const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_schoenflies(symbol,
			 lattice,
			 position,
			 types,
			 num_atom,
			 symprec);
}

int spgat_get_schoenflies(char symbol[10],
			  SPGCONST double lattice[3][3],
			  SPGCONST double position[][3],
			  const int types[],
			  const int num_atom,
			  const double symprec,
			  const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_schoenflies(symbol,
			 lattice,
			 position,
			 types,
			 num_atom,
			 symprec);
}

int spg_get_pointgroup(char symbol[6],
		       int trans_mat[3][3],
		       SPGCONST int rotations[][3][3],
		       const int num_rotations)
{
  int ptg_num;
  double tmp_trans_mat[3][3];
  Pointgroup ptgroup;

  ptg_num = ptg_get_pointgroup_number_by_rotations(rotations,
						   num_rotations);
  ptgroup = ptg_get_pointgroup(ptg_num);
  strcpy(symbol, ptgroup.symbol);
  ptg_get_transformation_matrix(tmp_trans_mat,
				rotations,
				num_rotations);
  mat_cast_matrix_3d_to_3i(trans_mat, tmp_trans_mat);
  return ptg_num + 1;
}

int spg_refine_cell(double lattice[3][3],
		    double position[][3],
		    int types[],
		    const int num_atom,
		    const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return refine_cell(lattice,
		     position,
		     types,
		     num_atom,
		     symprec);
}

int spgat_refine_cell(double lattice[3][3],
		      double position[][3],
		      int types[],
		      const int num_atom,
		      const double symprec,
		      const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return refine_cell(lattice,
		     position,
		     types,
		     num_atom,
		     symprec);
}

int spg_get_ir_kpoints(int map[],
		       SPGCONST double kpoints[][3],
		       const int num_kpoint,
		       SPGCONST double lattice[3][3],
		       SPGCONST double position[][3],
		       const int types[],
		       const int num_atom,
		       const int is_time_reversal,
		       const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_ir_kpoints(map,
			kpoints,
			num_kpoint,
			lattice,
			position,
			types,
			num_atom,
			is_time_reversal,
			symprec);
}

int spg_get_ir_reciprocal_mesh(int grid_point[][3],
			       int map[],
			       const int mesh[3],
			       const int is_shift[3],
			       const int is_time_reversal,
			       SPGCONST double lattice[3][3],
			       SPGCONST double position[][3],
			       const int types[],
			       const int num_atom,
			       const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_ir_reciprocal_mesh(grid_point,
				map,
				mesh,
				is_shift,
				is_time_reversal,
				lattice,
				position,
				types,
				num_atom,
				symprec);
}

int spg_get_stabilized_reciprocal_mesh(int grid_point[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const int num_rot,
				       SPGCONST int rotations[][3][3],
				       const int num_q,
				       SPGCONST double qpoints[][3])
{
  return get_stabilized_reciprocal_mesh(grid_point,
					map,
					mesh,
					is_shift,
					is_time_reversal,
					num_rot,
					rotations,
					num_q,
					qpoints);
}

SpglibTriplets * spg_get_triplets_reciprocal_mesh(const int mesh[3],
						  const int is_time_reversal,
						  const int num_rot,
						  SPGCONST int rotations[][3][3])
{
  return get_triplets_reciprocal_mesh(mesh,
				      is_time_reversal,
				      num_rot,
				      rotations);
}

void spg_free_triplets(SpglibTriplets * spg_triplets)
{
  free(spg_triplets->triplets);
  spg_triplets->triplets = NULL;
  free(spg_triplets->weights);
  spg_triplets->weights = NULL;
  free(spg_triplets);
  free(spg_triplets->mesh_points);
  spg_triplets->mesh_points = NULL;
  spg_triplets = NULL;
}

int spg_get_triplets_reciprocal_mesh_at_q(int weights[],
					  int grid_points[][3],
					  int third_q[],
					  const int grid_point,
					  const int mesh[3],
					  const int is_time_reversal,
					  const int num_rot,
					  SPGCONST int rotations[][3][3])
{
  return get_triplets_reciprocal_mesh_at_q(weights,
					   grid_points,
					   third_q,
					   grid_point,
					   mesh,
					   is_time_reversal,
					   num_rot,
					   rotations);
}

int spg_extract_triplets_reciprocal_mesh_at_q(int triplets_at_q[][3],
					      int weight_triplets_at_q[],
					      const int fixed_grid_number,
					      const int num_triplets,
					      SPGCONST int triplets[][3],
					      const int mesh[3],
					      const int is_time_reversal,
					      const int num_rot,
					      SPGCONST int rotations[][3][3])
{
  return extract_triplets_reciprocal_mesh_at_q(triplets_at_q,
					       weight_triplets_at_q,
					       fixed_grid_number,
					       num_triplets,
					       triplets,
					       mesh,
					       is_time_reversal,
					       num_rot,
					       rotations);
}


static SpglibDataset * get_dataset(SPGCONST double lattice[3][3],
				   SPGCONST double position[][3],
				   const int types[],
				   const int num_atom,
				   const double symprec)
{
  int i, j;
  int *mapping_table, *wyckoffs, *equiv_atoms, *equiv_atoms_prim;
  double tolerance;
  Spacegroup spacegroup;
  SpglibDataset *dataset;
  Cell *cell, *primitive;
  double inv_mat[3][3];
  Symmetry *symmetry;

  dataset = (SpglibDataset*) malloc(sizeof(SpglibDataset));

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);

  mapping_table = (int*) malloc(sizeof(int) * cell->size);
  primitive = prm_get_primitive_with_mapping_table(mapping_table,
						   cell,
						   symprec);
  tolerance = prm_get_current_tolerance();
  spacegroup = spa_get_spacegroup_with_primitive(primitive, tolerance);

  if (spacegroup.number > 0) {
    /* Spacegroup type, transformation matrix, origin shift */
    dataset->spacegroup_number = spacegroup.number;
    dataset->hall_number = spacegroup.hall_number;
    strcpy(dataset->international_symbol, spacegroup.international_short);
    strcpy(dataset->hall_symbol, spacegroup.hall_symbol);
    mat_inverse_matrix_d3(inv_mat, lattice, tolerance);
    mat_multiply_matrix_d3(dataset->transformation_matrix,
			   inv_mat,
			   spacegroup.bravais_lattice);
    mat_copy_vector_d3(dataset->origin_shift, spacegroup.origin_shift);

    /* Wyckoff positions */
    wyckoffs = (int*) malloc(sizeof(int) * primitive->size);
    equiv_atoms_prim = (int*) malloc(sizeof(int) * primitive->size);
    for (i = 0; i < primitive->size; i++) {
      wyckoffs[i] = -1;
      equiv_atoms_prim[i] = -1;
    }
    ref_get_Wyckoff_positions(wyckoffs, 
			      equiv_atoms_prim,
			      primitive,
			      &spacegroup,
			      tolerance);
    dataset->n_atoms = cell->size;
    dataset->wyckoffs = (int*) malloc(sizeof(int) * cell->size); 
    for (i = 0; i < cell->size; i++) {
      dataset->wyckoffs[i] = wyckoffs[mapping_table[i]];
    }
    
    free(wyckoffs);
    wyckoffs = NULL;

    /* Equivalent atoms */
    dataset->equivalent_atoms = (int*) malloc(sizeof(int) * cell->size);
    equiv_atoms = (int*) malloc(sizeof(int) * primitive->size);
    for (i = 0; i < primitive->size; i++) {
      for (j = 0; j < cell->size; j++) {
	if (mapping_table[j] == equiv_atoms_prim[i]) {
	  equiv_atoms[i] = j;
	  break;
	}
      }
    }
    for (i = 0; i < cell->size; i++) {
      dataset->equivalent_atoms[i] = equiv_atoms[mapping_table[i]];
    }
    free(equiv_atoms);
    equiv_atoms = NULL;

    free(equiv_atoms_prim);
    equiv_atoms_prim = NULL;

    /* Symmetry operations */
    symmetry = ref_get_refined_symmetry_operations(cell,
						   primitive,
						   &spacegroup,
						   tolerance);
    dataset->rotations = (int (*)[3][3]) malloc(sizeof(int[3][3]) * symmetry->size);
    dataset->translations = (double (*)[3]) malloc(sizeof(double[3]) * symmetry->size);
    dataset->n_operations = symmetry->size;
    for (i = 0; i < symmetry->size; i++) {
      mat_copy_matrix_i3(dataset->rotations[i], symmetry->rot[i]);
      mat_copy_vector_d3(dataset->translations[i], symmetry->trans[i]);
    }
    sym_free_symmetry(symmetry);

  } else {
    dataset->spacegroup_number = 0;
  }

  free(mapping_table);
  mapping_table = NULL;

  cel_free_cell(primitive);

  if (dataset->spacegroup_number == 0) {
    strcpy(dataset->international_symbol, "");
    strcpy(dataset->hall_symbol, "");
    dataset->origin_shift[0] = 0;
    dataset->origin_shift[1] = 0;
    dataset->origin_shift[2] = 0;
    dataset->n_atoms = 0;
    dataset->wyckoffs = NULL;
    dataset->equivalent_atoms = NULL;
    dataset->n_operations = 0;
    dataset->rotations = NULL;
    dataset->translations = NULL;
  }
  
  cel_free_cell(cell);
  return dataset;
}

static int get_symmetry(int rotation[][3][3],
			double translation[][3],
			const int max_size,
			SPGCONST double lattice[3][3],
			SPGCONST double position[][3],
			const int types[],
			const int num_atom,
			const double symprec)
{
  int i, j, size;
  Symmetry *symmetry;
  Cell *cell;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);
  symmetry = sym_get_operation(cell, symprec);

  if (symmetry->size > max_size) {
    fprintf(stderr, "spglib: Indicated max size(=%d) is less than number ", max_size);
    fprintf(stderr, "spglib: of symmetry operations(=%d).\n", symmetry->size);
    sym_free_symmetry(symmetry);
    return 0;
  }

  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotation[i], symmetry->rot[i]);
    for (j = 0; j < 3; j++) {
      translation[i][j] = symmetry->trans[i][j];
    }
  }

  size = symmetry->size;

  cel_free_cell(cell);
  sym_free_symmetry(symmetry);

  return size;
}

static int get_symmetry_with_collinear_spin(int rotation[][3][3],
					    double translation[][3],
					    const int max_size,
					    SPGCONST double lattice[3][3],
					    SPGCONST double position[][3],
					    const int types[],
					    const double spins[],
					    const int num_atom,
					    const double symprec)
{
  int i, j, size;
  Symmetry *symmetry;
  Cell *cell;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);
  symmetry = spn_get_collinear_operation(cell, spins, symprec);
  
  if (symmetry->size > max_size) {
    fprintf(stderr, "spglib: Indicated max size(=%d) is less than number ", max_size);
    fprintf(stderr, "spglib: of symmetry operations(=%d).\n", symmetry->size);
    sym_free_symmetry(symmetry);
    return 0;
  }

  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotation[i], symmetry->rot[i]);
    for (j = 0; j < 3; j++) {
      translation[i][j] = symmetry->trans[i][j];
    }
  }

  size = symmetry->size;

  cel_free_cell(cell);
  sym_free_symmetry(symmetry);

  return size;
}

static int get_multiplicity(SPGCONST double lattice[3][3],
			    SPGCONST double position[][3],
			    const int types[],
			    const int num_atom,
			    const double symprec)
{
  Symmetry *symmetry;
  Cell *cell;
  int size;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);
  symmetry = sym_get_operation(cell, symprec);

  size = symmetry->size;

  cel_free_cell(cell);
  sym_free_symmetry(symmetry);

  return size;
}

static int find_primitive(double lattice[3][3],
			  double position[][3],
			  int types[],
			  const int num_atom,
			  const double symprec)
{
  int i, j, num_prim_atom=0;
  Cell *cell, *primitive;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);

  /* find primitive cell */
  primitive = prm_get_primitive(cell, symprec);
  if (primitive->size == cell->size) { /* Already primitive */
    num_prim_atom = 0;
  } else { /* Primitive cell was found. */
    num_prim_atom = primitive->size;
    if (num_prim_atom < num_atom && num_prim_atom > 0 ) {
      mat_copy_matrix_d3(lattice, primitive->lattice);
      for (i = 0; i < primitive->size; i++) {
	types[i] = primitive->types[i];
	for (j=0; j<3; j++) {
	  position[i][j] = primitive->position[i][j];
	}
      }
    }
  }

  cel_free_cell(primitive);
  cel_free_cell(cell);
    
  return num_prim_atom;
}

static int get_international(char symbol[11],
			     SPGCONST double lattice[3][3],
			     SPGCONST double position[][3],
			     const int types[],
			     const int num_atom,
			     const double symprec)
{
  Cell *cell;
  Spacegroup spacegroup;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);
  spacegroup = spa_get_spacegroup(cell, symprec);
  if (spacegroup.number > 0) {
    strcpy(symbol, spacegroup.international_short);
  }

  cel_free_cell(cell);
  
  return spacegroup.number;
}

static int get_schoenflies(char symbol[10],
			   SPGCONST double lattice[3][3],
			   SPGCONST double position[][3],
			   const int types[],
			   const int num_atom,
			   const double symprec)
{
  Cell *cell;
  Spacegroup spacegroup;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);

  spacegroup = spa_get_spacegroup(cell, symprec);
  if (spacegroup.number > 0) {
    strcpy(symbol, spacegroup.schoenflies);
  }

  cel_free_cell(cell);

  return spacegroup.number;
}

static int refine_cell(double lattice[3][3],
		       double position[][3],
		       int types[],
		       const int num_atom,
		       const double symprec)
{
  int i, num_atom_bravais;
  Cell *cell, *bravais;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);

  bravais = ref_refine_cell(cell, symprec);
  cel_free_cell(cell);

  if (bravais->size > 0) {
    mat_copy_matrix_d3(lattice, bravais->lattice);
    num_atom_bravais = bravais->size;
    for (i = 0; i < bravais->size; i++) {
      types[i] = bravais->types[i];
      mat_copy_vector_d3(position[i], bravais->position[i]);
    }
  } else {
    num_atom_bravais = 0;
  }

  cel_free_cell(bravais);
  
  return num_atom_bravais;
}

static int get_ir_kpoints(int map[],
			  SPGCONST double kpoints[][3],
			  const int num_kpoint,
			  SPGCONST double lattice[3][3],
			  SPGCONST double position[][3],
			  const int types[],
			  const int num_atom,
			  const int is_time_reversal,
			  const double symprec)
{
  Symmetry *symmetry;
  Cell *cell;
  int num_ir_kpoint;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);
  symmetry = sym_get_operation(cell, symprec);

  num_ir_kpoint = kpt_get_irreducible_kpoints(map,
					      kpoints,
					      num_kpoint,
					      symmetry,
					      is_time_reversal,
					      symprec);


  cel_free_cell(cell);
  sym_free_symmetry(symmetry);

  return num_ir_kpoint;
}

static int get_ir_reciprocal_mesh(int grid_point[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  const int is_time_reversal,
				  SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec)
{
  Symmetry *symmetry;
  Cell *cell;
  int num_ir;

  cell = cel_alloc_cell(num_atom);
  cel_set_cell(cell, lattice, position, types);
  symmetry = sym_get_operation(cell, symprec);

  num_ir = kpt_get_irreducible_reciprocal_mesh(grid_point,
					       map,
					       mesh,
					       is_shift,
					       is_time_reversal,
					       symmetry);


  cel_free_cell(cell);
  sym_free_symmetry(symmetry);

  return num_ir;
}

static int get_stabilized_reciprocal_mesh(int grid_point[][3],
					  int map[],
					  const int mesh[3],
					  const int is_shift[3],
					  const int is_time_reversal,
					  const int num_rot,
					  SPGCONST int rotations[][3][3],
					  const int num_q,
					  SPGCONST double qpoints[][3])
{
  MatINT *rot_real;
  int i, num_ir;
  
  rot_real = mat_alloc_MatINT(num_rot);
  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot_real->mat[i], rotations[i]);
  }

  num_ir = kpt_get_stabilized_reciprocal_mesh(grid_point,
					      map,
					      mesh,
					      is_shift,
					      is_time_reversal,
					      rot_real,
					      num_q,
					      qpoints);

  mat_free_MatINT(rot_real);

  return num_ir;
}

static SpglibTriplets * get_triplets_reciprocal_mesh(const int mesh[3],
						     const int is_time_reversal,
						     const int num_rot,
						     SPGCONST int rotations[][3][3])
{
  int i, j, num_grid;
  MatINT *rot_real;
  Triplets *tps;
  SpglibTriplets *spg_triplets;
  
  num_grid = mesh[0] * mesh[1] * mesh[2];
  rot_real = mat_alloc_MatINT(num_rot);
  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot_real->mat[i], rotations[i]);
  }

  tps = kpt_get_triplets_reciprocal_mesh(mesh,
					 is_time_reversal,
					 rot_real);
  mat_free_MatINT(rot_real);

  spg_triplets = (SpglibTriplets*) malloc(sizeof(SpglibTriplets));
  spg_triplets->size = tps->size;
  spg_triplets->triplets = (int (*)[3]) malloc(sizeof(int[3]) * tps->size);
  spg_triplets->weights = (int*) malloc(sizeof(int) * tps->size);
  spg_triplets->mesh_points = (int (*)[3]) malloc(sizeof(int[3]) * num_grid);

  for (i = 0; i < 3; i++) {
    spg_triplets->mesh[i] = tps->mesh[i];
  }
  for (i = 0; i < num_grid; i++) {
    for (j = 0; j < 3; j++) {
      spg_triplets->mesh_points[i][j] = tps->mesh_points[i][j];
    }
  }

  for (i = 0; i < tps->size; i++) {
    for (j = 0; j < 3; j++) {
      spg_triplets->triplets[i][j] = tps->triplets[i][j];
    }
    spg_triplets->weights[i] = tps->weights[i];
  }
  kpt_free_triplets(tps);

  return spg_triplets;
}

static int get_triplets_reciprocal_mesh_at_q(int weights[],
					     int grid_points[][3],
					     int third_q[],
					     const int grid_point,
					     const int mesh[3],
					     const int is_time_reversal,
					     const int num_rot,
					     SPGCONST int rotations[][3][3])
{
  MatINT *rot_real;
  int i, num_ir;
  
  rot_real = mat_alloc_MatINT(num_rot);
  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot_real->mat[i], rotations[i]);
  }

  num_ir = kpt_get_ir_triplets_at_q(weights,
				    grid_points,
				    third_q,
				    grid_point,
				    mesh,
				    is_time_reversal,
				    rot_real);

  mat_free_MatINT(rot_real);

  return num_ir;
}

static int extract_triplets_reciprocal_mesh_at_q(int triplets_at_q[][3],
						 int weight_triplets_at_q[],
						 const int fixed_grid_number,
						 const int num_triplets,
						 SPGCONST int triplets[][3],
						 const int mesh[3],
						 const int is_time_reversal,
						 const int num_rot,
						 SPGCONST int rotations[][3][3])
{
  MatINT *rot_real;
  int i, num_ir;
  
  rot_real = mat_alloc_MatINT(num_rot);
  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot_real->mat[i], rotations[i]);
  }

  num_ir = kpt_extract_triplets_reciprocal_mesh_at_q(triplets_at_q,
						     weight_triplets_at_q,
						     fixed_grid_number,
						     num_triplets,
						     triplets,
						     mesh,
						     is_time_reversal,
						     rot_real);

  
  mat_free_MatINT(rot_real);

  return num_ir;
}
