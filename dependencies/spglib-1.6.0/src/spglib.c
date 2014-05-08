/* spglib.c */
/* Copyright (C) 2008 Atsushi Togo */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cell.h"
#include "kpoint.h"
#include "lattice.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "spglib.h"
#include "primitive.h"
#include "refinement.h"
#include "spacegroup.h"
#include "spg_database.h"
#include "spin.h"
#include "symmetry.h"
#include "tetrahedron_method.h"

/*---------*/
/* general */
/*---------*/
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
static int get_symmetry_from_dataset(int rotation[][3][3],
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

/*---------*/
/* kpoints */
/*---------*/
static int get_ir_reciprocal_mesh(int grid_address[][3],
				  int map[],
				  const int mesh[3],
				  const int is_shift[3],
				  const int is_time_reversal,
				  SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec);

static int get_stabilized_reciprocal_mesh(int grid_address[][3],
					  int map[],
					  const int mesh[3],
					  const int is_shift[3],
					  const int is_time_reversal,
					  const int num_rot,
					  SPGCONST int rotations[][3][3],
					  const int num_q,
					  SPGCONST double qpoints[][3]);
static int get_triplets_reciprocal_mesh_at_q(int weights[],
					     int grid_address[][3],
					     int third_q[],
					     const int grid_point,
					     const int mesh[3],
					     const int is_time_reversal,
					     const int num_rot,
					     SPGCONST int rotations[][3][3]);


/*========*/
/* global */
/*========*/

/*---------*/
/* general */
/*---------*/
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

  return get_symmetry_from_dataset(rotation,
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

  return get_symmetry_from_dataset(rotation,
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

SpglibSpacegroupType spg_get_spacegroup_type(int hall_number)
{
  SpglibSpacegroupType spglibtype;
  SpacegroupType spgtype;

  spgtype = spgdb_get_spacegroup_type(hall_number);
  spglibtype.number = spgtype.number;
  strcpy(spglibtype.schoenflies, spgtype.schoenflies);
  strcpy(spglibtype.hall_symbol, spgtype.hall_symbol);
  strcpy(spglibtype.international, spgtype.international);
  strcpy(spglibtype.international_full, spgtype.international_full);
  strcpy(spglibtype.international_short, spgtype.international_short);
  
  return spglibtype;
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

/*---------*/
/* kpoints */
/*---------*/
int spg_get_ir_reciprocal_mesh(int grid_address[][3],
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

  return get_ir_reciprocal_mesh(grid_address,
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

int spg_get_stabilized_reciprocal_mesh(int grid_address[][3],
				       int map[],
				       const int mesh[3],
				       const int is_shift[3],
				       const int is_time_reversal,
				       const int num_rot,
				       SPGCONST int rotations[][3][3],
				       const int num_q,
				       SPGCONST double qpoints[][3])
{
  return get_stabilized_reciprocal_mesh(grid_address,
					map,
					mesh,
					is_shift,
					is_time_reversal,
					num_rot,
					rotations,
					num_q,
					qpoints);
}

int spg_relocate_BZ_grid_address(int bz_grid_address[][3],
				 int bz_map[],
				 SPGCONST int grid_address[][3],
				 const int mesh[3],
				 SPGCONST double rec_lattice[3][3],
				 const int is_shift[3])
{
  return kpt_relocate_BZ_grid_address(bz_grid_address,
				      bz_map,
				      grid_address,
				      mesh,
				      rec_lattice,
				      is_shift);
}

int spg_get_triplets_reciprocal_mesh_at_q(int weights[],
					  int grid_address[][3],
					  int third_q[],
					  const int grid_point,
					  const int mesh[3],
					  const int is_time_reversal,
					  const int num_rot,
					  SPGCONST int rotations[][3][3])
{
  return get_triplets_reciprocal_mesh_at_q(weights,
					   grid_address,
					   third_q,
					   grid_point,
					   mesh,
					   is_time_reversal,
					   num_rot,
					   rotations);
}

int spg_get_BZ_triplets_at_q(int triplets[][3],
			     const int grid_point,
			     SPGCONST int bz_grid_address[][3],
			     const int bz_map[],
			     const int triplet_weights[],
			     const int mesh[3])

{
  return kpt_get_BZ_triplets_at_q(triplets,
				  grid_point,
				  bz_grid_address,
				  bz_map,
				  triplet_weights,
				  mesh);
}

void spg_get_neighboring_grid_points(int relative_grid_points[],
				     const int grid_point,
				     SPGCONST int relative_grid_address[][3],
				     const int num_relative_grid_address,
				     const int mesh[3],
				     SPGCONST int bz_grid_address[][3],
				     const int bz_map[])
{
  kpt_get_neighboring_grid_points(relative_grid_points,
				  grid_point,
				  relative_grid_address,
				  num_relative_grid_address,
				  mesh,
				  bz_grid_address,
				  bz_map);
}

/*--------------------*/
/* tetrahedron method */
/*--------------------*/
void
spg_get_tetrahedra_relative_grid_address(int relative_grid_address[24][4][3],
					 SPGCONST double rec_lattice[3][3])
{
  thm_get_relative_grid_address(relative_grid_address, rec_lattice);
}

double
spg_get_tetrahedra_integration_weight(const double omega,
				      SPGCONST double tetrahedra_omegas[24][4],
				      const char function)
{
  return thm_get_integration_weight(omega,
				    tetrahedra_omegas,
				    function);
}

void
spg_get_tetrahedra_integration_weight_at_omegas
(double integration_weights[],
 const int num_omegas,
 const double omegas[],
 SPGCONST double tetrahedra_omegas[24][4],
 const char function)
{
  thm_get_integration_weight_at_omegas(integration_weights,
				       num_omegas,
				       omegas,
				       tetrahedra_omegas,
				       function);
}


/*=======*/
/* local */
/*=======*/

/*---------*/
/* general */
/*---------*/
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
    dataset->rotations = (int (*)[3][3])
      malloc(sizeof(int[3][3]) * symmetry->size);
    dataset->translations = (double (*)[3])
      malloc(sizeof(double[3]) * symmetry->size);
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

static int get_symmetry_from_dataset(int rotation[][3][3],
				     double translation[][3],
				     const int max_size,
				     SPGCONST double lattice[3][3],
				     SPGCONST double position[][3],
				     const int types[],
				     const int num_atom,
				     const double symprec)
{
  int i, j, k, num_sym;
  SpglibDataset *dataset;

  dataset = get_dataset(lattice,
			position,
			types,
			num_atom,
			symprec);
  
  if (dataset->n_operations > max_size) {
    fprintf(stderr,
	    "spglib: Indicated max size(=%d) is less than number ", max_size);
    fprintf(stderr,
	    "spglib: of symmetry operations(=%d).\n", dataset->n_operations);
    num_sym = 0;
    goto ret;
  }

  num_sym = dataset->n_operations;
  for (i = 0; i < num_sym; i++) {
    for (j = 0; j < 3; j++) {
      translation[i][j] = dataset->translations[i][j];
      for (k = 0; k < 3; k++) {
	rotation[i][j][k] = dataset->rotations[i][j][k];
      }
    }	  
  }
  
 ret:
  spg_free_dataset(dataset);
  return num_sym;
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

/*---------*/
/* kpoints */
/*---------*/
static int get_ir_reciprocal_mesh(int grid_address[][3],
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
  SpglibDataset *dataset;
  int num_ir, i;
  MatINT *rotations;

  dataset = get_dataset(lattice,
			position,
			types,
			num_atom,
			symprec);
  rotations = mat_alloc_MatINT(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    mat_copy_matrix_i3(rotations->mat[i], dataset->rotations[i]);
  }
  num_ir = kpt_get_irreducible_reciprocal_mesh(grid_address,
					       map,
					       mesh,
					       is_shift,
					       is_time_reversal,
					       rotations);
  mat_free_MatINT(rotations);
  spg_free_dataset(dataset);
  return num_ir;
}

static int get_stabilized_reciprocal_mesh(int grid_address[][3],
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

  num_ir = kpt_get_stabilized_reciprocal_mesh(grid_address,
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

static int get_triplets_reciprocal_mesh_at_q(int weights[],
					     int grid_address[][3],
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
				    grid_address,
				    third_q,
				    grid_point,
				    mesh,
				    is_time_reversal,
				    rot_real);

  mat_free_MatINT(rot_real);

  return num_ir;
}

