/* Copyright (C) 2008 Atsushi Togo */
/* All rights reserved. */

/* This file is part of spglib. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the phonopy project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cell.h"
#include "debug.h"
#include "kgrid.h"
#include "kpoint.h"
#include "lattice.h"
#include "mathfunc.h"
#include "niggli.h"
#include "pointgroup.h"
#include "spglib.h"
#include "primitive.h"
#include "refinement.h"
#include "spacegroup.h"
#include "spg_database.h"
#include "spin.h"
#include "symmetry.h"
#include "version.h"

#define REDUCE_RATE 0.95

/*---------*/
/* general */
/*---------*/
static SpglibDataset * get_dataset(SPGCONST double lattice[3][3],
				   SPGCONST double position[][3],
				   const int types[],
				   const int num_atom,
				   const int hall_number,
				   const double symprec);
static int set_dataset(SpglibDataset * dataset,
		       SPGCONST Cell * cell,
		       SPGCONST Cell * primitive,
		       SPGCONST Spacegroup * spacegroup,
		       const int * mapping_table,
		       const double tolerance);
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
					    int equivalent_atoms[],
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
static int standardize_primitive(double lattice[3][3],
				 double position[][3],
				 int types[],
				 const int num_atom,
				 const double symprec);
static int standardize_cell(double lattice[3][3],
			    double position[][3],
			    int types[],
			    const int num_atom,
			    const double symprec);
static int get_standardized_cell(double lattice[3][3],
				 double position[][3],
				 int types[],
				 const int num_atom,
				 const int to_primitive,
				 const double symprec);
static Centering get_centering(int hall_number);
static void set_cell(double lattice[3][3],
		     double position[][3],
		     int types[],
		     Cell * cell);
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
static int get_symmetry_numerical(int rotation[][3][3],
				  double translation[][3],
				  const int max_size,
				  SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec);
static int is_rhombohedral(const int hall_number);

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

/*========*/
/* global */
/*========*/

/*--------------------------------------------*/
/* Version: spglib-[major].[minor].[micro] */
/*--------------------------------------------*/
int spg_get_major_version(void)
{
  return SPGLIB_MAJOR_VERSION;
}

int spg_get_minor_version(void)
{
  return SPGLIB_MINOR_VERSION;
}

int spg_get_micro_version(void)
{
  return SPGLIB_MICRO_VERSION;
}

/*---------*/
/* general */
/*---------*/
/* Return NULL if failed */
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
		     0,
		     symprec);
}

/* Return NULL if failed */
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
		     0,
		     symprec);
}

/* Return NULL if failed */
SpglibDataset * spg_get_dataset_with_hall_number(SPGCONST double lattice[3][3],
						 SPGCONST double position[][3],
						 const int types[],
						 const int num_atom,
						 const int hall_number,
						 const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_dataset(lattice,
		     position,
		     types,
		     num_atom,
		     hall_number,
		     symprec);
}

/* Return NULL if failed */
SpglibDataset *
spgat_get_dataset_with_hall_number(SPGCONST double lattice[3][3],
				   SPGCONST double position[][3],
				   const int types[],
				   const int num_atom,
				   const int hall_number,
				   const double symprec,
				   const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return get_dataset(lattice,
		     position,
		     types,
		     num_atom,
		     hall_number,
		     symprec);
}

void spg_free_dataset(SpglibDataset *dataset)
{
  if (dataset->n_operations > 0) {
    free(dataset->rotations);
    dataset->rotations = NULL;
    free(dataset->translations);
    dataset->translations = NULL;
    dataset->n_operations = 0;
  }

  if (dataset->n_atoms > 0) {
    free(dataset->wyckoffs);
    dataset->wyckoffs = NULL;
    free(dataset->equivalent_atoms);
    dataset->equivalent_atoms = NULL;
    dataset->n_atoms = 0;
  }

  if (dataset->n_std_atoms > 0) {
    free(dataset->std_positions);
    dataset->std_positions = NULL;
    free(dataset->std_types);
    dataset->std_types = NULL;
    dataset->n_std_atoms = 0;
  }

  dataset->spacegroup_number = 0;
  dataset->hall_number = 0;
  strcpy(dataset->international_symbol, "");
  strcpy(dataset->hall_symbol, "");
  strcpy(dataset->setting, "");
  
  free(dataset);
  dataset = NULL;
}

/* Return 0 if failed */
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

/* Return 0 if failed */
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

/* Return 0 if failed */
int spg_get_symmetry_numerical(int rotation[][3][3],
			       double translation[][3],
			       const int max_size,
			       SPGCONST double lattice[3][3],
			       SPGCONST double position[][3],
			       const int types[],
			       const int num_atom,
			       const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return get_symmetry_numerical(rotation,
				translation,
				max_size,
				lattice,
				position,
				types,
				num_atom,
				symprec);
}

/* Return 0 if failed */
int spgat_get_symmetry_numerical(int rotation[][3][3],
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

  return get_symmetry_numerical(rotation,
				translation,
				max_size,
				lattice,
				position,
				types,
				num_atom,
				symprec);
}

/* Return 0 if failed */
int spg_get_symmetry_with_collinear_spin(int rotation[][3][3],
					 double translation[][3],
					 int equivalent_atoms[],
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
					  equivalent_atoms,
					  max_size,
					  lattice,
					  position,
					  types,
					  spins,
					  num_atom,
					  symprec);
}

/* Return 0 if failed */
int spgat_get_symmetry_with_collinear_spin(int rotation[][3][3],
					   double translation[][3],
					   int equivalent_atoms[],
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
					  equivalent_atoms,
					  max_size,
					  lattice,
					  position,
					  types,
					  spins,
					  num_atom,
					  symprec);
}

/* Return 0 if failed */
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

/* Return 0 if failed */
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

/* Return 0 if failed */
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

/* Return 0 if failed */
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

/* Return 0 if failed */
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

/* Return 0 if failed */
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

/* Return 0 if failed */
int spg_get_pointgroup(char symbol[6],
		       int transform_mat[3][3],
		       SPGCONST int rotations[][3][3],
		       const int num_rotations)
{
  Pointgroup pointgroup;

  pointgroup = ptg_get_transformation_matrix(transform_mat,
					     rotations,
					     num_rotations);

  if (pointgroup.number == 0) {
    return 0;
  }

  strcpy(symbol, pointgroup.symbol);

  return pointgroup.number;
}

/* Return 0 if failed */
int spg_get_symmetry_from_database(int rotations[192][3][3],
				   double translations[192][3],
				   const int hall_number)
{
  int i, size;
  Symmetry *symmetry;

  symmetry = NULL;

  if ((symmetry = spgdb_get_spacegroup_operations(hall_number)) == NULL) {
    return 0;
  }

  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotations[i], symmetry->rot[i]);
    mat_copy_vector_d3(translations[i], symmetry->trans[i]);
  }
  size = symmetry->size;

  sym_free_symmetry(symmetry);

  return size;
}

/* Return spglibtype.number = 0 if failed */
SpglibSpacegroupType spg_get_spacegroup_type(const int hall_number)
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

/* Return 0 if failed */
int spg_standardize_cell(double lattice[3][3],
			 double position[][3],
			 int types[],
			 const int num_atom,
			 const int to_primitive,
			 const int no_idealize,
			 const double symprec)
{
  return spgat_standardize_cell(lattice,
				position,
				types,
				num_atom,
				to_primitive,
				no_idealize,
				symprec,
				-1.0);
}

/* Return 0 if failed */
int spgat_standardize_cell(double lattice[3][3],
			   double position[][3],
			   int types[],
			   const int num_atom,
			   const int to_primitive,
			   const int no_idealize,
			   const double symprec,
			   const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  if (to_primitive) {
    if (no_idealize) {
      return get_standardized_cell(lattice,
				   position,
				   types,
				   num_atom,
				   1,
				   symprec);
    } else {
      return standardize_primitive(lattice,
				   position,
				   types,
				   num_atom,
				   symprec);
    }
  } else {
    if (no_idealize) {
      return get_standardized_cell(lattice,
				   position,
				   types,
				   num_atom,
				   0,
				   symprec);
    } else {
      return standardize_cell(lattice,
			      position,
			      types,
			      num_atom,
			      symprec);
    }
  }
}

/* Return 0 if failed */
int spg_find_primitive(double lattice[3][3],
		       double position[][3],
		       int types[],
		       const int num_atom,
		       const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return standardize_primitive(lattice,
			       position,
			       types,
			       num_atom,
			       symprec);
}

/* Return 0 if failed */
int spgat_find_primitive(double lattice[3][3],
			 double position[][3],
			 int types[],
			 const int num_atom,
			 const double symprec,
			 const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return standardize_primitive(lattice,
			       position,
			       types,
			       num_atom,
			       symprec);
}

/* Return 0 if failed */
int spg_refine_cell(double lattice[3][3],
		    double position[][3],
		    int types[],
		    const int num_atom,
		    const double symprec)
{
  sym_set_angle_tolerance(-1.0);

  return standardize_cell(lattice,
			  position,
			  types,
			  num_atom,
			  symprec);
}

/* Return 0 if failed */
int spgat_refine_cell(double lattice[3][3],
		      double position[][3],
		      int types[],
		      const int num_atom,
		      const double symprec,
		      const double angle_tolerance)
{
  sym_set_angle_tolerance(angle_tolerance);

  return standardize_cell(lattice,
			  position,
			  types,
			  num_atom,
			  symprec);
}

/*---------*/
/* kpoints */
/*---------*/
int spg_get_grid_point_from_address(const int grid_address[3],
				    const int mesh[3])
{
  int address_double[3];
  int is_shift[3];

  is_shift[0] = 0;
  is_shift[1] = 0;
  is_shift[2] = 0;
  kgd_get_grid_address_double_mesh(address_double,
				   grid_address,
				   mesh,
				   is_shift);
  return kgd_get_grid_point_double_mesh(address_double, mesh);
}

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

int spg_get_grid_points_by_rotations(int rot_grid_points[],
				     const int address_orig[3],
				     const int num_rot,
				     SPGCONST int rot_reciprocal[][3][3],
				     const int mesh[3],
				     const int is_shift[3])
{
  int i;
  MatINT *rot;

  rot = NULL;

  if ((rot = mat_alloc_MatINT(num_rot)) == NULL) {
    return 0;
  }

  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot->mat[i], rot_reciprocal[i]);
  }
  kpt_get_grid_points_by_rotations(rot_grid_points,
				   address_orig,
				   rot,
				   mesh,
				   is_shift);
  mat_free_MatINT(rot);

  return 1;
}

int spg_get_BZ_grid_points_by_rotations(int rot_grid_points[],
					const int address_orig[3],
					const int num_rot,
					SPGCONST int rot_reciprocal[][3][3],
					const int mesh[3],
					const int is_shift[3],
					const int bz_map[])
{
  int i;
  MatINT *rot;

  rot = NULL;

  if ((rot = mat_alloc_MatINT(num_rot)) == NULL) {
    return 0;
  }

  for (i = 0; i < num_rot; i++) {
    mat_copy_matrix_i3(rot->mat[i], rot_reciprocal[i]);
  }
  kpt_get_BZ_grid_points_by_rotations(rot_grid_points,
				      address_orig,
				      rot,
				      mesh,
				      is_shift,
				      bz_map);
  mat_free_MatINT(rot);

  return 1;
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

/*--------*/
/* Niggli */
/*--------*/
/* Return 0 if failed */
int spg_niggli_reduce(double lattice[3][3], const double symprec)
{
  int i, j, succeeded;
  double vals[9];
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      vals[i * 3 + j] = lattice[i][j];
    }
  }

  succeeded = niggli_reduce(vals, symprec);

  if (succeeded) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
	lattice[i][j] = vals[i * 3 + j];
      }
    }
  }

  return succeeded;
}




/*=======*/
/* local */
/*=======*/

/*---------*/
/* general */
/*---------*/
/* Return NULL if failed */
static SpglibDataset * get_dataset(SPGCONST double lattice[3][3],
				   SPGCONST double position[][3],
				   const int types[],
				   const int num_atom,
				   const int hall_number,
				   const double symprec)
{
  Spacegroup spacegroup;
  SpacegroupType spacegroup_type;
  SpglibDataset *dataset;
  Cell *cell;
  Primitive *primitive;

  spacegroup.number = 0;
  dataset = NULL;
  cell = NULL;
  primitive = NULL;

  if ((dataset = (SpglibDataset*) malloc(sizeof(SpglibDataset))) == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    return NULL;
  }

  dataset->spacegroup_number = 0;
  dataset->hall_number = 0;
  strcpy(dataset->international_symbol, "");
  strcpy(dataset->hall_symbol, "");
  strcpy(dataset->setting, "");
  dataset->origin_shift[0] = 0;
  dataset->origin_shift[1] = 0;
  dataset->origin_shift[2] = 0;
  dataset->n_atoms = 0;
  dataset->wyckoffs = NULL;
  dataset->equivalent_atoms = NULL;
  dataset->n_operations = 0;
  dataset->rotations = NULL;
  dataset->translations = NULL;
  dataset->n_std_atoms = 0;
  dataset->std_positions = NULL;
  dataset->std_types = NULL;
  dataset->pointgroup_number = 0;
  strcpy(dataset->pointgroup_symbol, "");

  if ((cell = cel_alloc_cell(num_atom)) == NULL) {
    free(dataset);
    dataset = NULL;
    return NULL;
  }

  cel_set_cell(cell, lattice, position, types);

  primitive = spa_get_spacegroup(&spacegroup, cell, symprec);

  if ((spacegroup.number > 0) && (primitive != NULL)) {

    /* With hall_number > 0, specific choice is searched. */
    if (hall_number > 0) {
      spacegroup_type = spgdb_get_spacegroup_type(hall_number);
      if (spacegroup.number == spacegroup_type.number) {
	spacegroup = spa_get_spacegroup_with_hall_number(primitive,
							 hall_number);
      } else {
	goto err;
      }

      if (spacegroup.number == 0) {
	goto err;
      }
    }

    if (spacegroup.number > 0) {
      if ((set_dataset(dataset,
		       cell,
		       primitive->cell,
		       &spacegroup,
		       primitive->mapping_table,
		       primitive->tolerance)) == 0) {
	goto err;
      }
    }
  }

  cel_free_cell(cell);
  prm_free_primitive(primitive);

  return dataset;

 err:
  cel_free_cell(cell);
  prm_free_primitive(primitive);
  free(dataset);
  dataset = NULL;
  return NULL;
}

/* Return 0 if failed */
static int set_dataset(SpglibDataset * dataset,
		       SPGCONST Cell * cell,
		       SPGCONST Cell * primitive,
		       SPGCONST Spacegroup * spacegroup,
		       const int * mapping_table,
		       const double tolerance)
{
  int i;
  double inv_lat[3][3];
  Cell *bravais;
  Symmetry *symmetry;
  Pointgroup pointgroup;

  bravais = NULL;
  symmetry = NULL;

  /* Spacegroup type, transformation matrix, origin shift */
  dataset->n_atoms = cell->size;
  dataset->spacegroup_number = spacegroup->number;
  dataset->hall_number = spacegroup->hall_number;
  strcpy(dataset->international_symbol, spacegroup->international_short);
  strcpy(dataset->hall_symbol, spacegroup->hall_symbol);
  strcpy(dataset->setting, spacegroup->setting);
  mat_inverse_matrix_d3(inv_lat, spacegroup->bravais_lattice, 0);
  mat_multiply_matrix_d3(dataset->transformation_matrix,
			 inv_lat,
			 cell->lattice);
  mat_copy_vector_d3(dataset->origin_shift, spacegroup->origin_shift);

  /* Symmetry operations */
  if ((symmetry = ref_get_refined_symmetry_operations(cell,
						      primitive,
						      spacegroup,
						      tolerance)) == NULL) {
    return 0;
  }

  dataset->n_operations = symmetry->size;

  if ((dataset->rotations =
       (int (*)[3][3]) malloc(sizeof(int[3][3]) * dataset->n_operations))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((dataset->translations =
       (double (*)[3]) malloc(sizeof(double[3]) * dataset->n_operations))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(dataset->rotations[i], symmetry->rot[i]);
    mat_copy_vector_d3(dataset->translations[i], symmetry->trans[i]);
  }

  /* Wyckoff positions */
  if ((dataset->wyckoffs = (int*) malloc(sizeof(int) * dataset->n_atoms))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((dataset->equivalent_atoms =
       (int*) malloc(sizeof(int) * dataset->n_atoms))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((bravais = ref_get_Wyckoff_positions(dataset->wyckoffs, 
					   dataset->equivalent_atoms,
					   primitive,
					   cell,
					   spacegroup,
					   symmetry,
					   mapping_table,
					   tolerance)) == NULL) {
    goto err;
  }

  dataset->n_std_atoms = bravais->size;
  mat_copy_matrix_d3(dataset->std_lattice, bravais->lattice);

  if ((dataset->std_positions =
       (double (*)[3]) malloc(sizeof(double[3]) * dataset->n_std_atoms))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  if ((dataset->std_types = (int*) malloc(sizeof(int) * dataset->n_std_atoms))
      == NULL) {
    warning_print("spglib: Memory could not be allocated.");
    goto err;
  }

  for (i = 0; i < dataset->n_std_atoms; i++) {
    mat_copy_vector_d3(dataset->std_positions[i], bravais->position[i]);
    dataset->std_types[i] = bravais->types[i];
  }
  
  cel_free_cell(bravais);
  sym_free_symmetry(symmetry);

  dataset->pointgroup_number = spacegroup->pointgroup_number;
  pointgroup = ptg_get_pointgroup(spacegroup->pointgroup_number);
  strcpy(dataset->pointgroup_symbol, pointgroup.symbol);

  return 1;

 err:
  if (dataset->std_positions != NULL) {
    free(dataset->std_positions);
    dataset->std_positions = NULL;
  }
  if (bravais != NULL) {
    cel_free_cell(bravais);
  }
  if (dataset->equivalent_atoms != NULL) {
    free(dataset->equivalent_atoms);
    dataset->equivalent_atoms = NULL;
  }
  if (dataset->wyckoffs != NULL) {
    free(dataset->wyckoffs);
    dataset->wyckoffs = NULL;
  }
  if (dataset->translations != NULL) {
    free(dataset->translations);
    dataset->translations = NULL;
  }
  if (dataset->rotations != NULL) {
    free(dataset->rotations);
    dataset->rotations = NULL;
  }
  if (symmetry != NULL) {
    sym_free_symmetry(symmetry);
  }
  return 0;
}

/* Return 0 if failed */
static int get_symmetry_from_dataset(int rotation[][3][3],
				     double translation[][3],
				     const int max_size,
				     SPGCONST double lattice[3][3],
				     SPGCONST double position[][3],
				     const int types[],
				     const int num_atom,
				     const double symprec)
{
  int i, num_sym;
  SpglibDataset *dataset;

  num_sym = 0;
  dataset = NULL;

  if ((dataset = get_dataset(lattice,
			     position,
			     types,
			     num_atom,
			     0,
			     symprec)) == NULL) {
    return 0;
  }
  
  if (dataset->n_operations > max_size) {
    fprintf(stderr,
	    "spglib: Indicated max size(=%d) is less than number ", max_size);
    fprintf(stderr,
	    "spglib: of symmetry operations(=%d).\n", dataset->n_operations);
    goto ret;
  }

  num_sym = dataset->n_operations;
  for (i = 0; i < num_sym; i++) {
    mat_copy_matrix_i3(rotation[i], dataset->rotations[i]);
    mat_copy_vector_d3(translation[i], dataset->translations[i]);
  }
  
 ret:
  spg_free_dataset(dataset);
  return num_sym;
}

/* Return 0 if failed */
static int get_symmetry_with_collinear_spin(int rotation[][3][3],
					    double translation[][3],
					    int equivalent_atoms[],
					    const int max_size,
					    SPGCONST double lattice[3][3],
					    SPGCONST double position[][3],
					    const int types[],
					    const double spins[],
					    const int num_atom,
					    const double symprec)
{
  int i, size;
  Symmetry *symmetry, *sym_nonspin;
  Cell *cell;
  SpglibDataset *dataset;

  size = 0;
  symmetry = NULL;
  sym_nonspin = NULL;
  cell = NULL;
  dataset = NULL;

  if ((cell = cel_alloc_cell(num_atom)) == NULL) {
    goto err;
  }

  cel_set_cell(cell, lattice, position, types);

  if ((dataset = get_dataset(lattice,
			     position,
			     types,
			     num_atom,
			     0,
			     symprec)) == NULL) {
    cel_free_cell(cell);
    goto err;
  }

  if ((sym_nonspin = sym_alloc_symmetry(dataset->n_operations)) == NULL) {
    spg_free_dataset(dataset);
    cel_free_cell(cell);
    goto err;
  }

  for (i = 0; i < dataset->n_operations; i++) {
    mat_copy_matrix_i3(sym_nonspin->rot[i], dataset->rotations[i]);
    mat_copy_vector_d3(sym_nonspin->trans[i], dataset->translations[i]);
  }
  spg_free_dataset(dataset);

  if ((symmetry = spn_get_collinear_operations(equivalent_atoms,
					       sym_nonspin,
					       cell,
					       spins,
					       symprec)) == NULL) {
    sym_free_symmetry(sym_nonspin);
    cel_free_cell(cell);
    goto err;
  }

  sym_free_symmetry(sym_nonspin);
  
  if (symmetry->size > max_size) {
    fprintf(stderr, "spglib: Indicated max size(=%d) is less than number ",
	    max_size);
    fprintf(stderr, "spglib: of symmetry operations(=%d).\n", symmetry->size);
    goto ret;
  }

  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotation[i], symmetry->rot[i]);
    mat_copy_vector_d3(translation[i], symmetry->trans[i]);
  }

  size = symmetry->size;

 ret:
  sym_free_symmetry(symmetry);
  cel_free_cell(cell);

  return size;

 err:
  return 0;
}

/* Return 0 if failed */
static int get_multiplicity(SPGCONST double lattice[3][3],
			    SPGCONST double position[][3],
			    const int types[],
			    const int num_atom,
			    const double symprec)
{
  int size;
  SpglibDataset *dataset;

  size = 0;
  dataset = NULL;

  if ((dataset = get_dataset(lattice,
			     position,
			     types,
			     num_atom,
			     0,
			     symprec)) == NULL) {
    return 0;
  }

  size = dataset->n_operations;
  spg_free_dataset(dataset);

  return size;
}

static int standardize_primitive(double lattice[3][3],
				 double position[][3],
				 int types[],
				 const int num_atom,
				 const double symprec)
{
  int num_prim_atom;
  Centering centering;
  SpglibDataset *dataset;
  Cell *primitive, *bravais;

  double identity[3][3] = {{ 1, 0, 0 },
			   { 0, 1, 0 },
			   { 0, 0, 1 }};

  num_prim_atom = 0;
  dataset = NULL;
  primitive = NULL;
  bravais = NULL;

  if ((dataset = get_dataset(lattice,
			     position,
			     types,
			     num_atom,
			     0,
			     symprec)) == NULL) {
    return 0;
  }

  if ((centering = get_centering(dataset->hall_number)) == CENTERING_ERROR) {
    goto err;
  }

  if (is_rhombohedral(dataset->hall_number)) {
    centering = R_CENTER;
  }


  if ((bravais = cel_alloc_cell(dataset->n_std_atoms)) == NULL) {
    spg_free_dataset(dataset);
    return 0;
  }

  cel_set_cell(bravais,
	       dataset->std_lattice,
	       dataset->std_positions,
	       dataset->std_types);

  spg_free_dataset(dataset);

  primitive = spa_transform_to_primitive(bravais,
					 identity,
					 centering,
					 symprec);
  cel_free_cell(bravais);

  if (primitive == NULL) {
    goto err;
  }

  set_cell(lattice, position, types, primitive);
  num_prim_atom = primitive->size;

  cel_free_cell(primitive);

  return num_prim_atom;

 err:
  return 0;
}

static int standardize_cell(double lattice[3][3],
			    double position[][3],
			    int types[],
			    const int num_atom,
			    const double symprec)
{
  int i, n_std_atoms;
  SpglibDataset *dataset;

  n_std_atoms = 0;
  dataset = NULL;

  if ((dataset = get_dataset(lattice,
			     position,
			     types,
			     num_atom,
			     0,
			     symprec)) == NULL) {
    return 0;
  }

  n_std_atoms = dataset->n_std_atoms;

  mat_copy_matrix_d3(lattice, dataset->std_lattice);
  for (i = 0; i < dataset->n_std_atoms; i++) {
    types[i] = dataset->std_types[i];
    mat_copy_vector_d3(position[i], dataset->std_positions[i]);
  }

  spg_free_dataset(dataset);
  
  return n_std_atoms;
}

static int get_standardized_cell(double lattice[3][3],
				 double position[][3],
				 int types[],
				 const int num_atom,
				 const int to_primitive,
				 const double symprec)
{
  int num_std_atom;
  SpglibDataset *dataset;
  Cell *std_cell, *cell;
  Centering centering;

  num_std_atom = 0;
  dataset = NULL;
  std_cell = NULL;
  cell = NULL;

  if ((dataset = get_dataset(lattice,
			     position,
			     types,
			     num_atom,
			     0,
			     symprec)) == NULL) {
    goto err;
  }

  if (to_primitive) {
    if ((centering = get_centering(dataset->hall_number)) == CENTERING_ERROR) {
      goto err;
    }
    if (is_rhombohedral(dataset->hall_number)) {
      centering = R_CENTER;
    }
  } else {
    centering = PRIMITIVE;
  }

  if ((cell = cel_alloc_cell(num_atom)) == NULL) {
    spg_free_dataset(dataset);
    goto err;
  }
  
  cel_set_cell(cell, lattice, position, types);
  std_cell = spa_transform_to_primitive(cell,
					dataset->transformation_matrix,
					centering,
					symprec);
  spg_free_dataset(dataset);
  cel_free_cell(cell);

  if (std_cell == NULL) {
    goto err;
  }

  set_cell(lattice, position, types, std_cell);
  num_std_atom = std_cell->size;

  cel_free_cell(std_cell);

  return num_std_atom;

 err:
  return 0;
}

static void set_cell(double lattice[3][3],
		     double position[][3],
		     int types[],
		     Cell * cell)
{
  int i;

  mat_copy_matrix_d3(lattice, cell->lattice);
  for (i = 0; i < cell->size; i++) {
    types[i] = cell->types[i];
    mat_copy_vector_d3(position[i], cell->position[i]);
  }
}

static Centering get_centering(int hall_number)
{
  SpacegroupType spgtype;

  spgtype = spgdb_get_spacegroup_type(hall_number);

  return spgtype.centering;
}

static int get_international(char symbol[11],
			     SPGCONST double lattice[3][3],
			     SPGCONST double position[][3],
			     const int types[],
			     const int num_atom,
			     const double symprec)
{
  Cell *cell;
  Primitive *primitive;
  Spacegroup spacegroup;

  cell = NULL;
  primitive = NULL;
  spacegroup.number = 0;

  if ((cell = cel_alloc_cell(num_atom)) == NULL) {
    return 0;
  }

  cel_set_cell(cell, lattice, position, types);

  if ((primitive = spa_get_spacegroup(&spacegroup, cell, symprec)) != NULL) {
    prm_free_primitive(primitive);
    if (spacegroup.number > 0) {
      strcpy(symbol, spacegroup.international_short);
    }
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
  Primitive *primitive;
  Spacegroup spacegroup;

  cell = NULL;
  primitive = NULL;
  spacegroup.number = 0;

  if ((cell = cel_alloc_cell(num_atom)) == NULL) {
    return 0;
  }

  cel_set_cell(cell, lattice, position, types);

  if ((primitive = spa_get_spacegroup(&spacegroup, cell, symprec)) != NULL) {
    prm_free_primitive(primitive);
    if (spacegroup.number > 0) {
      strcpy(symbol, spacegroup.schoenflies);
    }
  }

  cel_free_cell(cell);

  return spacegroup.number;
}

/* Return 0 if failed */
static int get_symmetry_numerical(int rotation[][3][3],
				  double translation[][3],
				  const int max_size,
				  SPGCONST double lattice[3][3],
				  SPGCONST double position[][3],
				  const int types[],
				  const int num_atom,
				  const double symprec)
{
  int i, size;
  Cell *cell;
  Symmetry *symmetry;

  size = 0;
  cell = NULL;
  symmetry = NULL;

  if ((cell = cel_alloc_cell(num_atom)) == NULL) {
    return 0;
  }

  cel_set_cell(cell, lattice, position, types);

  if ((symmetry = sym_get_operation(cell, symprec)) == NULL) {
    cel_free_cell(cell);
    return 0;
  }

  if (symmetry->size > max_size) {
    fprintf(stderr, "spglib: Indicated max size(=%d) is less than number ",
	    max_size);
    fprintf(stderr, "spglib: of symmetry operations(=%d).\n", symmetry->size);
    goto ret;
  }

  for (i = 0; i < symmetry->size; i++) {
    mat_copy_matrix_i3(rotation[i], symmetry->rot[i]);
    mat_copy_vector_d3(translation[i], symmetry->trans[i]);
  }
  size = symmetry->size;

 ret:
  sym_free_symmetry(symmetry);
  cel_free_cell(cell);

  return size;
}

static int is_rhombohedral(const int hall_number)
{
  if (hall_number == 433 ||
      hall_number == 436 ||
      hall_number == 444 ||
      hall_number == 450 ||
      hall_number == 452 ||
      hall_number == 458 ||
      hall_number == 460) {
    return 1;
  } else {
    return 0;
  }
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
  MatINT *rotations, *rot_reciprocal;

  dataset = get_dataset(lattice,
			position,
			types,
			num_atom,
			0,
			symprec);

  if ((rotations = mat_alloc_MatINT(dataset->n_operations)) == NULL) {
    spg_free_dataset(dataset);
    return 0;
  }

  for (i = 0; i < dataset->n_operations; i++) {
    mat_copy_matrix_i3(rotations->mat[i], dataset->rotations[i]);
  }
  rot_reciprocal = kpt_get_point_group_reciprocal(rotations, is_time_reversal);
  num_ir = kpt_get_irreducible_reciprocal_mesh(grid_address,
					       map,
					       mesh,
					       is_shift,
					       rot_reciprocal);
  mat_free_MatINT(rot_reciprocal);
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
  
  rot_real = NULL;

  if ((rot_real = mat_alloc_MatINT(num_rot)) == NULL) {
    return 0;
  }

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
