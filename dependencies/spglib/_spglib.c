#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

static PyObject * get_dataset(PyObject *self, PyObject *args);
static PyObject * get_spacegroup(PyObject *self, PyObject *args);
static PyObject * get_pointgroup(PyObject *self, PyObject *args);
static PyObject * refine_cell(PyObject *self, PyObject *args);
static PyObject * get_symmetry(PyObject *self, PyObject *args);
static PyObject *
get_symmetry_with_collinear_spin(PyObject *self, PyObject *args);
static PyObject * find_primitive(PyObject *self, PyObject *args);
static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_grid_points_by_rotations(PyObject *self, PyObject *args);
static PyObject * get_BZ_grid_points_by_rotations(PyObject *self, PyObject *args);
static PyObject * relocate_BZ_grid_address(PyObject *self, PyObject *args);
static PyObject *
get_triplets_reciprocal_mesh_at_q(PyObject *self, PyObject *args);
static PyObject * get_BZ_triplets_at_q(PyObject *self, PyObject *args);
static PyObject * get_neighboring_grid_points(PyObject *self, PyObject *args);
static PyObject *
get_tetrahedra_relative_grid_address(PyObject *self, PyObject *args);
static PyObject *
get_all_tetrahedra_relative_grid_address(PyObject *self, PyObject *args);
static PyObject *
get_tetrahedra_integration_weight(PyObject *self, PyObject *args);
static PyObject *
get_tetrahedra_integration_weight_at_omegas(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"dataset", get_dataset, METH_VARARGS, "Dataset for crystal symmetry"},
  {"spacegroup", get_spacegroup, METH_VARARGS, "International symbol"},
  {"pointgroup", get_pointgroup, METH_VARARGS,
   "International symbol of pointgroup"},
  {"refine_cell", refine_cell, METH_VARARGS, "Refine cell"},
  {"symmetry", get_symmetry, METH_VARARGS, "Symmetry operations"},
  {"symmetry_with_collinear_spin", get_symmetry_with_collinear_spin,
   METH_VARARGS, "Symmetry operations with collinear spin magnetic moments"},
  {"primitive", find_primitive, METH_VARARGS,
   "Find primitive cell in the input cell"},
  {"ir_reciprocal_mesh", get_ir_reciprocal_mesh, METH_VARARGS,
   "Reciprocal mesh points with map"},
  {"stabilized_reciprocal_mesh", get_stabilized_reciprocal_mesh, METH_VARARGS,
   "Reciprocal mesh points with map"},
  {"grid_points_by_rotations", get_grid_points_by_rotations, METH_VARARGS,
   "Rotated grid points are returned"},
  {"BZ_grid_points_by_rotations", get_BZ_grid_points_by_rotations, METH_VARARGS,
   "Rotated grid points in BZ are returned"},
  {"BZ_grid_address", relocate_BZ_grid_address, METH_VARARGS,
   "Relocate grid addresses inside Brillouin zone"},
  {"triplets_reciprocal_mesh_at_q", get_triplets_reciprocal_mesh_at_q,
   METH_VARARGS, "Triplets on reciprocal mesh points at a specific q-point"},
  {"BZ_triplets_at_q", get_BZ_triplets_at_q,
   METH_VARARGS, "Triplets in reciprocal primitive lattice are transformed to those in BZ."},
  {"neighboring_grid_points", get_neighboring_grid_points,
   METH_VARARGS, "Neighboring grid points by relative grid addresses"},
  {"tetrahedra_relative_grid_address", get_tetrahedra_relative_grid_address,
   METH_VARARGS, "Relative grid addresses of vertices of 24 tetrahedra"},
  {"all_tetrahedra_relative_grid_address",
   get_all_tetrahedra_relative_grid_address,
   METH_VARARGS,
   "4 (all) sets of relative grid addresses of vertices of 24 tetrahedra"},
  {"tetrahedra_integration_weight", get_tetrahedra_integration_weight,
   METH_VARARGS, "Integration weight for tetrahedron method"},
  {"tetrahedra_integration_weight_at_omegas",
   get_tetrahedra_integration_weight_at_omegas,
   METH_VARARGS, "Integration weight for tetrahedron method at omegas"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_spglib(void)
{
  Py_InitModule3("_spglib", functions, "C-extension for spglib\n\n...\n");
  return;
}

static PyObject * get_dataset(PyObject *self, PyObject *args)
{
  int i, j, k;
  double symprec, angle_tolerance;
  SpglibDataset *dataset;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  PyObject* array, *vec, *mat, *rot, *trans, *wyckoffs, *equiv_atoms;
  
  if (!PyArg_ParseTuple(args, "OOOdd",
			&lattice,
			&position,
			&atom_type,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  const int num_atom = position->dimensions[0];
  const int* typat = (int*)atom_type->data;

  dataset = spgat_get_dataset(lat,
			      pos,
			      typat,
			      num_atom,
			      symprec,
			      angle_tolerance);

  array = PyList_New(9);

  /* Space group number, international symbol, hall symbol */
  PyList_SetItem(array, 0, PyInt_FromLong((long) dataset->spacegroup_number));
  PyList_SetItem(array, 1, PyString_FromString(dataset->international_symbol));
  PyList_SetItem(array, 2, PyString_FromString(dataset->hall_symbol));

  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }
  PyList_SetItem(array, 3, mat);

  /* Origin shift */
  vec = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(vec, i, PyFloat_FromDouble(dataset->origin_shift[i]));
  }
  PyList_SetItem(array, 4, vec);

  /* Rotation matrices */
  rot = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    mat = PyList_New(3);
    for (j = 0; j < 3; j++) {
      vec = PyList_New(3);
      for (k = 0; k < 3; k++) {
	PyList_SetItem(vec, k, PyInt_FromLong((long) dataset->rotations[i][j][k]));
      }
      PyList_SetItem(mat, j, vec);
    }
    PyList_SetItem(rot, i, mat);
  }
  PyList_SetItem(array, 5, rot);

  /* Translation vectors */
  trans = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->translations[i][j]));
    }
    PyList_SetItem(trans, i, vec);
  }
  PyList_SetItem(array, 6, trans);

  /* Wyckoff letters, Equivalent atoms */
  wyckoffs = PyList_New(dataset->n_atoms);
  equiv_atoms = PyList_New(dataset->n_atoms);
  for (i = 0; i < dataset->n_atoms; i++) {
    PyList_SetItem(wyckoffs, i, PyInt_FromLong((long) dataset->wyckoffs[i]));
    PyList_SetItem(equiv_atoms, i, PyInt_FromLong((long) dataset->equivalent_atoms[i]));
  }
  PyList_SetItem(array, 7, wyckoffs);
  PyList_SetItem(array, 8, equiv_atoms);
  spg_free_dataset(dataset);

  return array;
}

static PyObject * get_spacegroup(PyObject *self, PyObject *args)
{
  int i;
  double symprec, angle_tolerance;
  char symbol_with_number[17], spg_symbol[11];
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOdd",
			&lattice,
			&position,
			&atom_type,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  const int num_atom = position->dimensions[0];
  const int* typat = (int*)atom_type->data;

  const int num_spg = spgat_get_international(spg_symbol,
					      lat,
					      pos,
					      typat,
					      num_atom,
					      symprec,
					      angle_tolerance);
  for (i = 9; i > 0; i--) {
    if (! isspace(spg_symbol[i])) { break; }
  }
  spg_symbol[i + 1] = 0;
  sprintf(symbol_with_number, "%s (%d)", spg_symbol, num_spg);

  return PyString_FromString(symbol_with_number);
}

static PyObject * get_pointgroup(PyObject *self, PyObject *args)
{
  PyArrayObject* rotations;
  if (! PyArg_ParseTuple(args, "O", &rotations)) {
    return NULL;
  }

  int i, j;
  int trans_mat[3][3];
  char symbol[6];
  PyObject* array, * mat, * vec;

  SPGCONST int(*rot)[3][3] = (int(*)[3][3])rotations->data;
  const int num_rot = rotations->dimensions[0];
  const int ptg_num = spg_get_pointgroup(symbol, trans_mat, rot, num_rot);

  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyInt_FromLong((long)trans_mat[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }

  array = PyList_New(3);
  PyList_SetItem(array, 0, PyString_FromString(symbol));
  PyList_SetItem(array, 1, PyInt_FromLong((long) ptg_num));
  PyList_SetItem(array, 2, mat);

  return array;
}

static PyObject * refine_cell(PyObject *self, PyObject *args)
{
  int num_atom;
  double symprec, angle_tolerance;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOidd",
			&lattice,
			&position,
			&atom_type,
			&num_atom,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  int* typat = (int*)atom_type->data;

  int num_atom_brv = spgat_refine_cell(lat,
				       pos,
				       typat,
				       num_atom,
				       symprec,
				       angle_tolerance);

  return PyInt_FromLong((long) num_atom_brv);
}


static PyObject * find_primitive(PyObject *self, PyObject *args)
{
  double symprec, angle_tolerance;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOdd",
			&lattice,
			&position,
			&atom_type,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  double (*lat)[3] = (double(*)[3])lattice->data;
  double (*pos)[3] = (double(*)[3])position->data;
  int num_atom = position->dimensions[0];
  int* types = (int*)atom_type->data;

  int num_atom_prim = spgat_find_primitive(lat,
					   pos,
					   types,
					   num_atom,
					   symprec,
					   angle_tolerance);

  return PyInt_FromLong((long) num_atom_prim);
}

static PyObject * get_symmetry(PyObject *self, PyObject *args)
{
  double symprec, angle_tolerance;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* rotation;
  PyArrayObject* translation;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOOdd",
			&rotation,
			&translation,
			&lattice,
			&position,
			&atom_type,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  const int* types = (int*)atom_type->data;
  const int num_atom = position->dimensions[0];
  int (*rot)[3][3] = (int(*)[3][3])rotation->data;
  double (*trans)[3] = (double(*)[3])translation->data;
  const int num_sym_from_array_size = rotation->dimensions[0];
  
  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_sym = spgat_get_symmetry(rot,
					 trans,
					 num_sym_from_array_size,
					 lat,
					 pos,
					 types,
					 num_atom,
					 symprec,
					 angle_tolerance);
  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_symmetry_with_collinear_spin(PyObject *self,
						   PyObject *args)
{
  double symprec, angle_tolerance;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* rotation;
  PyArrayObject* translation;
  PyArrayObject* atom_type;
  PyArrayObject* magmom;
  if (!PyArg_ParseTuple(args, "OOOOOOdd",
			&rotation,
			&translation,
			&lattice,
			&position,
			&atom_type,
			&magmom,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  const double *spins = (double*) magmom->data;
  const int* types = (int*)atom_type->data;
  const int num_atom = position->dimensions[0];
  int (*rot)[3][3] = (int(*)[3][3])rotation->data;
  double (*trans)[3] = (double(*)[3])translation->data;
  const int num_sym_from_array_size = rotation->dimensions[0];

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_sym = 
    spgat_get_symmetry_with_collinear_spin(rot,
					   trans,
					   num_sym_from_array_size,
					   lat,
					   pos,
					   types,
					   spins,
					   num_atom,
					   symprec,
					   angle_tolerance);
  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args)
{
  double symprec;
  PyArrayObject* grid_address_py;
  PyArrayObject* map;
  PyArrayObject* mesh;
  PyArrayObject* is_shift;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOiOOOd",
			&grid_address_py,
			&map,
			&mesh,
			&is_shift,
			&is_time_reversal,
			&lattice,
			&position,
			&atom_type,
			&symprec))
    return NULL;

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  const int* types = (int*)atom_type->data;
  const int* mesh_int = (int*)mesh->data;
  const int* is_shift_int = (int*)is_shift->data;
  const int num_atom = position->dimensions[0];
  int (*grid_address)[3] = (int(*)[3])grid_address_py->data;
  int *map_int = (int*)map->data;

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_ir = spg_get_ir_reciprocal_mesh(grid_address,
						map_int,
						mesh_int,
						is_shift_int,
						is_time_reversal,
						lat,
						pos,
						types,
						num_atom,
						symprec);
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args)
{
  PyArrayObject* grid_address_py;
  PyArrayObject* map;
  PyArrayObject* mesh;
  PyArrayObject* is_shift;
  int is_time_reversal;
  PyArrayObject* rotations;
  PyArrayObject* qpoints;
  if (!PyArg_ParseTuple(args, "OOOOiOO",
			&grid_address_py,
			&map,
			&mesh,
			&is_shift,
			&is_time_reversal,
			&rotations,
			&qpoints)) {
    return NULL;
  }

  int (*grid_address)[3] = (int(*)[3])grid_address_py->data;
  int *map_int = (int*)map->data;
  const int* mesh_int = (int*)mesh->data;
  const int* is_shift_int = (int*)is_shift->data;
  SPGCONST int (*rot)[3][3] = (int(*)[3][3])rotations->data;
  const int num_rot = rotations->dimensions[0];
  SPGCONST double (*q)[3] = (double(*)[3])qpoints->data;
  const int num_q = qpoints->dimensions[0];

  const int num_ir = spg_get_stabilized_reciprocal_mesh(grid_address,
							map_int,
							mesh_int,
							is_shift_int,
							is_time_reversal,
							num_rot,
							rot,
							num_q,
							q);
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * get_grid_points_by_rotations(PyObject *self, PyObject *args)
{
  PyArrayObject* rot_grid_points_py;
  PyArrayObject* address_orig_py;
  PyArrayObject* rot_reciprocal_py;
  PyArrayObject* mesh_py;
  PyArrayObject* is_shift_py;
  if (!PyArg_ParseTuple(args, "OOOOO",
			&rot_grid_points_py,
			&address_orig_py,
			&rot_reciprocal_py,
			&mesh_py,
			&is_shift_py)) {
    return NULL;
  }

  int *rot_grid_points = (int*)rot_grid_points_py->data;
  const int *address_orig = (int*)address_orig_py->data;
  SPGCONST int (*rot_reciprocal)[3][3] = (int(*)[3][3])rot_reciprocal_py->data;
  const int num_rot = rot_reciprocal_py->dimensions[0];
  const int* mesh = (int*)mesh_py->data;
  const int* is_shift = (int*)is_shift_py->data;
  
  spg_get_grid_points_by_rotations(rot_grid_points,
				   address_orig,
				   num_rot,
				   rot_reciprocal,
				   mesh,
				   is_shift);
  Py_RETURN_NONE;
}

static PyObject * get_BZ_grid_points_by_rotations(PyObject *self, PyObject *args)
{
  PyArrayObject* rot_grid_points_py;
  PyArrayObject* address_orig_py;
  PyArrayObject* rot_reciprocal_py;
  PyArrayObject* mesh_py;
  PyArrayObject* is_shift_py;
  PyArrayObject* bz_map_py;
  if (!PyArg_ParseTuple(args, "OOOOOO",
			&rot_grid_points_py,
			&address_orig_py,
			&rot_reciprocal_py,
			&mesh_py,
			&is_shift_py,
			&bz_map_py)) {
    return NULL;
  }

  int *rot_grid_points = (int*)rot_grid_points_py->data;
  const int *address_orig = (int*)address_orig_py->data;
  SPGCONST int (*rot_reciprocal)[3][3] = (int(*)[3][3])rot_reciprocal_py->data;
  const int num_rot = rot_reciprocal_py->dimensions[0];
  const int* mesh = (int*)mesh_py->data;
  const int* is_shift = (int*)is_shift_py->data;
  const int* bz_map = (int*)bz_map_py->data;
  
  spg_get_BZ_grid_points_by_rotations(rot_grid_points,
				      address_orig,
				      num_rot,
				      rot_reciprocal,
				      mesh,
				      is_shift,
				      bz_map);
  Py_RETURN_NONE;
}

static PyObject * relocate_BZ_grid_address(PyObject *self, PyObject *args)
{
  PyArrayObject* bz_grid_address_py;
  PyArrayObject* bz_map_py;
  PyArrayObject* grid_address_py;
  PyArrayObject* mesh_py;
  PyArrayObject* is_shift_py;
  PyArrayObject* reciprocal_lattice_py;
  if (!PyArg_ParseTuple(args, "OOOOOO",
			&bz_grid_address_py,
			&bz_map_py,
			&grid_address_py,
			&mesh_py,
			&reciprocal_lattice_py,
			&is_shift_py)) {
    return NULL;
  }

  int (*bz_grid_address)[3] = (int(*)[3])bz_grid_address_py->data;
  int *bz_map = (int*)bz_map_py->data;
  SPGCONST int (*grid_address)[3] = (int(*)[3])grid_address_py->data;
  const int* mesh = (int*)mesh_py->data;
  const int* is_shift = (int*)is_shift_py->data;
  SPGCONST double (*reciprocal_lattice)[3]  =
    (double(*)[3])reciprocal_lattice_py->data;
  int num_ir_gp;
  
  num_ir_gp = spg_relocate_BZ_grid_address(bz_grid_address,
					   bz_map,
					   grid_address,
					   mesh,
					   reciprocal_lattice,
					   is_shift);

  return PyInt_FromLong((long) num_ir_gp);
}

static PyObject * get_triplets_reciprocal_mesh_at_q(PyObject *self, PyObject *args)
{
  PyArrayObject* map_triplets;
  PyArrayObject* grid_address_py;
  PyArrayObject* map_q;
  int fixed_grid_number;
  PyArrayObject* mesh;
  int is_time_reversal;
  PyArrayObject* rotations;
  if (!PyArg_ParseTuple(args, "OOOiOiO",
			&map_triplets,
			&map_q,
			&grid_address_py,
			&fixed_grid_number,
			&mesh,
			&is_time_reversal,
			&rotations)) {
    return NULL;
  }

  int (*grid_address)[3] = (int(*)[3])grid_address_py->data;
  int *map_triplets_int = (int*)map_triplets->data;
  int *map_q_int = (int*)map_q->data;

  const int* mesh_int = (int*)mesh->data;
  SPGCONST int (*rot)[3][3] = (int(*)[3][3])rotations->data;
  const int num_rot = rotations->dimensions[0];
  const int num_ir = 
    spg_get_triplets_reciprocal_mesh_at_q(map_triplets_int,
					  map_q_int,
					  grid_address,
					  fixed_grid_number,
					  mesh_int,
					  is_time_reversal,
					  num_rot,
					  rot);

  return PyInt_FromLong((long) num_ir);
}


static PyObject * get_BZ_triplets_at_q(PyObject *self, PyObject *args)
{
  PyArrayObject* triplets_py;
  PyArrayObject* bz_grid_address_py;
  PyArrayObject* bz_map_py;
  PyArrayObject* map_triplets_py;
  PyArrayObject* mesh_py;
  int grid_point;
  if (!PyArg_ParseTuple(args, "OiOOOO",
			&triplets_py,
			&grid_point,
			&bz_grid_address_py,
			&bz_map_py,
			&map_triplets_py,
			&mesh_py)) {
    return NULL;
  }

  int (*triplets)[3] = (int(*)[3])triplets_py->data;
  SPGCONST int (*bz_grid_address)[3] = (int(*)[3])bz_grid_address_py->data;
  const int *bz_map = (int*)bz_map_py->data;
  const int *map_triplets = (int*)map_triplets_py->data;
  const int num_map_triplets = (int)map_triplets_py->dimensions[0];
  const int *mesh = (int*)mesh_py->data;
  int num_ir;
  
  num_ir = spg_get_BZ_triplets_at_q(triplets,
				    grid_point,
				    bz_grid_address,
				    bz_map,
				    map_triplets,
				    num_map_triplets,
				    mesh);

  return PyInt_FromLong((long) num_ir);
}

static PyObject *get_neighboring_grid_points(PyObject *self, PyObject *args)
{
  PyArrayObject* relative_grid_points_py;
  PyArrayObject* relative_grid_address_py;
  PyArrayObject* mesh_py;
  PyArrayObject* bz_grid_address_py;
  PyArrayObject* bz_map_py;
  int grid_point;
  if (!PyArg_ParseTuple(args, "OiOOOO",
			&relative_grid_points_py,
			&grid_point,
			&relative_grid_address_py,
			&mesh_py,
			&bz_grid_address_py,
			&bz_map_py)) {
    return NULL;
  }

  int* relative_grid_points = (int*)relative_grid_points_py->data;
  SPGCONST int (*relative_grid_address)[3] =
    (int(*)[3])relative_grid_address_py->data;
  const int num_relative_grid_address = relative_grid_address_py->dimensions[0];
  const int *mesh = (int*)mesh_py->data;
  SPGCONST int (*bz_grid_address)[3] = (int(*)[3])bz_grid_address_py->data;
  const int *bz_map = (int*)bz_map_py->data;
  
  spg_get_neighboring_grid_points(relative_grid_points,
				  grid_point,
				  relative_grid_address,
				  num_relative_grid_address,
				  mesh,
				  bz_grid_address,
				  bz_map);
  Py_RETURN_NONE;
}

static PyObject *
get_tetrahedra_relative_grid_address(PyObject *self, PyObject *args)
{
  PyArrayObject* relative_grid_address_py;
  PyArrayObject* reciprocal_lattice_py;
  
  if (!PyArg_ParseTuple(args, "OO",
			&relative_grid_address_py,
			&reciprocal_lattice_py)) {
    return NULL;
  }

  int (*relative_grid_address)[4][3] =
    (int(*)[4][3])relative_grid_address_py->data;
  SPGCONST double (*reciprocal_lattice)[3] =
    (double(*)[3])reciprocal_lattice_py->data;

  spg_get_tetrahedra_relative_grid_address(relative_grid_address,
					   reciprocal_lattice);

  Py_RETURN_NONE;
}

static PyObject *
get_all_tetrahedra_relative_grid_address(PyObject *self, PyObject *args)
{
  PyArrayObject* relative_grid_address_py;
  
  if (!PyArg_ParseTuple(args, "O",
			&relative_grid_address_py)) {
    return NULL;
  }

  int (*relative_grid_address)[24][4][3] =
    (int(*)[24][4][3])relative_grid_address_py->data;

  spg_get_all_tetrahedra_relative_grid_address(relative_grid_address);

  Py_RETURN_NONE;
}

static PyObject *
get_tetrahedra_integration_weight(PyObject *self, PyObject *args)
{
  double omega;
  PyArrayObject* tetrahedra_omegas_py;
  char function;
  if (!PyArg_ParseTuple(args, "dOc",
			&omega,
			&tetrahedra_omegas_py,
			&function)) {
    return NULL;
  }

  SPGCONST double (*tetrahedra_omegas)[4] =
    (double(*)[4])tetrahedra_omegas_py->data;
  
  double iw = spg_get_tetrahedra_integration_weight(omega,
						    tetrahedra_omegas,
						    function);

  return PyFloat_FromDouble(iw);
}

static PyObject *
get_tetrahedra_integration_weight_at_omegas(PyObject *self, PyObject *args)
{
  PyArrayObject* integration_weights_py;
  PyArrayObject* omegas_py;
  PyArrayObject* tetrahedra_omegas_py;
  char function;
  if (!PyArg_ParseTuple(args, "OOOc",
			&integration_weights_py,
			&omegas_py,
			&tetrahedra_omegas_py,
			&function)) {
    return NULL;
  }

  const double *omegas = (double*)omegas_py->data;
  double *iw = (double*)integration_weights_py->data;
  const int num_omegas = (int)omegas_py->dimensions[0];
  SPGCONST double (*tetrahedra_omegas)[4] =
    (double(*)[4])tetrahedra_omegas_py->data;
  
  spg_get_tetrahedra_integration_weight_at_omegas(iw,
						  num_omegas,
						  omegas,
						  tetrahedra_omegas,
						  function);

  Py_RETURN_NONE;
}
