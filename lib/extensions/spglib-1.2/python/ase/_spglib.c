#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

static PyObject * get_dataset(PyObject *self, PyObject *args);
static PyObject * get_spacegroup(PyObject *self, PyObject *args);
static PyObject * get_pointgroup(PyObject *self, PyObject *args);
static PyObject * refine_cell(PyObject *self, PyObject *args);
static PyObject * get_symmetry(PyObject *self, PyObject *args);
static PyObject * get_symmetry_with_collinear_spin(PyObject *self, PyObject *args);
static PyObject * find_primitive(PyObject *self, PyObject *args);
static PyObject * get_ir_kpoints(PyObject *self, PyObject *args);
static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_triplets_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_triplets_reciprocal_mesh_at_q(PyObject *self, PyObject *args);
static PyObject * extract_triplets_reciprocal_mesh_at_q(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"dataset", get_dataset, METH_VARARGS,
   "Dataset for crystal symmetry"},
  {"spacegroup", get_spacegroup, METH_VARARGS,
   "International symbol"},
  {"pointgroup", get_pointgroup, METH_VARARGS,
   "International symbol of pointgroup"},
  {"refine_cell", refine_cell, METH_VARARGS,
   "Refine cell"},
  {"symmetry", get_symmetry, METH_VARARGS,
   "Symmetry operations"},
  {"symmetry_with_collinear_spin", get_symmetry_with_collinear_spin,
   METH_VARARGS, "Symmetry operations with collinear spin magnetic moments"},
  {"primitive", find_primitive, METH_VARARGS,
   "Find primitive cell in the input cell"},
  {"ir_kpoints", get_ir_kpoints, METH_VARARGS,
   "Irreducible k-points"},
  {"ir_reciprocal_mesh", get_ir_reciprocal_mesh, METH_VARARGS,
   "Reciprocal mesh points with map"},
  {"stabilized_reciprocal_mesh", get_stabilized_reciprocal_mesh, METH_VARARGS,
   "Reciprocal mesh points with map"},
  {"triplets_reciprocal_mesh", get_triplets_reciprocal_mesh, METH_VARARGS,
   "Triplets on reciprocal mesh points"},
  {"triplets_reciprocal_mesh_at_q", get_triplets_reciprocal_mesh_at_q,
   METH_VARARGS, "Triplets on reciprocal mesh points at a specific q-point"},
  {"triplets_reciprocal_mesh_at_q_from_triplets",
   extract_triplets_reciprocal_mesh_at_q, METH_VARARGS,
   "Triplets on reciprocal mesh points at a specific q-point extracted from full triplets"},
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
  const long* typat_long = (long*)atom_type->data;

  int typat[atom_type->dimensions[0]];
  for (i = 0; i < num_atom; i++) {
    typat[i] = (int)typat_long[i];
  }

  dataset = spgat_get_dataset(lat,
			      pos,
			      typat,
			      num_atom,
			      symprec,
			      angle_tolerance);

  array = PyList_New(0);

  /* Space group number, international symbol, hall symbol */
  PyList_Append(array, PyInt_FromLong((long) dataset->spacegroup_number));
  PyList_Append(array, PyString_FromString(dataset->international_symbol));
  PyList_Append(array, PyString_FromString(dataset->hall_symbol));

  /* Transformation matrix */
  mat = PyList_New(0);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(0);
    for (j = 0; j < 3; j++) {
      PyList_Append(vec, PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
    }
    PyList_Append(mat, vec);
  }
  PyList_Append(array, mat);

  /* Origin shift */
  vec = PyList_New(0);
  for (i = 0; i < 3; i++) {
    PyList_Append(vec, PyFloat_FromDouble(dataset->origin_shift[i]));
  }
  PyList_Append(array, vec);

  /* Rotation matrices */
  rot = PyList_New(0);
  for (i = 0; i < dataset->n_operations; i++) {
    mat = PyList_New(0);
    for (j = 0; j < 3; j++) {
      vec = PyList_New(0);
      for (k = 0; k < 3; k++) {
	PyList_Append(vec, PyInt_FromLong((long) dataset->rotations[i][j][k]));
      }
      PyList_Append(mat, vec);
    }
    PyList_Append(rot, mat);
  }
  PyList_Append(array, rot);

  /* Translation vectors */
  trans = PyList_New(0);
  for (i = 0; i < dataset->n_operations; i++) {
    vec = PyList_New(0);
    for (j = 0; j < 3; j++) {
      PyList_Append(vec, PyFloat_FromDouble(dataset->translations[i][j]));
    }
    PyList_Append(trans, vec);
  }
  PyList_Append(array, trans);

  /* Wyckoff letters, Equivalent atoms */
  wyckoffs = PyList_New(0);
  equiv_atoms = PyList_New(0);
  for (i = 0; i < dataset->n_atoms; i++) {
    PyList_Append(wyckoffs, PyInt_FromLong((long) dataset->wyckoffs[i]));
    PyList_Append(equiv_atoms, PyInt_FromLong((long) dataset->equivalent_atoms[i]));
  }
  PyList_Append(array, wyckoffs);
  PyList_Append(array, equiv_atoms);

  return array;
}

static PyObject * get_spacegroup(PyObject *self, PyObject *args)
{
  int i;
  double symprec, angle_tolerance;
  char symbol[26];
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
  const long* typat_long = (long*)atom_type->data;

  int typat[num_atom];
  for (i = 0; i < num_atom; i++)
    typat[i] = (int)typat_long[i];

  const int num_spg = spgat_get_international(symbol,
					      lat,
					      pos,
					      typat,
					      num_atom,
					      symprec,
					      angle_tolerance);
  sprintf(symbol, "%s (%d)", symbol, num_spg);

  return PyString_FromString(symbol);
}

static PyObject * get_pointgroup(PyObject *self, PyObject *args)
{
  PyArrayObject* rotations;
  if (! PyArg_ParseTuple(args, "O", &rotations)) {
    return NULL;
  }

  long *rot_long = (long*)rotations->data;

  int i, j, k;
  int trans_mat[3][3];
  char symbol[6];
  PyObject* array, * mat, * vec;
    
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for (i = 0; i < num_rot; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];
      }
    }
  }

  const int ptg_num = spg_get_pointgroup(symbol, trans_mat, rot, num_rot);

  /* Transformation matrix */
  mat = PyList_New(0);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(0);
    for (j = 0; j < 3; j++) {
      PyList_Append(vec, PyInt_FromLong((long)trans_mat[i][j]));
    }
    PyList_Append(mat, vec);
  }

  array = PyList_New(0);
  PyList_Append(array, PyString_FromString(symbol));
  PyList_Append(array, PyInt_FromLong((long) ptg_num));
  PyList_Append(array, mat);

  return array;
}

static PyObject * refine_cell(PyObject *self, PyObject *args)
{
  int i, num_atom;
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
  long* typat_long = (long*)atom_type->data;

  int typat[atom_type->dimensions[0]];
  for (i = 0; i < num_atom; i++) {
    typat[i] = (int)typat_long[i];
  }

  int num_atom_brv = spgat_refine_cell(lat,
				       pos,
				       typat,
				       num_atom,
				       symprec,
				       angle_tolerance);

  for (i = 0; i < num_atom_brv; i++) {
    typat_long[i] = typat[i];
  }

  return PyInt_FromLong((long) num_atom_brv);
}


static PyObject * find_primitive(PyObject *self, PyObject *args)
{
  int i;
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
  long* types_long = (long*)atom_type->data;

  int types[num_atom];
  for (i = 0; i < num_atom; i++) {
    types[i] = (int)types_long[i];
  }
  
  int num_atom_prim = spgat_find_primitive(lat,
					   pos,
					   types,
					   num_atom,
					   symprec,
					   angle_tolerance);

  for (i = 0; i < num_atom_prim; i++) {
    types_long[i] = types[i];
  }

  return PyInt_FromLong((long) num_atom_prim);
}

static PyObject * get_symmetry(PyObject *self, PyObject *args)
{
  int i, j, k;
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
  const long* types_long = (long*)atom_type->data;
  const int num_atom = position->dimensions[0];
  long *rot_long = (long*)rotation->data;
  double (*trans)[3] = (double(*)[3])translation->data;
  const int num_sym_from_array_size = rotation->dimensions[0];

  int rot[num_sym_from_array_size][3][3];
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++) {
    types[i] = (int)types_long[i];
  }
  
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
  for (i = 0; i < num_sym; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot_long[i*9+j*3+k] = (long)rot[i][j][k];
      }
    }
  }

  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_symmetry_with_collinear_spin(PyObject *self,
						   PyObject *args)
{
  int i, j, k;
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
  const long* types_long = (long*)atom_type->data;
  const int num_atom = position->dimensions[0];
  long *rot_long = (long*)rotation->data;
  double (*trans)[3] = (double(*)[3])translation->data;
  const int num_sym_from_array_size = rotation->dimensions[0];

  int rot[num_sym_from_array_size][3][3];
  int types[num_atom];
  for (i = 0; i < num_atom; i++) {
    types[i] = (int)types_long[i];
  }
  
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
  for (i = 0; i < num_sym; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot_long[i*9+j*3+k] = (long)rot[i][j][k];
      }
    }
  }

  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_ir_kpoints(PyObject *self, PyObject *args)
{
  int i;
  double symprec;
  int is_time_reversal;
  PyArrayObject* kpoint;
  PyArrayObject* kpoint_map;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOOid", &kpoint_map, &kpoint, &lattice, &position,
			&atom_type, &is_time_reversal, &symprec))
    return NULL;

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  SPGCONST double (*pos)[3] = (double(*)[3])position->data;
  SPGCONST double (*kpts)[3] = (double(*)[3])kpoint->data;
  const int num_kpoint = kpoint->dimensions[0];
  const long* types_long = (long*)atom_type->data;
  const int num_atom = position->dimensions[0];
  long *map_long = (long*)kpoint_map->data;
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];

  int map[num_kpoint];

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_ir_kpt = spg_get_ir_kpoints(map,
					    kpts,
					    num_kpoint,
					    lat,
					    pos,
					    types,
					    num_atom,
					    is_time_reversal,
					    symprec);

  for (i = 0; i < num_kpoint; i++)
    map_long[i] = (long) map[i];

  return PyInt_FromLong((long) num_ir_kpt);
}

static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args)
{
  int i, j;
  double symprec;
  PyArrayObject* grid_point;
  PyArrayObject* map;
  PyArrayObject* mesh;
  PyArrayObject* is_shift;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOiOOOd",
			&grid_point,
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
  const int num_grid = grid_point->dimensions[0];
  const long* types_long = (long*)atom_type->data;
  const long* mesh_long = (long*)mesh->data;
  const long* is_shift_long = (long*)is_shift->data;
  const int num_atom = position->dimensions[0];
  long *grid_long = (long*)grid_point->data;
  int grid_int[num_grid][3];
  long *map_long = (long*)map->data;
  int map_int[num_grid];
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++) {
    types[i] = (int)types_long[i];
  }

  int mesh_int[3];
  int is_shift_int[3];
  for (i = 0; i < 3; i++) {
    mesh_int[i] = (int) mesh_long[i];
    is_shift_int[i] = (int) is_shift_long[i];
  }  

  /* Check memory space */
  if (mesh_int[0]*mesh_int[1]*mesh_int[2] > num_grid) {
    return NULL;
  }

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_ir = spg_get_ir_reciprocal_mesh(grid_int,
						map_int,
						mesh_int,
						is_shift_int,
						is_time_reversal,
						lat,
						pos,
						types,
						num_atom,
						symprec);
  
  for (i = 0; i < mesh_int[0] * mesh_int[1] * mesh_int[2]; i++) {
    for (j = 0; j < 3; j++)
      grid_long[ i*3 + j ] = (long) grid_int[i][j];
    map_long[i] = (long) map_int[i];
  }
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args)
{
  int i, j, k;
  PyArrayObject* grid_point;
  PyArrayObject* map;
  PyArrayObject* mesh;
  PyArrayObject* is_shift;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* rotations;
  PyArrayObject* qpoints;
  double symprec;
  if (!PyArg_ParseTuple(args, "OOOOiOOOd",
			&grid_point,
			&map,
			&mesh,
			&is_shift,
			&is_time_reversal,
			&lattice,
			&rotations,
			&qpoints,
			&symprec))
    return NULL;

  long *grid_long = (long*)grid_point->data;
  const int num_grid = grid_point->dimensions[0];
  int grid_int[num_grid][3];

  long *map_long = (long*)map->data;
  int map_int[num_grid];

  int mesh_int[3];
  int is_shift_int[3];
  const long* mesh_long = (long*)mesh->data;
  const long* is_shift_long = (long*)is_shift->data;
  for (i = 0; i < 3; i++) {
    mesh_int[i] = (int) mesh_long[i];
    is_shift_int[i] = (int) is_shift_long[i];
  }  

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  const long* rot_long = (long*)rotations->data;
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for (i = 0; i < num_rot; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];
      }
    }
  }

  SPGCONST double (*q)[3] = (double(*)[3])qpoints->data;
  const int num_q = qpoints->dimensions[0];

  /* Check memory space */
  if (mesh_int[0]*mesh_int[1]*mesh_int[2] > num_grid) {
    return NULL;
  }

  const int num_ir = spg_get_stabilized_reciprocal_mesh(grid_int,
							map_int,
							mesh_int,
							is_shift_int,
							is_time_reversal,
							lat,
							num_rot,
							rot,
							num_q,
							q,
							symprec);
  
  for (i = 0; i < mesh_int[0] * mesh_int[1] * mesh_int[2]; i++) {
    for (j = 0; j < 3; j++) {
      grid_long[ i*3 + j ] = (long) grid_int[i][j];
    }
    map_long[i] = (long) map_int[i];
  }
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * get_triplets_reciprocal_mesh(PyObject *self, PyObject *args)
{
  PyArrayObject* mesh;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* rotations;
  double symprec;
  if (!PyArg_ParseTuple(args, "OiOOd",
			&mesh,
			&is_time_reversal,
			&lattice,
			&rotations,
			&symprec))
    return NULL;

  int i, j, k;
  PyObject * triplets, * weights, *tp, *ret_array, *mesh_points;

  int mesh_int[3];
  const long* mesh_long = (long*)mesh->data;
  for (i = 0; i < 3; i++) {
    mesh_int[i] = (int) mesh_long[i];
  }  
  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  const long* rot_long = (long*)rotations->data;
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for (i = 0; i < num_rot; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];
      }
    }
  }

  SpglibTriplets * spg_triplets = \
    spg_get_triplets_reciprocal_mesh(mesh_int,
				     is_time_reversal,
				     lat,
				     num_rot,
				     rot,
				     symprec);

  ret_array = PyList_New(0);
  triplets = PyList_New(0);
  weights = PyList_New(0);
  mesh_points = PyList_New(0);
  
  for (i = 0; i < spg_triplets->size; i++) {
    tp = PyList_New(0);
    for (j = 0; j < 3; j++) {
      PyList_Append(tp, PyInt_FromLong((long) spg_triplets->triplets[i][j]));
    }
    PyList_Append(triplets, tp);
    PyList_Append(weights, PyInt_FromLong((long) spg_triplets->weights[i]));
  }

  for (i = 0; i < mesh_int[0]*mesh_int[1]*mesh_int[2]; i++) {
    tp = PyList_New(0);
    for (j = 0; j < 3; j++) {
      PyList_Append(tp, PyInt_FromLong((long) spg_triplets->mesh_points[i][j]));
    }
    PyList_Append(mesh_points, tp);
  }

  PyList_Append(ret_array, triplets);
  PyList_Append(ret_array, weights);
  PyList_Append(ret_array, mesh_points);

  spg_free_triplets(spg_triplets);

  return ret_array;
}

static PyObject * get_triplets_reciprocal_mesh_at_q(PyObject *self, PyObject *args)
{
  PyArrayObject* weights;
  PyArrayObject* grid_points;
  PyArrayObject* third_q;
  int fixed_grid_number;
  PyArrayObject* mesh;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* rotations;
  double symprec;
  if (!PyArg_ParseTuple(args, "OOOiOiOOd",
			&weights,
			&grid_points,
			&third_q,
			&fixed_grid_number,
			&mesh,
			&is_time_reversal,
			&lattice,
			&rotations,
			&symprec))
    return NULL;

  int i, j, k;

  const int num_grid = grid_points->dimensions[0];
  long *grid_points_long = (long*)grid_points->data;
  int grid_points_int[num_grid][3];
  long *weights_long = (long*)weights->data;
  int weights_int[num_grid];
  long *third_q_long = (long*)third_q->data;
  int third_q_int[num_grid];

  int mesh_int[3];
  const long* mesh_long = (long*)mesh->data;
  for (i = 0; i < 3; i++) {
    mesh_int[i] = (int) mesh_long[i];
  }  

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  const long* rot_long = (long*)rotations->data;
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for (i = 0; i < num_rot; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];
      }
    }
  }

  const int num_ir = 
    spg_get_triplets_reciprocal_mesh_at_q(weights_int,
					  grid_points_int,
					  third_q_int,
					  fixed_grid_number,
					  mesh_int,
					  is_time_reversal,
					  lat,
					  num_rot,
					  rot,
					  symprec);

  for (i = 0; i < num_grid; i++) {
    weights_long[i] = (long) weights_int[i];
    third_q_long[i] = (long) third_q_int[i];
    for (j = 0; j < 3; j++) {
      grid_points_long[ i*3 + j ] = (long) grid_points_int[i][j];
    }
  }
  
  return PyInt_FromLong((long) num_ir);
}

static PyObject * extract_triplets_reciprocal_mesh_at_q(PyObject *self, PyObject *args)
{
  int i, j, k;
  PyArrayObject* triplets_at_q;
  PyArrayObject* weight_triplets_at_q;
  int fixed_grid_number;
  PyArrayObject* triplets;
  PyArrayObject* mesh;
  int is_time_reversal;
  PyArrayObject* lattice;
  PyArrayObject* rotations;
  double symprec;
  if (!PyArg_ParseTuple(args, "OOiOOiOOd",
			&triplets_at_q,
			&weight_triplets_at_q,
			&fixed_grid_number,
			&triplets,
			&mesh,
			&is_time_reversal,
			&lattice,
			&rotations,
			&symprec))
    return NULL;

  const long *triplets_long = (long*)triplets->data;
  const int num_triplets = triplets->dimensions[0];
  int triplets_int[num_triplets][3];

  for (i = 0; i < num_triplets; i++) {
    for (j = 0; j < 3; j++) {
      triplets_int[i][j] = (long) triplets_long[ i*3 + j ];
    }
  }

  long *triplets_at_q_long = (long*)triplets_at_q->data;
  int triplets_at_q_int[num_triplets][3];

  long *weight_triplets_at_q_long = (long*)weight_triplets_at_q->data;
  int weight_triplets_at_q_int[num_triplets];

  int mesh_int[3];
  const long* mesh_long = (long*)mesh->data;
  for (i = 0; i < 3; i++) {
    mesh_int[i] = (int) mesh_long[i];
  }  

  SPGCONST double (*lat)[3] = (double(*)[3])lattice->data;
  const long* rot_long = (long*)rotations->data;
  const int num_rot = rotations->dimensions[0];
  int rot[num_rot][3][3];
  for (i = 0; i < num_rot; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
	rot[i][j][k] = (int) rot_long[ i*9 + j*3 + k ];
      }
    }
  }

  const int num_triplets_at_q =
    spg_extract_triplets_reciprocal_mesh_at_q(triplets_at_q_int,
					      weight_triplets_at_q_int,
					      fixed_grid_number,
					      num_triplets,
					      triplets_int,
					      mesh_int,
					      is_time_reversal,
					      lat,
					      num_rot,
					      rot,
					      symprec);

  for (i = 0; i < num_triplets_at_q; i++) {
    weight_triplets_at_q_long[i] = (long) weight_triplets_at_q_int[i];
    for (j = 0; j < 3; j++) {
      triplets_at_q_long[ i*3 + j ] = (long) triplets_at_q_int[i][j];
    }
  }
  
  return PyInt_FromLong((long) num_triplets_at_q);
}


