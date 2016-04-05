/* Copyright (C) 2015 Atsushi Togo */
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

#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

#if (PY_MAJOR_VERSION < 3) && (PY_MINOR_VERSION < 6)
#define PYUNICODE_FROMSTRING PyString_FromString
#else
#define PYUNICODE_FROMSTRING PyUnicode_FromString
#endif

static PyObject * get_version(PyObject *self, PyObject *args);
static PyObject * get_dataset(PyObject *self, PyObject *args);
static PyObject * get_spacegroup_type(PyObject *self, PyObject *args);
static PyObject * get_pointgroup(PyObject *self, PyObject *args);
static PyObject * standardize_cell(PyObject *self, PyObject *args);
static PyObject * refine_cell(PyObject *self, PyObject *args);
static PyObject * get_symmetry(PyObject *self, PyObject *args);
static PyObject *
get_symmetry_with_collinear_spin(PyObject *self, PyObject *args);
static PyObject * find_primitive(PyObject *self, PyObject *args);
static PyObject * get_grid_point_from_address(PyObject *self, PyObject *args);
static PyObject * get_ir_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_stabilized_reciprocal_mesh(PyObject *self, PyObject *args);
static PyObject * get_grid_points_by_rotations(PyObject *self, PyObject *args);
static PyObject * get_BZ_grid_points_by_rotations(PyObject *self, PyObject *args);
static PyObject * relocate_BZ_grid_address(PyObject *self, PyObject *args);
static PyObject * get_symmetry_from_database(PyObject *self, PyObject *args);

struct module_state {
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject *
error_out(PyObject *m) {
  struct module_state *st = GETSTATE(m);
  PyErr_SetString(st->error, "something bad happened");
  return NULL;
}

static PyMethodDef _spglib_methods[] = {
  {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
  {"version", get_version, METH_VARARGS, "Spglib version"},
  {"dataset", get_dataset, METH_VARARGS, "Dataset for crystal symmetry"},
  {"spacegroup_type", get_spacegroup_type, METH_VARARGS, "Space-group type symbols"},
  {"symmetry_from_database", get_symmetry_from_database, METH_VARARGS,
   "Get symmetry operations from database"},
  {"pointgroup", get_pointgroup, METH_VARARGS,
   "International symbol of pointgroup"},
  {"standardize_cell", standardize_cell, METH_VARARGS, "Standardize cell"},
  {"refine_cell", refine_cell, METH_VARARGS, "Refine cell"},
  {"symmetry", get_symmetry, METH_VARARGS, "Symmetry operations"},
  {"symmetry_with_collinear_spin", get_symmetry_with_collinear_spin,
   METH_VARARGS, "Symmetry operations with collinear spin magnetic moments"},
  {"primitive", find_primitive, METH_VARARGS,
   "Find primitive cell in the input cell"},
  {"grid_point_from_address", get_grid_point_from_address, METH_VARARGS,
   "Translate grid adress to grid point index"},
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
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int _spglib_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int _spglib_clear(PyObject *m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_spglib",
  NULL,
  sizeof(struct module_state),
  _spglib_methods,
  NULL,
  _spglib_traverse,
  _spglib_clear,
  NULL
};

#define INITERROR return NULL

PyObject *
PyInit__spglib(void)

#else
#define INITERROR return

  void
  init_spglib(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
  PyObject *module = PyModule_Create(&moduledef);
#else
  PyObject *module = Py_InitModule("_spglib", _spglib_methods);
#endif

  if (module == NULL)
    INITERROR;
  struct module_state *st = GETSTATE(module);

  st->error = PyErr_NewException("_spglib.Error", NULL, NULL);
  if (st->error == NULL) {
    Py_DECREF(module);
    INITERROR;
  }

#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}

static PyObject * get_version(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  PyObject *array;
  int i;
  int version[3];

  version[0] = spg_get_major_version();
  version[1] = spg_get_minor_version();
  version[2] = spg_get_micro_version();

  array = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(array, i, PyLong_FromLong((long)version[i]));
  }

  return array;
}

static PyObject * get_dataset(PyObject *self, PyObject *args)
{
  int i, j, k, n;
  double symprec, angle_tolerance;
  SpglibDataset *dataset;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  PyObject *array, *vec, *mat, *rot, *trans, *wyckoffs, *equiv_atoms;
  PyObject *std_lattice, *std_types, *std_positions;

  if (!PyArg_ParseTuple(args, "OOOdd",
			&lattice,
			&position,
			&atom_type,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  const int num_atom = PyArray_DIMS(position)[0];
  const int* typat = (int*)PyArray_DATA(atom_type);

  dataset = spgat_get_dataset(lat,
			      pos,
			      typat,
			      num_atom,
			      symprec,
			      angle_tolerance);

  array = PyList_New(15);
  n = 0;

  /* Space group number, international symbol, hall symbol */
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->spacegroup_number));
  n++;
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->hall_number));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->international_symbol));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->hall_symbol));
  n++;

  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->transformation_matrix[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }
  PyList_SetItem(array, n, mat);
  n++;

  /* Origin shift */
  vec = PyList_New(3);
  for (i = 0; i < 3; i++) {
    PyList_SetItem(vec, i, PyFloat_FromDouble(dataset->origin_shift[i]));
  }
  PyList_SetItem(array, n, vec);
  n++;

  /* Rotation matrices */
  rot = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    mat = PyList_New(3);
    for (j = 0; j < 3; j++) {
      vec = PyList_New(3);
      for (k = 0; k < 3; k++) {
	PyList_SetItem(vec, k, PyLong_FromLong((long) dataset->rotations[i][j][k]));
      }
      PyList_SetItem(mat, j, vec);
    }
    PyList_SetItem(rot, i, mat);
  }
  PyList_SetItem(array, n, rot);
  n++;

  /* Translation vectors */
  trans = PyList_New(dataset->n_operations);
  for (i = 0; i < dataset->n_operations; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->translations[i][j]));
    }
    PyList_SetItem(trans, i, vec);
  }
  PyList_SetItem(array, n, trans);
  n++;

  /* Wyckoff letters, Equivalent atoms */
  wyckoffs = PyList_New(dataset->n_atoms);
  equiv_atoms = PyList_New(dataset->n_atoms);
  for (i = 0; i < dataset->n_atoms; i++) {
    PyList_SetItem(wyckoffs, i, PyLong_FromLong((long) dataset->wyckoffs[i]));
    PyList_SetItem(equiv_atoms, i, PyLong_FromLong((long) dataset->equivalent_atoms[i]));
  }
  PyList_SetItem(array, n, wyckoffs);
  n++;
  PyList_SetItem(array, n, equiv_atoms);
  n++;

  std_lattice = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->std_lattice[i][j]));
    }
    PyList_SetItem(std_lattice, i, vec);
  }
  PyList_SetItem(array, n, std_lattice);
  n++;

  /* Standardized unit cell */
  std_types = PyList_New(dataset->n_std_atoms);
  std_positions = PyList_New(dataset->n_std_atoms);
  for (i = 0; i < dataset->n_std_atoms; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyFloat_FromDouble(dataset->std_positions[i][j]));
    }
    PyList_SetItem(std_types, i, PyLong_FromLong((long) dataset->std_types[i]));
    PyList_SetItem(std_positions, i, vec);
  }
  PyList_SetItem(array, n, std_types);
  n++;
  PyList_SetItem(array, n, std_positions);
  n++;

  /* Point group */
  PyList_SetItem(array, n, PyLong_FromLong((long) dataset->pointgroup_number));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(dataset->pointgroup_symbol));
  n++;

  spg_free_dataset(dataset);

  return array;
}

static PyObject * get_symmetry_from_database(PyObject *self, PyObject *args)
{
  int hall_number;
  PyArrayObject* rotation;
  PyArrayObject* translation;
  if (!PyArg_ParseTuple(args, "OOi",
			&rotation,
			&translation,
			&hall_number)) {
    return NULL;
  }

  if (PyArray_DIMS(rotation)[0] < 192 || PyArray_DIMS(translation)[0] < 192) {
    Py_RETURN_NONE;
  }

  int (*rot)[3][3] = (int(*)[3][3])PyArray_DATA(rotation);
  double (*trans)[3] = (double(*)[3])PyArray_DATA(translation);


  const int num_sym = spg_get_symmetry_from_database(rot, trans, hall_number);

  return PyLong_FromLong((long) num_sym);
}

static PyObject * get_spacegroup_type(PyObject *self, PyObject *args)
{
  int n, hall_number;
  PyObject *array;
  SpglibSpacegroupType symbols;

  if (!PyArg_ParseTuple(args, "i",&hall_number)) {
    return NULL;
  }

  symbols = spg_get_spacegroup_type(hall_number);

  array = PyList_New(5);
  n = 0;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(symbols.schoenflies));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(symbols.hall_symbol));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(symbols.international));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(symbols.international_full));
  n++;
  PyList_SetItem(array, n, PYUNICODE_FROMSTRING(symbols.international_short));
  n++;

  return array;
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

  SPGCONST int(*rot)[3][3] = (int(*)[3][3])PyArray_DATA(rotations);
  const int num_rot = PyArray_DIMS(rotations)[0];
  const int ptg_num = spg_get_pointgroup(symbol, trans_mat, rot, num_rot);

  /* Transformation matrix */
  mat = PyList_New(3);
  for (i = 0; i < 3; i++) {
    vec = PyList_New(3);
    for (j = 0; j < 3; j++) {
      PyList_SetItem(vec, j, PyLong_FromLong((long)trans_mat[i][j]));
    }
    PyList_SetItem(mat, i, vec);
  }

  array = PyList_New(3);
  PyList_SetItem(array, 0, PYUNICODE_FROMSTRING(symbol));
  PyList_SetItem(array, 1, PyLong_FromLong((long) ptg_num));
  PyList_SetItem(array, 2, mat);

  return array;
}

static PyObject * standardize_cell(PyObject *self, PyObject *args)
{
  int num_atom, to_primitive, no_idealize;
  double symprec, angle_tolerance;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOiiidd",
			&lattice,
			&position,
			&atom_type,
			&num_atom,
			&to_primitive,
			&no_idealize,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  int* typat = (int*)PyArray_DATA(atom_type);

  int num_atom_std = spgat_standardize_cell(lat,
					    pos,
					    typat,
					    num_atom,
					    to_primitive,
					    no_idealize,
					    symprec,
					    angle_tolerance);

  return PyLong_FromLong((long) num_atom_std);
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

  double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  int* typat = (int*)PyArray_DATA(atom_type);

  int num_atom_std = spgat_refine_cell(lat,
				       pos,
				       typat,
				       num_atom,
				       symprec,
				       angle_tolerance);

  return PyLong_FromLong((long) num_atom_std);
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

  double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  int num_atom = PyArray_DIMS(position)[0];
  int* types = (int*)PyArray_DATA(atom_type);

  int num_atom_prim = spgat_find_primitive(lat,
					   pos,
					   types,
					   num_atom,
					   symprec,
					   angle_tolerance);

  return PyLong_FromLong((long) num_atom_prim);
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

  SPGCONST double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  const int* types = (int*)PyArray_DATA(atom_type);
  const int num_atom = PyArray_DIMS(position)[0];
  int (*rot)[3][3] = (int(*)[3][3])PyArray_DATA(rotation);
  double (*trans)[3] = (double(*)[3])PyArray_DATA(translation);
  const int num_sym_from_array_size = PyArray_DIMS(rotation)[0];

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
  return PyLong_FromLong((long) num_sym);
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
  PyArrayObject* equiv_atoms_py;

  if (!PyArg_ParseTuple(args, "OOOOOOOdd",
			&rotation,
			&translation,
			&equiv_atoms_py,
			&lattice,
			&position,
			&atom_type,
			&magmom,
			&symprec,
			&angle_tolerance)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  const double *spins = (double*)PyArray_DATA(magmom);
  const int *types = (int*)PyArray_DATA(atom_type);
  const int num_atom = PyArray_DIMS(position)[0];
  int (*rot)[3][3] = (int(*)[3][3])PyArray_DATA(rotation);
  double (*trans)[3] = (double(*)[3])PyArray_DATA(translation);
  int *equiv_atoms = (int*)PyArray_DATA(equiv_atoms_py);
  const int num_sym_from_array_size = PyArray_DIMS(rotation)[0];

  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_sym =
    spgat_get_symmetry_with_collinear_spin(rot,
					   trans,
					   equiv_atoms,
					   num_sym_from_array_size,
					   lat,
					   pos,
					   types,
					   spins,
					   num_atom,
					   symprec,
					   angle_tolerance);
  return PyLong_FromLong((long) num_sym);
}

static PyObject * get_grid_point_from_address(PyObject *self, PyObject *args)
{
  PyArrayObject* grid_address_py;
  PyArrayObject* mesh_py;
  if (!PyArg_ParseTuple(args, "OO",
			&grid_address_py,
			&mesh_py)) {
    return NULL;
  }

  const int* grid_address = (int*)PyArray_DATA(grid_address_py);
  const int* mesh = (int*)PyArray_DATA(mesh_py);

  const int gp = spg_get_grid_point_from_address(grid_address, mesh);

  return PyLong_FromLong((long) gp);
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
			&symprec)) {
    return NULL;
  }

  SPGCONST double (*lat)[3] = (double(*)[3])PyArray_DATA(lattice);
  SPGCONST double (*pos)[3] = (double(*)[3])PyArray_DATA(position);
  const int* types = (int*)PyArray_DATA(atom_type);
  const int* mesh_int = (int*)PyArray_DATA(mesh);
  const int* is_shift_int = (int*)PyArray_DATA(is_shift);
  const int num_atom = PyArray_DIMS(position)[0];
  int (*grid_address)[3] = (int(*)[3])PyArray_DATA(grid_address_py);
  int *map_int = (int*)PyArray_DATA(map);

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

  return PyLong_FromLong((long) num_ir);
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

  int (*grid_address)[3] = (int(*)[3])PyArray_DATA(grid_address_py);
  int *map_int = (int*)PyArray_DATA(map);
  const int* mesh_int = (int*)PyArray_DATA(mesh);
  const int* is_shift_int = (int*)PyArray_DATA(is_shift);
  SPGCONST int (*rot)[3][3] = (int(*)[3][3])PyArray_DATA(rotations);
  const int num_rot = PyArray_DIMS(rotations)[0];
  SPGCONST double (*q)[3] = (double(*)[3])PyArray_DATA(qpoints);
  const int num_q = PyArray_DIMS(qpoints)[0];

  const int num_ir = spg_get_stabilized_reciprocal_mesh(grid_address,
							map_int,
							mesh_int,
							is_shift_int,
							is_time_reversal,
							num_rot,
							rot,
							num_q,
							q);

  return PyLong_FromLong((long) num_ir);
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

  int *rot_grid_points = (int*)PyArray_DATA(rot_grid_points_py);
  const int *address_orig = (int*)PyArray_DATA(address_orig_py);
  SPGCONST int (*rot_reciprocal)[3][3] = (int(*)[3][3])PyArray_DATA(rot_reciprocal_py);
  const int num_rot = PyArray_DIMS(rot_reciprocal_py)[0];
  const int* mesh = (int*)PyArray_DATA(mesh_py);
  const int* is_shift = (int*)PyArray_DATA(is_shift_py);

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

  int *rot_grid_points = (int*)PyArray_DATA(rot_grid_points_py);
  const int *address_orig = (int*)PyArray_DATA(address_orig_py);
  SPGCONST int (*rot_reciprocal)[3][3] = (int(*)[3][3])PyArray_DATA(rot_reciprocal_py);
  const int num_rot = PyArray_DIMS(rot_reciprocal_py)[0];
  const int* mesh = (int*)PyArray_DATA(mesh_py);
  const int* is_shift = (int*)PyArray_DATA(is_shift_py);
  const int* bz_map = (int*)PyArray_DATA(bz_map_py);

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

  int (*bz_grid_address)[3] = (int(*)[3])PyArray_DATA(bz_grid_address_py);
  int *bz_map = (int*)PyArray_DATA(bz_map_py);
  SPGCONST int (*grid_address)[3] = (int(*)[3])PyArray_DATA(grid_address_py);
  const int* mesh = (int*)PyArray_DATA(mesh_py);
  const int* is_shift = (int*)PyArray_DATA(is_shift_py);
  SPGCONST double (*reciprocal_lattice)[3]  =
    (double(*)[3])PyArray_DATA(reciprocal_lattice_py);
  int num_ir_gp;

  num_ir_gp = spg_relocate_BZ_grid_address(bz_grid_address,
					   bz_map,
					   grid_address,
					   mesh,
					   reciprocal_lattice,
					   is_shift);

  return PyLong_FromLong((long) num_ir_gp);
}
