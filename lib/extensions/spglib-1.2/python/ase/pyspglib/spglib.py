"""
Spglib interface for ASE
"""

import pyspglib._spglib as spg
import numpy as np

def get_symmetry(bulk, symprec=1e-5, angle_tolerance=-1.0):
    """
    Return symmetry operations as hash.
    Hash key 'rotations' gives the numpy integer array
    of the rotation matrices for scaled positions
    Hash key 'translations' gives the numpy float64 array
    of the translation vectors in scaled positions
    """

    # Atomic positions have to be specified by scaled positions for spglib.
    positions = bulk.get_scaled_positions().copy()
    lattice = bulk.get_cell().T.copy()
    numbers = bulk.get_atomic_numbers()
  
    # Get number of symmetry operations and allocate symmetry operations
    # multi = spg.multiplicity(cell, positions, numbers, symprec)
    multi = 48 * bulk.get_number_of_atoms()
    rotation = np.zeros((multi, 3, 3), dtype=int)
    translation = np.zeros((multi, 3))
  
    # Get symmetry operations
    magmoms = bulk.get_magnetic_moments()
    if magmoms == None:
        num_sym = spg.symmetry(rotation,
                               translation,
                               lattice,
                               positions,
                               numbers,
                               symprec,
                               angle_tolerance)
    else:
        num_sym = spg.symmetry_with_collinear_spin(rotation,
                                                   translation,
                                                   lattice,
                                                   positions,
                                                   numbers,
                                                   magmoms,
                                                   symprec,
                                                   angle_tolerance)
  
    return {'rotations': rotation[:num_sym],
            'translations': translation[:num_sym]}

def get_symmetry_dataset(bulk, symprec=1e-5, angle_tolerance=-1.0):
    """
    number: International space group number
    international: International symbol
    hall: Hall symbol
    transformation_matrix:
      Transformation matrix from lattice of input cell to Bravais lattice
      L^bravais = L^original * Tmat
    origin shift: Origin shift in the setting of 'Bravais lattice'
    rotations, translations:
      Rotation matrices and translation vectors
      Space group operations are obtained by
        [(r,t) for r, t in zip(rotations, translations)]
    wyckoffs:
      Wyckoff letters
    """
    positions = bulk.get_scaled_positions().copy()
    lattice = bulk.get_cell().T.copy()
    numbers = bulk.get_atomic_numbers()
    keys = ('number',
            'international',
            'hall',
            'transformation_matrix',
            'origin_shift',
            'rotations',
            'translations',
            'wyckoffs',
            'equivalent_atoms')
    dataset = {}
    for key, data in zip(keys, spg.dataset(lattice,
                                           positions,
                                           numbers,
                                           symprec,
                                           angle_tolerance)):
        dataset[key] = data

    dataset['international'] = dataset['international'].strip()
    dataset['hall'] = dataset['hall'].strip()
    dataset['transformation_matrix'] = np.array(dataset['transformation_matrix'])
    dataset['origin_shift'] = np.array(dataset['origin_shift'])
    dataset['rotations'] = np.array(dataset['rotations'])
    dataset['translations'] = np.array(dataset['translations'])
    letters = "abcdefghijklmnopqrstuvwxyz"
    dataset['wyckoffs'] = [letters[x] for x in dataset['wyckoffs']]
    dataset['equivalent_atoms'] = np.array(dataset['equivalent_atoms'])

    return dataset

def get_spacegroup(bulk, symprec=1e-5, angle_tolerance=-1.0):
    """
    Return space group in international table symbol and number
    as a string.
    """
    # Atomic positions have to be specified by scaled positions for spglib.
    return spg.spacegroup(bulk.get_cell().T.copy(),
                          bulk.get_scaled_positions().copy(),
                          bulk.get_atomic_numbers(),
                          symprec,
                          angle_tolerance)

def get_pointgroup(rotations):
    """
    Return point group in international table symbol and number.
    The symbols are mapped to the numbers as follows:
    1   "1    "
    2   "-1   "
    3   "2    "
    4   "m    "
    5   "2/m  "
    6   "222  "
    7   "mm2  "
    8   "mmm  "
    9   "4    "
    10  "-4   "
    11  "4/m  "
    12  "422  "
    13  "4mm  "
    14  "-42m "
    15  "4/mmm"
    16  "3    "
    17  "-3   "
    18  "32   "
    19  "3m   "
    20  "-3m  "
    21  "6    "
    22  "-6   "
    23  "6/m  "
    24  "622  "
    25  "6mm  "
    26  "-62m "
    27  "6/mmm"
    28  "23   "
    29  "m-3  "
    30  "432  "
    31  "-43m "
    32  "m-3m "
    """

    # (symbol, pointgroup_number, transformation_matrix)
    return spg.pointgroup(rotations)

def refine_cell(bulk, symprec=1e-5, angle_tolerance=-1.0):
    """
    Return refined cell
    """
    # Atomic positions have to be specified by scaled positions for spglib.
    num_atom = bulk.get_number_of_atoms()
    lattice = bulk.get_cell().T.copy()
    pos = np.zeros((num_atom * 4, 3), dtype=float)
    pos[:num_atom] = bulk.get_scaled_positions()

    numbers = np.zeros(num_atom * 4, dtype=int)
    numbers[:num_atom] = bulk.get_atomic_numbers()
    num_atom_bravais = spg.refine_cell(lattice,
                                       pos,
                                       numbers,
                                       num_atom,
                                       symprec,
                                       angle_tolerance)

    return lattice.T.copy(), pos[:num_atom_bravais], numbers[:num_atom_bravais]

def find_primitive(bulk, symprec=1e-5, angle_tolerance=-1.0):
    """
    A primitive cell in the input cell is searched and returned
    as an object of Atoms class.
    If no primitive cell is found, (None, None, None) is returned.
    """

    # Atomic positions have to be specified by scaled positions for spglib.
    positions = bulk.get_scaled_positions().copy()
    lattice = bulk.get_cell().T.copy()
    numbers = bulk.get_atomic_numbers()

    # lattice is transposed with respect to the definition of Atoms class
    num_atom_prim = spg.primitive(lattice,
                                  positions,
                                  numbers,
                                  symprec,
                                  angle_tolerance)
    if num_atom_prim > 0:
        return lattice.T, positions[:num_atom_prim], numbers[:num_atom_prim]
    else:
        return None, None, None
  
def get_ir_kpoints(kpoint,
                   bulk,
                   is_time_reversal=True,
                   symprec=1e-5):
    """
    Retrun irreducible kpoints
    """
    mapping = np.zeros(kpoint.shape[0], dtype=int)
    spg.ir_kpoints(mapping,
                   kpoint,
                   bulk.get_cell().T.copy(),
                   bulk.get_scaled_positions().copy(),
                   bulk.get_atomic_numbers(),
                   is_time_reversal * 1,
                   symprec)
    return mapping
  
def get_ir_reciprocal_mesh(mesh,
                           bulk,
                           is_shift=np.zeros(3, dtype=int),
                           is_time_reversal=True,
                           symprec=1e-5):
    """
    Return k-points mesh and k-point map to the irreducible k-points
    The symmetry is serched from the input cell.
    is_shift=[0, 0, 0] gives Gamma center mesh.
    """

    mapping = np.zeros(np.prod(mesh), dtype=int)
    mesh_points = np.zeros((np.prod(mesh), 3), dtype=int)
    spg.ir_reciprocal_mesh(mesh_points,
                           mapping,
                           np.array(mesh),
                           np.array(is_shift),
                           is_time_reversal * 1,
                           bulk.get_cell().T.copy(),
                           bulk.get_scaled_positions().copy(),
                           bulk.get_atomic_numbers(), symprec)
  
    return mapping, mesh_points
  
def get_stabilized_reciprocal_mesh(mesh,
                                   lattice,
                                   rotations,
                                   is_shift=np.zeros(3, dtype=int),
                                   is_time_reversal=True,
                                   qpoints=np.array([], dtype=float),
                                   symprec=1e-5):
    """
    Return k-point map to the irreducible k-points and k-point grid points .

    The symmetry is searched from the input rotation matrices in real space.
    The convention of 'lattice' is:
       [[a_x, a_y, a_z],
        [b_x, b_y, b_z],
        [c_x, c_y, c_z]] (ASE convention)
    Since it has to be passed to the C extention in the following format:
       [[a_x, b_x, c_x],
        [a_y, b_y, c_y],
        [a_z, b_z, c_z]] (spglib convention)
    Therefore in this method, lattice is transposed.
    
    is_shift=[0, 0, 0] gives Gamma center mesh and the values 1 give
    half mesh distance shifts.
    """
    
    mapping = np.zeros(np.prod(mesh), dtype=int)
    mesh_points = np.zeros((np.prod(mesh), 3), dtype=int)
    qpoints = np.array(qpoints, dtype=float)
    if qpoints.shape == (3,):
        qpoints = np.array([qpoints])
    spg.stabilized_reciprocal_mesh(mesh_points,
                                   mapping,
                                   np.array(mesh, dtype=int),
                                   np.array(is_shift, dtype=int),
                                   is_time_reversal * 1,
                                   lattice.T.copy(),
                                   rotations.copy(),
                                   np.array(qpoints, dtype=float),
                                   symprec)
    
    return mapping, mesh_points

def get_triplets_reciprocal_mesh(mesh,
                                 lattice,
                                 pointgroup,
                                 is_time_reversal=True,
                                 symprec=1e-5):
    """
    Return symmetry reduced triplets (set of addresses) and
    k-point grid points corresponding to addresses.
    The k-point grid is accessed by mesh_points[address].

    The symmetry is searched from the input rotation matrices in real space.
    is_shift=[0, 0, 0] gives Gamma center mesh and the values 1 give
    half mesh distance shifts.
    """
    
    triplets, weights, mesh_points = \
        spg.triplets_reciprocal_mesh(np.array(mesh, dtype=int),
                                     is_time_reversal * 1,
                                     lattice.T.copy(),
                                     pointgroup.copy(),
                                     symprec)

    return np.array(triplets), np.array(weights), np.array(mesh_points)

def get_triplets_reciprocal_mesh_at_q(fixed_grid_number,
                                      mesh,
                                      lattice,
                                      rotations,
                                      is_time_reversal=True,
                                      symprec=1e-5):

    weights = np.zeros(np.prod(mesh), dtype=int)
    third_q = np.zeros(np.prod(mesh), dtype=int)
    mesh_points = np.zeros((np.prod(mesh), 3), dtype=int)
    

    spg.triplets_reciprocal_mesh_at_q(weights,
                                      mesh_points,
                                      third_q,
                                      fixed_grid_number,
                                      np.array(mesh, dtype=int),
                                      is_time_reversal * 1,
                                      lattice.T.copy(),
                                      rotations.copy(),
                                      symprec)

    return weights, third_q, mesh_points
        

def extract_triplets_reciprocal_mesh_at_q(fixed_grid_number,
                                          triplets,
                                          mesh,
                                          lattice,
                                          pointgroup,
                                          is_time_reversal=True,
                                          symprec=1e-5):

    triplets_with_q = np.zeros((len(triplets), 3), dtype=int)
    weights_with_q = np.zeros(len(triplets), dtype=int)

    num_triplets_with_q = \
        spg.triplets_reciprocal_mesh_at_q_from_triplets(triplets_with_q,
                                                        weights_with_q,
                                                        fixed_grid_number,
                                                        triplets,
                                                        np.array(mesh, dtype=int),
                                                        is_time_reversal * 1,
                                                        lattice.T.copy(),
                                                        pointgroup.copy(),
                                                        symprec)
    
    return \
        triplets_with_q[:num_triplets_with_q], \
        weights_with_q[:num_triplets_with_q]
