#!/usr/bin/env python

'''
An interface to the excellent spglib library by Atsushi Togo 
(http://spglib.sourceforge.net/) for pymatgen.

v1.0 - Now works with both ordered and disordered structure.

.. note::
    This is a *beta* version. Not all spglib functions are implemented. 
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"

import re

import numpy as np
import itertools

from pymatgen.core.structure import Structure
from pymatgen.symmetry.spacegroup import Spacegroup
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.core.operations import SymmOp

import pyspglib._spglib as spg


class SymmetryFinder(object):
    """
    Takes a pymatgen.core.structure.Structure object and a symprec.
    Uses pyspglib to perform various symmetry finding operations.
    """

    def __init__(self, structure, symprec=1e-5, angle_tolerance=5):
        """
        Args:
            structure:
                Structure object
            symprec:
                Tolerance for symmetry finding
            angle_tolerance:
                Angle tolerance for symmetry finding.
        """
        self._symprec = symprec
        self._angle_tol = angle_tolerance
        self._structure = structure
        self._lattice = structure.lattice.matrix
        self._positions = np.array([site.frac_coords for site in structure])
        unique_species = []
        zs = []

        for species, g in itertools.groupby(structure, key=lambda site: site.species_and_occu):
            try:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(g)))
            except ValueError:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(g)))

        self._unique_species = unique_species
        self._numbers = np.array(zs)

        self._spacegroup_data = spg.spacegroup(self._lattice.transpose().copy(), self._positions.copy(), self._numbers, self._symprec, self._angle_tol)

    def get_spacegroup(self):
        """
        Return space group in international table symbol and number
        as a string.
        """
        # Atomic positions have to be specified by scaled positions for spglib.    
        return Spacegroup(self.get_spacegroup_symbol(), self.get_spacegroup_number(), self.get_symmetry_operations())

    def get_spacegroup_symbol(self):
        return re.split("\s+", self._spacegroup_data)[0]

    def get_spacegroup_number(self):
        sgnum = re.split("\s+", self._spacegroup_data)[1]
        sgnum = int(re.sub("\D", "", sgnum))
        return sgnum

    def get_hall(self):
        ds = self.get_symmetry_dataset()
        return ds['hall']

    def get_pointgroup(self):
        ds = self.get_symmetry_dataset()
        return get_pointgroup(ds['rotations'])[0].strip()

    def get_crystal_system(self):
        n = self.get_spacegroup_number()

        f = lambda i, j: i <= n <= j
        cs = {}
        cs['triclinic'] = (1, 2)
        cs['monoclinic'] = (3, 15)
        cs['orthorhombic'] = (16, 74)
        cs['tetragonal'] = (75, 142)
        cs['trigonal'] = (143, 167)
        cs['hexagonal'] = (168, 194)
        cs['cubic'] = (195, 230)

        crystal_sytem = None

        for k, v in cs.items():
            if f(*v) == True:
                crystal_sytem = k
                break
        return crystal_sytem

    def get_symmetry_dataset(self):
        """
        Returns the symmetry dataset as a dict.
        
        Returns:
            number: 
                International space group number
            international: 
                International symbol
            hall: 
                Hall symbol
            transformation_matrix: 
                Transformation matrix from lattice of input cell to Bravais lattice
                L^bravais = L^original * Tmat
            origin shift: 
                Origin shift in the setting of 'Bravais lattice'
            rotations, translations:
                Rotation matrices and translation vectors. 
                Space group operations are obtained by
                [(r,t) for r, t in zip(rotations, translations)]
            wyckoffs:
                  Wyckoff letters
        """
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
        for key, data in zip(keys, spg.dataset(self._lattice.transpose().copy(), self._positions, self._numbers, self._symprec, self._angle_tol)):
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

    def get_symmetry(self):
        """
        Return symmetry operations as hash.
        Hash key 'rotations' gives the numpy integer array
        of the rotation matrices for scaled positions
        Hash key 'translations' gives the numpy float64 array
        of the translation vectors in scaled positions
        """

        # Get number of symmetry operations and allocate symmetry operations
        # multi = spg.multiplicity(cell, positions, numbers, symprec)
        multi = 48 * self._structure.num_sites
        rotation = np.zeros((multi, 3, 3), dtype=int)
        translation = np.zeros((multi, 3))

        num_sym = spg.symmetry(rotation, translation, self._lattice.transpose().copy(),
                                   self._positions, self._numbers, self._symprec, self._angle_tol)
        return (rotation[:num_sym], translation[:num_sym])

    def get_symmetry_operations(self, cartesian=False):
        """
        Return symmetry operations as a list of SymmOp objects.
        By default returns fractional coord symmops.
        But cartesian can be returned too.
        """
        (rotation, translation) = self.get_symmetry()
        symmops = []
        for rot, trans in zip(rotation, translation):
            if cartesian:
                rot = np.dot(self._structure.lattice.md2c, np.dot(rot, self._structure.lattice.mc2d))
                trans = np.dot(self._structure.lattice.md2c, trans)
            symmops.append(SymmOp.from_rotation_matrix_and_translation_vector(rot, trans))
        return symmops

    def get_symmetrized_structure(self):
        ds = self.get_symmetry_dataset()
        sg = Spacegroup(self.get_spacegroup_symbol(), self.get_spacegroup_number(), self.get_symmetry_operations())
        return SymmetrizedStructure(self.get_refined_structure(), sg, ds['equivalent_atoms'])

    def get_refined_structure(self):
        """
        Return refined Structure
        """
        # Atomic positions have to be specified by scaled positions for spglib.
        num_atom = self._structure.num_sites
        lattice = self._lattice.T.copy()
        pos = np.zeros((num_atom * 4, 3), dtype=float)
        pos[:num_atom] = self._positions.copy()

        numbers = np.zeros(num_atom * 4, dtype=int)
        numbers[:num_atom] = self._numbers.copy()
        num_atom_bravais = spg.refine_cell(lattice,
                                           pos,
                                           numbers,
                                           num_atom,
                                           self._symprec,
                                           self._angle_tol)
        zs = numbers[:num_atom_bravais]
        species = [self._unique_species[i - 1] for i in zs]
        return Structure(lattice.T.copy(), species, pos[:num_atom_bravais])


    def find_primitive(self):
        """
        A primitive cell in the input cell is searched and returned
        as an Structure object.
        If no primitive cell is found, (None, None, None) is returned.
        """

        # Atomic positions have to be specified by scaled positions for spglib.
        positions = self._positions.copy()
        lattice = self._lattice.T.copy()
        numbers = self._numbers.copy()
        # lattice is transposed with respect to the definition of Atoms class
        num_atom_prim = spg.primitive(lattice, positions, numbers, self._symprec, self._angle_tol)
        zs = numbers[:num_atom_prim]
        species = [self._unique_species[i - 1] for i in zs]

        if num_atom_prim > 0:
            return Structure(lattice.T, species, positions[:num_atom_prim])
        else:
            #Not sure if we should return None or just return the full structure.
            return None

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
    if rotations.dtype == 'float64':
        rotations = np.int_(rotations)
    # (symbol, pointgroup_number, transformation_matrix)
    return spg.pointgroup(rotations)



'''
Shyue
-----
TODO : The following methods have not yet been ported over.
 
def get_ir_kpoints(kpoint, bulk, is_time_reversal=True, symprec=1e-5):
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
'''
