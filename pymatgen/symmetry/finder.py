#!/usr/bin/env python

"""
An interface to the excellent spglib library by Atsushi Togo
(http://spglib.sourceforge.net/) for pymatgen.

v1.0 - Now works with both ordered and disordered structure.

.. note::
    Not all spglib functions are implemented.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"

import re
import itertools

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.spacegroup import Spacegroup
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.core.operations import SymmOp

try:
    import pymatgen._spglib as spg
except ImportError:
    try:
        import pyspglib._spglib as spg
    except ImportError:
        msg = "Spglib required. Please either run python setup.py install" + \
            " for pymatgen, or install pyspglib from spglib."
        raise ImportError(msg)


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

        #Spglib"s convention for the lattice definition is the transpose of the
        #pymatgen version.
        self._transposed_latt = structure.lattice.matrix.transpose()
        #Spglib requires numpy floats.
        self._transposed_latt = self._transposed_latt.astype(float)
        self._positions = np.array([site.frac_coords for site in structure])
        unique_species = []
        zs = []

        for species, g in itertools.groupby(structure,
                                            key=lambda s: s.species_and_occu):
            try:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(g)))
            except ValueError:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(g)))
        self._unique_species = unique_species
        self._numbers = np.array(zs)
        self._spacegroup_data = spg.spacegroup(self._transposed_latt.copy(),
                                               self._positions.copy(),
                                               self._numbers, self._symprec,
                                               self._angle_tol)

    def get_spacegroup(self):
        """
        Return space group in international table symbol and number
        as a string.
        """
        # Atomic positions have to be specified by scaled positions for spglib.
        return Spacegroup(self.get_spacegroup_symbol(),
                          self.get_spacegroup_number(),
                          self.get_symmetry_operations())

    def get_spacegroup_symbol(self):
        """
        Get the spacegroup symbol (e.g., Pnma) for structure.

        Returns:
            Spacegroup symbol for structure.
        """
        return self._spacegroup_data.split()[0]

    def get_spacegroup_number(self):
        """
        Get the international spacegroup number (e.g., 62) for structure.

        Returns:
            International spacegroup number for structure.
        """
        sgnum = self._spacegroup_data.split()[1]
        sgnum = int(re.sub("\D", "", sgnum))
        return sgnum

    def get_hall(self):
        """
        Returns Hall symbol for structure.
        """
        ds = self.get_symmetry_dataset()
        return ds["hall"]

    def get_pointgroup(self):
        """
        Get the point group associated with the structure.

        Returns:
            Pointgroup for structure.
        """
        ds = self.get_symmetry_dataset()
        return get_pointgroup(ds["rotations"])[0].strip()

    def get_crystal_system(self):
        """
        Get the crystal system for the structure, e.g., (triclinic,
        orthorhombic, cubic, etc.).

        Returns:
            Crystal system for structure.
        """
        n = self.get_spacegroup_number()

        f = lambda i, j: i <= n <= j
        cs = {}
        cs["triclinic"] = (1, 2)
        cs["monoclinic"] = (3, 15)
        cs["orthorhombic"] = (16, 74)
        cs["tetragonal"] = (75, 142)
        cs["trigonal"] = (143, 167)
        cs["hexagonal"] = (168, 194)
        cs["cubic"] = (195, 230)

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
                Transformation matrix from lattice of input cell to Bravais
                lattice
                L^bravais = L^original * Tmat
            origin shift:
                Origin shift in the setting of "Bravais lattice"
            rotations, translations:
                Rotation matrices and translation vectors.
                Space group operations are obtained by
                [(r,t) for r, t in zip(rotations, translations)]
            wyckoffs:
                  Wyckoff letters
        """
        keys = ("number",
                "international",
                "hall",
                "transformation_matrix",
                "origin_shift",
                "rotations",
                "translations",
                "wyckoffs",
                "equivalent_atoms")
        dataset = {}
        for key, data in zip(keys, spg.dataset(self._transposed_latt.copy(),
                                               self._positions, self._numbers,
                                               self._symprec,
                                               self._angle_tol)):
            dataset[key] = data

        dataset["international"] = dataset["international"].strip()
        dataset["hall"] = dataset["hall"].strip()
        dataset["transformation_matrix"] = \
            np.array(dataset["transformation_matrix"])
        dataset["origin_shift"] = np.array(dataset["origin_shift"])
        dataset["rotations"] = np.array(dataset["rotations"])
        dataset["translations"] = np.array(dataset["translations"])
        letters = "abcdefghijklmnopqrstuvwxyz"
        dataset["wyckoffs"] = [letters[x] for x in dataset["wyckoffs"]]
        dataset["equivalent_atoms"] = np.array(dataset["equivalent_atoms"])

        return dataset

    def _get_symmetry(self):
        """
        Get the symmetry operations associated with the structure.

        Returns:
            Symmetry operations as a tuple of two equal length sequences.
            (rotations, translations). "rotations" is the numpy integer array
            of the rotation matrices for scaled positions
            "translations" gives the numpy float64 array of the translation
            vectors in scaled positions.
        """

        # Get number of symmetry operations and allocate symmetry operations
        # multi = spg.multiplicity(cell, positions, numbers, symprec)
        multi = 48 * self._structure.num_sites
        rotation = np.zeros((multi, 3, 3), dtype=int)
        translation = np.zeros((multi, 3))

        num_sym = spg.symmetry(rotation, translation,
                               self._transposed_latt.copy(),
                               self._positions, self._numbers, self._symprec,
                               self._angle_tol)
        return (rotation[:num_sym], translation[:num_sym])

    def get_symmetry_operations(self, cartesian=False):
        """
        Return symmetry operations as a list of SymmOp objects.
        By default returns fractional coord symmops.
        But cartesian can be returned too.
        """
        (rotation, translation) = self._get_symmetry()
        symmops = []
        for rot, trans in zip(rotation, translation):
            if cartesian:
                rot = np.dot(self._structure.lattice.md2c, np.dot(rot,
                                                self._structure.lattice.mc2d))
                trans = np.dot(self._structure.lattice.md2c, trans)
            op = SymmOp.from_rotation_and_translation(rot, trans)
            symmops.append(op)
        return symmops

    def get_symmetrized_structure(self):
        """
        Get a symmetrized structure. A symmetrized structure is one where the
        sites have been grouped into symmetrically equivalent groups.

        Returns:
            pymatgen.symmetry.structure.SymmetrizedStructure object.
        """
        ds = self.get_symmetry_dataset()
        sg = Spacegroup(self.get_spacegroup_symbol(),
                        self.get_spacegroup_number(),
                        self.get_symmetry_operations())
        return SymmetrizedStructure(self._structure, sg,
                                    ds["equivalent_atoms"])

    def get_refined_structure(self):
        """
        Get the refined structure based on detected symmetry.

        Returns:
            Refined structure.
        """
        # Atomic positions have to be specified by scaled positions for spglib.
        num_atom = self._structure.num_sites
        lattice = self._transposed_latt.copy()
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
        Find a primitive version of the unit cell.

        Returns:
            A primitive cell in the input cell is searched and returned
            as an Structure object. If no primitive cell is found, None is
            returned.
        """

        # Atomic positions have to be specified by scaled positions for spglib.
        positions = self._positions.copy()
        lattice = self._transposed_latt.copy()
        numbers = self._numbers.copy()
        # lattice is transposed with respect to the definition of Atoms class
        num_atom_prim = spg.primitive(lattice, positions, numbers,
                                      self._symprec, self._angle_tol)
        zs = numbers[:num_atom_prim]
        species = [self._unique_species[i - 1] for i in zs]

        if num_atom_prim > 0:
            return Structure(lattice.T, species, positions[:num_atom_prim])
        else:
            #Not sure if we should return None or just return the full
            #structure.
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
    if rotations.dtype == "float64":
        rotations = np.int_(rotations)
    # (symbol, pointgroup_number, transformation_matrix)
    return spg.pointgroup(rotations)
