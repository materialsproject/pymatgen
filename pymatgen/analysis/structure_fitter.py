#!/usr/bin/env python

"""
This module provides classes to perform fitting of two structures.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Geoffroy Hautier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "$Sep 23, 2011M$"

import math
import itertools
import logging
import random

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor

logger = logging.getLogger(__name__)

class StructureFitter(object):
    """
    Class to perform fitting of two structures
    
    Attributes:
        fit_found:
            True if a fit was found
        mapping_op:
            Operation that maps the two structures onto one another.  
            None if no fit was found.
        structure_a:
            First structure
        structure_b:
            Second structure
    """
    FAST_FIT = 25000
    NORMAL_FIT = 50000
    ACCURATE_FIT = 100000
    EXTREME_FIT = 10000000

    def __init__(self, structure_a, structure_b, tolerance_cell_misfit=0.1,
                 tolerance_atomic_misfit=0.3, supercells_allowed=True,
                 anonymized=False, fitting_accuracy=FAST_FIT,
                 use_symmetry=False):
        """
        Fits two structures. All fitting parameters have been set with defaults 
        that should work in most cases. To use, initialize the structure fitter
        with parameters.
        E.g.,
        fitter = StructureFitter(a, b)
        print fitter.fit_found
        
        Args:
            structure_a : 
                First structure
            structure_b : 
                Second structure to try to match with first structure
            tolerance_cell_misfit : 
                Tolerance for cell misfit. Default = 0.1
            tolerance_atomic_misfit : 
                Tolerance for atomic misfit. Default = 0.3.
            supercells_allowed : 
                Whether supercell structures are allowed.  Default = True.
            anonymized : 
                Whether to attempt matching of different species.  Setting this
                to true will allow any two structures with the same framework,
                but different species to match to each other. Default = False.
            fitting_accuracy : 
                An integer setting for the fitting accuracy.  Corresponds to
                the max number of candidate rotations considered.  Use the
                static variables, 
                StructureFitter.FAST_FIT
                StructureFitter.NORMAL_FIT
                StructureFitter.ACCURATE_FIT
                to set the tradeoff between accuracy and speed.  The default, 
                FAST_FIT, should work reasonably well in most instances.
            use_symmetry:
                Whether to use pymatgen.symmetry to determine the spacegroup
                first. Eliminates most non-fits. Defaults to True.
        """

        self._tolerance_cell_misfit = tolerance_cell_misfit
        self._tolerance_atomic_misfit = tolerance_atomic_misfit
        self._supercells_allowed = supercells_allowed
        self._anonymized = anonymized
        self._max_rotations = fitting_accuracy
        #Sort structures first so that they have the same arrangement of species
        self._structure_a = structure_a.get_sorted_structure()
        self._structure_b = structure_b.get_sorted_structure()
        if use_symmetry:
            from pymatgen.symmetry.spglib_adaptor import SymmetryFinder
            finder_a = SymmetryFinder(self._structure_a, symprec=0.1)
            finder_b = SymmetryFinder(self._structure_b, symprec=0.1)
            same_sg = finder_a.get_spacegroup_number() == finder_b.get_spacegroup_number()

        if not use_symmetry or same_sg:
            self._mapping_op = None
            if not self._anonymized:
                self.fit(self._structure_a, self._structure_b)
                if self.fit_found:
                    self.el_mapping = {el:el for el in self._structure_a.composition.elements}
            else:
                comp_a = structure_a.composition
                comp_b = structure_b.composition
                if len(comp_a.elements) == len(comp_b.elements):
                    el_a = comp_a.elements
                    #Create permutations of the specie/elements in structure A 
                    for p in itertools.permutations(el_a):
                        # Create mapping of the specie/elements in structure B to that of A.
                        # Then create a modified structure with those elements and try to fit it.
                        el_mapping = dict(zip(comp_b.elements, p))
                        logger.debug("Using specie mapping " + str(el_mapping))
                        mod = StructureEditor(self._structure_b)
                        mod.replace_species(el_mapping)
                        self.fit(self._structure_a, mod.modified_structure)
                        if self._mapping_op != None:
                            #Store successful element mapping
                            self.el_mapping = {el_a:el_b for el_b, el_a in el_mapping.items()}
                            break
                else:
                    logger.debug("No. of elements in structures are unequal.")
        else:
            self._mapping_op = None
            logger.debug("Symmetry is different.")


    def fit(self, a, b):
        """
        Compares two structures and give the possible affine mapping that
        transforms one into another.
        """
        logger.debug("Structure a")
        logger.debug(str(a))
        logger.debug("Structure b")
        logger.debug(str(b))

        # Check composition first.
        #If compositions are not the same, do not need to fit further.
        if a.composition.reduced_formula != b.composition.reduced_formula or ((a.num_sites != b.num_sites) and not self._supercells_allowed):
            logger.debug('Compositions do not match')
            return None

        logger.debug('Compositions match')

        # Fitting is done by matching sites in one structure (to_fit)
        # to the other (fixed).  
        # We set the structure with more sites as fixed, 
        # and scale the structures to the same density
        (fixed, to_fit) = self._scale_structures(a, b)

        # Defines the atom misfit tolerance
        tol_atoms = self._tolerance_atomic_misfit
        logger.debug("Atomic misfit tolerance = %.4f" % (tol_atoms))

        max_sites = float('inf')
        # determine which type of sites to use for the mapping
        for sp in to_fit.species_and_occu:
            sp_sites = [site for site in to_fit if site.species_and_occu == sp]
            if len(sp_sites) < max_sites:
                fit_sites = sp_sites
                max_sites = len(sp_sites)

        # Set the arbitrary origin
        origin = fit_sites[0]
        logger.debug("Origin = " + str(origin))
        # now that candidate rotations have been found, shift origin of to_fit
        # because the symmetry operations will be applied from this origin

        oshift = SymmOp.from_rotation_matrix_and_translation_vector(np.eye(3), -origin.coords)
        shifted_to_fit = apply_operation(to_fit, oshift)
        found_map = False

        # This is cheating, but let's try a simple rotation first.  In many situations,
        # E.g., when a structure has been topotatically delithiated or substituted, you actually
        # can get a very fast answer without having to try all rotations.
        simple_rot = self._get_rot_matrix(fixed, to_fit)
        if simple_rot is not None:
            rot = SymmOp.from_rotation_matrix_and_translation_vector(simple_rot, np.array([0, 0, 0]))
            (found_map, mapping_op) = self._test_rotation(rot, origin, fixed, shifted_to_fit, tol_atoms)

        if not found_map: #If simple rotation matching does not work, we have to search and try all rotations.
            logger.debug("Identity matching failed. Finding candidate rotations.")
            #Get candidate rotations
            cand_rot = self._get_candidate_rotations(origin, fixed, shifted_to_fit)

            logger.debug(" FOUND {} candidate rotations ".format(len(cand_rot)))

            if len(cand_rot) == 0:
                logger.debug("No candidate rotations found, returning null. ")
                return None

            # sort the operations, the first ones are the ones with small shear
            # this assures us that we find the smallest cell misfit fits
            sorted_cand_rot = sorted(cand_rot.keys(), key=lambda r: cand_rot[r])

            for rot in sorted_cand_rot:
                (found_map, mapping_op) = self._test_rotation(rot, origin, fixed, shifted_to_fit, tol_atoms)

                if found_map:
                    break
        else:
            logger.debug("Identity matching found.")

        logger.debug("Done testing all candidate rotations")

        if mapping_op != None:
            rot = mapping_op.rotation_matrix # maps toFit to fixed
            p = sqrt_matrix(np.dot(rot.transpose(), rot))
            scale_matrix = np.eye(3) * self.scale
            newrot = np.dot(scale_matrix, rot)
            # we need to now make sure fitterdata.MappingOp maps b -> a and not
            # the other way around
            mshift = np.dot(rot, oshift.translation_vector)
            finaltranslation = mapping_op.translation_vector + mshift[0]
            composite_op = SymmOp.from_rotation_matrix_and_translation_vector(newrot, finaltranslation)
            self._mapping_op = composite_op if self.fixed_is_a else composite_op.inverse
            #self._mapping_op = mapping_op
            self._cell_misfit = shear_invariant(p)

    def _get_rot_matrix(self, fixed, to_fit):
        a1 = fixed.lattice.angles
        a2 = to_fit.lattice.angles
        for i in xrange(3):
            if abs(a1[i] - a2[i]) > 5:
                return None
        to_fit_unit_matrix = np.array([row / np.linalg.norm(row) for row in to_fit.lattice.matrix])
        fixed_unit_matrix = np.array([row / np.linalg.norm(row) for row in fixed.lattice.matrix])
        return np.dot(to_fit_unit_matrix.transpose(), np.linalg.inv(fixed_unit_matrix.transpose()))

    def _test_rotation(self, rot, origin, fixed, shifted_to_fit, tol_atoms):
        logger.debug("Trying candidate rotation : \n" + str(rot))
        for site in fixed:
            if site.species_and_occu == origin.species_and_occu:
                shift = site.coords
                op = SymmOp.from_rotation_matrix_and_translation_vector(rot.rotation_matrix, shift)
                logger.debug("Trying op : \n" + str(op))
                nstruct = apply_operation(shifted_to_fit, op)
                correspondance = compare_structures(fixed, nstruct, tol_atoms)
                if correspondance:
                    self.correspondance = correspondance
                    return (True, op)
        return (False, None)

    def __str__(self):
        output = ["Fitting structures"]
        output.append("\nStructure 1:")
        output.append(str(self._structure_a))
        output.append("\nStructure 2:")
        output.append(str(self._structure_b))
        output.append("\nFitting parameters:")
        output.append("\tTolerance cell misfit = " + str(self._tolerance_cell_misfit))
        output.append("\tSupercells allowed = " + str(self._supercells_allowed))
        output.append("\tAnonymized  = " + str(self._anonymized))
        output.append("\nFitting " + ("succeeded " if self._mapping_op != None else "failed"))
        if self._mapping_op != None:
            output.append("\tMapping op = " + str(self._mapping_op))
            output.append("\tCell misfit = " + str(self._cell_misfit))
            if self._anonymized:
                output.append("Element mapping")
                output.append("\n".join([str(k) + " -> " + str(v) for k, v in self.el_mapping.items()]))
        return "\n".join(output)

    def _scale_structures(self, a, b):
        # which structure do we want to fit to the other ?
        # assume that structure b has less sites and switch if needed
        self.fixed_is_a = a.num_sites > b.num_sites
        (fixed, to_fit) = (a, b) if self.fixed_is_a else (b, a)

        # scale the structures to the same density
        rho_a = fixed.num_sites / fixed.volume
        rho_b = to_fit.num_sites / to_fit.volume
        scale = (rho_b / rho_a) ** (1 / 3)
        self.scale = scale
        to_fit = Structure(Lattice(to_fit.lattice.matrix * scale), to_fit.species_and_occu, to_fit.frac_coords)
        return (fixed, to_fit)

    def _get_candidate_rotations(self, origin, fixed, to_fit):
        tol_shear = self._tolerance_cell_misfit
        fixed_basis = fixed.lattice.matrix.transpose()
        # need to generate candidate rotations ...
        lengths = fixed.lattice.abc

        shells = []
        for i in range(3):
            dr = lengths[i] * math.sqrt(tol_shear / 2)
            shell = to_fit.get_neighbors_in_shell(np.array([0, 0, 0]), lengths[i], dr)
            logger.debug("shell {} radius={} dr={}".format(i, lengths[i], dr))
            shell = filter(lambda x: x[0].species_and_occu == origin.species_and_occu, shell)
            shell = sorted(shell, key=lambda x: x[1])
            shells.append([site for (site, dist) in shell])
            logger.debug("No. in shell = {}".format(len(shells[-1])))
        # now generate candidate rotations
        cand_rot = {} #Dict of SymmOp : float
        a = len(shells[0])
        b = len(shells[1])
        c = len(shells[2])
        total_rots = a * b * c
        if total_rots < self._max_rotations:
            logger.debug("Total rots = {}. Using all rotations.".format(total_rots))
            test_rotations = itertools.product(*shells)
        else:
            logger.info("Total rots = {m} exceed max_rotations = {n}. Using {n} randomly selected rotations.".format(m=total_rots, n=self._max_rotations))
            def random_rot():
                considered_rots = []
                while len(considered_rots) < self._max_rotations:
                    (x, y, z) = [random.randint(0, i - 1) for i in [a, b, c]]
                    if (x, y, z) not in considered_rots:
                        considered_rots.append((x, y, z))
                        yield (shells[0][x], shells[1][y], shells[2][z])
            test_rotations = random_rot()

        for pool in test_rotations:
            # now, can a unitary transformation bring the cell vectors together
            cell_v = np.array([nn.coords for nn in pool]).transpose()
            det = np.linalg.det(cell_v)
            if abs(det) < 0.001 or abs(abs(det) - fixed.volume) > 0.01:
                continue
            rot = np.dot(fixed_basis, np.linalg.inv(cell_v))
            r = SymmOp.from_rotation_matrix_and_translation_vector(rot, np.array([0, 0, 0]))
            if r not in cand_rot:
                transf = r.rotation_matrix
                transf = np.dot(transf.transpose(), transf)
                transf = np.eye(3) if almost_identity(transf) else transf
                pbis = sqrt_matrix(transf)
                shear_inv = shear_invariant(pbis)
                if shear_inv < tol_shear:
                    cand_rot[r] = shear_inv
                else:
                    logging.debug("Shear {} exceeds tol of {}".format(shear_inv, tol_shear))
        return cand_rot

    @property
    def fit_found(self):
        """
        True if a fit was found.
        """
        return self._mapping_op != None

    @property
    def mapping_op(self):
        """
        Operation mapping the two structures.  None if no fit was found.
        """
        return self._mapping_op

    @property
    def structure_a(self):
        """
        First input structure
        """
        return self._structure_a

    @property
    def structure_b(self):
        """
        Second input structure
        """
        return self._structure_b

def apply_operation(structure, symmop):
    editor = StructureEditor(structure)
    editor.apply_operation(symmop)
    return editor.modified_structure

def sqrt_matrix(input_matrix):
    d, v = np.linalg.eig(input_matrix)
    diagonalbis = np.array([[d[0] ** 0.5, 0, 0], [0, d[1] ** 0.5, 0], [0, 0, d[2] ** 0.5]])
    temp = np.dot(diagonalbis, v.transpose())
    result = np.dot(v, temp)
    return result

def shear_invariant(matrix):
    return (matrix[0][0] - matrix[1][1]) ** 2 + (matrix[1][1] - matrix[2][2]) ** 2 + (matrix[0][0] - matrix[2][2]) ** 2 + 6 * (matrix[0][1] * matrix[0][1] + matrix[0][2] * matrix[0][2] + matrix[1][2] * matrix[1][2])

def almost_identity(mat):
    """
    This is to resolve a very very strange bug in numpy that causes random eigen vectors to be returned
    for matrices very very close to the identity matrix.  See test_eig for examples.
    """
    return (abs(mat - np.eye(3)) < 1e-10).all()


def compare_structures(larger_struct, smaller_struct, tol):
    assert (len(larger_struct) >= len(smaller_struct)), "Structures are not in the right order!"
    mapping = {}
    for site1 in smaller_struct:
        mapping[site1] = []
        for site2 in larger_struct:
            if site1.species_and_occu == site2.species_and_occu:
                #Map site in larger struct to smaller lattice
                fcoords = site1.lattice.get_fractional_coords(site2.coords)
                dist = site1.distance_and_image_from_frac_coords(fcoords)[0]
                if dist < tol:
                    mapping[site1].append(site2)
    num_map = [len(v) for v in mapping.values()]
    if all([i == num_map[0] for i in num_map]):
        return mapping
    return None
