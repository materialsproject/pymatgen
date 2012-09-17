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
__date__ = "Jun 29, 2012"

import math
import itertools
import logging
import random
import time
from collections import OrderedDict

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.symmetry.finder import SymmetryFinder

logger = logging.getLogger(__name__)


class StructureFitter(object):
    """
    Class to perform fitting of two structures.

    Attributes:
        fit_found:
            True if a fit was found
        mapping_op:
            Operation that maps the two structures onto one another. None if no
            fit was found.
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
                 tolerance_atomic_misfit=1.0, supercells_allowed=True,
                 anonymized=False, fitting_accuracy=FAST_FIT,
                 timeout=300, symmetry_tol=0):
        """
        Fits two structures. All fitting parameters have been set with defaults
        that should work in most cases.

        Args:
            structure_a:
                First structure
            structure_b:
                Second structure to try to match with first structure
            tolerance_cell_misfit:
                Tolerance for cell misfit. Default = 0.1
            tolerance_atomic_misfit:
                Tolerance for atomic misfit. Default = 1.0.
            supercells_allowed:
                Whether supercell structures are allowed.  Default = True.
            anonymized:
                Whether to attempt matching of different species.  Setting this
                to true will allow any two structures with the same framework,
                but different species to match to each other. Default = False.
            fitting_accuracy:
                An integer setting for the fitting accuracy.  Corresponds to
                the max number of candidate rotations considered.  Use the
                static variables FAST_FIT, NORMAL_FIT, ACCURATE_FIT and
                EXTREME_FIT to set the tradeoff between accuracy and speed. The
                default FAST_FIT should work reasonably well in most instances.
            timeout:
                Time out in seconds for generating and testing rotations. Note
                that the same timeout will apply to generating and testing
                separately. Defaults to 5 mins, which means a total of 10 mins.
            symmetry_tol:
                If > 0, symmetry checking is performed on the two structures
                based on that symmetry tolerance in angstrom. Structures with
                different symmetries are not fitted to each other. A good value
                is around 0.1. Defaults to 0.
        """

        self._tolerance_cell_misfit = tolerance_cell_misfit
        self._tolerance_atomic_misfit = tolerance_atomic_misfit
        self._supercells_allowed = supercells_allowed
        self._anonymized = anonymized
        self._max_rotations = fitting_accuracy
        self._timeout = timeout
        # Sort structures so that they have the same arrangement of species.
        self._structure_a = structure_a.get_sorted_structure()
        self._structure_b = structure_b.get_sorted_structure()

        self._start_time = time.time()

        if symmetry_tol:
            finder_a = SymmetryFinder(self._structure_a, symprec=symmetry_tol)
            finder_b = SymmetryFinder(self._structure_b, symprec=symmetry_tol)
            sg_a = finder_a.get_spacegroup_number()
            sg_b = finder_b.get_spacegroup_number()
            same_sg = sg_a == sg_b
            logger.debug("Spacegroup numbers: A-{}, B-{}".format(sg_a, sg_b))

        if not symmetry_tol or same_sg:
            self._mapping_op = None
            if not self._anonymized:
                self.fit(self._structure_a, self._structure_b)
                if self.fit_found:
                    comp = self._structure_a.composition
                    self.el_mapping = {el: el for el in comp.elements}
            else:
                comp_a = structure_a.composition
                comp_b = structure_b.composition
                if len(comp_a.elements) == len(comp_b.elements):
                    el_a = comp_a.elements
                    #Create permutations of the specie/elements in structure A
                    for p in itertools.permutations(el_a):
                        # Create mapping of the specie/elements in structure B
                        # to that of A.
                        # Then create a modified structure with those elements
                        # and try to fit it.
                        el_mapping = dict(zip(comp_b.elements, p))
                        logger.debug("Using specie mapping " + str(el_mapping))
                        mod = StructureEditor(self._structure_b)
                        mod.replace_species(el_mapping)
                        self.fit(self._structure_a, mod.modified_structure)
                        if self._mapping_op != None:
                            #Store successful element mapping
                            self.el_mapping = {el_a: el_b
                                               for el_b, el_a
                                               in el_mapping.items()}
                            break
                else:
                    logger.debug("No. of elements in structures are unequal. "
                                 "Cannot be fitted!")
        elif symmetry_tol:
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
        # Check composition first.  If compositions are not the same, do not
        # need to fit further.
        if a.composition.reduced_formula != b.composition.reduced_formula or \
           ((a.num_sites != b.num_sites) and not self._supercells_allowed):
            logger.debug('Compositions do not match')
            return None

        logger.debug('Compositions match')

        # This is cheating, but let's try an identity fit first.
        # In many situations, e.g., when a structure has been topotatically
        # delithiated or substituted, you actually can get a very fast answer
        # without having to try all rotations.
        (mapping_op, biggest_dist, oshift) = self.simple_fit(a, b)

        # If we can't do a simple fit, guess we'd have to do the painful
        # process of finding rotations.
        if not mapping_op:
            # Let's sanitize the structures for fitting. This results in more
            # orthogonal cells, and typically much faster fitting.
            a = a.copy(sanitize=True)
            b = b.copy(sanitize=True)

            # Fitting is done by matching sites in one structure (to_fit)
            # to the other (fixed).
            # We set the structure with fewer sites as fixed,
            # and scale the structures to the same density
            (fixed, to_fit) = self._scale_structures(a, b)

            # Defines the atom misfit tolerance
            tol_atoms = self._tolerance_atomic_misfit * \
                        (3 * 0.7405 * fixed.volume / \
                         (4 * math.pi * fixed.num_sites)) ** (1 / 3)
            logger.debug("Atomic misfit tolerance = %.4f" % (tol_atoms))

            max_sites = float('inf')
            # Determine which type of sites to use for the mapping
            for sp in to_fit.species_and_occu:
                sp_sites = [site for site in to_fit \
                            if site.species_and_occu == sp]
                if len(sp_sites) < max_sites:
                    fit_sites = sp_sites
                    max_sites = len(sp_sites)

            # Set the arbitrary origin
            origin = fit_sites[0]
            logger.debug("Origin = " + str(origin))

            oshift = SymmOp.from_rotation_and_translation(np.eye(3),
                                                          - origin.coords)
            shifted_to_fit = apply_operation(to_fit, oshift)

            #Get candidate rotations
            cand_rot = self._get_candidate_rotations(origin, fixed, to_fit)

            #Reset the clock.
            self._start_time = time.time()
            logger.debug(" Found {} candidate rotations".format(len(cand_rot)))
            if len(cand_rot) == 0:
                logger.debug("No candidate rotations found, returning null. ")
                return None

            # sort the operations, the first ones are the ones with small shear
            # this assures us that we find the smallest cell misfit fits
            sorted_cand_rot = sorted(cand_rot.keys(),
                                     key=lambda r: cand_rot[r])

            for rot in sorted_cand_rot:
                (found_map, mapping_op, biggest_dist) = \
                    self._test_rotation(rot, origin, fixed, shifted_to_fit,
                                        tol_atoms)
                if found_map:
                    break
        else:
            logger.debug("Identity fitting matched!")

        logger.debug("Done testing in {} secs".format(time.time() - \
                                                      self._start_time))

        self._atomic_misfit = biggest_dist / ((3 * 0.7405 * a.volume / \
                                (4 * math.pi * a.num_sites)) ** (1 / 3))

        if mapping_op != None:
            rot = mapping_op.rotation_matrix  # maps to_fit to fixed
            p = sqrt_matrix(np.dot(rot.transpose(), rot))
            scale_matrix = np.eye(3) * self.scale
            newrot = np.dot(scale_matrix, rot)
            # we need to now make sure fitterdata.MappingOp maps b -> a and not
            # the other way around
            mshift = np.dot(rot, oshift.translation_vector)
            finaltranslation = mapping_op.translation_vector + mshift[0]
            composite_op = SymmOp.from_rotation_and_translation(
                                                newrot, finaltranslation)
            self._mapping_op = composite_op if self.fixed_is_a \
                                            else composite_op.inverse
            #self._mapping_op = mapping_op
            self._cell_misfit = shear_invariant(p)

    def simple_fit(self, a, b):
        # Fitting is done by matching sites in one structure (to_fit)
        # to the other (fixed).
        # We set the structure with fewer sites as fixed,
        # and scale the structures to the same density
        (fixed, to_fit) = self._scale_structures(a, b)
        # Defines the atom misfit tolerance
        tol_atoms = self._tolerance_atomic_misfit * (3 * 0.7405 * fixed.volume\
                                  / (4 * math.pi * fixed.num_sites)) ** (1 / 3)
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

        oshift = SymmOp.from_rotation_and_translation(np.eye(3),
                                                               - origin.coords)
        shifted_to_fit = apply_operation(to_fit, oshift)

        simple_rots = self._get_simple_rotations(fixed, to_fit)
        for rot in simple_rots:
            rot_op = SymmOp.from_rotation_and_translation(rot,
                                                          np.array([0, 0, 0]))
            (found_map, mapping_op, biggest_dist) = self._test_rotation(rot_op,
                                    origin, fixed, shifted_to_fit, tol_atoms)
            if found_map:
                return (mapping_op, biggest_dist, oshift)

        return (None, None, None)

    def _test_rotation(self, rot, origin, fixed, to_fit, tol_atoms):
        tol_atoms_plus = 1.1 * tol_atoms
        found_map = False
        mapping_op = None
        biggest_dist = 0
        logger.debug("Trying candidate rotation : \n" + str(rot))
        for site in fixed:
            if time.time() - self._start_time > self._timeout:
                logger.debug("Timeout reached when testing rotations.")
                break
            if site.species_and_occu == origin.species_and_occu:
                shift = site.coords
                op = SymmOp.from_rotation_and_translation(
                                                    rot.rotation_matrix, shift)
                nstruct = apply_operation(to_fit, op)
                correspondance = OrderedDict()
                all_match = True
                biggest_dist = 0
                # check to see if transformed struct matches fixed structure
                for trans in nstruct:
                    cands = fixed.get_sites_in_sphere(trans.coords,
                                                      tol_atoms_plus)
                    if len(cands) == 0:
                        logger.debug("No candidates found.")
                        all_match = False
                        break
                    cands = sorted(cands, key=lambda a: a[1])
                    (closest, closest_dist) = cands[0]
                    if closest_dist > tol_atoms or \
                       closest.species_and_occu != trans.species_and_occu:
                        logger.debug("Closest dist too large! " + \
                                     "Closest dist = {}".format(closest_dist))
                        all_match = False
                        break
                    correspondance[trans] = closest

                    if closest_dist > biggest_dist:
                        biggest_dist = closest_dist

                if not all_match:
                    continue

                if not are_sites_unique(correspondance.values(), False):
                    all_match = False
                    continue
                else:
                    for k, v in correspondance.items():
                        logger.debug(str(k) + " fits on " + str(v))

                # now check to see if the converse is true -- do all of the
                # sites of fixed match up with a site in toFit
                # this used to not be here. This fixes a bug.
                logger.debug("Checking inverse mapping")
                inv_correspondance = OrderedDict()

                # it used to be fixed.getNumSites() != nStruct.getNumSites()
                # so only when the number of sites are different but it's
                # actually better to always check the reverse. This
                # elimininates weird situations where two atoms fit to one
                # (reduced in the unit cell)
                for fixed_site in fixed:
                    cands = nstruct.get_sites_in_sphere(fixed_site.coords,
                                                        tol_atoms_plus)
                    if len(cands) == 0:
                        logger.debug("Rejected because inverse mapping does "
                                     "not fit - Step 1")
                        all_match = False
                        break

                    cands = sorted(cands, key=lambda a: a[1])
                    (closest, closest_dist) = cands[0]

                    if closest_dist > tol_atoms or \
                       closest.species_and_occu != fixed_site.species_and_occu:
                        all_match = False
                        logger.debug("Rejected because inverse mapping does "
                                     "not fit - Step 2")
                        break
                    inv_correspondance[fixed_site] = closest

                if all_match:
                    if not are_sites_unique(inv_correspondance.values(),
                                            False):
                        all_match = False
                        logger.debug("Rejected because two atoms fit to the "
                                     "same site for the inverse")
                        continue

                    self.inv_correspondance = inv_correspondance
                    logger.debug("Correspondance for the inverse")
                    for k, v in inv_correspondance.items():
                        logger.debug("{} fits on {}".format(k, v))

                    # The smallest correspondance array shouldn't have any
                    # equivalent sites
                    if fixed.num_sites != to_fit.num_sites:
                        logger.debug("Testing sites unique")
                        if not are_sites_unique(correspondance.values()):
                            all_match = False
                            logger.debug("Rejected because the smallest "
                                         "correspondance array has equivalent"
                                         " sites.")
                            continue

                    found_map = True
                    mapping_op = op
                    self.correspondance = correspondance
                    break

        return (found_map, mapping_op, biggest_dist)

    def __str__(self):

        output = ["Fitting structures"]
        output.append("\nStructure 1:")
        output.append(str(self._structure_a))
        output.append("\nStructure 2:")
        output.append(str(self._structure_b))
        output.append("\nFitting parameters:")
        output.append("\tTolerance cell misfit = " + \
                      str(self._tolerance_cell_misfit))
        output.append("\tSupercells allowed = {}".format(self.
                                                         _supercells_allowed))
        output.append("\tAnonymized  = " + str(self._anonymized))
        output.append("\nFitting " + ("succeeded " if self._mapping_op != None\
                                                   else "failed"))
        if self._mapping_op != None:
            output.append("\tMapping op = " + str(self._mapping_op))
            output.append("\tCell misfit = " + str(self._cell_misfit))
            output.append("\tAtomic misfit = " + str(self._atomic_misfit))
            if self._anonymized:
                output.append("Element mapping")
                output.append("\n".join([str(k) + " -> " + str(v) \
                                         for k, v in self.el_mapping.items()]))
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
        to_fit = Structure(Lattice(to_fit.lattice.matrix * scale),
                           to_fit.species_and_occu, to_fit.frac_coords)
        return (fixed, to_fit)

    def _get_candidate_rotations(self, origin, fixed, to_fit):
        tol_shear = self._tolerance_cell_misfit
        fixed_basis = fixed.lattice.matrix.transpose()
        # need to generate candidate rotations ...
        lengths = fixed.lattice.abc

        shells = []
        for i in range(3):
            dr = lengths[i] * math.sqrt(tol_shear / 2)
            shell = to_fit.get_neighbors_in_shell(origin.coords, lengths[i],
                                                  dr)
            logger.debug("shell {} radius={} dr={}".format(i, lengths[i], dr))
            shells.append([site for (site, dist) in shell
                           if site.species_and_occu ==
                           origin.species_and_occu])
            logger.debug("No. in shell = {}".format(len(shells[-1])))
        # now generate candidate rotations
        cand_rot = {}  # Dict of SymmOp : float
        a = len(shells[0])
        b = len(shells[1])
        c = len(shells[2])
        total_rots = a * b * c
        if total_rots < self._max_rotations:
            logger.debug("Total rots = {}. Using all rotations.".format(
                                                                   total_rots))
            test_rotations = itertools.product(*shells)
        else:
            logger.warning("Total rots={m} exceed max_rotations={n}.".format(
                                        m=total_rots, n=self._max_rotations))

            def random_rot():
                considered_rots = []
                while len(considered_rots) < self._max_rotations:
                    (x, y, z) = [random.randint(0, i - 1) for i in [a, b, c]]
                    if (x, y, z) not in considered_rots:
                        considered_rots.append((x, y, z))
                        yield (shells[0][x], shells[1][y], shells[2][z])
            test_rotations = random_rot()

        for pool in test_rotations:
            if all([nn.species_and_occu == origin.species_and_occu \
                    for nn in pool]):
                # Can a unitary transformation bring the cell vectors together?
                coords = [nn.coords - origin.coords for nn in pool]
                cell_v = np.array(coords).transpose()
                det = np.linalg.det(cell_v)
                if abs(det) < 0.001 or abs(abs(det) - fixed.volume) > 0.01:
                    continue
                rot = np.dot(fixed_basis, np.linalg.inv(cell_v))
                r = SymmOp.from_rotation_and_translation(rot,
                                                         np.array([0, 0, 0]))

                if r not in cand_rot:
                    transf = r.rotation_matrix
                    transf = np.dot(transf.transpose(), transf)

                    # This resolves a very very strange bug in numpy that
                    # causes random eigenvectors to be returned for matrices
                    # very very close to the identity matrix.
                    transf = np.eye(3) if np.allclose(transf, np.eye(3)) \
                                       else transf
                    pbis = sqrt_matrix(transf)
                    if shear_invariant(pbis) < tol_shear:
                        cand_rot[r] = shear_invariant(pbis)
            if time.time() - self._start_time > self._timeout:
                logger.debug("Timeout reached when generating rotations.")
                break

        return cand_rot

    def _get_simple_rotations(self, fixed, to_fit):
        a1 = fixed.lattice.angles
        a2 = to_fit.lattice.angles
        for i in xrange(3):
            if abs(a1[i] - a2[i]) > 5:
                return []
        fixed_unit_matrix = np.array([row / np.linalg.norm(row) \
                                      for row in fixed.lattice.matrix])
        to_fit_unit_matrix = np.array([row / np.linalg.norm(row) \
                                       for row in to_fit.lattice.matrix])
        possible_rots = []
        inv_fixed = np.linalg.inv(fixed_unit_matrix.transpose())
        for p in itertools.permutations(to_fit_unit_matrix):
            possible_rots.append(np.dot(np.array(p).transpose(), inv_fixed))
        return possible_rots

    @property
    def fit_found(self):
        """
        True if a fit was found.
        """
        return self.mapping_op != None

    @property
    def mapping_op(self):
        """
        Operation mapping the two structures.  None if no fit was found.
        """
        return self._mapping_op

    @property
    def structure_a(self):
        """First input structure"""
        return self._structure_a

    @property
    def structure_b(self):
        """Second input structure"""
        return self._structure_b


def apply_operation(structure, symmop):
    """
    Applies a symmetry operation on a structure.

    Args:
        structure:
            Input structure
        symmop:
            SymmOp to apply

    Returns:
        Modified structure after applying symmop.
    """
    editor = StructureEditor(structure)
    editor.apply_operation(symmop)
    return editor.modified_structure


def sqrt_matrix(matrix):
    """
    Calculates sqrt matrix for input matrix.

    Args:
        matrix:
            Input matrix

    Returns:
        Sqrt matrix.
    """
    d, v = np.linalg.eig(matrix)
    diagonalbis = np.sqrt(d) * np.eye(3)
    return np.dot(v, np.dot(diagonalbis, v.transpose()))


def shear_invariant(matrix):
    """
    Calculates the shear invariant for a matrix.

    Args:
        matrix:
            Input matrix

    Returns:
        Shear invariant
    """
    ans = 0
    for i, j in itertools.combinations(xrange(3), 2):
        ans += (matrix[i][i] - matrix[j][j]) ** 2
        ans += 6 * matrix[i][j] ** 2
    return ans


def are_sites_unique(sites, allow_periodic_image=True):
    """
    Checks if the sites are unique.

    Args:
        sites:
            List of sites to check.
        allow_periodic_image:
            Whether to allow periodic images of sites to map onto one another.

    Returns:
        True if sites are unique, i.e., does not map onto each other.
    """
    for (site1, site2) in itertools.combinations(sites, 2):
        if allow_periodic_image and site1.is_periodic_image(site2):
            return False
        if site1.species_and_occu == site2.species_and_occu and \
           (abs(site1.coords - site2.coords) < 0.1).all():
            return False
    return True
