#!/usr/bin/env python

"""
This module implements a point group assigner for a molecule.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "5/8/13"

import logging
import itertools
from collections import defaultdict

import numpy as np
try:
    import scipy.cluster as spcluster
except ImportError:
    spcluster = None

from pymatgen.core.operations import SymmOp
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.util.decorators import requires

logger = logging.getLogger(__name__)


class PointGroup(list):
    """
    Defines a point group, which is essentially a sequence of symmetry
    operations.

    .. attribute:: sch_symbol

        Schoenflies symbol of the point group.
    """
    def __init__(self, sch_symbol, operations, tol=0.1):
        """
        Args:
            sch_symbol:
                The schoenflies symbol of the point group.
            operations:
                An initial set of symmetry operations. It is sufficient to
                provide only just enough operations to generate the full set
                of symmetries.
            tol:
                Tolerance to generate the full set of symmetry operations.
        """
        self.sch_symbol = sch_symbol
        super(PointGroup, self).__init__(
            generate_full_symmops(operations, tol))

    def __str__(self):
        return self.sch_symbol

    def __repr__(self):
        return self.__str__()


@requires(spcluster is not None, "Cannot import scipy. PointGroupAnalyzer "
                                 "requires scipy.cluster")
class PointGroupAnalyzer(object):
    """
    A class to analyze the point group of a molecule. The general outline of
    the algorithm is as follows:

    1. Center the molecule around its center of mass.
    2. Compute the inertia tensor and the eigenvalues and eigenvectors.
    3. Handle the symmetry detection based on eigenvalues.

        a. Linear molecules have one zero eigenvalue. Possible symmetry
           operations are C*v or D*v
        b. Asymetric top molecules have all different eigenvalues. The
           maximum rotational symmetry in such molecules is 2
        c. Symmetric top molecules have 1 unique eigenvalue, which gives a
           unique rotation axis.  All axial point groups are possible
           except the cubic groups (T & O) and I.
        d. Spherical top molecules have all three eigenvalues equal. They
           have the rare T, O or I point groups.

    .. attribute:: sch_symbol

        Schoenflies symbol of the detected point group.
    """
    inversion_op = SymmOp.inversion()

    def __init__(self, mol, tolerance=0.3, eigen_tolerance=0.01,
                 matrix_tol=0.1):
        """
        The default settings are usually sufficient.

        Args:
            mol:
                Molecule
            tolerance:
                Distance tolerance to consider sites as symmetrically
                equivalent. Defaults to 0.3 Angstrom.
            eigen_tolerance:
                Tolerance to compare eigen values of the inertia tensor.
                Defaults to 0.01.
            matrix_tol:
                Tolerance used to generate the full set of symmetry
                operations of the point group.
        """
        self.mol = mol
        self.centered_mol = mol.get_centered_molecule()
        self.tol = tolerance
        self.eig_tol = eigen_tolerance
        self.mat_tol = matrix_tol
        self._analyze()

    def _analyze(self):
        if len(self.centered_mol) == 1:
            self.sch_symbol = "Kh"
        else:
            inertia_tensor = np.zeros((3, 3))
            total_inertia = 0
            for site in self.mol:
                c = site.coords
                wt = site.species_and_occu.weight
                for i in xrange(3):
                    inertia_tensor[i, i] += wt * (c[(i + 1) % 3] ** 2
                                                  + c[(i + 2) % 3] ** 2)
                for i, j in itertools.combinations(xrange(3), 2):
                    inertia_tensor[i, j] += -wt * c[i] * c[j]
                    inertia_tensor[j, i] += -wt * c[j] * c[i]
                total_inertia += wt * np.dot(c, c)

            # Normalize the inertia tensor so that it does not scale with size
            # of the system.  This mitigates the problem of choosing a proper
            # comparison tolerance for the eigenvalues.
            inertia_tensor /= total_inertia
            eigvals, eigvecs = np.linalg.eig(inertia_tensor)
            self.principal_axes = eigvecs.T
            self.eigvals = eigvals
            v1, v2, v3 = eigvals
            eig_zero = abs(v1 * v2 * v3) < self.eig_tol ** 3
            eig_all_same = abs(v1 - v2) < self.eig_tol and abs(
                v1 - v3) < self.eig_tol
            eig_all_diff = abs(v1 - v2) > self.eig_tol and abs(
                v1 - v2) > self.eig_tol and abs(v2 - v3) > self.eig_tol

            self.rot_sym = []
            self.symmops = [SymmOp(np.eye(4))]

            if eig_zero:
                logger.debug("Linear molecule detected")
                self._proc_linear()
            elif eig_all_same:
                logger.debug("Spherical top molecule detected")
                self._proc_sph_top()
            elif eig_all_diff:
                logger.debug("Asymmetric top molecule detected")
                self._proc_asym_top()
            else:
                logger.debug("Symmetric top molecule detected")
                self._proc_sym_top()

    def _proc_linear(self):
        if self.is_valid_op(PointGroupAnalyzer.inversion_op):
            self.sch_symbol = "D*h"
            self.symmops.append(PointGroupAnalyzer.inversion_op)
        else:
            self.sch_symbol = "C*v"

    def _proc_asym_top(self):
        """
        Handles assymetric top molecules, which cannot contain rotational
        symmetry larger than 2.
        """
        self._check_R2_axes_asym()
        if len(self.rot_sym) == 0:
            logger.debug("No rotation symmetries detected.")
            self._proc_no_rot_sym()
        elif len(self.rot_sym) == 3:
            logger.debug("Dihedral group detected.")
            self._proc_dihedral()
        else:
            logger.debug("Cyclic group detected.")
            self._proc_cyclic()

    def _proc_sym_top(self):
        """
        Handles symetric top molecules which has one unique eigenvalue whose
        corresponding principal axis is a unique rotational axis.  More complex
        handling required to look for R2 axes perpendicular to this unique
        axis.
        """
        if abs(self.eigvals[0] - self.eigvals[1]) < self.eig_tol:
            ind = 2
        elif abs(self.eigvals[1] - self.eigvals[2]) < self.eig_tol:
            ind = 0
        else:
            ind = 1

        unique_axis = self.principal_axes[ind]
        self._check_rot_sym(unique_axis)
        if len(self.rot_sym) > 0:
            self._check_perpendicular_r2_axis(unique_axis)

        if len(self.rot_sym) >= 2:
            self._proc_dihedral()
        elif len(self.rot_sym) == 1:
            self._proc_cyclic()
        else:
            self._proc_no_rot_sym()

    def _proc_no_rot_sym(self):
        """
        Handles molecules with no rotational symmetry. Only possible point
        groups are C1, Cs and Ci.
        """
        self.sch_symbol = "C1"
        if self.is_valid_op(PointGroupAnalyzer.inversion_op):
            self.sch_symbol = "Ci"
            self.symmops.append(PointGroupAnalyzer.inversion_op)
        else:
            for v in self.principal_axes:
                mirror_type = self._find_mirror(v)
                if not mirror_type == "":
                    self.sch_symbol = "Cs"
                    break

    def _proc_cyclic(self):
        """
        Handles cyclic group molecules.
        """
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        self.sch_symbol = "C{}".format(rot)
        mirror_type = self._find_mirror(main_axis)
        if mirror_type == "h":
            self.sch_symbol += "h"
        elif mirror_type == "v":
            self.sch_symbol += "v"
        elif mirror_type == "":
            if self.is_valid_op(SymmOp.rotoreflection(main_axis,
                                                      angle=180 / rot)):
                self.sch_symbol = "S{}".format(2 * rot)

    def _proc_dihedral(self):
        """
        Handles dihedral group molecules, i.e those with intersecting R2 axes
        and a main axis.
        """
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        self.sch_symbol = "D{}".format(rot)
        mirror_type = self._find_mirror(main_axis)
        if mirror_type == "h":
            self.sch_symbol += "h"
        elif not mirror_type == "":
            self.sch_symbol += "d"

    def _check_R2_axes_asym(self):
        """
        Test for 2-fold rotation along the principal axes. Used to handle
        asymetric top molecules.
        """
        for v in self.principal_axes:
            op = SymmOp.from_axis_angle_and_translation(v, 180)
            if self.is_valid_op(op):
                self.symmops.append(op)
                self.rot_sym.append((v, 2))

    def _find_mirror(self, axis):
        """
        Looks for mirror symmetry of specified type about axis.  Possible
        types are "h" or "vd".  Horizontal (h) mirrors are perpendicular to
        the axis while vertical (v) or diagonal (d) mirrors are parallel.  v
        mirrors has atoms lying on the mirror plane while d mirrors do
        not.
        """
        mirror_type = ""

        #First test whether the axis itself is the normal to a mirror plane.
        if self.is_valid_op(SymmOp.reflection(axis)):
            self.symmops.append(SymmOp.reflection(axis))
            mirror_type = "h"
        else:
            # Iterate through all pairs of atoms to find mirror
            for s1, s2 in itertools.combinations(self.centered_mol, 2):
                if s1.species_and_occu == s2.species_and_occu:
                    normal = s1.coords - s2.coords
                    if np.dot(normal, axis) < self.tol:
                        op = SymmOp.reflection(normal)
                        if self.is_valid_op(op):
                            self.symmops.append(op)
                            if len(self.rot_sym) > 1:
                                mirror_type = "d"
                                for v, r in self.rot_sym:
                                    if not np.linalg.norm(v - axis) < self.tol:
                                        if np.dot(v, normal) < self.tol:
                                            mirror_type = "v"
                                            break
                            else:
                                mirror_type = "v"
                            break

        return mirror_type

    def _get_smallest_set_not_on_axis(self, axis):
        """
        Returns the smallest list of atoms with the same species and
        distance from origin AND does not lie on the specified axis.  This
        maximal set limits the possible rotational symmetry operations,
        since atoms lying on a test axis is irrelevant in testing rotational
        symmetryOperations.
        """
        def not_on_axis(site):
            v = np.cross(site.coords, axis)
            return np.linalg.norm(v) > self.tol

        valid_sets = []
        origin_site, dist_el_sites = cluster_sites(self.centered_mol, self.tol)
        for test_set in dist_el_sites.values():
            valid_set = filter(not_on_axis, test_set)
            if len(valid_set) > 0:
                valid_sets.append(valid_set)

        return min(valid_sets, key=lambda s: len(s))

    def _check_rot_sym(self, axis):
        """
        Determines the rotational symmetry about supplied axis.  Used only for
        symmetric top molecules which has possible rotational symmetry
        operations > 2.
        """
        min_set = self._get_smallest_set_not_on_axis(axis)
        max_sym = len(min_set)
        for i in xrange(max_sym, 0, -1):
            if max_sym % i != 0:
                continue
            op = SymmOp.from_axis_angle_and_translation(axis, 360 / i)
            rotvalid = self.is_valid_op(op)
            if rotvalid:
                self.symmops.append(op)
                self.rot_sym.append((axis, i))
                return i
        return 1

    def _check_perpendicular_r2_axis(self, axis):
        """
        Checks for R2 axes perpendicular to unique axis.  For handling
        symmetric top molecules.
        """
        min_set = self._get_smallest_set_not_on_axis(axis)
        for s1, s2 in itertools.combinations(min_set, 2):
            test_axis = np.cross(s1.coords - s2.coords, axis)
            if np.linalg.norm(test_axis) > self.tol:
                op = SymmOp.from_axis_angle_and_translation(test_axis, 180)
                r2present = self.is_valid_op(op)
                if r2present:
                    self.symmops.append(op)
                    self.rot_sym.append((test_axis, 2))
                    return True

    def _proc_sph_top(self):
        """
        Handles Sperhical Top Molecules, which belongs to the T, O or I point
        groups.
        """
        self._find_spherical_axes()
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        if len(self.rot_sym) == 0 or rot < 3:
            logger.debug("Accidental speherical top!")
            self._proc_sym_top()
        elif rot == 3:
            mirror_type = self._find_mirror(main_axis)
            if mirror_type != "":
                if self.is_valid_op(PointGroupAnalyzer.inversion_op):
                    self.symmops.append(PointGroupAnalyzer.inversion_op)
                    self.sch_symbol = "Th"
                else:
                    self.sch_symbol = "Td"
            else:
                self.sch_symbol = "T"
        elif rot == 4:
            if self.is_valid_op(PointGroupAnalyzer.inversion_op):
                self.symmops.append(PointGroupAnalyzer.inversion_op)
                self.sch_symbol = "Oh"
            else:
                self.sch_symbol = "O"
        elif rot == 5:
            if self.is_valid_op(PointGroupAnalyzer.inversion_op):
                self.symmops.append(PointGroupAnalyzer.inversion_op)
                self.sch_symbol = "Ih"
            else:
                self.sch_symbol = "I"

    def _find_spherical_axes(self):
        """
        Looks for R5, R4, R3 and R2 axes in speherical top molecules.  Point
        group T molecules have only one unique 3-fold and one unique 2-fold
        axis. O molecules have one unique 4, 3 and 2-fold axes. I molecules
        have a unique 5-fold axis.
        """
        rot_present = defaultdict(bool)
        origin_site, dist_el_sites = cluster_sites(self.centered_mol, self.tol)
        test_set = min(dist_el_sites.values(), key=lambda s: len(s))
        coords = [s.coords for s in test_set]
        for c1, c2, c3 in itertools.combinations(coords, 3):
            for cc1, cc2 in itertools.combinations([c1, c2, c3], 2):
                if not rot_present[2]:
                    test_axis = cc1 + cc2
                    if np.linalg.norm(test_axis) > self.tol:
                        op = SymmOp.from_axis_angle_and_translation(test_axis,
                                                                    180)
                        rot_present[2] = self.is_valid_op(op)
                        if rot_present[2]:
                            self.symmops.append(op)
                            self.rot_sym.append((test_axis, 2))

            test_axis = np.cross(c2 - c1, c3 - c1)
            if np.linalg.norm(test_axis) > self.tol:
                for r in (3, 4, 5):
                    if not rot_present[r]:
                        op = SymmOp.from_axis_angle_and_translation(
                            test_axis, 360/r)
                        rot_present[r] = self.is_valid_op(op)
                        if rot_present[r]:
                            self.symmops.append(op)
                            self.rot_sym.append((test_axis, r))
                            break
            if rot_present[2] and rot_present[3] and (
                    rot_present[4] or rot_present[5]):
                break

    def get_pointgroup(self):
        """
        Returns a PointGroup object for the molecule.
        """
        return PointGroup(self.sch_symbol, self.symmops, self.mat_tol)

    def is_valid_op(self, symmop):
        """
        Check if a particular symmetry operation is a valid symmetry operation
        for a molecule, i.e., the operation maps all atoms to another
        equivalent atom.

        Args:
            symmop:
                Symmetry op to test.
        """
        coords = self.centered_mol.cart_coords
        for site in self.centered_mol:
            coord = symmop.operate(site.coords)
            ind = find_in_coord_list(coords, coord, self.tol)
            if not (len(ind) == 1 and
                    self.centered_mol[ind[0]].species_and_occu
                    == site.species_and_occu):
                return False
        return True


@requires(spcluster is not None, "Cannot import scipy. cluster_sites require "
                                 "scipy.cluster.")
def cluster_sites(mol, tol):
    """
    Cluster sites based on distance and species type.

    Args:
        mol:
            Molecule (should be centered at center of mass).
        tol:
            Tolerance to use.

    Returns:
        (origin_site, clustered_sites). origin_site is a site at the center
        of mass (None if there are no origin atoms). clustered_sites is a
        dict of {(avg_dist, species_and_occu): [list of sites]}
    """
    # Cluster works for dim > 2 data. We just add a dummy 0 for second
    # coordinate.
    dists = [[np.linalg.norm(site.coords), 0] for site in mol]
    f = spcluster.hierarchy.fclusterdata(dists, tol, criterion='distance')
    clustered_dists = defaultdict(list)
    for i, site in enumerate(mol):
        clustered_dists[f[i]].append(dists[i])
    avg_dist = {label: np.mean(val) for label, val in clustered_dists.items()}
    clustered_sites = defaultdict(list)
    origin_site = None
    for i, site in enumerate(mol):
        if avg_dist[f[i]] < tol:
            origin_site = site
        else:
            clustered_sites[(avg_dist[f[i]],
                             site.species_and_occu)].append(site)
    return origin_site, clustered_sites


def generate_full_symmops(symmops, tol):
    """
    Recursive algorithm to permute through all possible combinations of the
    initially supplied symmetry operations to arrive at a complete set of
    operations mapping a single atom to all other equivalent atoms in the
    point group.  This assumes that the initial number already uniquely
    identifies all operations.

    Args:
        symmops:
            Initial set of symmetry operations.

    Returns:
        Full set of symmetry operations.
    """

    a = [o.affine_matrix for o in symmops]

    if len(symmops) > 300:
        logger.debug("Generation of symmetry operations in infinite loop.  " +
                     "Possible error in initial operations or tolerance too "
                     "low.")
    else:
        for op1, op2 in itertools.product(symmops, symmops):
            m = np.dot(op1.affine_matrix, op2.affine_matrix)
            d = np.abs(a - m) < tol
            if not np.any(np.all(np.all(d, axis=2), axis=1)):
                return generate_full_symmops(symmops + [SymmOp(m)], tol)

    return symmops
