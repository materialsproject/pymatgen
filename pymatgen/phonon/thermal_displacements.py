# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to handle thermal displacement matrices (anisotropic displacement parameters).
"""

import numpy as np
from monty.json import MSONable

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.cif import CifWriter

try:
    import phonopy
except ImportError as ex:
    print(ex)
    phonopy = None

__author__ = "J. George"
__copyright__ = "Copyright 2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "J. George"
__email__ = "janine.george@bam.de"
__status__ = "Testing"
__date__ = "August 09, 2022"


class ThermalDisplacementMatrices(MSONable):
    """
    Class to handle thermal displacement matrices
    This class stores thermal displacement matrices in Ucart format

    An earlier implementation based on Matlab can be found here:
    https://github.com/JaGeo/MolecularToolbox
    ( J. George, A. Wang, V. L. Deringer, R. Wang, R. Dronskowski, U. Englert, CrystEngComm, 2015, 17, 7414–7422.)
    """

    def __init__(self, thermal_displacement_matrix_cart, structure, temperature, thermal_displacement_matrix_cif=None):
        """
        Args:
            thermal_displacement_matrix: 2D numpy array including the thermal_displacement matrix Ucart
                1st dimension atom types, then compressed thermal displacement matrix will follow
                 U11, U22, U33, U23, U13, U12 (xx, yy, zz, yz, xz, xy)
                 convention similar to "thermal_displacement_matrices.yaml" in phonopy
            structure: A pymatgen Structure object
            temperature: temperature at which thermal displacement matrix was determined
            thermal_displacement_matrix: 2D numpy array including the thermal_displacement matrix Ucif format
                1st dimension atom types, then compressed thermal displacement matrix will follow
                 U11, U22, U33, U23, U13, U12 (xx, yy, zz, yz, xz, xy)
                 convention similar to "thermal_displacement_matrices.yaml" in phonopy
        """
        self.thermal_displacement_matrix_cart = np.array(thermal_displacement_matrix_cart)
        self.structure = structure
        self.temperature = temperature
        if thermal_displacement_matrix_cif is not None:
            self.thermal_displacement_matrix_cif = np.array(thermal_displacement_matrix_cif)
        else:
            self.thermal_displacement_matrix_cif = None

        # convert to matrix form easier handling of the data
        self.thermal_displacement_matrix_cart_matrixform = ThermalDisplacementMatrices.get_full_matrix(
            self.thermal_displacement_matrix_cart
        )

        if self.thermal_displacement_matrix_cif is not None:
            self.thermal_displacement_matrix_cif_matrixform = ThermalDisplacementMatrices.get_full_matrix(
                self.thermal_displacement_matrix_cif
            )

    @staticmethod
    def get_full_matrix(thermal_displacement):
        """
        transfers the reduced matrix to the full matrix (order of reduced matrix U11, U22, U33, U23, U13, U12)
        Args:
            thermal_displacement: 2d numpy array, first dimension are the atoms

        Returns:
            3d numpy array including thermal displacements, first dimensions are the atoms
        """
        matrixform = np.zeros((len(thermal_displacement), 3, 3))
        for imat, mat in enumerate(thermal_displacement):
            # xx, yy, zz, yz, xz, xy
            matrixform[imat][0][0] = mat[0]
            matrixform[imat][1][1] = mat[1]
            matrixform[imat][2][2] = mat[2]
            matrixform[imat][1][2] = mat[3]
            matrixform[imat][2][1] = mat[3]
            matrixform[imat][0][2] = mat[4]
            matrixform[imat][2][0] = mat[4]
            matrixform[imat][0][1] = mat[5]
            matrixform[imat][1][0] = mat[5]
        return matrixform

    @staticmethod
    def get_reduced_matrix(thermal_displacement):
        """
        transfers the full matrix to reduced matrix (order of reduced matrix U11, U22, U33, U23, U13, U12)
        Args:
            thermal_displacement: 2d numpy array, first dimension are the atoms

        Returns:
            3d numpy array including thermal displacements, first dimensions are the atoms
        """
        reduced_matrix = np.zeros((len(thermal_displacement), 6))
        for imat, mat in enumerate(thermal_displacement):
            # xx, yy, zz, yz, xz, xy
            reduced_matrix[imat][0] = mat[0][0]
            reduced_matrix[imat][1] = mat[1][1]
            reduced_matrix[imat][2] = mat[2][2]
            reduced_matrix[imat][3] = mat[1][2]
            reduced_matrix[imat][4] = mat[0][2]
            reduced_matrix[imat][5] = mat[0][1]
        return reduced_matrix

    @property
    def Ustar(self):
        """
        computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477–480.
        Returns: Ustar as a numpy array, first dimension are the atoms in the structure
        """
        A = self.structure.lattice.matrix.T
        Ainv = np.linalg.inv(A)
        Ustar = []
        for mat in self.thermal_displacement_matrix_cart_matrixform:
            mat_ustar = np.dot(np.dot(Ainv, mat), Ainv.T)
            Ustar.append(mat_ustar)
        return np.array(Ustar)

    @property
    def Ucif(self):
        """
        computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477–480.
        Returns: Ucif as a numpy array, first dimension are the atoms in the structure
        """
        if self.thermal_displacement_matrix_cif is None:
            # computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477–480.
            # will compute Ucif based on Ustar for each atom in the structure

            A = self.structure.lattice.matrix.T  # check this again?
            N = np.diag([np.linalg.norm(x) for x in np.linalg.inv(A)])
            Ninv = np.linalg.inv(N)
            Ucif = []
            Ustar = self.Ustar
            for mat in Ustar:
                mat_cif = np.dot(np.dot(Ninv, mat), Ninv.T)
                Ucif.append(mat_cif)
            return np.array(Ucif)
        else:
            return self.thermal_displacement_matrix_cif_matrixform

    @property
    def B(self):
        """
        computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477–480.
        Returns: B as a numpy array, first dimension are the atoms in the structure
        """

        B = []
        for mat in self.Ucif:
            mat_B = mat * 8 * np.pi**2
            B.append(mat_B)
        return np.array(B)

    @property
    def beta(self):
        """
        computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477–480.
        Returns: beta as a numpy array, first dimension are the atoms in the structure
        """
        # will compute beta based on Ustar
        beta = []
        for mat in self.Ustar:
            mat_beta = mat * 2 * np.pi**2
            beta.append(mat_beta)
        return beta

    @property
    def U1U2U3(self):
        """
        computation as described in R. W. Grosse-Kunstleve, P. D. Adams, J Appl Cryst 2002, 35, 477–480.
        Returns: numpy array of eigenvalues of Ucart,  first dimension are the atoms in the structure
        """
        U1U2U3 = []
        for mat in self.thermal_displacement_matrix_cart_matrixform:
            U1U2U3.append(np.linalg.eig(mat)[0])
        return U1U2U3

    def write_cif(self, filename):
        """
        writes a cif including thermal displacements
        Args:
            filename: name of the cif file
        """
        w = CifWriter(self.structure)
        w.write_file(filename)
        # This will simply append the thermal displacement part to the cif from the cifwriter
        # In the long run, cifwriter could be extended to handle thermal displacement matrices
        with open(filename, "a") as file:
            file.write("loop_ \n")
            file.write("_atom_site_aniso_label\n")
            file.write("_atom_site_aniso_U_11\n")
            file.write("_atom_site_aniso_U_22\n")
            file.write("_atom_site_aniso_U_33\n")
            file.write("_atom_site_aniso_U_23\n")
            file.write("_atom_site_aniso_U_13\n")
            file.write("_atom_site_aniso_U_12\n")
            file.write(f"# Additional Data for U_Aniso: {self.temperature}\n")

            count = 0
            for site, matrix in zip(self.structure, self.Ucif):
                file.write(
                    f"{site.specie.symbol}{count} {matrix[0][0]} {matrix[1][1]} {matrix[2][2]}"
                    f" {matrix[1][2]} {matrix[0][2]} {matrix[0][1]}\n"
                )
                count += 1

    @staticmethod
    def _angle_dot(a, b):
        dot_product = np.dot(a, b)
        prod_of_norms = np.linalg.norm(a) * np.linalg.norm(b)
        divided = dot_product / prod_of_norms
        angle_rad = np.arccos(np.round(divided, 10))
        angle = np.degrees(angle_rad)
        return angle

    def compute_directionality_quality_criterion(self, other):
        """
        Will compute directionality of prolate displacement ellipsoids as described in
        https://doi.org/10.1039/C9CE00794F with the earlier implementation: https://github.com/damMroz/Angle/
        Args:
            other: ThermalDisplacementMatrix
            please make sure that the order of the atoms in both objects that are compared
            is the same. Otherwise, this analysis will deliver wrong results

        Returns:
            will return a list including dicts for each atom that include "vector0"
            (largest principal axes of self object),
             "vector1" (largest principal axes of the other object), "angle" between both axes,
              These vectors can then, for example, be drawn into the structure with VESTA.
              Vectors are given in Cartesian coordinates
        """

        # compare the atoms string at least
        for spec1, spec2 in zip(self.structure.species, other.structure.species):
            if spec1 != spec2:
                raise ValueError(
                    "Species in both structures are not the same! "
                    "Please use structures that are similar to each other"
                )
        # check if structures match
        structure_match = StructureMatcher()
        if not structure_match.fit(struct1=self.structure, struct2=other.structure):
            raise ValueError("Structures have to be similar")

        results = []
        for self_Ucart, other_Ucart in zip(
            self.thermal_displacement_matrix_cart_matrixform, other.thermal_displacement_matrix_cart_matrixform
        ):
            result_dict = {}

            # determine eigenvalues and vectors for inverted Ucart
            invUcart_eig_self, invUcart_eigv_self = np.linalg.eig(np.linalg.inv(self_Ucart))
            invUcart_eig_other, invUcart_eigv_other = np.linalg.eig(np.linalg.inv(other_Ucart))

            argmin_self = np.argmin(invUcart_eig_self)
            vec_self = invUcart_eigv_self.transpose()[argmin_self]
            argmin_other = np.argmin(invUcart_eig_other)
            vec_other = invUcart_eigv_other.transpose()[argmin_other]
            # vector direction does not matter here, smallest angle should be given
            result_dict["angle"] = np.min(
                [self._angle_dot(vec_self, vec_other), self._angle_dot(vec_self, vec_other * -1)]
            )
            result_dict["vector0"] = vec_self
            result_dict["vector1"] = vec_other

            results.append(result_dict)

        return results

    @property
    def ratio_prolate(self):
        """
        This will compute ratio between largest eigenvalue of Ucart and smallest one
        Returns:
        """
        ratios = []
        for us in self.U1U2U3:
            ratios.append(np.max(us) / np.min(us))

        return np.array(ratios)

    @staticmethod
    def from_Ucif(thermal_displacement_matrix_cif, structure, temperature):
        """
        starting from a numpy array, it will convert Ucif values into Ucart values and initialize the class
        Args:
            thermal_displacement_matrix_cif: np.array,
                first dimension are the atoms,
                then reduced form of thermal displacement matrix will follow
                Order as above: U11, U22, U33, U23, U13, U12
            structure: Structure object
            temperature: float
                Corresponding temperature

        Returns:
            ThermalDisplacementMatrices
        """
        # get matrix form
        thermal_displacement_matrix_cif_matrix_form = ThermalDisplacementMatrices.get_full_matrix(
            thermal_displacement_matrix_cif
        )

        # convert the parameters for each atom
        A = structure.lattice.matrix.T
        N = np.diag([np.linalg.norm(x) for x in np.linalg.inv(A)])
        Ucart = []
        for mat in thermal_displacement_matrix_cif_matrix_form:
            mat_ustar = np.dot(np.dot(N, mat), N.T)
            mat_ucart = np.dot(np.dot(A, mat_ustar), A.T)
            Ucart.append(mat_ucart)

        thermal_displacement_matrix_cart = ThermalDisplacementMatrices.get_reduced_matrix(np.array(Ucart))

        # get ThermalDisplacementMatrices Object

        return ThermalDisplacementMatrices(
            thermal_displacement_matrix_cart=thermal_displacement_matrix_cart,
            thermal_displacement_matrix_cif=thermal_displacement_matrix_cif,
            structure=structure,
            temperature=temperature,
        )
