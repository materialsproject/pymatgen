# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to handle thermal displacement matrices (anisotropic displacement parameters).
"""

import numpy as np
from monty.json import MSONable

__author__ = "J. George"
__copyright__ = "Copyright 2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "J. George"
__email__ = "janine.george@bam.de"
__status__ = "Testing"
__date__ = "August 09, 2022"

try:
    import phonopy
except ImportError as ex:
    print(ex)
    phonopy = None
from pymatgen.io.cif import CifWriter


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
            thermal_displacement: 1d numpy array

        Returns:
            2d numpy array including thermal displacements
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
