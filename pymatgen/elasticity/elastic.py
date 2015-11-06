from __future__ import absolute_import
from pymatgen.elasticity import voigt_map
from pymatgen.elasticity.tensors import SQTensor
from pymatgen.elasticity.stress import Stress
from pymatgen.elasticity.strain import Strain
import numpy as np
import warnings
from six.moves import range

__author__ = "Maarten de Jong, Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 22, 2012"


class ElasticTensor(SQTensor):
    """
    This class extends SQTensor to describe the 6x6
    elastic tensor, C_{ij}, in Voigt-notation
    """

    def __new__(cls, input_matrix):
        """
        Create an ElasticTensor object.  The constructor throws an error if 
        the shape of the input_matrix argument is not 6x6, i. e. in Voigt-
        notation.  Also issues a warning if the input_matrix argument is
        not symmetric.  Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy
        ndarrays.

        Args:
            input_matrix (6x6 array-like): the Voigt-notation 6x6 array-like
                representing the elastic tensor
        """
        obj = SQTensor(input_matrix).view(cls)
        if obj.shape != (6,6):
            raise ValueError("Default elastic tensor constructor requires "
                             "input argument to be the Voigt-notation 6x6 "
                             "array.  To construct from a 3x3x3x3 array, use "
                             "ElasticTensor.from_full_tensor")
        if not obj.is_symmetric():
            warnings.warn("Elastic tensor input is not symmetric!")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "ElasticTensor({})".format(self.__str__())

    @property
    def compliance_tensor(self):
        """
        returns the compliance tensor, which is the matrix inverse of the
        Voigt-notation elastic tensor
        """
        return self.I

    @property
    def k_voigt(self):
        """
        returns the K_v bulk modulus
        """
        return self[:3, :3].mean()

    @property
    def g_voigt(self):
        """
        returns the G_v shear modulus 
        """
        return (2. * self[:3, :3].trace() - np.triu(self[:3, :3]).sum() +
                3 * self[3:, 3:].trace()) / 15.

    @property
    def k_reuss(self):
        """
        returns the K_r bulk modulus
        """
        return 1. / self.I[:3, :3].sum()

    @property
    def g_reuss(self):
        """
        returns the G_r shear modulus
        """
        return 15. / (8. * self.I[:3, :3].trace() -
                      4. * np.triu(self.I[:3, :3]).sum() +
                      3. * self.I[3:, 3:].trace())

    @property
    def k_vrh(self):
        """
        returns the K_vrh (Voigt-Reuss-Hill) average bulk modulus
        """
        return 0.5 * (self.k_voigt + self.k_reuss)

    @property
    def g_vrh(self):
        """
        returns the G_vrh (Voigt-Reuss-Hill) average shear modulus
        """
        return 0.5 * (self.g_voigt + self.g_reuss)

    @property
    def kg_average(self):
        """
        returns a list of Voigt, Reuss, and Voigt-Reuss-Hill averages of bulk
        and shear moduli similar to legacy behavior
        """
        return [self.k_voigt, self.g_voigt, self.k_reuss, self.g_reuss,
                self.k_vrh, self.g_vrh]

    @property
    def universal_anisotropy(self):
        """
        returns the universal anisotropy value
        """
        return 5. * self.g_voigt / self.g_reuss + \
            self.k_voigt / self.k_reuss - 6.

    @property
    def homogeneous_poisson(self):
        """
        returns the homogeneous poisson ratio
        """
        return (1. - 2. / 3. * self.g_vrh / self.k_vrh) / \
               (2. + 2. / 3. * self.g_vrh / self.k_vrh)

    @property
    def full_tensor(self):
        """
        Returns the tensor in standard notation (i. e. 
        a 4th order 3-dimensional tensor, C_{ijkl}), which
        is represented in a np.array with shape (3,3,3,3)
        """
        c = np.zeros((3, 3, 3, 3))
        for p in range(6):
            for q in range(6):
                i, j = voigt_map[p]
                k, l = voigt_map[q]
                c[i, j, k, l] = c[j, i, k, l] = c[i, j, l, k] = \
                    c[j, i, l, k] = c[k, l, i, j] = self[p, q]
        return c

    @classmethod
    def from_full_tensor(cls, c_ijkl, tol=1e-5):
        """
        Factory method to construct elastic tensor from fourth order
        tensor C_ijkl.  First tests for appropriate symmetries and then 
        constructs the 6x6 voigt notation tensor.

        Args:
            c_ijkl (3x3x3x3 array-like): fourth-order tensor corresponding
                to the full elastic tensor
            tol (float): tolerance for the symmetry test of the tensor
        """
        # Test symmetry of elastic tensor
        c_ijkl = np.array(c_ijkl)
        if not ((c_ijkl - np.transpose(c_ijkl, (1, 0, 2, 3)) < tol).all() and
                (c_ijkl - np.transpose(c_ijkl, (0, 1, 3, 2)) < tol).all() and
                (c_ijkl - np.transpose(c_ijkl, (1, 0, 3, 2)) < tol).all() and
                (c_ijkl - np.transpose(c_ijkl, (3, 2, 0, 1)) < tol).all()):
            raise ValueError("Input elasticity tensor does "
                             "not satisfy necessary symmetries")
        # Construct elastic tensor
        c_pq = np.zeros((6, 6))
        for p in range(6):
            for q in range(6):
                i, j = voigt_map[p]
                k, l = voigt_map[q]
                c_pq[p, q] = c_ijkl[i, j, k, l]

        return cls(c_pq)

    @classmethod
    def from_strain_stress_list(cls, strains, stresses):
        """
        Class method to fit an elastic tensor from stress/strain data.  Method
        uses Moore-Penrose pseudoinverse to invert the s = C*e equation with
        elastic tensor, stress, and strain in voigt notation

        Args:
            stresses (Nx3x3 array-like): list or array of stresses
            strains (Nx3x3 array-like): list or array of strains
        """
        # convert the stress/strain to Nx6 arrays of voigt-notation
        warnings.warn("Linear fitting of Strain/Stress lists may yield "
                      "questionable results from vasp data, use with caution.")
        stresses = np.array([Stress(stress).voigt for stress in stresses])
        with warnings.catch_warnings(record=True):
            strains = np.array([Strain(strain).voigt for strain in strains])

        return cls(np.transpose(np.dot(np.linalg.pinv(strains), stresses)))

    @classmethod
    def from_stress_dict(cls, stress_dict, tol=0.1, vasp=True):
        """
        Constructs the elastic tensor from IndependentStrain-Stress dictionary
        corresponding to legacy behavior of elasticity package.

        Args:
            stress_dict (dict): dictionary of stresses indexed by corresponding
                IndependentStrain objects.
        """
        inds = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
        c_ij = np.array([[np.polyfit([strain[ind1] for strain in list(stress_dict.keys())
                                     if (strain.i, strain.j) == ind1],
                                     [stress_dict[strain][ind2] for strain
                                     in list(stress_dict.keys())
                                     if (strain.i, strain.j) == ind1], 1)[0]
                         for ind1 in inds] for ind2 in inds])
        if vasp:
            c_ij *= -0.1  # Convert units/sign convention of vasp stress tensor
        c_ij[3:, 3:] = 0.5 * c_ij[3:, 3:]  # account for voigt doubling of e4,e5,e6
        c_ij = SQTensor(c_ij)
        c_ij = c_ij.zeroed(tol)
        c_ij = c_ij.symmetrized
        return cls(c_ij)
