from scipy.linalg import polar
import numpy as np

__author__ = "Maarten de Jong, Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 22, 2012"


class SQTensor(np.matrix):
    """
    Class for doing useful general operations on *square* second order tensors, 
    without restrictions on what type (stress, elastic, strain etc.).
    Error is thrown when the class is initialized with non-square matrix.
    """

    def __new__(cls, input_matrix):
        obj = np.asmatrix(input_matrix).view(cls)
        if obj.shape[0] != obj.shape[1]:
            raise ValueError("SQTensor only takes square arrays as input")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "SQTensor({})".format(self.__str__())

    @property
    def T(self):
        """
        shorthand for transpose on SQTensor, addresses issue in np.matrix
        """
        return SQTensor(np.transpose(self))

    @property
    def I(self):
        """
        shorthand for matrix inverse on SQTensor, addresses issue in np.matrix
        """
        return SQTensor(np.linalg.inv(self))

    @property
    def det(self):
        """
        shorthand for the determinant of the SQTensor
        """
        return np.linalg.det(self)

    def is_symmetric(self, tol=1e-5):
        """
        Test to see if tensor is symmetric to a user-defined tolerance.
        This is determined by subtracting the transpose; if any of the 
        resultant elements are above the specified tolerance, returns 
        False.  Otherwise returns true.

        Args:
            tol (float): tolerance to symmetry test
        """
        return (np.abs(self - self.T) < tol).all()

    def is_rotation(self, tol=1e-5):
        """
        Test to see if tensor is a valid rotation matrix, performs a 
        test to check whether the inverse is equal to the transpose
        and if the determinant is equal to one within the specified
        tolerance

        Args:
            tol (float): tolerance to both tests of whether the 
                the determinant is one and the inverse is equal
                to the transpose
        """

        return (np.abs(self.I - self.T) < tol).all() \
            and (np.linalg.det(self) - 1. < tol)

    @property
    def symmetrized(self):
        """
        Returns a symmetrized matrix from the input matrix, 
        calculated by taking the sum of the matrix and its
        transpose
        """
        return 0.5 * (self + self.T)

    def rotate(self, rotation):
        """
        Returns a rotated tensor based on input of a another
        rotation tensor.

        Args:
            rotation (3x3 array-like): rotation tensor, is tested
                for rotation properties and then operates on self
        """
        if self.shape() != (3, 3):
            raise NotImplementedError("Rotations are only implemented for "
                                      "3x3 tensors.")
        rotation = SQTensor(rotation)
        if not rotation.is_rotation():
            raise ValueError("Specified rotation matrix is invalid")
        return rotation * self * rotation.T

    def get_scaled(self, scale_factor):
        """
        Scales the tensor by a certain multiplicative scale factor
        """
        return SQTensor(self * scale_factor)

    @property
    def principal_invariants(self):
        """
        Returns a list of principal invariants for the tensor,
        which are the values of the coefficients of the characteristic 
        polynomial for the matrix
        """
        # TODO: JM asks whether this fulfills the necessary sign conventions
        return np.poly(self)[1:]

    @property
    def polar_decomposition(self, side='right'):
        """
        calculates matrices for polar decomposition
        """
        return polar(self, side=side)

    def zeroed(self, tol):
        """
        returns the matrix with all entries below a certain threshold
        (i.e. tol) set to zero
        """
        new_tensor = self.copy()
        new_tensor[new_tensor < tol] = 0
        return new_tensor


if __name__ == "__main__":
    import doctest

    doctest.testmod()

    eye = np.identity(3)
    sigma = SQTensor(np.random.randn(3, 3))
