import warnings
import sys
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
import os
__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 13, 2012"


class Strain(SQTensor):
    """
    Subclass of SQTensor that describes the strain matrix
    """
    # TODO: AJ says strain must subclass SQTensor
    # TODO: AJ says much of this class should be reimplemented from the standpoint of subclassing SQTensor
    # TODO: AJ says there should be a static from_strain(matrix) method that 
    #           constructs the object from a strain matrix rather than deformation matrix
    #           For an example of the above, see the various 'from_xxxxx' methods in 
    #           pymatgen.core.structure.Composition
    # TODO: JM says we might want to have the constructor use the strain matrix, and have a from_deformation_matrix
    #   method instead of a from_strain method
    
    def __init__(self, strain_matrix, dfm=None):
        self._dfm = None
        super(Strain,self).__init__(strain_matrix)

    @classmethod
    def from_deformation(cls, deformation_matrix):
        """
        constructor that returns a Strain object from a deformation
        matrix, rather than a strain matrix
        """
        return cls(0.5*(self._matrix*self._matrix.T - np.eye(3)),
                  dfm = deformation_matrix)

    # return Green-Lagrange strain matrix
    @property
    def strain(self):
        return self._strain

    @property
    def deformation_matrix(self):
        return self._dfm

# TODO: JM asks whether this method should be implemented
#    # construct def. matrix from indices and amount
#    @staticmethod
#    def from_ind_amt_dfm(matrixpos, amt):
#        F = np.identity(3)
#        F[matrixpos] = F[matrixpos] + amt
#        return Strain(F)

    def __eq__(self, other):
        df, df2 = self.deformation_matrix, other.deformation_matrix
        #TODO: AJ asks can't you just do return self.deformation_matrix == other.deformation_matrix?
        for i, row in enumerate(df):
            for j, item in enumerate(row):
                if np.any(df[i][j] !=df2[i][j]):
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        #for now, just use a sum of elements * row * column^2
        df = self.deformation_matrix.tolist()
        h_sum = 0.0
        for i, row in enumerate(df):
            for j, item in enumerate(row):
                h_sum += item * (i + 1) * (j + 1) * (j + 1) * (i + 2)**3
        return h_sum


class IndependentStrain(Strain):
   # todo: add polar decomposition method
   #

    def __init__(self, deformation,tol=0.00000001):

        super(IndependentStrain, self).__init__(deformation)
        (self._i, self._j) = self.check_F(tol)

    @staticmethod
    def from_ind_amt_dfm(matrixpos, amt):
        F = np.identity(3)
        F[matrixpos] = F[matrixpos] + amt
        return IndependentStrain(F)

    def check_F(self, tol=0.00001):
        df1 = self.deformation_matrix
        counter = 0
        checkmatrix = np.zeros((3,3))

        for c1 in range(0,3):
            for c2 in range(0,3):
                if c1 != c2:
                    if np.abs(df1[c1,c2]) > tol:
                        checkmatrix[c1,c2] = 1
                        counter = counter + 1
                else:
                    if np.abs(df1[c1,c2]-1) > tol:
                        checkmatrix[c1,c2] = 1
                        counter = counter + 1

        if counter != 1:
            raise ValueError("One independent deformation must be applied")

        return (checkmatrix.nonzero()[0][0], checkmatrix.nonzero()[1][0])


    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j


if __name__ == "__main__":

    mat = np.eye(3)
    mat[0,1] = 0.001
#    print mat

    my_strain = IndependentStrain(mat)
    my_strain.check_F()


#    print my_strain._strain
    
    
    
#    print type(mat)

#    print my_strain.deformation_matrix
#    print my_strain.strain

#    my_strain2 = IndependentStrain(mat)
#    print my_strain2.__dict__.keys()
#    print my_strain2.__hash__()

#    print my_strain2._j
#    print my_strain2.check_F()
#    my_strain2.checkF
#    print my_strain.__dict__.keys()
#    print my_strain.deformation_matrix
#    print my_strain.strain
#    my_strain.index
#    my_scaled_strain = my_strain.get_scaled(1.05)
#    print my_scaled_strain.deformation_matrix
#    print my_scaled_strain.strain
#    print my_strain == my_scaled_strain
#    mat2 = np.eye(3)
#    mat2[0,0] = 1.01
#    my_strain2 = Strain(mat)
#    print my_strain == my_strain2


