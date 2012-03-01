import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo') # (If one does not want to change $PYTHONPATH)
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

class Strain(object):
 
    def __init__(self, deformation_matrix):
        self._dfm = deformation_matrix
        self._strain = 0.5 * (np.matrix(self._dfm) * np.transpose(np.matrix(self._dfm)) - np.eye(3))

    # return a scaled version of this matrix
    def get_scaled(self, scale_factor):
        deformation_matrix = self._dfm * scale_factor
        return Strain(deformation_matrix)

    # return Green-Lagrange strain matrix
    @property
    def strain(self):
        return 0.5 * (np.matrix(self._dfm) * np.transpose(np.matrix(self._dfm)) - np.eye(3))

    @property
    def deformation_matrix(self):
        return self._dfm

#    # construct def. matrix from indices and amount
    @staticmethod
    def from_ind_amt_dfm(matrixpos, amt):
        F = np.identity(3)
        F[matrixpos] = F[matrixpos] + amt
        return Strain(F)

    def __eq__(self, other):
        df, df2 = self.deformation_matrix, other.deformation_matrix
        for i, row in enumerate(df):
            for j, item in enumerate(row):
                if df[i][j] != df2[i][j]:
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        #for now, just use a sum of elements * row * column^2
        df = self.deformation_matrix
        h_sum = 0.0
        for i, row in enumerate(df):
            for j, item in enumerate(row):
                h_sum += item * (i + 1) * (j + 1) * (j + 1)
        return h_sum


class IndependentStrain(Strain):
   #TODO: add polar decomposition method

    def __init__(self, deformation):
        
        super(IndependentStrain, self).__init__(deformation)
        (self._i, self._j) = self.check_F()

    def check_F(self):
        df1 = self.deformation_matrix
        sum1 = 0
        sum2 = 0
        for c1 in range(0,3):
            for c2 in range(0,3):
                if c1 != c2:
                    sum1 = sum1 + np.abs(df1[c1,c2])
                    if np.abs(df1[c1,c2]) > 0.00000001: 
                        sum2 = sum2 + 1

        if sum1<0.00000001: # if no shear components present
            if len(np.nonzero(df1-np.diag([1,1,1])>0.00000001)[0])>1: # check hom many diagonal components differ from unity
                raise ValueError("More than one normal mode was applied.")
            elif len(np.nonzero(df1-np.diag([1,1,1])>0.00000001)[0])==0: # if identity transformation
                return []
            else: # if proper transformation
                return np.nonzero(df1-np.diag([1,1,1])>0.00000001)[0][0], np.nonzero(df1-np.diag([1,1,1])>0.00000001)[1][0]
        
        else: # if shear components present
            if sum2 > 1: # if multiple shear components present
                raise ValueError("More than one shear mode was applied.")
            elif len(np.nonzero(np.abs(df1-np.diag([1,1,1]))>0.00000001)[0])>1: # if one shear def. present but also normal modes:
                raise ValueError("Shear and normal deformations were applied simultaneously.")
            else: # if proper transformation
                return (np.nonzero(np.abs(df1-np.diag([1,1,1]))>0.00000001)[0][0], np.nonzero(np.abs(df1-np.diag([1,1,1]))>0.00000001)[1][0])

    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j



  




if __name__ == "__main__":


    mat = np.eye(3)
    mat[2,1] = 1.05
#    mat[0,1] = 1.00001

    my_strain = Strain(mat)
#    print my_strain.deformation_matrix
#    print my_strain.strain

    my_strain2 = IndependentStrain(mat)
    print my_strain2._j


#    print my_strain2.check_F()

#    my_strain2.checkF

#    print my_strain.__dict__.keys()
#
#    print my_strain.deformation_matrix
#    print my_strain.strain
#    my_strain.index


#    my_scaled_strain = my_strain.get_scaled(1.05)
#    print my_scaled_strain.deformation_matrix
#    print my_scaled_strain.strain
#
#    print my_strain == my_scaled_strain
#
#    mat2 = np.eye(3)
#    mat2[0,0] = 1.01
#    my_strain2 = Strain(mat)
#    print my_strain == my_strain2



