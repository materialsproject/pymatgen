import warnings, sys, operator, os, unittest
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

__author__="Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Anubhav Jain, Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="March 20, 2012"


class CijTensor(object):
    
    """
    Class that takes a dict in the form {IndependentStrain:stress matrix} and fits the elastic tensor from this.

    Args:
        - strain_stress_dict: dict containing stress matrices
    """

    def __init__(self, strain_stress_dict):
        self._sig_eps = strain_stress_dict

    def _chain_stresses(self, stress, ind1, ind2):
        s = []
        for el in stress:
            s.append(el[ind1,ind2])
        return s

    def rsquared(self, yi, fi):
        ybar = 1.0/len(yi)*sum(yi)
        v1 = yi-fi
        v2 = yi-ybar
        
        sum1 = 0
        sum2 = 0

        for n in v1:
            sum1 = sum1 + n**2

        for n in v2:
            sum2 = sum2 + n**2

        r2 = 1.0-sum1/sum2
        return r2

    def fitCij(self, tol=0.1):
        stress_dict = self._sig_eps
        Cij = np.zeros((6, 6))

        inds = [[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]]
        
        for n1 in range(0, 6):
            strain = []
            stress = []
            count1 = 0

            for c in stress_dict:
                if c.i == inds[n1][0] and c.j== inds[n1][1]:
                    strain.append(c.strain[c.i, c.j])
                    stress.append(stress_dict[c].stress_matrix)

            for k in inds:
                true_data = self._chain_stresses(stress, k[0], k[1])
                p1 = np.polyfit(strain, true_data, 1)
                p2 = np.polyfit(strain, true_data, 2)
                p3 = np.polyfit(strain, true_data, 3)

                f1 = np.polyval(p1, strain)
                f2 = np.polyval(p2, strain)
                f3 = np.polyval(p3, strain)
                
                r2 = self.rsquared(true_data, f1)
                Cij[count1, n1] = -0.10*p1[0]
                count1 += 1                

        for n1 in range(0,6):
            for n2 in range(0,6):
                if np.abs(Cij[n1, n2]) < tol:
                    Cij[n1, n2] = 0

                if n2 > 2:
                    Cij[n1, n2] = Cij[n1, n2]*0.50

        
        return Cij

