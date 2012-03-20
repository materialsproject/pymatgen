import warnings, sys, operator, os, unittest
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/') # (If one does not want to change $PYTHONPATH)
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

    def sort_stress_strain(self, strain, stress):
        indices = [i for (i,j) in sorted(enumerate(strain), key=operator.itemgetter(1))]
        epsilon = sorted(strain)
        sigma = []

        for i in indices:
            sigma.append(-1.000*stress[i])
        return epsilon, sigma

    def chain_stresses(self, stress, ind1, ind2):
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

    def fitCij(self):
        
        stress_dict = self._sig_eps
        Cij = np.zeros((6, 6))
#        print Cij

        eps11 = []
        eps22 = []
        eps33 = []
        eps12 = []
        eps23 = []
        eps13 = []

        stress11 = []
        stress22 = []
        stress33 = []
        stress12 = []
        stress23 = []
        stress13 = []
        
        for c in stress_dict:

            if c.i==0 and c.j==0:
                eps11.append(c.strain[c.i, c.j])
                stress11.append(stress_dict[c])

            if c.i==1 and c.j==1:
                eps22.append(c.strain[c.i, c.j])
                stress22.append(stress_dict[c])

            if c.i==2 and c.j==2:
                eps33.append(c.strain[c.i, c.j])
                stress33.append(stress_dict[c])

            if c.i==0 and c.j==1:
                eps12.append(c.strain[c.i, c.j])
                stress12.append(stress_dict[c])

            if c.i==1 and c.j==2:
                eps23.append(c.strain[c.i, c.j])
                stress23.append(stress_dict[c])

            if c.i==0 and c.j==2:
                eps13.append(c.strain[c.i, c.j])
                stress13.append(stress_dict[c])
                
        [epsilon11, sigma11] = self.sort_stress_strain(eps11,stress11)
        [epsilon22, sigma22] = self.sort_stress_strain(eps22,stress22)
        [epsilon33, sigma33] = self.sort_stress_strain(eps33,stress33)
        [epsilon12, sigma12] = self.sort_stress_strain(eps12,stress12)
        [epsilon23, sigma23] = self.sort_stress_strain(eps23,stress23)
        [epsilon13, sigma13] = self.sort_stress_strain(eps13,stress13)

        inds = [[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]]

        count1 = 0
        for k in inds:

            true_data = self.chain_stresses(sigma11, k[0], k[1])

            p1 = np.polyfit(epsilon11, true_data, 1)
            p2 = np.polyfit(epsilon11, true_data, 2)
            p3 = np.polyfit(epsilon11, true_data, 3)

            f1 = np.polyval(p1, epsilon11)
            f2 = np.polyval(p2, epsilon11)
            f3 = np.polyval(p3, epsilon11)

            r2 = self.rsquared(true_data, f1)
            Cij[count1, 0] = 0.10*p1[0]
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = self.chain_stresses(sigma22, k[0], k[1])

            p1 = np.polyfit(epsilon22, true_data, 1)
            p2 = np.polyfit(epsilon22, true_data, 2)
            p3 = np.polyfit(epsilon22, true_data, 3)

            f1 = np.polyval(p1, epsilon22)
            f2 = np.polyval(p2, epsilon22)
            f3 = np.polyval(p3, epsilon22)

            r2 = self.rsquared(true_data, f1)
            Cij[count1, 1] = 0.10*p1[0]
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = self.chain_stresses(sigma33, k[0], k[1])

            p1 = np.polyfit(epsilon33, true_data, 1)
            p2 = np.polyfit(epsilon33, true_data, 2)
            p3 = np.polyfit(epsilon33, true_data, 3)

            f1 = np.polyval(p1, epsilon33)
            f2 = np.polyval(p2, epsilon33)
            f3 = np.polyval(p3, epsilon33)

            r2 = self.rsquared(true_data, f1)
            Cij[count1, 2] = 0.10*p1[0]
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = self.chain_stresses(sigma23, k[0], k[1])

            p1 = np.polyfit(epsilon23, true_data, 1)
            p2 = np.polyfit(epsilon23, true_data, 2)
            p3 = np.polyfit(epsilon23, true_data, 3)

            f1 = np.polyval(p1, epsilon23)
            f2 = np.polyval(p2, epsilon23)
            f3 = np.polyval(p3, epsilon23)

            r2 = self.rsquared(true_data, f1)
            Cij[count1, 3] = 0.10*p1[0]/2
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = self.chain_stresses(sigma13, k[0], k[1])

            p1 = np.polyfit(epsilon13, true_data, 1)
            p2 = np.polyfit(epsilon13, true_data, 2)
            p3 = np.polyfit(epsilon13, true_data, 3)

            f1 = np.polyval(p1, epsilon13)
            f2 = np.polyval(p2, epsilon13)
            f3 = np.polyval(p3, epsilon13)

            r2 = self.rsquared(true_data, f1)
            Cij[count1, 4] =  0.10*p1[0]/2
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = self.chain_stresses(sigma12, k[0], k[1])

            p1 = np.polyfit(epsilon12, true_data, 1)
            p2 = np.polyfit(epsilon12, true_data, 2)
            p3 = np.polyfit(epsilon12, true_data, 3)

            f1 = np.polyval(p1, epsilon12)
            f2 = np.polyval(p2, epsilon12)
            f3 = np.polyval(p3, epsilon12)

            r2 = self.rsquared(true_data, f1)
            Cij[count1, 5] =   0.10*p1[0]/2
            count1 += 1

        for n1 in range(0,6):
            for n2 in range(0,6):

                if np.abs(Cij[n1, n2]) < 1:

                    Cij[n1, n2] = 0

        
#        print Cij



