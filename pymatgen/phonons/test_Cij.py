import warnings
import sys
import operator
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/') # (If one does not want to change $PYTHONPATH)
import unittest
from StringIO import StringIO as sio
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
import random
from scipy import stats
import math

from pymatgen.phonons.strain import Strain
from pymatgen.phonons.strain import IndependentStrain
from pymatgen.phonons.genstrain import DeformGeometry
from pymatgen.phonons.Cij import CijTensor

np.set_printoptions(precision=3)

if __name__ == "__main__":

    f = open('stresses', 'r') 
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines]
    sig1 = np.matrix(np.loadtxt(sio('\n'.join(lines[0:3]))))

    stress_dict = {}

    struct = CifParser('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen/phonons/aluminum.cif').get_structures()[0]
    D = DeformGeometry(struct)
    D = D.deform(0.01, 0.004, 4, 4)

#   stress_dict: {IndependentStrain:stress matrix}

    strain_array1 = []
    strain_array2 = []
    strain_array3 = []
    strain_array4 = []
    strain_array5 = []
    strain_array6 = []

    count = 0
    
    numdefnormal = 0
    numdefshear = 0

    for c in D:
        if c.i==0 and c.j==0:
            strain_array1.append(c.strain[c.i, c.j])
            numdefnormal +=1

        elif c.i==1 and c.j==1:
            strain_array2.append(c.strain[c.i, c.j])

        elif c.i==2 and c.j==2:
            strain_array3.append(c.strain[c.i, c.j])

        elif c.i==1 and c.j==2:
            strain_array4.append(c.strain[c.i, c.j])
            numdefshear +=1

        elif c.i==0 and c.j==2:
            strain_array5.append(c.strain[c.i, c.j])

        elif c.i==0 and c.j==1:
            strain_array6.append(c.strain[c.i, c.j])
        
    count = 0
    
    strain_array1 = sorted(strain_array1)
    strain_array2 = sorted(strain_array2)
    strain_array3 = sorted(strain_array3)
    strain_array4 = sorted(strain_array4)
    strain_array5 = sorted(strain_array5)
    strain_array6 = sorted(strain_array6)

#   fill up 

    for n in strain_array1:
        for c in D:
            if c.i==0 and c.j==0 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1 

    for n in strain_array2:
        for c in D:
            if c.i==1 and c.j==1 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1

    for n in strain_array3:
        for c in D:
            if c.i==2 and c.j==2 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1

    for n in strain_array4:
        for c in D:
            if c.i==1 and c.j==2 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1

    for n in strain_array5:
        for c in D:
            if c.i==0 and c.j==2 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1

    for n in strain_array6:
        for c in D:
            if c.i==0 and c.j==1 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1


#    CijTensor(stress_dict).fitCij()
    Cij = CijTensor(stress_dict).fitCij()
    print Cij


#    print stress_dict

#   at this point, the stress_dict is created in such a way that fit_elas can operate on it, i.e. it is in the format stress_dict: {IndependentStrain:stress matrix}

    

    """
    def sort_stress_strain(strain, stress):
        indices = [i for (i,j) in sorted(enumerate(strain), key=operator.itemgetter(1))]
        epsilon = sorted(strain)
        sigma = []

        for i in indices:
            sigma.append(-1.000*stress[i])
        return epsilon, sigma

    def chain_stresses(stress, ind1, ind2):
        s = []
        for el in stress:
            s.append(el[ind1,ind2])
        return s

    def rsquared(yi, fi):

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
#        print math.isnan(1-sum1/sum2)

        return r2

    def fitCij(stress_dict):

        Cij = np.zeros((6, 6))

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
                

        [epsilon11, sigma11] = sort_stress_strain(eps11,stress11)
        [epsilon22, sigma22] = sort_stress_strain(eps22,stress22)
        [epsilon33, sigma33] = sort_stress_strain(eps33,stress33)
        [epsilon12, sigma12] = sort_stress_strain(eps12,stress12)
        [epsilon23, sigma23] = sort_stress_strain(eps23,stress23)
        [epsilon13, sigma13] = sort_stress_strain(eps13,stress13)

        inds = [[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]]

        count1 = 0
        for k in inds:

            true_data = chain_stresses(sigma11, k[0], k[1])

            p1 = np.polyfit(epsilon11, true_data, 1)
            p2 = np.polyfit(epsilon11, true_data, 2)
            p3 = np.polyfit(epsilon11, true_data, 3)

            f1 = np.polyval(p1, epsilon11)
            f2 = np.polyval(p2, epsilon11)
            f3 = np.polyval(p3, epsilon11)

            r2 = rsquared(true_data, f1)
            Cij[count1, 0] = 0.10*p1[0]
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = chain_stresses(sigma22, k[0], k[1])

            p1 = np.polyfit(epsilon22, true_data, 1)
            p2 = np.polyfit(epsilon22, true_data, 2)
            p3 = np.polyfit(epsilon22, true_data, 3)

            f1 = np.polyval(p1, epsilon22)
            f2 = np.polyval(p2, epsilon22)
            f3 = np.polyval(p3, epsilon22)

            r2 = rsquared(true_data, f1)
            Cij[count1, 1] = 0.10*p1[0]
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = chain_stresses(sigma33, k[0], k[1])

            p1 = np.polyfit(epsilon33, true_data, 1)
            p2 = np.polyfit(epsilon33, true_data, 2)
            p3 = np.polyfit(epsilon33, true_data, 3)

            f1 = np.polyval(p1, epsilon33)
            f2 = np.polyval(p2, epsilon33)
            f3 = np.polyval(p3, epsilon33)

            r2 = rsquared(true_data, f1)
            Cij[count1, 2] = 0.10*p1[0]
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = chain_stresses(sigma23, k[0], k[1])

            p1 = np.polyfit(epsilon23, true_data, 1)
            p2 = np.polyfit(epsilon23, true_data, 2)
            p3 = np.polyfit(epsilon23, true_data, 3)

            f1 = np.polyval(p1, epsilon23)
            f2 = np.polyval(p2, epsilon23)
            f3 = np.polyval(p3, epsilon23)

            r2 = rsquared(true_data, f1)
            Cij[count1, 3] = 0.10*p1[0]/2
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = chain_stresses(sigma13, k[0], k[1])

            p1 = np.polyfit(epsilon13, true_data, 1)
            p2 = np.polyfit(epsilon13, true_data, 2)
            p3 = np.polyfit(epsilon13, true_data, 3)

            f1 = np.polyval(p1, epsilon13)
            f2 = np.polyval(p2, epsilon13)
            f3 = np.polyval(p3, epsilon13)

            r2 = rsquared(true_data, f1)
            Cij[count1, 4] =  0.10*p1[0]/2
            count1 += 1

        count1 = 0

        for k in inds:

            true_data = chain_stresses(sigma12, k[0], k[1])

            p1 = np.polyfit(epsilon12, true_data, 1)
            p2 = np.polyfit(epsilon12, true_data, 2)
            p3 = np.polyfit(epsilon12, true_data, 3)

            f1 = np.polyval(p1, epsilon12)
            f2 = np.polyval(p2, epsilon12)
            f3 = np.polyval(p3, epsilon12)

            r2 = rsquared(true_data, f1)
            Cij[count1, 5] =   0.10*p1[0]/2
            count1 += 1

        for n1 in range(0,6):
            for n2 in range(0,6):

                if np.abs(Cij[n1, n2]) < 1:

                    Cij[n1, n2] = 0

        
        print Cij

    fitCij(stress_dict)

    """
    

