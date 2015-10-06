import warnings
import sys
import operator
import unittest
from StringIO import StringIO as sio
import pymatgen
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Vasprun
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
#from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
import random
#from scipy import stats
import math
from strain import generate_deformed_structures
from pymatgen.phonons.stress import Stress
from pymatgen.phonons.strain import Strain
from pymatgen.phonons.strain import IndependentStrain
from pymatgen.phonons.Cij import CijFitter
from pymatgen.phonons.tensors import SQTensor

np.set_printoptions(precision=3)

# example for rhenium, not aluminum!!

if __name__ == "__main__":

    f = open('stresses', 'r') 
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines]
    sig1 = np.matrix(np.loadtxt(sio('\n'.join(lines[0:3]))))

    stress_dict = {}

    struct = CifParser('aluminum.cif').get_structures(primitive=False)[0]
    
    D = generate_deformed_structures(struct, 0.01, 0.004, 4, 4)

#   stress_dict: {IndependentStrain object:Stress object}

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
            strain_array1.append(c[c.i, c.j])
            numdefnormal +=1

        elif c.i==1 and c.j==1:
            strain_array2.append(c[c.i, c.j])

        elif c.i==2 and c.j==2:
            strain_array3.append(c[c.i, c.j])

        elif c.i==1 and c.j==2:
            strain_array4.append(c[c.i, c.j])
            numdefshear +=1

        elif c.i==0 and c.j==2:
            strain_array5.append(c[c.i, c.j])

        elif c.i==0 and c.j==1:
            strain_array6.append(c[c.i, c.j])
        
    count = 0
    
    strain_array1 = sorted(strain_array1)
    strain_array2 = sorted(strain_array2)
    strain_array3 = sorted(strain_array3)
    strain_array4 = sorted(strain_array4)
    strain_array5 = sorted(strain_array5)
    strain_array6 = sorted(strain_array6)

    for n in strain_array1:
        for c in D:
            if c.i==0 and c.j==0 and np.abs(c[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = Stress(np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3])))))
                count+=1

    for n in strain_array2:
        for c in D:
            if c.i==1 and c.j==1 and np.abs(c[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = Stress(np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3])))))
                count+=1

    for n in strain_array3:
        for c in D:
            if c.i==2 and c.j==2 and np.abs(c[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = Stress(np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3])))))
                count+=1

    for n in strain_array4:
        for c in D:
            if c.i==1 and c.j==2 and np.abs(c[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = Stress(np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3])))))
                count+=1

    for n in strain_array5:
        for c in D:
            if c.i==0 and c.j==2 and np.abs(c[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = Stress(np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3])))))
                count+=1

    for n in strain_array6:
        for c in D:
            if c.i==0 and c.j==1 and np.abs(c[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = Stress(np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3])))))
                count+=1




#print stress_dict

#for c in D:

#    print c
#    print dir(stress_dict[c])
#    print stress_dict[c].MeanStress
#    print stress_dict[c].stress_matrix

Cij = CijFitter(stress_dict).fitCij()

print Cij

Cij2 = CijFitter(stress_dict).fitCij2()

print Cij2

Cij = 0.5*(Cij + np.transpose(Cij))

SQTensor(Cij).KG_average




