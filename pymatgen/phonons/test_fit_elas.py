import warnings
import sys
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
from pymatgen.phonons.strain import Strain
from pymatgen.phonons.strain import IndependentStrain
from pymatgen.phonons.genstrain import DeformGeometry

if __name__ == "__main__":

    f = open('stresses', 'r') 
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines]
    sig1 = np.matrix(np.loadtxt(sio('\n'.join(lines[0:3]))))

    stress_dict = {}

    struct = CifParser('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen/phonons/aluminum.cif').get_structures()[0]
    D = DeformGeometry(struct)
    D = D.deform()

#   stress_dict: {IndependentStrain:stress matrix}

    strain_array = []
    count = 0
    
    for c in D:
#        stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
        if c.i==0 and c.j==0:
            strain_array.append(c.strain[c.i, c.j])

    count = 0
    strain_array = sorted(strain_array)
    
    for n in strain_array:
        for c in D:
            if c.i==0 and c.j==0 and np.abs(c.strain[c.i, c.j] - n) < 0.00001:
                stress_dict[c] = np.matrix(np.loadtxt(sio('\n'.join(lines[3*count:3*count+3]))))
                count+=1 


    


#    print stress_dict
             

#    print strain_array


    for mode in stress_dict:

#        print stress_dict[mode]
        print mode._dfm


    
#       if mode.i==0 and mode.j==0:

#           print mode.strain

#        print mode.i, mode.j, mode.strain[mode.i, mode.j]       
#        print i
#        print stress_dict[i]
    
 
