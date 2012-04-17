from __future__ import division
import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen')
#sys.path.append('/home/MDEJONG1/local/lib/python2.7/site-packages/cctbx_sources')
#sys.path.append('/home/MDEJONG1/local/lib/python2.7/site-packages/boost_1_48_0/')
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
import os
from collections import defaultdict
import pymongo
#from cctbx import xray
#import cctbx
#import boost


#myCIF = CifParser('/home/MDEJONG1/pythonplayground/pymatgen/classes/7048.cif').get_structures()[0]
#w = Poscar(myCIF)
#w.write_file('file_test')


#s = StructureEditor(myCIF)

#newlat = np.zeros((3,3))
#newlat[0,0] = 3.141592
#newlat[1,1] = 5.670515
#newlat[2,2] = 6.141512

#print newlat

#newlat = Lattice(newlat)
#s.modify_lattice(newlat)

#print s.modified_structure

#print myCIF._lattice._matrix


trans = np.identity(3)
trans[1,2] = 0.1

#print type(myCIF._lattice._matrix)
#print type(np.matrix(myCIF._lattice._matrix))

#print np.matrix(myCIF._lattice._matrix)*trans



#s.apply_strain_transformation(trans)
#struct = s.modified_structure

#print struct
#print struct.__dict__.keys()
#print type(struct._sites)
#print len(struct._sites)
#print struct._lattice.__dict__.keys()
#print struct._lattice._matrix
#print newlat.coords
#print struct._sites[2].__dict__.keys()
#print struct._sites
#print type(struct._sites[2]._fcoords)
#print struct._sites[2]._coords

# test dictionaries

#d1 = dict([('sape', 4139), ('guido', 4127), ('jack', 4098)])
#print d1

#d2 = dict([(x, s.modified_structure) for x in (2, 4, 6)])
#print d2[2]

#l=[ [1, 'A', 'D'], [1, 'B'], [2, 'C'] ]
#print l[0][2]

#d1 = dict([  (1, [np.identity(3)])   ,   (2, [np.identity(3)])   ])
#d2 = ([np.identity(3)*4])

#d1[6] = d2

#print d1

#A = np.identity(3)*4

#for i in A:
#	print i

#print A
#if (A>3).any() == True:
#	warnings.warn('This is a warning message')

#print A
#string = "How are you?" 
#print "\033[34m" + string + "\033[0m"

#arr1 = []
#arr1.append([1,2,3])
#arr1.append([6,7,8])

#print arr1[0]

A=np.array([1, 2, 3, 4, 5])
B=np.array([2, 4, 6, 7, 11])

#print A
#print B

#print A*A

#print np.sum(A*A)






