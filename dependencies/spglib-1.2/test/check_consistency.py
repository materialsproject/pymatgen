#!/usr/bin/python
from numpy import *
import sys
import datetime
from phonopy.structure.spglib import *
from phonopy.structure.atoms import Atoms
from phonopy.interface.vasp import read_vasp

# lattice = array([[4.,0.,0.], [0.,4.,0.], [0.,0.,6.]])
# position = array([
# 	[0.0, 0.0, 0.0],
# 	[0.5, 0.5, 0.25],
# 	[0.3, 0.3, 0.0],
# 	[0.7, 0.7, 0.0],
# 	[0.2, 0.8, 0.25],
# 	[0.8, 0.2, 0.25],
# 	[0.0, 0.0, 0.5],
# 	[0.5, 0.5, 0.75],
# 	[0.3, 0.3, 0.5],
# 	[0.7, 0.7, 0.5],
# 	[0.2, 0.8, 0.75],
# 	[0.8, 0.2, 0.75]
# 	])
# atom_type = array([1,1,2,2,2,2,1,1,2,2,2,2])
# atoms = Atoms( cell=lattice, scaled_positions=position, numbers=atom_type )
atoms = read_vasp( sys.argv[1] )
# print datetime.datetime.now()
# print "get_symmetry_dataset"
dataset = get_symmetry_dataset(atoms, 1e-5)
# for i, (r,t) in enumerate( zip( dataset['rotations'], dataset['translations'] ) ):
#   print "--- %d ---" % (i+1)
#   print r
#   print t

# print
# print dataset['transformation_matrix']
# print dataset['origin_shift']
# print dataset['international'], dataset['number']
# print dataset['hall']

# print datetime.datetime.now()
# print "get_symmetry"
symmetry = get_symmetry(atoms, 1e-5)
# print datetime.datetime.now()

count = 0
# print "Check consistency"
for r1, t1 in zip( dataset['rotations'],
                   dataset['translations'] ):

  for r2, t2 in zip( symmetry['rotations'],
                     symmetry['translations'] ):

    
    if ( (r1 - r2) == 0 ).all():
      diff = t1 - t2
      diff -= diff.round()
      if ( abs(diff) < 1e-5 ).all():
        count += 1
        # print count
        # print r1
        # print t1
        break

if ( count==len(dataset['rotations']) ):
  print "OK"
else:
  print "BAD"
