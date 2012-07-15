#!/usr/bin/evn python

import sys
from optparse import OptionParser
from atoms import Atoms
from vasp import read_vasp
from pyspglib import spglib
import numpy as np

parser = OptionParser()
parser.set_defaults( mesh = None )
parser.add_option("-m", "--mesh", dest="mesh",
                  type="string",help="Mesh numbers")
(options, args) = parser.parse_args()

if options.mesh == None:
  mesh = [ 4, 4, 4 ]
else:
  mesh = [ int( x ) for x in options.mesh.split() ]
cell = read_vasp( args[0] )

# cell = Atoms( symbols=['Si']*2,
#               cell=[(0,2,2),
#                     (2,0,2),
#                     (2,2,0)],
#               scaled_positions=[(0, 0, 0),
#                                 (0.25, 0.25, 0.25)],
#                  pbc=True)
# mesh = [ 4, 4, 4 ]

dataset = spglib.get_symmetry_dataset( cell )
triplets, weights, grid_point = \
    spglib.get_triplets_reciprocal_mesh( mesh,
                                         cell.get_cell(),
                                         dataset['rotations'] )

for i, ( t, w ) in enumerate( zip( triplets, weights ) ):
  print i+1, t, w
print np.array(weights).sum()
