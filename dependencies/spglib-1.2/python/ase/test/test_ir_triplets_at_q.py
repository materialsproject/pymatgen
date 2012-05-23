#!/usr/bin/evn python

import sys
from optparse import OptionParser
from atoms import Atoms
from vasp import read_vasp
from pyspglib import spglib
import numpy as np

def parse_triplets( filename ):
    file = open(filename, 'r')
    triplets = []
    weights = []
    for line in file:
        if line.strip()[0] == "#":
            continue

        line_array = [ int(x) for x in line.split() ]
        triplets.append( line_array[:3] )
        weights.append( line_array[3] )

    return np.array( triplets ), np.array( weights )


parser = OptionParser()
parser.set_defaults( mesh = None,
                     grid_point = 0,
                     filename = None )
parser.add_option("-m", "--mesh", dest="mesh",
                  type="string",help="Mesh numbers")
parser.add_option("-f", dest="filename",
                  type="string",help="Filename of triplets at q")
parser.add_option("-g", dest="grid_point",
                  type="int",help="A grid point")
(options, args) = parser.parse_args()

if options.mesh == None:
  mesh = [ 4, 4, 4 ]
else:
  mesh = [ int( x ) for x in options.mesh.split() ]

cell = read_vasp( args[0] )

dataset = spglib.get_symmetry_dataset( cell )

weights, third_q, grid_points = \
    spglib.get_triplets_reciprocal_mesh_at_q( options.grid_point,
                                              mesh,
                                              cell.get_cell(),
                                              dataset['rotations'] )
if options.filename==None:
  for i, ( w, q ) in enumerate( zip ( weights, third_q ) ):
    if w > 0:
      print options.grid_point, i, q, w
else:
  triplets_in, weights_in = parse_triplets( options.filename )
  count = 0
  for i, ( w, q ) in enumerate( zip ( weights, third_q ) ):
    if w > 0:
      if triplets_in[count][0] == options.grid_point and \
            triplets_in[count][1] == i and \
            triplets_in[count][2] == q and \
            weights_in[count] == w:
        print options.grid_point, i, q, w, "OK"
      else:
        print "Don't match", options.grid_point, i, q, w
        sys.exit(0)
      count += 1

print weights.sum()
