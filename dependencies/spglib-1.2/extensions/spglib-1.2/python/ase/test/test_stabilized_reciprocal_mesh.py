#!/usr/bin/evn python

import sys
from optparse import OptionParser
from atoms import Atoms
from vasp import read_vasp
from pyspglib import spglib
import numpy as np

parser = OptionParser()
parser.set_defaults( mesh = None,
                     qpoints = None )
parser.add_option("-m", "--mesh", dest="mesh",
                  type="string",help="Mesh numbers")
parser.add_option("-q", "--qpoints", dest="qpoints",
                  type="string",help="Stabilizers")
(options, args) = parser.parse_args()

if options.mesh == None:
  mesh = [ 4, 4, 4 ]
else:
  mesh = [ int( x ) for x in options.mesh.split() ]

if options.qpoints == None:
  qpoints = np.array( [[ 0, 0, 0 ]], dtype=float )
else:
  qpoints = np.array( [ float( x ) for x in options.qpoints.split() ] ).reshape( -1, 3 )

cell = read_vasp( args[0] )

mapping, grid_point = \
    spglib.get_ir_reciprocal_mesh( mesh, cell )

print mapping

dataset = spglib.get_symmetry_dataset( cell )
mapping, grid_point = \
    spglib.get_stabilized_reciprocal_mesh( mesh,
                                           cell.get_cell(),
                                           dataset['rotations'],
                                           qpoints=qpoints )


print mapping
