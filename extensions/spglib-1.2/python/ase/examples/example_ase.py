#!/usr/bin/evn python

import sys
from ase import *
from pyspglib import spglib
import numpy as np

silicon = Atoms( symbols=['Si']*8,
                 cell=[(4,0,0),
                       (0,4,0),
                       (0,0,4)],
                 scaled_positions=[(0, 0, 0),
                                   (0, 0.5, 0.5),
                                   (0.5, 0, 0.5),
                                   (0.5, 0.5, 0),
                                   (0.25, 0.25, 0.25),
                                   (0.25, 0.75, 0.75),
                                   (0.75, 0.25, 0.75),
                                   (0.75, 0.75, 0.25)],
                 pbc=True)

silicon_prim = Atoms( symbols=['Si']*2,
                      cell=[(0,2,2),
                            (2,0,2),
                            (2,2,0)],
                      scaled_positions=[(0, 0, 0),
                                        (0.25, 0.25, 0.25)],
                 pbc=True)

rutile = Atoms( symbols=['Si']*2+['O']*4,
                cell=[ (4,0,0),
                       (0,4,0),
                       (0,0,3) ],
                scaled_positions=[(0, 0, 0),
                                  (0.5, 0.5, 0.5),
                                  (0.3, 0.3, 0.0),
                                  (0.7, 0.7, 0.0),
                                  (0.2, 0.8, 0.5),
                                  (0.8, 0.2, 0.5)],
                pbc=True )

# For VASP case
# import ase.io.vasp as vasp
# bulk = vasp.read_vasp(sys.argv[1])

print "[get_spacegroup]"
print "  Spacegroup of Silicon is ", spglib.get_spacegroup(silicon)
print ""
print "[get_spacegroup]"
print "  Spacegroup of Rutile is ", spglib.get_spacegroup(rutile)
print ""
print "[get_symmetry]"
print "  Symmetry operations of Rutile unitcell are:"
print ""
symmetry = spglib.get_symmetry( rutile )
for i in range(symmetry['rotations'].shape[0]):
  print "  --------------- %4d ---------------" % (i+1)
  rot = symmetry['rotations'][i]
  trans = symmetry['translations'][i]
  print "  rotation:"
  for x in rot:
    print "     [%2d %2d %2d]" % (x[0], x[1], x[2])
  print "  translation:"
  print "     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2])
print ""

dataset = spglib.get_symmetry_dataset( rutile )
print "[get_symmetry_dataset] ['international']"
print "  Spacegroup of Rutile is ", dataset['international']
print ""
print "[get_symmetry_dataset] ['wyckoffs']"
alphabet = "abcdefghijklmnopqrstuvwxyz"
print "  Wyckoff letters of Rutile are: ", dataset['wyckoffs']
print ""
print "[get_symmetry_dataset] ['equivalent_atoms']"
print "  Mapping to equivalent atoms of Rutile are: "
for i, x in enumerate( dataset['equivalent_atoms'] ):
  print "  %d -> %d" % ( i+1, x+1 )
print ""
print "[get_symmetry_dataset] ['rotations'], ['translations']"
print "  Symmetry operations of Rutile unitcell are:"
for i, (rot,trans) in enumerate( zip( dataset['rotations'], dataset['translations'] ) ):
  print "  --------------- %4d ---------------" % (i+1)
  print "  rotation:"
  for x in rot:
    print "     [%2d %2d %2d]" % (x[0], x[1], x[2])
  print "  translation:"
  print "     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2])
print ""

symmetry = spglib.get_symmetry(silicon)
print "[get_symmetry]"
print "  Number of symmetry operations of silicon convensional"
print "  unit cell is ", len( symmetry['rotations'] ), "(192)."
print ""

mapping, grid = spglib.get_ir_reciprocal_mesh( [ 11, 11, 11 ],
                                               silicon_prim,
                                               is_shift=[ 0, 0, 0 ] )
num_ir_kpt = len( np.unique( mapping ) )
print "[get_ir_reciprocal_mesh]"
print "  Number of irreducible k-points of primitive silicon with"
print "  11x11x11 Monkhorst-Pack mesh is ", num_ir_kpt, "(56)."
print ""

mapping, grid = spglib.get_ir_reciprocal_mesh( [ 8, 8, 8 ],
                                               rutile,
                                               is_shift=[ 1, 1, 1 ] )
num_ir_kpt = len( np.unique( mapping ) )
print "[get_ir_reciprocal_mesh]"
print "  Number of irreducible k-points of Rutile with"
print "  8x8x8 Monkhorst-Pack mesh is ", num_ir_kpt, "(40)."
print ""


mesh = np.array([8,8,8])
kpoints = []
for i in range(mesh[0]):
  for j in range(mesh[1]):
    for k in range(mesh[2]):
      kpoints.append([float(i)/mesh[0],
                      float(j)/mesh[1],
                      float(k)/mesh[2]])
kpoints = np.array(kpoints) + 0.5/mesh
kpoints = kpoints - 1 * ( kpoints > 0.5 ) 
mapping = spglib.get_ir_kpoints( kpoints, rutile )
num_ir_kpt = len( np.unique( mapping ) )
print "[get_ir_kpoints]"
print "  Number of irreducible k-points of Rutile with"
print "  8x8x8 Monkhorst-Pack mesh is ", num_ir_kpt, "(40)."
print ""

