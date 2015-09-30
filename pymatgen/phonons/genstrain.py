from __future__ import division
import warnings
import sys
sys.path.append('')
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
from pymatgen.phonons.strain import Strain
from pymatgen.phonons.strain import IndependentStrain

__author__="Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="Jan 24, 2012"

def deform_geometry(rlxd_str, nd=0.01, ns=0.08, m=4, n=4):
    '''
    function to deform the geometry of a structure.  Generates
    m x n deformed structures according to the supplied parameters.

    Args:
        rlxd_str (structure): structure to undergo deformation
        nd (float): maximum perturbation applied to deformation
        ns (float): maximum perturbation applied to shear deformation
        m (int): number of deformation structures to generate for 
            normal deformation
        n (int): number of deformation structures to generate for 
            shear deformation
    '''

    msteps = np.int(m)
    nsteps = np.int(n)

    if m%2 != 0:
        raise ValueError("m has to be even.")
    if n%2 != 0:
        raise ValueError("n has to be even.")

    # TODO: JM asks is this line necessary?
    # TODO: JM suggests we might want to not delete the center point.
    # mystrains = np.zeros((3, 3, np.int(m*3) + np.int(n*3)))   
    defs = np.linspace(-nd, nd, num=m+1)
    defs = np.delete(defs, np.int(m/2), 0)
    sheardef = np.linspace(-ns, ns, num=n+1)
    sheardef = np.delete(sheardef, np.int(n/2), 0)
    defstructures = {}

    # TODO: JM asks is this line necessary?
    # strcopy = rlxd_str.copy()

    # TODO: Make this section more pythonic
    # First apply normal deformations
    for i1 in range(0, 3):
        for i2 in range(0, len(defs)):
            s=rlxd_str.copy()
            F = np.identity(3)
            F[i1, i1] = F[i1, i1] + defs[i2] 
            StrainObject = IndependentStrain(F)
            s.apply_deformation_gradient(F)
            defstructures[StrainObject] = s

    # Now apply shear deformations #		
    F_index = [[0, 1], [0, 2], [1, 2]]
    for j1 in range(0, 3):
        for j2 in range(0, len(sheardef)):
            s=rlxd_str.copy()
            F = np.identity(3)
            F[F_index[j1][0], F_index[j1][1]] = F[F_index[j1][0], F_index[j1][1]] + sheardef[j2]
#           F = np.matrix(F)   # this needs to be checked carefully, might give problems in certain cases
            StrainObject = IndependentStrain(F)
            s.apply_deformation_gradient(F)
            defstructures[StrainObject] = s

    return defstructures

if __name__ == "__main__":

    struct = CifParser('aluminum.cif').get_structures(primitive=False)[0]
#    print struct
    D = DeformGeometry(struct)
    print D
#    print D


