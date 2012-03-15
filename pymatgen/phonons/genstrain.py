from __future__ import division
import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/') # (If one does not want to change $PYTHONPATH)
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
import os
import fit_elas
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

class DeformGeometry(object):

    """
    Class that takes a pymatgen core structure object and applies a range of deformation with the aim at determining the elastic tensor.

    Args:
        - rlxd_str: pymatgen core structure object (fully relaxed)
    """

    def __init__(self, rlxd_str):
        self._rlxd_str = rlxd_str

    def deform(self, nd=0.02, ns=0.02, m=4, n=4): 
        """
        This method applies deformations and stores the result in a dict with StrainObjects as keys

        Args:
            - nd, ns: maximum amount of normal and shear deformation, respectively
            - m,n: number of normal and shear deformation, respectively
        """	
        
        self.msteps = np.int(m)
        self.nsteps = np.int(n)

        if m%2!=0:
            raise ValueError("m has to be even.")
         
        if n%2!=0:
            raise ValueError("n has to be even.")
        
        mystrains = np.zeros((3, 3, np.int(m*3) + np.int(n*3)))
        
        defs = np.linspace(-nd, nd, num=m+1)
        defs = np.delete(defs, np.int(m/2), 0)
        sheardef = np.linspace(-ns, ns, num=n+1)
        sheardef = np.delete(sheardef, np.int(n/2), 0)
        defstructures = {}

        # First apply normal deformations
        for i1 in range(0, 3):

            for i2 in range(0, len(defs)):

                s = StructureEditor(self._rlxd_str)
                F = np.identity(3)
                F[i1, i1] = F[i1, i1] + defs[i2] 
                StrainObject = IndependentStrain(F)
                s.apply_strain_transformation(F)      
                defstructures[StrainObject] = s.modified_structure

        # Now apply shear deformations #		
        F_index = [[0, 1], [0, 2], [1, 2]]
        for j1 in range(0, 3):
		
            for j2 in range(0, len(sheardef)):

                s = StructureEditor(self._rlxd_str)
                F = np.identity(3)
                F[F_index[j1][0], F_index[j1][1]] = F[F_index[j1][0], F_index[j1][1]] + sheardef[j2]
#                F = np.matrix(F)   # this needs to be checked carefully
                StrainObject = IndependentStrain(F)
                s.apply_strain_transformation(F)   
                defstructures[StrainObject] = s.modified_structure          

        self.defstructures = defstructures
        
        return self.defstructures


if __name__ == "__main__":

    struct = CifParser('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen/phonons/aluminum.cif').get_structures()[0]
    D = DeformGeometry(struct)    
    D = D.deform()
    mynum = np.float64(35.91)
#    print type(D)
#    print D[mynum]
   


