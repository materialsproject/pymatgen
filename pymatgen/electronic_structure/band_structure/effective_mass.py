#!/usr/bin/env python

"""
This module provides classes and utilities to analyze effective masses in a band structure
"""

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "March 14, 2012"

class EffectiveMassAnalyzer():
    """
    this class is used to get effective masses data from a band structure
    """
    def __init(self, bs):
        self._bs = bs

    def get_effective_mass_tensor_symmetry(self, target):
        """
        this method analyze the symmetry of the target kpoint and gives what form should have the effective mass tensor
        (1 distinct eigenvalue, 2 or 3) as well as a few principal directions.
        """
        listUc = self._bs.get_pg_matrices_rec()
        list_return = []
        list_op = []
        for u in listUc:
            print u
            newkpoint = np.dot(u['matrix'], kpoint.T)
            if(np.linalg.norm(newkpoint - kpoint) < 0.001):
                list_op.append(u)
        fourfold = False
        threefold = False
        rot_higher_2 = False
        for u in list_op:
            if(u['Schoenflies'][-1] == 4):
                fourfold = True
            if(u['Schoenflies'][-1] == 3):
                threefold = True

        if(threefold == True and fourfold == True):
            return 1
        for u in list_op:
            if(u['Schoenflies'][-1] == 4 or u['Schoenflies'][-1] == 3 or u['Schoenflies'][-1] == 6):
                #2 eigenvalues, one princ axis along the rot axis another is perp.
                pass
