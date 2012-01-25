#!/usr/bin/env python

'''
This module provides input and output from the CSSR file format.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jan 24, 2012"


class Cssr(object):
    """
    Basic object for working with Cssr file. Right now, only conversion from 
    a Structure to a Cssr file is supported.
    """
    
    def __init__(self, structure):
        if not structure.is_ordered:
            raise ValueError("Cssr file can only be constructed from ordered structure")
        self._structure = structure
        
    def __str__(self):
        output = []
        output.append("{:.4f} {:.4f} {:.4f}".format(*self._structure.lattice.abc))
        output.append("{:.2f} {:.2f} {:.2f} SPGR =  1 P 1    OPT = 1".format(*self._structure.lattice.angles))
        output.append("{} 0".format(len(self._structure)))
        output.append("0 {}".format(self._structure.formula))
        for i, site in enumerate(self._structure.sites):
            output.append("{} {} {:.4f} {:.4f} {:.4f}".format(i+1, site.specie, site.a, site.b, site.c))
        return "\n".join(output)
    
    def write_file(self, filename):
        with open(filename, 'w') as f:
            f.write(str(self) + "\n")
        
