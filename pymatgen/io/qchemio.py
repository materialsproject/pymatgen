#!/usr/bin/env python

"""
This module implements input and output processing from QChem.
"""
from pymatgen.serializers.json_coders import MSONable

__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2013, The Electrolyte Genome Project"
__version__ = "0.1"
__maintainer__ = "Xiaohui Qu"
__email__ = "xhqu1981@gmail.com"
__date__ = "11/4/13"


class QcInput(MSONable):
    '''
        An object representing a QChem input file.
    '''
    def __init__(self, molecule=None, change=None, spin_multiplicity=None,
                 title=None, exchange="B3LYP", correlation=None,
                 basis_set="6-31+G*", rem_params=None, optional_params=None):
        '''
        Args:
            molecule:
        '''
        self.mol = molecule