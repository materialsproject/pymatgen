#!/usr/bin/env python

'''
This module implements abstract base classes for file io classes.  For seamless integration
with the rest of the code base, any classes providing to a file io function should extend this class.
'''

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import abc

class VaspInput(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod 
    def write_file(self, filename):
        '''
        Writes the contents to a file.
        '''
        return
    
