# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
#import re
#from monty.io import zopen
#from monty.re import regrep
from monty.json import MSONable

"""
Classes for reading/manipulating/writing QChem ouput files.
"""

__author__ = "Brandon Wood, Samuel Blau, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__email__ = "b.wood@berkeley.edu"

logger = logging.getLogger(__name__)


class QcJob(MSONable):
    """
    An object representing a qchem job
    """

    def __init__(self, molecule, rem, opt=None):
        self.molecule = molecule
        self.rem = rem
        self.opt = opt

    def __str__(self):
        combined_list = []
        # molecule section
        combined_list.append(molecule_template(self.molecule))
        combined_list.append("")
        # rem section
        combined_list.append(rem_template(self.rem))
        combined_list.append("")
        # opt section
        if self.opt:
            combined_list.append(opt_template(self.opt))
            combined_list.append("")
        return '\n'.join(combined_list)

    def write_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.__str__())
