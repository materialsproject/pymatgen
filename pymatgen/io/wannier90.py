# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Modules for working with wannier90 input and output.
"""

import numpy as np
from scipy.io import FortranFile, FortranEOFError

__author__ = "Mark Turiansky"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "4.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Jun 04, 2020"


class Unk:
    """
    Object representing the data in a UNK file.

    .. attribute:: 

        blah

    """

    def __init__(self, ik, data):
        """
        """
        self.ik = ik
        self.data = data

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        temp_val = np.array(value, dtype=np.complex128)
        if len(temp_val.shape) not in [4, 5]:
            raise ValueError('invalid data shape, must be (nbnd, ngx, ngy, ngz'
                             ') or (nbnd, 2, ngx, ngy, ngz) for noncollinear '
                             f'data, given {temp_val.shape}')
        if len(temp_val.shape) == 5 and temp_val.shape[1] != 2:
            raise ValueError('invalid noncollinear data, shape should be (nbnd'
                             f', 2, ngx, ngy, ngz), given {temp_val.shape}')
        self._data = temp_val
        self.is_noncollinear = (len(self.data.shape) == 5)
        self.ng = self.data.shape[-3:]

    @staticmethod
    def from_file(filename):
        """
        """
        data = []
        with FortranFile(filename, 'r') as f:
            *ng, ik, nbnd = f.read_ints()
            for _ in range(nbnd):
                data.append(f.read_record(np.complex128))
            try:
                for _ in range(nbnd):
                    data.append(f.read_record(np.complex128))
                is_noncollinear = True
            except FortranEOFError:
                is_noncollinear = False

        data = np.array(data, dtype=np.complex128)

        if is_noncollinear:
            data.shape = (nbnd, 2, *ng)
        else:
            data.shape = (nbnd, *ng)

        return Unk(ik, data)

    def write_file(self, filename):
        """
        """
        pass
