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
__version__ = "0.1"
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
        """
        """
        temp_val = np.array(value, dtype=np.complex128)
        if len(temp_val.shape) not in [4, 5]:
            raise ValueError('invalid data shape, must be (nbnd, ngx, ngy, ngz'
                             ') or (nbnd, 2, ngx, ngy, ngz) for noncollinear '
                             f'data, given {temp_val.shape}')
        if len(temp_val.shape) == 5 and temp_val.shape[1] != 2:
            raise ValueError('invalid noncollinear data, shape should be (nbnd'
                             f', 2, ngx, ngy, ngz), given {temp_val.shape}')
        self._data = temp_val

        # derived properties
        self.is_noncollinear = (len(self.data.shape) == 5)
        self.nbnd = self.data.shape[0]
        self.ng = self.data.shape[-3:]

    @staticmethod
    def from_file(filename):
        """
        """
        data = []
        with FortranFile(filename, 'r') as f:
            *ng, ik, nbnd = f.read_ints()
            for _ in range(nbnd):
                data.append(
                    # when reshaping need to specify ordering as fortran
                    f.read_record(np.complex128).reshape(ng, order='F')
                )
            try:
                for _ in range(nbnd):
                    data.append(
                        f.read_record(np.complex128).reshape(ng, order='F')
                    )
                is_noncollinear = True
            except FortranEOFError:
                is_noncollinear = False

        data = np.array(data, dtype=np.complex128)

        # spinors are interwoven, need to separate them
        if is_noncollinear:
            tdata = np.empty((nbnd, 2, *ng))
            tdata[:, 0, :, :, :] = data[::2, :, :, :]
            tdata[:, 1, :, :, :] = data[1::2, :, :, :]
            return Unk(ik, tdata)
        return Unk(ik, data)

    def write_file(self, filename):
        """
        """
        with FortranFile(filename, 'w') as f:
            f.write_record(
                np.array([*self.ng, self.ik, self.nbnd], dtype=np.int32)
            )
            for ib in range(self.nbnd):
                if self.is_noncollinear:
                    f.write_record(self.data[ib, 0].flatten('F'))
                    f.write_record(self.data[ib, 1].flatten('F'))
                else:
                    f.write_record(self.data[ib].flatten('F'))

    def __repr__(self):
        return f'<UNK ik={self.ik} nbnd={self.nbnd} ncl={self.is_noncollinear}' \
            + f' ngx={self.ng[0]} ngy={self.ng[1]} ngz={self.ng[2]}>'

    def __eq__(self, other):
        if not np.allclose(self.ng, other.ng):
            return False

        if self.ik != other.ik:
            return False

        if self.is_noncollinear != other.is_noncollinear:
            return False

        if self.nbnd != other.nbnd:
            return False

        for ib in range(self.nbnd):
            if self.is_noncollinear:
                if not (np.allclose(self.data[ib, 0], other.data[ib, 0])
                        and np.allclose(self.data[ib, 1], other.data[ib, 1])):
                    return False
            else:
                if not np.allclose(self.data[ib], other.data[ib]):
                    return False
        return True
