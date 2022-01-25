# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Modules for working with wannier90 input and output.
"""

from typing import Sequence

import numpy as np
from scipy.io import FortranEOFError, FortranFile

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

    .. attribute:: ik

        int index of kpoint for this file

    .. attribute:: data

        numpy.ndarray that contains the wavefunction data for in the UNK file.
        The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and
        (nbnd, 2, ngx, ngy, ngz) for noncollinear calculations.

    .. attribute:: is_noncollinear

        bool that specifies if data is from a noncollinear calculation

    .. attribute:: nbnd

        int number of bands in data

    .. attribute:: ng

        sequence of three integers that correspond to the grid size of the
        given data. The definition is ng = (ngx, ngy, ngz).

    """

    ik: int
    is_noncollinear: bool
    nbnd: int
    ng: Sequence[int]

    def __init__(self, ik: int, data: np.ndarray) -> None:
        """
        Initialize Unk class.

        Args:
            ik (int): index of the kpoint UNK file is for
            data (np.ndarray): data from the UNK file that has shape (nbnd,
                ngx, ngy, ngz) or (nbnd, 2, ngx, ngy, ngz) if noncollinear
        """
        self.ik = ik
        self.data = data

    @property
    def data(self) -> np.ndarray:
        """
        np.ndarray: contains the wavefunction data for in the UNK file.
        The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and
        (nbnd, 2, ngx, ngy, ngz) for noncollinear calculations.
        """
        return self._data

    @data.setter
    def data(self, value: np.ndarray) -> None:
        """
        Sets the value of data.

        Args:
            value (np.ndarray): data to replace stored data, must haveshape
                (nbnd, ngx, ngy, ngz) or (nbnd, 2, ngx, ngy, ngz) if
                noncollinear calculation
        """
        temp_val = np.array(value, dtype=np.complex128)
        if len(temp_val.shape) not in [4, 5]:
            raise ValueError(
                "invalid data shape, must be (nbnd, ngx, ngy, ngz"
                ") or (nbnd, 2, ngx, ngy, ngz) for noncollinear "
                f"data, given {temp_val.shape}"
            )
        if len(temp_val.shape) == 5 and temp_val.shape[1] != 2:
            raise ValueError(
                f"invalid noncollinear data, shape should be (nbnd, 2, ngx, ngy, ngz), given {temp_val.shape}"
            )
        self._data = temp_val

        # derived properties
        self.is_noncollinear = len(self.data.shape) == 5
        self.nbnd = self.data.shape[0]
        self.ng = self.data.shape[-3:]

    @staticmethod
    def from_file(filename: str) -> object:
        """
        Reads the UNK data from file.

        Args:
            filename (str): path to UNK file to read

        Returns:
            Unk object
        """
        input_data = []
        with FortranFile(filename, "r") as f:
            *ng, ik, nbnd = f.read_ints()
            for _ in range(nbnd):
                input_data.append(
                    # when reshaping need to specify ordering as fortran
                    f.read_record(np.complex128).reshape(ng, order="F")
                )
            try:
                for _ in range(nbnd):
                    input_data.append(f.read_record(np.complex128).reshape(ng, order="F"))
                is_noncollinear = True
            except FortranEOFError:
                is_noncollinear = False

        # mypy made me create an extra variable here >:(
        data = np.array(input_data, dtype=np.complex128)

        # spinors are interwoven, need to separate them
        if is_noncollinear:
            temp_data = np.empty((nbnd, 2, *ng), dtype=np.complex128)
            temp_data[:, 0, :, :, :] = data[::2, :, :, :]
            temp_data[:, 1, :, :, :] = data[1::2, :, :, :]
            return Unk(ik, temp_data)
        return Unk(ik, data)

    def write_file(self, filename: str) -> None:
        """
        Write the UNK file.

        Args:
            filename (str): path to UNK file to write, the name should have the
                form 'UNKXXXXX.YY' where XXXXX is the kpoint index (Unk.ik) and
                YY is 1 or 2 for the spin index or NC if noncollinear
        """
        with FortranFile(filename, "w") as f:
            f.write_record(np.array([*self.ng, self.ik, self.nbnd], dtype=np.int32))
            for ib in range(self.nbnd):
                if self.is_noncollinear:
                    f.write_record(self.data[ib, 0].flatten("F"))
                    f.write_record(self.data[ib, 1].flatten("F"))
                else:
                    f.write_record(self.data[ib].flatten("F"))

    def __repr__(self) -> str:
        return (
            f"<UNK ik={self.ik} nbnd={self.nbnd} ncl={self.is_noncollinear}"
            + f" ngx={self.ng[0]} ngy={self.ng[1]} ngz={self.ng[2]}>"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Unk):
            return NotImplemented

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
                if not (
                    np.allclose(self.data[ib, 0], other.data[ib, 0], atol=1e-4)
                    and np.allclose(self.data[ib, 1], other.data[ib, 1], atol=1e-4)
                ):
                    return False
            else:
                if not np.allclose(self.data[ib], other.data[ib], atol=1e-4):
                    return False
        return True
