"""Modules for working with wannier90 input and output."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.io import FortranEOFError, FortranFile

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

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

    Attributes:
        ik (int): Index of kpoint for this file.
        data (numpy.ndarray): Numpy array that contains the wavefunction data in the UNK file.
            The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and (nbnd, 2, ngx, ngy, ngz)
            for noncollinear calculations.
        is_noncollinear (bool): True if data is from a noncollinear calculation.
        nbnd (int): Number of bands in data.
        ng (tuple): Sequence of three integers that correspond to the grid size of the given data.
            The definition is ng = (ngx, ngy, ngz).
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
        np.ndarray: contains the wavefunction data in the UNK file.
        The shape should be (nbnd, ngx, ngy, ngz) for regular calculations and
        (nbnd, 2, ngx, ngy, ngz) for noncollinear calculations.
        """
        return self._data

    @data.setter
    def data(self, value: np.ndarray) -> None:
        """Set the value of data.

        Args:
            value (np.ndarray): data to replace stored data, must have shape
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

    @classmethod
    def from_file(cls, filename: str) -> Self:
        """
        Reads the UNK data from file.

        Args:
            filename (str): path to UNK file to read

        Returns:
            Unk object
        """
        input_data = []
        with FortranFile(filename, mode="r") as file:
            *ng, ik, nbnd = file.read_ints()
            for _ in range(nbnd):
                input_data.append(
                    # when reshaping need to specify ordering as fortran
                    file.read_record(np.complex128).reshape(ng, order="F")
                )
            try:
                for _ in range(nbnd):
                    input_data.append(file.read_record(np.complex128).reshape(ng, order="F"))
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
            return cls(ik, temp_data)
        return cls(ik, data)

    def write_file(self, filename: str) -> None:
        """Write the UNK file.

        Args:
            filename (str): path to UNK file to write, the name should have the
                form 'UNKXXXXX.YY' where XXXXX is the kpoint index (Unk.ik) and
                YY is 1 or 2 for the spin index or NC if noncollinear
        """
        with FortranFile(filename, mode="w") as file:
            file.write_record(np.array([*self.ng, self.ik, self.nbnd], dtype=np.int32))
            for ib in range(self.nbnd):
                if self.is_noncollinear:
                    file.write_record(self.data[ib, 0].flatten("F"))
                    file.write_record(self.data[ib, 1].flatten("F"))
                else:
                    file.write_record(self.data[ib].flatten("F"))

    def __repr__(self) -> str:
        ik, nbnd, ncl, ngx, ngy, ngz = (
            self.ik,
            self.nbnd,
            self.is_noncollinear,
            *self.ng,
        )
        return f"{(type(self).__name__)}({ik=}, {nbnd=}, {ncl=}, {ngx=}, {ngy=}, {ngz=})"

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
            elif not np.allclose(self.data[ib], other.data[ib], atol=1e-4):
                return False
        return True
