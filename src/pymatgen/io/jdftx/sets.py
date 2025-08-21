from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from pymatgen.io.core import InputSet
from pymatgen.io.jdftx.inputs import JDFTXInfile, JDFTXStructure

if TYPE_CHECKING:
    from pymatgen.core import Structure
    from pymatgen.util.typing import PathLike

FILE_NAMES = {"in": "init.in", "out": "jdftx.out"}


class JdftxInputSet(InputSet):
    """
    A class to represent a JDFTx input file as a JDFTx InputSet.

    Parameters
    ----------
    jdftxinput
        A JdftxInput object
    """

    def __init__(self, jdftxinput: JDFTXInfile, structure: Structure) -> None:
        self.structure = structure
        self.jdftxinput = jdftxinput

    def write_input(
        self,
        directory: str | Path,
        infile: PathLike = FILE_NAMES["in"],
        make_dir: bool = True,
        overwrite: bool = True,
    ) -> None:
        """Write JDFTx input file to a directory.

        Parameters
        ----------
        directory
            Directory to write input files to.
        make_dir
            Whether to create the directory if it does not already exist.
        overwrite
            Whether to overwrite an input file if it already exists.
        """
        directory = Path(directory)
        if make_dir:
            os.makedirs(directory, exist_ok=True)

        if not overwrite and (directory / infile).exists():
            raise FileExistsError(f"{directory / infile} already exists.")
        jdftx_structure = JDFTXStructure(structure=self.structure)
        jdftxinput = condense_jdftxinputs(self.jdftxinput, jdftx_structure)

        jdftxinput.write_file(filename=(directory / infile))

    @staticmethod
    def from_file(
        file: PathLike,
    ) -> JdftxInputSet:
        """Load a set of JDFTx inputs from a filename.

        Parameters
        ----------
        directory
            Input file to read JDFTx inputs from.
        """
        jdftxinput = JDFTXInfile.from_file(file)
        structure = jdftxinput.structure
        if structure is None:
            raise ValueError(f"Structure not defined in file {file}.")
        return JdftxInputSet(jdftxinput=jdftxinput, structure=structure)


def condense_jdftxinputs(jdftxinput: JDFTXInfile, jdftxstructure: JDFTXStructure) -> JDFTXInfile:
    """
    Combine JDFTXInfile and JDFTxStructure into complete JDFTXInfile.

    Function combines a JDFTXInfile class with calculation
    settings and a JDFTxStructure that defines the structure
    into one JDFTXInfile instance.

    Parameters
    ----------
        jdftxinput: JDFTXInfile
            A JDFTXInfile object with calculation settings.

        jdftxstructure: JDFTXStructure
            A JDFTXStructure object that defines the structure.

    Returns
    -------
        JDFTXInfile
            A JDFTXInfile that includes the calculation
            parameters and input structure.
    """
    # force Cartesian coordinates
    coords_type = jdftxinput.get("coords-type")
    return jdftxinput + JDFTXInfile.from_str(jdftxstructure.get_str(in_cart_coords=(coords_type == "Cartesian")))
