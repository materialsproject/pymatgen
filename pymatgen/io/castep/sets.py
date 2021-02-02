from dataclasses import dataclass
from pathlib import Path

from monty.json import MSONable

from pymatgen.io.castep.inputs import Param, Cell


@dataclass
class CastepInputSet(MSONable):
    """
    A base class defining the interface for a CastepInputSet
    """

    seedname: str
    param: Param
    cell: Cell

    def write_input(self, output_dir="."):
        """
        Write input files.

        Args:
            output_dir: destination path, by default will be current direcotry
        """
        path = Path(output_dir)
        self.cell.write_file(path / self.seedname)
        self.param.write_file(path / self.seedname)


class MPStaticSet(CastepInputSet):
    """
    An input set for running static, total-energy calculations, designed
    to be roughly analogous to the MPStaticSet(VaspInputSet).

    This input set has not been extensively tested and is not used for
    any production calculations. Modifications are welcome.
    """


