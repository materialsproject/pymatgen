"""
This module provides a pymatgen I/O interface to packmol.

This adopts the minimal core I/O interface (see pymatgen/io/core).
In this case, only a two classes are used. PackmolSet(InputSet) is the container
class that provides a run() method for running packmol locally.

PackmolBoxGen(InputGenerator) provides a recipe for packing molecules into a
box, which returns a PackmolSet object.

For the run() method to work, you need to install the packmol package
See http://m3g.iqm.unicamp.br/packmol or
http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml
for download and setup instructions. Note that packmol versions prior to 20.3.0
do not support paths with spaces.
After installation, you may need to manually add the path of the packmol
executable to the PATH environment variable.
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
from monty.os.path import which

from pymatgen.core import Molecule
from pymatgen.io.core import InputGenerator, InputSet

__author__ = "Tingzheng Hou, Ryan Kingsbury, Orion Cohen"
__version__ = "1.0"
__maintainer__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__date__ = "Nov 2021"


class PackmolSet(InputSet):
    """
    InputSet for the Packmol software. This class defines several attributes related
    to
    """

    def run(self, path: Union[str, Path], timeout=30):
        """Run packmol and write out the packed structure.
        Args:
            path: The path in which packmol input files are located.
            timeout: Timeout in seconds.
        Raises:
            ValueError if packmol does not succeed in packing the box.
            TimeoutExpiredError if packmold does not finish within the timeout.
        """
        wd = os.getcwd()
        if not which("packmol"):
            raise RuntimeError(
                "Running a PackmolSet requires the executable 'packmol' to be in "
                "the path. Please download packmol from "
                "https://github.com/leandromartinez98/packmol "
                "and follow the instructions in the README to compile. "
                "Don't forget to add the packmol binary to your path"
            )
        try:
            os.chdir(path)
            p = subprocess.run(
                f"packmol < '{self.inputfile}'",
                check=True,
                shell=True,
                timeout=timeout,
                capture_output=True,
            )
            # this workaround is needed because packmol can fail to find
            # a solution but still return a zero exit code
            # see https://github.com/m3g/packmol/issues/28
            if "ERROR" in p.stdout.decode():
                if "Could not open file." in p.stdout.decode():
                    raise ValueError(
                        "Your packmol might be too old to handle paths with spaces."
                        "Please try again with a newer version or use paths without spaces."
                    )
                msg = p.stdout.decode().split("ERROR")[-1]
                raise ValueError(f"Packmol failed with return code 0 and stdout: {msg}")
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Packmol failed with errorcode {e.returncode} and stderr: {e.stderr}") from e
        else:
            with open(Path(path, self.stdoutfile), "w") as out:
                out.write(p.stdout.decode())
        finally:
            os.chdir(wd)

    @classmethod
    def from_directory(cls, directory: Union[str, Path]):
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory: Directory to read input files from
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls}")


class PackmolBoxGen(InputGenerator):
    """
    Generator for a Packmol InputSet that packs one or more molecules into a rectangular
    simulation box.
    """

    def __init__(
        self,
        tolerance: float = 2.0,
        seed: int = 1,
        control_params: Optional[Dict] = None,
        inputfile: Union[str, Path] = "packmol.inp",
        outputfile: Union[str, Path] = "packmol_out.xyz",
        stdoutfile: Union[str, Path] = "packmol.stdout",
    ):
        """
        Instantiate a PackmolBoxGen class. The init method defines simulations parameters
        like filenames, random seed, tolerance, etc.

        Args:
            tolerance: Tolerance for packmol, in Å.
            seed: Random seed for packmol. Use a value of 1 (default) for deterministic
                output, or -1 to generate a new random seed from the current time.
            inputfile: Path to the input file. Default to 'packmol.inp'.
            outputfile: Path to the output file. Default to 'output.xyz'.
            stdoutfile: Path to the file where stdout will be recorded. Default to 'packmol.stdout'
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.stdoutfile = stdoutfile
        self.control_params = control_params if control_params else {}
        self.tolerance = tolerance
        self.seed = seed

    def get_input_set(  # type: ignore
        self,
        molecules: List[Dict],
        box: Optional[List[float]] = None,
    ):
        """
        Generate a Packmol InputSet for a set of molecules.

        Args:
            molecules: A list of dict containing information about molecules to pack
                into the box. Each dict requires three keys:
                    1. "name" - the structure name
                    2. "number" - the number of that molecule to pack into the box
                    3. "coords" - Coordinates in the form of either a Molecule object or
                        a path to a file.
                Example:
                    {"name": "water",
                     "number": 500,
                     "coords": "/path/to/input/file.xyz"}
            box: A list of box dimensions xlo, ylo, zlo, xhi, yhi, zhi, in Å. If set to None
                (default), pymatgen will estimate the required box size based on the volumes of
                the provided molecules.
        """
        mapping = {}
        file_contents = "# Packmol input generated by pymatgen.\n"
        file_contents += "# " + " + ".join(str(d["number"]) + " " + d["name"] for d in molecules) + "\n"

        for k, v in self.control_params.items():
            if isinstance(v, list):
                file_contents += f"{k} {' '.join(str(x) for x in v)}\n"
            else:
                file_contents += f"{k} {str(v)}\n"
        file_contents += f"seed {self.seed}\n"
        file_contents += f"tolerance {self.tolerance}\n\n"

        file_contents += "filetype xyz\n\n"
        if " " in str(self.outputfile):
            # NOTE - double quotes are deliberately used inside the f-string here, do not change
            # fmt: off
            file_contents += f'output "{self.outputfile}"\n\n'
            # fmt: on
        else:
            file_contents += f"output {self.outputfile}\n\n"

        if box:
            box_list = " ".join(str(i) for i in box)
        else:
            # estimate the total volume of all molecules in cubic Å
            net_volume = 0.0
            for d in molecules:
                if not isinstance(d["coords"], Molecule):
                    mol = Molecule.from_file(d["coords"])
                else:
                    mol = d["coords"]
                # pad the calculated length by an amount related to the tolerance parameter
                # the amount to add was determined arbitrarily
                length = (
                    max(np.max(mol.cart_coords[:, i]) - np.min(mol.cart_coords[:, i]) for i in range(3))
                    + self.tolerance
                )
                net_volume += (length**3.0) * float(d["number"])
            box_length = net_volume ** (1.0 / 3.0)
            print(f"Auto determined box size is {box_length:.1f} Å per side.")
            box_list = f"0.0 0.0 0.0 {box_length:.1f} {box_length:.1f} {box_length:.1f}"

        for d in molecules:
            if isinstance(d["coords"], str):
                mol = Molecule.from_file(d["coords"])
            elif isinstance(d["coords"], Path):
                mol = Molecule.from_file(str(d["coords"]))
            elif isinstance(d["coords"], Molecule):
                mol = d["coords"]
            fname = f"packmol_{d['name']}.xyz"
            mapping.update({fname: mol.to(fmt="xyz")})
            if " " in str(fname):
                # NOTE - double quotes are deliberately used inside the f-string here, do not change
                # fmt: off
                file_contents += f'structure "{fname}"\n'
                # fmt: on
            else:
                file_contents += f"structure {fname}\n"
            file_contents += f"  number {str(d['number'])}\n"
            file_contents += f"  inside box {box_list}\n"
            file_contents += "end structure\n\n"

        mapping.update({str(self.inputfile): file_contents})

        return PackmolSet(
            inputs=mapping,  # type: ignore
            seed=self.seed,
            inputfile=self.inputfile,
            outputfile=self.outputfile,
            stdoutfile=self.stdoutfile,
            control_params=self.control_params,
            tolerance=self.tolerance,
        )
