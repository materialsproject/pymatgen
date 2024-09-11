"""
This module provides a pymatgen I/O interface to PACKMOL.

- PackmolSet provides a "run" method to run PACKMOL locally.
- PackmolBoxGen provides "get_input_set" for packing molecules into a box,
which returns a PackmolSet object.

For the run() method to work, you need to install the PACKMOL package.
See https://m3g.iqm.unicamp.br/packmol for download and setup instructions.
After installation, you may need to add the path of the PACKMOL
executable to the PATH environment variable.

Note that PACKMOL versions prior to 20.3.0 do not support paths with spaces.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from shutil import which
from typing import TYPE_CHECKING

import numpy as np

from pymatgen.core import Molecule
from pymatgen.io.core import InputGenerator, InputSet

if TYPE_CHECKING:
    from pymatgen.util.typing import PathLike

__author__ = "Tingzheng Hou, Ryan Kingsbury, Orion Cohen"
__version__ = "1.0"
__maintainer__ = "Ryan Kingsbury"
__email__ = "RKingsbury@lbl.gov"
__date__ = "Nov 2021"


class PackmolSet(InputSet):
    """InputSet for the PACKMOL software. This class defines several attributes related to."""

    def run(self, path: PathLike, timeout: float = 30) -> None:
        """Run PACKMOL and write out the packed structure.

        Args:
            path (PathLike): The path in which packmol input files are located.
            timeout (float): Timeout in seconds.

        Raises:
            ValueError: if packmol does not succeed in packing the box.
            TimeoutExpiredError: if packmol does not finish within the timeout.
        """
        wd = os.getcwd()
        if not which("packmol"):
            raise RuntimeError(
                "Running a PackmolSet requires the executable 'packmol' to be in the path. Please "
                "download packmol from https://github.com/leandromartinez98/packmol and follow the "
                "instructions in the README to compile. Don't forget to add the packmol binary to your path"
            )
        try:
            os.chdir(path)
            with open(self.inputfile, encoding="utf-8") as infile:
                proc = subprocess.run(
                    ["packmol"],
                    stdin=infile,
                    check=True,
                    timeout=timeout,
                    capture_output=True,
                )
            # This workaround is needed because packmol can fail to find
            # a solution but still return a zero exit code.
            # See https://github.com/m3g/packmol/issues/28
            if "ERROR" in proc.stdout.decode():
                if "Could not open file." in proc.stdout.decode():
                    raise ValueError(
                        "Your packmol might be too old to handle paths with spaces."
                        "Please try again with a newer version or use paths without spaces."
                    )
                msg = proc.stdout.decode().split("ERROR")[-1]
                raise ValueError(f"Packmol failed with return code 0 and stdout: {msg}")

        except subprocess.CalledProcessError as exc:
            raise ValueError(f"Packmol failed with error code {exc.returncode} and stderr: {exc.stderr}") from exc
        else:
            with open(Path(path, self.stdoutfile), mode="w", encoding="utf-8") as out:
                out.write(proc.stdout.decode())
        finally:
            os.chdir(wd)

    @classmethod
    def from_directory(cls, directory: PathLike) -> None:
        """
        Construct an InputSet from a directory of one or more files.

        Args:
            directory (PathLike): Directory to read input files from.
        """
        raise NotImplementedError(f"from_directory has not been implemented in {cls.__name__}")


class PackmolBoxGen(InputGenerator):
    """
    Generator for a Packmol InputSet that packs one or more molecules into a rectangular
    simulation box.
    """

    def __init__(
        self,
        tolerance: float = 2.0,
        seed: int = 1,
        control_params: dict | None = None,
        inputfile: PathLike = "packmol.inp",
        outputfile: PathLike = "packmol_out.xyz",
        stdoutfile: PathLike = "packmol.stdout",
    ) -> None:
        """
        Instantiate a PackmolBoxGen class. The init method defines simulations parameters
        like filenames, random seed, tolerance, etc.

        Args:
            tolerance (float): Tolerance for packmol, in Å.
            seed (int): Random seed for packmol. Use 1 (default) for deterministic
                output, or -1 to generate a new random seed from the current time.
            inputfile (PathLike): Path to the input file. Default to "packmol.inp".
            outputfile (PathLike): Path to the output file. Default to "packmol_out.xyz".
            stdoutfile (PathLike): Path to the file where stdout will be recorded. Default to "packmol.stdout".
        """
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.stdoutfile = stdoutfile
        self.control_params = control_params or {}
        self.tolerance = tolerance
        self.seed = seed

    def get_input_set(
        self,
        molecules: list[dict],
        box: list[float] | None = None,
    ) -> PackmolSet:
        """Generate a Packmol InputSet for a set of molecules.

        Args:
            molecules (list[dict]): Information about molecules to pack
                into the box. Each dict requires three keys:
                    1. "name" - the structure name.
                    2. "number" - the number of that molecule to pack into the box.
                    3. "coords" - Coordinates in the form of either a Molecule
                        object or a path to a file.
                Example:
                    {
                        "name": "water",
                        "number": 500,
                        "coords": "/path/to/input/file.xyz",
                    }
            box (list[float]): Box dimensions xlo, ylo, zlo, xhi, yhi, zhi, in Å. If set to None
                (default), pymatgen will estimate the required box size based on the volumes of
                the provided molecules.

        Returns:
            PackmolSet
        """
        mapping: dict = {}
        file_contents: list[str] = [
            "# Packmol input generated by pymatgen.\n",
            f"# {' + '.join(str(d['number']) + ' ' + d['name'] for d in molecules)}",
        ]

        for key, val in self.control_params.items():
            if isinstance(val, list):
                file_contents.append(f"{key} {' '.join(str(x) for x in val)}")
            else:
                file_contents.append(f"{key} {val}")

        file_contents += [
            f"seed {self.seed}",
            f"tolerance {self.tolerance}\n",
            "filetype xyz\n",
        ]

        if " " in str(self.outputfile):
            # NOTE - double quotes are deliberately used inside the f-string here, do not change
            file_contents.append(f'output "{self.outputfile}"\n')
        else:
            file_contents.append(f"output {self.outputfile}\n")

        if box:
            box_list = " ".join(map(str, box))
        else:
            # Estimate the total volume of all molecules in cubic Å
            net_volume: float = 0.0
            for dct in molecules:
                mol = dct["coords"] if isinstance(dct["coords"], Molecule) else Molecule.from_file(dct["coords"])

                if mol is None:
                    raise ValueError("Molecule cannot be None.")
                # Pad the calculated length by an amount related to the tolerance parameter
                # the amount to add was determined arbitrarily
                length = (
                    max(np.max(mol.cart_coords[:, i]) - np.min(mol.cart_coords[:, i]) for i in range(3))
                    + self.tolerance
                )
                net_volume += (length**3.0) * float(dct["number"])
            box_length: float = net_volume ** (1 / 3)
            print(f"Auto determined box size is {box_length:.1f} Å per side.")
            box_list = f"0.0 0.0 0.0 {box_length:.1f} {box_length:.1f} {box_length:.1f}"

        for dct in molecules:
            if isinstance(dct["coords"], str | Path):
                mol = Molecule.from_file(dct["coords"])
            elif isinstance(dct["coords"], Molecule):
                mol = dct["coords"]
            else:
                raise TypeError("Molecule is not provided in supported format.")

            fname = f"packmol_{dct['name']}.xyz"
            if mol is None:
                raise ValueError("mol is None")
            mapping[fname] = mol.to(fmt="xyz")
            if " " in str(fname):
                file_contents.append(f"structure {fname!r}")
            else:
                file_contents.append(f"structure {fname}")

            file_contents.extend(
                (
                    f"  number {dct['number']}",
                    f"  inside box {box_list}",
                    "end structure\n\n",
                )
            )

        mapping |= {str(self.inputfile): "\n".join(file_contents)}

        return PackmolSet(
            inputs=mapping,
            seed=self.seed,
            inputfile=self.inputfile,
            outputfile=self.outputfile,
            stdoutfile=self.stdoutfile,
            control_params=self.control_params,
            tolerance=self.tolerance,
        )
