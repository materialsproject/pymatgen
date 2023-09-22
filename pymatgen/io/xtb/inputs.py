"""Classes for writing XTB input files."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

from monty.json import MSONable

if TYPE_CHECKING:
    from pymatgen.core import Molecule

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Alex Epstein"
__email__ = "aepstein@lbl.gov"
__credits__ = "Sam Blau, Evan Spotte-Smith"

logger = logging.getLogger(__name__)


class CRESTInput(MSONable):
    """
    An object representing  CREST input files.
    Because CREST is controlled through command line flags and external
    files, the CRESTInput class mainly consists of methods for containing
    and writing external files.
    """

    def __init__(
        self,
        molecule: Molecule,
        working_dir: str = ".",
        coords_filename: str | None = "crest_in.xyz",
        constraints: dict[str, list[int] | float] | None = None,
    ):
        """
        :param molecule (pymatgen Molecule object):
            Input molecule, the only required CREST input.
        :param working_dir (str):
            Location to write input files, defaults to current directory
        :param coords_filename (str):
            Name of input coordinates file
        :param constraints (dict):
            Dictionary of common editable parameters for .constrains file.
            {"atoms": [List of 1-indexed atoms to fix], "force_constant":
            float]
        """
        self.molecule = molecule
        self.coords_filename = coords_filename
        self.constraints = constraints
        self.working_dir = working_dir

    def write_input_files(self):
        """Write input files to working directory."""
        self.molecule.to(filename=os.path.join(self.working_dir, self.coords_filename))
        if self.constraints:
            constrains_string = self.constrains_template(
                molecule=self.molecule,
                reference_fnm=self.coords_filename,
                constraints=self.constraints,
            )
            with open(".constrains", "w") as f:
                f.write(constrains_string)

    @staticmethod
    def constrains_template(molecule, reference_fnm, constraints) -> str:
        """
        :param molecule (pymatgen Molecule):
            Molecule the constraints will be performed on
        :param reference_fnm:
            Name of file containing reference structure in same directory
        :param constraints:
            Dictionary of common editable parameters for .constrains file.
            {"atoms": [List of 1-indexed atoms to fix], "force_constant": float]

        Returns:
            str: for .constrains file
        """
        atoms_to_constrain = constraints["atoms"]
        force_constant = constraints["force_constant"]
        mol = molecule
        atoms_for_mtd = [idx for idx in range(1, len(mol) + 1) if idx not in atoms_to_constrain]
        # Write as 1-3,5 instead of 1,2,3,5
        interval_list = [atoms_for_mtd[0]]
        for i, v in enumerate(atoms_for_mtd):
            if v + 1 not in atoms_for_mtd:
                interval_list.append(v)
                if i != len(atoms_for_mtd) - 1:
                    interval_list.append(atoms_for_mtd[i + 1])
        allowed_mtd_string = ",".join(
            [f"{interval_list[i]}-{interval_list[i + 1]}" for i in range(len(interval_list)) if i % 2 == 0]
        )
        return (
            "$constrain\n"
            f"  atoms: {','.join(map(str, atoms_to_constrain))}\n"
            f"  force constant={force_constant}\n"
            f"  reference={reference_fnm}\n$metadyn\n"
            f"  atoms: {allowed_mtd_string}\n$end"
        )
