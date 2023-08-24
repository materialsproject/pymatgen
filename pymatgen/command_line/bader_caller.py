"""This module implements an interface to the Henkelmann et al.'s excellent
Fortran code for calculating a Bader charge analysis.

This module depends on a compiled bader executable available in the path.
Please download the library at http://theory.cm.utexas.edu/henkelman/code/bader/
and follow the instructions to compile the executable.

If you use this module, please cite:

G. Henkelman, A. Arnaldsson, and H. Jonsson, "A fast and robust algorithm for
Bader decomposition of charge density", Comput. Mater. Sci. 36, 254-360 (2006).
"""

from __future__ import annotations

import os
import subprocess
import warnings
from glob import glob
from shutil import which
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any

import numpy as np
from monty.shutil import decompress_file
from monty.tempfile import ScratchDir

from pymatgen.io.common import VolumetricData
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Chgcar

if TYPE_CHECKING:
    from pymatgen.core import Structure

__author__ = "shyuepingong"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Beta"
__date__ = "4/5/13"


BADER_EXE = which("bader") or which("bader.exe")


class BaderAnalysis:
    """Performs Bader analysis for Cube files and VASP outputs.

    Attributes:
        data (list[dict]): Atomic data parsed from bader analysis. Each dictionary in the list has the keys:
            "atomic_vol", "min_dist", "charge", "x", "y", "z".
        vacuum_volume (float): Vacuum volume of the Bader analysis.
        vacuum_charge (float): Vacuum charge of the Bader analysis.
        nelectrons (int): Number of electrons of the Bader analysis.
        chgcar (Chgcar): Chgcar object associated with input CHGCAR file.
        atomic_densities (list[dict]): List of charge densities for each atom centered on the atom.
            Excess 0's are removed from the array to reduce its size. Each dictionary has the keys:
            "data", "shift", "dim", where "data" is the charge density array,
            "shift" is the shift used to center the atomic charge density, and
            "dim" is the dimension of the original charge density map.
    """

    def __init__(
        self,
        chgcar_filename: str = "",
        potcar_filename: str = "",
        chgref_filename: str = "",
        parse_atomic_densities: bool = False,
        cube_filename: str = "",
        bader_exe_path: str | None = BADER_EXE,
    ) -> None:
        """Initializes the Bader caller.

        Args:
            chgcar_filename (str): The filename of the CHGCAR.
            potcar_filename (str): The filename of the POTCAR.
            chgref_filename (str): The filename of the reference charge density.
            parse_atomic_densities (bool, optional): Enable atomic partition of the charge density.
                Charge densities are atom centered. Defaults to False.
            cube_filename (str, optional): The filename of the cube file.
            bader_exe_path (str, optional): The path to the bader executable.
        """
        bader_exe = which(bader_exe_path or "")
        if bader_exe is None:
            raise RuntimeError(
                "BaderAnalysis requires bader or bader.exe to be in the PATH or the absolute path "
                f"to the binary to be specified via {bader_exe_path=}. Download the binary at "
                "https://theory.cm.utexas.edu/henkelman/code/bader."
            )

        if not (cube_filename or chgcar_filename):
            raise ValueError("You must provide either a cube file or a CHGCAR")
        if cube_filename and chgcar_filename:
            raise ValueError(
                f"You cannot parse a cube and a CHGCAR at the same time.\n{cube_filename=}\n{chgcar_filename=}"
            )
        self.parse_atomic_densities = parse_atomic_densities

        with ScratchDir("."):
            if chgcar_filename:
                self.is_vasp = True

                # decompress the file if compressed
                fpath = chgcar_fpath = decompress_file(filepath=chgcar_filename) or chgcar_filename
                self.chgcar = Chgcar.from_file(chgcar_fpath)
                self.structure = self.chgcar.structure

                self.potcar = Potcar.from_file(potcar_filename) if potcar_filename else None
                self.natoms = self.chgcar.poscar.natoms

                chgref_fpath = decompress_file(filepath=chgref_filename) or chgref_filename
                self.reference_used = bool(chgref_filename)

                # List of nelects for each atom from potcar
                potcar_indices = []
                for idx, val in enumerate(self.natoms):
                    potcar_indices += [idx] * val
                self.nelects = (
                    [self.potcar[potcar_indices[i]].nelectrons for i in range(len(self.structure))]
                    if self.potcar
                    else []
                )

            else:
                self.is_vasp = False
                fpath = cube_fpath = decompress_file(filepath=cube_filename) or cube_filename
                self.cube = VolumetricData.from_cube(cube_fpath)
                self.structure = self.cube.structure
                self.nelects = []
                chgref_fpath = decompress_file(filepath=chgref_filename) or chgref_filename
                self.reference_used = bool(chgref_filename)

            args = [bader_exe, fpath]

            if chgref_fpath:
                args += ["-ref", chgref_fpath]

            if parse_atomic_densities:
                args += ["-p", "all_atom"]

            with subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, close_fds=True) as proc:
                stdout, stderr = proc.communicate()
                if proc.returncode != 0:
                    raise RuntimeError(
                        f"{bader_exe} exit code: {proc.returncode}, error message: {stderr!s}.\nstdout: {stdout!s}"
                        "Please check your bader installation."
                    )

            try:
                self.version = float(stdout.split()[5])
            except ValueError:
                self.version = -1  # Unknown
            if self.version < 1.0:
                warnings.warn(
                    "Your installed version of Bader is outdated, calculation of vacuum charge may be incorrect.",
                    UserWarning,
                )

            data = []
            with open("ACF.dat") as f:
                lines = f.readlines()
                headers = ("x", "y", "z", "charge", "min_dist", "atomic_vol")
                lines.pop(0)
                lines.pop(0)
                while True:
                    line = lines.pop(0).strip()
                    if line.startswith("-"):
                        break
                    vals = map(float, line.split()[1:])
                    data.append(dict(zip(headers, vals)))
                for line in lines:
                    tokens = line.strip().split(":")
                    if tokens[0] == "VACUUM CHARGE":
                        self.vacuum_charge = float(tokens[1])
                    elif tokens[0] == "VACUUM VOLUME":
                        self.vacuum_volume = float(tokens[1])
                    elif tokens[0] == "NUMBER OF ELECTRONS":
                        self.nelectrons = float(tokens[1])
            self.data = data

            if self.parse_atomic_densities:
                # convert the charge density for each atom spit out by Bader into Chgcar objects for easy parsing
                atom_chgcars = [Chgcar.from_file(f"BvAt{idx + 1:04}.dat") for idx in range(len(self.chgcar.structure))]

                atomic_densities = []
                # For each atom in the structure
                for _, loc, chg in zip(
                    self.chgcar.structure,
                    self.chgcar.structure.frac_coords,
                    atom_chgcars,
                ):
                    # Find the index of the atom in the charge density atom
                    index = np.round(np.multiply(loc, chg.dim))

                    # Find the shift vector in the array
                    shift = (np.divide(chg.dim, 2) - index).astype(int)

                    # Shift the data so that the atomic charge density to the center for easier manipulation
                    shifted_data = np.roll(chg.data["total"], shift, axis=(0, 1, 2))

                    # Slices a central window from the data array
                    def slice_from_center(data, x_width, y_width, z_width):
                        x, y, z = data.shape
                        start_x = x // 2 - (x_width // 2)
                        start_y = y // 2 - (y_width // 2)
                        start_z = z // 2 - (z_width // 2)
                        return data[
                            start_x : start_x + x_width,
                            start_y : start_y + y_width,
                            start_z : start_z + z_width,
                        ]

                    def find_encompassing_vol(data: np.ndarray):
                        # Find the central encompassing volume which holds all the data within a precision
                        total = np.sum(data)
                        for idx in range(np.max(data.shape)):
                            sliced_data = slice_from_center(data, idx, idx, idx)
                            if total - np.sum(sliced_data) < 0.1:
                                return sliced_data
                        return None

                    dct = {
                        "data": find_encompassing_vol(shifted_data),
                        "shift": shift,
                        "dim": self.chgcar.dim,
                    }
                    atomic_densities.append(dct)
                self.atomic_densities = atomic_densities

    def get_charge(self, atom_index: int) -> float:
        """Convenience method to get the charge on a particular atom. This is the "raw"
        charge generated by the Bader program, not a partial atomic charge. If the cube file
        is a spin-density file, then this will return the spin density per atom with
        positive being spin up and negative being spin down.

        Args:
            atom_index (int): Index of atom.

        Returns:
            float: Charge associated with atom from the Bader analysis.
        """
        return self.data[atom_index]["charge"]

    def get_charge_transfer(self, atom_index: int, nelect: int | None = None) -> float:
        """Returns the charge transferred for a particular atom. A positive value means
        that the site has gained electron density (i.e. exhibits anionic character)
        whereas a negative value means the site has lost electron density (i.e. exhibits
        cationic character). If the arg nelect is not supplied, then POTCAR must be
        supplied to determine nelect.

        Args:
            atom_index (int): Index of atom.
            nelect (int): number of electrons associated with an isolated atom at this index.
                For most DFT codes this corresponds to the number of valence electrons
                associated with the pseudopotential (e.g. ZVAL for VASP).

        Returns:
            Charge transfer associated with atom from the Bader analysis.
            Given by bader charge on atom - nelect for associated atom.
        """
        if not self.nelects and nelect is None:
            raise ValueError("No NELECT info! Need POTCAR for VASP or nelect argument for cube file")
        return self.data[atom_index]["charge"] - (nelect if nelect is not None else self.nelects[atom_index])

    def get_partial_charge(self, atom_index: int, nelect: int | None = None) -> float:
        """Convenience method to get the partial charge on a particular atom. This is
        simply the negative value of the charge transferred. A positive value indicates
        that the atom has cationic character, whereas a negative value indicates the
        site has anionic character.

        Args:
            atom_index (int): Index of atom.
            nelect (int): number of electrons associated with an isolated atom at this index.
                For most DFT codes this corresponds to the number of valence electrons
                associated with the pseudopotential (e.g. ZVAL for VASP).

        Returns:
            Charge associated with atom from the Bader analysis.
        """
        return -self.get_charge_transfer(atom_index, nelect)

    def get_charge_decorated_structure(self) -> Structure:
        """Returns a charge decorated structure.

        Note, this assumes that the Bader analysis was correctly performed on a file
        with electron densities
        """
        charges = [-self.get_charge(i) for i in range(len(self.structure))]
        struct = self.structure.copy()
        struct.add_site_property("charge", charges)
        return struct

    def get_oxidation_state_decorated_structure(self, nelects: list[int] | None = None) -> Structure:
        """Returns an oxidation state decorated structure based on bader analysis results.
        Each site is assigned a charge based on the computed partial atomic charge from bader.

        Note, this assumes that the Bader analysis was correctly performed on a file
        with electron densities.

        Args:
            nelects (list[int]): number of electrons associated with an isolated atom at this index.

        Returns:
            Structure: with bader-analysis-based oxidation states.
        """
        charges = [self.get_partial_charge(i, None if not nelects else nelects[i]) for i in range(len(self.structure))]
        struct = self.structure.copy()
        struct.add_oxidation_state_by_site(charges)
        return struct

    def get_decorated_structure(self, property_name: str, average: bool = False) -> Structure:
        """Get a property-decorated structure from the Bader analysis.

        This is distinct from getting charge decorated structure, which assumes
        the "standard" Bader analysis of electron densities followed by converting
        electron count to charge. The expected way to use this is to call Bader on
        a non-charge density file such as a spin density file, electrostatic potential
        file, etc., while using the charge density file as the reference (chgref_filename)
        so that the partitioning is determined via the charge, but averaging or integrating
        is done for another property.

        User warning: Bader analysis cannot automatically determine what property is
        inside of the file. So if you want to use this for a non-conventional property
        like spin, you must ensure that you have the file is for the appropriate
        property and you have an appropriate reference file.

        Args:
            property_name (str): name of the property to assign to the structure, note that
                if name is "spin" this is handled as a special case, and the appropriate
                spin properties are set on the species in the structure
            average (bool): whether or not to return the average of this property, rather
                than the total, by dividing by the atomic volume.

        Returns:
            structure with site properties assigned via Bader Analysis
        """
        vals = np.array([self.get_charge(i) for i in range(len(self.structure))])
        struct = self.structure.copy()
        if average:
            vals = np.divide(vals, [d["atomic_vol"] for d in self.data])
        struct.add_site_property(property_name, vals)
        if property_name == "spin":
            struct.add_spin_by_site(vals)
        return struct

    @property
    def summary(self) -> dict[str, Any]:
        """Dict summary of key analysis, e.g., atomic volume, charge, etc."""
        summary = {
            "min_dist": [d["min_dist"] for d in self.data],
            "charge": [d["charge"] for d in self.data],
            "atomic_volume": [d["atomic_vol"] for d in self.data],
            "vacuum_charge": self.vacuum_charge,
            "vacuum_volume": self.vacuum_volume,
            "reference_used": self.reference_used,
            "bader_version": self.version,
        }

        if self.parse_atomic_densities:
            summary["charge_densities"] = self.atomic_densities

        if self.potcar:
            charge_transfer = [self.get_charge_transfer(i) for i in range(len(self.data))]
            summary["charge_transfer"] = charge_transfer

        return summary

    @classmethod
    def from_path(cls, path: str, suffix: str = "") -> BaderAnalysis:
        """Convenient constructor that takes in the path name of VASP run
        to perform Bader analysis.

        Args:
            path (str): Name of directory where VASP output files are
                stored.
            suffix (str): specific suffix to look for (e.g. '.relax1' for 'CHGCAR.relax1.gz').

        Returns:
            BaderAnalysis
        """

        def _get_filepath(filename):
            name_pattern = filename + suffix + "*" if filename != "POTCAR" else filename + "*"
            paths = glob(f"{path}/{name_pattern}")
            fpath = ""
            if len(paths) >= 1:
                # using reverse=True because, if multiple files are present, they likely
                # have suffixes 'static', 'relax', 'relax2', etc. and this would give
                # 'static' over 'relax2' over 'relax' however, better to use 'suffix'
                # kwarg to avoid this!
                paths.sort(reverse=True)
                if len(paths) > 1:
                    warnings.warn(f"Multiple files detected, using {paths[0]}")
                fpath = paths[0]
            else:
                msg = f"Could not find {filename!r}"
                if filename in ["AECCAR0", "AECCAR2"]:
                    msg += ", cannot calculate charge transfer."
                elif filename == "POTCAR":
                    msg += ", interpret Bader results with caution."
                warnings.warn(msg)
            return fpath

        chgcar_filename = _get_filepath("CHGCAR")
        if chgcar_filename is None:
            raise FileNotFoundError("Could not find CHGCAR!")
        potcar_filename = _get_filepath("POTCAR")
        aeccar0 = _get_filepath("AECCAR0")
        aeccar2 = _get_filepath("AECCAR2")
        if aeccar0 and aeccar2:
            # `chgsum.pl AECCAR0 AECCAR2` equivalent to obtain chgref_file
            chgref = Chgcar.from_file(aeccar0) + Chgcar.from_file(aeccar2)
            chgref_filename = f"{path}/CHGCAR_ref"
            chgref.write_file(chgref_filename)
        else:
            chgref_filename = ""
        return cls(
            chgcar_filename=chgcar_filename,
            potcar_filename=potcar_filename,
            chgref_filename=chgref_filename,
        )


def bader_analysis_from_path(path, suffix=""):
    """Convenience method to run Bader analysis on a folder containing
    typical VASP output files.

    This method will:

    1. Look for files CHGCAR, AECCAR0, AECCAR2, POTCAR or their gzipped
    counterparts.
    2. If AECCAR* files are present, constructs a temporary reference
    file as AECCAR0 + AECCAR2
    3. Runs Bader analysis twice: once for charge, and a second time
    for the charge difference (magnetization density).

    :param path: path to folder to search in
    :param suffix: specific suffix to look for (e.g. '.relax1' for 'CHGCAR.relax1.gz'
    Returns:
        summary dict
    """

    def _get_filepath(filename, warning, path=path, suffix=suffix):
        paths = glob(f"{path}/{filename}{suffix}*")
        if not paths:
            warnings.warn(warning)
            return None
        if len(paths) > 1:
            # using reverse=True because, if multiple files are present,
            # they likely have suffixes 'static', 'relax', 'relax2', etc.
            # and this would give 'static' over 'relax2' over 'relax'
            # however, better to use 'suffix' kwarg to avoid this!
            paths.sort(reverse=True)
            warnings.warn(f"Multiple files detected, using {os.path.basename(path)}")
        return paths[0]

    chgcar_path = _get_filepath("CHGCAR", "Could not find CHGCAR!")
    chgcar = Chgcar.from_file(chgcar_path)

    aeccar0_path = _get_filepath("AECCAR0", "Could not find AECCAR0, interpret Bader results with caution.")
    aeccar0 = Chgcar.from_file(aeccar0_path) if aeccar0_path else None

    aeccar2_path = _get_filepath("AECCAR2", "Could not find AECCAR2, interpret Bader results with caution.")
    aeccar2 = Chgcar.from_file(aeccar2_path) if aeccar2_path else None

    potcar_path = _get_filepath("POTCAR", "Could not find POTCAR, cannot calculate charge transfer.")
    potcar = Potcar.from_file(potcar_path) if potcar_path else None

    return bader_analysis_from_objects(chgcar, potcar, aeccar0, aeccar2)


def bader_analysis_from_objects(chgcar, potcar=None, aeccar0=None, aeccar2=None):
    """Convenience method to run Bader analysis from a set
    of pymatgen Chgcar and Potcar objects.

    This method will:

    1. If aeccar objects are present, constructs a temporary reference
    file as AECCAR0 + AECCAR2
    2. Runs Bader analysis twice: once for charge, and a second time
    for the charge difference (magnetization density).

    :param chgcar: Chgcar object
    :param potcar: (optional) Potcar object
    :param aeccar0: (optional) Chgcar object from aeccar0 file
    :param aeccar2: (optional) Chgcar object from aeccar2 file

    Returns:
        summary dict
    """
    orig_dir = os.getcwd()
    try:
        with TemporaryDirectory() as tmp_dir:
            os.chdir(tmp_dir)
            if aeccar0 and aeccar2:
                # construct reference file
                chgref = aeccar0.linear_add(aeccar2)
                chgref_path = f"{tmp_dir}/CHGCAR_ref"
                chgref.write_file(chgref_path)
            else:
                chgref_path = ""

            chgcar.write_file("CHGCAR")
            chgcar_path = f"{tmp_dir}/CHGCAR"

            if potcar:
                potcar.write_file("POTCAR")
                potcar_path = f"{tmp_dir}/POTCAR"
            else:
                potcar_path = ""

            ba = BaderAnalysis(
                chgcar_filename=chgcar_path,
                potcar_filename=potcar_path,
                chgref_filename=chgref_path,
            )

            summary = {
                "min_dist": [d["min_dist"] for d in ba.data],
                "charge": [d["charge"] for d in ba.data],
                "atomic_volume": [d["atomic_vol"] for d in ba.data],
                "vacuum_charge": ba.vacuum_charge,
                "vacuum_volume": ba.vacuum_volume,
                "reference_used": bool(chgref_path),
                "bader_version": ba.version,
            }

            if potcar:
                charge_transfer = [ba.get_charge_transfer(i) for i in range(len(ba.data))]
                summary["charge_transfer"] = charge_transfer

            if chgcar.is_spin_polarized:
                # write a CHGCAR containing magnetization density only
                chgcar.data["total"] = chgcar.data["diff"]
                chgcar.is_spin_polarized = False
                chgcar.write_file("CHGCAR_mag")

                chgcar_mag_path = f"{tmp_dir}/CHGCAR_mag"
                ba = BaderAnalysis(
                    chgcar_filename=chgcar_mag_path,
                    potcar_filename=potcar_path,
                    chgref_filename=chgref_path,
                )
                summary["magmom"] = [d["charge"] for d in ba.data]
    finally:
        os.chdir(orig_dir)

    return summary
