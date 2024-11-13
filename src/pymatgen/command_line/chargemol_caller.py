"""This module implements an interface to Thomas Manz's Chargemol code
https://sourceforge.net/projects/ddec for calculating DDEC3, DDEC6, and CM5 population analyses.

This module depends on a compiled chargemol executable being available in the path.
If you use this module, please cite the following based on which modules you use:

Chargemol:
(1) T. A. Manz and N. Gabaldon Limas, Chargemol program for performing DDEC analysis,
Version 3.5, 2017, ddec.sourceforge.net.

DDEC6 Charges:
(1) T. A. Manz and N. Gabaldon Limas, “Introducing DDEC6 atomic population analysis:
part 1. Charge partitioning theory and methodology,” RSC Adv., 6 (2016) 47771-47801.
(2) N. Gabaldon Limas and T. A. Manz, “Introducing DDEC6 atomic population analysis:
part 2. Computed results for a wide range of periodic and nonperiodic materials,”
(3) N. Gabaldon Limas and T. A. Manz, “Introducing DDEC6 atomic population analysis:
part 4. Efficient parallel computation of net atomic charges, atomic spin moments,
bond orders, and more,” RSC Adv., 8 (2018) 2678-2707.

CM5 Charges:
(1) A.V. Marenich, S.V. Jerome, C.J. Cramer, D.G. Truhlar, "Charge Model 5: An Extension
of Hirshfeld Population Analysis for the Accurate Description of Molecular Interactions
in Gaseous and Condensed Phases", J. Chem. Theory. Comput., 8 (2012) 527-541.

Spin Moments:
(1) T. A. Manz and D. S. Sholl, “Methods for Computing Accurate Atomic Spin Moments for
Collinear and Noncollinear Magnetism in Periodic and Nonperiodic Materials,”
J. Chem. Theory Comput. 7 (2011) 4146-4164.

Bond Orders:
(1) “Introducing DDEC6 atomic population analysis: part 3. Comprehensive method to compute
bond orders,” RSC Adv., 7 (2017) 45552-45581.

DDEC3 Charges:
(1) T. A. Manz and D. S. Sholl, “Improved Atoms-in-Molecule Charge Partitioning Functional
for Simultaneously Reproducing the Electrostatic Potential and Chemical States in Periodic
and Non-Periodic Materials,” J. Chem. Theory Comput. 8 (2012) 2844-2867.
(2) T. A. Manz and D. S. Sholl, “Chemically Meaningful Atomic Charges that Reproduce the
Electrostatic Potential in Periodic and Nonperiodic Materials,” J. Chem. Theory Comput. 6
(2010) 2455-2468.
"""

from __future__ import annotations

import os
import platform
import subprocess
import warnings
from glob import glob
from shutil import which
from typing import TYPE_CHECKING

import numpy as np
from monty.tempfile import ScratchDir

from pymatgen.core import Element
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Chgcar

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from pymatgen.core import Structure

__author__ = "Martin Siron, Andrew S. Rosen"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "01/18/21"

CHARGEMOL_EXE = (
    which("Chargemol_09_26_2017_linux_parallel") or which("Chargemol_09_26_2017_linux_serial") or which("chargemol")
)


class ChargemolAnalysis:
    """Chargemol analysis for DDEC3, DDEC6, and/or CM5 population analyses,
    including the calculation of partial atomic charges, atomic spin moments,
    bond orders, and related properties.
    """

    def __init__(
        self,
        path: str | Path | None = None,
        atomic_densities_path: str | Path | None = None,
        run_chargemol: bool = True,
    ) -> None:
        """Initialize the Chargemol Analysis.

        Args:
            path (str): Path to the CHGCAR, POTCAR, AECCAR0, and AECCAR files.
                The files can be gzipped or not. Default: None (current working directory).
            atomic_densities_path (str | None): Path to the atomic densities directory
                required by Chargemol. If None, Pymatgen assumes that this is
                defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment variable.
                Only used if run_chargemol is True. Default: None.
            run_chargemol (bool): Whether to run the Chargemol analysis. If False,
                the existing Chargemol output files will be read from path. Default: True.
        """
        path = path or os.getcwd()
        if run_chargemol and not CHARGEMOL_EXE:
            raise OSError(
                "ChargemolAnalysis requires the Chargemol executable to be in PATH."
                " Please download the library at https://sourceforge.net/projects/ddec/files"
                "and follow the instructions."
            )
        if atomic_densities_path == "":
            atomic_densities_path = os.getcwd()
        self._atomic_densities_path = atomic_densities_path

        self._chgcar_path = self._get_filepath(path, "CHGCAR")
        self._potcar_path = self._get_filepath(path, "POTCAR")
        self._aeccar0_path = self._get_filepath(path, "AECCAR0")
        self._aeccar2_path = self._get_filepath(path, "AECCAR2")

        if run_chargemol and not (
            self._chgcar_path and self._potcar_path and self._aeccar0_path and self._aeccar2_path
        ):
            raise FileNotFoundError("CHGCAR, AECCAR0, AECCAR2, and POTCAR are all needed for Chargemol.")

        if self._chgcar_path:
            self.chgcar: Chgcar | None = Chgcar.from_file(self._chgcar_path)
            self.structure: Structure | None = self.chgcar.structure
            self.natoms: list[int] | None = self.chgcar.poscar.natoms

        else:
            self.chgcar = self.structure = self.natoms = None
            warnings.warn("No CHGCAR found. Some properties may be unavailable.", UserWarning)

        if self._potcar_path:
            self.potcar = Potcar.from_file(self._potcar_path)
        else:
            warnings.warn("No POTCAR found. Some properties may be unavailable.", UserWarning)

        self.aeccar0 = Chgcar.from_file(self._aeccar0_path) if self._aeccar0_path else None
        self.aeccar2 = Chgcar.from_file(self._aeccar2_path) if self._aeccar2_path else None

        if run_chargemol:
            self._execute_chargemol()
        else:
            self._from_data_dir(chargemol_output_path=path)

    @staticmethod
    def _get_filepath(path, filename, suffix=""):
        """Get the full path to the filename in the path. Works even if the file has
        a .gz extension.

        Args:
            path (str): Path to the file.
            filename (str): Filename.
            suffix (str): Optional suffix at the end of the filename.

        Returns:
            str: Absolute path to the file.
        """
        name_pattern = f"{filename}{suffix}*" if filename != "POTCAR" else f"{filename}*"
        paths = glob(os.path.join(path, name_pattern))
        fpath = None
        if len(paths) >= 1:
            # using reverse=True because, if multiple files are present,
            # they likely have suffixes 'static', 'relax', 'relax2', etc.
            # and this would give 'static' over 'relax2' over 'relax'
            # however, better to use 'suffix' kwarg to avoid this!
            paths.sort(reverse=True)
            if len(paths) > 1:
                warnings.warn(f"Multiple files detected, using {os.path.basename(paths[0])}")
            fpath = paths[0]
        return fpath

    def _execute_chargemol(self, **job_control_kwargs):
        """Internal function to run Chargemol.

        Args:
            atomic_densities_path (str): Path to the atomic densities directory
            required by Chargemol. If None, Pymatgen assumes that this is
            defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment variable.
                Default: None.
            job_control_kwargs: Keyword arguments for _write_jobscript_for_chargemol.
        """
        with ScratchDir("."):
            try:
                os.symlink(self._chgcar_path, "./CHGCAR")
                os.symlink(self._potcar_path, "./POTCAR")
                os.symlink(self._aeccar0_path, "./AECCAR0")
                os.symlink(self._aeccar2_path, "./AECCAR2")
            except OSError as exc:
                print(f"Error creating symbolic link: {exc}")

            # write job_script file:
            self._write_jobscript_for_chargemol(**job_control_kwargs)

            # Run Chargemol
            with subprocess.Popen(
                CHARGEMOL_EXE,
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                close_fds=True,
            ) as rs:
                _stdout, stderr = rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError(
                    f"{CHARGEMOL_EXE} exit code: {rs.returncode}, error message: {stderr!s}. "
                    "Please check your Chargemol installation."
                )

            self._from_data_dir()

    def _from_data_dir(self, chargemol_output_path=None):
        """Internal command to parse Chargemol files from a directory.

        Args:
            chargemol_output_path (str): Path to the folder containing the
            Chargemol output files.
                Default: None (current working directory).
        """
        if chargemol_output_path is None:
            chargemol_output_path = "."

        charge_path = f"{chargemol_output_path}/DDEC6_even_tempered_net_atomic_charges.xyz"
        self.ddec_charges = self._get_data_from_xyz(charge_path)
        self.dipoles = self._get_dipole_info(charge_path)

        bond_order_path = f"{chargemol_output_path}/DDEC6_even_tempered_bond_orders.xyz"
        if os.path.isfile(bond_order_path):
            self.bond_order_sums = self._get_data_from_xyz(bond_order_path)
            self.bond_order_dict = self._get_bond_order_info(bond_order_path)
        else:
            self.bond_order_sums = self.bond_order_dict = None

        spin_moment_path = f"{chargemol_output_path}/DDEC6_even_tempered_atomic_spin_moments.xyz"
        if os.path.isfile(spin_moment_path):
            self.ddec_spin_moments = self._get_data_from_xyz(spin_moment_path)
        else:
            self.ddec_spin_moments = None

        rsquared_path = f"{chargemol_output_path}/DDEC_atomic_Rsquared_moments.xyz"
        if os.path.isfile(rsquared_path):
            self.ddec_rsquared_moments = self._get_data_from_xyz(rsquared_path)
        else:
            self.ddec_rsquared_moments = None

        rcubed_path = f"{chargemol_output_path}/DDEC_atomic_Rcubed_moments.xyz"
        if os.path.isfile(rcubed_path):
            self.ddec_rcubed_moments = self._get_data_from_xyz(rcubed_path)
        else:
            self.ddec_rcubed_moments = None

        rfourth_path = f"{chargemol_output_path}/DDEC_atomic_Rfourth_moments.xyz"
        if os.path.isfile(rfourth_path):
            self.ddec_rfourth_moments = self._get_data_from_xyz(rfourth_path)
        else:
            self.ddec_rfourth_moments = None

        ddec_analysis_path = f"{chargemol_output_path}/VASP_DDEC_analysis.output"
        if os.path.isfile(ddec_analysis_path):
            self.cm5_charges = self._get_cm5_data_from_output(ddec_analysis_path)
        else:
            self.cm5_charges = None

    def get_charge_transfer(self, atom_index, charge_type="ddec"):
        """Get the charge transferred for a particular atom. A positive value means
        that the site has gained electron density (i.e. exhibits anionic character)
        whereas a negative value means the site has lost electron density (i.e. exhibits
        cationic character). This is the same thing as the negative of the partial atomic
        charge.

        Args:
            atom_index (int): Index of atom to get charge transfer for.
            charge_type (str): Type of charge to use ("ddec" or "cm5").

        Returns:
            float: charge transferred at atom_index
        """
        if charge_type.lower() == "ddec":
            return -self.ddec_charges[atom_index]

        if charge_type.lower() == "cm5":
            return -self.cm5_charges[atom_index]

        raise ValueError(f"Invalid {charge_type=}")

    def get_charge(self, atom_index, nelect=None, charge_type: Literal["ddec", "cm5"] = "ddec"):
        """Convenience method to get the charge on a particular atom using the same
        sign convention as the BaderAnalysis. Note that this is *not* the partial
        atomic charge. This value is nelect (e.g. ZVAL from the POTCAR) + the
        charge transferred. If you want the partial atomic charge, use
        get_partial_charge().

        Args:
            atom_index (int): Index of atom to get charge for.
            nelect (int): number of electrons associated with an isolated atom at this index.
            For most DFT codes this corresponds to the number of valence electrons
            associated with the pseudopotential. If None, this value will be automatically
            obtained from the POTCAR (if present).
                Default: None.
            charge_type (str): Type of charge to use ("ddec" or "cm5").

        Returns:
            float: charge on atom_index
        """
        if nelect:
            charge = nelect + self.get_charge_transfer(atom_index, charge_type=charge_type)
        elif self.potcar and self.natoms:
            charge = None
            potcar_indices = []
            for idx, val in enumerate(self.natoms):
                potcar_indices += [idx] * val
            nelect = self.potcar[potcar_indices[atom_index]].nelectrons
            charge = nelect + self.get_charge_transfer(atom_index, charge_type=charge_type)
        else:
            charge = None
        return charge

    def get_partial_charge(self, atom_index, charge_type: Literal["ddec", "cm5"] = "ddec"):
        """Convenience method to get the partial atomic charge on a particular atom.
        This is the value printed in the Chargemol analysis.

        Args:
            atom_index (int): Index of atom to get charge for.
            charge_type (str): Type of charge to use ("ddec" or "cm5").
        """
        if charge_type.lower() == "ddec":
            return self.ddec_charges[atom_index]

        if charge_type.lower() == "cm5":
            return self.cm5_charges[atom_index]

        raise ValueError(f"Invalid charge_type: {charge_type}")

    def get_bond_order(self, index_from, index_to):
        """Convenience method to get the bond order between two atoms.

        Args:
            index_from (int): Index of atom to get bond order from.
            index_to (int): Index of atom to get bond order to.

        Returns:
            float: bond order between atoms
        """
        bonded_set = self.bond_order_dict[index_from]["bonded_to"]
        bond_orders = [v["bond_order"] for v in bonded_set if v["index"] == index_to]
        return 0.0 if bond_orders == [] else np.sum(bond_orders)

    def _write_jobscript_for_chargemol(
        self,
        net_charge=0.0,
        periodicity=(True, True, True),
        method="ddec6",
        compute_bond_orders=True,
    ):
        """Write job_script.txt for Chargemol execution.

        Args:
            net_charge (float): Net charge of the system.
                Defaults to 0.0.
            periodicity (tuple[bool]): Periodicity of the system.
                Default: (True, True, True).
            method (str): Method to use for the analysis. Options include "ddec6"
                and "ddec3". Default: "ddec6"
            compute_bond_orders (bool): Whether to compute bond orders. Default: True.
        """
        self.net_charge = net_charge
        self.periodicity = periodicity
        self.method = method

        lines = ""

        # Net Charge
        if net_charge:
            lines += f"<net charge>\n{net_charge}\n</net charge>\n"

        # Periodicity
        if periodicity:
            per_a = ".true." if periodicity[0] else ".false."
            per_b = ".true." if periodicity[1] else ".false."
            per_c = ".true." if periodicity[2] else ".false."
            lines += (
                f"<periodicity along A, B, and C vectors>\n{per_a}\n{per_b}\n{per_c}\n"
                "</periodicity along A, B, and C vectors>\n"
            )

        # atomic_densities dir
        atomic_densities_path = self._atomic_densities_path or os.getenv("DDEC6_ATOMIC_DENSITIES_DIR")
        if atomic_densities_path is None:
            raise OSError(
                "The DDEC6_ATOMIC_DENSITIES_DIR environment variable must be set or the atomic_densities_path must"
                " be specified"
            )
        if not os.path.isfile(atomic_densities_path):
            raise FileNotFoundError(f"{atomic_densities_path=} does not exist")

        # This is to fix a Chargemol filepath nuance
        if platform.system() == "Windows":
            if atomic_densities_path[-1] != "\\":
                atomic_densities_path += "\\"
        elif atomic_densities_path[-1] != "/":
            atomic_densities_path += "/"

        lines += (
            f"\n<atomic densities directory complete path>\n{atomic_densities_path}\n</atomic densities directory "
            "complete path>\n"
        )

        # Charge type
        lines += f"\n<charge type>\n{method.upper()}\n</charge type>\n"

        if compute_bond_orders:
            bo = ".true." if compute_bond_orders else ".false."
            lines += f"\n<compute BOs>\n{bo}\n</compute BOs>\n"

        with open("job_control.txt", mode="w") as file:
            file.write(lines)

    @staticmethod
    def _get_dipole_info(filepath):
        """Internal command to process dipoles.

        Args:
            filepath (str): The path to the DDEC6_even_tempered_net_atomic_charges.xyz file
        """
        idx = 0
        start = False
        dipoles = []
        with open(filepath) as file:
            for line in file:
                if "The following XYZ" in line:
                    start = True
                    idx += 1
                    continue
                if start and line.strip() == "":
                    break
                if idx >= 2:
                    dipoles.append(list(map(float, line.strip().split()[6:9])))
                if start:
                    idx += 1

        return dipoles

    @staticmethod
    def _get_bond_order_info(filename):
        """Internal command to process pairwise bond order information.

        Args:
            filename (str): The path to the DDEC6_even_tempered_bond_orders.xyz file
        """
        # Get where relevant info for each atom starts
        bond_order_info = {}

        with open(filename, encoding="utf-8") as file:
            start_idx = 0
            for line in file:
                split = line.strip().split()
                if "Printing BOs" in line:
                    start_idx = int(split[5]) - 1
                    start_el = Element(split[7])
                    bond_order_info[start_idx] = {"element": start_el, "bonded_to": []}

                elif "Bonded to the" in line:
                    direction = tuple(int(i.split(")")[0].split(",")[0]) for i in split[4:7])
                    end_idx = int(split[12]) - 1
                    end_el = Element(split[14])
                    bo = float(split[20])
                    spin_bo = float(split[-1])
                    bonded_to = {
                        "index": end_idx,
                        "element": end_el,
                        "bond_order": bo,
                        "direction": direction,
                        "spin_polarization": spin_bo,
                    }
                    bond_order_info[start_idx]["bonded_to"].append(bonded_to)

                elif "The sum of bond orders for this atom" in line:
                    bond_order_info[start_idx]["bond_order_sum"] = float(split[-1])

        return bond_order_info

    def get_property_decorated_structure(self):
        """Takes CHGCAR's structure object and updates it with properties
        from the Chargemol analysis.

        Returns:
            Pymatgen structure with site properties added
        """
        struct = self.structure.copy()
        struct.add_site_property("partial_charge_ddec6", self.ddec_charges)
        if self.dipoles:
            struct.add_site_property("dipole_ddec6", self.dipoles)
        if self.bond_order_sums:
            struct.add_site_property("bond_order_sum_ddec6", self.bond_order_sums)
        if self.ddec_spin_moments:
            struct.add_site_property("spin_moment_ddec6", self.ddec_spin_moments)
        if self.cm5_charges:
            struct.add_site_property("partial_charge_cm5", self.cm5_charges)
        return struct

    @property
    def summary(self):
        """A dictionary summary of the Chargemol analysis
        {
            "ddec": {
                "partial_charges": list[float],
                "spin_moments": list[float],
                "dipoles": list[float],
                "rsquared_moments": list[float],
                "rcubed_moments": list[float],
                "rfourth_moments": list[float],
                "bond_order_dict": dict
            },
            "cm5": {
                "partial_charges": list[float],
            }
        }.
        """
        summary = {}
        ddec_summary = {"partial_charges": self.ddec_charges}
        if self.bond_order_sums:
            ddec_summary["bond_order_sums"] = self.bond_order_sums
        if self.ddec_spin_moments:
            ddec_summary["spin_moments"] = self.ddec_spin_moments
        if self.dipoles:
            ddec_summary["dipoles"] = self.dipoles
        if self.ddec_rsquared_moments:
            ddec_summary["rsquared_moments"] = self.ddec_rsquared_moments
        if self.ddec_rcubed_moments:
            ddec_summary["rcubed_moments"] = self.ddec_rcubed_moments
        if self.ddec_rfourth_moments:
            ddec_summary["rfourth_moments"] = self.ddec_rfourth_moments
        if self.bond_order_dict:
            ddec_summary["bond_order_dict"] = self.bond_order_dict

        cm5_summary = {"partial_charges": self.cm5_charges} if self.cm5_charges else None

        summary["ddec"] = ddec_summary
        summary["cm5"] = cm5_summary

        return summary

    @staticmethod
    def _get_data_from_xyz(xyz_path) -> list[float]:
        """Internal command to process Chargemol XYZ files.

        Args:
            xyz_path (str): Path to XYZ file

        Returns:
            list[float]: site-specific properties
        """
        props = []
        if os.path.isfile(xyz_path):
            with open(xyz_path) as file:
                for idx, line in enumerate(file):
                    if idx <= 1:
                        continue
                    if line.strip() == "":
                        break
                    props.append(float(line.split()[-1]))
        else:
            raise FileNotFoundError(f"{xyz_path} not found")

        return props

    @staticmethod
    def _get_cm5_data_from_output(ddec_analysis_path) -> list[float]:
        """Internal command to process Chargemol CM5 data.

        Args:
            ddec_analysis_path (str): Path VASP_DDEC_analysis.output file

        Returns:
            list[float]: CM5 charges
        """
        props = []
        if os.path.isfile(ddec_analysis_path):
            start = False
            with open(ddec_analysis_path) as file:
                for line in file:
                    if "computed CM5" in line:
                        start = True
                        continue
                    if "Hirshfeld and CM5" in line:
                        break
                    if start:
                        vals = line.split()
                        props.extend([float(c) for c in [val.strip() for val in vals]])
        else:
            raise FileNotFoundError(f"{ddec_analysis_path} not found")
        return props
