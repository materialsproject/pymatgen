# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements an interface to the Thomas Manz's
Chargemol code (https://sourceforge.net/projects/ddec/files) for 
calculating DDEC3, DDEC6, and CM5 population analyses.
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
__author__ = "Martin Siron, Andrew S. Rosen"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "01/18/21"

import os
import subprocess
import shutil
import numpy as np
from pymatgen.core import Element

from monty.dev import requires
from monty.os.path import which
from monty.io import zopen
from monty.tempfile import ScratchDir
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Chgcar


CHARGEMOLEXE = (
    which("Chargemol_09_26_2017_linux_parallel") or which("Chargemol_09_26_2017_linux_serial") or which("chargemol")
)

# TODO: Add support for non-VASP output files, including non-zero charge
# TODO: Add support for CM5


class ChargemolAnalysis:
    """
    Chargemol analysis for DDEC3, DDEC6, and/or CM5 population analyses,
    including the calculation of partial atomic charges, atomic spin moments,
    bond orders, and related properties.
    """

    @requires(
        which("Chargemol_09_26_2017_linux_parallel")
        or which("Chargemol_09_26_2017_linux_serial")
        or which("chargemol"),
        "ChargemolAnalysis requires the Chargemol executable to be in the path."
        " Please download the library at https://sourceforge.net/projects/ddec/files"
        "and follow the instructions.",
    )
    def __init__(
        self,
        chgcar_filename=None,
        potcar_filename=None,
        aeccar0_filename=None,
        aeccar2_filename=None,
        atomic_densities_path=None,
        run_chargemol=True,
    ):
        """
        Initializes the Chargemol Analysis.

        Args:
            chgcar_filename (str): Filename of the CHGCAR file
            potcar_filename (str): Filename of the POTCAR file
            aeccar0_filename (str): Filename of the AECCAR0 file
            aeccar2_filename (str): Filename of the AECCAR2 file
            run_chargemol (bool): Whether to run the Chargemol analysis. If False,
                the existing Chargemol output files will be read.
                Default: True.
            atomic_densities_path (str): Path to the atomic densities directory
                required by Chargemol. If None, Pymatgen assumes that this is
                defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment variable.
                Only used if run_chargemol is True.
                Default: None.
        """

        for f in [chgcar_filename, potcar_filename, aeccar0_filename, aeccar2_filename]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"{f} not found")

        self.chgcar = Chgcar.from_file(chgcar_filename)
        self.aeccar0 = Chgcar.from_file(aeccar0_filename)
        self.aeccar2 = Chgcar.from_file(aeccar2_filename)
        self.structure = self.chgcar.structure
        self.potcar = Potcar.from_file(potcar_filename)
        self.natoms = self.chgcar.poscar.natoms

        # List of nelects for each atom from potcar
        potcar_indices = []
        for i, v in enumerate(self.natoms):
            potcar_indices += [i] * v
        self.nelects = (
            [self.potcar[potcar_indices[i]].nelectrons for i in range(len(self.structure))] if self.potcar else []
        )

        if run_chargemol:
            self._execute_chargemol(atomic_densities_path)
        else:
            self._from_data_dir()

    def _execute_chargemol(self, atomic_densities_path, **jobcontrol_kwargs):
        """
        Internal command to execute Chargemol, data analysis runs after
        :param atomic_densities_path: custom directory for atomic_densities. Otherwise
            gets value from "DDEC6_ATOMIC_DENSITIES_DIR" environment variable
        """

        with ScratchDir("."):
            with zopen(self.chgcarpath, "rt") as f_in:
                with open("CHGCAR", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self.potcarpath, "rt") as f_in:
                with open("POTCAR", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self.aeccar0path, "rt") as f_in:
                with open("AECCAR0", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self.aeccar2path, "rt") as f_in:
                with open("AECCAR2", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # write job_script file:
            self._write_jobscript_for_chargemol(atomic_densities_path=atomic_densities_path, **jobcontrol_kwargs)

            # Run Chargemol
            with subprocess.Popen(
                CHARGEMOLEXE,
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                close_fds=True,
            ) as rs:
                stdout, _ = rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError(
                    "Chargemol exited with return code %d. " "Please check your Chargemol installation." % rs.returncode
                )
            try:
                self.version = float(stdout.split()[7])
            except ValueError:
                self.version = -1  # Unknown

            self._from_data_dir()

    def _from_data_dir(self):
        """
        Internal command to load XYZ files and process post Chargemol executable
        """
        if os.path.exists("DDEC6_even_tempered_net_atomic_charges.xyz"):
            with open("DDEC6_even_tempered_net_atomic_charges.xyz", "r") as f:
                atomic_charges_lines = [line.strip() for line in f]
        else:
            raise FileNotFoundError(
                "DDEC6_even_tempered_net_atomic_charges.xyz not found. " "A Chargemol run was not completed."
            )
        if os.path.exists("DDEC6_even_tempered_bond_orders.xyz"):
            with open("DDEC6_even_tempered_bond_orders.xyz", "r") as f:
                bond_orders_lines = [line.strip() for line in f]
        else:
            bond_orders_lines = None
        if os.path.exists("DDEC6_even_tempered_atomic_spin_moments.xyz"):
            with open("DDEC6_even_tempered_atomic_spin_moments.xyz", "r") as f:
                spin_moment_lines = [line.strip() for line in f]
        else:
            spin_moment_lines = None
        if os.path.exists("overlap_populations.xyz"):
            with open("overlap_populations.xyz", "r") as f:
                overlap_populations_lines = [line.strip() for line in f]
        else:
            overlap_populations_lines = None
        if os.path.exists("DDEC_atomic_Rcubed_moments.xyz"):
            with open("DDEC_atomic_Rcubed_moments.xyz", "r") as f:
                Rcubed_moments = [line.strip() for line in f]
        else:
            Rcubed_moments = None
        if os.path.exists("DDEC_atomic_Rfourth_moments.xyz"):
            with open("DDEC_atomic_Rfourth_moments.xyz", "r") as f:
                Rfourth_moments = [line.strip() for line in f]
        else:
            Rfourth_moments = None
        if os.path.exists("DDEC_atomic_Rsquared_moments.xyz"):
            with open("DDEC_atomic_Rsquared_moments.xyz", "r") as f:
                Rsquared_moments = [line.strip() for line in f]
        else:
            Rsquared_moments = None
        if os.path.exists("VASP_DDEC_analysis.output"):
            with open("VASP_DDEC_analysis.output", "r") as f:
                chargemol_output_lines = [line.strip() for line in f]
        else:
            chargemol_output_lines = None

        self.raw_data = {
            "atomic_charges": atomic_charges_lines,
        }
        if spin_moment_lines:
            self.raw_data["atomic_spin_moments"] = spin_moment_lines
        if bond_orders_lines:
            self.raw_data["bond_orders"] = bond_orders_lines
        if overlap_populations_lines:
            self.raw_data["overlap_populations"] = overlap_populations_lines
        if Rfourth_moments:
            self.raw_data["r_fourth"] = Rfourth_moments
        if Rsquared_moments:
            self.raw_data["r_squared"] = Rsquared_moments
        if Rcubed_moments:
            self.raw_data["r_cubed"] = Rcubed_moments
        if chargemol_output_lines:
            self.raw_data["chargemol_output"] = chargemol_output_lines

        self._get_charge_and_dipole_info()
        self._get_bond_order_info()

    def get_charge_transfer(
        self,
        atom_index=None,
        element=None,
    ):
        """
        Get charge for a select index or average of charge for a Element
        :param index: (int) specie index
        :param element: (Element) Pymatgen element
        :return: atomic charge, or avg of atomic charge for an element
        """

        if atom_index is not None:
            return -self.atomic_charges[atom_index]
        elif element is not None:
            charges = []
            for c_element, c_charges in zip(self.species, self.atomic_charges):
                if c_element == element:
                    charges.append(c_charges)
            return np.average(charges)
        else:
            return -self.atomic_charges

    def get_charge(self, atom_index):
        """
        #     Calculates difference between the valence charge of an atomic specie
        #     and its DDEC6 calculated charge
        #     :param index: (int) specie index
        #     :return: charge transfer
        """
        potcar_indices = []
        for i, v in enumerate(self.chgcar.poscar.natoms):
            potcar_indices += [i] * v
        nelect = self.potcar[potcar_indices[atom_index]].nelectrons
        return nelect + self.get_charge_transfer(index=atom_index)

    def get_bond_order(self, index_from, index_to):
        """
        Returns bond order index of species connected to certain specie
        :param index_from: (int) specie originating
        :param index_to: (int) bonded to this specie
        :return: bond order index
        """
        if not self.bond_orders[index_from].get("all_bonds", False):
            return None
        elif not self.bond_orders[index_from].get("all_bonds", {}).get(index_to, False):
            return None
        else:
            return self.bond_orders[index_from].get("all_bonds", {}).get(index_to, {}).get("bond_order", None)

    def _get_info_from_xyz(self, raw_data_key, info_array):
        """
        Internal command to process XYZ files
        :param raw_data_key: key in raw_data collection
        :param info_array: key to analyze in XYZ header
        :return:
        """
        species_count = self.species_count

        # in all files
        all_info = {
            "coords": np.zeros([species_count, 3], dtype="float64"),
            "species": [],
        }

        for element in info_array:
            all_info[element] = np.zeros(species_count, dtype="float64")

        for line_number in range(len(all_info["coords"])):
            line = self.raw_data[raw_data_key][line_number + 2].split()
            all_info["species"].append(Element(line[0]))
            all_info["coords"][line_number][:] = line[1:4]
            for num, element in enumerate(info_array):
                all_info[element][line_number] = line[4 + num]

        return all_info

    def _write_jobscript_for_chargemol(
        self,
        atomic_densities_path=None,
        net_charge=0.0,
        periodicity=[True, True, True],
        charge_type="DDEC6",
        compute_bond_orders=True,
    ):
        """
        Writes job_script.txt for Chargemol execution
        :param atomic_densities_path: (str) atomic densities reference directory
        :param net_charge: (float|None) net charge of structure, 0.0 is default
        :param periodicity: (List[bool]) periodicity among a,b, and c
        :param charge_type: (str) charge type, DDEC6 is default
        """
        self.net_charge = net_charge
        self.periodicity = periodicity
        self.charge_type = charge_type

        lines = ""

        # Net Charge
        if net_charge:
            lines += f"<net charge>\n{net_charge}\n</net charge>\n"

        # Periodicity
        if periodicity:
            per_a = ".true." if periodicity[0] else ".false."
            per_b = ".true." if periodicity[1] else ".false."
            per_c = ".true." if periodicity[2] else ".false."
            lines += f"<periodicity along A, B, and C vectors>\n{per_a}\n{per_b}\n{per_c}\n</periodicity along A, B, and C vectors>\n"

        # atomic_densities dir
        atomic_densities_path = atomic_densities_path or os.environ.get("DDEC6_ATOMIC_DENSITIES_DIR", None)
        if atomic_densities_path is None:
            raise EnvironmentError(
                "The DDEC6_ATOMIC_DENSITIES_DIR environment variable must be set or the atomic_densities_path must be specified"
            )
        if not os.path.exists(atomic_densities_path):
            raise EnvironmentError(f"Cannot find the path to the atomic densities at {atomic_densities_path}")

        # This is to fix a Chargemol filepath nuance
        if os.name == "nt":
            if atomic_densities_path[-1] != "\\":
                atomic_densities_path += "\\"
        else:
            if atomic_densities_path[-1] != "/":
                atomic_densities_path += "/"

        lines += f"\n<atomic densities directory complete path>\n{atomic_densities_path}\n</atomic densities directory complete path>\n"

        # Charge type
        lines += f"\n<charge type>\n{charge_type}\n</charge type>\n"

        if compute_bond_orders:
            bo = ".true." if compute_bond_orders else ".false."
            lines += f"\n<compute BOs>\n{bo}\n</compute BOs>\n"

        with open("job_control.txt", "wt") as fh:
            fh.write(lines)

    def _get_charge_and_dipole_info(self):
        """
        Internal command to process atomic charges, species, and coordinates
        """
        self.species_count = int(self.raw_data["atomic_charges"][0])
        self.atomic_charges = []
        self.species = []
        self.coords = []
        self.dipoles = []
        for line in self.raw_data["atomic_charges"][2 : 2 + self.species_count]:
            self.atomic_charges.append(float(line.split()[-1]))
            self.species.append(Element(line.split()[0]))
            self.coords.append([float(c) for c in line.split()[1:-1]])

        n = [n for n, data in enumerate(self.raw_data["atomic_charges"]) if "The following XYZ" in data][0]
        dipole_data = [s.split()[6:9] for s in self.raw_data["atomic_charges"][n + 2 : n + 2 + self.species_count]]
        self.dipoles = dipole_data

    def _get_bond_order_info(self):
        """
        Internal command to process bond order information
        """
        # Get meta data
        # bo_xyz = self._get_info_from_xyz("bond_orders", ["bond_orders"])

        # Get where relevant info for each atom starts
        bond_order_info = {}
        for line_number, line_content in enumerate(self.raw_data["bond_orders"]):
            if "Printing" in line_content:
                species_index = line_content.split()[5]
                bond_order_info[int(species_index) - 1] = {"start": line_number}

        # combine all relevant info
        for atom in bond_order_info.keys():
            try:
                for bo_line in self.raw_data["bond_orders"][
                    bond_order_info[atom]["start"] + 2 : bond_order_info[atom + 1]["start"] - 4
                ]:

                    # Find total bond order
                    total_bo = float(bo_line.split()[-1])

                    # Find current info
                    c_bonded_to = int(bo_line.split()[12]) - 1
                    c_bonded_to_element = Element(bo_line.split()[14])
                    c_bonded_to_bo = float(bo_line.split()[20])
                    c_direction = (
                        int(bo_line.split()[4][:-1]),
                        int(bo_line.split()[5][:-1]),
                        int(bo_line.split()[6][:-1]),
                    )

                    c_bo_by_bond = {
                        c_bonded_to: {
                            "element": c_bonded_to_element,
                            "bond_order": c_bonded_to_bo,
                            "direction": c_direction,
                        }
                    }
                    bo_by_bond = c_bo_by_bond
                    if bond_order_info[atom].get("all_bonds"):
                        bo_by_bond = bond_order_info[atom].get("all_bonds")
                    bo_by_bond.update(c_bo_by_bond)

                    # update bondings, total_bo
                    bond_order_info[atom].update({"all_bonds": bo_by_bond, "total_bo": total_bo})

            except:
                for bo_line in self.raw_data["bond_orders"][bond_order_info[atom]["start"] + 2 : -3]:
                    # Find total bond order
                    total_bo = float(bo_line.split()[-1])

                    # Find current info
                    c_bonded_to = int(bo_line.split()[12]) - 1
                    c_bonded_to_element = Element(bo_line.split()[14])
                    c_bonded_to_bo = float(bo_line.split()[20])
                    c_direction = (
                        int(bo_line.split()[4][:-1]),
                        int(bo_line.split()[5][:-1]),
                        int(bo_line.split()[6][:-1]),
                    )

                    c_bo_by_bond = {
                        c_bonded_to: {
                            "element": c_bonded_to_element,
                            "bond_order": c_bonded_to_bo,
                            "direction": c_direction,
                        }
                    }
                    bo_by_bond = c_bo_by_bond
                    if bond_order_info[atom].get("all_bonds"):
                        bo_by_bond = bond_order_info[atom].get("all_bonds")
                    bo_by_bond.update(c_bo_by_bond)

                    # update bondings, total_bo
                    bond_order_info[atom].update({"all_bonds": bo_by_bond, "total_bo": total_bo})
        self.bond_orders = bond_order_info

    def get_property_decorated_structure(self):
        """
        Takes CHGCAR's structure object and updates it with atomic charges,
        and bond orders
        :return: updated structure
        """
        struc = self.structure.copy()
        struc.add_site_property("charge_ddec6", self.atomic_charges)
        struc.add_site_property("dipole_ddec6", self.dipoles)
        if self.bond_orders:
            struc.add_site_property("bond_order_sum_ddec6", self.bond_orders)
        if self.spin_moments:
            struc.add_site_property("spin_moment_ddec6", self.spin_moments)
        struc.add_site_property("charge_cm5", self.cm5_charges)
        return struc
