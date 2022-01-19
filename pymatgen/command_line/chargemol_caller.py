# coding: utf-8
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
__author__ = "mhsiron, arosen93"
__version__ = "0.1"
__maintainer__ = "Martin Siron"
__email__ = "mhsiron@lbl.gov"
__status__ = "Beta"
__date__ = "11/02/19"

import os
import subprocess
import shutil
import numpy as np
from pymatgen.core import Element

from monty.dev import requires
from monty.os.path import which
from monty.io import zopen
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Chgcar

EXECUTABLE_OPTIONS = [
    "Chargemol_09_26_2017_linux_parallel",
    "Chargemol_09_26_2017_linux_serial",
    "Chargemol",
]
CHARGEMOLEXE = "Chargemol"
for exe_option in EXECUTABLE_OPTIONS:
    if which(exe_option):
        CHARGEMOLEXE = exe_option
        break

# TODO: Scratch dir!!!
# TODO: Add support for CM5


class ChargemolAnalysis:
    """
    Chargemol Analysis
    """

    @requires(
        which(CHARGEMOLEXE),
        "ChargemolAnalysis requires the Chargemol executable to be in the path."
        " Please download the library at https://sourceforge.net/projects/ddec/files"
        "and follow the instructions.",
    )
    def __init__(
        self,
        chgcar_filename="CHGCAR",
        potcar_filename="POTCAR",
        aeccar0_filename="AECCAR0",
        aeccar2_filename="AECCAR2",
        run=True,
        atomic_densities_path=None,
        gzipped=True,
    ):
        """
        Initializes the Chargemol Analysis. Either runs Chargemol executable,
        or simply analyzes the files created by Chargemol
        :param chgcar_filename: (str) path to CHGCAR file
        :param potcar_filename: (str) path to POTCAR file
        :param aeccar0_filename: (list) path to AECCAR0 file
        :param aeccar2_filename: (list) path to AECCAR2 file
        :param run: (bool) whether or not to run Chargemol
        :param atomic_densities_path: custom directory for atomic_densities. Otherwise
            gets value from "DDEC6_ATOMIC_DENSITIES_DIR" environment variable
        :param gzipped: whether or not the files are gzipped,
            they will be unzipped if so
        """

        for f in [chgcar_filename, potcar_filename, aeccar0_filename, aeccar2_filename]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"File {f} not found")

        aeccar_filenames = [aeccar0_filename, aeccar2_filename]

        # Set properties
        self.species_count = 0
        self.atomic_charges = []
        self.species = []
        self.coords = []
        self.bond_orders = {}
        self.potcar = None

        self.chgcar = Chgcar.from_file(chgcar_filename)
        self.potcar = Potcar.from_file(potcar_filename)

        # Set paths
        self._chgcarpath = os.path.abspath(chgcar_filename)
        self._potcarpath = os.path.abspath(potcar_filename)
        self._aeccarpaths = [os.path.abspath(aeccar) for aeccar in aeccar_filenames]

        if run:
            self._execute_chargemol(atomic_densities_path, gzipped)
        else:
            self._from_data_dir()

    def _execute_chargemol(self, atomic_densities_path, gzipped, **jobcontrol_kwargs):
        """
        Internal command to execute Chargemol, data analysis runs after
        :param atomic_densities_path: custom directory for atomic_densities. Otherwise
            gets value from "DDEC6_ATOMIC_DENSITIES_DIR" environment variable
        :param gzipped: whether or not the files are gzipped,
            they will be unzipped if so
        """

        # Unzip if needed:
        if gzipped:
            with zopen(self._chgcarpath, "rt") as f_in:
                with open("CHGCAR", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self._potcarpath, "rt") as f_in:
                with open("POTCAR", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            for n, aeccarpath in enumerate(self._aeccarpaths):
                with zopen(aeccarpath, "rt") as f_in:
                    with open("AECCAR" + str(n), "wt") as f_out:
                        shutil.copyfileobj(f_in, f_out)

        # write job_script file:
        self._write_jobscript_for_chargemol(atomic_densities_path=atomic_densities_path, **jobcontrol_kwargs)

        # command arguments for Chargemol
        rs = subprocess.Popen(CHARGEMOLEXE, stdout=subprocess.PIPE, stdin=subprocess.PIPE, close_fds=True)
        rs.communicate()
        if rs.returncode != 0:
            raise RuntimeError(
                "Chargemol exited with return code %d. " "Please check your Chargemol installation." % rs.returncode
            )

        self._from_data_dir()

    def _from_data_dir(self):
        """
        Internal command to load XYZ files and process post Chargemol executable
        """
        with open("DDEC6_even_tempered_net_atomic_charges.xyz", "r") as f:
            atomic_charges_lines = [line.strip() for line in f]

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

        self._get_charge_and_dipole_info()
        self._get_bond_order_info()

    def get_charge_transfer(
        self,
        index=None,
        element=None,
    ):
        """
        Get charge for a select index or average of charge for a Element
        :param index: (int) specie index
        :param element: (Element) Pymatgen element
        :return: atomic charge, or avg of atomic charge for an element
        """

        if index is not None:
            return -self.atomic_charges[index]
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
        all_info = {"coords": np.zeros([species_count, 3], dtype="float64"), "species": []}

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
        atomic_densities_path = atomic_densities_path or os.environ.get("DDEC6_ATOMIC_DENSITIES_DIR", ".")
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
            self.coords.append(line.split()[1:-2])

        n = [n for n, data in enumerate(self.raw_data["atomic_charges"]) if "The " in data][0]
        self.dipoles = [s.split()[6:9] for s in self.raw_data["atomic_charges"][n + 2 : n + 2 + self.species_count]]

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

    def update_structure(self):
        """
        Takes CHGCAR's structure object and updates it with atomic charges,
        and bond orders
        :return: updated structure
        """
        structure = self.chgcar.structure
        structure.add_site_property("atomic_charges_ddec6", self.atomic_charges)
        structure.add_site_property("bond_orders_ddec6", self.bond_orders)
        structure.add_site_property("dipole_ddec6", self.dipoles)
        structure.add_site_property("spin_moments_ddec6", self.dipoles)
        return structure
