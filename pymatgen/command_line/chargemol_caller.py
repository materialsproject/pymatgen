# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements an interface to Thomas Manz's
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
__author__ = "Martin Siron, Andrew S. Rosen"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "01/18/21"

import os
import subprocess
import shutil
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
# TODO: Handle gzip'd automatically


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
        chgcar_filename="CHGCAR",
        potcar_filename="POTCAR",
        aeccar0_filename="AECCAR0",
        aeccar2_filename="AECCAR2",
        atomic_densities_path=None,
        run_chargemol=True,
        chargemol_output_path=None,
    ):
        """
        Initializes the Chargemol Analysis.

        Args:
            chgcar_filename (str): Filename of the CHGCAR file
            potcar_filename (str): Filename of the POTCAR file
            aeccar0_filename (str): Filename of the AECCAR0 file
            aeccar2_filename (str): Filename of the AECCAR2 file
            atomic_densities_path (str): Path to the atomic densities directory
                required by Chargemol. If None, Pymatgen assumes that this is
                defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment variable.
                Only used if run_chargemol is True.
                Default: None.
            run_chargemol (bool): Whether to run the Chargemol analysis. If False,
                the existing Chargemol output files will be read.
                Default: True.
            chargemol_data_path (str): Path to the Chargemol output files if
                run_chargemol is False. Default is the current working directory.
        """

        self._chgcarpath = os.path.abspath(chgcar_filename)
        self._potcarpath = os.path.abspath(potcar_filename)
        self._aeccar0path = os.path.abspath(aeccar0_filename)
        self._aeccar2path = os.path.abspath(aeccar2_filename)

        self.chgcar = Chgcar.from_file(self._chgcarpath)
        self.aeccar0 = Chgcar.from_file(self._aeccar0path)
        self.aeccar2 = Chgcar.from_file(self._aeccar2path)
        self.structure = self.chgcar.structure
        self.potcar = Potcar.from_file(self._potcarpath)
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
            self._from_data_dir(chargemol_output_path=chargemol_output_path)

    def _execute_chargemol(self, atomic_densities_path=None, **jobcontrol_kwargs):
        """
        Internal function to run Chargemol.

        Args:
            atomic_densities_path (str): Path to the atomic densities directory
            required by Chargemol. If None, Pymatgen assumes that this is
            defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment variable.
                Default: None.
            jobcontrol_kwargs: Keyword arguments for _write_jobscript_for_chargemol.
        """

        with ScratchDir("."):
            with zopen(self._chgcarpath, "rt") as f_in:
                with open("CHGCAR", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self._potcarpath, "rt") as f_in:
                with open("POTCAR", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self._aeccar0path, "rt") as f_in:
                with open("AECCAR0", "wt") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            with zopen(self._aeccar2path, "rt") as f_in:
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
                rs.communicate()
            if rs.returncode != 0:
                raise RuntimeError(
                    "Chargemol exited with return code %d. " "Please check your Chargemol installation." % rs.returncode
                )

            self._from_data_dir()

    def _from_data_dir(self, chargemol_output_path=None):
        """
        Internal command to parse Chargemol files from a directory.

        Args:
            chargemol_output_path (str): Path to the folder containing the
            Chargemol output files.
                Default: None (current working directory).
        """
        if chargemol_output_path is None:
            chargemol_output_path = ""

        charge_path = os.path.join(chargemol_output_path, "DDEC6_even_tempered_net_atomic_charges.xyz")
        if os.path.exists(charge_path):
            ddec_charges = self._get_data_from_xyz(charge_path)
            self.ddec_charges = ddec_charges
        else:
            raise FileNotFoundError(
                "DDEC6_even_tempered_net_atomic_charges.xyz not found. " "A Chargemol run was not completed."
            )

        bond_order_path = os.path.join(chargemol_output_path, "DDEC6_even_tempered_bond_orders.xyz")
        if os.path.exists(bond_order_path):
            self.bond_order_sums = self._get_data_from_xyz(bond_order_path)

        spin_moment_path = os.path.join(chargemol_output_path, "DDEC6_even_tempered_atomic_spin_moments.xyz")
        if os.path.exists(spin_moment_path):
            self.ddec_spin_moments = self._get_data_from_xyz(spin_moment_path)

        rsquared_path = os.path.join(chargemol_output_path, "DDEC_atomic_Rsquared_moments.xyz")
        if os.path.exists(rsquared_path):
            self.ddec_rsquared_moments = self._get_data_from_xyz(rsquared_path)

        rcubed_path = os.path.join(chargemol_output_path, "DDEC_atomic_Rcubed_moments.xyz")
        if os.path.exists(rcubed_path):
            self.ddec_rcubed_moments = self._get_data_from_xyz(rcubed_path)

        rfourth_path = os.path.join(chargemol_output_path, "DDEC_atomic_Rfourth_moments.xyz")
        if os.path.exists(rfourth_path):
            self.ddec_rfourth_moments = self._get_data_from_xyz(rfourth_path)

        ddec_analysis_path = os.path.join(chargemol_output_path, "VASP_DDEC_analysis.output")
        if os.path.exists(ddec_analysis_path):
            self.cm5_charges = self._get_cm5_data_from_output(ddec_analysis_path)

        # self._get_dipole_info()
        # self._get_bond_order_info()

    def get_charge_transfer(self, atom_index, charge_type="ddec"):
        """
        Returns the charge transferred for a particular atom. A positive value means
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
        elif charge_type.lower() == "cm5":
            return -self.cm5_charges[atom_index]
        else:
            raise ValueError("Invalid charge_type: %s" % charge_type)

    def get_charge(self, atom_index, charge_type="ddec"):
        """
        Convenience method to get the charge on a particular atom using the same
        sign convention as the BaderAnalysis. Note that this is *not* the partial
        atomic charge. This value is nelect (e.g. ZVAL from the POTCAR) + the
        charge transferred. If you want the partial atomic charge, use
        get_partial_charge().

        Args:
            atom_index (int): Index of atom to get charge for.
            charge_type (str): Type of charge to use ("ddec" or "cm5").

        Returns:
            float: charge on atom_index
        """
        potcar_indices = []
        for i, v in enumerate(self.chgcar.poscar.natoms):
            potcar_indices += [i] * v
        nelect = self.potcar[potcar_indices[atom_index]].nelectrons
        return nelect + self.get_charge_transfer(atom_index, charge_type=charge_type)

    def get_partial_charge(self, atom_index, charge_type="ddec"):
        """
        Convenience method to get the partial atomic charge on a particular atom.
        This is the value printed in the Chargemol analysis.

        Args:
            atom_index (int): Index of atom to get charge for.
            charge_type (str): Type of charge to use ("ddec" or "cm5").
        """
        if charge_type.lower() == "ddec":
            return self.ddec_charges[atom_index]
        elif charge_type.lower() == "cm5":
            return self.cm5_charges[atom_index]
        else:
            raise ValueError("Invalid charge_type: %s" % charge_type)

    def get_bond_order(self, index_from, index_to):
        """
        Convenience method to get the bond order between two atoms.

        Args:
            index_from (int): Index of atom to get bond order from.
            index_to (int): Index of atom to get bond order to.

        Returns:
            float: bond order between atoms
        """
        if not self.bond_orders[index_from].get("all_bonds", False):
            return None
        elif not self.bond_orders[index_from].get("all_bonds", {}).get(index_to, False):
            return None
        else:
            return self.bond_orders[index_from].get("all_bonds", {}).get(index_to, {}).get("bond_order", None)

    def _get_data_from_xyz(self, xyz_path):
        """
        Internal command to process Chargemol XYZ files

        Args:
            xyz_path (str): Path to XYZ file

        Returns:
            list[float]: site-specific properties
        """

        props = []
        if os.path.exists(xyz_path):
            with open(xyz_path, "r") as r:
                for i, line in enumerate(r):
                    if i <= 1:
                        continue
                    if line.strip() == "":
                        break
                    props.append(float(line.split()[-1]))
        else:
            raise FileNotFoundError(f"{xyz_path} not found")

        return props

    def _get_cm5_data_from_output(self, ddec_analysis_path):
        """
        Internal command to process Chargemol CM5 data

        Args:
            ddec_analysis_path (str): Path VASP_DDEC_analysis.output file

        Returns:
            list[float]: CM5 charges
        """
        props = []
        if os.path.exists(ddec_analysis_path):
            start = False
            with open(ddec_analysis_path, "r") as r:
                for line in r:
                    if "computed CM5" in line:
                        start = True
                        continue
                    elif "Hirshfeld and CM5" in line:
                        break
                    if start:
                        vals = line.split()
                        props.extend([float(c) for c in [val.strip() for val in vals]])
        else:
            raise FileNotFoundError(f"{ddec_analysis_path} not found")
        return props

    def _write_jobscript_for_chargemol(
        self,
        atomic_densities_path=None,
        net_charge=0.0,
        periodicity=[True, True, True],
        method="ddec6",
        compute_bond_orders=True,
    ):
        """
        Writes job_script.txt for Chargemol execution

        Args:
            atomic_densities_path (str): Path to the atomic densities directory
                required by Chargemol. If None, Pymatgen assumes that this is
                defined in a "DDEC6_ATOMIC_DENSITIES_DIR" environment variable.
                Only used if run_chargemol is True.
                Default: None.
            net_charge (float): Net charge of the system.
                Defaults to 0.0.
            periodicity (list[bool]): Periodicity of the system.
                Defaut: [True, True, True].
            method (str): Method to use for the analysis. Options include "ddec6"
            and "ddec3".
                Default: "ddec6"
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
        lines += f"\n<charge type>\n{method.upper()}\n</charge type>\n"

        if compute_bond_orders:
            bo = ".true." if compute_bond_orders else ".false."
            lines += f"\n<compute BOs>\n{bo}\n</compute BOs>\n"

        with open("job_control.txt", "wt") as fh:
            fh.write(lines)

    def _get_dipole_info(self):
        """
        Internal command to process dipoles
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
        dipole_data = [
            [float(d) for d in s.split()[6:9]]
            for s in self.raw_data["atomic_charges"][n + 2 : n + 2 + self.species_count]
        ]
        self.dipoles = dipole_data

    def _get_bond_order_info(self):
        """
        Internal command to process pairwise bond order information
        """
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
        Takes CHGCAR's structure object and updates it with properties
        from the Chargemol analysis.

        Returns
            Pymatgen structure with site properties added
        """
        struc = self.structure.copy()
        struc.add_site_property("partial_charge_ddec6", self.ddec_charges)
        if self.dipoles:
            struc.add_site_property("dipole_ddec6", self.dipoles)
        if self.bond_orders:
            struc.add_site_property("bond_order_sum_ddec6", self.bond_order_sums)
        if self.spin_moments:
            struc.add_site_property("spin_moment_ddec6", self.ddec_spin_moments)
        if self.cm5_charges:
            struc.add_site_property("partial_charge_cm5", self.cm5_charges)
        return struc
