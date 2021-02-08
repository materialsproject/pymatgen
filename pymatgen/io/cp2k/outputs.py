# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines the Cp2k output parser along with a few other functions for parsing cp2k-related
outputs.
"""

import glob
import logging
import os
import re
import warnings

import numpy as np
import pandas as pd
from monty.io import zopen
from monty.json import jsanitize
from monty.re import regrep

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.electronic_structure.dos import CompleteDos, Dos, add_densities
from pymatgen.io.cp2k.sets import Cp2kInput
from pymatgen.io.cp2k.utils import _postprocessor, natural_keys
from pymatgen.io.xyz import XYZ

__author__ = "Nicholas Winner"
__version__ = "0.3"
__status__ = "Development"

logger = logging.getLogger(__name__)

_hartree_to_ev_ = 2.72113838565563e01
_static_run_names_ = [
    "ENERGY",
    "ENERGY_FORCE",
    "WAVEFUNCTION_OPTIMIZATION",
    "WFN_OPT",
]


class Cp2kOutput:
    """
    Class for parsing output file from CP2K. The CP2K output file is very flexible in the way that it is returned.
    This class will automatically parse parameters that should always be present, but other parsing features may be
    called depending on the run type.
    """

    def __init__(self, filename, verbose=False, auto_load=False):
        """
        Initialize the Cp2kOutput object.

        Args:
            filename: (str) Name of the CP2K output file to parse
            verbose: (bool) Whether or not to parse with verbosity (will parse lots of data that may not be useful)
            auto_load (bool): Whether or not to automatically load basic info like energies and structures.
        """

        # IO Info
        self.filename = filename
        self.dir = os.path.dirname(filename)
        self.filenames = {}
        self.parse_files()
        self.data = {}

        # Material properties/results
        self.input = None
        self.initial_structure = None
        self.lattice = None
        self.final_structure = None
        self.composition = None
        self.efermi = None
        self.vbm = None
        self.cbm = None
        self.band_gap = None
        self.structures = []
        self.ionic_steps = []

        # parse the basic run parameters always
        self.parse_cp2k_params()
        self.parse_input()  # parse the input file
        self.parse_global_params()  # Always present, parse the global parameters, most important is what run type
        self.parse_dft_params()  # Present so long as a DFT calculation was performed
        self.parse_scf_params()
        self.parse_atomic_kind_info()

        # Auto-load will load the most crucial data into the data attribute
        if auto_load:
            self.ran_successfully()  # Only if job completed. No info about convergence etc.
            self.convergence()  # Checks to see if job converged

            self.parse_initial_structure()  # Get the initial structure by parsing lattice and then parsing coords
            self.parse_structures()  # collect all structures from the run

            self.parse_energies()  # get total energy for each ionic step
            self.parse_forces()  # get forces on all atoms (in order), if available
            self.parse_stresses()  # get stress tensor and total stress at each ionic step, if available
            self.parse_ionic_steps()  # collect energy, forces, and total stress into ionic steps variable
            self.parse_dos()

            self.parse_mo_eigenvalues()  # Get the eigenvalues of the MOs (for finding gaps, VBM, CBM)
            self.parse_homo_lumo()  # Get the HOMO LUMO gap as printed after the mo eigenvalues
            self.parse_timing()  # Get timing info (includes total CPU time consumed, but also much more)

            # TODO: Is this the best way to implement? Should there just be the option to select each individually?
            if verbose:
                self.parse_scf_opt()
                self.parse_opt_steps()
                self.parse_total_numbers()
                self.parse_mulliken()
                self.parse_hirshfeld()

    @property
    def cp2k_version(self):
        """
        The cp2k version used in the calculation
        """
        return self.data.get("cp2k_version", None)

    @property
    def completed(self):
        """
        Did the calculation complete
        """
        c = self.data.get("completed", False)
        if c:
            return c[0][0]
        return c

    @property
    def num_warnings(self):
        """
        How many warnings showed up during the run
        """
        return self.data.get("num_warnings", 0)

    @property
    def run_type(self):
        """
        What type of run (Energy, MD, etc.) was performed
        """
        return self.data.get("global").get("Run_type")

    @property
    def project_name(self):
        """
        What project name was used for this calculation
        """
        return self.data.get("global").get("project_name")

    @property
    def spin_polarized(self):
        """
        Was the calculation spin polarized
        """
        if ("UKS" or "UNRESTRICTED_KOHN_SHAM" or "LSD" or "SPIN_POLARIZED") in self.data["dft"].values():
            return True
        return False

    @property
    def is_metal(self):
        """
        Was a band gap found? i.e. is it a metal
        """
        if self.band_gap is None:
            return True
        if self.band_gap <= 0:
            return True
        return False

    def parse_files(self):
        """
        Identify files present in the directory with the cp2k output file. Looks for trajectories, dos, and cubes
        """
        pdos = glob.glob(os.path.join(self.dir, "*pdos*"))
        self.filenames["PDOS"] = []
        self.filenames["LDOS"] = []
        for p in pdos:
            if p.split("/")[-1].__contains__("list"):
                self.filenames["LDOS"].append(p)
            else:
                self.filenames["PDOS"].append(p)

        self.filenames["trajectory"] = glob.glob(os.path.join(self.dir, "*pos*.xyz*"))
        self.filenames["forces"] = glob.glob(os.path.join(self.dir, "*frc*.xyz*"))
        self.filenames["stress"] = glob.glob(os.path.join(self.dir, "*stress*"))
        self.filenames["cell"] = glob.glob(os.path.join(self.dir, "*.cell*"))
        self.filenames["electron_density"] = glob.glob(os.path.join(self.dir, "*ELECTRON_DENSITY*.cube*"))
        self.filenames["spin_density"] = glob.glob(os.path.join(self.dir, "*SPIN_DENSITY*.cube*"))
        self.filenames["v_hartree"] = glob.glob(os.path.join(self.dir, "*hartree*.cube*"))
        self.filenames["v_hartree"].sort(key=natural_keys)

        restart = glob.glob(os.path.join(self.dir, "*restart*"))
        self.filenames["restart.bak"] = []
        for r in restart:
            if r.split("/")[-1].__contains__("bak"):
                self.filenames["restart.bak"].append(r)
            else:
                self.filenames["restart"] = r
        wfn = glob.glob(os.path.join(self.dir, "*wfn*"))
        self.filenames["wfn.bak"] = []
        for w in wfn:
            if w.split("/")[-1].__contains__("bak"):
                self.filenames["wfn.bak"].append(w)
            else:
                self.filenames["wfn"] = w

    def parse_structures(self, trajectory_file=None, lattice_file=None):
        """
        Parses the structures from a cp2k calculation. Static calculations simply use the initial structure.
        For calculations with ionic motion, the function will look for the appropriate trajectory and lattice
        files based on naming convention. If no file is given, and no file is found, it is assumed
        that the lattice/structure remained constant, and the initial lattice/structure is used.
        Cp2k does not output the trajectory in the main output file by default, so non static calculations have to
        reference the trajectory file.
        """
        if lattice_file is None:
            if len(self.filenames["cell"]) == 0:
                lattice = self.parse_cell_params()
            elif len(self.filenames["cell"]) == 1:
                latfile = np.loadtxt(self.filenames["cell"][0])
                lattice = (
                    [l[2:11].reshape(3, 3) for l in latfile] if len(latfile.shape) > 1 else latfile[2:11].reshape(3, 3)
                )
                lattice.append(lattice[-1])  # TODO is this always needed? from re-eval at minimum
            else:
                raise FileNotFoundError("Unable to automatically determine lattice file. More than one exist.")
        else:
            latfile = np.loadtxt(lattice_file)
            lattice = [l[2:].reshape(3, 3) for l in latfile]

        if trajectory_file is None:
            if len(self.filenames["trajectory"]) == 0:
                self.structures = []
                self.structures.append(self.parse_initial_structure())
                self.final_structure = self.structures[-1]
            elif len(self.filenames["trajectory"]) == 1:
                mols = XYZ.from_file(self.filenames["trajectory"][0]).all_molecules
                self.structures = []
                for m, l in zip(mols, lattice):
                    self.structures.append(
                        Structure(
                            lattice=l,
                            coords=[s.coords for s in m.sites],
                            species=[s.specie for s in m.sites],
                            coords_are_cartesian=True,
                        )
                    )
                self.final_structure = self.structures[-1]
            else:
                raise FileNotFoundError("Unable to automatically determine trajectory file. More than one exist.")
        else:
            mols = XYZ.from_file(trajectory_file).all_molecules
            self.structures = []
            for m, l in zip(mols, lattice):
                self.structures.append(
                    Structure(
                        lattice=l,
                        coords=[s.coords for s in m.sites],
                        species=[s.specie for s in m.sites],
                        coords_are_cartesian=True,
                    )
                )
            self.final_structure = self.structures[-1]
            self.final_structure.set_charge(self.initial_structure.charge)

    def parse_initial_structure(self):
        """
        Parse the initial structure from the main cp2k output file
        """
        pattern = re.compile(r"- Atoms:\s+(\d+)")
        patterns = {"num_atoms": pattern}
        self.read_pattern(
            patterns=patterns,
            reverse=False,
            terminate_on_match=True,
            postprocess=int,
        )

        coord_table = []
        with zopen(self.filename, "rt") as f:
            while True:
                line = f.readline()
                if "Atom  Kind  Element       X           Y           Z          Z(eff)       Mass" in line:
                    f.readline()
                    for i in range(self.data["num_atoms"][0][0]):
                        coord_table.append(f.readline().split())
                    break

        lattice = self.parse_cell_params()
        gs = {}
        for k in self.data["atomic_kind_info"].values():
            if k["pseudo_potential"].upper() == "NONE":
                gs[k["kind_number"]] = True
            else:
                gs[k["kind_number"]] = False

        self.initial_structure = Structure(
            lattice[0],
            species=[i[2] for i in coord_table],
            coords=[[float(i[4]), float(i[5]), float(i[6])] for i in coord_table],
            coords_are_cartesian=True,
            site_properties={"ghost": [gs.get(int(i[1])) for i in coord_table]},
        )

        self.initial_structure.set_charge(self.input["FORCE_EVAL"]["DFT"].get("CHARGE", [0])[0])
        self.composition = self.initial_structure.composition
        return self.initial_structure

    def ran_successfully(self):
        """
        Sanity checks that the program ran successfully. Looks at the bottom of the CP2K output file
        for the "PROGRAM ENDED" line, which is printed when successfully ran. Also grabs the number
        of warnings issued.
        """
        program_ended_at = re.compile(r"PROGRAM ENDED AT\s+(\w+)")
        num_warnings = re.compile(r"The number of warnings for this run is : (\d+)")
        self.read_pattern(
            patterns={"completed": program_ended_at},
            reverse=True,
            terminate_on_match=True,
            postprocess=bool,
        )
        self.read_pattern(
            patterns={"num_warnings": num_warnings},
            reverse=True,
            terminate_on_match=True,
            postprocess=int,
        )

        if not self.completed:
            raise ValueError("The provided CP2K job did not finish running! Cannot parse the file reliably.")

    def convergence(self):
        """
        Check whether or not the SCF and geometry optimization cycles converged.
        """
        # SCF Loops
        uncoverged_inner_loop = re.compile(r"(Leaving inner SCF loop)")
        scf_converged = re.compile(r"(SCF run converged)|(SCF run NOT converged)")
        self.read_pattern(
            patterns={
                "uncoverged_inner_loop": uncoverged_inner_loop,
                "scf_converged": scf_converged,
            },
            reverse=True,
            terminate_on_match=False,
            postprocess=bool,
        )
        for i, x in enumerate(self.data["scf_converged"]):
            if x[0]:
                self.data["scf_converged"][i] = True
            else:
                self.data["scf_converged"][i] = False

        # GEO_OPT
        geo_opt_not_converged = re.compile(r"(MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED)")
        geo_opt_converged = re.compile(r"(GEOMETRY OPTIMIZATION COMPLETED)")
        self.read_pattern(
            patterns={
                "geo_opt_converged": geo_opt_converged,
                "geo_opt_not_converged": geo_opt_not_converged,
            },
            reverse=True,
            terminate_on_match=True,
            postprocess=bool,
        )

        if not all(self.data["scf_converged"]):
            warnings.warn(
                "There is at least one unconverged SCF cycle in the provided cp2k calculation",
                UserWarning,
            )
        if any(self.data["geo_opt_not_converged"]):
            warnings.warn("Geometry optimization did not converge", UserWarning)

    def parse_energies(self):
        """
        Get the total energy from the output file
        """
        toten_pattern = re.compile(r"Total FORCE_EVAL.*\s(-?\d+.\d+)")
        self.read_pattern(
            {"total_energy": toten_pattern},
            terminate_on_match=False,
            postprocess=float,
            reverse=False,
        )
        self.data["total_energy"] = np.multiply(self.data.get("total_energy", []), _hartree_to_ev_)
        self.final_energy = self.data.get("total_energy", [])[-1][-1]

    def parse_forces(self):
        """
        Get the forces from the output file
        """

        if len(self.filenames["forces"]) == 1:
            self.data["forces"] = [
                [list(atom.coords) for atom in step]
                for step in XYZ.from_file(self.filenames["forces"][0]).all_molecules
            ]
        else:
            header_pattern = r"ATOMIC FORCES.+Z"
            row_pattern = r"\s+\d+\s+\d+\s+\w+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
            footer_pattern = r"SUM OF ATOMIC FORCES"

            self.data["forces"] = self.read_table_pattern(
                header_pattern=header_pattern,
                row_pattern=row_pattern,
                footer_pattern=footer_pattern,
                postprocess=_postprocessor,
                last_one_only=False,
            )

    def parse_stresses(self):
        """
        Get the stresses from the output file.
        """
        if len(self.filenames["stress"]) == 1:
            dat = np.loadtxt(self.filenames["stress"][0], skiprows=1)
            self.data["stress_tensor"] = [[list(d[2:5]), list(d[5:8]), list(d[8:11])] for d in dat]
        else:
            header_pattern = r"STRESS TENSOR.+Z"
            row_pattern = r"\s+\w+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
            footer_pattern = r"^$"

            self.data["stress_tensor"] = self.read_table_pattern(
                header_pattern=header_pattern,
                row_pattern=row_pattern,
                footer_pattern=footer_pattern,
                postprocess=_postprocessor,
                last_one_only=False,
            )

            trace_pattern = re.compile(r"Trace\(stress tensor.+(-?\d+\.\d+E?-?\d+)")
            self.read_pattern(
                {"stress": trace_pattern},
                terminate_on_match=False,
                postprocess=float,
                reverse=False,
            )

    def parse_ionic_steps(self):
        """
        Parse the ionic step info
        """
        self.ionic_steps = []

        # TODO: find a better workaround. Currently when optimization is done there
        # is an extra scf step before the optimization starts causing size difference
        if len(self.structures) + 1 == len(self.data["total_energy"]):
            self.data["total_energy"] = self.data["total_energy"][1:]

        for i in range(len(self.data["total_energy"])):
            self.ionic_steps.append({})
            try:
                self.ionic_steps[i]["E"] = self.data["total_energy"][i][0]
            except (TypeError, IndexError):
                warnings.warn("No total energies identified! Check output file")
            try:
                self.ionic_steps[i]["forces"] = self.data["forces"][i]
            except (TypeError, IndexError):
                pass
            try:
                self.ionic_steps[i]["stress_tensor"] = self.data["stress_tensor"][i][0]
            except (TypeError, IndexError):
                pass
            try:
                self.ionic_steps[i]["structure"] = self.structures[i]
            except (TypeError, IndexError):
                warnings.warn("Structure corresponding to this ionic step was not found!")

    def parse_cp2k_params(self):
        """
        Parse the CP2K general parameters from CP2K output file into a dictionary.
        """
        version = re.compile(r"\s+CP2K\|.+(\d\.\d)")
        input_file = re.compile(r"\s+CP2K\|\s+Input file name\s+(.+)$")
        self.read_pattern(
            {"cp2k_version": version, "input_filename": input_file},
            terminate_on_match=True,
            reverse=False,
            postprocess=_postprocessor,
        )

    def parse_input(self):
        """
        Load in the input set from the input file (if it can be found)
        """
        if len(self.data["input_filename"]) == 0:
            return
        input_filename = self.data["input_filename"][0][0]
        for ext in ["", ".gz", ".GZ", ".z", ".Z", ".bz2", ".BZ2"]:
            if os.path.exists(os.path.join(self.dir, input_filename + ext)):
                self.input = Cp2kInput.from_file(os.path.join(self.dir, input_filename + ext))
                return
        warnings.warn("Original input file not found. Some info may be lost.")

    def parse_global_params(self):
        """
        Parse the GLOBAL section parameters from CP2K output file into a dictionary.
        """
        pat = re.compile(r"\s+GLOBAL\|\s+([\w+\s]*)\s+(\w+)")
        self.read_pattern({"global": pat}, terminate_on_match=False, reverse=False)
        for d in self.data["global"]:
            d[0], d[1] = _postprocessor(d[0]), str(d[1])
        self.data["global"] = dict(self.data["global"])

    def parse_dft_params(self):
        """
        Parse the DFT parameters (as well as functional, HF, vdW params)
        """
        pat = re.compile(r"\s+DFT\|\s+(\w.*)\s\s\s(.*)$")
        self.read_pattern(
            {"dft": pat},
            terminate_on_match=False,
            postprocess=_postprocessor,
            reverse=False,
        )
        self.data["dft"] = dict(self.data["dft"])

        self.data["dft"]["cutoffs"] = {}
        self.data["dft"]["cutoffs"]["density"] = self.data["dft"].pop("Cutoffs:_density", None)
        self.data["dft"]["cutoffs"]["gradient"] = self.data["dft"].pop("gradient", None)
        self.data["dft"]["cutoffs"]["tau"] = self.data["dft"].pop("tau", None)

        # Functional
        functional = re.compile(r"\s+FUNCTIONAL\|\s+(.+):")
        self.read_pattern(
            {"functional": functional},
            terminate_on_match=False,
            postprocess=_postprocessor,
            reverse=False,
        )
        self.data["dft"]["functional"] = [item for sublist in self.data.pop("functional", None) for item in sublist]

        # HF exchange info
        hfx = re.compile(r"\s+HFX_INFO\|\s+(.+):\s+(.*)$")
        self.read_pattern(
            {"hfx": hfx},
            terminate_on_match=False,
            postprocess=_postprocessor,
            reverse=False,
        )
        if len(self.data["hfx"]) > 0:
            self.data["dft"]["hfx"] = dict(self.data.pop("hfx"))

        # Van der waals correction
        vdw = re.compile(r"\s+vdW POTENTIAL\|\s+(DFT-D.)\s")
        self.read_pattern(
            {"vdw": vdw},
            terminate_on_match=False,
            postprocess=_postprocessor,
            reverse=False,
        )
        if len(self.data["vdw"]) > 0:
            self.data["dft"]["vdw"] = self.data.pop("vdw")[0][0]

    def parse_scf_params(self):
        """
        Retrieve the most import SCF parameters: the max number of scf cycles (max_scf),
        the convergence cutoff for scf (eps_scf),
        :return:
        """
        max_scf = re.compile(r"max_scf:\s+(\d+)")
        eps_scf = re.compile(r"eps_scf:\s+(\d+)")
        self.read_pattern(
            {"max_scf": max_scf, "eps_scf": eps_scf},
            terminate_on_match=True,
            reverse=False,
        )
        self.data["scf"] = {}
        self.data["scf"]["max_scf"] = self.data.pop("max_scf")[0][0] if self.data["max_scf"] else None
        self.data["scf"]["eps_scf"] = self.data.pop("eps_scf")[0][0] if self.data["eps_scf"] else None

    def parse_cell_params(self):
        """
        Parse the lattice parameters (initial) from the output file
        """
        cell_volume = re.compile(r"\s+CELL\|\sVolume.*\s(\d+\.\d+)")
        vectors = re.compile(r"\s+CELL\| Vector.*\s(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
        angles = re.compile(r"\s+CELL\| Angle.*\s(\d+\.\d+)")

        self.read_pattern(
            {"cell_volume": cell_volume, "lattice": vectors, "angles": angles},
            terminate_on_match=False,
            postprocess=float,
            reverse=False,
        )
        i = iter(self.data["lattice"])
        return list(zip(i, i, i))

    def parse_atomic_kind_info(self):
        """
        Parse info on what atomic kinds are present and what basis/pseudopotential is describing each of them.
        """
        kinds = re.compile(r"Atomic kind: (\w+)")
        orbital_basis_set = re.compile(r"Orbital Basis Set\s+(.+$)")
        potential_information = re.compile(r"(?:Potential information for\s+(.+$))|(?:atomic kind are GHOST atoms)")
        auxiliary_basis_set = re.compile(r"Auxiliary Fit Basis Set\s+(.+$)")
        core_electrons = re.compile(r"Total number of core electrons\s+(\d+)")
        valence_electrons = re.compile(r"Total number of valence electrons\s+(\d+)")
        pseudo_energy = re.compile(r"Total Pseudopotential Energy.+(-?\d+.\d+)")
        self.read_pattern(
            {
                "kinds": kinds,
                "orbital_basis_set": orbital_basis_set,
                "potential_info": potential_information,
                "auxiliary_basis_set": auxiliary_basis_set,
                "core_electrons": core_electrons,
                "valence_electrons": valence_electrons,
                "pseudo_energy": pseudo_energy,
            },
            terminate_on_match=True,
            postprocess=str,
            reverse=False,
        )
        atomic_kind_info = {}
        for i, kind in enumerate(self.data["kinds"]):
            atomic_kind_info[kind[0]] = {
                "orbital_basis_set": self.data.get("orbital_basis_set")[i][0],
                "pseudo_potential": self.data.get("potential_info")[i][0],
                "kind_number": i + 1,
            }
            try:
                atomic_kind_info[kind[0]]["valence_electrons"] = self.data.get("valence_electrons")[i][0]
            except (TypeError, IndexError):
                atomic_kind_info[kind[0]]["valence_electrons"] = None
            try:
                atomic_kind_info[kind[0]]["core_electrons"] = self.data.get("core_electrons")[i][0]
            except (TypeError, IndexError):
                atomic_kind_info[kind[0]]["core_electrons"] = None
            try:
                atomic_kind_info[kind[0]]["auxiliary_basis_set"] = self.data.get("auxiliary_basis_set")[i]
            except (TypeError, IndexError):
                atomic_kind_info[kind[0]]["auxiliary_basis_set"] = None
            try:
                atomic_kind_info[kind[0]]["total_pseudopotential_energy"] = (
                    self.data.get("total_pseudopotential_energy")[i][0] * _hartree_to_ev_
                )
            except (TypeError, IndexError):
                atomic_kind_info[kind[0]]["total_pseudopotential_energy"] = None
        self.data["atomic_kind_info"] = atomic_kind_info

    def parse_total_numbers(self):
        """
        Parse total numbers (not usually important)
        """
        atomic_kinds = r"- Atomic kinds:\s+(\d+)"
        atoms = r"- Atoms:\s+(\d+)"
        shell_sets = r"- Shell sets:\s+(\d+)"
        shells = r"- Shells:\s+(\d+)"
        primitive_funcs = r"- Primitive Cartesian functions:\s+(\d+)"
        cart_base_funcs = r"- Cartesian basis functions:\s+(\d+)"
        spher_base_funcs = r"- Spherical basis functions:\s+(\d+)"

        self.read_pattern(
            {
                "atomic_kinds": atomic_kinds,
                "atoms": atoms,
                "shell_sets": shell_sets,
                "shells": shells,
                "primitive_cartesian_functions": primitive_funcs,
                "cartesian_basis_functions": cart_base_funcs,
                "spherical_basis_functions": spher_base_funcs,
            },
            terminate_on_match=True,
        )

    def parse_scf_opt(self):
        """
        Parse the SCF cycles (not usually important)
        """
        header = r"Step\s+Update method\s+Time\s+Convergence\s+Total energy\s+Change" + r"\s+\-+"
        row = (
            r"(\d+)\s+(\S+\s?\S+)\s+(\d+\.\d+E\+\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)?"
            + r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+E[\+\-]?\d+)"
        )
        footer = r"^$"

        scfs = self.read_table_pattern(
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer,
            last_one_only=False,
        )

        self.data["electronic_steps"] = []
        self.data["convergence"] = []
        self.data["scf_time"] = []
        for i in scfs:
            self.data["scf_time"].append([float(j[-4]) for j in i])
            self.data["convergence"].append([float(j[-3]) for j in i if j[-3] != "None"])
            self.data["electronic_steps"].append([float(j[-2]) for j in i])

    def parse_timing(self):
        """
        Parse the timing info (how long did the run take).
        """
        header = (
            r"SUBROUTINE\s+CALLS\s+ASD\s+SELF TIME\s+TOTAL TIME" + r"\s+MAXIMUM\s+AVERAGE\s+MAXIMUM\s+AVERAGE\s+MAXIMUM"
        )
        row = r"(\w+)\s+(.+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
        footer = r"\-+"

        timing = self.read_table_pattern(
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer,
            last_one_only=True,
            postprocess=_postprocessor,
        )
        self.timing = {}
        for t in timing:
            self.timing[t[0]] = {
                "calls": {"max": t[1]},
                "asd": t[2],
                "self_time": {"average": t[3], "maximum": t[4]},
                "total_time": {"average": t[5], "maximum": t[6]},
            }

    def parse_opt_steps(self):
        """
        Parse the geometry optimization information
        """
        # "Informations at step =" Summary block (floating point terms)
        total_energy = re.compile(r"\s+Total Energy\s+=\s+(-?\d+.\d+)")
        real_energy_change = re.compile(r"\s+Real energy change\s+=\s+(-?\d+.\d+)")
        prediced_change_in_energy = re.compile(r"\s+Predicted change in energy\s+=\s+(-?\d+.\d+)")
        scaling_factor = re.compile(r"\s+Scaling factor\s+=\s+(-?\d+.\d+)")
        step_size = re.compile(r"\s+Step size\s+=\s+(-?\d+.\d+)")
        trust_radius = re.compile(r"\s+Trust radius\s+=\s+(-?\d+.\d+)")
        used_time = re.compile(r"\s+Used time\s+=\s+(-?\d+.\d+)")

        # For RUN_TYPE=CELL_OPT
        pressure_deviation = re.compile(r"\s+Pressure Deviation.*=\s+(-?\d+.\d+)")
        pressure_tolerance = re.compile(r"\s+Pressure Tolerance.*=\s+(-?\d+.\d+)")

        self.read_pattern(
            {
                "total_energy": total_energy,
                "real_energy_change": real_energy_change,
                "predicted_change_in_energy": prediced_change_in_energy,
                "scaling_factor": scaling_factor,
                "step_size": step_size,
                "trust_radius": trust_radius,
                "used_time": used_time,
                "pressure_deviation": pressure_deviation,
                "pressure_tolerance": pressure_tolerance,
            },
            terminate_on_match=False,
            postprocess=float,
        )

        # "Informations at step =" Summary block (bool terms)
        decrease_in_energy = re.compile(r"\s+Decrease in energy\s+=\s+(\w+)")
        converged_step_size = re.compile(r"\s+Convergence in step size\s+=\s+(\w+)")
        converged_rms_step = re.compile(r"\s+Convergence in RMS step\s+=\s+(\w+)")
        converged_in_grad = re.compile(r"\s+Conv\. in gradients\s+=\s+(\w+)")
        converged_in_rms_grad = re.compile(r"\s+Conv\. in RMS gradients\s+=\s+(\w+)")
        pressure_converged = re.compile(r"\s+Conv\. for  PRESSURE\s+=\s+(\w+)")

        self.read_pattern(
            {
                "decrease_in_energy": decrease_in_energy,
                "converged_step_size": converged_step_size,
                "converged_rms_step": converged_rms_step,
                "converged_in_grad": converged_in_grad,
                "converged_in_rms_grad": converged_in_rms_grad,
                "pressure_converged": pressure_converged,
            },
            terminate_on_match=False,
            postprocess=_postprocessor,
        )

    def parse_mulliken(self):
        """
        Parse the mulliken population analysis info for each step
        :return:
        """
        header = r"Mulliken Population Analysis.+Net charge"
        pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer = r".+Total charge"

        d = self.read_table_pattern(
            header_pattern=header,
            row_pattern=pattern,
            footer_pattern=footer,
            last_one_only=False,
        )
        if d:
            print("Found data, but not yet implemented!")

    def parse_hirshfeld(self):
        """
        parse the hirshfeld population analysis for each step
        """
        uks = self.spin_polarized
        header = r"Hirshfeld Charges.+Net charge"
        footer = r"^$"

        if not uks:
            pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
            d = self.read_table_pattern(
                header_pattern=header,
                row_pattern=pattern,
                footer_pattern=footer,
                last_one_only=False,
            )
            for i, ionic_step in enumerate(d):
                population = []
                net_charge = []
                for site in ionic_step:
                    population.append(site[4])
                    net_charge.append(site[5])
                hirshfeld = [{"population": population[j], "net_charge": net_charge[j]} for j in range(len(population))]
                self.structures[i].add_site_property("hirshfield", hirshfeld)
        else:
            pattern = (
                r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+"
                + r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
            )
            d = self.read_table_pattern(
                header_pattern=header,
                row_pattern=pattern,
                footer_pattern=footer,
                last_one_only=False,
            )
            for i, ionic_step in enumerate(d):
                population = []
                net_charge = []
                spin_moment = []
                for site in ionic_step:
                    population.append(tuple(site[4:5]))
                    spin_moment.append(site[6])
                    net_charge.append(site[7])
                hirshfeld = [
                    {
                        "population": population[j],
                        "net_charge": net_charge[j],
                        "spin_moment": spin_moment[j],
                    }
                    for j in range(len(population))
                ]
                self.structures[i].add_site_property("hirshfield", hirshfeld)

    def parse_mo_eigenvalues(self):
        """
        Parse the MO eigenvalues from the cp2k output file. Will get the eigenvalues (and band gap)
        at each ionic step (if more than one exist).

        Everything is decomposed by spin channel. If calculation was performed without spin polarization,
        then only Spin.up will be present, which represents the average of up and down.
        """
        eigenvalues = []
        band_gap = []
        efermi = []

        with zopen(self.filename, "rt") as f:
            lines = iter(f.readlines())
            for line in lines:
                try:
                    if line.__contains__(" occupied subspace spin"):
                        eigenvalues.append(
                            {
                                "occupied": {Spin.up: [], Spin.down: []},
                                "unoccupied": {Spin.up: [], Spin.down: []},
                            }
                        )
                        efermi.append({Spin.up: None, Spin.down: None})
                        next(lines)
                        while True:
                            line = next(lines)
                            if line.__contains__("Fermi"):
                                efermi[-1][Spin.up] = float(line.split()[-1])
                                break
                            eigenvalues[-1]["occupied"][Spin.up].extend(
                                [_hartree_to_ev_ * float(l) for l in line.split()]
                            )
                        next(lines)
                        line = next(lines)
                        if line.__contains__(" occupied subspace spin"):
                            next(lines)
                            while True:
                                line = next(lines)
                                if line.__contains__("Fermi"):
                                    efermi[-1][Spin.down] = float(line.split()[-1])
                                    break
                                eigenvalues[-1]["occupied"][Spin.down].extend(
                                    [_hartree_to_ev_ * float(l) for l in line.split()]
                                )
                    if line.__contains__(" unoccupied subspace spin"):
                        next(lines)
                        line = next(lines)
                        while True:
                            if line.__contains__("WARNING : did not converge"):
                                warnings.warn(
                                    "Convergence of eigenvalues for " "unoccupied subspace spin 1 did NOT converge"
                                )
                                next(lines)
                                next(lines)
                                next(lines)
                                line = next(lines)
                                eigenvalues[-1]["unoccupied"][Spin.up].extend(
                                    [_hartree_to_ev_ * float(l) for l in line.split()]
                                )
                                next(lines)
                                line = next(lines)
                                break

                            line = next(lines)
                            if "Eigenvalues" in line or "HOMO" in line:
                                break
                            eigenvalues[-1]["unoccupied"][Spin.up].extend(
                                [_hartree_to_ev_ * float(l) for l in line.split()]
                            )
                        if line.__contains__(" unoccupied subspace spin"):
                            next(lines)
                            line = next(lines)
                            while True:
                                if line.__contains__("WARNING : did not converge"):
                                    warnings.warn(
                                        "Convergence of eigenvalues for " "unoccupied subspace spin 2 did NOT converge"
                                    )
                                    next(lines)
                                    next(lines)
                                    next(lines)
                                    line = next(lines)
                                    eigenvalues[-1]["unoccupied"][Spin.down].extend(
                                        [_hartree_to_ev_ * float(l) for l in line.split()]
                                    )
                                    break

                                line = next(lines)
                                if line.__contains__("HOMO"):
                                    next(lines)
                                    break
                                try:
                                    eigenvalues[-1]["unoccupied"][Spin.down].extend(
                                        [_hartree_to_ev_ * float(l) for l in line.split()]
                                    )
                                except AttributeError:
                                    break

                except ValueError:
                    eigenvalues = [
                        {
                            "occupied": {Spin.up: None, Spin.down: None},
                            "unoccupied": {Spin.up: None, Spin.down: None},
                        }
                    ]
                    warnings.warn("Convergence of eigenvalues  for one or more subspaces did NOT converge")

        self.data["eigenvalues"] = eigenvalues
        self.data["band_gap"] = band_gap

        if len(eigenvalues) == 0:
            warnings.warn("No MO eigenvalues detected.")
            return

        # self.data will always contained the eigenvalues resolved by spin channel. The average vbm, cbm, gap,
        # and fermi are saved as class attributes, as there is (usually) no assymmetry in these values for
        # common materials
        if self.spin_polarized:
            self.data["vbm"] = {
                Spin.up: np.max(eigenvalues[-1]["occupied"][Spin.up]),
                Spin.down: np.max(eigenvalues[-1]["occupied"][Spin.down]),
            }
            self.data["cbm"] = {
                Spin.up: np.min(eigenvalues[-1]["unoccupied"][Spin.up]),
                Spin.down: np.min(eigenvalues[-1]["unoccupied"][Spin.down]),
            }
            self.vbm = (self.data["vbm"][Spin.up] + self.data["vbm"][Spin.down]) / 2
            self.cbm = (self.data["cbm"][Spin.up] + self.data["cbm"][Spin.down]) / 2
            self.efermi = (efermi[-1][Spin.up] + efermi[-1][Spin.down]) / 2
        else:
            self.data["vbm"] = {
                Spin.up: np.max(eigenvalues[-1]["occupied"][Spin.up]),
                Spin.down: None,
            }
            self.data["cbm"] = {
                Spin.up: np.min(eigenvalues[-1]["unoccupied"][Spin.up]),
                Spin.down: None,
            }
            self.vbm = self.data["vbm"][Spin.up]
            self.cbm = self.data["cbm"][Spin.up]
            self.efermi = efermi[-1][Spin.up]

    def parse_homo_lumo(self):
        """
        Find the HOMO - LUMO gap in [eV]. Returns the last value. For gaps/eigenvalues decomposed by
        spin up/spin down channel and over many ionic steps, see parse_mo_eigenvalues()
        """
        pattern = re.compile(r"HOMO.*-.*LUMO.*gap.*\s(-?\d+.\d+)")
        self.read_pattern(
            patterns={"band_gap": pattern},
            reverse=True,
            terminate_on_match=False,
            postprocess=float,
        )
        bg = {Spin.up: [], Spin.down: []}
        for i in range(len(self.data["band_gap"])):
            if self.spin_polarized:
                if i % 2:
                    bg[Spin.up].append(self.data["band_gap"][i][0])
                else:
                    bg[Spin.down].append(self.data["band_gap"][i][0])
            else:
                bg[Spin.up].append(self.data["band_gap"][i][0])
                bg[Spin.down].append(self.data["band_gap"][i][0])
        self.data["band_gap"] = bg
        self.band_gap = (bg[Spin.up][-1] + bg[Spin.down][-1]) / 2 if bg[Spin.up] and bg[Spin.down] else None

    def parse_dos(self, pdos_files=None, ldos_files=None, sigma=0):
        """
        Parse the pdos_ALPHA files created by cp2k, and assimilate them into a CompleteDos object.
        Either provide a list of PDOS file paths, or use glob to find the .pdos_ALPHA extension in
        the calculation directory.

        Args:
            pdos_files (list): list of pdos file paths, otherwise they will be inferred
            ldos_Files (list): list of ldos file paths, otherwise they will be inferred
            sigma (float): Gaussian smearing parameter, if desired. Because cp2k is generally
                used as a gamma-point only code, this is often needed to get smooth DOS that
                are comparable to k-point averaged DOS
        """
        if pdos_files is None:
            pdos_files = self.filenames["PDOS"]

        if ldos_files is None:
            ldos_files = self.filenames["LDOS"]

        # Parse specie projected dos
        tdos, pdoss, ldoss = None, {}, {}
        for pdos_file in pdos_files:
            _pdos, _tdos = parse_dos(pdos_file, total=True, sigma=sigma)
            for k in _pdos:
                if k in pdoss:
                    for orbital in _pdos[k]:
                        pdoss[k][orbital].densities.update(_pdos[k][orbital].densities)
                else:
                    pdoss.update(_pdos)
            if not tdos:
                tdos = _tdos
            else:
                if not all([_tdos.densities.keys() == tdos.densities.keys()]):
                    tdos.densities.update(_tdos.densities)
                else:
                    tdos.densities = add_densities(density1=_tdos.densities, density2=tdos.densities)

        # parse any site-projected dos
        for ldos_file in ldos_files:
            _pdos = parse_dos(ldos_file, sigma=sigma)
            for k in _pdos:
                if k in ldoss:
                    for orbital in _pdos[k]:
                        ldoss[k][orbital].densities.update(_pdos[k][orbital].densities)
                else:
                    ldoss.update(_pdos)

        self.data["pdos"] = jsanitize(pdoss, strict=True)
        self.data["ldos"] = jsanitize(ldoss, strict=True)
        self.data["tdos"] = jsanitize(tdos, strict=True)

        # If number of site-projected dos == number of sites, assume they are bijective
        # and create the CompleteDos object
        _ldoss = {}
        if len(ldoss) == len(self.initial_structure):
            for k in self.data["ldos"]:
                _ldoss[self.initial_structure[int(k) - 1]] = self.data["ldos"][k]
            self.data["cdos"] = CompleteDos(self.final_structure, total_dos=tdos, pdoss=_ldoss)

    @staticmethod
    def _gauss_smear(densities, energies, npts, width):

        if not width:
            return densities

        """Return a gaussian smeared DOS"""
        d = np.zeros(npts)
        e_s = np.linspace(min(energies), max(energies), npts)

        for e, _pd in zip(energies, densities):
            weight = np.exp(-(((e_s - e) / width) ** 2)) / (np.sqrt(np.pi) * width)
            d += _pd * weight

        return d

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False, postprocess=str):
        r"""
        This function originally comes from pymatgen.io.vasp.outputs Outcar class

        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.,
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            reverse (bool): Read files in reverse. Defaults to false. Useful for
                large files, esp OUTCARs, especially when used with
                terminate_on_match.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"} will set the
            value of self.data["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned values
            are lists of lists, because you can grep multiple items on one line.
        """
        matches = regrep(
            self.filename,
            patterns,
            reverse=reverse,
            terminate_on_match=terminate_on_match,
            postprocess=postprocess,
        )
        for k in patterns.keys():
            self.data[k] = [i[0] for i in matches.get(k, [])]

    def read_table_pattern(
        self,
        header_pattern,
        row_pattern,
        footer_pattern,
        postprocess=str,
        attribute_name=None,
        last_one_only=True,
    ):
        r"""
        This function originally comes from pymatgen.io.vasp.outputs Outcar class

        Parse table-like data. A table composes of three parts: header,
        main body, footer. All the data matches "row pattern" in the main body
        will be returned.

        Args:
            header_pattern (str): The regular expression pattern matches the
                table header. This pattern should match all the text
                immediately before the main body of the table. For multiple
                sections table match the text until the section of
                interest. MULTILINE and DOTALL options are enforced, as a
                result, the "." meta-character will also match "\n" in this
                section.
            row_pattern (str): The regular expression matches a single line in
                the table. Capture interested field using regular expression
                groups.
            footer_pattern (str): The regular expression matches the end of the
                table. E.g. a long dash line.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.
            attribute_name (str): Name of this table. If present the parsed data
                will be attached to "data. e.g. self.data["efg"] = [...]
            last_one_only (bool): All the tables will be parsed, if this option
                is set to True, only the last table will be returned. The
                enclosing list will be removed. i.e. Only a single table will
                be returned. Default to be True.

        Returns:
            List of tables. 1) A table is a list of rows. 2) A row if either a list of
            attribute values in case the the capturing group is defined without name in
            row_pattern, or a dict in case that named capturing groups are defined by
            row_pattern.
        """
        with zopen(self.filename, "rt") as f:
            text = f.read()
        table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + row_pattern + r")+)\s+" + footer_pattern
        table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        tables = []
        for mt in table_pattern.finditer(text):
            table_body_text = mt.group("table_body")
            table_contents = []
            for line in table_body_text.split("\n"):
                ml = rp.search(line)
                d = ml.groupdict()
                if len(d) > 0:
                    processed_line = {k: postprocess(v) for k, v in d.items()}
                else:
                    processed_line = [postprocess(v) for v in ml.groups()]
                table_contents.append(processed_line)
            tables.append(table_contents)
        if last_one_only:
            retained_data = tables[-1]
        else:
            retained_data = tables
        if attribute_name is not None:
            self.data[attribute_name] = retained_data
        return retained_data

    def as_dict(self):
        """
        Return dictionary representation of the output
        """
        d = {"input": {}, "output": {}}
        d["total_time"] = self.timing["CP2K"]["total_time"]["maximum"]
        d["run_type"] = self.run_type
        d["input"]["global"] = self.data.get("global")
        d["input"]["dft"] = self.data.get("dft", None)
        d["input"]["scf"] = self.data.get("scf", None)
        d["input"]["structure"] = self.initial_structure.as_dict()
        d["input"]["atomic_kind_info"] = self.data.get("atomic_kind_info", None)
        d["input"]["cp2k_input"] = self.input
        d["ran_successfully"] = self.completed
        d["cp2k_version"] = self.cp2k_version
        d["output"]["structure"] = self.final_structure.as_dict()
        d["output"]["ionic_steps"] = self.ionic_steps
        d["composition"] = self.composition.as_dict()
        d["output"]["energy"] = self.final_energy
        d["output"]["energy_per_atom"] = self.final_energy / self.composition.num_atoms
        d["output"]["bandgap"] = self.band_gap
        d["output"]["cbm"] = self.cbm
        d["output"]["vbm"] = self.vbm
        d["output"]["efermi"] = self.efermi
        d["output"]["is_metal"] = self.is_metal
        return d


def parse_energy_file(energy_file):
    """
    Parses energy file for calculations with multiple ionic steps.
    """
    columns = [
        "step",
        "kinetic_energy",
        "temp",
        "potential_energy",
        "conserved_quantity",
        "used_time",
    ]
    df = pd.read_table(energy_file, skiprows=1, names=columns, sep=r"\s+")
    df["kinetic_energy"] = df["kinetic_energy"] * _hartree_to_ev_
    df["potential_energy"] = df["potential_energy"] * _hartree_to_ev_
    df["conserved_quantity"] = df["conserved_quantity"] * _hartree_to_ev_
    df.astype(float)
    d = {c: df[c].values for c in columns}
    return d


def parse_dos(dos_file=None, spin_channel=None, total=False, sigma=0):
    """
    Parse a single DOS file created by cp2k. Must contain one PDOS snapshot. i.e. you cannot
    use this cannot deal with multiple concatenated dos files.

    Args:
        dos_file (list): list of pdos_ALPHA file paths
        spin_channel (int): Which spin channel the file corresponds to. By default, CP2K will
            write the file with ALPHA or BETA in the filename (for spin up or down), but
            you can specify this here, in case you have a manual file name.
            spin_channel == 1 --> spin up, spin_channel == -1 --> spin down.
        total (bool): Whether to grab the total occupations, or the orbital decomposed ones.
        sigma (float): width for gaussian smearing, if desired

    Returns:
        Everything necessary to create a dos object, in dict format:
            (1) orbital decomposed DOS dict:
                i.e. pdoss = {specie: {orbital.s: {Spin.up: ... }, orbital.px: {Spin.up: ... } ...}}
            (2) energy levels of this dos file
            (3) fermi energy (in eV).
        DOS object is not created here

    """
    if spin_channel:
        spin = Spin(spin_channel)
    else:
        spin = Spin.down if os.path.split(dos_file)[-1].__contains__("BETA") else Spin.up

    with zopen(dos_file, "rt") as f:
        lines = f.readlines()
        kind = re.search(r"atomic kind\s(.*)\sat iter", lines[0]) or re.search(r"list\s(\d+)\s(.*)\sat iter", lines[0])
        kind = kind.groups()[0]

        efermi = float(lines[0].split()[-2]) * _hartree_to_ev_
        header = re.split(r"\s{2,}", lines[1].replace("#", "").strip())[2:]
        dat = np.loadtxt(dos_file)

        def cp2k_to_pmg_labels(x):
            if x == "p":
                return "px"
            if x == "d":
                return "dxy"
            if x == "f":
                return "f_3"
            if x == "d-2":
                return "dxy"
            if x == "d-1":
                return "dyz"
            if x == "d0":
                return "dz2"
            if x == "d+1":
                return "dxz"
            if x == "d+2":
                return "dx2"
            if x == "f-3":
                return "f_3"
            if x == "f-2":
                return "f_2"
            if x == "f-1":
                return "f_1"
            if x == "f0":
                return "f0"
            if x == "f+1":
                return "f1"
            if x == "f+2":
                return "f2"
            if x == "f+3":
                return "f3"
            return x

        header = [cp2k_to_pmg_labels(h) for h in header]

        data = dat[:, 1:]
        data[:, 0] *= _hartree_to_ev_
        energies = data[:, 0] * _hartree_to_ev_
        data = gauss_smear(data, sigma)

        pdos = {
            kind: {
                getattr(Orbital, h): Dos(efermi=efermi, energies=energies, densities={spin: data[:, i + 2]})
                for i, h in enumerate(header)
            }
        }
        if total:
            tdos = Dos(
                efermi=efermi,
                energies=energies,
                densities={spin: np.sum(data[:, 2:], axis=1)},
            )
            return pdos, tdos
        return pdos


def gauss_smear(data, width):
    """Return a gaussian smeared DOS"""

    if not width:
        return data

    npts, nOrbitals = data.shape

    e_s = np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), data.shape[0])
    grid = np.multiply(np.ones((npts, npts)), e_s).T

    def smear(d):
        return np.sum(
            np.multiply(
                np.exp(-((np.subtract(grid, data[:, 0]) / width) ** 2)) / (np.sqrt(np.pi) * width),
                d,
            ),
            axis=1,
        )

    return np.array([smear(data[:, i]) for i in range(1, nOrbitals)]).T
