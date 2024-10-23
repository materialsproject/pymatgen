"""
This module defines the CP2K output parser along with a few other functions
for parsing CP2K-related outputs.
"""

from __future__ import annotations

import os
import re
import warnings
from glob import glob
from itertools import chain

import numpy as np
import pandas as pd
from monty.io import zopen
from monty.json import MSONable, jsanitize
from monty.re import regrep

from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.units import Ha_to_eV
from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.io.cp2k.inputs import Keyword
from pymatgen.io.cp2k.sets import Cp2kInput
from pymatgen.io.cp2k.utils import natural_keys, postprocessor
from pymatgen.io.xyz import XYZ

__author__ = "Nicholas Winner"
__version__ = "2.0"
__status__ = "Production"


class Cp2kOutput:
    """Parse output file from CP2K. The CP2K output file is very flexible in the way that
    it is returned. This class will automatically parse parameters that should always be present,
    but other parsing features may be called depending on the run type.
    """

    def __init__(self, filename, verbose=False, auto_load=False):
        """
        Initialize the Cp2kOutput object.

        Args:
            filename: (str) Name of the CP2K output file to parse
            verbose: (bool) Whether or not to parse with verbosity (will parse lots of data that
                may not be useful)
            auto_load (bool): Whether or not to automatically load basic info like energies
                and structures.
        """
        # IO Info
        self.filename = filename
        self.dir = os.path.dirname(filename)
        self.filenames: dict = {}
        self.parse_files()
        self.data: dict = {}

        # Material properties/results
        self.input = self.initial_structure = self.lattice = self.final_structure = self.composition = None
        self.efermi = self.vbm = self.cbm = self.band_gap = None
        self.structures: list = []
        self.ionic_steps: list = []

        # parse the basic run parameters always
        self.parse_cp2k_params()
        self.parse_input()
        self.parse_global_params()
        self.parse_atomic_kind_info()
        self.parse_dft_params()
        self.parse_scf_params()

        # Auto-load will load the most crucial data into the data attribute
        if auto_load:
            self.ran_successfully()
            self.convergence()

            self.parse_structures()
            self.parse_energies()
            self.parse_forces()
            self.parse_stresses()
            self.parse_ionic_steps()

            self.parse_dos()
            self.parse_bandstructure()
            if not self.band_gap:
                self.parse_homo_lumo()
            if not self.vbm or not self.cbm:
                self.parse_mo_eigenvalues()
            self.parse_timing()

    @property
    def cp2k_version(self):
        """The CP2K version used in the calculation."""
        return self.data.get("cp2k_version")[0][0]

    @property
    def completed(self):
        """Did the calculation complete."""
        if c := self.data.get("completed", False):
            return c[0][0]
        return c

    @property
    def num_warnings(self):
        """How many warnings showed up during the run."""
        return self.data.get("num_warnings", 0)

    @property
    def run_type(self):
        """What type of run (Energy, MD, etc.) was performed."""
        return self.data.get("global").get("Run_type")

    @property
    def calculation_type(self):
        """The calculation type (what io.vasp.outputs calls run_type)."""
        LDA_TYPES = {
            "LDA",
            "PADE",
            "BECKE88",
            "BECKE88_LR",
            "BECKE88_LR_ADIABATIC",
            "BECKE97",
        }

        GGA_TYPES = {"PBE", "PW92"}

        HYBRID_TYPES = {"BLYP", "B3LYP"}

        METAGGA_TYPES = {
            "TPSS": "TPSS",
            "RTPSS": "revTPSS",
            "M06L": "M06-L",
            "MBJ": "modified Becke-Johnson",
            "SCAN": "SCAN",
            "MS0": "MadeSimple0",
            "MS1": "MadeSimple1",
            "MS2": "MadeSimple2",
        }

        functional = self.data.get("dft", {}).get("functional", [None])
        ip = self.data.get("dft", {}).get("hfx", {}).get("Interaction_Potential")
        frac = self.data.get("dft", {}).get("hfx", {}).get("FRACTION")

        if len(functional) > 1:
            rt = "Mixed: " + ", ".join(functional)
            functional = " ".join(functional)
            if "HYP" in functional or (ip and frac) or (functional in HYBRID_TYPES):
                rt = "Hybrid"
        else:
            functional = functional[0]

            if functional is None:
                rt = "None"
            elif "HYP" in functional or (ip and frac) or (functional) in HYBRID_TYPES:
                rt = "Hybrid"
            elif "MGGA" in functional or functional in METAGGA_TYPES:
                rt = "METAGGA"
            elif "GGA" in functional or functional in GGA_TYPES:
                rt = "GGA"
            elif "LDA" in functional or functional in LDA_TYPES:
                rt = "LDA"
            else:
                rt = "Unknown"

        if self.is_hubbard:
            rt += "+U"
        if self.data.get("dft").get("vdw"):
            rt += "+VDW"

        return rt

    @property
    def project_name(self) -> str:
        """What project name was used for this calculation."""
        return self.data.get("global", {}).get("project_name")

    @property
    def spin_polarized(self) -> bool:
        """Was the calculation spin polarized."""
        keys = ("UKS", "UNRESTRICTED_KOHN_SHAM", "LSD", "SPIN_POLARIZED")
        return any(key in self.data["dft"].values() for key in keys)

    @property
    def charge(self) -> float:
        """Charge from the input file."""
        return self.input["FORCE_EVAL"]["DFT"].get("CHARGE", Keyword("", 0)).values[0]

    @property
    def multiplicity(self) -> int:
        """The spin multiplicity from input file."""
        return self.input["FORCE_EVAL"]["DFT"].get("Multiplicity", Keyword("")).values[0]

    @property
    def is_molecule(self) -> bool:
        """
        True if the CP2K output was generated for a molecule (i.e.
        no periodicity in the cell).
        """
        return self.data.get("poisson_periodicity", [[""]])[0][0].upper() == "NONE"

    @property
    def is_metal(self) -> bool:
        """Was a band gap found? i.e. is it a metal."""
        return True if self.band_gap is None else self.band_gap <= 0

    @property
    def is_hubbard(self) -> bool:
        """True if hubbard +U correction was used."""
        for val in self.data.get("atomic_kind_info", {}).values():
            if val.get("DFT_PLUS_U", {}).get("U_MINUS_J", 0) > 0:
                return True
        return False

    def parse_files(self):
        """
        Identify files present in the directory with the CP2K output file. Looks for trajectories,
        dos, and cubes.
        """
        self.filenames["DOS"] = glob(os.path.join(self.dir, "*.dos*"))
        pdos = glob(os.path.join(self.dir, "*pdos*"))
        self.filenames["PDOS"] = []
        self.filenames["LDOS"] = []
        for p in pdos:
            if "list" in p.split("/")[-1]:
                self.filenames["LDOS"].append(p)
            else:
                self.filenames["PDOS"].append(p)
        self.filenames["band_structure"] = glob(os.path.join(self.dir, "*BAND.bs*"))
        self.filenames["trajectory"] = glob(os.path.join(self.dir, "*pos*.xyz*"))
        self.filenames["forces"] = glob(os.path.join(self.dir, "*frc*.xyz*"))
        self.filenames["stress"] = glob(os.path.join(self.dir, "*stress*"))
        self.filenames["cell"] = glob(os.path.join(self.dir, "*.cell*"))
        self.filenames["ener"] = glob(os.path.join(self.dir, "*.ener*"))
        self.filenames["electron_density"] = glob(os.path.join(self.dir, "*ELECTRON_DENSITY*.cube*"))
        self.filenames["spin_density"] = glob(os.path.join(self.dir, "*SPIN_DENSITY*.cube*"))
        self.filenames["v_hartree"] = glob(os.path.join(self.dir, "*hartree*.cube*"))
        self.filenames["hyperfine_tensor"] = glob(os.path.join(self.dir, "*HYPERFINE*eprhyp*"))
        self.filenames["g_tensor"] = glob(os.path.join(self.dir, "*GTENSOR*data*"))
        self.filenames["spinspin_tensor"] = glob(os.path.join(self.dir, "*K*data*"))
        self.filenames["chi_tensor"] = glob(os.path.join(self.dir, "*CHI*data*"))
        self.filenames["nmr_shift"] = glob(os.path.join(self.dir, "*SHIFT*data*"))
        self.filenames["raman"] = glob(os.path.join(self.dir, "*raman*data*"))
        restart = glob(os.path.join(self.dir, "*restart*"))
        self.filenames["restart.bak"] = []
        self.filenames["restart"] = []
        for r in restart:
            if "bak" in r.split("/")[-1]:
                self.filenames["restart.bak"].append(r)
            else:
                self.filenames["restart"].append(r)

        wfn = glob(os.path.join(self.dir, "*.wfn*")) + glob(os.path.join(self.dir, "*.kp*"))
        self.filenames["wfn.bak"] = []
        for w in wfn:
            if "bak" in w.split("/")[-1]:
                self.filenames["wfn.bak"].append(w)
            else:
                self.filenames["wfn"] = w
        for filename in self.filenames.values():
            if hasattr(filename, "sort"):
                filename.sort(key=natural_keys)

    def parse_structures(self, trajectory_file=None, lattice_file=None):
        """
        Parse the structures from a CP2K calculation. Static calculations simply use the initial
        structure. For calculations with ionic motion, the function will look for the appropriate
        trajectory and lattice files based on naming convention. If no file is given, and no file
        is found, it is assumed that the lattice/structure remained constant, and the initial
        lattice/structure is used. CP2K does not output the trajectory in the main output file by
        default, so non static calculations have to reference the trajectory file.
        """
        self.parse_initial_structure()
        trajectory_file = trajectory_file or self.filenames.get("trajectory")
        if isinstance(trajectory_file, list):
            if len(trajectory_file) == 1:
                trajectory_file = trajectory_file[0]
            elif len(trajectory_file) > 1:
                raise FileNotFoundError("Unable to automatically determine trajectory file. More than one exist.")

        if lattice_file is None:
            if len(self.filenames["cell"]) == 0:
                lattices = self.parse_cell_params()
            elif len(self.filenames["cell"]) == 1:
                latt_file = np.loadtxt(self.filenames["cell"][0])
                lattices = (
                    [latt[2:11].reshape(3, 3) for latt in latt_file]
                    if len(latt_file.shape) > 1
                    else [latt_file[2:11].reshape(3, 3)]
                )
            else:
                raise FileNotFoundError("Unable to automatically determine lattice file. More than one exist.")
        else:
            latt_file = np.loadtxt(lattice_file)
            lattices = [latt[2:].reshape(3, 3) for latt in latt_file]

        if not trajectory_file:
            self.structures = [self.initial_structure]
            self.final_structure = self.structures[-1]
        else:
            mols = XYZ.from_file(trajectory_file).all_molecules
            for mol in mols:
                mol.set_charge_and_spin(charge=self.charge, spin_multiplicity=self.multiplicity)
            self.structures = []
            gs = self.initial_structure.site_properties.get("ghost")
            if not self.is_molecule:
                for mol, latt in zip(mols, lattices, strict=False):
                    self.structures.append(
                        Structure(
                            lattice=latt,
                            coords=[s.coords for s in mol],
                            species=[s.specie for s in mol],
                            coords_are_cartesian=True,
                            site_properties={"ghost": gs} if gs else {},
                            charge=self.charge,
                        )
                    )
            else:
                self.structures = mols
            self.final_structure = self.structures[-1]

    def parse_initial_structure(self):
        """Parse the initial structure from the main CP2K output file."""
        patterns = {"num_atoms": re.compile(r"- Atoms:\s+(\d+)")}
        self.read_pattern(
            patterns=patterns,
            reverse=False,
            terminate_on_match=True,
            postprocess=int,
        )

        coord_table = []
        with zopen(self.filename, mode="rt") as file:
            while True:
                line = file.readline()
                if re.search(r"Atom\s+Kind\s+Element\s+X\s+Y\s+Z\s+Z\(eff\)\s+Mass", line):
                    for _ in range(self.data["num_atoms"][0][0]):
                        line = file.readline().split()
                        if line == []:
                            line = file.readline().split()
                        coord_table.append(line)
                    break

        lattice = self.parse_cell_params()
        ghost_atoms = {}
        self.data["atomic_kind_list"] = []
        for val in self.data["atomic_kind_info"].values():
            ghost_atoms[val["kind_number"]] = val["pseudo_potential"].upper() == "NONE"

        for coord in coord_table:
            for key, val in self.data["atomic_kind_info"].items():
                if int(val["kind_number"]) == int(coord[1]):
                    val["element"] = coord[2]
                    self.data["atomic_kind_list"].append(key)
                    break

        if self.is_molecule:
            self.initial_structure = Molecule(
                species=[coord[2] for coord in coord_table],
                coords=[[float(i[4]), float(i[5]), float(i[6])] for i in coord_table],
                site_properties={"ghost": [ghost_atoms.get(int(i[1])) for i in coord_table]},
                charge=self.charge,
                spin_multiplicity=self.multiplicity,
            )
        else:
            self.initial_structure = Structure(
                lattice,
                species=[coord[2] for coord in coord_table],
                coords=[[float(i[4]), float(i[5]), float(i[6])] for i in coord_table],
                coords_are_cartesian=True,
                site_properties={"ghost": [ghost_atoms.get(int(i[1])) for i in coord_table]},
                charge=self.charge,
            )

        self.composition = self.initial_structure.composition
        return self.initial_structure

    def ran_successfully(self):
        """Sanity checks that the program ran successfully. Looks at the bottom of the CP2K output
        file for the "PROGRAM ENDED" line, which is printed when successfully ran. Also grabs
        the number of warnings issued.
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
        """Check whether or not the SCF and geometry optimization cycles converged."""
        # SCF Loops
        unconverged_inner_loop = re.compile(r"(Leaving inner SCF loop)")
        scf_converged = re.compile(r"(SCF run converged)|(SCF run NOT converged)")
        self.read_pattern(
            patterns={
                "unconverged_inner_loop": unconverged_inner_loop,
                "scf_converged": scf_converged,
            },
            reverse=True,
            terminate_on_match=False,
            postprocess=bool,
        )
        for idx, val in enumerate(self.data["scf_converged"]):
            self.data["scf_converged"][idx] = bool(val[0])

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
                "There is at least one unconverged SCF cycle in the provided CP2K calculation",
                UserWarning,
            )
        if any(self.data["geo_opt_not_converged"]):
            warnings.warn("Geometry optimization did not converge", UserWarning)

    def parse_energies(self):
        """Get the total energy from a CP2K calculation. Presently, the energy reported in the
        trajectory (pos.xyz) file takes precedence over the energy reported in the main output
        file. This is because the trajectory file keeps track of energies in between restarts,
        while the main output file may or may not depending on whether a particular machine
        overwrites or appends it.
        """
        if self.filenames.get("trajectory"):
            toten_pattern = r".*E\s+\=\s+(-?\d+.\d+)"
            matches = regrep(
                self.filenames["trajectory"][-1],
                {"total_energy": toten_pattern},
                postprocess=float,
            )
            self.data["total_energy"] = list(
                chain.from_iterable(np.multiply([i[0] for i in matches.get("total_energy", [[]])], Ha_to_eV))
            )
        else:
            toten_pattern = re.compile(r"Total FORCE_EVAL.*\s(-?\d+.\d+)")
            self.read_pattern(
                {"total_energy": toten_pattern},
                terminate_on_match=False,
                postprocess=float,
                reverse=False,
            )
            self.data["total_energy"] = list(
                chain.from_iterable(np.multiply(self.data.get("total_energy", []), Ha_to_eV))
            )
        self.final_energy = self.data.get("total_energy", [])[-1]

    def parse_forces(self):
        """Get the forces from the forces file, or from the main output file."""
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
                postprocess=postprocessor,
                last_one_only=False,
            )

    def parse_stresses(self):
        """Get the stresses from stress file, or from the main output file."""
        if len(self.filenames["stress"]) == 1:
            dat = np.genfromtxt(self.filenames["stress"][0], skip_header=1)
            dat = [dat] if len(np.shape(dat)) == 1 else dat
            self.data["stress_tensor"] = [[list(d[2:5]), list(d[5:8]), list(d[8:11])] for d in dat]
        else:
            header_pattern = r"STRESS\|\s+x\s+y\s+z"
            row_pattern = (
                r"STRESS\|\s+[?:x|y|z]\s+(-?\d+\.\d+E?[-|\+]?\d+)\s+"
                r"(-?\d+\.\d+E?[-|\+]?\d+)\s+(-?\d+\.\d+E?[-|\+]?\d+).*$"
            )
            footer_pattern = r"^$"
            d = self.read_table_pattern(
                header_pattern=header_pattern,
                row_pattern=row_pattern,
                footer_pattern=footer_pattern,
                postprocess=postprocessor,
                last_one_only=False,
            )

            def chunks(lst, n):
                """Yield successive n-sized chunks from lst."""
                for i in range(0, len(lst), n):
                    if i % 2 == 0:
                        yield lst[i : i + n]

            if d:
                self.data["stress_tensor"] = list(chunks(d[0], 3))

    def parse_ionic_steps(self):
        """Parse the ionic step info. If already parsed, this will just assimilate."""
        if not self.structures:
            self.parse_structures()
        if not self.data.get("total_energy"):
            self.parse_energies()
        if not self.data.get("forces"):
            self.parse_forces()
        if not self.data.get("stress_tensor"):
            self.parse_stresses()

        for i, (structure, energy) in enumerate(zip(self.structures, self.data.get("total_energy"), strict=False)):
            self.ionic_steps.append(
                {
                    "structure": structure,
                    "E": energy,
                    "forces": (self.data["forces"][i] if self.data.get("forces") else None),
                    "stress_tensor": (self.data["stress_tensor"][i] if self.data.get("stress_tensor") else None),
                }
            )

        return self.ionic_steps

    def parse_cp2k_params(self):
        """Parse the CP2K general parameters from CP2K output file into a dictionary."""
        version = re.compile(r"\s+CP2K\|.+version\s+(.+)")
        input_file = re.compile(r"\s+CP2K\|\s+Input file name\s+(.+)$")
        self.read_pattern(
            {"cp2k_version": version, "input_filename": input_file},
            terminate_on_match=True,
            reverse=False,
            postprocess=str,
        )

    def parse_plus_u_params(self):
        """Parse the DFT+U params."""
        method = re.compile(r"\s+DFT\+U\|\s+Method\s+()$")
        self.read_pattern(
            {"dft_plus_u_method": method},
            terminate_on_match=True,
            reverse=False,
            postprocess=postprocessor,
        )

    def parse_input(self):
        """Load in the input set from the input file (if it can be found)."""
        if len(self.data["input_filename"]) == 0:
            return
        input_filename = self.data["input_filename"][0][0]
        for ext in ["", ".gz", ".GZ", ".z", ".Z", ".bz2", ".BZ2"]:
            if os.path.isfile(os.path.join(self.dir, input_filename + ext)):
                self.input = Cp2kInput.from_file(os.path.join(self.dir, input_filename + ext))
                return
        warnings.warn("Original input file not found. Some info may be lost.")

    def parse_global_params(self):
        """Parse the GLOBAL section parameters from CP2K output file into a dictionary."""
        pat = re.compile(r"\s+GLOBAL\|\s+([\w+\s]*)\s+(\w+)")
        self.read_pattern({"global": pat}, terminate_on_match=False, reverse=False)
        for d in self.data["global"]:
            d[0], d[1] = postprocessor(d[0]), str(d[1])
        self.data["global"] = dict(self.data["global"])

    def parse_dft_params(self):
        """Parse the DFT parameters (as well as functional, HF, vdW params)."""
        pat = re.compile(r"\s+DFT\|\s+(\w.*)\s\s\s(.*)$")
        self.read_pattern(
            {"dft": pat},
            terminate_on_match=False,
            postprocess=postprocessor,
            reverse=False,
        )
        self.data["dft"] = dict(self.data["dft"])

        self.data["dft"]["cutoffs"] = {}
        self.data["dft"]["cutoffs"]["density"] = self.data["dft"].pop("Cutoffs:_density", None)
        self.data["dft"]["cutoffs"]["gradient"] = self.data["dft"].pop("gradient", None)
        self.data["dft"]["cutoffs"]["tau"] = self.data["dft"].pop("tau", None)

        # Functional
        if self.input and self.input.check("FORCE_EVAL/DFT/XC/XC_FUNCTIONAL"):
            if xc_funcs := list(self.input["force_eval"]["dft"]["xc"]["xc_functional"].subsections):
                self.data["dft"]["functional"] = xc_funcs
            else:
                for v in self.input["force_eval"]["dft"]["xc"].subsections.values():
                    if v.name.upper() == "XC_FUNCTIONAL":
                        self.data["dft"]["functional"] = v.section_parameters
        else:
            functional = re.compile(r"\s+FUNCTIONAL\|\s+(.+):")
            self.read_pattern(
                {"functional": functional},
                terminate_on_match=False,
                postprocess=postprocessor,
                reverse=False,
            )
            self.data["dft"]["functional"] = [item for sublist in self.data.pop("functional", None) for item in sublist]

        # DFT+U
        self.data["dft"]["dft_plus_u"] = self.is_hubbard

        # HF exchange info
        hfx = re.compile(r"\s+HFX_INFO\|\s+(.+):\s+(.*)$")
        self.read_pattern(
            {"hfx": hfx},
            terminate_on_match=False,
            postprocess=postprocessor,
            reverse=False,
        )
        self.data["dft"]["hfx"] = dict(self.data.pop("hfx"))

        # Van der waals correction
        vdw = re.compile(r"\s+vdW POTENTIAL\|(.+)$")
        self.read_pattern(
            {"vdw": vdw},
            terminate_on_match=False,
            postprocess=lambda x: str(x).strip(),
            reverse=False,
        )
        if self.data.get("vdw"):
            found = False
            suffix = ""
            for ll in self.data.get("vdw"):
                for _possible, _name in zip(
                    ("RVV10", "LMKLL", "DRSLL", "DFT-D3", "DFT-D2"),
                    ("RVV10", "LMKLL", "DRSLL", "D3", "D2"),
                    strict=True,
                ):
                    if _possible in ll[0]:
                        found = _name
                    if "BJ" in ll[0]:
                        suffix = "(BJ)"

            self.data["dft"]["vdw"] = found + suffix if found else self.data.pop("vdw")[0][0]

        poisson_periodic = {"poisson_periodicity": re.compile(r"POISSON\| Periodicity\s+(\w+)")}
        self.read_pattern(poisson_periodic, terminate_on_match=True)

    def parse_qs_params(self):
        """Parse the DFT parameters (as well as functional, HF, vdW params)."""
        pat = re.compile(r"\s+QS\|\s+(\w.*)\s\s\s(.*)$")
        self.read_pattern(
            {"QS": pat},
            terminate_on_match=False,
            postprocess=postprocessor,
            reverse=False,
        )
        self.data["QS"] = dict(self.data["QS"])
        tmp = {}
        i = 1
        for k in list(self.data["QS"]):
            if "grid_level" in str(k) and "Number" not in str(k):
                tmp[i] = self.data["QS"].pop(k)
                i += 1
        self.data["QS"]["Multi_grid_cutoffs_[a.u.]"] = tmp

    def parse_overlap_condition(self):
        """Retrieve the overlap condition number."""
        overlap_condition = re.compile(r"\|A\|\*\|A\^-1\|.+=\s+(-?\d+\.\d+E[+\-]?\d+)\s+Log")
        self.read_pattern(
            {"overlap_condition_number": overlap_condition},
            terminate_on_match=True,
            reverse=False,
            postprocess=postprocessor,
        )

    def parse_scf_params(self):
        """
        Retrieve the most import SCF parameters: the max number of scf cycles (max_scf),
        the convergence cutoff for scf (eps_scf),.
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
        """Parse the lattice parameters (initial) from the output file."""
        if self.input and self.input.check("force_eval/subsys/cell"):
            cell = self.input["force_eval"]["subsys"]["cell"]
            if cell.get("abc"):
                return [
                    [cell["abc"].values[0], 0, 0],
                    [0, cell["abc"].values[1], 0],
                    [0, 0, cell["abc"].values[2]],
                ]
            return [
                list(cell.get("A").values),
                list(cell.get("B").values),
                list(cell.get("C").values),
            ]

        warnings.warn(
            "Input file lost. Reading cell params from summary at top of output. Precision errors may result."
        )
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
        lattices = list(zip(i, i, i, strict=True))
        return lattices[0]

    def parse_atomic_kind_info(self):
        """Parse info on what atomic kinds are present and what basis/pseudopotential is describing
        each of them.
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
        _kinds = []
        for _ in list(chain.from_iterable(self.data["kinds"])):
            if _ not in _kinds:
                _kinds.append(_)
        for i, kind in enumerate(_kinds):
            atomic_kind_info[kind] = {
                "orbital_basis_set": self.data.get("orbital_basis_set")[i][0],
                "pseudo_potential": self.data.get("potential_info")[i][0],
                "kind_number": i + 1,
            }
            try:
                if self.data.get("valence_electrons"):
                    tmp = self.data.get("valence_electrons")[i][0]
                elif self.data.get("potential_info")[i][0].upper() == "NONE":
                    tmp = 0
                else:
                    tmp = self.data.get("potential_info")[i][0].split("q")[-1]
                atomic_kind_info[kind]["valence_electrons"] = int(tmp)
            except (TypeError, IndexError, ValueError):
                atomic_kind_info[kind]["valence_electrons"] = None
            try:
                atomic_kind_info[kind]["core_electrons"] = int(self.data.get("core_electrons")[i][0])
            except (TypeError, IndexError, ValueError):
                atomic_kind_info[kind]["core_electrons"] = None
            try:
                atomic_kind_info[kind]["auxiliary_basis_set"] = self.data.get("auxiliary_basis_set")[i]
            except (TypeError, IndexError):
                atomic_kind_info[kind]["auxiliary_basis_set"] = None
            try:
                atomic_kind_info[kind]["total_pseudopotential_energy"] = float(
                    self.data.get("total_pseudopotential_energy")[i][0] * Ha_to_eV
                )
            except (TypeError, IndexError, ValueError):
                atomic_kind_info[kind]["total_pseudopotential_energy"] = None

        with zopen(self.filename, mode="rt") as file:
            j = -1
            lines = file.readlines()
            for k, line in enumerate(lines):
                if "MOLECULE KIND INFORMATION" in line:
                    break
                if "Atomic kind" in line:
                    j += 1
                if "DFT+U correction" in line:
                    atomic_kind_info[_kinds[j]]["DFT_PLUS_U"] = {
                        "L": int(lines[k + 1].split()[-1]),
                        "U_MINUS_J": float(lines[k + 2].split()[-1]),
                    }
                k += 1

        self.data["atomic_kind_info"] = atomic_kind_info

    def parse_total_numbers(self):
        """Parse total numbers (not usually important)."""
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
        """Parse the SCF cycles (not usually important)."""
        header = r"Step\s+Update method\s+Time\s+Convergence\s+Total energy\s+Change\s+\-+"
        row = (
            r"(\d+)"
            r"\s+([A-Za-z\./_]+\s?[A-Za-z\./]+)"
            r"\s+(-?\d+\.\d+(?:[eE][+\-]?\d+)?)"
            r"\s+(-?\d+\.\d+(?:[eE][+\-]?\d+)?)"
            r"(\s+-?\d+\.\d+(?:[eE][+\-]?\d+)?)?"
            r"\s+(-?\d+\.\d+(?:[eE][+\-]?\d+)?)"
            r"(\s+-?\d+\.\d+(?:[eE][+\-]?\d+)?)?"
        )

        footer = r"^$"

        scfs = self.read_table_pattern(
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer,
            last_one_only=False,
            strip=["HFX_MEM_INFO", "*"],
        )

        self.data["electronic_steps"] = []
        self.data["convergence"] = []
        self.data["scf_time"] = []
        for i in scfs:
            self.data["scf_time"].append([float(j[-4]) for j in i])
            self.data["convergence"].append([float(j[-3]) for j in i if j[-3] != "None"])
            self.data["electronic_steps"].append([float(j[-2]) for j in i])

    def parse_timing(self):
        """Parse the timing info (how long did the run take)."""
        header = r"SUBROUTINE\s+CALLS\s+ASD\s+SELF TIME\s+TOTAL TIME\s+MAXIMUM\s+AVERAGE\s+MAXIMUM\s+AVERAGE\s+MAXIMUM"
        row = r"(\w+)\s+(.+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
        footer = r"\-+"

        timing = self.read_table_pattern(
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer,
            last_one_only=True,
            postprocess=postprocessor,
        )
        self.timing = {}
        for time in timing:
            self.timing[time[0]] = {
                "calls": {"max": time[1]},
                "asd": time[2],
                "self_time": {"average": time[3], "maximum": time[4]},
                "total_time": {"average": time[5], "maximum": time[6]},
            }

    def parse_opt_steps(self):
        """Parse the geometry optimization information."""
        # "Information at step =" Summary block (floating point terms)
        total_energy = re.compile(r"\s+Total Energy\s+=\s+(-?\d+.\d+)")
        real_energy_change = re.compile(r"\s+Real energy change\s+=\s+(-?\d+.\d+)")
        predicted_change_in_energy = re.compile(r"\s+Predicted change in energy\s+=\s+(-?\d+.\d+)")
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
                "predicted_change_in_energy": predicted_change_in_energy,
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

        # "Information at step =" Summary block (bool terms)
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
            postprocess=postprocessor,
        )

    def parse_mulliken(self):
        """Parse the mulliken population analysis info for each step."""
        header = r"Mulliken Population Analysis.+Net charge"
        pattern = r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        footer = r".+Total charge"

        if self.read_table_pattern(
            header_pattern=header,
            row_pattern=pattern,
            footer_pattern=footer,
            last_one_only=False,
        ):
            print("Found data, but not yet implemented!")

    def parse_hirshfeld(self):
        """Parse the Hirshfeld population analysis for each step."""
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
                self.structures[i].add_site_property("hirshfeld", hirshfeld)
        else:
            pattern = (
                r"\s+(\d)\s+(\w+)\s+(\d+)\s+(-?\d+\.\d+)\s+"
                r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
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
                self.structures[i].add_site_property("hirshfeld", hirshfeld)

    def parse_mo_eigenvalues(self):
        """Parse the MO eigenvalues from the CP2K output file. Will get the eigenvalues (and band gap)
        at each ionic step (if more than one exist).

        Everything is decomposed by spin channel. If calculation was performed without spin
        polarization, then only Spin.up will be present, which represents the average of up and
        down.
        """
        eigenvalues = []
        efermi = []

        with zopen(self.filename, mode="rt") as file:
            lines = iter(file.readlines())
            for line in lines:
                try:
                    if " occupied subspace spin" in line:
                        eigenvalues.append(
                            {
                                "occupied": {Spin.up: [], Spin.down: []},
                                "unoccupied": {Spin.up: [], Spin.down: []},
                            }
                        )
                        efermi.append({Spin.up: np.nan, Spin.down: np.nan})
                        next(lines)
                        while True:
                            line = next(lines)
                            if "Fermi" in line:
                                efermi[-1][Spin.up] = float(line.split()[-1])
                                break
                            eigenvalues[-1]["occupied"][Spin.up] += [Ha_to_eV * float(val) for val in line.split()]
                        next(lines)
                        line = next(lines)
                        if " occupied subspace spin" in line:
                            next(lines)
                            while True:
                                line = next(lines)
                                if "Fermi" in line:
                                    efermi[-1][Spin.down] = float(line.split()[-1])
                                    break
                                eigenvalues[-1]["occupied"][Spin.down] += [
                                    Ha_to_eV * float(val) for val in line.split()
                                ]
                    if " unoccupied subspace spin" in line:
                        next(lines)
                        line = next(lines)
                        while True:
                            if "WARNING : did not converge" in line:
                                warnings.warn(
                                    "Convergence of eigenvalues for unoccupied subspace spin 1 did NOT converge"
                                )
                                next(lines)
                                next(lines)
                                next(lines)
                                line = next(lines)
                                eigenvalues[-1]["unoccupied"][Spin.up] += [
                                    Ha_to_eV * float(line) for line in line.split()
                                ]
                                next(lines)
                                line = next(lines)
                                break

                            if "convergence" in line:
                                line = next(lines)

                            if "eigenvalues" in line.lower() or "HOMO" in line or "|" in line:
                                break
                            eigenvalues[-1]["unoccupied"][Spin.up] += [Ha_to_eV * float(val) for val in line.split()]
                            line = next(lines)

                        if " unoccupied subspace spin" in line:
                            next(lines)
                            line = next(lines)
                            while True:
                                if "WARNING : did not converge" in line:
                                    warnings.warn(
                                        "Convergence of eigenvalues for unoccupied subspace spin 2 did NOT converge"
                                    )
                                    next(lines)
                                    next(lines)
                                    next(lines)
                                    line = next(lines)
                                    eigenvalues[-1]["unoccupied"][Spin.down] += [
                                        Ha_to_eV * float(line) for line in line.split()
                                    ]
                                    break

                                if "convergence" in line:
                                    line = next(lines)

                                if "HOMO" in line or "|" in line:
                                    next(lines)
                                    break
                                try:
                                    eigenvalues[-1]["unoccupied"][Spin.down] += [
                                        Ha_to_eV * float(val) for val in line.split()
                                    ]
                                except AttributeError:
                                    break
                                line = next(lines)

                except ValueError:
                    eigenvalues = [
                        {
                            "occupied": {Spin.up: None, Spin.down: None},
                            "unoccupied": {Spin.up: None, Spin.down: None},
                        }
                    ]
                    warnings.warn("Convergence of eigenvalues for one or more subspaces did NOT converge")

        self.data["eigenvalues"] = eigenvalues

        if len(eigenvalues) == 0:
            return

        if self.spin_polarized:
            self.data["vbm"] = {
                Spin.up: np.max(eigenvalues[-1]["occupied"][Spin.up]),
                Spin.down: np.max(eigenvalues[-1]["occupied"][Spin.down]),
            }
            self.data["cbm"] = {
                Spin.up: np.nanmin(eigenvalues[-1]["unoccupied"][Spin.up] or np.nan),
                Spin.down: np.nanmin(eigenvalues[-1]["unoccupied"][Spin.down] or np.nan),
            }
            self.vbm = np.nanmean(list(self.data["vbm"].values()))
            self.cbm = np.nanmean(list(self.data["cbm"].values()))
            self.efermi = np.nanmean(list(efermi[-1].values()))
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

        n_occ = len(eigenvalues[-1]["occupied"][Spin.up])
        n_unocc = len(eigenvalues[-1]["unoccupied"][Spin.up])
        self.data["tdos"] = Dos(
            efermi=self.vbm + 1e-6,
            energies=list(eigenvalues[-1]["occupied"][Spin.up]) + list(eigenvalues[-1]["unoccupied"][Spin.down]),
            densities={
                Spin.up: [1 for _ in range(n_occ)] + [0 for _ in range(n_unocc)],
                Spin.down: [1 for _ in range(n_occ)] + [0 for _ in range(n_unocc)],
            },
        )

    def parse_homo_lumo(self):
        """Find the HOMO - LUMO gap in [eV]. Returns the last value. For gaps/eigenvalues decomposed
        by spin up/spin down channel and over many ionic steps, see parse_mo_eigenvalues().
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

    def parse_dos(self, dos_file=None, pdos_files=None, ldos_files=None):
        """Parse the dos files produced by CP2K calculation. CP2K produces different files based
        on the input file rather than assimilating them all into one file.

        One file type is the overall DOS file, which is used for k-point calculations. For
        non-kpoint calculation, the overall DOS is generally not calculated, but the
        element-projected pDOS is. Separate files are created for each spin channel and each
        atom kind. If requested, CP2K can also do site/local projected dos (ldos). Each site
        requested will have a separate file for each spin channel (if spin polarized calculation
        is performed).

        If possible, this function will assimilate the ldos files into a CompleteDos object.
        Either provide a list of PDOS file paths, or use glob to find the .pdos_ALPHA extension
        in the calculation directory.

        Args:
            dos_file (str): Name of the dos file, otherwise will be inferred
            pdos_files (list): list of pdos file paths, otherwise they will be inferred
            ldos_files (list): list of ldos file paths, otherwise they will be inferred
        """
        if dos_file is None:
            dos_file = self.filenames["DOS"][0] if self.filenames["DOS"] else None

        if pdos_files is None:
            pdos_files = self.filenames["PDOS"]

        if ldos_files is None:
            ldos_files = self.filenames["LDOS"]

        tdos, pdoss, ldoss = None, {}, {}
        # Parse specie projected dos
        for pdos_file in pdos_files:
            _pdos, _tdos = parse_pdos(pdos_file, total=True)
            for k in _pdos:
                if k in pdoss:
                    for orbital in _pdos[k]:
                        pdoss[k][orbital].densities.update(_pdos[k][orbital].densities)
                else:
                    pdoss.update(_pdos)
            if not tdos:
                tdos = _tdos
            else:
                for k, v in _tdos.densities.copy().items():
                    if k not in tdos.densities:
                        tdos.densities[Spin(int(k))] = [0] * len(v)
                    tdos.densities[k] = np.array(tdos.densities[k]) + np.array(_tdos.densities[k])

        # parse any site-projected dos
        for ldos_file in ldos_files:
            _pdos = parse_pdos(ldos_file)
            for k in _pdos:
                if k in ldoss:
                    for orbital in _pdos[k]:
                        ldoss[k][orbital].densities.update(_pdos[k][orbital].densities)
                else:
                    ldoss.update(_pdos)

        self.data["pdos"] = jsanitize(pdoss, strict=True)
        self.data["ldos"] = jsanitize(ldoss, strict=True)

        self.data["tdos"] = parse_dos(dos_file) if dos_file else tdos

        if self.data.get("tdos"):
            self.band_gap = self.data["tdos"].get_gap()
            self.cbm, self.vbm = self.data["tdos"].get_cbm_vbm()
            self.efermi = (self.cbm + self.vbm) / 2

        # If number of site-projected dos == number of sites, assume they are bijective
        # and create the CompleteDos object
        _ldoss = {}
        if self.initial_structure and len(ldoss) == len(self.initial_structure):
            for k, lds in ldoss.items():
                _ldoss[self.initial_structure[int(k) - 1]] = {Orbital(orb): lds[orb].densities for orb in lds}

            self.data["cdos"] = CompleteDos(self.final_structure, total_dos=tdos, pdoss=_ldoss)

    @property
    def complete_dos(self) -> CompleteDos | None:
        """Complete dos object if it has been parsed."""
        return self.data.get("cdos")

    @property
    def band_structure(self) -> BandStructure | None:
        """Band structure object if it has been parsed."""
        return self.data.get("band_structure")

    def parse_bandstructure(self, bandstructure_filename=None) -> None:
        """Parse a CP2K bandstructure file.

        Args:
            bandstructure_filename: Filename containing bandstructure info. If
            not provided, then the pmg name of "BAND.bs" will be assumed by
            the filename parser.
        """
        if not bandstructure_filename:
            if self.filenames["band_structure"]:
                bandstructure_filename = self.filenames["band_structure"][0]
            else:
                return

        with open(bandstructure_filename, encoding="utf-8") as file:
            lines = file.read().split("\n")

        bands_data = np.loadtxt(bandstructure_filename)
        nkpts = int(lines[0].split()[6])
        nbands = int(lines[0].split()[-2])
        rec_lat = (
            self.final_structure.lattice.reciprocal_lattice
            if self.final_structure
            else self.initial_structure.lattice.reciprocal_lattice
        )

        labels = {}
        kpts = []
        nkpts = 0
        for line in lines:
            if not line.startswith("#"):
                continue
            if line.split()[1] == "Set":
                nkpts += int(lines[0].split()[6])
            elif line.split()[1] == "Point":
                kpts.append(list(map(float, line.split()[-4:-1])))
            elif line.split()[1] == "Special":
                splt = line.split()
                label = splt[7]
                if label.upper() == "GAMMA":
                    label = "\\Gamma"
                kpt = np.array(splt[4:7]).astype(float).tolist()
                if label.upper() != "NONE":
                    labels[label] = kpt

        if self.spin_polarized:
            kpts = kpts[::2]

        eigenvals = {}
        if self.spin_polarized:
            up = bands_data.reshape(-1, nbands * 2, bands_data.shape[1])[:, :nbands].reshape(-1, bands_data.shape[1])
            down = bands_data.reshape(-1, nbands * 2, bands_data.shape[1])[:, nbands:].reshape(-1, bands_data.shape[1])
            eigenvals = {
                Spin.up: up[:, 1].reshape((nkpts, nbands)).T.tolist(),
                Spin.down: down[:, 1].reshape((nkpts, nbands)).T.tolist(),
            }
        else:
            eigenvals = {Spin.up: bands_data.reshape((nbands, nkpts))}

        occ = bands_data[:, 1][bands_data[:, -1] != 0.0]
        homo = np.max(occ)
        unocc = bands_data[:, 1][bands_data[:, -1] == 0.0]
        lumo = np.min(unocc)
        efermi = (lumo + homo) / 2
        self.efermi = efermi

        self.data["band_structure"] = BandStructureSymmLine(
            kpoints=kpts,
            eigenvals=eigenvals,
            lattice=rec_lat,
            efermi=efermi,
            labels_dict=labels,
            structure=self.final_structure,
            projections=None,  # not implemented in CP2K
        )

        self.band_gap = self.data["band_structure"].get_band_gap().get("energy")
        self.vbm = self.data["band_structure"].get_vbm().get("energy")
        self.cbm = self.data["band_structure"].get_cbm().get("energy")

    def parse_hyperfine(self, hyperfine_filename=None):
        """Parse a file containing hyperfine coupling tensors for each atomic site."""
        if not hyperfine_filename:
            if self.filenames["hyperfine_tensor"]:
                hyperfine_filename = self.filenames["hyperfine_tensor"][0]
            else:
                return None

        with zopen(hyperfine_filename, mode="rt") as file:
            lines = [line for line in file.read().split("\n") if line]

        hyperfine = [[] for _ in self.ionic_steps]
        for i in range(2, len(lines), 5):
            x = list(map(float, lines[i + 2].split()[-3:]))
            y = list(map(float, lines[i + 3].split()[-3:]))
            z = list(map(float, lines[i + 4].split()[-3:]))
            hyperfine[-1].append([x, y, z])

        self.data["hyperfine_tensor"] = hyperfine
        return hyperfine

    def parse_gtensor(self, gtensor_filename=None):
        """Parse a file containing g tensor."""
        if not gtensor_filename:
            if self.filenames["g_tensor"]:
                gtensor_filename = self.filenames["g_tensor"][0]
            else:
                return None

        with zopen(gtensor_filename, mode="rt") as file:
            lines = [line for line in file.read().split("\n") if line]

        data = {}
        data["gmatrix_zke"] = []
        data["gmatrix_so"] = []
        data["gmatrix_soo"] = []
        data["gmatrix_total"] = []
        data["gtensor_total"] = []
        data["delta_g"] = []
        ionic = -1
        dat = None
        for line in lines:
            first = line.strip()
            if first == "G tensor":
                ionic += 1
                for d in data.values():
                    d.append([])
            elif first in data:
                dat = first
            elif first.startswith("delta_g"):
                dat = "delta_g"
            else:
                splt = [postprocessor(s) for s in line.split()]
                splt = [s for s in splt if isinstance(s, float)]
                data[dat][ionic].append(list(map(float, splt[-3:])))
        self.data.update(data)
        return data["gtensor_total"][-1]

    def parse_chi_tensor(self, chi_filename=None):
        """Parse the magnetic susceptibility tensor."""
        if not chi_filename:
            if self.filenames["chi_tensor"]:
                chi_filename = self.filenames["chi_tensor"][0]
            else:
                return None

        with zopen(chi_filename, mode="rt") as file:
            lines = [line for line in file.read().split("\n") if line]

        data = {}
        data["chi_soft"] = []
        data["chi_local"] = []
        data["chi_total"] = []
        data["chi_total_ppm_cgs"] = []
        data["PV1"] = []
        data["PV2"] = []
        data["PV3"] = []
        data["ISO"] = []
        data["ANISO"] = []
        ionic = -1
        dat = None
        for line in lines:
            first = line.strip()
            if first == "Magnetic Susceptibility Tensor":
                ionic += 1
                for d in data.values():
                    d.append([])
            elif first in data:
                dat = first
            elif "SOFT" in first:
                dat = "chi_soft"
            elif "LOCAL" in first:
                dat = "chi_local"
            elif "Total" in first:
                dat = "chi_total_ppm_cgs" if "ppm" in first else "chi_total"
            elif first.startswith("PV1"):
                splt = [postprocessor(s) for s in line.split()]
                splt = [s for s in splt if isinstance(s, float)]
                data["PV1"][ionic] = splt[0]
                data["PV2"][ionic] = splt[1]
                data["PV3"][ionic] = splt[2]
            elif first.startswith("ISO"):
                splt = [postprocessor(s) for s in line.split()]
                splt = [s for s in splt if isinstance(s, float)]
                data["ISO"][ionic] = splt[0]
                data["ANISO"][ionic] = splt[1]
            else:
                splt = [postprocessor(s) for s in line.split()]
                splt = [s for s in splt if isinstance(s, float)]
                data[dat][ionic].append(list(map(float, splt)))
        self.data.update(data)
        return data["chi_total"][-1]

    def parse_nmr_shift(self):
        """Parse NMR calculation."""
        raise NotImplementedError("NMR Parsing not yet implemented")

    def parse_tddfpt(self):
        """Parse TDDFPT calculation."""
        raise NotImplementedError("TDDFPT excited states parsing not yet implemented")

    def parse_raman(self):
        """Parse raman calculation."""
        raise NotImplementedError("Raman parsing not yet implemented")

    @staticmethod
    def _gauss_smear(densities, energies, npts, width):
        """Return a gaussian smeared DOS."""
        if not width:
            return densities

        dct = np.zeros(npts)
        e_s = np.linspace(min(energies), max(energies), npts)

        for e, _pd in zip(energies, densities, strict=False):
            weight = np.exp(-(((e_s - e) / width) ** 2)) / (np.sqrt(np.pi) * width)
            dct += _pd * weight

        return dct

    def read_pattern(self, patterns, reverse=False, terminate_on_match=False, postprocess=str):
        r"""Originally from pymatgen.io.vasp.outputs.Outcar.

        General pattern reading. Uses monty's regrep method. Takes the same
        arguments.

        Args:
            patterns (dict): A dict of patterns, e.g.
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
        for k in patterns:
            self.data[k] = [i[0] for i in matches.get(k, [])]

    def read_table_pattern(
        self,
        header_pattern,
        row_pattern,
        footer_pattern,
        postprocess=str,
        attribute_name=None,
        last_one_only=True,
        strip=None,
    ):
        r"""This function originated in pymatgen.io.vasp.outputs.Outcar.

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
            strip (list): Whether or not to strip contents out of the file before
                reading for a table pattern. This is mainly used by parse_scf_opt(),
                to strip HFX info out of the SCF loop start or DFT+U warnings out
                of the SCF loop iterations.

        Returns:
            List of tables. 1) A table is a list of rows. 2) A row if either a list of
            attribute values in case the the capturing group is defined without name in
            row_pattern, or a dict in case that named capturing groups are defined by
            row_pattern.
        """
        with zopen(self.filename, mode="rt") as file:
            if strip:
                lines = file.readlines()
                text = "".join(
                    [
                        lines[i]
                        for i in range(1, len(lines) - 1)
                        if all(
                            not lines[i].strip().startswith(c)
                            and not lines[i - 1].strip().startswith(c)
                            and not lines[i + 1].strip().startswith(c)
                            for c in strip
                        )
                    ]
                )
            else:
                text = file.read()

        table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + row_pattern + r")+)\s+" + footer_pattern
        table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
        rp = re.compile(row_pattern)
        tables = []
        for mt in table_pattern.finditer(text):
            table_body_text = mt.group("table_body")
            table_contents = []
            for line in table_body_text.split("\n"):
                ml = rp.search(line)
                if ml is None:
                    continue
                d = ml.groupdict()
                if len(d) > 0:
                    processed_line = {k: postprocess(v) for k, v in d.items()}
                else:
                    processed_line = [postprocess(v) for v in ml.groups()]
                table_contents.append(processed_line)
            tables.append(table_contents)
        retained_data = tables[-1] if last_one_only else tables
        if attribute_name is not None:
            self.data[attribute_name] = retained_data
        return retained_data

    def as_dict(self):
        """Return dictionary representation of the output."""
        dct = {"input": {}, "output": {}}
        dct["total_time"] = self.timing["CP2K"]["total_time"]["maximum"]
        dct["run_type"] = self.run_type
        dct["input"]["global"] = self.data.get("global")
        dct["input"]["dft"] = self.data.get("dft")
        dct["input"]["scf"] = self.data.get("scf")
        dct["input"]["qs"] = self.data.get("QS")
        dct["input"]["run_type"] = self.run_type
        dct["input"]["calculation_type"] = self.calculation_type
        dct["input"]["structure"] = self.initial_structure.as_dict()
        dct["input"]["atomic_kind_info"] = self.data.get("atomic_kind_info")
        dct["input"]["cp2k_input"] = self.input.as_dict()
        dct["ran_successfully"] = self.completed
        dct["cp2k_version"] = self.cp2k_version
        dct["output"]["structure"] = self.final_structure.as_dict()
        dct["output"]["forces"] = self.data.get("forces", [None])[-1]
        dct["output"]["stress"] = self.data.get("stress_tensor", [None])[-1]
        dct["output"]["ionic_steps"] = [
            {k: v.as_dict() if isinstance(v, MSONable) else v for k, v in step.items()} for step in self.ionic_steps
        ]
        dct["composition"] = self.composition.as_dict()
        dct["output"]["energy"] = self.final_energy
        dct["output"]["energy_per_atom"] = self.final_energy / self.composition.num_atoms
        dct["output"]["bandgap"] = self.cbm - self.vbm if self.cbm and self.vbm else None
        dct["output"]["cbm"] = self.cbm
        dct["output"]["vbm"] = self.vbm
        dct["output"]["efermi"] = self.efermi
        dct["output"]["is_metal"] = self.is_metal
        return dct


# TODO should store as pandas? Maybe it should be stored as a dict so it's python native
def parse_energy_file(energy_file):
    """Parse energy file for calculations with multiple ionic steps."""
    columns = [
        "step",
        "kinetic_energy",
        "temp",
        "potential_energy",
        "conserved_quantity",
        "used_time",
    ]
    df_energies = pd.read_csv(energy_file, skiprows=1, names=columns, sep=r"\s+")
    df_energies["kinetic_energy"] *= Ha_to_eV
    df_energies["potential_energy"] *= Ha_to_eV
    df_energies["conserved_quantity"] *= Ha_to_eV
    df_energies = df_energies.astype(float)
    return {c: df_energies[c].to_numpy() for c in columns}


# TODO: The DOS file that CP2K outputs as of 2022.1 seems to have a lot of problems.
def parse_dos(dos_file=None):
    """Parse a dos file. This format is different from the pdos files."""
    data = np.loadtxt(dos_file)
    data[:, 0] *= Ha_to_eV
    energies = data[:, 0]
    vbm_top = None
    for idx, val in enumerate(data[:, 1]):
        if val == 0:
            break
        vbm_top = idx

    efermi = energies[vbm_top] + 1e-6
    densities = {Spin.up: data[:, 1]}

    if data.shape[1] > 3:
        densities[Spin.down] = data[:, 3]
    return Dos(efermi=efermi, energies=energies, densities=densities)


def parse_pdos(dos_file=None, spin_channel=None, total=False):
    """
    Parse a single DOS file created by CP2K. Must contain one PDOS snapshot. i.e. you cannot
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
    spin = Spin(spin_channel) if spin_channel else Spin.down if "BETA" in os.path.split(dos_file)[-1] else Spin.up

    with zopen(dos_file, mode="rt") as file:
        lines = file.readlines()
        kind = re.search(r"atomic kind\s(.*)\sat iter", lines[0]) or re.search(r"list\s(\d+)\s(.*)\sat iter", lines[0])
        kind = kind.groups()[0]

        header = re.split(r"\s{2,}", lines[1].replace("#", "").strip())[2:]
        dat = np.loadtxt(dos_file)

        def cp2k_to_pmg_labels(label: str) -> str:
            if label == "p":
                return "px"
            if label == "d":
                return "dxy"
            if label == "f":
                return "f_3"
            if label == "d-2":
                return "dxy"
            if label == "d-1":
                return "dyz"
            if label == "d0":
                return "dz2"
            if label == "d+1":
                return "dxz"
            if label == "d+2":
                return "dx2"
            if label == "f-3":
                return "f_3"
            if label == "f-2":
                return "f_2"
            if label == "f-1":
                return "f_1"
            if label == "f0":
                return "f0"
            if label == "f+1":
                return "f1"
            if label == "f+2":
                return "f2"
            if label == "f+3":
                return "f3"
            return label

        header = [cp2k_to_pmg_labels(h) for h in header]

        data = np.delete(dat, 0, 1)
        occupations = data[:, 1]
        data = np.delete(data, 1, 1)
        data[:, 0] *= Ha_to_eV
        energies = data[:, 0]
        vbm_top = None
        for idx, occu in enumerate(occupations):
            if occu == 0:
                break
            vbm_top = idx

        # set Fermi level to be vbm plus tolerance for
        # PMG compatibility
        # *not* middle of the gap, which pdos might report
        efermi = energies[vbm_top] + 1e-6

        # for pymatgen's dos class. VASP creates an evenly spaced grid of energy states, which
        # leads to 0 density states in the band gap. CP2K does not do this. PMG's Dos class was
        # created with VASP in mind so the way it searches for vbm and cbm relies on grid points
        # in between VBM and CBM, so here we introduce trivial ones
        energies = np.insert(
            energies,
            vbm_top + 1,
            np.linspace(energies[vbm_top] + 1e-6, energies[vbm_top + 1] - 1e-6, 2),
        )
        data = np.insert(data, vbm_top + 1, np.zeros((2, data.shape[1])), axis=0)

        pdos = {
            kind: {
                getattr(Orbital, h): Dos(efermi=efermi, energies=energies, densities={spin: data[:, i + 1]})
                for i, h in enumerate(header)
            }
        }
        if total:
            tdos = Dos(
                efermi=efermi,
                energies=energies,
                densities={spin: np.sum(data[:, 1:], axis=1)},
            )
            return pdos, tdos
        return pdos
