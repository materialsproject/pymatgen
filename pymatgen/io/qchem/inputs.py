"""Classes for reading/manipulating/writing QChem input files."""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING, Literal

import numpy as np
from monty.io import zopen

from pymatgen.core import Molecule
from pymatgen.io.core import InputFile

from .utils import lower_and_check_unique, read_pattern, read_table_pattern

if TYPE_CHECKING:
    from pathlib import Path

__author__ = "Brandon Wood, Samuel Blau, Shyam Dwaraknath, Julian Self, Evan Spotte-Smith, Ryan Kingsbury"
__copyright__ = "Copyright 2018-2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__credits__ = "Xiaohui Qu"

logger = logging.getLogger(__name__)


class QCInput(InputFile):
    """
    An object representing a QChem input file. QCInput attributes represent different sections of a QChem input file.
    To add a new section one needs to modify __init__, __str__, from_sting and add static methods
    to read and write the new section i.e. section_template and read_section. By design, there is very little (or no)
    checking that input parameters conform to the appropriate QChem format, this responsible lands on the user or a
    separate error handling software.
    """

    def __init__(
        self,
        molecule: Molecule | list[Molecule] | Literal["read"],
        rem: dict,
        opt: dict[str, list] | None = None,
        pcm: dict | None = None,
        solvent: dict | None = None,
        smx: dict | None = None,
        scan: dict[str, list] | None = None,
        van_der_waals: dict[str, float] | None = None,
        vdw_mode: str = "atomic",
        plots: dict | None = None,
        nbo: dict | None = None,
        geom_opt: dict | None = None,
        cdft: list[list[dict]] | None = None,
        almo_coupling: list[list[tuple[int, int]]] | None = None,
        svp: dict | None = None,
        pcm_nonels: dict | None = None,
    ):
        """
        Args:
            molecule (pymatgen Molecule object, list of Molecule objects, or "read"):
                Input molecule(s). molecule can be set as a pymatgen Molecule object, a list of such
                Molecule objects, or as the string "read". "read" can be used in multi_job QChem input
                files where the molecule is read in from the previous calculation.
            rem (dict):
                A dictionary of all the input parameters for the rem section of QChem input file.
                Ex. rem = {'method': 'rimp2', 'basis': '6-31*G++' ... }
            opt (dict of lists):
                A dictionary of opt sections, where each opt section is a key and the corresponding
                values are a list of strings. Strings must be formatted as instructed by the QChem manual.
                The different opt sections are: CONSTRAINT, FIXED, DUMMY, and CONNECT
                Ex. opt = {"CONSTRAINT": ["tors 2 3 4 5 25.0", "tors 2 5 7 9 80.0"], "FIXED": ["2 XY"]}
            pcm (dict):
                A dictionary of the PCM section, defining behavior for use of the polarizable continuum model.
                Ex: pcm = {"theory": "cpcm", "hpoints": 194}
            solvent (dict):
                A dictionary defining the solvent parameters used with PCM.
                Ex: solvent = {"dielectric": 78.39, "temperature": 298.15}
            smx (dict):
                A dictionary defining solvent parameters used with the SMD method, a solvent method that adds
                short-range terms to PCM.
                Ex: smx = {"solvent": "water"}
            scan (dict of lists):
                A dictionary of scan variables. Because two constraints of the same type are allowed (for instance, two
                torsions or two bond stretches), each TYPE of variable (stre, bend, tors) should be its own key in the
                dict, rather than each variable. Note that the total number of variable (sum of lengths of all lists)
                CANNOT be
                more than two.
                Ex. scan = {"stre": ["3 6 1.5 1.9 0.1"], "tors": ["1 2 3 4 -180 180 15"]}
            van_der_waals (dict):
                A dictionary of custom van der Waals radii to be used when constructing cavities for the PCM
                model or when computing, e.g. Mulliken charges. They keys are strs whose meaning depends on
                the value of vdw_mode, and the values are the custom radii in angstroms.
            vdw_mode (str): Method of specifying custom van der Waals radii - 'atomic' or 'sequential'.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., 12 = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
            plots (dict):
                    A dictionary of all the input parameters for the plots section of the QChem input file.
            nbo (dict):
                    A dictionary of all the input parameters for the nbo section of the QChem input file.
            geom_opt (dict):
                    A dictionary of input parameters for the geom_opt section of the QChem input file.
                    This section is required when using the new libopt3 geometry optimizer.
            cdft (list of lists of dicts):
                    A list of lists of dictionaries, where each dictionary represents a charge constraint in the
                    cdft section of the QChem input file.

                    Each entry in the main list represents one state (allowing for multi-configuration calculations
                    using constrained density functional theory - configuration interaction (CDFT-CI).
                    Each state is represented by a list, which itself contains some number of constraints
                    (dictionaries).

                    Ex:

                    1. For a single-state calculation with two constraints:
                     cdft=[[
                        {"value": 1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [2], "types": [None]},
                        {"value": 2.0, "coefficients": [1.0, -1.0], "first_atoms": [1, 17], "last_atoms": [3, 19],
                            "types": ["s"]}
                    ]]

                    Note that a type of None will default to a charge constraint (which can also be accessed by
                    requesting a type of "c" or "charge".

                    2. For a multi-reference calculation:
                    cdft=[
                        [
                            {"value": 1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27],
                                "types": ["c"]},
                            {"value": 0.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27],
                                "types": ["s"]},
                        ],
                        [
                            {"value": 0.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27],
                                "types": ["c"]},
                            {"value": -1.0, "coefficients": [1.0], "first_atoms": [1], "last_atoms": [27],
                                "types": ["s"]},
                        ]
                    ]
            almo_coupling (list of lists of int 2-tuples):
                A list of lists of int 2-tuples used for calculations of diabatization and state coupling calculations
                    relying on the absolutely localized molecular orbitals (ALMO) methodology. Each entry in the main
                    list represents a single state (two states are included in an ALMO calculation). Within a single
                    state, each 2-tuple represents the charge and spin multiplicity of a single fragment.
                ex: almo=[
                            [
                                (1, 2),
                                (0, 1)
                            ],
                            [
                                (0, 1),
                                (1, 2)
                            ]
                        ]
        """
        self.molecule = molecule
        self.rem = lower_and_check_unique(rem)
        self.opt = opt
        self.pcm = lower_and_check_unique(pcm)
        self.solvent = lower_and_check_unique(solvent)
        self.smx = lower_and_check_unique(smx)
        self.scan = lower_and_check_unique(scan)
        self.van_der_waals = lower_and_check_unique(van_der_waals)
        self.vdw_mode = vdw_mode
        self.plots = lower_and_check_unique(plots)
        self.nbo = lower_and_check_unique(nbo)
        self.geom_opt = lower_and_check_unique(geom_opt)
        self.cdft = cdft
        self.almo_coupling = almo_coupling
        self.svp = lower_and_check_unique(svp)
        self.pcm_nonels = lower_and_check_unique(pcm_nonels)

        # Make sure rem is valid:
        #   - Has a basis
        #   - Has a method or DFT exchange functional
        #   - Has a valid job_type or jobtype

        valid_job_types = [
            "opt",
            "optimization",
            "sp",
            "freq",
            "frequency",
            "force",
            "nmr",
            "ts",
            "pes_scan",
        ]

        if "basis" not in self.rem:
            raise ValueError("The rem dictionary must contain a 'basis' entry")
        if "method" not in self.rem and "exchange" not in self.rem:
            raise ValueError("The rem dictionary must contain either a 'method' entry or an 'exchange' entry")
        if "job_type" not in self.rem:
            raise ValueError("The rem dictionary must contain a 'job_type' entry")
        if self.rem.get("job_type").lower() not in valid_job_types:
            raise ValueError("The rem dictionary must contain a valid 'job_type' entry")

        # Still to do:
        #   - Check that the method or functional is valid
        #   - Check that basis is valid
        #   - Check that basis is defined for all species in the molecule
        #   - Validity checks specific to job type?
        #   - Check OPT and PCM sections?

    @np.deprecate(message="Use get_str instead")
    def get_string(self, *args, **kwargs) -> str:
        return self.get_str(*args, **kwargs)

    def get_str(self) -> str:
        """Return a string representation of an entire input file."""
        return str(self)

    def __str__(self):
        combined_list = []
        # molecule section
        combined_list.append(self.molecule_template(self.molecule))
        combined_list.append("")
        # rem section
        combined_list.append(self.rem_template(self.rem))
        combined_list.append("")
        # opt section
        if self.opt:
            combined_list.append(self.opt_template(self.opt))
            combined_list.append("")
        # pcm section
        if self.pcm:
            combined_list.append(self.pcm_template(self.pcm))
            combined_list.append("")
        # solvent section
        if self.solvent:
            combined_list.append(self.solvent_template(self.solvent))
            combined_list.append("")
        if self.smx:
            combined_list.append(self.smx_template(self.smx))
            combined_list.append("")
        # section for pes_scan
        if self.scan:
            combined_list.append(self.scan_template(self.scan))
            combined_list.append("")
        # section for van_der_waals radii
        if self.van_der_waals:
            combined_list.append(self.van_der_waals_template(self.van_der_waals, self.vdw_mode))
            combined_list.append("")
        # plots section
        if self.plots:
            combined_list.append(self.plots_template(self.plots))
            combined_list.append("")
        # nbo section
        if self.nbo is not None:
            combined_list.append(self.nbo_template(self.nbo))
            combined_list.append("")
        # geom_opt section
        if self.geom_opt is not None:
            combined_list.append(self.geom_opt_template(self.geom_opt))
            combined_list.append("")
        # cdft section
        if self.cdft is not None:
            combined_list.append(self.cdft_template(self.cdft))
            combined_list.append("")
        # almo section
        if self.almo_coupling is not None:
            combined_list.append(self.almo_template(self.almo_coupling))
            combined_list.append("")
        # svp section
        if self.svp:
            combined_list.append(self.svp_template(self.svp))
            combined_list.append("")
        # pcm_nonels section
        if self.pcm_nonels:
            combined_list.append(self.pcm_nonels_template(self.pcm_nonels))
        return "\n".join(combined_list)

    @staticmethod
    def multi_job_string(job_list: list[QCInput]) -> str:
        """
        Args:
            job_list (): List of jobs.

        Returns:
            (str) String representation of multi job input file.
        """
        multi_job_string = ""
        for i, job_i in enumerate(job_list):
            if i < len(job_list) - 1:
                multi_job_string += str(job_i) + "\n@@@\n\n"
            else:
                multi_job_string += str(job_i)
        return multi_job_string

    @classmethod
    def from_str(cls, string: str) -> QCInput:
        """
        Read QcInput from string.

        Args:
            string (str): String input.

        Returns:
            QcInput
        """
        sections = cls.find_sections(string)
        molecule = cls.read_molecule(string)
        rem = cls.read_rem(string)
        # only molecule and rem are necessary everything else is checked
        opt = pcm = solvent = smx = scan = vdw = None
        vdw_mode = "atomic"
        plots = nbo = geom_opt = cdft = almo_coupling = svp = pcm_nonels = None
        if "opt" in sections:
            opt = cls.read_opt(string)
        if "pcm" in sections:
            pcm = cls.read_pcm(string)
        if "solvent" in sections:
            solvent = cls.read_solvent(string)
        if "smx" in sections:
            smx = cls.read_smx(string)
        if "scan" in sections:
            scan = cls.read_scan(string)
        if "van_der_waals" in sections:
            vdw_mode, vdw = cls.read_vdw(string)
        if "plots" in sections:
            plots = cls.read_plots(string)
        if "nbo" in sections:
            nbo = cls.read_nbo(string)
        if "geom_opt" in sections:
            geom_opt = cls.read_geom_opt(string)
        if "cdft" in sections:
            cdft = cls.read_cdft(string)
        if "almo_coupling" in sections:
            almo_coupling = cls.read_almo(string)
        if "svp" in sections:
            svp = cls.read_svp(string)
        if "pcm_nonels" in sections:
            pcm_nonels = cls.read_pcm_nonels(string)
        return cls(
            molecule,
            rem,
            opt=opt,
            solvent=solvent,
            pcm=pcm,
            smx=smx,
            scan=scan,
            van_der_waals=vdw,
            vdw_mode=vdw_mode,
            plots=plots,
            nbo=nbo,
            geom_opt=geom_opt,
            cdft=cdft,
            almo_coupling=almo_coupling,
            svp=svp,
            pcm_nonels=pcm_nonels,
        )

    @staticmethod
    def write_multi_job_file(job_list: list[QCInput], filename: str):
        """
        Write a multijob file.

        Args:
            job_list (): List of jobs.
            filename (): Filename
        """
        with zopen(filename, "wt") as f:
            f.write(QCInput.multi_job_string(job_list))

    @staticmethod
    def from_file(filename: str | Path) -> QCInput:
        """
        Create QcInput from file.

        Args:
            filename (str): Filename

        Returns:
            QcInput
        """
        with zopen(filename, "rt") as f:
            return QCInput.from_str(f.read())

    @classmethod
    def from_multi_jobs_file(cls, filename: str) -> list[QCInput]:
        """
        Create list of QcInput from a file.

        Args:
            filename (str): Filename

        Returns:
            List of QCInput objects
        """
        with zopen(filename, "rt") as f:
            # the delimiter between QChem jobs is @@@
            multi_job_strings = f.read().split("@@@")
            # list of individual QChem jobs
            return [cls.from_str(i) for i in multi_job_strings]

    @staticmethod
    def molecule_template(molecule: Molecule | list[Molecule] | Literal["read"]) -> str:
        """
        Args:
            molecule (Molecule, list of Molecules, or "read").

        Returns:
            (str) Molecule template.

        """
        # TODO: add ghost atoms
        mol_list = []
        mol_list.append("$molecule")

        # Edge case; can't express molecule as fragments with only one fragment
        if isinstance(molecule, list) and len(molecule) == 1:
            molecule = molecule[0]

        if isinstance(molecule, str):
            if molecule == "read":
                mol_list.append(" read")
            else:
                raise ValueError('The only acceptable text value for molecule is "read"')
        elif isinstance(molecule, Molecule):
            mol_list.append(f" {int(molecule.charge)} {molecule.spin_multiplicity}")
            for site in molecule.sites:
                mol_list.append(f" {site.species_string}     {site.x: .10f}     {site.y: .10f}     {site.z: .10f}")
        else:
            overall_charge = sum(x.charge for x in molecule)
            unpaired_electrons = sum(x.spin_multiplicity - 1 for x in molecule)
            overall_spin = unpaired_electrons + 1

            mol_list.append(f" {int(overall_charge)} {int(overall_spin)}")

            for fragment in molecule:
                mol_list.append("--")
                mol_list.append(f" {int(fragment.charge)} {fragment.spin_multiplicity}")
                for site in fragment.sites:
                    mol_list.append(f" {site.species_string}     {site.x: .10f}     {site.y: .10f}     {site.z: .10f}")

        mol_list.append("$end")
        return "\n".join(mol_list)

    @staticmethod
    def rem_template(rem: dict) -> str:
        """
        Args:
            rem ():

        Returns:
            (str)
        """
        rem_list = []
        rem_list.append("$rem")
        for key, value in rem.items():
            rem_list.append(f"   {key} = {value}")
        rem_list.append("$end")
        return "\n".join(rem_list)

    @staticmethod
    def opt_template(opt: dict[str, list]) -> str:
        """
        Optimization template.

        Args:
            opt ():

        Returns:
            (str)
        """
        opt_list = []
        opt_list.append("$opt")
        # loops over all opt sections
        for key, value in opt.items():
            opt_list.append(f"{key}")
            # loops over all values within the section
            for i in value:
                opt_list.append(f"   {i}")
            opt_list.append(f"END{key}")
            opt_list.append("")
        # this deletes the empty space after the last section
        del opt_list[-1]
        opt_list.append("$end")
        return "\n".join(opt_list)

    @staticmethod
    def pcm_template(pcm: dict) -> str:
        """
        Pcm run template.

        Args:
            pcm ():

        Returns:
            (str)
        """
        pcm_list = []
        pcm_list.append("$pcm")
        for key, value in pcm.items():
            pcm_list.append(f"   {key} {value}")
        pcm_list.append("$end")
        return "\n".join(pcm_list)

    @staticmethod
    def solvent_template(solvent: dict) -> str:
        """
        Solvent template.

        Args:
            solvent ():

        Returns:
            (str)
        """
        solvent_list = []
        solvent_list.append("$solvent")
        for key, value in solvent.items():
            solvent_list.append(f"   {key} {value}")
        solvent_list.append("$end")
        return "\n".join(solvent_list)

    @staticmethod
    def smx_template(smx: dict) -> str:
        """
        Args:
            smx ():

        Returns:
            (str)
        """
        smx_list = []
        smx_list.append("$smx")
        for key, value in smx.items():
            if value == "tetrahydrofuran":
                smx_list.append(f"   {key} thf")
            # Q-Chem bug, see https://talk.q-chem.com/t/smd-unrecognized-solvent/204
            elif value == "dimethyl sulfoxide":
                smx_list.append(f"   {key} dmso")
            else:
                smx_list.append(f"   {key} {value}")
        smx_list.append("$end")
        return "\n".join(smx_list)

    @staticmethod
    def scan_template(scan: dict[str, list]) -> str:
        """
        Args:
            scan (dict): Dictionary with scan section information.
                Ex: {"stre": ["3 6 1.5 1.9 0.1"], "tors": ["1 2 3 4 -180 180 15"]}.

        Returns:
            String representing Q-Chem input format for scan section
        """
        scan_list = []
        scan_list.append("$scan")
        total_vars = sum(len(v) for v in scan.values())
        if total_vars > 2:
            raise ValueError("Q-Chem only supports PES_SCAN with two or less variables.")
        for var_type, variables in scan.items():
            if variables not in [None, []]:
                for var in variables:
                    scan_list.append(f"   {var_type} {var}")
        scan_list.append("$end")
        return "\n".join(scan_list)

    @staticmethod
    def van_der_waals_template(radii: dict[str, float], mode: str = "atomic") -> str:
        """
        Args:
            radii (dict): Dictionary with custom van der Waals radii, in
                Angstroms, keyed by either atomic number or sequential
                atom number (see 'mode' kwarg).
                Ex: {1: 1.20, 12: 1.70}
            mode: 'atomic' or 'sequential'. In 'atomic' mode (default), dict keys
                represent the atomic number associated with each radius (e.g., '12' = carbon).
                In 'sequential' mode, dict keys represent the sequential position of
                a single specific atom in the input structure.
                **NOTE: keys must be given as strings even though they are numbers!**.

        Returns:
            String representing Q-Chem input format for van_der_waals section
        """
        vdw_list = []
        vdw_list.append("$van_der_waals")
        if mode == "atomic":
            vdw_list.append("1")
        elif mode == "sequential":
            vdw_list.append("2")
        else:
            raise ValueError(f"Invalid {mode=}, must be 'atomic' or 'sequential'")

        for num, radius in radii.items():
            vdw_list.append(f"   {num} {radius}")
        vdw_list.append("$end")
        return "\n".join(vdw_list)

    @staticmethod
    def plots_template(plots: dict) -> str:
        """
        Args:
            plots ():

        Returns:
            (str)
        """
        plots_list = []
        plots_list.append("$plots")
        for key, value in plots.items():
            plots_list.append(f"   {key} {value}")
        plots_list.append("$end")
        return "\n".join(plots_list)

    @staticmethod
    def nbo_template(nbo: dict) -> str:
        """
        Args:
            nbo ():

        Returns:
            (str)
        """
        nbo_list = []
        nbo_list.append("$nbo")
        for key, value in nbo.items():
            nbo_list.append(f"   {key} = {value}")
        nbo_list.append("$end")
        return "\n".join(nbo_list)

    @staticmethod
    def svp_template(svp: dict) -> str:
        """
        Template for the $svp section.

        Args:
            svp: dict of SVP parameters, e.g.
            {"rhoiso": "0.001", "nptleb": "1202", "itrngr": "2", "irotgr": "2"}

        Returns:
            str: the $svp section. Note that all parameters will be concatenated onto
                 a single line formatted as a FORTRAN namelist. This is necessary
                 because the isodensity SS(V)PE model in Q-Chem calls a secondary code.
        """
        svp_list = []
        svp_list.append("$svp")
        param_list = [f"{_key}={value}" for _key, value in svp.items()]
        svp_list.append(", ".join(param_list))
        svp_list.append("$end")
        return "\n".join(svp_list)

    @staticmethod
    def geom_opt_template(geom_opt: dict) -> str:
        """
        Args:
            geom_opt ():

        Returns:
            (str) geom_opt parameters.
        """
        geom_opt_list = []
        geom_opt_list.append("$geom_opt")
        for key, value in geom_opt.items():
            geom_opt_list.append(f"   {key} = {value}")
        geom_opt_list.append("$end")
        return "\n".join(geom_opt_list)

    @staticmethod
    def cdft_template(cdft: list[list[dict]]) -> str:
        """
        Args:
            cdft: list of lists of dicts.

        Returns:
            (str)
        """
        cdft_list = []
        cdft_list.append("$cdft")
        for ii, state in enumerate(cdft):
            for constraint in state:
                types = constraint["types"]
                cdft_list.append(f"   {constraint['value']}")

                type_strings = []
                for typ in types:
                    if typ is None or typ.lower() in ["c", "charge"]:
                        type_strings.append("")
                    elif typ.lower() in ["s", "spin"]:
                        type_strings.append("s")
                    else:
                        raise ValueError("Invalid CDFT constraint type!")

                for coef, first, last, type_string in zip(
                    constraint["coefficients"], constraint["first_atoms"], constraint["last_atoms"], type_strings
                ):
                    if type_string != "":
                        cdft_list.append(f"   {coef} {first} {last} {type_string}")
                    else:
                        cdft_list.append(f"   {coef} {first} {last}")
            if len(cdft) != 1 and ii + 1 < len(state):
                cdft_list.append("--------------")

        # Ensure that you don't have a line indicating a state that doesn't exist
        if cdft_list[-1] == "--------------":
            del cdft_list[-1]

        cdft_list.append("$end")
        return "\n".join(cdft_list)

    @staticmethod
    def almo_template(almo_coupling: list[list[tuple[int, int]]]) -> str:
        """
        Args:
            almo: list of lists of int 2-tuples.

        Returns:
            (str)
        """
        almo_list = []
        almo_list.append("$almo_coupling")

        # ALMO coupling calculations always involve 2 states
        if len(almo_coupling) != 2:
            raise ValueError("ALMO coupling calculations require exactly two states!")

        state_1 = almo_coupling[0]
        state_2 = almo_coupling[1]

        for frag in state_1:
            # Casting to int probably unnecessary, given type hint
            # Doesn't hurt, though
            almo_list.append(f"   {int(frag[0])} {int(frag[1])}")
        almo_list.append("   --")
        for frag in state_2:
            almo_list.append(f"   {int(frag[0])} {int(frag[1])}")

        almo_list.append("$end")
        return "\n".join(almo_list)

    @staticmethod
    def pcm_nonels_template(pcm_nonels: dict) -> str:
        """
        Template for the $pcm_nonels section.

        Arg
            pcm_nonels: dict of CMIRS parameters, e.g.
            {
                "a": "-0.006736",
                "b": "0.032698",
                "c": "-1249.6",
                "d": "-21.405",
                "gamma": "3.7",
                "solvrho": "0.05",
                "delta": 7,
                "gaulag_n": 40,
            }

        Returns:
            (str)
        """
        pcm_nonels_list = []
        pcm_nonels_list.append("$pcm_nonels")
        for key, value in pcm_nonels.items():
            # if the value is None, don't write it to output
            if value is not None:
                pcm_nonels_list.append(f"   {key} {value}")
        pcm_nonels_list.append("$end")
        return "\n".join(pcm_nonels_list)

    @staticmethod
    def find_sections(string: str) -> list:
        """
        Find sections in the string.

        Args:
            string (str): String

        Returns:
            List of sections.
        """
        patterns = {"sections": r"^\s*?\$([a-z_]+)", "multiple_jobs": r"(@@@)"}
        matches = read_pattern(string, patterns)
        # list of the sections present
        sections = [val[0] for val in matches["sections"]]
        # remove end from sections
        sections = [sec for sec in sections if sec != "end"]
        # this error should be replaced by a multi job read function when it is added
        if "multiple_jobs" in matches:
            raise ValueError("Output file contains multiple qchem jobs please parse separately")
        if "molecule" not in sections:
            raise ValueError("Output file does not contain a molecule section")
        if "rem" not in sections:
            raise ValueError("Output file does not contain a rem section")
        return sections

    @staticmethod
    def read_molecule(string: str) -> Molecule | list[Molecule] | Literal["read"]:
        """
        Read molecule from string.

        Args:
            string (str): String

        Returns:
            Molecule
        """
        charge = spin_mult = None
        patterns = {
            "read": r"^\s*\$molecule\n\s*(read)",
            "charge": r"^\s*\$molecule\n\s*((?:\-)*\d+)\s+\d+",
            "spin_mult": r"^\s*\$molecule\n\s(?:\-)*\d+\s*((?:\-)*\d+)",
            "fragment": r"^\s*\$molecule\n\s*(?:\-)*\d+\s+\d+\s*\n\s*(\-\-)",
        }
        matches = read_pattern(string, patterns)
        if "read" in matches:
            return "read"
        if "charge" in matches:
            charge = float(matches["charge"][0][0])
        if "spin_mult" in matches:
            spin_mult = int(matches["spin_mult"][0][0])
        multi_mol = "fragment" in matches

        if not multi_mol:
            header = r"^\s*\$molecule\n\s*(?:\-)*\d+\s+(?:\-)*\d+"
            row = r"\s*([A-Za-z]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer = r"^\$end"
            mol_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
            species = [val[0] for val in mol_table[0]]
            coords = [[float(val[1]), float(val[2]), float(val[3])] for val in mol_table[0]]
            if charge is None:
                mol = Molecule(species=species, coords=coords)
            else:
                mol = Molecule(species=species, coords=coords, charge=charge, spin_multiplicity=spin_mult)
            return mol

        header = r"\s*(?:\-)*\d+\s+(?:\-)*\d+"
        row = r"\s*([A-Za-z]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer = r"(:?(:?\-\-)|(:?\$end))"

        molecules = []

        patterns = {"charge_spin": r"\s*\-\-\s*([\-0-9]+)\s+([\-0-9]+)"}
        matches = read_pattern(string, patterns)

        mol_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        for match, table in zip(matches.get("charge_spin"), mol_table):
            charge = int(match[0])
            spin = int(match[1])
            species = [val[0] for val in table]
            coords = [[float(val[1]), float(val[2]), float(val[3])] for val in table]
            mol = Molecule(species=species, coords=coords, charge=charge, spin_multiplicity=spin)
            molecules.append(mol)

        return molecules

    @staticmethod
    def read_rem(string: str) -> dict:
        """
        Parse rem from string.

        Args:
            string (str): String

        Returns:
            (dict) rem
        """
        header = r"^\s*\$rem"
        row = r"\s*([a-zA-Z\_\d]+)\s*=?\s*(\S+)"
        footer = r"^\s*\$end"
        rem_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        return dict(rem_table[0])

    @staticmethod
    def read_opt(string: str) -> dict[str, list]:
        """
        Read opt section from string.

        Args:
            string (str): String

        Returns:
            (dict) Opt section
        """
        patterns = {
            "CONSTRAINT": r"^\s*CONSTRAINT",
            "FIXED": r"^\s*FIXED",
            "DUMMY": r"^\s*DUMMY",
            "CONNECT": r"^\s*CONNECT",
        }
        opt_matches = read_pattern(string, patterns)
        opt_sections = list(opt_matches)
        opt = {}
        if "CONSTRAINT" in opt_sections:
            c_header = r"^\s*CONSTRAINT\n"
            c_row = r"(\w.*)\n"
            c_footer = r"^\s*ENDCONSTRAINT\n"
            c_table = read_table_pattern(string, header_pattern=c_header, row_pattern=c_row, footer_pattern=c_footer)
            opt["CONSTRAINT"] = [val[0] for val in c_table[0]]
        if "FIXED" in opt_sections:
            f_header = r"^\s*FIXED\n"
            f_row = r"(\w.*)\n"
            f_footer = r"^\s*ENDFIXED\n"
            f_table = read_table_pattern(
                string,
                header_pattern=f_header,
                row_pattern=f_row,
                footer_pattern=f_footer,
            )
            opt["FIXED"] = [val[0] for val in f_table[0]]
        if "DUMMY" in opt_sections:
            d_header = r"^\s*DUMMY\n"
            d_row = r"(\w.*)\n"
            d_footer = r"^\s*ENDDUMMY\n"
            d_table = read_table_pattern(
                string,
                header_pattern=d_header,
                row_pattern=d_row,
                footer_pattern=d_footer,
            )
            opt["DUMMY"] = [val[0] for val in d_table[0]]
        if "CONNECT" in opt_sections:
            cc_header = r"^\s*CONNECT\n"
            cc_row = r"(\w.*)\n"
            cc_footer = r"^\s*ENDCONNECT\n"
            cc_table = read_table_pattern(
                string,
                header_pattern=cc_header,
                row_pattern=cc_row,
                footer_pattern=cc_footer,
            )
            opt["CONNECT"] = [val[0] for val in cc_table[0]]
        return opt

    @staticmethod
    def read_pcm(string: str) -> dict:
        """
        Read pcm parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) PCM parameters
        """
        header = r"^\s*\$pcm"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        pcm_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if not pcm_table:
            print("No valid PCM inputs found. Note that there should be no '=' characters in PCM input lines.")
            return {}

        return dict(pcm_table[0])

    @staticmethod
    def read_vdw(string: str) -> tuple[str, dict]:
        """
        Read van der Waals parameters from string.

        Args:
            string (str): String

        Returns:
            (str, dict) vdW mode ('atomic' or 'sequential') and dict of van der Waals radii.
        """
        header = r"^\s*\$van_der_waals"
        row = r"[^\d]*(\d+).?(\d+.\d+)?.*"
        footer = r"^\s*\$end"
        vdw_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if not vdw_table:
            print("No valid vdW inputs found. Note that there should be no '=' characters in vdW input lines.")
            return "", {}

        mode = "sequential" if vdw_table[0][0][0] == 2 else "atomic"

        return mode, dict(vdw_table[0][1:])

    @staticmethod
    def read_solvent(string: str) -> dict:
        """
        Read solvent parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) Solvent parameters
        """
        header = r"^\s*\$solvent"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        solvent_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if not solvent_table:
            print("No valid solvent inputs found. Note that there should be no '=' characters in solvent input lines.")
            return {}

        return dict(solvent_table[0])

    @staticmethod
    def read_smx(string: str) -> dict:
        """
        Read smx parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) SMX parameters.
        """
        header = r"^\s*\$smx"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        smx_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if not smx_table:
            print("No valid smx inputs found. Note that there should be no '=' characters in smx input lines.")
            return {}
        smx = {}
        for key, val in smx_table[0]:
            smx[key] = val
        if smx["solvent"] == "tetrahydrofuran":
            smx["solvent"] = "thf"
        # Q-Chem bug, see https://talk.q-chem.com/t/smd-unrecognized-solvent/204
        elif smx["solvent"] == "dimethyl sulfoxide":
            smx["solvent"] = "dmso"
        return smx

    @staticmethod
    def read_scan(string: str) -> dict[str, list]:
        """
        Read scan section from a string.

        Args:
            string: String to be parsed

        Returns:
            Dict representing Q-Chem scan section
        """
        header = r"^\s*\$scan"
        row = r"\s*(stre|bend|tors|STRE|BEND|TORS)\s+((?:[\-\.0-9]+\s*)+)"
        footer = r"^\s*\$end"
        scan_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if scan_table == []:
            print("No valid scan inputs found. Note that there should be no '=' characters in scan input lines.")
            return {}

        stre = []
        bend = []
        tors = []
        for row in scan_table[0]:
            if row[0].lower() == "stre":
                stre.append(row[1].replace("\n", "").rstrip())
            elif row[0].lower() == "bend":
                bend.append(row[1].replace("\n", "").rstrip())
            elif row[0].lower() == "tors":
                tors.append(row[1].replace("\n", "").rstrip())

        if len(stre) + len(bend) + len(tors) > 2:
            raise ValueError("No more than two variables are allows in the scan section!")

        return {"stre": stre, "bend": bend, "tors": tors}

    @staticmethod
    def read_plots(string: str) -> dict:
        """
        Read plots parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) plots parameters.
        """
        header = r"^\s*\$plots"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        plots_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if plots_table == []:
            print("No valid plots inputs found. Note that there should be no '=' characters in plots input lines.")
            return {}
        plots = {}
        for key, val in plots_table[0]:
            plots[key] = val
        return plots

    @staticmethod
    def read_nbo(string: str) -> dict:
        """
        Read nbo parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) nbo parameters.
        """
        header = r"^\s*\$nbo"
        row = r"\s*([a-zA-Z\_\d]+)\s*=?\s*(\S+)"
        footer = r"^\s*\$end"
        nbo_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if nbo_table == []:
            print("No valid nbo inputs found.")
            return {}
        nbo = {}
        for key, val in nbo_table[0]:
            nbo[key] = val
        return nbo

    @staticmethod
    def read_geom_opt(string: str) -> dict:
        """
        Read geom_opt parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) geom_opt parameters.
        """
        header = r"^\s*\$geom_opt"
        row = r"\s*([a-zA-Z\_]+)\s*=?\s*(\S+)"
        footer = r"^\s*\$end"
        geom_opt_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if geom_opt_table == []:
            print("No valid geom_opt inputs found.")
            return {}
        geom_opt = {}
        for key, val in geom_opt_table[0]:
            geom_opt[key] = val
        return geom_opt

    @staticmethod
    def read_cdft(string: str) -> list[list[dict]]:
        """
        Read cdft parameters from string.

        Args:
            string (str): String

        Returns:
            list[list[dict]]: cdft parameters
        """
        pattern_sec = {
            "full_section": r"\$cdft((:?(:?\s*[0-9\.\-]+\s+[0-9]+\s+[0-9]+(:?\s+[A-Za-z]+)?\s*\n)+|"
            r"(:?\s*[0-9\.\-]+\s*\n)|(:?\s*\-+\s*\n))+)\$end"
        }

        pattern_const = {
            "constraint": r"\s*([\-\.0-9]+)\s*\n((?:\s*(?:[\-\.0-9]+)\s+(?:\d+)\s+(?:\d+)(?:\s+[A-Za-z]+)?\s*)+)"
        }

        section = read_pattern(string, pattern_sec)["full_section"]
        if len(section) == 0:
            print("No valid cdft inputs found.")
            return []

        cdft = []
        section = section[0][0]
        states = re.split(r"\-{2,25}", section)
        for state in states:
            state_list = []
            const_out = list(read_pattern(state, pattern_const).get("constraint"))
            if len(const_out) == 0:
                continue
            for const in const_out:
                const_dict = {
                    "value": float(const[0]),
                    "coefficients": [],
                    "first_atoms": [],
                    "last_atoms": [],
                    "types": [],
                }  # type: ignore
                subconsts = const[1].strip().split("\n")
                for subconst in subconsts:
                    tokens = subconst.split()
                    const_dict["coefficients"].append(float(tokens[0]))  # type: ignore
                    const_dict["first_atoms"].append(int(tokens[1]))  # type: ignore
                    const_dict["last_atoms"].append(int(tokens[2]))  # type: ignore
                    if len(tokens) > 3:
                        const_dict["types"].append(tokens[3])  # type: ignore
                    else:
                        const_dict["types"].append(None)  # type: ignore

                state_list.append(const_dict)

            cdft.append(state_list)

        return cdft

    @staticmethod
    def read_almo(string: str) -> list[list[tuple[int, int]]]:
        """
        Read ALMO coupling parameters from string.

        Args:
            string (str): String

        Returns:
            (list of lists of int 2-tuples) almo_coupling parameters
        """
        pattern = {
            "key": r"\$almo_coupling\s*\n((?:\s*[\-0-9]+\s+[\-0-9]+\s*\n)+)\s*\-\-"
            r"((?:\s*[\-0-9]+\s+[\-0-9]+\s*\n)+)\s*\$end"
        }

        section = read_pattern(string, pattern)["key"]

        if len(section) == 0:
            print("No valid almo inputs found.")
            return []

        section = section[0]

        almo_coupling = [[], []]  # type: ignore

        state_1 = section[0]
        for line in state_1.strip().split("\n"):
            contents = line.split()
            almo_coupling[0].append((int(contents[0]), int(contents[1])))

        state_2 = section[1]
        for line in state_2.strip().split("\n"):
            contents = line.split()
            almo_coupling[1].append((int(contents[0]), int(contents[1])))

        return almo_coupling

    @staticmethod
    def read_svp(string: str) -> dict:
        """Read svp parameters from string."""
        header = r"^\s*\$svp"
        row = r"(\w.*)\n"
        footer = r"^\s*\$end"
        svp_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if svp_table == []:
            print("No valid svp inputs found.")
            return {}
        svp_list = svp_table[0][0][0].split(", ")
        svp_dict = {}
        for s in svp_list:
            svp_dict[s.split("=")[0]] = s.split("=")[1]
        return svp_dict

    @staticmethod
    def read_pcm_nonels(string: str) -> dict:
        """
        Read pcm_nonels parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) PCM parameters
        """
        header = r"^\s*\$pcm_nonels"
        row = r"\s*([a-zA-Z\_]+)\s+(.+)"
        footer = r"^\s*\$end"
        pcm_nonels_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if not pcm_nonels_table:
            print(
                "No valid $pcm_nonels inputs found. Note that there should be no '=' "
                "characters in $pcm_nonels input lines."
            )
            return {}

        return dict(pcm_nonels_table[0])
