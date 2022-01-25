# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Classes for reading/manipulating/writing QChem input files.
"""
import logging
from typing import Dict, List, Literal, Optional, Tuple, Union

from monty.io import zopen

from pymatgen.core import Molecule
from pymatgen.io.core import InputFile

from .utils import lower_and_check_unique, read_pattern, read_table_pattern

__author__ = "Brandon Wood, Samuel Blau, Shyam Dwaraknath, Julian Self, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__email__ = "b.wood@berkeley.edu"
__credits__ = "Xiaohui Qu"

logger = logging.getLogger(__name__)


class QCInput(InputFile):
    """
    An object representing a QChem input file. QCInput attributes represent different sections of a QChem input file.
    To add a new section one needs to modify __init__, __str__, from_sting and add staticmethods
    to read and write the new section i.e. section_template and read_section. By design, there is very little (or no)
    checking that input parameters conform to the appropriate QChem format, this responsible lands on the user or a
    separate error handling software.
    """

    def __init__(
        self,
        molecule: Union[Molecule, Literal["read"]],
        rem: Dict,
        opt: Optional[Dict[str, List]] = None,
        pcm: Optional[Dict] = None,
        solvent: Optional[Dict] = None,
        smx: Optional[Dict] = None,
        scan: Optional[Dict[str, List]] = None,
        van_der_waals: Optional[Dict[str, float]] = None,
        vdw_mode: str = "atomic",
        plots: Optional[Dict] = None,
        nbo: Optional[Dict] = None,
    ):
        """
        Args:
            molecule (pymatgen Molecule object or "read"):
                Input molecule. molecule can be set as either a pymatgen Molecule object or as the str "read".
                "read" can be used in multi_job QChem input files where the molecule is read in from the
                previous calculation.
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
                    A dictionary of all the input parameters for the plots section of QChem input file.
            nbo (dict):
                    A dictionary of all the input parameters for the nbo section of QChem input file.

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
        if "method" not in self.rem:
            if "exchange" not in self.rem:
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

    def get_string(self):
        """
        Return a string representation of an entire input file.
        """
        return self.__str__()

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
        return "\n".join(combined_list)

    @staticmethod
    def multi_job_string(job_list: List["QCInput"]) -> str:
        """
        Args:
            job_list (): List of jobs

        Returns:
            (str) String representation of multi job input file.
        """
        multi_job_string = ""
        for i, job_i in enumerate(job_list):
            if i < len(job_list) - 1:
                multi_job_string += job_i.__str__() + "\n@@@\n\n"
            else:
                multi_job_string += job_i.__str__()
        return multi_job_string

    @classmethod
    def from_string(cls, string: str) -> "QCInput":
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
        opt = None
        pcm = None
        solvent = None
        smx = None
        scan = None
        vdw = None
        vdw_mode = "atomic"
        plots = None
        nbo = None
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
        )

    @staticmethod
    def write_multi_job_file(job_list: List["QCInput"], filename: str):
        """
        Write a multijob file.

        Args:
            job_list (): List of jobs.
            filename (): Filename
        """
        with zopen(filename, "wt") as f:
            f.write(QCInput.multi_job_string(job_list))

    @classmethod
    def from_multi_jobs_file(cls, filename: str) -> List["QCInput"]:
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
            input_list = [cls.from_string(i) for i in multi_job_strings]
            return input_list

    @staticmethod
    def molecule_template(molecule: Union[Molecule, Literal["read"]]) -> str:
        """
        Args:
            molecule (Molecule): molecule

        Returns:
            (str) Molecule template.
        """
        # todo: add ghost atoms
        mol_list = []
        mol_list.append("$molecule")
        if isinstance(molecule, str):
            if molecule == "read":
                mol_list.append(" read")
            else:
                raise ValueError('The only acceptable text value for molecule is "read"')
        else:
            mol_list.append(f" {int(molecule.charge)} {molecule.spin_multiplicity}")
            for site in molecule.sites:
                mol_list.append(
                    " {atom}     {x: .10f}     {y: .10f}     {z: .10f}".format(
                        atom=site.species_string, x=site.x, y=site.y, z=site.z
                    )
                )
        mol_list.append("$end")
        return "\n".join(mol_list)

    @staticmethod
    def rem_template(rem: Dict) -> str:
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
    def opt_template(opt: Dict[str, List]) -> str:
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
    def pcm_template(pcm: Dict) -> str:
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
    def solvent_template(solvent: Dict) -> str:
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
    def smx_template(smx: Dict) -> str:
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
            else:
                smx_list.append(f"   {key} {value}")
        smx_list.append("$end")
        return "\n".join(smx_list)

    @staticmethod
    def scan_template(scan: Dict[str, List]) -> str:
        """
        Args:
            scan (dict): Dictionary with scan section information.
                Ex: {"stre": ["3 6 1.5 1.9 0.1"], "tors": ["1 2 3 4 -180 180 15"]}

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
    def van_der_waals_template(radii: Dict[str, float], mode: str = "atomic") -> str:
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
                **NOTE: keys must be given as strings even though they are numbers!**

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
            raise ValueError(f"Invalid value {mode} given for 'mode' kwarg.")

        for num, radius in radii.items():
            vdw_list.append(f"   {num} {radius}")
        vdw_list.append("$end")
        return "\n".join(vdw_list)

    @staticmethod
    def plots_template(plots: Dict) -> str:
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
    def nbo_template(nbo: Dict) -> str:
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
    def find_sections(string: str) -> List:
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
        if "multiple_jobs" in matches.keys():
            raise ValueError("Output file contains multiple qchem jobs please parse separately")
        if "molecule" not in sections:
            raise ValueError("Output file does not contain a molecule section")
        if "rem" not in sections:
            raise ValueError("Output file does not contain a rem section")
        return sections

    @staticmethod
    def read_molecule(string: str) -> Union[Molecule, Literal["read"]]:
        """
        Read molecule from string.

        Args:
            string (str): String

        Returns:
            Molecule
        """
        charge = None
        spin_mult = None
        patterns = {
            "read": r"^\s*\$molecule\n\s*(read)",
            "charge": r"^\s*\$molecule\n\s*((?:\-)*\d+)\s+\d",
            "spin_mult": r"^\s*\$molecule\n\s(?:\-)*\d+\s*(\d)",
        }
        matches = read_pattern(string, patterns)
        if "read" in matches.keys():
            return "read"
        if "charge" in matches.keys():
            charge = float(matches["charge"][0][0])
        if "spin_mult" in matches.keys():
            spin_mult = int(matches["spin_mult"][0][0])
        header = r"^\s*\$molecule\n\s*(?:\-)*\d+\s*\d"
        row = r"\s*((?i)[a-z]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer = r"^\$end"
        mol_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        species = [val[0] for val in mol_table[0]]
        coords = [[float(val[1]), float(val[2]), float(val[3])] for val in mol_table[0]]
        if charge is None:
            mol = Molecule(species=species, coords=coords)
        else:
            mol = Molecule(species=species, coords=coords, charge=charge, spin_multiplicity=spin_mult)
        return mol

    @staticmethod
    def read_rem(string: str) -> Dict:
        """
        Parse rem from string.

        Args:
            string (str): String

        Returns:
            (dict) rem
        """
        header = r"^\s*\$rem"
        row = r"\s*([a-zA-Z\_]+)\s*=?\s*(\S+)"
        footer = r"^\s*\$end"
        rem_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        return dict(rem_table[0])

    @staticmethod
    def read_opt(string: str) -> Dict[str, List]:
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
        opt_sections = list(opt_matches.keys())
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
    def read_pcm(string: str) -> Dict:
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
    def read_vdw(string: str) -> Tuple[str, Dict]:
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

        if vdw_table[0][0][0] == 2:
            mode = "sequential"
        else:
            mode = "atomic"

        return mode, dict(vdw_table[0][1:])

    @staticmethod
    def read_solvent(string: str) -> Dict:
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
    def read_smx(string: str) -> Dict:
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
        return smx

    @staticmethod
    def read_scan(string: str) -> Dict[str, List]:
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
    def read_plots(string: str) -> Dict:
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
    def read_nbo(string: str) -> Dict:
        """
        Read nbo parameters from string.

        Args:
            string (str): String

        Returns:
            (dict) nbo parameters.
        """
        header = r"^\s*\$nbo"
        row = r"\s*([a-zA-Z\_]+)\s*=?\s*(\S+)"
        footer = r"^\s*\$end"
        nbo_table = read_table_pattern(string, header_pattern=header, row_pattern=row, footer_pattern=footer)
        if nbo_table == []:
            print("No valid nbo inputs found.")
            return {}
        nbo = {}
        for key, val in nbo_table[0]:
            nbo[key] = val
        return nbo
