# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import logging
from monty.json import MSONable
from monty.io import zopen
from pymatgen.core import Molecule
from .utils import read_table_pattern, read_pattern, lower_and_check_unique

# Classes for reading/manipulating/writing QChem ouput files.


__author__ = "Brandon Wood, Samuel Blau, Shyam Dwaraknath, Julian Self"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__email__ = "b.wood@berkeley.edu"
__credits__ = "Xiaohui Qu"

logger = logging.getLogger(__name__)


class QCInput(MSONable):
    """
    An object representing a QChem input file. QCInput attributes represent different sections of a QChem input file.
    To add a new section one needs to modify __init__, __str__, from_sting and add staticmethods
    to read and write the new section i.e. section_template and read_section. By design, there is very little (or no)
    checking that input parameters conform to the appropriate QChem format, this responsible lands on the user or a
    separate error handling software.

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
            values are a list of strings. Stings must be formatted as instructed by the QChem manual.
            The different opt sections are: CONSTRAINT, FIXED, DUMMY, and CONNECT
            Ex. opt = {"CONSTRAINT": ["tors 2 3 4 5 25.0", "tors 2 5 7 9 80.0"], "FIXED": ["2 XY"]}
    """

    def __init__(self, molecule, rem, opt=None, pcm=None, solvent=None, smx=None):
        self.molecule = molecule
        self.rem = lower_and_check_unique(rem)
        self.opt = opt
        self.pcm = lower_and_check_unique(pcm)
        self.solvent = lower_and_check_unique(solvent)
        self.smx = lower_and_check_unique(smx)

        # Make sure molecule is valid: either the string "read" or a pymatgen molecule object

        if isinstance(self.molecule, str):
            self.molecule = self.molecule.lower()
            if self.molecule != "read":
                raise ValueError(
                    'The only acceptable text value for molecule is "read"')
        elif not isinstance(self.molecule, Molecule):
            raise ValueError(
                "The molecule must either be the string 'read' or be a pymatgen Molecule object"
            )

        # Make sure rem is valid:
        #   - Has a basis
        #   - Has a method or DFT exchange functional
        #   - Has a valid job_type or jobtype

        valid_job_types = [
            "opt", "optimization", "sp", "freq", "frequency", "nmr"
        ]

        if "basis" not in self.rem:
            raise ValueError("The rem dictionary must contain a 'basis' entry")
        if "method" not in self.rem:
            if "exchange" not in self.rem:
                raise ValueError(
                    "The rem dictionary must contain either a 'method' entry or an 'exchange' entry"
                )
        if "job_type" not in self.rem:
            raise ValueError(
                "The rem dictionary must contain a 'job_type' entry")
        if self.rem.get("job_type").lower() not in valid_job_types:
            raise ValueError(
                "The rem dictionary must contain a valid 'job_type' entry")

        # Still to do:
        #   - Check that the method or functional is valid
        #   - Check that basis is valid
        #   - Check that basis is defined for all species in the molecule
        #   - Validity checks specific to job type?
        #   - Check OPT and PCM sections?

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
        return '\n'.join(combined_list)

    @staticmethod
    def multi_job_string(job_list):
        multi_job_string = str()
        for i, job_i in enumerate(job_list):
            if i < len(job_list) - 1:
                multi_job_string += job_i.__str__() + "\n@@@\n\n"
            else:
                multi_job_string += job_i.__str__()
        return multi_job_string

    @classmethod
    def from_string(cls, string):
        sections = cls.find_sections(string)
        molecule = cls.read_molecule(string)
        rem = cls.read_rem(string)
        # only molecule and rem are necessary everything else is checked
        opt = None
        pcm = None
        solvent = None
        smx=None
        if "opt" in sections:
            opt = cls.read_opt(string)
        if "pcm" in sections:
            pcm = cls.read_pcm(string)
        if "solvent" in sections:
            solvent = cls.read_solvent(string)
        if "smx" in sections:
            smx = cls.read_smx(string)
        return cls(molecule, rem, opt=opt, pcm=pcm, solvent=solvent, smx=smx)

    def write_file(self, filename):
        with zopen(filename, 'wt') as f:
            f.write(self.__str__())

    @staticmethod
    def write_multi_job_file(job_list, filename):
        with zopen(filename, 'wt') as f:
            f.write(QCInput.multi_job_string(job_list))

    @staticmethod
    def from_file(filename):
        with zopen(filename, 'rt') as f:
            return QCInput.from_string(f.read())

    @classmethod
    def from_multi_jobs_file(cls, filename):
        # returns a list of QCInput objects
        with zopen(filename, 'rt') as f:
            # the delimiter between QChem jobs is @@@
            multi_job_strings = f.read().split("@@@")
            # list of individual QChem jobs
            input_list = [cls.from_string(i) for i in multi_job_strings]
            return input_list

    @staticmethod
    def molecule_template(molecule):
        # todo: add ghost atoms
        mol_list = []
        mol_list.append("$molecule")
        if isinstance(molecule, str):
            if molecule == "read":
                mol_list.append(" read")
            else:
                raise ValueError('The only acceptable text value for molecule is "read"')
        else:
            mol_list.append(" {charge} {spin_mult}".format(
                charge=int(molecule.charge),
                spin_mult=molecule.spin_multiplicity))
            for site in molecule.sites:
                mol_list.append(
                    " {atom}     {x: .10f}     {y: .10f}     {z: .10f}".format(
                        atom=site.species_string, x=site.x, y=site.y,
                        z=site.z))
        mol_list.append("$end")
        return '\n'.join(mol_list)

    @staticmethod
    def rem_template(rem):
        rem_list = []
        rem_list.append("$rem")
        for key, value in rem.items():
            rem_list.append("   {key} = {value}".format(key=key, value=value))
        rem_list.append("$end")
        return '\n'.join(rem_list)

    @staticmethod
    def opt_template(opt):
        opt_list = []
        opt_list.append("$opt")
        # loops over all opt sections
        for key, value in opt.items():
            opt_list.append("{section}".format(section=key))
            # loops over all values within the section
            for i in value:
                opt_list.append("   {val}".format(val=i))
            opt_list.append("END{section}".format(section=key))
            opt_list.append("")
        # this deletes the empty space after the last section
        del opt_list[-1]
        opt_list.append("$end")
        return '\n'.join(opt_list)

    @staticmethod
    def pcm_template(pcm):
        pcm_list = []
        pcm_list.append("$pcm")
        for key, value in pcm.items():
            pcm_list.append("   {key} {value}".format(key=key, value=value))
        pcm_list.append("$end")
        return '\n'.join(pcm_list)

    @staticmethod
    def solvent_template(solvent):
        solvent_list = []
        solvent_list.append("$solvent")
        for key, value in solvent.items():
            solvent_list.append("   {key} {value}".format(
                key=key, value=value))
        solvent_list.append("$end")
        return '\n'.join(solvent_list)

    @staticmethod
    def smx_template(smx):
        smx_list = []
        smx_list.append("$smx")
        for key, value in smx.items():
            smx_list.append("   {key} {value}".format(
                key=key, value=value))
        smx_list.append("$end")
        return '\n'.join(smx_list)

    @staticmethod
    def find_sections(string):
        patterns = {"sections": r"^\s*?\$([a-z]+)", "multiple_jobs": r"(@@@)"}
        matches = read_pattern(string, patterns)
        # list of the sections present
        sections = [val[0] for val in matches["sections"]]
        # remove end from sections
        sections = [sec for sec in sections if sec != 'end']
        # this error should be replaced by a multi job read function when it is added
        if "multiple_jobs" in matches.keys():
            raise ValueError(
                "Output file contains multiple qchem jobs please parse separately"
            )
        if "molecule" not in sections:
            raise ValueError("Output file does not contain a molecule section")
        if "rem" not in sections:
            raise ValueError("Output file does not contain a rem section")
        return sections

    @staticmethod
    def read_molecule(string):
        charge = None
        spin_mult = None
        patterns = {
            "read": r"^\s*\$molecule\n\s*(read)",
            "charge": r"^\s*\$molecule\n\s*((?:\-)*\d+)\s+\d",
            "spin_mult": r"^\s*\$molecule\n\s*\d+\s*(\d)"
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
        mol_table = read_table_pattern(
            string,
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer)
        species = [val[0] for val in mol_table[0]]
        coords = [[float(val[1]), float(val[2]),
                   float(val[3])] for val in mol_table[0]]
        mol = Molecule(
            species=species,
            coords=coords,
            charge=charge,
            spin_multiplicity=spin_mult)
        return mol

    @staticmethod
    def read_rem(string):
        header = r"^\s*\$rem"
        row = r"\s*([a-zA-Z\_]+)\s*=?\s*(\S+)"
        footer = r"^\s*\$end"
        rem_table = read_table_pattern(
            string,
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer)
        rem = {key: val for key, val in rem_table[0]}
        return rem

    @staticmethod
    def read_opt(string):
        patterns = {
            "CONSTRAINT": r"^\s*CONSTRAINT",
            "FIXED": r"^\s*FIXED",
            "DUMMY": r"^\s*DUMMY",
            "CONNECT": r"^\s*CONNECT"
        }
        opt_matches = read_pattern(string, patterns)
        opt_sections = [key for key in opt_matches.keys()]
        opt = {}
        if "CONSTRAINT" in opt_sections:
            c_header = r"^\s*CONSTRAINT\n"
            c_row = r"(\w.*)\n"
            c_footer = r"^\s*ENDCONSTRAINT\n"
            c_table = read_table_pattern(
                string,
                header_pattern=c_header,
                row_pattern=c_row,
                footer_pattern=c_footer)
            opt["CONSTRAINT"] = [val[0] for val in c_table[0]]
        if "FIXED" in opt_sections:
            f_header = r"^\s*FIXED\n"
            f_row = r"(\w.*)\n"
            f_footer = r"^\s*ENDFIXED\n"
            f_table = read_table_pattern(
                string,
                header_pattern=f_header,
                row_pattern=f_row,
                footer_pattern=f_footer)
            opt["FIXED"] = [val[0] for val in f_table[0]]
        if "DUMMY" in opt_sections:
            d_header = r"^\s*DUMMY\n"
            d_row = r"(\w.*)\n"
            d_footer = r"^\s*ENDDUMMY\n"
            d_table = read_table_pattern(
                string,
                header_pattern=d_header,
                row_pattern=d_row,
                footer_pattern=d_footer)
            opt["DUMMY"] = [val[0] for val in d_table[0]]
        if "CONNECT" in opt_sections:
            cc_header = r"^\s*CONNECT\n"
            cc_row = r"(\w.*)\n"
            cc_footer = r"^\s*ENDCONNECT\n"
            cc_table = read_table_pattern(
                string,
                header_pattern=cc_header,
                row_pattern=cc_row,
                footer_pattern=cc_footer)
            opt["CONNECT"] = [val[0] for val in cc_table[0]]
        return opt

    @staticmethod
    def read_pcm(string):
        header = r"^\s*\$pcm"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        pcm_table = read_table_pattern(
            string,
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer)
        if pcm_table == []:
            print(
                "No valid PCM inputs found. Note that there should be no '=' chracters in PCM input lines."
            )
            return {}
        else:
            pcm = {key: val for key, val in pcm_table[0]}
            return pcm

    @staticmethod
    def read_solvent(string):
        header = r"^\s*\$solvent"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        solvent_table = read_table_pattern(
            string,
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer)
        if solvent_table == []:
            print(
                "No valid solvent inputs found. Note that there should be no '=' chracters in solvent input lines."
            )
            return {}
        else:
            solvent = {key: val for key, val in solvent_table[0]}
            return solvent

    @staticmethod
    def read_smx(string):
        header = r"^\s*\$smx"
        row = r"\s*([a-zA-Z\_]+)\s+(\S+)"
        footer = r"^\s*\$end"
        smx_table = read_table_pattern(
            string,
            header_pattern=header,
            row_pattern=row,
            footer_pattern=footer)
        if smx_table == []:
            print(
                "No valid smx inputs found. Note that there should be no '=' chracters in smx input lines."
            )
            return {}
        else:
            smx = {key: val for key, val in smx_table[0]}
            return smx