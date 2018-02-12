# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re
import logging
import os

from monty.io import zopen
from monty.json import MSONable

from .utils import read_table_pattern, read_pattern
"""
Classes for reading/manipulating/writing QChem ouput files.
"""

__author__ = "Samuel Blau, Brandon Woods, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class QCOutput(MSONable):
    """
    Data in a single QChem Calculations

    Args:
        filename (str): OUTCAR filename to parse.
    """

    def __init__(self, filename):
        self.filename = filename
        self.data = {}
        self.text = ""
        with zopen(filename, 'rt') as f:
            self.text = f.read()

        # Check if output file contains multiple output files. If so, print an error message and exit
        self.data["multiple_outputs"] = read_pattern(
            self.text, {
                "key": r"Job\s+\d+\s+of\s+(\d+)\s+"
            }, terminate_on_match=True).get('key')
        if not (self.data.get('multiple_outputs') == None or self.data.get('multiple_outputs') == [['1']]):
            print("ERROR: multiple calculation outputs found in file " + filename +
                  ". Please instead call QCOutput.mulitple_outputs_from_file(QCOutput,'" + filename + "')")
            print("Exiting...")
            exit()

        # Check if calculation finished. If not, proceed with caution
        self.data["completion"] = read_pattern(self.text, {
            "key": r"Thank you very much for using Q-Chem.\s+Have a nice day."
        }).get('key')
        # if not self.data.get('completion'):
        #     print("WARNING: calculation did not reach successful completion")

        # Check if calculation is unrestricted
        self.data["unrestricted"] = read_pattern(
            self.text, {
                "key": r"A(?:n)*\sunrestricted[\s\w\-]+SCF\scalculation\swill\sbe"
            }, terminate_on_match=True).get('key')

        # Check if calculation uses GEN_SCFMAN
        self.data["using_GEN_SCFMAN"] = read_pattern(
            self.text, {
                "key": r"\s+GEN_SCFMAN: A general SCF calculation manager"
            }, terminate_on_match=True).get('key')

        # Parse the SCF
        if self.data.get('using_GEN_SCFMAN', []):
            self._read_GEN_SCFMAN()
        else:
            self._read_SCF()

        # Parse the Mulliken charges
        if self.data.get('unrestricted', []):
            self._read_unrestricted_mulliken()
        else:
            self._read_restricted_mulliken()

        # Parse the final energy
        self.data["final_energy"] = read_pattern(self.text, {"key": r"Final\senergy\sis\s+([\d\-\.]+)"}).get('key')

        # Parse the S2 values in the case of an unrestricted calculation
        if self.data.get('unrestricted', []):
            self.data["S2"] = read_pattern(self.text, {"key": r"<S\^2>\s=\s+([\d\-\.]+)"}).get('key')

        # Check if the calculation is a geometry optimization. If so, parse the relevant output
        self.data["optimization"] = read_pattern(self.text, {"key": r"(?i)\s*job(?:_)*type\s+=\s+opt"}).get('key')
        if self.data.get('optimization', []):
            self.data["energy_trajectory"] = read_pattern(self.text, {"key": r"\sEnergy\sis\s+([\d\-\.]+)"}).get('key')
            self._read_optimized_geometry()

        # Check if the calculation is a frequency analysis. If so, parse the relevant output
        self.data["frequency_job"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s+=\s+freq"
            }, terminate_on_match=True).get('key')
        if self.data.get('frequency_job', []):
            temp_dict = read_pattern(
                self.text, {
                    "frequencies": r"\s*Frequency:\s+([\d\-\.]+)(?:\s+([\d\-\.]+)(?:\s+([\d\-\.]+))*)*",
                    "enthalpy": r"\s*Total Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                    "entropy": r"\s*Total Entropy:\s+([\d\-\.]+)\s+cal/mol\.K"
                })
            for key in temp_dict:
                self.data[key] = temp_dict.get(key)

    @staticmethod
    def multiple_outputs_from_file(cls, filename, keep_sub_files=True):
        """
            Parses a QChem output file with multiple calculations
            1.) Seperates the output into sub-files
                e.g. qcout -> qcout.0, qcout.1, qcout.2 ... qcout.N
                a.) Find delimeter for multiple calcualtions
                b.) Make seperate output sub-files
            2.) Creates seperate QCCalcs for each one from the sub-files
        """
        to_return = []
        with zopen(filename, 'rt') as f:
            text = re.split('\s*(?:Running\s+)*Job\s+\d+\s+of\s+\d+\s+', f.read())
        if text[0] == '':
            text = text[1:]
        for i, sub_text in enumerate(text):
            temp = open(filename + '.' + str(i), 'w')
            temp.write(sub_text)
            temp.close()
            tempOutput = cls(filename + '.' + str(i))
            to_return.append(tempOutput)
            if not keep_sub_files:
                os.remove(filename + '.' + str(i))
        return to_return

    def _read_GEN_SCFMAN(self):
        """
        Parses all GEN_SCFMANs
        """
        header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*(?:RMS Gradient)*\s+\-+(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n"
        table_pattern = r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s+([\d\-\.]+)\s+([\d\-\.]+)e([\d\-\.\+]+)(?:\s+Convergence criterion met)*(?:\s+Preconditoned Steepest Descent)*(?:\s+Roothaan Step)*(?:\s+(?:Normal\s+)*BFGS [Ss]tep)*(?:\s+LineSearch Step)*(?:\s+Line search: overstep)*(?:\s+Descent step)*"
        footer_pattern = r"^\s*\-+\n"

        self.data["GEN_SCFMAN"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)

    def _read_SCF(self):
        """
        Parses all old-style SCFs. Starts by checking if the SCF failed to converge and setting the footer accordingly.
        """
        self.data["SCF_failed_to_converge"] = read_pattern(
            self.text, {
                "key": r"SCF failed to converge"
            }, terminate_on_match=True).get('key')
        if self.data.get("SCF_failed_to_converge", []):
            footer_pattern = r"^\s*\d+\s*[\d\-\.]+\s+[\d\-\.]+E[\d\-\.]+\s+Convergence\s+failure\n"
        else:
            footer_pattern = r"^\s*\-+\n"
        header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+DIIS Error\s+\-+\n"
        table_pattern = r"\s*\d+\s*([\d\-\.]+)\s+([\d\-\.]+)E([\d\-\.\+]+)(?:\s*\n\s*cpu\s+[\d\-\.]+\swall\s+[\d\-\.]+)*(?:\nin dftxc\.C, eleTot sum is:[\d\-\.]+, tauTot is\:[\d\-\.]+)*(?:\s+Convergence criterion met)*(?:\s+Done RCA\. Switching to DIIS)*(?:\n\s*Warning: not using a symmetric Q)*(?:\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+(?:\s*\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+)*)*"

        self.data["SCF"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)

    def _read_restricted_mulliken(self):
        """
        Parses Mulliken charges given a restricted SCF.
        """
        header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
        table_pattern = r"\s+\d+\s(\w+)\s+([\d\-\.]+)"
        footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

        self.data["restricted_Mulliken"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)

    def _read_unrestricted_mulliken(self):
        """
        Parses Mulliken charges and spins given an unrestricted SCF.
        """
        header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+Spin\s\(a\.u\.\)\s+\-+"
        table_pattern = r"\s+\d+\s(\w+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

        self.data["unrestricted_Mulliken"] = read_table_pattern(self.text, header_pattern, table_pattern,
                                                                footer_pattern)

    def _read_optimized_geometry(self):
        """
        Parses optimized XYZ coordinates. If not present, parses optimized Z-matrix.
        """
        header_pattern = r"\*+\s+OPTIMIZATION\s+CONVERGED\s+\*+\s+\*+\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+(\w+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Z-matrix Print:"

        self.data["optimized_geometry"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)

        if self.data.get('optimized_geometry') == []:
            header_pattern = r"^\s+\*+\s+OPTIMIZATION CONVERGED\s+\*+\s+\*+\s+Z-matrix\s+Print:\s+\$molecule\s+[\d\-]+\s+[\d\-]+\n"
            table_pattern = r"\s*(\w+)(?:\s+(\d+)\s+([\d\-\.]+)(?:\s+(\d+)\s+([\d\-\.]+)(?:\s+(\d+)\s+([\d\-\.]+))*)*)*(?:\s+0)*"
            footer_pattern = r"^\$end\n"

            self.data["optimized_zmat"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
