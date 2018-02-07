# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re
import logging
import os

from monty.io import zopen
from monty.json import MSONable
from monty.re import regrep

from temp_utils import new_read_table_pattern

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
        self.text=""
        with zopen(filename, 'rt') as f:
            self.text = f.read()

        # Check if output file contains multiple output files. If so, print an error message and exit
        self.read_pattern({"multiple_outputs": r"Job\s+\d+\s+of\s+(\d+)\s+"}, terminate_on_match=True)
        if not (self.data.get('multiple_outputs') == [] or self.data.get('multiple_outputs') == [['1']]):
            print "ERROR: multiple calculation outputs found in file "+filename+". Please instead call QCOutput.mulitple_outputs_from_file(QCOutput,'"+filename+"')"
            print "Exiting..."
            exit()

        # Check if calculation finished. If not, proceed with caution
        self.read_pattern({"completion": r"Thank you very much for using Q-Chem.\s+Have a nice day."})
        if not self.data.get('completion'):
            print "WARNING: calculation did not reach successful completion"

        # Check if calculation is unrestricted
        self.read_pattern({"unrestricted": r"A(?:n)*\sunrestricted[\s\w\-]+SCF\scalculation\swill\sbe"}, terminate_on_match=True)

        # Check if calculation uses GEN_SCFMAN
        self.read_pattern({"using_GEN_SCFMAN": r"\s+GEN_SCFMAN: A general SCF calculation manager"}, terminate_on_match=True)

        # Parse the SCF
        if self.data.get('using_GEN_SCFMAN',[]):
            self._read_GEN_SCFMAN()
        else:
            self._read_SCF()

        # Parse the Mulliken charges
        if self.data.get('unrestricted',[]):
            self._read_unrestricted_mulliken()
        else:
            self._read_restricted_mulliken()

        # Parse the final energy
        self.read_pattern({"final_energy": r"Final\senergy\sis\s+([\d\-\.]+)"})

        # Parse the S2 values in the case of an unrestricted calculation
        if self.data.get('unrestricted',[]):
            self.read_pattern({"S2": r"<S\^2>\s=\s+([\d\-\.]+)"})

        # Check if the calculation is a geometry optimization. If so, parse the relevant output
        self.read_pattern({"optimization": r"(?i)\s*job(?:_)*type\s+=\s+opt"}, terminate_on_match=True)
        if self.data.get('optimization',[]):
            self.read_pattern({'energy_trajectory': r"\sEnergy\sis\s+([\d\-\.]+)"})
            self._read_optimized_geometry()

        # Check if the calculation is a frequency analysis. If so, parse the relevant output
        self.read_pattern({"frequency_job": r"(?i)\s*job(?:_)*type\s+=\s+freq"}, terminate_on_match=True)
        if self.data.get('frequency_job',[]):
            self.read_pattern({"frequencies": r"\s*Frequency:\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)",
                               "enthalpy": r"\s*Total Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                               "entropy": r"\s*Total Entropy:\s+([\d\-\.]+)\s+cal/mol\.K"})


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
            text = re.split('Job\s+\d+\s+of\s+\d+\s+',f.read())
        for ii in range(len(text)):
            temp = open(filename+'.'+str(ii), 'w')
            temp.write(text[ii])
            temp.close()
            tempOutput = cls(filename+'.'+str(ii))
            to_return.append(tempOutput)
            if not keep_sub_files:
                os.remove(filename+'.'+str(ii))
        return to_return
        

    def _read_GEN_SCFMAN(self):
        """
        Parses all GEN_SCFMANs
        """
        header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*(?:RMS Gradient)*\s+\-+(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+\-+)*\n"
        table_pattern = r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s+([\d\-\.]+)\s+([\d\-\.]+)e([\d\-\.\+]+)(?:\s+Convergence criterion met)*(?:\s+Preconditoned Steepest Descent)*(?:\s+Roothaan Step)*(?:\s+(?:Normal\s+)*BFGS [Ss]tep)*(?:\s+LineSearch Step)*(?:\s+Line search: overstep)*(?:\s+Descent step)*" 
        footer_pattern = r"^\s*\-+\n"
        self.data["GEN_SCFMAN"] = new_read_table_pattern(self.text,header_pattern,table_pattern,footer_pattern)



    def _read_SCF(self):
        """
        Parses all old-style SCFs
        """
        #Check if the SCF failed to converge and set the footer accordingly
        self.read_pattern({"failed_to_converge": r"SCF failed to converge"}, terminate_on_match=True)
        if self.data.get("failed_to_converge",[]):
            footer_pattern = r"^\s*\d+\s+[\d\-\.]+\s+[\d\-\.]+E[\d\-\.]+\s+Convergence\s+failure\n"
        else:
            footer_pattern = r"^\s*\-+\n"
        header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+DIIS Error\s+\-+\n"
        table_pattern = r"\s*\d+\s+([\d\-\.]+)\s+([\d\-\.]+)E([\d\-\.]+)(?:\s*\n\s*cpu\s+[\d\-\.]+\swall\s+[\d\-\.]+)*(?:\nin dftxc\.C, eleTot sum is:[\d\-\.]+, tauTot is\:[\d\-\.]+)*(?:\s+Convergence criterion met)*(?:\s+Done RCA\. Switching to DIIS)*"
        
        self.data["SCF"] = new_read_table_pattern(self.text,header_pattern,table_pattern,footer_pattern)


    def _read_restricted_mulliken(self):
        """
        Parses Mulliken charges given a restricted SCF.
        """
        header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
        table_pattern = r"\s+\d+\s(\w+)\s+([\d\-\.]+)"
        footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

        self.data["restricted_Mulliken"] = new_read_table_pattern(self.text,header_pattern,table_pattern,footer_pattern)


    def _read_unrestricted_mulliken(self):
        """
        Parses Mulliken charges and spins given an unrestricted SCF. 
        """
        header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+Spin\s\(a\.u\.\)\s+\-+"
        table_pattern = r"\s+\d+\s(\w+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

        self.data["unrestricted_Mulliken"] = new_read_table_pattern(self.text,header_pattern,table_pattern,footer_pattern)


    def _read_optimized_geometry(self):
        """
        Parses optimized XYZ coordinates
        """
        header_pattern = r"\*+\s+OPTIMIZATION\s+CONVERGED\s+\*+\s+\*+\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+(\w+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Z-matrix Print:"

        self.data["optimized_geometry"] = new_read_table_pattern(self.text,header_pattern,table_pattern,footer_pattern)


    def read_pattern(self, patterns, reverse=False, terminate_on_match=False, postprocess=str):
        """
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
        matches = regrep(self.filename, patterns, reverse=reverse, terminate_on_match=terminate_on_match, postprocess=postprocess)
        # print matches
        for k in patterns.keys():
            self.data[k] = [i[0] for i in matches.get(k, [])]
