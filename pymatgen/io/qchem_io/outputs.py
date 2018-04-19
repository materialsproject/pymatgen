# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re
import logging
import os
import numpy as np

from monty.io import zopen
from monty.json import MSONable
from pymatgen.core import Molecule

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
        self.data["errors"] = []
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

        # Parse the charge and multiplicity
        self._read_charge_and_multiplicity()

        # Check if calculation finished
        self.data["completion"] = read_pattern(
            self.text, {
                "key": r"Thank you very much for using Q-Chem.\s+Have a nice day."
            }, terminate_on_match=True).get('key')

        # Check if calculation is unrestricted
        self.data["unrestricted"] = read_pattern(
            self.text, {
                "key": r"A(?:n)*\sunrestricted[\s\w\-]+SCF\scalculation\swill\sbe"
            }, terminate_on_match=True).get('key')

        # Check if calculation uses GEN_SCFMAN, multiple potential output formats
        self.data["using_GEN_SCFMAN"] = read_pattern(
            self.text, {
                "key": r"\s+GEN_SCFMAN: A general SCF calculation manager"
            }, terminate_on_match=True).get('key')
        if not self.data["using_GEN_SCFMAN"]:
            self.data["using_GEN_SCFMAN"] = read_pattern(
                self.text, {
                    "key": r"\s+General SCF calculation program by"
                }, terminate_on_match=True).get('key')

        # Check if the SCF failed to converge
        if read_pattern(self.text, {"key": r"SCF failed to converge"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["SCF_failed_to_converge"]

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
            # Then, if no optimized geometry or z-matrix is found, and no errors have been previously
            # idenfied, check to see if the optimization failed to converge or if Lambda wasn't able
            # to be determined. Also, if there is an energy trajectory, read the last geometry in the
            # optimization trajectory for use in the next input file.
            if self.data.get("errors") == [] and self.data.get('optimized_geometry') == [] and self.data.get('optimized_zmat') == []:
                self._check_optimization_errors()
                if self.data.get('energy_trajectory') != None:
                    self._read_last_geometry()

        # Check if the calculation contains a constraint in an $opt section.
        self.data["opt_constraint"] = read_pattern(self.text, {"key": r"\$opt\s+CONSTRAINT"}).get('key')
        if self.data.get('opt_constraint'):
            self.data["dihedral_constraint"] = read_pattern(
                self.text, {
                    "key": r"Constraints and their Current Values\s+Value\s+Constraint\s+Dihedral\:\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
                }).get('key')

        # Check if the calculation is a frequency analysis. If so, parse the relevant output
        self.data["frequency_job"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s+=\s+freq"
            }, terminate_on_match=True).get('key')
        if self.data.get('frequency_job', []):
            self._read_frequency_data()

        # If the calculation did not finish and no errors have been identified yet, check for other errors
        if not self.data.get('completion',[]) and self.data.get("errors") == []:
            self._check_completion_errors()


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

    def _read_charge_and_multiplicity(self):
        temp_charge = read_pattern(self.text, {"key": r"\$molecule\s+([\-\d])+\s+\d"}, terminate_on_match=True).get('key')
        if temp_charge != None:
            self.data["charge"] = int(temp_charge[0][0])
        else:
            temp_charge = read_pattern(self.text, {"key": r"Sum of atomic charges \=\s+([\d\-\.\+]+)"}, terminate_on_match=True).get('key')
            if temp_charge != None:
                self.data["charge"] = int(float(temp_charge[0][0]))
            else:
                self.data["charge"] = None

        temp_multiplicity = read_pattern(self.text, {"key": r"\$molecule\s+[\-\d]+\s+(\d)"}, terminate_on_match=True).get('key')
        if temp_multiplicity != None:
            self.data["multiplicity"] = int(temp_multiplicity[0][0])
        else:
            temp_multiplicity = read_pattern(self.text, {"key": r"Sum of spin\s+charges \=\s+([\d\-\.\+]+)"}, terminate_on_match=True).get('key')
            if temp_multiplicity == None:
                self.data["multiplicity"] = 1
            else:
                self.data["multiplicity"] = int(float(temp_multiplicity[0][0]))+1

    def _read_GEN_SCFMAN(self):
        """
        Parses all GEN_SCFMANs.
        """
        if "SCF_failed_to_converge" in self.data.get("errors"):
            footer_pattern = r"^\s*gen_scfman_exception: SCF failed to converge"
        else:
            footer_pattern = r"^\s*\-+\n"
        header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*(?:RMS Gradient)*\s+\-+(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n"
        table_pattern = r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s+([\d\-\.]+)\s+([\d\-\.]+)e([\d\-\.\+]+)(?:\s+Convergence criterion met)*(?:\s+Preconditoned Steepest Descent)*(?:\s+Roothaan Step)*(?:\s+(?:Normal\s+)*BFGS [Ss]tep)*(?:\s+LineSearch Step)*(?:\s+Line search: overstep)*(?:\s+Descent step)*"

        self.data["GEN_SCFMAN"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)

    def _read_SCF(self):
        """
        Parses all old-style SCFs.
        """
        if "SCF_failed_to_converge" in self.data.get("errors"):
            footer_pattern = r"^\s*\d+\s*[\d\-\.]+\s+[\d\-\.]+E[\d\-\.]+\s+Convergence\s+failure\n"
        else:
            footer_pattern = r"^\s*\-+\n"
        header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+DIIS Error\s+\-+\n"
        table_pattern = r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s*([\d\-\.]+)\s+([\d\-\.]+)E([\d\-\.\+]+)(?:\s*\n\s*cpu\s+[\d\-\.]+\swall\s+[\d\-\.]+)*(?:\nin dftxc\.C, eleTot sum is:[\d\-\.]+, tauTot is\:[\d\-\.]+)*(?:\s+Convergence criterion met)*(?:\s+Done RCA\. Switching to DIIS)*(?:\n\s*Warning: not using a symmetric Q)*(?:\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+(?:\s*\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+)*)*"

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

        self.data["unrestricted_Mulliken"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)

    def _make_geometry_into_molecule(self, geometry):
        """
        Takes a parsed geometry and makes it into a pymatgen Molecule object for storage.
        Makes it easier to use a geometry simply from other modules.
        """
        coords = []
        species = []
        for entry in enumerate(geometry):
            temp_coords = []
            species += [entry[1][0]]
            for jj in range(3):
                temp_coords += [float(entry[1][jj+1])]
            coords += [temp_coords]
        return Molecule(species=species, coords=coords, charge=self.data.get('charge'), spin_multiplicity=self.data.get('multiplicity'))

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
        elif self.data.get('charge') != None:
            self.data["molecule_from_optimized_geometry"] = self._make_geometry_into_molecule(self.data.get('optimized_geometry')[0])

    def _read_last_geometry(self):
        """
        Parses the last geometry from an optimization trajectory for use in a new input file.
        """
        header_pattern = r"\s+Optimization\sCycle:\s+"+str(len(self.data.get("energy_trajectory")))+"\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+(\w+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Point Group\:\s+[\d\w]+\s+Number of degrees of freedom\:\s+\d+"

        self.data["last_geometry"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        if self.data.get('last_geometry') != [] and self.data.get('charge') != None:
            self.data["molecule_from_last_geometry"] = self._make_geometry_into_molecule(self.data.get('last_geometry')[0])

    def _read_frequency_data(self):
        """
        Parses frequencies, enthalpy, entropy, mode vectors, and the geometry of the molecule.
        """
        temp_dict = read_pattern(self.text, {
                    "frequencies": r"\s*Frequency:\s+([\d\-\.]+)(?:\s+([\d\-\.]+)(?:\s+([\d\-\.]+))*)*",
                    "enthalpy": r"\s*Total Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                    "entropy": r"\s*Total Entropy:\s+([\d\-\.]+)\s+cal/mol\.K"})
        for key in temp_dict:
            self.data[key] = temp_dict.get(key)

        header_pattern = r"\s*Raman Active:\s+[YESNO]+\s+(?:[YESNO]+\s+)*X\s+Y\s+Z\s+(?:X\s+Y\s+Z\s+)*"
        table_pattern = r"\s*([a-zA-Z][a-zA-Z\s])\s*([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*(?:([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*(?:([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+))*)*"
        footer_pattern = r"TransDip\s+[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+\s*(?:[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+\s*)*"
        self.data["frequency_mode_vectors"] = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        if self.data.get('frequencies') != None:
            if float(self.data.get('frequencies')[0][0]) < 0.0:
                all_vecs = self.data.get('frequency_mode_vectors')
                temp_vecs = np.zeros(shape=(len(all_vecs[0]),3),dtype=float)
                for ii, entry in enumerate(all_vecs[0]):
                    for jj in range(1,4):
                        temp_vecs[ii,jj-1] = float(entry[jj])
                self.data["negative_freq_vecs"] = temp_vecs

        header_pattern = r"Standard Nuclear Orientation \(Angstroms\)\s+I\s+Atom\s+X\s+Y\s+Z\s+-+"
        table_pattern = r"\s*\d+\s+([a-zA-Z]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*"
        footer_pattern = r"\s*-+"
        tmp_freq_geom = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)[0]
        freq_species = []
        freq_geometry = np.zeros(shape=(len(tmp_freq_geom),3),dtype=float)
        for ii, entry in enumerate(tmp_freq_geom):
            freq_species += [entry[0]]
            for jj in range(3):
                freq_geometry[ii,jj] = float(entry[jj+1])
        self.data["freq_species"] = freq_species
        self.data["freq_geometry"] = freq_geometry

    def _check_optimization_errors(self):
        """
        Parses three potential optimization errors: failing to converge within the allowed number
        of optimization cycles, failure to determine the lamda needed to continue, and inconsistent
        size of MO files due to a linear dependence in the AO basis.
        """
        if read_pattern(self.text, {"key": r"MAXIMUM OPTIMIZATION CYCLES REACHED"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["out_of_opt_cycles"]
        elif read_pattern(self.text, {"key": r"UNABLE TO DETERMINE Lamda IN FormD"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["unable_to_determine_lamda"]
        elif read_pattern(self.text, {"key": r"Inconsistent size for SCF MO coefficient file"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["linear_dependent_basis"]

    def _check_completion_errors(self):
        """
        Parses four potential errors that can cause jobs to crash: inability to transform
        coordinates due to a bad symmetric specification, an input file that fails to pass
        inspection, and errors reading and writing files.
        """
        if read_pattern(self.text, {"key": r"Coordinates do not transform within specified threshold"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["failed_to_transform_coords"]
        elif read_pattern(self.text, {"key": r"The Q\-Chem input file has failed to pass inspection"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["input_file_error"]
        elif read_pattern(self.text, {"key": r"Error opening input stream"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["failed_to_read_input"]
        elif read_pattern(self.text, {"key": r"FileMan error: End of file reached prematurely"}, terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["IO_error"]
        else:
            self.data["errors"] += ["unknown_error"]




