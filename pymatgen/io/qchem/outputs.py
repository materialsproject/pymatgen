# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re
import logging
import os
import numpy as np
import math
import copy

from monty.io import zopen
from monty.json import jsanitize
from monty.json import MSONable
from pymatgen.core import Molecule

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import metal_edge_extender
import networkx as nx

try:
    import openbabel as ob

    have_babel = True
except ImportError:
    ob = None
    have_babel = False

from .utils import read_table_pattern, read_pattern, process_parsed_coords

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class QCOutput(MSONable):
    """
    Class to parse QChem output files.
    """

    def __init__(self, filename):
        """
        Args:
            filename (str): Filename to parse
        """
        self.filename = filename
        self.data = {}
        self.data["errors"] = []
        self.data["warnings"] = {}
        self.text = ""
        with zopen(filename, 'rt') as f:
            self.text = f.read()

        # Check if output file contains multiple output files. If so, print an error message and exit
        self.data["multiple_outputs"] = read_pattern(
            self.text, {
                "key": r"Job\s+\d+\s+of\s+(\d+)\s+"
            },
            terminate_on_match=True).get('key')
        if not (self.data.get('multiple_outputs') is None
                or self.data.get('multiple_outputs') == [['1']]):
            print(
                "ERROR: multiple calculation outputs found in file " +
                filename +
                ". Please instead call QCOutput.mulitple_outputs_from_file(QCOutput,'"
                + filename + "')")
            print("Exiting...")
            exit()

        # Parse the molecular details: charge, multiplicity,
        # species, and initial geometry.
        self._read_charge_and_multiplicity()
        if read_pattern(
                self.text, {
                    "key": r"Nuclear Repulsion Energy"
                },
                terminate_on_match=True).get('key') == [[]]:
            self._read_species_and_inital_geometry()

        # Check if calculation finished
        self.data["completion"] = read_pattern(
            self.text, {
                "key":
                    r"Thank you very much for using Q-Chem.\s+Have a nice day."
            },
            terminate_on_match=True).get('key')

        # If the calculation finished, parse the job time.
        if self.data.get('completion', []):
            temp_timings = read_pattern(
                self.text, {
                    "key":
                        r"Total job time\:\s*([\d\-\.]+)s\(wall\)\,\s*([\d\-\.]+)s\(cpu\)"
                }).get('key')
            if temp_timings is not None:
                self.data["walltime"] = float(temp_timings[0][0])
                self.data["cputime"] = float(temp_timings[0][1])
            else:
                self.data["walltime"] = None
                self.data["cputime"] = None

        # Check if calculation is unrestricted
        self.data["unrestricted"] = read_pattern(
            self.text, {
                "key":
                    r"A(?:n)*\sunrestricted[\s\w\-]+SCF\scalculation\swill\sbe"
            },
            terminate_on_match=True).get('key')

        # Check if calculation uses GEN_SCFMAN, multiple potential output formats
        self.data["using_GEN_SCFMAN"] = read_pattern(
            self.text, {
                "key": r"\s+GEN_SCFMAN: A general SCF calculation manager"
            },
            terminate_on_match=True).get('key')
        if not self.data["using_GEN_SCFMAN"]:
            self.data["using_GEN_SCFMAN"] = read_pattern(
                self.text, {
                    "key": r"\s+General SCF calculation program by"
                },
                terminate_on_match=True).get('key')

        # Check if the SCF failed to converge
        if read_pattern(
                self.text, {
                    "key": r"SCF failed to converge"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["SCF_failed_to_converge"]

        # Parse the SCF
        self._read_SCF()

        # Parse the Mulliken/ESP/RESP charges
        self._read_charges()

        # Check for various warnings
        self._detect_general_warnings()

        # Check to see if PCM or SMD are present
        self.data["solvent_method"] = None
        self.data["solvent_data"] = None
        if read_pattern(
                self.text, {
                    "key": r"solvent_method\s*=?\s*pcm"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["solvent_method"] = "PCM"
        if read_pattern(
                self.text, {
                    "key": r"solvent_method\s*=?\s*smd"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["solvent_method"] = "SMD"

        # Parse information specific to a solvent model
        if self.data["solvent_method"] == "PCM":
            self.data["solvent_data"] = {}
            temp_dielectric = read_pattern(
                self.text, {
                    "key": r"dielectric\s*([\d\-\.]+)"
                },
                terminate_on_match=True).get('key')
            self.data["solvent_data"]["PCM_dielectric"] = float(temp_dielectric[0][0])
            self._read_pcm_information()
        elif self.data["solvent_method"] == "SMD":
            if read_pattern(
                    self.text, {
                        "key": r"Unrecognized solvent"
                    },
                    terminate_on_match=True).get('key') == [[]]:
                if not self.data.get('completion', []):
                    self.data["errors"] += ["unrecognized_solvent"]
                else:
                    self.data["warnings"]["unrecognized_solvent"] = True
            self.data["solvent_data"] = {}
            temp_solvent = read_pattern(
                self.text, {
                    "key": r"\s[Ss]olvent:? ([a-zA-Z]+)"
                }).get('key')
            for val in temp_solvent:
                if val[0] != temp_solvent[0][0]:
                    if val[0] != "for":
                        self.data["warnings"]["SMD_two_solvents"] = str(temp_solvent[0][0]) + " and " + str(val[0])
                    else:
                        if "unrecognized_solvent" not in self.data["errors"] and "unrecognized_solvent" not in \
                                self.data["warnings"]:
                            self.data["warnings"]["questionable_SMD_parsing"] = True
            self.data["solvent_data"]["SMD_solvent"] = temp_solvent[0][0]
            self._read_smd_information()

        # Parse the final energy
        temp_final_energy = read_pattern(
            self.text, {
                "key": r"Final\senergy\sis\s+([\d\-\.]+)"
            }).get('key')
        if temp_final_energy is None:
            self.data["final_energy"] = None
        else:
            self.data["final_energy"] = float(temp_final_energy[0][0])

        # Check if calculation is using dft_d and parse relevant info if so
        self.data["using_dft_d3"] = read_pattern(
            self.text, {
                "key":
                    r"dft_d\s*= d3"
            },
            terminate_on_match=True).get('key')
        if self.data.get('using_dft_d3', []):
            temp_d3 = read_pattern(self.text, {
                "key": r"\-D3 energy without 3body term =\s*([\d\.\-]+) hartrees"
            }).get('key')
            real_d3 = np.zeros(len(temp_d3))
            if temp_d3 is None:
                self.data["dft_d3"] = None
            elif len(temp_d3) == 1:
                self.data["dft_d3"] = float(temp_d3[0][0])
            else:
                for ii, entry in enumerate(temp_d3):
                    real_d3[ii] = float(entry[0])
                self.data["dft_d3"] = real_d3

        # Parse the S2 values in the case of an unrestricted calculation
        if self.data.get('unrestricted', []):
            correct_s2 = 0.5 * (self.data["multiplicity"] - 1) * (0.5 * (self.data["multiplicity"] - 1) + 1)
            temp_S2 = read_pattern(self.text, {
                "key": r"<S\^2>\s=\s+([\d\-\.]+)"
            }).get('key')
            if temp_S2 is None:
                self.data["S2"] = None
            elif len(temp_S2) == 1:
                self.data["S2"] = float(temp_S2[0][0])
                if abs(correct_s2 - self.data["S2"]) > 0.01:
                    self.data["warnings"]["spin_contamination"] = abs(correct_s2 - self.data["S2"])
            else:
                real_S2 = np.zeros(len(temp_S2))
                have_spin_contamination = False
                for ii, entry in enumerate(temp_S2):
                    real_S2[ii] = float(entry[0])
                    if abs(correct_s2 - real_S2[ii]) > 0.01:
                        have_spin_contamination = True
                self.data["S2"] = real_S2
                if have_spin_contamination:
                    spin_contamination = np.zeros(len(self.data["S2"]))
                    for ii, entry in enumerate(self.data["S2"]):
                        spin_contamination[ii] = abs(correct_s2 - entry)
                    self.data["warnings"]["spin_contamination"] = spin_contamination

        # Check if the calculation is a geometry optimization. If so, parse the relevant output
        self.data["optimization"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*opt"
            }).get('key')
        if self.data.get('optimization', []):
            self._read_optimization_data()

        # Check if the calculation contains a constraint in an $opt section.
        self.data["opt_constraint"] = read_pattern(self.text, {
            "key": r"\$opt\s+CONSTRAINT"
        }).get('key')
        if self.data.get('opt_constraint'):
            temp_constraint = read_pattern(
                self.text, {
                    "key": r"Constraints and their Current Values\s+Value\s+Constraint\s+(\w+)\:\s+([\d\-\.]+)\s+"
                           r"([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
                }).get('key')
            if temp_constraint is not None:
                self.data["opt_constraint"] = temp_constraint[0]
                if float(self.data.get('opt_constraint')[5]) != float(
                        self.data.get('opt_constraint')[6]):
                    if abs(float(self.data.get('opt_constraint')[5])) != abs(
                            float(self.data.get('opt_constraint')[6])):
                        raise ValueError(
                            "ERROR: Opt section value and constraint should be the same!"
                        )
                    elif abs(float(
                            self.data.get('opt_constraint')[5])) not in [
                        0.0, 180.0
                    ]:
                        raise ValueError(
                            "ERROR: Opt section value and constraint can only differ by a sign at 0.0 and 180.0!"
                        )

        # Check if the calculation is a frequency analysis. If so, parse the relevant output
        self.data["frequency_job"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*freq"
            },
            terminate_on_match=True).get('key')
        if self.data.get('frequency_job', []):
            self._read_frequency_data()

        self.data["single_point_job"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*sp"
            },
            terminate_on_match=True).get("key")
        if self.data.get("single_point_job", []):
            self._read_single_point_data()

        self.data["force_job"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*force"
            },
            terminate_on_match=True).get("key")
        if self.data.get("force", []):
            self._read_force_data()

        # If the calculation did not finish and no errors have been identified yet, check for other errors
        if not self.data.get('completion',
                             []) and self.data.get("errors") == []:
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
            text = re.split(r'\s*(?:Running\s+)*Job\s+\d+\s+of\s+\d+\s+',
                            f.read())
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
        """
        Parses charge and multiplicity.
        """
        temp_charge = read_pattern(
            self.text, {
                "key": r"\$molecule\s+([\-\d]+)\s+\d"
            },
            terminate_on_match=True).get('key')
        if temp_charge is not None:
            self.data["charge"] = int(temp_charge[0][0])
        else:
            temp_charge = read_pattern(
                self.text, {
                    "key": r"Sum of atomic charges \=\s+([\d\-\.\+]+)"
                },
                terminate_on_match=True).get('key')
            if temp_charge is None:
                self.data["charge"] = None
            else:
                self.data["charge"] = int(float(temp_charge[0][0]))

        temp_multiplicity = read_pattern(
            self.text, {
                "key": r"\$molecule\s+[\-\d]+\s+(\d)"
            },
            terminate_on_match=True).get('key')
        if temp_multiplicity is not None:
            self.data["multiplicity"] = int(temp_multiplicity[0][0])
        else:
            temp_multiplicity = read_pattern(
                self.text, {
                    "key": r"Sum of spin\s+charges \=\s+([\d\-\.\+]+)"
                },
                terminate_on_match=True).get('key')
            if temp_multiplicity is None:
                self.data["multiplicity"] = 1
            else:
                self.data["multiplicity"] = int(
                    float(temp_multiplicity[0][0])) + 1

    def _read_species_and_inital_geometry(self):
        """
        Parses species and initial geometry.
        """
        header_pattern = r"Standard Nuclear Orientation \(Angstroms\)\s+I\s+Atom\s+X\s+Y\s+Z\s+-+"
        table_pattern = r"\s*\d+\s+([a-zA-Z]+)\s*([\d\-\.\*]+)\s*([\d\-\.\*]+)\s*([\d\-\.\*]+)\s*"
        footer_pattern = r"\s*-+"
        temp_geom = read_table_pattern(self.text, header_pattern,
                                       table_pattern, footer_pattern)
        if temp_geom is None or len(temp_geom) == 0:
            self.data["species"] = None
            self.data["initial_geometry"] = None
            self.data["initial_molecule"] = None
            self.data["point_group"] = None
        else:
            temp_point_group = read_pattern(
                self.text, {
                    "key":
                        r"Molecular Point Group\s+([A-Za-z\d\*]+)"
                },
                terminate_on_match=True).get('key')
            if temp_point_group is not None:
                self.data["point_group"] = temp_point_group[0][0]
            else:
                self.data["point_group"] = None
            temp_geom = temp_geom[0]
            species = []
            geometry = np.zeros(shape=(len(temp_geom), 3), dtype=float)
            for ii, entry in enumerate(temp_geom):
                species += [entry[0]]
                for jj in range(3):
                    if "*" in entry[jj + 1]:
                        geometry[ii, jj] = 10000000000.0
                    else:
                        geometry[ii, jj] = float(entry[jj + 1])
            self.data["species"] = species
            self.data["initial_geometry"] = geometry
            if self.data["charge"] is not None and self.data["multiplicity"] is not None:
                self.data["initial_molecule"] = Molecule(
                    species=species,
                    coords=geometry,
                    charge=self.data.get('charge'),
                    spin_multiplicity=self.data.get('multiplicity'))
            else:
                self.data["initial_molecule"] = None

    def _read_SCF(self):
        """
        Parses both old and new SCFs.
        """
        if self.data.get('using_GEN_SCFMAN', []):
            if "SCF_failed_to_converge" in self.data.get("errors"):
                footer_pattern = r"^\s*gen_scfman_exception: SCF failed to converge"
            else:
                footer_pattern = r"^\s*\-+\n\s+SCF time"
            header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*(?:RMS Gradient)*\s+\-+" \
                             r"(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+" \
                             r"(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, " \
                             r"Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n"
            table_pattern = r"(?:\n[a-zA-Z_\s/]+(?:\.C)*::(?:WARNING energy changes are now smaller than effective " \
                            r"accuracy\.)*(?:\s+calculation will continue, but THRESH should be increased)*" \
                            r"(?:\s+or SCF_CONVERGENCE decreased\. )*(?:\s+effective_thresh = [\d\-\.]+e[\d\-]+)*)*" \
                            r"(?:\s*Nonlocal correlation = [\d\-\.]+e[\d\-]+)*" \
                            r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+" \
                            r"Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s+" \
                            r"([\d\-\.]+)\s+([\d\-\.]+)e([\d\-\.\+]+)(?:\s+Convergence criterion met)*" \
                            r"(?:\s+Preconditoned Steepest Descent)*(?:\s+Roothaan Step)*(?:\s+" \
                            r"(?:Normal\s+)*BFGS [Ss]tep)*(?:\s+LineSearch Step)*(?:\s+Line search: overstep)*" \
                            r"(?:\s+Dog-leg BFGS step)*(?:\s+Line search: understep)*" \
                            r"(?:\s+Descent step)*(?:\s+Done DIIS. Switching to GDM)*" \
                            r"(?:\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*" \
                            r"(?:RMS Gradient)*\s+\-+(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+" \
                            r"(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, " \
                            r"Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n)*"
        else:
            if "SCF_failed_to_converge" in self.data.get("errors"):
                footer_pattern = r"^\s*\d+\s*[\d\-\.]+\s+[\d\-\.]+E[\d\-\.]+\s+Convergence\s+failure\n"
            else:
                footer_pattern = r"^\s*\-+\n"
            header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+DIIS Error\s+\-+\n"
            table_pattern = r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+" \
                            r"Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s*" \
                            r"([\d\-\.]+)\s+([\d\-\.]+)E([\d\-\.\+]+)(?:\s*\n\s*cpu\s+[\d\-\.]+\swall\s+[\d\-\.]+)*" \
                            r"(?:\nin dftxc\.C, eleTot sum is:[\d\-\.]+, tauTot is\:[\d\-\.]+)*" \
                            r"(?:\s+Convergence criterion met)*(?:\s+Done RCA\. Switching to DIIS)*" \
                            r"(?:\n\s*Warning: not using a symmetric Q)*" \
                            r"(?:\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+" \
                            r"(?:\s*\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+)*)*"

        temp_scf = read_table_pattern(self.text, header_pattern, table_pattern,
                                      footer_pattern)
        real_scf = []
        for one_scf in temp_scf:
            temp = np.zeros(shape=(len(one_scf), 2))
            for ii, entry in enumerate(one_scf):
                temp[ii, 0] = float(entry[0])
                temp[ii, 1] = float(entry[1]) * 10 ** float(entry[2])
            real_scf += [temp]

        self.data["SCF"] = real_scf

        temp_thresh_warning = read_pattern(self.text, {
            "key": r"\n[a-zA-Z_\s/]+\.C::WARNING energy changes are now smaller than effective accuracy"
                   r"\.\n[a-zA-Z_\s/]+\.C::\s+calculation will continue, but THRESH should be increased\n"
                   r"[a-zA-Z_\s/]+\.C::\s+or SCF_CONVERGENCE decreased\. \n"
                   r"[a-zA-Z_\s/]+\.C::\s+effective_thresh = ([\d\-\.]+e[\d\-]+)"
        }).get('key')
        if temp_thresh_warning is not None:
            if len(temp_thresh_warning) == 1:
                self.data["warnings"]["thresh"] = float(temp_thresh_warning[0][0])
            else:
                thresh_warning = np.zeros(len(temp_thresh_warning))
                for ii, entry in enumerate(temp_thresh_warning):
                    thresh_warning[ii] = float(entry[0])
                self.data["warnings"]["thresh"] = thresh_warning

        temp_SCF_energy = read_pattern(self.text, {
            "key": r"SCF   energy in the final basis set =\s*([\d\-\.]+)"
        }).get('key')
        if temp_SCF_energy is not None:
            if len(temp_SCF_energy) == 1:
                self.data["SCF_energy_in_the_final_basis_set"] = float(temp_SCF_energy[0][0])
            else:
                SCF_energy = np.zeros(len(temp_SCF_energy))
                for ii, val in enumerate(temp_SCF_energy):
                    SCF_energy[ii] = float(val[0])
                self.data["SCF_energy_in_the_final_basis_set"] = SCF_energy

        temp_Total_energy = read_pattern(self.text, {
            "key": r"Total energy in the final basis set =\s*([\d\-\.]+)"
        }).get('key')
        if temp_Total_energy is not None:
            if len(temp_Total_energy) == 1:
                self.data["Total_energy_in_the_final_basis_set"] = float(temp_Total_energy[0][0])
            else:
                Total_energy = np.zeros(len(temp_Total_energy))
                for ii, val in enumerate(temp_Total_energy):
                    Total_energy[ii] = float(val[0])
                self.data["Total_energy_in_the_final_basis_set"] = Total_energy

    def _read_charges(self):
        """
        Parses Mulliken/ESP/RESP charges. Also parses spins given an unrestricted SCF.
        """
        if self.data.get('unrestricted', []):
            header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+" \
                             r"Spin\s\(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"
        else:
            header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

        temp_mulliken = read_table_pattern(self.text, header_pattern,
                                           table_pattern, footer_pattern)
        real_mulliken = []
        for one_mulliken in temp_mulliken:
            if self.data.get('unrestricted', []):
                temp = np.zeros(shape=(len(one_mulliken), 2))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii, 0] = float(entry[0])
                    temp[ii, 1] = float(entry[1])
            else:
                temp = np.zeros(len(one_mulliken))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii] = float(entry[0])
            real_mulliken += [temp]

        self.data["Mulliken"] = real_mulliken

        # Check for ESP/RESP charges
        esp_or_resp = read_pattern(
            self.text, {
                "key": r"Merz-Kollman (R?ESP) Net Atomic Charges"
            }).get('key')
        if esp_or_resp is not None:
            header_pattern = r"Merz-Kollman (R?ESP) Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

            temp_esp_or_resp = read_table_pattern(self.text, header_pattern,
                                                  table_pattern, footer_pattern)
            real_esp_or_resp = []
            for one_entry in temp_esp_or_resp:
                temp = np.zeros(len(one_entry))
                for ii, entry in enumerate(one_entry):
                    temp[ii] = float(entry[0])
                real_esp_or_resp += [temp]
            self.data[esp_or_resp[0][0]] = real_esp_or_resp

    def _detect_general_warnings(self):
        # Check for inaccurate integrated density
        temp_inac_integ = read_pattern(
            self.text, {
                "key": r"Inaccurate integrated density:\n\s+Number of electrons\s+=\s+([\d\-\.]+)\n\s+"
                       r"Numerical integral\s+=\s+([\d\-\.]+)\n\s+Relative error\s+=\s+([\d\-\.]+)\s+\%\n"
            }).get('key')
        if temp_inac_integ is not None:
            inaccurate_integrated_density = np.zeros(shape=(len(temp_inac_integ), 3))
            for ii, entry in enumerate(temp_inac_integ):
                for jj, val in enumerate(entry):
                    inaccurate_integrated_density[ii][jj] = float(val)
            self.data["warnings"]["inaccurate_integrated_density"] = inaccurate_integrated_density

        # Check for an MKL error
        if read_pattern(
                self.text, {
                    "key": r"Intel MKL ERROR"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["mkl"] = True

        # Check if the job is being hindered by a lack of analytical derivatives
        if read_pattern(
                self.text, {
                    "key": r"Starting finite difference calculation for IDERIV"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["missing_analytical_derivates"] = True

        # Check if the job is complaining about MO files of inconsistent size
        if read_pattern(
                self.text, {
                    "key": r"Inconsistent size for SCF MO coefficient file"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["inconsistent_size"] = True

        # Check for AO linear depend
        if read_pattern(
                self.text, {
                    "key": r"Linear dependence detected in AO basis"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["linear_dependence"] = True

        # Check for Hessian without desired local structure
        if read_pattern(
                self.text, {
                    "key": r"\*\*WARNING\*\* Hessian does not have the Desired Local Structure"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["hessian_local_structure"] = True

        # Check if GetCART cycle iterations ever exceeded
        if read_pattern(
                self.text, {
                    "key": r"\*\*\*ERROR\*\*\* Exceeded allowed number of iterative cycles in GetCART"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["GetCART_cycles"] = True

        # Check for problems with internal coordinates
        if read_pattern(
                self.text, {
                    "key": r"\*\*WARNING\*\* Problems with Internal Coordinates"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["internal_coordinates"] = True

        # Check for problem with eigenvalue magnitude
        if read_pattern(
                self.text, {
                    "key": r"\*\*WARNING\*\* Magnitude of eigenvalue"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["eigenvalue_magnitude"] = True

        # Check for problem with hereditary postivive definiteness
        if read_pattern(
                self.text, {
                    "key": r"\*\*WARNING\*\* Hereditary positive definiteness endangered"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["positive_definiteness_endangered"] = True

        # Check if there were problems with a colinear bend
        if read_pattern(
                self.text, {
                    "key": r"\*\*\*ERROR\*\*\* Angle[\s\d]+is near\-linear\s+"
                           r"But No atom available to define colinear bend"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["colinear_bend"] = True

        # Check if there were problems diagonalizing B*B(t)
        if read_pattern(
                self.text, {
                    "key": r"\*\*\*ERROR\*\*\* Unable to Diagonalize B\*B\(t\) in <MakeNIC>"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["warnings"]["diagonalizing_BBt"] = True

        # Check for bad Roothaan step
        for scf in self.data["SCF"]:
            if abs(scf[0][0] - scf[1][0]) > 10.0:
                self.data["warnings"]["bad_roothaan"] = True
                if abs(scf[0][0] - scf[1][0]) > 100.0:
                    self.data["warnings"]["very_bad_roothaan"] = True

    def _read_geometries(self):
        """
        Parses all geometries from an optimization trajectory.
        """
        geoms = []
        header_pattern = r"\s+Optimization\sCycle:\s+\d+\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+\w+\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Point Group\:\s+[\d\w\*]+\s+Number of degrees of freedom\:\s+\d+"

        parsed_geometries = read_table_pattern(
            self.text, header_pattern, table_pattern, footer_pattern)
        for ii, parsed_geometry in enumerate(parsed_geometries):
            if parsed_geometry == [] or None:
                geoms.append(None)
            else:
                geoms.append(process_parsed_coords(parsed_geometry))
        self.data["geometries"] = geoms
        self.data["last_geometry"] = geoms[-1]
        if self.data.get('charge') is not None:
            self.data["molecule_from_last_geometry"] = Molecule(
                species=self.data.get('species'),
                coords=self.data.get('last_geometry'),
                charge=self.data.get('charge'),
                spin_multiplicity=self.data.get('multiplicity'))

        # Parses optimized XYZ coordinates. If not present, parses optimized Z-matrix.
        header_pattern = r"\*+\s+OPTIMIZATION\s+CONVERGED\s+\*+\s+\*+\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+\w+\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Z-matrix Print:"

        parsed_optimized_geometry = read_table_pattern(
            self.text, header_pattern, table_pattern, footer_pattern)
        if parsed_optimized_geometry == [] or None:
            self.data["optimized_geometry"] = None
            header_pattern = r"^\s+\*+\s+OPTIMIZATION CONVERGED\s+\*+\s+\*+\s+Z-matrix\s+" \
                             r"Print:\s+\$molecule\s+[\d\-]+\s+[\d\-]+\n"
            table_pattern = r"\s*(\w+)(?:\s+(\d+)\s+([\d\-\.]+)(?:\s+(\d+)\s+([\d\-\.]+)" \
                            r"(?:\s+(\d+)\s+([\d\-\.]+))*)*)*(?:\s+0)*"
            footer_pattern = r"^\$end\n"

            self.data["optimized_zmat"] = read_table_pattern(
                self.text, header_pattern, table_pattern, footer_pattern)
        else:
            self.data["optimized_geometry"] = process_parsed_coords(
                parsed_optimized_geometry[0])
            if self.data.get('charge') is not None:
                self.data["molecule_from_optimized_geometry"] = Molecule(
                    species=self.data.get('species'),
                    coords=self.data.get('optimized_geometry'),
                    charge=self.data.get('charge'),
                    spin_multiplicity=self.data.get('multiplicity'))

    def _get_grad_format_length(self, header):
        """
        Determines the maximum number of gradient entries printed on a line,
        which changes for different versions of Q-Chem
        """
        found_end = False
        index = 1
        pattern = header
        while not found_end:
            if read_pattern(
                    self.text, {
                        "key": pattern
                    },
                    terminate_on_match=True).get('key') != [[]]:
                found_end = True
            else:
                pattern = pattern + r"\s+" + str(index)
                index += 1
        return index - 2

    def _read_gradients(self):
        """
        Parses all gradients obtained during an optimization trajectory
        """

        grad_header_pattern = r"Gradient of (?:SCF)?(?:MP2)? Energy(?: \(in au\.\))?"
        footer_pattern = r"(?:Max gradient component|Gradient time)"

        grad_format_length = self._get_grad_format_length(grad_header_pattern)
        grad_table_pattern = r"(?:\s+\d+(?:\s+\d+)?(?:\s+\d+)?(?:\s+\d+)?(?:\s+\d+)?(?:\s+\d+)?)?\n\s\s\s\s[1-3]\s*" \
                             r"(\-?[\d\.]{9,12})"
        if grad_format_length > 1:
            for ii in range(1, grad_format_length):
                grad_table_pattern = grad_table_pattern + r"(?:\s*(\-?[\d\.]{9,12}))?"

        parsed_gradients = read_table_pattern(
            self.text, grad_header_pattern, grad_table_pattern, footer_pattern)
        sorted_gradients = np.zeros(shape=(len(parsed_gradients), len(self.data["initial_molecule"]), 3))
        for ii, grad in enumerate(parsed_gradients):
            for jj in range(int(len(grad) / 3)):
                for kk in range(grad_format_length):
                    if grad[jj * 3][kk] != 'None':
                        sorted_gradients[ii][jj * grad_format_length + kk][0] = grad[jj * 3][kk]
                        sorted_gradients[ii][jj * grad_format_length + kk][1] = grad[jj * 3 + 1][kk]
                        sorted_gradients[ii][jj * grad_format_length + kk][2] = grad[jj * 3 + 2][kk]

        self.data["gradients"] = sorted_gradients

        if self.data["solvent_method"] is not None:
            header_pattern = r"total gradient after adding PCM contribution --\s+-+\s+Atom\s+X\s+Y\s+Z\s+-+"
            table_pattern = r"\s+\d+\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s"
            footer_pattern = r"-+"

            parsed_gradients = read_table_pattern(
                self.text, header_pattern, table_pattern, footer_pattern)

            pcm_gradients = np.zeros(shape=(len(parsed_gradients), len(self.data["initial_molecule"]), 3))
            for ii, grad in enumerate(parsed_gradients):
                for jj, entry in enumerate(grad):
                    for kk, val in enumerate(entry):
                        pcm_gradients[ii][jj][kk] = float(val)

            self.data["pcm_gradients"] = pcm_gradients
        else:
            self.data["pcm_gradients"] = None

        if read_pattern(self.text, {
            "key": r"Gradient of CDS energy"
        },
                        terminate_on_match=True).get('key') == [[]]:
            header_pattern = r"Gradient of CDS energy"

            parsed_gradients = read_table_pattern(
                self.text, header_pattern, grad_table_pattern, grad_header_pattern)

            sorted_gradients = np.zeros(shape=(len(parsed_gradients), len(self.data["initial_molecule"]), 3))
            for ii, grad in enumerate(parsed_gradients):
                for jj in range(int(len(grad) / 3)):
                    for kk in range(grad_format_length):
                        if grad[jj * 3][kk] != 'None':
                            sorted_gradients[ii][jj * grad_format_length + kk][0] = grad[jj * 3][kk]
                            sorted_gradients[ii][jj * grad_format_length + kk][1] = grad[jj * 3 + 1][kk]
                            sorted_gradients[ii][jj * grad_format_length + kk][2] = grad[jj * 3 + 2][kk]

            self.data["CDS_gradients"] = sorted_gradients
        else:
            self.data["CDS_gradients"] = None

    def _read_optimization_data(self):
        temp_energy_trajectory = read_pattern(
            self.text, {
                "key": r"\sEnergy\sis\s+([\d\-\.]+)"
            }).get('key')
        if temp_energy_trajectory is None:
            self.data["energy_trajectory"] = []
        else:
            real_energy_trajectory = np.zeros(len(temp_energy_trajectory))
            for ii, entry in enumerate(temp_energy_trajectory):
                real_energy_trajectory[ii] = float(entry[0])
            self.data["energy_trajectory"] = real_energy_trajectory
            self._read_geometries()
            if have_babel:
                self.data["structure_change"] = check_for_structure_changes(
                    self.data["initial_molecule"],
                    self.data["molecule_from_last_geometry"])
            self._read_gradients()
            # Then, if no optimized geometry or z-matrix is found, and no errors have been previously
            # idenfied, check to see if the optimization failed to converge or if Lambda wasn't able
            # to be determined.
            if len(self.data.get("errors")) == 0 and self.data.get('optimized_geometry') is None \
                    and len(self.data.get('optimized_zmat')) == 0:
                if read_pattern(
                        self.text, {
                            "key": r"MAXIMUM OPTIMIZATION CYCLES REACHED"
                        },
                        terminate_on_match=True).get('key') == [[]]:
                    self.data["errors"] += ["out_of_opt_cycles"]
                elif read_pattern(
                        self.text, {
                            "key": r"UNABLE TO DETERMINE Lamda IN FormD"
                        },
                        terminate_on_match=True).get('key') == [[]]:
                    self.data["errors"] += ["unable_to_determine_lamda"]

    def _read_frequency_data(self):
        """
        Parses frequencies, enthalpy, entropy, and mode vectors.
        """
        temp_dict = read_pattern(
            self.text, {
                "frequencies":
                    r"\s*Frequency:\s+(\-?[\d\.\*]+)(?:\s+(\-?[\d\.\*]+)(?:\s+(\-?[\d\.\*]+))*)*",
                "trans_dip":
                    r"TransDip\s+(\-?[\d\.]{5,7}|\*{5,7})\s*(\-?[\d\.]{5,7}|\*{5,7})\s*(\-?[\d\.]{5,7}|\*{5,7})\s*"
                    r"(?:(\-?[\d\.]{5,7}|\*{5,7})\s*(\-?[\d\.]{5,7}|\*{5,7})\s*(\-?[\d\.]{5,7}|\*{5,7})\s*"
                    r"(?:(\-?[\d\.]{5,7}|\*{5,7})\s*(\-?[\d\.]{5,7}|\*{5,7})\s*(\-?[\d\.]{5,7}|\*{5,7}))*)*",
                "IR_intens":
                    r"\s*IR Intens:\s*(\-?[\d\.\*]+)(?:\s+(\-?[\d\.\*]+)(?:\s+(\-?[\d\.\*]+))*)*",
                "IR_active":
                    r"\s*IR Active:\s+([YESNO]+)(?:\s+([YESNO]+)(?:\s+([YESNO]+))*)*",
                "ZPE":
                    r"\s*Zero point vibrational energy:\s+([\d\-\.]+)\s+kcal/mol",
                "trans_enthalpy":
                    r"\s*Translational Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                "rot_enthalpy":
                    r"\s*Rotational Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                "vib_enthalpy":
                    r"\s*Vibrational Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                "gas_constant":
                    r"\s*gas constant \(RT\):\s+([\d\-\.]+)\s+kcal/mol",
                "trans_entropy":
                    r"\s*Translational Entropy:\s+([\d\-\.]+)\s+cal/mol\.K",
                "rot_entropy":
                    r"\s*Rotational Entropy:\s+([\d\-\.]+)\s+cal/mol\.K",
                "vib_entropy":
                    r"\s*Vibrational Entropy:\s+([\d\-\.]+)\s+cal/mol\.K",
                "total_enthalpy":
                    r"\s*Total Enthalpy:\s+([\d\-\.]+)\s+kcal/mol",
                "total_entropy":
                    r"\s*Total Entropy:\s+([\d\-\.]+)\s+cal/mol\.K"
            })

        keys = ["ZPE", "trans_enthalpy", "rot_enthalpy", "vib_enthalpy", "gas_constant", "trans_entropy", "rot_entropy",
                "vib_entropy", "total_enthalpy", "total_entropy"]

        for key in keys:
            if temp_dict.get(key) is None:
                self.data[key] = None
            else:
                self.data[key] = float(temp_dict.get(key)[0][0])

        if temp_dict.get('frequencies') is None:
            self.data['frequencies'] = None
            self.data['IR_intens'] = None
            self.data['IR_active'] = None
            self.data['trans_dip'] = None
        else:
            temp_freqs = [
                value for entry in temp_dict.get('frequencies')
                for value in entry
            ]
            temp_intens = [
                value for entry in temp_dict.get('IR_intens')
                for value in entry
            ]
            active = [
                value for entry in temp_dict.get('IR_active')
                for value in entry
            ]
            temp_trans_dip = [
                value for entry in temp_dict.get('trans_dip')
                for value in entry
            ]
            self.data['IR_active'] = active

            trans_dip = np.zeros(shape=(int((len(temp_trans_dip) - temp_trans_dip.count('None')) / 3), 3))
            for ii, entry in enumerate(temp_trans_dip):
                if entry != 'None':
                    if "*" in entry:
                        trans_dip[int(ii / 3)][ii % 3] = float("inf")
                    else:
                        trans_dip[int(ii / 3)][ii % 3] = float(entry)
            self.data['trans_dip'] = trans_dip

            freqs = np.zeros(len(temp_freqs) - temp_freqs.count('None'))
            for ii, entry in enumerate(temp_freqs):
                if entry != 'None':
                    if "*" in entry:
                        if ii == 0:
                            freqs[ii] = -float("inf")
                        elif ii == len(freqs) - 1:
                            freqs[ii] = float("inf")
                        elif freqs[ii - 1] == -float("inf"):
                            freqs[ii] = -float("inf")
                        elif "*" in temp_freqs[ii + 1]:
                            freqs[ii] = float("inf")
                        else:
                            raise RuntimeError(
                                "ERROR: Encountered an undefined frequency not at the beginning or end of the "
                                "frequency list, which makes no sense! Exiting...")
                        if not self.data.get('completion', []):
                            if "undefined_frequency" not in self.data["errors"]:
                                self.data["errors"] += ["undefined_frequency"]
                        else:
                            if "undefined_frequency" not in self.data["warnings"]:
                                self.data["warnings"]["undefined_frequency"] = True
                    else:
                        freqs[ii] = float(entry)
            self.data['frequencies'] = freqs

            intens = np.zeros(len(temp_intens) - temp_intens.count('None'))
            for ii, entry in enumerate(temp_intens):
                if entry != 'None':
                    if "*" in entry:
                        intens[ii] = float("inf")
                    else:
                        intens[ii] = float(entry)
            self.data['IR_intens'] = intens

            header_pattern = r"\s*Raman Active:\s+[YESNO]+\s+(?:[YESNO]+\s+)*X\s+Y\s+Z\s+(?:X\s+Y\s+Z\s+)*"
            table_pattern = r"\s*[a-zA-Z][a-zA-Z\s]\s*([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*(?:([\d\-\.]+)\s*" \
                            r"([\d\-\.]+)\s*([\d\-\.]+)\s*(?:([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+))*)*"
            footer_pattern = r"TransDip\s+\-?[\d\.\*]+\s*\-?[\d\.\*]+\s*\-?[\d\.\*]+\s*(?:\-?[\d\.\*]+\s*\-?" \
                             r"[\d\.\*]+\s*\-?[\d\.\*]+\s*)*"
            temp_freq_mode_vecs = read_table_pattern(
                self.text, header_pattern, table_pattern, footer_pattern)
            freq_mode_vecs = np.zeros(
                shape=(len(freqs), len(temp_freq_mode_vecs[0]), 3))

            for ii, triple_FMV in enumerate(temp_freq_mode_vecs):
                for jj, line in enumerate(triple_FMV):
                    for kk, entry in enumerate(line):
                        if entry != 'None':
                            freq_mode_vecs[int(ii * 3 + math.floor(kk / 3)),
                                           jj, kk % 3] = float(entry)

            self.data["frequency_mode_vectors"] = freq_mode_vecs
            freq_length = len(self.data["frequencies"])
            if len(self.data["frequency_mode_vectors"]) != freq_length or len(
                    self.data["IR_intens"]) != freq_length or len(self.data["IR_active"]) != freq_length:
                self.data["warnings"]["frequency_length_inconsistency"] = True

    def _read_single_point_data(self):
        """
        Parses final free energy information from single-point calculations.
        """
        temp_dict = read_pattern(
            self.text, {
                "final_energy":
                    r"\s*SCF\s+energy in the final basis set\s+=\s*([\d\-\.]+)"
            })

        if temp_dict.get('final_energy') is None:
            self.data['final_energy'] = None
        else:
            # -1 in case of pcm
            # Two lines will match the above; we want final calculation
            self.data['final_energy'] = float(temp_dict.get('final_energy')[-1][0])

    def _read_force_data(self):
        self._read_gradients()

    def _read_pcm_information(self):
        """
        Parses information from PCM solvent calculations.
        """

        temp_dict = read_pattern(
            self.text, {
                "g_electrostatic": r"\s*G_electrostatic\s+=\s+([\d\-\.]+)\s+hartree\s+=\s+([\d\-\.]+)\s+kcal/mol\s*",
                "g_cavitation": r"\s*G_cavitation\s+=\s+([\d\-\.]+)\s+hartree\s+=\s+([\d\-\.]+)\s+kcal/mol\s*",
                "g_dispersion": r"\s*G_dispersion\s+=\s+([\d\-\.]+)\s+hartree\s+=\s+([\d\-\.]+)\s+kcal/mol\s*",
                "g_repulsion": r"\s*G_repulsion\s+=\s+([\d\-\.]+)\s+hartree\s+=\s+([\d\-\.]+)\s+kcal/mol\s*",
                "total_contribution_pcm": r"\s*Total\s+=\s+([\d\-\.]+)\s+hartree\s+=\s+([\d\-\.]+)\s+kcal/mol\s*",
                "solute_internal_energy": r"Solute Internal Energy \(H0\)\s*=\s*([\d\-\.]+)"
            }
        )

        for key in temp_dict:
            if temp_dict.get(key) is None:
                self.data["solvent_data"][key] = None
            elif len(temp_dict.get(key)) == 1:
                self.data["solvent_data"][key] = float(temp_dict.get(key)[0][0])
            else:
                temp_result = np.zeros(len(temp_dict.get(key)))
                for ii, entry in enumerate(temp_dict.get(key)):
                    temp_result[ii] = float(entry[0])
                self.data["solvent_data"][key] = temp_result

        smd_keys = ["smd0", "smd3", "smd4", "smd6", "smd9"]
        for key in smd_keys:
            self.data["solvent_data"][key] = None

    def _read_smd_information(self):
        """
        Parses information from SMD solvent calculations.
        """

        temp_dict = read_pattern(
            self.text, {
                "smd0": r"E-EN\(g\) gas\-phase elect\-nuc energy\s*([\d\-\.]+) a\.u\.",
                "smd3": r"G\-ENP\(liq\) elect\-nuc\-pol free energy of system\s*([\d\-\.]+) a\.u\.",
                "smd4": r"G\-CDS\(liq\) cavity\-dispersion\-solvent structure\s*free energy\s*([\d\-\.]+) kcal\/mol",
                "smd6": r"G\-S\(liq\) free energy of system\s*([\d\-\.]+) a\.u\.",
                "smd9": r"DeltaG\-S\(liq\) free energy of\s*solvation\s*\(9\) = \(6\) \- \(0\)\s*([\d\-\.]+) kcal\/mol"
            }
        )

        for key in temp_dict:
            if temp_dict.get(key) is None:
                self.data["solvent_data"][key] = None
            elif len(temp_dict.get(key)) == 1:
                self.data["solvent_data"][key] = float(temp_dict.get(key)[0][0])
            else:
                temp_result = np.zeros(len(temp_dict.get(key)))
                for ii, entry in enumerate(temp_dict.get(key)):
                    temp_result[ii] = float(entry[0])
                self.data["solvent_data"][key] = temp_result

        pcm_keys = ["g_electrostatic", "g_cavitation", "g_dispersion", "g_repulsion", "total_contribution_pcm",
                    "solute_internal_energy"]
        for key in pcm_keys:
            self.data["solvent_data"][key] = None

    def _check_completion_errors(self):
        """
        Parses potential errors that can cause jobs to crash
        """
        if read_pattern(
                self.text, {
                    "key":
                        r"Coordinates do not transform within specified threshold"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["failed_to_transform_coords"]
        elif read_pattern(
                self.text,
                {
                    "key": r"The Q\-Chem input file has failed to pass inspection"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["input_file_error"]
        elif read_pattern(
                self.text, {
                    "key": r"Error opening input stream"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["failed_to_read_input"]
        elif read_pattern(
                self.text, {
                    "key": r"FileMan error: End of file reached prematurely"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["premature_end_FileMan_error"]
        elif read_pattern(
                self.text, {
                    "key": r"method not available"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["method_not_available"]
        elif read_pattern(
                self.text, {
                    "key": r"Could not find \$molecule section in ParseQInput"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["read_molecule_error"]
        elif read_pattern(
                self.text, {
                    "key": r"Welcome to Q-Chem"
                },
                terminate_on_match=True).get('key') != [[]]:
            self.data["errors"] += ["never_called_qchem"]
        elif read_pattern(
                self.text, {
                    "key": r"\*\*\*ERROR\*\*\* Hessian Appears to have all zero or negative eigenvalues"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["hessian_eigenvalue_error"]
        elif read_pattern(
                self.text, {
                    "key": r"FlexNet Licensing error"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["licensing_error"]
        elif read_pattern(
                self.text, {
                    "key": r"Unable to validate license"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["licensing_error"]
        elif read_pattern(
                self.text, {
                    "key": r"Could not open driver file in ReadDriverFromDisk"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["driver_error"]
        elif read_pattern(
                self.text, {
                    "key": r"gen_scfman_exception:  GDM:: Zero or negative preconditioner scaling factor"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["gdm_neg_precon_error"]
        else:
            tmp_failed_line_searches = read_pattern(
                self.text, {
                    "key": r"\d+\s+failed line searches\.\s+Resetting"
                },
                terminate_on_match=False).get('key')
            if tmp_failed_line_searches is not None:
                if len(tmp_failed_line_searches) > 10:
                    self.data["errors"] += ["SCF_failed_to_converge"]
        if self.data.get("errors") == []:
            self.data["errors"] += ["unknown_error"]

    def as_dict(self):
        d = {}
        d["data"] = self.data
        d["text"] = self.text
        d["filename"] = self.filename
        return jsanitize(d, strict=True)


def check_for_structure_changes(mol1, mol2):
    # special_elements = ["Li", "Na", "Mg", "Ca", "Zn"]
    # mol_list = [copy.deepcopy(mol1), copy.deepcopy(mol2)]

    if mol1.composition != mol2.composition:
        raise RuntimeError("Molecules have different compositions! Exiting...")

    # for ii, site in enumerate(mol1):
    #     if site.specie.symbol != mol2[ii].specie.symbol:
    #         print(
    #             "WARNING: Comparing molecules with different atom ordering! Turning off special treatment for "
    #             "coordinating metals.")
    #         special_elements = []

    # special_sites = [[], []]
    # for ii, mol in enumerate(mol_list):
    #     for jj, site in enumerate(mol):
    #         if site.specie.symbol in special_elements:
    #             distances = [[kk, site.distance(other_site)] for kk, other_site in enumerate(mol)]
    #             special_sites[ii].append([jj, site, distances])
    #     for jj, site in enumerate(mol):
    #         if site.specie.symbol in special_elements:
    #             mol.__delitem__(jj)

    # Can add logic to check the distances in the future if desired

    initial_mol_graph = MoleculeGraph.with_local_env_strategy(mol1,
                                                              OpenBabelNN(),
                                                              reorder=False,
                                                              extend_structure=False)
    initial_mol_graph = metal_edge_extender(initial_mol_graph)
    initial_graph = initial_mol_graph.graph
    last_mol_graph = MoleculeGraph.with_local_env_strategy(mol2,
                                                           OpenBabelNN(),
                                                           reorder=False,
                                                           extend_structure=False)
    last_mol_graph = metal_edge_extender(last_mol_graph)
    last_graph = last_mol_graph.graph
    if initial_mol_graph.isomorphic_to(last_mol_graph):
        return "no_change"
    else:
        if nx.is_connected(initial_graph.to_undirected()) and not nx.is_connected(last_graph.to_undirected()):
            return "unconnected_fragments"
        elif last_graph.number_of_edges() < initial_graph.number_of_edges():
            return "fewer_bonds"
        elif last_graph.number_of_edges() > initial_graph.number_of_edges():
            return "more_bonds"
        else:
            return "bond_change"
