# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import re
import logging
import os
import numpy as np
import math

from monty.io import zopen
from monty.json import jsanitize
from monty.json import MSONable
from pymatgen.core import Molecule

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
import networkx as nx
from pymatgen.io.babel import BabelMolAdaptor
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
        self.text = ""
        with zopen(filename, 'rt') as f:
            self.text = f.read()

        # Check if output file contains multiple output files. If so, print an error message and exit
        self.data["multiple_outputs"] = read_pattern(
            self.text, {
                "key": r"Job\s+\d+\s+of\s+(\d+)\s+"
            },
            terminate_on_match=True).get('key')
        if not (self.data.get('multiple_outputs') == None
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
        if self.data.get('charge') is not None:
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

        # Parse the Mulliken charges
        self._read_mulliken()

        # Parse PCM information
        self._read_pcm_information()

        # Parse the final energy
        temp_final_energy = read_pattern(
            self.text, {
                "key": r"Final\senergy\sis\s+([\d\-\.]+)"
            }).get('key')
        if temp_final_energy == None:
            self.data["final_energy"] = None
        else:
            self.data["final_energy"] = float(temp_final_energy[0][0])

        # Parse the S2 values in the case of an unrestricted calculation
        if self.data.get('unrestricted', []):
            temp_S2 = read_pattern(self.text, {
                "key": r"<S\^2>\s=\s+([\d\-\.]+)"
            }).get('key')
            if temp_S2 == None:
                self.data["S2"] = None
            elif len(temp_S2) == 1:
                self.data["S2"] = float(temp_S2[0][0])
            else:
                real_S2 = np.zeros(len(temp_S2))
                for ii, entry in enumerate(temp_S2):
                    real_S2[ii] = float(entry[0])
                self.data["S2"] = real_S2

        # Check if the calculation is a geometry optimization. If so, parse the relevant output
        self.data["optimization"] = read_pattern(
            self.text, {
                "key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*opt"
            }).get('key')
        if self.data.get('optimization', []):
            temp_energy_trajectory = read_pattern(
                self.text, {
                    "key": r"\sEnergy\sis\s+([\d\-\.]+)"
                }).get('key')
            if temp_energy_trajectory == None:
                self.data["energy_trajectory"] = []
            else:
                real_energy_trajectory = np.zeros(len(temp_energy_trajectory))
                for ii, entry in enumerate(temp_energy_trajectory):
                    real_energy_trajectory[ii] = float(entry[0])
                self.data["energy_trajectory"] = real_energy_trajectory
                self._read_last_geometry()
                if have_babel:
                    self._check_for_structure_changes()
                self._read_optimized_geometry()
                # Then, if no optimized geometry or z-matrix is found, and no errors have been previously
                # idenfied, check to see if the optimization failed to converge or if Lambda wasn't able
                # to be determined.
                if len(self.data.get("errors")) == 0 and self.data.get('optimized_geometry') is None \
                        and len(self.data.get('optimized_zmat')) == 0:
                    self._check_optimization_errors()

        # Check if the calculation contains a constraint in an $opt section.
        self.data["opt_constraint"] = read_pattern(self.text, {
            "key": r"\$opt\s+CONSTRAINT"
        }).get('key')
        if self.data.get('opt_constraint'):
            temp_constraint = read_pattern(
                self.text, {
                    "key":
                    r"Constraints and their Current Values\s+Value\s+Constraint\s+(\w+)\:\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
                }).get('key')
            if temp_constraint != None:
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
        if temp_charge != None:
            self.data["charge"] = int(temp_charge[0][0])
        else:
            temp_charge = read_pattern(
                self.text, {
                    "key": r"Sum of atomic charges \=\s+([\d\-\.\+]+)"
                },
                terminate_on_match=True).get('key')
            if temp_charge == None:
                self.data["charge"] = None
            else:
                self.data["charge"] = int(float(temp_charge[0][0]))

        temp_multiplicity = read_pattern(
            self.text, {
                "key": r"\$molecule\s+[\-\d]+\s+(\d)"
            },
            terminate_on_match=True).get('key')
        if temp_multiplicity != None:
            self.data["multiplicity"] = int(temp_multiplicity[0][0])
        else:
            temp_multiplicity = read_pattern(
                self.text, {
                    "key": r"Sum of spin\s+charges \=\s+([\d\-\.\+]+)"
                },
                terminate_on_match=True).get('key')
            if temp_multiplicity == None:
                self.data["multiplicity"] = 1
            else:
                self.data["multiplicity"] = int(
                    float(temp_multiplicity[0][0])) + 1

    def _read_species_and_inital_geometry(self):
        """
        Parses species and initial geometry.
        """
        header_pattern = r"Standard Nuclear Orientation \(Angstroms\)\s+I\s+Atom\s+X\s+Y\s+Z\s+-+"
        table_pattern = r"\s*\d+\s+([a-zA-Z]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*"
        footer_pattern = r"\s*-+"
        temp_geom = read_table_pattern(self.text, header_pattern,
                                       table_pattern, footer_pattern)
        if temp_geom == None or len(temp_geom) == 0:
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
            if temp_point_group != None:
                self.data["point_group"] = temp_point_group[0][0]
            else:
                self.data["point_group"] = None
            temp_geom = temp_geom[0]
            species = []
            geometry = np.zeros(shape=(len(temp_geom), 3), dtype=float)
            for ii, entry in enumerate(temp_geom):
                species += [entry[0]]
                for jj in range(3):
                    geometry[ii, jj] = float(entry[jj + 1])
            self.data["species"] = species
            self.data["initial_geometry"] = geometry
            self.data["initial_molecule"] = Molecule(
                species=species,
                coords=geometry,
                charge=self.data.get('charge'),
                spin_multiplicity=self.data.get('multiplicity'))

    def _read_SCF(self):
        """
        Parses both old and new SCFs.
        """
        if self.data.get('using_GEN_SCFMAN', []):
            if "SCF_failed_to_converge" in self.data.get("errors"):
                footer_pattern = r"^\s*gen_scfman_exception: SCF failed to converge"
            else:
                footer_pattern = r"^\s*\-+\n\s+SCF time"
            header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*(?:RMS Gradient)*\s+\-+(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n"
            table_pattern = r"(?:\s*Nonlocal correlation = [\d\-\.]+e[\d\-]+)*(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s+([\d\-\.]+)\s+([\d\-\.]+)e([\d\-\.\+]+)(?:\s+Convergence criterion met)*(?:\s+Preconditoned Steepest Descent)*(?:\s+Roothaan Step)*(?:\s+(?:Normal\s+)*BFGS [Ss]tep)*(?:\s+LineSearch Step)*(?:\s+Line search: overstep)*(?:\s+Dog-leg BFGS step)*(?:\s+Line search: understep)*(?:\s+Descent step)*"
        else:
            if "SCF_failed_to_converge" in self.data.get("errors"):
                footer_pattern = r"^\s*\d+\s*[\d\-\.]+\s+[\d\-\.]+E[\d\-\.]+\s+Convergence\s+failure\n"
            else:
                footer_pattern = r"^\s*\-+\n"
            header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+DIIS Error\s+\-+\n"
            table_pattern = r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s*([\d\-\.]+)\s+([\d\-\.]+)E([\d\-\.\+]+)(?:\s*\n\s*cpu\s+[\d\-\.]+\swall\s+[\d\-\.]+)*(?:\nin dftxc\.C, eleTot sum is:[\d\-\.]+, tauTot is\:[\d\-\.]+)*(?:\s+Convergence criterion met)*(?:\s+Done RCA\. Switching to DIIS)*(?:\n\s*Warning: not using a symmetric Q)*(?:\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+(?:\s*\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+)*)*"

        temp_scf = read_table_pattern(self.text, header_pattern, table_pattern,
                                      footer_pattern)
        real_scf = []
        for one_scf in temp_scf:
            temp = np.zeros(shape=(len(one_scf), 2))
            for ii, entry in enumerate(one_scf):
                temp[ii, 0] = float(entry[0])
                temp[ii, 1] = float(entry[1]) * 10**float(entry[2])
            real_scf += [temp]

        self.data["SCF"] = real_scf

    def _read_mulliken(self):
        """
        Parses Mulliken charges. Also parses spins given an unrestricted SCF.
        """
        if self.data.get('unrestricted', []):
            header_pattern = r"\-+\s+Ground-State Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+Spin\s\(a\.u\.\)\s+\-+"
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

    def _read_optimized_geometry(self):
        """
        Parses optimized XYZ coordinates. If not present, parses optimized Z-matrix.
        """
        header_pattern = r"\*+\s+OPTIMIZATION\s+CONVERGED\s+\*+\s+\*+\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+\w+\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Z-matrix Print:"

        parsed_optimized_geometry = read_table_pattern(
            self.text, header_pattern, table_pattern, footer_pattern)
        if parsed_optimized_geometry == [] or None:
            self.data["optimized_geometry"] = None
            header_pattern = r"^\s+\*+\s+OPTIMIZATION CONVERGED\s+\*+\s+\*+\s+Z-matrix\s+Print:\s+\$molecule\s+[\d\-]+\s+[\d\-]+\n"
            table_pattern = r"\s*(\w+)(?:\s+(\d+)\s+([\d\-\.]+)(?:\s+(\d+)\s+([\d\-\.]+)(?:\s+(\d+)\s+([\d\-\.]+))*)*)*(?:\s+0)*"
            footer_pattern = r"^\$end\n"

            self.data["optimized_zmat"] = read_table_pattern(
                self.text, header_pattern, table_pattern, footer_pattern)
        else:
            self.data["optimized_geometry"] = process_parsed_coords(
                parsed_optimized_geometry[0])
            if self.data.get('charge') != None:
                self.data["molecule_from_optimized_geometry"] = Molecule(
                    species=self.data.get('species'),
                    coords=self.data.get('optimized_geometry'),
                    charge=self.data.get('charge'),
                    spin_multiplicity=self.data.get('multiplicity'))

    def _read_last_geometry(self):
        """
        Parses the last geometry from an optimization trajectory for use in a new input file.
        """
        header_pattern = r"\s+Optimization\sCycle:\s+" + \
            str(len(self.data.get("energy_trajectory"))) + \
            r"\s+Coordinates \(Angstroms\)\s+ATOM\s+X\s+Y\s+Z"
        table_pattern = r"\s+\d+\s+\w+\s+([\d\-\.]+)\s+([\d\-\.]+)\s+([\d\-\.]+)"
        footer_pattern = r"\s+Point Group\:\s+[\d\w\*]+\s+Number of degrees of freedom\:\s+\d+"

        parsed_last_geometry = read_table_pattern(
            self.text, header_pattern, table_pattern, footer_pattern)
        if parsed_last_geometry == [] or None:
            self.data["last_geometry"] = None
        else:
            self.data["last_geometry"] = process_parsed_coords(
                parsed_last_geometry[0])
            if self.data.get('charge') != None:
                self.data["molecule_from_last_geometry"] = Molecule(
                    species=self.data.get('species'),
                    coords=self.data.get('last_geometry'),
                    charge=self.data.get('charge'),
                    spin_multiplicity=self.data.get('multiplicity'))

    def _check_for_structure_changes(self):
        initial_mol_graph = MoleculeGraph.with_local_env_strategy(self.data["initial_molecule"],
                                                                  OpenBabelNN(),
                                                                  reorder=False,
                                                                  extend_structure=False)
        initial_graph = initial_mol_graph.graph
        last_mol_graph = MoleculeGraph.with_local_env_strategy(self.data["molecule_from_last_geometry"],
                                                               OpenBabelNN(),
                                                               reorder=False,
                                                               extend_structure=False)
        last_graph = last_mol_graph.graph
        if initial_mol_graph.isomorphic_to(last_mol_graph):
            self.data["structure_change"] = "no_change"
        else:
            if nx.is_connected(initial_graph.to_undirected()) and not nx.is_connected(last_graph.to_undirected()):
                self.data["structure_change"] = "unconnected_fragments"
            elif last_graph.number_of_edges() < initial_graph.number_of_edges():
                self.data["structure_change"] = "fewer_bonds"
            elif last_graph.number_of_edges() > initial_graph.number_of_edges():
                self.data["structure_change"] = "more_bonds"
            else:
                self.data["structure_change"] = "bond_change"

    def _read_frequency_data(self):
        """
        Parses frequencies, enthalpy, entropy, and mode vectors.
        """
        temp_dict = read_pattern(
            self.text, {
                "frequencies":
                r"\s*Frequency:\s+([\d\-\.]+)(?:\s+([\d\-\.]+)(?:\s+([\d\-\.]+))*)*",
                "IR_intens":
                r"\s*IR Intens:\s+([\d\-\.]+)(?:\s+([\d\-\.]+)(?:\s+([\d\-\.]+))*)*",
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

        keys = ["ZPE", "trans_enthalpy", "rot_enthalpy", "vib_enthalpy", "gas_constant", "trans_entropy", "rot_entropy", "vib_entropy", "total_enthalpy", "total_entropy"]

        for key in keys:
            if temp_dict.get(key) == None:
                self.data[key] = None
            else:
                self.data[key] = float(temp_dict.get(key)[0][0])

        if temp_dict.get('frequencies') == None:
            self.data['frequencies'] = None
            self.data['IR_intens'] = None
            self.data['IR_active'] = None
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
            self.data['IR_active'] = active

            freqs = np.zeros(len(temp_freqs) - temp_freqs.count('None'))
            for ii, entry in enumerate(temp_freqs):
                if entry != 'None':
                    freqs[ii] = float(entry)
            self.data['frequencies'] = freqs

            intens = np.zeros(len(temp_intens) - temp_intens.count('None'))
            for ii, entry in enumerate(temp_intens):
                if entry != 'None':
                    intens[ii] = float(entry)
            self.data['IR_intens'] = intens

            header_pattern = r"\s*Raman Active:\s+[YESNO]+\s+(?:[YESNO]+\s+)*X\s+Y\s+Z\s+(?:X\s+Y\s+Z\s+)*"
            table_pattern = r"\s*[a-zA-Z][a-zA-Z\s]\s*([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*(?:([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+)\s*(?:([\d\-\.]+)\s*([\d\-\.]+)\s*([\d\-\.]+))*)*"
            footer_pattern = r"TransDip\s+[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+\s*(?:[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+\s*)*"
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

    def _read_single_point_data(self):
        """
        Parses final free energy information from single-point calculations.
        """
        temp_dict = read_pattern(
            self.text, {
                "final_energy":
                    r"\s*SCF\s+energy in the final basis set\s+=\s*([\d\-\.]+)"
            })

        if temp_dict.get('final_energy') == None:
            self.data['final_energy'] = None
        else:
            # -1 in case of pcm
            # Two lines will match the above; we want final calculation
            self.data['final_energy'] = float(temp_dict.get('final_energy')[-1][0])

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
            }
        )

        if temp_dict.get("g_electrostatic") is None:
            self.data["g_electrostatic"] = None
        else:
            self.data["g_electrostatic"] = float(temp_dict.get("g_electrostatic")[0][0])

        if temp_dict.get("g_cavitation") is None:
            self.data["g_cavitation"] = None
        else:
            self.data["g_cavitation"] = float(temp_dict.get("g_cavitation")[0][0])

        if temp_dict.get("g_dispersion") is None:
            self.data["g_dispersion"] = None
        else:
            self.data["g_dispersion"] = float(temp_dict.get("g_dispersion")[0][0])

        if temp_dict.get("g_repulsion") is None:
            self.data["g_repulsion"] = None
        else:
            self.data["g_repulsion"] = float(temp_dict.get("g_repulsion")[0][0])

        if temp_dict.get("total_contribution_pcm") is None:
            self.data["total_contribution_pcm"] = []
        else:
            self.data["total_contribution_pcm"] = float(temp_dict.get("total_contribution_pcm")[0][0])

    def _check_optimization_errors(self):
        """
        Parses three potential optimization errors: failing to converge within the allowed number
        of optimization cycles, failure to determine the lamda needed to continue, and inconsistent
        size of MO files due to a linear dependence in the AO basis.
        """
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
        elif read_pattern(
                self.text, {
                    "key": r"Inconsistent size for SCF MO coefficient file"
                },
                terminate_on_match=True).get('key') == [[]]:
            self.data["errors"] += ["linear_dependent_basis"]

    def _check_completion_errors(self):
        """
        Parses four potential errors that can cause jobs to crash: inability to transform
        coordinates due to a bad symmetric specification, an input file that fails to pass
        inspection, and errors reading and writing files.
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
            self.data["errors"] += ["IO_error"]
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
        else:
            self.data["errors"] += ["unknown_error"]

    def as_dict(self):
        d = {}
        d["data"] = self.data
        d["text"] = self.text
        d["filename"] = self.filename
        return jsanitize(d, strict=True)