# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
import argparse
import os

import numpy as np
import re

from monty.json import MSONable
from monty.io import zopen

from pymatgen.core import Molecule
from pymatgen.io.xyz import XYZ
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.qchem.utils import (read_table_pattern,
                                     read_pattern,
                                     lower_and_check_unique)

# Classes for reading/manipulating/writing QChem output files.

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"

logger = logging.getLogger(__name__)


class XTBOutput(MSONable):
    """
    Class to parse XTB outputs
    """

    def __init__(self, path, output_filename, namespace=None):
        """
        Args:
            path (str): Path to directory including output_filename and all
                other xtb output files (xtbopt.log, vibspectrum, etc.)
            output_filename (str): Filename to parse
            namespace (str): If the namespace is not None (default), this will
                be prepended to all filenames other than output_filename
        """

        self.path = path
        self.filename = output_filename
        self.namespace = namespace
        self.data = dict()

        self.data["errors"] = list()
        self.data["warnings"] = dict()

        self.text = ""
        with zopen(output_filename, "rt") as f:
            self.text = f.read()

        self.data["input"] = dict()
        self.data["output"] = dict()

        self._read_setup()
        self._read_parameters()


        # Need some condition to see if this is an opt job

    def _read_setup(self):
        """
        Read calculation setup to extract, for instance, the command used
        and the namespace

        Note that this does not attempt to extract all information from the
            setup.
        """

        parser = argparse.ArgumentParser()
        parser.add_argument("geometry")
        parser.add_argument("-c", "--chrg", type=int, nargs=1)
        parser.add_argument("-u", "--uhf", type=int, nargs=1)
        parser.add_argument("-a", "--acc", type=float, nargs=1)
        parser.add_argument("--gfn", type=int, nargs=1)
        parser.add_argument("-g", "--gbsa", nargs="+")
        parser.add_argument("--scc", "--sp", action="store_true")
        parser.add_argument("--grad", action="store_true")
        parser.add_argument("-o", "--opt", nargs="?", default="normal")
        parser.add_argument("--hess", action="store_true")
        parser.add_argument("--pop", action="store_true")
        parser.add_argument("--ohess", nargs="?", default="normal")
        parser.add_argument("-I", "--input", nargs=1, type=str)
        parser.add_argument("--namespace", nargs=1, type=str)
        parser.add_argument("-P", "--parallel", nargs=1, type=int)

        temp_dict = read_pattern(
        self.text, {
            "program_call": r"\s+program call\s+:\s+xtb (.+)\n",
            "coordinate file": r"\s+coordinate file\s+:\s+([:\\\.\-/_0-9a-zA-Z]+)\n",
            "charge": r"\s+charge\s+:\s+([\-0-9]+)\n",
            "spin": r"\s+spin\s+:\s+([\.0-9]+)\n"
        })

        if temp_dict.get("program_call") is None:
            self.data["setup"]["input"]["program_call"] = None
        else:
            namespace = parser.parse_args(temp_dict.get("program_call")[0][0].split())

            self.data["input"]["accuracy"] = namespace.acc


            if namespace.gbsa is None:
                self.data["input"]["gbsa"] = None
            else:
                self.data["input"]["gbsa"] = dict()
                if len(namespace.gbsa) == 1:
                    self.data["input"]["gbsa"]["solvent"] = namespace.gbsa[0]
                elif len(namespace.gbsa) == 2:
                    self.data["input"]["gbsa"]["solvent"] = namespace.gbsa[0]
                    self.data["input"]["gbsa"]["grid"] = namespace.gbsa[1]
                else:
                    self.data["input"]["gbsa"]["solvent"] = namespace.gbsa[0]
                    self.data["input"]["gbsa"]["grid"] = namespace.gbsa[1]
                    self.data["input"]["gbsa"]["other"] = namespace.gbsa[2:]







class VibspectrumOutput(MSONable):
    """
    Class to parse vibspectrum files from xTB frequency output.
    """
    pass


class OptLogOutput(MSONable):
    """
    Class to parse optimization log files from xTB optimization output.
    """
    pass


class FakeGaussOutput(MSONable):
    """
    Class to parse "g98" files from xTB frequency output.
    """
    pass

class CRESTOutput(MSONable):

    def __init__(self, path, output_filename):
        """
        Currently assumes runtype is iMTD-GC [default]
        Args:
            path (str): Path to directory including output_filename and all
                other xtb output files (xtbopt.log, vibspectrum, etc.)
            output_filename (str): Filename to parse
        """

        self.path = path
        self.filename = output_filename
        self.cmd_options = dict()
        self.data = dict()
        self.data = dict()
        self.sorted_structures_energies = []

        self._parse_crest_output()

    def _parse_crest_output(self):
        """
        Parse output file to extract all command line inputs and output files
        Needs functionality for incomplete jobs, error processing
        """
        output_filepath = os.path.join(self.path, self.filename)

        crest_cmd = None
        with open(output_filepath, 'r') as xtbout_file:
            for line in xtbout_file:
                if '> crest' in line:
                    crest_cmd = line.strip()[8:]
                    break

        split_cmd = crest_cmd.split(' ')
        self.coord_file = os.path.join(self.path, split_cmd[0])
        self.input_structure = Molecule.from_file(filename=self.coord_file)
        for i, entry in enumerate(split_cmd):
            value = None
            option = None
            if entry:
                if '-' in entry:
                    option = entry[1:]
                    if i + 1 < len(split_cmd):
                        if '-' not in split_cmd[i + 1]:
                            value = split_cmd[i + 1]
                    self.cmd_options[option] = value

        with open(output_filepath, 'rb+') as xtbout_file:
            xtbout_file.seek(-2, 2)
            while xtbout_file.read(1) != b"\n":
                xtbout_file.seek(-2, 1)
            end_bstring = xtbout_file.read()
            if b'CREST terminated normally.' in end_bstring:
                self.properly_terminated = True

        if self.properly_terminated:
            self.lowest_energy_structure = Molecule.from_file(self.path+'/'+'crest_best.xyz')

            rotamer_structures = XYZ.from_file(os.path.join(self.path, 'crest_conformers.xyz')).all_molecules

            conformer_pattern = re.compile(
                r"\s+\d+\s+(?P<Erel>\d*\.\d*)\s+(?P<Etot>-*\d+\.\d+)\s+(?P<weight>-*\d+\.\d+)\s+(?P<conformer>-*\d+\.\d+)\s+(?P<set>\d+)\s+(?P<degen>\d+)")
            rotamer_pattern = re.compile(
                r"\s+\d+\s+(?P<Erel>\d*\.\d*)\s+(?P<Etot>-*\d+\.\d+)\s+(?P<weight>-*\d+\.\d+)\s+\w+\n")
            conformer_degeneracies = []
            energies = []
            with open(output_filepath, 'r') as xtbout_file:
                for line in xtbout_file:
                    conformer_match = conformer_pattern.match(line)
                    rotamer_match = rotamer_pattern.match(line)
                    if conformer_match:
                        conformer_degeneracies.append(int(conformer_match['degen']))
                        energies.append(conformer_match['Etot'])
                    elif rotamer_match:
                        energies.append(rotamer_match['Etot'])
            start = 0
            for n, d in enumerate(conformer_degeneracies):
                self.sorted_structures_energies.append([])
                for i in range(start,start+d):
                   self.sorted_structures_energies[n].append([rotamer_structures[i], energies[i]])
                start = i



 # def parse_xtb_output(file_path="xtb.out"):
#     """
#     Things we need to parse:
#     - Final energy
#     - Final enthalpy
#     - Final entropy
#     - Final free energy
#     - Heat capacity
#     - Charge?
#
#     :param file_path:
#     :return:
#     """
#
#     # In Hartree
#     energy_pattern = re.compile(r"\s+\| TOTAL ENERGY\s+(?P<energy>[\-\.0-9]+) Eh")
#
#     # In hartree
#     tot_enthalpy = re.compile(r"\s+\| TOTAL ENTHALPY\s+(?P<total_enthalpy>[\-\.0-9]+) Eh")
#
#     # In hartree
#     tot_gibbs = re.compile(r"\s+\| TOTAL FREE ENERGY\s+(?P<total_free_energy>[\-\.0-9]+) Eh")
#
#     # In cal/mol, cal/mol-K, cal/mol-K, respectively
#     tot_pattern = re.compile(r"TOT\s+(?P<enthalpy>[\-\.0-9]+)\s+(?P<heat_capacity>[\-\.0-9]+)\s+(?P<entropy>[\-\.0-9]+)")
#
#     with open(file_path, 'r') as xtbout_file:
#         contents = xtbout_file.read()
#
#         total_energy = float(energy_pattern.search(contents).group("energy"))
#         total_enthalpy = float(tot_enthalpy.search(contents).group("total_enthalpy"))
#         total_gibbs = float(tot_gibbs.search(contents).group("total_free_energy"))
#
#         tot_match = tot_pattern.search(contents)
#         enthalpy = float(tot_match.group("enthalpy")) / 1000
#         heat_capacity = float(tot_match.group("heat_capacity"))
#         entropy = float(tot_match.group("entropy"))
#
#         return {"energy": total_energy,
#                 "total_h": total_enthalpy,
#                 "total_g": total_gibbs,
#                 "enthalpy": enthalpy,
#                 "cp": heat_capacity,
#                 "entropy": entropy}