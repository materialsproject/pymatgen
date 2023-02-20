# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Parsers for XTB output files and directories
"""
from __future__ import annotations

import logging
import os
import re

from monty.json import MSONable

from pymatgen.core import Molecule
from pymatgen.io.xyz import XYZ

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Alex Epstein"
__email__ = "aepstein@lbl.gov"
__credits__ = "Sam Blau, Evan Spotte-Smith"

logger = logging.getLogger(__name__)


class CRESTOutput(MSONable):
    """
    Class to parse CREST output files
    """

    def __init__(self, output_filename, path="."):
        """
        Assumes runtype is iMTD-GC [default]

        Args:
            output_filename (str): Filename to parse
            path (str): Path to directory including output_filename and all
                other xtb output files (crest_best.xyz, etc.)
        """
        self.path = path
        self.filename = output_filename

        self.cmd_options = {}
        self.sorted_structures_energies = []
        self.properly_terminated = False
        self._parse_crest_output()

    def _parse_crest_output(self):
        """
        Parse output file and directory to extract all command line inputs
            and output files.
        Sets the attributes:
            cmd_options: Dict of type {flag: value}
            sorted_structrues_energies: n x m x 2 list, for n conformers,
                m rotamers per conformer, and tuple of
                [Molecule, energy]
            properly_terminated: True or False if run properly terminated
        """
        output_filepath = os.path.join(self.path, self.filename)

        # Get CREST command
        crest_cmd = None
        with open(output_filepath) as xtbout_file:
            for line in xtbout_file:
                if "> crest" in line:
                    crest_cmd = line.strip()[8:]
                    break

        split_cmd = crest_cmd.split(" ")

        # Get input file if present
        try:
            self.coord_file = os.path.join(self.path, split_cmd[0])
            self.input_structure = Molecule.from_file(filename=self.coord_file)
        except FileNotFoundError:
            print(f"Input file {split_cmd[0]} not found")

        # Get CREST input flags
        for i, entry in enumerate(split_cmd):
            value = None
            if entry and "-" in entry:
                option = entry[1:]
                if i + 1 < len(split_cmd) and "-" not in split_cmd[i + 1]:
                    value = split_cmd[i + 1]
                self.cmd_options[option] = value
        # Get input charge for decorating parsed molecules
        chg = 0
        if "chrg" in self.cmd_options:
            str_chg = self.cmd_options["chrg"]
            chg = int(str_chg) if "-" in str_chg else int(str_chg[-1])
        elif "c" in self.cmd_options:
            str_chg = self.cmd_options["c"]
            chg = int(str_chg) if "-" in str_chg else int(str_chg[-1])

        # Check for proper termination
        with open(output_filepath, "rb+") as xtbout_file:
            xtbout_file.seek(-2, 2)
            while xtbout_file.read(1) != b"\n":
                xtbout_file.seek(-2, 1)
            end_bstring = xtbout_file.read()
            if b"CREST terminated normally." in end_bstring:
                self.properly_terminated = True

        if self.properly_terminated:
            # Parse for number of conformers and rotamers
            conformer_pattern = re.compile(
                r"\s+\d+\s+(?P<Erel>\d*\.\d*)\s+(?P<Etot>-*\d+\.\d+)\s+"
                r"(?P<weight>-*\d+\.\d+)\s+"
                r"(?P<conformer>-*\d+\.\d+)\s+(?P<set>\d+)\s+(?P<degen>\d+)\s+"
                r"(?P<origin>\w+)\n"
            )
            rotamer_pattern = re.compile(
                r"\s+\d+\s+(?P<Erel>\d*\.\d*)\s+(?P<Etot>-*\d+\.\d+)\s+"
                r"(?P<weight>-*\d+\.\d+)\s+"
                r"(?P<origin>\w+)\n"
            )
            conformer_degeneracies = []
            energies = []
            with open(output_filepath) as xtbout_file:
                for line in xtbout_file:
                    conformer_match = conformer_pattern.match(line)
                    rotamer_match = rotamer_pattern.match(line)
                    if conformer_match:
                        conformer_degeneracies.append(int(conformer_match["degen"]))
                        energies.append(conformer_match["Etot"])
                    elif rotamer_match:
                        energies.append(rotamer_match["Etot"])
            # Get final rotamers file and read in all molecules,
            # sorted by conformer type and energy
            if "crest_rotamers.xyz" in os.listdir(self.path):
                final_rotamer_filename = "crest_rotamers.xyz"
            else:
                n_rot_files = []
                for f in os.listdir(self.path):
                    if "crest_rotamers" in f:
                        n_rot_file = int(os.path.splitext(f)[0].split("_")[2])
                        n_rot_files.append(n_rot_file)
                if len(n_rot_files) > 0:
                    final_rotamer_filename = f"crest_rotamers_{max(n_rot_files)}.xyz"
            try:
                rotamers_path = os.path.join(self.path, final_rotamer_filename)
                rotamer_structures = XYZ.from_file(rotamers_path).all_molecules
                for r in rotamer_structures:
                    r.set_charge_and_spin(charge=chg)
                start = 0
                for n, d in enumerate(conformer_degeneracies):
                    self.sorted_structures_energies.append([])
                    i = 0
                    for i in range(start, start + d):
                        self.sorted_structures_energies[n].append([rotamer_structures[i], energies[i]])
                    start = i + 1
            except FileNotFoundError:
                print(f"{final_rotamer_filename} not found, no rotamer list processed")

            # Get lowest energy conformer from 'crest_best.xyz'
            crestbest_path = os.path.join(self.path, "crest_best.xyz")
            try:
                lowest_e_struct = Molecule.from_file(crestbest_path)
                lowest_e_struct.set_charge_and_spin(charge=chg)
                self.lowest_energy_structure = lowest_e_struct
            except FileNotFoundError:
                print(f"{crestbest_path} not found")

        else:
            crestbest_path = os.path.join(self.path, "crest_best.xyz")
            try:
                lowest_e_struct = Molecule.from_file(crestbest_path)
                lowest_e_struct.set_charge_and_spin(charge=chg)
                self.lowest_energy_structure = lowest_e_struct
            except FileNotFoundError:
                print(f"{crestbest_path} not found")
