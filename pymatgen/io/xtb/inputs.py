# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import logging
import re
from monty.json import MSONable
from monty.io import zopen

# Classes for reading/manipulating/writing XTB Input and Output files

__author__ = "Evan Spotte-Smith, Shyam Dwaraknath"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__email__ = "ewcspottesmith@lbl.gov"

logger = logging.getLogger(__name__)


def xcontrol_eval(string):
    """
    Convert possible xTB input values into appropriate Python types

    Args:
        string: str to be evaluated
    return:
        Depends on type of value held in string

    """

    if "true" in string.lower():
        return True

    if "false" in string.lower():
        return False

    for val_type in (int, float):
        try:
            return val_type(string)
        except:
            pass

    return string


class XTBInput(MSONable):
    """
    An object representing a XTB Input file.
    """

    def __init__(self, blocks):
        """
        Args:
            blocks: dict mapping block names to blocks of input data
        """
        self.blocks = blocks

    def __repr__(self):
        data = []
        for block_name, block in self.blocks.items():
            if isinstance(block, dict):
                # For structured, multi-value blocks
                data.append(f"${block_name.lower()}")
                for key, val in block.items():
                    data.append(f"   {key.lower()} = {str(val)}")
            else:
                # For single-value blocks
                data.append(f"${block_name.lower()} {str(block)}")
        data.append("$end")

        return "\n".join(data)

    @classmethod
    def from_file(cls, filename):
        """
        Builds an XTBInput object from a file

        Args:
            filename: path str for the xTB input file.
        """
        lines = list()
        with zopen(filename) as f:
            lines = f.readlines()

        block_regex = re.compile(r"\s*\$(\w*)")
        data_regex = re.compile(r"\s*(\w.*)=\s*(.*)\s*\n")
        blocks = [
            (line_number, block_regex.match(line).group(1))
            for line_number, line in enumerate(lines)
            if "$" in line
        ]

        block_names = [b[1] for b in blocks]
        block_beginnings = [b[0] + 1 for b in blocks][:-1]
        block_endings = [b[0] for b in blocks][1:]

        block_data = {}
        for name, start, end in zip(block_names, block_beginnings, block_endings):
            block_data[name] = {}
            for line_number in range(start, end, 1):
                match = data_regex.match(lines[line_number])
                key = match.group(1)
                val = match.group(2)
                block_data[name][key] = xcontrol_eval(val.title())

        return cls(blocks=block_data)

    @classmethod
    def from_defaults(cls):
        """
        Construct an XTBInput object based on some basic default values
        """

        blocks = dict()

        blocks["gfn"] = {"method": 2,
                         "scc": True,
                         "periodic": False}
        blocks["scc"] = {"maxiterations": 250,
                         "temp": 298.15,
                         "broydamp": 0.40}
        blocks["thermo"] = {"temp": 298.15}
        blocks["chrg"] = 0
        blocks["write"] = {"density": False,
                           "charges": True,
                           "mulliken": True,
                           "distances": True,
                           "angles": True,
                           "torsions": True}

        return cls(blocks)

    def to_file(self, path, filename, charge=None, unpaired=None):
        """
        Write xTB inputs to file.

        Args:
            path: path str to a directory in which to print input files
            filename: str for the main xTB input file
            charge: If not None (default), then this will be the
            unpaired: If not None (default)
        returns:
            None
        """
        contents = self.__repr__()

        with open(os.path.join(path, filename), 'w') as file:
            file.write(contents)

        if charge is not None:
            if isinstance(charge, int) or isinstance(charge, float):
                with open(os.path.join(path, ".CHRG"), 'w') as file:
                    file.write(str(int(charge)))
            else:
                raise TypeError("unpaired must be an integer (preferred) or float!")

        if unpaired is not None:
            if isinstance(unpaired, int) or isinstance(charge, float):
                if unpaired > 0:
                    with open(os.path.join(path, ".UHF"), 'w') as file:
                        file.write(str(int(unpaired)))
                else:
                    raise ValueError("unpaired must be > 0!")
            else:
                raise TypeError("unpaired must be an integer (preferred) or float!")