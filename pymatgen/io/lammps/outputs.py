# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals,\
    absolute_import

import re
import glob
from io import StringIO

import numpy as np
import pandas as pd

from monty.json import MSONable
from monty.io import zopen

from pymatgen.io.lammps.data import LammpsBox

"""
This module implements classes and methods for processing LAMMPS output
files (log and dump).

"""

__author__ = "Kiran Mathew, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"


class LammpsDump(MSONable):

    """
    Object for representing dump data for a single snapshot.

    """

    def __init__(self, timestep, natoms, box, data):
        """
        Base constructor.

        Args:
            timestep (int): Current timestep.
            natoms (int): Total number of atoms in the box.
            box (LammpsBox): Simulation box.
            data (pd.DataFrame): Dumped atomic data.

        """
        self.timestep = timestep
        self.natoms = natoms
        self.box = box
        self.data = data

    @classmethod
    def from_string(cls, string):
        """
        Constructor from string parsing.

        Args:
            string (str): Input string.

        """
        lines = string.split("\n")
        timestep = int(lines[1])
        natoms = int(lines[3])
        box_arr = np.loadtxt(StringIO("\n".join(lines[5:8])))
        bounds = box_arr[:, :2]
        tilt = None
        if "xy xz yz" in lines[4]:
            tilt = box_arr[:, 2]
            x = (0, tilt[0], tilt[1], tilt[0] + tilt[1])
            y = (0, tilt[2])
            bounds -= np.array([[min(x), max(x)], [min(y), max(y)], [0, 0]])
        box = LammpsBox(bounds, tilt)
        data_head = lines[8].replace("ITEM: ATOMS", "").split()
        data = pd.read_csv(StringIO("\n".join(lines[9:])), names=data_head,
                           delim_whitespace=True)
        return cls(timestep, natoms, box, data)

    @classmethod
    def from_dict(cls, d):
        items = {"timestep": d["timestep"], "natoms": d["natoms"]}
        items["box"] = LammpsBox.from_dict(d["box"])
        items["data"] = pd.read_json(d["data"], orient="split")
        return cls(**items)

    def as_dict(self):
        d = dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["timestep"] = self.timestep
        d["natoms"] = self.natoms
        d["box"] = self.box.as_dict()
        d["data"] = self.data.to_json(orient="split")
        return d


def parse_lammps_dumps(file_pattern):
    """
    Generator that parses dump file(s).

    Args:
        file_pattern (str): Filename to parse. The timestep wildcard
            (e.g., dump.atom.'*') is supported and the files are parsed
            in the sequence of timestep.

    Yields:
        LammpsDump for each available snapshot.

    """
    files = glob.glob(file_pattern)
    if len(files) > 1:
        pattern = r"%s" % file_pattern.replace("*", "([0-9]+)")
        pattern = pattern.replace("\\", "\\\\")
        files = sorted(files,
                       key=lambda f: int(re.match(pattern, f).group(1)))

    for fname in files:
        with zopen(fname, "rt") as f:
            dump_cache = []
            for line in f:
                if line.startswith("ITEM: TIMESTEP"):
                    if len(dump_cache) > 0:
                        yield LammpsDump.from_string("".join(dump_cache))
                    dump_cache = [line]
                else:
                    dump_cache.append(line)
            yield LammpsDump.from_string("".join(dump_cache))


def parse_lammps_log(filename="log.lammps"):
    """
    Parses log file with focus on thermo data. Both one and multi line
    formats are supported. Any incomplete runs (no "Loop time" marker)
    will not be parsed.

    Notes:
        SHAKE stats printed with thermo data are not supported yet.
        They are ignored in multi line format, while they may cause
        issues with dataframe parsing in one line format.

    Args:
        filename (str): Filename to parse.

    Returns:
        [pd.DataFrame] containing thermo data for each completed run.

    """
    with open(filename) as f:
        lines = f.readlines()
    begin_flag = ("Memory usage per processor =",
                  "Per MPI rank memory allocation (min/avg/max) =")
    end_flag = "Loop time of"
    begins, ends = [], []
    for i, l in enumerate(lines):
        if l.startswith(begin_flag):
            begins.append(i)
        elif l.startswith(end_flag):
            ends.append(i)

    def _parse_thermo(lines):
        multi_pattern = r"-+\s+Step\s+([0-9]+)\s+-+"
        # multi line thermo data
        if re.match(multi_pattern, lines[0]):
            timestep_marks = [i for i, l in enumerate(lines)
                              if re.match(multi_pattern, l)]
            timesteps = np.split(lines, timestep_marks)[1:]
            dicts = []
            kv_pattern = r"([0-9A-Za-z_\[\]]+)\s+=\s+([0-9eE\.+-]+)"
            for ts in timesteps:
                data = {}
                data["Step"] = int(re.match(multi_pattern, ts[0]).group(1))
                data.update({k: float(v) for k, v
                             in re.findall(kv_pattern, "".join(ts[1:]))})
                dicts.append(data)
            df = pd.DataFrame(dicts)
            # rearrange the sequence of columns
            columns = ["Step"] + [k for k, v in
                                  re.findall(kv_pattern,
                                             "".join(timesteps[0][1:]))]
            df = df[columns]
        # one line thermo data
        else:
            df = pd.read_csv(StringIO("".join(lines)), delim_whitespace=True)
        return df

    runs = []
    for b, e in zip(begins, ends):
        runs.append(_parse_thermo(lines[b + 1:e]))
    return runs
