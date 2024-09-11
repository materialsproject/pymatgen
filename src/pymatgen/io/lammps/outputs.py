"""
This module implements classes and methods for processing LAMMPS output
files (log and dump).
"""

from __future__ import annotations

import re
from glob import glob
from io import StringIO
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.lammps.data import LammpsBox

if TYPE_CHECKING:
    from typing import Any

    from typing_extensions import Self

__author__ = "Kiran Mathew, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"


class LammpsDump(MSONable):
    """Object for representing dump data for a single snapshot."""

    def __init__(self, timestep: int, natoms: int, box: LammpsBox, data: pd.DataFrame) -> None:
        """
        Base constructor.

        Args:
            timestep (int): Current time step.
            natoms (int): Total number of atoms in the box.
            box (LammpsBox): Simulation box.
            data (pd.DataFrame): Dumped atomic data.
        """
        self.timestep = timestep
        self.natoms = natoms
        self.box = box
        self.data = data

    @classmethod
    def from_str(cls, string: str) -> Self:
        """
        Constructor from string parsing.

        Args:
            string (str): Input string.
        """
        lines = string.split("\n")
        time_step = int(lines[1])
        n_atoms = int(lines[3])
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
        data = pd.read_csv(StringIO("\n".join(lines[9:])), names=data_head, sep=r"\s+")
        return cls(time_step, n_atoms, box, data)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            LammpsDump
        """
        items = {"timestep": dct["timestep"], "natoms": dct["natoms"]}
        items["box"] = LammpsBox.from_dict(dct["box"])
        items["data"] = pd.read_json(dct["data"], orient="split")
        return cls(**items)

    def as_dict(self) -> dict[str, Any]:
        """Get MSONable dict."""
        dct: dict[str, Any] = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["timestep"] = self.timestep
        dct["natoms"] = self.natoms
        dct["box"] = self.box.as_dict()
        dct["data"] = self.data.to_json(orient="split")
        return dct


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
    files = glob(file_pattern)
    if len(files) > 1:
        pattern = file_pattern.replace("*", "([0-9]+)").replace("\\", "\\\\")
        files = sorted(files, key=lambda f: int(re.match(pattern, f)[1]))

    for filename in files:
        with zopen(filename, mode="rt") as file:
            dump_cache = []
            for line in file:
                if line.startswith("ITEM: TIMESTEP"):
                    if len(dump_cache) > 0:
                        yield LammpsDump.from_str("".join(dump_cache))
                    dump_cache = [line]
                else:
                    dump_cache.append(line)
            yield LammpsDump.from_str("".join(dump_cache))


def parse_lammps_log(filename: str = "log.lammps") -> list[pd.DataFrame]:
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
    with zopen(filename, mode="rt") as file:
        lines = file.readlines()
    begin_flag = (
        "Memory usage per processor =",
        "Per MPI rank memory allocation (min/avg/max) =",
    )
    end_flag = "Loop time of"
    begins, ends = [], []
    for idx, line in enumerate(lines):
        if line.startswith(begin_flag):
            begins.append(idx)
        elif line.startswith(end_flag):
            ends.append(idx)

    def _parse_thermo(lines: list[str]) -> pd.DataFrame:
        multi_pattern = r"-+\s+Step\s+([0-9]+)\s+-+"
        # multi line thermo data
        if re.match(multi_pattern, lines[0]):
            timestep_marks = [idx for idx, line in enumerate(lines) if re.match(multi_pattern, line)]
            time_steps = np.split(lines, timestep_marks)[1:]
            dicts = []
            kv_pattern = r"([0-9A-Za-z_\[\]]+)\s+=\s+([0-9eE\.+-]+)"
            for ts in time_steps:
                data = {}
                step = re.match(multi_pattern, ts[0])
                if step is None:
                    raise ValueError("step is None")
                data["Step"] = int(step[1])
                data |= {k: float(v) for k, v in re.findall(kv_pattern, "".join(ts[1:]))}
                dicts.append(data)
            df_thermo = pd.DataFrame(dicts)
            # rearrange the sequence of columns
            columns = ["Step"] + [k for k, v in re.findall(kv_pattern, "".join(time_steps[0][1:]))]
            df_thermo = df_thermo[columns]
        # one line thermo data
        else:
            df_thermo = pd.read_csv(StringIO("".join(lines)), sep=r"\s+")
        return df_thermo

    runs = []
    for b, e in zip(begins, ends, strict=True):
        runs.append(_parse_thermo(lines[b + 1 : e]))
    return runs
