# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements classes and methods for processing LAMMPS output
files (log and dump).
"""

import glob
import re
from io import StringIO

import numpy as np
import pandas as pd
from monty.io import zopen
from monty.json import MSONable
from scipy.linalg import norm

from pymatgen.io.lammps.data import LammpsBox
from pymatgen.util.coord import pbc_diff

__author__ = "Kiran Mathew, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"


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
        pattern = file_pattern.replace("*", "([0-9]+)").replace("\\", "\\\\")
        files = sorted(files, key=lambda f: int(re.match(pattern, f).group(1)))

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
    with zopen(filename, "rt") as f:
        lines = f.readlines()
    begin_flag = (
        "Memory usage per processor =",
        "Per MPI rank memory allocation (min/avg/max) =",
    )
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
            timestep_marks = [i for i, l in enumerate(lines) if re.match(multi_pattern, l)]
            timesteps = np.split(lines, timestep_marks)[1:]
            dicts = []
            kv_pattern = r"([0-9A-Za-z_\[\]]+)\s+=\s+([0-9eE\.+-]+)"
            for ts in timesteps:
                data = {}
                data["Step"] = int(re.match(multi_pattern, ts[0]).group(1))
                data.update({k: float(v) for k, v in re.findall(kv_pattern, "".join(ts[1:]))})
                dicts.append(data)
            df = pd.DataFrame(dicts)
            # rearrange the sequence of columns
            columns = ["Step"] + [k for k, v in re.findall(kv_pattern, "".join(timesteps[0][1:]))]
            df = df[columns]
        # one line thermo data
        else:
            df = pd.read_csv(StringIO("".join(lines)), delim_whitespace=True)
        return df

    runs = []
    for b, e in zip(begins, ends):
        runs.append(_parse_thermo(lines[b + 1 : e]))
    return runs


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
        if "id" in data.columns:
            data.set_index("id", inplace=True)
            data.sort_index(inplace=True)
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
        data = pd.read_csv(StringIO("\n".join(lines[9:])), names=data_head, delim_whitespace=True)
        return cls(timestep, natoms, box, data)

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            LammpsDump
        """
        items = {
            "timestep": d["timestep"],
            "natoms": d["natoms"],
            "box": LammpsBox.from_dict(d["box"]),
            "data": pd.read_json(d["data"], orient="split"),
        }
        items["data"].index.name = "id"
        return cls(**items)

    def as_dict(self):
        """
        Returns: MSONable dict
        """
        d = {}
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        d["timestep"] = self.timestep
        d["natoms"] = self.natoms
        d["box"] = self.box.as_dict()
        d["data"] = self.data.to_json(orient="split")
        return d


class LammpsTrajectory(MSONable):
    """
    Object for representing a trajectory of atoms using dump files obtained with Lammps.
    The timestep is considered to be constant for the time being.

    Args:
        trajectory: dictionary with timesteps as keys and LammpsDump objects as values
    """

    def __init__(self, trajectory):
        self.trajectory = trajectory

        # Getting the timestep (assumed to be constant)
        timestep = np.diff(list(trajectory.keys()))[0]

        self.timestep = timestep

        if self.timestep <= 0:
            raise ValueError("The timestep cannot be <= 0.")

    @classmethod
    def from_file_pattern(cls, file_pattern, timestep_min=None, timestep_max=None):
        """
        Creates a LammpsTrajectory object from a file pattern for the dump files.
        The created trajectory object contains timesteps between timestep_min and
        timestep_max.

        Args:
            file_pattern (str): Filename to parse. The timestep wildcard
                (e.g., dump.atom.'*') is supported and the files are parsed
                in the sequence of timestep.
            timestep_min: minimal timestep to include in the trajectory
            timestep_max: maximal timestep to include in the trajectory
        """
        if timestep_min is None:
            timestep_min = 0
        if timestep_max is None:
            timestep_max = 1e20

        trajectory = {}
        for dump in parse_lammps_dumps(file_pattern):
            if timestep_min <= dump.timestep <= timestep_max:
                trajectory[dump.timestep] = dump

        return cls(trajectory)

    @classmethod
    def from_dict(cls, d):
        trajectory = {}
        for timestep, dump in d.items():
            if not (isinstance(timestep, int) or isinstance(timestep, float)):
                raise ValueError("The provided keys should be timesteps.")
            if not isinstance(dump, LammpsDump):
                raise ValueError("The provided values should be LammpsDump objects.")

            trajectory[timestep] = dump

        return LammpsTrajectory(trajectory)

    def get_msd(self, atom_type=None, time_average=True):
        """
        Compute the mean-square displacement for a trajectory.
        MSD = (1/N) sum_i^N < | R_i(t' + t) - R_i(t') |^2 >
        The <> denote the time average over possible starting times t'
        allowed by the elapsed time t.
        This computation assumes that the lattice does not change with the timestep.

        Args:
            atom_type: int or array of ints containing the types of atoms
                       for which the MSD will be computed
            time_average: True if an average is to be performed on the possible
                          starting times t' (better statistical average).

        Returns:
            Array containing the MSD for each timestep in the trajectory.
        """
        # First and last steps of the trajectory
        t0 = np.min(list(self.trajectory.keys()))
        tmax = np.max(list(self.trajectory.keys()))

        if atom_type is None:
            # Set the default atom_type to all if not specified
            atom_type = np.unique(self.trajectory[t0].data["type"])
        atom_type = np.atleast_1d(atom_type)

        # Initialize MSD [n_atom_type, n_timesteps, 4]
        # The last 4 components are the MSD along x, y, z, and total
        MSD = np.zeros([len(atom_type), len(self.trajectory.keys()), 4])

        # Loop over timesteps and atom types to get their MSD(t)
        # The average over initial times is done if time_average is True

        # The MSD for the first time-step is zero
        trajectory_after_t0 = self.trajectory.copy()
        # We discard it for the loop
        trajectory_after_t0.pop(t0)
        for time in trajectory_after_t0.keys():
            # Get the index in the array of the timesteps
            # instead of its value
            time_index = int(np.round(time / self.timestep))

            # Loop over atom types
            for ityp, typ in enumerate(atom_type):
                # Now we perform the time average
                if time_average:
                    # Get the possible starting times t' for a given time spent t and maximum tmax
                    tprime = [key for key in self.trajectory.keys() if key <= tmax - time]
                else:
                    tprime = [0]

                for ipstep in tprime:
                    # X(t')
                    Xtp = self.trajectory[ipstep].data
                    Xtp = Xtp.loc[Xtp["type"] == typ]
                    Xtp = Xtp[["xs", "ys", "zs"]].values
                    Xtp = np.atleast_2d(Xtp)

                    # X(t' + t)
                    Xtpt = self.trajectory[ipstep + time].data
                    Xtpt = Xtpt.loc[Xtpt["type"] == typ]
                    Xtpt = Xtpt[["xs", "ys", "zs"]].values
                    Xtpt = np.atleast_2d(Xtpt)

                    # Assumes that the lattice does not change (it should not)
                    diff = pbc_diff(Xtp, Xtpt)
                    diff = self.trajectory[t0].box.to_lattice().get_cartesian_coords(diff)
                    MSD[ityp, time_index, 0:3] += np.average(diff**2, axis=0)
                    MSD[ityp, time_index, 3] += np.average(norm(diff, axis=1) ** 2)

                MSD[ityp, time_index] = MSD[ityp, time_index] / len(tprime)

        return MSD

    def get_diffusion_coefficient(self, MSD=None, start=None, stop=None):
        """
        Compute the diffusion coefficient from the MSD

        Args:
            MSD: array containing (one or multiple components of) the mean-square displacement.
            start: index along MSD determining the start of the diffusion computation.
                   This can be used to select the linear regime of the MSD.
                   By default, removes 10% of the data at the start of the MSD.
            stop: index along MSD determining the stop of the diffusion computation.
                  By default, removes 10% of the data at the end of the MSD.

        Returns:
            Tuple (D, time) with D the diffusion coefficient for the atoms for which the MSD is provided,
            and timesteps given in time. The time dimension of D and time is 1 less than that of
            the MSD. The last dimension of D corresponds to [Dx, Dy, Dz, Dtot] based on the same
            corresponding MSD.
        """
        # First we get the MSD if not present
        if MSD is None:
            MSD = self.get_msd()
        # Then we determine the start and end indexes of the time
        # dimension of the MSD, that should correspond to the linear part.
        # If they are not specified, we remove the first and final 10% of the data.
        if start is None:
            start = int(np.ceil(len(MSD[0, :, 0]) * 0.1))
        if stop is None:
            stop = int(np.floor(len(MSD[0, :, 0]) * 0.9))

        MSD = np.atleast_3d(MSD[:, start:stop, :])
        time = np.linspace(
            (start + 1) * self.timestep,  # + 1 because D starts one timestep after MSD
            (start + len(MSD[0, :, 0]) - 1) * self.timestep,
            len(MSD[0, :, 0]) - 1,
        )
        time_div = np.atleast_3d(time)
        D = np.zeros([MSD.shape[0], MSD.shape[1] - 1, MSD.shape[2]])
        for i in range(MSD.shape[0]):
            for j in range(MSD.shape[1] - 1):
                D[i] += MSD[i, j + 1, :] - MSD[i, 0, :]
        D = D / (6 * (time_div - time_div[0, 0, 0] + self.timestep))

        return D, time
