"""
Module to call mcsqs, distributed with AT-AT
https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
"""

import os
import tempfile
import warnings
from pathlib import Path
from subprocess import Popen, TimeoutExpired
from collections import namedtuple
from typing import Dict, List, Optional, Union

from monty.dev import requires
from monty.os.path import which

from pymatgen.core.structure import Structure


Sqs = namedtuple("Sqs", "bestsqs objective_function allsqs clusters directory")
"""
Return type for run_mcsqs.
bestsqs: Structure
objective_function: Union[float, str]
allsqs: List
clusters: List
directory: str
"""


@requires(
    which("mcsqs") and which("str2cif"),
    "run_mcsqs requires first installing AT-AT, " "see https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/",
)
def run_mcsqs(
    structure: Structure,
    clusters: Dict[int, float],
    scaling: Union[int, List[int]] = 1,
    search_time: float = 60,
    directory: Optional[str] = None,
    instances: Optional[int] = None,
    temperature: Union[int, float] = 1,
    wr: float = 1,
    wn: float = 1,
    wd: float = 0.5,
    tol: float = 1e-3,
) -> Sqs:
    """
    Helper function for calling mcsqs with different arguments
    Args:
        structure (Structure): Disordered pymatgen Structure object
        clusters (dict): Dictionary of cluster interactions with entries in the form
            number of atoms: cutoff in angstroms
        scaling (int or list): Scaling factor to determine supercell. Two options are possible:
                a. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,
                   scaling=4 would lead to a 32 atom supercell
                b. A sequence of three scaling factors, e.g., [2, 1, 1], which
                   specifies that the supercell should have dimensions 2a x b x c
            Defaults to 1.
        search_time (float): Time spent looking for the ideal SQS in minutes (default: 60)
        directory (str): Directory to run mcsqs calculation and store files (default: None
            runs calculations in a temp directory)
        instances (int): Specifies the number of parallel instances of mcsqs to run
            (default: number of cpu cores detected by Python)
        temperature (int or float): Monte Carlo temperature (default: 1), "T" in atat code
        wr (int or float): Weight assigned to range of perfect correlation match in objective
            function (default = 1)
        wn (int or float): Multiplicative decrease in weight per additional point in cluster (default: 1)
        wd (int or float): Exponent of decay in weight as function of cluster diameter (default: 0.5)
        tol (int or float): Tolerance for matching correlations (default: 1e-3)

    Returns:
        Tuple of Pymatgen structure SQS of the input structure, the mcsqs objective function,
            list of all SQS structures, and the directory where calculations are run
    """

    num_atoms = len(structure)

    if structure.is_ordered:
        raise ValueError("Pick a disordered structure")

    if instances is None:
        # os.cpu_count() can return None if detection fails
        instances = os.cpu_count()

    original_directory = os.getcwd()
    if not directory:
        directory = tempfile.mkdtemp()
    os.chdir(directory)

    if isinstance(scaling, (int, float)):

        if scaling % 1:
            raise ValueError("Scaling should be an integer, not {}".format(scaling))
        mcsqs_find_sqs_cmd = ["mcsqs", "-n {}".format(scaling * num_atoms)]

    else:

        # Set supercell to identity (will make supercell with pymatgen)
        with open("sqscell.out", "w") as f:
            f.write("1\n1 0 0\n0 1 0\n0 0 1\n")
        structure = structure * scaling
        mcsqs_find_sqs_cmd = ["mcsqs", "-rc", "-n {}".format(num_atoms)]

    structure.to(filename="rndstr.in")

    # Generate clusters
    mcsqs_generate_clusters_cmd = ["mcsqs"]
    for num in clusters:
        mcsqs_generate_clusters_cmd.append("-" + str(num) + "=" + str(clusters[num]))

    # Run mcsqs to find clusters
    p = Popen(mcsqs_generate_clusters_cmd)
    p.communicate()

    # Generate SQS structures
    add_ons = [
        "-T {}".format(temperature),
        "-wr {}".format(wr),
        "-wn {}".format(wn),
        "-wd {}".format(wd),
        "-tol {}".format(tol),
    ]

    mcsqs_find_sqs_processes = []
    if instances and instances > 1:
        # if multiple instances, run a range of commands using "-ip"
        for i in range(instances):
            instance_cmd = ["-ip {}".format(i + 1)]
            cmd = mcsqs_find_sqs_cmd + add_ons + instance_cmd
            p = Popen(cmd)
            mcsqs_find_sqs_processes.append(p)
    else:
        # run normal mcsqs command
        cmd = mcsqs_find_sqs_cmd + add_ons
        p = Popen(cmd)
        mcsqs_find_sqs_processes.append(p)

    try:
        for idx, p in enumerate(mcsqs_find_sqs_processes):
            p.communicate(timeout=search_time * 60)

        if instances and instances > 1:
            p = Popen(["mcsqs", "-best"])
            p.communicate()

        if os.path.exists("bestsqs.out") and os.path.exists("bestcorr.out"):
            return _parse_sqs_path(".")

        raise RuntimeError("mcsqs exited before timeout reached")

    except TimeoutExpired:
        for p in mcsqs_find_sqs_processes:
            p.kill()
            p.communicate()

        # Find the best sqs structures
        if instances and instances > 1:

            if not os.path.exists("bestcorr1.out"):
                raise RuntimeError(
                    "mcsqs did not generate output files, "
                    "is search_time sufficient or are number of instances too high?"
                )

            p = Popen(["mcsqs", "-best"])
            p.communicate()

        if os.path.exists("bestsqs.out") and os.path.exists("bestcorr.out"):
            sqs = _parse_sqs_path(".")
            return sqs

        os.chdir(original_directory)
        raise TimeoutError("Cluster expansion took too long.")


def _parse_sqs_path(path) -> Sqs:
    """
    Private function to parse mcsqs output directory
    Args:
        path: directory to perform parsing

    Returns:
        Tuple of Pymatgen structure SQS of the input structure, the mcsqs objective function,
            list of all SQS structures, and the directory where calculations are run
    """

    path = Path(path)

    # detected instances will be 0 if mcsqs was run in series, or number of instances
    detected_instances = len(list(path.glob("bestsqs*[0-9]*.out")))

    # Convert best SQS structure to cif file and pymatgen Structure
    p = Popen("str2cif < bestsqs.out > bestsqs.cif", shell=True, cwd=path)
    p.communicate()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bestsqs = Structure.from_file(path / "bestsqs.out")

    # Get best SQS objective function
    with open(path / "bestcorr.out", "r") as f:
        lines = f.readlines()

    objective_function_str = lines[-1].split("=")[-1].strip()
    objective_function: Union[float, str]
    if objective_function_str != "Perfect_match":
        objective_function = float(objective_function_str)
    else:
        objective_function = "Perfect_match"

    # Get all SQS structures and objective functions
    allsqs = []

    for i in range(detected_instances):
        sqs_out = "bestsqs{}.out".format(i + 1)
        sqs_cif = "bestsqs{}.cif".format(i + 1)
        corr_out = "bestcorr{}.out".format(i + 1)
        p = Popen("str2cif < {} > {}".format(sqs_out, sqs_cif), shell=True, cwd=path)
        p.communicate()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sqs = Structure.from_file(path / sqs_out)
        with open(path / corr_out, "r") as f:
            lines = f.readlines()

        objective_function_str = lines[-1].split("=")[-1].strip()
        obj: Union[float, str]
        if objective_function_str != "Perfect_match":
            obj = float(objective_function_str)
        else:
            obj = "Perfect_match"
        allsqs.append({"structure": sqs, "objective_function": obj})

    clusters = _parse_clusters(path / "clusters.out")

    return Sqs(
        bestsqs=bestsqs,
        objective_function=objective_function,
        allsqs=allsqs,
        directory=str(path.resolve()),
        clusters=clusters,
    )


def _parse_clusters(filename):
    """
    Private function to parse clusters.out file
    Args:
        path: directory to perform parsing

    Returns:
        List of dicts
    """

    with open(filename, "r") as f:
        lines = f.readlines()

    clusters = []
    cluster_block = []
    for line in lines:
        line = line.split("\n")[0]
        if line == "":
            clusters.append(cluster_block)
            cluster_block = []
        else:
            cluster_block.append(line)

    cluster_dicts = []
    for cluster in clusters:
        cluster_dict = {
            "multiplicity": int(cluster[0]),
            "longest_pair_length": float(cluster[1]),
            "num_points_in_cluster": int(cluster[2]),
        }
        points = []
        for point in range(cluster_dict["num_points_in_cluster"]):
            line = cluster[3 + point].split(" ")
            point_dict = {}
            point_dict["coordinates"] = [float(line) for line in line[0:3]]
            point_dict["num_possible_species"] = int(line[3]) + 2  # see ATAT manual for why +2
            point_dict["cluster_function"] = float(line[4])  # see ATAT manual for what "function" is
            points.append(point_dict)

        cluster_dict["coordinates"] = points
        cluster_dicts.append(cluster_dict)

    return cluster_dicts
