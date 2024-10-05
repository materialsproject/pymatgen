"""Module to call mcsqs, distributed with AT-AT
https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/.
"""

from __future__ import annotations

import os
import tempfile
import warnings
from pathlib import Path
from shutil import which
from subprocess import Popen, TimeoutExpired
from typing import TYPE_CHECKING, NamedTuple

from monty.dev import requires

from pymatgen.core.structure import Structure

if TYPE_CHECKING:
    from pymatgen.core.structure import IStructure


class Sqs(NamedTuple):
    """Return type for run_mcsqs."""

    bestsqs: Structure | IStructure
    objective_function: float | str
    allsqs: list
    clusters: list | str
    directory: str


@requires(
    which("mcsqs") and which("str2cif"),
    "run_mcsqs requires first installing AT-AT, see https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/",
)
def run_mcsqs(
    structure: Structure,
    clusters: dict[int, float],
    scaling: int | list[int] = 1,
    search_time: float = 60,
    directory: str | None = None,
    instances: int | None = None,
    temperature: float = 1,
    wr: float = 1,
    wn: float = 1,
    wd: float = 0.5,
    tol: float = 1e-3,
) -> Sqs:
    """Helper function for calling mcsqs with different arguments
    Args:
        structure (Structure): Disordered pymatgen Structure object
        clusters (dict): Dictionary of cluster interactions with entries in the form
            number of atoms: cutoff in angstroms
        scaling (int or list): Scaling factor to determine supercell. Two options are possible:
                a. (preferred) Scales number of atoms, e.g. for a structure with 8 atoms,
                   scaling=4 would lead to a 32 atom supercell
                b. A sequence of three scaling factors, e.g. [2, 1, 1], which
                   specifies that the supercell should have dimensions 2a x b x c
            Defaults to 1.
        search_time (float): Time spent looking for the ideal SQS in minutes (default: 60)
        directory (str): Directory to run mcsqs calculation and store files (default: None
            runs calculations in a temp directory)
        instances (int): Specifies the number of parallel instances of mcsqs to run
            (default: number of cpu cores detected by Python)
        temperature (float): Monte Carlo temperature (default: 1), "T" in atat code
        wr (float): Weight assigned to range of perfect correlation match in objective
            function (default = 1)
        wn (float): Multiplicative decrease in weight per additional point in cluster (default: 1)
        wd (float): Exponent of decay in weight as function of cluster diameter (default: 0.5)
        tol (float): Tolerance for matching correlations (default: 1e-3).

    Returns:
        tuple: Pymatgen structure SQS of the input structure, the mcsqs objective function,
            list of all SQS structures, and the directory where calculations are run
    """
    n_atoms = len(structure)

    if structure.is_ordered:
        raise ValueError("Pick a disordered structure")

    if instances is None:
        # os.cpu_count() can return None if detection fails
        instances = os.cpu_count()

    original_directory = os.getcwd()
    directory = directory or tempfile.mkdtemp()
    os.chdir(directory)

    if isinstance(scaling, int | float):
        if scaling % 1 != 0:
            raise ValueError(f"{scaling=} should be an integer")
        mcsqs_find_sqs_cmd = ["mcsqs", f"-n {scaling * n_atoms}"]

    else:
        # Set supercell to identity (will make supercell with pymatgen)
        with open("sqscell.out", mode="w") as file:
            file.write("1\n1 0 0\n0 1 0\n0 0 1\n")
        structure *= scaling
        mcsqs_find_sqs_cmd = ["mcsqs", "-rc", f"-n {n_atoms}"]

    structure.to(filename="rndstr.in")

    # Generate clusters
    mcsqs_generate_clusters_cmd = ["mcsqs"]
    for num, cutoff in clusters.items():
        mcsqs_generate_clusters_cmd.append(f"-{num}={cutoff}")

    # Run mcsqs to find clusters
    with Popen(mcsqs_generate_clusters_cmd) as process:
        process.communicate()

    # Generate SQS structures
    add_ons = [f"-T {temperature}", f"-wr {wr}", f"-wn {wn}", f"-wd {wd}", f"-tol {tol}"]

    mcsqs_find_sqs_processes = []
    if instances and instances > 1:
        # if multiple instances, run a range of commands using "-ip"
        for i in range(instances):
            instance_cmd = [f"-ip {i + 1}"]
            cmd = mcsqs_find_sqs_cmd + add_ons + instance_cmd
            process = Popen(cmd)
            mcsqs_find_sqs_processes.append(process)
    else:
        # run normal mcsqs command
        cmd = mcsqs_find_sqs_cmd + add_ons
        process = Popen(cmd)
        mcsqs_find_sqs_processes.append(process)

    try:
        for process in mcsqs_find_sqs_processes:
            process.communicate(timeout=search_time * 60)

        if instances and instances > 1:
            process = Popen(["mcsqs", "-best"])
            process.communicate()

        if os.path.isfile("bestsqs.out") and os.path.isfile("bestcorr.out"):
            return _parse_sqs_path(".")

        raise RuntimeError("mcsqs exited before timeout reached")

    except TimeoutExpired:
        for process in mcsqs_find_sqs_processes:
            process.kill()
            process.communicate()

        # Find the best sqs structures
        if instances and instances > 1:
            if not os.path.isfile("bestcorr1.out"):
                raise RuntimeError(
                    "mcsqs did not generate output files, "
                    "is search_time sufficient or are number of instances too high?"
                )

            process = Popen(["mcsqs", "-best"])
            process.communicate()

        if os.path.isfile("bestsqs.out") and os.path.isfile("bestcorr.out"):
            return _parse_sqs_path(".")

        os.chdir(original_directory)
        raise TimeoutError("Cluster expansion took too long.")


def _parse_sqs_path(path) -> Sqs:
    """Private function to parse mcsqs output directory
    Args:
        path: directory to perform parsing.

    Returns:
        tuple: Pymatgen structure SQS of the input structure, the mcsqs objective function,
            list of all SQS structures, and the directory where calculations are run
    """
    path = Path(path)

    # detected instances will be 0 if mcsqs was run in series, or number of instances
    detected_instances = len(list(path.glob("bestsqs*[0-9]*.out")))

    # Convert best SQS structure to CIF file and pymatgen Structure
    with (
        open(os.path.join(path, "bestsqs.out")) as input_file,
        open(os.path.join(path, "bestsqs.cif"), "w") as output_file,
    ):
        process = Popen(["str2cif"], stdin=input_file, stdout=output_file, cwd=path)
        process.communicate()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        best_sqs = Structure.from_file(path / "bestsqs.out")

    # Get best SQS objective function
    with open(path / "bestcorr.out") as file:
        lines = file.readlines()

    objective_function_str = lines[-1].split("=")[-1].strip()
    objective_function: float | str
    objective_function = float(objective_function_str) if objective_function_str != "Perfect_match" else "Perfect_match"

    # Get all SQS structures and objective functions
    all_sqs = []

    for idx in range(detected_instances):
        sqs_out = os.path.join(path, f"bestsqs{idx + 1}.out")
        sqs_cif = os.path.join(path, f"bestsqs{idx + 1}.cif")

        with open(sqs_out) as input_file, open(sqs_cif, "w") as output_file:
            process = Popen(["str2cif"], stdin=input_file, stdout=output_file, cwd=path)
            process.communicate()
        sqs = Structure.from_file(path / sqs_out)

        corr_out = f"bestcorr{idx + 1}.out"
        with open(path / corr_out) as file:
            lines = file.readlines()

        objective_function_str = lines[-1].split("=")[-1].strip()
        obj: float | str
        obj = float(objective_function_str) if objective_function_str != "Perfect_match" else "Perfect_match"
        all_sqs.append({"structure": sqs, "objective_function": obj})

    clusters = _parse_clusters(path / "clusters.out")

    return Sqs(
        bestsqs=best_sqs,
        objective_function=objective_function,
        allsqs=all_sqs,
        directory=str(path.resolve()),
        clusters=clusters,
    )


def _parse_clusters(filename):
    """Private function to parse clusters.out file
    Args:
        path: directory to perform parsing.

    Returns:
        list[dict]: List of cluster dictionaries with keys:
            multiplicity: int
            longest_pair_length: float
            num_points_in_cluster: int
            coordinates: list[dict] of points with keys:
                coordinates: list[float]
                num_possible_species: int
                cluster_function: float
    """
    with open(filename) as file:
        lines = file.readlines()

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
            point_dict["coordinates"] = [float(line) for line in line[:3]]
            point_dict["num_possible_species"] = int(line[3]) + 2  # see ATAT manual for why +2
            point_dict["cluster_function"] = float(line[4])  # see ATAT manual for what "function" is
            points.append(point_dict)

        cluster_dict["coordinates"] = points
        cluster_dicts.append(cluster_dict)

    return cluster_dicts
