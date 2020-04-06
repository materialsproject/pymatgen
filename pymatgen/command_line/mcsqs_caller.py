"""
Module to call mcsqs, distributed with AT-AT
https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
"""

import os
from subprocess import Popen, PIPE, TimeoutExpired
from typing import Dict, Union, List, NamedTuple, Optional, AnyStr

from monty.tempfile import ScratchDir
from monty.dev import requires
from monty.os.path import which
import tempfile

from pymatgen import Structure


class Sqs(NamedTuple):
    """
    Return type for run_mcsqs.
    """
    bestsqs: Structure
    objective_function: float
    allsqs: List
    directory: AnyStr

@requires(which("mcsqs") and which("str2cif"),
          "run_mcsqs requires first installing AT-AT, "
          "see https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/")
def run_mcsqs(
    structure: Structure,
    clusters: Dict[int, float],
    scaling: Union[int, List[int]],
    search_time: float = 0.01,
    directory: Optional[str] = None,
    instances: Optional[int] = None,
    temperature: Union[int, float] = 1,
    wr: float = 1,
    wn: float = 1,
    wd: float = 0,
    tol: float = 1e-3
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
    Keyword Args:
        search_time (int): Time spent looking for the ideal SQS in minutes (default: 0.1)
        directory (str): Directory to reun mcsqs calculation and store files (default: None
            runs calculations in a temp directory)
        instances (int): Specifies the number of parallel instances of mcsqs to run
            (default: number of cpu cores detected by Python)
        temperature (int or float): Monte Carlo temperature (default: 1), "T" in atat code
        wr (int or float): Weight assigned to range of perfect correlation match in objective
            function (default = 1)
        wn (int or float): Multiplicative decrease in weight per additional point in cluster (default: 1)
        wd (int or float): Exponent of decay in weight as function of cluster diameter (default: 0)
        tol (int or float): Tolerance for matching correlations (default: 1e-3)

    Returns:
        Tuple of Pymatgen structure SQS of the input structure, the mcsqs objective function,
            list of all SQS structures, and the directory where calculations are run
    """

    num_atoms = len(structure)

    if structure.is_ordered:
        raise ValueError("Pick a disordered structure")

    if instances is None:
        instances = os.cpu_count()

    # TODO: figure out how to handle instances=1 (doesn't generate bestsqs)

    original_directory = os.getcwd()
    if not directory:
        directory = tempfile.mkdtemp()
    os.chdir(directory)

    if isinstance(scaling, (int, float)):

        if scaling % 1:
            raise ValueError("Scaling should be an integer, not {}".format(scaling))
        mcsqs_find_sqs_cmd = ["mcsqs", "-n {}".format(scaling*num_atoms)]

    else:

        # Set supercell to identity (will make supercell with pymatgen)
        with open("sqscell.out", "w") as f:
            f.write("1\n"
                    "1 0 0\n"
                    "0 1 0\n"
                    "0 0 1\n")
        structure = structure*scaling
        mcsqs_find_sqs_cmd = ["mcsqs", "-rc", "-n {}".format(num_atoms)]

    structure.to(filename="rndstr.in")

    # Generate clusters
    mcsqs_generate_clusters_cmd = ["mcsqs"]
    for num in clusters:
        mcsqs_generate_clusters_cmd.append("-" + str(num) + "=" + str(clusters[num]))

    # Run mcsqs to find clusters
    p = Popen(
        mcsqs_generate_clusters_cmd
    )
    p.communicate()

    # Generate SQS structures
    add_ons = [
                  "-T {}".format(temperature),
                  "-wr {}".format(wr),
                  "-wn {}".format(wn),
                  "-wd {}".format(wd),
                  "-tol {}".format(tol)
              ]

    mcsqs_find_sqs_processes = []
    for i in range(instances):
        instance_cmd = ["-ip {}".format(i + 1)]
        cmd = mcsqs_find_sqs_cmd + add_ons + instance_cmd
        p = Popen(cmd)
        mcsqs_find_sqs_processes.append(p)

    try:
        for idx, p in enumerate(mcsqs_find_sqs_processes):
            # print("Starting comm with {}".format(idx))
            p.communicate(timeout=search_time * 60)
            # print("Finish comm with {}".format(idx))
            os.chdir(original_directory)
            raise Exception("mcsqs exited before timeout reached")
    except TimeoutExpired:
        for p in mcsqs_find_sqs_processes:
            p.kill()
            p.communicate()

        # Find the best sqs structures
        if instances > 1:
            p = Popen(
                ["mcsqs", "-best"]
            )
            p.communicate()

        if os.path.exists("bestsqs.out") and os.path.exists("bestcorr.out"):

            # Convert best SQS structure to cif file and pymatgen Structure
            p = Popen(
                "str2cif < bestsqs.out > bestsqs.cif",
                shell=True
            )
            p.communicate()
            bestsqs = Structure.from_file("bestsqs.cif")

            # Get best SQS objective function
            with open('bestcorr.out', 'r') as f:
                lines = f.readlines()
            objective_function = float(lines[-1].split('=')[-1].strip())

            # Get all SQS structures and objective functions
            allsqs = []
            for i in range(instances):
                sqs_out = 'bestsqs' + str(i+1) + '.out'
                sqs_cif = 'bestsqs' + str(i+1) + '.cif'
                corr_out = 'bestcorr' + str(i+1) + '.out'
                p = Popen(
                    "str2cif <" + sqs_out + ">" + sqs_cif,
                    shell=True
                )
                p.communicate()
                sqs = Structure.from_file(sqs_cif)
                with open(corr_out, 'r') as f:
                    lines = f.readlines()
                obj = float(lines[-1].split('=')[-1].strip())
                allsqs.append({'structure': sqs, 'objective_function': obj})

            # Change back to original directory
            os.chdir(original_directory)

            return Sqs(
                bestsqs=bestsqs,
                objective_function=objective_function,
                allsqs=allsqs,
                directory=directory
            )

        else:
            os.chdir(original_directory)
            raise TimeoutError("Cluster expansion took too long.")

