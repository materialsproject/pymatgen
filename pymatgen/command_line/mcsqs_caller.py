"""
Module to call mcsqs, distributed with AT-AT
https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
"""

import os
from subprocess import Popen, PIPE, TimeoutExpired
from typing import Dict, Union, List, NamedTuple

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


@requires(which("mcsqs") and which("str2cif"),
          "run_mcsqs requires first installing AT-AT, "
          "see https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/")
def run_mcsqs(
    structure: Structure,
    clusters: Dict[int, float],
    scaling: Union[int, List[int]],
    search_time: float = 0.01,
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
        search_time (int): Time spent looking for the ideal SQS in minutes
        temperature (int or float): Monte Carlo temperature (default: 1), "T" in atat code
        wr (int or float): Weight assigned to range of perfect correlation match in objective
            function (default = 1)
        wn (int or float): Multiplicative decrease in weight per additional point in cluster (default: 1)
        wd (int or float): Exponent of decay in weight as function of cluster diameter (default: 0)
        tolerance (int or float): Tolerance for matching correlations (default: 1e-3)

    Returns:
        Tuple of Pymatgen structure, which is an SQS of the input structure, and the
            mcsqs objective function
    """

    num_atoms = len(structure)

    if structure.is_ordered:
        raise ValueError("Pick a disordered structure")

    tempdir = tempfile.mkdtemp()
    os.chdir(tempdir)

    if isinstance(scaling, (int, float)):

        if scaling % 1:
            raise ValueError("Scaling should be an integer, not {}".format(scaling))

        mcsqs_find_sqs_cmd = ["mcsqs", "-n {}".format(scaling*num_atoms),
                              "-T {}".format(temperature),
                              "-wr {}".format(wr),
                              "-wn {}".format(wn),
                              "-wd {}".format(wd),
                              "-tol {}".format(tol)]

    else:

        # Set supercell to identity (will make supercell with pymatgen)
        with open("sqscell.out", "w") as f:
            f.write("1\n"
                    "1 0 0\n"
                    "0 1 0\n"
                    "0 0 1\n")

        structure = structure*scaling
        mcsqs_find_sqs_cmd = ["mcsqs", "-rc", "-n {}".format(num_atoms),
                              "-T {}".format(temperature),
                              "-wr {}".format(wr),
                              "-wn {}".format(wn),
                              "-wd {}".format(wd),
                              "-tol {}".format(tol)]

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

    # Run mcsqs to find sqs structure
    p = Popen(
        mcsqs_find_sqs_cmd
    )

    try:
        p.communicate(timeout=search_time * 60)
        raise Exception("mcsqs exited before timeout reached")
    except TimeoutExpired:
        p.kill()
        p.communicate()

        if os.path.exists("bestsqs.out") and os.path.exists("bestcorr.out"):

            # Convert output sqs structure to cif file
            p = Popen(
                "str2cif < bestsqs.out > bestsqs.cif",
                shell=True
            )
            p.communicate()

            # Get objective function
            with open('bestcorr.out', 'r') as f:
                lines = f.readlines()
            objective_function = float(lines[-1].split('=')[-1].strip())

            return Sqs(
                bestsqs=Structure.from_file("bestsqs.cif"),
                objective_function=objective_function
            )

        else:
            raise TimeoutError("Cluster expansion took too long.")

