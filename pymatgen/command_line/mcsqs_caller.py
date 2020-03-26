"""
Module to call mcsqs, distributed with AT-AT
https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
"""

import os
import subprocess
from typing import Dict, Union, List, NamedTuple, Optional

import numpy as np
from monty.tempfile import ScratchDir
from monty.dev import requires
from monty.os.path import which

from pymatgen import Structure


class Sqs(NamedTuple):
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
) -> Optional[Sqs]:
    """
    Helper function for calling mcsqs with different arguments
    Args:
        structure (Structure): disordered pymatgen Structure object
        clusters (dict): dictionary of cluster interactions with entries in the form
            number of atoms: cutoff in angstroms
        scaling (int or list): scaling factor to determine supercell. Two options are possible:
                a. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,
                   scaling=4 would lead to a 32 atom supercell
                b. An sequence of three scaling factors, e.g., [2, 1, 1], which
                   specifies that the supercell should have dimensions 2a x b x c
        search_time (int): The time spent looking for the ideal SQS in minutes

    Returns:
        Tuple of Pymatgen structure, which is an SQS of the input structure, and the
            mcsqs objective function
    """

    num_atoms = len(structure)

    if structure.is_ordered:
        raise ValueError("Pick a disordered structure")

    with ScratchDir("."):

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
        p = subprocess.Popen(
            mcsqs_generate_clusters_cmd
        )
        p.communicate()

        # Run mcsqs to find sqs structure
        p = subprocess.Popen(
            mcsqs_find_sqs_cmd
        )

        try:
            p.communicate(timeout=search_time * 60)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
            if os.path.exists("bestsqs.out") and os.path.exists("bestcorr.out"):

                # Convert output sqs structure to cif file
                p = subprocess.Popen(
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
