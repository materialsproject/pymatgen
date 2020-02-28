import subprocess
import os
import numpy as np

from pymatgen.io import atat

from monty.tempfile import ScratchDir


def run_mcsqs(structure, clusters, supercell=None, total_atoms=None, search_time=0.01,
              keep_artifacts=False):
    """
    Helper function for calling mcsqs with different arguments

    Args:
        clusters (dict): dictionary of cluster interactions with entries in the form
        # atoms: cutoff in angstroms
        supercell (list): dimensions of the supercell in units of the original unit cell
        total_atoms(int): total number of atoms in the final SQS. Choose either
        this OR supercell
        search_time (int): The time spent looking for the ideal SQS in minutes
        keep_artifacts (bool): If True, retain mcsqs input and output files

    Returns:
        Pymatgen structure which is an SQS of the input structure
    """
    num_atoms = len(structure)
    total_atoms = total_atoms or num_atoms

    if supercell and total_atoms != num_atoms:
        raise ValueError("Pick supercell OR number of atoms")




    with ScratchDir('.', copy_to_current_on_exit=keep_artifacts):
        # Set supercell
        cell = np.eye(3)
        text_file = open("sqscell.out", "w")
        text_file.write("1\n")
        for i in range(len(cell)):
            text_file.write("\n")
            for j in range(len(cell[i])):
                text_file.write(str(cell[i][j]) + " ")
        text_file.close()
        struccopy = structure.copy()

        if supercell:
            struccopy.make_supercell(supercell)
            mcsqs_command = ["mcsqs", "-rc", "-n {}".format(len(structure))]
        else:
            mcsqs_command = ["mcsqs", "-n {}".format(total_atoms)]

        struc = atat.Mcsqs(struccopy)
        text_file = open("rndstr.in", "w")
        text_file.write(struc.to_string())
        text_file.close()

        # Generate Clusters
        command = ["mcsqs"]
        for num in clusters:
            command.append("-" + str(num) + "=" + str(clusters[num]))

        p = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
        )
        p.communicate()

        p = subprocess.Popen(
            mcsqs_command,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
        )
        try:
            p.communicate(timeout=search_time * 60)
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
            if os.path.exists("bestsqs.out"):
                text_file = open("bestsqs.out", "r")
                bestsqs = text_file.read()
                text_file.close()

                return atat.Mcsqs.structure_from_string(bestsqs)
            else:
                raise TimeoutError("Cluster expansion took too long.")
