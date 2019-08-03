import subprocess
import os 
import numpy as np
import time
from threading import Timer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io import atat

def run_mcsqs(structure, clusters, supercell = None, total_atoms = None, max_time = None):
    """
    Helper function for calling mcsqs with different arguments
    max_time (int): max time for running mcsqs in minutes
    """
    num_atoms = len(structure)
    if total_atoms == None:
        total_atoms = num_atoms

    if supercell is not None and total_atoms != num_atoms:
        print("pick supercell OR number of atoms")
        return

    ## Set supercell

    cell = np.eye(3)
    text_file = open("sqscell.out", "w")
    text_file.write('1\n')
    for i in range(len(cell)):
        text_file.write('\n')
        for j in range(len(cell[i])):
            text_file.write(str(cell[i][j]) +' ')
    text_file.close()
    struccopy = structure.copy()

    if supercell is not None:
        struccopy.make_supercell(supercell)
        struc = atat.Mcsqs(struccopy)
        text_file = open("rndstr.in", "w")
        text_file.write(struc.to_string())
        text_file.close()

        
        ## Generate Clusters
        command = ['mcsqs']
        for num in clusters:
            command.append('-'+ str(num) + '=' + str(clusters[num]))

        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        ## Create stop file
        text_file = open("stopsqs", "w")
        text_file.close()
        ## Run mcsqs
        time.sleep(0.01)
        command = ['mcsqs',  '-rc', "-n {}".format(len(structure))]
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        if max_time:
            timed_out = False
            timer = Timer(max_time*60, lambda p: p.kill(), [p])

            try:
                timer.start()
                output = p.communicate()[0].decode("utf-8")
            finally:
                if not timer.is_alive():
                    timed_out = True
                timer.cancel()

            if timed_out:
                raise TimeoutError('Cluster expansion took too long.')
    else:
        struc = atat.Mcsqs(struccopy)
        text_file = open("rndstr.in", "w")
        text_file.write(struc.to_string())
        text_file.close()

        
        ## Generate Clusters
        command = ['mcsqs']
        for num in clusters:
            command.append('-'+ str(num) + '=' + str(clusters[num]))

        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        ## Create stop file
        text_file = open("stopsqs", "w")
        text_file.close()
        ## Run mcsqs

        time.sleep(0.01)
        command = ['mcsqs', "-n {}".format(total_atoms)]
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        if max_time:
            timed_out = False
            timer = Timer(max_time*60, lambda p: p.kill(), [p])

            try:
                timer.start()
                output = p.communicate()[0].decode("utf-8")
            finally:
                if not timer.is_alive():
                    timed_out = True
                timer.cancel()

            if timed_out:
                raise TimeoutError('Cluster expansion took too long.')
        
    time.sleep(0.01)
    text_file = open("bestsqs.out", "r")
    bestsqs = text_file.read()
    text_file.close()

    return atat.Mcsqs.structure_from_string(bestsqs)
