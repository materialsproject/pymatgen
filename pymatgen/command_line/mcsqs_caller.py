import subprocess
import os 
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io import atat

def run_mcsqs(structure, clusters, supercell = None, total_atoms = None, max_time = None):
    """
    Helper function for calling aconvasp with different arguments
    """
    num_atoms = len(structure)
    if total_atoms == None:
        total_atoms = num_atoms

    if supercell != None and total_atoms != num_atoms:
        print("pick supercell or number of atoms")
        return

    ## Set supercell
    if supercell == None:
        supercell = np.eye(3)
    else:
        supercell = np.array(supercell)
        text_file = open("sqscell.out", "w")
        text_file.write('1\n')
        for i in range(len(supercell)):
            text_file.write('\n')
            for j in range(len(supercell[i])):
                text_file.write(str(supercell[i][j]) +' ')
        text_file.close()
        
    struc = atat.Mcsqs(structure)
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
    command = ['mcsqs -n ' + str(np.prod(supercell)*total_atoms)]
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    os.remove('sqscell.out')
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
            raise TimeoutError('Enumeration took too long.')
            
    text_file = open("bestsqs.out", "r")
    bestsqs = text_file.read()

    return atat.Mcsqs.structure_from_string(bestsqs)
