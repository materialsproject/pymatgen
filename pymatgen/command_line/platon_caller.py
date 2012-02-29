'''
Interface with command line platon
http://www.cryst.chem.uu.nl/spek/xraysoft/
Only tested on Linux
inspired by Shyue's qhull_caller
WARNING: you need to have a platon in your path for this to work
'''

__author__="Geoffroy Hautier"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy.hautier@uclouvain.be"
__status__ = "Production"
__date__ ="$February 24, 2012M$"

import subprocess
import pymatgen.io.cifio
import os


def run_platon_command(command, structure):
    """
    Helper function for calling platon with different arguments
    I know it's uggly to actually write a file and call platon but I did not manage to make it work in another way
    """
    writer=pymatgen.io.cifio.CifWriter(structure)
    #incommand=str(writer)
    #print incommand
    writer.write_file("/tmp/tmp.cif")
    command.append("/tmp/tmp.cif")
    p = subprocess.Popen(command,shell = False, stdout = subprocess.PIPE, stdin = subprocess.PIPE)
    output = p.communicate()
    os.remove("/tmp/tmp.cif")
    return output

def get_space_group(structure):
    output=run_platon_command(['platon', '-o', '-c'], structure)
    dictio={}
    for line in output[0].split("\n"):
        #   print line
        if(line.find("Space Group")!=-1):
            list_tmp=line.split()
            #       print list_tmp
            for i in range(len(list_tmp)):
                if(list_tmp[i]=='Group'):
                    dictio['SG_HM']=list_tmp[i+1]
                if(list_tmp[i]=='No:'):
                    dictio['SG_NB']=list_tmp[i+1]
    return dictio
