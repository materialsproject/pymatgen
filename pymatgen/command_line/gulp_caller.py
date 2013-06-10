#!/usr/bin/env python

"""
Interface with command line GULP. http://projects.ivec.org/gulp/help/manuals.html
WARNING: you need to have GULP in your path for this to work
"""

__author__ = "Wenhao Sun"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Wenhao Sun"
__email__ = "wenhao@mit.edu"
__status__ = "Production"
__date__ = "$Jan 22, 2013M$"

import subprocess
import os
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.command_line.aconvasp_caller import run_aconvasp_command
from pymatgen.core.periodic_table import Element


_anion_set = {
              # Chalcogen group
              Element("O"), Element("S"), 
              # Halogen group
              Element("F"), Element("Cl"), Element("Br"),
              # Pnictogen group
              Element("N"), Element("P")
              }
_cation_set = {
               # alkali metals
               Element("Li"), Element("Na"), Element("K"), 
               # alkaline metals
               Element("Be"), Element("Mg"), Element("Ca"),                
               # other most probably used elements 
               Element("Fe"), Element("Co"), Element("Ni"),
               Element("Cr"), Element("Mn"), Element("Zn"), 
               Element("Ti"), Element("Al"), Element("Zr"),
               Element("Cu"), Element("Ag"), Element("Au"),
               Element("Mn")
               }

def _gulp_structure(structure, cation_shell_flag=False):
    """
    Generates GULP input string for pymatgen structure 
    """
    #identify cations
    #generate two lines with core and shell if cation_shell_flag
    #For anions, generate two lines with core and shell
    gin = ""
    gin += "cell\n"
    lat = structure.lattice
    gin += str(lat.a)+" "+str(lat.b)+" "+str(lat.c)+" "
    gin += str(lat.alpha)+" "+str(lat.beta)+" "+str(lat.gamma)+"\n"
    gin += "frac\n"
    for site in structure.sites: 
        site_desc = site.specie.symbol+" core "+str(site.frac_coords)+"\n"
        specie = site.specie
        if specie in _anion_set or (specie in _cation_set and 
                                         cation_shell_flag):
            site_desc += site.specie.symbol+" shel "+str(site.frac_coords)+"\n"
        elif specie in _cation_set and not cation_shell_flag:
            pass
        else:
            raise GulpError(site.specie.symbol+
                    ": Don't know if element is cation or anion")
        gin += site_desc
    return gin

def run_gulp_input(ginfile):
    """
    Given a complete gulp.in file, run GULP
    """
    command=["gulp"]
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                     stdin=subprocess.PIPE,
                     stderr=subprocess.PIPE)
    output = p.communicate(ginfile)
    print "----output_0---------"
    print output[0]
    print "----End of output_0------"
    print "----output_1--------"
    print output[1]
    print "----End of output_1------"
    if "Error" in output[1] or "error" in output[1]:
        raise GulpError(output[1])

    # We may not need this
    if "ERROR" in output[0]:
        raise GulpExecutionError()

    # Mott_Littleton method may fail to reach proper convergence
    conv_err_string = "Conditions for a minimum have not been satisfied"
    if conv_err_string in output:
        raise GulpConvergenceError()

    gout_string = ""
    for line in output[0].split("\n"):
        gout_string = gout_string + line + "\n"
    return gout_string

def _gulp_energy(gulpout):
    energy = None
    for line in gulpout.split("\n"):
        if "Total lattice energy" in line and "eV" in line:
            energy = line.split()
    if energy:
        return energy[4]
    else:
        print gulpout
        raise GulpError("Energy not found in Gulp output")

def binaryoxide_tersoff_gulpinput(structure):
    '''
    Gets a GULP input for an oxide structure
    CURRENTLY ONLY WORKS FOR BINARY OXIDES WITH A SINGLE OXIDATION STATE
    Calculates oxidation state from the formula
    '''
    
    #gin=get_gulpinput(structure)
    gin=_gulp_structure(structure)
    gin="static noelectrostatics \n "+gin
    specs=structure.sites
    
    comp=structure.composition.get_el_amt_dict()

    ginput=""
    c=0
    for line in gin.split("\n"):
        c=c+1
        if c != 2 and c != 3 and c != 8:
            if c==1:
                ginput = ginput+line   
            elif c>=10 and c < 10+structure.num_sites:
                d=c-10
                ginput = ginput + "\n" + str(specs[d].specie) +" core "+line
            else:
                ginput = ginput+"\n"+line 
    ginput = ginput + "species \n"
    
    lastspec=[]
    metaloxi=[]
    endstring=""
    for ii in specs:
        if lastspec != ii.specie:
            lastspec=ii.specie
            
            '''ONE DAY UPDATE THIS WITH BVANALYZER'''
            '''The only challenge is, I don't know how to deal with cation-cation potential'''
            
            if str(lastspec) != "O":
                nummet=comp[str(lastspec)]
                numoxi=comp["O"]
                oxistate=numoxi*2/nummet
                if oxistate%1 != 0:
                    raise SystemError("Oxide has mixed valence on metal")
                oxidationstring=str(lastspec)+" core "+str(oxistate)
                metaloxi.append(lastspec)
            else:
                oxidationstring=str(lastspec)+" core -2"
                
            ginput=ginput+oxidationstring+ "\n"
            endstring=endstring+"qerfc \n"+str(lastspec)+" "+str(lastspec)+"  0.6000 10.0000 \n" 
    
    ginput=ginput+"# noelectrostatics \n Morse \n"
    
    for metal in metaloxi:
        metal=str(metal)+"("+str(int(oxistate))+")"
        MetOxiTers=Tersoff_pot().data[metal]
        ginput=ginput+MetOxiTers
    
    ginput=ginput+endstring
    
    return ginput


def gulpduplicatecheck(structure):
    '''
    Gets a GULP input for any structure
    Uses only an Au potential for all atoms.
    Not to actually get energies - just used to attribute some arbitrary energy to the structure.
    Identical structures will have the same 'energy', which can be used to pre-screen structurematcher for large structures. 
    
    '''
    
    #gin=get_gulpinput(structure)
    gin=_gulp_structure(structure)
    gin="single \n "+gin
    
    comp=structure.composition.get_el_amt_dict()

    ginput=""
    c=0
    for line in gin.split("\n"):
        c=c+1
        if c != 2 and c != 3 and c != 8:
            if c==1:
                ginput = ginput+line   
            elif c>=10 and c < 10+structure.num_sites:
                d=c-10
                ginput = ginput + "\n Au "+line
            else:
                ginput = ginput+"\n"+line 
    ginput = ginput + """lennard 12 6
Au core Au core 214180.2000 625.482 40.000 0 0"""
    output2=run_gulp_input(ginput)
    return float(_gulp_energy(output2))

def get_binoxi_gulp_energy(structure):
    output=binaryoxide_tersoff_gulpinput(structure)
    output2=run_gulp_input(output)
    return float(_gulp_energy(output2))


class GulpError(Exception):
    """
    Exception class for GULP.
    Raised when the GULP gives an error
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return "GulpError : " + self.msg

class GulpConvergenceError(Exception):
    """
    Exception class for GULP.
    Raised when proper convergence is not reached in Mott-Littleton
    defect energy optimisation procedure in GULP
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return "GulpConvergenceError : " + self.msg


class Tersoff_pot(object):
    def __init__(self, verbose=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fid = open(os.path.join(module_dir, "OxideTersoffPotentials"), "rU")
        data = dict()
        for row in fid:
            metaloxi=row.split()[0]
            line=row.split(")")
            data[metaloxi]=line[1]
        fid.close()
        self.data=data
        

def get_gulpinput(structure):
    """
    From a structure, get a starting gulp.in file using aconvasp
    """
    output=run_aconvasp_command(["aconvasp", "--gulp"], structure)
    gin_string = ""
    for line in output[0].split("\n"):
        if gin_string=="":
            gin_string = gin_string + line
        else:
            gin_string = gin_string +  "\n" +line
    
    return gin_string

