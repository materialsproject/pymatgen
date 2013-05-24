#!/usr/bin/env python

"""
Interface with command line GULP.
http://projects.ivec.org/gulp/help/manuals.html
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

from pymatgen.command_line.aconvasp_caller import run_aconvasp_command


def get_gulpinput(structure):
    """
    From a structure, get a starting gulp.in file using aconvasp
    """
    output = run_aconvasp_command(["aconvasp", "--gulp"], structure)
    gin_string = ""
    for line in output[0].split("\n"):
        if gin_string == "":
            gin_string = gin_string + line
        else:
            gin_string = gin_string + "\n" +line

    return gin_string


def run_gulp_input(ginfile):
    """
    Given a complete gulp.in file, run GULP
    """
    command = ["gulp"]
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate(ginfile)
    if "Error" in output[1]:
        raise GulpError(output[1])

    gout_string = ""
    for line in output[0].split("\n"):
        gout_string = gout_string + line + "\n"
    return gout_string


def _gulp_energy(gulpout):
    for line in gulpout.split("\n"):
        if "Total lattice energy" in line and "eV" in line:
            energy = line.split()
    return energy[4]


def binaryoxide_tersoff_gulpinput(structure):
    """
    Gets a GULP input for an oxide structure
    CURRENTLY ONLY WORKS FOR BINARY OXIDES WITH A SINGLE OXIDATION STATE
    Calculates oxidation state from the formula
    """

    gin = get_gulpinput(structure)
    gin = "static noelectrostatics \n " + gin
    specs = structure.sites

    comp = structure.composition.get_el_amt_dict()

    ginput = ""
    c = 0
    for line in gin.split("\n"):
        c += 1
        if c != 2 and c != 3 and c != 8:
            if c == 1:
                ginput += line
            elif 10 <= c < 10 + structure.num_sites:
                d = c - 10
                ginput += "\n" + str(specs[d].specie) + " core " + line
            else:
                ginput += "\n" + line
    ginput += "species \n"

    lastspec = []
    metaloxi = []
    endstring = ""
    for ii in specs:
        if lastspec != ii.specie:
            lastspec = ii.specie

            '''ONE DAY UPDATE THIS WITH BVANALYZER'''
            '''The only challenge is, I don't know how to deal with cation-cation potential'''

            if str(lastspec) != "O":
                nummet = comp[str(lastspec)]
                numoxi = comp["O"]
                oxistate = numoxi*2/nummet
                if oxistate % 1 != 0:
                    raise SystemError("Oxide has mixed valence on metal")
                oxidationstring = str(lastspec)+" core "+str(oxistate)
                metaloxi.append(lastspec)
            else:
                oxidationstring = str(lastspec)+" core -2"

            ginput += oxidationstring+ "\n"
            endstring += "qerfc \n" + str(lastspec) + " " + str(lastspec) + \
                         "  0.6000 10.0000 \n"

    ginput += "# noelectrostatics \n Morse \n"

    for metal in metaloxi:
        metal = str(metal) + "(" + str(int(oxistate)) + ")"
        MetOxiTers = Tersoff_pot().data[metal]
        ginput += MetOxiTers

    ginput += endstring

    return ginput


def gulpduplicatecheck(structure):

    gin = get_gulpinput(structure)
    gin = "single \n "+gin

    ginput = ""
    c = 0
    for line in gin.split("\n"):
        c += 1
        if c != 2 and c != 3 and c != 8:
            if c == 1:
                ginput += line
            elif 10 <= c < 10+structure.num_sites:
                ginput += "\n" + "Au core "+line
            else:
                ginput += "\n" + line
    ginput += "lennard 12 6 \n"
    ginput += "Au core Au core 214180.2000 625.482 40.000 0 0"

    output2 = run_gulp_input(ginput)
    return float(_gulp_energy(output2))


def get_binoxi_gulp_energy(structure):
    output = binaryoxide_tersoff_gulpinput(structure)
    output2 = run_gulp_input(output)
    return float(_gulp_energy(output2))


class GulpError(Exception):
    """
    Exception class for aconvasp.
    Raised when the aconvasp gives an error
    """
    pass


class Tersoff_pot(object):
    def __init__(self, verbose=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fid = open(os.path.join(module_dir, "OxideTersoffPotentials"), "rU")
        data = dict()
        for row in fid:
            metaloxi = row.split()[0]
            line = row.split(")")
            data[metaloxi] = line[1]
        fid.close()
        self.data = data
