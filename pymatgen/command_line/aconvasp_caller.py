# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

"""
Interface with command line aconvasp. http://aflowlib.org/
Only tested on Linux. Inspired by Shyue"s qhull_caller
WARNING: you need to have a convasp in your path for this to work
"""

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy.hautier@uclouvain.be"
__status__ = "Production"
__date__ = "$Nov 22, 2011M$"

import subprocess

from pymatgen.io.vasp.inputs import Poscar


def run_aconvasp_command(command, structure):
    """
    Helper function for calling aconvasp with different arguments
    """
    poscar = Poscar(structure)
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output = p.communicate(input=poscar.get_string())
    return output


def get_num_division_kpoints(structure, kppa):
    """
    Get kpoint divisions for a given k-point density (per reciprocal-atom):
    kppa and a given structure
    """
    output = run_aconvasp_command(["aconvasp", "--kpoints", str(kppa)],
                                  structure)
    tmp = output[0].rsplit("\n")[6].rsplit(" ")
    return [int(tmp[5]), int(tmp[6]), int(tmp[7])]


def get_minkowski_red(structure):
    """
    Get a minkowski reduced structure
    """
    output = run_aconvasp_command(["aconvasp", "--kpath"], structure)
    started = False
    poscar_string = ""
    if "ERROR" in output[1]:
        raise AconvaspError(output[1])
    for line in output[0].split("\n"):
        if started or line.find("KPOINTS TO RUN") != -1:
            poscar_string = poscar_string + line + "\n"
        if line.find("STRUCTURE TO RUN") != -1:
            started = True
        if line.find("KPOINTS TO RUN") != -1:
            started = False
    return Poscar.from_string(poscar_string).structure


def get_conv_struct(structure):
    """
    Get a minkowski reduced structure
    """
    output = run_aconvasp_command(["aconvasp", "--std_conv"], structure)
    if "ERROR" in output[1]:
        raise AconvaspError(output[1])
    tmp = Poscar.from_string(output[0])
    return {'struct': tmp.structure, 'comm': tmp.comment}


def get_prim_struct(structure):
    """
    Get standard primitive
    """
    output = run_aconvasp_command(["aconvasp", "--std_prim"], structure)
    if "ERROR" in output[1]:
        raise AconvaspError(output[1])
    tmp = Poscar.from_string(output[0])
    return {'struct': tmp.structure, 'comm': tmp.comment}


def get_vasp_kpoint_file_sym(structure):
    """
    get a kpoint file ready to be ran in VASP along the symmetry lines of the
    Brillouin Zone
    """
    output = run_aconvasp_command(["aconvasp", "--kpath"], structure)
    if "ERROR" in output[1]:
        raise AconvaspError(output[1])
    started = False
    kpoints_string = ""
    for line in output[0].split("\n"):
        #print line
        if started or line.find("END") != -1:
            kpoints_string = kpoints_string + line + "\n"
        if line.find("KPOINTS TO RUN") != -1:
            started = True
        if line.find("END") != -1:
            started = False
    return kpoints_string


class AconvaspError(Exception):
    """
    Exception class for aconvasp.
    Raised when the aconvasp gives an error
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return "AconvaspError : " + self.msg
