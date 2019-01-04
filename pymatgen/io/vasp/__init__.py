# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This package implements modules for input and output to and from VASP. It
imports the key classes form both vasp_input and vasp_output to allow most
classes to be simply called as pymatgen.io.vasp.Incar for example, to retain
backwards compatibility.
"""

from .inputs import *
from .outputs import *
