"""
This package implements modules for input and output to and from VASP. It
imports the key classes form both vasp_input and vasp_output to allow most
classes to be simply called as pymatgen.io.vaspio.Incar for example, to retain
backwards compatibility.
"""

from .vasp_input import *
from .vasp_output import *
